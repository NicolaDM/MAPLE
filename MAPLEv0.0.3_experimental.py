import sys
from math import log
import argparse
from time import time
import os.path

#©EMBL-European Bioinformatics Institute, 2021

# MAPLE code to estimate a tree by maximum likelihood from a MAPLE format input.



#things that it makes sense to try soon:

#TODO develop more nuanced exploration of topology space: reattempt to re-place only subtrees that have been affected in the last round of topology search, and their neighbours.


#things one could maybe do in future:
#TODO given an efficient way to optimize placement branch length, for example with an analytical derivative solution, we could directly find the likelihood 
# for the best position within a branch without having to re-traverse these branches and trying explicitly different branch lengths?
#TODO would it be convenient to remove the genome position element from entries of type ACGTO ? Sounds like a lot of changes, and might make things slower and less readable?
#TODO if the model is reversible, no need to consider complicated case of root position, so the genome vector lists can be simpler and I can avoid additional calculations?
#TODO does it make sense to have createFurtherMidNodes() or maybe rather save on the memory and create these genome lists anew each time? Maybe rather find optimal values for the threshold.
#TODO re-write updateBLen(), updatePartialsFromTop(), and updatePartialsFromBottom() as non-recursive function(s)?
#TODO try to replace appendProb() with appendProbNode()? it might be slightly slower, but it would be one fewer function in the code.



parser = argparse.ArgumentParser(description='Estimate a tree from a diff format and using iterative approximate maximum likelihood sample placement.')
parser.add_argument('--input',default="/Users/demaio/Desktop/GISAID-hCoV-19-phylogeny-2021-03-12/phylogenetic_inference/2021-03-31_unmasked_differences_reduced.txt_consensus-based.txt", help='input diff file name.')
parser.add_argument('--reference',default="/Users/demaio/Desktop/GISAID-hCoV-19-phylogeny-2021-03-12/phylogenetic_inference/2021-03-31_unmasked_differences_reduced.txt_consensus.fa", help='input reference file name found in folder specified by --path.')
parser.add_argument('--output',default="/Users/demaio/Desktop/GISAID-hCoV-19-phylogeny-2021-03-12/phylogenetic_inference/MAPLE", help='output path and identifier to be used for newick output file.')
parser.add_argument("--onlyNambiguities", help="Treat all ambiguities as N (total missing information).", action="store_true")
parser.add_argument("--thresholdProb",help="relative probability threshold used to ignore possible states with very low probabilities.",  type=float, default=0.0000001)
parser.add_argument("--thresholdLogLK",help="logLK difference threshold to consider a logLk close to optimal.",  type=float, default=60.0)
parser.add_argument("--thresholdLogLKtopology",help="logLK difference threshold to consider a logLk close to optimal when looking for topology improvements.",  type=float, default=60.0)
parser.add_argument("--allowedFails",help="Number of times one can go down the tree without inclreasing placement likelihood before the tree traversal is stopped (only applies to non-0 branch lengths).",  type=int, default=3)
parser.add_argument("--allowedFailsTopology",help="Number of times one can crawl along the tree decreasing placement likelihood before the tree traversal is stopped during topology search (only applies to non-0 branch lengths).",  type=int, default=3)
parser.add_argument("--bLenAdjustment",help="If >1, try placing also with a longer bLen than the standard bLen.",  type=int, default=1)
parser.add_argument("--verbose", help="Print to screen a lot of stuff.", action="store_true")
parser.add_argument("--model", help="Which substitution model should be used. Allowed models so far are JC, GTR (default) or UNREST.", default="GTR")
parser.add_argument("--overwrite", help="Overwrite previous results if already present.", action="store_true")
parser.add_argument("--binaryTree", help="Write output tree as binary.", action="store_true")
parser.add_argument("--numTopologyImprovements",help="Number of times we traverse the tree looking for topological improvements.",  type=int, default=3)
parser.add_argument("--thresholdTopologyPlacement",help="Don't try to re-place nodes that have current appending logLK cost above this threshold.",  type=float, default=-0.2)
parser.add_argument("--minBLenForMidNode",help="Don't try to place samples mid-branch if the considered branch is below this length in number of mutations genome-wide.",  type=float, default=4.1)
parser.add_argument("--updateSubstMatrixEveryThisSamples",help="How many new samples to place before each update of the substitution rate matrix.",  type=int, default=100)
parser.add_argument("--nonStrictInitialStopRules", help="If specified, then during the initial placement stage, a slower, non-strict rule for stopping the placement search is applied: the search is stopped if enough many consencutive LK worsening are observed, AND if LK is below the considered threshold.", action="store_true")
parser.add_argument("--nonStrictTopologyStopRules", help="If specified, then during the topological improvement stage, a slower, non-strict rule for stopping the SPR search is applied: the search is stopped if enough many consencutive LK worsening are observed, AND if LK is below the considered threshold.", action="store_true")
parser.add_argument("--thresholdDiffForUpdate",help="Consider the probability of a new partial changed if the difference between old and new if above this threshold.",  type=float, default=0.01)
parser.add_argument("--thresholdFoldChangeUpdate",help="Consider the probability of a new partial changed, if the fold difference between old and new if above this threshold.",  type=float, default=1.5)
parser.add_argument("--minBLen",help="Don't attempt branch lengths below this number of mutations (genome-wide): just consider 0 length instead.",  type=float, default=0.2)
parser.add_argument("--maxBLen",help="Don't attempt branch lengths above this number of mutations (genome-wide).",  type=float, default=40.0)
parser.add_argument("--thresholdLogLKconsecutivePlacement",help="logLK difference threshold to consider something as a significant decrease in log-LK when considering consecutive likelihood decreases.",  type=float, default=1.0)
args = parser.parse_args()

onlyNambiguities=args.onlyNambiguities
thresholdProb=args.thresholdProb
verbose=args.verbose
inputFile=args.input
outputFile=args.output
refFile=args.reference
allowedFails=args.allowedFails
allowedFailsTopology=args.allowedFailsTopology
model=args.model
bLenAdjustment=args.bLenAdjustment
if bLenAdjustment>1.0+thresholdProb:
	tryOtherBLen=True
else:
	tryOtherBLen=False
thresholdLogLK=args.thresholdLogLK
thresholdLogLKtopology=args.thresholdLogLKtopology
overwrite=args.overwrite
binaryTree=args.binaryTree
numTopologyImprovements=args.numTopologyImprovements
thresholdTopologyPlacement=args.thresholdTopologyPlacement
minBLenForMidNode=args.minBLenForMidNode
updateSubstMatrixEveryThisSamples=args.updateSubstMatrixEveryThisSamples
strictInitialStopRules=(not args.nonStrictInitialStopRules)
strictTopologyStopRules=(not args.nonStrictTopologyStopRules)
thresholdDiffForUpdate=args.thresholdDiffForUpdate
thresholdFoldChangeUpdate=args.thresholdFoldChangeUpdate
minBLen=args.minBLen
maxBLen=args.maxBLen
thresholdLogLKconsecutivePlacement=args.thresholdLogLKconsecutivePlacement
example=False

minimumCarryOver=sys.float_info.min*(1e50)

if os.path.isfile(outputFile+"_tree.tree")  and (not overwrite):
	print("File "+outputFile+"_tree.tree already exists, quitting fastLK tree inference. Use option --overwrite if you want to overwirte previous inference.")
	exit()
if not os.path.isfile(inputFile):
	print("Input file in Maple format "+inputFile+" not found, quitting MAPLE tree inference. Use option --input to specify a valid input file.")
	exit()
if not os.path.isfile(refFile):
	print("Input reference fasta file "+refFile+" not found, quitting MAPLE tree inference. Use option --reference to specify a valid input reference file.")
	exit()
lineSplit=outputFile.split("/")
lineSplit[-1]=""
outFolder="/".join(lineSplit)
if not os.path.isdir(outFolder):
	print("Path to output file "+outFolder+" does not exist, quitting MAPLE tree inference. Use option --output to specify a valid output file path and output file name.")
	exit()

#Class defininf nodes of the tree
class Tree(object):
	def __init__(self, name='', children=None, dist=1.0):
		if name!='':
			self.name = name
		self.dist = dist
		self.children = []
		self.up=None
		if children is not None:
			for child in children:
				self.add_child(child)
	def __repr__(self):
		try:
			return self.name
		except AttributeError:
			return ""
	def add_child(self, node):
		assert isinstance(node, Tree)
		self.children.append(node)


alleles={"A":0,"C":1,"G":2,"T":3}
allelesList=["A","C","G","T"]
allelesLow={"a":0,"c":1,"g":2,"t":3}
allelesUpOrLow={"a":0,"c":1,"g":2,"t":3,"A":0,"C":1,"G":2,"T":3}
allelesListLow=["a","c","g","t"]
ambiguities={"y":[0.0,0.5,0.0,0.5],"r":[0.5,0.0,0.5,0.0],"w":[0.5,0.0,0.0,0.5],"s":[0.0,0.5,0.5,0.0],"k":[0.0,0.0,0.5,0.5],"m":[0.5,0.5,0.0,0.0],"d":[1.0/3,0.0,1.0/3,1.0/3],"v":[1.0/3,1.0/3,1.0/3,0.0],"h":[1.0/3,1.0/3,0.0,1.0/3],"b":[0.0,1.0/3,1.0/3,1.0/3]}

#collect reference and calculate vector of cumulative numbers of nucleotides and calculate root frequencies.
def collectReference(fileName):
	file=open(fileName)
	line=file.readline()
	ref=""
	while line!="":
		line=file.readline()
		ref+=line.replace("\n","")
	ref=ref.lower()
	lRef=len(ref)
	#print("Ref genome length: "+str(lRef))
	file.close()
	#vector to count how many bases of each type are cumulatively in the reference genome up to a certain position
	cumulativeBases=[[0,0,0,0]]
	for i in range(lRef):
		cumulativeBases.append(list(cumulativeBases[i]))
		cumulativeBases[i+1][allelesUpOrLow[ref[i]]]+=1
	rootFreqs=[0.0,0.0,0.0,0.0]
	rootFreqsLog=[0.0,0.0,0.0,0.0]
	for i in range(4):
		rootFreqs[i]=cumulativeBases[-1][i]/float(lRef)
		rootFreqsLog[i]=log(rootFreqs[i])
	return ref, cumulativeBases, rootFreqs, rootFreqsLog

ref, cumulativeBases, rootFreqs, rootFreqsLog = collectReference(refFile)
lRef=len(ref)
print("Length of reference genome: "+str(lRef))
refIndeces=[]
for i in range(lRef):
	refIndeces.append(allelesLow[ref[i]])
if model=="JC":
	rootFreqs=[0.25,0.25,0.25,0.25]
	rootFreqsLog=[log(0.25),log(0.25),log(0.25),log(0.25)]
oneMutBLen=1.0/lRef
minBLen=minBLen*oneMutBLen
maxBLen=maxBLen*oneMutBLen
minBLenForMidNode=minBLenForMidNode*oneMutBLen



#Read input file
def readConciseAlignment(fileName):
	fileI=open(fileName)
	line=fileI.readline()
	nSeqs=0
	data=[]
	while line!="" and line!="\n":
		seqList=[]
		#name=line.replace(">","").replace("\n","")
		line=fileI.readline()
		while line!="" and line!="\n" and line[0]!=">":
			linelist=line.split()
			if len(linelist)>2:
				entry=(linelist[0].lower(),int(linelist[1]),int(linelist[2]))
			else:
				entry=(linelist[0].lower(),int(linelist[1]))
			if ref[entry[1]-1]==entry[0]:
				print("Mutation observed into reference nucleotide at position "+str(entry[1])+" , nucleotide "+entry[0]+". Wrong reference and/or diff file?")
				exit()
			seqList.append(entry)
			line=fileI.readline()
		#data[name]=seqList
		data.append(seqList)
		nSeqs+=1
	fileI.close()
	print(str(nSeqs)+" sequences in diff file.")
	return data

runOnlyExample=False
if not runOnlyExample:	
	#read sequence data from file
	data=readConciseAlignment(inputFile)

range4=range(4)
#stricter thresholds than the input one for probabilities.
thresholdProb2=thresholdProb*thresholdProb
thresholdProb4=thresholdProb2*thresholdProb2

#IMPORTANT definitions of genome list entries structure.
#The first element of a genome list entry represnts its type: 0="A", 1="C", 2="G", 3="T", 4="R", 5="N", 6="O"
#the second element represents the last position of the stretch of genome that they represent; 
# a value of 1 refers to the first position of the genome.
#For entries of type "O"/6 the last element is always a normalized vector of likelihoods (they always have to sum to 1.0).
#If present, anoter element after the position one represents the evolutionary distance since the considered type was observed (this is always missing with type "N"/5);
#If the element is not not present, this distance is assumed to be 0.0.
#If the distance (branch length) value is present, another distance value can also be present for entries of type <5 ; 
#this is to account for the fact that the observation might have occurred on the other side of the phylogeny with respect to the root; 
#when this additional branch length is also present, for example if the entry is (1,234,0.0001,0.0002) then this means that a "C" was observed at genome position 234, and the observation
#is separated from the root by a distance of 0.0001, while a distance of 0.0002 separates the root from the current node (or position along a branch) considered.

#if probability mass is concentrated in one nucleotide, simplify the entry from "O" to another type.
def simplfy(vec,refA):
	maxP=0.0
	maxI=0
	numA=0
	for i in range4:
		if vec[i]>maxP:
			maxP=vec[i]
			maxI=i
		if vec[i]>thresholdProb:
			numA+=1
	if maxP<thresholdProb4:
		print("Inside simplify(), all values in vector are too small - something wrong?")
		print(vec)
		exit()
	if numA==1:
		if maxI==refA:
			return 4
		else:
			return maxI
	else:
		return 6

#Shorten genome list by merging together R entries that are mergeable
def shorten(vec):
	entryOld=vec[0]
	index=0
	while index<len(vec)-1:
		newVec=vec[index+1]
		if newVec[0]==4 and entryOld[0]==4 and len(newVec)==len(entryOld):
			if len(newVec)==2:
				vec.pop(index)
			elif abs(newVec[2]-entryOld[2])>thresholdProb:
				index+=1
				entryOld=vec[index]
			elif len(newVec)==3:
				vec.pop(index)
			elif abs(newVec[3]-entryOld[3])<thresholdProb:
				vec.pop(index)
			else:
				index+=1
				entryOld=vec[index]
		else:
			index+=1
			entryOld=vec[index]
	return

#Like shorten(), but specific to lower likelihoods that don't need to consider root frequencies. 
#UNNECESSARY, one could always use shorten() instead of it.
def shortenLower(vec):
	entryOld=vec[0]
	index=0
	while index<len(vec)-1:
		newVec=vec[index+1]
		if newVec[0]==4 and entryOld[0]==4 and len(newVec)==len(entryOld):
			if len(newVec)==2:
				vec.pop(index)
			elif abs(newVec[2]-entryOld[2])>thresholdProb:
				index+=1
				entryOld=vec[index]
			else:
				vec.pop(index)
		else:
			index+=1
			entryOld=vec[index]
	return



#define partial likelihood vector for a sample given its data
def probVectTerminalNode(diffs):
	if diffs is None:
		probVect=[(5,lRef)]
		return probVect
	pos=1
	probVect=[]
	for m in diffs:
			currPos=m[1]
			if currPos>pos:
				#region where the node with branch length bLen is identical to the ref.
				probVect.append((4,currPos-1))
				pos=currPos
			if m[0]=="n" or m[0]=="-":
				if len(m)>2:
					length=m[2]
				else:
					length=1
				#region with no info, store last position and length.
				probVect.append((5,currPos+length-1))
				pos=currPos+length
			elif m[0] in allelesLow:
				#position at which node allele is sure but is different from the reference.
				probVect.append((allelesLow[m[0]],currPos))
				pos=currPos+1
			else:
				# non-"n" ambiguity character; for now interpret this as ambiguity instead of as a polymorphism.
				if onlyNambiguities:
					# if user asks to, to make things easier, interpret any ambiguity as an "n".
					probVect.append((5,currPos))
				else:
					#otherwise, store as an "other" scenario, where each nucleotide has its own partial likelihood.
					probVect.append((6,currPos,ambiguities[m[0]]))
				pos=currPos+1
	if pos<=lRef:
		probVect.append((4,lRef))
	return probVect




#Update and normalize the mutation rate matrix, given new mutation counts
def updateSubMatrix(pseudoMutCounts,model,oldMutMatrix): #,mutMatrixBLen,mutMatrixBLenAdjusted
	mutMatrix=[[0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0]]
	if model=="UNREST":
		for i in range4:
			tot=0.0
			for j in range4:
				if j!=i:
					mutMatrix[i][j]=pseudoMutCounts[i][j]/rootFreqs[i]
					tot+=mutMatrix[i][j]
			mutMatrix[i][i]=-tot
	elif model=="GTR":
		for i in range4:
			tot=0.0
			for j in range4:
				if j!=i:
					mutMatrix[i][j]=(pseudoMutCounts[i][j]+pseudoMutCounts[j][i])/rootFreqs[i]
					tot+=mutMatrix[i][j]
			mutMatrix[i][i]=-tot
	else:
		print("Substitution model not recognised! Exiting")
		exit()
	totRate=-(rootFreqs[0]*mutMatrix[0][0]+ rootFreqs[1]*mutMatrix[1][1]+ rootFreqs[2]*mutMatrix[2][2]+ rootFreqs[3]*mutMatrix[3][3] )
	for i in range4:
		for j in range4:
			mutMatrix[i][j]=mutMatrix[i][j]/totRate
	matChange=0.0
	for i in range4:
		for j in range4:
			if j!=i:
				matChange+=abs(mutMatrix[i][j]-oldMutMatrix[i][j])
	if matChange > 0.001: #the matrix has changed significantly, update it
		for i in range4:
			for j in range4:
				oldMutMatrix[i][j]=mutMatrix[i][j]
		#the matrix has changed significantly, return True so that cumulative rates can be updated as well.
		return True
	else:
		return False

#Initialize the mutation rate matrix
#If a fixed rate matrix is needed for SARS-CoV-2, this is the nucleotide mutation rate matrix from De Maio et al 2021
#mutMatrix=[[0.0,0.039,0.310,0.123],[0.140,0.0,0.022,3.028],[0.747,0.113,0.0,2.953],[0.056,0.261,0.036,0.0]]
if model=="JC":
	mutMatrix=[[-1.0,1.0/3,1.0/3,1.0/3],[1.0/3,-1.0,1.0/3,1.0/3],[1.0/3,1.0/3,-1.0,1.0/3],[1.0/3,1.0/3,1.0/3,-1.0]]
else:
	mutMatrix=[[0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0]]
	pseudoMutCounts=[[0.0,1.0,5.0,2.0],[2.0,0.0,1.0,40.0],[5.0,2.0,0.0,20.0],[2.0,3.0,1.0,0.0]]
	updateSubMatrix(pseudoMutCounts,model,mutMatrix)

nonMutRates=[0,0,0,0]
for i in range4:
	nonMutRates[i]=mutMatrix[i][i]

cumulativeRate=[0.0]
for i in range(lRef):
	ind=refIndeces[i]
	cumulativeRate.append(cumulativeRate[-1]+nonMutRates[ind])





#Sort samples based on distance from reference, but punishing more isolated N's and ambiguity characters.
#more ambiguous sequences are placed last this way - this is useful since ambiguous sequences are harder to place and are more likely to be less informative (and so be removed from the analysis altogether)
#def distancesFromRefPunishNs(data,samples):
def distancesFromRefPunishNs(data):
	sampleDistances=[]
	#for sample in samples:
	for diffIndex in range(len(data)):
		diffs=data[diffIndex]
		pos=1
		comparisons=0
		diffNum=0
		for m in diffs:
			currPos=m[1]
			if currPos>pos:
				#region identical to the ref.
				comparisons+=currPos-pos
				pos=currPos
			if m[0]=="n" or m[0]=="-":
				if len(m)>2:
					pos=currPos+m[2]
				else:
					pos=currPos+1
				diffNum+=1
			elif m[0] in allelesLow:
				comparisons+=1
				diffNum+=1
				pos=currPos+1
			else:
				pos=currPos+1
				diffNum+=1
		if pos<=lRef:
			comparisons+=lRef+1-pos
		#sampleDistances.append((diffNum*1000+lRef-comparisons,sample))
		sampleDistances.append((diffNum*1000+lRef-comparisons,diffIndex))

	from operator import itemgetter
	print("Now doing sorting")
	sampleDistances.sort(reverse=True,key=itemgetter(0))
	return sampleDistances



#function to check is one sequence is less informative than another;
#returns 0 when the 2 sequences are not comparable, otherwise returns 1 if the first is more informative or if they are identical, and 2 otherwise.
def isMinorSequence(probVect1,probVect2):
	indexEntry1, indexEntry2, pos = 0, 0, 0
	entry1=probVect1[indexEntry1]
	entry2=probVect2[indexEntry2]
	found1bigger=False
	found2bigger=False
	while True:
		if entry1[0]!=entry2[0]:
			if entry1[0]==5:
				pos=min(entry1[1],entry2[1])
				found2bigger=True
			elif entry2[0]==5:
				pos=min(entry1[1],entry2[1])
				found1bigger=True
			elif entry1[0]==6:
				if entry2[0]==4:
					i2=refIndeces[pos]
				else:
					i2=entry2[0]
				if entry1[-1][i2]>0.1:
					found2bigger=True
				else:
					return 0
				pos+=1
			elif entry2[0]==6:
				if entry1[0]==4:
					i1=refIndeces[pos]
				else:
					i1=entry1[0]
				if entry2[-1][i1]>0.1:
					found1bigger=True
				else:
					return 0
				pos+=1
			else:
				return 0
		elif entry1[0]==6:
			for j in range4:
				if entry2[-1][j]>0.1 and entry1[-1][j]<0.1:
					found1bigger=True
				elif entry1[-1][j]>0.1 and entry2[-1][j]<0.1:
					found2bigger=True
			pos+=1
		else:
			pos=min(entry1[1],entry2[1])
		if found1bigger and found2bigger:
			return 0
		if pos==lRef:
			break
		if pos==entry1[1]:
			indexEntry1+=1
			entry1=probVect1[indexEntry1]
		if pos==entry2[1]:
			indexEntry2+=1
			entry2=probVect2[indexEntry2]

	if found1bigger:
		if found2bigger:
			return 0
		else: 
			return 1
	else:
		if found2bigger:
			return 2
		else:
			return 1






#function to calculate likelihood cost of appending sample to parent node
#UNNECESSARY this can be considered an optimized version of function appendProbNode().
#as such, when rewriting the code, one may get rid of this function altogether and use just appendProbNode() instead.
def appendProb(probVectP,probVectC,bLen,mutMatrix):
	if not bLen:
		bLen=0.0
	Lkcost, indexEntry1, indexEntry2, totalFactor, pos = 0.0, 0, 0, 1.0, 0
	entry1=probVectP[indexEntry1]
	entry2=probVectC[indexEntry2]
	end=min(entry1[1],entry2[1])
	contribLength=bLen
	while True:
		if entry2[0]==5: # case N
			pos=min(entry1[1],entry2[1])
		elif entry1[0]==5: # case N
			#if parent node is type "N", in theory we might have to calculate the contribution of root nucleotides; 
			# however, if this node is "N" then every other node in the current tree is "N", so we can ignore this since this contribution cancels out in relative terms.
			pos=min(entry1[1],entry2[1])
		elif entry1[0]==4: # case entry1 is R
			if entry2[0]==4:
				end=min(entry1[1],entry2[1])
				if len(entry1)==2:
					Lkcost+=bLen*(cumulativeRate[end]-cumulativeRate[pos])
				else:
					contribLength=bLen+entry1[2]
					if len(entry1)==3:
						Lkcost+=contribLength*(cumulativeRate[end]-cumulativeRate[pos])
					else:
						#here contribution from root frequency gets added and subtracted so it's ignored
						Lkcost+=(contribLength+entry1[3])*(cumulativeRate[end]-cumulativeRate[pos])
				pos=end
			elif entry2[0]==6:
				i1=refIndeces[pos]
				if len(entry1)==4:
					contribLength=bLen+entry1[3]
					if entry2[2][i1]>0.1:
						contribLength+=entry1[2]
						#here contribution from root frequency can also be also ignored
						Lkcost+=nonMutRates[i1]*contribLength 
					else:
						tot=0.0
						for i in range4:
							if i1==i:
								tot2=rootFreqs[i]*(1.0+nonMutRates[i]*entry1[2])
							else:
								tot2=rootFreqs[i]*mutMatrix[i][i1]*entry1[2]
							tot3=0.0
							for j in range4:
								if entry2[2][j]>0.1:
									tot3+=mutMatrix[i][j]
							tot3*=contribLength
							if entry2[2][i]>0.1:
								tot3+=1.0
							tot+=tot2*tot3
						totalFactor*=(tot/rootFreqs[i1])
				else:
					if entry2[2][i1]>0.1:
						if len(entry1)==3:
							Lkcost+=nonMutRates[i1]*(bLen+entry1[2])
						else:
							Lkcost+=nonMutRates[i1]*bLen
					else:
						tot=0.0
						for j in range4:
							if entry2[2][j]>0.1:
								tot+=mutMatrix[i1][j]
						if len(entry1)==3:
							totalFactor*=tot*(bLen+entry1[2])
						else:
							totalFactor*=tot*bLen
				pos+=1

			else: #entry1 is reference and entry2 is a different but single nucleotide
				if len(entry1)==2:
					totalFactor*=mutMatrix[refIndeces[pos]][entry2[0]]*bLen
				elif len(entry1)==3:
					totalFactor*=mutMatrix[refIndeces[pos]][entry2[0]]*(bLen+entry1[2])
				else:
					i1=refIndeces[pos]
					i2=entry2[0]
					totalFactor*=((rootFreqs[i1]*mutMatrix[i1][i2]*(bLen+entry1[3])*(1.0+nonMutRates[i1]*entry1[2])+rootFreqs[i2]*mutMatrix[i2][i1]*entry1[2]*(1.0+nonMutRates[i2]*(bLen+entry1[3])))/rootFreqs[i1])
				pos+=1

		# entry1 is of type "O"
		elif entry1[0]==6:
			if len(entry1)==3:
				bLen13=bLen
			else:
				bLen13=bLen+entry1[2]
			if entry2[0]==6:
				tot=0.0
				for j in range4:
					tot2=0.0
					for j2 in range4:
						if entry2[2][j2]>0.1:
							tot2+=mutMatrix[j][j2]
					tot2*=bLen13
					if entry2[2][j]>0.1:
						tot2+=1.0
					tot+=tot2*entry1[-1][j]
				totalFactor*=tot
			else:
				if entry2[0]==4:
					i2=refIndeces[pos]
				else:
					i2=entry2[0]
				totalFactor*=(entry1[-1][i2]+bLen13*(entry1[-1][0]*mutMatrix[0][i2]+entry1[-1][1]*mutMatrix[1][i2]+entry1[-1][2]*mutMatrix[2][i2]+entry1[-1][3]*mutMatrix[3][i2]))
			pos+=1

		else: #entry1 is a non-ref nuc
			i1=entry1[0]
			if entry2[0]==i1:
				if len(entry1)==2:
					Lkcost+=nonMutRates[i1]*bLen
				elif len(entry1)==3:
					Lkcost+=nonMutRates[i1]*(bLen+entry1[2])
				else:
					Lkcost+=nonMutRates[i1]*(bLen+entry1[2]+entry1[3])
			else: #entry1 and entry2 are of different types
				if entry2[0]==6:
					if len(entry1)==4:
						bLen15=bLen+entry1[3]
						if entry2[2][i1]>0.1:
							Lkcost+=nonMutRates[i1]*(bLen15+entry1[2])
						else:
							tot=0.0
							for i in range4:
								if i1==i:
									tot2=rootFreqs[i]*(1.0+nonMutRates[i1]*entry1[2])
								else:
									tot2=rootFreqs[i]*mutMatrix[i][i1]*entry1[2]
								tot3=0.0
								for j in range4:
									if entry2[2][j]>0.1:
										tot3+=mutMatrix[i][j]
								if entry2[2][i]>0.1:
									tot+=tot2*(1.0+bLen15*tot3)
								else:
									tot+=tot2*bLen15*tot3
							totalFactor*=(tot/rootFreqs[i1])
					else:
						if entry2[2][i1]>0.1:
							if len(entry1)==2:
								Lkcost+=nonMutRates[i1]*bLen
							else:
								Lkcost+=nonMutRates[i1]*(bLen+entry1[2])
						else:
							tot=0.0
							for j in range4:
								if entry2[2][j]>0.1:
									tot+=mutMatrix[i1][j]
							if len(entry1)==2:
								totalFactor*=tot*bLen
							else:
								totalFactor*=tot*(bLen+entry1[2])
				#entry2 is a nucleotide type (like entry1)		
				else:
					if entry2[0]==4:
						i2=refIndeces[pos]
					else:
						i2=entry2[0]
					if len(entry1)==2:
						totalFactor*=mutMatrix[i1][i2]*bLen
					elif len(entry1)==3:
						totalFactor*=mutMatrix[i1][i2]*(bLen+entry1[2])
					else:
						#here we ignore contribution of non-parsimonious mutational histories
						totalFactor*=((rootFreqs[i1]*mutMatrix[i1][i2]*(bLen+entry1[3])*(1.0+nonMutRates[i1]*entry1[2]) + rootFreqs[i2]*mutMatrix[i2][i1]*entry1[2]*(1.0+nonMutRates[i2]*(bLen+entry1[3])) )/rootFreqs[i1])						
			pos+=1

		if totalFactor<=minimumCarryOver:
			if totalFactor<sys.float_info.min:
				return float("-inf")
			Lkcost+=log(totalFactor)
			totalFactor=1.0

		if pos==lRef:
			break
		if pos==entry1[1]:
			indexEntry1+=1
			entry1=probVectP[indexEntry1]
		if pos==entry2[1]:
			indexEntry2+=1
			entry2=probVectC[indexEntry2]
	return Lkcost+log(totalFactor)





#number of samples that could have been placed as major of another sample but weren't due to sample placement order
totalMissedMinors=[0]

#function to find the best node in the tree where to append the new sample; traverses the tree and tries to append the sample at each node and mid-branch nodes.
def findBestParent(t1,diffs,sample,mutMatrix):
	bestNodeSoFar=t1
	bestLKdiff=float("-inf")
	bestIsMidNode=False
	bestUpLK=float("-inf")
	bestDownLK=float("-inf")
	bestDownNode=None
	nodesToVisit=[(t1,float("-inf"),0)]
	adjustBLen=False
	while nodesToVisit:
		t1,parentLK,failedPasses=nodesToVisit.pop()
		if not t1.children: #check if the new leaf is strictly less informative than already placed leaf
			comparison=isMinorSequence(t1.probVect,diffs)
			if comparison==1:
				t1.minorSequences.append(sample)
				return t1, 1.0, False, float("-inf"), float("-inf"), None, False
			elif comparison==2:
				totalMissedMinors[0]+=1
	
		adjusted=False
		if t1.dist and t1.up!=None: # try first placing as a descendant of the mid-branch point of the branch above the current node.
			LKdiff2=appendProb(t1.probVectTotUp,diffs,oneMutBLen,mutMatrix)
			if tryOtherBLen: #try also placing with a longer new terminal branch, which can be useful if the sample has many new mutations.
				# this may become totally unnecessary if one will develop a derivative-based analytical branch length optimization.
				newLKdiff2=appendProb(t1.probVectTotUp,diffs,oneMutBLen*bLenAdjustment,mutMatrix)
				if newLKdiff2>LKdiff2:
					LKdiff2=newLKdiff2
					adjusted=True
				else:
					adjusted=False
			if LKdiff2>bestLKdiff:
				adjustBLen=adjusted
				bestLKdiff=LKdiff2
				bestNodeSoFar=t1
				failedPasses=0
				bestIsMidNode=True
		else:
			LKdiff2=float("-inf")

		if t1.dist: #now, try to place as descendant of the current node (this is skipped if the node has top branch length 0 and so is part of a polytomy).
			LKdiff=appendProb(t1.probVectTot,diffs,oneMutBLen,mutMatrix)
			if tryOtherBLen:
				newLKdiff=appendProb(t1.probVectTot,diffs,oneMutBLen*bLenAdjustment,mutMatrix)
				if newLKdiff>LKdiff:
					LKdiff=newLKdiff
					adjusted=True
				else:
					adjusted=False
			if LKdiff>bestLKdiff:
				adjustBLen=adjusted
				bestLKdiff=LKdiff
				bestNodeSoFar=t1
				failedPasses=0
				bestIsMidNode=False
				bestUpLK=LKdiff2
			elif LKdiff2>=(bestLKdiff-thresholdProb):
				bestUpLK=parentLK
				bestDownLK=LKdiff
				bestDownNode=t1
			elif LKdiff<(parentLK-thresholdLogLKconsecutivePlacement): #placement at current node is considered failed if placement likelihood is not improved by a certain margin compared to best placement so far for the nodes above it.
				failedPasses+=1
		else:
			LKdiff=parentLK

		#keep trying to place at children nodes, unless placement has failed too many times already or the likelihood cost suboptimality is above a certain threshold
		if strictInitialStopRules:
			if failedPasses<=allowedFails and LKdiff>(bestLKdiff-thresholdLogLK): 
				for c in t1.children:
					nodesToVisit.append((c,LKdiff,failedPasses))
		else:
			if failedPasses<=allowedFails or LKdiff>(bestLKdiff-thresholdLogLK): 
				for c in t1.children:
					nodesToVisit.append((c,LKdiff,failedPasses))

	#exploration of the tree is finished, and we are left with the node found so far with the best appending likelihood cost.
	#Now we explore placement just below this node for more fine-grained placement within its descendant branches.
	bestOfChildLK=float("-inf")
	bestChild=None
	if bestIsMidNode:
		return bestNodeSoFar , bestLKdiff , bestIsMidNode, bestUpLK, bestDownLK, bestDownNode, adjustBLen
	else:
		#current node might be part of a polytomy (represented by 0 branch lengths) so we want to explore all the children of the current node
		#to find out if the best placement is actually in any of the branches below the current node.
		nodesToVisit=[]
		for c in bestNodeSoFar.children:
			nodesToVisit.append(c)
		while nodesToVisit:
			t1=nodesToVisit.pop()
			if not t1.dist:
				for c in t1.children:
					nodesToVisit.append(c)
			else:
				#now try to place on the current branch below the best node, at an height above the mid-branch.
				newBLen2=t1.dist/2
				bestLKdiff2=float("-inf")
				furtherNode=-1
				newProbVect2=t1.probVectTotUp
				while True:
					newLKdiff2=appendProb(newProbVect2,diffs,oneMutBLen,mutMatrix)
					if tryOtherBLen:
						newLKdiff3=appendProb(newProbVect2,diffs,oneMutBLen*bLenAdjustment,mutMatrix)
						if newLKdiff3>newLKdiff2:
							newLKdiff2=newLKdiff3
					if newLKdiff2>bestLKdiff2:
						bestLKdiff2=newLKdiff2
					else:
						break
					newBLen2=newBLen2/2
					if newBLen2<=(minBLenForMidNode/2):
						break
					furtherNode+=1
					newProbVect2=t1.furtherMidNodes[furtherNode]
				
				if bestLKdiff2>bestOfChildLK:
					bestOfChildLK=bestLKdiff2
					bestChild=t1
		#pass on the best child found to the next function, which will place the new sample somewhere between the best node and the best child.
		return bestNodeSoFar , bestLKdiff , bestIsMidNode, bestUpLK, bestOfChildLK, bestChild, adjustBLen








#for the root, take lower likelihood genome list probVect, and create an overall likelihood (or upper right or upper left) genome list by multiplying likelihoods by root frequencies.
def rootVector(probVect,bLen,mutMatrix):
	newProbVect=[]
	for entry in probVect:
		if entry[0]==5:
			newProbVect.append(entry)
		elif entry[0]==6:
			if len(entry)==4:
				totBLen=entry[2]
				if bLen:
					totBLen+=bLen
			else:
				totBLen=bLen
			newVect=[]
			if totBLen:
				for i in range4:
					tot=0.0
					for j in range4:
						tot+=mutMatrix[i][j]*entry[-1][j]
					tot*=totBLen
					tot+=entry[-1][i]
					newVect.append(tot*rootFreqs[i])
				totSum=sum(newVect)
				for i in range4:
					newVect[i]/=totSum
				newProbVect.append((6,entry[1],newVect))
			else:
				for i in range4:
					newVect.append(entry[-1][i]*rootFreqs[i])
				totSum=sum(newVect)
				for i in range4:
					newVect[i]/=totSum
				newProbVect.append((6,entry[1],newVect))
		else:
			if len(entry)==3:
				if bLen:
					newProbVect.append((entry[0],entry[1],entry[2]+bLen,0.0))
				else:
					newProbVect.append((entry[0],entry[1],entry[2],0.0))
			else:
				if bLen:
					newProbVect.append((entry[0],entry[1],bLen,0.0))
				else:
					newProbVect.append((entry[0],entry[1]))
	return newProbVect







#merge two partial likelihood vectors, one from above, probVect1, and one from below, probVect2
#unlike appendProb(), this function is not used on a large part of the tree at each placement, but only in a small neighbourhood;
def mergeVectorsUpDown(probVect1,bLenUp,probVect2,bLenDown,mutMatrix):
	indexEntry1, indexEntry2, pos, end = 0, 0, 0, 0
	probVect=[]
	entry1=probVect1[indexEntry1]
	entry2=probVect2[indexEntry2]
	while True:
		if entry1[0]==5:
			if entry2[0]==5:
				pos=min(entry1[1],entry2[1])
				probVect.append((5,pos))
				
			elif entry2[0]<5:
				pos=min(entry1[1],entry2[1])
				if len(entry2)==3:
					if bLenDown:
						probVect.append((entry2[0],pos,entry2[2]+bLenDown,0.0))
					else:
						probVect.append((entry2[0],pos,entry2[2],0.0))
				else:
					if bLenDown:
						probVect.append((entry2[0],pos,bLenDown,0.0))
					else:
						probVect.append((entry2[0],pos))
			else: # entry2 case "O", entry 1 is "N"
				pos+=1
				if len(entry2)==4:
					totBLen=entry2[2]
					if bLenDown:
						totBLen+=bLenDown
				else:
					totBLen=bLenDown
				newVect=[]
				if totBLen:
					for i in range4:
						tot=0.0
						for j in range4:
							tot+=mutMatrix[i][j]*entry2[-1][j]
						tot*=totBLen
						tot+=entry2[-1][i]
						newVect.append(tot*rootFreqs[i])
					totSum=sum(newVect)
					for i in range4:
						newVect[i]/=totSum
					probVect.append((6,pos,newVect))
				else:
					for i in range4:
						newVect.append(entry2[-1][i]*rootFreqs[i])
					totSum=sum(newVect)
					for i in range4:
						newVect[i]/=totSum
					probVect.append((6,pos,newVect))
		elif entry2[0]==5: #entry2 is N
			if entry1[0]<5:
				pos=min(entry1[1],entry2[1])
				if len(entry1)==2:
					if bLenUp:
						probVect.append((entry1[0],pos,bLenUp))
					else:
						probVect.append((entry1[0],pos))
				elif len(entry1)==3:
					if bLenUp:
						probVect.append((entry1[0],pos,entry1[2]+bLenUp))
					else:
						probVect.append((entry1[0],pos,entry1[2]))
				else:
					if bLenUp:
						probVect.append((entry1[0],pos,entry1[2],entry1[3]+bLenUp))
					else:
						probVect.append((entry1[0],pos,entry1[2],entry1[3]))

			else: #entry1 is "O", entry 2 is "N"
				pos+=1
				if len(entry1)==4:
					totBLen=entry1[2]
					if bLenUp:
						totBLen+=bLenUp
				elif bLenUp:
					totBLen=bLenUp
				else:
					totBLen=False
				newVect=[]
				if totBLen:
					for i in range4:
						tot=0.0
						for j in range4:
							tot+=entry1[-1][j]*mutMatrix[j][i]
						tot*=totBLen
						tot+=entry1[-1][i]
						newVect.append(tot)
					totSum=sum(newVect)
					for i in range4:
						newVect[i]/=totSum
					probVect.append((6,pos,newVect))
				else:
					probVect.append((6,pos,entry1[-1]))
		elif entry2[0]==entry1[0] and entry1[0]<5:
			pos=min(entry1[1],entry2[1])
			probVect.append((entry2[0],pos))
		else: #cases where the new genome list entry will likely be of type "O"
			if entry1[0]<5:
				if len(entry1)==2:
					totLen1=bLenUp
				else:
					totLen1=entry1[2]
					if bLenUp:
						totLen1+=bLenUp
					if len(entry1)==4:
						totLen1+=entry1[3]
			else:
				if len(entry1)==3:
					totLen1=bLenUp
				else:
					totLen1=entry1[2]
					if bLenUp:
						totLen1+=bLenUp
			
			if entry2[0]<5:
				if len(entry2)==2:
					totLen2=bLenDown
				else:
					totLen2=entry2[2]
					if bLenDown:
						totLen2+=bLenDown
			else:
				if len(entry2)==3:
					totLen2=bLenDown
				else:
					totLen2=entry2[2]
					if bLenDown:
						totLen2+=bLenDown
			if entry2[0]<5 and (not totLen2): #due to 0 distance, the entry will be of same type as entry2
				if (not totLen1) and entry1[0]<5:
					return None
				pos=min(entry1[1],entry2[1])
				probVect.append((entry2[0],pos))
			elif entry1[0]<5 and (not totLen1): #due to 0 distance, the entry will be of same type as entry1
				pos=min(entry1[1],entry2[1])
				probVect.append((entry1[0],pos))
			elif entry1[0]<5:
				if entry1[0]==4:
					i1=refIndeces[pos]
				else:
					i1=entry1[0]
				newVec=[]
				if len(entry1)==4:
					rootVec=list(rootFreqs)
					for i in range4:
						if i==i1:
							rootVec[i]*=(1.0+nonMutRates[i1]*(entry1[2]))
						else:
							rootVec[i]*=mutMatrix[i][i1]*(entry1[2])
					if bLenUp:
						lenToRoot=entry1[3]+bLenUp
					else:
						lenToRoot=entry1[3]
					for j in range4:
						tot=0.0
						for i in range4:
							tot+=mutMatrix[i][j]*rootVec[i]
						tot*=lenToRoot
						tot+=rootVec[j]
						newVec.append(tot)
				else:
					if totLen1:
						for i in range4:
							if i==i1:
								newVec.append(1.0+nonMutRates[i]*totLen1)
							else:
								newVec.append(mutMatrix[i1][i]*totLen1)
					else:
						for i in range4:
							if i==i1:
								newVec.append(1.0)
							else:
								newVec.append(0.0)
				if entry2[0]==6: #entry 2 is "O" and entry1 is a nucleotide
					for j in range4:
						tot=0.0
						if totLen2:
							for i in range4:
								tot+=mutMatrix[j][i]*entry2[-1][i]
							tot*=totLen2
						tot+=entry2[-1][j]
						newVec[j]*=tot
					sumV=sum(newVec)
					for i in range4:
						newVec[i]=newVec[i]/sumV
					state =simplfy(newVec,refIndeces[pos])
					pos+=1
					if state==6:
						probVect.append((6,pos,newVec))
					else:
						probVect.append((state,pos))
				else: #entry 1 and entry 2 are different nucleotides (possibly reference)
					if entry2[0]==4:
						i2=refIndeces[pos]
					else:
						i2=entry2[0]
					if totLen2:
						for i in range4:
							if i==i2:
								newVec[i]*=1.0+nonMutRates[i]*totLen2
							else:
								newVec[i]*=mutMatrix[i][i2]*totLen2
					else:
						for i in range4:
							if i!=i2:
								newVec[i]=0
					sumV=sum(newVec)
					for i in range4:
						newVec[i]=newVec[i]/sumV
					pos+=1
					probVect.append((6,pos,newVec))
			else: #entry1[0]==6:
				if totLen1:
					newVec=[]
					for i in range4:
						tot=0.0
						for j in range4:
							tot+=mutMatrix[j][i]*entry1[-1][j]
						tot*=totLen1
						tot+=entry1[-1][i]
						newVec.append(tot)
				else:
					newVec=list(entry1[-1])

				if entry2[0]==6:
					if totLen2:
						for i in range4:
							tot=0.0
							for j in range4:
								tot+=mutMatrix[i][j]*entry2[-1][j]
							tot*=totLen2
							tot+=entry2[-1][i]
							newVec[i]*=tot
					else:
						for i in range4:
							newVec[i]*=entry2[-1][i]
				else: #entry1 is "O" and entry2 is a nucleotide
					if entry2[0]==4:
						i2=refIndeces[pos]
					else:
						i2=entry2[0]
					if totLen2:
						for i in range4:
							if i==i2:
								newVec[i]*=(1.0+nonMutRates[i]*totLen2)
							else:
								newVec[i]*=mutMatrix[i][i2]*totLen2
					else:
						for i in range4:
							if i!=i2:
								newVec[i]=0.0
				sumV=sum(newVec)
				if not sumV:
					return None
				for i in range4:
					newVec[i]=newVec[i]/sumV
				state =simplfy(newVec,refIndeces[pos])
				pos+=1
				if state==6:
					probVect.append((6,pos,newVec))
				else:
					probVect.append((state,pos))

		if pos==lRef:
			break
		if pos==entry1[1]:
			indexEntry1+=1
			entry1=probVect1[indexEntry1]
		if pos==entry2[1]:
			indexEntry2+=1
			entry2=probVect2[indexEntry2]

	if verbose:
		print("Merged up-down ")
		print(probVect)
	#check if the final  probVect can be simplified by merging consecutive entries
	shorten(probVect)
	if verbose:
		print("Shortened up-down merging ")
		print(probVect)
	return probVect







#merge two lower child partial likelihood vectors to create a new one 
#(and also calculate the logLk of the merging if necessary, which is currently only useful for the root but could also be useful to calculate the overall total likelihood of the tree).
def mergeVectors(probVect1,bLen1,probVect2,bLen2,mutMatrix,returnLK=False):
	indexEntry1, indexEntry2, pos = 0, 0, 0
	probVect=[]
	cumulPartLk=0.0
	entry1=probVect1[indexEntry1]
	entry2=probVect2[indexEntry2]
	end=min(entry1[1],entry2[1])
	while True:
		if entry1[0]==5:
			if entry2[0]==5:
				pos=min(entry1[1],entry2[1])
				probVect.append((5,pos))
			elif entry2[0]<5:
				pos=min(entry1[1],entry2[1])
				if len(entry2)==2:
					if bLen2:
						probVect.append((entry2[0],pos,bLen2))
					else:
						probVect.append((entry2[0],pos))
				else:
					if bLen2:
						probVect.append((entry2[0],pos,entry2[2]+bLen2))
					else:
						probVect.append((entry2[0],pos,entry2[2]))		
			else: # case entry2 is "O" and entry1 is "N"
				pos+=1
				if len(entry2)==3:
					if bLen2:
						probVect.append((6,pos,bLen2,entry2[-1]))
					else:
						probVect.append((6,pos,entry2[-1]))
				else:
					if bLen2:
						probVect.append((6,pos,entry2[2]+bLen2,entry2[-1]))
					else:
						probVect.append((6,pos,entry2[2],entry2[-1]))
		elif entry2[0]==5: #entry2 is N
			if entry1[0]<5:
				pos=min(entry1[1],entry2[1])
				if len(entry1)==2:
					if bLen1:
						probVect.append((entry1[0],pos,bLen1))
					else:
						probVect.append((entry1[0],pos))
				else:
					if bLen1:
						probVect.append((entry1[0],pos,entry1[2]+bLen1))
					else:
						probVect.append((entry1[0],pos,entry1[2]))
			else: #entry1 is "O" and entry2 is "N"
				pos+=1
				if len(entry1)==3:
					if bLen1:
						probVect.append((6,pos,bLen1,entry1[-1]))
					else:
						probVect.append((6,pos,entry1[-1]))
				else:
					if bLen1:
						probVect.append((6,pos,entry1[2]+bLen1,entry1[-1]))
					else:
						probVect.append((6,pos,entry1[2],entry1[-1]))
					
		else: #entry1 and entry2 are not "N"
			if entry1[0]<5:
				if len(entry1)==2:
					totLen1=bLen1
				else:
					totLen1=entry1[2]
					if bLen1:
						totLen1+=bLen1
			else:
				if len(entry1)==3:
					totLen1=bLen1
				else:
					totLen1=entry1[2]
					if bLen1:
						totLen1+=bLen1
			
			if entry2[0]<5:
				if len(entry2)==2:
					totLen2=bLen2
				else:
					totLen2=entry2[2]
					if bLen2:
						totLen2+=bLen2
			else:
				if len(entry2)==3:
					totLen2=bLen2
				else:
					totLen2=entry2[2]
					if bLen2:
						totLen2+=bLen2

			if entry2[0]==entry1[0] and entry2[0]<5: #entry1 and entry2 are two identical nucleotides
				end=min(entry1[1],entry2[1])
				probVect.append((entry2[0],end))
				if returnLK:
					if entry2[0]==4:
						cumulPartLk+=(totLen1+totLen2)*(cumulativeRate[end]-cumulativeRate[pos])
					else:
						cumulPartLk+=nonMutRates[entry1[0]]*(totLen1+totLen2)
				pos=end
			elif (not totLen1) and (not totLen2) and entry1[0]<5 and entry2[0]<5: #0 distance between different nucleotides: merge is not possible
				if returnLK:
					return None, float("-inf")
				else:
					return None
			elif entry1[0]<5: #entry1 is a nucleotide
				if entry1[0]==4:
					i1=refIndeces[pos]
				else:
					i1=entry1[0]
				newVec=[]
				if totLen1:
					for i in range4:
						if i==i1:
							newVec.append(1.0+nonMutRates[i]*totLen1)
						else:
							newVec.append(mutMatrix[i][i1]*totLen1)
				else:
					newVec=[0.0,0.0,0.0,0.0]
					newVec[i1]=1.0

				if entry2[0]==6: #entry1 is a nucleotide and entry2 is "O"
					if totLen2:
						for j in range4:
							tot=0.0
							for i in range4:
								tot+=mutMatrix[j][i]*entry2[-1][i]
							tot*=totLen2
							tot+=entry2[-1][j]
							newVec[j]*=tot
					else:
						for j in range4:
							newVec[j]*=entry2[-1][j]
					sumV=sum(newVec)
					if not sumV:
						if returnLK:
							return None, float("-inf")
						else:
							return None
					for i in range4:
						newVec[i]=newVec[i]/sumV
					state =simplfy(newVec,refIndeces[pos])
					pos+=1
					if state==6:
						probVect.append((6,pos,newVec))
					else:
						probVect.append((state,pos))
					if returnLK:
						cumulPartLk+=log(sumV)
				else: #entry1 and entry2 are nucleotides
					if entry2[0]==4:
						i2=refIndeces[pos]
					else:
						i2=entry2[0]
					if totLen2:
						for i in range4:
							if i==i2:
								newVec[i]*=1.0+nonMutRates[i]*totLen2
							else:
								newVec[i]*=mutMatrix[i][i2]*totLen2
						sumV=sum(newVec)
						for i in range4:
							newVec[i]=newVec[i]/sumV
						pos+=1
						probVect.append((6,pos,newVec))
						if returnLK:
							cumulPartLk+=log(sumV)
					else:
						pos+=1
						probVect.append((entry2[0],pos))
						if returnLK:
							cumulPartLk+=log(newVec[i2])
				
			else: #entry1[0]==6:
				if totLen1:
					newVec=[]
					for i in range4:
						tot=0.0
						for j in range4:
							tot+=mutMatrix[i][j]*entry1[-1][j]
						tot*=totLen1
						tot+=entry1[-1][i]
						newVec.append(tot)
				else:
					newVec=list(entry1[-1])

				if entry2[0]==6:
					if totLen2:
						for i in range4:
							tot=0.0
							for j in range4:
								tot+=mutMatrix[i][j]*entry2[-1][j]
							tot*=totLen2
							tot+=entry2[-1][i]
							newVec[i]*=tot
					else:
						for i in range4:
							newVec[i]*=entry2[-1][i]
					sumV=sum(newVec)
					if not sumV:
						if returnLK:
							return None, float("-inf")
						else:
							return None
					for i in range4:
						newVec[i]=newVec[i]/sumV
					state =simplfy(newVec,refIndeces[pos])
					pos+=1
					if state==6:
						probVect.append((6,pos,newVec))
					else:
						probVect.append((state,pos))
					if returnLK:
						cumulPartLk+=log(sumV)
				else: #entry2 is a nucleotide and entry1 is "O"
					if entry2[0]==4:
						i2=refIndeces[pos]
					else:
						i2=entry2[0]
					if totLen2:
						for i in range4:
							if i==i2:
								newVec[i]*=(1.0+nonMutRates[i]*totLen2)
							else:
								newVec[i]*=mutMatrix[i][i2]*totLen2
						sumV=sum(newVec)
						for i in range4:
							newVec[i]=newVec[i]/sumV
						state =simplfy(newVec,refIndeces[pos])
						pos+=1
						if state==6:
							probVect.append((6,pos,newVec))
						else:
							probVect.append((state,pos))
						if returnLK:
							cumulPartLk+=log(sumV)
					else:
						if not newVec[i2]:
							if returnLK:
								return None, float("-inf")
							else:
								return None
						pos+=1
						probVect.append((entry2[0],pos))
						if returnLK:
							cumulPartLk+=log(newVec[i2])

		if pos==lRef:
			break
		if pos==entry1[1]:
			indexEntry1+=1
			entry1=probVect1[indexEntry1]
		if pos==entry2[1]:
			indexEntry2+=1
			entry2=probVect2[indexEntry2]

	if verbose:
		print("Merged vector at root ")
		print(probVect)
	#check if the final  probVect can be simplified by merging consecutive entries
	shorten(probVect)
	if verbose:
		print("Shortened root merging ")
		print(probVect)
	if returnLK:
		return probVect, cumulPartLk
	else:
		return probVect











#calculate the probability that results from combining a lower likelihood genome list of the root with root frequencies.
#Could maybe be combined with function rootVector() ?
def findProbRoot(probVect):
	logLK=0.0
	logFactor=1.0
	pos=0
	for entry in probVect:
		if entry[0]==4:
			for i in range4:
				logLK+=rootFreqsLog[i]*(cumulativeBases[entry[1]][i]-cumulativeBases[pos][i])
		elif entry[0]<4:
			logLK+=rootFreqsLog[entry[0]]
		elif entry[0]==6:
			tot=0.0
			for i in range4:
				tot+=rootFreqs[i]*entry[-1][i]
			logFactor*=tot
		pos=entry[1]
	logLK+=log(logFactor)
	return logLK








#Check if two genome lists represent the same partial likelihoods or not.
#This is used to traverse the tree and only update genome lists if they changed after a change in the tree.
def areVectorsDifferent(probVect1,probVect2):
	if probVect2==None:
		return True
	indexEntry1, indexEntry2, pos = 0, 0, 0
	entry1=probVect1[indexEntry1]
	entry2=probVect2[indexEntry2]
	while True:
		if entry1[0]!=entry2[0]:	
			return True
		if len(entry1)!=len(entry2):
			return True
		if entry1[0]<5:
			if len(entry1)>2:
				if abs(entry1[2] - entry2[2])>thresholdProb:
					return True
				if len(entry1)==4:
					if abs(entry1[3] - entry2[3])>thresholdProb:
						return True
		if entry1[0]==6:
			if len(entry1)==4:
				if abs(entry1[2] - entry2[2])>thresholdProb:
					return True
			for i in range4:
				diffVal=abs(entry1[-1][i] - entry2[-1][i])
				#example thresholdDiffForUpdate=0.01
				#example thresholdFoldChangeUpdate=1.5
				if diffVal:
					if (not entry1[-1][i]) or (not entry2[-1][i]):
						return True
					if diffVal>thresholdDiffForUpdate or (diffVal>thresholdProb and ( (diffVal/entry1[-1][i]>thresholdFoldChangeUpdate)  or  (diffVal/entry2[-1][i]>thresholdFoldChangeUpdate))):
						return True
		#update pos, end, etc
		pos=min(entry1[1],entry2[1])
		if pos==lRef:
			break
		if pos==entry1[1]:
			indexEntry1+=1
			entry1=probVect1[indexEntry1]
		if pos==entry2[1]:
			indexEntry2+=1
			entry2=probVect2[indexEntry2]
	return False







#create further mid-nodes for longer branches, so to make it faster to calculate appending probabilities.
#This is useful when we have located the best (internal) node for appending a new samples, and we want to know which child branch might be best for the appending.
#These node.furtherMidNodes are not used very often (not as much as probVectTot and probVectTotUp), and so represent a compromise between memory and time, and it might be convenient to get rid of them altogether.
def createFurtherMidNodes(node,vectUp):
	node.furtherMidNodes=[]
	newBLen2=node.dist/4
	while newBLen2>=minBLenForMidNode/2:
		newProbVect2=mergeVectorsUpDown(vectUp,newBLen2,node.probVect,node.dist-newBLen2,mutMatrix)
		node.furtherMidNodes.append(newProbVect2)
		newBLen2=newBLen2/2







#if updating genome lists creates an inconsistency, this function can increase the length of a 0-length branch to resolve the inconsistency.
def updateBLen(node,childNum,mutMatrix):
	cNode=node.children[childNum]
	if childNum:
		vectUp=node.probVectUpLeft
	else:
		vectUp=node.probVectUpRight
	bestLK=appendProbNode(vectUp,cNode.probVect,oneMutBLen,mutMatrix)
	bestLen=oneMutBLen
	while bestLen>minBLen:
		newBLen=bestLen/2
		newLK=appendProbNode(vectUp,cNode.probVect,newBLen,mutMatrix)
		if newLK>bestLK:
			bestLK=newLK
			bestLen=newBLen
		else:
			break
	if bestLen>0.7*oneMutBLen:
		while bestLen<maxBLen:
			newBLen=bestLen*2
			newLK=appendProbNode(vectUp,cNode.probVect,newBLen,mutMatrix)
			if newLK>bestLK:
				bestLK=newLK
				bestLen=newBLen
			else:
				break
	cNode.dist=bestLen
	updatePartialsFromTop(cNode,vectUp,mutMatrix)
	updatePartialsFromBottom(node,cNode.probVect,childNum,cNode,mutMatrix)









#Update partial likelihood vectors in the tree after the addition of a new tip.
#This function traverses the tree downward.
def updatePartialsFromTop(node,probVectUp,mutMatrix):
	finalBool=False
	if node.dist : #if necessary, update the total probabilities at the mid node.
		newTot=mergeVectorsUpDown(probVectUp,node.dist/2,node.probVect,node.dist/2,mutMatrix)
		if newTot==None:
			if node.up.children[0]==node:
				updateBLen(node.up,0,mutMatrix)
			else:
				updateBLen(node.up,1,mutMatrix)
			return True
		node.probVectTotUp=newTot
		if node.dist>=2*minBLenForMidNode:
			createFurtherMidNodes(node,probVectUp)
	if len(node.children)==0 and node.dist: #if necessary, update the total probability vector at the terminal node.
		newTot=mergeVectorsUpDown(probVectUp,node.dist,node.probVect,False,mutMatrix)
		if newTot==None:
			if node.up.children[0]==node:
				updateBLen(node.up,0,mutMatrix)
			else:
				updateBLen(node.up,1,mutMatrix)
			return True
		node.probVectTot=newTot
		node.dirty=True
	elif len(node.children)>0: #at valid internal node, update upLeft and upRight, and if necessary pass the function on to children.
		child0Vect=node.children[0].probVect
		child1Vect=node.children[1].probVect
		dist0=node.children[0].dist
		dist1=node.children[1].dist
		newUpRight=mergeVectorsUpDown(probVectUp,node.dist,child1Vect,dist1,mutMatrix)
		newUpLeft=mergeVectorsUpDown(probVectUp,node.dist,child0Vect,dist0,mutMatrix)
		if newUpRight==None:
			if not node.dist :
				if node.up.children[0]==node:
					updateBLen(node.up,0,mutMatrix)
				else:
					updateBLen(node.up,1,mutMatrix)
			else:
				updateBLen(node,1,mutMatrix)
			return True
		if newUpLeft==None:
			if not node.dist :
				if node.up.children[0]==node:
					updateBLen(node.up,0,mutMatrix)
				else:
					updateBLen(node.up,1,mutMatrix)
			else:
				updateBLen(node,0,mutMatrix)
			return True
		updatedTot=False
		if areVectorsDifferent(node.probVectUpRight,newUpRight):
			node.probVectUpRight=newUpRight
			node.dirty=True
			if node.dist:
				newTot=mergeVectorsUpDown(newUpRight,False,child0Vect,dist0,mutMatrix)
				if newTot==None:
					updateBLen(node,0,mutMatrix)
					return True
				node.probVectTot=newTot
				updatedTot=True
			if updatePartialsFromTop(node.children[0],newUpRight,mutMatrix):
				finalBool=True
		if areVectorsDifferent(node.probVectUpLeft,newUpLeft):
			node.probVectUpLeft=newUpLeft
			node.dirty=True
			if (not updatedTot) and node.dist:
				newTot=mergeVectorsUpDown(newUpLeft,False,child1Vect,dist1,mutMatrix)
				if newTot==None:
					updateBLen(node,1,mutMatrix)
					return True
				node.probVectTot=newTot
			if updatePartialsFromTop(node.children[1],newUpLeft,mutMatrix):
				finalBool=True
	return finalBool









#Update partial likelihood vectors in the tree after the addition of a new tip.
#This function traverses the tree upward (but also branches downward).
def updatePartialsFromBottom(node,probVectDown,childNum,childNode,mutMatrix):
	if childNum==-1:
		if node.children[0]==childNode:
			childNum=0
		else:
			childNum=1
	otherChildNum=1-childNum
	childDist=node.children[childNum].dist
	otherChildDist=node.children[otherChildNum].dist
	otherChildVect=node.children[otherChildNum].probVect
	newVect=mergeVectors(otherChildVect,otherChildDist,probVectDown,childDist,mutMatrix)
	finalBool=False
	if newVect==None:
		if not childDist:
			updateBLen(node,childNum,mutMatrix)
		else:
			updateBLen(node,otherChildNum,mutMatrix)
		return True
	if verbose:
		print("Updating vectors from bottom, new vec")
		print(newVect)
		print("Old vec:")
		print(node.probVect)
	updatedTot=False
	if node.dist:
		if childNum==0:
			newTot=mergeVectorsUpDown(node.probVectUpRight,False,probVectDown,childDist,mutMatrix)
			if newTot==None:
				updateBLen(node,0,mutMatrix)
				return True
		else:
			newTot=mergeVectorsUpDown(node.probVectUpLeft,False,probVectDown,childDist,mutMatrix)
			if newTot==None:
				updateBLen(node,1,mutMatrix)
				return True
		node.probVectTot=newTot
		if verbose:
			print("new tot vect")
			print(node.probVectTot)
		updatedTot=True
	oldProbVect=node.probVect
	node.probVect=newVect
	if areVectorsDifferent(node.probVect,oldProbVect):
		node.dirty=True
		if node.up != None:
			if updatePartialsFromBottom(node.up,newVect,-1,node,mutMatrix):
				finalBool=True
		if verbose:
			print("Old upRight and UpLeft:")
			print(node.probVectUpRight)
			print(node.probVectUpLeft)
	if node.up != None:
		if verbose:
			print("Trying to move up while updating partials")
		if node.up.children[0]==node:
			vectUp=node.up.probVectUpRight
		else:
			vectUp=node.up.probVectUpLeft
		if childNum==0:
			newUpLeftVect=mergeVectorsUpDown(vectUp,node.dist,probVectDown,childDist,mutMatrix)
			if newUpLeftVect==None:
				updateBLen(node,0,mutMatrix)
				return True
			if areVectorsDifferent(node.probVectUpLeft,newUpLeftVect):
				node.dirty=True
				if verbose:
					print("Moving up the update from left; new UpLeft and old UpLeft")
					print(newUpLeftVect)
					print(node.probVectUpLeft)
				node.probVectUpLeft=newUpLeftVect
				if (not updatedTot) and node.dist:
					if verbose:
						print("Updating tot: old tot and new tot")
						print(node.probVectTot)
					newTot=mergeVectorsUpDown(newUpLeftVect,False,otherChildVect,otherChildDist,mutMatrix)
					if newTot==None:
						updateBLen(node,1,mutMatrix)
						return True
					node.probVectTot=newTot
					if verbose:
						print(newTot)
					updatedTot=True
				if updatePartialsFromTop(node.children[1],newUpLeftVect,mutMatrix):
					finalBool=True
		else:
			newUpRightVect=mergeVectorsUpDown(vectUp,node.dist,probVectDown,childDist,mutMatrix)
			if newUpRightVect==None:
				updateBLen(node,1,mutMatrix)
				return True
			if areVectorsDifferent(node.probVectUpRight,newUpRightVect):
				node.dirty=True
				if verbose:
					print("Moving up the update from right")
				node.probVectUpRight=newUpRightVect
				if (not updatedTot) and node.dist:
					newTot=mergeVectorsUpDown(newUpRightVect,False,otherChildVect,otherChildDist,mutMatrix)
					if newTot==None:
						updateBLen(node,0,mutMatrix)
						return True
					node.probVectTot=newTot
					updatedTot=True
				if updatePartialsFromTop(node.children[0],newUpRightVect,mutMatrix):
					finalBool=True
		if updatedTot:
			node.probVectTotUp=mergeVectorsUpDown(vectUp,node.dist/2,node.probVect,node.dist/2,mutMatrix)
			if node.dist>=2*minBLenForMidNode:
				createFurtherMidNodes(node,vectUp)
	#case node is root
	else:
		if verbose:
			print("reached root while moving up")
		if childNum==0:
			newUpLeftVect=rootVector(probVectDown,childDist,mutMatrix)
			if areVectorsDifferent(node.probVectUpLeft,newUpLeftVect):
				node.probVectUpLeft=newUpLeftVect
				if not updatedTot:
					node.dirty=True
					#node.probVectTot=mergeVectorsUpDown(newUpLeftVect,False,otherChildVect,otherChildDist,mutMatrix)
					newTot=mergeVectorsUpDown(newUpLeftVect,False,otherChildVect,otherChildDist,mutMatrix)
					if newTot==None:
						updateBLen(node,1,mutMatrix)
						return True
					node.probVectTot=newTot
				updatePartialsFromTop(node.children[1],newUpLeftVect,mutMatrix)

		else:
			newUpRightVect=rootVector(probVectDown,childDist,mutMatrix)
			if areVectorsDifferent(node.probVectUpRight,newUpRightVect):
				node.probVectUpRight=newUpRightVect
				if not updatedTot:
					node.dirty=True
					#node.probVectTot=mergeVectorsUpDown(newUpRightVect,False,otherChildVect,otherChildDist,mutMatrix)
					newTot=mergeVectorsUpDown(newUpRightVect,False,otherChildVect,otherChildDist,mutMatrix)
					if newTot==None:
						updateBLen(node,0,mutMatrix)
						return True
					node.probVectTot=newTot
				updatePartialsFromTop(node.children[0],newUpRightVect,mutMatrix)
	if verbose:
		print("New upRight, UpLeft and tot:")
		print(node.probVectUpRight)
		print(node.probVectUpLeft)
		print(node.probVectTot)
	return finalBool










#function to add new mutation events in new sample to the pre-exisitng pseudocounts so to improve the estimate of the substitution rates.
# probVect1 is the genome list at the node where the appending happens; probVect2 is the genome list for the new sample.
def updatePesudoCounts(probVect1,probVect2,pseudoMutCounts):
	if model!="JC":
		indexEntry1, indexEntry2, pos = 0, 0, 0
		entry1=probVect1[indexEntry1]
		entry2=probVect2[indexEntry2]
		while True:
			if entry1[0]!=entry2[0] and entry1[0]<5 and entry2[0]<5:
				if entry1[0]==4:
					pseudoMutCounts[refIndeces[pos]][entry2[0]]+=1		
				elif entry2[0]==4:
					pseudoMutCounts[entry1[0]][refIndeces[pos]]+=1
				else:
					pseudoMutCounts[entry1[0]][entry2[0]]+=1
				pos+=1
			else:
				pos=min(entry1[1],entry2[1])

			if pos==lRef:
				break
			if pos==entry1[1]:
				indexEntry1+=1
				entry1=probVect1[indexEntry1]
			if pos==entry2[1]:
				indexEntry2+=1
				entry2=probVect2[indexEntry2]







#we know that sample "sample", with partials "newPartials", is best placed near a node resulting in logLK contribution of newChildLK
# now explore exactly which position near node is best for placement (direct child, sibling or sibling of existing child), and which are the best branch lengths,
#then add the sample at that position of the tree, and update all the internal probability vectors.
def placeSampleOnTreeNew(node,newPartials,sample,newChildLK,isMidNode, bestUpLK, bestDownLK, bestDownNode,mutMatrix,pseudoMutCounts, adjustBLen):
	if adjustBLen:
		factor=float(bLenAdjustment)
		bLen=oneMutBLen*bLenAdjustment
	else:
		factor=1.0
		bLen=oneMutBLen

	if isMidNode:
		if node==node.up.children[0]:
			child=0
			vectUp=node.up.probVectUpRight
		else:
			child=1
			vectUp=node.up.probVectUpLeft
		bestSplit=0.5
		bestSplitLK=newChildLK
		childBestVect=node.probVectTotUp
		newSplit=0.25
		#try different positions on the existing branch
		while newSplit*node.dist>minBLen:
			probVectParentNew=mergeVectorsUpDown(vectUp,node.dist*newSplit,node.probVect,node.dist*(1.0-newSplit),mutMatrix)
			probChild=appendProb(probVectParentNew,newPartials,bLen,mutMatrix)
			if probChild>bestSplitLK:
				bestSplitLK=probChild
				bestSplit=newSplit
				childBestVect=probVectParentNew
			else:
				break
			newSplit=bestSplit/2
		if bestSplit>0.49:
			newSplit=0.25
			while newSplit*node.dist>minBLen:
				probVectParentNew=mergeVectorsUpDown(vectUp,node.dist*(1.0-newSplit),node.probVect,node.dist*newSplit,mutMatrix)
				probChild=appendProb(probVectParentNew,newPartials,bLen,mutMatrix)
				if probChild>bestSplitLK:
					bestSplitLK=probChild
					bestSplit=newSplit
					childBestVect=probVectParentNew
				else:
					bestSplit=1.0-bestSplit
					break
				newSplit=bestSplit/2
		#now try different lengths for the new branch
		LK1=bestSplitLK
		bestLen=oneMutBLen*factor
		while bestLen>minBLen:
			newBLen=bestLen/2
			probChild=appendProb(childBestVect,newPartials,newBLen,mutMatrix)
			if probChild>LK1:
				LK1=probChild
				bestLen=newBLen
			else:
				break
		if bestLen>0.7*oneMutBLen*factor:
			while bestLen<maxBLen:
				newBLen=bestLen*2
				probChild=appendProb(childBestVect,newPartials,newBLen,mutMatrix)
				if probChild>LK1:
					LK1=probChild
					bestLen=newBLen
				else:
					break
		if bestLen<minBLen:
			LK0=appendProb(childBestVect,newPartials,False,mutMatrix)
			if LK0>LK1:
				bestLen=False
		#now create new internal node and append child to it
		distTop=node.dist*bestSplit
		newInternalNode=Tree()
		node.up.children[child]=newInternalNode
		newInternalNode.up=node.up
		distBottom=node.dist*(1.0-bestSplit)
		newInternalNode.add_child(node)
		node.up=newInternalNode
		node.dist=distBottom
		newNode=Tree(name=sample,dist=bestLen)
		newNode.minorSequences=[]
		newNode.up=newInternalNode
		newInternalNode.add_child(newNode)
		newInternalNode.dist=distTop
		newInternalNode.children[1].probVect=newPartials
		newInternalNode.probVectUpRight=mergeVectorsUpDown(vectUp,distTop,newPartials,bestLen,mutMatrix)
		newInternalNode.probVectUpLeft=childBestVect
		newInternalNode.probVect=mergeVectors(node.probVect,node.dist,newPartials,bestLen,mutMatrix)
		newInternalNode.probVectTotUp=mergeVectorsUpDown(vectUp,distTop/2,newInternalNode.probVect,distTop/2,mutMatrix)
		newInternalNode.probVectTot=mergeVectorsUpDown(childBestVect,False,newPartials,bestLen,mutMatrix)
		if newInternalNode.probVectTot==None:
			print("Problem, None vector when placing sample, below node")
			print(childBestVect)
			print(vectUp)
			print(newPartials)
			print(bestLen)
			print(distTop)
			print(distBottom)
		if distTop>=2*minBLenForMidNode:
			createFurtherMidNodes(newInternalNode,vectUp)
		if bestLen:
			newNode.probVectTot=mergeVectorsUpDown(childBestVect,bestLen,newPartials,False,mutMatrix)
			newInternalNode.children[1].probVectTotUp=mergeVectorsUpDown(childBestVect,bestLen/2,newPartials,bestLen/2,mutMatrix)
			if bestLen>=2*minBLenForMidNode:
				createFurtherMidNodes(newInternalNode.children[1],childBestVect)
		updatePesudoCounts(childBestVect,newPartials,pseudoMutCounts)
		if verbose:
			print("new internal node added to tree")
			print(newInternalNode.probVect)
			print(newInternalNode.probVectUpRight)
			print(newInternalNode.probVectUpLeft)
			print(newInternalNode.probVectTot)
		updatePartialsFromTop(node,newInternalNode.probVectUpRight,mutMatrix)
		updatePartialsFromBottom(newInternalNode.up,newInternalNode.probVect,child,newInternalNode,mutMatrix)

	#best lk so far is for appending directly to existing node
	else:
		if bestDownNode!=None:
			if bestDownNode==bestDownNode.up.children[0]:
				child=0
				vectUp=bestDownNode.up.probVectUpRight
			else:
				child=1
				vectUp=bestDownNode.up.probVectUpLeft
			bestSplit=0.5
			bestSplitLK=bestDownLK
			childBestVect=bestDownNode.probVectTotUp
			newSplit=0.25
			while newSplit*bestDownNode.dist>minBLen:
				probVectParentNew=mergeVectorsUpDown(vectUp,bestDownNode.dist*newSplit,bestDownNode.probVect,bestDownNode.dist*(1.0-newSplit),mutMatrix)
				probChild=appendProb(probVectParentNew,newPartials,bLen,mutMatrix)
				if probChild>bestSplitLK:
					bestSplitLK=probChild
					bestSplit=newSplit
					childBestVect=probVectParentNew
				else:
					break
				newSplit=bestSplit/2
			bestChildLK=bestSplitLK
			bestChildSplit=bestSplit
		else:
			bestChildLK=float("-inf")

		#if node is root, try to place as sibling of the current root.
		if node.up==None:
			probOldRoot = findProbRoot(node.probVect)
			probVectRoot,probRoot = mergeVectors(node.probVect,oneMutBLen,newPartials,oneMutBLen*factor,mutMatrix,returnLK=True)
			probRoot+= findProbRoot(probVectRoot)
			parentLKdiff=probRoot-probOldRoot
			bestRootBL=oneMutBLen
			parentBestVect=probVectRoot
			newBL=0.5*oneMutBLen
			while newBL>minBLen:
				probVectRoot,probRoot = mergeVectors(node.probVect,newBL,newPartials,oneMutBLen*factor,mutMatrix,returnLK=True)
				probRoot+= findProbRoot(probVectRoot)
				newDiff=probRoot-probOldRoot
				if newDiff>parentLKdiff:
					parentLKdiff=newDiff
					bestRootBL=newBL
					parentBestVect=probVectRoot
				else:
					break
				newBL=bestRootBL/2

		else: #node is not root
			if node==node.up.children[0]:
				child=0
				vectUp=node.up.probVectUpRight
			else:
				child=1
				vectUp=node.up.probVectUpLeft
			bestSplit=0.5
			bestSplitLK=bestUpLK
			parentBestVect=node.probVectTotUp
			newSplit=0.25
			while newSplit*node.dist>minBLen:
				probVectParentNew=mergeVectorsUpDown(vectUp,node.dist*(1.0-newSplit),node.probVect,node.dist*newSplit,mutMatrix)
				probChild=appendProb(probVectParentNew,newPartials,bLen,mutMatrix)
				if probChild>bestSplitLK:
					bestSplitLK=probChild
					bestSplit=newSplit
					parentBestVect=probVectParentNew
				else:
					break
				newSplit=bestSplit/2
			parentLKdiff=bestSplitLK
			bestParentSplit=bestSplit
		
		#Best placement is below node: add internal node below "node"
		if bestChildLK>=parentLKdiff and bestChildLK>=newChildLK:
			if bestDownNode==bestDownNode.up.children[0]:
				child=0
				vectUp=bestDownNode.up.probVectUpRight
			else:
				child=1
				vectUp=bestDownNode.up.probVectUpLeft

			LK1=bestChildLK
			bestLen=oneMutBLen*factor
			while bestLen>minBLen:
				newBLen=bestLen/2
				probChild=appendProb(childBestVect,newPartials,newBLen,mutMatrix)
				if probChild>LK1:
					LK1=probChild
					bestLen=newBLen
				else:
					break
			if bestLen>0.7*oneMutBLen*factor:
				while bestLen<maxBLen:
					newBLen=bestLen*2
					probChild=appendProb(childBestVect,newPartials,newBLen,mutMatrix)
					if probChild>LK1:
						LK1=probChild
						bestLen=newBLen
					else:
						break
			if bestLen<minBLen:
				LK0=appendProb(childBestVect,newPartials,False,mutMatrix)
				if LK0>LK1:
					bestLen=False
			#now create new internal node and append child to it
			newInternalNode=Tree()
			bestDownNode.up.children[child]=newInternalNode
			newInternalNode.up=bestDownNode.up
			distBottom=bestDownNode.dist*(1.0-bestChildSplit)
			distTop=bestDownNode.dist*bestChildSplit
			bestDownNode.up=newInternalNode
			bestDownNode.dist=distBottom
			newInternalNode.add_child(bestDownNode)
			newNode=Tree(name=sample,dist=bestLen)
			newNode.minorSequences=[]
			newNode.up=newInternalNode
			newInternalNode.add_child(newNode)
			newInternalNode.dist=distTop
			newInternalNode.children[1].probVect=newPartials
			newInternalNode.probVectUpRight=mergeVectorsUpDown(vectUp,newInternalNode.dist,newPartials,bestLen,mutMatrix)
			newInternalNode.probVectUpLeft=childBestVect
			newInternalNode.probVect=mergeVectors(bestDownNode.probVect,bestDownNode.dist,newPartials,bestLen,mutMatrix)
			newInternalNode.probVectTotUp=mergeVectorsUpDown(vectUp,distTop/2,newInternalNode.probVect,distTop/2,mutMatrix)
			newInternalNode.probVectTot=mergeVectorsUpDown(childBestVect,False,newPartials,bestLen,mutMatrix)
			if newInternalNode.probVectTot==None:
					print("Problem, None vector when placing sample, below node")
					print(childBestVect)
					print(vectUp)
					print(newPartials)
					print(bestLen)
					print(distTop)
					print(distBottom)
			if distTop>=2*minBLenForMidNode:
				createFurtherMidNodes(newInternalNode,vectUp)
			if bestLen:
				newNode.probVectTot=mergeVectorsUpDown(childBestVect,bestLen,newPartials,False,mutMatrix)
				newNode.probVectTotUp=mergeVectorsUpDown(childBestVect,bestLen/2,newPartials,bestLen/2,mutMatrix)
				if bestLen>=2*minBLenForMidNode:
					createFurtherMidNodes(newNode,childBestVect)
			updatePesudoCounts(childBestVect,newPartials,pseudoMutCounts)
			if verbose:
				print("new internal node added to tree")
				print(newInternalNode.probVect)
				print(newInternalNode.probVectUpRight)
				print(newInternalNode.probVectUpLeft)
				print(newInternalNode.probVectTot)
			updatePartialsFromTop(bestDownNode,newInternalNode.probVectUpRight,mutMatrix)
			updatePartialsFromBottom(newInternalNode.up,newInternalNode.probVect,child,newInternalNode,mutMatrix)
		
		#add new parent to "node"
		else:
			#new parent is actually part of a polytomy since best placement is exactly at the node
			if newChildLK>=parentLKdiff:
				bestRootBL=False
				bestParentSplit=False
				parentLKdiff=newChildLK
				parentBestVect=node.probVectTot
				if node.up==None:
					probVectRoot,probRoot = mergeVectors(node.probVect,False,newPartials,oneMutBLen*factor,mutMatrix,returnLK=True)
					parentBestVect=probVectRoot

			#add parent to the root
			if node.up==None:
				#now try different lengths for right branch
				bestLen2=oneMutBLen*factor
				while bestLen2>minBLen:
					newBLen=bestLen2/2
					newProbVectRoot,newProbRoot = mergeVectors(node.probVect,bestRootBL,newPartials,newBLen,mutMatrix,returnLK=True)
					newProbRoot+= findProbRoot(newProbVectRoot)
					LKdiffRoot=newProbRoot-probOldRoot
					if LKdiffRoot>parentLKdiff:
						parentLKdiff=LKdiffRoot
						bestLen2=newBLen
						probVectRoot=newProbVectRoot
					else:
						break
				if bestLen2>0.7*oneMutBLen*factor:
					while bestLen2<maxBLen:
						newBLen=bestLen2*2
						newProbVectRoot,newProbRoot = mergeVectors(node.probVect,bestRootBL,newPartials,newBLen,mutMatrix,returnLK=True)
						newProbRoot+= findProbRoot(newProbVectRoot)
						LKdiffRoot=newProbRoot-probOldRoot
						if LKdiffRoot>parentLKdiff:
							parentLKdiff=LKdiffRoot
							bestLen2=newBLen
							probVectRoot=newProbVectRoot
						else:
							break
				if bestLen2<minBLen:
					newProbVectRoot,newProbRoot = mergeVectors(node.probVect,bestRootBL,newPartials,False,mutMatrix,returnLK=True)
					newProbRoot+= findProbRoot(newProbVectRoot)
					LK0=newProbRoot-probOldRoot
					if LK0>parentLKdiff:
						bestLen2=False
						parentLKdiff=LK0
						probVectRoot=newProbVectRoot

				newRoot=Tree()
				newRoot.probVect=probVectRoot
				newRoot.probVectTot=rootVector(probVectRoot,False,mutMatrix)
				newRoot.probVectUpRight=rootVector(newPartials,bestLen2,mutMatrix)
				newRoot.probVectUpLeft=rootVector(node.probVect,bestRootBL,mutMatrix)
				if newRoot.probVectTot==None:
					print("Problem, None vector when placing sample, new root")
					print(probVectRoot)
					print(node.probVect)
					print(newPartials)
					print(bestLen2)
					print(bestRootBL)
				node.up=newRoot
				node.dist=bestRootBL
				if not bestRootBL:
					node.probVectTot=None
					node.probVectTotUp=None
					node.furtherMidNodes=None
				newRoot.add_child(node)
				newNode=Tree(name=sample,dist=bestLen2)
				newNode.minorSequences=[]
				newNode.up=newRoot
				newRoot.add_child(newNode)
				newNode.probVect=newPartials
				if bestLen2:
					newNode.probVectTot=mergeVectorsUpDown(newRoot.probVectUpLeft,bestLen2,newPartials,False,mutMatrix)
					newNode.probVectTotUp=mergeVectorsUpDown(newRoot.probVectUpLeft,bestLen2/2,newPartials,bestLen2/2,mutMatrix)
					if bestLen2>=2*minBLenForMidNode:
						createFurtherMidNodes(newRoot.children[1],newRoot.probVectUpLeft)
				if verbose:
					print("new root added to tree")
					print(newRoot.probVect)
					print(newRoot.children[0].probVect)
					print(newNode.probVect)
				updatePartialsFromTop(node,newRoot.probVectUpRight,mutMatrix)
				return newRoot

			#add parent to non-root node
			else:
				if node==node.up.children[0]:
					child=0
					vectUp=node.up.probVectUpRight
				else:
					child=1
					vectUp=node.up.probVectUpLeft
				
				#now try different lengths for the new branch
				LK1=parentLKdiff
				bestLen=oneMutBLen*factor
				while bestLen>minBLen:
					newBLen=bestLen/2
					probChild=appendProb(parentBestVect,newPartials,newBLen,mutMatrix)
					if probChild>LK1:
						LK1=probChild
						bestLen=newBLen
					else:
						break
				if bestLen>0.7*oneMutBLen*factor:
					while bestLen<maxBLen:
						newBLen=bestLen*2
						probChild=appendProb(parentBestVect,newPartials,newBLen,mutMatrix)
						if probChild>LK1:
							LK1=probChild
							bestLen=newBLen
						else:
							break
				if bestLen<minBLen:
					LK0=appendProb(parentBestVect,newPartials,False,mutMatrix)
					if LK0>LK1:
						bestLen=False
				#now create new internal node and append child to it
				newInternalNode=Tree()
				node.up.children[child]=newInternalNode
				newInternalNode.up=node.up
				if bestParentSplit:
					distBottom=node.dist*bestParentSplit
					distTop=node.dist*(1.0-bestParentSplit)
				else:
					distBottom=False
					distTop=node.dist
					node.probVectTot=None
					node.probVectTotUp=None
					node.furtherMidNodes=None
				node.dist=distBottom
				node.up=newInternalNode
				newInternalNode.add_child(node)
				newNode=Tree(name=sample,dist=bestLen)
				newNode.minorSequences=[]
				newNode.up=newInternalNode
				newInternalNode.add_child(newNode)
				newInternalNode.dist=distTop
				newInternalNode.children[1].probVect=newPartials
				newInternalNode.probVectUpRight=mergeVectorsUpDown(vectUp,newInternalNode.dist,newPartials,bestLen,mutMatrix)
				newInternalNode.probVectUpLeft=parentBestVect
				newInternalNode.probVect=mergeVectors(node.probVect,node.dist,newPartials,bestLen,mutMatrix)
				newInternalNode.probVectTotUp=mergeVectorsUpDown(vectUp,distTop/2,newInternalNode.probVect,distTop/2,mutMatrix)
				newInternalNode.probVectTot=mergeVectorsUpDown(parentBestVect,False,newPartials,bestLen,mutMatrix)
				if newInternalNode.probVectTot==None:
					print("Problem, None vector when placing sample, new parent")
					print(parentBestVect)
					print(vectUp)
					print(newPartials)
					print(bestLen)
					print(distTop)
					print(distBottom)
				if distTop>=2*minBLenForMidNode:
					createFurtherMidNodes(newInternalNode,vectUp)
				if bestLen:
					newNode.probVectTot=mergeVectorsUpDown(parentBestVect,bestLen,newPartials,False,mutMatrix)
					newNode.probVectTotUp=mergeVectorsUpDown(parentBestVect,bestLen/2,newPartials,bestLen/2,mutMatrix)
					if bestLen>=2*minBLenForMidNode:
						createFurtherMidNodes(newNode,parentBestVect)
				updatePesudoCounts(parentBestVect,newPartials,pseudoMutCounts)
				if verbose:
					print("new internal node added to tree")
					print(newInternalNode.probVect)
					print(newInternalNode.probVectUpRight)
					print(newInternalNode.probVectUpLeft)
					print(newInternalNode.probVectTot)
				updatePartialsFromTop(node,newInternalNode.probVectUpRight,mutMatrix)
				updatePartialsFromBottom(newInternalNode.up,newInternalNode.probVect,child,newInternalNode,mutMatrix)

	return None







#set all descendant nodes to dirty.
#So far thses flags are used to prevent traversing the same part of the tree multiple times.
def setAllDirty(node):
	nextLeaves=[node]
	while nextLeaves:
		nextNode=nextLeaves.pop()
		nextNode.dirty=True
		for c in nextNode.children:
			nextLeaves.append(c)








#function to calculate likelihood cost of appending node to parent node 
#differently from appendProb, this allows the bottom node to be internal, not just a sample.
def appendProbNode(probVectP,probVectC,bLen,mutMatrix):
	Lkcost, indexEntry1, indexEntry2, totalFactor, pos = 0.0, 0, 0, 1.0, 0
	entry1=probVectP[indexEntry1]
	entry2=probVectC[indexEntry2]
	end=min(entry1[1],entry2[1])
	contribLength=bLen
	while True:
		if entry2[0]==5: # case entry1 is N
			pos=min(entry1[1],entry2[1])
			pass
		elif entry1[0]==5: # case entry2 is N
			#if parent node is type "N", in theory we might have to calculate the contribution of root nucleotides; 
			# however, if this node is "N" then every other node in the current tree is "N", so we can ignore this since this contribution cancels out in relative terms.
			pos=min(entry1[1],entry2[1])
			pass
		else:
			#contribLength will be here the total length from the root or from the upper node, down to the down node.
			if entry1[0]<5:
				if len(entry1)==2:
					contribLength=bLen
				elif len(entry1)==3:
					contribLength=entry1[2]
					if bLen:
						contribLength+=bLen
				else:
					contribLength=entry1[3]
					if bLen:
						contribLength+=bLen
			else:
				if len(entry1)==3:
					contribLength=bLen
				else:
					contribLength=entry1[2]
					if bLen:
						contribLength+=bLen
			if entry2[0]<5:
				if len(entry2)==3:
					if contribLength:
						contribLength+=entry2[2]
					else:
						contribLength=entry2[2]
			else:
				if len(entry2)==4:
					if contribLength:
						contribLength+=entry2[2]
					else:
						contribLength=entry2[2]

			if entry1[0]==4: # case entry1 is R	
				if entry2[0]==4:
					if len(entry1)==4:
						end=min(entry1[1],entry2[1])
						contribLength+=entry1[2]
						Lkcost+=contribLength*(cumulativeRate[end]-cumulativeRate[pos])
						pos=end
					else:
						if contribLength:
							end=min(entry1[1],entry2[1])
							Lkcost+=contribLength*(cumulativeRate[end]-cumulativeRate[pos])
							pos=end
						else:
							pos=min(entry1[1],entry2[1])

				#entry1 is reference and entry2 is of type "O"
				elif entry2[0]==6:
					i1=refIndeces[pos]
					if len(entry1)==4:
						tot=0.0
						for i in range4:
							if i1==i:
								tot2=rootFreqs[i]*(1.0+nonMutRates[i]*entry1[2])
							else:
								tot2=rootFreqs[i]*mutMatrix[i][i1]*entry1[2]
							if contribLength:
								tot3=0.0
								for j in range4:
									tot3+=mutMatrix[i][j]*entry2[-1][j]
								tot+=tot2*(entry2[-1][i] + contribLength*tot3)
							else:
								tot+=tot2*entry2[-1][i]
					else:
						if contribLength:
							tot=0.0
							for j in range4:
								tot+=mutMatrix[i1][j]*entry2[-1][j]
							tot*=contribLength
							tot+=entry2[-1][i1]
						else:
							tot=entry2[-1][i1]
					totalFactor*=tot
					pos+=1

				else: #entry1 is R and entry2 is a different but single nucleotide
					if len(entry1)==4:
						i1=refIndeces[pos]
						i2=entry2[0]
						if contribLength:
							totalFactor*=((rootFreqs[i1]*mutMatrix[i1][i2]*contribLength*(1.0+nonMutRates[i1]*entry1[2])+rootFreqs[i2]*mutMatrix[i2][i1]*entry1[2]*(1.0+nonMutRates[i2]*contribLength))/rootFreqs[i1])
						else:
							totalFactor*=((rootFreqs[i2]*mutMatrix[i2][i1]*entry1[2])/rootFreqs[i1])
					else:
						if contribLength:
							totalFactor*=mutMatrix[refIndeces[pos]][entry2[0]]*contribLength
						else:
							return float("-inf")
					pos+=1

			# entry1 is of type "O"
			elif entry1[0]==6:
				if entry2[0]==6:
					if contribLength:
						tot=0.0
						for j in range4:
							tot+=entry1[-1][j]*(entry2[-1][j] + contribLength*(mutMatrix[j][0]*entry2[-1][0]+mutMatrix[j][1]*entry2[-1][1]+mutMatrix[j][2]*entry2[-1][2]+mutMatrix[j][3]*entry2[-1][3]) )
						totalFactor*=tot
					else:
						tot=0.0
						for j in range4:
							tot+=entry1[-1][j]*entry2[-1][j]
						totalFactor*=tot
				else: #entry1 is "O" and entry2 is a nucleotide
					if entry2[0]==4:
						i2=refIndeces[pos]
					else:
						i2=entry2[0]
					if contribLength:
						totalFactor*=(entry1[-1][i2]+contribLength*(entry1[-1][0]*mutMatrix[0][i2]+entry1[-1][1]*mutMatrix[1][i2]+entry1[-1][2]*mutMatrix[2][i2]+entry1[-1][3]*mutMatrix[3][i2]))
					else:
						totalFactor*=entry1[-1][i2]
				pos+=1

			else: #entry1 is a non-ref nuc
				if entry2[0]==entry1[0]:
					if len(entry1)==4:
						contribLength+=entry1[2]
					if contribLength:
						Lkcost+=nonMutRates[entry1[0]]*contribLength
				else: #entry1 is a nucleotide and entry2 is not the same as entry1
					i1=entry1[0]
					if entry2[0]<5: #entry2 is a nucleotide
						if entry2[0]==4:
							i2=refIndeces[pos]
						else:
							i2=entry2[0]
						if len(entry1)==4:
							if contribLength:
								totalFactor*=(( rootFreqs[i1]*mutMatrix[i1][i2]*contribLength*(1.0+nonMutRates[i1]*entry1[2]) + rootFreqs[i2]*mutMatrix[i2][i1]*entry1[2]*(1.0+nonMutRates[i2]*contribLength) )/rootFreqs[i1])
							else:
								totalFactor*=((rootFreqs[i2]*mutMatrix[i2][i1]*entry1[2] )/rootFreqs[i1])
						else:
							if contribLength:
								totalFactor*=mutMatrix[i1][i2]*contribLength
							else:
								return float("-inf")
					else: #entry1 is a nucleotide and entry2 is of type "O"
						if len(entry1)==4:
							tot=0.0
							for i in range4:
								if i1==i:
									tot2=rootFreqs[i]*(1.0+nonMutRates[i]*entry1[2])
								else:
									tot2=rootFreqs[i]*mutMatrix[i][i1]*entry1[2]
								tot3=0.0
								for j in range4:
									tot3+=mutMatrix[i][j]*entry2[-1][j]
								tot+=tot2*(entry2[-1][i] + contribLength*tot3)
							totalFactor*=(tot/rootFreqs[i1])
						else:
							tot=0.0
							for j in range4:
								tot+=mutMatrix[i1][j]*entry2[-1][j]
							tot*=contribLength
							tot+=entry2[-1][i1]
							totalFactor*=tot
				pos+=1

		if totalFactor<=minimumCarryOver:
			if totalFactor<sys.float_info.min:
				return float("-inf")
			Lkcost+=log(totalFactor)
			totalFactor=1.0
		if pos==lRef:
			break
		if pos==entry1[1]:
			indexEntry1+=1
			entry1=probVectP[indexEntry1]
		if pos==entry2[1]:
			indexEntry2+=1
			entry2=probVectC[indexEntry2]
	return Lkcost+log(totalFactor)









#function traversing the tree to find the best node in the tree where to re-append the given subtree (rooted at node.children[child]) to improve the topology of the current tree.
# bestLKdiff is the best likelihood cost found for the current placement (after optimizing the branch length).
# removedBLen is such branch length that optimizes the current placement - it will be used to place the subtree attached at other nodes of the tree.
def findBestParentTopology(node,child,bestLKdiff,removedBLen,mutMatrix):
	bestNodeSoFar=node
	bestIsMidNode=False
	nodesToVisit=[]
	removedPartials=node.children[child].probVect
	if node.up!=None:
		if node.up.children[0]==node:
			childUp=1
			vectUpUp=node.up.probVectUpRight
		else:
			childUp=2
			vectUpUp=node.up.probVectUpLeft
		#the list nodesToVisit keeps trak of the nodes of the tree we wil need to traverse, (first element of each entry),
		# the direction from where we are visitng them (0=from parent, and 1,2=from one of the children),
		# an updated genome list from the direction where we come from (taking into account the removal of the given subtree),
		# a branch length value separating the node from this updated genome list (useful for the fact that the removal of the subtree changes the branch length at the removal node),
		# a flag that says if the updated genome list passed needs still updating, or if it has become identical to the pre-existing genome list in the tree (which usually happens after a while),
		# the likelihood cost of appending at the last node encountered in this direction,
		# a number of consecutively falied traversal steps since the last improvement found (if this number goes beyond a threshold, traversal in the considered direction might be halted).
		nodesToVisit.append((node.up,childUp,node.children[1-child].probVect,node.children[1-child].dist+node.dist,True,bestLKdiff,0))
		nodesToVisit.append((node.children[1-child],0,vectUpUp,node.children[1-child].dist+node.dist,True,bestLKdiff,0))
	else:
		# case node is root
		if node.children[1-child].children: # case there is only one sample outside of the subtree doesn't need to be considered
			child1=node.children[1-child].children[0]
			child2=node.children[1-child].children[1]
			vectUp1=rootVector(child2.probVect,child2.dist,mutMatrix)
			nodesToVisit.append((child1,0,vectUp1,child1.dist,True,bestLKdiff,0))
			vectUp2=rootVector(child1.probVect,child1.dist,mutMatrix)
			nodesToVisit.append((child2,0,vectUp2,child2.dist,True,bestLKdiff,0))

	while nodesToVisit:
		t1,direction,passedPartials,distance,needsUpdating,lastLK,failedPasses=nodesToVisit.pop()
		if direction==0:
			#consider the case we are moving from a parent to a child
			if t1.dist:
				if not (t1.up==node or t1.up==None): #try to append mid-branch
					if needsUpdating:
						midTot=mergeVectorsUpDown(passedPartials,distance/2,t1.probVect,distance/2,mutMatrix)
					else:
						midTot=t1.probVectTotUp
					midProb=appendProbNode(midTot,removedPartials,removedBLen,mutMatrix)
					if midProb>bestLKdiff:
						bestNodeSoFar=t1
						bestLKdiff=midProb
						bestIsMidNode=True
						failedPasses=0
				#now try appending exactly at node
				if needsUpdating:
					nodeTot=mergeVectorsUpDown(passedPartials,distance,t1.probVect,False,mutMatrix)
					if not areVectorsDifferent(nodeTot,t1.probVectTot):
						needsUpdating=False
				else:
					nodeTot=t1.probVectTot
				nodeProb=appendProbNode(nodeTot,removedPartials,removedBLen,mutMatrix)
				if nodeProb<(lastLK-thresholdLogLKconsecutivePlacement): #placement at current node is considered failed if placement likelihood is not improved by a certain margin compared to best placement so far for the nodes above it.
					failedPasses+=1
				elif nodeProb>bestLKdiff:
					bestNodeSoFar=t1
					bestLKdiff=nodeProb
					bestIsMidNode=False
					failedPasses=0
				#TODO check if this commented out below is any useful?
				#elif nodeProb>lastLK+1.0 and failedPasses>0:
				#	failedPasses-=1
			else:
				nodeProb=lastLK
			
			#keep crawling down into children nodes unless the stop criteria for the traversal are satisfied.
			traverseChildren=False
			if strictTopologyStopRules:
				if failedPasses<=allowedFailsTopology and nodeProb>(bestLKdiff-thresholdLogLKtopology) and t1.children:
					traverseChildren=True
			else:
				if failedPasses<=allowedFailsTopology or nodeProb>(bestLKdiff-thresholdLogLKtopology):
					if t1.children:
						traverseChildren=True
			if traverseChildren:
				#if (failedPasses<=allowedFailsTopology or nodeProb>(bestLKdiff-thresholdLogLKtopology) ) and len(t1.children)==2:
				child=t1.children[0]
				otherChild=t1.children[1]
				if needsUpdating:
					vectUpRight=mergeVectorsUpDown(passedPartials,distance,otherChild.probVect,otherChild.dist,mutMatrix)
				else:
					vectUpRight=t1.probVectUpRight
				if vectUpRight!=None:
					nodesToVisit.append((child,0,vectUpRight,child.dist,needsUpdating,nodeProb,failedPasses))
				child=t1.children[1]
				otherChild=t1.children[0]
				if needsUpdating:
					vectUpLeft=mergeVectorsUpDown(passedPartials,distance,otherChild.probVect,otherChild.dist,mutMatrix)
				else:
					vectUpLeft=t1.probVectUpLeft
				if vectUpLeft!=None:
					nodesToVisit.append((child,0,vectUpLeft,child.dist,needsUpdating,nodeProb,failedPasses))

		else: #case when crawling up from child to parent
			if t1.dist or t1.up==None: #append directly at the node
				if needsUpdating:
					if direction==1:
						nodeTot=mergeVectorsUpDown(t1.probVectUpRight,False,passedPartials,distance,mutMatrix)
					else:
						nodeTot=mergeVectorsUpDown(t1.probVectUpLeft,False,passedPartials,distance,mutMatrix)
					if nodeTot==None:
						#print("Removing a node created an inconsistency while moving up.")
						continue
					elif not areVectorsDifferent(nodeTot,t1.probVectTot):
						needsUpdating=False
				else:
					nodeTot=t1.probVectTot
				nodeProb=appendProbNode(nodeTot,removedPartials,removedBLen,mutMatrix)
				if nodeProb<(lastLK-thresholdLogLKconsecutivePlacement): #placement at current node is considered failed if placement likelihood is not improved by a certain margin compared to best placement so far for the nodes above it.
					failedPasses+=1
				elif nodeProb>bestLKdiff:
					bestNodeSoFar=t1
					bestLKdiff=nodeProb
					bestIsMidNode=False
					failedPasses=0
				#TODO check if this commented out below is any useful?
				#elif nodeProb>lastLK+1.0 and failedPasses>0:
				#	failedPasses-=1
			else:
				nodeProb=lastLK

			otherChild=t1.children[2-direction]
			midBottom=None
			if t1.dist and t1.up!=None: #try appending mid-branch
				if needsUpdating:
					midBottom=mergeVectors(otherChild.probVect,otherChild.dist,passedPartials,distance,mutMatrix)
					if midBottom==None:
						continue
					if t1==t1.up.children[0]:
						vectUp=t1.up.probVectUpRight
					else:
						vectUp=t1.up.probVectUpLeft
					midTot=mergeVectorsUpDown(vectUp,t1.dist/2,midBottom,t1.dist/2,mutMatrix)
				else:
					midTot=t1.probVectTotUp
				midProb=appendProbNode(midTot,removedPartials,removedBLen,mutMatrix)
				if midProb>bestLKdiff:
					bestNodeSoFar=t1
					bestLKdiff=midProb
					bestIsMidNode=True
					failedPasses=0
			
			#testing for the stop rule of the traversal process
			keepTraversing=False
			if strictTopologyStopRules:
				if failedPasses<=allowedFailsTopology and nodeProb>(bestLKdiff-thresholdLogLKtopology):
					keepTraversing=True
			else:
				if failedPasses<=allowedFailsTopology or nodeProb>(bestLKdiff-thresholdLogLKtopology):
					keepTraversing=True
			if keepTraversing:
				#if failedPasses<=allowedFailsTopology or nodeProb>(bestLKdiff-thresholdLogLKtopology) :
				# keep crawling up into parent and sibling node
				if t1.up!=None: #case the node is not the root
					#first pass the crawling down the other child
					if t1==t1.up.children[0]:
						upChild=0
						if needsUpdating:
							vectUpUp=t1.up.probVectUpRight
					else:
						upChild=1
						if needsUpdating:
							vectUpUp=t1.up.probVectUpLeft
					if needsUpdating:
						vectUp=mergeVectorsUpDown(vectUpUp,t1.dist,passedPartials,distance,mutMatrix)
					else:
						if direction==1:
							vectUp=t1.probVectUpLeft
						else:
							vectUp=t1.probVectUpRight

					if vectUp==None:
						continue
					else:
						nodesToVisit.append((otherChild,0,vectUp,otherChild.dist,needsUpdating,nodeProb,failedPasses))
					#now pass the crawling up to the parent node
					if needsUpdating:
						if midBottom==None:
							midBottom=mergeVectors(otherChild.probVect,otherChild.dist,passedPartials,distance,mutMatrix)
							if midBottom==None:
								continue
					else:
						midBottom=t1.probVect
					nodesToVisit.append((t1.up,upChild+1,midBottom,t1.dist,needsUpdating,nodeProb,failedPasses))
				#now consider case of root node
				else:
					if needsUpdating:
						vectUp=rootVector(passedPartials,distance,mutMatrix)
					else:
						if direction==1:
							vectUp=t1.probVectUpLeft
						else:
							vectUp=t1.probVectUpRight
					nodesToVisit.append((otherChild,0,vectUp,otherChild.dist,needsUpdating,nodeProb,failedPasses))
	return bestNodeSoFar , bestLKdiff , bestIsMidNode









#we know that sample "appendedNode", with partials "newPartials", is best placed as child of "node" resulting in logLK contribution of newChildLK
# now explore exactly which position near node is best for placement (direct child, sibling or sibling of existing child), and which are the best branch lengths,
#then add the sample at that position of the tree, and update all the internal probability vectors.
def placeSampleOnTreeTopology(node,newPartials,appendedNode,newBranchL,newChildLK,isMidNode,mutMatrix):
	if not node: #in case of a polytomy, reach first the top of the polytomy, which is the only node at which appending is allowed.
		while (not node.dist) and node.up!=None:
			node=node.up
	if isMidNode and node.up!=None: #best appending location found so far was mid-branch above node
		if node==node.up.children[0]:
			child=0
			vectUp=node.up.probVectUpRight
		else:
			child=1
			vectUp=node.up.probVectUpLeft
		bestSplit=0.5
		bestSplitLK=newChildLK
		childBestVect=node.probVectTotUp
		newSplit=0.25
		#print("Trying to append mid-node with LK "+str(bestSplitLK)+" and bL "+str(newBranchL)+" at node with current length "+str(node.dist))
		#try different positions on the existing branch
		while newSplit*node.dist>minBLen:
			probVectParentNew=mergeVectorsUpDown(vectUp,node.dist*newSplit,node.probVect,node.dist*(1.0-newSplit),mutMatrix)
			probChild=appendProbNode(probVectParentNew,newPartials,newBranchL,mutMatrix)
			if probChild>bestSplitLK:
				bestSplitLK=probChild
				bestSplit=newSplit
				childBestVect=probVectParentNew
			else:
				break
			newSplit=bestSplit/2
		if bestSplit>0.49:
			newSplit=0.25
			#print("Now trying the reverse direction")
			while newSplit*node.dist>minBLen:
				probVectParentNew=mergeVectorsUpDown(vectUp,node.dist*(1.0-newSplit),node.probVect,node.dist*newSplit,mutMatrix)
				probChild=appendProbNode(probVectParentNew,newPartials,newBranchL,mutMatrix)
				if probChild>bestSplitLK:
					bestSplitLK=probChild
					bestSplit=newSplit
					childBestVect=probVectParentNew
				else:
					bestSplit=1.0-bestSplit
					break
				newSplit=bestSplit/2
		#now try different lengths for the new branch
		LK1=bestSplitLK
		bestLen=newBranchL
		#print("Initial bLen: "+str(bestLen)+" and LK: "+str(LK1))
		if not bestLen:
			#print("weird, best is mid-branch node, but node is part of a polytomy.")
			bestLen=minBLen
			LK1=appendProbNode(childBestVect,newPartials,bestLen,mutMatrix)
		while bestLen>minBLen:
			newBLen=bestLen/2
			probChild=appendProbNode(childBestVect,newPartials,newBLen,mutMatrix)
			if probChild>LK1:
				LK1=probChild
				bestLen=newBLen
			else:
				break
		if (not newBranchL) or bestLen>0.7*newBranchL:
			while bestLen<maxBLen:
				newBLen=bestLen*2
				probChild=appendProbNode(childBestVect,newPartials,newBLen,mutMatrix)
				if probChild>LK1:
					LK1=probChild
					bestLen=newBLen
				else:
					break
		if bestLen<2*minBLen:
			LK0=appendProbNode(childBestVect,newPartials,False,mutMatrix)
			if LK0>LK1:
				bestLen=False

		#now create new internal node and append child to it
		distTop=node.dist*bestSplit
		newInternalNode=Tree()
		newInternalNode.dirty=True
		node.up.children[child]=newInternalNode
		newInternalNode.up=node.up
		distBottom=node.dist*(1.0-bestSplit)
		newInternalNode.add_child(node)
		node.up=newInternalNode
		node.dist=distBottom
		appendedNode.dist=bestLen
		appendedNode.up=newInternalNode
		newInternalNode.add_child(appendedNode)
		newInternalNode.dist=distTop
		newInternalNode.probVectUpRight=mergeVectorsUpDown(vectUp,newInternalNode.dist,newPartials,bestLen,mutMatrix)
		newInternalNode.probVectUpLeft=childBestVect
		newInternalNode.probVect=mergeVectors(node.probVect,node.dist,newPartials,bestLen,mutMatrix)
		newInternalNode.probVectTotUp=mergeVectorsUpDown(vectUp,distTop/2,newInternalNode.probVect,distTop/2,mutMatrix)
		newVect=mergeVectorsUpDown(childBestVect,False,newPartials,bestLen,mutMatrix)
		if newVect==None:
			print("Problem, None vector when re-placing sample, mid-branch")
			print(newVect)
			print(childBestVect)
			print(newPartials)
			print(bestLen)
		newInternalNode.probVectTot=newVect
		if verbose:
			print("new internal node added to tree")
			print(newInternalNode.probVect)
			print(newInternalNode.probVectUpRight)
			print(newInternalNode.probVectUpLeft)
			print(newInternalNode.probVectTot)
		if distTop>=2*minBLenForMidNode:
			createFurtherMidNodes(newInternalNode,vectUp)
		updatePartialsFromTop(node,newInternalNode.probVectUpRight,mutMatrix)
		updatePartialsFromTop(appendedNode,newInternalNode.probVectUpLeft,mutMatrix)
		updatePartialsFromBottom(newInternalNode.up,newInternalNode.probVect,child,newInternalNode,mutMatrix)

	#best lk so far is for appending directly to existing node
	else:
		#We first explore placement just below the best placement node for more fine-grained placement within its descendant branches (accounting for polytomies).
		bestDownLK=float("-inf")
		bestDownNode=None
		#current node might be part of a polytomy (represented by 0 branch lengths) so we want to explore all the children of the current node
		#to find out if the best placement is actually in any of the branches below the current node.
		nodesToVisit=[]
		bestSplit=0.5
		for c in node.children:
			nodesToVisit.append(c)
		while nodesToVisit:
			t1=nodesToVisit.pop()
			if not t1.dist:
				for c in t1.children:
					nodesToVisit.append(c)
			else:
				newSplit=0.5
				newBestSplit=0.5
				#now try to place on the current branch below the best node, at an height above or equal to the mid-branch.
				newBLen2=t1.dist*newSplit
				bestLKdiff2=float("-inf")
				furtherNode=-1
				newProbVect2=t1.probVectTotUp
				while True:
					newLKdiff2=appendProbNode(newProbVect2,newPartials,newBranchL,mutMatrix)
					if newLKdiff2>bestLKdiff2:
						bestLKdiff2=newLKdiff2
						newBestSplit=newSplit
					else:
						break
					newSplit=newSplit/2
					newBLen2=t1.dist*newSplit
					if newBLen2<=minBLenForMidNode/2:
						break
					furtherNode+=1
					newProbVect2=t1.furtherMidNodes[furtherNode]
				
				if bestLKdiff2>bestDownLK:
					bestDownLK=bestLKdiff2
					bestDownNode=t1
					bestSplit=newBestSplit

		if bestDownNode!=None: #now that the best descendant branch is found, find more precisely at which height within the branch the best placement is.
			if bestDownNode==bestDownNode.up.children[0]:
				child=0
				vectUp=bestDownNode.up.probVectUpRight
			else:
				child=1
				vectUp=bestDownNode.up.probVectUpLeft
			bestSplitLK=bestDownLK
			childBestVect=bestDownNode.probVectTotUp
			newSplit=bestSplit/2
			while newSplit*bestDownNode.dist>minBLen:
				probVectParentNew=mergeVectorsUpDown(vectUp,bestDownNode.dist*newSplit,bestDownNode.probVect,bestDownNode.dist*(1.0-newSplit),mutMatrix)
				probChild=appendProbNode(probVectParentNew,newPartials,newBranchL,mutMatrix)
				if probChild>bestSplitLK:
					bestSplitLK=probChild
					bestSplit=newSplit
					childBestVect=probVectParentNew
				else:
					break
				newSplit=bestSplit/2
			bestChildLK=bestSplitLK
			bestChildSplit=bestSplit
		else:
			bestChildLK=float("-inf")

		#if node is root, try to place as sibling of the current root.
		if node.up==None:
			probOldRoot = findProbRoot(node.probVect)
			probVectRoot,probRoot = mergeVectors(node.probVect,oneMutBLen,newPartials,newBranchL,mutMatrix,returnLK=True)
			probRoot+= findProbRoot(probVectRoot)
			parentLKdiff=probRoot-probOldRoot
			bestRootBL=oneMutBLen
			parentBestVect=probVectRoot
			newBL=0.5*oneMutBLen
			while newBL>minBLen:
				probVectRoot,probRoot = mergeVectors(node.probVect,newBL,newPartials,newBranchL,mutMatrix,returnLK=True)
				probRoot+= findProbRoot(probVectRoot)
				newDiff=probRoot-probOldRoot
				if newDiff>parentLKdiff:
					parentLKdiff=newDiff
					bestRootBL=newBL
					parentBestVect=probVectRoot
				else:
					break
				newBL=bestRootBL/2

		else: #try to append just above node
			if node==node.up.children[0]:
				child=0
				vectUp=node.up.probVectUpRight
			else:
				child=1
				vectUp=node.up.probVectUpLeft
			bestSplit=0.5
			parentBestVect=node.probVectTotUp
			bestSplitLK=appendProbNode(parentBestVect,newPartials,newBranchL,mutMatrix)
			newSplit=0.25
			while newSplit*node.dist>minBLen:
				probVectParentNew=mergeVectorsUpDown(vectUp,node.dist*(1.0-newSplit),node.probVect,node.dist*newSplit,mutMatrix)
				probChild=appendProbNode(probVectParentNew,newPartials,newBranchL,mutMatrix)
				if probChild>bestSplitLK:
					bestSplitLK=probChild
					bestSplit=newSplit
					parentBestVect=probVectParentNew
				else:
					break
				newSplit=bestSplit/2
			parentLKdiff=bestSplitLK
			bestParentSplit=bestSplit
		#now we have three likelihood costs,
		# bestChildLK is the likelihood score of appending below node;
		# parentLKdiff is the likelihood score of appending above node;
		# newChildLK is the likelihood cost of appending exactly at node.
		if bestChildLK>=parentLKdiff and bestChildLK>=newChildLK: #Best placement is below node: add internal node below "node"
			if bestDownNode==bestDownNode.up.children[0]:
				child=0
				vectUp=bestDownNode.up.probVectUpRight
			else:
				child=1
				vectUp=bestDownNode.up.probVectUpLeft

			LK1=bestChildLK
			bestLen=newBranchL
			if not bestLen:
				bestLen=minBLen
				LK1=appendProbNode(childBestVect,newPartials,bestLen,mutMatrix)
			while bestLen>minBLen:
				newBLen=bestLen/2
				probChild=appendProbNode(childBestVect,newPartials,newBLen,mutMatrix)
				if probChild>LK1:
					LK1=probChild
					bestLen=newBLen
				else:
					break
			if (not newBranchL) or bestLen>0.7*newBranchL:
				while bestLen<maxBLen:
					newBLen=bestLen*2
					probChild=appendProbNode(childBestVect,newPartials,newBLen,mutMatrix)
					if probChild>LK1:
						LK1=probChild
						bestLen=newBLen
					else:
						break
			if bestLen<2*minBLen:
				LK0=appendProbNode(childBestVect,newPartials,False,mutMatrix)
				if LK0>LK1:
					bestLen=False
			#now create new internal node and append child to it
			newInternalNode=Tree()
			newInternalNode.dirty=True
			bestDownNode.up.children[child]=newInternalNode
			newInternalNode.up=bestDownNode.up
			distBottom=bestDownNode.dist*(1.0-bestChildSplit)
			distTop=bestDownNode.dist*bestChildSplit
			bestDownNode.up=newInternalNode
			bestDownNode.dist=distBottom
			newInternalNode.add_child(bestDownNode)
			appendedNode.dist=bestLen
			appendedNode.up=newInternalNode
			newInternalNode.add_child(appendedNode)
			newInternalNode.dist=distTop
			newInternalNode.probVectUpRight=mergeVectorsUpDown(vectUp,newInternalNode.dist,newPartials,bestLen,mutMatrix)
			newInternalNode.probVectUpLeft=childBestVect
			newInternalNode.probVect=mergeVectors(bestDownNode.probVect,bestDownNode.dist,newPartials,bestLen,mutMatrix)
			newInternalNode.probVectTotUp=mergeVectorsUpDown(vectUp,distTop/2,newInternalNode.probVect,distTop/2,mutMatrix)
			newInternalNode.probVectTot=mergeVectorsUpDown(childBestVect,False,newPartials,bestLen,mutMatrix)
			if newInternalNode.probVectTot==None:
				print("Problem, None vector when re-placing sample, position 2")
				print(newInternalNode.probVectTot)
				print(childBestVect)
				print(newPartials)
				print(bestLen)
			if distTop>=2*minBLenForMidNode:
				createFurtherMidNodes(newInternalNode,vectUp)
			if verbose:
				print("new internal node added to tree")
				print(newInternalNode.probVect)
				print(newInternalNode.probVectUpRight)
				print(newInternalNode.probVectUpLeft)
				print(newInternalNode.probVectTot)
			updatePartialsFromTop(bestDownNode,newInternalNode.probVectUpRight,mutMatrix)
			updatePartialsFromTop(appendedNode,newInternalNode.probVectUpLeft,mutMatrix)
			updatePartialsFromBottom(newInternalNode.up,newInternalNode.probVect,child,newInternalNode,mutMatrix)
		
		else: # the placement will be exactly at node, by creating (or extending) a polytomy, or above node.
			if newChildLK>=parentLKdiff: #new parent will be part of a polytomy since best placement is exactly at the node
				bestRootBL=False
				bestParentSplit=False
				parentLKdiff=newChildLK
				parentBestVect=node.probVectTot
				if node.up==None:
					probVectRoot,probRoot = mergeVectors(node.probVect,False,newPartials,newBranchL,mutMatrix,returnLK=True)
					parentBestVect=probVectRoot

			#add parent to the root
			if node.up==None:
				#now try different lengths for right branch
				bestLen2=newBranchL
				if not bestLen2:
					bestLen2=minBLen
					newProbVectRoot,newProbRoot = mergeVectors(node.probVect,bestRootBL,newPartials,bestLen2,mutMatrix,returnLK=True)
					newProbRoot+= findProbRoot(newProbVectRoot)
					parentLKdiff=newProbRoot-probOldRoot
				while bestLen2>minBLen:
					newBLen=bestLen2/2
					newProbVectRoot,newProbRoot = mergeVectors(node.probVect,bestRootBL,newPartials,newBLen,mutMatrix,returnLK=True)
					newProbRoot+= findProbRoot(newProbVectRoot)
					LKdiffRoot=newProbRoot-probOldRoot
					if LKdiffRoot>parentLKdiff:
						parentLKdiff=LKdiffRoot
						bestLen2=newBLen
						probVectRoot=newProbVectRoot
					else:
						break
				if (not newBranchL) or bestLen2>0.7*newBranchL:
					while bestLen2<maxBLen:
						newBLen=bestLen2*2
						newProbVectRoot,newProbRoot = mergeVectors(node.probVect,bestRootBL,newPartials,newBLen,mutMatrix,returnLK=True)
						newProbRoot+= findProbRoot(newProbVectRoot)
						LKdiffRoot=newProbRoot-probOldRoot
						if LKdiffRoot>parentLKdiff:
							parentLKdiff=LKdiffRoot
							bestLen2=newBLen
							probVectRoot=newProbVectRoot
						else:
							break
				if bestLen2<2*minBLen:
					newProbVectRoot,newProbRoot = mergeVectors(node.probVect,bestRootBL,newPartials,False,mutMatrix,returnLK=True)
					newProbRoot+= findProbRoot(newProbVectRoot)
					LK0=newProbRoot-probOldRoot
					if LK0>parentLKdiff:
						bestLen2=False
						parentLKdiff=LK0
						probVectRoot=newProbVectRoot

				newRoot=Tree()
				newRoot.dirty=True
				newRoot.probVect=probVectRoot
				newRoot.probVectTot=rootVector(probVectRoot,False,mutMatrix)
				newRoot.probVectUpRight=rootVector(newPartials,bestLen2,mutMatrix)
				newRoot.probVectUpLeft=rootVector(node.probVect,bestRootBL,mutMatrix)
				if newRoot.probVectTot==None:
					print("Problem, None vector when re-placing sample, position root")
					print(newRoot.probVectTot)
					print(newPartials)
					print(bestLen2)
					print(bestRootBL)
					print(node.probVect)
				node.up=newRoot
				node.dist=bestRootBL
				newRoot.add_child(node)
				appendedNode.up=newRoot
				newRoot.add_child(appendedNode)
				appendedNode.dist=bestLen2
				if verbose:
					print("new root added to tree")
					print(newRoot.probVect)
					print(newRoot.children[0].probVect)
					print(newRoot.children[1].probVect)
				updatePartialsFromTop(node,newRoot.probVectUpRight,mutMatrix)
				updatePartialsFromTop(appendedNode,newRoot.probVectUpLeft,mutMatrix)
				return newRoot

			#add parent to non-root node
			else:
				if node==node.up.children[0]:
					child=0
					vectUp=node.up.probVectUpRight
				else:
					child=1
					vectUp=node.up.probVectUpLeft
				
				#now try different lengths for the new branch
				LK1=parentLKdiff
				bestLen=newBranchL
				if not bestLen:
					bestLen=minBLen
					LK1=appendProbNode(parentBestVect,newPartials,bestLen,mutMatrix)
				while bestLen>minBLen:
					newBLen=bestLen/2
					probChild=appendProbNode(parentBestVect,newPartials,newBLen,mutMatrix)
					if probChild>LK1:
						LK1=probChild
						bestLen=newBLen
					else:
						break
				if (not newBranchL) or bestLen>0.7*newBranchL:
					while bestLen<10*newBranchL:
						newBLen=bestLen*2
						probChild=appendProbNode(parentBestVect,newPartials,newBLen,mutMatrix)
						if probChild>LK1:
							LK1=probChild
							bestLen=newBLen
						else:
							break
				#print("final branch length above node "+str(bestLen)+" with LK "+str(LK1))
				if bestLen<2*minBLen:
					LK0=appendProbNode(parentBestVect,newPartials,False,mutMatrix)
					if LK0>LK1:
						bestLen=False
				newInternalNode=Tree()
				newInternalNode.dirty=True
				node.up.children[child]=newInternalNode
				newInternalNode.up=node.up
				if bestParentSplit:
					distBottom=node.dist*bestParentSplit
					distTop=node.dist*(1.0-bestParentSplit)
				else:
					distBottom=False
					distTop=node.dist
					node.probVectTot=None
					node.probVectTotUp=None
					node.furtherMidNodes=None
				node.dist=distBottom
				node.up=newInternalNode
				newInternalNode.add_child(node)
				appendedNode.up=newInternalNode
				newInternalNode.add_child(appendedNode)
				appendedNode.dist=bestLen
				newInternalNode.dist=distTop
				newInternalNode.probVectUpLeft=parentBestVect
				newVect=mergeVectors(node.probVect,node.dist,newPartials,bestLen,mutMatrix)
				if newVect==None:
					bestLen=minBLen
					appendedNode.dist=bestLen
					newVect=mergeVectors(node.probVect,node.dist,newPartials,bestLen,mutMatrix)
				newInternalNode.probVect=newVect
				newInternalNode.probVectUpRight=mergeVectorsUpDown(vectUp,newInternalNode.dist,newPartials,bestLen,mutMatrix)
				newInternalNode.probVectTotUp=mergeVectorsUpDown(vectUp,distTop/2,newInternalNode.probVect,distTop/2,mutMatrix)
				newInternalNode.probVectTot=mergeVectorsUpDown(parentBestVect,False,newPartials,bestLen,mutMatrix)
				if newInternalNode.probVectTot==None:
					print("Problem, None vector when re-placing sample, position 3")
					print(newPartials)
					print(bestLen)
					print(parentBestVect)
					print(vectUp)
					print(distTop)
					print(distBottom)
				if distTop>=2*minBLenForMidNode:
					createFurtherMidNodes(newInternalNode,vectUp)
				if verbose:
					print("new internal node added to tree")
					print(newInternalNode.probVect)
					print(newInternalNode.probVectUpRight)
					print(newInternalNode.probVectUpLeft)
					print(newInternalNode.probVectTot)
				updatePartialsFromTop(node,newInternalNode.probVectUpRight,mutMatrix)
				updatePartialsFromTop(appendedNode,newInternalNode.probVectUpLeft,mutMatrix)
				updatePartialsFromBottom(newInternalNode.up,newInternalNode.probVect,child,newInternalNode,mutMatrix)

	return None







#remove node from the current position in the tree and re-attach it at a new given place new bestNode.
# First remove node from the tree, then update the genome lists;
# then find the exact best reattachment of node and update the genome lists again using function placeSampleOnTreeTopology().
def cutAndPasteNode(node,bestNode,isMidNode,attachmentBLen,bestLK,mutMatrix):
	#remove node from the tree
	parentNode=node.up
	if node==parentNode.children[0]:
		sibling=parentNode.children[1]
	else:
		sibling=parentNode.children[0]
	if parentNode.up!=None:
		if parentNode==parentNode.up.children[0]:
			childP=0
			vectUp=parentNode.up.probVectUpRight
		else:
			childP=1
			vectUp=parentNode.up.probVectUpLeft
		parentNode.up.children[childP]=sibling
	sibling.up=parentNode.up
	if sibling.dist:
		if parentNode.dist:
			sibling.dist+=parentNode.dist
	else:
		sibling.dist=parentNode.dist
	#update likelihood lists after node removal
	if sibling.up==None:
		sibling.probVectTot=rootVector(sibling.probVect,False,mutMatrix)
		if sibling.children:
			sibling.probVectUpRight=rootVector(sibling.children[1].probVect,sibling.children[1].dist,mutMatrix)
			sibling.probVectUpLeft=rootVector(sibling.children[0].probVect,sibling.children[0].dist,mutMatrix)
			updatePartialsFromTop(sibling.children[0],sibling.probVectUpRight,mutMatrix)
			updatePartialsFromTop(sibling.children[1],sibling.probVectUpLeft,mutMatrix)
	else:
		updatePartialsFromTop(sibling,vectUp,mutMatrix)
		updatePartialsFromBottom(sibling.up,sibling.probVect,childP,sibling,mutMatrix)
	#re-place the node and re-update the vector lists
	newRoot = placeSampleOnTreeTopology(bestNode,node.probVect,node,attachmentBLen,bestLK,isMidNode, mutMatrix)
	#if the root of the tree has changed, return the new root
	if sibling.up==None:
		return sibling
	else:
		return newRoot


# try to find a re-placement of a dirty node of the tree to improve the topology.
# Cut out the subtree at this node, and look for somwhere else in the tree where to attach it (an SPR move).
# To find the best location of the new re-attachment, we traverse the tree starting at the current attachment node, and for each node visited we evaluate the reattachment,
# as done by the findBestParentTopology() function.
# To avoid traversing the whole tree for each SPR move, we use stopping conditions similar to those used in the initial sample placement process.
# After we find the best SPR move for the given node, we execute it using the cutAndPasteNode() function.
def traverseTreeForTopologyUpdate(node,mutMatrix):
	#track if the root has changed, so that the new root node can be returned.
	newRoot=None
	# has the branch length been updated so far? If we change the branch length above the current node, even if we don't perform any SPR, we still need to update the genome lists in the tree.
	bLenChanged=False
	totalImprovement=0.0
	# we avoid the root node since it cannot be re-placed with SPR moves
	if node.up!=None:
		#evaluate current placement
		parentNode=node.up
		if parentNode.children[0]==node:
			child=0
			vectUp=parentNode.probVectUpRight
		else:
			child=1
			vectUp=parentNode.probVectUpLeft
		#score of current tree
		bestCurrenBLen=node.dist
		bestCurrentLK=appendProbNode(vectUp,node.probVect,bestCurrenBLen,mutMatrix)
		if bestCurrentLK<thresholdTopologyPlacement:
			#try different branch lengths for the current node placement (just in case branch length can be improved, in which case it counts both as tree improvment and
			# better strategy to find a new placement).
			if not node.dist:
				originalLK=bestCurrentLK
				bestCurrenBLen=minBLen
				bestCurrentLK=appendProbNode(vectUp,node.probVect,bestCurrenBLen,mutMatrix)
			bestSplit=1.0
			newSplit=0.5
			while newSplit*bestCurrenBLen>minBLen:
				newLK=appendProbNode(vectUp,node.probVect,newSplit*bestCurrenBLen,mutMatrix)
				if newLK>bestCurrentLK:
					bestCurrentLK=newLK
					bestSplit=newSplit
					bLenChanged=True
				else:
					break
				newSplit=bestSplit/2
			if bestSplit>0.7:
				newSplit=2.0
				while newSplit*bestCurrenBLen<maxBLen:
					newLK=appendProbNode(vectUp,node.probVect,newSplit*bestCurrenBLen,mutMatrix)
					if newLK>bestCurrentLK:
						bestCurrentLK=newLK
						bestSplit=newSplit
						bLenChanged=True
					else:
						break
					newSplit=bestSplit*2
			bestCurrenBLen=bestCurrenBLen*bestSplit
			if not node.dist:
				if originalLK>bestCurrentLK:
					bestCurrentLK=originalLK

		if bestCurrentLK<thresholdTopologyPlacement:
			#now find the best place on the tree where to re-attach the subtree rooted ar "node"
			#but to do that we need to consider new vector probabilities after removing the node that we want to replace
			# this is done using findBestParentTopology().
			topologyUpdated=False
			bestNodeSoFar , bestLKdiff , bestIsMidNode = findBestParentTopology(parentNode,child,bestCurrentLK,bestCurrenBLen,mutMatrix)
			if bestLKdiff>thresholdProb2:
				print("Strange, LK cost is positive")
				exit()
			elif bestLKdiff<-1e50:
				print("Error: found likelihood cost is very heavy, this might mean that the reference used is not the same used to generate the input diff file")
				exit()
			if bestLKdiff+thresholdTopologyPlacement>bestCurrentLK:
				if bestNodeSoFar==parentNode:
					print("Strange, re-placement is at same node")
				elif bestNodeSoFar==parentNode.children[1-child] and bestIsMidNode:
					print("Re-placement is above sibling node")
				else:
					topNode1=bestNodeSoFar
					while (not topNode1.dist) and (topNode1.up!=None): #reach the top of a multifurcation, which is the only place in a multifurcation where placement is allowed.
						topNode1=topNode1.up
					if topNode1!=bestNodeSoFar:
						print("Strange, placement node not at top of polytomy")
					topNode2=parentNode
					while (not topNode2.dist) and (topNode2.up!=None):
						topNode2=topNode2.up
					if topNode2==topNode1 and (not bestIsMidNode):
						print("Re-placement at same polytomy, not going forward with move")
					else:
						totalImprovement=(bestLKdiff-bestCurrentLK)
						newRoot = cutAndPasteNode(node,bestNodeSoFar,bestIsMidNode,bestCurrenBLen,bestLKdiff,mutMatrix)
						topologyUpdated=True
				if (not topologyUpdated) and bLenChanged:
					node.dist=bestCurrenBLen
					updatePartialsFromTop(node,vectUp,mutMatrix)
					updatePartialsFromBottom(node.up,node.probVect,child,node,mutMatrix)
			elif bLenChanged:
				node.dist=bestCurrenBLen
				updatePartialsFromTop(node,vectUp,mutMatrix)
				updatePartialsFromBottom(node.up,node.probVect,child,node,mutMatrix)
		elif bLenChanged:
			node.dist=bestCurrenBLen
			updatePartialsFromTop(node,vectUp,mutMatrix)
			updatePartialsFromBottom(node.up,node.probVect,child,node,mutMatrix)
	return newRoot,totalImprovement


#traverse the tree (here the input "node" will usually be the root), and for each dirty node ancountered, call traverseTreeForTopologyUpdate() 
# to attempt an SPR move by cutting the subtree rooted at this dirty node and trying to re-append it elsewhere.
def startTopologyUpdates(node,mutMatrix):
	nodesToVisit=[node]
	totalImprovement=0.0
	newRoot=None
	numNodes=0
	while nodesToVisit:
		newNode=nodesToVisit.pop()
		for c in newNode.children:
			#if c.dirty:
			nodesToVisit.append(c)
		if newNode.dirty:
			newNode.dirty=False
			newRoot2,improvement=traverseTreeForTopologyUpdate(newNode,mutMatrix)
			totalImprovement+=improvement
			if newRoot2!=None:
				newRoot=newRoot2
			numNodes+=1
			if (numNodes%1000)==0:
				print("Processed topology for "+str(numNodes)+" nodes.")
	return newRoot,totalImprovement





def createNewick(node):
	nextNode=node
	stringList=[]
	direction=0
	while nextNode!=None:
		if nextNode.children:
			if direction==0:
				stringList.append("(")
				nextNode=nextNode.children[0]
			elif direction==1:
				stringList.append(",")
				nextNode=nextNode.children[1]
				direction=0
			else:
				stringList.append("):"+str(nextNode.dist))
				if nextNode.up!=None:
					if nextNode.up.children[0]==nextNode:
						direction=1
					else:
						direction=2
				nextNode=nextNode.up
		else:
			if len(nextNode.minorSequences)>0:
				stringList.append("("+nextNode.name+":0")
				for s2 in nextNode.minorSequences:
					stringList.append(","+s2+":0")
				stringList.append("):"+str(nextNode.dist))
			else:
				stringList.append(nextNode.name+":"+str(nextNode.dist))
			if nextNode.up!=None:
				if nextNode.up.children[0]==nextNode:
					direction=1
				else:
					direction=2
			nextNode=nextNode.up
	stringList.append(";")
	return "".join(stringList)






def createBinaryNewick(node):
	nextNode=node
	stringList=[]
	direction=0
	while nextNode!=None:
		if nextNode.children:
			if direction==0:
				stringList.append("(")
				nextNode=nextNode.children[0]
			elif direction==1:
				stringList.append(",")
				nextNode=nextNode.children[1]
				direction=0
			else:
				if nextNode.dist>thresholdProb:
					stringList.append("):"+str(nextNode.dist))
				else:
					stringList.append("):"+str(thresholdProb))
				if nextNode.up!=None:
					if nextNode.up.children[0]==nextNode:
						direction=1
					else:
						direction=2
				nextNode=nextNode.up
		else:
			if len(nextNode.minorSequences)>0:
				for i in nextNode.minorSequences:
					stringList.append("(")
				stringList.append(nextNode.name+":")
				for s2 in nextNode.minorSequences:
					stringList.append(str(thresholdProb)+","+s2+":"+str(thresholdProb)+"):")
				if nextNode.dist>thresholdProb:
					stringList.append(str(nextNode.dist))
				else:
					stringList.append(str(thresholdProb))
			else:
				if nextNode.dist>thresholdProb:
					stringList.append(nextNode.name+":"+str(nextNode.dist))
				else:
					stringList.append(nextNode.name+":"+str(thresholdProb))
			if nextNode.up!=None:
				if nextNode.up.children[0]==nextNode:
					direction=1
				else:
					direction=2
			nextNode=nextNode.up
	stringList.append(";")
	return "".join(stringList)



numNodes=[0,0,0,0,0]

#Given a tree, and a final substitution rate matrix, re-calculate all genome lists within the tree according to this matrix.
# this is useful once the matrix estimation has finished, to make sure all genome lists replect this matrix. 
def reCalculateAllGenomeLists(root,mutMatrix):
	#first pass to update all lower likelihoods.
	node=root
	#direction 0 means from parent, direction 1 means from a child
	lastNode=None
	direction=0
	while node!=None:
		if direction==0:
			if node.children:
				node=node.children[0]
			else:
				numNodes[0]+=1
				for entry in node.probVect:
					if entry[0]<4:
						numNodes[1]+=1
					elif entry[0]==4:
						numNodes[2]+=1
					elif entry[0]==5:
						numNodes[3]+=1
					else:
						numNodes[4]+=1

				lastNode=node
				node=node.up
				direction=1
		else :
			if lastNode==node.children[0]:
				node=node.children[1]
				direction=0
			else:
				node.probVect=mergeVectors(node.children[0].probVect,node.children[0].dist,node.children[1].probVect,node.children[1].dist,mutMatrix,returnLK=False)

				numNodes[0]+=1
				for entry in node.probVect:
					if entry[0]<4:
						numNodes[1]+=1
					elif entry[0]==4:
						numNodes[2]+=1
					elif entry[0]==5:
						numNodes[3]+=1
					else:
						numNodes[4]+=1

				lastNode=node
				node=node.up
				direction=1

				

	#now update the other genome lists for the root
	node=root
	node.probVectTot=rootVector(node.probVect,False,mutMatrix)
	if node.children:
		node.probVectUpRight=rootVector(node.children[1].probVect,node.children[1].dist,mutMatrix)
		node.probVectUpLeft=rootVector(node.children[0].probVect,node.children[0].dist,mutMatrix)
		
		#now traverse the tree downward and update the non-lower genome lists for all other nodes of the tree.
		lastNode=None
		node=node.children[0]
		direction=0
		while node!=None:
			if direction==0:
				if node==node.up.children[0]:
					vectUp=node.up.probVectUpRight
				else:
					vectUp=node.up.probVectUpLeft
				if node.dist:
					newVect=mergeVectorsUpDown(vectUp,node.dist,node.probVect,False,mutMatrix)
					#if areVectorsDifferent(newVect,node.probVectTot):
					#	updatedGenomeLists+=1
					node.probVectTot=newVect
					if node.probVectTot==None:
						if node==node.up.children[0]:
							updateBLen(node.up,0,mutMatrix)
						else:
							updateBLen(node.up,1,mutMatrix)
					else:
						node.probVectTotUp=mergeVectorsUpDown(vectUp,node.dist/2,node.probVect,node.dist/2,mutMatrix)
						if node.dist>=2*minBLenForMidNode:
							createFurtherMidNodes(node,vectUp)
				if node.children:
					node.probVectUpRight=mergeVectorsUpDown(vectUp,node.dist,node.children[1].probVect,node.children[1].dist,mutMatrix)
					node.probVectUpLeft=mergeVectorsUpDown(vectUp,node.dist,node.children[0].probVect,node.children[0].dist,mutMatrix)
					node=node.children[0]
				else:
					lastNode=node
					node=node.up
					direction=1
			else:
				if lastNode==node.children[0]:
					node=node.children[1]
					direction=0
				else:
					lastNode=node
					node=node.up
					direction=1






if runOnlyExample:
	data={}
	data["refSeq"]=[]
	data["Seq1"]=[('t',1),('-',10,1)]
	data["Seq1bis"]=[('t',1),('-',11,1)]
	data["Seq2"]=[('t',1),('a',3)]
	data["Seq3"]=[('t',1),('a',3),('t',5)]
	data["Seq4"]=[('a',2)]
	data["Seq5"]=[('a',2),('t',4)]
	data["Seq6"]=[('a',2),('t',6)]
	print(data)
	print(ref[:10])
	samples=data.keys()

#distances=distancesFromRefPunishNs(data,samples)
distances=distancesFromRefPunishNs(data)
print("Distances from the reference calculated")
#extract root genome among those closest to the reference but not empty
firstSample=distances.pop()

#initialize tree to just the initial root sample
t1=Tree(name=firstSample[1])
#t1.probVect=probVectTerminalNode(data.pop(firstSample[1]))
t1.probVect=probVectTerminalNode(data[firstSample[1]])
data[firstSample[1]]=None
t1.probVectTot=rootVector(t1.probVect,False,mutMatrix)
t1.minorSequences=[]

timeFinding=0.0
timePlacing=0.0
numSamples=0
while distances:
	d=distances.pop()
	numSamples+=1
	sample=d[1]
	newPartials=probVectTerminalNode(data[sample])
	data[sample]=None
	if (numSamples%updateSubstMatrixEveryThisSamples)==0:
		if updateSubMatrix(pseudoMutCounts,model,mutMatrix):
			for i in range(lRef):
				cumulativeRate[i+1]=cumulativeRate[i]+nonMutRates[refIndeces[i]]
			for i in range4:
				nonMutRates[i]=mutMatrix[i][i]
	if (numSamples%1000)==0:
		print("Sample num "+str(numSamples))
	start=time()
	node , bestNewLK, isMidNode, bestUpLK, bestDownLK, bestDownNode, adjustBLen=findBestParent(t1,newPartials,sample,mutMatrix)
	timeFinding+=(time()-start)
	if bestNewLK<0.5:
		start=time()
		newRoot=placeSampleOnTreeNew(node,newPartials,sample,bestNewLK,isMidNode, bestUpLK, bestDownLK, bestDownNode,mutMatrix,pseudoMutCounts, adjustBLen)
		if newRoot!=None:
			t1=newRoot
		timePlacing+=(time()-start)

if runOnlyExample:
	print("Tree after initial placement:")
	newickString=createNewick(t1)
	print(newickString)
	print("Now making change to the tree to create imperfection. Moving")
	nodeToReplace=t1.children[1].children[0].children[0].children[1]
	destination=t1.children[0].children[1].children[1]
	newickString=createNewick(nodeToReplace)
	print(newickString)
	print("to")
	newickString=createNewick(destination)
	print(newickString)
	cutAndPasteNode(nodeToReplace,destination,False,0.0001,-100.0,mutMatrix)
	newickString=createNewick(t1)
	print(newickString)

data.clear()



#recalculate all genome lists according to the final substitution model
if numTopologyImprovements>0:
	start=time()
	reCalculateAllGenomeLists(t1,mutMatrix)
	timeRecalculation=time()-start
	print("Time to recalculate all genome lists: "+str(timeRecalculation))

print("Number of nodes: "+str(numNodes[0]))
print("Os per node: "+str(float(numNodes[4])/numNodes[0]))
print("Nucs per node: "+str(float(numNodes[1])/numNodes[0]))
print("Ns per node: "+str(float(numNodes[3])/numNodes[0]))
print("R per node: "+str(float(numNodes[2])/numNodes[0]))
print("Non-O per node: "+str(float(numNodes[1]+numNodes[2]+numNodes[3])/numNodes[0]))


#now run topological improvements
timeTopology=0.0
for i in range(numTopologyImprovements):
	print("Starting topological impromevement attempt traversing number "+str(i+1))
	start=time()
	setAllDirty(t1)
	newRoot,improvement=startTopologyUpdates(t1,mutMatrix)
	if newRoot!=None:
		t1=newRoot
	timeForUpdatingTopology=(time()-start)
	print("Time for this traversal of the tree to update the topology or branch length: "+str(timeForUpdatingTopology))
	timeTopology+=timeForUpdatingTopology
	print("LK improvement apparently brought: "+str(improvement))

	if runOnlyExample:
		print("Tree after improvement "+str(i))
		newickString=createNewick(t1)
		print(newickString)

	if improvement<thresholdLogLKtopology:
		print("Small improvement, stopping topological search.")
		break
	
	#run improvements only on the nodes that have been affected by some changes in the last round, and so on
	start=time()
	subRound=0
	while True:
		print("Topological subround "+str(subRound+1))
		newRoot,improvement=startTopologyUpdates(t1,mutMatrix)
		if newRoot!=None:
			t1=newRoot
		print("LK improvement apparently brought: "+str(improvement))
		if improvement<thresholdLogLKtopology:
			break
	timeForUpdatingTopology=(time()-start)
	print("Time for the subrounds of this traversal of the tree: "+str(timeForUpdatingTopology))
	timeTopology+=timeForUpdatingTopology


#Free space by deleting the genome lists in the tree.
nextLeaves=[t1]
while nextLeaves:
	node=nextLeaves.pop()
	node.probVect=None
	if node.dist:
		node.probVectTot=None
		if node.up!=None:
			node.probVectTotUp=None
			node.furtherMidNodes=None
	if node.children:
		node.probVectUpRight=None
		node.probVectUpLeft=None
		for c in node.children:
			nextLeaves.append(c)
print("Deleted genome lists.")

#Read input file to collect all sample names
fileI=open(inputFile)
line=fileI.readline()
sampleNames=[]
while line!="" and line!="\n":
	name=line.replace(">","").replace("\n","")
	line=fileI.readline()
	while line!="" and line!="\n" and line[0]!=">":
		line=fileI.readline()
	sampleNames.append(name)
fileI.close()
print("Sample names read.")

#put sample names in the tree
nextLeaves=[t1]
name=""
while nextLeaves:
	node=nextLeaves.pop()
	if not node.children:
		name=sampleNames[node.name]
		sampleNames[node.name]=None
		node.name=name
		for m in range(len(node.minorSequences)):
			name=sampleNames[node.minorSequences[m]]
			sampleNames[node.minorSequences[m]]=None
			node.minorSequences[m]=name
	else:
		for c in node.children:
			nextLeaves.append(c)
print("Sample names assigned.")


if binaryTree:
	newickString=createBinaryNewick(t1)
else:
	newickString=createNewick(t1)

if runOnlyExample:
	print("Final tree:")
	print(newickString)
	exit()

file=open(outputFile+"_tree.tree","w")
file.write(newickString)
file.close()
file=open(outputFile+"_subs.txt","w")
for i in range4:
	for j in range4:
		file.write(str(mutMatrix[i][j])+"\t")
	file.write("\n")
file.close()
print("Missed minor samples: "+str(totalMissedMinors[0]))
print("Final Substitution matrix:")
print(mutMatrix)
print("Time spent finding placement nodes: "+str(timeFinding))
print("Time spent placing samples on the tree: "+str(timePlacing))
print("Time spent in total updating the topology and branch lengths: "+str(timeTopology))
exit()





