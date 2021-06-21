import sys
import math
import numpy as np
import argparse
from ete3 import Tree

#©EMBL-European Bioinformatics Institues, 2021

#Estimate a tree using fastLK from a diff format and using iterative sample placement.

parser = argparse.ArgumentParser(description='Run estimation of mutation rates and synonymous site selection from SARS-CoV-2 data, and generate plots.')
parser.add_argument('--path',default="", help='path where to find files and plot results.')
parser.add_argument('--input',default="", help='input diff file name found in folder specified by --path.')
parser.add_argument('--reference',default="EPI_ISL_402124_lowercase.fasta", help='input reference file name found in folder specified by --path.')
parser.add_argument('--output',default="iterativeFastLK.tree", help='output newick file name to be written in the folder --path.')
parser.add_argument("--onlyNambiguities", help="Treat all ambiguities as N (total missing information).", action="store_true")
parser.add_argument("--useLogs", help="Calculate logarithms of non-mutation probabilities, otherwise approximate them.", action="store_true")
parser.add_argument("--thresholdProb",help="relative probability threshold used to ignore possible states with very low probabilities.",  type=float, default=0.0000001)
parser.add_argument("--allowedFails",help="Number of times one can go down the tree without inclreasing placement likelihood before the tree traversal is stopped (only applies to non-0 branch lengths).",  type=int, default=2)
parser.add_argument("--verbose", help="Print to screen a lot of stuff.", action="store_true")
parser.add_argument("--useJC", help="Use JC model.", action="store_true")
args = parser.parse_args()

onlyNambiguities=args.onlyNambiguities
useLogs=args.useLogs
thresholdProb=args.thresholdProb
verbose=args.verbose
pathSimu=args.path
inputFile=pathSimu+args.input
outputFile=pathSimu+args.output
refFile=pathSimu+args.reference
allowedFails=args.allowedFails
useJC=args.useJC

alleles={"A":0,"C":1,"G":2,"T":3}
allelesList=["A","C","G","T"]
allelesLow={"a":0,"c":1,"g":2,"t":3}
allelesListLow=["a","c","g","t"]
ambiguities={"y":[0.0,1.0,0.0,1.0],"r":[1.0,0.0,1.0,0.0],"w":[1.0,0.0,0.0,1.0],"s":[0.0,1.0,1.0,0.0],"k":[0.0,0.0,1.0,1.0],"m":[1.0,1.0,0.0,0.0],"d":[1.0,0.0,1.0,1.0],"v":[1.0,1.0,1.0,0.0],"h":[1.0,1.0,0.0,1.0],"b":[0.0,1.0,1.0,1.0]}

#collect reference
def collectReference(fileName):
	file=open(fileName)
	line=file.readline()
	ref=""
	while line!="":
		line=file.readline()
		ref+=line.replace("\n","")
	lRef=len(ref)
	print("Ref genome length: "+str(lRef))
	file.close()
	#vector to count how many bases of each type are cumulatively in the reference genome up to a certain position
	cumulativeBases=np.zeros((lRef+1,4))
	for i in range(lRef):
		for k in range(4):
			cumulativeBases[i+1][k]=cumulativeBases[i][k]
		cumulativeBases[i+1][allelesLow[ref[i]]]+=1
	print("cumulativeBases")
	print(cumulativeBases[-1])
	rootFreqs=np.zeros(4)
	rootFreqsLog=np.zeros(4)
	for i in range(4):
		rootFreqs[i]=cumulativeBases[-1][i]/float(lRef)
		rootFreqsLog[i]=math.log(rootFreqs[i])
	print("ref base frequencies and log")
	print(rootFreqs)
	print(rootFreqsLog)
	return ref, cumulativeBases, rootFreqs, rootFreqsLog

ref, cumulativeBases, rootFreqs, rootFreqsLog = collectReference(refFile)
lRef=len(ref)
if useJC:
	rootFreqs=[0.25,0.25,0.25,0.25]
	rootFreqsLog=[math.log(0.25),math.log(0.25),math.log(0.25),math.log(0.25)]




def readConciseAlignment(fileName):
	#start = time.time()
	fileI=open(fileName)
	line=fileI.readline()
	nSeqs=0
	data={}
	while line!="" and line!="\n":
		nSeqs+=1
		seqList=[]
		name=line.replace(">","").replace("\n","")
		line=fileI.readline()
		while line!="" and line!="\n" and line[0]!=">":
			linelist=line.split()
			if len(linelist)>2:
				entry=(linelist[0],int(linelist[1]),int(linelist[2]))
			else:
				entry=(linelist[0],int(linelist[1]))
			seqList.append(entry)
			line=fileI.readline()
		data[name]=seqList
	fileI.close()
	#time2 = time.time() - start
	#print("Time to read DNA reduced data file: "+str(time2))
	print(str(nSeqs)+" sequences in diff file.")
	return data

#read sequence data from file
data=readConciseAlignment(inputFile)


samples=data.keys()
print(str(len(samples))+" sequences in the concise DNA data file")

range4=range(4)
thresholdProb2=thresholdProb*thresholdProb
thresholdProb4=thresholdProb2*thresholdProb2

#if probability mass is concentrated in one nucleotide, simply the entry from O to another type.
def simplfy(vec,refA):
	maxP=0.0
	maxI=0
	for i in range4:
		if vec[i]>maxP:
			maxP=vec[i]
			maxI=i
	if maxP<thresholdProb4:
		print("Inside simplify(), all values in vector are too small - something wrong?")
		print(vec)
		exit()
	newVec=vec/maxP
	numA=0
	for i in range4:
		if newVec[i]>thresholdProb:
			numA+=1
	if numA==1:
		if allelesListLow[maxI]==refA:
			return "R", newVec, maxP
		else:
			return allelesList[maxI], newVec, maxP
	else:
		return "O", newVec, maxP

#Shorten genome list by mergeing together R entries that are mergeable
def shorten(vec):
	newVec=[]
	entryOld=vec[0]
	for i in range(len(vec)-1):
		if entryOld[0]!="R":
			newVec.append(entryOld)
			entryOld=vec[i+1]
		else:
			entryNew=vec[i+1]
			if entryNew[0]!="R":
				newVec.append(entryOld)
				entryOld=vec[i+1]
			else:
				if abs(entryOld[3]-entryNew[3])<thresholdProb:
					entryOld[2]+=entryNew[2]
				else:
					newVec.append(entryOld)
					entryOld=vec[i+1]
	newVec.append(entryOld)
	return newVec

#define partial likelihood vector for a sample
def probVectTerminalNode(diffs,bLen):
	if diffs is None:
		probVect=[["N",1,lRef,bLen]]
		return probVect
	pos=1
	probVect=[]
	for m in diffs:
			currPos=m[1]
			if currPos>pos:
				#region where the node with branch length bLen is identical to the ref.
				probVect.append(["R",pos,currPos-pos,bLen,False])
				pos=currPos
			if m[0]=="n" or m[0]=="-":
				length=m[2]
				#region with no info, store first position and length.
				probVect.append(["N",currPos,length,bLen,False])
				pos=currPos+length
			elif m[0] in allelesLow:
				#position at which node allele is sure but is different from the reference.
				probVect.append([m[0].upper(),currPos,1,bLen,False])
				pos=currPos+1
			else:
				# non-"n" ambiguity character; for now interpret this as ambiguity instead of as a polymorphism.
				if onlyNambiguities:
					# if user asks to, to make things easier, interpret any ambiguity as an "n".
					probVect.append(["N",currPos,1,bLen,False])
				else:
					#otherwise, store as an "other" scenario, where each nucleotide has its own partial likelihood.
					probVect.append(["O",currPos,1,bLen,ambiguities[m[0]]])
				pos=currPos+1
	if pos<=lRef:
		probVect.append(["R",pos,lRef+1-pos,bLen,False])
	return probVect

#preliminary nuc mutation rate matrix, from De Maio et al 2021
mutMatrix=[[0.0,0.039,0.310,0.123],[0.140,0.0,0.022,3.028],[0.747,0.113,0.0,2.953],[0.056,0.261,0.036,0.0]]
if useJC:
	mutMatrix=[[0.0,0.333,0.333,0.333],[0.333,0.0,0.333,0.333],[0.333,0.333,0.0,0.333],[0.333,0.333,0.333,0.0]]
mutMatrix[0][0]=-(mutMatrix[0][1]+mutMatrix[0][2]+mutMatrix[0][3])
mutMatrix[1][1]=-(mutMatrix[1][0]+mutMatrix[1][2]+mutMatrix[1][3])
mutMatrix[2][2]=-(mutMatrix[2][0]+mutMatrix[2][1]+mutMatrix[2][3])
mutMatrix[3][3]=-(mutMatrix[3][0]+mutMatrix[3][1]+mutMatrix[3][2])


#Sort samples based on distance from reference, but punishing more isolated N's and ambiguity codes.
def distancesFromRefPunishNs(data):
	sampleDistances=[]
	for sample in samples:
		diffs=data[sample]
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
				pos=currPos+m[2]
				#if m[2]==1:
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
		sampleDistances.append((diffNum*1000+lRef-comparisons,sample))

	from operator import itemgetter
	print("Now doing sorting")
	sampleDistances.sort(key=itemgetter(0))
	return sampleDistances



#function to check is one sequence is inferior to another
#returns 0 when the 2 sequences are not comparable, otherwise returns 1 if the first is more informative or if they are identical, and 2 otherwise.
def isMinorSequence(probVect1,probVect2):
			indexEntry1=0
			indexEntry2=0
			pos=1
			entry1=probVect1[indexEntry1]
			pos1=entry1[1]
			if entry1[0]!="N" and entry1[0]!="R":
				end1=pos1
			else:
				end1=pos1+entry1[2]-1
			entry2=probVect2[indexEntry2]
			pos2=entry2[1]
			if entry2[0]!="N" and entry2[0]!="R":
				end2=pos2
			else:
				end2=pos2+entry2[2]-1
			end=end1
			if end2<end:
				end=end2
			length=end+1-pos
			found1bigger=False
			found2bigger=False
			while True:
				if entry1[0]!=entry2[0]:
					if entry1[0]=="N":
						found2bigger=True
					elif entry2[0]=="N":
						found1bigger=True
					elif entry1[0]=="O":
						found2bigger=True
					elif entry2[0]=="O":
						found1bigger=True
					else:
						found2bigger=True
						found1bigger=True
				elif entry1[0]=="O":
					for j in range4:
						if entry2[4][j]>entry1[4][j]+0.01:
							found1bigger=True
						elif entry1[4][j]>entry2[4][j]+0.01:
							found2bigger=True

				pos+=length
				if pos>lRef:
					break
				if found1bigger and found2bigger:
					break
				if pos>end1:
					indexEntry1+=1
					entry1=probVect1[indexEntry1]
					pos1=entry1[1]
					if pos1<pos:
						print("Error with posititons 1")
						print(entry1)
						exit()
					if entry1[0]!="N" and entry1[0]!="R":
						end1=pos1
					else:
						end1=pos1+entry1[2]-1
				if pos>end2:
					indexEntry2+=1
					entry2=probVect2[indexEntry2]
					pos2=entry2[1]
					if pos2<pos:
						print("Error with posititons 2")
						print(entry2)
						exit()
					if entry2[0]!="N" and entry2[0]!="R":
						end2=pos2
					else:
						end2=pos2+entry2[2]-1
				end=end1
				if end2<end:
					end=end2
				length=end+1-pos

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





#function to calculate likelihood cost of appending child to parent
def appendProb(probVectP,probVectC,bLen):
	Lkcost=0.0
	LkcostOriginal=0.0
	indexEntry1=0
	indexEntry2=0
	pos=1
	entry1=probVectP[indexEntry1]
	pos1=entry1[1]
	if entry1[0]!="N" and entry1[0]!="R":
		end1=pos1
	else:
		end1=pos1+entry1[2]-1
	entry2=probVectC[indexEntry2]
	pos2=entry2[1]
	if entry2[0]!="N" and entry2[0]!="R":
		end2=pos2
	else:
		end2=pos2+entry2[2]-1
	end=end1
	if end2<end:
		end=end2
	length=end+1-pos
	while True:
		if entry1[0]=="N":
			if not entry1[4]:
				print("Warning, should have been True")
			if entry2[0]=="N":
				pass
			elif entry2[0]=="R":
				for i in range4:
					Lkcost+=rootFreqsLog[i]*(cumulativeBases[end][i]-cumulativeBases[pos-1][i])
			elif entry2[0]=="O":
				tot=0.0
				for i in range4:
					tot+=rootFreqs[i]*entry2[4][i]
				Lkcost+=math.log(tot)
			else:
				i2=alleles[entry2[0]]
				Lkcost+=rootFreqsLog[i2]
		elif entry2[0]=="N":
			pass
		elif entry1[0]=="R":
			if entry2[0]=="R":
					for i in range4:
						Lkcost+=mutMatrix[i][i]*bLen*(cumulativeBases[end][i]-cumulativeBases[pos-1][i])
			elif entry2[0]=="O":
				i1=allelesLow[ref[pos-1]]
				tot=0.0
				for j in range4:
					if i1!=j:
						tot+=mutMatrix[i1][j]*bLen*entry2[4][j]
					else:
						tot+=(1.0+mutMatrix[i1][j]*bLen)*entry2[4][j]
				if tot>sys.float_info.min:
					Lkcost+=math.log(tot)
				else:
					return float("-inf")
			else:
				if bLen<thresholdProb2:
					return float("-inf")
				i2=alleles[entry2[0]]
				i1=allelesLow[ref[pos-1]]
				Lkcost+=math.log(mutMatrix[i1][i2]*bLen)
		elif entry1[0]=="O":
			tot1=0.0
			for i in range4:
				tot1+=entry1[4][i]
			LkcostOriginal+=math.log(tot1)
			if entry2[0]=="O":
				tot=0.0
				for i in range4:
					for j in range4:
						if i==j:
							tot+=(1.0+mutMatrix[i][i]*bLen)*entry1[4][i]*entry2[4][j]
						else:
							tot+=mutMatrix[i][j]*bLen*entry1[4][i]*entry2[4][j]
				if tot>sys.float_info.min:
					Lkcost+=math.log(tot)
				else:
					return float("-inf")
			else:
				if entry2[0]=="R":
					i2=allelesLow[ref[pos-1]]
				else:
					i2=alleles[entry2[0]]
				tot=0.0
				for i in range4:
					if i==i2:
						tot+=entry1[4][i]*(1.0+mutMatrix[i][i]*bLen)
					else:
						tot+=entry1[4][i]*mutMatrix[i][i2]*bLen
				if tot>sys.float_info.min:
					Lkcost+=math.log(tot)
				else:
					return float("-inf")
		#entry1 is a non-ref nuc
		else:
			i1=alleles[entry1[0]]
			if entry2[0]==entry1[0]:
				Lkcost+=mutMatrix[i1][i1]*bLen
			else:
				if entry2[0]=="O":
						tot=0.0
						for j in range4:
							if i1==j:
								tot+=(1.0+mutMatrix[i1][i1]*bLen)*entry2[4][j]
							else:
								tot+=mutMatrix[i1][j]*bLen*entry2[4][j]
						if tot>sys.float_info.min:
							Lkcost+=math.log(tot)
						else:
							return float("-inf")

				else:
					if bLen<thresholdProb2:
						return float("-inf")
					if entry2[0]=="R":
						i2=allelesLow[ref[pos-1]]
					else:
						i2=alleles[entry2[0]]
					Lkcost+=math.log(mutMatrix[i1][i2]*bLen)

		pos+=length
		if pos>lRef:
			break
		if pos>end1:
			indexEntry1+=1
			entry1=probVectP[indexEntry1]
			pos1=entry1[1]
			if entry1[0]!="N" and entry1[0]!="R":
				end1=pos1
			else:
				end1=pos1+entry1[2]-1
		if pos>end2:
			indexEntry2+=1
			entry2=probVectC[indexEntry2]
			pos2=entry2[1]
			if entry2[0]!="N" and entry2[0]!="R":
				end2=pos2
			else:
				end2=pos2+entry2[2]-1
		end=end1
		if end2<end:
			end=end2
		length=end+1-pos
	return Lkcost-LkcostOriginal

#number of samples that could have been placed as major of another sample, but weren't due to sample placement order
totalMissedMinors=[0]
#function to find the best node in the tree where to append the new sample, based on the node's tot posterior state probabilities.
def findBestParent(t1,diffs,sample,bLen,bestLKdiff,bestNodeSoFar,failedPasses,bestIsMidNode):
	if t1.is_leaf():
		probVect=t1.probVect
		comparison=isMinorSequence(probVect,diffs)
		if comparison==1:
			t1.minorSequences.append(sample)
			return t1, 1.0, False
		elif comparison==2:
			if verbose:
				print("Second sequence is more informative than first, this could be done better, but for now this fact is ignored and the two sequences are not merged.")
				print(t1.name)
				print(sample)
				print(probVect)
				print(diffs)
			totalMissedMinors[0]+=1
	
	if t1.dist>thresholdProb2 and t1.up!=None:
		probVect=t1.probVectTot
		LKdiff=appendProb(probVect,diffs,bLen)
		isMidNode=False
		if verbose:
			print("Trying a new parent")
			if t1.is_leaf():
				print(t1.name)
			print(probVect)
			#print(diffs)
			print("Failed passes before this step: "+str(failedPasses)+" logLK score: "+str(LKdiff))
		#print(t1.write())
		#print(t1.dist)
		probVect2=t1.probVectTotUp
		LKdiff2=appendProb(probVect2,diffs,bLen)
		if LKdiff2>LKdiff:
			LKdiff=LKdiff2
			probVect=probVect2
			isMidNode=True
		if LKdiff>bestLKdiff:
			if verbose:
				print("New best LK")
			bestLKdiff=LKdiff
			bestNodeSoFar=t1
			failedPasses=0
			bestIsMidNode=isMidNode
		elif LKdiff<bestLKdiff-0.1:
			failedPasses+=1
			if failedPasses>allowedFails:
				return bestNodeSoFar , bestLKdiff , bestIsMidNode

	for c in t1.children:
		node, score , childIsMidNode = findBestParent(c,diffs,sample,bLen,bestLKdiff,bestNodeSoFar,failedPasses,bestIsMidNode)
		if score>0.5:
			return node, score , False
		if score>bestLKdiff:
			bestLKdiff=score
			bestNodeSoFar=node
			bestIsMidNode=childIsMidNode
	
	return bestNodeSoFar , bestLKdiff , bestIsMidNode



#for the root, return a vector taking into account the root frequency of each state
def rootVector(probVect,rootFreqs,bLen):
	newProbVect=[]
	for entry in probVect:
		newEntry=list(entry)
		if entry[0]=="N":
			newEntry[4]=True
		elif entry[0]=="O":
			newEntry[4]=[]
			for i in range4:
					tot=0.0
					for j in range4:
						if i==j:
							tot+=rootFreqs[i]*(1.0+mutMatrix[i][i]*(entry[3]+bLen))*entry[4][j]
						else:
							tot+=rootFreqs[i]*mutMatrix[i][j]*(entry[3]+bLen)*entry[4][j]
					newEntry[4].append(tot)
			newEntry[3]=0.0
		else:
			newEntry[4]=True
			#distance from the root of the upward branch
			newEntry[3]+=bLen
			#distance from the root on the downward branch
			newEntry.append(0.0)
		newProbVect.append(newEntry)
	return newProbVect


#for the root, return a vector taking into account the root frequency of each state (same as above, but assuming the branch length separating from the root is 0)
def rootVector0(probVect,rootFreqs):
	newProbVect=[]
	for entry in probVect:
		newEntry=list(entry)
		if entry[0]=="N":
			newEntry[4]=True
		elif entry[0]=="O":
			newEntry[4]=[]
			for i in range4:
				newEntry[4].append(entry[4][i]*rootFreqs[i])
			newEntry[3]=0.0
		else:
			newEntry[4]=True
			#distance from the root on the downward branch
			newEntry.append(0.0)
		newProbVect.append(newEntry)
	return newProbVect





#merge two partial likelihood vectors, one from above and one from below
def mergeVectorsUpDown(probVectUp,bLenUp,probVectDown,bLenDown,mutMatrix):
			probVect1= probVectUp
			probVect2=probVectDown
			probVect=[]
			indexEntry1=0
			indexEntry2=0
			pos=1
			entry1=probVect1[indexEntry1]
			pos1=entry1[1]
			if entry1[0]!="N" and entry1[0]!="R":
				end1=pos1
			else:
				end1=pos1+entry1[2]-1
			entry2=probVect2[indexEntry2]
			pos2=entry2[1]
			if entry2[0]!="N" and entry2[0]!="R":
				end2=pos2
			else:
				end2=pos2+entry2[2]-1
			end=end1
			if end2<end:
				end=end2
			length=end+1-pos
			while True:
				if entry1[0]=="N":
					probVect.append(list(entry2))
					if entry2[0]=="N" or entry2[0]=="R": 
						probVect[-1][2]=length
						probVect[-1][1]=pos
					if entry2[0]!="O":
						probVect[-1][4]=True
						probVect[-1].append(0.0)
						probVect[-1][3]+=bLenDown
					else:
						probVect[-1][4]=[]
						for i in range4:
							tot=0.0
							for j in range4:
								if i==j:
									tot+=rootFreqs[i]*(1.0+mutMatrix[i][i]*(entry2[3]+bLenDown))*entry2[4][j]
								else:
									tot+=rootFreqs[i]*mutMatrix[i][j]*(entry2[3]+bLenDown)*entry2[4][j]
							probVect[-1][4].append(tot)
						probVect[-1][3]=0
				elif entry2[0]=="N":
					probVect.append(list(entry1))
					if entry1[0]=="N" or entry1[0]=="R": 
						probVect[-1][2]=length
						probVect[-1][1]=pos
					if entry1[0]=="O":
						probVect[-1][4]=list(entry1[4])
						probVect[-1][3]+=bLenUp
					elif probVect[-1][4]:
						probVect[-1][5]+=bLenUp
					else:
						probVect[-1][3]+=bLenUp
				else:
					totLen1=bLenUp+entry1[3]
					totLen2=bLenDown+entry2[3]
					if entry1[0]!="O" and entry1[4]:
						totLen1+=entry1[5]
					if entry2[0]!="O" and entry2[4]:
						totLen2+=entry2[5]
					if totLen1+totLen2<thresholdProb2 and entry1[0]!=entry2[0] and entry1[0]!="O" and entry2[0]!="O":
						print("Warning, inside mergeVectorsUpDown, branch lengths are 0, but the two vectors are different")
						print(probVectUp)
						print(bLenUp)
						print(probVectDown)
						print(bLenDown)
						print(entry1)
						print(entry2)
						return None
					elif entry2[0]=="R" and (totLen2<thresholdProb2):
						probVect.append(["R",pos,length,0.0,False])
					elif (totLen2<thresholdProb2) and entry2[0]!="O":
						probVect.append([entry2[0],pos,1,0.0,False])
					elif entry1[0]=="R":
						if entry2[0]=="R" or (totLen1<thresholdProb2):
							probVect.append(["R",pos,length,0.0,False])
						else:
							i1=allelesLow[ref[pos-1]]
							newVec=np.ones(4)
							if entry1[4]:
								rootVec=np.ones(4)
								for i in range4:
									rootVec[i]=rootFreqs[i]
									if i==i1:
										rootVec[i]*=(1.0+mutMatrix[i][i1]*(entry1[3]))
									else:
										rootVec[i]*=mutMatrix[i][i1]*(entry1[3])
								for j in range4:
									tot=0.0
									for i in range4:
										if j==i:
											tot+=(1.0+mutMatrix[i][j]*(entry1[5]+bLenUp))*rootVec[i]
										else:
											tot+=mutMatrix[i][j]*(entry1[5]+bLenUp)*rootVec[i]
									newVec[j]=tot
							else:
								for i in range4:
									if i==i1:
										newVec[i]=1.0+mutMatrix[i][i]*(entry1[3]+bLenUp)
									else:
										newVec[i]=mutMatrix[i1][i]*(entry1[3]+bLenUp)
							if entry2[0]=="O":
								for j in range4:
									tot=0.0
									for i in range4:
										if j==i:
											tot+=(1.0+mutMatrix[j][i]*(entry2[3]+bLenDown))*entry2[4][i]
										else:
											tot+=mutMatrix[j][i]*(entry2[3]+bLenDown)*entry2[4][i]
									newVec[j]*=tot
								state, newVec, maxP =simplfy(newVec,ref[pos-1])
								if state=="O":
									probVect.append([state,pos,1,0.0,newVec])
								else:
									probVect.append([state,pos,1,0.0,False])
							else:
								i2=alleles[entry2[0]]
								for i in range4:
									if i==i2:
										newVec[i]*=1.0+mutMatrix[i][i]*(entry2[3]+bLenDown)
									else:
										newVec[i]*=mutMatrix[i][i2]*(entry2[3]+bLenDown)
								probVect.append(["O",pos,1,0.0,newVec])
					elif entry1[0]=="O":
						newVec=np.ones(4)
						if entry2[0]=="O":
							for i in range4:
								tot1=0.0
								tot2=0.0
								for j in range4:
									if i==j:
										tot1+=(1.0+mutMatrix[i][i]*(entry1[3]+bLenUp))*entry1[4][j]
										tot2+=(1.0+mutMatrix[i][i]*(entry2[3]+bLenDown))*entry2[4][j]
									else:
										tot1+=mutMatrix[j][i]*(entry1[3]+bLenUp)*entry1[4][j]
										tot2+=mutMatrix[i][j]*(entry2[3]+bLenDown)*entry2[4][j]
								newVec[i]*=tot1*tot2
						else:
							for i in range4:
								tot1=0.0
								for j in range4:
									if i==j:
										tot1+=(1.0+mutMatrix[i][i]*(entry1[3]+bLenUp))*entry1[4][j]
									else:
										tot1+=mutMatrix[j][i]*(entry1[3]+bLenUp)*entry1[4][j]
								newVec[i]=tot1
							if entry2[0]=="R":
								i2=allelesLow[ref[pos-1]]
								for i in range4:
									if i==i2:
										newVec[i]*=1.0+mutMatrix[i][i]*(entry2[3]+bLenDown)
									else:
										newVec[i]*=mutMatrix[i][i2]*(entry2[3]+bLenDown)
							else:
								i2=alleles[entry2[0]]
								for i in range4:
									if i==i2:
										newVec[i]*=1.0+mutMatrix[i][i]*(entry2[3]+bLenDown)
									else:
										newVec[i]*=mutMatrix[i][i2]*(entry2[3]+bLenDown)
						state, newVec, maxP =simplfy(newVec,ref[pos-1])
						if state=="O":
							probVect.append([state,pos,1,0.0,newVec])
						else:
							probVect.append([state,pos,1,0.0,False])
					#entry1 is a non-ref nuc
					else:
						if entry2[0]==entry1[0] or (totLen1<thresholdProb2):
							probVect.append([entry1[0],pos,1,0.0,False])
						else:
							i1=alleles[entry1[0]]
							newVec=np.ones(4)
							if entry1[4]:
								rootVec=np.ones(4)
								for i in range4:
									rootVec[i]=rootFreqs[i]
									if i==i1:
										rootVec[i]*=(1.0+mutMatrix[i][i1]*(entry1[3]))
									else:
										rootVec[i]*=mutMatrix[i][i1]*(entry1[3])
								for j in range4:
									tot=0.0
									for i in range4:
										if j==i:
											tot+=(1.0+mutMatrix[i][j]*(entry1[5]+bLenUp))*rootVec[i]
										else:
											tot+=mutMatrix[i][j]*(entry1[5]+bLenUp)*rootVec[i]
									newVec[j]=tot
							else:
								for i in range4:
									if i==i1:
										newVec[i]=1.0+mutMatrix[i][i]*(entry1[3]+bLenUp)
									else:
										newVec[i]=mutMatrix[i1][i]*(entry1[3]+bLenUp)
							if entry2[0]=="O":
								for i in range4:
									tot=0.0
									for j in range4:
										if i==j:
											tot+=(1.0+mutMatrix[i][i]*(entry2[3]+bLenDown))*entry2[4][j]
										else:
											tot+=mutMatrix[i][j]*(entry2[3]+bLenDown)*entry2[4][j]
									newVec[i]*=tot
								state, newVec, maxP =simplfy(newVec,ref[pos-1])
								if state=="O":
									probVect.append([state,pos,1,0.0,newVec])
								else:
									probVect.append([state,pos,1,0.0,False])
							else:
								if entry2[0]=="R":
									i2=allelesLow[ref[pos-1]]
								else:
									i2=alleles[entry2[0]]
								for i in range4:
									if i==i2:
										newVec[i]*=1.0+mutMatrix[i][i]*(entry2[3]+bLenDown)
									else:
										newVec[i]*=mutMatrix[i][i2]*(entry2[3]+bLenDown)
								probVect.append(["O",pos,1,0.0,newVec])

				#update pos, end, etc
				pos+=length
				if pos>lRef:
					break
				if pos>end1:
					indexEntry1+=1
					entry1=probVect1[indexEntry1]
					pos1=entry1[1]
					if entry1[0]!="N" and entry1[0]!="R":
						end1=pos1
					else:
						end1=pos1+entry1[2]-1
				if pos>end2:
					indexEntry2+=1
					entry2=probVect2[indexEntry2]
					pos2=entry2[1]
					if entry2[0]!="N" and entry2[0]!="R":
						end2=pos2
					else:
						end2=pos2+entry2[2]-1
				end=end1
				if end2<end:
					end=end2
				length=end+1-pos

			if verbose:
				print("Merged up-down ")
				print(probVect)
			probVect1=probVect

			#check if the final  probVect can be simplified by merging consecutive entries
			probVect1 =shorten(probVect1)
			if verbose:
				print("Shortened up-down merging ")
				print(probVect1)
			return probVect1




#merge two childn partial likelihood vectors at the root and calculate the logLk, and return vector and logLK.
def mergeVectorsRoot(probVect1,bLen1,probVect2,bLen2,mutMatrix):
			cumulPartLk=0.0
			probVect=[]
			indexEntry1=0
			indexEntry2=0
			pos=1
			entry1=probVect1[indexEntry1]
			pos1=entry1[1]
			if entry1[0]!="N" and entry1[0]!="R":
				end1=pos1
			else:
				end1=pos1+entry1[2]-1
			entry2=probVect2[indexEntry2]
			pos2=entry2[1]
			if entry2[0]!="N" and entry2[0]!="R":
				end2=pos2
			else:
				end2=pos2+entry2[2]-1
			end=end1
			if end2<end:
				end=end2
			length=end+1-pos
			while True:
				totLen=entry1[3]+entry2[3]+bLen1+bLen2
				if entry1[0]=="N":
					probVect.append(list(entry2))
					if entry2[0]=="N" or entry2[0]=="R": 
						probVect[-1][2]=length
						probVect[-1][1]=pos
					probVect[-1][3]+=bLen2
					if entry2[0]=="O":
						probVect[-1][4]=list(entry2[4])
				elif entry2[0]=="N":
					probVect.append(list(entry1))
					if entry1[0]=="N" or entry1[0]=="R": 
						probVect[-1][2]=length
						probVect[-1][1]=pos
					probVect[-1][3]+=bLen1
					if entry1[0]=="O":
						probVect[-1][4]=list(entry1[4])
				elif entry1[0]=="R":
					if entry2[0]=="R":
						probVect.append(["R",pos,length,0.0,False])
						#Now update partial likelihood
						if useLogs:
							for i in range4:
								cumulPartLk+=math.log(1.0+mutMatrix[i][i]*totLen)*(cumulativeBases[end][i]-cumulativeBases[pos-1][i])
						else:
							for i in range4:
								cumulPartLk+=mutMatrix[i][i]*totLen*(cumulativeBases[end][i]-cumulativeBases[pos-1][i])
					elif entry2[0]=="O":
						i1=allelesLow[ref[pos-1]]
						newVec=np.ones(4)
						for i in range4:
							if i==i1:
								newVec[i]=1.0+mutMatrix[i][i]*(entry1[3]+bLen1)
							else:
								newVec[i]=mutMatrix[i][i1]*(entry1[3]+bLen1)
							tot=0.0
							for j in range4:
								if i==j:
									tot+=(1.0+mutMatrix[i][i]*(entry2[3]+bLen2))*entry2[4][j]
								else:
									tot+=mutMatrix[i][j]*(entry2[3]+bLen2)*entry2[4][j]
							newVec[i]*=tot
						state, newVec, maxP =simplfy(newVec,ref[pos-1])
						probVect.append([state,pos,1,0.0,newVec])
						cumulPartLk+=math.log(maxP)
					else:
						i2=alleles[entry2[0]]
						i1=allelesLow[ref[pos-1]]
						newVec=np.ones(4)
						for i in range4:
							if i==i2:
								newVec[i]=1.0+mutMatrix[i][i]*(entry2[3]+bLen2)
							else:
								newVec[i]=mutMatrix[i][i2]*(entry2[3]+bLen2)
							if i==i1:
								newVec[i]*=1.0+mutMatrix[i][i]*(entry1[3]+bLen1)
							else:
								newVec[i]*=mutMatrix[i][i1]*(entry1[3]+bLen1)
						probVect.append(["O",pos,1,0.0,newVec])
				elif entry1[0]=="O":
					newVec=np.ones(4)
					if entry2[3]+entry1[3]<thresholdProb2:
							entry2[3]=thresholdProb
							entry1[3]=thresholdProb
					if entry2[0]=="O":
						for i in range4:
							tot1=0.0
							tot2=0.0
							for j in range4:
								if i==j:
									tot1+=(1.0+mutMatrix[i][i]*(entry1[3]+bLen1))*entry1[4][j]
									tot2+=(1.0+mutMatrix[i][i]*(entry2[3]+bLen2))*entry2[4][j]
								else:
									tot1+=mutMatrix[i][j]*(entry1[3]+bLen1)*entry1[4][j]
									tot2+=mutMatrix[i][j]*(entry2[3]+bLen2)*entry2[4][j]
							newVec[i]*=tot1*tot2
					else:
						for i in range4:
							tot1=0.0
							for j in range4:
								if i==j:
									tot1+=(1.0+mutMatrix[i][i]*(entry1[3]+bLen1))*entry1[4][j]
								else:
									tot1+=mutMatrix[i][j]*(entry1[3]+bLen1)*entry1[4][j]
							newVec[i]=tot1
						if entry2[0]=="R":
							i2=allelesLow[ref[pos-1]]
							for i in range4:
								if i==i2:
									newVec[i]*=1.0+mutMatrix[i][i]*(entry2[3]+bLen2)
								else:
									newVec[i]*=mutMatrix[i][i2]*(entry2[3]+bLen2)
						else:
							i2=alleles[entry2[0]]
							for i in range4:
								if i==i2:
									newVec[i]*=1.0+mutMatrix[i][i]*(entry2[3]+bLen2)
								else:
									newVec[i]*=mutMatrix[i][i2]*(entry2[3]+bLen2)
					state, newVec, maxP =simplfy(newVec,ref[pos-1])
					probVect.append([state,pos,1,0.0,newVec])
					cumulPartLk+=math.log(maxP)
				#entry1 is a non-ref nuc
				else:
					if entry2[0]==entry1[0]:
						i2=alleles[entry2[0]]
						probVect.append([entry1[0],pos,1,0.0,False])
						#Now update partial likelihood
						if useLogs:
							cumulPartLk+=math.log(1.0+mutMatrix[i2][i2]*totLen)
						else:
							cumulPartLk+=mutMatrix[i2][i2]*totLen
					else:
						i1=alleles[entry1[0]]
						newVec=np.ones(4)
						for i in range4:
							if i==i1:
								newVec[i]=1.0+mutMatrix[i][i]*(entry1[3]+bLen1)
							else:
								newVec[i]=mutMatrix[i][i1]*(entry1[3]+bLen1)
						if entry2[0]=="O":
							for i in range4:
								tot=0.0
								for j in range4:
									if i==j:
										tot+=(1.0+mutMatrix[i][i]*(entry2[3]+bLen2))*entry2[4][j]
									else:
										tot+=mutMatrix[i][j]*(entry2[3]+bLen2)*entry2[4][j]
								newVec[i]*=tot
							state, newVec, maxP =simplfy(newVec,ref[pos-1])
							probVect.append([state,pos,1,0.0,newVec])
							cumulPartLk+=math.log(maxP)
						else:
							if entry2[0]=="R":
								i2=allelesLow[ref[pos-1]]
							else:
								i2=alleles[entry2[0]]
							for i in range4:
								if i==i2:
									newVec[i]*=1.0+mutMatrix[i][i]*(entry2[3]+bLen2)
								else:
									newVec[i]*=mutMatrix[i][i2]*(entry2[3]+bLen2)
							probVect.append(["O",pos,1,0.0,newVec])

				#update pos, end, etc
				pos+=length
				if pos>lRef:
					break
				if pos>end1:
					indexEntry1+=1
					entry1=probVect1[indexEntry1]
					pos1=entry1[1]
					if entry1[0]!="N" and entry1[0]!="R":
						end1=pos1
					else:
						end1=pos1+entry1[2]-1
				if pos>end2:
					indexEntry2+=1
					entry2=probVect2[indexEntry2]
					pos2=entry2[1]
					if entry2[0]!="N" and entry2[0]!="R":
						end2=pos2
					else:
						end2=pos2+entry2[2]-1
				end=end1
				if end2<end:
					end=end2
				length=end+1-pos
			return probVect, cumulPartLk




#merge two child partial likelihood vectors at an internal node.
def mergeVectors(probVect1,bLen1,probVect2,bLen2,mutMatrix):
			probVect=[]
			indexEntry1=0
			indexEntry2=0
			pos=1
			entry1=probVect1[indexEntry1]
			pos1=entry1[1]
			if entry1[0]!="N" and entry1[0]!="R":
				end1=pos1
			else:
				end1=pos1+entry1[2]-1
			entry2=probVect2[indexEntry2]
			pos2=entry2[1]
			if entry2[0]!="N" and entry2[0]!="R":
				end2=pos2
			else:
				end2=pos2+entry2[2]-1
			end=end1
			if end2<end:
				end=end2
			length=end+1-pos
			while True:
				if entry1[0]=="N":
					probVect.append(list(entry2))
					if entry2[0]=="N" or entry2[0]=="R": 
						probVect[-1][2]=length
						probVect[-1][1]=pos
					probVect[-1][3]+=bLen2
					if entry2[0]=="O":
						probVect[-1][4]=list(entry2[4])
				elif entry2[0]=="N":
					probVect.append(list(entry1))
					if entry1[0]=="N" or entry1[0]=="R": 
						probVect[-1][2]=length
						probVect[-1][1]=pos
					probVect[-1][3]+=bLen1
					if entry1[0]=="O":
						probVect[-1][4]=list(entry1[4])
				elif bLen1+bLen2+entry1[3]+entry2[3]<thresholdProb2 and entry1[0]!=entry2[0] and entry1[0]!="O" and entry2[0]!="O":
					print("Warning, inside mergeVectors, branch lengths are 0, but the two vectors are different")
					print(probVect1)
					print(bLen1)
					print(probVect2)
					print(bLen2)
					print(entry1)
					print(entry2)
					return None
				elif entry2[0]=="R" and (entry2[3]+bLen2<thresholdProb2):
					probVect.append(["R",pos,length,0.0,False])
				elif (entry2[3]+bLen2<thresholdProb2) and entry2[0]!="O":
					probVect.append([entry2[0],pos,1,0.0,False])
				elif entry1[0]=="R":
					if entry2[0]=="R" or (entry1[3]+bLen1<thresholdProb2):
						probVect.append(["R",pos,length,0.0,False])
					elif entry2[0]=="O":
						i1=allelesLow[ref[pos-1]]
						newVec=np.ones(4)
						for i in range4:
							if i==i1:
								newVec[i]=1.0+mutMatrix[i][i]*(entry1[3]+bLen1)
							else:
								newVec[i]=mutMatrix[i][i1]*(entry1[3]+bLen1)
							tot=0.0
							for j in range4:
								if i==j:
									tot+=(1.0+mutMatrix[i][i]*(entry2[3]+bLen2))*entry2[4][j]
								else:
									tot+=mutMatrix[i][j]*(entry2[3]+bLen2)*entry2[4][j]
							newVec[i]*=tot
						state, newVec, maxP =simplfy(newVec,ref[pos-1])
						if state=="O":
							probVect.append([state,pos,1,0.0,newVec])
						else:
							probVect.append([state,pos,1,0.0,False])
					else:
						i2=alleles[entry2[0]]
						i1=allelesLow[ref[pos-1]]
						newVec=np.ones(4)
						for i in range4:
							if i==i2:
								newVec[i]=1.0+mutMatrix[i][i]*(entry2[3]+bLen2)
							else:
								newVec[i]=mutMatrix[i][i2]*(entry2[3]+bLen2)
							if i==i1:
								newVec[i]*=1.0+mutMatrix[i][i]*(entry1[3]+bLen1)
							else:
								newVec[i]*=mutMatrix[i][i1]*(entry1[3]+bLen1)
						probVect.append(["O",pos,1,0.0,newVec])
				elif entry1[0]=="O":
					newVec=np.ones(4)
					if entry2[0]=="O":
						#print("O")
						for i in range4:
							tot1=0.0
							tot2=0.0
							for j in range4:
								if i==j:
									tot1+=(1.0+mutMatrix[i][i]*(entry1[3]+bLen1))*entry1[4][j]
									tot2+=(1.0+mutMatrix[i][i]*(entry2[3]+bLen2))*entry2[4][j]
								else:
									tot1+=mutMatrix[i][j]*(entry1[3]+bLen1)*entry1[4][j]
									tot2+=mutMatrix[i][j]*(entry2[3]+bLen2)*entry2[4][j]
							newVec[i]*=tot1*tot2
					else:
						for i in range4:
							tot1=0.0
							for j in range4:
								if i==j:
									tot1+=(1.0+mutMatrix[i][i]*(entry1[3]+bLen1))*entry1[4][j]
								else:
									tot1+=mutMatrix[i][j]*(entry1[3]+bLen1)*entry1[4][j]
							newVec[i]=tot1
						if entry2[0]=="R":
							i2=allelesLow[ref[pos-1]]
							for i in range4:
								if i==i2:
									newVec[i]*=1.0+mutMatrix[i][i]*(entry2[3]+bLen2)
								else:
									newVec[i]*=mutMatrix[i][i2]*(entry2[3]+bLen2)
						else:
							i2=alleles[entry2[0]]
							for i in range4:
								if i==i2:
									newVec[i]*=1.0+mutMatrix[i][i]*(entry2[3]+bLen2)
								else:
									newVec[i]*=mutMatrix[i][i2]*(entry2[3]+bLen2)
					state, newVec, maxP =simplfy(newVec,ref[pos-1])
					if state=="O":
						probVect.append([state,pos,1,0.0,newVec])
					else:
						probVect.append([state,pos,1,0.0,False])
				#entry1 is a non-ref nuc
				else:
					if entry2[0]==entry1[0] or (entry1[3]+bLen1<thresholdProb2):
						probVect.append([entry1[0],pos,1,0.0,False])
					else:
						i1=alleles[entry1[0]]
						newVec=np.ones(4)
						for i in range4:
							if i==i1:
								newVec[i]=1.0+mutMatrix[i][i]*(entry1[3]+bLen1)
							else:
								newVec[i]=mutMatrix[i][i1]*(entry1[3]+bLen1)
						if entry2[0]=="O":
							for i in range4:
								tot=0.0
								for j in range4:
									if i==j:
										tot+=(1.0+mutMatrix[i][i]*(entry2[3]+bLen2))*entry2[4][j]
									else:
										tot+=mutMatrix[i][j]*(entry2[3]+bLen2)*entry2[4][j]
								newVec[i]*=tot
							state, newVec, maxP =simplfy(newVec,ref[pos-1])
							if state=="O":
								probVect.append([state,pos,1,0.0,newVec])
							else:
								probVect.append([state,pos,1,0.0,False])
						else:
							if entry2[0]=="R":
								i2=allelesLow[ref[pos-1]]
							else:
								i2=alleles[entry2[0]]
							for i in range4:
								if i==i2:
									newVec[i]*=1.0+mutMatrix[i][i]*(entry2[3]+bLen2)
								else:
									newVec[i]*=mutMatrix[i][i2]*(entry2[3]+bLen2)
							probVect.append(["O",pos,1,0.0,newVec])

				#update pos, end, etc
				pos+=length
				if pos>lRef:
					break
				if pos>end1:
					indexEntry1+=1
					entry1=probVect1[indexEntry1]
					pos1=entry1[1]
					if entry1[0]!="N" and entry1[0]!="R":
						end1=pos1
					else:
						end1=pos1+entry1[2]-1
				if pos>end2:
					indexEntry2+=1
					entry2=probVect2[indexEntry2]
					pos2=entry2[1]
					if entry2[0]!="N" and entry2[0]!="R":
						end2=pos2
					else:
						end2=pos2+entry2[2]-1
				end=end1
				if end2<end:
					end=end2
				length=end+1-pos

			if verbose:
				print("Merged children vectors ")
				print(probVect)
			probVect1=probVect

			#check if the final  probVect can be simplified by merging consecutive entries
			probVect1 =shorten(probVect1)
			if verbose:
				print("Shortened after merging children vectors ")
				print(probVect1)
			return probVect1





#Total probability of partial likelihood vector of the root after merging with root frequencies
def findProbRoot(probVect):
	logLK=0.0
	for entry in probVect:
			if entry[0]=="R":
				pos=entry[1]
				end=entry[2]+pos-1
				for i in range(4):
					logLK+=rootFreqsLog[i]*(cumulativeBases[end][i]-cumulativeBases[pos-1][i])
			elif entry[0]=="N":
				pass
			elif entry[0]=="O":
				tot=0.0
				for i in range(4):
					tot+=rootFreqs[i]*entry[4][i]
				logLK+=math.log(tot)
			else:
				i1=alleles[entry[0]]
				logLK+=rootFreqsLog[i1]
	return logLK




#Check if two genome lists represent the same partial likelihoods or not
def areVectorsDifferent(probVect1,probVect2):
	indexEntry1=0
	indexEntry2=0
	pos=1
	entry1=probVect1[indexEntry1]
	pos1=entry1[1]
	if entry1[0]!="N" and entry1[0]!="R":
		end1=pos1
	else:
		end1=pos1+entry1[2]-1
	entry2=probVect2[indexEntry2]
	pos2=entry2[1]
	if entry2[0]!="N" and entry2[0]!="R":
		end2=pos2
	else:
		end2=pos2+entry2[2]-1
	end=end1
	if end2<end:
		end=end2
	length=end+1-pos
	while True:
		if entry1[0]!=entry2[0]:
			if entry1[0]!="O":
				if entry2[0]!="O":
					return True
				else:
					if entry1[0]=="N":
						return True
					vSum=np.sum(entry2[4])
					vect=entry2[4]/vSum
					if not vSum>thresholdProb4:
						print("Problem? sum of partial likelihoods is very low or nan")
						print(probVect1)
						print(probVect2)
						exit()
					if entry1[0]=="R":
						i1=allelesLow[ref[pos-1]]
					else:
						i1=alleles[entry1[0]]
					if vect[i1]<0.999:
						return True
			else:
				if entry2[0]=="N":
					return True
				vSum=np.sum(entry1[4])
				if not vSum>thresholdProb4:
					print("Problem? sum of partial likelihoods is very low or nan")
					print(probVect1)
					print(probVect2)
					exit()
				vect=entry1[4]/vSum
				if entry2[0]=="R":
					i2=allelesLow[ref[pos-1]]
				else:
					i2=alleles[entry2[0]]
				if vect[i2]<0.999:
					return True

		elif entry1[0]=="O":
			vect1=entry1[4]
			vect2=entry2[4]
			maxP1=-1.0
			maxP2=-1.0
			for i in range4:
				if vect1[i]>maxP1:
					maxP1=vect1[i]
				if vect2[i]>maxP2:
					maxP2=vect2[i]
			for i in range4:
				if abs((vect1[i]/maxP1) - (vect2[i]/maxP2))>0.001:
					return True

		#update pos, end, etc
		pos+=length
		if pos>lRef:
			break
		if pos>end1:
			indexEntry1+=1
			entry1=probVect1[indexEntry1]
			pos1=entry1[1]
			if entry1[0]!="N" and entry1[0]!="R":
				end1=pos1
			else:
				end1=pos1+entry1[2]-1
		if pos>end2:
			indexEntry2+=1
			entry2=probVect2[indexEntry2]
			pos2=entry2[1]
			if entry2[0]!="N" and entry2[0]!="R":
				end2=pos2
			else:
				end2=pos2+entry2[2]-1
		end=end1
		if end2<end:
			end=end2
		length=end+1-pos
	return False


#if updating lk creates an inconsistency, this function can increase the bLen from 0 so that the inconsistency is resolved.
def updateBLen(node,childNum,bLen):
	cNode=node.children[childNum]
	if childNum==0:
		vectUp=node.probVectUpRight
	else:
		vectUp=node.probVectUpLeft
	bestLK=appendProb(vectUp,cNode.probVect,bLen)
	bestLen=bLen
	while bestLen>0.1*bLen:
		newBLen=bestLen/2
		newLK=appendProb(vectUp,cNode.probVect,newBLen)
		if newLK>bestLK:
			bestLK=newLK
			bestLen=newBLen
		else:
			break
	if bestLen>0.7*bLen:
		while bestLen<10*bLen:
			newBLen=bestLen*2
			newLK=appendProb(vectUp,cNode.probVect,newBLen)
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
	if node.dist>thresholdProb2:
		newTot=mergeVectorsUpDown(probVectUp,node.dist/2,node.probVect,node.dist/2,mutMatrix)
		if newTot==None:
			if node.up.children[0]==node:
				updateBLen(node.up,0,bLen)
			else:
				updateBLen(node.up,1,bLen)
			return
		newTot=shorten(newTot)
		node.probVectTotUp=newTot
	if len(node.children)==0:
		newTot=mergeVectorsUpDown(probVectUp,node.dist,node.probVect,0.0,mutMatrix)
		if newTot==None:
			if node.up.children[0]==node:
				updateBLen(node.up,0,bLen)
			else:
				updateBLen(node.up,1,bLen)
			return
		newTot=shorten(newTot)
		node.probVectTot=newTot
	else:
		child0Vect=node.children[0].probVect
		child1Vect=node.children[1].probVect
		dist0=node.children[0].dist
		dist1=node.children[1].dist
		newUpRight=mergeVectorsUpDown(probVectUp,node.dist,child1Vect,dist1,mutMatrix)
		newUpLeft=mergeVectorsUpDown(probVectUp,node.dist,child0Vect,dist0,mutMatrix)
		updatedTot=False
		if areVectorsDifferent(node.probVectUpRight,newUpRight):
			newUpRight =shorten(newUpRight)
			node.probVectUpRight=newUpRight
			newTot=mergeVectorsUpDown(newUpRight,0.0,child0Vect,dist0,mutMatrix)
			if newTot==None:
				updateBLen(node,0,bLen)
				return
			newTot =shorten(newTot)
			node.probVectTot=newTot
			updatedTot=True
			updatePartialsFromTop(node.children[0],newUpRight,mutMatrix)
		if areVectorsDifferent(node.probVectUpLeft,newUpLeft):
			newUpLeft =shorten(newUpLeft)
			node.probVectUpLeft=newUpLeft
			if not updatedTot:
				newTot=mergeVectorsUpDown(newUpLeft,0.0,child1Vect,dist1,mutMatrix)
				if newTot==None:
					updateBLen(node,1,bLen)
					return
				newTot =shorten(newTot)
				node.probVectTot=newTot
			updatePartialsFromTop(node.children[1],newUpLeft,mutMatrix)


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
	if verbose:
		print("Updating vectors from bottom, new vec")
		print(newVect)
		print("Old vec:")
		print(node.probVect)
	updatedTot=False
	if areVectorsDifferent(node.probVect,newVect):
		if verbose:
			print("Vectors are different")
		newVect =shorten(newVect)
		node.probVect=newVect
		if childNum==0:
			newTot=mergeVectorsUpDown(node.probVectUpRight,0.0,probVectDown,childDist,mutMatrix)
			if newTot==None:
				updateBLen(node,0,bLen)
				return
		else:
			newTot=mergeVectorsUpDown(node.probVectUpLeft,0.0,probVectDown,childDist,mutMatrix)
			if newTot==None:
				updateBLen(node,1,bLen)
				return
		newTot=shorten(newTot)
		node.probVectTot=newTot
		if verbose:
			print("new tot vect")
			print(node.probVectTot)
		updatedTot=True
		if node.up != None:
			updatePartialsFromBottom(node.up,newVect,-1,node,mutMatrix)
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
			if areVectorsDifferent(node.probVectUpLeft,newUpLeftVect):
				if verbose:
					print("Moving up the update from left; new UpLeft and old UpLeft")
					print(newUpLeftVect)
					print(node.probVectUpLeft)
				newUpLeftVect =shorten(newUpLeftVect)
				node.probVectUpLeft=newUpLeftVect
				if not updatedTot:
					if verbose:
						print("Updating tot: old tot and new tot")
						print(node.probVectTot)
					newTot=mergeVectorsUpDown(newUpLeftVect,0.0,otherChildVect,otherChildDist,mutMatrix)
					if newTot==None:
						updateBLen(node,1,bLen)
						return
					newTot=shorten(newTot)
					node.probVectTot=newTot
					if verbose:
						print(newTot)
				
				updatePartialsFromTop(node.children[1],newUpLeftVect,mutMatrix)
		else:
			newUpRightVect=mergeVectorsUpDown(vectUp,node.dist,probVectDown,childDist,mutMatrix)
			if areVectorsDifferent(node.probVectUpRight,newUpRightVect):
				if verbose:
					print("Moving up the update from right")
				newUpRightVect =shorten(newUpRightVect)
				node.probVectUpRight=newUpRightVect
				if not updatedTot:
					newTot=mergeVectorsUpDown(newUpRightVect,0.0,otherChildVect,otherChildDist,mutMatrix)
					if newTot==None:
						updateBLen(node,0,bLen)
						return
					newTot=shorten(newTot)
					node.probVectTot=newTot
				updatePartialsFromTop(node.children[0],newUpRightVect,mutMatrix)
		if updatedTot:
			newTot=mergeVectorsUpDown(vectUp,node.dist/2,node.probVect,node.dist/2,mutMatrix)
			newTot=shorten(newTot)
			node.probVectTotUp=newTot
	#case node is root
	else:
		if verbose:
			print("reached root while moving up")
		if childNum==0:
			newUpLeftVect=rootVector(probVectDown,rootFreqs,childDist)
			if areVectorsDifferent(node.probVectUpLeft,newUpLeftVect):
				newUpLeftVect =shorten(newUpLeftVect)
				node.probVectUpLeft=newUpLeftVect
				if not updatedTot:
					newTot=mergeVectorsUpDown(newUpLeftVect,0.0,otherChildVect,otherChildDist,mutMatrix)
					newTot=shorten(newTot)
					node.probVectTot=newTot
				updatePartialsFromTop(node.children[1],newUpLeftVect,mutMatrix)

		else:
			newUpRightVect=rootVector(probVectDown,rootFreqs,childDist)
			if areVectorsDifferent(node.probVectUpRight,newUpRightVect):
				newUpRightVect =shorten(newUpRightVect)
				node.probVectUpRight=newUpRightVect
				if not updatedTot:
					newTot=mergeVectorsUpDown(newUpRightVect,0.0,otherChildVect,otherChildDist,mutMatrix)
					node.probVectTot=newTot
					newTot=shorten(newTot)
				updatePartialsFromTop(node.children[0],newUpRightVect,mutMatrix)
	if verbose:
		print("New upRight, UpLeft and tot:")
		print(node.probVectUpRight)
		print(node.probVectUpLeft)
		print(node.probVectTot)
	return

#try all children branches (if meets a 0-length branch, pass the function to child node) to find the best branch where to place (mid-length) a new sample.
#This basically takes into account polytomies (here politomies are given by 0-length branches in a binary tree).
def findBestChildBranch(node,newPartials,sample,bLen):
	bestLKSoFar=float("-inf")
	bestNodeSoFar=None
	bestVectSoFar=None
	for c in range(len(node.children)):
		child=node.children[c]
		if child.dist>thresholdProb2:
			if c==0:
				probVectChild=mergeVectorsUpDown(node.probVectUpRight,child.dist/2,child.probVect,child.dist/2,mutMatrix)
			else:
				probVectChild=mergeVectorsUpDown(node.probVectUpLeft,child.dist/2,child.probVect,child.dist/2,mutMatrix)
			probChild=appendProb(probVectChild,newPartials,bLen)
			if probChild>bestLKSoFar:
				bestLKSoFar=probChild
				bestNodeSoFar=child
				bestVectSoFar=probVectChild
		else:
			childBestNode, childBestLK, childBestVect=findBestChildBranch(child,newPartials,sample,bLen)
			if childBestLK>bestLKSoFar:
				bestLKSoFar=childBestLK
				bestNodeSoFar=childBestNode
				bestVectSoFar=childBestVect
	return bestNodeSoFar, bestLKSoFar, bestVectSoFar
	


#we know that sample "sample", with partials "newPartials", is best placed as child of node resulting in logLK contribution of newChildLK
# now explore exactly which position near node is best for placement (direct child, sibling or sibling of existing child), and which are the best branch lengths,
#then add the sample at that position of the tree, and update all the internal probability vectors.
#NEW VERSION: now not going through the children or parent, but just attaching to the best place found so far.
def placeSampleOnTreeNew(node,newPartials,sample,bLen,newChildLK,isMidNode):
	#if node is root, try to place as sibling of the current root.
	if node.is_root():
		probOldRoot = findProbRoot(node.probVect)
		probVectRoot,probRoot = mergeVectorsRoot(node.probVect,bLen,newPartials,bLen,mutMatrix)
		if verbose:
			print("parent is root, trying to add new root; node.up, node.probvect, root old prob, new root vect, new root prob")
			print(node.up)
			print(node.probVect)
			print(probOldRoot)
			print(probVectRoot)
			print(probRoot)
		probRoot+= findProbRoot(probVectRoot)
		parentLKdiff=probRoot-probOldRoot
		if verbose:
			print(" new root prob after adding root freqs; difference between old root prob and new")
			print(probRoot)
			print(parentLKdiff)
	else:
		parentLKdiff=float("-inf")
	
	#append as direct descendant of node
	if newChildLK>=parentLKdiff and (not isMidNode):
		#now try different lengths for the new branch
		LK1=newChildLK
		bestLen=bLen
		while bestLen>0.1*bLen:
			newBLen=bestLen/2
			probChild=appendProb(node.probVectTot,newPartials,newBLen)
			if probChild>LK1:
				LK1=probChild
				bestLen=newBLen
			else:
				break
		if bestLen>0.7*bLen:
			while bestLen<10*bLen:
				newBLen=bestLen*2
				probChild=appendProb(node.probVectTot,newPartials,newBLen)
				if probChild>LK1:
					LK1=probChild
					bestLen=newBLen
				else:
					break
		if bestLen<0.2*bLen:
			LK0=appendProb(node.probVectTot,newPartials,0.0)
			if LK0>LK1:
				bestLen=0.0

		#add parent to the root, but coinciding with the current root
		if node.is_root():
			newRoot=Tree()
			newRoot.add_child(child=node,dist=0.0)
			node.up=newRoot
			newRoot.add_child(dist=bestLen, name=sample)
			newRoot.children[1].up=newRoot
			newRoot.probVectUpLeft=node.probVectTot
			newRoot.probVectUpRight=rootVector(newPartials,rootFreqs,bestLen)
			newRoot.probVect=mergeVectors(node.probVect,0.0,newPartials,bestLen,mutMatrix)
			newVect=rootVector0(newRoot.probVect,rootFreqs)
			newVect=shorten(newVect)
			newRoot.probVectTot=newVect
			newVect=mergeVectorsUpDown(newRoot.probVectUpLeft,bestLen,newPartials,0.0,mutMatrix)
			newVect=shorten(newVect)
			newRoot.children[1].probVectTot=newVect
			newRoot.children[1].minorSequences=[]
			newRoot.children[1].probVect=newPartials
			newVect=mergeVectorsUpDown(newRoot.probVectUpLeft,bestLen/2,newPartials,bestLen/2,mutMatrix)
			newVect=shorten(newVect)
			newRoot.children[1].probVectTotUp=newVect
			if verbose:
				print("new root added to tree, identical to previous root")
				print(newRoot.probVect)
				print(newRoot.children[0].probVect)
				print(newRoot.children[1].probVect)
			updatePartialsFromTop(node,newRoot.probVectUpRight,mutMatrix)
			return newRoot

		#add parent to node, but coinciding with it
		else:
			if node==node.up.children[0]:
				child=0
				vectUp=node.up.probVectUpRight
			else:
				child=1
				vectUp=node.up.probVectUpLeft
			#now create new internal node and append child to it
			newInternalNode=Tree()
			node.up.children[child]=newInternalNode
			newInternalNode.up=node.up
			newInternalNode.dist=node.dist
			newInternalNode.add_child(dist=0.0,child=node)
			node.up=newInternalNode
			newInternalNode.add_child(dist=bestLen, name=sample)
			newInternalNode.children[1].minorSequences=[]
			newInternalNode.children[1].up=newInternalNode
			newInternalNode.children[1].probVect=newPartials
			newVect=mergeVectorsUpDown(vectUp,newInternalNode.dist,newPartials,bestLen,mutMatrix)
			#newVect=shorten(newVect)
			newInternalNode.probVectUpRight=newVect
			newInternalNode.probVectUpLeft=node.probVectTot
			newVect=mergeVectors(node.probVect,0.0,newPartials,bestLen,mutMatrix)
			#newVect=shorten(newVect)
			newInternalNode.probVect=newVect
			newVect=mergeVectorsUpDown(newInternalNode.probVectUpLeft,0.0,newPartials,bestLen,mutMatrix)
			newVect=shorten(newVect)
			newInternalNode.probVectTot=newVect
			newVect=mergeVectorsUpDown(vectUp,newInternalNode.dist/2,newInternalNode.probVect,newInternalNode.dist/2,mutMatrix)
			newVect=shorten(newVect)
			newInternalNode.probVectTotUp=newVect
			newVect=mergeVectorsUpDown(newInternalNode.probVectUpLeft,bestLen,newPartials,0.0,mutMatrix)
			newVect=shorten(newVect)
			newInternalNode.children[1].probVectTot=newVect
			newVect=mergeVectorsUpDown(newInternalNode.probVectUpLeft,bestLen/2,newPartials,bestLen/2,mutMatrix)
			newVect=shorten(newVect)
			newInternalNode.children[1].probVectTotUp=newVect
			if verbose:
				print("new internal node and child node added to tree")
				print(newInternalNode.probVect)
				print(newInternalNode.probVectUpRight)
				print(newInternalNode.probVectUpLeft)
				print(newInternalNode.probVectTot)
			updatePartialsFromTop(node,newInternalNode.probVectUpRight,mutMatrix)
			updatePartialsFromBottom(newInternalNode.up,newInternalNode.probVect,child,newInternalNode,mutMatrix)

	else: #add a new parent to node
		#add a new root!
		if node.is_root():
			bestLen1=bLen
			bestLK1=parentLKdiff
			rootBestVect=probVectRoot
			while bestLen1>0.1*bLen:
				newBLen=bestLen1/2
				newProbVectRoot,newProbRoot = mergeVectorsRoot(node.probVect,newBLen,newPartials,bLen,mutMatrix)
				newProbRoot+= findProbRoot(newProbVectRoot)
				LKdiffRoot=newProbRoot-probOldRoot
				if LKdiffRoot>bestLK1:
					bestLK1=LKdiffRoot
					bestLen1=newBLen
					rootBestVect=newProbVectRoot
				else:
					break
			if bestLen1>0.7*bLen:
				while bestLen1<10*bLen:
					newBLen=bestLen1*2
					newProbVectRoot,newProbRoot = mergeVectorsRoot(node.probVect,newBLen,newPartials,bLen,mutMatrix)
					newProbRoot+= findProbRoot(newProbVectRoot)
					LKdiffRoot=newProbRoot-probOldRoot
					if LKdiffRoot>bestLK1:
						bestLK1=LKdiffRoot
						bestLen1=newBLen
						rootBestVect=newProbVectRoot
					else:
						break
			#now try different lengths for right branch
			bestLen2=bLen
			while bestLen2>0.1*bLen:
				newBLen=bestLen2/2
				newProbVectRoot,newProbRoot = mergeVectorsRoot(node.probVect,bestLen1,newPartials,newBLen,mutMatrix)
				newProbRoot+= findProbRoot(newProbVectRoot)
				LKdiffRoot=newProbRoot-probOldRoot
				if LKdiffRoot>bestLK1:
					bestLK1=LKdiffRoot
					bestLen2=newBLen
					rootBestVect=newProbVectRoot
				else:
					break
			if bestLen2>0.7*bLen:
				while bestLen2<10*bLen:
					newBLen=bestLen2*2
					newProbVectRoot,newProbRoot = mergeVectorsRoot(node.probVect,bestLen1,newPartials,newBLen,mutMatrix)
					newProbRoot+= findProbRoot(newProbVectRoot)
					LKdiffRoot=newProbRoot-probOldRoot
					if LKdiffRoot>bestLK1:
						bestLK1=LKdiffRoot
						bestLen2=newBLen
						rootBestVect=newProbVectRoot
					else:
						break
			if bestLen2<0.2*bLen:
				newProbVectRoot,newProbRoot = mergeVectorsRoot(node.probVect,bestLen1,newPartials,0.0,mutMatrix)
				newProbRoot+= findProbRoot(newProbVectRoot)
				LK0=newProbRoot-probOldRoot
				if LK0>bestLK1:
					bestLen2=0.0
					bestLK1=LK0
					rootBestVect=newProbVectRoot
			#now creating new root
			newRoot=Tree()
			newRoot.probVect=rootBestVect
			newVect=rootVector0(rootBestVect,rootFreqs)
			newVect=shorten(newVect)
			newRoot.probVectTot=newVect
			newRoot.probVectUpRight=rootVector(newPartials,rootFreqs,bestLen2)
			newRoot.probVectUpLeft=rootVector(node.probVect,rootFreqs,bestLen1)
			newRoot.add_child(child=node,dist=bestLen1)
			node.up=newRoot
			newRoot.add_child(dist=bestLen2, name=sample)
			newRoot.children[1].probVect=newPartials
			newVect=mergeVectorsUpDown(newRoot.probVectUpLeft,bestLen2,newPartials,0.0,mutMatrix)
			newVect=shorten(newVect)
			newRoot.children[1].probVectTot=newVect
			newVect=mergeVectorsUpDown(newRoot.probVectUpLeft,bestLen2/2,newPartials,bestLen2/2,mutMatrix)
			newVect=shorten(newVect)
			newRoot.children[1].probVectTotUp=newVect
			newRoot.children[1].minorSequences=[]
			newRoot.children[1].up=newRoot
			if verbose:
				print("new root added to tree")
				print(newRoot.probVect)
				print(newRoot.children[0].probVect)
				print(newRoot.children[1].probVect)
			updatePartialsFromTop(node,newRoot.probVectUpRight,mutMatrix)
			return newRoot
		
		#add parent to node, knowing that node is not the root
		else:
			if node==node.up.children[0]:
				child=0
				vectUp=node.up.probVectUpRight
			else:
				child=1
				vectUp=node.up.probVectUpLeft
			bestSplit=0.5
			bestSplitLK=parentLKdiff
			childBestVect=node.probVectTotUp
			#try different positions on the existing branch
			while bestSplit>0.05:
				newSplit=bestSplit/2
				probVectParentNew=mergeVectorsUpDown(vectUp,node.dist*newSplit,node.probVect,node.dist*(1.0-newSplit),mutMatrix)
				probChild=appendProb(probVectParentNew,newPartials,bLen)
				if probChild>bestSplitLK:
					bestSplitLK=probChild
					bestSplit=newSplit
					childBestVect=probVectParentNew
				else:
					break
			if bestSplit>0.49:
				#print("Now trying the reverse direction")
				while bestSplit>0.05:
					newSplit=bestSplit/2
					probVectParentNew=mergeVectorsUpDown(vectUp,node.dist*(1.0-newSplit),node.probVect,node.dist*newSplit,mutMatrix)
					probChild=appendProb(probVectParentNew,newPartials,bLen)
					if probChild>bestSplitLK:
						bestSplitLK=probChild
						bestSplit=newSplit
						childBestVect=probVectParentNew
					else:
						bestSplit=1.0-bestSplit
						break
			#now try different lengths for the new branch
			childBestVect=shorten(childBestVect)
			LK1=bestSplitLK
			bestLen=bLen
			while bestLen>0.1*bLen:
				newBLen=bestLen/2
				probChild=appendProb(childBestVect,newPartials,newBLen)
				if probChild>LK1:
					LK1=probChild
					bestLen=newBLen
				else:
					break
			if bestLen>0.7*bLen:
				while bestLen<10*bLen:
					newBLen=bestLen*2
					probChild=appendProb(childBestVect,newPartials,newBLen)
					if probChild>LK1:
						LK1=probChild
						bestLen=newBLen
					else:
						break
			if bestLen<0.2*bLen:
				LK0=appendProb(childBestVect,newPartials,0.0)
				if LK0>LK1:
					bestLen=0.0
			#now create new internal node and append child to it
			newInternalNode=Tree()
			node.up.children[child]=newInternalNode
			newInternalNode.up=node.up
			distBottom=node.dist*(1.0-bestSplit)
			distTop=node.dist*bestSplit
			newInternalNode.add_child(dist=distBottom,child=node)
			node.up=newInternalNode
			newInternalNode.add_child(dist=bestLen, name=sample)
			newInternalNode.children[1].minorSequences=[]
			newInternalNode.children[1].up=newInternalNode
			newInternalNode.dist=distTop
			newInternalNode.children[1].probVect=newPartials
			newVect=mergeVectorsUpDown(vectUp,newInternalNode.dist,newPartials,bestLen,mutMatrix)
			#newVect=shorten(newVect)
			newInternalNode.probVectUpRight=newVect
			newInternalNode.probVectUpLeft=childBestVect
			newVect=mergeVectors(node.probVect,node.dist,newPartials,bestLen,mutMatrix)
			#newVect=shorten(newVect)
			newInternalNode.probVect=newVect
			newVect=mergeVectorsUpDown(vectUp,distTop/2,newInternalNode.probVect,distTop/2,mutMatrix)
			newVect=shorten(newVect)
			newInternalNode.probVectTotUp=newVect
			newVect=mergeVectorsUpDown(childBestVect,0.0,newPartials,bestLen,mutMatrix)
			newVect=shorten(newVect)
			newInternalNode.probVectTot=newVect
			newVect=mergeVectorsUpDown(childBestVect,bestLen,newPartials,0.0,mutMatrix)
			newVect=shorten(newVect)
			newInternalNode.children[1].probVectTot=newVect
			newVect=mergeVectorsUpDown(childBestVect,bestLen/2,newPartials,bestLen/2,mutMatrix)
			newVect=shorten(newVect)
			newInternalNode.children[1].probVectTotUp=newVect
			if verbose:
				print("new internal node added to tree")
				print(newInternalNode.probVect)
				print(newInternalNode.probVectUpRight)
				print(newInternalNode.probVectUpLeft)
				print(newInternalNode.probVectTot)
			updatePartialsFromTop(node,newInternalNode.probVectUpRight,mutMatrix)
			updatePartialsFromBottom(newInternalNode.up,newInternalNode.probVect,child,newInternalNode,mutMatrix)

	return None





#we know that sample "sample", with partials "newPartials", is best placed as child of node resulting in logLK contribution of newChildLK
# now explore exactly which position near node is best for placement (direct child, sibling or sibling of existing child), and which are the best branch lengths,
#then add the sample at that position of the tree, and update all the internal probability vectors.
def placeSampleOnTree(node,newPartials,sample,bLen,newChildLK):
	if len(node.children)>0:
		childBestNode, childBestLK, childBestVect=findBestChildBranch(node,newPartials,sample,bLen)
	else:
		childBestLK=float("-inf")
		childBestNode=None
		childBestVect=None
	#if node is root, try to place as sibling of the current root.
	if node.is_root():
		probOldRoot = findProbRoot(node.probVect)
		probVectRoot,probRoot = mergeVectorsRoot(node.probVect,bLen,newPartials,bLen,mutMatrix)
		if verbose:
			print("parent is root, trying to add new root; node.up, node.probvect, root old prob, new root vect, new root prob")
			print(node.up)
			print(node.probVect)
			print(probOldRoot)
			print(probVectRoot)
			print(probRoot)
		probRoot+= findProbRoot(probVectRoot)
		parentLKdiff=probRoot-probOldRoot
		if verbose:
			print(" new root prob after adding root freqs; difference between old root prob and new")
			print(probRoot)
			print(parentLKdiff)
	else:
		if node==node.up.children[0]:
			probVectParent=mergeVectorsUpDown(node.up.probVectUpRight,node.dist/2,node.probVect,node.dist/2,mutMatrix)
		else:
			probVectParent=mergeVectorsUpDown(node.up.probVectUpLeft,node.dist/2,node.probVect,node.dist/2,mutMatrix)
		parentLKdiff=appendProb(probVectParent,newPartials,bLen)
		if verbose:
			print(" non-root parent placement prob and vector:")
			print(probVectParent)
			print(parentLKdiff)
	
	#append as direct descendant of node
	if newChildLK>=parentLKdiff and newChildLK>=childBestLK:
		#now try different lengths for the new branch
		LK1=newChildLK
		bestLen=bLen
		while bestLen>0.1*bLen:
			newBLen=bestLen/2
			probChild=appendProb(node.probVectTot,newPartials,newBLen)
			if probChild>LK1:
				LK1=probChild
				bestLen=newBLen
			else:
				break
		if bestLen>0.7*bLen:
			while bestLen<10*bLen:
				newBLen=bestLen*2
				probChild=appendProb(node.probVectTot,newPartials,newBLen)
				if probChild>LK1:
					LK1=probChild
					bestLen=newBLen
				else:
					break
		LK0=appendProb(node.probVectTot,newPartials,0.0)
		if LK0>LK1:
			bestLen=0.0

		#add parent to the root, but coinciding with the current root
		if node.is_root():
			newRoot=Tree()
			newRoot.add_child(child=node,dist=0.0)
			node.up=newRoot
			newRoot.add_child(dist=bestLen, name=sample)
			newRoot.children[1].up=newRoot
			newRoot.probVectUpLeft=node.probVectTot
			newRoot.probVectUpRight=rootVector(newPartials,rootFreqs,bestLen)
			newRoot.probVect=mergeVectors(node.probVect,0.0,newPartials,bestLen,mutMatrix)
			newRoot.probVectTot=rootVector0(newRoot.probVect,rootFreqs)
			newRoot.children[1].probVectTot=mergeVectorsUpDown(newRoot.probVectUpLeft,bestLen,newPartials,0.0,mutMatrix)
			newRoot.children[1].minorSequences=[]
			newRoot.children[1].probVect=newPartials
			if verbose:
				print("new root added to tree, identical to previous root")
				print(newRoot.probVect)
				print(newRoot.children[0].probVect)
				print(newRoot.children[1].probVect)
			updatePartialsFromTop(node,newRoot.probVectUpRight,mutMatrix)
			return newRoot

		#add parent to node, but coinciding with it
		else:
			if node==node.up.children[0]:
				child=0
				vectUp=node.up.probVectUpRight
			else:
				child=1
				vectUp=node.up.probVectUpLeft
			#now create new internal node and append child to it
			newInternalNode=Tree()
			node.up.children[child]=newInternalNode
			newInternalNode.up=node.up
			newInternalNode.dist=node.dist
			newInternalNode.add_child(dist=0.0,child=node)
			node.up=newInternalNode
			newInternalNode.add_child(dist=bestLen, name=sample)
			newInternalNode.children[1].minorSequences=[]
			newInternalNode.children[1].up=newInternalNode
			newInternalNode.children[1].probVect=newPartials
			newVect=mergeVectorsUpDown(vectUp,newInternalNode.dist,newPartials,bestLen,mutMatrix)
			newVect=shorten(newVect)
			newInternalNode.probVectUpRight=newVect
			newInternalNode.probVectUpLeft=node.probVectTot
			newVect=mergeVectors(node.probVect,0.0,newPartials,bestLen,mutMatrix)
			newVect=shorten(newVect)
			newInternalNode.probVect=newVect
			newVect=mergeVectorsUpDown(newInternalNode.probVectUpLeft,0.0,newPartials,bestLen,mutMatrix)
			newVect=shorten(newVect)
			newInternalNode.probVectTot=newVect
			newVect=mergeVectorsUpDown(newInternalNode.probVectUpLeft,bestLen,newPartials,0.0,mutMatrix)
			newVect=shorten(newVect)
			newInternalNode.children[1].probVectTot=newVect
			if verbose:
				print("new internal node and child node added to tree")
				print(newInternalNode.probVect)
				print(newInternalNode.probVectUpRight)
				print(newInternalNode.probVectUpLeft)
				print(newInternalNode.probVectTot)
			updatePartialsFromTop(node,newInternalNode.probVectUpRight,mutMatrix)
			updatePartialsFromBottom(newInternalNode.up,newInternalNode.probVect,child,newInternalNode,mutMatrix)

	#put new sample as descendant of a new internal node on one of the children branches
	elif childBestLK>=parentLKdiff: 
		if childBestNode==childBestNode.up.children[0]:
			child=0
			vectUp=childBestNode.up.probVectUpRight
		else:
			child=1
			vectUp=childBestNode.up.probVectUpLeft
		bestSplit=0.5
		bestSplitLK=childBestLK
		#try different positions on the existing branch
		while bestSplit>0.01:
			newSplit=bestSplit/2
			probVectChildNew=mergeVectorsUpDown(vectUp,childBestNode.dist*newSplit,childBestNode.probVect,childBestNode.dist*(1.0-newSplit),mutMatrix)
			probChild=appendProb(probVectChildNew,newPartials,bLen)
			if probChild>bestSplitLK:
				bestSplitLK=probChild
				bestSplit=newSplit
				childBestVect=probVectChildNew
			else:
				break
		if bestSplit>0.49:
			while bestSplit>0.01:
				newSplit=bestSplit/2
				probVectChildNew=mergeVectorsUpDown(vectUp,childBestNode.dist*(1.0-newSplit),childBestNode.probVect,childBestNode.dist*newSplit,mutMatrix)
				probChild=appendProb(probVectChildNew,newPartials,bLen)
				if probChild>bestSplitLK:
					bestSplitLK=probChild
					bestSplit=newSplit
					childBestVect=probVectChildNew
				else:
					bestSplit=1.0-bestSplit
					break
		#now try different lengths for the new branch
		childBestVect=shorten(childBestVect)
		LK1=bestSplitLK
		bestLen=bLen
		while bestLen>0.1*bLen:
			newBLen=bestLen/2
			probChild=appendProb(childBestVect,newPartials,newBLen)
			if probChild>LK1:
				LK1=probChild
				bestLen=newBLen
			else:
				break
		if bestLen>0.7*bLen:
			while bestLen<10*bLen:
				newBLen=bestLen*2
				probChild=appendProb(childBestVect,newPartials,newBLen)
				if probChild>LK1:
					LK1=probChild
					bestLen=newBLen
				else:
					break
		LK0=appendProb(childBestVect,newPartials,0.0)
		if LK0>LK1:
			bestLen=0.0
		#now create new internal node and append child to it
		newInternalNode=Tree()
		childBestNode.up.children[child]=newInternalNode
		newInternalNode.up=childBestNode.up
		distBottom=childBestNode.dist*(1.0-bestSplit)
		distTop=childBestNode.dist*bestSplit
		newInternalNode.add_child(dist=distBottom,child=childBestNode)
		childBestNode.up=newInternalNode
		newInternalNode.add_child(dist=bestLen, name=sample)
		newInternalNode.children[1].minorSequences=[]
		newInternalNode.children[1].up=newInternalNode
		newInternalNode.dist=distTop
		newInternalNode.children[1].probVect=newPartials
		newVect=mergeVectorsUpDown(vectUp,newInternalNode.dist,newPartials,bestLen,mutMatrix)
		newVect=shorten(newVect)
		newInternalNode.probVectUpRight=newVect
		newInternalNode.probVectUpLeft=childBestVect
		newVect=mergeVectors(childBestNode.probVect,childBestNode.dist,newPartials,bestLen,mutMatrix)
		newVect=shorten(newVect)
		newInternalNode.probVect=newVect
		newVect=mergeVectorsUpDown(childBestVect,0.0,newPartials,bestLen,mutMatrix)
		newVect=shorten(newVect)
		newInternalNode.probVectTot=newVect
		newVect=mergeVectorsUpDown(childBestVect,bestLen,newPartials,0.0,mutMatrix)
		newVect=shorten(newVect)
		newInternalNode.children[1].probVectTot=newVect
		if verbose:
			print("new internal node added to tree")
			print(newInternalNode.probVect)
			print(newInternalNode.probVectUpRight)
			print(newInternalNode.probVectUpLeft)
			print(newInternalNode.probVectTot)
		updatePartialsFromTop(childBestNode,newInternalNode.probVectUpRight,mutMatrix)
		updatePartialsFromBottom(newInternalNode.up,newInternalNode.probVect,child,newInternalNode,mutMatrix)
	
	else: #add a new parent to node
		#add a new root!
		if node.is_root():
			bestLen1=bLen
			bestLK1=parentLKdiff
			rootBestVect=probVectRoot
			while bestLen1>0.1*bLen:
				newBLen=bestLen1/2
				newProbVectRoot,newProbRoot = mergeVectorsRoot(node.probVect,newBLen,newPartials,bLen,mutMatrix)
				newProbRoot+= findProbRoot(newProbVectRoot)
				LKdiffRoot=newProbRoot-probOldRoot
				if LKdiffRoot>bestLK1:
					bestLK1=LKdiffRoot
					bestLen1=newBLen
					rootBestVect=newProbVectRoot
				else:
					break
			if bestLen1>0.7*bLen:
				while bestLen1<10*bLen:
					newBLen=bestLen1*2
					newProbVectRoot,newProbRoot = mergeVectorsRoot(node.probVect,newBLen,newPartials,bLen,mutMatrix)
					newProbRoot+= findProbRoot(newProbVectRoot)
					LKdiffRoot=newProbRoot-probOldRoot
					if LKdiffRoot>bestLK1:
						bestLK1=LKdiffRoot
						bestLen1=newBLen
						rootBestVect=newProbVectRoot
					else:
						break
			#now try different lengths for right branch
			bestLen2=bLen
			while bestLen2>0.1*bLen:
				newBLen=bestLen2/2
				newProbVectRoot,newProbRoot = mergeVectorsRoot(node.probVect,bestLen1,newPartials,newBLen,mutMatrix)
				newProbRoot+= findProbRoot(newProbVectRoot)
				LKdiffRoot=newProbRoot-probOldRoot
				if LKdiffRoot>bestLK1:
					bestLK1=LKdiffRoot
					bestLen2=newBLen
					rootBestVect=newProbVectRoot
				else:
					break
			if bestLen2>0.7*bLen:
				while bestLen2<10*bLen:
					newBLen=bestLen2*2
					newProbVectRoot,newProbRoot = mergeVectorsRoot(node.probVect,bestLen1,newPartials,newBLen,mutMatrix)
					newProbRoot+= findProbRoot(newProbVectRoot)
					LKdiffRoot=newProbRoot-probOldRoot
					if LKdiffRoot>bestLK1:
						bestLK1=LKdiffRoot
						bestLen2=newBLen
						rootBestVect=newProbVectRoot
					else:
						break
			newProbVectRoot,newProbRoot = mergeVectorsRoot(node.probVect,bestLen1,newPartials,0.0,mutMatrix)
			newProbRoot+= findProbRoot(newProbVectRoot)
			LK0=newProbRoot-probOldRoot
			if LK0>bestLK1:
				bestLen2=0.0
				bestLK1=LK0
				rootBestVect=newProbVectRoot
			#now creating new root
			newRoot=Tree()
			newRoot.probVect=rootBestVect
			newRoot.probVectTot=rootVector0(rootBestVect,rootFreqs)
			newRoot.probVectUpRight=rootVector(newPartials,rootFreqs,bestLen2)
			newRoot.probVectUpLeft=rootVector(node.probVect,rootFreqs,bestLen1)
			newRoot.add_child(child=node,dist=bestLen1)
			node.up=newRoot
			newRoot.add_child(dist=bestLen2, name=sample)
			newRoot.children[1].probVect=newPartials
			newRoot.children[1].probVectTot=mergeVectorsUpDown(newRoot.probVectUpLeft,bestLen2,newPartials,0.0,mutMatrix)
			newRoot.children[1].minorSequences=[]
			newRoot.children[1].up=newRoot
			if verbose:
				print("new root added to tree")
				print(newRoot.probVect)
				print(newRoot.children[0].probVect)
				print(newRoot.children[1].probVect)
			updatePartialsFromTop(node,newRoot.probVectUpRight,mutMatrix)
			return newRoot
		
		#add parent to node, knowing that node is not the root
		else:
			if node==node.up.children[0]:
				child=0
				vectUp=node.up.probVectUpRight
			else:
				child=1
				vectUp=node.up.probVectUpLeft
			bestSplit=0.5
			bestSplitLK=parentLKdiff
			childBestVect=probVectParent
			#try different positions on the existing branch
			while bestSplit>0.01:
				newSplit=bestSplit/2
				probVectParentNew=mergeVectorsUpDown(vectUp,node.dist*newSplit,node.probVect,node.dist*(1.0-newSplit),mutMatrix)
				probChild=appendProb(probVectParentNew,newPartials,bLen)
				if probChild>bestSplitLK:
					bestSplitLK=probChild
					bestSplit=newSplit
					childBestVect=probVectParentNew
				else:
					break
			if bestSplit>0.49:
				#print("Now trying the reverse direction")
				while bestSplit>0.01:
					newSplit=bestSplit/2
					probVectParentNew=mergeVectorsUpDown(vectUp,node.dist*(1.0-newSplit),node.probVect,node.dist*newSplit,mutMatrix)
					probChild=appendProb(probVectParentNew,newPartials,bLen)
					if probChild>bestSplitLK:
						bestSplitLK=probChild
						bestSplit=newSplit
						childBestVect=probVectParentNew
					else:
						bestSplit=1.0-bestSplit
						break
			#now try different lengths for the new branch
			childBestVect=shorten(childBestVect)
			LK1=bestSplitLK
			bestLen=bLen
			while bestLen>0.1*bLen:
				newBLen=bestLen/2
				probChild=appendProb(childBestVect,newPartials,newBLen)
				if probChild>LK1:
					LK1=probChild
					bestLen=newBLen
				else:
					break
			if bestLen>0.7*bLen:
				while bestLen<10*bLen:
					newBLen=bestLen*2
					probChild=appendProb(childBestVect,newPartials,newBLen)
					if probChild>LK1:
						LK1=probChild
						bestLen=newBLen
					else:
						break
			LK0=appendProb(childBestVect,newPartials,0.0)
			if LK0>LK1:
				bestLen=0.0
			#now create new internal node and append child to it
			newInternalNode=Tree()
			node.up.children[child]=newInternalNode
			newInternalNode.up=node.up
			distBottom=node.dist*(1.0-bestSplit)
			distTop=node.dist*bestSplit
			newInternalNode.add_child(dist=distBottom,child=node)
			node.up=newInternalNode
			newInternalNode.add_child(dist=bestLen, name=sample)
			newInternalNode.children[1].minorSequences=[]
			newInternalNode.children[1].up=newInternalNode
			newInternalNode.dist=distTop
			newInternalNode.children[1].probVect=newPartials
			newVect=mergeVectorsUpDown(vectUp,newInternalNode.dist,newPartials,bestLen,mutMatrix)
			newVect=shorten(newVect)
			newInternalNode.probVectUpRight=newVect
			newInternalNode.probVectUpLeft=childBestVect
			newVect=mergeVectors(node.probVect,node.dist,newPartials,bestLen,mutMatrix)
			newVect=shorten(newVect)
			newInternalNode.probVect=newVect
			newVect=mergeVectorsUpDown(childBestVect,0.0,newPartials,bestLen,mutMatrix)
			newVect=shorten(newVect)
			newInternalNode.probVectTot=newVect
			newVect=mergeVectorsUpDown(childBestVect,bestLen,newPartials,0.0,mutMatrix)
			newVect=shorten(newVect)
			newInternalNode.children[1].probVectTot=newVect
			if verbose:
				print("new internal node added to tree")
				print(newInternalNode.probVect)
				print(newInternalNode.probVectUpRight)
				print(newInternalNode.probVectUpLeft)
				print(newInternalNode.probVectTot)
			updatePartialsFromTop(node,newInternalNode.probVectUpRight,mutMatrix)
			updatePartialsFromBottom(newInternalNode.up,newInternalNode.probVect,child,newInternalNode,mutMatrix)

	return None








distances=distancesFromRefPunishNs(data)
print("Distances from the reference calculated")
#extract root genome among those closest to the reference but not empty
root=distances.pop(0)
#print("initial root:")
#print(root)
#print(data[root[1]])

#default initial length for branches
bLen=1.0/lRef

#initialize tree to just the initial root sample
t1=Tree()
t1.name=root[1]
t1.probVect=probVectTerminalNode(data[root[1]],0.0)
t1.probVectTot=rootVector0(t1.probVect,rootFreqs)
t1.minorSequences=[]

numSamples=0
for d in distances:
	numSamples+=1
	sample=d[1]
	newPartials=probVectTerminalNode(data[sample],0.0)
	#print("\n")
	#print(sample)
	if (numSamples%500)==0:
		print("Sample num "+str(numSamples))
		#print(sample)
		#print(newPartials)
	node , bestNewLK, isMidNode=findBestParent(t1,newPartials,sample,bLen,float('-inf'),t1,0,False)
	if bestNewLK<0.5:
		newRoot=placeSampleOnTreeNew(node,newPartials,sample,bLen,bestNewLK,isMidNode)
		if newRoot!=None:
			t1=newRoot

print("Tree finalized")
treeString=(t1.write()).replace(")1:","):")
totMinors=0
for s in t1.get_leaves():
	if len(s.minorSequences)>0:
		totMinors+=len(s.minorSequences)
		newString="("+s.name+":0"
		for s2 in s.minorSequences:
			newString+=","+s2+":0"
		newString+="):"
		treeString=treeString.replace(s.name+":",newString)
#print("leaves in tree: ")
#print(len(t1.get_leaves()))
#print("tot minors:")
#print(totMinors)
		
#print(treeString)
file=open(outputFile,"w")
file.write(treeString)
file.close()
print("Missed minor samples: "+str(totalMissedMinors[0]))
exit()

#1)ALLOW SAMPLE PLACEMENT IN THE MIDDLE OF LONG BRANCHES.
#2)ESTIMATE SUBSTITUTION MODEL ON THE GO.


#estimate susbstitution rates?
#optimize branch lengths?
#Optimize root position?
#Calculate pairwise distances?
#try to improve topology?
















