import sys
#import math
from math import log
#import numpy as np
import argparse
from time import time
#from ete3 import Tree
import os.path

#©EMBL-European Bioinformatics Institute, 2021

#Estimate a tree using fastLK from a diff format and using iterative sample placement.

parser = argparse.ArgumentParser(description='Run estimation of mutation rates and synonymous site selection from SARS-CoV-2 data, and generate plots.')
#parser.add_argument('--path',default="", help='path where to find files and plot results.')
parser.add_argument('--input',default="/Users/demaio/Desktop/GISAID-hCoV-19-phylogeny-2021-03-12/phylogenetic_inference/2021-03-31_unmasked_differences_reduced.txt_consensus-based.txt", help='input diff file name.')
parser.add_argument('--reference',default="/Users/demaio/Desktop/GISAID-hCoV-19-phylogeny-2021-03-12/phylogenetic_inference/2021-03-31_unmasked_differences_reduced.txt_consensus.fa", help='input reference file name found in folder specified by --path.')
parser.add_argument('--output',default="/Users/demaio/Desktop/GISAID-hCoV-19-phylogeny-2021-03-12/phylogenetic_inference/iterativeFastLK", help='output newick file name to be written in the folder --path.')
parser.add_argument("--onlyNambiguities", help="Treat all ambiguities as N (total missing information).", action="store_true")
parser.add_argument("--useLogs", help="Calculate logarithms of non-mutation probabilities, otherwise approximate them.", action="store_true")
parser.add_argument("--thresholdProb",help="relative probability threshold used to ignore possible states with very low probabilities.",  type=float, default=0.0000001)
parser.add_argument("--thresholdLogLK",help="logLK difference threshold to consider a logLk close to optimal.",  type=float, default=200.0)
parser.add_argument("--allowedFails",help="Number of times one can go down the tree without inclreasing placement likelihood before the tree traversal is stopped (only applies to non-0 branch lengths).",  type=int, default=2)
parser.add_argument("--bLenAdjustment",help="If >1, try placing also with a longer bLen than the standard bLen.",  type=int, default=1)
parser.add_argument("--verbose", help="Print to screen a lot of stuff.", action="store_true")
parser.add_argument("--model", help="Which substitution model should be used. Allowed models so far are JC, GTR (default) or UNREST.", default="GTR")
parser.add_argument("--bLenFactor",help="split branch length by this factor when looking for best child of best node.",  type=float, default=2.0)
parser.add_argument("--overwrite", help="Overwrite previous results if already present.", action="store_true")
args = parser.parse_args()

onlyNambiguities=args.onlyNambiguities
useLogs=args.useLogs
thresholdProb=args.thresholdProb
verbose=args.verbose
#pathSimu=args.path
#inputFile=pathSimu+args.input
#outputFile=pathSimu+args.output
#refFile=pathSimu+args.reference
inputFile=args.input
outputFile=args.output
refFile=args.reference
allowedFails=args.allowedFails
model=args.model
bLenAdjustment=args.bLenAdjustment
thresholdLogLK=args.thresholdLogLK
bLenFactor=args.bLenFactor
overwrite=args.overwrite
example=False

if os.path.isfile(outputFile+"_tree.tree")  and (not overwrite):
	print("File "+outputFile+"_tree.tree already exists, quitting fastLK tree inference. Use option --overwrite if you want to overwirte previous inference.")
	exit()



class Tree(object):
	def __init__(self, name='', children=None, dist=1.0):
		self.name = name
		self.dist = dist
		self.children = []
		self.up=None
		if children is not None:
			for child in children:
				self.add_child(child)
	def __repr__(self):
		return self.name
	def add_child(self, node):
		assert isinstance(node, Tree)
		self.children.append(node)

# class Tree:
# 	dist = 1.0
# 	children=[]
# 	up=None
# 	name=""
	#minorSequences=[]
	# def add_child(self,name="",child=None,dist=1.0):
	# 	if child==None:
	# 		newChild=Tree()
	# 		newChild.name=name
	# 		newChild.dist=dist
	# 		#newChild.up=self
	# 		(self.children).append(newChild)
	# 	else:
	# 		child.dist=dist
	# 		#child.up=self
	# 		#newChild=child
	# 		(self.children).append(child)
	# def __init__(self, name="", up=None, dist=1.0):
	# 	self.name = name
	# 	self.up = up
	# 	self.dist=dist

# root=Tree()
# root.name="root"
# newNode=Tree()
# newNode.name="sample"
# newNode.dist=0.1
# newNode.up=root
# root.add_child(newNode)
# #(root.children).append(newNode)
# print(root.name)
# print(root.children)
# print(root.children[0].name)
# print(root.children[0].children)
# print((root.children[0]).name)
# print((root.children[0]).children)
# print(newNode.name)
# print(newNode.children)
# print(newNode.children[0].name)
# print(newNode.children[0].children)
# print((newNode.children[0]).name)
# print((newNode.children[0]).children)
# exit()



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
	#print("Ref genome length: "+str(lRef))
	file.close()
	#vector to count how many bases of each type are cumulatively in the reference genome up to a certain position
	#cumulativeBases=np.zeros((lRef+1,4))
	cumulativeBases=[[0,0,0,0]]
	for i in range(lRef):
		cumulativeBases.append(list(cumulativeBases[i]))
		#for k in range(4):
		#	cumulativeBases[i+1][k]=cumulativeBases[i][k]
		cumulativeBases[i+1][allelesLow[ref[i]]]+=1
	#print("cumulativeBases")
	#print(cumulativeBases[-1])
	rootFreqs=[0.0,0.0,0.0,0.0]
	#rootFreqs=np.zeros(4)
	rootFreqsLog=[0.0,0.0,0.0,0.0]
	#rootFreqsLog=np.zeros(4)
	for i in range(4):
		rootFreqs[i]=cumulativeBases[-1][i]/float(lRef)
		rootFreqsLog[i]=log(rootFreqs[i])
	#print("ref base frequencies and log")
	#print(rootFreqs)
	#print(rootFreqsLog)
	return ref, cumulativeBases, rootFreqs, rootFreqsLog

ref, cumulativeBases, rootFreqs, rootFreqsLog = collectReference(refFile)
lRef=len(ref)
if model=="JC":
	rootFreqs=[0.25,0.25,0.25,0.25]
	rootFreqsLog=[log(0.25),log(0.25),log(0.25),log(0.25)]

# if example:
# 	ref="aaccggtt"
# 	lRef=len(ref)
# 	rootFreqs=[0.25,0.25,0.25,0.25]
# 	rootFreqsLog=[log(0.25),log(0.25),log(0.25),log(0.25)]
# 	cumulativeBases=np.zeros((lRef+1,4))
# 	for i in range(lRef):
# 		for k in range(4):
# 			cumulativeBases[i+1][k]=cumulativeBases[i][k]
# 		cumulativeBases[i+1][allelesLow[ref[i]]]+=1

#TO DO
#recreate a simple scenario of a recombinant sample and with just a few samples and small genome, and see if there is still a bias in placing at the root, and why that is.



def readConciseAlignment(fileName):
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
	print(str(nSeqs)+" sequences in diff file.")
	return data

#read sequence data from file
data=readConciseAlignment(inputFile)

# if example:
# 	data={}
# 	data["sampleCG"]=[]
# 	data["sampleTG"]=[("t",3)]
# 	data["sampleCT"]=[("t",5)]
# 	data["sampleTT"]=[("t",3),("t",5)]

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
	newVec=[a/maxP for a in vec]
	#newVec=vec/maxP
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
			elif entryNew[4]!=entryOld[4]:
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

def updateSubMatrix(pseudoMutCounts,model):
	#mutMatrix=np.zeros((4,4))
	mutMatrix=[[0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0]]
	if model=="UNREST":
		for i in range4:
			for j in range4:
				mutMatrix[i][j]=pseudoMutCounts[i][j]/rootFreqs[i]
	elif model=="GTR":
		for i in range4:
			for j in range4:
				mutMatrix[i][j]=(pseudoMutCounts[i][j]+pseudoMutCounts[j][i])/rootFreqs[i]
	else:
		print("Substitution model not recognised! Exiting")
		exit()
	mutMatrix[0][0]=-(mutMatrix[0][1]+mutMatrix[0][2]+mutMatrix[0][3])
	mutMatrix[1][1]=-(mutMatrix[1][0]+mutMatrix[1][2]+mutMatrix[1][3])
	mutMatrix[2][2]=-(mutMatrix[2][0]+mutMatrix[2][1]+mutMatrix[2][3])
	mutMatrix[3][3]=-(mutMatrix[3][0]+mutMatrix[3][1]+mutMatrix[3][2])
	totRate=-(rootFreqs[0]*mutMatrix[0][0]+ rootFreqs[1]*mutMatrix[1][1]+ rootFreqs[2]*mutMatrix[2][2]+ rootFreqs[3]*mutMatrix[3][3] )
	for i in range4:
		for j in range4:
			mutMatrix[i][j]=mutMatrix[i][j]/totRate
	return mutMatrix

#preliminary nuc mutation rate matrix, from De Maio et al 2021
mutMatrix=[[0.0,0.039,0.310,0.123],[0.140,0.0,0.022,3.028],[0.747,0.113,0.0,2.953],[0.056,0.261,0.036,0.0]]
if model=="JC":
	mutMatrix=[[0.0,0.333,0.333,0.333],[0.333,0.0,0.333,0.333],[0.333,0.333,0.0,0.333],[0.333,0.333,0.333,0.0]]
else:
	pseudoMutCounts=[[0.0,1.0,5.0,2.0],[2.0,0.0,1.0,40.0],[5.0,2.0,0.0,20.0],[2.0,3.0,1.0,0.0]]
	mutMatrix=updateSubMatrix(pseudoMutCounts,model)

cumulativeRate=[0.0]
for i in range(lRef):
	ind=allelesLow[ref[i]]
	cumulativeRate.append(cumulativeRate[-1]+mutMatrix[ind][ind])



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
						if entry2[0]=="R":
							i2=allelesLow[ref[pos-1]]
						else:
							i2=alleles[entry2[0]]
						if entry1[4][i2]>0.5:
							found2bigger=True
						else:
							return 0
					elif entry2[0]=="O":
						if entry1[0]=="R":
							i1=allelesLow[ref[pos-1]]
						else:
							i1=alleles[entry1[0]]
						if entry2[4][i1]>0.5:
							found1bigger=True
						else:
							return 0
					else:
						return 0
				elif entry1[0]=="O":
					for j in range4:
						if entry2[4][j]>entry1[4][j]+0.01:
							found1bigger=True
						elif entry1[4][j]>entry2[4][j]+0.01:
							found2bigger=True
				if found1bigger and found2bigger:
					return 0
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
def appendProb(probVectP,probVectC,bLen,mutMatrix):
	Lkcost, indexEntry1, indexEntry2, totalFactor, pos = 0.0, 0, 0, 1.0, 1
	#indexEntry1=0
	#indexEntry2=0
	#totalFactor=1.0
	#pos=1
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
	while 1:
		if entry2[0]=="N":
			pass
		elif entry1[0]=="N":
			pass
		elif entry1[0]=="R":
			if entry2[0]=="R":
				#cB1=cumulativeBases[end]
				#cB2=cumulativeBases[pos-1]
				if entry1[4]:
					#Lkcost+=(bLen+entry1[3]+entry1[5])*(mutMatrix[0][0]*(cumulativeBases[end][0]-cumulativeBases[posm1][0])+mutMatrix[1][1]*(cumulativeBases[end][1]-cumulativeBases[posm1][1])+mutMatrix[2][2]*(cumulativeBases[end][2]-cumulativeBases[posm1][2])+mutMatrix[3][3]*(cumulativeBases[end][3]-cumulativeBases[posm1][3]))
					#Lkcost+=(bLen+entry1[3]+entry1[5])*(mutMatrix[0][0]*(cB1[0]-cB2[0])+mutMatrix[1][1]*(cB1[1]-cB2[1])+mutMatrix[2][2]*(cB1[2]-cB2[2])+mutMatrix[3][3]*(cB1[3]-cB2[3]))
					Lkcost+=(bLen+entry1[3]+entry1[5])*(cumulativeRate[end]-cumulativeRate[pos-1])
					#tot=0.0
					#for i in range4:
					#	tot+=mutMatrix[i][i]*(cumulativeBases[end][i]-cumulativeBases[posm1][i])
					#tot*=(bLen+entry1[3]+entry1[5])
					#Lkcost+=tot
				else:
					#Lkcost+=(bLen+entry1[3])*(mutMatrix[0][0]*(cumulativeBases[end][0]-cumulativeBases[posm1][0])+mutMatrix[1][1]*(cumulativeBases[end][1]-cumulativeBases[posm1][1])+mutMatrix[2][2]*(cumulativeBases[end][2]-cumulativeBases[posm1][2])+mutMatrix[3][3]*(cumulativeBases[end][3]-cumulativeBases[posm1][3]))
					#Lkcost+=(bLen+entry1[3])*(mutMatrix[0][0]*(cB1[0]-cB2[0])+mutMatrix[1][1]*(cB1[1]-cB2[1])+mutMatrix[2][2]*(cB1[2]-cB2[2])+mutMatrix[3][3]*(cB1[3]-cB2[3]))
					Lkcost+=(bLen+entry1[3])*(cumulativeRate[end]-cumulativeRate[pos-1])
					#tot=0.0
					#for i in range4:
					#	tot+=mutMatrix[i][i]*(cumulativeBases[end][i]-cumulativeBases[posm1][i])
					#tot*=(bLen+entry1[3])
					#Lkcost+=tot
			elif entry2[0]=="O":
				i1=allelesLow[ref[pos-1]]
				if entry1[4]:
					if entry2[4][i1]>0.5:
						Lkcost+=mutMatrix[i1][i1]*(bLen+entry1[3]+entry1[5])
					else:
						bLen15=bLen+entry1[5]
						tot=0.0
						for i in range4:
							if i1==i:
								tot2=rootFreqs[i]*(1.0+mutMatrix[i][i1]*entry1[3])
							else:
								tot2=rootFreqs[i]*mutMatrix[i][i1]*entry1[3]

							tot+=tot2*(entry2[4][i] + bLen15*(mutMatrix[i][0]*entry2[4][0]+mutMatrix[i][1]*entry2[4][1]+mutMatrix[i][2]*entry2[4][2]+mutMatrix[i][3]*entry2[4][3]))
							#tot3=0.0
							#	if i!=j:
							#		tot3+=mutMatrix[i][j]*bLen15*entry2[4][j]
							#	else:
							#		tot3+=(1.0+mutMatrix[i][j]*bLen15)*entry2[4][j]
							#tot2=tot2*tot3
							#tot+=tot2
						#if tot>sys.float_info.min:
						totalFactor*=(tot/rootFreqs[i1])
							#Lkcost+=log(tot/rootFreqs[i1])
						#else:
						#	return float("-inf")
				else:
					if entry2[4][i1]>0.5:
						Lkcost+=mutMatrix[i1][i1]*(bLen+entry1[3])
					else:
						#totalFactor*=(entry2[4][i1]+(bLen+entry1[3])*(mutMatrix[i1][0]*entry2[4][0]+mutMatrix[i1][1]*entry2[4][1]+mutMatrix[i1][2]*entry2[4][2]+mutMatrix[i1][3]*entry2[4][3]))
						tot=0.0
						for j in range4:
							if entry2[4][j]>0.5:
								tot+=mutMatrix[i1][j]
						tot*=(bLen+entry1[3])
						tot+=entry2[4][i1]
						#bLen13=bLen+entry1[3]
						#	if i1!=j:
						#		tot+=mutMatrix[i1][j]*bLen13*entry2[4][j]
						#	else:
						#		tot+=(1.0+mutMatrix[i1][j]*bLen13)*entry2[4][j]
						totalFactor*=tot
			else:
				if entry1[4]:
					i2=alleles[entry2[0]]
					i1=allelesLow[ref[pos-1]]
					#tot=rootFreqs[i1]*mutMatrix[i1][i2]*(bLen+entry1[5]) + rootFreqs[i2]*mutMatrix[i2][i1]*entry1[3] + entry1[3]*(bLen+entry1[5])*(rootFreqs[0]*mutMatrix[0][i1]*mutMatrix[0][i2]+rootFreqs[1]*mutMatrix[1][i1]*mutMatrix[1][i2]+rootFreqs[2]*mutMatrix[2][i1]*mutMatrix[2][i2]+rootFreqs[3]*mutMatrix[3][i1]*mutMatrix[3][i2]) 
					tot=0.0
					bLen15=bLen+entry1[5]
					for i in range4:
						if i1==i:
							tot+=rootFreqs[i]*(1.0+mutMatrix[i][i1]*entry1[3])*mutMatrix[i][i2]*bLen15
						elif i2==i:
							tot+=rootFreqs[i]*mutMatrix[i][i1]*entry1[3]*(1.0+mutMatrix[i][i2]*bLen15)
						else:
							tot+=rootFreqs[i]*mutMatrix[i][i1]*entry1[3]*mutMatrix[i][i2]*bLen15
					totalFactor*=(tot/rootFreqs[i1])
					#Lkcost+=log(tot/rootFreqs[i1])
				else:
					#if (bLen+entry1[3])<thresholdProb2:
					#	return float("-inf")
					#i2=alleles[entry2[0]]
					#i1=allelesLow[ref[pos-1]]
					totalFactor*=(mutMatrix[allelesLow[ref[pos-1]]][alleles[entry2[0]]]*(bLen+entry1[3]))
					#Lkcost+=log(mutMatrix[allelesLow[ref[pos-1]]][alleles[entry2[0]]]*bLen13)
		elif entry1[0]=="O":
			bLen13=bLen+entry1[3]
			#tot1=np.sum(entry1[4])
			#LkcostOriginal+=log(tot1)
			if entry2[0]=="O":
				tot=0.0
				for j in range4:
					if entry2[4][j]>0.5:
						tot+=entry1[4][j] + bLen13*(mutMatrix[0][j]*entry1[4][0]+mutMatrix[1][j]*entry1[4][1]+mutMatrix[2][j]*entry1[4][2]+mutMatrix[3][j]*entry1[4][3])
						#for i in range4:
						#	if i==j:
						#		tot+=(1.0+mutMatrix[i][i]*bLen13)*entry1[4][i]
						#	else:
						#		tot+=mutMatrix[i][j]*bLen13*entry1[4][i]
				#totalFactor*=(tot/np.sum(entry1[4]))
				totalFactor*=(tot/sum(entry1[4]))
				#if tot>sys.float_info.min:
				#	Lkcost+=log(tot/tot1)
				#else:
				#	return float("-inf")
			else:
				if entry2[0]=="R":
					i2=allelesLow[ref[pos-1]]
				else:
					i2=alleles[entry2[0]]
				#tot=entry1[4][i2]+bLen13*(entry1[4][0]*mutMatrix[0][i2]+entry1[4][1]*mutMatrix[1][i2]+entry1[4][2]*mutMatrix[2][i2]+entry1[4][3]*mutMatrix[3][i2])
				totalFactor*=((entry1[4][i2]+bLen13*(entry1[4][0]*mutMatrix[0][i2]+entry1[4][1]*mutMatrix[1][i2]+entry1[4][2]*mutMatrix[2][i2]+entry1[4][3]*mutMatrix[3][i2]))/sum(entry1[4]))
				#tot=0.0
				#for i in range4:
				#	if i==i2:
				#		tot+=entry1[4][i]*(1.0+mutMatrix[i][i]*bLen13)
				#	else:
				#		tot+=entry1[4][i]*mutMatrix[i][i2]*bLen13
				#totalFactor*=(tot/np.sum(entry1[4]))
				#totalFactor*=(tot/sum(entry1[4]))
				#if tot>sys.float_info.min:
				#	Lkcost+=log(tot/tot1)
				#else:
				#	return float("-inf")
		#entry1 is a non-ref nuc
		else:
			i1=alleles[entry1[0]]
			if entry2[0]==entry1[0]:
				if entry1[4]:
					Lkcost+=mutMatrix[i1][i1]*(bLen+entry1[3]+entry1[5])
				else:
					Lkcost+=mutMatrix[i1][i1]*(bLen+entry1[3])
			else:
				if entry2[0]=="O":
					if entry1[4]:
						if entry2[4][i1]>0.5:
							Lkcost+=mutMatrix[i1][i1]*(bLen+entry1[3]+entry1[5])
						else:
							bLen15=bLen+entry1[5]
							tot=0.0
							for i in range4:
								if i1==i:
									tot2=rootFreqs[i]*(1.0+mutMatrix[i][i1]*entry1[3])
								else:
									tot2=rootFreqs[i]*mutMatrix[i][i1]*entry1[3]
								tot+=tot2*(entry2[4][i] + bLen15*(mutMatrix[i][0]*entry2[4][0]+mutMatrix[i][1]*entry2[4][1]+mutMatrix[i][2]*entry2[4][2]+mutMatrix[i][3]*entry2[4][3]))
								#tot3=0.0
								#for j in range4:
								#	if entry2[4][j]>0.5:
								#		if i!=j:
								#			tot3+=mutMatrix[i][j]*bLen15
								#		else:
								#			tot3+=(1.0+mutMatrix[i][j]*bLen15)
								#tot+=tot2*tot3
							totalFactor*=(tot/rootFreqs[i1])
							#if tot>sys.float_info.min:
							#	Lkcost+=log(tot/rootFreqs[i1])
							#else:
							#	return float("-inf")
					else:
						if entry2[4][i1]>0.5:
							Lkcost+=mutMatrix[i1][i1]*(bLen+entry1[3])
						else:
							#totalFactor*=(entry2[4][i1]+(bLen+entry1[3])*(mutMatrix[i1][0]*entry2[4][0]+mutMatrix[i1][1]*entry2[4][1]+mutMatrix[i1][2]*entry2[4][2]+mutMatrix[i1][3]*entry2[4][3]))
							tot=0.0
							#bLen13=bLen+entry1[3]
							for j in range4:
								if entry2[4][j]>0.5:
							#		if i1==j:
							#			tot+=(1.0+mutMatrix[i1][i1]*bLen13)
							#		else:
										tot+=mutMatrix[i1][j]
							totalFactor*=(tot*(bLen+entry1[3]))

				else:
					if entry2[0]=="R":
						i2=allelesLow[ref[pos-1]]
					else:
						i2=alleles[entry2[0]]
					if entry1[4]:
						#tot=rootFreqs[i1]*mutMatrix[i1][i2]*(bLen+entry1[5]) + rootFreqs[i2]*mutMatrix[i2][i1]*entry1[3] + entry1[3]*(bLen+entry1[5])*(rootFreqs[0]*mutMatrix[0][i1]*mutMatrix[0][i2]+rootFreqs[1]*mutMatrix[1][i1]*mutMatrix[1][i2]+rootFreqs[2]*mutMatrix[2][i1]*mutMatrix[2][i2]+rootFreqs[3]*mutMatrix[3][i1]*mutMatrix[3][i2]) 
						totalFactor*=((rootFreqs[i1]*mutMatrix[i1][i2]*(bLen+entry1[5]) + rootFreqs[i2]*mutMatrix[i2][i1]*entry1[3] + entry1[3]*(bLen+entry1[5])*(rootFreqs[0]*mutMatrix[0][i1]*mutMatrix[0][i2]+rootFreqs[1]*mutMatrix[1][i1]*mutMatrix[1][i2]+rootFreqs[2]*mutMatrix[2][i1]*mutMatrix[2][i2]+rootFreqs[3]*mutMatrix[3][i1]*mutMatrix[3][i2]) )/rootFreqs[i1])
						#tot=0.0
						#bLen15=bLen+entry1[5]
						#for i in range4:
						#	if i1==i:
						#		tot+=rootFreqs[i]*(1.0+mutMatrix[i][i1]*entry1[3])*mutMatrix[i][i2]*bLen15
						#	elif i2==i:
						#		tot+=rootFreqs[i]*mutMatrix[i][i1]*entry1[3]*(1.0+mutMatrix[i][i2]*bLen15)
						#	else:
						#		tot+=rootFreqs[i]*mutMatrix[i][i1]*entry1[3]*mutMatrix[i][i2]*bLen15
						#totalFactor*=(tot/rootFreqs[i1])
					else:
						#if bLen13<thresholdProb2:
						#	return float("-inf")
						#if entry2[0]=="R":
						#	i2=allelesLow[ref[pos-1]]
						#else:
						#	i2=alleles[entry2[0]]
						totalFactor*=(mutMatrix[i1][i2]*(bLen+entry1[3]))
						#Lkcost+=log(mutMatrix[i1][i2]*bLen13)

		pos+=length
		if pos>lRef:
			break
		if pos>end1:
			indexEntry1+=1
			entry1=probVectP[indexEntry1]
			pos1=entry1[1]
			#if entry1[0]!="N" and entry1[0]!="R":
			#	end1=pos1
			#else:
			end1=pos1+entry1[2]-1
		if pos>end2:
			indexEntry2+=1
			entry2=probVectC[indexEntry2]
			pos2=entry2[1]
			#if entry2[0]!="N" and entry2[0]!="R":
			#	end2=pos2
			#else:
			end2=pos2+entry2[2]-1
		end=min(end1,end2)
		#if end2<end:
		#	end=end2
		length=end+1-pos
	if totalFactor>sys.float_info.min:
		return Lkcost+log(totalFactor)
	else:
		return float("-inf")
	#return Lkcost+log(totalFactor)
	#-LkcostOriginal



#number of samples that could have been placed as major of another sample, but weren't due to sample placement order
totalMissedMinors=[0]
#function to find the best node in the tree where to append the new sample, based on the node's tot posterior state probabilities.
def findBestParent(t1,diffs,sample,bLen,mutMatrix):
	# print("findBestParent(), Node")
	# print(t1)
	# print("With probVectTot: ")
	# print(t1.probVectTot)
	bestNodeSoFar=t1
	bestLKdiff=float("-inf")
	bestIsMidNode=False
	bestUpLK=float("-inf")
	bestDownLK=float("-inf")
	bestDownNode=None
	nodesToVisit=[(t1,float("-inf"),0)]
	adjustBLen=False

	#while len(nodesToVisit)>0:
	while nodesToVisit:
		t1,parentLK,failedPasses=nodesToVisit.pop()
		#element=nodesToVisit.pop()
		#t1=element[0]
		#parentLK=element[1]
		#failedPasses=element[2]

		#if t1.is_leaf():
		if not t1.children:
			comparison=isMinorSequence(t1.probVect,diffs)
			if comparison==1:
				t1.minorSequences.append(sample)
				return t1, 1.0, False, float("-inf"), float("-inf"), None, False
			elif comparison==2:
				#if verbose:
				#	print("Second sequence is more informative than first, but for now the two sequences are not merged.")
				#	print(t1.name)
				#	print(sample)
				totalMissedMinors[0]+=1
	
		adjusted=False
		if t1.dist>thresholdProb2 and t1.up!=None:
			LKdiff2=appendProb(t1.probVectTotUp,diffs,bLen,mutMatrix)
			#if LKdiff2<(bestLKdiff-thresholdLogLK):
			#	return None , float("-inf") , bestIsMidNode, bestUpLK, bestDownLK, bestDownNode, LKdiff2, t1, adjustBLen
			#else:
			if bLenAdjustment>1:
				newLKdiff2=appendProb(t1.probVectTotUp,diffs,bLen*bLenAdjustment,mutMatrix)
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

		if t1.dist>thresholdProb2:
			LKdiff=appendProb(t1.probVectTot,diffs,bLen,mutMatrix)
			#if LKdiff<(bestLKdiff-thresholdLogLK):
			#	return None , float("-inf") , bestIsMidNode, bestUpLK, bestDownLK, bestDownNode, LKdiff2, t1, adjustBLen
			if bLenAdjustment>1:
				newLKdiff=appendProb(t1.probVectTot,diffs,bLen*bLenAdjustment,mutMatrix)
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
			elif LKdiff<(parentLK-1.0):
				failedPasses+=1
				#if failedPasses>allowedFails and LKdiff<(bestLKdiff-thresholdLogLK):
					#return bestNodeSoFar , bestLKdiff , bestIsMidNode, bestUpLK, bestDownLK, bestDownNode, LKdiff2, t1, adjustBLen
				#	return None , float("-inf") , bestIsMidNode, bestUpLK, bestDownLK, bestDownNode, LKdiff2, t1, adjustBLen
		else:
			#LKdiff=float("-inf")
			LKdiff=parentLK
		# print("t1.dist "+str(t1.dist))
		# print("LKdiff "+str(LKdiff))
		# print("LKdiff2 "+str(LKdiff2))

		if failedPasses<=allowedFails or LKdiff>(bestLKdiff-thresholdLogLK):
			for c in t1.children:
				nodesToVisit.append((c,LKdiff,failedPasses))

	bestOfChildLK=float("-inf")
	bestChild=None
	if bestIsMidNode:
		return bestNodeSoFar , bestLKdiff , bestIsMidNode, bestUpLK, bestDownLK, bestDownNode, adjustBLen
	else:
		nodesToVisit=[]
		for c in bestNodeSoFar.children:
			nodesToVisit.append(c)
		while nodesToVisit:
			t1=nodesToVisit.pop()
			# print("Before child loop, best child")
			# print(bestChild)
			# print(bestOfChildLK)
			if t1.dist<=thresholdProb2:
				for c in t1.children:
					nodesToVisit.append(c)
			else:
				newBLen2=t1.dist/2
				bestLKdiff2=float("-inf")
				furtherNode=-1
				newProbVect2=t1.probVectTotUp
				while True:
					newLKdiff2=appendProb(newProbVect2,diffs,bLen,mutMatrix)
					if bLenAdjustment>1:
						newLKdiff3=appendProb(newProbVect2,diffs,bLen*bLenAdjustment,mutMatrix)
						if newLKdiff3>newLKdiff2:
							newLKdiff2=newLKdiff3
					if newLKdiff2>bestLKdiff2:
						bestLKdiff2=newLKdiff2
					else:
						break
					newBLen2=newBLen2/2
					if newBLen2<=bLen/(bLenFactor+thresholdProb):
						break
					furtherNode+=1
					newProbVect2=t1.furtherMidNodes[furtherNode]
				
				if bestLKdiff2>bestOfChildLK:
					bestOfChildLK=bestLKdiff2
					bestChild=t1
		return bestNodeSoFar , bestLKdiff , bestIsMidNode, bestUpLK, bestOfChildLK, bestChild, adjustBLen








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
			if newEntry[3]>thresholdProb4:
				for i in range4:
					tot=0.0
					for j in range4:
						if i==j:
							tot+=rootFreqs[i]*(1.0+mutMatrix[i][i]*(entry[3]))*entry[4][j]
						else:
							tot+=rootFreqs[i]*mutMatrix[i][j]*(entry[3])*entry[4][j]
					newEntry[4].append(tot)
			else:
				for i in range4:
					newEntry[4].append(entry[4][i]*rootFreqs[i])
		else:
			newEntry[4]=True
			#newEntry[4]=False
			#distance from the root on the downward branch
			newEntry.append(0.0)
			#if newEntry[2]>1000 and newEntry[5]>thresholdProb2:
			#	print("Problem with rootVector0")
			#	exit()
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
					if entry2[0]!="O" and entry2[0]!="N":
						probVect[-1][4]=True
						probVect[-1].append(0.0)
						probVect[-1][3]+=bLenDown
					elif entry2[0]=="O":
						probVect[-1][4]=[]
						for i in range4:
							tot=0.0
							for j in range4:
								if i==j:
									tot+=rootFreqs[i]*(1.0+mutMatrix[i][i]*(entry2[3]+bLenDown))*entry2[4][j]
								else:
									tot+=rootFreqs[i]*mutMatrix[i][j]*(entry2[3]+bLenDown)*entry2[4][j]
							probVect[-1][4].append(tot)
						probVect[-1][3]=0.0
				elif entry2[0]=="N":
					probVect.append(list(entry1))
					if entry1[0]=="R": 
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
					#print("totLens")
					#print(totLen1)
					#print(totLen2)
					if entry1[0]!="O" and entry1[4]:
						totLen1+=entry1[5]
					#if entry2[0]!="O" and entry2[4]:
					#	totLen2+=entry2[5]
					if totLen1+totLen2<thresholdProb2 and entry1[0]!=entry2[0] and entry1[0]!="O" and entry2[0]!="O":
						# print("Warning, inside mergeVectorsUpDown, branch lengths are 0, but the two vectors are different")
						# print(probVectUp)
						# print(bLenUp)
						# print(probVectDown)
						# print(bLenDown)
						# print(entry1)
						# print(entry2)
						return None
					elif entry2[0]=="R" and (totLen2<thresholdProb2):
						probVect.append(["R",pos,length,0.0,False])
					elif (totLen2<thresholdProb2) and entry2[0]!="O":
						probVect.append([entry2[0],pos,1,0.0,False])
					elif entry1[0]=="R":
						if entry2[0]=="R" or (totLen1<thresholdProb2):
							probVect.append(["R",pos,length,0.0,False])
						else:
							#print("Inside big loop")
							i1=allelesLow[ref[pos-1]]
							#newVec=np.ones(4)
							newVec=[1.0,1.0,1.0,1.0]
							if entry1[4]:
								rootVec=[1.0,1.0,1.0,1.0]
								#rootVec=np.ones(4)
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
						newVec=[1.0,1.0,1.0,1.0]
						#newVec=np.ones(4)
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
							#newVec=np.ones(4)
							newVec=[1.0,1.0,1.0,1.0]
							if entry1[4]:
								#rootVec=np.ones(4)
								rootVec=[1.0,1.0,1.0,1.0]
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

				# if probVect[-1][0]=="R" and probVect[-1][2]>20000 and probVect[-1][4]==True:
				# 	print("Strange vector with True inside mergeVectorsUpDown")
				# 	print(length)
				# 	print(pos)
				# 	print(bLenUp)
				# 	print(bLenDown)
				# 	print(entry1)
				# 	print(entry2)
				# 	print(probVect[-1])
				# 	exit()

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




#merge two child partial likelihood vectors at the root and calculate the logLk, and return vector and logLK.
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
								cumulPartLk+=log(1.0+mutMatrix[i][i]*totLen)*(cumulativeBases[end][i]-cumulativeBases[pos-1][i])
						else:
							for i in range4:
								cumulPartLk+=mutMatrix[i][i]*totLen*(cumulativeBases[end][i]-cumulativeBases[pos-1][i])
					elif entry2[0]=="O":
						i1=allelesLow[ref[pos-1]]
						newVec=[1.0,1.0,1.0,1.0]
						#newVec=np.ones(4)
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
						cumulPartLk+=log(maxP)
					else:
						i2=alleles[entry2[0]]
						i1=allelesLow[ref[pos-1]]
						#newVec=np.ones(4)
						newVec=[1.0,1.0,1.0,1.0]
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
					newVec=[1.0,1.0,1.0,1.0]
					#newVec=np.ones(4)
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
					if state=="O":
						probVect.append([state,pos,1,0.0,newVec])
					else:
						probVect.append([state,pos,1,0.0,False])
					cumulPartLk+=log(maxP)
				#entry1 is a non-ref nuc
				else:
					if entry2[0]==entry1[0]:
						i2=alleles[entry2[0]]
						probVect.append([entry1[0],pos,1,0.0,False])
						#Now update partial likelihood
						if useLogs:
							cumulPartLk+=log(1.0+mutMatrix[i2][i2]*totLen)
						else:
							cumulPartLk+=mutMatrix[i2][i2]*totLen
					else:
						i1=alleles[entry1[0]]
						#newVec=np.ones(4)
						newVec=[1.0,1.0,1.0,1.0]
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
							cumulPartLk+=log(maxP)
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
				elif bLen1+bLen2+entry1[3]+entry2[3]<thresholdProb2 and entry1[0]!=entry2[0] and entry1[0]!="O" and entry2[0]!="O" :
					#and (not (entry1[4] and entry1[5]>0))
					print("Warning, inside mergeVectors, branch lengths are 0, but the two vectors are different")
					print(probVect1)
					print(bLen1)
					print(probVect2)
					print(bLen2)
					print(entry1)
					print(entry2)
					#exit()
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
						#newVec=np.ones(4)
						newVec=[1.0,1.0,1.0,1.0]
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
						#newVec=np.ones(4)
						newVec=[1.0,1.0,1.0,1.0]
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
					#newVec=np.ones(4)
					newVec=[1.0,1.0,1.0,1.0]
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
						#newVec=np.ones(4)
						newVec=[1.0,1.0,1.0,1.0]
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
				logLK+=log(tot)
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
					vSum=sum(entry2[4])
					#vect=entry2[4]/vSum
					vect=[a/vSum for a in entry2[4]]
					if not vSum>thresholdProb4:
						print("Problem? sum of partial likelihoods is very low or nan")
						print(probVect1)
						print(probVect2)
						exit()
					if entry1[0]=="R":
						i1=allelesLow[ref[pos-1]]
					else:
						i1=alleles[entry1[0]]
					if vect[i1]<(1.0-thresholdProb):	
					#if vect[i1]<0.9999:
						return True
			else:
				if entry2[0]=="N":
					return True
				vSum=sum(entry1[4])
				if not vSum>thresholdProb4:
					print("Problem? sum of partial likelihoods is very low or nan")
					print(probVect1)
					print(probVect2)
					exit()
				#vect=entry1[4]/vSum
				vect=[a/vSum for a in entry1[4]]
				if entry2[0]=="R":
					i2=allelesLow[ref[pos-1]]
				else:
					i2=alleles[entry2[0]]
				#if vect[i2]<0.9999:
				if vect[i2]<(1.0-thresholdProb):	
					return True

		elif entry1[0]=="O":
			su1=sum(entry1[4])
			vect1=[a/su1 for a in entry1[4]]
			su2=sum(entry2[4])
			vect2=[a/su2 for a in entry2[4]]
			#vect1=entry1[4]/sum(entry1[4])
			#vect2=entry2[4]/sum(entry2[4])
			# maxP1=-1.0
			# maxP2=-1.0
			# for i in range4:
			# 	if vect1[i]>maxP1:
			# 		maxP1=vect1[i]
			# 	if vect2[i]>maxP2:
			# 		maxP2=vect2[i]
			# for i in range4:
			# 	if abs((vect1[i]/maxP1) - (vect2[i]/maxP2))>0.001:
			# 		return True
			for i in range4:
				diffVal=abs(vect1[i] - vect2[i])
				if diffVal>0.001 or (diffVal>thresholdProb and ((vect1[i]>thresholdProb2 and (diffVal/vect1[i])>1.5) or (vect2[i]>thresholdProb2 and (diffVal/vect2[i])>1.5))):
					return True
		elif entry1[0]!="N":
			if entry1[4]!=entry2[4]:
				return True
			diffVal=abs(entry1[3] - entry2[3])
			if diffVal>bLen/100:
				return True
			if entry1[4]:
				diffVal=abs(entry1[5] - entry2[5])
				if diffVal>bLen/100:
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



#create further mid-nodes for longer branches, so to make it faster to calculate appending probabilities
def createFurtherMidNodes(node,vectUp,bLen):
			node.furtherMidNodes=[]
			newBLen2=node.dist/4
			while newBLen2>bLen/(bLenFactor+thresholdProb):
				newProbVect2=mergeVectorsUpDown(vectUp,newBLen2,node.probVect,node.dist-newBLen2,mutMatrix)
				node.furtherMidNodes.append(newProbVect2)
				newBLen2=newBLen2/2


#if updating lk creates an inconsistency, this function can increase the bLen from 0 so that the inconsistency is resolved.
def updateBLen(node,childNum,bLen,mutMatrix):
	cNode=node.children[childNum]
	if childNum==0:
		vectUp=node.probVectUpRight
	else:
		vectUp=node.probVectUpLeft
	bestLK=appendProb(vectUp,cNode.probVect,bLen,mutMatrix)
	bestLen=bLen
	while bestLen>0.1*bLen:
		newBLen=bestLen/2
		newLK=appendProb(vectUp,cNode.probVect,newBLen,mutMatrix)
		if newLK>bestLK:
			bestLK=newLK
			bestLen=newBLen
		else:
			break
	if bestLen>0.7*bLen:
		while bestLen<10*bLen:
			newBLen=bestLen*2
			newLK=appendProb(vectUp,cNode.probVect,newBLen,mutMatrix)
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
	# print("Updating partials from top")
	# print(node.children)
	# print(node.name)
	# print(node.dist)
	# print("Name of children")
	# print(node.children[0].name)
	# print(node.children[1].name)
	if node.dist>thresholdProb2:
		newTot=mergeVectorsUpDown(probVectUp,node.dist/2,node.probVect,node.dist/2,mutMatrix)
		if newTot==None:
			if node.up.children[0]==node:
				updateBLen(node.up,0,bLen,mutMatrix)
			else:
				updateBLen(node.up,1,bLen,mutMatrix)
			return
		newTot=shorten(newTot)
		node.probVectTotUp=newTot
		#if node.dist>bLen/2:
		if node.dist>4*bLen/(bLenFactor+thresholdProb):
			createFurtherMidNodes(node,probVectUp,bLen)
	if len(node.children)==0 and node.dist>thresholdProb2:
		newTot=mergeVectorsUpDown(probVectUp,node.dist,node.probVect,0.0,mutMatrix)
		if newTot==None:
			if node.up.children[0]==node:
				updateBLen(node.up,0,bLen,mutMatrix)
			else:
				updateBLen(node.up,1,bLen,mutMatrix)
			return
		newTot=shorten(newTot)
		node.probVectTot=newTot
	elif len(node.children)>0:
		child0Vect=node.children[0].probVect
		child1Vect=node.children[1].probVect
		dist0=node.children[0].dist
		dist1=node.children[1].dist
		newUpRight=mergeVectorsUpDown(probVectUp,node.dist,child1Vect,dist1,mutMatrix)
		newUpLeft=mergeVectorsUpDown(probVectUp,node.dist,child0Vect,dist0,mutMatrix)
		if newUpRight==None or newUpLeft==None:
			if node.up.children[0]==node:
				updateBLen(node.up,0,bLen,mutMatrix)
			else:
				updateBLen(node.up,1,bLen,mutMatrix)
			return
		# print("Updating from top")
		# print(node.probVectUpRight)
		# print(newUpRight)
		# print(node.probVectUpLeft)
		# print(newUpLeft)
		#exit()
		updatedTot=False
		if areVectorsDifferent(node.probVectUpRight,newUpRight):
			# print("Vectors up right are different")
			# print(node.probVectUpRight)
			# print(newUpRight)
			# exit()
			newUpRight =shorten(newUpRight)
			node.probVectUpRight=newUpRight
			if node.dist>thresholdProb2:
				newTot=mergeVectorsUpDown(newUpRight,0.0,child0Vect,dist0,mutMatrix)
				if newTot==None:
					updateBLen(node,0,bLen,mutMatrix)
					return
				newTot =shorten(newTot)
				node.probVectTot=newTot
				updatedTot=True
			updatePartialsFromTop(node.children[0],newUpRight,mutMatrix)
		if areVectorsDifferent(node.probVectUpLeft,newUpLeft):
			# print("Vectors up left are different")
			# print(node.probVectUpLeft)
			# print(newUpLeft)
			# exit()
			newUpLeft =shorten(newUpLeft)
			node.probVectUpLeft=newUpLeft
			if (not updatedTot) and node.dist>thresholdProb2:
				newTot=mergeVectorsUpDown(newUpLeft,0.0,child1Vect,dist1,mutMatrix)
				if newTot==None:
					updateBLen(node,1,bLen,mutMatrix)
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
	if newVect==None:
		updateBLen(node,childNum,bLen,mutMatrix)
		return
	if verbose:
		print("Updating vectors from bottom, new vec")
		print(newVect)
		print("Old vec:")
		print(node.probVect)
	updatedTot=False
	if areVectorsDifferent(node.probVect,newVect):
		newVect =shorten(newVect)
		node.probVect=newVect
		if node.dist>thresholdProb2:
			if childNum==0:
				newTot=mergeVectorsUpDown(node.probVectUpRight,0.0,probVectDown,childDist,mutMatrix)
				if newTot==None:
					updateBLen(node,0,bLen,mutMatrix)
					return
			else:
				newTot=mergeVectorsUpDown(node.probVectUpLeft,0.0,probVectDown,childDist,mutMatrix)
				if newTot==None:
					updateBLen(node,1,bLen,mutMatrix)
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
			if newUpLeftVect==None:
				updateBLen(node,1,bLen,mutMatrix)
				return
			if areVectorsDifferent(node.probVectUpLeft,newUpLeftVect):
				if verbose:
					print("Moving up the update from left; new UpLeft and old UpLeft")
					print(newUpLeftVect)
					print(node.probVectUpLeft)
				newUpLeftVect =shorten(newUpLeftVect)
				node.probVectUpLeft=newUpLeftVect
				if (not updatedTot) and node.dist>thresholdProb2:
					if verbose:
						print("Updating tot: old tot and new tot")
						print(node.probVectTot)
					newTot=mergeVectorsUpDown(newUpLeftVect,0.0,otherChildVect,otherChildDist,mutMatrix)
					if newTot==None:
						updateBLen(node,1,bLen,mutMatrix)
						return
					newTot=shorten(newTot)
					node.probVectTot=newTot
					if verbose:
						print(newTot)
					updatedTot=True
				updatePartialsFromTop(node.children[1],newUpLeftVect,mutMatrix)
		else:
			newUpRightVect=mergeVectorsUpDown(vectUp,node.dist,probVectDown,childDist,mutMatrix)
			if newUpRightVect==None:
				updateBLen(node,1,bLen,mutMatrix)
				return
			if areVectorsDifferent(node.probVectUpRight,newUpRightVect):
				if verbose:
					print("Moving up the update from right")
				newUpRightVect =shorten(newUpRightVect)
				node.probVectUpRight=newUpRightVect
				if (not updatedTot) and node.dist>thresholdProb2:
					newTot=mergeVectorsUpDown(newUpRightVect,0.0,otherChildVect,otherChildDist,mutMatrix)
					if newTot==None:
						updateBLen(node,0,bLen,mutMatrix)
						return
					newTot=shorten(newTot)
					node.probVectTot=newTot
					updatedTot=True
				updatePartialsFromTop(node.children[0],newUpRightVect,mutMatrix)
		if updatedTot:
			newTot=mergeVectorsUpDown(vectUp,node.dist/2,node.probVect,node.dist/2,mutMatrix)
			newTot=shorten(newTot)
			node.probVectTotUp=newTot
			if node.dist>4*bLen/(bLenFactor+thresholdProb):
			#if node.dist>bLen/2:
				createFurtherMidNodes(node,vectUp,bLen)
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
					newTot=shorten(newTot)
					node.probVectTot=newTot
				updatePartialsFromTop(node.children[0],newUpRightVect,mutMatrix)
	if verbose:
		print("New upRight, UpLeft and tot:")
		print(node.probVectUpRight)
		print(node.probVectUpLeft)
		print(node.probVectTot)
	return


#function to add new mutation events in new sample to the pre-exisitng pseudocounts so to improve the estimate of the substitution rates
def updatePesudoCounts(oldVector,newPartials,pseudoMutCounts):
	if model!="JC":
			probVect1=oldVector
			probVect2=newPartials
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
				if entry1[0]=="R":
						if entry2[0]!="R" and entry2[0]!="N" and entry2[0]!="O":
							i1=allelesLow[ref[pos-1]]
							i2=alleles[entry2[0]]
							pseudoMutCounts[i1][i2]+=1
								
				elif entry1[0]!="O" and entry1[0]!="N":
					if entry2[0]=="R":
							i1=alleles[entry1[0]]
							i2=allelesLow[ref[pos-1]]
							pseudoMutCounts[i1][i2]+=1
					elif entry2[0]!="N" and entry2[0]!="O":
						i1=alleles[entry1[0]]
						i2=alleles[entry2[0]]
						pseudoMutCounts[i1][i2]+=1

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




#we know that sample "sample", with partials "newPartials", is best placed as child of node resulting in logLK contribution of newChildLK
# now explore exactly which position near node is best for placement (direct child, sibling or sibling of existing child), and which are the best branch lengths,
#then add the sample at that position of the tree, and update all the internal probability vectors.
#NEW VERSION: now not going through the children or parent, but just attaching to the best place found so far.
def placeSampleOnTreeNew(node,newPartials,sample,bLen,newChildLK,isMidNode, bestUpLK, bestDownLK, bestDownNode,mutMatrix,pseudoMutCounts, adjustBLen):
	# print("\n")
	# print("placeSampleOnTreeNew on node")
	# print(node)
	# print(" with tot:")
	# print(node.probVectTot)

	if adjustBLen:
		factor=float(bLenAdjustment)
	else:
		factor=1.0

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
		while newSplit*node.dist>0.1*bLen:
			probVectParentNew=mergeVectorsUpDown(vectUp,node.dist*newSplit,node.probVect,node.dist*(1.0-newSplit),mutMatrix)
			probChild=appendProb(probVectParentNew,newPartials,bLen*factor,mutMatrix)
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
			while newSplit*node.dist>0.1*bLen:
				probVectParentNew=mergeVectorsUpDown(vectUp,node.dist*(1.0-newSplit),node.probVect,node.dist*newSplit,mutMatrix)
				probChild=appendProb(probVectParentNew,newPartials,bLen*factor,mutMatrix)
				if probChild>bestSplitLK:
					bestSplitLK=probChild
					bestSplit=newSplit
					childBestVect=probVectParentNew
				else:
					bestSplit=1.0-bestSplit
					break
				newSplit=bestSplit/2
		#now try different lengths for the new branch
		childBestVect=shorten(childBestVect)
		LK1=bestSplitLK
		bestLen=bLen*factor
		while bestLen>0.1*bLen:
			newBLen=bestLen/2
			probChild=appendProb(childBestVect,newPartials,newBLen,mutMatrix)
			if probChild>LK1:
				LK1=probChild
				bestLen=newBLen
			else:
				break
		if bestLen>0.7*bLen*factor:
			while bestLen<10*bLen*factor:
				newBLen=bestLen*2
				probChild=appendProb(childBestVect,newPartials,newBLen,mutMatrix)
				if probChild>LK1:
					LK1=probChild
					bestLen=newBLen
				else:
					break
		if bestLen<0.2*bLen:
			LK0=appendProb(childBestVect,newPartials,0.0,mutMatrix)
			if LK0>LK1:
				bestLen=0.0
		#now create new internal node and append child to it

		distTop=node.dist*bestSplit
		newInternalNode=Tree()
		node.up.children[child]=newInternalNode
		newInternalNode.up=node.up
		distBottom=node.dist*(1.0-bestSplit)
		# print("Inside placeSampleOnTreeNew, add node midBranch ")
		# print(sample)
		# print(distBottom)
		# print(bestLen)
		newInternalNode.add_child(node)
		node.up=newInternalNode
		node.dist=distBottom
		#newInternalNode.children.append(node)
		#newInternalNode.add_child(dist=bestLen, name=sample)
		newNode=Tree(name=sample,dist=bestLen)
		#newNode.name=sample
		#newNode.dist=bestLen
		newNode.minorSequences=[]
		newNode.up=newInternalNode
		newInternalNode.add_child(newNode)
		#newInternalNode.children.append(newNode)
		#newInternalNode.children[1].minorSequences=[]
		#newInternalNode.children[1].up=newInternalNode
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
		if distTop>4*bLen/(bLenFactor+thresholdProb):
		#if distTop>bLen/2:
			createFurtherMidNodes(newInternalNode,vectUp,bLen)
		newVect=mergeVectorsUpDown(childBestVect,bestLen,newPartials,0.0,mutMatrix)
		newVect=shorten(newVect)
		newInternalNode.children[1].probVectTot=newVect
		if bestLen>thresholdProb:
			newVect=mergeVectorsUpDown(childBestVect,bestLen/2,newPartials,bestLen/2,mutMatrix)
			newVect=shorten(newVect)
			newInternalNode.children[1].probVectTotUp=newVect
			if bestLen>4*bLen/(bLenFactor+thresholdProb):
			#if bestLen>bLen/2:
				createFurtherMidNodes(newInternalNode.children[1],childBestVect,bLen)
		updatePesudoCounts(childBestVect,newPartials,pseudoMutCounts)
		if verbose:
			print("new internal node added to tree")
			print(newInternalNode.probVect)
			print(newInternalNode.probVectUpRight)
			print(newInternalNode.probVectUpLeft)
			print(newInternalNode.probVectTot)
		updatePartialsFromTop(node,newInternalNode.probVectUpRight,mutMatrix)
		updatePartialsFromBottom(newInternalNode.up,newInternalNode.probVect,child,newInternalNode,mutMatrix)

	#best lk so far is for uppending directly to existing node
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
			while newSplit*bestDownNode.dist>0.1*bLen:
				probVectParentNew=mergeVectorsUpDown(vectUp,bestDownNode.dist*newSplit,bestDownNode.probVect,bestDownNode.dist*(1.0-newSplit),mutMatrix)
				probChild=appendProb(probVectParentNew,newPartials,bLen*factor,mutMatrix)
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
		# print("best down node:")
		# print(bestDownNode)
		# print("bestChildLK: "+str(bestChildLK))


		#if node is root, try to place as sibling of the current root.
		if node.up==None:
			probOldRoot = findProbRoot(node.probVect)
			probVectRoot,probRoot = mergeVectorsRoot(node.probVect,bLen,newPartials,bLen*factor,mutMatrix)
			probRoot+= findProbRoot(probVectRoot)
			parentLKdiff=probRoot-probOldRoot
			bestRootBL=bLen
			parentBestVect=probVectRoot
			newBL=0.5*bLen
			while newBL>0.1*bLen:
				probVectRoot,probRoot = mergeVectorsRoot(node.probVect,newBL,newPartials,bLen*factor,mutMatrix)
				probRoot+= findProbRoot(probVectRoot)
				newDiff=probRoot-probOldRoot
				if newDiff>parentLKdiff:
					parentLKdiff=newDiff
					bestRootBL=newBL
					parentBestVect=probVectRoot
				else:
					break
				newBL=bestRootBL/2
			#print("node is root, parentLKdiff: "+str(parentLKdiff))

		else:
			if node==node.up.children[0]:
				child=0
				vectUp=node.up.probVectUpRight
			else:
				child=1
				vectUp=node.up.probVectUpLeft
			#print("node distance from parent: "+str(node.dist))
			#print("bestUpLK: "+str(bestUpLK))
			bestSplit=0.5
			bestSplitLK=bestUpLK
			parentBestVect=node.probVectTotUp
			newSplit=0.25
			while newSplit*node.dist>0.1*bLen:
				probVectParentNew=mergeVectorsUpDown(vectUp,node.dist*(1.0-newSplit),node.probVect,node.dist*newSplit,mutMatrix)
				probChild=appendProb(probVectParentNew,newPartials,bLen*factor,mutMatrix)
				if probChild>bestSplitLK:
					bestSplitLK=probChild
					bestSplit=newSplit
					parentBestVect=probVectParentNew
				else:
					break
				newSplit=bestSplit/2
			parentLKdiff=bestSplitLK
			bestParentSplit=bestSplit
			#print("node is not root, parentLKdiff: "+str(parentLKdiff))
			#print("newChildLK: "+str(newChildLK))
		
		#add internal node below "node"
		if bestChildLK>=parentLKdiff and bestChildLK>=newChildLK:
			if bestDownNode==bestDownNode.up.children[0]:
				child=0
				vectUp=bestDownNode.up.probVectUpRight
			else:
				child=1
				vectUp=bestDownNode.up.probVectUpLeft
			childBestVect=shorten(childBestVect)

			LK1=bestChildLK
			bestLen=bLen*factor
			while bestLen>0.1*bLen:
				newBLen=bestLen/2
				probChild=appendProb(childBestVect,newPartials,newBLen,mutMatrix)
				if probChild>LK1:
					LK1=probChild
					bestLen=newBLen
				else:
					break
			if bestLen>0.7*bLen*factor:
				while bestLen<10*bLen*factor:
					newBLen=bestLen*2
					probChild=appendProb(childBestVect,newPartials,newBLen,mutMatrix)
					if probChild>LK1:
						LK1=probChild
						bestLen=newBLen
					else:
						break
			if bestLen<0.2*bLen:
				LK0=appendProb(childBestVect,newPartials,0.0,mutMatrix)
				if LK0>LK1:
					bestLen=0.0
			#now create new internal node and append child to it
			newInternalNode=Tree()
			bestDownNode.up.children[child]=newInternalNode
			newInternalNode.up=bestDownNode.up
			distBottom=bestDownNode.dist*(1.0-bestChildSplit)
			distTop=bestDownNode.dist*bestChildSplit
			# print("Inside placeSampleOnTreeNew, add node below node ")
			# print(sample)
			# print(distBottom)
			# print(bestLen)
			#newInternalNode.add_child(dist=distBottom,child=bestDownNode)
			bestDownNode.up=newInternalNode
			bestDownNode.dist=distBottom
			newInternalNode.add_child(bestDownNode)
			#newInternalNode.children.append(bestDownNode)

			#newInternalNode.add_child(dist=bestLen, name=sample)
			newNode=Tree(name=sample,dist=bestLen)
			#newNode.name=sample
			#newNode.dist=bestLen
			newNode.minorSequences=[]
			newNode.up=newInternalNode
			#newInternalNode.children.append(newNode)
			newInternalNode.add_child(newNode)
			#newInternalNode.children[1].minorSequences=[]
			#newInternalNode.children[1].up=newInternalNode
			newInternalNode.dist=distTop
			newInternalNode.children[1].probVect=newPartials
			newVect=mergeVectorsUpDown(vectUp,newInternalNode.dist,newPartials,bestLen,mutMatrix)
			#newVect=shorten(newVect)
			newInternalNode.probVectUpRight=newVect
			newInternalNode.probVectUpLeft=childBestVect
			newVect=mergeVectors(bestDownNode.probVect,bestDownNode.dist,newPartials,bestLen,mutMatrix)
			#newVect=shorten(newVect)
			newInternalNode.probVect=newVect
			newVect=mergeVectorsUpDown(vectUp,distTop/2,newInternalNode.probVect,distTop/2,mutMatrix)
			newVect=shorten(newVect)
			newInternalNode.probVectTotUp=newVect
			newVect=mergeVectorsUpDown(childBestVect,0.0,newPartials,bestLen,mutMatrix)
			newVect=shorten(newVect)
			newInternalNode.probVectTot=newVect
			if distTop>4*bLen/(bLenFactor+thresholdProb):
			#if distTop>bLen/2:
				createFurtherMidNodes(newInternalNode,vectUp,bLen)
			newVect=mergeVectorsUpDown(childBestVect,bestLen,newPartials,0.0,mutMatrix)
			newVect=shorten(newVect)
			newInternalNode.children[1].probVectTot=newVect
			if bestLen>thresholdProb4:
				newVect=mergeVectorsUpDown(childBestVect,bestLen/2,newPartials,bestLen/2,mutMatrix)
				newVect=shorten(newVect)
				newInternalNode.children[1].probVectTotUp=newVect
				if bestLen>4*bLen/(bLenFactor+thresholdProb):
				#if bestLen>bLen/2:
					createFurtherMidNodes(newInternalNode.children[1],childBestVect,bLen)
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
			#new parent is actually part of a polytomy
			if newChildLK>=parentLKdiff:
				bestRootBL=0.0
				bestParentSplit=0.0
				#print("Appending directly at node with tot")
				#print(node.probVectTot)
				#print("with loglk: "+str(newChildLK))
				parentLKdiff=newChildLK
				parentBestVect=node.probVectTot
				if node.up==None:
					probVectRoot,probRoot = mergeVectorsRoot(node.probVect,0.0,newPartials,bLen*factor,mutMatrix)
					parentBestVect=probVectRoot

			#add parent to the root
			if node.up==None:

				#now try different lengths for right branch
				bestLen2=bLen*factor
				while bestLen2>0.1*bLen:
					newBLen=bestLen2/2
					newProbVectRoot,newProbRoot = mergeVectorsRoot(node.probVect,bestRootBL,newPartials,newBLen,mutMatrix)
					newProbRoot+= findProbRoot(newProbVectRoot)
					LKdiffRoot=newProbRoot-probOldRoot
					if LKdiffRoot>parentLKdiff:
						parentLKdiff=LKdiffRoot
						bestLen2=newBLen
						probVectRoot=newProbVectRoot
					else:
						break
				if bestLen2>0.7*bLen*factor:
					while bestLen2<10*bLen*factor:
						newBLen=bestLen2*2
						newProbVectRoot,newProbRoot = mergeVectorsRoot(node.probVect,bestRootBL,newPartials,newBLen,mutMatrix)
						newProbRoot+= findProbRoot(newProbVectRoot)
						LKdiffRoot=newProbRoot-probOldRoot
						if LKdiffRoot>parentLKdiff:
							parentLKdiff=LKdiffRoot
							bestLen2=newBLen
							probVectRoot=newProbVectRoot
						else:
							break
				if bestLen2<0.2*bLen:
					newProbVectRoot,newProbRoot = mergeVectorsRoot(node.probVect,bestRootBL,newPartials,0.0,mutMatrix)
					newProbRoot+= findProbRoot(newProbVectRoot)
					LK0=newProbRoot-probOldRoot
					if LK0>parentLKdiff:
						bestLen2=0.0
						parentLKdiff=LK0
						probVectRoot=newProbVectRoot

				#print("new root to be added to tree")
				#print(node.name)
				#print(sample)
				#print(node.probVect)
				#print(node.children)
				#print(node.probVectUpRight)
				#print(newPartials)
				#print(bestRootBL)
				#print(bestLen2)
				newRoot=Tree()
				newRoot.name="newRoot"
				# print("Inside placeSampleOnTreeNew, add node above root ")
				# print(sample)
				# print(bestRootBL)
				# print(bestLen2)
				newRoot.probVect=probVectRoot
				newVect=rootVector0(probVectRoot,rootFreqs)
				newVect=shorten(newVect)
				newRoot.probVectTot=newVect
				newRoot.probVectUpRight=rootVector(newPartials,rootFreqs,bestLen2)
				newRoot.probVectUpLeft=rootVector(node.probVect,rootFreqs,bestRootBL)
				#newRoot.add_child(child=node,dist=bestRootBL)
				node.up=newRoot
				node.dist=bestRootBL
				newRoot.add_child(node)
				#(newRoot.children).append(node)
				#newRoot.add_child(dist=bestLen2, name=sample)
				newNode=Tree(name=sample,dist=bestLen2)
				#newNode.name=sample
				#newNode.dist=bestLen2
				newNode.minorSequences=[]
				newNode.up=newRoot
				#(newRoot.children).append(newNode)
				newRoot.add_child(newNode)
				newRoot.children[1].probVect=newPartials
				newVect=mergeVectorsUpDown(newRoot.probVectUpLeft,bestLen2,newPartials,0.0,mutMatrix)
				newVect=shorten(newVect)
				newRoot.children[1].probVectTot=newVect
				if bestLen2>thresholdProb4:
					newVect=mergeVectorsUpDown(newRoot.probVectUpLeft,bestLen2/2,newPartials,bestLen2/2,mutMatrix)
					newVect=shorten(newVect)
					newRoot.children[1].probVectTotUp=newVect
					if bestLen2>4*bLen/(bLenFactor+thresholdProb):
					#if bestLen2>bLen/2:
						createFurtherMidNodes(newRoot.children[1],newRoot.probVectUpLeft,bLen)
				#newRoot.children[1].minorSequences=[]
				#newRoot.children[1].up=newRoot
				if verbose:
					print("new root added to tree")
					print(newRoot.probVect)
					print(newRoot.children[0].probVect)
					print(newRoot.children[1].probVect)
				# print("new root added to tree")
				# print(newRoot.name)
				# print(newRoot.children[0].name)
				# print(newRoot.children[1].name)
				# print(newRoot.dist)
				# print(newRoot.children[0].dist)
				# print(newRoot.children[1].dist)
				# print(newRoot.probVect)
				# print(newRoot.children[0].probVect)
				# print(newRoot.children[1].probVect)
				# print("updating partials from top using probVectUpRight")
				# print("children of child")
				# print((newRoot.children[1]).children)
				# print(newRoot.probVectUpRight)
				# print("name of node: "+node.name)
				updatePartialsFromTop(node,newRoot.probVectUpRight,mutMatrix)
				# if sample=="EPI_ISL_1113446":
				# 	print(newRoot.probVectTot)
				# 	print(newRoot.probVectUpRight)
				# 	print(newRoot.probVectUpLeft)
				# 	print(newRoot.probVect)
				# 	print(bestRootBL)
				# 	print(node.probVect)
				# 	print(bestLen2)
				# 	exit()
				return newRoot

			#add parent to node
			else:
				if node==node.up.children[0]:
					child=0
					vectUp=node.up.probVectUpRight
				else:
					child=1
					vectUp=node.up.probVectUpLeft
				
				#now try different lengths for the new branch
				LK1=parentLKdiff
				bestLen=bLen*factor
				while bestLen>0.1*bLen:
					newBLen=bestLen/2
					probChild=appendProb(parentBestVect,newPartials,newBLen,mutMatrix)
					if probChild>LK1:
						LK1=probChild
						bestLen=newBLen
					else:
						break
				if bestLen>0.7*bLen*factor:
					while bestLen<10*bLen*factor:
						newBLen=bestLen*2
						probChild=appendProb(parentBestVect,newPartials,newBLen,mutMatrix)
						if probChild>LK1:
							LK1=probChild
							bestLen=newBLen
						else:
							break
				if bestLen<0.2*bLen:
					LK0=appendProb(parentBestVect,newPartials,0.0,mutMatrix)
					if LK0>LK1:
						bestLen=0.0
						#print("Appending with child branch 0, with loglk LK0: "+str(LK0))
						#print(parentBestVect)
						#print(newPartials)
				#now create new internal node and append child to it
				# print("Inside placeSampleOnTreeNew, add node above node ")
				# print(sample)
				# print(bestParentSplit)
				# print(bestLen)
				# print(parentBestVect)
				# print("\n")
				# print("node.up.probVectTot")
				# print(node.up.probVectTot)
				# print(node.up.is_root())
				# print(node.up.up==None)
				# print("node.up.probVectUpRight and Left")
				# print(node.up.probVectUpRight)
				# print(node.up.probVectUpLeft)
				# print("node.probVect and sibling")
				# print(node.up.children[0].probVect)
				# print(node.up.children[1].probVect)
				# if node.up.up!=None:
				# 	print("node.up.up.probVectUpRight and Left")
				# 	print(node.up.up.probVectUpRight)
				# 	print(node.up.up.probVectUpLeft)
				newInternalNode=Tree()
				node.up.children[child]=newInternalNode
				newInternalNode.up=node.up
				distBottom=node.dist*bestParentSplit
				distTop=node.dist*(1.0-bestParentSplit)
				#print(node.dist)
				#print(bestParentSplit)
				#newInternalNode.add_child(dist=distBottom,child=node)
				node.dist=distBottom
				#newInternalNode.children.append(node)
				node.up=newInternalNode
				newInternalNode.add_child(node)
				#newInternalNode.add_child(dist=bestLen, name=sample)
				newNode=Tree(name=sample,dist=bestLen)
				#newNode.name=sample
				#newNode.dist=bestLen
				newNode.minorSequences=[]
				newNode.up=newInternalNode
				newInternalNode.add_child(newNode)
				#newInternalNode.children.append(newNode)
				#print(node.dist)
				#print(bestLen)
				#newInternalNode.children[1].minorSequences=[]
				#newInternalNode.children[1].up=newInternalNode
				newInternalNode.dist=distTop
				newInternalNode.children[1].probVect=newPartials
				newVect=mergeVectorsUpDown(vectUp,newInternalNode.dist,newPartials,bestLen,mutMatrix)
				#newVect=shorten(newVect)
				newInternalNode.probVectUpRight=newVect
				newInternalNode.probVectUpLeft=parentBestVect
				# print("node.probVectTot")
				# print(node.probVectTot)
				# print("node.probVect")
				# print(node.probVect)
				# print("newPartials")
				# print(newPartials)
				# print("vectUp")
				# print(vectUp)
				newVect=mergeVectors(node.probVect,node.dist,newPartials,bestLen,mutMatrix)
				#newVect=shorten(newVect)
				newInternalNode.probVect=newVect
				newVect=mergeVectorsUpDown(vectUp,distTop/2,newInternalNode.probVect,distTop/2,mutMatrix)
				newVect=shorten(newVect)
				newInternalNode.probVectTotUp=newVect
				newVect=mergeVectorsUpDown(parentBestVect,0.0,newPartials,bestLen,mutMatrix)
				newVect=shorten(newVect)
				newInternalNode.probVectTot=newVect
				if distTop>4*bLen/(bLenFactor+thresholdProb):
				#if distTop>bLen/2:
					createFurtherMidNodes(newInternalNode,vectUp,bLen)
				newVect=mergeVectorsUpDown(parentBestVect,bestLen,newPartials,0.0,mutMatrix)
				newVect=shorten(newVect)
				newInternalNode.children[1].probVectTot=newVect
				if bestLen>thresholdProb4:
					newVect=mergeVectorsUpDown(parentBestVect,bestLen/2,newPartials,bestLen/2,mutMatrix)
					newVect=shorten(newVect)
					newInternalNode.children[1].probVectTotUp=newVect
					if bestLen>4*bLen/(bLenFactor+thresholdProb):
					#if bestLen>bLen/2:
						createFurtherMidNodes(newInternalNode.children[1],parentBestVect,bLen)
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






def updateCumulativeNonMutationProb(cumulativeRate):
	#cumulativeRate=[0.0]
	for i in range(lRef):
		ind=allelesLow[ref[i]]
		cumulativeRate[i+1]=cumulativeRate[i]+mutMatrix[ind][ind]
		#cumulativeRate.append(mutMatrix[ind][ind])





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

timeFinding=0.0
timePlacing=0.0

numSamples=0
for d in distances:
	numSamples+=1
	sample=d[1]
	newPartials=probVectTerminalNode(data[sample],0.0)
	#print("\n\n\n\n\n")
	#print(sample)
	if (numSamples%40)==0:
		oldMatrix=mutMatrix
		mutMatrix=updateSubMatrix(pseudoMutCounts,model)
		if abs(oldMatrix[0][0]-mutMatrix[0][0])>0.001 or abs(oldMatrix[1][1]-mutMatrix[1][1])>0.001 or abs(oldMatrix[2][2]-mutMatrix[2][2])>0.001 or abs(oldMatrix[3][3]-mutMatrix[3][3])>0.001:
			updateCumulativeNonMutationProb(cumulativeRate)

		#print(mutMatrix)
	if (numSamples%500)==0:
		print("Sample num "+str(numSamples))
		#print()
		#print(sample)
		#print(newPartials)
	start=time()
	#node , bestNewLK, isMidNode, bestUpLK, bestDownLK, bestDownNode, LKdiff2, directChildNode, adjustBLen=findBestParent(t1,newPartials,sample,bLen,float('-inf'),t1,0,False,float('-inf'),float('-inf'),None,float('-inf'),mutMatrix, False)
	node , bestNewLK, isMidNode, bestUpLK, bestDownLK, bestDownNode, adjustBLen=findBestParent(t1,newPartials,sample,bLen,mutMatrix)
	timeFinding+=(time()-start)
	if bestNewLK<0.5:
		start=time()
		newRoot=placeSampleOnTreeNew(node,newPartials,sample,bLen,bestNewLK,isMidNode, bestUpLK, bestDownLK, bestDownNode,mutMatrix,pseudoMutCounts, adjustBLen)
		if newRoot!=None:
			t1=newRoot
		timePlacing+=(time()-start)
	#print("Updated tree")
	#print(t1)
	#print("\n")

	#if sample=="EPI_ISL_1058154":
	#	break

def createNewick(node):
	if node.children:
		return "("+createNewick(node.children[0])+","+createNewick(node.children[1])+"):"+str(node.dist)
	else:
		if len(node.minorSequences)>0:
			newList=["(",node.name,":0"]
			for s2 in node.minorSequences:
				newList.append(","+s2+":0")
			newList.append("):"+str(node.dist))
			return "".join(newList)
			#treeString=treeString.replace(s.name+":",newString)
		else:
			return node.name+":"+str(node.dist)

# print("Tree finalized")
# treeString=(t1.write()).replace(")1:","):")
# totMinors=0
# for s in t1.get_leaves():
# 	if len(s.minorSequences)>0:
# 		totMinors+=len(s.minorSequences)
# 		newString="("+s.name+":0"
# 		for s2 in s.minorSequences:
# 			newString+=","+s2+":0"
# 		newString+="):"
# 		treeString=treeString.replace(s.name+":",newString)
#print("leaves in tree: ")
#print(len(t1.get_leaves()))
#print("tot minors:")
#print(totMinors)

newickString=createNewick(t1)

#print(treeString)
file=open(outputFile+"_tree.tree","w")
file.write(newickString+";")
file.close()
file=open(outputFile+"_subs.txt","w")
for i in range4:
	for j in range4:
		file.write(str(mutMatrix[i][j])+"\t")
	file.write("\n")
file.close()
print("Missed minor samples: "+str(totalMissedMinors[0]))
print("Final Substitution matrix:")
print(print(mutMatrix))
print("Time spent finding placement nodes: "+str(timeFinding))
print("Time spent updating the tree: "+str(timePlacing))
exit()



#optimize branch lengths?
#Optimize root position?
#Calculate pairwise distances?
#try to improve topology?
















