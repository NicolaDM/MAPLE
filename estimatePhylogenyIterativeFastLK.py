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
parser.add_argument("--model", help="Which substitution model should be used. Allowed models so far are JC, GTR (default) or UNREST.", default="GTR")
args = parser.parse_args()

bLenAdjustment=10.0

onlyNambiguities=args.onlyNambiguities
useLogs=args.useLogs
thresholdProb=args.thresholdProb
verbose=args.verbose
pathSimu=args.path
inputFile=pathSimu+args.input
outputFile=pathSimu+args.output
refFile=pathSimu+args.reference
allowedFails=args.allowedFails
model=args.model

example=False

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
	cumulativeBases=np.zeros((lRef+1,4))
	for i in range(lRef):
		for k in range(4):
			cumulativeBases[i+1][k]=cumulativeBases[i][k]
		cumulativeBases[i+1][allelesLow[ref[i]]]+=1
	#print("cumulativeBases")
	#print(cumulativeBases[-1])
	rootFreqs=np.zeros(4)
	rootFreqsLog=np.zeros(4)
	for i in range(4):
		rootFreqs[i]=cumulativeBases[-1][i]/float(lRef)
		rootFreqsLog[i]=math.log(rootFreqs[i])
	#print("ref base frequencies and log")
	#print(rootFreqs)
	#print(rootFreqsLog)
	return ref, cumulativeBases, rootFreqs, rootFreqsLog

ref, cumulativeBases, rootFreqs, rootFreqsLog = collectReference(refFile)
lRef=len(ref)
if model=="JC":
	rootFreqs=[0.25,0.25,0.25,0.25]
	rootFreqsLog=[math.log(0.25),math.log(0.25),math.log(0.25),math.log(0.25)]

if example:
	ref="aaccggtt"
	lRef=len(ref)
	rootFreqs=[0.25,0.25,0.25,0.25]
	rootFreqsLog=[math.log(0.25),math.log(0.25),math.log(0.25),math.log(0.25)]
	cumulativeBases=np.zeros((lRef+1,4))
	for i in range(lRef):
		for k in range(4):
			cumulativeBases[i+1][k]=cumulativeBases[i][k]
		cumulativeBases[i+1][allelesLow[ref[i]]]+=1

#TO DO
#recreate a simple scenario of a recombinant sample and with just a few samples and small genome, and see if there is still a bias in placing at the root, and why that is.



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

if example:
	data={}
	data["sampleCG"]=[]
	data["sampleTG"]=[("t",3)]
	data["sampleCT"]=[("t",5)]
	data["sampleTT"]=[("t",3),("t",5)]

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

def updateSubMatrix(pseudoMutCounts,model):
	mutMatrix=np.zeros((4,4))
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
def appendProb(probVectP,probVectC,bLen,mutMatrix):
	Lkcost=0.0
	#LkcostOriginal=0.0
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
				if entry1[4]:
					for i in range4:
						Lkcost+=mutMatrix[i][i]*(bLen+entry1[3]+entry1[5])*(cumulativeBases[end][i]-cumulativeBases[pos-1][i])
				else:
					for i in range4:
						Lkcost+=mutMatrix[i][i]*(bLen+entry1[3])*(cumulativeBases[end][i]-cumulativeBases[pos-1][i])
			elif entry2[0]=="O":
				i1=allelesLow[ref[pos-1]]
				tot=0.0
				if entry1[4]:
					for i in range4:
						if i1==i:
							tot2=rootFreqs[i]*(1.0+mutMatrix[i][i1]*entry1[3])
						else:
							tot2=rootFreqs[i]*mutMatrix[i][i1]*entry1[3]
						tot3=0.0
						for j in range4:
							if i!=j:
								tot3+=mutMatrix[i][j]*(bLen+entry1[5])*entry2[4][j]
							else:
								tot3+=(1.0+mutMatrix[i][j]*(bLen+entry1[5]))*entry2[4][j]
						tot2=tot2*tot3
						tot+=tot2
					if tot>sys.float_info.min:
						Lkcost+=math.log(tot/rootFreqs[i1])
					else:
						return float("-inf")
				else:
					for j in range4:
						if i1!=j:
							tot+=mutMatrix[i1][j]*(bLen+entry1[3])*entry2[4][j]
						else:
							tot+=(1.0+mutMatrix[i1][j]*(bLen+entry1[3]))*entry2[4][j]
					if tot>sys.float_info.min:
						Lkcost+=math.log(tot)
					else:
						return float("-inf")
			else:
				if entry1[4]:
					if (bLen+entry1[3]+entry1[5])<thresholdProb2:
						return float("-inf")
					i2=alleles[entry2[0]]
					i1=allelesLow[ref[pos-1]]
					tot=0.0
					for i in range4:
						if i1==i:
							tot+=rootFreqs[i]*(1.0+mutMatrix[i][i1]*entry1[3])*mutMatrix[i][i2]*(entry1[5]+bLen)
						elif i2==i:
							tot+=rootFreqs[i]*mutMatrix[i][i1]*entry1[3]*(1.0+mutMatrix[i][i2]*(entry1[5]+bLen))
						else:
							tot+=rootFreqs[i]*mutMatrix[i][i1]*entry1[3]*mutMatrix[i][i2]*(entry1[5]+bLen)
					Lkcost+=math.log(tot/rootFreqs[i1])
				else:
					if (bLen+entry1[3])<thresholdProb2:
						return float("-inf")
					i2=alleles[entry2[0]]
					i1=allelesLow[ref[pos-1]]
					Lkcost+=math.log(mutMatrix[i1][i2]*(bLen+entry1[3]))
		elif entry1[0]=="O":
			tot1=0.0
			for i in range4:
				tot1+=entry1[4][i]
			#LkcostOriginal+=math.log(tot1)
			if entry2[0]=="O":
				tot=0.0
				for i in range4:
					for j in range4:
						if i==j:
							tot+=(1.0+mutMatrix[i][i]*(bLen+entry1[3]))*entry1[4][i]*entry2[4][j]
						else:
							tot+=mutMatrix[i][j]*(bLen+entry1[3])*entry1[4][i]*entry2[4][j]
				if tot>sys.float_info.min:
					Lkcost+=math.log(tot/tot1)
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
						tot+=entry1[4][i]*(1.0+mutMatrix[i][i]*(bLen+entry1[3]))
					else:
						tot+=entry1[4][i]*mutMatrix[i][i2]*(bLen+entry1[3])
				if tot>sys.float_info.min:
					Lkcost+=math.log(tot/tot1)
				else:
					return float("-inf")
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
					tot=0.0
					if entry1[4]:
						for i in range4:
							if i1==i:
								tot2=rootFreqs[i]*(1.0+mutMatrix[i][i1]*entry1[3])
							else:
								tot2=rootFreqs[i]*mutMatrix[i][i1]*entry1[3]
							tot3=0.0
							for j in range4:
								if i!=j:
									tot3+=mutMatrix[i][j]*(bLen+entry1[5])*entry2[4][j]
								else:
									tot3+=(1.0+mutMatrix[i][j]*(bLen+entry1[5]))*entry2[4][j]
							tot2=tot2*tot3
							tot+=tot2
						if tot>sys.float_info.min:
							Lkcost+=math.log(tot/rootFreqs[i1])
						else:
							return float("-inf")
					else:
						for j in range4:
							if i1==j:
								tot+=(1.0+mutMatrix[i1][i1]*(bLen+entry1[3]))*entry2[4][j]
							else:
								tot+=mutMatrix[i1][j]*(bLen+entry1[3])*entry2[4][j]
						if tot>sys.float_info.min:
							Lkcost+=math.log(tot)
						else:
							return float("-inf")

				else:
					if entry1[4]:
						if (bLen+entry1[3]+entry1[5])<thresholdProb2:
							return float("-inf")
						if entry2[0]=="R":
							i2=allelesLow[ref[pos-1]]
						else:
							i2=alleles[entry2[0]]
						tot=0.0
						for i in range4:
							if i1==i:
								tot+=rootFreqs[i]*(1.0+mutMatrix[i][i1]*entry1[3])*mutMatrix[i][i2]*(entry1[5]+bLen)
							elif i2==i:
								tot+=rootFreqs[i]*mutMatrix[i][i1]*entry1[3]*(1.0+mutMatrix[i][i2]*(entry1[5]+bLen))
							else:
								tot+=rootFreqs[i]*mutMatrix[i][i1]*entry1[3]*mutMatrix[i][i2]*(entry1[5]+bLen)
						Lkcost+=math.log(tot/rootFreqs[i1])
					else:
						if (bLen+entry1[3])<thresholdProb2:
							return float("-inf")
						if entry2[0]=="R":
							i2=allelesLow[ref[pos-1]]
						else:
							i2=alleles[entry2[0]]
						Lkcost+=math.log(mutMatrix[i1][i2]*(bLen+entry1[3]))

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
	return Lkcost 
	#-LkcostOriginal

#function to calculate likelihood cost of appending child to parent
def appendProbOld(probVectP,probVectC,bLen,mutMatrix):
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
def findBestParent(t1,diffs,sample,bLen,bestLKdiff,bestNodeSoFar,failedPasses,bestIsMidNode,bestUpLK,bestDownLK,bestDownNode,parentLK,mutMatrix,adjustBLen):
	# print("findBestParent(), Node")
	# print(t1)
	# print("With probVectTot: ")
	# print(t1.probVectTot)
	if t1.is_leaf():
		probVect=t1.probVect
		comparison=isMinorSequence(probVect,diffs)
		if comparison==1:
			t1.minorSequences.append(sample)
			return t1, 1.0, False, float("-inf"), float("-inf"), None, float("-inf"), None, False
		elif comparison==2:
			if verbose:
				print("Second sequence is more informative than first, but for now the two sequences are not merged.")
				print(t1.name)
				print(sample)
				#print(probVect)
				#print(diffs)
			totalMissedMinors[0]+=1
	
	if t1.dist>thresholdProb2 and t1.up!=None:
		probVect2=t1.probVectTotUp
		LKdiff2=appendProb(probVect2,diffs,bLen,mutMatrix)
		newLKdiff2=appendProb(probVect2,diffs,bLen*bLenAdjustment,mutMatrix)
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
		probVect=t1.probVectTot
		LKdiff=appendProb(probVect,diffs,bLen,mutMatrix)
		if verbose:
			print("Trying a new parent")
			if t1.is_leaf():
				print(t1.name)
			print(probVect)
			print("Failed passes before this step: "+str(failedPasses)+" logLK score: "+str(LKdiff))
		newLKdiff=appendProb(probVect,diffs,bLen*bLenAdjustment,mutMatrix)
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
		elif LKdiff2>=bestLKdiff-thresholdProb2:
			bestUpLK=parentLK
			bestDownLK=LKdiff
			bestDownNode=t1
		elif LKdiff<bestLKdiff-0.1:
			failedPasses+=1
			if failedPasses>allowedFails:
				return bestNodeSoFar , bestLKdiff , bestIsMidNode, bestUpLK, bestDownLK, bestDownNode, LKdiff2, t1, adjustBLen
	else:
		LKdiff=float("-inf")
	# print("t1.dist "+str(t1.dist))
	# print("LKdiff "+str(LKdiff))
	# print("LKdiff2 "+str(LKdiff2))
	

	bestOfChildLK=float("-inf")
	bestChild=None
	# print("Before child loop, best child")
	# print(bestChild)
	# print(bestOfChildLK)
	for c in t1.children:
		node, score , childIsMidNode, childUpLK, childDownLK, childBestDownNode, currentBestChildLK, currentBestChild , adjustBLenChild = findBestParent(c,diffs,sample,bLen,bestLKdiff,bestNodeSoFar,failedPasses,bestIsMidNode,bestUpLK,bestDownLK,bestDownNode,LKdiff,mutMatrix,adjustBLen)
		if score>0.5:
			return node, score , False, float("-inf"), float("-inf"), None, float("-inf"), None, False
		if score>bestLKdiff:
			bestLKdiff=score
			bestNodeSoFar=node
			bestIsMidNode=childIsMidNode
			bestUpLK=childUpLK
			bestDownLK=childDownLK
			bestDownNode=childBestDownNode
			adjustBLen=adjustBLenChild
		if currentBestChildLK>bestOfChildLK:
			bestOfChildLK=currentBestChildLK
			bestChild=currentBestChild
		# print("Inside child loop, best child")
		# print(bestChild)
		# print(bestOfChildLK)
	# print("End of child loop for node")
	# print(t1)
	# print(t1.dist)
	# print("best child")
	# print(bestChild)
	# if bestChild!=None:
	# 	print(bestChild.dist)
	# print(bestOfChildLK)
	if LKdiff>=bestLKdiff-thresholdProb2:
		bestDownLK=bestOfChildLK
		bestDownNode=bestChild
	# print("best node in absolute")
	# print(bestNodeSoFar)
	# print(bestIsMidNode)
	# print("with lk "+str(bestLKdiff))
	# print("LKdiff: "+str(LKdiff))
	# print("with best child ")
	# print(bestDownNode)
	# print("with lk "+str(bestDownLK))
	# print("while current best child")
	# print(bestChild)
	# print("with lk "+str(bestOfChildLK))
	# print("\n\n")
	if t1.dist>thresholdProb2 and t1.up!=None:
		if t1.up.children[0]==t1:
			vectUp=t1.up.probVectUpRight
		else:
			vectUp=t1.up.probVectUpLeft
		if t1.dist>bLen/2:
			newBLen2=t1.dist/4
			bestLKdiff2=LKdiff2
			while newBLen2>bLen/4:
				newProbVect2=mergeVectorsUpDown(vectUp,newBLen2,t1.probVect,t1.dist-newBLen2,mutMatrix)
				newLKdiff2=appendProb(newProbVect2,diffs,bLen,mutMatrix)
				newLKdiff3=appendProb(newProbVect2,diffs,bLen*bLenAdjustment,mutMatrix)
				if newLKdiff3>newLKdiff2:
					newLKdiff2=newLKdiff3
					adjusted=True
				else:
					adjusted=False
				if newLKdiff2>bestLKdiff2:
					bestLKdiff2=newLKdiff2
					adjustBLen=adjusted
				else:
					break
				newBLen2=newBLen2/2
			return bestNodeSoFar , bestLKdiff , bestIsMidNode, bestUpLK, bestDownLK, bestDownNode, bestLKdiff2, t1, adjustBLen
		else:	
			return bestNodeSoFar , bestLKdiff , bestIsMidNode, bestUpLK, bestDownLK, bestDownNode, LKdiff2, t1, adjustBLen
	else:
		return bestNodeSoFar , bestLKdiff , bestIsMidNode, bestUpLK, bestDownLK, bestDownNode, bestOfChildLK, bestChild, adjustBLen



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
					#if entry2[0]!="O" and entry2[4]:
					#	totLen2+=entry2[5]
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
						if state=="O":
							probVect.append([state,pos,1,0.0,newVec])
						else:
							probVect.append([state,pos,1,0.0,False])
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
					if state=="O":
						probVect.append([state,pos,1,0.0,newVec])
					else:
						probVect.append([state,pos,1,0.0,False])
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
							if state=="O":
								probVect.append([state,pos,1,0.0,newVec])
							else:
								probVect.append([state,pos,1,0.0,False])
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
					if vect[i1]<0.9999:
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
				if vect[i2]<0.9999:
					return True

		elif entry1[0]=="O":
			vect1=entry1[4]/np.sum(entry1[4])
			vect2=entry2[4]/np.sum(entry2[4])
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
			if diffVal>bLen/10:
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
	if len(node.children)==0:
		newTot=mergeVectorsUpDown(probVectUp,node.dist,node.probVect,0.0,mutMatrix)
		if newTot==None:
			if node.up.children[0]==node:
				updateBLen(node.up,0,bLen,mutMatrix)
			else:
				updateBLen(node.up,1,bLen,mutMatrix)
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
			if not updatedTot:
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
						updateBLen(node,1,bLen,mutMatrix)
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
						updateBLen(node,0,bLen,mutMatrix)
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
		factor=bLenAdjustment
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
		#now create new internal node and append child to it
		newInternalNode=Tree()
		node.up.children[child]=newInternalNode
		newInternalNode.up=node.up
		distBottom=node.dist*(1.0-bestSplit)
		distTop=node.dist*bestSplit
		# print("Inside placeSampleOnTreeNew, add node midBranch ")
		# print(sample)
		# print(distBottom)
		# print(bestLen)
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
		if node.is_root():
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
			newInternalNode.add_child(dist=distBottom,child=bestDownNode)
			bestDownNode.up=newInternalNode
			newInternalNode.add_child(dist=bestLen, name=sample)
			newInternalNode.children[1].minorSequences=[]
			newInternalNode.children[1].up=newInternalNode
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
			newVect=mergeVectorsUpDown(childBestVect,bestLen,newPartials,0.0,mutMatrix)
			newVect=shorten(newVect)
			newInternalNode.children[1].probVectTot=newVect
			if bestLen>thresholdProb4:
				newVect=mergeVectorsUpDown(childBestVect,bestLen/2,newPartials,bestLen/2,mutMatrix)
				newVect=shorten(newVect)
				newInternalNode.children[1].probVectTotUp=newVect
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
				parentLKdiff=newChildLK
				parentBestVect=node.probVectTot
				if node.is_root():
					probVectRoot,probRoot = mergeVectorsRoot(node.probVect,0.0,newPartials,bLen*factor,mutMatrix)
					parentBestVect=probVectRoot

			#add parent to the root
			if node.is_root():

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

				newRoot=Tree()
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
				newRoot.add_child(child=node,dist=bestRootBL)
				node.up=newRoot
				newRoot.add_child(dist=bestLen2, name=sample)
				newRoot.children[1].probVect=newPartials
				newVect=mergeVectorsUpDown(newRoot.probVectUpLeft,bestLen2,newPartials,0.0,mutMatrix)
				newVect=shorten(newVect)
				newRoot.children[1].probVectTot=newVect
				if bestLen2>thresholdProb4:
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
				#now create new internal node and append child to it
				# print("Inside placeSampleOnTreeNew, add node above node ")
				# print(sample)
				# print(bestParentSplit)
				# print(bestLen)
				# print(parentBestVect)
				# print("\n")
				newInternalNode=Tree()
				node.up.children[child]=newInternalNode
				newInternalNode.up=node.up
				distBottom=node.dist*bestParentSplit
				distTop=node.dist*(1.0-bestParentSplit)
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
				newInternalNode.probVectUpLeft=parentBestVect
				newVect=mergeVectors(node.probVect,node.dist,newPartials,bestLen,mutMatrix)
				#newVect=shorten(newVect)
				newInternalNode.probVect=newVect
				newVect=mergeVectorsUpDown(vectUp,distTop/2,newInternalNode.probVect,distTop/2,mutMatrix)
				newVect=shorten(newVect)
				newInternalNode.probVectTotUp=newVect
				newVect=mergeVectorsUpDown(parentBestVect,0.0,newPartials,bestLen,mutMatrix)
				newVect=shorten(newVect)
				newInternalNode.probVectTot=newVect
				newVect=mergeVectorsUpDown(parentBestVect,bestLen,newPartials,0.0,mutMatrix)
				newVect=shorten(newVect)
				newInternalNode.children[1].probVectTot=newVect
				if bestLen>thresholdProb4:
					newVect=mergeVectorsUpDown(parentBestVect,bestLen/2,newPartials,bestLen/2,mutMatrix)
					newVect=shorten(newVect)
					newInternalNode.children[1].probVectTotUp=newVect
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
	if (numSamples%5)==0:
		mutMatrix=updateSubMatrix(pseudoMutCounts,model)
		#print(mutMatrix)
	if (numSamples%500)==0:
		print("Sample num "+str(numSamples))
		#print()
		#print(sample)
		#print(newPartials)
	node , bestNewLK, isMidNode, bestUpLK, bestDownLK, bestDownNode, LKdiff2, directChildNode, adjustBLen=findBestParent(t1,newPartials,sample,bLen,float('-inf'),t1,0,False,float('-inf'),float('-inf'),None,float('-inf'),mutMatrix, False)
	if bestNewLK<0.5:
		newRoot=placeSampleOnTreeNew(node,newPartials,sample,bLen,bestNewLK,isMidNode, bestUpLK, bestDownLK, bestDownNode,mutMatrix,pseudoMutCounts, adjustBLen)
		if newRoot!=None:
			t1=newRoot
	#print("Updated tree")
	#print(t1)
	#print("\n")

	#if sample=="EPI_ISL_1058154":
	#	break

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
print("Final Substitution matrix:")
print(print(mutMatrix))
exit()


#estimate susbstitution rates?
#optimize branch lengths?
#Optimize root position?
#Calculate pairwise distances?
#try to improve topology?
















