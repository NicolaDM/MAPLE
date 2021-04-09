import sys
import os
import math
import numpy as np
import random
#import os.path
from os import path
import argparse
#import matplotlib
#import matplotlib.pyplot as plt
#from scipy.stats import ttest_ind, ttest_ind_from_stats, chi2_contingency
#from scipy.special import stdtr
#from discreteMarkovChain import markovChain
from Bio.Data import CodonTable
table = CodonTable.ambiguous_dna_by_id[1]
from Bio.Seq import _translate_str
import time

from ete3 import Tree

#Work on efficiently calculating phylogenetic likelihoods.
#Example run command line: python3 /Users/demaio/Desktop/GISAID-hCoV-19-phylogeny-2021-03-12/createDiffsFile.py --path /Users/demaio/Desktop/GISAID-hCoV-19-phylogeny-2021-03-12/

parser = argparse.ArgumentParser(description='Run estimation of mutation rates and synonymous site selection from SARS-CoV-2 data, and generate plots.')
parser.add_argument('--path',default="", help='path where to find files and plot results.')
parser.add_argument("--createDiffFile", help="create a concise file replacing a MSA.", action="store_true")
parser.add_argument("--reduceDiffFile", help="Reduce the size of a concise file by only including samples in the tree.", action="store_true")
parser.add_argument("--onlyNambiguities", help="Treat all ambiguities as N (total missing information).", action="store_true")
parser.add_argument("--useLogs", help="Calculate logarithms of non-mutation probabilities, otherwise approximate them.", action="store_true")
parser.add_argument("--thresholdProb",help="relative probability threshold used to ignore possible states with very low probabilities.",  type=float, default=0.000001)
parser.add_argument("--verbose", help="Print to screen a lot of stuff.", action="store_true")
parser.add_argument("--runAllData", help="Run the fast likelihood on the huge tree (takes about 2 minutes).", action="store_true")
parser.add_argument("--subsample",help="Sample this number of samples from the original dataset and run likelihood comparison with PhyML.",  type=int, default=0)
parser.add_argument("--nRepeats",help="When subsampling the tree, this is the number of repeats to be done.",  type=int, default=1)
parser.add_argument("--seed",help="Seed for random number generator.",  type=int, default=1)
args = parser.parse_args()

createDiffFile=args.createDiffFile
reduceDiffFile=args.reduceDiffFile
onlyNambiguities=args.onlyNambiguities
useLogs=args.useLogs
thresholdProb=args.thresholdProb
verbose=args.verbose
runAllData=args.runAllData
subsample=args.subsample
nRepeats=args.nRepeats
seed=args.seed

random.seed(a=seed)

pathSimu=args.path

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
	#print(cumulativeBases)
	print(cumulativeBases[-1])
	rootFreqs=np.zeros(4)
	rootFreqsLog=np.zeros(4)
	for i in range(4):
		rootFreqs[i]=cumulativeBases[-1][i]/float(lRef)
		rootFreqsLog[i]=math.log(rootFreqs[i])
	print(rootFreqs)
	print(rootFreqsLog)
	return ref, cumulativeBases, rootFreqs, rootFreqsLog


ref, cumulativeBases, rootFreqs, rootFreqsLog = collectReference(pathSimu+"EPI_ISL_402124_lowercase.fasta")
lRef=len(ref)

if createDiffFile:
	start = time.time()
	#collect alignment and translate into diff file
	fileI=open(pathSimu+"2021-03-31_unmasked.fa")
	fileO=open(pathSimu+"2021-03-31_unmasked_differences.txt","w")
	line=fileI.readline()
	nSeqs=0
	while line!="" and line!="\n":
		nSeqs+=1
		seq=""
		name=line.replace(">","").replace("\n","")
		fileO.write(line)
		line=fileI.readline()
		while line!="" and line!="\n" and line[0]!=">":
			seq+=line.replace("\n","")
			line=fileI.readline()
		if len(seq)!=lRef:
			print("Seq "+name+" has length "+str(len(seq))+" while ref is "+str(lRef))
			exit()
		# state 0=ref; 1=N; 2=-; 
		state=0
		seqList=[]
		length=0
		for i in range(lRef):
			if state==1:
				if seq[i]=="n":
					length+=1
				else:
					seqList.append(("n",i+1-length,length))
					length=0
					state=0
			elif state==2:
				if seq[i]=="-":
					length+=1
				else:
					seqList.append(("-",i+1-length,length))
					length=0
					state=0
			elif seq[i]=="n" and state!=1:
				length=1
				state=1
			elif seq[i]=="-" and state!=2:
				length=1
				state=2
			if seq[i]!=ref[i] and seq[i]!="-" and seq[i]!="n":
				seqList.append((seq[i],i+1))
		if state==1:
				seqList.append(("n",lRef+1-length,length))
		elif state==2:
				seqList.append(("-",lRef+1-length,length))
		for m in seqList:
			if len(m)==2:
				fileO.write(m[0]+"\t"+str(m[1])+"\n")
			else:
				fileO.write(m[0]+"\t"+str(m[1])+"\t"+str(m[2])+"\n")
		if (nSeqs%1000)==0:
			print(nSeqs)
	fileI.close()
	fileO.close()

	time2 = time.time() - start
	print("Time to convert alignment file: "+str(time2))
	print(str(nSeqs)+" sequences converted.")
	#Time take to convert alignment file: 15593.673340320587
	#915508 sequences converted.



#Read tree
start = time.time()
phylo = Tree(pathSimu+"global.tree")
time2 = time.time() - start
print("Time to read newick tree: "+str(time2))
print(str(len(phylo))+" sequences in the tree.")




def readConciseAlignment(fileName):
	start = time.time()
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
	time2 = time.time() - start
	print("Time to read DNA reduced data file: "+str(time2))
	print(str(nSeqs)+" sequences in file.")
	return data


if reduceDiffFile:
	#read sequence data from file
	data=readConciseAlignment(pathSimu+"2021-03-31_unmasked_differences.txt")

	dataReduced={}
	found=[0]
	nFound=[0]
	def traverseTree(node,found,nFound):
		if len(node.children)<1:
			if node.name in data:
				dataReduced[node.name]=data[node.name]
				found[0]+=1
			else:
				nFound[0]+=1

		else:
			for c in node.children:
				traverseTree(c,found,nFound)
	traverseTree(phylo,found,nFound)
	print(found)
	print(nFound)
	
	fileO=open(pathSimu+"2021-03-31_unmasked_differences_reduced.txt","w")
	for k in dataReduced.keys():
		fileO.write(">"+k+"\n")
		for m in dataReduced[k]:
			if len(m)==2:
				fileO.write(m[0]+"\t"+str(m[1])+"\n")
			else:
				fileO.write(m[0]+"\t"+str(m[1])+"\t"+str(m[2])+"\n")
	fileO.close()
	data=dataReduced
	#540520 sequences in the tree found
	#971 sequences not found
else:
	#read sequence data from file
	data=readConciseAlignment(pathSimu+"2021-03-31_unmasked_differences_reduced.txt")


samples=data.keys()
print(str(len(samples))+" sequences in the concise DNA data file")


#exit()
range4=range(4)

def simplfy(vec,refA):
	maxP=0.0
	maxI=0
	for i in range4:
		if vec[i]>maxP:
			maxP=vec[i]
			maxI=i
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




#preliminary nuc mutation rate matrix, from De Maio et al 2021
mutMatrix=[[0.0,0.039,0.310,0.123],[0.140,0.0,0.022,3.028],[0.747,0.113,0.0,2.953],[0.056,0.261,0.036,0.0]]
mutMatrix[0][0]=-(mutMatrix[0][1]+mutMatrix[0][2]+mutMatrix[0][3])
mutMatrix[1][1]=-(mutMatrix[1][0]+mutMatrix[1][2]+mutMatrix[1][3])
mutMatrix[2][2]=-(mutMatrix[2][0]+mutMatrix[2][1]+mutMatrix[2][3])
mutMatrix[3][3]=-(mutMatrix[3][0]+mutMatrix[3][1]+mutMatrix[3][2])
thresholdProb2=thresholdProb*thresholdProb
#Calculate likelihood?
def pruningFast(node,mutMatrix,data,ref, cumulativeBases):
	children=node.children
	bLen=node.dist
	#cumulative contribution of the subtree to the tree likelihood (excluding probs in the vector probVect)
	cumulPartLk=0.0
	#tree tip
	if len(children)==0:
		if not (node.name in data):
			probVect=[["N",1,lRef,bLen]]
			if verbose:
				print("\n"+"Terminal node "+node.name)
				print(probVect)
			return cumulPartLk, probVect
		diffs=data[node.name]
		pos=1
		#vector of states and probabilities to be passed on up the tree
		probVect=[]
		for m in diffs:
			currPos=m[1]
			if currPos>pos:
				#region where the node with branch length bLen is identical to the ref.
				probVect.append(["R",pos,currPos-pos,bLen])
				pos=currPos
			if m[0]=="n" or m[0]=="-":
				length=m[2]
				#region with no info, store first position and length.
				probVect.append(["N",currPos,length,bLen])
				pos=currPos+length
			elif m[0] in allelesLow:
				#position at which node allele is sure but is different from the reference.
				probVect.append([m[0].upper(),currPos,1,bLen])
				pos=currPos+1
			else:
				# non-"n" ambiguity character; for now interpret this as ambiguity instead of as a polymorphism.
				if onlyNambiguities:
					# if user asks to, to make things easier, interpret any ambiguity as an "n".
					probVect.append(["N",currPos,1,bLen])
				else:
					#otherwise, store as an "other" scenario, where each nucleotide has its own partial likelihood.
					probVect.append(["O",currPos,1,bLen,ambiguities[m[0]]])
				pos=currPos+1
		if pos<=lRef:
			probVect.append(["R",pos,lRef+1-pos,bLen])
		node.cumulPartLk=cumulPartLk
		node.probVect=probVect
		if verbose:
			print("\n"+"Terminal node "+node.name+" "+str(node.dist))
			print(diffs)
			print(probVect)
		#print("\n"+"Terminal node "+node.name+" "+str(node.dist))
		#print(probVect)
		return cumulPartLk, probVect
	
	#internal node
	else:
		cumulPartLk1, probVect1= pruningFast(children[0],mutMatrix,data,ref, cumulativeBases)
		cumulPartLk+=cumulPartLk1
		if verbose:
			print("\n"+"Internal node "+node.name)
			print("First child "+children[0].name)
			print(probVect1)
			print(cumulPartLk1)
		for c in range(len(children)-1):
			newChild=children[c+1]
			cumulPartLk2, probVect2= pruningFast(newChild,mutMatrix,data,ref, cumulativeBases)
			cumulPartLk+=cumulPartLk2
			if verbose:
				print("Next child "+newChild.name)
				print(probVect2)
				print(cumulPartLk2)

			#now traverse the two genome vectors and create a new probVect while also updating the cumulPartLk
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
			totLen=entry1[3]+entry2[3]
			#print("totLen")
			#print(totLen)
			#print(entry1[3])
			#print(entry2[3])
			#print("")
			while True:
				totLen=entry1[3]+entry2[3]
				if verbose:
					print("pos "+str(pos)+" length "+str(length)+" end "+str(end)+" indexEntry1 "+str(indexEntry1)+" indexEntry2 "+str(indexEntry2))
				if length<0:
					print("Issue with length")
					exit()
				if entry1[0]=="N":
					probVect.append(list(entry2))
					if entry2[0]=="N" or entry2[0]=="R": 
						probVect[-1][2]=length
						probVect[-1][1]=pos
				elif entry2[0]=="N":
					probVect.append(list(entry1))
					if entry1[0]=="N" or entry1[0]=="R": 
						probVect[-1][2]=length
						probVect[-1][1]=pos
				elif entry1[0]=="R":
					if entry2[0]=="R":
						probVect.append(["R",pos,length,0.0])
						#Now update partial likelihood
						if useLogs:
							for i in range4:
								cumulPartLk+=math.log(1.0+mutMatrix[i][i]*totLen)*(cumulativeBases[end][i]-cumulativeBases[pos-1][i])
						else:
							for i in range4:
								#print(totLen)
								#print(mutMatrix[i][i])
								#print(cumulativeBases[end][i]-cumulativeBases[pos-1][i])
								cumulPartLk+=mutMatrix[i][i]*totLen*(cumulativeBases[end][i]-cumulativeBases[pos-1][i])
								#print(cumulPartLk)
								#print("")
					elif entry2[0]=="O":
						i1=allelesLow[ref[pos-1]]
						newVec=np.ones(4)
						if entry2[3]+entry1[3]<thresholdProb2:
							if verbose:
								print("Warning: there might be a substitution on a 0-length cherry. Pretending branch lengths are not exactly zero.")
								print(entry1)
								print(entry2)
								print(i1)
							entry2[3]=thresholdProb
							entry1[3]=thresholdProb
						for i in range4:
							if i==i1:
								newVec[i]=1.0+mutMatrix[i][i]*entry1[3]
							else:
								newVec[i]=mutMatrix[i][i1]*entry1[3]
							tot=0.0
							for j in range4:
								if i==j:
									tot+=(1.0+mutMatrix[i][i]*entry2[3])*entry2[4][j]
								else:
									tot+=mutMatrix[i][j]*entry2[3]*entry2[4][j]
							newVec[i]*=tot
						state, newVec, maxP =simplfy(newVec,ref[pos-1])
						probVect.append([state,pos,1,0.0,newVec])
						cumulPartLk+=math.log(maxP)
						if np.sum(newVec)<thresholdProb2:
								if verbose:
									print("Issue merging R and O")
									print(entry1)
									print(entry2)
									print(probVect[-1])
								#exit()
					else:
						i2=alleles[entry2[0]]
						i1=allelesLow[ref[pos-1]]
						newVec=np.ones(4)
						if entry2[3]+entry1[3]<thresholdProb2:
							if verbose:
								print("Warning: there is a substitution on a 0-length cherry. Pretending branch lengths are not exactly zero.")
								print(entry1)
								print(entry2)
							entry2[3]=thresholdProb
							entry1[3]=thresholdProb
							#exit()
						for i in range4:
							if i==i2:
								newVec[i]=1.0+mutMatrix[i][i]*entry2[3]
							else:
								newVec[i]=mutMatrix[i][i2]*entry2[3]
							if i==i1:
								newVec[i]*=1.0+mutMatrix[i][i]*entry1[3]
							else:
								newVec[i]*=mutMatrix[i][i1]*entry1[3]
						if entry2[3]<thresholdProb2 or entry1[3]<thresholdProb2 :
							state, newVec, maxP =simplfy(newVec,ref[pos-1])
							probVect.append([state,pos,1,0.0,newVec])
						else:
							probVect.append(["O",pos,1,0.0,newVec])
						if np.sum(newVec)<thresholdProb2:
								if verbose:
									print("Issue merging R and non-R")
									print(entry1)
									print(entry2)
									print(probVect[-1])
								#exit()
				elif entry1[0]=="O":
					newVec=np.ones(4)
					if entry2[3]+entry1[3]<thresholdProb2:
							if verbose:
								print("Warning: there might be a substitution on this branch which has 0 length - Pretending branch lengths are not exactly zero.")
								print(entry1)
								print(entry2)
							entry2[3]=thresholdProb
							entry1[3]=thresholdProb
					if entry2[0]=="O":
						#print("O")
						for i in range4:
							tot1=0.0
							tot2=0.0
							for j in range4:
								if i==j:
									tot1+=(1.0+mutMatrix[i][i]*entry1[3])*entry1[4][j]
									tot2+=(1.0+mutMatrix[i][i]*entry2[3])*entry2[4][j]
								else:
									tot1+=mutMatrix[i][j]*entry1[3]*entry1[4][j]
									tot2+=mutMatrix[i][j]*entry2[3]*entry2[4][j]
							newVec[i]*=tot1*tot2
					else:
						for i in range4:
							tot1=0.0
							for j in range4:
								if i==j:
									tot1+=(1.0+mutMatrix[i][i]*entry1[3])*entry1[4][j]
								else:
									tot1+=mutMatrix[i][j]*entry1[3]*entry1[4][j]
							newVec[i]=tot1
						if verbose:
							print("Half-way newVec")
							print(entry1)
							print(newVec)
						if entry2[0]=="R":
							i2=allelesLow[ref[pos-1]]
							for i in range4:
								if i==i2:
									newVec[i]*=1.0+mutMatrix[i][i]*entry2[3]
								else:
									newVec[i]*=mutMatrix[i][i2]*entry2[3]
						else:
							i2=alleles[entry2[0]]
							for i in range4:
								if i==i2:
									newVec[i]*=1.0+mutMatrix[i][i]*entry2[3]
								else:
									newVec[i]*=mutMatrix[i][i2]*entry2[3]
					state, newVec, maxP =simplfy(newVec,ref[pos-1])
					probVect.append([state,pos,1,0.0,newVec])
					if verbose:
						print(newVec)
						print(state)
						print(maxP)
						if np.sum(newVec)<thresholdProb2:
								print("Issue merging O and something else")
								print(entry1)
								print(entry2)
								print(probVect[-1])
								#exit()
					cumulPartLk+=math.log(maxP)
				#entry1 is a non-ref nuc
				else:
					if entry2[0]==entry1[0]:
						i2=alleles[entry2[0]]
						probVect.append([entry1[0],pos,1,0.0])
						#Now update partial likelihood
						if useLogs:
							cumulPartLk+=math.log(1.0+mutMatrix[i2][i2]*totLen)
						else:
							cumulPartLk+=mutMatrix[i2][i2]*totLen
					else:
						if entry2[3]+entry1[3]<thresholdProb2:
							if verbose:
								print("Warning: there might be a substitution on this branch which has 0 length - Pretending branch lengths are not exactly zero.")
								print(entry1)
								print(entry2)
							entry2[3]=thresholdProb
							entry1[3]=thresholdProb
						i1=alleles[entry1[0]]
						newVec=np.ones(4)
						for i in range4:
							if i==i1:
								newVec[i]=1.0+mutMatrix[i][i]*entry1[3]
							else:
								newVec[i]=mutMatrix[i][i1]*entry1[3]
						if entry2[0]=="O":
							for i in range4:
								tot=0.0
								for j in range4:
									if i==j:
										tot+=(1.0+mutMatrix[i][i]*entry2[3])*entry2[4][j]
									else:
										tot+=mutMatrix[i][j]*entry2[3]*entry2[4][j]
								newVec[i]*=tot
							state, newVec, maxP =simplfy(newVec,ref[pos-1])
							probVect.append([state,pos,1,0.0,newVec])
							cumulPartLk+=math.log(maxP)
							if np.sum(newVec)<thresholdProb2:
								if verbose:
									print("Issue merging non-R and O")
									print(entry1)
									print(entry2)
									print(probVect[-1])
								#exit()
						else:
							if entry2[0]=="R":
								i2=allelesLow[ref[pos-1]]
							else:
								i2=alleles[entry2[0]]
							for i in range4:
								if i==i2:
									newVec[i]*=1.0+mutMatrix[i][i]*entry2[3]
								else:
									newVec[i]*=mutMatrix[i][i2]*entry2[3]
							probVect.append(["O",pos,1,0.0,newVec])
							if np.sum(newVec)<thresholdProb2:
								if verbose:
									print("Issue merging non-R and (R or non-R)")
									print(entry1)
									print(entry2)
									print(probVect[-1])
								#exit()

				#update pos, end, etc
				if verbose:
					print(pos)
					print(length)
				pos+=length
				if verbose:
					print(entry1)
					print(entry2)
					print(probVect1)
					print(probVect2)
					print("New values:")
					print(pos)
				if pos>lRef:
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
				if verbose:
					print(entry1)
					print(pos1)
					print(end1)
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
				if verbose:
					print(entry2)
					print(pos2)
					print(end2)
				end=end1
				if end2<end:
					end=end2
				length=end+1-pos
				if verbose:
					print(end)
					print("New length:")
					print(length)

			if verbose:
				print("Updated genome vector for internal node "+node.name)
				print(probVect)
			probVect1=probVect

			#check if the final  probVect can be simplified by merging consecutive entries
			probVect1 =shorten(probVect1)
			if verbose:
				print("Shortened genome vector for internal node "+node.name)
				print(probVect1)

		#update probVect to include the impact of its own branch
		for entry in probVect1:
			entry[3]+=node.dist
		
		node.cumulPartLk=cumulPartLk
		node.probVect=probVect1
		if verbose:
			print("\n"+"Internal node "+node.name+", final genome vector and likelihood contribution:")
			print(probVect1)
			print(cumulPartLk)
		#print("\n"+"Internal node with top branch length "+str(node.dist))
		#print(probVect1)
		return cumulPartLk, probVect1

#once reached the root, finalize the likelihood
def finalLK(probVect,ref,rootFreqsLog,rootFreqs):
	logLK=0.0
	for entry in probVect:
		if entry[3]<thresholdProb:
			#approximate assuming the distance is 0.0
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
		else:
			if entry[0]=="R":
				pos=entry[1]
				end=entry[2]+pos-1
				for j in range(4):
					prob=0.0
					for i in range(4):
						if i==j:
							prob+=rootFreqs[i]*(1.0+mutMatrix[i][j]*entry[3])
						else:
							prob+=rootFreqs[i]*entry[3]*mutMatrix[i][j]
					logLK+=math.log(prob)*(cumulativeBases[end][j]-cumulativeBases[pos-1][j])
			elif entry[0]=="N":
				pass
			elif entry[0]=="O":
				tot=0.0
				for i in range(4):
					for j in range(4):
						if i==j:
							tot+=rootFreqs[i]*entry[4][j]*(1.0+mutMatrix[i][j]*entry[3])
						else:
							tot+=rootFreqs[i]*entry[4][j]*mutMatrix[i][j]*entry[3]
				logLK+=math.log(tot)
			else:
				j=alleles[entry[0]]
				prob=0.0
				for i in range(4):
						if i==j:
							prob+=rootFreqs[i]*(1.0+mutMatrix[i][j]*entry[3])
						else:
							prob+=rootFreqs[i]*entry[3]*mutMatrix[i][j]
				logLK+=math.log(prob)
	return logLK



if runAllData:		
	#calculate likelihood of phylogeny
	start = time.time()
	cumulPartLk, probVect = pruningFast(phylo,mutMatrix,data,ref, cumulativeBases)
	time2 = time.time() - start
	print("Fast likelihood finished")
	print(cumulPartLk)
	print(probVect)
	#-8526086.483727645
	#[['R', 1, 14, 0.0], ['T', 15, 1, 0.0, array([0., 0., 0., 1.])], ['T', 16, 1, 0.0], ['R', 17, 8765, 0.0], ['T', 8782, 1, 0.0], ['R', 8783, 19361, 0.0], ['C', 28144, 1, 0.0], ['R', 28145, 1709, 0.0], ['G', 29854, 1, 0.0, array([0., 0., 1., 0.])], ['R', 29855, 1, 0.0], ['A', 29856, 1, 0.0, array([1., 0., 0., 0.])], ['G', 29857, 1, 0.0, array([0., 0., 1., 0.])], ['R', 29858, 3, 0.0], ['C', 29861, 1, 0.0, array([0., 1., 0., 0.])], ['R', 29862, 30, 0.0]]
	#now use root frequencies to finalize the likelihood.
	cumulPartLk+=finalLK(probVect,ref,rootFreqsLog,rootFreqs)
	print("Final likelihood after finalization:")
	print(cumulPartLk)
	print("Time taken "+str(time2))








#function to run likelihood calculation in a subsampled tree to compare it to full likelihood calculation in other methods.
def countLeaves(phylo,newSamples):
		if len(phylo.children)==0:
			if phylo.name in newSamples:
				#print(phylo.name)
				#print(samples)
				#print(phylo.name in samples)
				#exit()
				phylo.descendants=1
			else:
				phylo.descendants=0
		else:
			tot=0
			for c in phylo.children:
				countLeaves(c,newSamples)
				tot+=c.descendants
			phylo.descendants=tot

def extractSubtree(phylo,t1,minBlen=0.00001):
		if phylo.descendants>0:
			nChildren=0
			if len(phylo.children)==0:
				t1.name=phylo.name
				t1.dist+=phylo.dist
			newChildren=[]
			for c in phylo.children:
				if c.descendants>0:
					nChildren+=1
					newChildren.append(c)
			if nChildren==1:
				t1.dist+=phylo.dist
				extractSubtree(newChildren[0],t1,minBlen=minBlen)
			if nChildren>1:
				for i in range(nChildren-1):
					t1.add_child(dist=minBlen)
					extractSubtree(newChildren[i],t1.children[0],minBlen=minBlen)
					if i!=(nChildren-2):
						t1.add_child(dist=minBlen)
						t1=t1.children[1]
				t1.add_child(dist=minBlen)
				extractSubtree(newChildren[-1],t1.children[1],minBlen=minBlen)
						
if subsample>0:
	#sample "subsample" many leaves of the tree
	#samples=data.keys()
	for repeat in range(nRepeats):
		newSamples=random.sample(samples,subsample)
		#print(newSamples)

		#extract subtree containing these samples and replace 0 lengths with non-zero lengths and multifurcations with bifurcations.
		countLeaves(phylo,newSamples)
		t1 = Tree()	
		extractSubtree(phylo,t1)
		print(str(len(t1))+" sequences in the new tree.")
		#exit()
		#print(t1.write())
		t1.write(outfile=pathSimu+"2021-03-31_subsample"+str(subsample)+"_repeat"+str(repeat)+".tree")

		#possibly remove parts of the genome
		LKlimit=10.3
		genomeEnd=29970
		for s in newSamples:
			for m in range(len(data[s])):
				if data[s][m][1]>genomeEnd:
					data[s]=data[s][:m]
					break
		genomeStart=-10
		for s in newSamples:
			for m in range(len(data[s])):
				if data[s][-(m+1)][1]<genomeStart:
					if m==0:
						data[s]=[]
					else:
						data[s]=data[s][-(m):]
					break

		#write data files for subsample
		#concise file
		fileO=open(pathSimu+"2021-03-31_subsample"+str(subsample)+"_repeat"+str(repeat)+".txt","w")
		for s in newSamples:
			#data[s]=[["t",1],["g",2],["m",4],["-",5,29880],["r",29891]]
			#data[s]=[]
			#if s==newSamples[1]:
			#	data[s]=[["c",1],["y",3],["-",5,29884],["y",29891]]
			#elif s==newSamples[2]:
			#	data[s]=[["g",1],["g",2],["g",3],["-",5,29882],["m",29891]]
				#data[s]=[["-",1,29890]]
			fileO.write(">"+s+"\n")
			for m in data[s]:
				if len(m)==2:
						fileO.write(m[0]+"\t"+str(m[1])+"\n")
				else:
						fileO.write(m[0]+"\t"+str(m[1])+"\t"+str(m[2])+"\n")

		fileO.close()
		#alignment file
		fileO=open(pathSimu+"2021-03-31_subsample"+str(subsample)+"_repeat"+str(repeat)+".phylip","w")
		fileO.write(str(subsample)+"\t"+str(lRef)+"\n")
		for s in newSamples:
			fileO.write(s+" ")
			refList=list(ref)
			for m in data[s]:
				if len(m)==2:
					refList[m[1]-1]=m[0]
				else:
					for i in range(m[2]):
						refList[m[1]+i-1]=m[0]
			#print(refList[:100])
			#print(data[s])
			fileO.write("".join(refList)+"\n")
		fileO.close()

		#Run fast likelihood
		print("Starting new Felsenstein on subsample of "+str(subsample)+" sequences, repeat "+str(repeat))
		start = time.time()
		dataSample=readConciseAlignment(pathSimu+"2021-03-31_subsample"+str(subsample)+"_repeat"+str(repeat)+".txt")
		#print(dataSample)
		ref, cumulativeBases, rootFreqs, rootFreqsLog = collectReference(pathSimu+"EPI_ISL_402124_lowercase.fasta")
		lRef=len(ref)
		phylo2 = Tree(pathSimu+"2021-03-31_subsample"+str(subsample)+"_repeat"+str(repeat)+".tree")
		mutMatrix=[[0.0,0.04,0.30,0.10],[0.04,0.0,0.02,1.0],[0.30,0.02,0.0,1.0],[0.10,1.0,1.0,0.0]]
		tot=0.0
		for i in range(4):
			for j in range(4):
				tot+=mutMatrix[i][j]
		tot=tot/4
		for i in range(4):
			for j in range(4):
				mutMatrix[i][j]=mutMatrix[i][j]/tot
		for i in range(4):
			tot=0.0
			for j in range(4):
				tot+=mutMatrix[i][j]
			mutMatrix[i][i]=-tot
		print("Substitution rate matrix:")
		print(mutMatrix)
		rootFreqs=[0.25,0.25,0.25,0.25]
		rootFreqsLog=[math.log(0.25),math.log(0.25),math.log(0.25),math.log(0.25)]
		cumulPartLk, probVect = pruningFast(phylo2,mutMatrix,dataSample,ref, cumulativeBases)
		print("Fast likelihood finished")
		print(cumulPartLk)
		print(probVect)
		cumulPartLk+=finalLK(probVect,ref,rootFreqsLog,rootFreqs)
		time2 = time.time() - start
		print("Final log-likelihood after finalization:")
		print(cumulPartLk)
		print("Time taken "+str(time2))

		#Run PhyML
		start = time.time()
		os.system("/Applications/phyml-master/src/phyml --no_memory_check --leave_duplicates --run_id subsample"+str(subsample)+"_repeat"+str(repeat)+"_phyml --print_trace -o n -c 1 -b 0 -m 0.04,0.3,0.1,0.02,1.0,1.0 -f 0.25,0.25,0.25,0.25 -i "+pathSimu+"2021-03-31_subsample"+str(subsample)+"_repeat"+str(repeat)+".phylip -u "+pathSimu+"2021-03-31_subsample"+str(subsample)+"_repeat"+str(repeat)+".tree")
		time2 = time.time() - start
		print("Time taken by PhyML "+str(time2))
		filePhy=open(pathSimu+"2021-03-31_subsample"+str(subsample)+"_repeat"+str(repeat)+".phylip_phyml_stats_subsample"+str(subsample)+"_repeat"+str(repeat)+"_phyml.txt")
		line=filePhy.readline()
		linelist=line.split()
		while len(linelist)<2 or linelist[1]!="Log-likelihood:":
			line=filePhy.readline()
			linelist=line.split()
		phymlLK=float(linelist[2])
		if abs(phymlLK-cumulPartLk)>LKlimit:
			print("Found discrepancy in likelihood")
			exit()



			
					
					







#estimate susbstitution rates?
#optimize branch lengths?
#Optimize root position?
#Calculate pairwise distances?
#reduce tree by removing sequences that are basically identical to others?
#try to improve topology?

exit()






















