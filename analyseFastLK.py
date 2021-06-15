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

#©EMBL-European Bioinformatics Institues, 2021

#Run a bunch of analyses regarding efficiently calculating phylogenetic likelihoods.
#Also includes experimental code for estimating an initial tree (not finished yet).
#Example run command line: python3 /Users/demaio/Desktop/GISAID-hCoV-19-phylogeny-2021-03-12/createDiffsFile.py --path /Users/demaio/Desktop/GISAID-hCoV-19-phylogeny-2021-03-12/

parser = argparse.ArgumentParser(description='Run estimation of mutation rates and synonymous site selection from SARS-CoV-2 data, and generate plots.')
parser.add_argument('--path',default="", help='path where to find files and plot results.')
parser.add_argument("--createDiffFile", help="create a concise file replacing a MSA.", action="store_true")
parser.add_argument("--reduceDiffFile", help="Reduce the size of a concise file by only including samples in the tree.", action="store_true")
parser.add_argument("--onlyNambiguities", help="Treat all ambiguities as N (total missing information).", action="store_true")
parser.add_argument("--useLogs", help="Calculate logarithms of non-mutation probabilities, otherwise approximate them.", action="store_true")
parser.add_argument("--thresholdProb",help="relative probability threshold used to ignore possible states with very low probabilities.",  type=float, default=0.0000001)
parser.add_argument("--verbose", help="Print to screen a lot of stuff.", action="store_true")
parser.add_argument("--runAllData", help="Run the fast likelihood on the huge tree (takes about 2 minutes).", action="store_true")
parser.add_argument("--subsample",help="Sample this number of samples from the original dataset and run likelihood comparison with PhyML.",  type=int, default=0)
parser.add_argument("--nRepeats",help="When subsampling the tree, this is the number of repeats to be done.",  type=int, default=1)
parser.add_argument("--seed",help="Seed for random number generator.",  type=int, default=1)
parser.add_argument("--debug", help="Run continuously on random samples until it finds a discrepancy in log-LK with PAML, at which point investigate which part of the genome is the cause.", action="store_true")
parser.add_argument("--newPhylogeny", help="Estimate a brand new phylogeny from the data.", action="store_true")
parser.add_argument("--runPhyml", help="Compare likelihood and running time with PhyML when subsampling.", action="store_true")
parser.add_argument("--runIQtree", help="Compare likelihood and running time with IQ-TREE when subsampling.", action="store_true")
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
debug=args.debug
newPhylogeny=args.newPhylogeny
runPhyml=args.runPhyml
runIQtree=args.runIQtree

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

#print(ref[21985:21989])
#print(ref[21960:22010])

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


if reduceDiffFile or runAllData or subsample>0:
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

checkD142G=False
if checkD142G:
	equalRef=0
	diffRef=0
	tot21991=0
	tot21988=0
	tot21990=0
	tot21993=0
	tot21987a=0
	tot21987t=0
	tot21985t=0
	tot21986a=0
	for s in samples:
		found=False
		found21991=False
		found21988=False
		found21990=False
		found21993=False
		found21987a=False
		found21986a=False
		found21987t=False
		found21985t=False
		for d in data[s]:
			if (d[0]=="-" or d[0]=="n") and d[1]==21991:
				found21991=True
				tot21991+=1
			elif (d[0]=="-" or d[0]=="n") and d[1]==21988:
				found21988=True
				tot21988+=1
			elif (d[0]=="-" or d[0]=="n") and d[1]==21990:
				found21990=True
				tot21990+=1
			elif (d[0]=="t") and d[1]==21993:
				found21993=True
				tot21993+=1
			elif d[1]==21987 and d[0]=="a":
				found21987a=True
				tot21987a+=1
			elif d[1]==21986 and d[0]=="a":
				found21986a=True
				tot21986a+=1
			elif d[1]==21987 and d[0]=="t":
				found21987t=True
				tot21987t+=1
			elif d[1]==21985 and d[0]=="t":
				found21985t=True
				tot21985t+=1
			elif d[1]>21983 and d[1]<21993 and d[0]!="n" and d[0]!="k" and d[0]!="r" and d[0]!="y":
				print(s)
				print(d)
				print(found21991)
				print(found21988)
				print(found21990)
				print(found21993)
				print(found21987a)
				print(found21987t)
				print(found21985t)
				print(found21986a)
				print("")
				found=True
		if found:
			diffRef+=1
		else:
			equalRef+=1
	print("Different from ref: "+str(diffRef)+"  with 21991 del: "+str(tot21991)+"  with 21988 del: "+str(tot21988)+"  with 21990 del: "+str(tot21990)+"  t at 21993: "+str(tot21993)+"  a at 21987 "+str(tot21987a)+"  t at 21987 "+str(tot21987t)+"   t at 21985 "+str(tot21985t)+"   a at 21986 "+str(tot21986a)+"  rest: "+str(equalRef))
	exit()




#exit()
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
mutMatrix[0][0]=-(mutMatrix[0][1]+mutMatrix[0][2]+mutMatrix[0][3])
mutMatrix[1][1]=-(mutMatrix[1][0]+mutMatrix[1][2]+mutMatrix[1][3])
mutMatrix[2][2]=-(mutMatrix[2][0]+mutMatrix[2][1]+mutMatrix[2][3])
mutMatrix[3][3]=-(mutMatrix[3][0]+mutMatrix[3][1]+mutMatrix[3][2])
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
		# probVect=[]
		# for m in diffs:
		# 	currPos=m[1]
		# 	if currPos>pos:
		# 		#region where the node with branch length bLen is identical to the ref.
		# 		probVect.append(["R",pos,currPos-pos,bLen])
		# 		pos=currPos
		# 	if m[0]=="n" or m[0]=="-":
		# 		length=m[2]
		# 		#region with no info, store first position and length.
		# 		probVect.append(["N",currPos,length,bLen])
		# 		pos=currPos+length
		# 	elif m[0] in allelesLow:
		# 		#position at which node allele is sure but is different from the reference.
		# 		probVect.append([m[0].upper(),currPos,1,bLen])
		# 		pos=currPos+1
		# 	else:
		# 		# non-"n" ambiguity character; for now interpret this as ambiguity instead of as a polymorphism.
		# 		if onlyNambiguities:
		# 			# if user asks to, to make things easier, interpret any ambiguity as an "n".
		# 			probVect.append(["N",currPos,1,bLen])
		# 		else:
		# 			#otherwise, store as an "other" scenario, where each nucleotide has its own partial likelihood.
		# 			probVect.append(["O",currPos,1,bLen,ambiguities[m[0]]])
		# 		pos=currPos+1
		# if pos<=lRef:
		# 	probVect.append(["R",pos,lRef+1-pos,bLen])
		probVect=probVectTerminalNode(diffs,bLen)
		node.cumulPartLk=cumulPartLk
		node.probVect=probVect
		if verbose:
			print("\n"+"Terminal node "+node.name+" "+str(node.dist))
			print(diffs)
			print(probVect)
		if debug:
			print("\n"+"Terminal node "+node.name+" "+str(node.dist))
			print(probVect)
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
			#totLen=entry1[3]+entry2[3]
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
						print("Simplified O")
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
					#print(entry1)
					#print(entry2)
					#print(probVect1)
					#print(probVect2)
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
		if debug:
			print("\n"+"Internal node with top branch length "+str(node.dist))
			print(probVect1)
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





if newPhylogeny:
	#read metadata to find sample times
	#they are not exact times, just approximations to sort the samples based on collection times.
	file=open(pathSimu+"metadata.csv")
	line=file.readline()
	line=file.readline()
	dates={}
	while line!="" and line!="\n":
		linelist=line.split(",")
		name=linelist[0]
		dateList=(linelist[-1].replace("\n","")).split("-")
		dates[name]=(int(dateList[0])-2019)*366
		if len(dateList)>1:
			dates[name]+=31*int(dateList[1])
		else:
			dates[name]+=31*12
		if len(dateList)>2:
			dates[name]+=int(dateList[2])
		else:
			dates[name]+=31
		line=file.readline()
	file.close()






#sort sequences according to their sampling date. Not used as one ends up placing highly ambiguous sequences first, which is not optimal.
def sortByDate(dates):
	#start = time.time()
	sampleDistances=[]
	for sample in samples:
		date=dates[sample]
		sampleDistances.append((date,sample))

	from operator import itemgetter
	print("Now doing sorting")
	sampleDistances.sort(key=itemgetter(0))
	return sampleDistances


#Same as distancesFromRef() but punishing more isolated N's and ambiguity codes.
def distancesFromRefPunishNs(data):
	#start = time.time()
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
		#print(sample)
		#print(diffs)
		#print(str(comparisons)+" "+str(diffNum)+"\n")
		#if comparisons==0:
		#	sampleDistances.append((lRef,sample))
		#else:
		sampleDistances.append((diffNum*1000+lRef-comparisons,sample))

	from operator import itemgetter
	print("Now doing sorting")
	sampleDistances.sort(key=itemgetter(0))
	#print(sampleDistances)
	#time2 = time.time() - start
	#print("Time to create a sorted list of distances from reference: "+str(time2))
	return sampleDistances


#Calculate distances of sequences from the reference, and return sorted list, with closest sequences first
def distancesFromRef(data):
	#start = time.time()
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
			elif m[0] in allelesLow:
				comparisons+=1
				diffNum+=1
				pos=currPos+1
			else:
				pos=currPos+1
		if pos<=lRef:
			comparisons+=lRef+1-pos
		#print(sample)
		#print(diffs)
		#print(str(comparisons)+" "+str(diffNum)+"\n")
		#if comparisons==0:
		#	sampleDistances.append((lRef,sample))
		#else:
		sampleDistances.append((diffNum*(lRef+1)+lRef-comparisons,sample))

	from operator import itemgetter
	print("Now doing sorting")
	sampleDistances.sort(key=itemgetter(0))
	#print(sampleDistances)
	#time2 = time.time() - start
	#print("Time to create a sorted list of distances from reference: "+str(time2))
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
			#print("totLen")
			#print(totLen)
			#print(entry1[3])
			#print(entry2[3])
			#print("")
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





#function to calculate likelihood of appending child to parent
#CHECKED!
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
	#implement parent-child likelihoood cost calculation!!!
	while True:
		#print(entry1)
		#print(entry2)
		#print("new cost updated:")
		#print(Lkcost)
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
			#print("LkcostOriginal updated "+str(LkcostOriginal))
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

totalMissedMinors=[0]
allowedFails=3
def findBestParent(t1,diffs,sample,bLen,bestLKdiff,bestNodeSoFar,failedPasses):
	if t1.is_leaf():
		#print("Trying to place on leaf "+str(t1.dist))
		probVect=t1.probVect
		comparison=isMinorSequence(probVect,diffs)
		if comparison==1:
			t1.minorSequences.append(sample)
			#print(diffs)
			#print("appended to ")
			#print(probVect)
			return t1, 1.0
		elif comparison==2:
			print("Second sequence is more informative than first, this could be done better, but for now this fact is ignored and the two sequences are not merged.")
			print(t1.name)
			print(sample)
			print(probVect)
			print(diffs)
			totalMissedMinors[0]+=1
			#exit()
		#return bestNodeSoFar, bestLKdiff
	#if internal node, and if distant enough from its parent, try to append child here.
	#print("Trying to place "+str(t1.dist)+" "+str(len(t1.children)))
	#print(t1.probVect)
	if t1.dist>0.000000001:
		probVect=t1.probVectTot
		LKdiff=appendProb(probVect,diffs,bLen)
		#if LKdiff>-40:
		#	print(LKdiff)
		#	print(probVect)
		#	print(t1.dist)
		if verbose:
			print("Trying a new parent")
			print(probVect)
			print(diffs)
			print("Failed passes before this step: "+str(failedPasses)+" logLK score: "+str(LKdiff))
		foundBetter=False
		if LKdiff>bestLKdiff:
			if verbose:
				print("New best LK")
			bestLKdiff=LKdiff
			bestNodeSoFar=t1
			foundBetter=True
			failedPasses=0
		elif LKdiff<bestLKdiff-0.1:
			failedPasses+=1
			if failedPasses>allowedFails:
				#if verbose:
				#print("failed too many times, not passing to children")
				return bestNodeSoFar , bestLKdiff
	for c in t1.children:
		node, score = findBestParent(c,diffs,sample,bLen,bestLKdiff,bestNodeSoFar,failedPasses)
		if score>bestLKdiff:
			bestLKdiff=score
			bestNodeSoFar=node
	return bestNodeSoFar , bestLKdiff



#for the root, return a vector taking into account the root frequency of each state
#CHECKED
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
#CHECKED
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





#in order to insert a new child and calculate its cost, merge to vectors, one from above and one from below, at an intermediate position 
#CHECKED!
def mergeVectorsUpDown(probVectUp,bLenUp,probVectDown,bLenDown,mutMatrix):
			probVect1= probVectUp
			probVect2=probVectDown
			#now traverse the two genome vectors and create a new probVect
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
				#totLen=bLenUp+bLenDown+entry2[3]+entry1[3]
				if entry1[0]=="N":
					probVect.append(list(entry2))
					if entry2[0]=="N" or entry2[0]=="R": 
						probVect[-1][2]=length
						probVect[-1][1]=pos
					#probVect[-1][3]+=bLenDown
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
					#probVect[-1][3]+=bLenUp
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
							#print("O")
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
						#print("non-ref nuc")
						if entry2[0]==entry1[0] or (totLen1<thresholdProb2):
							probVect.append([entry1[0],pos,1,0.0,False])
						else:
							i1=alleles[entry1[0]]
							newVec=np.ones(4)
							if entry1[4]:
								#print("root stuff")
								#print(entry1[3])
								#print(entry1[5]+bLenUp)
								#print(mutMatrix[1][3])
								#print(mutMatrix[3][1])
								rootVec=np.ones(4)
								for i in range4:
									rootVec[i]=rootFreqs[i]
									if i==i1:
										rootVec[i]*=(1.0+mutMatrix[i][i1]*(entry1[3]))
									else:
										rootVec[i]*=mutMatrix[i][i1]*(entry1[3])
								#print(rootVec)
								for j in range4:
									tot=0.0
									for i in range4:
										if j==i:
											tot+=(1.0+mutMatrix[i][j]*(entry1[5]+bLenUp))*rootVec[i]
										else:
											tot+=mutMatrix[i][j]*(entry1[5]+bLenUp)*rootVec[i]
									newVec[j]=tot
								#print("newVec:")
								#print(newVec)
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
								#print("newVec again")
								#print(newVec)
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




#merge two child vectors at the root and calculate the logLk, and return vector and logLK.
#CHECKED!
def mergeVectorsRoot(probVect1,bLen1,probVect2,bLen2,mutMatrix):
			#now traverse the two genome vectors and create a new probVect while also updating the cumulPartLk
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




#merge two child vectors at an internal node.
#CHECKED!
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
				#totLen=entry1[3]+entry2[3]+bLen1+bLen2
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
					if not vSum>.000000000000000000000001:
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
				if not vSum>.000000000000000000000001:
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
	#print("Updating branch length following an inconsistency")
	cNode=node.children[childNum]
	if childNum==0:
		vectUp=node.probVectUpRight
	else:
		vectUp=node.probVectUpLeft
	bestLK=appendProb(vectUp,cNode.probVect,bLen)
	bestLen=bLen
	#print("First len "+str(bestLen)+" gives lk "+str(bestLK))
	while bestLen>0.1*bLen:
		newBLen=bestLen/2
		newLK=appendProb(vectUp,cNode.probVect,newBLen)
		#print("Attempting new len "+str(newBLen)+" gives lk "+str(newLK))
		if newLK>bestLK:
			bestLK=newLK
			bestLen=newBLen
		else:
			break
	if bestLen>0.7*bLen:
		while bestLen<10*bLen:
			newBLen=bestLen*2
			newLK=appendProb(vectUp,cNode.probVect,newBLen)
			#print("Attempting new len "+str(newBLen)+" gives lk "+str(newLK))
			if newLK>bestLK:
				bestLK=newLK
				bestLen=newBLen
			else:
				break
	cNode.dist=bestLen
	#print("Found best len "+str(bestLen)+" with lk "+str(bestLK))
	updatePartialsFromTop(cNode,vectUp,mutMatrix)
	updatePartialsFromBottom(node,cNode.probVect,childNum,cNode,mutMatrix)
	#print("Finished updating bLen")
	#exit()



#Update partial likelihood vectors in the tree after the addition of a new tip.
#This function traverses the tree downward.
def updatePartialsFromTop(node,probVectUp,mutMatrix):
	#print("\n\n Inside updatePartialsFromTop()")
	#print(probVectUp)
	#print(node.probVect)
	#print(node.dist)
	if len(node.children)==0:
		#print("\n updatePartialsFromTop Leaf node. Old tot:")
		#print(node.probVectTot)
		#print("Old patials")
		#print(node.probVect)
		#print("New vector from top")
		#print(probVectUp)
		#print(node.dist)
		newTot=mergeVectorsUpDown(probVectUp,node.dist,node.probVect,0.0,mutMatrix)
		if newTot==None:
			if node.up.children[0]==node:
				updateBLen(node.up,0,bLen)
			else:
				updateBLen(node.up,1,bLen)
			return
		node.probVectTot=newTot
		#print("New tot:")
		#print(node.probVectTot)
	else:
		child0Vect=node.children[0].probVect
		child1Vect=node.children[1].probVect
		dist0=node.children[0].dist
		dist1=node.children[1].dist
		#print("updatePartialsFromTop internal node. Old upRight:")
		#print(node.probVectUpRight)
		#print("probVectUp")
		#print(probVectUp)
		#print(node.dist)
		#print("child vector")
		#print(child1Vect)
		#print(dist1)
		newUpRight=mergeVectorsUpDown(probVectUp,node.dist,child1Vect,dist1,mutMatrix)
		#print("New upRight:")
		#print(newUpRight)
		#print("Updating upLeft")
		#print(probVectUp)
		#print(node.dist)
		#print(child0Vect)
		#print(dist0)
		#print("updatePartialsFromTop internal node. Old upLeft:")
		#print(node.probVectUpRight)
		newUpLeft=mergeVectorsUpDown(probVectUp,node.dist,child0Vect,dist0,mutMatrix)
		#print("New upLeft:")
		#print(newUpLeft)
		updatedTot=False
		#print("Inside updatePartialsFromTop")
		#print(node.probVectUpRight)
		#print(newUpRight)
		#print(node.probVectUpLeft)
		#print(newUpLeft)
		if areVectorsDifferent(node.probVectUpRight,newUpRight):
			#print("upRight is different, making the update.")
			#print(node.probVectUpRight)
			#print(newUpRight)
			newUpRight =shorten(newUpRight)
			#print("shortened")
			#print(newUpRight)
			node.probVectUpRight=newUpRight
			#print("updatePartialsFromTop internal node. Old probVectTot:")
			#print(node.probVectTot)
			#print("child vector")
			#print(child0Vect)
			newTot=mergeVectorsUpDown(newUpRight,0.0,child0Vect,dist0,mutMatrix)
			if newTot==None:
				updateBLen(node,0,bLen)
				return
			#print("New probVectTot:")
			#print(newTot)
			newTot =shorten(newTot)
			node.probVectTot=newTot
			updatedTot=True
			updatePartialsFromTop(node.children[0],newUpRight,mutMatrix)
		if areVectorsDifferent(node.probVectUpLeft,newUpLeft):
			#print("upLeft is different, making the update.")
			#print(node.probVectUpLeft)
			#print(newUpLeft)
			newUpLeft =shorten(newUpLeft)
			#print("shortened")
			#print(newUpLeft)
			node.probVectUpLeft=newUpLeft
			if not updatedTot:
				#print("updatePartialsFromTop internal node. Old probVectTot:")
				#print(node.probVectTot)
				#print("child vector")
				#print(child1Vect)
				newTot=mergeVectorsUpDown(newUpLeft,0.0,child1Vect,dist1,mutMatrix)
				if newTot==None:
					updateBLen(node,1,bLen)
					return
				newTot =shorten(newTot)
				#print("New probVectTot:")
				#print(newTot)
				node.probVectTot=newTot
			updatePartialsFromTop(node.children[1],newUpLeft,mutMatrix)
	#print("End of updatePartialsFromTop")
	#print(logLK)
	#print("\n")


#Update partial likelihood vectors in the tree after the addition of a new tip.
#This function traverses the tree upward.
def updatePartialsFromBottom(node,probVectDown,childNum,childNode,mutMatrix):
	#print("\n\n Inside updatePartialsFromBottom()")
	#print(node.probVect)
	#print(node.probVectUpRight)
	#print(node.probVectUpLeft)
	#print(node.probVectTot)
	#print(childNum)
	if childNum==-1:
		if node.children[0]==childNode:
			childNum=0
		else:
			childNum=1
	#print(childNum)
	otherChildNum=1-childNum
	childDist=node.children[childNum].dist
	otherChildDist=node.children[otherChildNum].dist
	otherChildVect=node.children[otherChildNum].probVect
	#print(childDist)
	#print(otherChildDist)
	#print(node.dist)
	#print(node.children[childNum].probVect)
	#print(probVectDown)
	#print(otherChildVect)
	newVect=mergeVectors(otherChildVect,otherChildDist,probVectDown,childDist,mutMatrix)
	#print("new children merged vect:")
	#print(newVect)
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
			#print("Inside updatePartialsFromBottom, child 0 old probVectTot:")
			#print(node.probVectTot)
			newTot=mergeVectorsUpDown(node.probVectUpRight,0.0,probVectDown,childDist,mutMatrix)
			if newTot==None:
				updateBLen(node,0,bLen)
				return
			#print("new probVectTot:")
			#print(newTot)
		else:
			#print("Inside updatePartialsFromBottom, child 1 old probVectTot:")
			#print(node.probVectTot)
			newTot=mergeVectorsUpDown(node.probVectUpLeft,0.0,probVectDown,childDist,mutMatrix)
			if newTot==None:
				updateBLen(node,1,bLen)
				return
			#print("new probVectTot:")
			#print(newTot)
		newTot=shorten(newTot)
		node.probVectTot=newTot
		if verbose:
			print("new tot vect")
			print(node.probVectTot)
		updatedTot=True
		if node.up != None:
			#print("updatePartialsFromBottom from the function itself")
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
			#print("Inside updatePartialsFromBottom, child 0 old upLeft:")
			#print(node.probVectUpLeft)
			newUpLeftVect=mergeVectorsUpDown(vectUp,node.dist,probVectDown,childDist,mutMatrix)
			#print("new upLeft:")
			#print(newUpLeftVect)
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
					#print("Inside updatePartialsFromBottom, child 0 old probVectTot:")
					#print(node.probVectTot)
					newTot=mergeVectorsUpDown(newUpLeftVect,0.0,otherChildVect,otherChildDist,mutMatrix)
					if newTot==None:
						updateBLen(node,1,bLen)
						return
					#print("new probVectTot:")
					#print(newTot)
					newTot=shorten(newTot)
					node.probVectTot=newTot
					if verbose:
						print(newTot)
				
				updatePartialsFromTop(node.children[1],newUpLeftVect,mutMatrix)
		else:
			#print("Inside updatePartialsFromBottom, child 1 old upRight:")
			#print(node.probVectUpRight)
			newUpRightVect=mergeVectorsUpDown(vectUp,node.dist,probVectDown,childDist,mutMatrix)
			#print("new upRight:")
			#print(newUpRightVect)
			if areVectorsDifferent(node.probVectUpRight,newUpRightVect):
				if verbose:
					print("Moving up the update from right")
				newUpRightVect =shorten(newUpRightVect)
				node.probVectUpRight=newUpRightVect
				if not updatedTot:
					#print("Inside updatePartialsFromBottom, child 1 old probVectTot:")
					#print(node.probVectTot)
					newTot=mergeVectorsUpDown(newUpRightVect,0.0,otherChildVect,otherChildDist,mutMatrix)
					if newTot==None:
						updateBLen(node,0,bLen)
						return
					#print("new probVectTot:")
					#print(newTot)
					newTot=shorten(newTot)
					node.probVectTot=newTot
				updatePartialsFromTop(node.children[0],newUpRightVect,mutMatrix)
	#case node is root
	else:
		if verbose:
			print("reached root while moving up")
		if childNum==0:
			#newUpLeftVect=finalLKwithVector(probVectDown,ref,rootFreqsLog,rootFreqs,bLen)
			newUpLeftVect=rootVector(probVectDown,rootFreqs,childDist)
			if areVectorsDifferent(node.probVectUpLeft,newUpLeftVect):
				newUpLeftVect =shorten(newUpLeftVect)
				node.probVectUpLeft=newUpLeftVect
				if not updatedTot:
					#print("Inside updatePartialsFromBottom, root, child 0 old probVectTot:")
					#print(node.probVectTot)
					newTot=mergeVectorsUpDown(newUpLeftVect,0.0,otherChildVect,otherChildDist,mutMatrix)
					#print("new probVectTot:")
					#print(newTot)
					newTot=shorten(newTot)
					node.probVectTot=newTot
				updatePartialsFromTop(node.children[1],newUpLeftVect,mutMatrix)

		else:
			#newUpRightVect=finalLKwithVector(probVectDown,ref,rootFreqsLog,rootFreqs,bLen)
			newUpRightVect=rootVector(probVectDown,rootFreqs,childDist)
			if areVectorsDifferent(node.probVectUpRight,newUpRightVect):
				newUpRightVect =shorten(newUpRightVect)
				node.probVectUpRight=newUpRightVect
				if not updatedTot:
					#print("Inside updatePartialsFromBottom, root, child 1 old probVectTot:")
					#print(node.probVectTot)
					newTot=mergeVectorsUpDown(newUpRightVect,0.0,otherChildVect,otherChildDist,mutMatrix)
					#print("new probVectTot:")
					#print(newTot)
					node.probVectTot=newTot
					newTot=shorten(newTot)
				updatePartialsFromTop(node.children[0],newUpRightVect,mutMatrix)
	if verbose:
		print("New upRight, UpLeft and tot:")
		print(node.probVectUpRight)
		print(node.probVectUpLeft)
		print(node.probVectTot)
	#print("End of updatePartialsFromBottom")
	#print(node.probVectUpRight)
	#print(node.probVectUpLeft)
	#print(node.probVectTot)
	#print("\n")
	return

#try all children branches (if meets a 0-length branch, pass the function to clid node) to find the best branch where to place (mid-length) a child node.
#This basically takes into account polytomies (here politomies are given by 0-length branches in a binary tree).
def findBestChildBranch(node,newPartials,sample,bLen):
	bestLKSoFar=float("-inf")
	bestNodeSoFar=None
	bestVectSoFar=None
	for c in range(len(node.children)):
		child=node.children[c]
		if child.dist>thresholdProb2:
			if c==0:
				#print("Inside findBestChildBranch 0")
				probVectChild=mergeVectorsUpDown(node.probVectUpRight,child.dist/2,child.probVect,child.dist/2,mutMatrix)
			else:
				#print("Inside findBestChildBranch 1")
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
# now explore exactly which position near node is best for placement, and which are the best branch lengths,
#then add the sample at that position of the tree, and update all the internal probability vectors.
def placeSampleOnTree(node,newPartials,sample,bLen,newChildLK):
	#place as sibling of either child
	#print("\n\n\n Running placeSampleOnTree(). node.tot:")
	#print(node.probVectTot)
	#print("New partials: ")
	#print(newPartials)
	#print("LK of initial placement as child: "+str(newChildLK))
	if len(node.children)>0:
		#print("Placing sample at internal node")
		childBestNode, childBestLK, childBestVect=findBestChildBranch(node,newPartials,sample,bLen)
		#print("Found best child node, with tot LK")
		# if childBestNode!=None:
		# 	print(childBestNode.probVectTot)
		# 	print("With placement lk: "+str(childBestLK)+" and vector: ")
		# 	print(childBestVect)
		# else:
		# 	print("No suitable children (with dist>0) found.")
	else:
		#print("Placing sample at tip.")
		childBestLK=float("-inf")
		childBestNode=None
		childBestVect=None
	#if node is root, try to place as sibling of the current root.
	#print(node.up)
	if node.is_root():
		probOldRoot = findProbRoot(node.probVect)
		#print("node is root. LK of old root: "+str(probOldRoot))
		probVectRoot,probRoot = mergeVectorsRoot(node.probVect,bLen,newPartials,bLen,mutMatrix)
		#print("If new root created, vector")
		#print(probVectRoot)
		if verbose:
			print("parent is root, trying to add new root; node.up, node.probvect, root old prob, new root vect, new root prob")
			print(node.up)
			print(node.probVect)
			print(probOldRoot)
			print(probVectRoot)
			print(probRoot)
		probRoot+= findProbRoot(probVectRoot)
		parentLKdiff=probRoot-probOldRoot
		#print("and lk would be "+str(probRoot)+" , for a difference of "+str(parentLKdiff))
		if verbose:
			print(" new root prob after adding root freqs; difference between old root prob and new")
			print(probRoot)
			print(parentLKdiff)
	else:
		if node==node.up.children[0]:
			#print("inside placeSampleOnTree 0")
			probVectParent=mergeVectorsUpDown(node.up.probVectUpRight,node.dist/2,node.probVect,node.dist/2,mutMatrix)
		else:
			#print("inside placeSampleOnTree 1")
			probVectParent=mergeVectorsUpDown(node.up.probVectUpLeft,node.dist/2,node.probVect,node.dist/2,mutMatrix)
		#print("If inserted as parent, it would create vector")
		#print(probVectParent)
		parentLKdiff=appendProb(probVectParent,newPartials,bLen)
		#print("With lk contribution "+str(parentLKdiff))
		if verbose:
			print(" non-root parent placement prob and vector:")
			print(probVectParent)
			print(parentLKdiff)
	
	#append as direct descendant of node
	if newChildLK>=parentLKdiff and newChildLK>=childBestLK:
		#now try different lengths for the new branch
		#print("Going for direct descendant of node")
		LK1=newChildLK
		bestLen=bLen
		while bestLen>0.1*bLen:
			newBLen=bestLen/2
			probChild=appendProb(node.probVectTot,newPartials,newBLen)
			#print("New len "+str(newBLen)+" give lk "+str(probChild))
			if probChild>LK1:
				LK1=probChild
				bestLen=newBLen
			else:
				break
		if bestLen>0.7*bLen:
			while bestLen<10*bLen:
				newBLen=bestLen*2
				probChild=appendProb(node.probVectTot,newPartials,newBLen)
				#print("New len "+str(newBLen)+" give lk "+str(probChild))
				if probChild>LK1:
					LK1=probChild
					bestLen=newBLen
				else:
					break
		LK0=appendProb(node.probVectTot,newPartials,0.0)
		#print("len 0 give lk "+str(LK0))
		if LK0>LK1:
			bestLen=0.0
		#print("Bestlen "+str(bestLen))

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
			#print("new tot root parent")
			#print(newRoot.probVectTot)
			newRoot.children[1].probVectTot=mergeVectorsUpDown(newRoot.probVectUpLeft,bestLen,newPartials,0.0,mutMatrix)
			#print("new child tot")
			#print(newRoot.children[1].probVectTot)
			newRoot.children[1].minorSequences=[]
			newRoot.children[1].probVect=newPartials
			if verbose:
				print("new root added to tree, identical to previous root")
				print(newRoot.probVect)
				print(newRoot.children[0].probVect)
				print(newRoot.children[1].probVect)
			#CALL FUNCTIONS TO UPDATE PROBABILITY VECTORS IN DESCENDANTS
			updatePartialsFromTop(node,newRoot.probVectUpRight,mutMatrix)
			return newRoot

		#add parent to node, but coinciding with it
		else:
			#print("Inside placeSampleOnTree, direct child of existing node, bLen "+str(bestLen))
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
			#print("new upRight vector, direct descendant non-root")
			newVect=mergeVectorsUpDown(vectUp,newInternalNode.dist,newPartials,bestLen,mutMatrix)
			newVect=shorten(newVect)
			newInternalNode.probVectUpRight=newVect
			#print(newInternalNode.probVectUpRight)
			newInternalNode.probVectUpLeft=node.probVectTot
			newVect=mergeVectors(node.probVect,0.0,newPartials,bestLen,mutMatrix)
			newVect=shorten(newVect)
			newInternalNode.probVect=newVect
			#print("new probVect")
			#print(newInternalNode.probVect)
			newVect=mergeVectorsUpDown(newInternalNode.probVectUpLeft,0.0,newPartials,bestLen,mutMatrix)
			newVect=shorten(newVect)
			newInternalNode.probVectTot=newVect
			#print("new tot vector")
			#print(newInternalNode.probVectTot)
			#print("new child tot")
			newVect=mergeVectorsUpDown(newInternalNode.probVectUpLeft,bestLen,newPartials,0.0,mutMatrix)
			newVect=shorten(newVect)
			newInternalNode.children[1].probVectTot=newVect
			#print("new child tot vector")
			#print(newInternalNode.children[1].probVectTot)
			if verbose:
				print("new internal node and child node added to tree")
				print(newInternalNode.probVect)
				print(newInternalNode.probVectUpRight)
				print(newInternalNode.probVectUpLeft)
				print(newInternalNode.probVectTot)
			#CALL FUNCTIONS TO UPDATE PROBABILITY VECTORS IN DESCENDANTS and ANCESTORS
			updatePartialsFromTop(node,newInternalNode.probVectUpRight,mutMatrix)
			#print("updatePartialsFromBottom from internal node")
			updatePartialsFromBottom(newInternalNode.up,newInternalNode.probVect,child,newInternalNode,mutMatrix)

	#put new sample as descendant of a new internal node on one of the children branches
	elif childBestLK>=parentLKdiff: 
		#print("Going for sibling of a descendant of node")
		if childBestNode==childBestNode.up.children[0]:
			child=0
			vectUp=childBestNode.up.probVectUpRight
		else:
			child=1
			vectUp=childBestNode.up.probVectUpLeft
		bestSplit=0.5
		bestSplitLK=childBestLK
		#print("Best likelihood so far "+str(bestSplitLK))
		#try different positions on the existing branch
		while bestSplit>0.01:
			newSplit=bestSplit/2
			probVectChildNew=mergeVectorsUpDown(vectUp,childBestNode.dist*newSplit,childBestNode.probVect,childBestNode.dist*(1.0-newSplit),mutMatrix)
			probChild=appendProb(probVectChildNew,newPartials,bLen)
			#print("new split "+str(newSplit)+" leads to lk "+str(probChild)+"  top dist "+str(childBestNode.dist*newSplit)+"  bottom dist "+str(childBestNode.dist*(1.0-newSplit))+" and vector:")
			#print(probVectChildNew)
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
				#print("new split "+str(newSplit)+" leads to lk "+str(probChild)+"  bottom dist "+str(childBestNode.dist*newSplit)+"  top dist "+str(childBestNode.dist*(1.0-newSplit))+" and vector:")
				#print(probVectChildNew)
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
			#print("new bLen "+str(newBLen)+" leads to lk "+str(probChild)+" ")
			if probChild>LK1:
				LK1=probChild
				bestLen=newBLen
			else:
				break
		if bestLen>0.7*bLen:
			while bestLen<10*bLen:
				newBLen=bestLen*2
				probChild=appendProb(childBestVect,newPartials,newBLen)
				#print("new bLen "+str(newBLen)+" leads to lk "+str(probChild)+" ")
				if probChild>LK1:
					LK1=probChild
					bestLen=newBLen
				else:
					break
		LK0=appendProb(childBestVect,newPartials,0.0)
		#print("new bLen 0.0 leads to lk "+str(LK0)+" ")
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
		#print("Inside placeSampleOnTree, internal node, new upRight")
		newVect=mergeVectorsUpDown(vectUp,newInternalNode.dist,newPartials,bestLen,mutMatrix)
		newVect=shorten(newVect)
		newInternalNode.probVectUpRight=newVect
		#print("new upRight vector")
		#print(newInternalNode.probVectUpRight)
		newInternalNode.probVectUpLeft=childBestVect
		newVect=mergeVectors(childBestNode.probVect,childBestNode.dist,newPartials,bestLen,mutMatrix)
		newVect=shorten(newVect)
		newInternalNode.probVect=newVect
		#print("Inside placeSampleOnTree, internal node, new tot")
		newVect=mergeVectorsUpDown(childBestVect,0.0,newPartials,bestLen,mutMatrix)
		newVect=shorten(newVect)
		newInternalNode.probVectTot=newVect
		#print("new tot vector")
		#print(newInternalNode.probVectTot)
		#print("Inside placeSampleOnTree, internal node, new child tot")
		newVect=mergeVectorsUpDown(childBestVect,bestLen,newPartials,0.0,mutMatrix)
		newVect=shorten(newVect)
		newInternalNode.children[1].probVectTot=newVect
		#print("new child tot vector")
		#print(newInternalNode.children[1].probVectTot)
		if verbose:
			print("new internal node added to tree")
			print(newInternalNode.probVect)
			print(newInternalNode.probVectUpRight)
			print(newInternalNode.probVectUpLeft)
			print(newInternalNode.probVectTot)
		#CALL FUNCTIONS TO UPDATE PROBABILITY VECTORS IN DESCENDANTS and ANCESTORS
		updatePartialsFromTop(childBestNode,newInternalNode.probVectUpRight,mutMatrix)
		#print("updatePartialsFromBottom from internal node")
		updatePartialsFromBottom(newInternalNode.up,newInternalNode.probVect,child,newInternalNode,mutMatrix)
	
	else: #add a new parent to node
		#print("Going for sibling of node")
		#add a new root!
		if node.is_root():
			#probOldRoot = findProbRoot(node.probVect)
			#probVectRoot,probRoot = mergeVectorsRoot(node.probVect,bLen,newPartials,bLen,mutMatrix)
			#parentLKdiff=probRoot-probOldRoot
			bestLen1=bLen
			bestLK1=parentLKdiff
			rootBestVect=probVectRoot
			#print("best lk to start with: "+str(bestLK1)+" with vect")
			#print(rootBestVect)
			#try different lengths of left branch
			while bestLen1>0.1*bLen:
				newBLen=bestLen1/2
				newProbVectRoot,newProbRoot = mergeVectorsRoot(node.probVect,newBLen,newPartials,bLen,mutMatrix)
				newProbRoot+= findProbRoot(newProbVectRoot)
				LKdiffRoot=newProbRoot-probOldRoot
				#print("new bLen1 "+str(newBLen)+" leads to lk diff "+str(LKdiffRoot)+" and vector:")
				#print(newProbVectRoot)
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
					#print("new bLen1 "+str(newBLen)+" leads to lk diff "+str(LKdiffRoot)+" and vector:")
					#print(newProbVectRoot)
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
				#print("new bLen2 "+str(newBLen)+" leads to lk diff "+str(LKdiffRoot)+" and vector:")
				#print(newProbVectRoot)
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
					#print("new bLen2 "+str(newBLen)+" leads to lk diff "+str(LKdiffRoot)+" and vector:")
					#print(newProbVectRoot)
					if LKdiffRoot>bestLK1:
						bestLK1=LKdiffRoot
						bestLen2=newBLen
						rootBestVect=newProbVectRoot
					else:
						break
			newProbVectRoot,newProbRoot = mergeVectorsRoot(node.probVect,bestLen1,newPartials,0.0,mutMatrix)
			newProbRoot+= findProbRoot(newProbVectRoot)
			LK0=newProbRoot-probOldRoot
			#print("if bLen2 0.0 leads to lk diff "+str(LK0)+" and vector:")
			#print(newProbVectRoot)
			if LK0>bestLK1:
				bestLen2=0.0
				bestLK1=LK0
				rootBestVect=newProbVectRoot
			#now creating new root
			newRoot=Tree()
			newRoot.probVect=rootBestVect
			#print("Final bLen1 "+str(bestLen1)+" final bLen2 "+str(bestLen2))
			#print(rootBestVect)
			newRoot.probVectTot=rootVector0(rootBestVect,rootFreqs)
			#print("new tot")
			#print(newRoot.probVectTot)
			newRoot.probVectUpRight=rootVector(newPartials,rootFreqs,bestLen2)
			newRoot.probVectUpLeft=rootVector(node.probVect,rootFreqs,bestLen1)
			newRoot.add_child(child=node,dist=bestLen1)
			node.up=newRoot
			newRoot.add_child(dist=bestLen2, name=sample)
			newRoot.children[1].probVect=newPartials
			#print("new child tot")
			newRoot.children[1].probVectTot=mergeVectorsUpDown(newRoot.probVectUpLeft,bestLen2,newPartials,0.0,mutMatrix)
			#print("new child tot")
			#print(newRoot.children[1].probVectTot)
			newRoot.children[1].minorSequences=[]
			newRoot.children[1].up=newRoot
			if verbose:
				print("new root added to tree")
				print(newRoot.probVect)
				print(newRoot.children[0].probVect)
				print(newRoot.children[1].probVect)
			#CALL FUNCTIONS TO UPDATE PROBABILITY VECTORS IN DESCENDANTS
			updatePartialsFromTop(node,newRoot.probVectUpRight,mutMatrix)
			return newRoot
		
		#add parent to node, knowing that node is not the root
		else:
			#print("not root, vectUp:")
			if node==node.up.children[0]:
				child=0
				vectUp=node.up.probVectUpRight
			else:
				child=1
				vectUp=node.up.probVectUpLeft
			#print(vectUp)
			bestSplit=0.5
			bestSplitLK=parentLKdiff
			childBestVect=probVectParent
			#try different positions on the existing branch
			while bestSplit>0.01:
				newSplit=bestSplit/2
				probVectParentNew=mergeVectorsUpDown(vectUp,node.dist*newSplit,node.probVect,node.dist*(1.0-newSplit),mutMatrix)
				probChild=appendProb(probVectParentNew,newPartials,bLen)
				#print("new split "+str(newSplit)+" leads to lk "+str(probChild)+"  top dist "+str(node.dist*newSplit)+"  bottom dist "+str(node.dist*(1.0-newSplit))+" and vector:")
				#print(probVectParentNew)
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
					#print("new split "+str(newSplit)+" leads to lk "+str(probChild)+"  top dist "+str(node.dist*(1.0-newSplit))+"  bottom dist "+str(node.dist*newSplit)+" and vector:")
					#print(probVectParentNew)
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
				#print("new length "+str(newBLen)+" leads to lk "+str(probChild))
				if probChild>LK1:
					LK1=probChild
					bestLen=newBLen
				else:
					break
			if bestLen>0.7*bLen:
				while bestLen<10*bLen:
					newBLen=bestLen*2
					probChild=appendProb(childBestVect,newPartials,newBLen)
					#print("new length "+str(newBLen)+" leads to lk "+str(probChild))
					if probChild>LK1:
						LK1=probChild
						bestLen=newBLen
					else:
						break
			LK0=appendProb(childBestVect,newPartials,0.0)
			#print("0 length leads to lk "+str(LK0))
			if LK0>LK1:
				bestLen=0.0
			#print("Best length: "+str(bestLen))
			#now create new internal node and append child to it
			newInternalNode=Tree()
			node.up.children[child]=newInternalNode
			newInternalNode.up=node.up
			distBottom=node.dist*(1.0-bestSplit)
			distTop=node.dist*bestSplit
			newInternalNode.add_child(dist=distBottom,child=node)
			node.up=newInternalNode
			newInternalNode.add_child(dist=bestLen, name=sample)
			#print("top dist "+str(distTop)+" bottom dist "+str(distBottom)+" new bLen "+str(bestLen))
			newInternalNode.children[1].minorSequences=[]
			newInternalNode.children[1].up=newInternalNode
			newInternalNode.dist=distTop
			newInternalNode.children[1].probVect=newPartials
			#print("Inside placeSampleOnTree, internal node, new upRight")
			newVect=mergeVectorsUpDown(vectUp,newInternalNode.dist,newPartials,bestLen,mutMatrix)
			newVect=shorten(newVect)
			newInternalNode.probVectUpRight=newVect
			#print("new upRight vector")
			#print(newInternalNode.probVectUpRight)
			newInternalNode.probVectUpLeft=childBestVect
			newVect=mergeVectors(node.probVect,node.dist,newPartials,bestLen,mutMatrix)
			newVect=shorten(newVect)
			newInternalNode.probVect=newVect
			#print("Inside placeSampleOnTree, internal node, new tot")
			newVect=mergeVectorsUpDown(childBestVect,0.0,newPartials,bestLen,mutMatrix)
			newVect=shorten(newVect)
			newInternalNode.probVectTot=newVect
			#print("new tot vector")
			#print(newInternalNode.probVectTot)
			#print("Inside placeSampleOnTree, internal node, new child tot")
			newVect=mergeVectorsUpDown(childBestVect,bestLen,newPartials,0.0,mutMatrix)
			newVect=shorten(newVect)
			newInternalNode.children[1].probVectTot=newVect
			#print("new child tot vector")
			#print(newInternalNode.children[1].probVectTot)
			if verbose:
				print("new internal node added to tree")
				print(newInternalNode.probVect)
				print(newInternalNode.probVectUpRight)
				print(newInternalNode.probVectUpLeft)
				print(newInternalNode.probVectTot)
			#CALL FUNCTIONS TO UPDATE PROBABILITY VECTORS IN DESCENDANTS and ANCESTORS
			updatePartialsFromTop(node,newInternalNode.probVectUpRight,mutMatrix)
			#print("updatePartialsFromBottom from internal node")
			updatePartialsFromBottom(newInternalNode.up,newInternalNode.probVect,child,newInternalNode,mutMatrix)

	return None





#IMPROVEMENTS THAT NEED TO BE DONE:
#2) Use 0-length branches to speed up calculations of vector merges?
#2) TEST EACH FUNCTION ONE AT THE TIME BY PRINTING ITS INPUT OUTPUT ON THE TOTAL RUN


stopAt=2000
sampleOfInterest="EPI_ISL_425829"
sampleOfInterest="EPI_ISL_422130"
sampleOfInterest=""
if newPhylogeny:
	#calculate list of distances from ref
	start = time.time()
	print("Creating list of distances from ref")
	#distances=distancesFromRef(data)
	#distances=sortByDate(dates)
	distances=distancesFromRefPunishNs(data)
	time2 = time.time() - start
	#print(distances)
	print("Distance calculation finished in time "+str(time2))
	#distances[i][0] is name, distances[i][1] is number of non-ambiguous bp, distances[i][2] is number of differences from the reference.

	#extract root genome among those closest to the reference but not empty
	root=distances.pop(0)
	print("initial root:")
	print(root)
	print(data[root[1]])
	print("root freqs")
	print(rootFreqs)

	#we only consider one possible bLen!=0 for now
	bLen=1.0/lRef

	#initialize tree to just the initial root sample
	t1=Tree()
	t1.name=root[1]
	t1.probVect=probVectTerminalNode(data[root[1]],0.0)
	t1.probVectTot=rootVector0(t1.probVect,rootFreqs)
	t1.minorSequences=[]
	#LKtot, probVectTot = finalLKwithVector(t1.probVect,ref,rootFreqsLog,rootFreqs,0.0)
	#probVectTot = rootVector(t1.probVect,rootFreqs,0.0)
	#t1.probVectTot=probVectTot
	numSamples=0
	for d in distances:
		numSamples+=1
		sample=d[1]
		newPartials=probVectTerminalNode(data[sample],0.0)
		if (numSamples%100)==0:
			print("Sample num "+str(numSamples))
			print(sample)
			print(newPartials)
		if sample==sampleOfInterest:
			print("Sample of interest")
			verbose=True
			print(t1)
		node , bestNewLK=findBestParent(t1,newPartials,sample,bLen,float('-inf'),t1,0)
		if bestNewLK<0.5:
			#print("Major sequence not found, now placing sequence with partials")
			newRoot=placeSampleOnTree(node,newPartials,sample,bLen,bestNewLK)
			if newRoot!=None:
				t1=newRoot
			#print("Root vector:")
			#print(t1.probVect)
		if numSamples>=stopAt or sample==sampleOfInterest:
			print("Final tree:")
			treeString=(t1.write()).replace(")1:","):")
			print(treeString)
			file=open(pathSimu+str(stopAt)+"samplesPhylogeny_iterativeFastLK.tree","w")
			file.write(treeString)
			file.close()
			#t1.write(outfile=pathSimu+str(stopAt)+"samplesPhylogeny_iterativeFastLK.tree")
			#print(t1)
			print(totalMissedMinors[0])
			exit()


	#t1.up
	






























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
			t1.dist+=phylo.dist
			if len(phylo.children)==0:
				t1.name=phylo.name
				#print(t1.name)
				#print("New tip branch length: "+str(t1.dist))
				#print("Old tip branch length: "+str(phylo.dist))
				#t1.dist+=phylo.dist
			newChildren=[]
			for c in phylo.children:
				if c.descendants>0:
					nChildren+=1
					newChildren.append(c)
			#if nChildren>1:
			#	print("Number of children: "+str(nChildren))
			if nChildren==1:
				#t1.dist+=phylo.dist
				extractSubtree(newChildren[0],t1,minBlen=minBlen)
			if nChildren>1:
				for i in range(nChildren-1):
					t1.add_child(dist=minBlen)
					extractSubtree(newChildren[i],t1.children[-1],minBlen=minBlen)
					#print("Branch length of new child: "+str(t1.children[-1].dist))
					#print("Its old bLen: "+str(newChildren[i].dist))
					if i!=(nChildren-2):
						#print("splitting polytomy")
						t1.add_child(dist=minBlen)
						#print(t1.dist)
						t1=t1.children[1]
				t1.add_child(dist=minBlen)
				extractSubtree(newChildren[-1],t1.children[-1],minBlen=minBlen)
				#print("Branch length of second child: "+str(t1.children[-1].dist))
				#print("Its old bLen: "+str(newChildren[-1].dist))




def dataSubgenome(genomeStart,genomeEnd,data,newSamples):
	newData={}
	for s in newSamples:
		newData[s]=[]
		for m in range(len(data[s])):
			if data[s][m][1]<=genomeEnd and data[s][m][1]>=genomeStart:
				newData[s].append(data[s][m])
	return newData

def logLKdiff(data,newSamples,string,ref):
		#write data files for subsample
		#concise file
		subsample=len(newSamples)
		#if debug:
		#	fileO=open(pathSimu+"2021-03-31_debugging.txt","w")
		#else:
		fileO=open(pathSimu+"2021-03-31_"+string+".txt","w")
		for s in newSamples:
			fileO.write(">"+s+"\n")
			for m in data[s]:
				if len(m)==2:
						fileO.write(m[0]+"\t"+str(m[1])+"\n")
				else:
						fileO.write(m[0]+"\t"+str(m[1])+"\t"+str(m[2])+"\n")

		fileO.close()
		#alignment file
		#if debug:
		#	fileO=open(pathSimu+"2021-03-31_debugging.phylip","w")
		#else:
		fileO=open(pathSimu+"2021-03-31_"+string+".phylip","w")
		lRef=len(ref)
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
		print("Starting new Felsenstein on  "+string+" ")
		start = time.time()
		#if debug:
		#	dataSample=readConciseAlignment(pathSimu+"2021-03-31_debugging.txt")
		#else:
		dataSample=readConciseAlignment(pathSimu+"2021-03-31_"+string+".txt")
		#print(dataSample)
		ref, cumulativeBases, rootFreqs, rootFreqsLog = collectReference(pathSimu+"EPI_ISL_402124_lowercase.fasta")
		lRef=len(ref)
		#if debug:
		#	phylo2 = Tree(pathSimu+"2021-03-31_debugging.tree")
		#else:
		phylo2 = Tree(pathSimu+"2021-03-31_"+string+".tree")
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
		time2 = time.time() - start
		print("Time taken to prepare for likelihood calculation "+str(time2))
		start = time.time()
		cumulPartLk, probVect = pruningFast(phylo2,mutMatrix,dataSample,ref, cumulativeBases)
		print("Fast likelihood finished")
		print(cumulPartLk)
		print(probVect)
		cumulPartLk+=finalLK(probVect,ref,rootFreqsLog,rootFreqs)
		time3 = time.time() - start
		print("Final log-likelihood after finalization:")
		print(cumulPartLk)
		print("Time for likelihood: "+str(time3))
		print("Time taken "+str(time2+time3))

		#Run PhyML
		if runPhyml:
			start = time.time()
			#if debug:
			#	os.system("/Applications/phyml-master/src/phyml --no_memory_check --leave_duplicates --run_id debugging_phyml --print_trace -o n -c 1 -b 0 -m 0.04,0.3,0.1,0.02,1.0,1.0 -f 0.25,0.25,0.25,0.25 -i "+pathSimu+"2021-03-31_debugging.phylip -u "+pathSimu+"2021-03-31_debugging.tree")
			#else:
			os.system("/Applications/phyml-master/src/phyml --no_memory_check --leave_duplicates --run_id "+string+"_phyml --print_trace -o n -c 1 -b 0 -m 0.04,0.3,0.1,0.02,1.0,1.0 -f 0.25,0.25,0.25,0.25 -i "+pathSimu+"2021-03-31_"+string+".phylip -u "+pathSimu+"2021-03-31_"+string+".tree")
			time2 = time.time() - start
			print("Time taken by PhyML "+str(time2))



			#if debug:
			#	filePhy=open(pathSimu+"2021-03-31_debugging.phylip_phyml_stats_debugging_phyml.txt")
			#else:
			filePhy=open(pathSimu+"2021-03-31_"+string+".phylip_phyml_stats_"+string+"_phyml.txt")
			line=filePhy.readline()
			linelist=line.split()
			while len(linelist)<2 or linelist[1]!="Log-likelihood:":
				line=filePhy.readline()
				linelist=line.split()
			phymlLK=float(linelist[2])

			if runIQtree:
				print("I can only compare with one of IQtree or PhyML at the time, this time skipping IQtree")
		elif runIQtree:
			start = time.time()
			print("Starting IQtree")
			os.system("/Applications/iqtree-1.6.12-MacOSX/bin/iqtree -s "+pathSimu+"2021-03-31_"+string+".phylip -st DNA -te "+pathSimu+"2021-03-31_"+string+".tree -pre "+string+"_iqtree1 -blfix -keep-ident -m GTR\{0.04,0.3,0.1,0.02,1.0,1.0\}+FQ -quiet -redo  -nt 1 ")
			time2 = time.time() - start
			print("Time taken by IQtree "+str(time2))

			start = time.time()
			print("Starting IQtree with 4 threads")
			#os.system("/Applications/iqtree-1.6.12-MacOSX/bin/iqtree -s "+pathSimu+"2021-03-31_"+string+".phylip -st DNA -te "+pathSimu+"2021-03-31_"+string+".tree -pre "+string+"_iqtree1_4threads -blfix -keep-ident -m GTR\{0.04,0.3,0.1,0.02,1.0,1.0\}+FQ -quiet -redo -nt 4 ")
			time2 = time.time() - start
			print("Time taken by IQtree with 4 threads "+str(time2))

			start = time.time()
			print("Starting IQtree with no partialLK memory")
			#os.system("/Applications/iqtree-1.6.12-MacOSX/bin/iqtree -s "+pathSimu+"2021-03-31_"+string+".phylip -st DNA -te "+pathSimu+"2021-03-31_"+string+".tree -pre "+string+"_iqtree1_noMem -blfix -keep-ident -m GTR\{0.04,0.3,0.1,0.02,1.0,1.0\}+FQ -quiet -redo -mem 0  -nt 1 ")
			time2 = time.time() - start
			print("Time taken by IQtree with no partialLK memory "+str(time2))

			start = time.time()
			print("Starting IQtree with 4 threads, no partialLK memory")
			#os.system("/Applications/iqtree-1.6.12-MacOSX/bin/iqtree -s "+pathSimu+"2021-03-31_"+string+".phylip -st DNA -te "+pathSimu+"2021-03-31_"+string+".tree -pre "+string+"_iqtree1_4threads_noMem -blfix -keep-ident -m GTR\{0.04,0.3,0.1,0.02,1.0,1.0\}+FQ -quiet -redo -mem 0 -nt 4 ")
			time2 = time.time() - start
			print("Time taken by IQtree with 4 threads no partialLK memory "+str(time2))

			start = time.time()
			print("Starting IQtree, estimate bLengths")
			#os.system("/Applications/iqtree-1.6.12-MacOSX/bin/iqtree -s "+pathSimu+"2021-03-31_"+string+".phylip -st DNA -te "+pathSimu+"2021-03-31_"+string+".tree -pre "+string+"_iqtree1_bLen -nt 1 -keep-ident -m GTR\{0.04,0.3,0.1,0.02,1.0,1.0\}+FQ -quiet -redo ")
			time2 = time.time() - start
			print("Time taken by IQtree, also bLengths "+str(time2))

			start = time.time()
			print("Starting IQtree with no partialLK memory, estimate bLengths")
			#os.system("/Applications/iqtree-1.6.12-MacOSX/bin/iqtree -s "+pathSimu+"2021-03-31_"+string+".phylip -st DNA -te "+pathSimu+"2021-03-31_"+string+".tree -pre "+string+"_iqtree1_bLen_noMem -nt 1 -keep-ident -m GTR\{0.04,0.3,0.1,0.02,1.0,1.0\}+FQ -quiet -redo -mem 0 ")
			time2 = time.time() - start
			print("Time taken by IQtree with no partialLK memory, also bLengths "+str(time2))

		else:
			print("For comparinson to work, you have to choose one between PhyML or IQtree in the options.")

		return phymlLK-cumulPartLk


						
if subsample>0:
	#sample "subsample" many leaves of the tree
	#samples=data.keys()
	for repeat in range(nRepeats):
		if debug:
			print("Seed "+str(seed+repeat))
			random.seed(seed+repeat)
		newSamples=random.sample(samples,subsample)
		#print(newSamples)
		#phylo=Tree("((((seq1:0.001,seq2:0.001):0.001,seq3:0.001):0.001,(seq4:0.001,seq5:0.001):0.001):0.001,(((seq6:0.001,seq7:0.001):0.001,seq8:0.001):0.001,(seq9:0.001,seq10:0.001):0.001):0.001):0.001;")
		#newSamples=["seq1","seq3","seq5","seq6","seq9"]

		#extract subtree containing these samples and replace 0 lengths with non-zero lengths and multifurcations with bifurcations.
		countLeaves(phylo,newSamples)
		t1 = Tree()	
		extractSubtree(phylo,t1)
		print(str(len(t1))+" sequences in the new tree.")
		#exit()
		#print(t1.write())
		#exit()
		if debug:
			t1.write(outfile=pathSimu+"2021-03-31_debugging.tree")
		else:
			t1.write(outfile=pathSimu+"2021-03-31_subsample"+str(subsample)+"_repeat"+str(repeat)+".tree")

		
		if debug:
			LKdiff=logLKdiff(data,newSamples,"debugging",ref)
		else:
			LKdiff=logLKdiff(data,newSamples,"subsample"+str(subsample)+"_repeat"+str(repeat),ref)
			print("likelihood difference is: "+str(LKdiff))

		LKlimit=0.3
		if debug and abs(LKdiff)>LKlimit:
			print("Found discrepancy in likelihood")
			genomeEnd=29981
			genomeStart=0
			finished=False
			while (genomeEnd - genomeStart >10) and (not finished):
				print("\n"+"Genome start "+str(genomeStart)+" genome end "+str(genomeEnd))
				newGenomeEnd=genomeStart+(genomeEnd-genomeStart)/2
				print("\n"+"New genome end  "+str(newGenomeEnd))
				if logLKdiff(dataSubgenome(genomeStart,newGenomeEnd,data,newSamples),newSamples,"debugging",ref)<LKlimit:
					genomeStart=newGenomeEnd
					print("\n"+"New genome start  "+str(genomeStart))
					if logLKdiff(dataSubgenome(genomeStart,genomeEnd,data,newSamples),newSamples,"debugging",ref)<LKlimit:
						print("There is not a single site with worrying LK contribution")
						finished=True
						#exit()
				else:
					genomeEnd=newGenomeEnd
			if (not finished):
				print("Found portion of the genome with significant LK contribution")
				exit()





			
					
					







#estimate susbstitution rates?
#optimize branch lengths?
#Optimize root position?
#Calculate pairwise distances?
#reduce tree by removing sequences that are basically identical to others?
#try to improve topology?

exit()








#for the root, finalize the partial likelihood by multiplying by root frequencies, and return new vector (only affects O states)
#OBSOLETE - now using rootVector0 which uses the new format and tracks of which entries need to be multiplied by root frequencies.
# def multipyByRootProbs(probVect,rootFreqs):
# 	newProbVect=[]
# 	for entry in probVect:
# 		newEntry=list(entry)
# 		if entry[0]=="O":
# 			newEntry[4]=np.ones(4)
# 			for i in range4:
# 				newEntry[4][i]=entry[4][i]*rootFreqs[i]
# 		newProbVect.append(newEntry)
# 	return newProbVect


#for the root, finalize the likelihood and also return a vector of the root probability of each state
#OBSOLETE!
# def finalLKwithVector(probVect,ref,rootFreqsLog,rootFreqs,bLen):
# 	logLK=0.0
# 	newProbVect=[]
# 	for entry in probVect:
# 		#if entry[3]<thresholdProb:
# 			#approximate assuming the distance is 0.0
# 			if entry[0]=="R":
# 				pos=entry[1]
# 				end=entry[2]+pos-1
# 				for i in range(4):
# 					logLK+=rootFreqsLog[i]*(cumulativeBases[end][i]-cumulativeBases[pos-1][i])
# 				newProbVect.append(("R",pos,entry[2],entry[3]+bLen))
# 			elif entry[0]=="N":
# 				newProbVect.append(("N",entry[1],entry[2],0.0))
# 			elif entry[0]=="O":
# 				#tot=0.0
# 				newProbVect.append(("O",entry[1],entry[2],0.0,[]))
# 				for i in range(4):
# 					#tot+=rootFreqs[i]*entry[4][i]
# 					tot=0.0
# 					for j in range4:
# 						if i==j:
# 							tot+=rootFreqs[i]*(1.0+mutMatrix[i][i]*(entry[3]+bLen))*entry[4][j]
# 						else:
# 							tot+=rootFreqs[i]*mutMatrix[i][j]*(entry[3]+bLen)*entry[4][j]
# 					#newProbVect[-1][4].append(rootFreqs[i]*entry[4][i])
# 					newProbVect[-1][4].append(tot)
# 				#logLK+=math.log(tot)
# 			else:
# 				i1=alleles[entry[0]]
# 				newProbVect.append((entry[0],entry[1],entry[2],entry[3]+bLen))
# 				logLK+=rootFreqsLog[i1]
# 	return logLK, newProbVect


#we know that sample "sample", with partials "newPartials", is best placed as child of node resulting in logLK contribution of newChildLK
# now explore exactly which position near node is best for placement, and which are the best branch lengths,
#then add the sample at that position of the tree, and update all the internal probability vectors.
#OBSOLETE!
def placeSampleOnTreeOld(node,newPartials,sample,bLen,newChildLK):
	#place as sibling of either child
	if len(node.children)==2:
		if verbose:
			print("Trying to place sample at new internal node; current UpRight, child0, UpLef, child1.")
			print(node.probVectUpRight)
			print(node.children[0].probVect)
			print(node.probVectUpLeft)
			print(node.children[1].probVect)
		print("Inside placeSampleOnTree, child 1")
		probVectChid1=mergeVectorsUpDown(node.probVectUpRight,bLen,node.children[0].probVect,bLen,mutMatrix)
		print("child 1 vector")
		print(probVectChid1)
		probChild1=appendProb(probVectChid1,newPartials,bLen)
		print("Inside placeSampleOnTree, child 1")
		probVectChid2=mergeVectorsUpDown(node.probVectUpLeft,bLen,node.children[1].probVect,bLen,mutMatrix)
		print("child 2 vector")
		print(probVectChid2)
		probChild2=appendProb(probVectChid2,newPartials,bLen)
		if verbose:
			print("Vector 1 to append to, LK1, vector 2 to append to, LK2")
			print(probVectChid1)
			print(probChild1)
			print(probVectChid2)
			print(probChild2)
		if probChild2>probChild1:
			probChild=probChild2
			probVectChid=probVectChid2
			vectUp=node.probVectUpLeft
			child=1
		else:
			probChild=probChild1
			probVectChid=probVectChid1
			vectUp=node.probVectUpRight
			child=0
		if verbose:
			print("probChild")
			print(probChild)
	addChild=True
	#if node is root, try to place as sibling of the current root.
	#print(node.up)
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
		if verbose:
			print(" new root prob after adding root freqs; difference between old root prob and new")
			print(probRoot)
			print(probRoot-probOldRoot)
		if (len(node.children)<1) or (probRoot-probOldRoot>probChild):
			#CREATE NEW ROOT
			if verbose:
				print("Creating new root")
			newRoot=Tree()
			newRoot.probVect=probVectRoot
			#LKtot, probVectTot = finalLKwithVector(probVectRoot,ref,rootFreqsLog,rootFreqs,0.0)
			probVectTot = rootVector0(probVectRoot,rootFreqs)
			#probVectTot = rootVector(probVectRoot,rootFreqs,0.0)
			newRoot.probVectTot=probVectTot
			#LKUpR, newRoot.probVectUpRight=finalLKwithVector(newPartials,ref,rootFreqsLog,rootFreqs,bLen)
			#LKUpL, newRoot.probVectUpLeft=finalLKwithVector(node.probVect,ref,rootFreqsLog,rootFreqs,bLen)
			newRoot.probVectUpRight=rootVector(newPartials,rootFreqs,bLen)
			newRoot.probVectUpLeft=rootVector(node.probVect,rootFreqs,bLen)
			newRoot.add_child(child=node,dist=bLen)
			node.up=newRoot
			newRoot.add_child(dist=bLen, name=sample)
			newRoot.children[1].probVect=newPartials
			print("Inside placeSampleOnTree, new tot")
			newRoot.children[1].probVectTot=mergeVectorsUpDown(newRoot.probVectUpLeft,bLen,newPartials,0.0,mutMatrix)
			print("new tot")
			print(newRoot.children[1].probVectTot)
			newRoot.children[1].minorSequences=[]
			newRoot.children[1].up=newRoot
			if verbose:
				print("new root added to tree")
				print(newRoot.probVect)
				print(newRoot.children[0].probVect)
				print(newRoot.children[1].probVect)
			#CALL FUNCTIONS TO UPDATE PROBABILITY VECTORS IN DESCENDANTS
			updatePartialsFromTop(node,newRoot.probVectUpRight,mutMatrix)
			addChild=False
			return newRoot
		else:
			if verbose:
				print("Adding child node")
			
	#CREATE NEW INTERNAL NODE AND ADD CHILD INSTEAD
	if addChild:
		if len(node.children)==2:
			newInternalNode=Tree()
			newInternalNode.add_child(dist=bLen,child=node.children[child])
			node.children[child].up=newInternalNode
			newInternalNode.add_child(dist=bLen, name=sample)
			newInternalNode.children[1].minorSequences=[]
			newInternalNode.children[1].up=newInternalNode
			node.children[child]=newInternalNode
			newInternalNode.up=node
			newInternalNode.dist=bLen
			newInternalNode.children[1].probVect=newPartials
			print("Inside placeSampleOnTree, internal node, new upRight")
			newInternalNode.probVectUpRight=mergeVectorsUpDown(vectUp,bLen,newPartials,bLen,mutMatrix)
			print("new vector")
			print(newInternalNode.probVectUpRight)
			newInternalNode.probVectUpLeft=probVectChid
			newInternalNode.probVect=mergeVectors(newInternalNode.children[0].probVect,bLen,newPartials,bLen,mutMatrix)
			print("Inside placeSampleOnTree, internal node, new tot")
			newInternalNode.probVectTot=mergeVectorsUpDown(newInternalNode.probVectUpLeft,0.0,newPartials,bLen,mutMatrix)
			print("new vector")
			print(newInternalNode.probVectTot)
			print("Inside placeSampleOnTree, internal node, new child tot")
			newInternalNode.children[1].probVectTot=mergeVectorsUpDown(newInternalNode.probVectUpLeft,bLen,newPartials,0.0,mutMatrix)
			print("new vector")
			print(newInternalNode.children[1].probVectTot)
			if verbose:
				print("new internal node added to tree")
				print(newInternalNode.probVect)
				print(newInternalNode.probVectUpRight)
				print(newInternalNode.probVectUpLeft)
				print(newInternalNode.probVectTot)
			#CALL FUNCTIONS TO UPDATE PROBABILITY VECTORS IN DESCENDANTS and ANCESTORS
			updatePartialsFromTop(newInternalNode.children[0],newInternalNode.probVectUpRight,mutMatrix)
			#print("updatePartialsFromBottom from internal node")
			updatePartialsFromBottom(node,newInternalNode.probVect,child,newInternalNode,mutMatrix)
		#adding child of leaf
		else:
			newInternalNode=Tree()
			if node.up.children[0]==node:
				child=0
				vectUp=node.up.probVectUpRight
			else:
				child=1
				vectUp=node.up.probVectUpLeft
			node.up.children[child]=newInternalNode
			newInternalNode.up=node.up
			newInternalNode.dist=node.dist
			newInternalNode.add_child(dist=0.0,child=node)
			node.dist=0.0
			node.up=newInternalNode
			newInternalNode.add_child(dist=bLen, name=sample)
			newInternalNode.children[1].minorSequences=[]
			newInternalNode.children[1].up=newInternalNode
			newInternalNode.children[1].probVect=newPartials
			print("Inside placeSampleOnTree, leaf node, new upRight")
			newInternalNode.probVectUpRight=mergeVectorsUpDown(vectUp,newInternalNode.dist,newPartials,bLen,mutMatrix)
			print("new vector")
			print(newInternalNode.probVectUpRight)
			newInternalNode.probVectUpLeft=node.probVectTot
			newInternalNode.probVect=mergeVectors(newInternalNode.children[0].probVect,0.0,newPartials,bLen,mutMatrix)
			print("Inside placeSampleOnTree, leaf node, new tot")
			newInternalNode.probVectTot=mergeVectorsUpDown(newInternalNode.probVectUpLeft,0.0,newPartials,bLen,mutMatrix)
			print("new vector")
			print(newInternalNode.probVectTot)
			print("Inside placeSampleOnTree, leaf node, new child tot")
			newInternalNode.children[1].probVectTot=mergeVectorsUpDown(newInternalNode.probVectUpLeft,bLen,newPartials,0.0,mutMatrix)
			print("new vector")
			print(newInternalNode.children[1].probVectTot)
			if verbose:
				print("new internal node added to tree")
				print(newInternalNode.probVect)
				print(newInternalNode.probVectUpRight)
				print(newInternalNode.probVectUpLeft)
				print(newInternalNode.probVectTot)
			#CALL FUNCTIONS TO UPDATE PROBABILITY VECTORS IN DESCENDANTS and ANCESTORS
			updatePartialsFromTop(node,newInternalNode.probVectUpRight,mutMatrix)
			#print("updatePartialsFromBottom from leaf")
			updatePartialsFromBottom(newInternalNode.up,newInternalNode.probVect,child,None,mutMatrix)
		return None











