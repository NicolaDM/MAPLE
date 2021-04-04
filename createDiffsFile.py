import sys
import os
import math
import numpy as np
import os.path
from os import path
import argparse
import matplotlib
import matplotlib.pyplot as plt
from scipy.stats import ttest_ind, ttest_ind_from_stats, chi2_contingency
from scipy.special import stdtr
from discreteMarkovChain import markovChain
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
args = parser.parse_args()

createDiffFile=args.createDiffFile
reduceDiffFile=args.reduceDiffFile
onlyNambiguities=args.onlyNambiguities
useLogs=args.useLogs
thresholdProb=args.thresholdProb

pathSimu=args.path

alleles={"A":0,"C":1,"G":2,"T":3}
allelesList=["A","C","G","T"]
allelesLow={"a":0,"c":1,"g":2,"t":3}
allelesListLow=["a","c","g","t"]
ambiguities={"y":[0.0,1.0,0.0,1.0],"r":[1.0,0.0,1.0,0.0],"w":[1.0,0.0,0.0,1.0],"s":[0.0,1.0,1.0,0.0],"k":[0.0,0.0,1.0,1.0],"m":[1.0,1.0,0.0,0.0],"d":[1.0,0.0,1.0,1.0],"v":[1.0,1.0,1.0,0.0],"h":[1.0,1.0,0.0,1.0],"b":[0.0,1.0,1.0,1.0]}

#collect reference
file=open(pathSimu+"EPI_ISL_402124_lowercase.fasta")
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



if reduceDiffFile:
	#read sequence data from file
	start = time.time()
	fileI=open(pathSimu+"2021-03-31_unmasked_differences.txt")
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
	print("Time to read DNA data file: "+str(time2))
	print(str(nSeqs)+" sequences in file.")

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
	start = time.time()
	fileI=open(pathSimu+"2021-03-31_unmasked_differences_reduced.txt")
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


print(str(len(data.keys()))+" sequences in the concise DNA data file")


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
				entryOld[2]+=entryNew[2]
	newVec.append(entryOld)
	return newVec




#preliminary nuc mutation rate matrix, from De Maio et al 2021
mutMatrix=[[0.0,0.039,0.310,0.123],[0.140,0.0,0.022,3.028],[0.747,0.113,0.0,2.953],[0.056,0.261,0.036,0.0]]
mutMatrix[0][0]=-(mutMatrix[0][1]+mutMatrix[0][2]+mutMatrix[0][3])
mutMatrix[1][1]=-(mutMatrix[1][0]+mutMatrix[1][2]+mutMatrix[1][3])
mutMatrix[2][2]=-(mutMatrix[2][0]+mutMatrix[2][1]+mutMatrix[2][3])
mutMatrix[3][3]=-(mutMatrix[3][0]+mutMatrix[3][1]+mutMatrix[3][2])
#Calculate likelihood?
def pruningFast(node,mutMatrix):
	children=node.children
	bLen=node.dist
	#cumulative contribution of the subtree to the tree likelihood (excluding probs in the vector probVect)
	cumulPartLk=0.0
	#tree tip
	if len(children)==0:
		if not (node.name in data):
			probVect=[["N",1,refL,bLen]]
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
				print("Ambiguity in sample")
				print(probVect)
				exit()
		node.cumulPartLk=cumulPartLk
		node.probVect=probVect
		print("\n"+"Terminal node "+node.name+" "+str(node.dist))
		print(diffs)
		print(probVect)
		return cumulPartLk, probVect
	
	#internal node
	else:
		cumulPartLk1, probVect1= pruningFast(children[0],mutMatrix)
		cumulPartLk+=cumulPartLk1
		print("\n"+"Internal node "+node.name)
		print("First child "+children[0].name)
		print(probVect1)
		print(cumulPartLk1)
		for c in range(len(children)-1):
			newChild=children[c+1]
			cumulPartLk2, probVect2= pruningFast(newChild,mutMatrix)
			cumulPartLk+=cumulPartLk2
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
			#nBasesOdd=np.zeros(4,dtype=int)
			while True:
				print("pos "+str(pos)+" length "+str(length)+" end "+str(end)+" indexEntry1 "+str(indexEntry1)+" indexEntry2 "+str(indexEntry2))
				if length<0:
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
								cumulPartLk+=mutMatrix[i][i]*totLen*(cumulativeBases[end][i]-cumulativeBases[pos-1][i])
					elif entry2[0]=="O":
						i1=allelesLow[ref[pos-1]]
						newVec=np.ones(4)
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
						if np.sum(newVec)<thresholdProb:
								print("Issue merging R and O")
								print(entry1)
								print(entry2)
								print(probVect[-1])
								exit()
					else:
						i2=alleles[entry2[0]]
						i1=allelesLow[ref[pos-1]]
						newVec=np.ones(4)
						if entry2[3]+entry1[3]<thresholdProb:
							print("Warning: there is a substitution on a 0-length cherry. Pretending branch lengths are not exactly zero.")
							print(entry1)
							print(entry2)
							entry2[3]=thresholdProb
							entry1[3]=thresholdProb
							exit()
						for i in range4:
							if i==i2:
								newVec[i]=1.0+mutMatrix[i][i]*entry2[3]
							else:
								newVec[i]=mutMatrix[i][i2]*entry2[3]
							if i==i1:
								newVec[i]*=1.0+mutMatrix[i][i]*entry1[3]
							else:
								newVec[i]*=mutMatrix[i][i1]*entry1[3]
						probVect.append(["O",pos,1,0.0,newVec])
						if np.sum(newVec)<thresholdProb:
								print("Issue merging R and non-R")
								print(entry1)
								print(entry2)
								print(probVect[-1])
								exit()
				elif entry1[0]=="O":
					newVec=np.ones(4)
					if entry2[0]=="O":
						print("O")
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
						print("Half-way newVec")
						print(entry1)
						print(newVec)
						if entry2[0]=="R":
							print("R")
							i2=allelesLow[ref[pos-1]]
							for i in range4:
								if i==i2:
									newVec[i]*=1.0+mutMatrix[i][i]*entry2[3]
								else:
									newVec[i]*=mutMatrix[i][i2]*entry2[3]
						else:
							print("Other")
							i2=alleles[entry2[0]]
							for i in range4:
								if i==i2:
									newVec[i]*=1.0+mutMatrix[i][i]*entry2[3]
								else:
									newVec[i]*=mutMatrix[i][i2]*entry2[3]
					print(newVec)
					state, newVec, maxP =simplfy(newVec,ref[pos-1])
					probVect.append([state,pos,1,0.0,newVec])
					print(newVec)
					print(state)
					print(maxP)
					if np.sum(newVec)<thresholdProb:
								print("Issue merging O and something else")
								print(entry1)
								print(entry2)
								print(probVect[-1])
								exit()
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
							if np.sum(newVec)<thresholdProb:
								print("Issue merging non-R and O")
								print(entry1)
								print(entry2)
								print(probVect[-1])
								exit()
						else:
							if entry2[3]+entry1[3]<thresholdProb:
								print("Warning: there is a substitution on a 0-length branches junction. Pretending branch lengths are not exactly zero, but assigning a threshold value.")
								entry1[3]=thresholdProb
								entry2[3]=thresholdProb
								for i in range4:
									if i==i1:
										newVec[i]=1.0+mutMatrix[i][i]*entry1[3]
									else:
										newVec[i]=mutMatrix[i][i1]*entry1[3]
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
							if np.sum(newVec)<thresholdProb:
								print("Issue merging non-R and (R or non-R)")
								print(entry1)
								print(entry2)
								print(probVect[-1])
								exit()

				#update pos, end, etc
				print(pos)
				print(length)
				pos+=length
				print(pos)
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
				print(pos1)
				print(end1)
				if pos>end2:
					indexEntry2+=1
					entry2=probVect2[indexEntry2]
					pos2=entry2[1]
					if entry2[0]!="N" and entry2[0]!="R":
						end2=pos2
					else:
						end2=pos2+entry2[2]-1
				print(pos2)
				print(end2)
				end=end1
				if end2<end:
					end=end2
				print(end)
				length=end+1-pos
				print("New length:")
				print(length)

			print("Updated genome vector for internal node "+node.name)
			print(probVect)
			probVect1=probVect

			#check if the final  probVect can be simplified by merging consecutive entries
			probVect1 =shorten(probVect1)
			print("Shortened genome vector for internal node "+node.name)
			print(probVect1)

			#exit()

		#update probVect to include the impact of its own branch
		for entry in probVect1:
			entry[3]+=node.dist
		
		node.cumulPartLk=cumulPartLk
		node.probVect=probVect1
		print("\n"+"Internal node "+node.name+", final genome vector and likelihood contribution:")
		print(probVect1)
		print(cumulPartLk)
		exit()
		return cumulPartLk, probVect

						

				
#calculate likelihood of phylogeny
pruningFast(phylo,mutMatrix)

#now use root frequencies to finalize the likelihood, if necessary.
						
				
			
					
					







#estimate susbstitution rates?
#optimize branch lengths?
#Optimize root position?
#Calculate pairwise distances?
#reduce tree by removing sequences that are basically identical to others?
#try to improve topology?

exit()






















#collect masked sites
file=open(pathSimu+"problematic_sites_sarsCov2_22_12_2020.vcf")
line=file.readline()
masked=[]
while line[0]=="#":
	line=file.readline()
while line!="" and line!="\n":
	linelist=line.split()
	pos=int(linelist[1])
	recc=linelist[6]
	if pos>55 and pos<29804 and recc=="mask":
		masked.append(pos)
	line=file.readline()
file.close()

#genome-wide gene annotation
geneEnds=[[266,13468],[13468,21555],[21563,25384],[25393,26220],[26245,26472],[26523,27191],[27202,27387],[27394,27759],[27756,27887],[27894,28259],[28274,29533],[29558,29674]]
geneEndsNames=["ORF1a","ORF1b","S","ORF3a","E","M","ORF6","ORF7a","ORF7b","ORF8","N","ORF10"]
#counts of opportunities for mutation events: intergenic, synonymous,nonsynonymous,nonsense; normal and C->T context-dependency
countsPossibleMuts=[np.zeros((4,4)),[],np.zeros((4,4)),[],np.zeros((4,4)),[],np.zeros((4,4))]
fourFolds=[0,0,0,0]
fourFoldMuts=np.zeros((2,4,4))
contextRange=10
countsContextPossibilities=np.zeros((2,contextRange*2,4))
countsContext=np.zeros((2,contextRange*2,4))
countsContextRates=np.zeros((2,contextRange*2,4))
#NOW WITH VARIANT COUNTS
countsContextSites=np.zeros((2,contextRange*2,4))
countsContextRatesSites=np.zeros((2,contextRange*2,4))
#new vectors taking into account immediate context (for CpG and APOBEC target).
# first dimension: type of mutation (noncoding, synonymous, nonsynonymous, 4-fold degenerate site, nonsense).
#dimensions 1-3: reference bases at positions -1, 0 and 1 respectively. 
#dimension 4: destination allele.
countsPossibleMutsContext=np.zeros((5,4,4,4,4))
#no context
countsPossibleMutsTree=np.zeros((5,4,4),dtype=int)
#new vector for U selection variation across genes.
# first dimension: gene.
#dimensions 2: reference base. 
#dimension 3: destination allele.
countsPossibleMutsContextU=np.zeros((len(geneEnds),4,4))

#count possible mutation events of different types along the genome
gene=0
intoStop=np.zeros((4,4))
for pos in range(len(ref)):
	if pos<100 or pos>29802:
		continue
	if (pos+1) in masked:
		continue
	p=pos
	typeOfPos=4 #0=intergenic, 2=genic, 3=start codon, 4=stop codon
	if gene<len(geneEnds) and p>=geneEnds[gene][1]:
		gene+=1
	if gene>=len(geneEnds) or p<(geneEnds[gene][0]-1):
		typeOfPos=0
	if gene<len(geneEnds):
		if p>=(geneEnds[gene][0]-1) and p<geneEnds[gene][1]:
			if p<=(geneEnds[gene][0]+1):
				typeOfPos=3
			elif p>=(geneEnds[gene][1]-3):
				typeOfPos=4
			else:
				typeOfPos=2
	if typeOfPos==0:
		for a in range(4):
			if a!=alleles[ref[pos]]:
				countsPossibleMuts[0][alleles[ref[pos]]][a]+=1
				countsPossibleMutsContext[0,alleles[ref[pos-1]],alleles[ref[pos]],alleles[ref[pos+1]],a]+=1
				countsPossibleMutsTree[0,alleles[ref[pos]],a]+=1
	elif typeOfPos==2:
		rest=(pos+1-geneEnds[gene][0])%3
		cod1=pos-rest
		refCodon=ref[cod1:cod1+3]
		refA=_translate_str(refCodon, table)
		if refA=="*":
			print("error")
			exit()
		fourFoldCount=0
		for a in range(4):
			if a!=alleles[ref[pos]]:
				altCod=list(refCodon)
				altCod[rest]=allelesList[a]
				altCod="".join(altCod)
				altA=_translate_str(altCod, table)
				if altA=="*":	
					intoStop[alleles[ref[pos]],a]+=1
					countsPossibleMuts[6][alleles[ref[pos]]][a]+=1
					countsPossibleMutsContext[4,alleles[ref[pos-1]],alleles[ref[pos]],alleles[ref[pos+1]],a]+=1
					countsPossibleMutsTree[4,alleles[ref[pos]],a]+=1
				if altA!=refA and altA!="*":
					countsPossibleMuts[4][alleles[ref[pos]]][a]+=1
					countsPossibleMutsContext[2,alleles[ref[pos-1]],alleles[ref[pos]],alleles[ref[pos+1]],a]+=1
					countsPossibleMutsTree[2,alleles[ref[pos]],a]+=1
				elif altA==refA:
					fourFoldCount+=1
					countsPossibleMuts[2][alleles[ref[pos]]][a]+=1
					countsPossibleMutsContext[1,alleles[ref[pos-1]],alleles[ref[pos]],alleles[ref[pos+1]],a]+=1
					countsPossibleMutsTree[1,alleles[ref[pos]],a]+=1
					countsPossibleMutsContextU[gene,alleles[ref[pos]],a]+=1
					#count contexts
					if alleles[ref[pos]]==1 and a==3:
						for ip in range(contextRange):
							countsContextPossibilities[0][ip][alleles[ref[pos-contextRange+ip]]]+=1
							countsContextPossibilities[0][ip+contextRange][alleles[ref[pos+(ip+1)]]]+=1
					elif alleles[ref[pos]]==2 and a==3:
						for ip in range(contextRange):
							countsContextPossibilities[1][ip][alleles[ref[pos-contextRange+ip]]]+=1
							countsContextPossibilities[1][ip+contextRange][alleles[ref[pos+(ip+1)]]]+=1
		if fourFoldCount==3:
			fourFolds[alleles[ref[pos]]]+=1
			for a in range(4):
				if a!=alleles[ref[pos]]:
					countsPossibleMutsContext[3,alleles[ref[pos-1]],alleles[ref[pos]],alleles[ref[pos+1]],a]+=1
					countsPossibleMutsTree[3,alleles[ref[pos]],a]+=1
print("Possible mutation events that would result into a stop codon")
print(intoStop)
#print("Numbers of opportunities for different mutation events (intergenic, synonymous, nonsynonymous; normal, C->U context)")		
#print(countsPossibleMuts)
print("Numbers of 4-fold degenerate sites")
print(fourFolds)
sum4=0.0
for a in range(4):
	sum4+=fourFolds[a]
for a in range(4):
	fourFolds[a]=fourFolds[a]/sum4
print("Frequencies at 4-fold degenerate sites")
print(fourFolds)
#print(countsContextPossibilities)
#print("New enocmpassing vectors for context and mutation possibilities:")
#print(countsPossibleMutsContext)



#plot a barplot for mutation rates and similar
def barplot(valuesLists,axisLabels,plotFileName,colors=[]):
	#plot the rates in a barplot
	labels = ['A>C', 'A>G', 'A>U', 'C>A', 'C>G','C>U','G>A','G>C','G>U','U>A','U>C','U>G']
	values=[]
	for k in range(len(axisLabels)):
		values.append([])
		for i in range(4):
			for j in range(4):
				if j!=i:
					values[k].append(valuesLists[k][i][j])
	x = np.arange(len(labels))
	#width = 0.35  # the width of the bars
	width=0.8/len(axisLabels)
	fig, ax = plt.subplots()
	reacts=[]
	for k in range(len(axisLabels)):
		if len(axisLabels)%2==0:
			if len(colors)==0:
				reacts.append(ax.bar(x +k*width - (len(axisLabels)-1)*width/2, values[k], width, label=axisLabels[k]))
			else:
				reacts.append(ax.bar(x +k*width - (len(axisLabels)-1)*width/2, values[k], width, color=colors[k], label=axisLabels[k]))
		else:
			if len(colors)==0:
				reacts.append(ax.bar(x +k*width -int(len(axisLabels)/2)*width, values[k], width, label=axisLabels[k]))
			else:
				reacts.append(ax.bar(x +k*width -int(len(axisLabels)/2)*width, values[k], width, color=colors[k], label=axisLabels[k]))
	ax.set_xlabel('Mutation types')
	ax.set_xticks(x)
	ax.set_xticklabels(labels)
	ax.legend()
	fig.tight_layout()
	#plt.show()
	fig.savefig(plotFileName)
	
#plot a barplot for mutation rates and similar when considering the near mutational context
def barplotContext(valuesLists,axisLabels,plotFileName,xName='Mutation context',noSave=False,colors=[],labels=['A_A', 'A_C', 'A_G', 'A_U', 'C_A', 'C_C', 'C_G','C_U','G_A','G_C', 'G_G','G_U','U_A','U_C','U_G', 'U_U']):
	#plot the rates in a barplot
	values=[]
	for k in range(len(axisLabels)):
		values.append([])
		if len(labels)==16:
			for i in range(4):
				for j in range(4):
				#if j!=i:
					values[k].append(valuesLists[k][i][j])
		else:
			for i in range(len(labels)):
				values[k].append(valuesLists[k][i])
	x = np.arange(len(labels))
	#width = 0.35  # the width of the bars
	width=0.8/len(axisLabels)
	fig, ax = plt.subplots()
	reacts=[]
	for k in range(len(axisLabels)):
		if len(axisLabels)%2==0:
			if len(colors)==0:
				reacts.append(ax.bar(x +k*width - (len(axisLabels)-1)*width/2, values[k], width, label=axisLabels[k]))
			else:
				reacts.append(ax.bar(x +k*width - (len(axisLabels)-1)*width/2, values[k], width, color=colors[k], label=axisLabels[k]))
		else:
			if len(colors)==0:
				reacts.append(ax.bar(x +k*width -int(len(axisLabels)/2)*width, values[k], width, label=axisLabels[k]))
			else:
				reacts.append(ax.bar(x +k*width -int(len(axisLabels)/2)*width, values[k], width, color=colors[k], label=axisLabels[k]))
	ax.set_xlabel(xName)
	ax.set_xticks(x)
	ax.set_xticklabels(labels)
	ax.legend()
	fig.tight_layout()
	#plt.show()
	if not noSave:
		fig.savefig(plotFileName)
	else:
		return fig

#plot a barplot for U selection across genes
def barplotGene(valuesLists,axisLabels,plotFileName,colors=[]):
	#plot the rates in a barplot
	labels = geneEndsNames
	values=[]
	for k in range(len(axisLabels)):
		values.append([])
		for g in range(len(geneEnds)):
			values[k].append(valuesLists[g][k])
	x = np.arange(len(labels))
	#width = 0.35  # the width of the bars
	width=0.8/len(axisLabels)
	fig, ax = plt.subplots()
	reacts=[]
	for k in range(len(axisLabels)):
		if len(axisLabels)%2==0:
			if len(colors)==0:
				reacts.append(ax.bar(x +k*width - (len(axisLabels)-1)*width/2, values[k], width, label=axisLabels[k]))
			else:
				reacts.append(ax.bar(x +k*width - (len(axisLabels)-1)*width/2, values[k], width, color=colors[k], label=axisLabels[k]))
		else:
			if len(colors)==0:
				reacts.append(ax.bar(x +k*width -int(len(axisLabels)/2)*width, values[k], width, label=axisLabels[k]))
			else:
				reacts.append(ax.bar(x +k*width -int(len(axisLabels)/2)*width, values[k], width, color=colors[k], label=axisLabels[k]))
	ax.set_xlabel('Gene')
	ax.set_xticks(x)
	ax.set_xticklabels(labels)
	ax.legend()
	fig.tight_layout()
	#plt.show()
	fig.savefig(plotFileName)


print("New estimation of mutation rates, no approximations")
file=open(pathSimu+"mutations_13_11_20.tsv")
dict_mutations={}
line=file.readline()
mutBranchSpectrum=np.zeros(20,dtype=int)
while line!="" and line!="\n":
	linelist=line.split()
	if len(linelist)>1:
		linelist2=linelist[1].split(",")
		mutBranchSpectrum[len(linelist2)]+=1
		mutList=[]
		for m in linelist2:
			mutList.append([m[0],m[-1],int(m[1:-1])])
		dict_mutations[linelist[0]]=mutList
	else:
		mutBranchSpectrum[0]+=1
	line=file.readline()
print("Branches with mutations:")
print(len(dict_mutations.keys()))
print(mutBranchSpectrum)
print(sum(mutBranchSpectrum)-(mutBranchSpectrum[0]+mutBranchSpectrum[1]))
meanPerBranch=0.0
for i in range(20):
	meanPerBranch+=mutBranchSpectrum[i]*i
meanPerBranch=meanPerBranch/sum(mutBranchSpectrum)
print("mean mutations per branch: "+str(meanPerBranch))

t = Tree(pathSimu+"gisaid-13-11-20.nwk",format=1)

#vector that says which codon mutations are synonymous and which are not (and which are nonsense)
isSynonymous=np.zeros((4,4,4,3,4),dtype=int)
for i1 in range(4):
	for i2 in range(4):
		for i3 in range(4):
			refCodonList=[allelesList[i1],allelesList[i2],allelesList[i3]]
			refCodon=allelesList[i1]+allelesList[i2]+allelesList[i3]
			refA=_translate_str(refCodon, table)
			for i in range(3):
				if i==0:
					index=i1
				elif i==1:
					index=i2
				else:
					index=i3
				altCodList=list(refCodon)
				for j in range(4):
					if j!=index:
						altCodList[i]=allelesList[j]
						altCod="".join(altCodList)
						altA=_translate_str(altCod, table)
						if refA==altA:
							isSynonymous[i1,i2,i3,i,j]=1
						elif altA=="*":
							isSynonymous[i1,i2,i3,i,j]=4
						elif refA!="*":
							isSynonymous[i1,i2,i3,i,j]=2
						#print(refCodon+" "+refA+" "+altCod+" "+altA+" "+str(isSynonymous[i1,i2,i3,i,j]))
#exit()
#print("synonymicity vector:")
#print(isSynonymous)
#print(isSynonymous[3,3,3,2,0])
#print(isSynonymous[3,3,3,2,1])
#print(isSynonymous[3,3,3,2,2])
#exit()

#new vectors taking into account immediate context (for CpG and APOBEC target).
# first dimension: type of mutation (noncoding, synonymous, nonsynonymous, 4-fold degenerate site, nonsense).
#dimensions 2-3: fromA, toA 
countsPossibleMuts_treewise=np.zeros((5,4,4),dtype=int)
#use this vector the record the difference between the baseline possibilities and those of the current ancestral genome
possibleMutsDifferences=np.zeros((5,4,4),dtype=int)
#dictionary of mutations occurred from root up to this branch
dict_current_mutations={}
allMutations=[]
#for i in range(len(ref)):
#	allMutations[i+1]={}
nMutations=[0]
totBLen=[0]
#function that traverses the tree recording information on each mutation event,
#like num of descendants, context of the mutation, and recording mutation possibilities.
def traverseTree(t,dict_current_mutations,possibleMutsDifferences,countsPossibleMuts_treewise,nMutations,totBLen):
	if t.name in dict_mutations:
		bLen=int(len(dict_mutations[t.name]))
	else:
		bLen=int(0)
	if bLen>0.001:
		totBLen[0]+=bLen
		removedMutations={}
		increment=np.zeros((5,4,4),dtype=int)
		#print(countsPossibleMuts_treewise)
		#print(possibleMutsDifferences)
		#print(bLen)
		#print(possibleMutsDifferences*bLen)
		countsPossibleMuts_treewise+=possibleMutsDifferences*bLen
		mutationAdded=False
		#print("Branch length "+str(t.dist))
		for m in dict_mutations[t.name]:
			#print(m)
			pos=m[2]-1
			if pos<100 or pos>29802:
				continue
			if (pos+1) in masked:
				continue
			p=pos
			gene=0
			while gene<len(geneEnds) and p>=geneEnds[gene][1]:
				gene+=1
			if gene>=len(geneEnds) or p<(geneEnds[gene][0]-1):
				typeOfPos=0
			if gene<len(geneEnds):
				if p>=(geneEnds[gene][0]-1) and p<geneEnds[gene][1]:
					if p<=(geneEnds[gene][0]+1):
						typeOfPos=3
					elif p>=(geneEnds[gene][1]-3):
						typeOfPos=4
					else:
						typeOfPos=2
			if typeOfPos==0:
				for a in range(4):
					if a!=alleles[m[0]]:
						increment[0,alleles[m[0]],a]-=1
					if a!=alleles[m[1]]:
						increment[0,alleles[m[1]],a]+=1
				allMutations.append([0,m[0],m[1],m[2],0])
				nMutations[0]+=1
				#print("Noncoding mutation")
				mutationAdded=True
			elif typeOfPos==2:
				rest=(pos+1-geneEnds[gene][0])%3
				cod1=pos-rest
				refCodonList=[ref[cod1],ref[cod1+1],ref[cod1+2]]
				for i in range(3):
					if (cod1+i+1) in dict_current_mutations:
						refCodonList[i]=dict_current_mutations[cod1+i+1][0]
				refCodon="".join(refCodonList)
				refA=_translate_str(refCodon, table)
				if refA=="*":
					#print("mutated stop codon "+refA+" "+refCodon+" "+str(rest)+" "+str(cod1))
					#print(m)
					#print([ref[cod1],ref[cod1+1],ref[cod1+2]])
					continue
				#exit()
				else:
					nMutations[0]+=1
					mutationAdded=True
					altCodList=list(refCodon)
					altCodList[rest]=m[1]
					altCod="".join(altCodList)
					altA=_translate_str(altCod, table)
					indeces=[alleles[refCodonList[0]],alleles[refCodonList[1]],alleles[refCodonList[2]]]
					if altA==refA:
						#print("Synonymous mutation")
						allMutations.append([1,m[0],m[1],m[2],0])
						for a in range(4):
							syn=isSynonymous[indeces[0],indeces[1],indeces[2],rest,a]
							if a!=alleles[m[0]]:
								increment[syn,alleles[m[0]],a]-=1
							if a!=alleles[m[1]]:
								increment[syn,alleles[m[1]],a]+=1

					elif altA=="*":
						#print("Nonsense mutation")
						allMutations.append([4,m[0],m[1],m[2],0])
						for i in range(3):
							for a in range(4):
								if a!=indeces[i]:
									syn=isSynonymous[indeces[0],indeces[1],indeces[2],i,a]
									increment[syn,indeces[i],a]-=1
					else:
						#print("NonSynonymous mutation")
						indeces2=[alleles[altCodList[0]],alleles[altCodList[1]],alleles[altCodList[2]]]
						allMutations.append([2,m[0],m[1],m[2],0])
						for i in range(3):
							for a in range(4):
								if a!=indeces[i]:
									syn=isSynonymous[indeces[0],indeces[1],indeces[2],i,a]
									increment[syn,indeces[i],a]-=1
								if a!=indeces2[i]:
									syn=isSynonymous[indeces2[0],indeces2[1],indeces2[2],i,a]
									increment[syn,indeces2[i],a]+=1

			if mutationAdded:
				#add mutations to the current mutation dictionary
				#print(nMutations)
				if m[2] in dict_current_mutations:
					if dict_current_mutations[m[2]][0]==m[0]:
						removedMutations[m[2]]=dict_current_mutations[m[2]]
						dict_current_mutations[m[2]]=[m[1],nMutations[0]-1]
					else:
						print("Error, at position ",str(m[2])," mutation dictionary entry is "+dict_current_mutations[m[2]]+" but the mutation is from "+m[0])
						exit()
				else:
					if m[0]==ref[m[2]-1]:
						dict_current_mutations[m[2]]=[m[1],nMutations[0]-1]
					else:
						print("Error, at position "+str(m[2])+" reference is "+ref[m[2]-1]+" but the mutation is from "+m[0])
						print(nMutations)
						#exit()

		possibleMutsDifferences+=increment
	if len(t.children)==0:
		for pos in dict_current_mutations.keys():
			allMutations[dict_current_mutations[pos][1]][-1]+=1
	for c in t.children:
		traverseTree(c,dict_current_mutations,possibleMutsDifferences,countsPossibleMuts_treewise,nMutations,totBLen)
	#reverse mutations to return to original state of the current mutation dictionary
	if bLen>0.001:
		possibleMutsDifferences-=increment
		for m in dict_mutations[t.name]:
			if m[2] in dict_current_mutations:
				del dict_current_mutations[m[2]]
		for pos in removedMutations.keys():
			dict_current_mutations[pos]=removedMutations[pos]

traverseTree(t,dict_current_mutations,possibleMutsDifferences,countsPossibleMuts_treewise,nMutations,totBLen)
#print(allMutations)
print("total mutations: "+str(nMutations))
print("total tree length: "+str(totBLen))

#print("mutation possibilities from reference: ")
#print(countsPossibleMutsTree)
countsPossibleMuts_treewise+=countsPossibleMutsTree*totBLen[0]
print("mutation possibilities from traversing the tree: ")
print(countsPossibleMuts_treewise)

print("mutation events from traversing the tree: ")
mutationsHappened=np.zeros((5,4,4,3))
for m in allMutations:
	#print(m)
	if m[4]>4:
		index=2
	elif m[4]>1:
		index=1
	else:
		index=0
	mutationsHappened[m[0],alleles[m[1]],alleles[m[2]],index]+=1
print(mutationsHappened)

mutationRatesTree=np.zeros((5,4,4,3))
for i in range(5):
	for i1 in range(4):
		for i2 in range(4):
			for fre in range(3):
				if countsPossibleMuts_treewise[i,i1,i2]>0:
					mutationRatesTree[i,i1,i2,fre] = mutationsHappened[i,i1,i2,fre]/countsPossibleMuts_treewise[i,i1,i2]
print(mutationRatesTree)

#exit()


#collect information regarding mutations
file=open(pathSimu+"find_parsimonious_assignments.13-11-20.out") #.collapsed
line=file.readline()
moreThan2=np.zeros((4,4))
mutations1000=0
mutCounts=[np.zeros((4,4)),[],np.zeros((4,4)),[],np.zeros((4,4)),[],np.zeros((4,4))]
siteCounts=[np.zeros((4,4)),[],np.zeros((4,4)),[],np.zeros((4,4)),[],np.zeros((4,4))]
mutRates=[np.zeros((4,4)),[],np.zeros((4,4)),[],np.zeros((4,4)),[],np.zeros((4,4))]
siteRates=[np.zeros((4,4)),[],np.zeros((4,4)),[],np.zeros((4,4)),[],np.zeros((4,4))]
mutationSpectrum=np.zeros((3,4,4,200),dtype=int)
mutationSpectrumRates=np.zeros((3,4,4,200))
numMutsCounts=[]
for k in range(3):
	numMutsCounts.append([])
	for i in range(4):
		numMutsCounts[k].append([])
		for j in range(4):
			numMutsCounts[k][i].append([])
#new vectors taking into account immediate context (for CpG and APOBEC target).
# first dimension: type of mutation (noncoding, synonymous, nonsynonymous, 4-fold degenerate sites, and nonsense mutations).
#dimensions 1-3: reference bases at positions -1, 0 and 1 respectively. 
#dimension 4: destination allele.
#dimension 5: variant allele count.
countsMutsContext=np.zeros((5,4,4,4,4,100))
countsMutsContextSites=np.zeros((5,4,4,4,4,100))
#new vectors for gene-wie U selection.
# first dimension: gene number
#dimensions 2: reference base. 
#dimension 3: destination allele.
#dimension 4: variant allele count.
countsMutsContextU=np.zeros((len(geneEnds),4,4,100))
#countsPossibleMutsContext=np.zeros((4,4,4,4,4))
logo=[]
countZeros=0
countZeros2=0
countStop=0
countStop2=0
countStop3=0
highHomoplMuts=[]
totMutsCounted=[0,0]
totDiffs=0
#print("\n Now plotting the most recurrent mutations of the genome, with particular emphasis on those that are not C->U or G->U")
while line!="" and line!="\n":
	line=file.readline()
	linelist=line.split()
	if len(linelist)==0:
		continue
	refA=linelist[0][0]
	linelist2=linelist[0].split(",")
	site=int(linelist2[0][1:-1])

	if site<100 or site>29802:
		continue
	if (site) in masked:
		continue

	p=site
	pos=site
	typeOfPos=4 #0=intergenic, 2=genic, 3=start codon, 4=stop codon
	gene=0
	while gene<len(geneEnds) and p>=geneEnds[gene][1]:
		gene+=1
	if gene>=len(geneEnds) or p<(geneEnds[gene][0]-1):
		typeOfPos=0
	if gene<len(geneEnds):
		if p>=(geneEnds[gene][0]-1) and p<geneEnds[gene][1]:
			if p<=(geneEnds[gene][0]+1):
				typeOfPos=3
			elif p>=(geneEnds[gene][1]-3):
				typeOfPos=4
			else:
				typeOfPos=2

	if refA!=ref[site-1]:
		print(line)
		print(ref[site-1])
		print(str(site))
		print("problem with reference?")
		exit()
	
	#create list of alternative alleles with counts
	altA=[]
	lineIndex=1
	for s in linelist2:
		a=s[-1]
		if a in allelesList:
			n=int(linelist[lineIndex].split("=")[1])
			a=linelist[lineIndex][0]
			if n>0:
				altA.append([a,n])
			lineIndex+=1
	#print(altA)

	#create list of mutation events
	found=False
	mutations=[]
	mutationsAtSite=np.zeros((4,4),dtype=int)
	for s in linelist:
		if s.find("mutation_clade_sizes") !=-1:
			Afrom=s[0]
			Ato=s[2]
			sizes=s.split("=")[1].split(",")
			for i in range(len(sizes)):
				sizes[i]=int(sizes[i])
			if len(sizes)>60:
				#print(linelist[0])
				#print(s)
				#print(linelist[0]+"\t"+str(len(sizes)))
				highHomoplMuts.append([linelist[0],len(sizes),Afrom,Ato,ref[site-11:site+10],sum(sizes)])
				#if ((not (Afrom=="C" and Ato=="T")) and (not (Afrom=="G" and Ato=="T"))) or len(sizes)>100:
				#	print(linelist[0])
				#	print(s)
			for i in range(len(sizes)):
				#sizes[i]=int(sizes[i])
				if sizes[i]>lowestCount:
					moreThan2[alleles[Afrom],alleles[Ato]]+=1
					mutationsAtSite[alleles[Afrom],alleles[Ato]]+=1
				if sizes[i]>1000:
					found=True
			mutations.append([Afrom,Ato,sizes])
			if Ato!=ref[pos-1]:
				for s in sizes:
					totDiffs+=s
			if Afrom==ref[pos-1]:
				totMutsCounted[0]+=len(sizes)
			else:
				totMutsCounted[1]+=len(sizes)
	#if np.sum(mutationsAtSite)>2:
	#	print(linelist)
	#	print(ref[pos-1])
	#	print(mutationsAtSite)
	#	print("\n")
	if found:
		#print(linelist[0])
		#print(mutations)
		#print("\n")
		mutations1000+=1
	
	if alleles[ref[pos-1]]==1:
		countContext=0
		if alleles[ref[pos-2]]==0 or alleles[ref[pos-2]]==3:
			countContext+=1
		if alleles[ref[pos]]==0 or alleles[ref[pos]]==3:
			countContext+=1

	if 	typeOfPos==0:
		for j in range(4):
			if j!=alleles[ref[pos-1]]:
				numMutsCounts[0][alleles[ref[pos-1]]][j].append(mutationsAtSite[alleles[ref[pos-1]]][j])
				mutationSpectrum[0][alleles[ref[pos-1]]][j][mutationsAtSite[alleles[ref[pos-1]]][j]]+=1
		
		for m in mutations:
			if m[0]==ref[pos-1]:
				for s in m[2]:
					if s>=100:
						countsMutsContext[0,alleles[ref[pos-2]],alleles[ref[pos-1]],alleles[ref[pos]],alleles[m[1]],99]+=1
					else:
						countsMutsContext[0,alleles[ref[pos-2]],alleles[ref[pos-1]],alleles[ref[pos]],alleles[m[1]],s-1]+=1
					if s>lowestCount:
					#if m[0]==ref[pos-1]:
						mutCounts[0][alleles[m[0]]][alleles[m[1]]]+=1
						#if m[0]=="C":
						#	mutCounts[1][countContext][alleles[m[1]]]+=1
				
		for s in altA:
			if s[1]>=100:
				countsMutsContextSites[0,alleles[ref[pos-2]],alleles[ref[pos-1]],alleles[ref[pos]],alleles[s[0]],99]+=1
			else:
				countsMutsContextSites[0,alleles[ref[pos-2]],alleles[ref[pos-1]],alleles[ref[pos]],alleles[s[0]],s[1]-1]+=1
			if s[1]>lowestCount:
				siteCounts[0][alleles[refA]][alleles[s[0]]]+=1
				#if refA=="C":
				#	siteCounts[1][countContext][alleles[s[0]]]+=1
				if refA=="G" and s[0]=="T":
					logo.append(ref[pos-6:pos+5])
	elif typeOfPos==2:
		rest=(pos-geneEnds[gene][0])%3
		cod1=pos-1-rest
		refCodon=ref[cod1:cod1+3]
		refAA=_translate_str(refCodon, table)
		if refAA=="*":
			print("error")
			exit()
		
		fourFoldCount=0
		for j in range(4):
			if j!=alleles[ref[pos-1]]:
				altCod=list(refCodon)
				altCod[rest]=allelesList[j]
				altCod="".join(altCod)
				altAA=_translate_str(altCod, table)
				if mutationsAtSite[alleles[ref[pos-1]]][j]>60:
					#print(refAA+" "+refCodon+" -> "+altAA+" "+altCod)
					highHomoplMuts[-1].append(refAA+" "+refCodon+" -> "+altAA+" "+altCod)
				if altAA!=refAA and altAA!="*":
					numMutsCounts[2][alleles[ref[pos-1]]][j].append(mutationsAtSite[alleles[ref[pos-1]]][j])
					mutationSpectrum[2][alleles[ref[pos-1]]][j][mutationsAtSite[alleles[ref[pos-1]]][j]]+=1
				elif altAA==refAA:
					numMutsCounts[1][alleles[ref[pos-1]]][j].append(mutationsAtSite[alleles[ref[pos-1]]][j])
					mutationSpectrum[1][alleles[ref[pos-1]]][j][mutationsAtSite[alleles[ref[pos-1]]][j]]+=1
					fourFoldCount+=1
				elif altAA=="*" and mutationsAtSite[alleles[ref[pos-1]]][j]>0:
					countStop3+=1

			
		for m in mutations:
			altCod=list(refCodon)
			altCod[rest]=m[1]
			altCod="".join(altCod)
			altAA=_translate_str(altCod, table)
			if m[0]==ref[pos-1]:
				for s in m[2]:
						if altAA!=refAA and altAA!="*":
							if s>=100:
								countsMutsContext[2,alleles[ref[pos-2]],alleles[ref[pos-1]],alleles[ref[pos]],alleles[m[1]],99]+=1
							else:
								countsMutsContext[2,alleles[ref[pos-2]],alleles[ref[pos-1]],alleles[ref[pos]],alleles[m[1]],s-1]+=1
							if s>lowestCount:
								mutCounts[4][alleles[m[0]]][alleles[m[1]]]+=1
								#if m[0]=="C":
								#	mutCounts[5][countContext][alleles[m[1]]]+=1
						elif altAA==refAA:
							if s>=100:
								countsMutsContext[1,alleles[ref[pos-2]],alleles[ref[pos-1]],alleles[ref[pos]],alleles[m[1]],99]+=1
								countsMutsContextU[gene,alleles[ref[pos-1]],alleles[m[1]],99]+=1
								if fourFoldCount==3:
									countsMutsContext[3,alleles[ref[pos-2]],alleles[ref[pos-1]],alleles[ref[pos]],alleles[m[1]],99]+=1
							else:
								countsMutsContext[1,alleles[ref[pos-2]],alleles[ref[pos-1]],alleles[ref[pos]],alleles[m[1]],s-1]+=1
								countsMutsContextU[gene,alleles[ref[pos-1]],alleles[m[1]],s-1]+=1
								if fourFoldCount==3:
									countsMutsContext[3,alleles[ref[pos-2]],alleles[ref[pos-1]],alleles[ref[pos]],alleles[m[1]],s-1]+=1
							if s>lowestCount:
								mutCounts[2][alleles[m[0]]][alleles[m[1]]]+=1
								if m[0]=="C":
									#mutCounts[3][countContext][alleles[m[1]]]+=1
									if m[1]=="T":
										for ip in range(contextRange):
											countsContext[0][ip][alleles[ref[pos-contextRange+ip-1]]]+=1
											countsContext[0][ip+contextRange][alleles[ref[pos+ip]]]+=1
								if m[0]=="G" and m[1]=="T":
										for ip in range(contextRange):
											countsContext[1][ip][alleles[ref[pos-contextRange+ip-1]]]+=1
											countsContext[1][ip+contextRange][alleles[ref[pos+ip]]]+=1
								if fourFoldCount==3:
									fourFoldMuts[0][alleles[m[0]]][alleles[m[1]]]+=1
							
						else:
							if s>=100:
								countsMutsContext[4,alleles[ref[pos-2]],alleles[ref[pos-1]],alleles[ref[pos]],alleles[m[1]],99]+=1
							else:
								countsMutsContext[4,alleles[ref[pos-2]],alleles[ref[pos-1]],alleles[ref[pos]],alleles[m[1]],s-1]+=1
							if s>lowestCount:
								mutCounts[6][alleles[m[0]]][alleles[m[1]]]+=1
								countStop2+=1
								#print(pos)
								#print(linelist[0])
		for s in altA:
				altCod=list(refCodon)
				altCod[rest]=s[0]
				altCod="".join(altCod)
				altAA=_translate_str(altCod, table)
				#if s[1]>lowestCount:
				if altAA!=refAA and altAA!="*":
					if s[1]>=100:
						countsMutsContextSites[2,alleles[ref[pos-2]],alleles[ref[pos-1]],alleles[ref[pos]],alleles[s[0]],99]+=1
					else:
						countsMutsContextSites[2,alleles[ref[pos-2]],alleles[ref[pos-1]],alleles[ref[pos]],alleles[s[0]],s[1]-1]+=1
					if s[1]>lowestCount:
						siteCounts[4][alleles[refA]][alleles[s[0]]]+=1
						#if refA=="C":
						#	siteCounts[5][countContext][alleles[s[0]]]+=1
				elif altAA==refAA:
					if s[1]>=100:
						countsMutsContextSites[1,alleles[ref[pos-2]],alleles[ref[pos-1]],alleles[ref[pos]],alleles[s[0]],99]+=1
						if fourFoldCount==3:
							countsMutsContextSites[3,alleles[ref[pos-2]],alleles[ref[pos-1]],alleles[ref[pos]],alleles[s[0]],99]+=1
					else:
						countsMutsContextSites[1,alleles[ref[pos-2]],alleles[ref[pos-1]],alleles[ref[pos]],alleles[s[0]],s[1]-1]+=1
						if fourFoldCount==3:
							countsMutsContextSites[3,alleles[ref[pos-2]],alleles[ref[pos-1]],alleles[ref[pos]],alleles[s[0]],s[1]-1]+=1
					if s[1]>lowestCount:
						siteCounts[2][alleles[refA]][alleles[s[0]]]+=1
						if refA=="C":
							#siteCounts[3][countContext][alleles[s[0]]]+=1
							if s[0]=="T":
								for ip in range(contextRange):
									countsContextSites[0][ip][alleles[ref[pos-contextRange+ip-1]]]+=1
									countsContextSites[0][ip+contextRange][alleles[ref[pos+ip]]]+=1
						if refA=="G" and s[0]=="T":
								for ip in range(contextRange):
									countsContextSites[1][ip][alleles[ref[pos-contextRange+ip-1]]]+=1
									countsContextSites[1][ip+contextRange][alleles[ref[pos+ip]]]+=1
						if fourFoldCount==3:
							fourFoldMuts[1][alleles[refA]][alleles[s[0]]]+=1
				else:
					if s[1]>=100:
						countsMutsContextSites[4,alleles[ref[pos-2]],alleles[ref[pos-1]],alleles[ref[pos]],alleles[s[0]],99]+=1
					else:
						countsMutsContextSites[4,alleles[ref[pos-2]],alleles[ref[pos-1]],alleles[ref[pos]],alleles[s[0]],s[1]-1]+=1
					if s[1]>lowestCount:
						siteCounts[6][alleles[refA]][alleles[s[0]]]+=1
						countStop+=1
				if refA=="G" and s[0]=="T":
					if s[1]>lowestCount:
						logo.append(ref[pos-6:pos+5])
file.close()
						
for i in range(len(mutCounts)):
	for j in range(len(mutCounts[i])):
		for k in range(len(mutCounts[i][j])):
			mutRates[i][j][k]=mutCounts[i][j][k]/countsPossibleMuts[i][j][k]
			siteRates[i][j][k]=siteCounts[i][j][k]/countsPossibleMuts[i][j][k]
			
#sort most homoplasic mutations by number of mutation events
highHomoplMutsSorted=highHomoplMuts.sort(key=lambda x: x[1], reverse=True)
print("Most homoplasic mutations:")
for si in highHomoplMuts:
	#print(si)
	if (si[3]!="T" or (si[2]!="C" and si[2]!="G")):
		print("non C->T or G->T")
	print(si)
print(len(highHomoplMuts))

print("\n"+"Counted and not counted mutations:")
print(totMutsCounted)
print("\n"+"Total differences wrt the reference:")
print(totDiffs)
print("Per genome:")
print(float(totDiffs)/147137)
pM=(float(totDiffs)/147137)/len(ref)
pM=12.0/len(ref)
pconstant=(1.0-pM)*(1.0-pM)
pchange=1.0-pconstant
print(pconstant,pchange)
#exit()

#print("Mutations of different types and rates:")
#print("Possibilities:")
#print(countsPossibleMuts)
#print("Counts of mutation events:")
#print(mutCounts)
#print("Rates from numbers of mutation events:")
#print(mutRates)
#print("Counts of variable sites:")
#print(siteCounts)
#print("Rates from numbers of variable sites:")
#print(siteRates)

#print("4-fold degenerate site mutations, using mutation event counts:  ")
#print(fourFoldMuts[0])
#print("Using variant site counts:")
#print(fourFoldMuts[1])
for k in range(2):
	total=0.0
	for i in range(4):
		for j in range(4):
			fourFoldMuts[k][i][j]=fourFoldMuts[k][i][j]/float(fourFolds[i])
			total+=fourFoldMuts[k][i][j]*float(fourFolds[i])
	for i in range(4):
		for j in range(4):
			fourFoldMuts[k][i][j]=fourFoldMuts[k][i][j]/total
#print("Normalized 4-fold degenrate site rates, using mutation event counts: ")
#print(fourFoldMuts[0])
#print("Using variant site counts:")
#print(fourFoldMuts[1])
#barplot([fourFoldMuts[0],fourFoldMuts[1]],['mutation events','variant alleles'],pathSimu+"barplot_rates4fold"+str(lowestCount+1)+".pdf")

#print("recurring mutations:")
maxSpectrum=0
for k in range(3):
	for i in range(4):
		for j in range(4):
			if i!=j:
				mutationSpectrum[k][i][j][0]+=(countsPossibleMuts[k*2][i][j]-np.sum(mutationSpectrum[k][i][j]))
			for l in range(100):
				if mutationSpectrum[k][i][j][l]>1 and l>maxSpectrum:
					maxSpectrum=l
				if i!=j:
					mutationSpectrumRates[k][i][j][l]=mutationSpectrum[k][i][j][l]/countsPossibleMuts[k*2][i][j]
			#print(str(k)+" "+str(i)+" "+str(j)+" ")
			#print(countsPossibleMuts[k*2][i][j])
			#print(np.sum(mutationSpectrum[k][i][j]))
					
countsContextRelative=np.zeros((2,2*contextRange,4))
countsContextRelativeSites=np.zeros((2,2*contextRange,4))
#print("context of CU and GU mutations")
#print(countsContextPossibilities)
#print(countsContext)
#print(countsContextSites)
for k in range(2):
	for ip in range(2*contextRange):
		for j in range(4):
			countsContextRates[k][ip][j]=countsContext[k][ip][j]/countsContextPossibilities[k][ip][j]
			countsContextRelative[k][ip][j]=(countsContextRates[k][ip][j]-mutRates[2][1+k][3])/mutRates[2][1+k][3]
			countsContextRatesSites[k][ip][j]=countsContextSites[k][ip][j]/countsContextPossibilities[k][ip][j]
			countsContextRelativeSites[k][ip][j]=(countsContextRatesSites[k][ip][j]-siteRates[2][1+k][3])/siteRates[2][1+k][3]
#print(countsContextRates)
#print(countsContextRatesSites)
#print("Context rates relative to baseline")
#print(countsContextRelative)
#print(countsContextRelativeSites)
print("Now printing context plots for context up to 5bp upstream and downstream.")
allelesListRNA=["A","C","G","U"]
#barplot for context
namesMut=["C->U","G->U"]
cm = plt.get_cmap('tab20')
for k in range(2):
	#print("k="+str(k)+namesMut[k])
	NUM_COLORS=4
	space=0.17
	
	fig, ax = plt.subplots()
	#print([cm(1.*i/NUM_COLORS) for i in range(NUM_COLORS)])
	ax.set_prop_cycle('color',[cm(1.*i/NUM_COLORS) for i in range(NUM_COLORS)])
	reacts=[]
	#for i in range(contextRange):
	xRangeArray=np.concatenate((np.arange(contextRange)-contextRange,np.arange(contextRange) +1))
	#print(xRangeArray)
	for j in range(4):
		step=(j-2+0.5)*space
		reacts.append(ax.bar(xRangeArray +step, countsContextSites[k,:,j], space, label=allelesListRNA[j]))
		#reacts.append(ax.bar(np.arange(contextRange) +step+1, countsContextRelative[k,contextRange:,j], space, label=allelesList[j]))
	ax.legend()
	fig.savefig(pathSimu+"histogram_context_countsSites_"+namesMut[k]+str(lowestCount)+".pdf", scale=16)
	
	fig, ax = plt.subplots()
	ax.set_prop_cycle('color',[cm(1.*i/NUM_COLORS) for i in range(NUM_COLORS)])
	reacts=[]
	#for i in range(contextRange):
	xRangeArray=np.concatenate((np.arange(contextRange)-contextRange,np.arange(contextRange) +1))
	#print(xRangeArray)
	for j in range(4):
		step=(j-2+0.5)*space
		reacts.append(ax.bar(xRangeArray +step, countsContext[k,:,j], space, label=allelesListRNA[j]))
		#reacts.append(ax.bar(np.arange(contextRange) +step+1, countsContextRelative[k,contextRange:,j], space, label=allelesList[j]))
	ax.legend()
	fig.savefig(pathSimu+"histogram_context_counts_"+namesMut[k]+str(lowestCount)+".pdf", scale=16)
	
	fig, ax = plt.subplots()
	ax.set_prop_cycle('color',[cm(1.*i/NUM_COLORS) for i in range(NUM_COLORS)])
	#ax.set_color_cycle([cm(1.*i/NUM_COLORS) for i in range(NUM_COLORS)])
	reacts=[]
	#for i in range(contextRange):
	xRangeArray=np.concatenate((np.arange(contextRange)-contextRange,np.arange(contextRange) +1))
	#print(xRangeArray)
	for j in range(4):
		step=(j-2+0.5)*space
		reacts.append(ax.bar(xRangeArray +step, countsContextRelative[k,:,j], space, label=allelesListRNA[j]))
		#reacts.append(ax.bar(np.arange(contextRange) +step+1, countsContextRelative[k,contextRange:,j], space, label=allelesList[j]))
	ax.legend()
	fig.savefig(pathSimu+"histogram_context_"+namesMut[k]+str(lowestCount)+".pdf", scale=16)
	
	fig, ax = plt.subplots()
	ax.set_prop_cycle('color',[cm(1.*i/NUM_COLORS) for i in range(NUM_COLORS)])
	reacts=[]
	#for i in range(contextRange):
	for j in range(4):
		step=(j-2+0.5)*space
		reacts.append(ax.bar(xRangeArray +step, countsContextPossibilities[k,:,j], space, label=allelesListRNA[j]))
		#reacts.append(ax.bar(np.arange(contextRange) +step+1, countsContextPossibilities[k,contextRange:,j], space, label=allelesList[j]))
	ax.legend()
	fig.savefig(pathSimu+"histogram_context_possibilities_"+namesMut[k]+str(lowestCount)+".pdf", scale=16)
	
	fig, ax = plt.subplots()
	ax.set_prop_cycle('color',[cm(1.*i/NUM_COLORS) for i in range(NUM_COLORS)])
	reacts=[]
	#for i in range(contextRange):
	for j in range(4):
		step=(j-2+0.5)*space
		reacts.append(ax.bar(xRangeArray +step, countsContextRelativeSites[k,:,j], space, label=allelesListRNA[j]))
		#reacts.append(ax.bar(np.arange(contextRange) +step+1, countsContextRatesSites[k,contextRange:,j], space, label=allelesList[j]))
	ax.legend()
	fig.savefig(pathSimu+"histogram_context_sites_"+namesMut[k]+str(lowestCount)+".pdf", scale=16)
					
#maxSpectrum=int(np.max(mutationSpectrum))
namesMut=["noncoding","synonymous","nonsyn"]
cm = plt.get_cmap('tab20')
for k in range(3):
	#print("k="+str(k)+namesMut[k])
	NUM_COLORS=12
	fig, ax = plt.subplots()
	ax.set_prop_cycle('color',[cm(1.*i/NUM_COLORS) for i in range(NUM_COLORS)])
	reacts=[]
	space=0.07
	for i in range(4):
		for j in range(4):
			if i!=j:
				#if k==1 or k==2:
				#	print(allelesList[i]+allelesList[j])
				#	print(mutationSpectrum[k][i][j][:maxSpectrum+1])
				#print(np.mean(numMutsCounts[k][i][j]))
				if i>j:
					step=i*3*space+j*space - 5*space
				else:
					step=i*3*space+(j-1)*space - 5*space
				reacts.append(ax.bar(np.arange(maxSpectrum) +step+1, mutationSpectrum[k][i][j][1:maxSpectrum+1], space, label=allelesListRNA[i]+"->"+allelesListRNA[j]))
	ax.legend()
	fig.savefig(pathSimu+"histogram_"+namesMut[k]+str(lowestCount)+".pdf", scale=16)
	
	#print("Rates")
	NUM_COLORS=12
	fig, ax = plt.subplots()
	ax.set_prop_cycle('color',[cm(1.*i/NUM_COLORS) for i in range(NUM_COLORS)])
	reacts=[]
	space=0.07
	for i in range(4):
		for j in range(4):
			if i!=j:
				#if k==1 or k==2:
				#	print(allelesList[i]+allelesList[j])
				#	print(mutationSpectrumRates[k][i][j][:maxSpectrum+1])
				#print(np.mean(numMutsCounts[k][i][j]))
				if i>j:
					step=i*3*space+j*space - 5*space
				else:
					step=i*3*space+(j-1)*space - 5*space
				reacts.append(ax.bar(np.arange(maxSpectrum) +step+1, mutationSpectrumRates[k][i][j][1:maxSpectrum+1], space, label=allelesListRNA[i]+"->"+allelesListRNA[j]))
	ax.legend()
	fig.savefig(pathSimu+"histogram_rates_"+namesMut[k]+str(lowestCount)+".pdf", scale=16)
	plt.close()
	
	maxPrinted=5
	#print("Rates")
	NUM_COLORS=12
	fig, ax = plt.subplots()
	ax.set_prop_cycle('color',[cm(1.*i/NUM_COLORS) for i in range(NUM_COLORS)])
	xnames=[]
	for k2 in range(maxPrinted+1):
		xnames.append(str(k2))
	xnames.append(">"+str(maxPrinted))
	plt.xticks(range(maxPrinted+2),xnames)
	#print(range(maxPrinted+2))
	#print(xnames)
	reacts=[]
	space=0.05
	for i in range(4):
		for j in range(4):
			if i!=j:
				#if k==1 or k==2:
				#	print(allelesList[i]+allelesList[j])
				#	print(mutationSpectrumRates[k][i][j][:maxSpectrum+1])
				#print(np.mean(numMutsCounts[k][i][j]))
				if i>j:
					step=i*3*space+j*space - 5*space
				else:
					step=i*3*space+(j-1)*space - 5*space
				newList=[]
				for k2 in range(maxPrinted+1):
					newList.append(mutationSpectrumRates[k][i][j][k2])
				sumBig=np.sum(mutationSpectrumRates[k][i][j][maxPrinted+1:])
				newList.append(sumBig)
				#print(newList)
				reacts.append(ax.bar(np.arange(maxPrinted+2) +step, newList, space, label=allelesListRNA[i]+"->"+allelesListRNA[j]))
	ax.legend()
	fig.savefig(pathSimu+"histogram_rates_new_"+namesMut[k]+str(lowestCount)+".pdf", scale=16)
	plt.close()

#print("Summing non-coding and non-synonymous")
for k in range(2):
	tot=0.0
	totSites=0.0
	for i in range(len(mutCounts[k])):
		for j in range(len(mutCounts[k][i])):
			mutCounts[k][i][j]+=mutCounts[k+2][i][j]
			siteCounts[k][i][j]+=siteCounts[k+2][i][j]
			countsPossibleMuts[k][i][j]+=countsPossibleMuts[k+2][i][j]
			mutRates[k][i][j]=mutCounts[k][i][j]/countsPossibleMuts[k][i][j]
			siteRates[k][i][j]=siteCounts[k][i][j]/countsPossibleMuts[k][i][j]
			if not math.isnan(mutRates[k][i][j]):
				tot+=mutRates[k][i][j]
			if not math.isnan(siteRates[k][i][j]):
				totSites+=siteRates[k][i][j]
	#print(mutRates[k])
			
#barplot([mutCounts[0],countsPossibleMuts[0]],['mutation events','possibility counts'],pathSimu+"barplot_counts"+str(lowestCount+1)+".pdf")

file=open(pathSimu+"mutRates"+str(lowestCount+1)+".txt","w")
file.write("mutRates\n")
for i in range(len(mutRates)):
	for j in range(len(mutRates[i])):
		for k in range(len(mutRates[i][j])):
			file.write(str(mutRates[i][j][k])+"\t")
		file.write("\n")
	file.write("\n")
file.write("\n")
file.write("siteRates\n")
for i in range(len(mutRates)):
	for j in range(len(mutRates[i])):
		for k in range(len(mutRates[i][j])):
			file.write(str(siteRates[i][j][k])+"\t")
		file.write("\n")
	file.write("\n")
	
#barplot([mutRates[0],siteRates[0]],['mutation events','site counts'],pathSimu+"barplot"+str(lowestCount+1)+".pdf")			

#print("\n Rates of synonymous-noncoding mutations combined")
#print(mutRates[0])
#print(siteRates[0])
#print(mutRates[1])
#print(siteRates[1])
siteRatesNorm=np.copy(mutRates[0])
siteRatesNorm2=np.copy(siteRates[0])
counts=[ref.count("A"),ref.count("C"),ref.count("G"),ref.count("T")]
tot=sum(counts)
for i in range(4):
	counts[i]=float(counts[i])/tot
#print("Reference frequencies:")
#print(counts)
totSum=0.0
totSum2=0.0
for i in range(4):
	for j in range(4):
		if i!=j:
			totSum+=siteRatesNorm[i][j]*counts[i]
			totSum2+=siteRatesNorm2[i][j]*counts[i]
for i in range(4):
	for j in range(4):
		if i!=j:
			siteRatesNorm[i][j]=siteRatesNorm[i][j]/totSum
			siteRatesNorm2[i][j]=siteRatesNorm2[i][j]/totSum2
#print("\n Rates normalized for phylo inference")
#print(siteRatesNorm)
#print("\n Site Rates normalized for phylo inference")
#print(siteRatesNorm2)

omega=0.0
omegaPoss=0.0
omegaSites=0.0
for i in range(len(mutCounts[4])):
	for j in range(len(mutCounts[4][i])):
			omega+=mutCounts[4][i][j]
			omegaSites+=siteCounts[4][i][j]
			omegaPoss+=countsPossibleMuts[4][i][j]
omegaS=0.0
omegaPossS=0.0
omegaSitesS=0.0
for i in range(len(mutCounts[2])):
	for j in range(len(mutCounts[2][i])):
			omegaS+=mutCounts[2][i][j]
			omegaSitesS+=siteCounts[2][i][j]
			omegaPossS+=countsPossibleMuts[2][i][j]
omega=(omega*omegaPossS)/(omegaPoss*omegaS)
omegaSites=(omegaSites*omegaPossS)/(omegaPoss*omegaSitesS)
#print("Omegas:")
#print(omega)
#print(omegaSites)

file.write("omega\n")
file.write(str(omega)+"\n")
file.write("omegaSites\n")
file.write(str(omegaSites)+"\n")

print("Nonsense mutation counts:")
print(countStop3)
#print(countStop2)
#print(countStop)

#print("new  vectors:")
#print(countsPossibleMutsContext)
#print(countsMutsContext)
#print(countsMutsContextSites)

basicCounts=np.zeros((5,4,4,3),dtype=int)
basicCountsSites=np.zeros((5,4,4,3),dtype=int)
basicCountsOpportunities=np.zeros((5,4,4))
basicRates=np.zeros((5,4,4,3))
basicRatesSites=np.zeros((5,4,4,3))
for k in range(5):
	for i1 in range(4):
		for i2 in range(4):
			for i3 in range(4):
				for j in range(4):
					basicCountsOpportunities[k,i2,j]+=countsPossibleMutsContext[k,i1,i2,i3,j]
					for s in range(100):
						if s==0:
							basicCounts[k,i2,j,0]+=countsMutsContext[k,i1,i2,i3,j,s]
							basicCountsSites[k,i2,j,0]+=countsMutsContextSites[k,i1,i2,i3,j,s]
						elif s>0 and s<4:
							basicCounts[k,i2,j,1]+=countsMutsContext[k,i1,i2,i3,j,s]
							basicCountsSites[k,i2,j,1]+=countsMutsContextSites[k,i1,i2,i3,j,s]
						else:
							basicCounts[k,i2,j,2]+=countsMutsContext[k,i1,i2,i3,j,s]
							basicCountsSites[k,i2,j,2]+=countsMutsContextSites[k,i1,i2,i3,j,s]

for k in range(4):
	for i2 in range(4):
		for j in range(4):
			for s in range(3):
				basicRates[k,i2,j,s]=basicCounts[k,i2,j,s]/basicCountsOpportunities[k,i2,j]
				basicRatesSites[k,i2,j,s]=basicCountsSites[k,i2,j,s]/basicCountsOpportunities[k,i2,j]


normRates1=np.copy(basicRates[3,:,:,0])
normRates2=np.copy(basicRates[3,:,:,1])
normRates3=np.copy(basicRates[3,:,:,2])
normRates4=np.copy(basicRates[3,:,:,2])
normRates4=normRates4+basicRates[3,:,:,1]
tot1=0.0
tot2=0.0
tot3=0.0
tot4=0.0
freqs=np.zeros(4)
for aa in range(4):
	bb=(aa+1)%4
	freqs[aa]=basicCountsOpportunities[3,aa,bb]
freqs=freqs/np.sum(freqs)
print("4-fold degenerate nucleotide freqs")
print(freqs)
for aa in range(4):
	for bb in range(4):
		if bb!=aa:
			tot1+=basicRates[3,aa,bb,0]*freqs[aa]
			tot2+=basicRates[3,aa,bb,1]*freqs[aa]
			tot3+=basicRates[3,aa,bb,2]*freqs[aa]
			tot4+=normRates4[aa,bb]*freqs[aa]
#print(tot1)
#print(normRates1)
normRates1=normRates1/tot1
normRates2=normRates2/tot2
normRates3=normRates3/tot3
normRates4=normRates4/tot4
#print("4-fold degenerate rates normalized, single low and high freqs")
#print(normRates1)
#print(normRates2)
#print(normRates3)
#print("4-fold degenerate rates normalized, low+high freqs")
#print(normRates4)

#countsMutsContext=np.zeros((4,4,4,4,4,100))
#countsMutsContextSites=np.zeros((4,4,4,4,4,100))
#countsPossibleMutsContext=np.zeros((4,4,4,4,4))

#Plot one barplot with multiple (for different thresholds) counts and and one with rates
names=["nonCoding","synonymous","nonSynonymous","4fold","nonsense"]
colors=["green","red","orange","yellow","darkblue","blue","lightblue"]
colors2=["green","darkblue","blue","lightblue"]
for k in range(4):
	barplot([basicCountsOpportunities[k,:,:],basicCounts[k,:,:,0],basicCounts[k,:,:,1],basicCounts[k,:,:,2],basicCountsSites[k,:,:,0]+basicCountsSites[k,:,:,1]+basicCountsSites[k,:,:,2],basicCountsSites[k,:,:,1]+basicCountsSites[k,:,:,2],basicCountsSites[k,:,:,2]],['possible mut sites','1-descendant muts','2-4 descendants muts','>4 descendants muts','variants observed','non-singleton variants','frequent variants'],pathSimu+"barplot_counts_all_"+names[k]+".pdf",colors=colors)
	barplot([basicCountsOpportunities[k,:,:],basicCounts[k,:,:,0],basicCounts[k,:,:,1],basicCounts[k,:,:,2]],['possible mut sites','1-descendant muts','2-4 descendants muts','>4 descendants muts'],pathSimu+"barplot_counts_mutations_"+names[k]+".pdf",colors=colors)
	barplot([basicCountsOpportunities[k,:,:],basicCountsSites[k,:,:,0]+basicCountsSites[k,:,:,1]+basicCountsSites[k,:,:,2],basicCountsSites[k,:,:,1]+basicCountsSites[k,:,:,2],basicCountsSites[k,:,:,2]],['possible mut sites','sites with vars','sites with non-singleton vars','sites with frequent vars'],pathSimu+"barplot_counts_sites_"+names[k]+".pdf",colors=colors2)
	plt.close()
barplot([basicCountsOpportunities[0,:,:]+basicCountsOpportunities[1,:,:]+basicCountsOpportunities[2,:,:],basicCounts[0,:,:,0]+basicCounts[1,:,:,0]+basicCounts[2,:,:,0],basicCounts[0,:,:,1]+basicCounts[1,:,:,1]+basicCounts[2,:,:,1],basicCounts[0,:,:,2]+basicCounts[1,:,:,2]+basicCounts[2,:,:,2],basicCountsSites[0,:,:,0]+basicCountsSites[0,:,:,1]+basicCountsSites[0,:,:,2]+basicCountsSites[1,:,:,0]+basicCountsSites[1,:,:,1]+basicCountsSites[1,:,:,2]+basicCountsSites[2,:,:,0]+basicCountsSites[2,:,:,1]+basicCountsSites[2,:,:,2],basicCountsSites[0,:,:,1]+basicCountsSites[0,:,:,2]+basicCountsSites[1,:,:,1]+basicCountsSites[1,:,:,2]+basicCountsSites[2,:,:,1]+basicCountsSites[2,:,:,2],basicCountsSites[0,:,:,2]+basicCountsSites[1,:,:,2]+basicCountsSites[2,:,:,2]],['possible mut sites','1-descendant muts','2-4 descendants muts','>4 descendants muts','variants observed','non-singleton variants','frequent variants'],pathSimu+"barplot_counts_all_wholeGenome.pdf",colors=colors)
barplot([basicCountsOpportunities[0,:,:]+basicCountsOpportunities[1,:,:]+basicCountsOpportunities[2,:,:],basicCounts[0,:,:,0]+basicCounts[1,:,:,0]+basicCounts[2,:,:,0],basicCounts[0,:,:,1]+basicCounts[1,:,:,1]+basicCounts[2,:,:,1],basicCounts[0,:,:,2]+basicCounts[1,:,:,2]+basicCounts[2,:,:,2]],['possible mut sites','1-descendant muts','2-4 descendants muts','>4 descendants muts'],pathSimu+"barplot_counts_mutations_wholeGenome.pdf",colors=colors)
barplot([basicCountsOpportunities[0,:,:]+basicCountsOpportunities[1,:,:]+basicCountsOpportunities[2,:,:],basicCountsSites[0,:,:,0]+basicCountsSites[0,:,:,1]+basicCountsSites[0,:,:,2]+basicCountsSites[1,:,:,0]+basicCountsSites[1,:,:,1]+basicCountsSites[1,:,:,2]+basicCountsSites[2,:,:,0]+basicCountsSites[2,:,:,1]+basicCountsSites[2,:,:,2],basicCountsSites[0,:,:,1]+basicCountsSites[0,:,:,2]+basicCountsSites[1,:,:,1]+basicCountsSites[1,:,:,2]+basicCountsSites[2,:,:,1]+basicCountsSites[2,:,:,2],basicCountsSites[0,:,:,2]+basicCountsSites[1,:,:,2]+basicCountsSites[2,:,:,2]],['possible mut sites','sites with vars','sites with non-singleton vars','sites with frequent vars'],pathSimu+"barplot_counts_sites_wholeGenome.pdf",colors=colors2)	
plt.close()
#barplot([basicCountsOpportunities[1,:,:],basicCounts[1,:,:,0],basicCounts[1,:,:,1],basicCounts[1,:,:,2],basicCountsSites[1,:,:,0]+basicCountsSites[1,:,:,1]+basicCountsSites[1,:,:,2],basicCountsSites[1,:,:,1]+basicCountsSites[1,:,:,2],basicCountsSites[1,:,:,2]],['possible mut sites','1-descendant muts','few descendants muts','many descendants muts','variants obbserved','non-singleton variants','frequent variants'],pathSimu+"barplot_counts_all_non-coding.pdf")

#rates
colors=["red","orange","yellow","darkblue","blue","lightblue"]
colors2=["darkblue","blue","lightblue"]
colors3=["red","darkgreen","orange","green","yellow","lightgreen"]
for k in range(4):
	barplot([basicCounts[k,:,:,0]/basicCountsOpportunities[k,:,:],basicCounts[k,:,:,1]/basicCountsOpportunities[k,:,:],basicCounts[k,:,:,2]/basicCountsOpportunities[k,:,:],(basicCountsSites[k,:,:,0]+basicCountsSites[k,:,:,1]+basicCountsSites[k,:,:,2])/basicCountsOpportunities[k,:,:],(basicCountsSites[k,:,:,1]+basicCountsSites[k,:,:,2])/basicCountsOpportunities[k,:,:],basicCountsSites[k,:,:,2]/basicCountsOpportunities[k,:,:]],['1-descendant muts','2-4 descendants muts','>4 descendants muts','variants observed','non-singleton variants','frequent variants'],pathSimu+"barplot_rates_all_"+names[k]+".pdf",colors=colors)
	barplot([basicCounts[k,:,:,0]/basicCountsOpportunities[k,:,:],basicCounts[k,:,:,1]/basicCountsOpportunities[k,:,:],basicCounts[k,:,:,2]/basicCountsOpportunities[k,:,:]],['1-descendant muts','2-4 descendants muts','>4 descendants muts'],pathSimu+"barplot_rates_mutations_"+names[k]+".pdf",colors=colors)
	barplot([(basicCountsSites[k,:,:,0]+basicCountsSites[k,:,:,1]+basicCountsSites[k,:,:,2])/basicCountsOpportunities[k,:,:],(basicCountsSites[k,:,:,1]+basicCountsSites[k,:,:,2])/basicCountsOpportunities[k,:,:],basicCountsSites[k,:,:,2]/basicCountsOpportunities[k,:,:]],['variants observed','non-singleton variants','frequent variants'],pathSimu+"barplot_rates_sites_"+names[k]+".pdf",colors=colors2)	
	plt.close()
barplot([(basicCounts[0,:,:,0]+basicCounts[1,:,:,0]+basicCounts[2,:,:,0])/(basicCountsOpportunities[0,:,:]+basicCountsOpportunities[1,:,:]+basicCountsOpportunities[2,:,:]),(basicCounts[0,:,:,1]+basicCounts[1,:,:,1]+basicCounts[2,:,:,1])/(basicCountsOpportunities[0,:,:]+basicCountsOpportunities[1,:,:]+basicCountsOpportunities[2,:,:]),(basicCounts[0,:,:,2]+basicCounts[1,:,:,2]+basicCounts[2,:,:,2])/(basicCountsOpportunities[0,:,:]+basicCountsOpportunities[1,:,:]+basicCountsOpportunities[2,:,:]),(basicCountsSites[0,:,:,0]+basicCountsSites[0,:,:,1]+basicCountsSites[0,:,:,2]+basicCountsSites[1,:,:,0]+basicCountsSites[1,:,:,1]+basicCountsSites[1,:,:,2]+basicCountsSites[2,:,:,0]+basicCountsSites[2,:,:,1]+basicCountsSites[2,:,:,2])/(basicCountsOpportunities[0,:,:]+basicCountsOpportunities[1,:,:]+basicCountsOpportunities[2,:,:]),(basicCountsSites[0,:,:,1]+basicCountsSites[0,:,:,2]+basicCountsSites[1,:,:,1]+basicCountsSites[1,:,:,2]+basicCountsSites[2,:,:,1]+basicCountsSites[2,:,:,2])/(basicCountsOpportunities[0,:,:]+basicCountsOpportunities[1,:,:]+basicCountsOpportunities[2,:,:]),(basicCountsSites[0,:,:,2]+basicCountsSites[1,:,:,2]+basicCountsSites[2,:,:,2])/(basicCountsOpportunities[0,:,:]+basicCountsOpportunities[1,:,:]+basicCountsOpportunities[2,:,:])],['1-descendant muts','2-4 descendants muts','>4 descendants muts','variants observed','non-singleton variants','frequent variants'],pathSimu+"barplot_rates_all_wholeGenome.pdf",colors=colors)
barplot([(basicCounts[0,:,:,0]+basicCounts[1,:,:,0]+basicCounts[2,:,:,0])/(basicCountsOpportunities[0,:,:]+basicCountsOpportunities[1,:,:]+basicCountsOpportunities[2,:,:]),(basicCounts[0,:,:,1]+basicCounts[1,:,:,1]+basicCounts[2,:,:,1])/(basicCountsOpportunities[0,:,:]+basicCountsOpportunities[1,:,:]+basicCountsOpportunities[2,:,:]),(basicCounts[0,:,:,2]+basicCounts[1,:,:,2]+basicCounts[2,:,:,2])/(basicCountsOpportunities[0,:,:]+basicCountsOpportunities[1,:,:]+basicCountsOpportunities[2,:,:])],['1-descendant muts','2-4 descendants muts','>4 descendants muts'],pathSimu+"barplot_rates_mutations_wholeGenome.pdf",colors=colors)
barplot([(basicCountsSites[0,:,:,0]+basicCountsSites[0,:,:,1]+basicCountsSites[0,:,:,2]+basicCountsSites[1,:,:,0]+basicCountsSites[1,:,:,1]+basicCountsSites[1,:,:,2]+basicCountsSites[2,:,:,0]+basicCountsSites[2,:,:,1]+basicCountsSites[2,:,:,2])/(basicCountsOpportunities[0,:,:]+basicCountsOpportunities[1,:,:]+basicCountsOpportunities[2,:,:]),(basicCountsSites[0,:,:,1]+basicCountsSites[0,:,:,2]+basicCountsSites[1,:,:,1]+basicCountsSites[1,:,:,2]+basicCountsSites[2,:,:,1]+basicCountsSites[2,:,:,2])/(basicCountsOpportunities[0,:,:]+basicCountsOpportunities[1,:,:]+basicCountsOpportunities[2,:,:]),(basicCountsSites[0,:,:,2]+basicCountsSites[1,:,:,2]+basicCountsSites[2,:,:,2])/(basicCountsOpportunities[0,:,:]+basicCountsOpportunities[1,:,:]+basicCountsOpportunities[2,:,:])],['sites with vars','sites with non-singleton vars','sites with frequent vars'],pathSimu+"barplot_rates_sites_wholeGenome.pdf",colors=colors2)	
#nonsense
barplot([basicCounts[4,:,:,0]/basicCountsOpportunities[4,:,:],basicCounts[4,:,:,1]/basicCountsOpportunities[4,:,:],basicCounts[4,:,:,2]/basicCountsOpportunities[4,:,:],(basicCountsSites[4,:,:,0]+basicCountsSites[4,:,:,1]+basicCountsSites[4,:,:,2])/basicCountsOpportunities[4,:,:],(basicCountsSites[4,:,:,1]+basicCountsSites[4,:,:,2])/basicCountsOpportunities[4,:,:],basicCountsSites[4,:,:,2]/basicCountsOpportunities[4,:,:]],['1-descendant muts','2-4 descendants muts','>4 descendants muts','variants observed','non-singleton variants','frequent variants'],pathSimu+"barplot_rates_all_"+names[4]+".pdf",colors=colors)
plt.close()
#normalized synonymous rates
unnormalizedRates=np.zeros((3,4,4))
normalizedRates=np.zeros((3,4,4))
for i in range(3):
	tot=0.0
	for aa in range(4):
		for bb in range(4):
			if bb!=aa:
				unnormalizedRates[i,aa,bb]=basicCounts[1,aa,bb,i]/basicCountsOpportunities[1,aa,bb]
				tot+=unnormalizedRates[i,aa,bb]
	for aa in range(4):
		for bb in range(4):
			if bb!=aa:
				normalizedRates[i,aa,bb]=unnormalizedRates[i,aa,bb]/tot
barplot([normalizedRates[0,:,:],normalizedRates[1,:,:],normalizedRates[2,:,:]],['1-descendant muts','2-4 descendants muts','>4 descendants muts'],pathSimu+"barplot_normRates_mutations_"+names[k]+".pdf",colors=colors)

#new plots with tree-wise estimates
for k in range(3):
	unnormalizedRates=np.zeros((3,4,4))
	normalizedRates=np.zeros((3,4,4))
	unnormalizedRates2=np.zeros((3,4,4))
	normalizedRates2=np.zeros((3,4,4))
	for i in range(3):
		tot=0.0
		tot2=0.0
		for aa in range(4):
			for bb in range(4):
				if bb!=aa:
					unnormalizedRates[i,aa,bb]=basicCounts[k,aa,bb,i]/basicCountsOpportunities[k,aa,bb]
					tot+=unnormalizedRates[i,aa,bb]
					unnormalizedRates2[i,aa,bb]=mutationRatesTree[k,aa,bb,i]
					tot2+=unnormalizedRates2[i,aa,bb]
		for aa in range(4):
			for bb in range(4):
				if bb!=aa:
					normalizedRates[i,aa,bb]=unnormalizedRates[i,aa,bb]/tot
					normalizedRates2[i,aa,bb]=unnormalizedRates2[i,aa,bb]/tot2
	barplot([normalizedRates[0,:,:],normalizedRates2[0,:,:],normalizedRates[1,:,:],normalizedRates2[1,:,:],normalizedRates[2,:,:],normalizedRates2[2,:,:]],['1-descendant muts','1-descendant muts tree','2-4 descendants muts','2-4 descendants muts tree','>4 descendants muts','>4 descendants muts tree'],pathSimu+"barplot_normRates_mutations_tree_"+names[k]+".pdf",colors=colors3)
plt.close()
				

#PRINT to screen mutation rates at 4-fold degenerate sites, then decrease the G->T mutation rate 9-fold;
#estimate equilibrium frequencies in both cases, and compare to observed frequencies at 4-fold degenerate sites.
Rates4fold=np.zeros((4,4))
Rates4fold2=np.zeros((4,4))
Rates4fold3=np.zeros((4,4))
for aa in range(4):
	for bb in range(4):
		if aa!=bb:
			Rates4fold[aa,bb]=basicCounts[3,aa,bb,1]/basicCountsOpportunities[3,aa,bb]
			Rates4fold2[aa,bb]=basicCounts[3,aa,bb,2]/basicCountsOpportunities[3,aa,bb]
			Rates4fold3[aa,bb]=Rates4fold2[aa,bb]+Rates4fold[aa,bb]
Rates4foldDecreasedGU=np.copy(Rates4fold)
Rates4fold2DecreasedGU=np.copy(Rates4fold2)
Rates4fold3DecreasedGU=np.copy(Rates4fold3)
Rates4foldDecreasedGU[2,3]=Rates4foldDecreasedGU[2,3]/9
Rates4fold2DecreasedGU[2,3]=Rates4fold2DecreasedGU[2,3]/9
Rates4fold3DecreasedGU[2,3]=Rates4fold3DecreasedGU[2,3]/9
RiceEtAl=np.array([[0.0,0.016,0.129,0.049],[0.058,0.0,0.016,0.541],[0.239,0.051,0.0,0.522],[0.025,0.115,0.013,0.0]])
RiceEtAlDecreasedGU=np.copy(RiceEtAl)
RiceEtAlDecreasedGU[2,3]=RiceEtAlDecreasedGU[2,3]/9
#print(RiceEtAl)
for aa in range(4):
	Rates4fold[aa,aa]=-np.sum(Rates4fold[aa,:])
	Rates4fold2[aa,aa]=-np.sum(Rates4fold2[aa,:])
	Rates4fold3[aa,aa]=-np.sum(Rates4fold3[aa,:])
	Rates4foldDecreasedGU[aa,aa]=-np.sum(Rates4foldDecreasedGU[aa,:])
	Rates4fold2DecreasedGU[aa,aa]=-np.sum(Rates4fold2DecreasedGU[aa,:])
	Rates4fold3DecreasedGU[aa,aa]=-np.sum(Rates4fold3DecreasedGU[aa,:])
	RiceEtAl[aa,aa]=-np.sum(RiceEtAl[aa,:])
	RiceEtAlDecreasedGU[aa,aa]=-np.sum(RiceEtAlDecreasedGU[aa,:])
print("Mutation rates 4-fold from Rice et al 2020:")
print(RiceEtAl)
print("Rates at 4-fold degenerate sites: from low-frequency variants")
print(Rates4fold)
print("from high-frequency variants")
print(Rates4fold2)
print("from both combined:")
print(Rates4fold3)
	
print("equilibrium frequencies low freq:")
#P = np.array([[0.5,0.5],[0.6,0.4]])
mc = markovChain(Rates4fold)
mc.computePi('linear') #We can also use 'power', 'krylov' or 'eigen'
print(mc.pi)
print("equilibrium frequencies high freq:")
mc = markovChain(Rates4fold2)
mc.computePi('linear')
print(mc.pi)
print("equilibrium frequencies low and high freq muts combined:")
mc = markovChain(Rates4fold3)
mc.computePi('linear')
print(mc.pi)
print("equilibrium frequencies low freq reduced GU:")
mc = markovChain(Rates4foldDecreasedGU)
mc.computePi('linear')
print(mc.pi)
print("equilibrium frequencies high freq reduced GU:")
mc = markovChain(Rates4fold2DecreasedGU)
mc.computePi('linear')
print(mc.pi)
print("equilibrium frequencies from low and high freq muts combined, reduced GU:")
mc = markovChain(Rates4fold3DecreasedGU)
mc.computePi('linear')
print(mc.pi)
print("equilibrium frequencies in Rice et al 4-fold sites:")
mc = markovChain(RiceEtAl)
mc.computePi('linear')
print(mc.pi)
print("equilibrium frequencies in Rice et al 4-fold sites with reduced GU rate:")
mc = markovChain(RiceEtAlDecreasedGU)
mc.computePi('linear')
print(mc.pi)

#observed frequencies
ObservedFreq4fold=np.zeros(4)
ObservedFreq4fold[0]=basicCountsOpportunities[3,0,1]
ObservedFreq4fold[1]=basicCountsOpportunities[3,1,0]
ObservedFreq4fold[2]=basicCountsOpportunities[3,2,1]
ObservedFreq4fold[3]=basicCountsOpportunities[3,3,1]
ObservedFreq4fold=ObservedFreq4fold/np.sum(ObservedFreq4fold)
print("observed frequencies:")
print(ObservedFreq4fold)


#Plot rates for triplets separately ACA->ATA, ACC->ATC,... AGA->ATA, AGC->ATC, ... and maybe stack the 3 types of mutation on top of each other for mutation counts? 
#It would be also nice to add a measure of uncertainty. By approximating the Poisson distributions as normal distributions, and representing the estimated rate as \hat{\lambda}, then an estimate of 95% confidence interval is \hat{\lambda}+- 1.96*math.sqrt(\hat{\lambda}/n)    where n is the number of possible mutatable sites. 
contextCountsC=np.zeros((4,4,4,3),dtype=int)
contextCountsSitesC=np.zeros((4,4,4,3),dtype=int)
contextCountsOpportunitiesC=np.zeros((4,4,4))
contextRatesC=np.zeros((4,4,4,3))
contextRatesSitesC=np.zeros((4,4,4,3))
contextCountsG=np.zeros((4,4,4,3),dtype=int)
contextCountsSitesG=np.zeros((4,4,4,3),dtype=int)
contextCountsOpportunitiesG=np.zeros((4,4,4))
contextRatesG=np.zeros((4,4,4,3))
contextRatesSitesG=np.zeros((4,4,4,3))
for k in range(4):
	for i1 in range(4):
		#for i2 in range(4):
			for i3 in range(4):
				#for j in range(4):
					contextCountsOpportunitiesC[k,i1,i3]+=countsPossibleMutsContext[k,i1,1,i3,3]
					contextCountsOpportunitiesG[k,i1,i3]+=countsPossibleMutsContext[k,i1,2,i3,3]
					for s in range(100):
						if s==0:
							contextCountsC[k,i1,i3,0]+=countsMutsContext[k,i1,1,i3,3,s]
							contextCountsSitesC[k,i1,i3,0]+=countsMutsContextSites[k,i1,1,i3,3,s]
							contextCountsG[k,i1,i3,0]+=countsMutsContext[k,i1,2,i3,3,s]
							contextCountsSitesG[k,i1,i3,0]+=countsMutsContextSites[k,i1,2,i3,3,s]
						elif s>0 and s<4:
							contextCountsC[k,i1,i3,1]+=countsMutsContext[k,i1,1,i3,3,s]
							contextCountsSitesC[k,i1,i3,1]+=countsMutsContextSites[k,i1,1,i3,3,s]
							contextCountsG[k,i1,i3,1]+=countsMutsContext[k,i1,2,i3,3,s]
							contextCountsSitesG[k,i1,i3,1]+=countsMutsContextSites[k,i1,2,i3,3,s]
						else:
							contextCountsC[k,i1,i3,2]+=countsMutsContext[k,i1,1,i3,3,s]
							contextCountsSitesC[k,i1,i3,2]+=countsMutsContextSites[k,i1,1,i3,3,s]
							contextCountsG[k,i1,i3,2]+=countsMutsContext[k,i1,2,i3,3,s]
							contextCountsSitesG[k,i1,i3,2]+=countsMutsContextSites[k,i1,2,i3,3,s]
							
#context counts
names=["nonCoding","synonymous","nonSynonymous","4fold"]
colors=["green","red","orange","yellow","darkblue","blue","lightblue"]
for k in range(4):
	barplotContext([contextCountsOpportunitiesC[k,:,:],contextCountsC[k,:,:,0],contextCountsC[k,:,:,1],contextCountsC[k,:,:,2],contextCountsSitesC[k,:,:,0]+contextCountsSitesC[k,:,:,1]+contextCountsSitesC[k,:,:,2],contextCountsSitesC[k,:,:,1]+contextCountsSitesC[k,:,:,2],contextCountsSitesC[k,:,:,2]],['possible mut sites','1-descendant muts','2-4 descendants muts','>4 descendants muts','variants observed','non-singleton variants','frequent variants'],pathSimu+"barplot_contextC_counts_all_"+names[k]+".pdf",colors=colors)	
	plt.close()
barplotContext([contextCountsOpportunitiesC[0,:,:]+contextCountsOpportunitiesC[1,:,:]+contextCountsOpportunitiesC[2,:,:],contextCountsC[0,:,:,0]+contextCountsC[1,:,:,0]+contextCountsC[2,:,:,0],contextCountsC[0,:,:,1]+contextCountsC[1,:,:,1]+contextCountsC[2,:,:,1],contextCountsC[0,:,:,2]+contextCountsC[1,:,:,2]+contextCountsC[2,:,:,2],contextCountsSitesC[0,:,:,0]+contextCountsSitesC[0,:,:,1]+contextCountsSitesC[0,:,:,2]+contextCountsSitesC[1,:,:,0]+contextCountsSitesC[1,:,:,1]+contextCountsSitesC[1,:,:,2]+contextCountsSitesC[2,:,:,0]+contextCountsSitesC[2,:,:,1]+contextCountsSitesC[2,:,:,2],contextCountsSitesC[0,:,:,1]+contextCountsSitesC[0,:,:,2]+contextCountsSitesC[1,:,:,1]+contextCountsSitesC[1,:,:,2]+contextCountsSitesC[2,:,:,1]+contextCountsSitesC[2,:,:,2],contextCountsSitesC[0,:,:,2]+contextCountsSitesC[1,:,:,2]+contextCountsSitesC[2,:,:,2]],['possible mut sites','1-descendant muts','2-4 descendants muts','>4 descendants muts','variants observed','non-singleton variants','frequent variants'],pathSimu+"barplot_contextC_counts_all_wholeGenome.pdf",colors=colors)	
plt.close()

for k in range(4):
	barplotContext([contextCountsOpportunitiesG[k,:,:],contextCountsG[k,:,:,0],contextCountsG[k,:,:,1],contextCountsG[k,:,:,2],contextCountsSitesG[k,:,:,0]+contextCountsSitesG[k,:,:,1]+contextCountsSitesG[k,:,:,2],contextCountsSitesG[k,:,:,1]+contextCountsSitesG[k,:,:,2],contextCountsSitesG[k,:,:,2]],['possible mut sites','1-descendant muts','2-4 descendants muts','>4 descendants muts','variants observed','non-singleton variants','frequent variants'],pathSimu+"barplot_contextG_counts_all_"+names[k]+".pdf",colors=colors)	
	plt.close()
barplotContext([contextCountsOpportunitiesG[0,:,:]+contextCountsOpportunitiesG[1,:,:]+contextCountsOpportunitiesG[2,:,:],contextCountsG[0,:,:,0]+contextCountsG[1,:,:,0]+contextCountsG[2,:,:,0],contextCountsG[0,:,:,1]+contextCountsG[1,:,:,1]+contextCountsG[2,:,:,1],contextCountsG[0,:,:,2]+contextCountsG[1,:,:,2]+contextCountsG[2,:,:,2],contextCountsSitesG[0,:,:,0]+contextCountsSitesG[0,:,:,1]+contextCountsSitesG[0,:,:,2]+contextCountsSitesG[1,:,:,0]+contextCountsSitesG[1,:,:,1]+contextCountsSitesG[1,:,:,2]+contextCountsSitesG[2,:,:,0]+contextCountsSitesG[2,:,:,1]+contextCountsSitesG[2,:,:,2],contextCountsSitesG[0,:,:,1]+contextCountsSitesG[0,:,:,2]+contextCountsSitesG[1,:,:,1]+contextCountsSitesG[1,:,:,2]+contextCountsSitesG[2,:,:,1]+contextCountsSitesG[2,:,:,2],contextCountsSitesG[0,:,:,2]+contextCountsSitesG[1,:,:,2]+contextCountsSitesG[2,:,:,2]],['possible mut sites','1-descendant muts','2-4 descendants muts','>4 descendants muts','variants observed','non-singleton variants','frequent variants'],pathSimu+"barplot_contextG_counts_all_wholeGenome.pdf",colors=colors)	
plt.close()

#context rates
colors=["red","orange","yellow","darkblue","blue","lightblue"]
for k in range(4):
	barplotContext([contextCountsC[k,:,:,0]/contextCountsOpportunitiesC[k,:,:],contextCountsC[k,:,:,1]/contextCountsOpportunitiesC[k,:,:],contextCountsC[k,:,:,2]/contextCountsOpportunitiesC[k,:,:],(contextCountsSitesC[k,:,:,0]+contextCountsSitesC[k,:,:,1]+contextCountsSitesC[k,:,:,2])/contextCountsOpportunitiesC[k,:,:],(contextCountsSitesC[k,:,:,1]+contextCountsSitesC[k,:,:,2])/contextCountsOpportunitiesC[k,:,:],contextCountsSitesC[k,:,:,2]/contextCountsOpportunitiesC[k,:,:]],['1-descendant muts','2-4 descendants muts','>4 descendants muts','variants observed','non-singleton variants','frequent variants'],pathSimu+"barplot_contextC_rates_all_"+names[k]+".pdf",colors=colors)	
	plt.close()
barplotContext([(contextCountsC[0,:,:,0]+contextCountsC[1,:,:,0]+contextCountsC[2,:,:,0])/(contextCountsOpportunitiesC[0,:,:]+contextCountsOpportunitiesC[1,:,:]+contextCountsOpportunitiesC[2,:,:]),(contextCountsC[0,:,:,1]+contextCountsC[1,:,:,1]+contextCountsC[2,:,:,1])/(contextCountsOpportunitiesC[0,:,:]+contextCountsOpportunitiesC[1,:,:]+contextCountsOpportunitiesC[2,:,:]),(contextCountsC[0,:,:,2]+contextCountsC[1,:,:,2]+contextCountsC[2,:,:,2])/(contextCountsOpportunitiesC[0,:,:]+contextCountsOpportunitiesC[1,:,:]+contextCountsOpportunitiesC[2,:,:]),(contextCountsSitesC[0,:,:,0]+contextCountsSitesC[0,:,:,1]+contextCountsSitesC[0,:,:,2]+contextCountsSitesC[1,:,:,0]+contextCountsSitesC[1,:,:,1]+contextCountsSitesC[1,:,:,2]+contextCountsSitesC[2,:,:,0]+contextCountsSitesC[2,:,:,1]+contextCountsSitesC[2,:,:,2])/(contextCountsOpportunitiesC[0,:,:]+contextCountsOpportunitiesC[1,:,:]+contextCountsOpportunitiesC[2,:,:]),(contextCountsSitesC[0,:,:,1]+contextCountsSitesC[0,:,:,2]+contextCountsSitesC[1,:,:,1]+contextCountsSitesC[1,:,:,2]+contextCountsSitesC[2,:,:,1]+contextCountsSitesC[2,:,:,2])/(contextCountsOpportunitiesC[0,:,:]+contextCountsOpportunitiesC[1,:,:]+contextCountsOpportunitiesC[2,:,:]),(contextCountsSitesC[0,:,:,2]+contextCountsSitesC[1,:,:,2]+contextCountsSitesC[2,:,:,2])/(contextCountsOpportunitiesC[0,:,:]+contextCountsOpportunitiesC[1,:,:]+contextCountsOpportunitiesC[2,:,:])],['1-descendant muts','2-4 descendants muts','>4 descendants muts','variants observed','non-singleton variants','frequent variants'],pathSimu+"barplot_contextC_rates_all_wholeGenome.pdf",colors=colors)	
barplotContext([contextCountsC[1,:,:,0]/contextCountsOpportunitiesC[1,:,:],contextCountsC[1,:,:,1]/contextCountsOpportunitiesC[1,:,:],contextCountsC[1,:,:,2]/contextCountsOpportunitiesC[1,:,:]],['1-descendant muts','2-4 descendants muts','>4 descendants muts'],pathSimu+"barplot_contextC_rates_mutations_synonymous.pdf",colors=colors)	

plt.close()

for k in range(4):
	barplotContext([contextCountsG[k,:,:,0]/contextCountsOpportunitiesG[k,:,:],contextCountsG[k,:,:,1]/contextCountsOpportunitiesG[k,:,:],contextCountsG[k,:,:,2]/contextCountsOpportunitiesG[k,:,:],(contextCountsSitesG[k,:,:,0]+contextCountsSitesG[k,:,:,1]+contextCountsSitesG[k,:,:,2])/contextCountsOpportunitiesG[k,:,:],(contextCountsSitesG[k,:,:,1]+contextCountsSitesG[k,:,:,2])/contextCountsOpportunitiesG[k,:,:],contextCountsSitesG[k,:,:,2]/contextCountsOpportunitiesG[k,:,:]],['1-descendant muts','2-4 descendants muts','>4 descendants muts','variants observed','non-singleton variants','frequent variants'],pathSimu+"barplot_contextG_rates_all_"+names[k]+".pdf",colors=colors)	
	plt.close()
barplotContext([(contextCountsG[0,:,:,0]+contextCountsG[1,:,:,0]+contextCountsG[2,:,:,0])/(contextCountsOpportunitiesG[0,:,:]+contextCountsOpportunitiesG[1,:,:]+contextCountsOpportunitiesG[2,:,:]),(contextCountsG[0,:,:,1]+contextCountsG[1,:,:,1]+contextCountsG[2,:,:,1])/(contextCountsOpportunitiesG[0,:,:]+contextCountsOpportunitiesG[1,:,:]+contextCountsOpportunitiesG[2,:,:]),(contextCountsG[0,:,:,2]+contextCountsG[1,:,:,2]+contextCountsG[2,:,:,2])/(contextCountsOpportunitiesG[0,:,:]+contextCountsOpportunitiesG[1,:,:]+contextCountsOpportunitiesG[2,:,:]),(contextCountsSitesG[0,:,:,0]+contextCountsSitesG[0,:,:,1]+contextCountsSitesG[0,:,:,2]+contextCountsSitesG[1,:,:,0]+contextCountsSitesG[1,:,:,1]+contextCountsSitesG[1,:,:,2]+contextCountsSitesG[2,:,:,0]+contextCountsSitesG[2,:,:,1]+contextCountsSitesG[2,:,:,2])/(contextCountsOpportunitiesG[0,:,:]+contextCountsOpportunitiesG[1,:,:]+contextCountsOpportunitiesG[2,:,:]),(contextCountsSitesG[0,:,:,1]+contextCountsSitesG[0,:,:,2]+contextCountsSitesG[1,:,:,1]+contextCountsSitesG[1,:,:,2]+contextCountsSitesG[2,:,:,1]+contextCountsSitesG[2,:,:,2])/(contextCountsOpportunitiesG[0,:,:]+contextCountsOpportunitiesG[1,:,:]+contextCountsOpportunitiesG[2,:,:]),(contextCountsSitesG[0,:,:,2]+contextCountsSitesG[1,:,:,2]+contextCountsSitesG[2,:,:,2])/(contextCountsOpportunitiesG[0,:,:]+contextCountsOpportunitiesG[1,:,:]+contextCountsOpportunitiesG[2,:,:])],['1-descendant muts','2-4 descendants muts','>4 descendants muts','variants observed','non-singleton variants','frequent variants'],pathSimu+"barplot_contextG_rates_all_wholeGenome.pdf",colors=colors)	
barplotContext([contextCountsG[1,:,:,0]/contextCountsOpportunitiesG[1,:,:],contextCountsG[1,:,:,1]/contextCountsOpportunitiesG[1,:,:],contextCountsG[1,:,:,2]/contextCountsOpportunitiesG[1,:,:]],['1-descendant muts','2-4 descendants muts','>4 descendants muts'],pathSimu+"barplot_contextG_rates_mutations_synonymous.pdf",colors=colors)	
plt.close()

print("\n Now plotting selection tests and its significance")
#Now plot CpG and all-v-U mutation and selection effects.
#decreasing, neutral, increasing ; singleton, low freq, high freq
countsCpGsyn=np.zeros((3,3))
countsCpGsynOpportunities=np.zeros(3,dtype=int)
countsGCsyn=np.zeros((3,3))
countsGCsynOpportunities=np.zeros(3,dtype=int)
countsUsyn=np.zeros((3,3))
countsUsynOpportunities=np.zeros(3,dtype=int)
countsUsynGene=np.zeros((len(geneEnds),3,3))
countsUsynOpportunitiesGene=np.zeros((len(geneEnds),3),dtype=int)
#for k in range(4):
for i2 in range(4):
				for j in range(4):
					if i2!=j:
						if i2==3:
							U=0
						elif j==3:
							U=2
						else:
							U=1
						for g in range(len(geneEnds)):
							countsUsynOpportunitiesGene[g,U]+=countsPossibleMutsContextU[g,i2,j]
							for s in range(100):
								if s==0:
									countsUsynGene[g,U,0]+=countsMutsContextU[g,i2,j,s]
								elif s>0 and s<4:
									countsUsynGene[g,U,1]+=countsMutsContextU[g,i2,j,s]
								else:
									countsUsynGene[g,U,2]+=countsMutsContextU[g,i2,j,s]
for i1 in range(4):
		for i2 in range(4):
			for i3 in range(4):
				for j in range(4):
					if i2!=j:
						if i2==3:
							U=0
						elif j==3:
							U=2
						else:
							U=1
						if i2==1 and i3==2 and (j!=2 or i1!=1):
							CpG=0
						elif j==1 and i3==2 and (i2!=2 or i1!=1):
							CpG=2
						elif i1==1 and i2==2 and (j!=1 or i3!=2):
							CpG=0
						elif i1==1 and j==2 and (i2!=1 or i3!=2):
							CpG=2
						else:
							CpG=1
						if (i2==0 or i2==3) and (j==1 or j==2):
							GC=2
						elif (i2==1 or i2==2) and (j==0 or j==3):
							GC=0
						else:
							GC=1
						
						countsCpGsynOpportunities[CpG]+=countsPossibleMutsContext[1,i1,i2,i3,j]
						countsUsynOpportunities[U]+=countsPossibleMutsContext[1,i1,i2,i3,j]
						countsGCsynOpportunities[GC]+=countsPossibleMutsContext[1,i1,i2,i3,j]
						for s in range(100):
							if s==0:
								countsCpGsyn[CpG,0]+=countsMutsContext[1,i1,i2,i3,j,s]
								countsUsyn[U,0]+=countsMutsContext[1,i1,i2,i3,j,s]
								countsGCsyn[GC,0]+=countsMutsContext[1,i1,i2,i3,j,s]
							elif s>0 and s<4:
								countsCpGsyn[CpG,1]+=countsMutsContext[1,i1,i2,i3,j,s]
								countsUsyn[U,1]+=countsMutsContext[1,i1,i2,i3,j,s]
								countsGCsyn[GC,1]+=countsMutsContext[1,i1,i2,i3,j,s]
							else:
								countsCpGsyn[CpG,2]+=countsMutsContext[1,i1,i2,i3,j,s]
								countsUsyn[U,2]+=countsMutsContext[1,i1,i2,i3,j,s]
								countsGCsyn[GC,2]+=countsMutsContext[1,i1,i2,i3,j,s]
#print(countsCpGsynOpportunities)
#print(countsUsynOpportunities)
#print(countsCpGsyn[:,1])
#print(countsCpGsyn[:,0])
#print(countsCpGsyn[:,1]/countsCpGsyn[:,0])				
#context counts
#names=["nonCoding","synonymous","nonSynonymous","4fold"]
#labels=['<CpG','=CpG','>CpG']
print("\n Tests on CpG selection")
labels2=['<CpG','=CpG','>CpG']
colors=["green","red","orange","yellow"]
barplotContext([countsCpGsynOpportunities,countsCpGsyn[:,0],countsCpGsyn[:,1],countsCpGsyn[:,2]],['possible mut sites','1-descendant muts','2-4 descendants muts','>4 descendants muts'],pathSimu+"barplot_CpG_counts_all_synonymous.pdf",xName='Mutation type',colors=colors,labels=labels2)	
plt.close()
#REVERT ORDER - use colors for CpG effect, and have 3 sets of bars one for each comparison
colors=["purple","brown","black"]
labels=['low vs single','high vs low','high vs single']
comparison=np.zeros((3,3))
for i in range(3):
	comparison[0,i]=countsCpGsyn[i,1]/countsCpGsyn[i,0]
	comparison[1,i]=countsCpGsyn[i,2]/countsCpGsyn[i,1]
	comparison[2,i]=countsCpGsyn[i,2]/countsCpGsyn[i,0]
for i in range(3):
	if i==0:
		i1=1
		i2=0
	elif i==1:
		i1=2
		i2=1
	else:
		i1=2
		i2=0
	print(labels2[i1]+" vs "+labels2[i2])
	for j in range(3):
		print(labels[j])
		print(j)
		if j==0:
			j1=1
			j2=0
		elif j==1:
			j1=2
			j2=1
		else:
			j1=2
			j2=0
		print(chi2_contingency([[countsCpGsyn[i1,j1],  countsCpGsyn[i1,j2]],[countsCpGsyn[i2,j1], countsCpGsyn[i2,j2]]]))
		print(str(countsCpGsyn[i1,j1])+" "+str(countsCpGsyn[i1,j2])+" "+str(countsCpGsyn[i2,j1])+" "+str(countsCpGsyn[i2,j2]))
		#population1 = np.random.binomial(1, countsCpGsyn[i1,j1]/(countsCpGsyn[i1,j2]+countsCpGsyn[i1,j1]), int(countsCpGsyn[i1,j2]+countsCpGsyn[i1,j1]))
		#population2 = np.random.binomial(1, countsCpGsyn[i2,j1]/(countsCpGsyn[i2,j2]+countsCpGsyn[i2,j1]), int(countsCpGsyn[i2,j2]+countsCpGsyn[i2,j1]))
		#print(sm.stats.ttest_ind(population1, population2))
fig=barplotContext([comparison[:,0],comparison[:,1],comparison[:,2]],['<CpG','=CpG','>CpG'],pathSimu+"barplot_CpG_ratios_all_synonymous.pdf",noSave=True,xName='Comparison',colors=colors,labels=labels)	
#plt.text(1.0, 0.53, "p=0.014", ha='center', va='bottom', color='black')
#plt.plot([0.75, 0.75, 1.25, 1.25], [0.515, 0.525, 0.525, 0.515], lw=1.5, c='black')
plt.text(1.125, 0.57, "p=0.073", ha='center', va='bottom', color='black')
plt.plot([1.0, 1.0, 1.25, 1.25], [0.56, 0.565, 0.565, 0.56], lw=1.5, c='black')
fig.savefig(pathSimu+"barplot_CpG_ratios_all_synonymous.pdf")
plt.close()

print("\n Tests on U selection")
labels2=['<U','=U','>U']
colors=["green","red","orange","yellow"]
barplotContext([countsUsynOpportunities,countsUsyn[:,0],countsUsyn[:,1],countsUsyn[:,2]],['possible mut sites','1-descendant muts','2-4 descendants muts','>4 descendants muts'],pathSimu+"barplot_U_counts_all_synonymous.pdf",xName='Mutation type',colors=colors,labels=labels2)	
plt.close()
colors=["purple","brown","black"]
labels=['low vs single','high vs low','high vs single']
comparison=np.zeros((3,3))
for i in range(3):
	comparison[0,i]=countsUsyn[i,1]/countsUsyn[i,0]
	comparison[1,i]=countsUsyn[i,2]/countsUsyn[i,1]
	comparison[2,i]=countsUsyn[i,2]/countsUsyn[i,0]
for i in range(3):
	if i==0:
		i1=1
		i2=0
	elif i==1:
		i1=2
		i2=1
	else:
		i1=2
		i2=0
	print(labels2[i1]+" vs "+labels2[i2])
	for j in range(3):
		print(labels[j])
		if j==0:
			j1=1
			j2=0
		elif j==1:
			j1=2
			j2=1
		else:
			j1=2
			j2=0
		print(chi2_contingency([[countsUsyn[i1,j1],  countsUsyn[i1,j2]],[countsUsyn[i2,j1], countsUsyn[i2,j2]]]))
		print(str(countsUsyn[i1,j1])+" "+str(countsUsyn[i1,j2])+" "+str(countsUsyn[i2,j1])+" "+str(countsUsyn[i2,j2]))
		#population1 = np.random.binomial(1, countsUsyn[i1,j1]/(countsUsyn[i1,j2]+countsUsyn[i1,j1]), int(countsUsyn[i1,j2]+countsUsyn[i1,j1]))
		#population2 = np.random.binomial(1, countsUsyn[i2,j1]/(countsUsyn[i2,j2]+countsUsyn[i2,j1]), int(countsUsyn[i2,j2]+countsUsyn[i2,j1]))
		#print(sm.stats.ttest_ind(population1, population2))
fig=barplotContext([comparison[:,0],comparison[:,1],comparison[:,2]],['<U','=U','>U'],pathSimu+"barplot_U_ratios_all_synonymous.pdf",noSave=True,xName='Comparison',colors=colors,labels=labels)	
plt.text(2.0, 0.28, "p=0.0083", ha='center', va='bottom', color='black')
plt.plot([1.75, 1.75, 2.25, 2.25], [0.27, 0.275, 0.275, 0.27], lw=1.5, c='black')
plt.text(0.0, 0.41, "p=0.051", ha='center', va='bottom', color='black')
plt.plot([-0.25, -0.25, 0.25, 0.25], [0.40, 0.405, 0.405, 0.40], lw=1.5, c='black')
plt.text(2.125, 0.24, "p=1.6e-05", ha='center', va='bottom', color='black')
plt.plot([2.0, 2.0, 2.25, 2.25], [0.23, 0.235, 0.235, 0.23], lw=1.5, c='black')
plt.text(1.125, 0.59, "p=0.0037", ha='center', va='bottom', color='black')
plt.plot([1.0, 1.0, 1.25, 1.25], [0.58, 0.585, 0.585, 0.58], lw=1.5, c='black')
fig.savefig(pathSimu+"barplot_U_ratios_all_synonymous.pdf")
#plt.hlines(countsUsyn[1,1]/countsUsyn[1,0], 0, 100, colors=colors[0], linestyles='dashed')
plt.close()

print("\n Tests on U selection, gene-wise")
labels2=['<U','=U','>U']
labels=['low vs single','high vs low','high vs single']
comparison=np.zeros((len(geneEnds),3,3))
for g in range(len(geneEnds)):
	colors=["green","red","orange","yellow"]
	barplotContext([countsUsynOpportunitiesGene[g],countsUsynGene[g,:,0],countsUsynGene[g,:,1],countsUsynGene[g,:,2]],['possible mut sites','1-descendant muts','2-4 descendants muts','>4 descendants muts'],pathSimu+"barplot_U_counts_all_synonymous"+geneEndsNames[g]+".pdf",xName='Mutation type',colors=colors,labels=labels2)	
	plt.close()
	colors=["purple","brown","black"]
	print("\n"+geneEndsNames[g])
	for i in range(3):
		comparison[g,0,i]=countsUsynGene[g,i,1]/countsUsynGene[g,i,0]
		comparison[g,1,i]=countsUsynGene[g,i,2]/countsUsynGene[g,i,1]
		comparison[g,2,i]=countsUsynGene[g,i,2]/countsUsynGene[g,i,0]
	for i in range(3):
		if i==0:
			i1=1
			i2=0
		elif i==1:
			i1=2
			i2=1
		else:
			i1=2
			i2=0
		print(labels2[i1]+" vs "+labels2[i2])
		for j in range(3):
			print(labels[j])
			if j==0:
				j1=1
				j2=0
			elif j==1:
				j1=2
				j2=1
			else:
				j1=2
				j2=0
			print(chi2_contingency([[countsUsynGene[g,i1,j1],  countsUsynGene[g,i1,j2]],[countsUsynGene[g,i2,j1], countsUsynGene[g,i2,j2]]]))
			print(str(countsUsynGene[g,i1,j1])+" "+str(countsUsynGene[g,i1,j2])+" "+str(countsUsynGene[g,i2,j1])+" "+str(countsUsynGene[g,i2,j2]))
			#population1 = np.random.binomial(1, countsUsyn[i1,j1]/(countsUsyn[i1,j2]+countsUsyn[i1,j1]), int(countsUsyn[i1,j2]+countsUsyn[i1,j1]))
			#population2 = np.random.binomial(1, countsUsyn[i2,j1]/(countsUsyn[i2,j2]+countsUsyn[i2,j1]), int(countsUsyn[i2,j2]+countsUsyn[i2,j1]))
			#print(sm.stats.ttest_ind(population1, population2))
	fig=barplotContext([comparison[g,:,0],comparison[g,:,1],comparison[g,:,2]],['<U','=U','>U'],pathSimu+"barplot_U_ratios_all_synonymous_gene"+geneEndsNames[g]+".pdf",noSave=True,xName='Comparison',colors=colors,labels=labels)	
	if g==9:
		plt.text(0.125, 0.63, "p=0.049", ha='center', va='bottom', color='black')
		plt.plot([0.0, 0.0, 0.25, 0.25], [0.62, 0.625, 0.625, 0.62], lw=1.5, c='black')
	if g==2:
		plt.text(2.0, 0.28, "p=0.042", ha='center', va='bottom', color='black')
		plt.plot([1.75, 1.75, 2.25, 2.25], [0.26, 0.265, 0.265, 0.26], lw=1.5, c='black')
	if g==1:
		plt.text(2.125, 0.22, "p=0.046", ha='center', va='bottom', color='black')
		plt.plot([2.0, 2.0, 2.25, 2.25], [0.21, 0.215, 0.215, 0.21], lw=1.5, c='black')
	if g==0:
		plt.text(2.125, 0.25, "p=0.00027", ha='center', va='bottom', color='black')
		plt.plot([2.0, 2.0, 2.25, 2.25], [0.24, 0.245, 0.245, 0.24], lw=1.5, c='black')
		plt.text(1.125, 0.58, "p=0.016", ha='center', va='bottom', color='black')
		plt.plot([1.0, 1.0, 1.25, 1.25], [0.57, 0.575, 0.575, 0.57], lw=1.5, c='black')
		plt.text(1.875, 0.22, "p=0.020", ha='center', va='bottom', color='black')
		plt.plot([1.75, 1.75, 2.0, 2.0], [0.21, 0.215, 0.215, 0.21], lw=1.5, c='black')
	fig.savefig(pathSimu+"barplot_U_ratios_all_synonymous_gene"+geneEndsNames[g]+".pdf")
	#plt.hlines(countsUsyn[1,1]/countsUsyn[1,0], 0, 100, colors=colors[0], linestyles='dashed')
	plt.close()
	#barplotGene()

print("\n Tests on GC selection")
labels2=['<GC','=GC','>GC']
colors=["green","red","orange","yellow"]
barplotContext([countsGCsynOpportunities,countsGCsyn[:,0],countsGCsyn[:,1],countsGCsyn[:,2]],['possible mut sites','1-descendant muts','2-4 descendants muts','>4 descendants muts'],pathSimu+"barplot_GC_counts_all_synonymous.pdf",xName='Mutation type',colors=colors,labels=labels2)	
plt.close()
#REVERT ORDER - use colors for CpG effect, and have 3 sets of bars one for each comparison
colors=["purple","brown","black"]
labels=['low vs single','high vs low','high vs single']
comparison=np.zeros((3,3))
for i in range(3):
	comparison[0,i]=countsGCsyn[i,1]/countsGCsyn[i,0]
	comparison[1,i]=countsGCsyn[i,2]/countsGCsyn[i,1]
	comparison[2,i]=countsGCsyn[i,2]/countsGCsyn[i,0]
for i in range(3):
	if i==0:
		i1=1
		i2=0
	elif i==1:
		i1=2
		i2=1
	else:
		i1=2
		i2=0
	print(labels2[i1]+" vs "+labels2[i2])
	for j in range(3):
		print(labels[j])
		print(j)
		if j==0:
			j1=1
			j2=0
		elif j==1:
			j1=2
			j2=1
		else:
			j1=2
			j2=0
		print(chi2_contingency([[countsGCsyn[i1,j1],  countsGCsyn[i1,j2]],[countsGCsyn[i2,j1], countsGCsyn[i2,j2]]]))
		print(str(countsGCsyn[i1,j1])+" "+str(countsGCsyn[i1,j2])+" "+str(countsGCsyn[i2,j1])+" "+str(countsGCsyn[i2,j2]))
		#population1 = np.random.binomial(1, countsCpGsyn[i1,j1]/(countsCpGsyn[i1,j2]+countsCpGsyn[i1,j1]), int(countsCpGsyn[i1,j2]+countsCpGsyn[i1,j1]))
		#population2 = np.random.binomial(1, countsCpGsyn[i2,j1]/(countsCpGsyn[i2,j2]+countsCpGsyn[i2,j1]), int(countsCpGsyn[i2,j2]+countsCpGsyn[i2,j1]))
		#print(sm.stats.ttest_ind(population1, population2))
fig=barplotContext([comparison[:,0],comparison[:,1],comparison[:,2]],['<GC','=GC','>GC'],pathSimu+"barplot_GC_ratios_all_synonymous.pdf",noSave=True,xName='Comparison',colors=colors,labels=labels)	
plt.text(1.875, 0.24, "p=0.072", ha='center', va='bottom', color='black')
plt.plot([1.75, 1.75, 2.0, 2.0], [0.23, 0.235, 0.235, 0.23], lw=1.5, c='black')
plt.text(2.0, 0.29, "p=0.0032", ha='center', va='bottom', color='black')
plt.plot([1.75, 1.75, 2.25, 2.25], [0.28, 0.285, 0.285, 0.28], lw=1.5, c='black')
fig.savefig(pathSimu+"barplot_GC_ratios_all_synonymous.pdf")
plt.close()

#Synonymous-nonsynonymous

#Now plot CpG and all-v-U mutation and selection effects.
#decreasing, neutral, increasing ; singleton, low freq, high freq
countsSyn=np.zeros((2,3))
countsSynOpportunities=np.zeros(2,dtype=int)
#for k in range(4):
for i1 in range(4):
		for i2 in range(4):
			for i3 in range(4):
				for j in range(4):
					if i2!=j:
						countsSynOpportunities[0]+=countsPossibleMutsContext[1,i1,i2,i3,j]
						countsSynOpportunities[1]+=countsPossibleMutsContext[2,i1,i2,i3,j]
						for s in range(100):
							if s==0:
								countsSyn[0,0]+=countsMutsContext[1,i1,i2,i3,j,s]
								countsSyn[1,0]+=countsMutsContext[2,i1,i2,i3,j,s]
							elif s>0 and s<4:
								countsSyn[0,1]+=countsMutsContext[1,i1,i2,i3,j,s]
								countsSyn[1,1]+=countsMutsContext[2,i1,i2,i3,j,s]
							else:
								countsSyn[0,2]+=countsMutsContext[1,i1,i2,i3,j,s]
								countsSyn[1,2]+=countsMutsContext[2,i1,i2,i3,j,s]
#print(countsSynOpportunities)
#print(countsSyn)			

print("\n Tests on syn-nonsyn selection")
labels2=['Synonymous','Nonsynonymous']
colors=["green","red","orange","yellow"]
barplotContext([countsSynOpportunities,countsSyn[:,0],countsSyn[:,1],countsSyn[:,2]],['possible mut sites','1-descendant muts','2-4 descendants muts','>4 descendants muts'],pathSimu+"barplot_counts_all_syn-nonsyn.pdf",xName='Mutation type',colors=colors,labels=labels2)	
plt.close()
colors=["purple","brown","black"]
labels=['low vs single','high vs low','high vs single']
comparison=np.zeros((3,2))
for i in range(2):
	comparison[0,i]=countsSyn[i,1]/countsSyn[i,0]
	comparison[1,i]=countsSyn[i,2]/countsSyn[i,1]
	comparison[2,i]=countsSyn[i,2]/countsSyn[i,0]
for i in range(1):
	i1=1
	i2=0
	print(labels2[i1]+" vs "+labels2[i2])
	for j in range(3):
		print(labels[j])
		if j==0:
			j1=1
			j2=0
		elif j==1:
			j1=2
			j2=1
		else:
			j1=2
			j2=0
		print(chi2_contingency([[countsSyn[i1,j1],  countsSyn[i1,j2]],[countsSyn[i2,j1], countsSyn[i2,j2]]]))
		print(str(countsSyn[i1,j1])+" "+str(countsSyn[i1,j2])+" "+str(countsSyn[i2,j1])+" "+str(countsSyn[i2,j2]))
fig=barplotContext([comparison[:,0],comparison[:,1]],['Synonymous','Nonsynonymous'],pathSimu+"barplot_ratios_all_syn-nonsyn.pdf",noSave=True,xName='Comparison',colors=colors,labels=labels)
plt.text(0.0, 0.40, "p=0.00011", ha='center', va='bottom', color='black')
plt.plot([-0.2, -0.2, 0.2, 0.2], [0.39, 0.395, 0.395, 0.39], lw=1.5, c='black')
plt.text(1.0, 0.56, "p=0.029", ha='center', va='bottom', color='black')
plt.plot([0.8, 0.8, 1.2, 1.2], [0.55, 0.555, 0.555, 0.55], lw=1.5, c='black')
plt.text(2.0, 0.23, "p=3.6e-08", ha='center', va='bottom', color='black')
plt.plot([1.8, 1.8, 2.2, 2.2], [0.22, 0.225, 0.225, 0.22], lw=1.5, c='black')
fig.savefig(pathSimu+"barplot_ratios_all_syn-nonsyn.pdf")
#plt.hlines(countsUsyn[1,1]/countsUsyn[1,0], 0, 100, colors=colors[0], linestyles='dashed')
plt.close()


exit()



