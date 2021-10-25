import sys
import os
import math
import numpy as np
import random
from os import path
import argparse
from Bio.Data import CodonTable
table = CodonTable.ambiguous_dna_by_id[1]
from Bio.Seq import _translate_str
import time
from ete3 import Tree

#©EMBL-European Bioinformatics Institues, 2021

#Efficiently calculating phylogenetic likelihoods.
#Example run command line: python3 calculateFastLK.py --path /pathToFolder/ --reference EPI_ISL_402124_lowercase.fasta

parser = argparse.ArgumentParser(description='Run estimation of mutation rates and synonymous site selection from SARS-CoV-2 data, and generate plots.')
parser.add_argument('--path',default="", help='path where to find files and plot results.')
parser.add_argument('--reference',default="EPI_ISL_402124_lowercase.fasta", help='name of the reference sequence file within the --path.')
parser.add_argument("--diffFile",default="2021-03-31_unmasked_differences_reduced.txt", help="input diff (concise alignment) file.")
parser.add_argument("--treeFile",default="global.txt", help="input tree (newick) file.")
parser.add_argument("--onlyNambiguities", help="Treat all ambiguities as N (total missing information).", action="store_true")
parser.add_argument("--thresholdProb",help="relative probability threshold used to ignore possible states with very low probabilities.",  type=float, default=0.000001)
parser.add_argument("--verbose", help="Print to screen a lot of stuff.", action="store_true")
args = parser.parse_args()

pathSimu=args.path
reference=args.reference
diffFile=args.diffFile
treeFile=args.treeFile
onlyNambiguities=args.onlyNambiguities
thresholdProb=args.thresholdProb
verbose=args.verbose

alleles={"A":0,"C":1,"G":2,"T":3}
allelesList=["A","C","G","T"]
allelesLow={"a":0,"c":1,"g":2,"t":3}
allelesListLow=["a","c","g","t"]
ambiguities={"y":[0.0,1.0,0.0,1.0],"r":[1.0,0.0,1.0,0.0],"w":[1.0,0.0,0.0,1.0],"s":[0.0,1.0,1.0,0.0],"k":[0.0,0.0,1.0,1.0],"m":[1.0,1.0,0.0,0.0],"d":[1.0,0.0,1.0,1.0],"v":[1.0,1.0,1.0,0.0],"h":[1.0,1.0,0.0,1.0],"b":[0.0,1.0,1.0,1.0]}

#collect reference
def collectReferenceAndFreqs(fileName):
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
	print(cumulativeBases[-1])
	#collect base frequencies at the reference, posssibly to be used as root frequencies
	rootFreqs=np.zeros(4)
	rootFreqsLog=np.zeros(4)
	for i in range(4):
		rootFreqs[i]=cumulativeBases[-1][i]/float(lRef)
		rootFreqsLog[i]=math.log(rootFreqs[i])
	print(rootFreqs)
	print(rootFreqsLog)
	return ref, cumulativeBases, rootFreqs, rootFreqsLog


ref, cumulativeBases, rootFreqs, rootFreqsLog = collectReferenceAndFreqs(pathSimu+reference)
lRef=len(ref)


#read the input concise alignment (diff) file
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


#read sequence data from file
data=readConciseAlignment(pathSimu+diffFile)

#list of sample/taxon names
samples=data.keys()
print(str(len(samples))+" sequences in the concise DNA data file")

start = time.time()
phylo = Tree(pathSimu+treeFile)
time2 = time.time() - start
print("Time to read newick tree: "+str(time2))
print(str(len(phylo))+" sequences in the tree.")


range4=range(4)
#function to normalize the partial likelihood in an "O" genome list entry, and, 
# in case the partial likelihoods are all condensed in one entry,
# redefine the entry as non-O 
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

#Shorten genome list by merging consecutive entries if they are both of R type 
# and if their time value is (almost) identical
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

#define genome list (partial likelihoods) for a leaf node, 
# given its data and given its branch length
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
	return probVect

mutMatrix=[[0.0,0.04,0.3,0.1],[0.04,0.0,0.02,1.0],[0.3,0.02,0.0,1.0],[0.1,1.0,1.0,0.0]]
rootFreqs=[0.25,0.25,0.25,0.25]
rootFreqsLog=[math.log(0.25),math.log(0.25),math.log(0.25),math.log(0.25)]
#0.04,0.3,0.1,0.02,1.0,1.0
#preliminary nuc mutation rate matrix, taken from De Maio et al 2021 GBE
#mutMatrix=[[0.0,0.039,0.310,0.123],[0.140,0.0,0.022,3.028],[0.747,0.113,0.0,2.953],[0.056,0.261,0.036,0.0]]
mutMatrix[0][0]=-(mutMatrix[0][1]+mutMatrix[0][2]+mutMatrix[0][3])
mutMatrix[1][1]=-(mutMatrix[1][0]+mutMatrix[1][2]+mutMatrix[1][3])
mutMatrix[2][2]=-(mutMatrix[2][0]+mutMatrix[2][1]+mutMatrix[2][3])
mutMatrix[3][3]=-(mutMatrix[3][0]+mutMatrix[3][1]+mutMatrix[3][2])
thresholdProb2=thresholdProb*thresholdProb


#Calculate likelihood with fast approximate approach.
#This function is called recursively (parents call it on children nodes).
#It could probably be made faster and less memory-intensive by making it non-recursive.
def pruningFast(node,mutMatrix,data,ref, cumulativeBases):
	children=node.children
	bLen=node.dist
	#cumulative contribution of the subtree to the tree likelihood (excluding probs in the vector probVect)
	cumulPartLk=0.0
	#case in which the considered node is a leaf node
	if len(children)==0:
		#if the leaf name is not in the diff file, then assign an empty sequence to it
		if not (node.name in data):
			probVect=[["N",1,lRef,bLen]]
			if verbose:
				print("\n"+"Terminal node not in diff file, name: "+node.name)
			return cumulPartLk, probVect
		diffs=data[node.name]
		pos=1
		probVect=probVectTerminalNode(diffs,bLen)
		node.cumulPartLk=cumulPartLk
		node.probVect=probVect
		if verbose:
			print("\n"+"Terminal node "+node.name+" "+str(node.dist))
			print(diffs)
			print(probVect)
		return cumulPartLk, probVect
	
	#internal node
	else:
		#extracting values for first child
		cumulPartLk1, probVect1= pruningFast(children[0],mutMatrix,data,ref, cumulativeBases)
		cumulPartLk+=cumulPartLk1
		if verbose:
			node.name=""
			for c in range(len(node.children)):
				if c!=0:
					node.name+="-"
				node.name+=node.children[c].name
			print("\n"+"Internal node "+node.name)
			print("First child "+children[0].name)
			print(probVect1)
			print(cumulPartLk1)
		#now iterating over node children. Here I don't assume a bifurcating tree, so one can have a polytomy, 
		# in which case I merge with one child at the time.
		for c in range(len(children)-1):
			#extract values for next child
			newChild=children[c+1]
			cumulPartLk2, probVect2= pruningFast(newChild,mutMatrix,data,ref, cumulativeBases)
			cumulPartLk+=cumulPartLk2
			if verbose:
				print("Next child "+newChild.name)
				print(probVect2)
				print(cumulPartLk2)

			#now traverse the two genome lists and create a new probVect while also updating the cumulPartLk
			#new genome list
			probVect=[]
			#Here I initialize the indeces for traversing the two genome lists (for the two children)
			# simultaneously. "pos" and "end" refer to the first and last position of the intersection
			# of the two entries (from the two genome lists) currently considered.
			# "indexEntry1" refers to the index of the current entry considered within the first genome list 
			# (the one from the first child). "pos1" and "end1" refer to start and end position of the 
			# currently considered entry of the first genome list. "pos2", "end2" and "indexEntry2" 
			# do the same for the genome list of the second child.
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
				#Here I consider a single entry of the new genome list, originating from
				# the intersection of "entry1" and "entry2".
				totLen=entry1[3]+entry2[3]
				if verbose:
					print("pos "+str(pos)+" length "+str(length)+" end "+str(end)+" indexEntry1 "+str(indexEntry1)+" indexEntry2 "+str(indexEntry2))
				#If one or the other entry is "N", then only the length term needs updating (in which case the update is done later, not here).
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
				#case the first child is in state "R"
				elif entry1[0]=="R":
					if entry2[0]=="R":
						probVect.append(["R",pos,length,0.0])
						#Now update partial likelihood
						for i in range4:
							cumulPartLk+=mutMatrix[i][i]*totLen*(cumulativeBases[end][i]-cumulativeBases[pos-1][i])
					elif entry2[0]=="O":
						i1=allelesLow[ref[pos-1]]
						newVec=np.ones(4)
						if entry2[3]+entry1[3]<thresholdProb2:
							#I am doing this bit here because in my example there are some discordances between
							# phylogeny and data, since I was not the one who estimated the phylogeny, and since
							# who estimated the phylogeny masked some of the data.
							# For this reason I am assigning a minimum brnach length.
							# However, in normal circumstances one would skip this step and just consider 0-length branches as such.
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
						# Check if we can convert the "O" entry into something else 
						# (in this case only "R" would be possible, so maybe one could do this check more efficiently)
						state, newVec, maxP =simplfy(newVec,ref[pos-1])
						probVect.append([state,pos,1,0.0,newVec])
						cumulPartLk+=math.log(maxP)
					else:
						#Case in which the second child is in A, C, G or T state.
						i2=alleles[entry2[0]]
						i1=allelesLow[ref[pos-1]]
						newVec=np.ones(4)
						if entry2[3]+entry1[3]<thresholdProb2:
							#This is optional, same as before, it's just to avoid potential errors which
							# however normally shouldn't happen, only if the data and the tree 
							# don't fit each other
							entry2[3]=thresholdProb
							entry1[3]=thresholdProb
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
				#Fist child is in state "O"
				elif entry1[0]=="O":
					newVec=np.ones(4)
					if entry2[3]+entry1[3]<thresholdProb2:
							#Again, optional
							entry2[3]=thresholdProb
							entry1[3]=thresholdProb
					if entry2[0]=="O":
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
					#Case in which the seond child is either in R, A, C, or T state
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
					#Normalize partial likelihoods and possibly replace
					# "O" state with R, A, C, G or T
					state, newVec, maxP =simplfy(newVec,ref[pos-1])
					probVect.append([state,pos,1,0.0,newVec])
					if verbose:
						print("Simplified O")
						print(newVec)
						print(state)
						print(maxP)
					cumulPartLk+=math.log(maxP)
				#entry1 is A, C, G, or T
				else:
					if entry2[0]==entry1[0]:
						i2=alleles[entry2[0]]
						probVect.append([entry1[0],pos,1,0.0])
						#Now update partial likelihood
						cumulPartLk+=mutMatrix[i2][i2]*totLen
					else:
						if entry2[3]+entry1[3]<thresholdProb2:
							#Again, possibly not needed/wanted in normal circumstances
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
						#Case second child is R, A, C, G or T, and second child state is not equal to first child state
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

				#update pos, end, etc so that we can move to the next intersection between 
				# entries of the genome lists of the two children.
				if verbose:
					print(pos)
					print(length)
				pos+=length
				if verbose:
					print("New values:")
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
		return cumulPartLk, probVect1

#once reached the root, finalize the likelihood by using the root partial likelihoods and the root frequencies estimated from the reference genome
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
		# I am doing this second case for the scenario in which 
		# mutation rates are not stationary, and so the frequencies might change over time.
		# This is however probably just a waste of effort, and it would make sense to just
		# approximate every case entry[3]==0.0
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


	
#calculate likelihood of phylogeny
start = time.time()
cumulPartLk, probVect = pruningFast(phylo,mutMatrix,data,ref, cumulativeBases)
time2 = time.time() - start
print("Fast likelihood finished")
print(cumulPartLk)
print(probVect)
#now use root frequencies to finalize the likelihood.
cumulPartLk+=finalLK(probVect,ref,rootFreqsLog,rootFreqs)
print("Final likelihood after finalization:")
print(cumulPartLk)
print("Time taken "+str(time2))

exit()


