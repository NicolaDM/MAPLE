import sys
import os
import math
import numpy as np
import time

#Some functions that can be useful in maximum likelihood inference or in estimating an initial tree.
#For example, one is for comparing two samples and deciding if one of them can be removed as 
# strictly less informative than the other, and attached to the other sample after maximum
# likelihood inference. This is a generalization of the approach of removing sequences identical to
# other sequences in the alignment, and can lead to substantial time saving with SARS-CoV-2.

alleles={"A":0,"C":1,"G":2,"T":3}
allelesList=["A","C","G","T"]
allelesLow={"a":0,"c":1,"g":2,"t":3}
allelesListLow=["a","c","g","t"]
ambiguities={"y":[0.0,1.0,0.0,1.0],"r":[1.0,0.0,1.0,0.0],"w":[1.0,0.0,0.0,1.0],"s":[0.0,1.0,1.0,0.0],"k":[0.0,0.0,1.0,1.0],"m":[1.0,1.0,0.0,0.0],"d":[1.0,0.0,1.0,1.0],"v":[1.0,1.0,1.0,0.0],"h":[1.0,1.0,0.0,1.0],"b":[0.0,1.0,1.0,1.0]}


range4=range(4)


#create list with one entry per sample, and sort it based on the
# distance of the reference from the sample. Distance is calculated so that genomes with fewer 
# unsequenced positions have smaller distances.
def distancesFromRef(data):
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
		sampleDistances.append((diffNum*(lRef+1)+lRef-comparisons,sample))

	from operator import itemgetter
	print("Now doing sorting")
	sampleDistances.sort(key=itemgetter(0))
	return sampleDistances

#function to check is one sequence is inferior to another.
#Inferior sequences can be removed when doing maximum likelihood inference, 
# and then attached as children of their superior samples.
# For each comparison, returns 0 when the 2 sequences are not comparable, 
# otherwise returns 1 if the first is more informative (superior) 
# or if they are identical, and 2 otherwise.
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

exit()






















