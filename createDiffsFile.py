import sys
import os
import math
import numpy as np
from os import path
import argparse
from Bio.Data import CodonTable
table = CodonTable.ambiguous_dna_by_id[1]
from Bio.Seq import _translate_str
import time

#©EMBL-European Bioinformatics Institues, 2021

#Translate fasta alignment into diff file.
#Example run command line: python3 createDiffsFile.py --path /pathToFolder/ --reference EPI_ISL_402124_lowercase.fasta --fasta 2021-03-31_unmasked.fa --output 2021-03-31_unmasked_differences.txt

parser = argparse.ArgumentParser(description='Translate fasta alignment into a diff file.')
parser.add_argument('--path',default="", help='path where to find and write files.')
parser.add_argument('--reference',default="EPI_ISL_402124_lowercase.fasta", help='name of the reference sequence file within the --path.')
parser.add_argument("--fasta",default="2021-03-31_unmasked.fa", help="name of the input fasta alignment file.")
parser.add_argument("--output",default="2021-03-31_unmasked_differences.txt", help="name of the output diff file.")
args = parser.parse_args()

pathSimu=args.path
reference=args.reference
fasta=args.fasta
output=args.output

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
	return ref


ref = collectReference(pathSimu+reference)
lRef=len(ref)


start = time.time()
#collect alignment and translate into diff file
fileI=open(pathSimu+fasta)
fileO=open(pathSimu+output,"w")
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
	if (nSeqs%10000)==0:
		print("Processes "+str(nSeqs)+" sequences")
fileI.close()
fileO.close()

time2 = time.time() - start
print("Time to convert alignment file: "+str(time2))
print(str(nSeqs)+" sequences converted.")

exit()



