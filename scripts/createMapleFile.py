import sys
import os
from os import path
import argparse
import time

#©EMBL-European Bioinformatics Institues, 2021

#Translate fasta alignment into a MAPLE format file.
#Example run command line: python3 createDiffsFile.py --path /pathToFolder/ --reference EPI_ISL_402124_lowercase.fasta --fasta 2021-03-31_unmasked.fa --output 2021-03-31_unmasked_differences.txt

#TODO use consensus sequence when no reference is provided.
#TODO create new MAPLE format output with reference first.

parser = argparse.ArgumentParser(description='Translate fasta alignment into a MAPLE file.')
parser.add_argument('--path',default="", help='path where to find and write files.')
parser.add_argument('--reference',default="", help='name of the reference sequence file within the --path. By default creates a new reference from the input alignment consensus.')
parser.add_argument("--fasta",default="2021-03-31_unmasked.fa", help="name of the input fasta alignment file.")
parser.add_argument("--output",default="2021-03-31_unmasked_differences.txt", help="name of the output diff file.")
parser.add_argument("--overwrite", help="Overwrite previous MAPLE file with the same output name if already present.", action="store_true")
args = parser.parse_args()

pathSimu=args.path
reference=args.reference
fasta=args.fasta
output=args.output
overwrite=args.overwrite

if not os.path.isdir(pathSimu):
	print("ERROR path "+pathSimu+" does not exist, quitting createMapleFile.py . Use option --path to specify a valid path for input and output files location.")
	exit()
if not os.path.isfile(pathSimu+fasta):
	print("ERROR input file in fasta format "+pathSimu+fasta+" not found, quitting createMapleFile.py . Use option --fasta to specify a valid input fasta file name to be found within the folder spcefified by option --path .")
	exit()
if reference!="" and ( not os.path.isfile(pathSimu+reference)):
	print("ERROR input reference fasta file "+pathSimu+reference+" not found, quitting createMapleFile.py . You can skip option --reference and this script will use the alignmnt consensus as a reference instead .")
	exit()
if os.path.isfile(pathSimu+output)  and (not overwrite):
	print("ERROR file "+pathSimu+output+" already exists, quitting createMapleFile.py . Use option --overwrite if you want to overwirte files previously created.")
	exit()



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
	ref=ref.lower()
	return ref

#extract consensus sequence from an alignment
def extractConsensus(fileName):
	file=open(fileName)
	line=file.readline()
	counts=[]
	while line!="":
		seq=""
		name=line.replace(">","").replace("\n","")
		line=file.readline()
		while line!="" and line!="\n" and line[0]!=">":
			seq+=line.replace("\n","").lower()
			line=file.readline()
		if not counts:
			for i in range(len(seq)):
				counts.append([0,0,0,0])
		if len(seq)!=len(counts):
			print("ERROR sequence of sample "+name+" has length "+str(len(seq))+" instead of "+str(len(counts))+" of the first sequence in the file. Exiting createMapleFile.py .")
			exit()
		for i in range(len(seq)):
			if seq[i] in allelesListLow:
				counts[i][allelesLow[seq[i]]]+=1
		while line=="\n":
			line=file.readline()
	consensus=""
	for i in range(len(counts)):
		maxI=0
		maxV=0
		for j in range(4):
			if counts[i][j]>maxV:
				maxI=j
				maxV=counts[i][j]
		if maxV>0:
			consensus+=allelesListLow[maxI]
		else:
			print("WARNING no nucleotide observed at position "+str(i+1)+" of the alignment. Consensus is assigned as n, which can create problems down the line if more sequence will be analysed with the same reference." )
			consensus+="n"
	return consensus

if reference!="":
	ref = collectReference(pathSimu+reference)
else:
	ref = extractConsensus(pathSimu+fasta)
lRef=len(ref)


start = time.time()
#collect alignment and translate into diff file
fileI=open(pathSimu+fasta)
fileO=open(pathSimu+output,"w")
fileO.write(">reference\n"+ref+"\n")
line=fileI.readline()
nSeqs=0
while line!="":
	while line=="\n":
		line=fileI.readline()
	if line=="":
		break
	nSeqs+=1
	seq=""
	name=line.replace(">","").replace("\n","")
	fileO.write(line)
	line=fileI.readline()
	while line!="" and line!="\n" and line[0]!=">":
		seq+=line.replace("\n","")
		line=fileI.readline()
	if len(seq)!=lRef:
		print("Seq "+name+" has length "+str(len(seq))+" while reference is "+str(lRef))
		exit()
	# state 0=ref; 1=N; 2=-; 
	state=0
	seqList=[]
	length=0
	seq=seq.lower()
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
		if seq[i]=="n" and state!=1:
			length=1
			state=1
		elif seq[i]=="-" and state!=2:
			length=1
			state=2
		elif seq[i]!=ref[i] and seq[i]!="-" and seq[i]!="n":
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




