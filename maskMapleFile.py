import sys
import os
from os import path
import argparse
import time

#©EMBL-European Bioinformatics Institues, 2023

#Mask entries of a MAPLE file filter out possibly error-prone positions.
#Example run command line: python3 maskMapleFile.py --maskFile problematic_sites_sarsCov2.vcf --input alignment.maple --output maskedAlignment.maple --minimumPos 266 --maximumPos 29673

#TODO allow option to remove masked positions from the alignment so to reduce runtime

parser = argparse.ArgumentParser(description='Translate fasta alignment into a MAPLE file.')
parser.add_argument('--maskFile',default="problematic_sites_sarsCov2.vcf", help='file containing the masking recommendations.')
parser.add_argument("--minimumPos",help="First genome position not to be masked - all position before this one will be masked.",  type=int, default=1)
parser.add_argument("--maximumPos",help="Last genome position not to be masked - all position after this one will be masked.",  type=int, default=float('inf'))
parser.add_argument("--input",default="", help="name of the input maple alignment file - it must have a reference at the beginning of the file.")
parser.add_argument("--output",default="", help="name of the output masked maple file.")
parser.add_argument("--overwrite", help="Overwrite previous MAPLE file with the same output name if already present.", action="store_true")
parser.add_argument("--reduceAlignment", help="Remove masked positions from the output MAPLE file - this makes inference from the alignment faster, but will change the length of the alignment and the interpretability of the alignment positions.", action="store_true")
args = parser.parse_args()

maskFile=args.maskFile
input=args.input
minimumPos=args.minimumPos
maximumPos=args.maximumPos
output=args.output
overwrite=args.overwrite
reduceAlignment=args.reduceAlignment

if not os.path.isfile(maskFile):
	print("ERROR "+maskFile+" does not exist, quitting maskMapleFile.py . Use option --maskFile to specify a valid file for masking recommendations.")
	exit()
if not os.path.isfile(input):
	print("ERROR input maple alignment file "+input+" not found, quitting maskMapleFile.py . Use option --input to specify a valid input maple file name to mask.")
	exit()
if os.path.isfile(output)  and (not overwrite):
	print("ERROR file "+output+" already exists, quitting maskMapleFile.py . Use option --overwrite if you want to overwirte files previously created.")
	exit()

alleles={"A":0,"C":1,"G":2,"T":3}
allelesList=["A","C","G","T"]
allelesLow={"a":0,"c":1,"g":2,"t":3}
allelesListLow=["a","c","g","t"]
ambiguities={"y":[0.0,1.0,0.0,1.0],"r":[1.0,0.0,1.0,0.0],"w":[1.0,0.0,0.0,1.0],"s":[0.0,1.0,1.0,0.0],"k":[0.0,0.0,1.0,1.0],"m":[1.0,1.0,0.0,0.0],"d":[1.0,0.0,1.0,1.0],"v":[1.0,1.0,1.0,0.0],"h":[1.0,1.0,0.0,1.0],"b":[0.0,1.0,1.0,1.0]}

start = time.time()

#read masking file
masks=[]
if minimumPos>1:
	masks.append((1,minimumPos-1))
fileM=open(maskFile)
line=fileM.readline()
lastLine=line
while line[0]=="#":
	lastLine=line
	line=fileM.readline()
linelist=lastLine.split()
posIndex=-1
filterIndex=-1
for i in range(len(linelist)):
	if linelist[i]=="POS":
		posIndex=i
	elif linelist[i]=="FILTER":
		filterIndex=i
if filterIndex==-1:
	print(lastLine)
	print("error, not found column FILTER in the input VCF masking file. Please follow the format in https://github.com/W-L/ProblematicSites_SARS-CoV2/blob/master/problematic_sites_sarsCov2.vcf ")
	exit()
if posIndex==-1:
	print(lastLine)
	print("error, not found column POS in the input VCF masking file. Please follow the format in https://github.com/W-L/ProblematicSites_SARS-CoV2/blob/master/problematic_sites_sarsCov2.vcf ")
	exit()
while line!="" and line[0]!="\n":
	linelist=line.split()
	pos=int(linelist[posIndex])
	if linelist[filterIndex]=="mask":
		if pos>=minimumPos and pos<=maximumPos:
			masks.append((pos,pos))
	line=fileM.readline()
fileM.close()

#Read and write reference
fileI=open(input)
fileO=open(output,"w")
line=fileI.readline()
fileO.write(line)
line=fileI.readline()
ref=""
while line!="" and line[0]!=">":
	ref+=line.replace("\n","")
	line=fileI.readline()
ref=ref.lower()
lRef=len(ref)

if lRef>maximumPos:
	masks.append((maximumPos+1,lRef))
masks.append((lRef+1,lRef+1))

if reduceAlignment:
	pos=1
	newRef=""
	for mask in masks:
		if mask[0]>pos and pos<=lRef:
			newRef+=ref[pos-1:mask[0]-1]
			pos=mask[1]+1
		else:
			if (mask[1]+1)>pos:
				pos=mask[1]+1
	if pos<=lRef:
		newRef+=ref[pos-1:]
	fileO.write(newRef+"\n")
else:
	fileO.write(ref+"\n")

#Read sequences and mask them
nSeqs=0
while line!="" and line!="\n":
	seqList=[]
	fileO.write(line)
	line=fileI.readline()
	pos=0
	indexMask=0
	lastMaskPos=0
	if reduceAlignment:
		cumulativeCountMasked=0
	readNewLine=True
	while line!="" and line!="\n" and line[0]!=">":
		linelist=line.split()
		if len(linelist)>2:
			entry=(linelist[0].lower(),int(linelist[1]),int(linelist[2]))
		elif len(linelist)<2:
			print("In input file "+input+" found line with only one column: \n"+line+"ERROR Please check for errors in the alignment format.")
			raise Exception("exit")
		else:
			entry=(linelist[0].lower(),int(linelist[1]))
		if ref[entry[1]-1]==entry[0] and entry[0]!="n" and entry[0]!="-":
			print("Mutation observed into reference nucleotide at position "+str(entry[1])+" , nucleotide "+entry[0]+". Wrong reference and/or diff file?")
			raise Exception("exit")
		
		if reduceAlignment:
			if entry[1]>pos:
				pos=entry[1]
		else:
			pos=entry[1]
		duration=1
		if len(entry)>2:
			duration=entry[2]
		lastPos=entry[1]+duration-1
		if reduceAlignment:
			if len(entry)>2:
				duration=entry[2]+(entry[1]-pos)
		#print("Entry:")
		#print(entry)
		#print("Mask:")
		#print(masks[indexMask])
		#print("cumulativeCountMasked "+str(cumulativeCountMasked)+" lastMaskPos "+str(lastMaskPos)+" lastPos "+str(lastPos)+" pos "+str(pos))
		while masks[indexMask][1]<pos or masks[indexMask][0]<lastMaskPos:
			#print("138 "+str(indexMask)+" "+str(lastMaskPos)+" "+str(pos))
			if masks[indexMask][1]>lastMaskPos:
				firstMaskPos=max(lastMaskPos+1,masks[indexMask][0])
				lastMaskPos=masks[indexMask][1]
				if reduceAlignment:
					cumulativeCountMasked+=(lastMaskPos+1-firstMaskPos)
				else:
					fileO.write("n"+"\t"+str(firstMaskPos)+"\t"+str(lastMaskPos+1-firstMaskPos)+"\n")
			indexMask+=1
			#print("148 "+str(indexMask)+" "+str(lastMaskPos)+" "+str(pos)+" "+str(firstMaskPos))	
			
		if lastPos>lastMaskPos:
			if pos<masks[indexMask][0]:
				if lastPos<masks[indexMask][0]:
					if pos>lastMaskPos:
						if len(entry)==2:
							if reduceAlignment:
								fileO.write(entry[0]+"\t"+str(entry[1]-cumulativeCountMasked)+"\n")
								#print("Writing "+entry[0]+"\t"+str(entry[1]-cumulativeCountMasked))
							else:
								fileO.write(entry[0]+"\t"+str(entry[1])+"\n")
						else:
							if reduceAlignment:
								fileO.write(entry[0]+"\t"+str(pos-cumulativeCountMasked)+"\t"+str(duration)+"\n")
								#print("Writing "+entry[0]+"\t"+str(pos-cumulativeCountMasked)+"\t"+str(duration))
							else:
								fileO.write(entry[0]+"\t"+str(entry[1])+"\t"+str(entry[2])+"\n")
					else:
						if reduceAlignment:
							firstPosToPrint=lastMaskPos+1-cumulativeCountMasked
							if (1+lastPos-firstPosToPrint)>0:
								fileO.write(entry[0]+"\t"+str(firstPosToPrint)+"\t"+str(1+lastPos-firstPosToPrint)+"\n")
								#print("Writing "+entry[0]+"\t"+str(lastMaskPos+1-cumulativeCountMasked)+"\t"+str(duration))
						else:
							firstMaskPos=lastMaskPos+1
							lastMaskPos=lastPos
							fileO.write("n"+"\t"+str(firstMaskPos)+"\t"+str(lastMaskPos+1-firstMaskPos)+"\n")
				else:
					if reduceAlignment:
						if pos>lastMaskPos:
							firstPosToPrint=pos-cumulativeCountMasked
						else:
							firstPosToPrint=lastMaskPos+1-cumulativeCountMasked
						while masks[indexMask][1]<lastPos:
							cumulativeCountMasked+=(masks[indexMask][1]+1-masks[indexMask][0])
							lastMaskPos=masks[indexMask][1]
							indexMask+=1
							#print("173 "+str(indexMask)+" "+str(masks[indexMask][0])+" "+str(masks[indexMask][1]))
						strange=False
						if masks[indexMask][0]<=lastPos:
							lastMaskPos=masks[indexMask][1]
							cumulativeCountMasked+=(lastPos+1-masks[indexMask][0])
						fileO.write("n"+"\t"+str(firstPosToPrint)+"\t"+str(lastPos+1-(cumulativeCountMasked+firstPosToPrint))+"\n")
						#print("Writing "+"n"+"\t"+str(firstPosToPrint)+"\t"+str(lastPos+1-(cumulativeCountMasked+firstPosToPrint)))
						if masks[indexMask][0]<=lastPos:
							cumulativeCountMasked+=(masks[indexMask][1]-lastPos)
							indexMask+=1
							#print("181 "+str(indexMask)+" "+str(masks[indexMask][0])+" "+str(masks[indexMask][1]))
					else:
						firstMaskPos=min(pos,masks[indexMask][0])
						firstMaskPos=max(firstMaskPos,lastMaskPos+1)
						lastMaskPos=max(masks[indexMask][1],lastPos)
						fileO.write("n"+"\t"+str(firstMaskPos)+"\t"+str(lastMaskPos+1-firstMaskPos)+"\n")
						indexMask+=1
						#print("188 "+str(indexMask)+" "+str(masks[indexMask][0])+" "+str(masks[indexMask][1]))
			else:
				firstMaskPos=max(masks[indexMask][0],lastMaskPos+1)
				if reduceAlignment:
					cumulativeCountMasked+=(masks[indexMask][1]+1-firstMaskPos)
					lastMaskPos=masks[indexMask][1]
					readNewLine=False
					pos=lastMaskPos+1
				else:
					lastMaskPos=max(masks[indexMask][1],lastPos)
					fileO.write("n"+"\t"+str(firstMaskPos)+"\t"+str(lastMaskPos+1-firstMaskPos)+"\n")
				indexMask+=1
				#print("200 "+str(indexMask)+" "+str(masks[indexMask][0])+" "+str(masks[indexMask][1]))
		if readNewLine:
			line=fileI.readline()
		else:
			readNewLine=True
	
	while masks[indexMask][1]<=lRef:
		if masks[indexMask][1]>lastMaskPos:
			firstMaskPos=max(masks[indexMask][0],lastMaskPos+1)
			lastMaskPos=masks[indexMask][1]
			if reduceAlignment:
				cumulativeCountMasked+=(lastMaskPos+1-firstMaskPos)
			else:
				fileO.write("n"+"\t"+str(firstMaskPos)+"\t"+str(lastMaskPos+1-firstMaskPos)+"\n")
		indexMask+=1
		#print("214 "+str(indexMask)+" "+str(masks[indexMask][0])+" "+str(masks[indexMask][1]))
	nSeqs+=1
fileI.close()
fileO.close()


time2 = time.time() - start
print("Time to mask maple file: "+str(time2))
print(str(nSeqs)+" sequences masked.")

exit()




