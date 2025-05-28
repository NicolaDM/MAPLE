import sys
import os
import math
import argparse
import os.path
import gzip

#©EMBL-European Bioinformatics Institute, 2025


parser = argparse.ArgumentParser(description='Prepare Viridian genomes for running in MAPLE.')
parser.add_argument("--maskAlignment", help="Mask problematic positions from the alignments", action="store_true")
parser.add_argument("--removeOnlyShortCommonDeletions", help="Remove deletion entries from the MAPLE alignment.", action="store_true")
parser.add_argument("--summarizeQCfiles", help="Look at Viridian QC files and create a summary listing positions with heterozygosity and coverage variation along the genome. Use this option to specify how many cores have been used to parallelize this task",  type=int, default=0)
parser.add_argument("--coreQC",help="When option --summarizeQCfiles is used, this is the number of the particular core assigned to this run of the script.",  type=int, default=1)
parser.add_argument("--analyseQCsummaryFiles", help="Calculate statistics from QC summary files.", action="store_true")

args = parser.parse_args()
maskAlignment=args.maskAlignment
removeOnlyShortCommonDeletions=args.removeOnlyShortCommonDeletions
summarizeQCfiles=args.summarizeQCfiles
coreQC=args.coreQC
analyseQCsummaryFiles=args.analyseQCsummaryFiles



#positions to be masked
maskedPoss={}
maskedPoss[25202]=True
maskedPoss[21987]=True
maskedPoss[27507]=True
maskedPoss[8835]=True
maskedPoss[15521]=True
maskedPoss[26766]=True
maskedPoss[8008]=True
maskedPoss[8012]=True
maskedPoss[15510]=True
maskedPoss[17259]=True
maskedPoss[19413]=True
maskedPoss[22786]=True
maskedPoss[22882]=True
maskedPoss[23948]=True
maskedPoss[8826]=True
maskedPoss[8829]=True
maskedPoss[15854]=True
maskedPoss[19672]=True
maskedPoss[21650]=True
maskedPoss[23118]=True
maskedPoss[25296]=True
maskedPoss[25324]=True
maskedPoss[25336]=True
maskedPoss[29687]=True
maskedPoss[22026]=True
maskedPoss[22027]=True
maskedPoss[22028]=True
maskedPoss[22029]=True
maskedPoss[22030]=True
maskedPoss[22031]=True
maskedPoss[22032]=True
maskedPoss[22033]=True
maskedPoss[22034]=True
maskedPoss[22195]=True
maskedPoss[22197]=True
maskedPoss[22198]=True
maskedPoss[22202]=True
maskedPoss[22204]=True
maskedPoss[274]=True
maskedPoss[4321]=True
maskedPoss[26530]=True
maskedPoss[28245]=True
maskedPoss[28247]=True
maskedPoss[28249]=True
maskedPoss[28253]=True
maskedPoss[28251]=True
maskedPoss[28254]=True




# Read all Viridian QC files, and create summaries about the coverage level and heterozygosity in all the samples.
if summarizeQCfiles:
	createBashForQCprocess=False
	if createBashForQCprocess:
		file=open("QCsummaryBash.sh","w")
		file.write("for i in $(seq 1 "+str(summarizeQCfiles)+")\n do\n\n")
		file.write("sbatch -J VirQC -t 100:00:00 --mem=10G -o QCsummaryBash\"$i\"_console_output.txt -e QCsummaryBash\"$i\"_console_error.txt "
								+" --wrap=\"pypy3.10 processMartinHuntData.py --summarizeQCfiles "+str(summarizeQCfiles)+" --coreQC \"$i\" \" \n ")
		file.write("done\n\n")
		file.close()
		exit()

	alFile=open("alignment.maple")
	sampleNum=1
	lineAl=alFile.readline()
	lineAl=alFile.readline()
	lineAl=alFile.readline()
	oFile=open("QCsummary"+str(coreQC)+".txt","w")
	while lineAl!="" and lineAl!="\n":
		if lineAl[0]==">":
			name=lineAl[1:-1]
			if ((sampleNum%summarizeQCfiles)+1)==coreQC:
				if name[0]=="E":
					try:
						file=gzip.open("Vdn_all_ena/Reads/"+name[0]+"/"+name[1:6]+"/"+name[6:8]+"/"+name[8:10]+"/vdn.v1.0.0/qc.tsv.gz",'rt')
					except:
						try:
							file=gzip.open("Vdn_all_ena/Reads/"+name[0]+"/"+name[1:7]+"/"+name[7:9]+"/"+name[9:11]+"/vdn.v1.0.0/qc.tsv.gz",'rt')
						except:
							print(name)
							lineAl=alFile.readline()
							continue
				elif name[0]=="S":
					file=gzip.open("Vdn_all_ena/Reads/"+name[0]+"/"+name[1:7]+"/"+name[7:9]+"/"+name[9:11]+"/vdn.v1.0.0/qc.tsv.gz",'rt')
				elif name[0]=="D":
					file=gzip.open("Vdn_all_ena/Reads/"+name[0]+"/"+name[1:5]+"/"+name[5:7]+"/"+name[7:9]+"/vdn.v1.0.0/qc.tsv.gz",'rt')
				else:
					print(name)
					lineAl=alFile.readline()
					continue
				oFile.write(lineAl)
				line=file.readline()
				line=file.readline()
				covState=-1
				numPosPrinted=0
				while line!="" and line!="\n":
					linelist=line.split()
					try:
						pos=int(linelist[0])
					except:
						print(line)
					try:
						cov=int(linelist[9])
					except:
						cov=0
					if cov<20:
						newCovState=0
					elif cov<100:
						newCovState=1
					else:
						newCovState=2
					if newCovState!=covState:
						oFile.write(str(pos)+"\t")
						if cov<20:
							oFile.write("cov<20\n")
						elif cov<100:
							oFile.write("cov<100\n")
						else:
							oFile.write("cov>=100\n")
						covState=newCovState
					if newCovState:
						numNucs=0
						numNucs2=0
						counts=[0,0,0,0]
						for i in range(4):
							counts[i]+=int(linelist[11+2*i])
							counts[i]+=int(linelist[12+2*i])
							if counts[i]>9 and (float(counts[i])/cov)>0.05:
								numNucs+=1
								if counts[i]>29 and (float(counts[i])/cov)>0.15:
									numNucs2+=1
						if numNucs>1:
							numPosPrinted+=1
							oFile.write(line)
						if numPosPrinted>50:
							break
					line=file.readline()
				
				file.close()

			if (sampleNum%1000)==0:
				print(sampleNum)
			sampleNum+=1
		lineAl=alFile.readline()
	alFile.close()
	oFile.close()




# Go through the QC summary files and find overall heterozygosity and coverage distributions so to find possible thresholds to remove samples that might be contaminated or low coverage.
# filter out samples with :
# >2 positions with Het>20%;
# >7 positions with Het>10%;
# >30 positions with Het>5%;
# >1500 positions with cov<=20;
# >2500 positions with cov<=100;
if analyseQCsummaryFiles:
	thresholdFreqs=[0.05,0.1,0.2]
	thresholdHetCov=[9,19,39]
	numPositionsThresholds=[1500,2500]
	maxNumHet=[30,7,2]
	barplotHet=[[0]*52,[0]*52,[0]*52]
	barplotCov=[[0]*102,[0]*102]
	filteredOutSamples={}
	numMaskedForHet=[0,0,0]
	numMaskedForCov=[0,0]
	for iFile in range(100):
		print("File QCsummary"+str(iFile+1)+".txt")
		file=open("QCsummary"+str(iFile+1)+".txt")
		line=file.readline()
		while line!="" and line!="\n":
			numHet=[0,0,0]
			numPosBelow=[0,0]
			currentPos=1
			currentCov=0
			name=line[1:-1]
			line=file.readline()
			while line!="" and line!="\n" and line[0]!=">":
				linelist=line.split()
				pos=int(linelist[0])
				if len(linelist[1])>1:
					pos=int(linelist[0])
					if currentCov<2:
						numPosBelow[currentCov]+=pos-currentPos
						if currentCov<1:
							numPosBelow[1]+=pos-currentPos
					if linelist[1]=="cov<20":
						currentCov=0
					elif linelist[1]=="cov>=100":
						currentCov=2
					else:
						currentCov=1
					currentPos=pos
				else:
					if currentCov:
						if not (pos in maskedPoss):
							try:
								cov=int(linelist[9])
							except:
								cov=0
							if cov:
								numNucs=[0,0,0]
								counts=[0,0,0,0]
								for i in range(4):
									counts[i]+=int(linelist[11+2*i])
									counts[i]+=int(linelist[12+2*i])
									for j in range(len(thresholdFreqs)):
										if counts[i]>thresholdHetCov[j] and (float(counts[i])/cov)>thresholdFreqs[j]:
											numNucs[j]+=1
								for j in range(len(thresholdFreqs)):
									if numNucs[j]>1:
										numHet[j]+=1
				line=file.readline()
			if currentCov<2:
				numPosBelow[currentCov]+=29904-currentPos
			if numHet[0]>50:
				barplotHet[0][51]+=1
				filteredOutSamples[name]=True
			else:
				for j in range(len(thresholdFreqs)):
					barplotHet[j][numHet[j]]+=1
					if numHet[j]>maxNumHet[j]:
						filteredOutSamples[name]=True
						numMaskedForHet[j]+=1
				for j in range(len(numPosBelow)):
					if numPosBelow[j]>numPositionsThresholds[j]:
						filteredOutSamples[name]=True
						numMaskedForCov[j]+=1
					barplotCov[j][int(numPosBelow[j]/300)]+=1

		file.close()
	print("barplots Het:")
	for j in range(len(thresholdFreqs)):
		print(barplotHet[j])
	print("barplots Cov:")
	for j in range(len(numPosBelow)):
		print(barplotCov[j])
	print("Number of samples removed for Het")
	print(numMaskedForHet)
	print("Number of samples removed for Cov")
	print(numMaskedForCov)

	#now filter out bad samples
	file=open("alignment.maple")
	fileO=open("alignment_filtered.maple","w")
	line=file.readline()
	fileO.write(line)
	line=file.readline()
	fileO.write(line)
	line=file.readline()
	while line!="" and line!="\n":
		name=line[1:-1]
		if not (name in filteredOutSamples):
			fileO.write(line)
			line=file.readline()
			while line!="" and line!="\n" and line[0]!=">":
				fileO.write(line)
				line=file.readline()
		else:
			line=file.readline()
			while line!="" and line!="\n" and line[0]!=">":
				line=file.readline()
	file.close()
	fileO.close()
	exit()







#create alignment without deletions - deletions seem to be negatively affecting phylogenetic inference due to sparse samples with errors, causing these sites to have heavy uncertainty and artificial ancestral mutations.
if removeOnlyShortCommonDeletions:
	file=open("alignment_filtered.maple")
	createSortedListOfDeletionAbundances=False
	if createSortedListOfDeletionAbundances:
		line=file.readline()
		line=file.readline()
		line=file.readline()
		deletionsDict={}
		while line!="" and line!="\n":
			if line[0]=="-":
				linelist=line.split()
				if len(linelist)>2:
					lineCode=linelist[1]+"+"+linelist[2]
				else:
					lineCode=linelist[1]+"+1"
				if lineCode in deletionsDict:
					deletionsDict[lineCode]+=1
				else:
					deletionsDict[lineCode]=1
			line=file.readline()
		sortedDeletions=sorted(deletionsDict.items(), key=lambda item: item[1])
		print(sortedDeletions[-1000:])
		exit()

	fileO=open("alignment_filtered_noShortDeletions.maple","w")
	line=file.readline()
	fileO.write(line)
	line=file.readline()
	fileO.write(line)
	line=file.readline()
	while line!="" and line!="\n":
		if line[0]!="-":
			fileO.write(line)
		else:
			linelist=line.split()
			if (len(linelist)<3) or (int(linelist[2])>30):
				fileO.write(line)
		line=file.readline()
	file.close()
	fileO.close()
	exit()






#mask positions of recurrent errors
if maskAlignment:
	alleles={"A":0,"C":1,"G":2,"T":3,"a":0,"c":1,"g":2,"t":3}
	ambiguities={"y":[0.0,1.0,0.0,1.0],"r":[1.0,0.0,1.0,0.0],"w":[1.0,0.0,0.0,1.0],"s":[0.0,1.0,1.0,0.0],"k":[0.0,0.0,1.0,1.0],"m":[1.0,1.0,0.0,0.0],"d":[1.0,0.0,1.0,1.0],"v":[1.0,1.0,1.0,0.0],"h":[1.0,1.0,0.0,1.0],"b":[0.0,1.0,1.0,1.0]}
	samples={}
	file=open("alignment_filtered_noShortDeletions.maple")
	line=file.readline()
	line=file.readline()
	line=file.readline()
	name=line.replace(">","").replace("\n","")
	samples[name]=True
	line=file.readline()
	while line!="":
		while line[0]!=">":
			line=file.readline()
			if line=="" or line=="\n":
				break
		if line!="" and line!="\n":
			name=line.replace(">","").replace("\n","")
			samples[name]=True
		line=file.readline()
	file.close()
	
	file=open("alignment_filtered_noShortDeletions.maple")
	fileO=open("alignment_filtered_noShortDeletions_masked.maple","w")
	line=file.readline()
	fileO.write(line)
	line=file.readline()
	fileO.write(line)
	line=file.readline()
	name=line.replace(">","").replace("\n","")
	print(name)
	nameInSamples=False
	if name in samples:
		fileO.write(line)
		nameInSamples=True
	line=file.readline()
	while line!="":
		while line[0]!=">":
			if nameInSamples:
				linelist=line.split()
				pos=int(linelist[1])
				if len(linelist)==3:
					numNs=int(linelist[2])
					if (pos in maskedPoss) and numNs==1:
						pass
					else:
						fileO.write(line)
				else:
					if not (pos in maskedPoss):
						fileO.write(line)
			line=file.readline()
			if line=="" or line=="\n":
				break
		if line!="" and line!="\n":
			name=line.replace(">","").replace("\n","")
			if name in samples:
				fileO.write(line)
				nameInSamples=True
			else:
				nameInSamples=False
		line=file.readline()
	file.close()
	fileO.close()
	exit()










		

exit()



