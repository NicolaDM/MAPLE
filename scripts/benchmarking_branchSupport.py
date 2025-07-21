import sys
import os.path
from time import time
import argparse
import random
import time

#©EMBL-European Bioinformatics Institute, 2023-2025

# create data:
# pypy3 benchmarking_branchSupport.py --createBashScript
# sh createSubsampleInputFiles.sh

#submit MAPLE (ignoring UShER for now)
# pypy3 benchmarking_branchSupport.py --createBashScript
# sh submitMAPLE.sh
# sh submitSPRTA.sh

#submit IQtree support measures using MAPLE trees as starting topologies
# pypy3 benchmarking_branchSupport.py --createBashScript
# sh submitIQTREE2.sh

#transform IQTREE trees into mat
# pypy3 benchmarking_branchSupport.py --createBashScript
# sh submitMatEstimation.sh

# collect results about branch support
# pypy3 benchmarking_branchSupport.py --compareMats --numSamples 1000
# etc

# then run figure generation on local computer:
# python3 benchmarking_branchSupport.py --runFigureGeneration
# python3 benchmarking_branchSupport.py --runFigureGeneration --scenario 0
# source sklearn-env_p3/bin/activate
# python3.12 benchmarking_branchSupport.py --runFigureGeneration --scenario 1

# check with an example which nodes are incorrect, why, and what support they have:
# pypy3 benchmarking_branchSupport.py --compareMats --compareMats_example

parser = argparse.ArgumentParser(description='Create files for the benchmarking of MAPLE, both input files and cluster shell scripts.')
parser.add_argument('--pathToSimulationFolder',default="", help='path to the folder for the new simulations.')
parser.add_argument('--pathToErrorSimulationFolder',default="", help='path to the folder for the error estimation simulations.')
parser.add_argument('--inputRealTree',default="public-latest.all.nwk", help='path to the (real) input tree for phastSim simulations.')
parser.add_argument('--inputSimulationReference',default="MN908947.3.fasta", help='path to the reference genome to be used for phastSim simulations.')

parser.add_argument("--createBashScript", help="Create bash script to run the file creation on the cluster in parallel.", action="store_true")

parser.add_argument("--createInputData", help="Subsample real and simulated data to create input datasets for benchmarking.", action="store_true")
parser.add_argument("--subSampleNum",help="Number of subsamples to extract.",  type=int, default=20)
parser.add_argument("--repeat",help="Which repeat to simulate, typically between 1-10. Only considered for option --createInputData .",  type=int, default=1)
parser.add_argument("--numSamples",help="Number of samples to consider when running different estimation methods, or when checking the results; typically between 1000 and 1000000. if 0 (default) consider the whole range.",  type=int, default=0)
parser.add_argument("--scenario",help="Which scenario to subsample, between 0 (no Ns) or 1 (with Ns).",  type=int, default=0)

parser.add_argument("--compareMats", help="Look into branch support results by comparing MATs.", action="store_true")
parser.add_argument("--compareMats_example", help="Look into branch support results by comparing MATs, this time just comparing the first estimated tree with the simulated one and exiting.", action="store_true")
parser.add_argument("--runFigureGeneration", help="Create the figures from the benchmarking experiment - it needs to be used from somewhere with matplotlib.", action="store_true")
parser.add_argument("--onlyNePlots", help="Look into branch support results by comparing MATs.", action="store_true")
args = parser.parse_args()

#subsample datasets
subSampleNums=[1000, 2000, 5000,10000,20000,50000,100000,200000]
subSampleFasta=[1000, 2000, 5000,10000,20000]

#define binary tree structure for quickly finding samples 
class BinarySearchTree:
	def __init__(self):
		self.root = None
	def __iter__(self):
		return self.root.__iter__()
	def put(self,key):
		if self.root:
			self._put(key,self.root)
		else:
			self.root = TreeNodeSearch(key)
	def _put(self,key,currentNode):
		if key < currentNode.key:
			if currentNode.hasLeftChild():
				self._put(key,currentNode.leftChild)
			else:
				currentNode.leftChild = TreeNodeSearch(key)
		else:
			if currentNode.hasRightChild():
				self._put(key,currentNode.rightChild)
			else:
				currentNode.rightChild = TreeNodeSearch(key)
	def get(self,key):
		if self.root:
			return self._get(key,self.root)
		else:
			return False

	def _get(self,key,currentNode):
		if not currentNode:
			return False
		elif currentNode.key == key:
			return True
		elif key < currentNode.key:
			return self._get(key,currentNode.leftChild)
		else:
			return self._get(key,currentNode.rightChild)

	def __getitem__(self,key):
		return self.get(key)

class TreeNodeSearch:
	def __init__(self,key,left=None,right=None):
		self.key = key
		self.leftChild = left
		self.rightChild = right
	def hasLeftChild(self):
		return self.leftChild
	def hasRightChild(self):
		return self.rightChild

#Class defining nodes of the tree
class Tree(object):
	def __init__(self, name='', children=None, dist=1.0):
		if name!='':
			self.name = name
		self.dist = dist
		self.children = []
		self.up=None
		self.dirty=True
		if children is not None:
			for child in children:
				self.add_child(child)
	def __repr__(self):
		try:
			return str(self.name)
		except AttributeError:
			try:
				return self.name
			except AttributeError:
				return ""
	def add_child(self, node):
		assert isinstance(node, Tree)
		self.children.append(node)


#function to read input mutation-annotated newick string and create a dictionary of mutation events
def readNewick(nwFile):
	phyloFile=open(nwFile)
	line=phyloFile.readline()
	while line=="\n":
		line=phyloFile.readline()
	if line[0]=="#":
		while line!="begin trees;\n":
			line=phyloFile.readline()
		line=phyloFile.readline()
		nwString=line.split()[-1]
	else:
		nwString=line.replace("\n","")
	#print(nwFile)
	index=0
	node=Tree()
	name=""
	distStr=""
	mutStr=""
	stateStr=""
	supportStr=""
	finished=False
	while index<len(nwString):
		if nwString[index]=="(":
			newNode=Tree()
			node.add_child(newNode)
			newNode.up=node
			node=newNode
			index+=1
		elif nwString[index]==";":
			root=node
			if mutStr!="":
				node.mutations=mutStr.split(",")
			else:
				node.mutations=[]
			if supportStr!="":
				#print(";")
				#print(supportStr)
				node.support=float(supportStr)
			if stateStr!="":
				node.rootState=stateStr.split(",")
			else:
				node.rootState=[]
			finished=True
			break
		elif nwString[index]=="[":
			index+=1
			while nwString[index]!="]":
				currentFeature=""
				if nwString[index]==",":
					index+=1
				while index<len(nwString) and nwString[index]!="]" and nwString[index]!="=":
					currentFeature+=nwString[index]
					index+=1
				if index>=len(nwString):
					print(currentFeature)
					print("issue within features")
					exit()
				if currentFeature=="&mutations" or currentFeature=="mutations" or currentFeature=="mutationsInf" or currentFeature=="&mutationsInf":
					while nwString[index]!="{":
						index+=1
					index+=1
					while nwString[index]!="}":
						mutStr+=nwString[index]
						index+=1
					index+=1
				elif currentFeature=="rootState" or currentFeature=="&rootState":
					while nwString[index]!="{":
						index+=1
					index+=1
					while nwString[index]!="}":
						stateStr+=nwString[index]
						index+=1
					index+=1
				elif currentFeature=="support" or currentFeature=="&support" or currentFeature=="IQsupport" or currentFeature=="&IQsupport":
					while nwString[index]!="=":
						index+=1
					index+=1
					while nwString[index]!="," and nwString[index]!="]":
						supportStr+=nwString[index]
						index+=1
					#index+=1
				elif currentFeature!="":
					index+=1
					if nwString[index]=="{":
						while nwString[index]!="}":
							index+=1
						index+=1
					else:
						while nwString[index]!="," and nwString[index]!="]":
							index+=1
				
			index+=1
		elif nwString[index]==":":
			index+=1
			while nwString[index]!="," and nwString[index]!=")" and nwString[index]!=";":
				distStr+=nwString[index]
				index+=1
		elif nwString[index]==",":
			if name!="":
				node.name=name
				name=""
			if distStr!="":
				node.dist=float(distStr)
				distStr=""
			if mutStr!="":
				node.mutations=mutStr.split(",")
			else:
				node.mutations=[]
			mutStr=""
			if supportStr!="":
				#print(",")
				#print(supportStr)
				node.support=float(supportStr)
			supportStr=""
			newNode=Tree()
			node=node.up
			node.add_child(newNode)
			newNode.up=node
			node=newNode
			index+=1
		elif nwString[index]==")":
			if name!="":
				node.name=name
				name=""
			if distStr!="":
				node.dist=float(distStr)
				distStr=""
			if mutStr!="":
				node.mutations=mutStr.split(",")
			else:
				node.mutations=[]
			mutStr=""
			if supportStr!="":
				#print(")")
				#print(supportStr)
				node.support=float(supportStr)
			supportStr=""
			node=node.up
			index+=1
		else:
			name+=nwString[index]
			index+=1
	if not finished:
		print("Error, final character ; not found in newick string in file "+nwFile+".")
		exit()

	phyloFile.close()
	return root


#function that changes multifurcating tree structure into a binary tree by adding 0-length branches/nodes
def makeTreeBinary(root):
	nodesToVisit=[root]
	while nodesToVisit:
		node=nodesToVisit.pop()
		if node.children:
			while len(node.children)>2:
				child2=node.children.pop()
				child1=node.children.pop()
				newParent=Tree(dist=False)
				newParent.add_child(child1)
				newParent.add_child(child2)
				child1.up=newParent
				child2.up=newParent
				newParent.up=node
				node.children.append(newParent)
			nodesToVisit.append(node.children[0])
			nodesToVisit.append(node.children[1])



#create newick string for a tree 
def createNewick(node,writeMutations=False):
	nextNode=node
	stringList=[]
	direction=0
	lastNode=None
	while nextNode!=None:
		if nextNode.children:
			if direction==0:
				#print("Will go into first child")
				stringList.append("(")
				lastNode=nextNode
				nextNode=nextNode.children[0]
			elif direction==1 and lastNode!=nextNode.children[-1]:
				#print("Will go into another child")
				stringList.append(",")
				childNum=0
				#print("Num children: "+str(len(nextNode.children)))
				while nextNode.children[childNum]!=lastNode:
					#print("Not child: "+str(childNum))
					childNum+=1
				lastNode=nextNode
				nextNode=nextNode.children[childNum+1]
				direction=0
			else:
				#print("Coming from last child")
				stringList.append(")")
				if writeMutations:
					stringList.append("[&mutations={")
					for iNode in range(len(nextNode.mutations)):
						stringList.append(nextNode.mutations[iNode])
						if iNode<(len(nextNode.mutations)-1):
							stringList.append(",")
					stringList.append("}]")

				if nextNode.dist:
					stringList.append(":"+str(nextNode.dist))
				else:
					stringList.append(":"+str(0.0))
				direction=1
				lastNode=nextNode
				nextNode=nextNode.up
		else:
			#print("Terminal node "+nextNode.name)
			stringList.append(nextNode.name)
			if writeMutations:
				stringList.append("[&mutations={")
				for iNode in range(len(nextNode.mutations)):
					stringList.append(nextNode.mutations[iNode])
					if iNode<(len(nextNode.mutations)-1):
						stringList.append(",")
				stringList.append("}]")
			if nextNode.dist:
				stringList.append(":"+str(nextNode.dist))
			else:
				stringList.append(":"+str(0.0))
			direction=1
			lastNode=nextNode
			nextNode=nextNode.up
	stringList.append(";")
	return "".join(stringList)

#create postorder traversal list of nodes for an input tree
def postorderList(phylo):
	numNodes=0
	if not phylo.children:
		return [phylo]
	nodeList=[]
	lastNode=phylo
	node=phylo.children[0]
	while node!=phylo or lastNode!=phylo.children[-1]:
		if lastNode==node.up:
			numNodes+=1
			if node.children:
				lastNode=node
				node=node.children[0]
			else:
				nodeList.append(node)
				lastNode=node
				node=node.up
		else:
			childIndex=0
			while lastNode!=node.children[childIndex]:
				childIndex+=1
			if childIndex==len(node.children)-1:
				lastNode=node
				nodeList.append(node)
				node=node.up
			else:
				lastNode=node
				node=node.children[childIndex+1]
	nodeList.append(phylo)
	print("Numbers of nodes:")
	print(numNodes+1)
	print(len(nodeList))
	return nodeList

def readConciseAlignment(fileName,numbersFirst=True,shift01pos=False):
	start = time.time()
	fileI=open(fileName)
	line=fileI.readline()
	nSeqs=0
	data={}
	addendum=0
	if shift01pos:
		addendum=1
	while line!="" and line!="\n":
		nSeqs+=1
		seqList=[]
		name=line.replace(">","").replace("\n","")
		line=fileI.readline()
		while line!="" and line!="\n" and line[0]!=">":
			linelist=line.split()
			if len(linelist)>2:
				print("Format problem")
				exit()
				#entry=(linelist[0],int(linelist[1]),int(linelist[2]))
			else:
				if numbersFirst:
					entry=(linelist[1],int(linelist[0])+addendum)
				else:
					entry=(linelist[0],int(linelist[1])+addendum)
			seqList.append(entry)
			line=fileI.readline()
		data[name]=seqList
	fileI.close()
	time2 = time.time() - start
	print("Time to read DNA reduced data file: "+str(time2))
	print(str(nSeqs)+" sequences in file.")
	return data

def readConciseAlignmentReal(fileName):
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
			if len(linelist)==3 and int(linelist[2])>1:
				entry=("N",int(linelist[1]),int(linelist[2]))
			else:
				if linelist[0]=="-":
					entry=("N",int(linelist[1]))
				else:
					entry=(linelist[0].upper(),int(linelist[1]))
			seqList.append(entry)
			line=fileI.readline()
		data[name]=seqList
	fileI.close()
	time2 = time.time() - start
	print("Time to read DNA reduced data file: "+str(time2))
	print(str(nSeqs)+" sequences in file.")
	return data


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
		if line=="":
			break
		if line[0]==">":
			break
	#lRef=len(ref)
	#print("Ref genome length: "+str(lRef))
	file.close()
	return ref
	

	
def mergeMutationsIntoList(newMutations,existingList,posStart=0,keepReversions=False):
	index1=0
	index2=0
	while index1<len(newMutations):
		mutation1=newMutations[index1]
		#pos1=int(mutation1[1:-1])
		pos1=int(mutation1[posStart:-1])
		if index2<len(existingList):
			mutation2=existingList[index2]
			#pos2=int(mutation2[1:-1])
			pos2=int(mutation2[posStart:-1])
			
			if pos1>pos2:
				index2+=1
			elif pos2>pos1:
				if mutation1[-1]!=ref[pos1-1]:
					#existingList.insert(index2,ref[pos1-1]+mutation1[1:])
					existingList.insert(index2,mutation1)
				elif keepReversions:
					existingList.insert(index2,mutation1)
				index1+=1
				index2=0
			else:
				if mutation1[-1]!=ref[pos1-1]:
					#existingList[index2]=ref[pos1-1]+mutation1[1:]
					existingList[index2]=mutation1
				elif keepReversions and mutation1[0]!=ref[pos1-1]:
					existingList[index2]=mutation1
				else:
					existingList.pop(index2)
				index1+=1
				index2=0
		else:
			if mutation1[-1]!=ref[pos1-1]:
				existingList.append(mutation1)
			elif keepReversions:
				existingList.append(mutation1)
			index1+=1
			index2=0
		

#creating folders for output files
createFolders=False
if createFolders:
	scenario=args.scenario
	folders=["simulationsBranchSupportSubsamples_noNs","simulationsBranchSupportSubsamples"]
	folder=folders[scenario]
	folderNameSimu=args.pathToSimulationFolder+folder+'/'
	if not os.path.isdir(folderNameSimu):
		os.mkdir(folderNameSimu)
	for j in [1000, 2000, 5000, 10000, 20000, 50000, 100000, 200000]:
		folderNameSimu=args.pathToSimulationFolder+folder+'/'+str(j)+"subsamples/"
		if not os.path.isdir(folderNameSimu):
			os.mkdir(folderNameSimu)
		for i in range(20):
			folderNameSimu=args.pathToSimulationFolder+folder+'/'+str(j)+"subsamples/output_repl"+str(i+1)+"/"
			if not os.path.isdir(folderNameSimu):
				os.mkdir(folderNameSimu)
	exit()


#create subsample files 
if args.createInputData:

	#creating folders for output files
	#folders=["simulationsBranchSupportSubsamples"]
	scenario=args.scenario
	folders=["simulationsBranchSupportSubsamples_noNs","simulationsBranchSupportSubsamples"]
	folder=folders[scenario]


	subsampleTreeInference=args.subSampleNum
	alignments=[args.pathToErrorSimulationFolder+"phastSim_genomes_alpha02.txt",args.pathToErrorSimulationFolder+"phastSim_genomes_alpha02_Ns.txt"]
	phylos=[args.pathToErrorSimulationFolder+"phastSim_genomes_alpha02.tree",args.pathToErrorSimulationFolder+"phastSim_genomes_alpha02.tree"]
	seed=args.repeat
	repeat=seed

	ref = collectReference(args.inputSimulationReference)
	pathToRepeat=args.pathToSimulationFolder+folder+'/'+str(subsampleTreeInference)+"subsamples/repeat"+str(repeat)+"_"+str(subsampleTreeInference)+"samples_"+folders[scenario]
	if scenario==0:
		data=readConciseAlignment(alignments[scenario],shift01pos=True)
	else:
		data=readConciseAlignmentReal(alignments[scenario])
	samples=data.keys()
	print(str(len(samples))+" sequences in the total concise DNA data file")

	subsample=subsampleTreeInference
	random.seed(seed)
	newSamples=random.sample(samples,subsampleTreeInference)

	#create subsampled MAPLE file
	diffFile=pathToRepeat+".txt"
	fileO=open(diffFile,"w")
	for s in newSamples:
		fileO.write(">"+s+"\n")
		for m in data[s]:
			if m[0]!="R":
				if len(m)==2:
					fileO.write(m[0]+"\t"+str(m[1])+"\n")
				else:
					fileO.write(m[0]+"\t"+str(m[1])+"\t"+str(m[2])+"\n")
		data[s]=None
	fileO.close()
	print("MAPLE file created")
	del data
	del samples
	
	#create fasta and phylip files for subsample
	data=readConciseAlignmentReal(diffFile)
	newSamples=data.keys()
	phylipFile=pathToRepeat+".phy"
	fastaFile=pathToRepeat+".fa"
	fileFa=open(fastaFile,"w")
	fileO=open(phylipFile,"w")
	lRef=len(ref)
	fileO.write(str(subsample)+"\t"+str(lRef)+"\n")
	for s in newSamples:
		fileO.write(s+" ")
		fileFa.write(">"+s+"\n")
		refList=list(ref)
		for m in data[s]:
			if m[0]!="R":
				if len(m)==2:
					refList[m[1]-1]=m[0]
				else:
					for i in range(m[2]):
						refList[m[1]+i-1]=m[0]
		data[s]=None
		seq="".join(refList)
		fileO.write(seq+"\n")
		fileFa.write(seq+"\n")
	fileO.close()
	fileFa.close()
	print("phylip file created")

	#create VCF file for USheER
	newFastaFile=pathToRepeat+"_withRef.fa"
	vcfFile=pathToRepeat+".vcf"
	os.system("cat "+args.inputSimulationReference+" "+fastaFile+" > "+newFastaFile+"\n")
	os.system("faToVcf "+newFastaFile+" "+vcfFile+"\n")
	os.system("rm "+newFastaFile+"\n")
	os.system("rm "+fastaFile+"\n")
	print("VCF file created")


	#extract subtree of the larger simulated tree
	#create search tree for current sample names:
	print("Total number of leaves to be found: "+str(len(newSamples)))
	binTree=BinarySearchTree()
	for s in newSamples:
		binTree.put(s)

	phylo = readNewick(phylos[scenario])
	print("Phylogeny read")
	nodeList=postorderList(phylo)

	#extract subtree for subsample
	numNode=0
	for node in nodeList:
		numNode+=1
		#print("numNode "+str(numNode))
		if len(node.children)==0:
			if binTree[node.name]:
				node.subtree=Tree(name=node.name,dist=float(node.dist)/len(ref))
				node.subtree.mutations=list(node.mutations)
			else:
				node.subtree=None
		else:
			numDesc=0
			for c in node.children:
				if c.subtree!=None:
					numDesc+=1
					child=c
			if numDesc==0:
				node.subtree=None
			elif numDesc==1:
				node.subtree=child.subtree
				node.subtree.dist+=float(node.dist)/len(ref)
				mergeMutationsIntoList(node.mutations,node.subtree.mutations,posStart=1,keepReversions=True)
				#node.subtree.mutations.extend(node.mutations)
			else:
				node.subtree=Tree(dist=float(node.dist)/len(ref))
				node.subtree.mutations=list(node.mutations)
				for c in node.children:
					if c.subtree:
						node.subtree.add_child(c.subtree)
						c.subtree.up=node.subtree
	newickString=createNewick(phylo.subtree,writeMutations=True)
	phylo=None
	del data
	phyloFile=open(pathToRepeat+"_MAT.nw","w")
	phyloFile.write(newickString+"\n")
	phyloFile.close()
	print("Subtree extracted and written to file "+pathToRepeat+"_MAT.nw")







if args.createBashScript:
	if args.numSamples==0:
		numSamples=[1000, 2000, 5000, 10000, 20000, 50000, 100000, 200000]
	else:
		numSamples=[args.numSamples]

	MAPLEversions=["0.3.2"]
	version=MAPLEversions[0]
	versionSPRTA="0.6.8"
	scenario=args.scenario
	folders=["simulationsBranchSupportSubsamples_noNs","simulationsBranchSupportSubsamples"]
	folder=folders[scenario]

	#creating bash script to run this python script in parallel on the cluster to create subsampled data
	file=open(args.pathToSimulationFolder+"createSubsampleInputFiles.sh","w")
	for j in numSamples:
		file.write("for i in $(seq 1 20)\n"+"do \n\t")
		#for scenario in range(len(folders)):
		folderNameSimu=args.pathToSimulationFolder+folders[scenario]+'/'+str(j)+"subsamples/"
		file.write("bsub -M "+str(int(15000+j/10))+" -o "+folderNameSimu+"repl\"$i\"_fileCreation_console_output.txt -e "+folderNameSimu+"repl\"$i\"_fileCreation_console_error.txt"
		+" pypy3 benchmarking_branchSupport.py --createInputData --scenario "+str(scenario)+" --subSampleNum "+str(j)+" --repeat \"$i\" \n\n\t")
		file.write("done\n\n")
	file.close()
	print("Created bash script "+args.pathToSimulationFolder+"createSubsampleInputFiles.sh")

	#creating bash scripts for running all methods
	methods=["matOptimize","MAPLE","IQTREE2","SPRTA"]
	#IQtree options:
	# UFBoot, without --fast and --tree-fix, ML tree with weird support written to ".treefile"
	# -B, --ufboot NUM     Replicates for ultrafast bootstrap (>=1000)

	# Felsenstein's bootsrap, without --tree-fix, ML tree with weird support written to ".treefile"
	# -b, --boot NUM       Replicates for bootstrap + ML tree + consensus tree
	# TBE, without --tree-fix, ML tree with support written to ".tbe.tree"
	# --tbe                Transfer bootstrap expectation

	# works with all options, including --tree-fix, ML tree with support written to ".treefile"
	# --alrt NUM           Replicates for SH approximate likelihood ratio test (Guindon et al 2010) New Algorithms and Methods to Estimate Maximum-Likelihood Phylogenies: Assessing the Performance of PhyML 3.0
	# works with all options, including --tree-fix, ML tree with support written to ".treefile"
	# --alrt 0             Parametric aLRT test (Anisimova and Gascuel 2006)
	# works with all options, including --tree-fix, ML tree with support written to ".treefile"
	# --abayes             approximate Bayes test (Anisimova et al. 2011)
	# works with all options, including --tree-fix, ML tree with support written to ".treefile"
	# --lbp NUM            Replicates for fast local bootstrap probabilities (Adachi and Hasegawa, 1996b) https://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.476.8552
	methodOptions=[[""],[""],[" -B 1000"," --fast -b 100"," --fast -b 100 --tbe"," --fast --alrt 1000"," --fast --alrt 0"," --fast --abayes"," --fast --lbp 1000"],[""]]
	methodOptionsNames=[[""],[""],["_UFBoot","_bootstrap","_tbe","_aLRT-SH","_aLRT","_aBayes","_LBP"],[""]]
	maxNumSamples=[[200000],[200000],[5000,5000,5000,20000,20000,20000,20000],[200000]]
	fileMatConv=open(args.pathToSimulationFolder+"submitMatOptimizeConversion.sh","w")
	fileMatEstimation=open(args.pathToSimulationFolder+"submitMatEstimation.sh","w")

	for method in range(len(methods)):
		methodName=methods[method]
		file=open(args.pathToSimulationFolder+"submit"+methodName+".sh","w")
		if method==0:
			file.write("export PATH=\"/miniconda3/bin:$PATH\" ; source /miniconda3/etc/profile.d/conda.sh ; conda activate usher-env \n")
			fileMatConv.write("bgadd /matConv"+str(numSamples[0])+"\n")
			fileMatConv.write("export PATH=\"/miniconda3/bin:$PATH\" ; source /miniconda3/etc/profile.d/conda.sh ; conda activate usher-env \n")
		for option in range(len(methodOptions[method])):
			for j in numSamples:
				if j<=maxNumSamples[method][option]:
					#file.write("bgadd /"+methodName+str(j)+"\n")
					file.write("for i in $(seq 1 20)\n do\n\n")
					if method==0:
						fileMatConv.write("for i in $(seq 1 20)\n do\n\n")
					if method==2 or method==0:
						fileMatEstimation.write("for i in $(seq 1 20)\n do\n\n")
					#for scenario in range(len(folders)):
					refFile=args.inputSimulationReference
					folderNameSimu=args.pathToSimulationFolder+folder+'/'+str(j)+"subsamples/output_repl\"$i\"/"
					pathToFile=args.pathToSimulationFolder+folder+'/'+str(j)+"subsamples/repeat\"$i\"_"+str(j)+"samples_"+folder

					#remove existing result if already exists - this helps in case of re-running an analysis, and in case the new version fails, to wrongly use the old tree as estimate of the new version.
					fileNamePre=folderNameSimu+methodName+methodOptionsNames[method][option]
					rateVarNames=[""]
					outputNames=["_console_output.txt","_console_error.txt"]
					for rateVarName in rateVarNames:
						for outputName in outputNames:
							file.write("rm -f "+fileNamePre+rateVarName+outputName+" || true\n")
					if method==0:
						treeFiles=[folderNameSimu+"matOptimize_output.nh"]
					elif method==1:
						treeFiles=[folderNameSimu+"MAPLE"+version+"_tree.tree"]
					elif method==2:
						treeFiles=[fileNamePre+".treefile"]

					for treeFile in treeFiles:
						file.write("rm -f "+treeFile+" || true\n")
								
					#run methods

					#matOptimize
					if method==0:
						file.write("bsub -g /"+methodName+str(j)+" -M "+str(int(400+j/5))+" -o "+fileNamePre+"_console_output.txt -e "+fileNamePre+"_console_error.txt "
						+"matOptimize -n -v "+pathToFile+".vcf -o "+folderNameSimu+"matOptimize_output.pb -T 1 -t "+folderNameSimu+"MAPLE"+version+"_tree.tree -E "+folderNameSimu+"matOptimize_support.tree \n")
						fileMatConv.write("bsub -g /matConv"+str(j)+" -M "+str(int(400+j/5))+" -o "+fileNamePre+"_conv_console_output.txt -e "+fileNamePre+"_conv_console_error.txt "
						+"matUtils extract -i "+folderNameSimu+"matOptimize_output.pb -T 1 -d "+folderNameSimu+" -t matOptimize_output.nh \n")
						fileMatEstimation.write("bsub -g /"+"MAT"+str(j)+" -M "+str(int(400+j/20))+" -o "+fileNamePre+"_MAT_console_output.txt -e "+fileNamePre+"_MAT_console_error.txt "
						+"pypy3 MAPLEv"+version+".py --reference "+refFile+" --input "
						+pathToFile+".txt --overwrite --calculateLKfinalTree --model UNREST --rateVariation --estimateMAT --keepInputUsherSupports --output "+folderNameSimu+methodName+methodOptionsNames[method][option]+"_MAT --noFastTopologyInitialSearch --numTopologyImprovements 0 --inputTree "+folderNameSimu+"matOptimize_support.tree \n")
					
					#MAPLE
					elif method==1:
						
						file.write("sbatch --parsable -J MAPLE"+str(j)+" -t 100:00:00 --mem=10G -o "+fileNamePre+"_console_output.txt -e "+fileNamePre+"_console_error.txt "
						+" --wrap=\"pypy3.10 MAPLEv"+version+".py --reference "+refFile+" --input "
						+pathToFile+".txt --overwrite --calculateLKfinalTree --model UNREST --rateVariation --networkOutput --estimateMAT --aBayesPlus --output "+folderNameSimu+"MAPLE"+version+"\" > "+fileNamePre+"_console_report.txt \n")

					#IQtree
					elif method==2:
						inputTreeExtension=".treefile"
						if option==2:
							inputTreeExtension=".tbe.tree"
						file.write("cd "+folderNameSimu+"\n\t"+"bsub -g /"+methodName+str(j)+" -M "+str(int(j*3.0))+" -o "+fileNamePre+"_console_output.txt -e "+fileNamePre+"_console_error.txt "
						+"iqtree2 -s "+pathToFile+".phy --seqtype DNA --prefix "+methodName+methodOptionsNames[method][option]+" -t "+folderNameSimu+"MAPLE"+version+"_tree.tree --seed 1 -m GTR+F+G4 --quiet --redo -nt 1 "+methodOptions[method][option]+" \n")
						
						fileMatEstimation.write("sbatch --parsable -J MAT"+str(j)+" -t 10:00:00 --mem=10G -o "+fileNamePre+"_MAT_console_output.txt -e "+fileNamePre+"MAT_console_error.txt "
						+" --wrap=\"pypy3.10 MAPLEv"+version+".py --reference "+refFile+" --input "
						+pathToFile+".txt --overwrite --calculateLKfinalTree --model UNREST --rateVariation --estimateMAT --keepInputIQtreeSupports --output "+folderNameSimu+methodName+methodOptionsNames[method][option]+"_MAT --noFastTopologyInitialSearch --numTopologyImprovements 0 --inputTree "+folderNameSimu+methodName+methodOptionsNames[method][option]+inputTreeExtension+" \" > "+fileNamePre+"_MAT_console_report.txt \n")

					#SPRTA
					elif method==3:

						file.write("sbatch --parsable -J SPRTA"+str(j)+" -t 10:00:00 --mem=10G -o "+fileNamePre+"_console_output.txt -e "+fileNamePre+"_console_error.txt "
						+" --wrap=\"pypy3.10 MAPLEv"+versionSPRTA+".py --reference "+refFile+" --input "
						+pathToFile+".txt --overwrite --model UNREST --rateVariation --networkOutput --estimateMAT --aBayesPlus --numTopologyImprovements 0 --doNotImproveTopology --inputTree "+folderNameSimu+"MAPLE"+version+"_tree.tree --output "+folderNameSimu+"SPRTA"+versionSPRTA+"\" > "+fileNamePre+"_console_report.txt \n")

					file.write("done\n\n")
					if method==0:
						fileMatConv.write("done\n\n")
					if method==2 or method==0:
						fileMatEstimation.write("done\n\n")
		file.close()
		print("Created "+methodName+" bash script "+args.pathToSimulationFolder+"submit"+methodName+".sh")
	fileMatConv.close()
	fileMatEstimation.close()
	print("Created matConvert bash script "+args.pathToSimulationFolder+"submitMatOptimizeConversion.sh")
	print("Created fileMatEstimation bash script "+args.pathToSimulationFolder+"submitMatEstimation.sh")






def highProbMutations(muts,threshold=0.9):
	highProbMuts=[]
	for mut in muts:
		mutList=mut.split(":")
		if len(mutList)==1:
			#highProbMuts.append(mut)
			highProbMuts.append(mut[1:])
		else:
			supp=float(mutList[1])
			if supp>threshold:
				#highProbMuts.append(mutList[0])
				highProbMuts.append(mutList[0][1:])
	return highProbMuts

# when .support is available for the branch/node, store the support of each mutation into the dictionary
def traversTreeMutations(node,currentMuts,mutDict,alsoReverse=True,supportDict={},registerSupport=False,compareTo={}):
	nodesToVisit=[(node,currentMuts)]
	while nodesToVisit:
		node,muts=nodesToVisit.pop()
		for child in node.children:
			#TODO add mutations and their reverse to dictionary
			key=",".join(muts)
			childMuts=highProbMutations(child.mutations)
			if	compareTo: 
				if not (key in compareTo):
					print("Parent genome "+key+" not found for node ")
					if hasattr(child,"name"):
						print(child.name)
					else:
						print(child)
					if hasattr(child,"support"):
						print("Support: "+str(child.support))
					print("with mutations")
					print(childMuts)
					print("")
				else:
					notFoundMu=[]
					for mu in childMuts:
						if not (mu in compareTo[key]):
							notFoundMu.append(mu)
					if notFoundMu:
						print("Some mutations not found in parent genome "+key+" for child node ")
						if hasattr(child,"name"):
							print(child.name)
						else:
							print(child)
						if hasattr(child,"support"):
							print("Support: "+str(child.support))
						print(notFoundMu)
						print("Out of total")
						print(childMuts)
						print("")
			if key in mutDict:
				mutDict[key].extend(childMuts)
			else:
				mutDict[key]=childMuts
			if registerSupport: 
				if hasattr(child, 'support') or hasattr(child, 'IQsupport'):
					supp=child.support
					if not (key in supportDict):
						supportDict[key]={}
					for mut in childMuts:
						supportDict[key][mut]=supp
			newMuts=list(muts)
			mergeMutationsIntoList(childMuts,newMuts)
			nodesToVisit.append((child,newMuts))

def ancestralGenome(muts,threshold=0.9):
	key=[]
	for mut in muts:
		if mut[1]<threshold:
			if mut[1]>1.0-threshold:
				return False
		else:
			key.append(mut[0])
	if key:
		return ",".join(key)
	else:
		return ""
	

#account for mutation probabilities, remove mutations with probability below threshold
def mergeUncertainMutationsIntoList(newMutations,existingList,posStart=0,threshold=0.9):
	index1=0
	index2=0
	while index1<len(newMutations):
		mut=newMutations[index1]
		mutList=mut.split(":")
		supp=float(mutList[1])
		mutation1=mutList[0]
		pos1=int(mutation1[1:-1])
		if index2<len(existingList):
			mut2=existingList[index2]
			mutation2=mut2[0]
			supp2=mut2[1]
			pos2=int(mutation2[posStart:-1])
			
			if pos1>pos2:
				if supp2<=1.0-threshold:
					existingList.pop(index2)
				else:
					index2+=1
			elif pos2>pos1:
				if mutation1[-1].lower()!=ref[pos1-1].lower() and supp>(1.0-threshold):
					existingList.insert(index2,(mutation1[1:],supp))
				index1+=1
				index2=0
			else:
				if mutation1[-1].lower()==mutation2[-1].lower():
					existingList[index2]=(mutation2,supp+supp2)
				elif mutation1[-1].lower()!=ref[pos1-1].lower():
					if mutation1[0].lower()==mutation2[-1].lower():
						if supp>supp2-supp:
							if supp>1.0-threshold:
								existingList[index2]=(mutation1[1:],supp)
							else:
								existingList.pop(index2)
						else:
							if supp2-supp>1.0-threshold:
								existingList[index2]=(mutation2,supp2-supp)
							else:
								existingList.pop(index2)
					else:
						if supp>supp2:
							if supp>1.0-threshold:
								existingList[index2]=(mutation1[1:],supp)
							else:
								existingList.pop(index2)
						else:
							if supp2>1.0-threshold:
								existingList[index2]=(mutation2,supp2-supp)
							else:
								existingList.pop(index2)

				else:
					if supp2-supp>1.0-threshold:
						existingList[index2]=(mutation2,supp2-supp)
					else:
						existingList.pop(index2)
				index1+=1
				index2=0
		else:
			if mutation1[-1].lower()!=ref[pos1-1].lower() and supp>(1.0-threshold):
				existingList.append((mutation1[1:],supp))
			index1+=1
			index2=0



# only consider correctness for nodes whose parent node genome is not uncertain
def traversConfidentTreeMutations(node,currentMuts,mutDict,alsoReverse=True,supportDict={},registerSupport=False,compareTo={}):
	nodesToVisit=[(node,currentMuts)]
	while nodesToVisit:
		node,muts=nodesToVisit.pop()
		for child in node.children:
			childMuts=highProbMutations(child.mutations)
			key=ancestralGenome(muts)
			if	compareTo: 
				if not (key in compareTo):
					print("Parent genome "+key+" not found for node ")
					if hasattr(child,"name"):
						print(child.name)
					else:
						print(child)
					if hasattr(child,"support"):
						print("Support: "+str(child.support))
					print("with mutations")
					print(childMuts)
					print("")
				else:
					notFoundMu=[]
					for mu in childMuts:
						if not (mu in compareTo[key]):
							notFoundMu.append(mu)
					if notFoundMu:
						print("Some mutations not found in parent genome "+key+" for child node ")
						if hasattr(child,"name"):
							print(child.name)
						else:
							print(child)
						if hasattr(child,"support"):
							print("Support: "+str(child.support))
						print(notFoundMu)
						print("Out of total")
						print(childMuts)
						print("")
			if key!=False:
				if key in mutDict:
					mutDict[key].extend(childMuts)
				else:
					mutDict[key]=childMuts
				if registerSupport: 
					if hasattr(child, 'support') or hasattr(child, 'IQsupport'):
						supp=child.support
						if not (key in supportDict):
							supportDict[key]={}
						for mut in childMuts:
							supportDict[key][mut]=supp

			newMuts=list(muts)
			mergeUncertainMutationsIntoList(child.mutations,newMuts)
			nodesToVisit.append((child,newMuts))


# when .support is available for the branch/node, store the support of each mutation into the dictionary
def createMutDict(rootNode,isTrueTree=True,alsoReverse=True,registerSupport=False,supportDict={},compareTo={},threshold=0.9):
	mutDict={}
	if isTrueTree:
		muts=[]
		mergeMutationsIntoList(highProbMutations(rootNode.mutations),muts)
		traversTreeMutations(rootNode,muts,mutDict,alsoReverse=alsoReverse)
	else:
		muts=[]
		for state in rootNode.rootState:
			stateList=state.split(":")
			supp=float(stateList[1])
			if supp>1.0-threshold:
				state=stateList[0]
				pos=int(state[1:])
				if ref[pos-1].lower()!=state[0].lower():
					muts.append((str(pos)+state[0],supp))
		if not alsoReverse:
			traversConfidentTreeMutations(rootNode,muts,mutDict,alsoReverse=alsoReverse,supportDict=supportDict,registerSupport=registerSupport,compareTo=compareTo)
		else:
			print("reversing mutations not supported anymore")
			exit()
	return mutDict


# measure accuracy and computational demand, and write measures to file.
if args.compareMats:
	scenario=args.scenario
	folders=["simulationsBranchSupportSubsamples_noNs","simulationsBranchSupportSubsamples"]
	folder=folders[scenario]
	ref = collectReference(args.inputSimulationReference)
	if args.numSamples==0:
		numSamples=[1000, 2000, 5000, 10000, 20000, 50000, 100000, 200000]
	else:
		numSamples=[args.numSamples]

	MAPLEversions=["0.3.2"]
	version=MAPLEversions[0]
	versionSPRTA="0.6.8"
	#creating bash scripts for running all methods
	methods=["matOptimize","MAPLE","IQTREE2","SPRTA"]
	methodOptions=[[""],[""],[" -B 1000"," --fast -b 100"," --fast -b 100 --tbe"," --fast --alrt 1000"," --fast --alrt 0"," --fast --abayes"," --fast --lbp 1000"],[""]]
	methodOptionsNames=[[""],[""],["_UFBoot","_bootstrap","_tbe","_aLRT-SH","_aLRT","_aBayes","_LBP"],[""]]
	maxNumSamples=[[200000],[200000],[5000,2000,2000,20000,20000,20000,20000],[200000]]
	
	for j in numSamples:
		refFile=args.inputSimulationReference

		supportsCorrect=[]
		supportsWrong=[]
		correctMutations=[]
		wrongMutations=[]
		times=[]
		memories=[]
		maxMems=[]
		for method in range(len(methods)):
			supportsCorrect.append([])
			supportsWrong.append([])
			correctMutations.append([])
			wrongMutations.append([])
			times.append([])
			memories.append([])
			maxMems.append([])
			for option in range(len(methodOptions[method])):
				#if j<=maxNumSamples[method][option]:
				supportsCorrect[method].append([])
				supportsWrong[method].append([])
				correctMutations[method].append(0)
				wrongMutations[method].append(0)
				times[method].append([])
				memories[method].append([])
				maxMems[method].append([])
				for nRepeat in range(1,21):
					supportsCorrect[method][option].append([])
					supportsWrong[method][option].append([])

		for nRepeat in range(1,21):
			#print("Repeat "+str(nRepeat))
			folderNameSimu=args.pathToSimulationFolder+folder+'/'+str(j)+"subsamples/output_repl"+str(nRepeat)+"/"
			pathToFile=args.pathToSimulationFolder+folder+'/'+str(j)+"subsamples/repeat"+str(nRepeat)+"_"+str(j)+"samples_"+folder
			pathToRepeat=args.pathToSimulationFolder+folder+'/'+str(j)+"subsamples/repeat"+str(nRepeat)+"_"+str(j)+"samples_"+folders[scenario]
			trueTreeFile=pathToRepeat+"_MAT.nw"

			trueTree=readNewick(trueTreeFile)
			# create true mutation dictionary
			#trueMutDict=createMutDict(trueTree,isTrueTree=True)
			trueMutDict=createMutDict(trueTree,isTrueTree=True,alsoReverse=False)

			for method in range(len(methods)):
				for option in range(len(methodOptions[method])):
					if j<=maxNumSamples[method][option]:

						#tree file names
						fileNamePre=folderNameSimu+methods[method]+methodOptionsNames[method][option]

						if method==0:
							treeFiles=[]
							matFiles=[]
						elif method==1:
							treeFiles=[folderNameSimu+"MAPLE"+version+"_nexusTree.tree"]
							matFiles=[folderNameSimu+"MAPLE"+version+"_nexusTree.tree"]
						elif method==2:
							treeFiles=[fileNamePre+".treefile"]
							if option==2:
								treeFiles=[fileNamePre+".tbe.tree"]
							matFiles=[fileNamePre+"_MAT_nexusTree.tree"]
						elif method==3:
							treeFiles=[folderNameSimu+"SPRTA"+versionSPRTA+"_nexusTree.tree"]
							matFiles=[folderNameSimu+"SPRTA"+versionSPRTA+"_nexusTree.tree"]
						
						for matFile in matFiles:
							#print("Exploring MAT estimated by "+methods[method]+methodOptionsNames[method][option])
							if not os.path.exists(matFile):
								print("Problem with opening this MAT, maybe it doesn't exists?")
								print(matFile)
								continue
							else:
								inferredTree=readNewick(matFile)

							supportDict={}
							if args.compareMats_example and method==1 and j==1000 and nRepeat==1:
								estimatedMutDict=createMutDict(inferredTree,isTrueTree=False,alsoReverse=False,registerSupport=True,supportDict=supportDict,compareTo=trueMutDict)
								exit()
							else:
								estimatedMutDict=createMutDict(inferredTree,isTrueTree=False,alsoReverse=False,registerSupport=True,supportDict=supportDict)
							estKeys=list(estimatedMutDict.keys())
							#requirePresenceAlsoInReverse=True
							
							#foundSupports=[]
							#notFoundSupports=[]
							for key in estKeys:
								if key in trueMutDict:
									#keyFound+=1
									for mut in estimatedMutDict[key]:
										if mut in trueMutDict[key]:
											#found+=1
											correctMutations[method][option]+=1
											if (key in supportDict) and (mut in supportDict[key]):
												#foundSupports.append(supportDict[key][mut])
												supportsCorrect[method][option][nRepeat-1].append(supportDict[key][mut])
										else:
											#notFound+=1
											wrongMutations[method][option]+=1
											if (key in supportDict) and (mut in supportDict[key]):
												#notFoundSupports.append(supportDict[key][mut])
												supportsWrong[method][option][nRepeat-1].append(supportDict[key][mut])
								else:
									#keyNotFound+=1
									for mut in estimatedMutDict[key]:
										#notFound+=1
										wrongMutations[method][option]+=1
										if (key in supportDict) and (mut in supportDict[key]):
											#notFoundSupports.append(supportDict[key][mut])
											supportsWrong[method][option][nRepeat-1].append(supportDict[key][mut])

							if os.path.isfile(fileNamePre+"_console_output.txt"):

								if (method==1 or method==3) and scenario==1:
									file=open(fileNamePre+"_console_report.txt")
									line=file.readline()
									jobNumber=line.replace("\n","")
									file.close()
									os.system("jobinfo "+jobNumber+" > "+fileNamePre+"_console_report2.txt")
									file=open(fileNamePre+"_console_report2.txt")
									for indexA in range(8):
										line=file.readline()
									linelist=line.split()
									if len(linelist)>1 and linelist[2]=="COMPLETED":
										for indexA in range(7):
											line=file.readline()
										linelist=line.split()
										if len(linelist)>1:
											#print(linelist)
											nDaysList=linelist[3].split("-")
											if len(nDaysList)>1:
												nDays=int(nDaysList[0])
												timeList=nDaysList[1].split(":")
											else:
												nDays=0
												timeList=nDaysList[0].split(":")
											#timeList=linelist[3].split(":")
											timeV=int(timeList[2])+int(timeList[1])*60+int(timeList[0])*3600+int(nDays)*24*3600
											times[method][option].append(timeV)
										else:
											times[method][option].append(float("nan"))
										for indexA in range(5):
											line=file.readline()
										linelist=line.split()
										if len(linelist)>1:
											if linelist[4][-1]=="M":
												memV=float(linelist[4].replace("M",""))
											elif linelist[4][-1]=="K":
												memV=float(linelist[4].replace("K",""))
											else:
												memV=float(linelist[4].replace("G",""))*1000
											if memV<80:
												memV=80.0
											maxMems[method][option].append(memV)
										else:
											maxMems[method][option].append(float("nan"))
									else:
										times[method][option].append(float("nan"))
										maxMems[method][option].append(float("nan"))
									memories[method][option].append(float("nan"))
									file.close()
								
								else:

									file=open(fileNamePre+"_console_output.txt")
									line=file.readline()
									while line!="Resource usage summary:\n" and line!="":
										line=file.readline()
									if line!="":
										line=file.readline()
										line=file.readline()
										times[method][option].append(float(line.split()[3]))
										line=file.readline()
										if line.split()[3]=="-":
											pass
										else:
											maxMems[method][option].append(float(line.split()[3]))
											line=file.readline()
											if line.split()[3]=="-":
												pass
											else:
												memories[method][option].append(float(line.split()[3]))
									else:
										pass
										print("time-memory info not found in "+fileNamePre+"_console_output.txt")
									file.close()
							else:
								print("Problem with opening this cluster output file, it doesn't exists?")
								print(fileNamePre+"_console_output.txt")



		for method in range(len(methods)):
			for option in range(len(methodOptions[method])):
				if j<=maxNumSamples[method][option]:
					print("MATs estimated by "+methods[method]+methodOptionsNames[method][option])
					if times[method][option]:
						print("Mean time: "+str(sum(times[method][option])/len(times[method][option]))+ " from "+str(len(times[method][option]))+" values")
						pathToFileT=args.pathToSimulationFolder+folder+"/times_"+str(j)+"subsamples_"+methods[method]+methodOptionsNames[method][option]+".txt"
						file=open(pathToFileT,"w")
						for supp in times[method][option]:
							file.write(str(supp)+"\t")
						file.close()
					if memories[method][option]:
						print("Mean memory: "+str(sum(memories[method][option])/len(memories[method][option]))+ " from "+str(len(memories[method][option]))+" values")
						pathToFileT=args.pathToSimulationFolder+folder+"/memories_"+str(j)+"subsamples_"+methods[method]+methodOptionsNames[method][option]+".txt"
						file=open(pathToFileT,"w")
						for supp in memories[method][option]:
							file.write(str(supp)+"\t")
						file.close()
					if maxMems[method][option]:
						print("Mean maxMemory: "+str(sum(maxMems[method][option])/len(maxMems[method][option]))+ " from "+str(len(maxMems[method][option]))+" values")
						pathToFileT=args.pathToSimulationFolder+folder+"/maxMems_"+str(j)+"subsamples_"+methods[method]+methodOptionsNames[method][option]+".txt"
						file=open(pathToFileT,"w")
						for supp in maxMems[method][option]:
							file.write(str(supp)+"\t")
						file.close()

					print("Correctly inferred mutations: "+str(correctMutations[method][option])+"   , not found: "+str(wrongMutations[method][option]))
					allSupportsCorrect=[]
					allSupportsWrong=[]
					for nRepeat in range(20):
						allSupportsCorrect.extend(supportsCorrect[method][option][nRepeat])
						allSupportsWrong.extend(supportsWrong[method][option][nRepeat])
					
					#print("Correctly inferred ancestral genomes: "+str(keyFound)+"   , not found: "+str(keyNotFound))
					
					if allSupportsCorrect:
						print("Mean support of correctly inferred mutations: "+str(sum(allSupportsCorrect)/len(allSupportsCorrect))+ " from "+str(len(allSupportsCorrect))+" values")
						pathToResultsCorrect=args.pathToSimulationFolder+folder+"/supports_correct_"+str(j)+"subsamples_"+methods[method]+methodOptionsNames[method][option]+".txt"
						file=open(pathToResultsCorrect,"w")
						for nRepeat in range(20):
							for supp in supportsCorrect[method][option][nRepeat]:
								file.write(str(supp)+"\t")
							file.write("\n")
						file.close()
						
					if allSupportsWrong:
						print("Mean support of incorrectly inferred mutations: "+str(sum(allSupportsWrong)/len(allSupportsWrong))+ " from "+str(len(allSupportsWrong))+" values")
						pathToResultsWrong=args.pathToSimulationFolder+folder+"/supports_wrong_"+str(j)+"subsamples_"+methods[method]+methodOptionsNames[method][option]+".txt"
						file=open(pathToResultsWrong,"w")
						for nRepeat in range(20):
							#if method==1:
							#	print(sum(supportsWrong[method][option][nRepeat])/len(supportsWrong[method][option][nRepeat]))
							for supp in supportsWrong[method][option][nRepeat]:
								file.write(str(supp)+"\t")
							file.write("\n")
						file.close()






def printMinMeanMax(valueList):
	print(valueList)
	print(str(min(valueList))+" \t"+str(sum(valueList)/len(valueList))+" \t"+str(max(valueList)))








# create plot comparing support scores for mutations in the real data vs simulations
plotSuppDistribution=False
if plotSuppDistribution:
	#import numpy as np
	#import seaborn as sns
	import matplotlib.pyplot as plt
	folders=["simulationsResults_noNs/","simulationsResults/"]
	supps=[[],[]]
	for folderN in range(len(folders)):
		print("supports_correct_200000subsamples_SPRTA.txt")
		folder=folders[folderN]
		file=open(folder+"supports_correct_200000subsamples_SPRTA.txt")
		for nRepl in range(20):
			line=file.readline().split()
			for l in line:
				supps[folderN].append(float(l))
		file.close()
		print("supports_wrong_200000subsamples_SPRTA.txt")
		file=open(folder+"supports_wrong_200000subsamples_SPRTA.txt")
		for nRepl in range(20):
			line=file.readline().split()
			for l in line:
				supps[folderN].append(float(l))
		file.close()
		print(len(supps[folderN]))
		print(sum(supps[folderN])/len(supps[folderN]))
	realSupps=[]
	ref = collectReference("Viridian_2Malignment_filtered_masked_noShortDeletions.maple")
	print("Viridian_2M_noShortDel_deeper_SPRTA_nexus_tree.tree")
	inferredTree=readNewick("Viridian_2M_noShortDel_deeper_SPRTA_nexus_tree.tree")
	supportDict={}
	print("creating dictionary")
	estimatedMutDict=createMutDict(inferredTree,isTrueTree=False,alsoReverse=False,registerSupport=True,supportDict=supportDict)
	print("creating list")
	genomes=list(supportDict.keys())
	print(len(genomes))
	for genome in genomes:
		muts=list(supportDict[genome].keys())
		for mut in muts:
			realSupps.append(supportDict[genome][mut])
	print(len(realSupps))
	print(sum(realSupps)/len(realSupps))

	#fig, ax = plt.subplots()
	fig = plt.figure(figsize=(15, 9))
	ax1 = fig.add_subplot(111)
	ax1.set_axisbelow(True)
	ax1.grid(False)
	ax1.spines[['right', 'top']].set_visible(False)
	for axis in ['top','bottom','left','right']:
		ax1.spines[axis].set_linewidth(2.5)
	ax1.tick_params(width=2.2)
	plt.xticks(fontsize=28)
	plt.yticks(fontsize=28)
	ax1.set_xlabel('Support score', fontsize=34)
	ax1.set_ylabel('Cumulative distribution', fontsize=34)
	ax1.ecdf(realSupps,color='blue',label='Real data, mean=0.926',linewidth=3)
	ax1.ecdf(supps[1],color='red',label='Simulations Ns, mean=0.957',linewidth=3)
	ax1.ecdf(supps[0],color='green',label='Simulations no Ns, mean=0.974',linewidth=3)

	handles, labels = plt.gca().get_legend_handles_labels()
	plt.legend(labels=['Real data', 'Simulations no Ns', 'Simulations Ns'])
	leg = plt.legend(handles, labels,loc="upper left", prop={'size': 25})
	for i in range(len(leg.get_lines())):
		leg.get_lines()[i].set_linewidth(3)
	plt.savefig('distributions_supportScores.pdf')
	plt.close('all')

	exit()


#Actually create the figures.
#this has to be executed somewhere where there is matpolotlib installed - so not on the cluster
if args.runFigureGeneration:
	minV=-0.001
	minY=0.01
	maxY=8000

	leavesRange=[1000, 2000, 5000, 10000, 20000, 50000, 100000, 200000]
	nLeaves=leavesRange
	nLeavesText=["1000", "2000", "5000","10000","20000","50000","100000","200000"]
	nLeavesTextTicks=["1000", "2000",'5000',"10⁴","2x10⁴","5x10⁴","10⁵","2x10⁵"]
	#methods=["UShER","matOptimize","MAPLE","IQTREE2"]
	methods=["MAPLE","IQTREE2"]
	methods=["SPRTA","IQTREE2"]
	#methodOptions=[[""],[""],[""],[" -B 1000"," --fast -b 100"," --fast -b 100 --tbe"," --fast --alrt 1000"," --fast --alrt 0"," --fast --abayes"," --fast --lbp 1000"]]
	#methodOptionsNames=[[""],[""],[""],["_UFBoot","_bootstrap","_tbe","_aLRT-SH","_aLRT","_aBayes","_LBP"]]
	methodOptionsNames=[[""],["_UFBoot","_bootstrap","_tbe","_aLRT-SH","_aLRT","_aBayes","_LBP"]]
	methodOptionsNamesLegend=[["SPRTA"],["UFBoot","FB","TBE","aLRT-SH","aLRT","aBayes","LBP"]]
	methodOptionsNames=[[""],["_aLRT","_aBayes","_aLRT-SH","_LBP","_UFBoot","_tbe","_bootstrap"]]
	methodOptionsNamesLegend=[["SPRTA"],["aLRT","aBayes","aLRT-SH","LBP","UFBoot","TBE","FB"]]
	#maxNumSamples=[[200000],[200000],[200000],[5000,2000,2000,20000,20000,20000,20000]]
	maxNumSamples=[[200000],[5000,2000,2000,20000,20000,20000,20000]]
	maxNumSamples=[[200000],[20000,20000,20000,20000,5000,2000,2000]]

	valueNames=["mean_correct","mean_wrong","all_correct","all_wrong","mean_correct-wrong","sum_correct-wrong","num_correct","num_wrong","wrong-correct_num_ratio","time","memory","maxMemory"]
	valueNames=["mean_correct","mean_wrong","all_correct","all_wrong","mean_correct-wrong","sum_correct-wrong","num_correct","num_wrong","wrong-correct_num_ratio","time","maxMemory"]
	yAxisLabels=["Mean supp. correct mut.","Mean supp. wrong mut.","Supports of all correctly inferred mutations","Supports of all wrongly inferred mutations","Difference of the means of correct and wrong mutation support","Sum of supports of correct mutations minus supports of wrong mutations","Number of correctly inferred mutations with support","Number of wrongly inferred mutations with support","Proportion of wrongly inferred mutations with support","Time (seconds)","Average memory (Mb)","Maximum memory (Mb)"]
	yAxisLabels=["Mean supp. correct mut.","Mean supp. wrong mut.","Supports of all correctly inferred mutations","Supports of all wrongly inferred mutations","Diff. correct vs wrong support","Sum of supports of correct mutations minus supports of wrong mutations","Number of correctly inferred mutations with support","Number of wrongly inferred mutations with support","Proportion of wrongly inferred mutations with support","Time (seconds)","Maximum memory (Mb)"]

	#logYlist=[True,True,True,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False]
	logYlist=[False,False,False,False,False,False,True,True,False,True,True,True]
	#compareToMaxList=[False,False,False,True,True,True,True,False,False,True,False,False,False,False,False,False,False,False,False,False]

	#TODO hardcoded for now
	if args.scenario==0:
		figureFolder="/simulationsResults_noNs/"
	elif args.scenario==1:
		figureFolder="/simulationsResults/"
	summaryfolder=figureFolder

	import numpy as np
	import matplotlib.pyplot as plt
	import matplotlib.patches as mpatches
	import matplotlib.ticker as ticker
	import math

	def SuperScriptinate(number):
		return number.replace('0','⁰').replace('1','¹').replace('2','²').replace('3','³').replace('4','⁴').replace('5','⁵').replace('6','⁶').replace('7','⁷').replace('8','⁸').replace('9','⁹').replace('-','⁻')

	def sci_notation(number, sig_fig=2):
		ret_string = "{0:.{1:d}e}".format(number, sig_fig)
		a,b = ret_string.split("e")
		b = int(b)         # removed leading "+" and strips leading zeros too.
		return a + "x10^" + SuperScriptinate(str(b))

	#BOXPlot drawing
	def boxplot(valuesLists,axisLabels,plotFileName,labels,colors, xLabel,topPlot,degreeSkew=45):

		fig, ax1 = plt.subplots(figsize=(10, 6))
		fig.canvas.set_window_title('Simulation Running times')
		fig.subplots_adjust(left=0.075, right=0.95, top=0.9, bottom=0.25)
		
		data=np.zeros((len(valuesLists)*len(valuesLists[0]),len(valuesLists[0][0])))
		print((len(valuesLists)*len(valuesLists[0]),len(valuesLists[0][0])))
		positions=[]
		for i in range(len(valuesLists)):
			for j in range(len(valuesLists[i])):
				c=colors[j]
				positions.append(i+(j-len(valuesLists[i])/2.0)*0.5/len(valuesLists[i]))
				for k in range(len(valuesLists[i][j])):
					data[i*len(valuesLists[0])+j][k]=valuesLists[i][j][k]
				position=[i+(j-(len(valuesLists[i])-1)/2.0)*0.5/len(valuesLists[i])]
				plt.boxplot(data[i*len(valuesLists[0])+j], positions=position, notch=False, patch_artist=True, widths=0.5/len(colors), manage_ticks=False, 
					boxprops=dict(facecolor=c, color=c),
					capprops=dict(color=c),
					whiskerprops=dict(color=c),
					medianprops=dict(color=c),
					flierprops = dict(marker='o', markerfacecolor=c, markersize=1.5,
						linestyle='none', markeredgecolor=c)
				)
		
		# Add a horizontal grid to the plot, but make it very light in color
		# so we can use it for reading data values but not be distracting
		ax1.yaxis.grid(True, linestyle='-', which='major', color='lightgrey',   alpha=0.5)

		# Hide these grid behind plot objects
		ax1.set_axisbelow(True)
		ax1.set_title('Comparison of Simulation Running times')
		ax1.set_xlabel(xLabel)
		ax1.set_ylabel('Time (seconds)')

		# Set the axes ranges and axes labels
		ax1.set_xlim(-0.7, len(valuesLists) -0.5)
		top = topPlot
		bottom = -1
		ax1.set_ylim(bottom, top)
		ax1.set_xticks(np.arange(len(valuesLists)))
		ax1.set_xticklabels(labels, rotation=degreeSkew, fontsize=11)

		# Finally, add a basic legend
		for i in range(len(colors)):
			fig.text(0.10, 0.85-0.045*i, axisLabels[i], backgroundcolor=colors[i], color='black', weight='roman', size='medium')
		
		fig.savefig(plotFileName)
		plt.close()

	def errplot(times,labels,plotFileName,n_leaves,colors,topPlot,linestyles,title="Comparison of Tree Estimation Running Times",yAxisLabel="Time (seconds)",logY=True,ymin=None,ymax=None,violin=True,legendSize=14,tickSize=12,labelSize=18,drawLegend=True,legendLoc='upper left',valuesXticks=None,baseLocator=None,linewidth=3,linewidthViolin=False,thickAxis=True,ncolLegend=1):
		values=[]
		mean_times = []
		errors = []
		for i in range(len(times[0])):
			values.append([])
			for j in range(len(times)):
				if len(times[j][i])>1:
					values[-1].append([])
					for k in range(len(times[j][i])):
						values[-1][-1].append(times[j][i][k])
				else: #does it work?
					values[-1].append([])
		for t1 in times:
			mean_times.append([])
			errors.append([])
			for t2 in t1:
				if len(t2)>1:
					mean_times[-1].append(np.mean([float(t) for t in t2]))
					if violin:
						errors[-1].append(0.0)
					else:
						errors[-1].append(np.std([float(t) for t in t2]))
				else:
					mean_times[-1].append(float("nan"))
					errors[-1].append(0.0)

		mean_times = np.array(mean_times).T
		errors = np.array(errors).T    
		
		if violin:
			patches=[]
			for i in range(len(labels)):
				patches.append(mpatches.Patch(color=colors[i]))

		x = n_leaves
		wids=[]
		for xi in x:
			wids.append(xi/4)
		#y = range(100,200)
		fig = plt.figure(figsize=(15, 9))
		ax1 = fig.add_subplot(111)
		ax1.set_axisbelow(True)
		ax1.grid(False)
		ax1.set_title(title, fontsize=24)
		ax1.set_xlabel(topPlot, fontsize=labelSize)
		ax1.set_ylabel(yAxisLabel, fontsize=labelSize)
		if valuesXticks==None:
			plt.xticks(fontsize=tickSize)
		plt.yticks(fontsize=tickSize)
		if ymin!=None:
			ax1.set_ylim([ymin, ymax])
		ax1.set_xscale('log')
		if logY:
			ax1.set_yscale('log')
		ax1.yaxis.grid(True, linestyle='-', which='major', color='lightgrey',   alpha=0.5)
		if baseLocator!=None:
			loc = ticker.MultipleLocator(base=baseLocator) # this locator puts ticks at regular intervals
			ax1.yaxis.set_major_locator(loc)

		for i in range(len(labels)):
			#ax1.errorbar(x, mean_times[i, :], yerr=errors[i, :], marker='x', c=colors[i], label=labels[i], capsize=2)
			if violin:
				if len(values[i])>0:
					#print(values[i])
					positionsPlot=[]
					widsPlot=[]
					valuesPlot=[]
					for j in range(len(values[i])):
						if len(values[i][j])>0:
							positionsPlot.append(x[j])
							widsPlot.append(wids[j])
							valuesPlot.append(values[i][j])
					if len(positionsPlot)>0:
						if linewidthViolin:
							violin_parts = ax1.violinplot(valuesPlot, positions=positionsPlot, vert=True, widths=widsPlot, showmeans=False, showextrema=False, showmedians=False, quantiles=None, points=100, bw_method=None, data=None)
						else:
							violin_parts = ax1.violinplot(valuesPlot, positions=positionsPlot, vert=True, widths=widsPlot, showmeans=False, showextrema=False, showmedians=False, quantiles=None, points=100, bw_method=None, data=None)
					
						for pc in violin_parts['bodies']:
							pc.set_facecolor(colors[i])
							pc.set_edgecolor(colors[i])
						#ax1.errorbar(x, mean_times[i, :], yerr=errors[i, :], marker='.', c=colors[i], label=labels[i], capsize=2, linestyle=linestyles[i])
						if linestyles!=None:
							ax1.plot(x, mean_times[i, :], marker='.', markersize=13, c=colors[i], label=labels[i], linestyle=linestyles[i],linewidth=linewidth)
						else:
							ax1.plot(x, mean_times[i, :], marker='.', markersize=13, c=colors[i], label=labels[i],linewidth=linewidth)
			else:
				if linestyles!=None:
					ax1.errorbar(x, mean_times[i, :], yerr=errors[i, :], marker='x', c=colors[i], label=labels[i], capsize=2, linestyle=linestyles[i],linewidth=linewidth)
				else:
					ax1.errorbar(x, mean_times[i, :], yerr=errors[i, :], marker='x', c=colors[i], label=labels[i], capsize=2,linewidth=linewidth)
		if drawLegend:
			# reordering the labels
			handles, labels = plt.gca().get_legend_handles_labels()
			handles.reverse()
			labels.reverse()
			leg = plt.legend(handles, labels,loc=legendLoc, prop={'size': legendSize}, ncol=ncolLegend, borderpad=0.1, labelspacing=0.2, handlelength=1.5, handletextpad=0.2, borderaxespad=0.1, columnspacing=1.0)
			for i in range(len(leg.get_lines())):
				leg.get_lines()[i].set_linewidth(linewidth)
	
		if valuesXticks!=None:
			plt.xticks(x,valuesXticks[:len(x)],fontsize=tickSize)
		ax1.grid(False)
		ax1.spines[['right', 'top']].set_visible(False)

		if thickAxis:
			for axis in ['top','bottom','left','right']:
				ax1.spines[axis].set_linewidth(2.5)
			# increase tick width
			ax1.tick_params(width=2.2)

		fig.savefig(plotFileName,bbox_inches='tight')
		plt.close('all')

	def scatterplot(values,valuesErr,lRef,plotFileName,topPlot,colors=["blue","red"],annotations=[],annotations2=[],title="Estimated rates along the genome",yAxisLabel="Rate",logY=True,ymin=0.000001,ymax=10000.0,legendSize=14,tickSize=12,labelSize=18,drawLegend=True,legendLoc='upper left',valuesXticks=None,baseLocator=None,linewidth=3,labels=["Substitution rate","Error probability"]): 
		x = range(1,lRef+1)
		fig = plt.figure(figsize=(13, 7))
		ax1 = fig.add_subplot(111)
		ax1.set_axisbelow(True)
		#ax1.grid(False)
		ax1.set_title(title, fontsize=24)
		ax1.set_xlabel(topPlot, fontsize=labelSize)
		ax1.set_ylabel(yAxisLabel, fontsize=labelSize)
		ax1.set_ylim([ymin, ymax])
		if logY:
			ax1.set_yscale('log')

		ax1.plot(x,values,color=colors[0],linestyle="None",marker=".",label=labels[0])
		ax1.plot(x,valuesErr,color=colors[1],linestyle="None",marker=".",label=labels[1])
		for annotation in annotations:
			ax1.annotate(annotation[2],(annotation[0],annotation[1]), ha='left', rotation=45,color=colors[0])
		for annotation in annotations2:
			ax1.annotate(annotation[2],(annotation[0],annotation[1]), ha='left', rotation=45,color=colors[1])

		if drawLegend:
			# reordering the labels
			handles, labels = plt.gca().get_legend_handles_labels()
			#handles.reverse()
			#labels.reverse()
			leg = plt.legend(handles, labels,loc=legendLoc, prop={'size': legendSize})
			for i in range(len(leg.get_lines())):
				leg.get_lines()[i].set_linewidth(linewidth)
	
		#if valuesXticks!=None:
		#	plt.xticks(x,valuesXticks[:len(x)],fontsize=tickSize)
		ax1.grid(False)
		ax1.spines[['right', 'top']].set_visible(False)
		fig.savefig(plotFileName,bbox_inches='tight')
		plt.close('all')

	def coverageplot(times,labels,plotFileName,colors,topPlot,linestyles,title='Coverage analysis',yAxisLabel="Percentage correct",legendSize=14,tickSize=12,labelSize=18,drawLegend=True,legendLoc='upper left',linewidth=3,thickAxis=True):
		values=times
		grain=len(times[0])
		widthBin=100/float(grain)
		x = list(range(grain))
		for xi in x:
			x[xi]=x[xi]*widthBin+widthBin/2

		wids=[]
		for xi in x:
			wids.append(xi/4)
		fig = plt.figure(figsize=(15, 9))
		ax1 = fig.add_subplot(111)
		ax1.set_axisbelow(True)
		ax1.grid(False)
		ax1.set_title(title, fontsize=24)
		ax1.set_xlabel(topPlot, fontsize=labelSize)
		ax1.set_ylabel(yAxisLabel, fontsize=labelSize)
		#if valuesXticks==None:
		plt.xticks(fontsize=tickSize)
		plt.yticks(fontsize=tickSize)
		ax1.set_ylim([0, 100])
		ax1.set_xlim([0, 100])

		for i in range(len(labels)):
			if linestyles!=None:
				ax1.plot(x, values[i], marker='.', c=colors[i], label=labels[i], linestyle=linestyles[i],linewidth=linewidth)
			else:
				ax1.plot(x, values[i], marker='.', c=colors[i], label=labels[i],linewidth=linewidth)

		if drawLegend:
			# reordering the labels
			handles, labels = plt.gca().get_legend_handles_labels()
			handles.reverse()
			labels.reverse()
			leg = plt.legend(handles, labels,loc=legendLoc, prop={'size': legendSize})
			for i in range(len(leg.get_lines())):
				leg.get_lines()[i].set_linewidth(linewidth)
	
		ax1.grid(False)
		ax1.spines[['right', 'top']].set_visible(False)

		if thickAxis:
			for axis in ['top','bottom','left','right']:
				ax1.spines[axis].set_linewidth(2.5)
			# increase tick width
			ax1.tick_params(width=2.2)

		fig.savefig(plotFileName,bbox_inches='tight')
		plt.close('all')



	#make plots
	plotTypes=valueNames
	titles=valueNames
	nans=[]
	for i in range(20):
		nans.append(float("NaN"))
	speeds=methodOptionsNames
	colors=[["blue"],["red","darkorange","olive","purple","green","black","brown"]]
	linestyles=[   [(0, ())]   ,  [(0, (1, 1))   ,   (0, (1, 10))    ,    (5, (10, 3))   ,   (0, (5, 10))    ,   (0, (5, 5))    ,    (0, (5, 1))    ,   (0, (3, 10, 1, 10)) ]   ] # ,    (0, (3, 5, 1, 5))    ,    (0, (3, 5, 1, 5, 1, 5)) 
	linestyles=[   [(0, ())]   ,  [(0, ()) , (0, ()) , (0, ()) , (0, ()) , (0, ()) , (0, ()) , (0, ())]    ] 
	indecesVeryShort=[True,False,True,True,True,False,True,False]

	#create lists of method names etc for the plotting
	namesAll=[]
	colorsAll=[]
	linestylesAll=[]
	for method in range(len(methods)):
		for option in range(len(methodOptionsNames[method])):
			namesAll.append(methodOptionsNamesLegend[method][option])
			colorsAll.append(colors[method][option])
			linestylesAll.append(linestyles[method][option])

	onlyAUROC=False
	if onlyAUROC:
		#AUROC analysis plots
		from sklearn.metrics import roc_curve
		from sklearn.metrics import roc_auc_score
		from sklearn.metrics import precision_recall_curve
		from sklearn.metrics import auc
		for i in range(len(leavesRange)):
			j=leavesRange[i]
			folder=figureFolder
			probList=[]
			truthList=[]
			for method in range(len(methods)):
				probList.append([])
				truthList.append([])
				for option in range(len(methodOptionsNames[method])):
					probList[method].append([])
					truthList[method].append([])
					if j<=maxNumSamples[method][option]:
						file=open(folder+"supports_correct_"+str(j)+"subsamples_"+methods[method]+methodOptionsNames[method][option]+".txt")
						for nRepl in range(20):
							line=file.readline().split()
							for l in line:
								if method:
									#probList[method][option].append(1.0-float(l))
									probList[method][option].append(float(l))
								else:
									probList[method][option].append(float(l))
								truthList[method][option].append(True)
						file.close()
						file=open(folder+"supports_wrong_"+str(j)+"subsamples_"+methods[method]+methodOptionsNames[method][option]+".txt")
						for nRepl in range(20):
							line=file.readline().split()
							for l in line:
								if method:
									#probList[method][option].append(1.0-float(l))
									probList[method][option].append(float(l))
								else:
									probList[method][option].append(float(l))
								truthList[method][option].append(False)
						file.close()

			fig = plt.figure(figsize=(15, 9))
			ax1 = fig.add_subplot(111)
			ax1.set_axisbelow(True)
			ax1.grid(False)
			ax1.spines[['right', 'top']].set_visible(False)
			for axis in ['top','bottom','left','right']:
				ax1.spines[axis].set_linewidth(2.5)
			ax1.tick_params(width=2.2)
			plt.xticks(fontsize=28)
			plt.yticks(fontsize=28)
			ax1.set_xlabel('False Positive Rate', fontsize=34)
			ax1.set_ylabel('True Positive Rate', fontsize=34)

			for method in range(len(methods)):
				for option in range(len(methodOptionsNames[method])):
					if j<=maxNumSamples[method][option]:
						fpr, tpr, thresholds = roc_curve(truthList[method][option], probList[method][option])
						auc_score=roc_auc_score(truthList[method][option], probList[method][option])
						labelAUC=methodOptionsNamesLegend[method][option]+" "+"{0:.2g}".format(auc_score)
						ax1.plot(fpr, tpr, c=colors[method][option], label=labelAUC,linewidth=3)

			handles, labels = plt.gca().get_legend_handles_labels()
			handles.reverse()
			labels.reverse()
			leg = plt.legend(handles, labels,loc="lower right", prop={'size': 25})
			for i in range(len(leg.get_lines())):
				leg.get_lines()[i].set_linewidth(3)
			plt.savefig(figureFolder+"figure_ROC_"+str(j)+"leaves.pdf",bbox_inches='tight')
			plt.close('all')

			fig2 = plt.figure(figsize=(15, 9))
			ax2 = fig2.add_subplot(111)
			ax2.set_axisbelow(True)
			ax2.grid(False)
			ax2.spines[['right', 'top']].set_visible(False)
			for axis in ['top','bottom','left','right']:
				ax2.spines[axis].set_linewidth(2.5)
			ax2.tick_params(width=2.2)
			plt.xticks(fontsize=28)
			plt.yticks(fontsize=28)
			ax2.set_xlabel('Recall', fontsize=34)
			ax2.set_ylabel('Precision', fontsize=34)
			ax2.set_ylim([0.4, 1.0])

			for method in range(len(methods)):
				for option in range(len(methodOptionsNames[method])):
					if j<=maxNumSamples[method][option]:
						precision, recall, thresholds = precision_recall_curve(truthList[method][option], probList[method][option])
						auc_score = auc(recall, precision)
						labelAUC=methodOptionsNamesLegend[method][option]+" "+"{0:.2g}".format(auc_score)
						ax2.plot(recall, precision, c=colors[method][option], label=labelAUC,linewidth=3)

			handles, labels = plt.gca().get_legend_handles_labels()
			handles.reverse()
			labels.reverse()
			leg = plt.legend(handles, labels,loc="lower right", prop={'size': 25})
			for i in range(len(leg.get_lines())):
				leg.get_lines()[i].set_linewidth(3)
			plt.savefig(figureFolder+"figure_PRC_"+str(j)+"leaves.pdf",bbox_inches='tight')
			plt.close('all')
		exit()		

	if not args.onlyNePlots:
	
		for value in range(len(valueNames)):
			#focusing only on plots in the manuscript
			if value==2 or value==3 or value==5 or value==6 or value==7 or value==8 :
				continue
			plotValues=[]
			for i in range(len(leavesRange)):
				j=leavesRange[i]
				plotValues.append([])
				folder=figureFolder
				for method in range(len(methods)):
					for option in range(len(methodOptionsNames[method])):
						if j<=maxNumSamples[method][option]:
							valueList=[]
							if value==0:
								file=open(folder+"supports_correct_"+str(j)+"subsamples_"+methods[method]+methodOptionsNames[method][option]+".txt")
								for nRepl in range(20):
									line=file.readline().split()
									supps=[]
									for l in line:
										supps.append(float(l))
									if supps:
										valueList.append(sum(supps)/len(supps))
								file.close()
							elif value==1:
								file=open(folder+"supports_wrong_"+str(j)+"subsamples_"+methods[method]+methodOptionsNames[method][option]+".txt")
								for nRepl in range(20):
									line=file.readline().split()
									supps=[]
									for l in line:
										supps.append(float(l))
									if supps:
										valueList.append(sum(supps)/len(supps))
								file.close()
							elif value==2:
								file=open(folder+"supports_correct_"+str(j)+"subsamples_"+methods[method]+methodOptionsNames[method][option]+".txt")
								for nRepl in range(20):
									line=file.readline().split()
									for l in line:
										valueList.append(float(l))
								file.close()
							elif value==3:
								file=open(folder+"supports_wrong_"+str(j)+"subsamples_"+methods[method]+methodOptionsNames[method][option]+".txt")
								for nRepl in range(20):
									line=file.readline().split()
									for l in line:
										valueList.append(float(l))
								file.close()
							elif value==4:
								file=open(folder+"supports_correct_"+str(j)+"subsamples_"+methods[method]+methodOptionsNames[method][option]+".txt")
								for nRepl in range(20):
									line=file.readline().split()
									supps=[]
									for l in line:
										supps.append(float(l))
									if supps:
										valueList.append(sum(supps)/len(supps))
								file.close()
								file=open(folder+"supports_wrong_"+str(j)+"subsamples_"+methods[method]+methodOptionsNames[method][option]+".txt")
								index=0
								for nRepl in range(20):
									line=file.readline().split()
									supps=[]
									for l in line:
										supps.append(float(l))
									if supps:
										valueList[index]-=(sum(supps)/len(supps))
										index+=1
								file.close()
							elif value==5:
								file=open(folder+"supports_correct_"+str(j)+"subsamples_"+methods[method]+methodOptionsNames[method][option]+".txt")
								for nRepl in range(20):
									line=file.readline().split()
									for l in line:
										valueList.append(float(l))
								file.close()
								file=open(folder+"supports_wrong_"+str(j)+"subsamples_"+methods[method]+methodOptionsNames[method][option]+".txt")
								for nRepl in range(20):
									line=file.readline().split()
									for l in line:
										valueList.append(1.0-float(l))
								file.close()
							elif value==6:
								file=open(folder+"supports_correct_"+str(j)+"subsamples_"+methods[method]+methodOptionsNames[method][option]+".txt")
								for nRepl in range(20):
									line=file.readline().split()
									if line:
										valueList.append(len(line))
								file.close()
							elif value==7:
								file=open(folder+"supports_wrong_"+str(j)+"subsamples_"+methods[method]+methodOptionsNames[method][option]+".txt")
								for nRepl in range(20):
									line=file.readline().split()
									if line:
										valueList.append(len(line))
								file.close()
							elif value==8:
								file1=open(folder+"supports_correct_"+str(j)+"subsamples_"+methods[method]+methodOptionsNames[method][option]+".txt")
								file=open(folder+"supports_wrong_"+str(j)+"subsamples_"+methods[method]+methodOptionsNames[method][option]+".txt")
								for nRepl in range(20):
									line1=file1.readline().split()
									line=file.readline().split()
									if line:
										valueList.append(float(len(line))/float(len(line)+len(line1)))
								file1.close()
								file.close()
							elif value==9:
								file=open(folder+"times_"+str(j)+"subsamples_"+methods[method]+methodOptionsNames[method][option]+".txt")
								line=file.readline().split()
								if line:
									for valueT in line:
										valueList.append(float(valueT))
								file.close()
							elif value==10:
								file=open(folder+"maxMems_"+str(j)+"subsamples_"+methods[method]+methodOptionsNames[method][option]+".txt")
								line=file.readline().split()
								if line:
									for valueT in line:
										valueList.append(float(valueT))
								file.close()

							plotValues[i].append(valueList)
						else:
							plotValues[i].append([])
							#plotValues[i].append([])
			legendLoc='lower right'
			legendSize=28
			tickSize=28
			labelSize=34
			drawLegends=[1,0,0,0,1,1,0,0,0,1,0,0]
			if value==9:
				ymin=10
				ymax=150000
				ncolLegend=2
				legendLoc='upper right'
				#legendSize=26
			elif value==10:
				ymin=50
				ymax=100000
				ncolLegend=1
			elif value==0 or value==1:
				if args.scenario==0:
					ymin=0.75
				else:
					ymin=0.75
					#ymin=0.82
				ymax=1.01
				ncolLegend=1
			else:
				ymin=None
				ymax=None
				ncolLegend=1
			errplot(plotValues,namesAll,figureFolder+"figure_supports_"+valueNames[value]+".pdf",nLeaves,colorsAll,'Number of genomes',linestylesAll,title='',yAxisLabel=yAxisLabels[value],legendSize=legendSize,ymin=ymin,ymax=ymax,tickSize=tickSize,labelSize=labelSize,logY=logYlist[value],drawLegend=drawLegends[value],legendLoc=legendLoc,valuesXticks=nLeavesTextTicks,linewidth=4,thickAxis=True,ncolLegend=ncolLegend)

	exit()

	#coverage analysis plots
	numLeaves=20000
	grain=10
	plotValues=[]
	j=numLeaves
	folder=figureFolder
	namesCoverage=[]
	colorsCoverage=[]
	linestylesCoverage=[]
	for method in range(len(methods)):
		for option in range(len(methodOptionsNames[method])):
			if j<=maxNumSamples[method][option]:
				namesCoverage.append(methodOptionsNamesLegend[method][option])
				colorsCoverage.append(colors[method][option])
				linestylesCoverage.append(linestyles[method][option])
				correctValueList=[0]*grain
				file=open(folder+"supports_correct_"+str(j)+"subsamples_"+methods[method]+methodOptionsNames[method][option]+".txt")
				for nRepl in range(20):
					line=file.readline().split()
					for l in line:
						correctValueList[math.ceil(float(l)*grain)-1]+=1
				file.close()
				wrongValueList=[0]*grain
				file=open(folder+"supports_wrong_"+str(j)+"subsamples_"+methods[method]+methodOptionsNames[method][option]+".txt")
				for nRepl in range(20):
					line=file.readline().split()
					for l in line:
						wrongValueList[math.ceil(float(l)*grain)-1]+=1
				file.close()
				proportionList=[]
				for i in range(grain):
					if correctValueList[i]+wrongValueList[i]>0:
						proportionList.append(100*float(correctValueList[i])/(correctValueList[i]+wrongValueList[i]))
					else:
						proportionList.append(math.nan)
				
				plotValues.append(proportionList)
				print(namesCoverage[-1])
				print(correctValueList)
				print(wrongValueList)
				print(proportionList)

	legendLoc='upper left'
	legendSize=28
	tickSize=28
	labelSize=34
	coverageplot(plotValues,namesCoverage,figureFolder+"figure_coverage_"+str(numLeaves)+"leaves.pdf",colorsCoverage,'Branch support',linestylesCoverage,title='',yAxisLabel="Percentage correct",legendSize=legendSize,tickSize=tickSize,labelSize=labelSize,drawLegend=True,legendLoc=legendLoc,linewidth=4)
				
				



exit()


















