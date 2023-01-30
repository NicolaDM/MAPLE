import sys
import os.path
from time import time
import argparse
import random
import time

#©EMBL-European Bioinformatics Institute, 2021
folderName="/Users/demaio/Desktop/GISAID-hCoV-19-phylogeny-2021-03-12"
folderCluster="/nfs/research/goldman/demaio"
pypy3Path="/hps/software/users/goldman/pypy3/pypy3.7-v7.3.5-linux64/bin/pypy3"
iqTreePath="/hps/software/users/goldman/iqtree-2.1.3-Linux/bin/iqtree2"
minicondaPath="/hps/software/users/goldman/nicola/miniconda3/"
fastTreePath="/hps/software/users/goldman/FastTree2/FastTree"
# Running MAPLE simulations with single script:
# downloaded tree from http://hgdownload.soe.ucsc.edu/goldenPath/wuhCor1/UShER_SARS-CoV-2/
# Run phastSim, once for each simulation scenario (so far 3, normal, 4 catrgories, and gamma-distributed rate variation):
# cd folderName/simulations/phastSim
# python3 bin/phastSim --outpath folderName/simulations/ --treeFile folderName/simulations/public-latest.all.nwk --scale 0.00003344 --reference MN908947.3.fasta --categoryProbs 0.25 0.25 0.25 0.25 --categoryRates 0.1 0.5 1.0 2.0 --outputFile phastSim_genomes_4categories_UNREST --createNewick 
# python3 bin/phastSim --outpath folderName/simulations/ --treeFile folderName/simulations/public-latest.all.nwk --scale 0.00003344 --reference MN908947.3.fasta --outputFile phastSim_genomes_UNREST --createNewick 
# python3 bin/phastSim --outpath folderName/simulations/ --treeFile folderName/simulations/public-latest.all.nwk --scale 0.00003344 --reference MN908947.3.fasta --alpha 0.1 --outputFile phastSim_genomes_alpha_UNREST --createNewick 
# Then move files to the cluster;

# Now run the script for creating the simulations inputs:
# scp folderName/*.py  folderCluster/fastLK/code
# bsub -M 20000 pypy3Path MAPLE_benchmarking.py --createTotalData
# pypy3Path folderCluster/fastLK/code/MAPLE_benchmarking.py --createBashScript
# sh folderCluster/fastLK/simulationsNew/createSubsampleInputFiles.sh

# Run MAPLE version and option testing: first modify the flags in appropriately, then run 
# pypy3Path folderCluster/fastLK/code/MAPLE_benchmarking.py --createBashScript
# sh folderCluster/fastLK/simulationsNew/submitMAPLE0.1.8.sh
# etc
# sh folderCluster/fastLK/simulationsNew/submitRF_MAPLE0.1.8.sh
# etc
# sh folderCluster/fastLK/simulationsNew/submitIQtreeLK_MAPLE0.1.4.sh
# etc
# sh folderCluster/fastLK/simulationsNew/submitLK_MAPLE0.1.5.sh
# etc
# pypy3Path folderCluster/fastLK/code/MAPLE_benchmarking.py --collectResults

# Running proper benchmarking:
# first run bash script creation for each sample number with:
# pypy3Path folderCluster/fastLK/code/MAPLE_benchmarking.py --createBashScript --numSamples 100000
# and then for each number of samples first run the methods:
# sh folderCluster/fastLK/simulationsNew/submitUShER.sh ; sh folderCluster/fastLK/simulationsNew/submitFastTree.sh ; sh folderCluster/fastLK/simulationsNew/submitIQtree.sh ; sh folderCluster/fastLK/simulationsNew/submitRAxML.sh ; sh folderCluster/fastLK/simulationsNew/submitRAxML-NG.sh ; sh folderCluster/fastLK/simulationsNew/submitMaple.sh

# when usher is finished also run matOptimize on the UShER output:
# sh folderCluster/fastLK/simulationsNew/submitmatOptimize.sh 
# sh folderCluster/fastLK/simulationsNew/submitMatOptimizeConversion.sh

# when methods are finished on a certain number of samples, run the post-run processing (remembering to first modify the bash scripts to inclde the right number of samples with MAPLE_benchmarking.py):
# sh folderCluster/fastLK/simulationsNew/submitIQtreeLK_UShER.sh ; sh folderCluster/fastLK/simulationsNew/submitMapleLK_UShER.sh ; sh folderCluster/fastLK/simulationsNew/submitRF_UShER.sh ; sh folderCluster/fastLK/simulationsNew/submitParsimony_UShER.sh
# sh folderCluster/fastLK/simulationsNew/submitIQtreeLK_matOptimize.sh ; sh folderCluster/fastLK/simulationsNew/submitMapleLK_matOptimize.sh ; sh folderCluster/fastLK/simulationsNew/submitRF_matOptimize.sh ; sh folderCluster/fastLK/simulationsNew/submitParsimony_matOptimize.sh
# sh folderCluster/fastLK/simulationsNew/submitIQtreeLK_IQtree.sh ; sh folderCluster/fastLK/simulationsNew/submitMapleLK_IQtree.sh 
# sh folderCluster/fastLK/simulationsNew/submitRF_IQtree.sh ; sh folderCluster/fastLK/simulationsNew/submitParsimony_IQtree.sh
# sh folderCluster/fastLK/simulationsNew/submitIQtreeLK_FastTree.sh ; sh folderCluster/fastLK/simulationsNew/submitMapleLK_FastTree.sh
# sh folderCluster/fastLK/simulationsNew/submitRF_FastTree.sh ; sh folderCluster/fastLK/simulationsNew/submitParsimony_FastTree.sh
# sh folderCluster/fastLK/simulationsNew/submitIQtreeLK_RAxML.sh ; sh folderCluster/fastLK/simulationsNew/submitMapleLK_RAxML.sh ; sh folderCluster/fastLK/simulationsNew/submitRF_RAxML.sh ; sh folderCluster/fastLK/simulationsNew/submitParsimony_RAxML.sh
# sh folderCluster/fastLK/simulationsNew/submitIQtreeLK_RAxML-NG.sh ; sh folderCluster/fastLK/simulationsNew/submitMapleLK_RAxML-NG.sh ; sh folderCluster/fastLK/simulationsNew/submitRF_RAxML-NG.sh ; sh folderCluster/fastLK/simulationsNew/submitParsimony_RAxML-NG.sh
# sh folderCluster/fastLK/simulationsNew/submitIQtreeLK_Maple.sh
# sh folderCluster/fastLK/simulationsNew/submitMapleLK_Maple.sh
# sh folderCluster/fastLK/simulationsNew/submitRF_Maple.sh ; sh folderCluster/fastLK/simulationsNew/submitParsimony_Maple.sh

# checking which estimations are still running:
# bjobs -g /IQtree2000
# bjobs -g /RAxML-NG2000
# bjobs -g /RAxML2000
# bjobs -g /FastTree20000
# bjobs -g /Maple100000
# bjobs -g /UShER100000
# bjobs -g /matOptimize10000
# bjobs -g /matConv20000

# To summarize results (still to be done with 2000 - ):
# pypy3Path folderCluster/fastLK/code/MAPLE_benchmarking.py --collectResults

parser = argparse.ArgumentParser(description='Create files for the benchmarking of MAPLE, both input files and cluster shell scripts.')
parser.add_argument('--inputRealData',default=folderCluster+"/fastLK/realData/2021-03-31_unmasked_differences_reduced.txt_consensus-based.txt", help='input real data file; should contain the difference of all samples with respet to the reference.')
parser.add_argument('--inputRealDataReference',default=folderCluster+"/fastLK/realData/2021-03-31_unmasked_differences_reduced.txt_consensus.fa", help='input real data file with the reference genome used.')
parser.add_argument('--pathToSimulationFolder',default=folderCluster+"/fastLK/simulationsNew/", help='path to the phastSim python script.')
parser.add_argument('--inputRealTree',default=folderCluster+"/fastLK/simulationsNew/public-latest.all.nwk", help='path to the (real) input tree for phastSim simulations.')
parser.add_argument('--inputSimulationReference',default=folderCluster+"/fastLK/simulationsNew/MN908947.3.fasta", help='path to the reference genome to be used for phastSim simulations.')

parser.add_argument("--createTotalData", help="Extract realized trees from global ones, and mask simulated alignments similarly to the real one.", action="store_true")

parser.add_argument("--createBashScript", help="Create bash script to run the file creation on the cluster in parallel.", action="store_true")

parser.add_argument("--createInputData", help="Subsample real and simulated data to create input datasets for benchmarking.", action="store_true")
parser.add_argument("--subSampleNum",help="Number of subsamples to extract.",  type=int, default=10)
parser.add_argument("--repeat",help="Which repeat to simulate, typically between 1-10. Only considered for option --createInputData .",  type=int, default=1)
#subsamples of 0) real data; 1) simulations basic scenario; 2) simulated with rate variation; 3) simulated with high (alpha) rate variation; 4) simulations with N's.
parser.add_argument("--scenario",help="Which scenario to subsample, typically between 0-4. Only considered for option --createInputData .",  type=int, default=0)

parser.add_argument("--numSamples",help="Number of samples to consider when running different estimation methods, or when checking the results; typically between 1000 and 1000000. if 0 (default) consider the whole range.",  type=int, default=0)

parser.add_argument("--collectResults", help="Collect the results from the benchmarking.", action="store_true")
parser.add_argument("--methodToFocusOn",help="Choose which method to investigate during the collection of the results (when using option --collectResults). By default, -1, means look at all methods. Otherwise choose method number from 0 to 6.",  type=int, default=-1)

parser.add_argument("--createFigures", help="Create the result files to generate figures from the benchmarking experiment.", action="store_true")

parser.add_argument("--runFigureGeneration", help="Create the figures from the benchmarking experiment - it needs to be used from somewhere with matplotlib.", action="store_true")

parser.add_argument("--runExtraAnalysesDivergence", help="Create files etc just for the extra analysis involving the effect of divergence.", action="store_true")
args = parser.parse_args()


#subsample datasets
subSampleNums=[1000,2000,5000,10000,20000,50000,100000,200000,500000]
subSampleFasta=[1000,2000,5000,10000,20000]

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

#function to read input newick string
def readNewick(nwFile,multipleTrees=False,dirtiness=True):
	phyloFile=open(nwFile)
	trees=[]
	line=phyloFile.readline()
	while line!="":
		while line=="\n":
			line=phyloFile.readline()
		if line=="":
			break
		nwString=line.replace("\n","")
	
		index=0
		node=Tree()
		node.dirty=dirtiness
		name=""
		distStr=""
		finished=False
		while index<len(nwString):
			if nwString[index]=="(":
				newNode=Tree()
				newNode.minorSequences=[]
				newNode.dirty=dirtiness
				node.add_child(newNode)
				newNode.up=node
				node=newNode
				index+=1
			elif nwString[index]==";":
				trees.append(node)
				finished=True
				break
				#return node
			elif nwString[index]=="[":
				while nwString[index]!="]":
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
				newNode=Tree()
				newNode.minorSequences=[]
				newNode.dirty=dirtiness
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
				node=node.up
				index+=1
			else:
				name+=nwString[index]
				index+=1
		if not finished:
			print("Error, final character ; not found in newick string in file "+nwFile+".")
			exit()

		if not multipleTrees:
			break
		line=phyloFile.readline()

	phyloFile.close()
	return trees

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
def createNewick(node):
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
				if nextNode.dist:
					stringList.append("):"+str(nextNode.dist))
				else:
					stringList.append("):"+str(0.0))
				#if nextNode.up!=None:
				direction=1
				lastNode=nextNode
					# if nextNode.up.children[0]==nextNode:
					# 	direction=1
					# else:
					# 	direction=2
				nextNode=nextNode.up
		else:
			#print("Terminal node "+nextNode.name)
			if nextNode.dist:
				stringList.append(nextNode.name+":"+str(nextNode.dist))
			else:
				stringList.append(nextNode.name+":"+str(0.0))
			#if nextNode.up!=None:
			direction=1
			lastNode=nextNode
				# if nextNode.up.children[0]==nextNode:
				# 	direction=1
				# else:
				# 	direction=2
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

#Robinson-Foulds distance (1981) using a simplification of the algorithm from Day 1985.
#this function prepare the data to compare trees to a reference one t1.
#I split in two functions so that I don't have to repeat these steps for the reference tree if I compare to the same reference tree multiple times.
def prepareTreeComparison(t1,rooted=False,minimumBLen=0.000006):
	#dictionary of values given to sequence names
	leafNameDict={}
	#list of sequence names sorted according to value
	leafNameDictReverse=[]
	#table containing clusters in the tree
	nodeTable=[]
	#if comparing as unrooted trees, calculate tot num of leaves (using a postorder traversal), which will become useful later
	if not rooted:
		nLeaves=0
		node=t1
		movingFrom=0
		while node!=t1.up:
			if movingFrom==0: #0 means reaching node from parent, 1 means coming back from a child
				if len(node.children)==0:
					nLeaves+=1
					nextNode=node.up
					movingFrom=1
					nodeTable.append([0,0])
				else:
					nextNode=node.children[0]
					movingFrom=0
					node.exploredChildren=0
			else:
				nChildren=len(node.children)
				node.exploredChildren+=1
				if node.exploredChildren==nChildren:
					nextNode=node.up
					movingFrom=1
				else:
					nextNode=node.children[node.exploredChildren]
					movingFrom=0
			node=nextNode
			
	#implementing a non-recursive postorder traversal to assign values to internal nodes to fill nodeTable
	leafCount=0
	node=t1
	movingFrom=0
	lastL=float("inf")
	lastR=float("-inf")
	lastDesc=0
	numBranches=0
	while node!=t1.up:
		if movingFrom==0: #0 means reaching node from parent, 1 means coming back from a child
			if len(node.children)==0:
				node.name=(node.name).replace("?","_").replace("&","_")
				leafNameDict[node.name]=leafCount
				leafNameDictReverse.append(node.name)
				if rooted:
					nodeTable.append([0,0])
				lastL=leafCount
				lastR=leafCount
				lastDesc=1
				leafCount+=1
				nextNode=node.up
				movingFrom=1
			else:
				node.exploredChildren=0
				node.maxSoFar=float("-inf")
				node.minSoFar=float("inf")
				node.nDescendants=0
				nextNode=node.children[0]
				movingFrom=0
		else:
			nChildren=len(node.children)
			node.exploredChildren+=1
			if lastL<node.minSoFar:
				node.minSoFar=lastL
			if lastR>node.maxSoFar:
				node.maxSoFar=lastR
			node.nDescendants+=lastDesc
			if node.exploredChildren==nChildren:
				nextNode=node.up
				movingFrom=1
				lastL=node.minSoFar
				lastR=node.maxSoFar
				lastDesc=node.nDescendants
				if node==t1:
					nodeTable[lastR][0]=lastL
					nodeTable[lastR][1]=lastR
				else:
					if node.dist>minimumBLen:
						numBranches+=1
						if rooted or lastL>0:
							if node==node.up.children[-1]:
								nodeTable[lastL][0]=lastL
								nodeTable[lastL][1]=lastR
							else:
								nodeTable[lastR][0]=lastL
								nodeTable[lastR][1]=lastR
						else: # re-root at leaf 0, so flip the values for the current branch if it contains leaf 0.
								flippedL=lastR+1
								flippedR=nLeaves-1
								nodeTable[flippedL][0]=flippedL
								nodeTable[flippedL][1]=flippedR
			else:
				nextNode=node.children[node.exploredChildren]
				movingFrom=0
		node=nextNode
	return leafNameDict, nodeTable, leafCount, numBranches

#Robinson-Foulds distance (1981) using a simplification of the algorithm from Day 1985.
#this function compares the current tree t2 to a previous one for which prepareTreeComparison() was run.
def RobinsonFouldsWithDay1985(t2,leafNameDict,nodeTable,leafCount,numBranches,rooted=False,minimumBLen=0.000006):
	#implementing a non-recursive postorder traversal to check branch existance in the reference tree
	node=t2
	#branches in reference tree that are also in t2
	foundBranches=0
	#branches in t2 that are not found in the reference
	missedBranches=0
	movingFrom=0
	lastL=float("inf")
	lastR=float("-inf")
	lastDesc=0
	visitedLeaves=0
	while node!=t2.up:
		if movingFrom==0: #0 means reaching node from parent, 1 means coming back from a child
			if len(node.children)==0:
				node.name=(node.name).replace("?","_").replace("&","_")
				if node.name in leafNameDict:
					leafNum=leafNameDict[node.name]
				else:
					print(node.name+" not in reference tree - aborting RF distance")
					return None, None, None, None, None, None
				lastL=leafNum
				lastR=leafNum
				lastDesc=1
				nextNode=node.up
				movingFrom=1
				visitedLeaves+=1
			else:
				node.exploredChildren=0
				node.maxSoFar=float("-inf")
				node.minSoFar=float("inf")
				node.nDescendants=0
				nextNode=node.children[0]
				movingFrom=0
		else:
			nChildren=len(node.children)
			node.exploredChildren+=1
			if lastL<node.minSoFar:
				node.minSoFar=lastL
			if lastR>node.maxSoFar:
				node.maxSoFar=lastR
			node.nDescendants+=lastDesc
			if node.exploredChildren==nChildren:
				nextNode=node.up
				movingFrom=1
				lastL=node.minSoFar
				lastR=node.maxSoFar
				lastDesc=node.nDescendants
				if node!=t2:
					if node.dist>minimumBLen:
						if (lastR+1-lastL)==lastDesc:
							if rooted or lastL>0:
								if nodeTable[lastL][0]==lastL and nodeTable[lastL][1]==lastR:
									foundBranches+=1
								elif nodeTable[lastR][0]==lastL and nodeTable[lastR][1]==lastR:
									foundBranches+=1
								else:
									missedBranches+=1
							else: # re-root at leaf 0, so flip the values for the current branch if it contains leaf 0.
								flippedL=lastR+1
								flippedR=leafCount-1
								if nodeTable[flippedL][0]==flippedL and nodeTable[flippedL][1]==flippedR:
									foundBranches+=1
								elif nodeTable[flippedR][0]==flippedL and nodeTable[flippedR][1]==flippedR:
									foundBranches+=1
								else:
									missedBranches+=1
						else:
							missedBranches+=1
			else:
				nextNode=node.children[node.exploredChildren]
				movingFrom=0
		node=nextNode
	if visitedLeaves<leafCount:
		print("There are leaves in the reference that have not been found in this new tree")
		return None, None, None, None, None, None
	#first value is number of differences, second value is max number of differences just in case one wants the normalized values; 
	#the other values are there just in case on wants more detail.
	numDiffs=((numBranches-foundBranches)+missedBranches)
	return numDiffs, float(numDiffs)/(2*(leafCount-3)), leafCount, foundBranches, missedBranches, (numBranches-foundBranches)



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

#assign ambiguities from leafAmb to tip "leaf"
def applyAmbiguities(leaf,leafAmb):
	#count number of isolated ambiguities and create subvector of leafAmb with just long stretches of "N"s.
	onlyAmb=[]
	numAmb=0
	lastPos=1
	for entry in leafAmb:
		if (len(entry)==2 or entry[2]==1) and (entry[0]=="N" or (not entry[0] in allelesList)):
			numAmb+=1
		elif entry[0]=="N" and len(entry)>2 and entry[2]>1:
			if entry[1]>lastPos:
				onlyAmb.append(["R",lastPos,entry[1]-lastPos])
			onlyAmb.append(entry)
			lastPos=entry[1]+entry[2]
	if lastPos<=len(ref)-12:
		onlyAmb.append(["R",lastPos,len(ref)-12+1-lastPos])
	#count number of non-reference entries in "leaf", and sample those to me masked.
	numDiff=0
	lastPos=1
	leafR=[]
	for entry in leaf:
		if entry[1]<=len(ref)-12:
			if (len(entry)<3 or entry[2]==1) and entry[0]!="N":
				numDiff+=1
			if entry[1]>lastPos:
				leafR.append(["R",lastPos,entry[1]-lastPos])
			leafR.append(entry)
			if len(entry)==3:
				lastPos=entry[1]+entry[2]
			else:
				lastPos=entry[1]+1
	if lastPos<=len(ref)-12:
		leafR.append(["R",lastPos,len(ref)-12+1-lastPos])
	if numAmb>=numDiff:
		masked=range(numAmb)
	else:
		masked=random.sample(range(numDiff), numAmb)

	#now create "newLeaf" by masking "leaf"
	indexEntry1=0
	indexEntry2=0
	pos=1
	entry1=leafR[indexEntry1]
	pos1=entry1[1]
	if entry1[0]!="N" and entry1[0]!="R":
		end1=pos1
	else:
		end1=pos1+entry1[2]-1

	entry2=onlyAmb[indexEntry2]
	pos2=entry2[1]
	end2=pos2+entry2[2]-1
	
	end=min(end1,end2)
	length=end+1-pos
	newLeaf=[]
	diffIndex=-1
	while True:
		if entry1[0]!="N" and entry1[0]!="O" and entry1[0]!="R":
			diffIndex+=1
			if diffIndex in masked or entry2[0]=="N":
				newLeaf.append(["N",pos,length])
			else:
				newLeaf.append([entry1[0],pos,length])
		elif entry2[0]=="N":
			newLeaf.append(["N",pos,length])
		else:
			newLeaf.append([entry1[0],pos,length])

		pos+=length
		if pos>lRef-12:
			break
		if pos>end1:
			indexEntry1+=1
			entry1=leafR[indexEntry1]
			pos1=entry1[1]
			if entry1[0]!="N" and entry1[0]!="R":
				end1=pos1
			else:
				end1=pos1+entry1[2]-1
		if pos>end2:
			indexEntry2+=1
			entry2=onlyAmb[indexEntry2]
			pos2=entry2[1]
			if entry2[0]!="N" and entry2[0]!="R":
				end2=pos2
			else:
				end2=pos2+entry2[2]-1
		end=min(end1,end2)
		length=end+1-pos

	return newLeaf

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
	#lRef=len(ref)
	#print("Ref genome length: "+str(lRef))
	file.close()
	return ref

if not args.runFigureGeneration:
	ref = collectReference(args.inputRealDataReference)
	lRef=len(ref)
	simuRef = collectReference(args.inputSimulationReference)
	lSimuRef=len(simuRef)

#create the input data for the benchmarking analysis
if args.createTotalData:
	random.seed(a=1)
	runExtraAnalysesDivergence=True
	if args.runExtraAnalysesDivergence:
		divergences=["01X","02X","05X","1X","2X","5X","10X","20X","50X","100X","200X","500X","1000X"]
		alignments=[]
		for div in divergences:
			alignments.append(args.pathToSimulationFolder+"phastSim_genomes_"+div+".txt")

	else:
		alignments=[args.pathToSimulationFolder+"phastSim_genomes_UNREST.txt",args.pathToSimulationFolder+"phastSim_genomes_4categories_UNREST.txt",args.pathToSimulationFolder+"phastSim_genomes_alpha_UNREST.txt"]
	
	#create realized trees, collapsing branches without mutations
	for alignment in alignments:
		print("Preparing mutation-informed simulated global tree from "+alignment)
		file=open(alignment.replace("txt","tree"))
		treeLine=file.readline()
		treeLine=treeLine.replace("&mutations={}","")
		charList=[]
		character=treeLine[0]
		index=0
		while character!="\n":
			if character=="[":
				index+=1
				if treeLine[index]=="]":
					while treeLine[index]!="," and treeLine[index]!=")" and treeLine[index]!=";":
						index+=1
					character=treeLine[index]
					charList.append(":0.0")
				else:
					while treeLine[index]!="]":
						index+=1
					index+=1
					character=treeLine[index]
			else:
				charList.append(character)
				index+=1
				character=treeLine[index]
		charList.append("\n")
		newTreeString="".join(charList)
		file=open(alignment.replace(".txt","_realizedTree.tree"),"w")
		file.write(newTreeString)
		file.close()
		print("Extracted realized tree")

		#apply N's to the simulated data 
		diffFile=(alignments[0].replace(".txt",""))+"_Ns.txt"
		data=readConciseAlignment(alignments[0],shift01pos=True)
		dataReal=readConciseAlignmentReal(args.inputRealData)
		samples=data.keys()
		print(str(len(samples))+" sequences in the concise DNA data file; creating simulated alignment with N's.")
		#Now create version of the data with ambiguous characters
		keys2=list(dataReal.keys())
		count=0
		fileO=open(diffFile,"w")
		for name1 in samples:
			fileO.write(">"+name1+"\n")
			i2=random.randint(0,len(keys2)-1)
			name2=keys2[i2]
			newLeaf=applyAmbiguities(data[name1],dataReal[name2])
			for m in newLeaf:
				if m[0]!="R":
					if len(m)==2:
						fileO.write(m[0]+"\t"+str(m[1])+"\n")
					else:
						fileO.write(m[0]+"\t"+str(m[1])+"\t"+str(m[2])+"\n")
			count+=1
			if count==1:
				print(data[name1])
				print(dataReal[name2])
				print(newLeaf)
			data[name1]=None
			if (count%100000)==0:
				print(count)
		fileO.close()
		del dataReal
		del data


#create folders
if args.createBashScript:
	if args.runExtraAnalysesDivergence:
		folderNameSimu=args.pathToSimulationFolder+'divergence/'
		if not os.path.isdir(folderNameSimu):
			os.mkdir(folderNameSimu)
		divergences=["01X","02X","05X","1X","2X","5X","10X","20X","50X","100X","200X","500X","1000X"]
		for div in divergences:
			folderNameSimu=args.pathToSimulationFolder+'divergence/'+div+"/"
			if not os.path.isdir(folderNameSimu):
				os.mkdir(folderNameSimu)
			for i in range(10):
				folderNameSimu=args.pathToSimulationFolder+'divergence/'+div+"/output_repl"+str(i+1)+"/"
				if not os.path.isdir(folderNameSimu):
					os.mkdir(folderNameSimu)
	else:
		folders=["realDataSubsamples","simulationsSubsamples","simulations4catSubsamples","simulationsAlphaSubsamples","simulationsNsSubsamples"]
		#creating folders for output files
		for scenario in range(len(folders)):
			folder=folders[scenario]
			for j in [1000, 2000, 5000, 10000, 20000, 50000, 100000, 200000, 500000]:
				for i in range(10):
					folderNameSimu=args.pathToSimulationFolder+folder+'/'+str(j)+"subsamples/output_repl"+str(i+1)+"/"
					if not os.path.isdir(folderNameSimu):
						os.mkdir(folderNameSimu)

	




#create subsample files 
if args.createInputData:

	if args.runExtraAnalysesDivergence:
		ref = collectReference(args.inputSimulationReference)
		scenario=args.scenario
		subsampleTreeInference=args.subSampleNum
		divergences=["01X","02X","05X","1X","2X","5X","10X","20X","50X","100X","200X","500X","1000X"]
		folder="divergence"
		seed=args.repeat
		repeat=seed

		#alignments=[args.inputRealData, args.pathToSimulationFolder+"phastSim_genomes_UNREST.txt", args.pathToSimulationFolder+"phastSim_genomes_4categories_UNREST.txt", args.pathToSimulationFolder+"phastSim_genomes_alpha_UNREST.txt", args.pathToSimulationFolder+"phastSim_genomes_UNREST_Ns.txt"]
		#phylos=["",args.pathToSimulationFolder+"phastSim_genomes_UNREST_realizedTree.tree", args.pathToSimulationFolder+"phastSim_genomes_4categories_UNREST_realizedTree.tree", args.pathToSimulationFolder+"phastSim_genomes_alpha_UNREST_realizedTree.tree", args.pathToSimulationFolder+"phastSim_genomes_UNREST_realizedTree.tree"]
		phylos=[args.pathToSimulationFolder+"phastSim_genomes_"+divergences[scenario]+"_realizedTree.tree"]

		pathToRepeat=args.pathToSimulationFolder+folder+'/'+divergences[scenario]+"/repeat"+str(repeat)+"_"+divergences[scenario]+"_"+str(subsampleTreeInference)+"samples"
		print("Considering scenario "+divergences[scenario])
		if not os.path.isdir(args.pathToSimulationFolder+folder):
			os.mkdir(args.pathToSimulationFolder+folder)
		data=readConciseAlignment(args.pathToSimulationFolder+"phastSim_genomes_"+divergences[scenario]+".txt",shift01pos=True)
		samples=data.keys()
		print(str(len(samples))+" sequences in the total concise DNA data file")

		print("Creating subsample files for "+divergences[scenario])
		subsample=args.subSampleNum
		random.seed(seed)
		newSamples=random.sample(samples,subsampleTreeInference)

		#create subsampled MAPLE file
		diffFile=pathToRepeat+".txt"
		fileO=open(diffFile,"w")
		numSeq=0
		for s in newSamples:
			numSeq+=1
			if numSeq==1:
				name1=s
				print(name1)
			if numSeq==2:
				name2=s
				print(name2)
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
		print("Fasta file created")

		#create initial (dummy) tree for USheER
		file=open(pathToRepeat+"_initialTree.tree","w")
		file.write("("+name1+":10,"+name2+":10):1;\n")
		file.close()
		#create VCF file for USheER
		newFastaFile=pathToRepeat+"_withRef.fa"
		vcfFile=pathToRepeat+".vcf"
		os.system("cat "+args.inputSimulationReference+" "+fastaFile+" > "+newFastaFile+"\n")
		os.system(""+minicondaPath+"envs/usher-env/bin/faToVcf "+newFastaFile+" "+vcfFile+"\n")
		print("VCF file created")

		#extract subtree of the larger simulated tree
		# first create search tree for current sample names:
		print("Total number of leaves to be found: "+str(len(newSamples)))
		binTree=BinarySearchTree()
		for s in newSamples:
			binTree.put(s)
		
		phylo = readNewick(phylos[0])[0]
		nodeList=postorderList(phylo)

		#extract subtree for subsample
		numNode=0
		numDescList=[0,0,0,0,0]
		for node in nodeList:
			numNode+=1
			#print("numNode "+str(numNode))
			if len(node.children)==0:
				numDescList[0]+=1
				#print(node.name)
				if binTree[node.name]:
					node.subtree=Tree(name=node.name,dist=node.dist)
				else:
					node.subtree=None
			else:
				if len(node.children)<4:
					numDescList[len(node.children)]+=1
				else:
					numDescList[-1]+=1
				#print(len(node.children))
				numDesc=0
				for c in node.children:
					if c.subtree!=None:
						numDesc+=1
						child=c
				if numDesc==0:
					node.subtree=None
				elif numDesc==1:
					node.subtree=child.subtree
					node.subtree.dist+=node.dist
				else:
					node.subtree=Tree(dist=node.dist)
					for c in node.children:
						if c.subtree:
							node.subtree.add_child(c.subtree)
							c.subtree.up=node.subtree
							
		newickString=createNewick(phylo.subtree)
		phylo=None
		del data
		phyloFile=open(pathToRepeat+"_realized.nw","w")
		phyloFile.write(newickString+"\n")
		phyloFile.close()
		print("Subtree extracted and written to file "+pathToRepeat+"_realized.nw")

		

	else:
		scenario=args.scenario
		subsampleTreeInference=args.subSampleNum
		#subsamples of 1) real data; 2) simulations basic scenario; 3) simulated with rate variation; 4) simulated with high (alpha) rate variation; 5) simulations with N's.
		folders=["realDataSubsamples","simulationsSubsamples","simulations4catSubsamples","simulationsAlphaSubsamples","simulationsNsSubsamples"]
		alignments=[args.inputRealData, args.pathToSimulationFolder+"phastSim_genomes_UNREST.txt", args.pathToSimulationFolder+"phastSim_genomes_4categories_UNREST.txt", args.pathToSimulationFolder+"phastSim_genomes_alpha_UNREST.txt", args.pathToSimulationFolder+"phastSim_genomes_UNREST_Ns.txt"]
		phylos=["",args.pathToSimulationFolder+"phastSim_genomes_UNREST_realizedTree.tree", args.pathToSimulationFolder+"phastSim_genomes_4categories_UNREST_realizedTree.tree", args.pathToSimulationFolder+"phastSim_genomes_alpha_UNREST_realizedTree.tree", args.pathToSimulationFolder+"phastSim_genomes_UNREST_realizedTree.tree"]
		seed=args.repeat
		repeat=seed

		if scenario==0:
			ref = collectReference(args.inputRealDataReference)
		else:
			ref = collectReference(args.inputSimulationReference)
		folder=folders[scenario]
		pathToRepeat=args.pathToSimulationFolder+folder+'/'+str(subsampleTreeInference)+"subsamples/repeat"+str(repeat)+"_"+str(subsampleTreeInference)+"samples_"+folders[scenario]
		print("Considering scenario "+folder)
		if not os.path.isdir(args.pathToSimulationFolder+folder):
			os.mkdir(args.pathToSimulationFolder+folder)
		if scenario==0 or scenario==4:
			data=readConciseAlignmentReal(alignments[scenario])
		else:
			data=readConciseAlignment(alignments[scenario],shift01pos=True)
		samples=data.keys()
		print(str(len(samples))+" sequences in the total concise DNA data file")

		print("Creating subsample files for "+str(subsampleTreeInference)+" subsamples")
		if not os.path.isdir(args.pathToSimulationFolder+folder+'/'+str(subsampleTreeInference)+"subsamples/"):
			os.mkdir(args.pathToSimulationFolder+folder+'/'+str(subsampleTreeInference)+"subsamples/")
		subsample=subsampleTreeInference
		random.seed(seed)
		newSamples=random.sample(samples,subsampleTreeInference)

		#create subsampled MAPLE file
		diffFile=pathToRepeat+".txt"
		fileO=open(diffFile,"w")
		numSeq=0
		for s in newSamples:
			numSeq+=1
			if numSeq==1:
				name1=s
				print(name1)
			if numSeq==2:
				name2=s
				print(name2)
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
		#if subsampleTreeInference in subSampleFasta:
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
		print("Fasta file created")
		#del data

		#create initial (dummy) tree for USheER
		file=open(pathToRepeat+"_initialTree.tree","w")
		file.write("("+name1+":10,"+name2+":10):1;\n")
		file.close()
		#create VCF file for USheER
		newFastaFile=pathToRepeat+"_withRef.fa"
		vcfFile=pathToRepeat+".vcf"
		if scenario==0:
			os.system("cat "+args.inputRealDataReference+" "+fastaFile+" > "+newFastaFile+"\n")
		else:
			os.system("cat "+args.inputSimulationReference+" "+fastaFile+" > "+newFastaFile+"\n")
		os.system(""+minicondaPath+"envs/usher-env/bin/faToVcf "+newFastaFile+" "+vcfFile+"\n")
		print("VCF file created")

		#extract subtree of the larger simulated tree
		if scenario>0:
			#create search tree for current sample names:
			print("Total number of leaves to be found: "+str(len(newSamples)))
			binTree=BinarySearchTree()
			for s in newSamples:
				binTree.put(s)

			phylo = readNewick(phylos[scenario])[0]
			nodeList=postorderList(phylo)

			#extract subtree for subsample
			numNode=0
			numDescList=[0,0,0,0,0]
			for node in nodeList:
				numNode+=1
				#print("numNode "+str(numNode))
				if len(node.children)==0:
					numDescList[0]+=1
					#print(node.name)
					if binTree[node.name]:
						node.subtree=Tree(name=node.name,dist=node.dist)
						#if len(newSamples)<15:
						#	print("found "+node.name)
					else:
						node.subtree=None
				else:
					if len(node.children)<4:
						numDescList[len(node.children)]+=1
					else:
						numDescList[-1]+=1
					#print(len(node.children))
					numDesc=0
					for c in node.children:
						if c.subtree!=None:
							numDesc+=1
							child=c
					if numDesc==0:
						node.subtree=None
					elif numDesc==1:
						node.subtree=child.subtree
						node.subtree.dist+=node.dist
					else:
						node.subtree=Tree(dist=node.dist)
						for c in node.children:
							if c.subtree:
								node.subtree.add_child(c.subtree)
								c.subtree.up=node.subtree
								#node.subtree.children[-1].dist=c.subtree.dist

			#print("Distribution of number of descendants:")
			#print(numDescList)
			#print("num nodes explored: "+str(numNode))
			#print("root "+str(numDesc))
			#print("Replicate "+str(repeat))
			#print("num samples "+str(subsampleTreeInference))
			# if len(newSamples)<15:
			# 	print(newSamples)
			# 	print(len(phylo.subtree.children))
			# 	nodeList=postorderList(phylo.subtree)
			# 	print(nodeList[0].name)
			# 	print(nodeList[1].name)
			# 	print(len(nodeList))
			# 	print(pathToRepeat+"_realized.nw")
			newickString=createNewick(phylo.subtree)
			#if len(newSamples)<15:
			#	print(newickString)
			phylo=None
			#del newSamples
			del data
			phyloFile=open(pathToRepeat+"_realized.nw","w")
			phyloFile.write(newickString+"\n")
			phyloFile.close()
			print("Subtree extracted and written to file "+pathToRepeat+"_realized.nw")













if args.createBashScript:
	if args.runExtraAnalysesDivergence:
		ref = collectReference(args.inputSimulationReference)
		scenario=args.scenario
		subsampleTreeInference=args.subSampleNum
		divergences=["01X","02X","05X","1X","2X","5X","10X","20X","50X","100X","200X","500X","1000X"]
		folder="divergence"
		seed=args.repeat
		repeat=seed

		phylos=[args.pathToSimulationFolder+"phastSim_genomes_"+divergences[scenario]+"_realizedTree.tree"]
		pathToRepeat=args.pathToSimulationFolder+folder+'/'+divergences[scenario]+"/repeat"+str(repeat)+"_"+divergences[scenario]+"_"+str(subsampleTreeInference)+"samples"
		
		if args.numSamples==0:
			numSamples=[2000, 5000]
			numSamplesIQtreeLK=[2000, 5000]
		else:
			numSamples=[args.numSamples]
			numSamplesIQtreeLK=[args.numSamples]

		#MAPLEversions=["0.1.4","0.1.5","0.1.6","0.1.7","0.1.9"]
		MAPLEversions=["0.2.0"]
		#folders=["realDataSubsamples","simulationsSubsamples","simulations4catSubsamples","simulationsAlphaSubsamples","simulationsNsSubsamples"]

		MAPLEoptions=[""]
		MAPLEoptionsNames=[""]
		testOptions=False
		
		parameterValues=[""]
		parameterNames=[""]
		testParameters=False

		#creating bash script to run this python script in parallel on the cluster to create subsampled data
		file=open(args.pathToSimulationFolder+"createSubsampleInputFiles_divergence.sh","w")
		for j in numSamples:
			file.write("for i in $(seq 1 10)\n"+"do \n\t")
			for scenario in range(len(divergences)):
				folderNameSimu=args.pathToSimulationFolder+folder+'/'+divergences[scenario]+"/"
				file.write("bsub -M "+str(int(25000))+" -o "+folderNameSimu+"repl\"$i\"_"+str(j)+"samples_fileCreation_console_output.txt -e "+folderNameSimu+"repl\"$i\"_"+str(j)+"samples_fileCreation_console_error.txt"
				+" "+pypy3Path+" "+folderCluster+"/fastLK/code/MAPLE_benchmarking.py --createInputData --runExtraAnalysesDivergence --scenario "+str(scenario)+" --subSampleNum "+str(j)+" --repeat \"$i\" \n\n\t")
			file.write("done\n\n")
		file.close()
		print("Created bash script "+args.pathToSimulationFolder+"createSubsampleInputFiles_divergence.sh")

		version=MAPLEversions[0]
		versionForRF="0.2.0"
		versionForLK="0.2.0"

		#creating bash scripts for running all methods
		#foldersForTreeFile=["","simulationsSubsamples","simulations4catSubsamples","simulationsAlphaSubsamples","simulationsSubsamples"]
		methods=["UShER","matOptimize","IQtree","FastTree","RAxML-NG",'Maple']
		methodOptions=[[""],[""],[" -fast"],[" -fastest"],[" --blmin 0.000005 --tree pars\{1\}"],[""]]
		methodOptionsNames=[[""],[""],["_fast"],["_fastest"],["_fast"],[""]]
		#maxNumSamples=[[200000],[200000],[5000,10000,20000],[20000,20000,20000,20000],[5000,5000],[5000,5000],[500000,500000,500000,500000,500000]]
		fileMatConv=open(args.pathToSimulationFolder+"submitMatOptimizeConversion_divergence.sh","w")
		for method in range(len(methods)):
			methodName=methods[method]
			file=open(args.pathToSimulationFolder+"submit"+methodName+"_divergence.sh","w")
			file.write("bgadd /"+methodName+divergences[scenario]+str(numSamples[0])+"\n")
			fileIQtreeEvaluation=open(args.pathToSimulationFolder+"submitIQtreeLK_"+methodName+"_divergence.sh","w")
			fileMapleEvaluation=open(args.pathToSimulationFolder+"submitMapleLK_"+methodName+"_divergence.sh","w")
			fileRF=open(args.pathToSimulationFolder+"submitRF_"+methodName+"_divergence.sh","w")
			fileParsimony=open(args.pathToSimulationFolder+"submitParsimony_"+methodName+"_divergence.sh","w")
			fileParsimony.write("export PATH=\""+minicondaPath+"bin:$PATH\" ; source "+minicondaPath+"etc/profile.d/conda.sh ; conda activate usher-env \n")
			#if method==4:
			#	file.write("module load raxml-8.2.11-gcc-9.3.0-mjwrm3x \n")
			if method==4:
				file.write("module load raxml-ng-1.0.2-gcc-9.3.0-uicuzej \n")
			elif method==0 or method==1:
				file.write("export PATH=\""+minicondaPath+"bin:$PATH\" ; source "+minicondaPath+"etc/profile.d/conda.sh ; conda activate usher-env \n")
			if method==1:
				fileMatConv.write("bgadd /matConv"+str(numSamples[0])+"\n")
				fileMatConv.write("export PATH=\""+minicondaPath+"bin:$PATH\" ; source "+minicondaPath+"etc/profile.d/conda.sh ; conda activate usher-env \n")
			for option in range(len(methodOptions[method])):
				for j in numSamples:
					#if j<=maxNumSamples[method][option]:
						file.write("for i in $(seq 1 10)\n do\n\n")
						fileIQtreeEvaluation.write("for i in $(seq 1 10)\n do\n\n")
						fileMapleEvaluation.write("for i in $(seq 1 10)\n do\n\n")
						fileRF.write("for i in $(seq 1 10)\n do\n\n")
						fileParsimony.write("for i in $(seq 1 10)\n do\n\n")
						if method==1:
							fileMatConv.write("for i in $(seq 1 10)\n do\n\n")
						for scenario in range(len(divergences)):
							refFile=args.inputSimulationReference
							#folder=folders[scenario]
							#folderNameSimu=args.pathToSimulationFolder+folder+'/'+str(j)+"subsamples/output_repl\"$i\"/"
							folderNameSimu=args.pathToSimulationFolder+folder+'/'+divergences[scenario]+"/output_repl\"$i\"/"
							#pathToFile=args.pathToSimulationFolder+folder+'/'+str(j)+"subsamples/repeat\"$i\"_"+str(j)+"samples_"+folder
							#pathToTreeFile=args.pathToSimulationFolder+folder+'/'+str(j)+"subsamples/repeat\"$i\"_"+str(j)+"samples_"+folder
							pathToFile=args.pathToSimulationFolder+folder+'/'+divergences[scenario]+"/repeat\"$i\"_"+divergences[scenario]+"_"+str(j)+"samples"

							#remove existing result if already exists - this helps in case of re-running an analysis, and in case the new version fails, to wrongly use the old tree as estimate of the new version.
							fileNamePre=folderNameSimu+methodName+methodOptionsNames[method][option]+"_"+str(j)+"samples"
							rateVarNames=[""]
							outputNames=["_console_output.txt","_console_error.txt"]
							IQtreeLKoptionNames=["_IQtreeLK"]
							mapleLKoptionNames=["_MapleLK"]
							for rateVarName in rateVarNames:
								for outputName in outputNames:
									file.write("rm -f "+fileNamePre+rateVarName+outputName+" || true\n")
									for analysesName in IQtreeLKoptionNames:
										fileIQtreeEvaluation.write("rm -f "+fileNamePre+rateVarName+analysesName+outputName+" || true\n")
										fileIQtreeEvaluation.write("rm -f "+fileNamePre+rateVarName+analysesName+"_tree.tree || true\n")
									for analysesName in mapleLKoptionNames:
										fileMapleEvaluation.write("rm -f "+fileNamePre+rateVarName+analysesName+outputName+" || true\n")
									fileRF.write("rm -f "+fileNamePre+rateVarName+"_RF"+outputName+" || true\n")
									fileParsimony.write("rm -f "+fileNamePre+rateVarName+"_parsimony"+outputName+" || true\n")
							if method==0:
								treeFiles=[folderNameSimu+str(j)+"samples_UShER_output/final-tree.nh"]
							elif method==1:
								treeFiles=[folderNameSimu+str(j)+"samples_matOptimize_output.nh"]
							elif method==2:
								treeFiles=[fileNamePre+".treefile"]
							elif method==3:
								treeFiles=[fileNamePre+"_console_output.txt"]	
							elif method==4:
								treeFiles=[fileNamePre+"_output.raxml.bestTree"]
							elif method==5:
								treeFiles=[fileNamePre+"_tree.tree"]
							for treeFile in treeFiles:
								file.write("rm -f "+treeFile+" || true\n")
										
							#run method
							#UShER
							if method==0:
								file.write("bsub -g /"+methodName+divergences[scenario]+str(numSamples[0])+" -M "+str(int(10000+j*3))+" -o "+fileNamePre+"_console_output.txt -e "+fileNamePre+"_console_error.txt "
								+"usher -v "+pathToFile+".vcf -d "+folderNameSimu+str(j)+"samples_UShER_output/ -T 1 -t "+pathToFile+"_initialTree.tree \n")
							
							#matOptimize
							elif method==1:
								file.write("bsub -g /"+methodName+divergences[scenario]+str(numSamples[0])+" -M "+str(int(10000+j*3))+" -o "+fileNamePre+"_console_output.txt -e "+fileNamePre+"_console_error.txt "
								+"matOptimize -n -v "+pathToFile+".vcf -o "+folderNameSimu+str(j)+"samples_matOptimize_output.pb -T 1 -t "+folderNameSimu+str(j)+"samples_UShER_output/final-tree.nh \n")
								fileMatConv.write("bsub -g /matConv"+divergences[scenario]+str(numSamples[0])+" -M "+str(int(10000+j*3))+" -o "+fileNamePre+"_conv_console_output.txt -e "+fileNamePre+"_conv_console_error.txt "
								+"matUtils extract -i "+folderNameSimu+str(j)+"samples_matOptimize_output.pb -T 1 -d "+folderNameSimu+" -t "+str(j)+"samples_matOptimize_output.nh \n")
							
							#IQtree
							elif method==2:
								file.write("cd "+folderNameSimu+"\n\t"+"bsub -g /"+methodName+divergences[scenario]+str(numSamples[0])+" -M "+str(int(10000+j*3))+" -o "+fileNamePre+"_console_output.txt -e "+fileNamePre+"_console_error.txt "
								+""+iqTreePath+" -s "+pathToFile+".phy -st DNA -pre "+methodName+methodOptionsNames[method][option]+"_"+str(j)+"samples"+" -m GTR -quiet -redo -nt 1 "+methodOptions[method][option]+" \n")
								
							#FastTree
							elif method==3:
								file.write("cd "+folderNameSimu+"\n\t"+"bsub -g /"+methodName+divergences[scenario]+str(numSamples[0])+" -M "+str(int(10000+j*3))+" -o "+fileNamePre+"_console_output.txt -e "+fileNamePre+"_console_error.txt "
								+""+fastTreePath+" -quiet -nosupport -nt -gtr -nocat "+methodOptions[method][option]+" "+pathToFile+".fa \n")

							#RAxML-NG
							elif method==4:
								file.write("bsub -g /"+methodName+divergences[scenario]+str(numSamples[0])+" -M "+str(int(j*3+10000))+" -o "+fileNamePre+"_console_output.txt -e "+fileNamePre+"_console_error.txt "
								+"raxml-ng --search "+methodOptions[method][option]+" --msa "+pathToFile+".fa --msa-format FASTA --data-type DNA --redo --prefix "+fileNamePre+"_output --model GTR --threads 1 \n")
								
							#MAPLE
							elif method==5:
								file.write("bsub -g /"+methodName+divergences[scenario]+str(numSamples[0])+" -M "+str(int(25000))+" -o "+fileNamePre+"_console_output.txt -e "+fileNamePre+"_console_error.txt "
								+""+pypy3Path+" "+folderCluster+"/fastLK/code/MAPLEv"+version+".py --reference "+refFile+methodOptions[method][option]+" --input "
								+pathToFile+".txt --overwrite --output "+fileNamePre+"\n")

							#run IQtree LK evaluation for all methods
							fileIQtreeEvaluation.write("cd "+folderNameSimu+"\n\t"+"bsub -g /"+methodName+divergences[scenario]+str(numSamples[0])+"IQLK -M "+str(int(10000+j*3))+" -o "+fileNamePre+"_IQtreeLK_console_output.txt -e "+fileNamePre+"_IQtreeLK_console_error.txt "
							+""+iqTreePath+" -s "+pathToFile+".phy -st DNA -te "+treeFiles[0]+" -pre "+methodName+methodOptionsNames[method][option]+"_"+str(j)+"samples"+"_IQtreeLK -m GTR -quiet -nt 1 -keep-ident -redo -blmin 0.000000005 \n")
							
							#Maple LK evaluation for all methods
							mapleLKoptionNames=["_MapleLK"]
							mapleLKoptions=[""]
							bLenScaleOption=""
							if method<2:
								bLenScaleOption=" --normalizeInputBLen 0.000033 "
							for mapleLKoptionIndex in range(len(mapleLKoptionNames)):
								fileMapleEvaluation.write("bsub -g /"+methodName+divergences[scenario]+str(numSamples[0])+"MLK"+str(mapleLKoptionIndex)+" -M "+str(int(10000+j*3))+" -o "+fileNamePre+mapleLKoptionNames[mapleLKoptionIndex]+"_console_output.txt -e "+fileNamePre+mapleLKoptionNames[mapleLKoptionIndex]+"_console_error.txt "
								+""+pypy3Path+" "+folderCluster+"/fastLK/code/MAPLEv"+versionForLK+".py --reference "+refFile+mapleLKoptions[mapleLKoptionIndex]+bLenScaleOption+" --inputTree "+treeFiles[0]+" --input "
								+pathToFile+".txt --calculateLKfinalTree --numTopologyImprovements 0 --noFastTopologyInitialSearch  --overwrite --output "+fileNamePre+mapleLKoptionNames[mapleLKoptionIndex]+"\n")
								
							#RF calculations
							fileRF.write("bsub -g /"+methodName+divergences[scenario]+str(numSamples[0])+"RF -M "+str(int(4000+j/20))+" -o "+fileNamePre+"_RF_console_output.txt -e "+fileNamePre+"_RF_console_error.txt "
							+""+pypy3Path+" "+folderCluster+"/fastLK/code/MAPLEv"+versionForRF+".py --inputRFtrees "+treeFiles[0]+" --inputTree "+pathToFile+"_realized.nw "
							+" --overwrite --output "+fileNamePre+"_RF\n")

							#parsimony calculations
							fileParsimony.write("bsub -g /"+methodName+divergences[scenario]+str(numSamples[0])+"par -M "+str(int(1000+j*3))+" -o "+fileNamePre+"_parsimony_console_output.txt -e "+fileNamePre+"_parsimony_console_error.txt "
							+"usher -v "+pathToFile+".vcf -t "+treeFiles[0]+" -o "+fileNamePre+"_parsimony.mat \n")

						file.write("done\n\n")
						fileIQtreeEvaluation.write("done\n\n")
						fileMapleEvaluation.write("done\n\n")
						fileRF.write("done\n\n")
						fileParsimony.write("done\n\n")
						if method==1:
							fileMatConv.write("done\n\n")
			file.close()
			fileIQtreeEvaluation.close()
			fileMapleEvaluation.close()
			fileRF.close()
			fileParsimony.close()
			print("Created "+methodName+" bash script "+args.pathToSimulationFolder+"submit"+methodName+"_divergence.sh")
			print("Created "+methodName+" bash script "+args.pathToSimulationFolder+"submitIQtreeLK_"+methodName+"_divergence.sh")
			print("Created "+methodName+" bash script "+args.pathToSimulationFolder+"submitMapleLK_"+methodName+"_divergence.sh")
			print("Created "+methodName+" bash script "+args.pathToSimulationFolder+"submitRF_"+methodName+"_divergence.sh")
			print("Created "+methodName+" bash script "+args.pathToSimulationFolder+"submitParsimony_"+methodName+"_divergence.sh")
		fileMatConv.close()
		print("Created matConvert bash script "+args.pathToSimulationFolder+"submitMatOptimizeConversion_divergence.sh")



	else:
		numSamplesIQtreeLK=[1000, 2000, 5000, 10000, 20000]
		if args.numSamples==0:
			numSamples=[1000, 2000, 5000, 10000, 20000, 50000, 100000, 200000, 500000]
			numSamplesIQtreeLK=[1000, 2000, 5000, 10000, 20000]
		else:
			numSamples=[args.numSamples]

		#MAPLEversions=["0.1.4","0.1.5","0.1.6","0.1.7","0.1.9"]
		MAPLEversions=["0.1.9"]
		folders=["realDataSubsamples","simulationsSubsamples","simulations4catSubsamples","simulationsAlphaSubsamples","simulationsNsSubsamples"]

		MAPLEoptions=[""," --fast"]
		MAPLEoptionsNames=["","_fast"]
		testOptions=False
		if testOptions:
			#MAPLEoptions=[""," --fast"," --rateVariation"," --fast --rateVariation"," --model UNREST"," --fast --model UNREST"," --model UNREST --rateVariation"," --fast --model UNREST --rateVariation"]
			#MAPLEoptionsNames=["","_fast","_rateVar","_fast_rateVar","_unrest","_fast_unrest","_unrest_rateVar","_fast_unrest_rateVar"]
			MAPLEoptions=[""," --fast"," --rateVariation"," --model UNREST"," --model UNREST --rateVariation"]
			MAPLEoptionsNames=["","_fast","_rateVar","_unrest","_unrest_rateVar"]
		
		parameterValues=[""]
		parameterNames=[""]
		testParameters=False
		if testParameters:
			#parameterValues=[" --minBLenSensitivity 0.1"," --minBLenSensitivity 0.01"," --minBLenSensitivity 0.001"," --minBLenSensitivity 0.0001"," --minBLenSensitivity 0.00001"," --minBLenSensitivity 0.000001"," --factorOptimizePlacementLKvsSearchLK 0.004"," --factorOptimizePlacementLKvsSearchLK 0.01"," --factorOptimizePlacementLKvsSearchLK 0.04"," --factorOptimizePlacementLKvsSearchLK 0.1"," --factorOptimizePlacementLKvsSearchLK 0.2"," --factorOptimizePlacementLKvsSearchLK 0.4"]
			parameterValues=[" --minBLenSensitivity 0.01"," --minBLenSensitivity 0.001"," --minBLenSensitivity 0.0001"," --thresholdLogLKoptimization 18.0"," --thresholdLogLKoptimization 12.0"," --thresholdLogLKoptimization 6.0"," --thresholdLogLKoptimization 1.0"," --thresholdLogLKoptimization 0.1"," --thresholdLogLKoptimizationTopology 18.0"," --thresholdLogLKoptimizationTopology 12.0"," --thresholdLogLKoptimizationTopology 6.0"," --thresholdLogLKoptimizationTopology 1.0"," --thresholdLogLKoptimizationTopology 0.1"]
			parameterNames=["_minBLen01","_minBLen001","_minBLen0001","_LK18","_LK12","_LK6","_LK1","_LK01","_LKtopology18","_LKtopology12","_LKtopology6","_LKtopology1","_LKtopology01"]
			#parameterNames=["_minBLen1","_minBLen01","_minBLen001","_minBLen0001","_minBLen00001","_minBLen000001","_factor004","_factor01","_factor04","_factor1","_factor2","_factor4"]
			folders=["realDataSubsamples"]

		#creating bash script to run this python script in parallel on the cluster to create subsampled data
		file=open(args.pathToSimulationFolder+"createSubsampleInputFiles.sh","w")
		for j in numSamples:
			file.write("for i in $(seq 1 10)\n"+"do \n\t")
			for scenario in range(len(folders)):
				folderNameSimu=args.pathToSimulationFolder+folders[scenario]+'/'+str(j)+"subsamples/"
				file.write("bsub -M "+str(int(15000+j/20))+" -o "+folderNameSimu+"repl\"$i\"_fileCreation_console_output.txt -e "+folderNameSimu+"repl\"$i\"_fileCreation_console_error.txt"
				+" "+pypy3Path+" "+folderCluster+"/fastLK/code/MAPLE_benchmarking.py --createInputData --scenario "+str(scenario)+" --subSampleNum "+str(j)+" --repeat \"$i\" \n\n\t")
			file.write("done\n\n")
		file.close()
		print("Created bash script "+args.pathToSimulationFolder+"createSubsampleInputFiles.sh")

		version=MAPLEversions[0]
		versionForRF="0.1.9"
		versionForLK="0.1.9"
		if testOptions or testParameters or len(MAPLEversions)>1:
			#creating bash scripts for MAPLE
			#foldersForTreeFile=["","simulationsSubsamples","simulations4catSubsamples","simulationsAlphaSubsamples","simulationsSubsamples"]
			for versNum in range(len(MAPLEversions)):
				version=MAPLEversions[versNum]
				file=open(args.pathToSimulationFolder+"submitMAPLE"+version+".sh","w")
				for j in numSamples:
					file.write("for i in $(seq 1 10)\n do\n\n")
					for scenario in range(len(folders)):
						if scenario==0:
							refFile=args.inputRealDataReference
						else:
							refFile=args.inputSimulationReference
						folder=folders[scenario]
						folderNameSimu=args.pathToSimulationFolder+folder+'/'+str(j)+"subsamples/output_repl\"$i\"/"
						pathToFile=args.pathToSimulationFolder+folder+'/'+str(j)+"subsamples/repeat\"$i\"_"+str(j)+"samples_"+folder
						for option in range(len(MAPLEoptions)):
							for paramNum in range(len(parameterValues)):
								#remove existing result if already exists - this helps in case of re-running an analysis, and in case the new version fails, to wrongly use the old tree as estimate of the new version.
								file.write("rm -f "+folderNameSimu+"MAPLE"+version+MAPLEoptionsNames[option]+parameterNames[paramNum]+"_tree.tree || true\n")
								file.write("rm -f "+folderNameSimu+"MAPLE"+version+MAPLEoptionsNames[option]+parameterNames[paramNum]+"_console_output.txt || true\n")
								file.write("rm -f "+folderNameSimu+"MAPLE"+version+MAPLEoptionsNames[option]+parameterNames[paramNum]+"_console_error.txt || true\n")
								#run MAPLE
								file.write("bsub -M "+str(int(400+j/20))+" -o "+folderNameSimu+"MAPLE"+version+MAPLEoptionsNames[option]+parameterNames[paramNum]+"_console_output.txt -e "+folderNameSimu+"MAPLE"+version+MAPLEoptionsNames[option]+parameterNames[paramNum]+"_console_error.txt "
								+""+pypy3Path+" "+folderCluster+"/fastLK/code/MAPLEv"+version+".py --reference "+refFile+MAPLEoptions[option]+parameterValues[paramNum]+" --input "
								+pathToFile+".txt --overwrite --output "+folderNameSimu+"MAPLE"+version+MAPLEoptionsNames[option]+parameterNames[paramNum]+"\n")
					file.write("done\n\n")
				file.close()
				print("Created MAPLE bash script "+args.pathToSimulationFolder+"submitMAPLE"+version+".sh")

				#create bash script for running robinson-foulds distance estimations
				file=open(args.pathToSimulationFolder+"submitRF_MAPLE"+version+".sh","w")
				for j in numSamples:
					file.write("for i in $(seq 1 10)\n do\n\n")
					for scenario in range(len(folders)):
						if scenario>0:
							refFile=args.inputSimulationReference
							folder=folders[scenario]
							folderNameSimu=args.pathToSimulationFolder+folder+'/'+str(j)+"subsamples/output_repl\"$i\"/"
							pathToFile=args.pathToSimulationFolder+folder+'/'+str(j)+"subsamples/repeat\"$i\"_"+str(j)+"samples_"+folder
							pathToTreeFile=args.pathToSimulationFolder+folder+'/'+str(j)+"subsamples/repeat\"$i\"_"+str(j)+"samples_"+folder
							#pathToTreeFile=args.pathToSimulationFolder+folder+'/'+str(j)+"subsamples/repeat\"$i\"_"+str(j)+"samples_"+foldersForTreeFile[scenario]
							#"+folderCluster+"/fastLK/simulationsNew/simulationsNsSubsamples/1000subsamples/repeat6_1000samples_simulationsSubsamples_realized.nw
							for option in range(len(MAPLEoptions)):
								for paramNum in range(len(parameterValues)):
									#remove existing result if already exists - this helps in case of re-running an analysis, and in case the new version fails, to wrongly use the old tree as estimate of the new version.
									file.write("rm -f "+folderNameSimu+"MAPLE"+version+MAPLEoptionsNames[option]+parameterNames[paramNum]+"_RFdistances.txt || true\n")
									file.write("rm -f "+folderNameSimu+"MAPLE"+version+MAPLEoptionsNames[option]+parameterNames[paramNum]+"_RF_console_output.txt || true\n")
									file.write("rm -f "+folderNameSimu+"MAPLE"+version+MAPLEoptionsNames[option]+parameterNames[paramNum]+"_RF_console_error.txt || true\n")
									#run MAPLE RF alculation
									file.write("bsub -M "+str(int(400+j/20))+" -o "+folderNameSimu+"MAPLE"+version+MAPLEoptionsNames[option]+parameterNames[paramNum]+"_RF_console_output.txt -e "+folderNameSimu+"MAPLE"+version+MAPLEoptionsNames[option]+parameterNames[paramNum]+"_RF_console_error.txt "
									+""+pypy3Path+" "+folderCluster+"/fastLK/code/MAPLEv"+versionForRF+".py --inputRFtrees "+folderNameSimu+"MAPLE"+version+MAPLEoptionsNames[option]+parameterNames[paramNum]+"_tree.tree --inputTree "+pathToTreeFile+"_realized.nw "
									+" --overwrite --output "+folderNameSimu+"MAPLE"+version+MAPLEoptionsNames[option]+parameterNames[paramNum]+"\n")
					file.write("done\n\n")
				file.close()
				print("Created RF bash script "+args.pathToSimulationFolder+"submitRF_MAPLE"+version+".sh")

				#create bash script for running likelihood evaluations (in MAPLE)
				furtherLKoptions=[""," --model UNREST"," --rateVariation"," --model UNREST --rateVariation"]
				furtherLKoptionsNames=["","_unrest","_rateVar","_unrest_rateVar"]
				file=open(args.pathToSimulationFolder+"submitLK_MAPLE"+version+".sh","w")
				for j in numSamples:
					file.write("for i in $(seq 1 10)\n do\n\n")
					for scenario in range(len(folders)):
						if scenario==0:
							refFile=args.inputRealDataReference
						else:
							refFile=args.inputSimulationReference
						folder=folders[scenario]
						folderNameSimu=args.pathToSimulationFolder+folder+'/'+str(j)+"subsamples/output_repl\"$i\"/"
						pathToFile=args.pathToSimulationFolder+folder+'/'+str(j)+"subsamples/repeat\"$i\"_"+str(j)+"samples_"+folder

						for option in range(len(MAPLEoptions)):
							for paramNum in range(len(parameterValues)):
								for furtherOptionNum in range(len(furtherLKoptions)):
									#remove existing result if already exists - this helps in case of re-running an analysis, and in case the new version fails, to wrongly use the old tree as estimate of the new version.
									file.write("rm -f "+folderNameSimu+"MAPLE"+version+MAPLEoptionsNames[option]+parameterNames[paramNum]+"_LK"+furtherLKoptionsNames[furtherOptionNum]+"_LK.txt || true\n")
									file.write("rm -f "+folderNameSimu+"MAPLE"+version+MAPLEoptionsNames[option]+parameterNames[paramNum]+"_LK"+furtherLKoptionsNames[furtherOptionNum]+"_console_output.txt || true\n")
									file.write("rm -f "+folderNameSimu+"MAPLE"+version+MAPLEoptionsNames[option]+parameterNames[paramNum]+"_LK"+furtherLKoptionsNames[furtherOptionNum]+"_console_error.txt || true\n")
									#run MAPLE LK calculation
									file.write("bsub -M "+str(int(400+j/20))+" -o "+folderNameSimu+"MAPLE"+version+MAPLEoptionsNames[option]+parameterNames[paramNum]+"_LK"+furtherLKoptionsNames[furtherOptionNum]+"_console_output.txt -e "+folderNameSimu+"MAPLE"+version+MAPLEoptionsNames[option]+parameterNames[paramNum]+"_LK"+furtherLKoptionsNames[furtherOptionNum]+"_console_error.txt "
									+""+pypy3Path+" "+folderCluster+"/fastLK/code/MAPLEv"+versionForLK+".py --reference "+refFile+furtherLKoptions[furtherOptionNum]+" --inputTree "+folderNameSimu+"MAPLE"+version+MAPLEoptionsNames[option]+parameterNames[paramNum]+"_tree.tree --input "
									+pathToFile+".txt --calculateLKfinalTree --numTopologyImprovements 0 --noFastTopologyInitialSearch  --overwrite --output "+folderNameSimu+"MAPLE"+version+MAPLEoptionsNames[option]+parameterNames[paramNum]+"_LK"+furtherLKoptionsNames[furtherOptionNum]+"\n")

					file.write("done\n\n")
				file.close()
				print("Created LK bash script "+args.pathToSimulationFolder+"submitLK_MAPLE"+version+".sh")

				#create bash script for running likelihood evaluations (in IQtree)
				#folders=["realDataSubsamples","simulationsSubsamples","simulations4catSubsamples","simulationsAlphaSubsamples","simulationsNsSubsamples"]
				file=open(args.pathToSimulationFolder+"submitIQtreeLK_MAPLE"+version+".sh","w")
				for j in numSamplesIQtreeLK:
					file.write("for i in $(seq 1 10)\n do\n\n")
					for scenario in range(len(folders)):
						folder=folders[scenario]
						folderNameSimu=args.pathToSimulationFolder+folder+'/'+str(j)+"subsamples/output_repl\"$i\"/"
						pathToFile=args.pathToSimulationFolder+folder+'/'+str(j)+"subsamples/repeat\"$i\"_"+str(j)+"samples_"+folder

						for option in range(len(MAPLEoptions)):
							file.write("rm -f "+folderNameSimu+"MAPLE"+version+MAPLEoptionsNames[option]+"_IQtreeLK_console_output.txt || true\n")
							file.write("rm -f "+folderNameSimu+"MAPLE"+version+MAPLEoptionsNames[option]+"_IQtreeLK_console_error.txt || true\n")
							file.write("cd "+folderNameSimu+"\n\t"+"bsub -M "+str(int(500+j*1.7))+" -o "+folderNameSimu+"MAPLE"+version+MAPLEoptionsNames[option]+"_IQtreeLK_console_output.txt -e "+folderNameSimu+"MAPLE"+version+MAPLEoptionsNames[option]+"_IQtreeLK_console_error.txt "
							+""+iqTreePath+" -s "+pathToFile+".phy -st DNA -te "+folderNameSimu+"MAPLE"+version+MAPLEoptionsNames[option]+"_tree.tree -pre MAPLE"+version+MAPLEoptionsNames[option]+"_IQtreeLK -m GTR -quiet -nt 1 -keep-ident -redo -blmin 0.000000005 \n")
							if scenario==0 or scenario==2 or scenario==3:
								file.write("rm -f "+folderNameSimu+"MAPLE"+version+MAPLEoptionsNames[option]+"_IQtreeLK_rateVar_console_output.txt || true\n")
								file.write("rm -f "+folderNameSimu+"MAPLE"+version+MAPLEoptionsNames[option]+"_IQtreeLK_rateVar_console_error.txt || true\n")
								file.write("cd "+folderNameSimu+"\n\t"+"bsub -M "+str(int(500+j*3.0))+" -o "+folderNameSimu+"MAPLE"+version+MAPLEoptionsNames[option]+"_IQtreeLK_rateVar_console_output.txt -e "+folderNameSimu+"MAPLE"+version+MAPLEoptionsNames[option]+"_IQtreeLK_rateVar_console_error.txt "
								+""+iqTreePath+" -s "+pathToFile+".phy -st DNA -te "+folderNameSimu+"MAPLE"+version+MAPLEoptionsNames[option]+"_tree.tree -pre MAPLE"+version+MAPLEoptionsNames[option]+"_IQtreeLK_rateVar -m GTR+G -quiet -nt 1 -keep-ident -redo -blmin 0.000000005 \n")

					file.write("done\n\n")
				file.close()
				print("Created IQtreeLK bash script "+args.pathToSimulationFolder+"submitIQtreeLK_MAPLE"+version+".sh")






		#creating bash scripts for running all methods
		#foldersForTreeFile=["","simulationsSubsamples","simulations4catSubsamples","simulationsAlphaSubsamples","simulationsSubsamples"]
		methods=["UShER","matOptimize","IQtree","FastTree","RAxML","RAxML-NG",'Maple']
		methodOptions=[[""],[""],[" -pers 0.1 -nstop 500 -blmin 0.000000005"," -blmin 0.000000005"," -fast"],[" -spr 4 -mlacc 2 -slownni"," -spr 4",""," -fastest"],[""," -D"],[" --blmin 0.000000005 --tree pars\{3\}"," --blmin 0.000005 --tree pars\{1\}"],[""," --fast"," --rateVariation"," --model UNREST"," --model UNREST --rateVariation"]]
		methodOptionsNames=[[""],[""],["_slow","_medium","_fast"],["_slow","_medium","_fast","_fastest"],["_slow","_fast"],["_slow","_fast"],["","_fast","_rateVar","_unrest","_unrest_rateVar"]]
		maxNumSamples=[[200000],[200000],[5000,10000,20000],[20000,20000,20000,20000],[5000,5000],[5000,5000],[500000,500000,500000,500000,500000]]
		fileMatConv=open(args.pathToSimulationFolder+"submitMatOptimizeConversion.sh","w")
		for method in range(len(methods)):
			methodName=methods[method]
			file=open(args.pathToSimulationFolder+"submit"+methodName+".sh","w")
			file.write("bgadd /"+methodName+str(numSamples[0])+"\n")
			fileIQtreeEvaluation=open(args.pathToSimulationFolder+"submitIQtreeLK_"+methodName+".sh","w")
			fileMapleEvaluation=open(args.pathToSimulationFolder+"submitMapleLK_"+methodName+".sh","w")
			fileRF=open(args.pathToSimulationFolder+"submitRF_"+methodName+".sh","w")
			fileParsimony=open(args.pathToSimulationFolder+"submitParsimony_"+methodName+".sh","w")
			fileParsimony.write("export PATH=\""+minicondaPath+"bin:$PATH\" ; source "+minicondaPath+"etc/profile.d/conda.sh ; conda activate usher-env \n")
			if method==4:
				file.write("module load raxml-8.2.11-gcc-9.3.0-mjwrm3x \n")
			elif method==5:
				file.write("module load raxml-ng-1.0.2-gcc-9.3.0-uicuzej \n")
			elif method==0 or method==1:
				file.write("export PATH=\""+minicondaPath+"bin:$PATH\" ; source "+minicondaPath+"etc/profile.d/conda.sh ; conda activate usher-env \n")
			if method==1:
				fileMatConv.write("bgadd /matConv"+str(numSamples[0])+"\n")
				fileMatConv.write("export PATH=\""+minicondaPath+"bin:$PATH\" ; source "+minicondaPath+"etc/profile.d/conda.sh ; conda activate usher-env \n")
			for option in range(len(methodOptions[method])):
				for j in numSamples:
					if j<=maxNumSamples[method][option]:
						file.write("for i in $(seq 1 10)\n do\n\n")
						fileIQtreeEvaluation.write("for i in $(seq 1 10)\n do\n\n")
						fileMapleEvaluation.write("for i in $(seq 1 10)\n do\n\n")
						fileRF.write("for i in $(seq 1 10)\n do\n\n")
						fileParsimony.write("for i in $(seq 1 10)\n do\n\n")
						if method==1:
							fileMatConv.write("for i in $(seq 1 10)\n do\n\n")
						for scenario in range(len(folders)):
							if scenario==0:
								refFile=args.inputRealDataReference
							else:
								refFile=args.inputSimulationReference
							folder=folders[scenario]
							folderNameSimu=args.pathToSimulationFolder+folder+'/'+str(j)+"subsamples/output_repl\"$i\"/"
							pathToFile=args.pathToSimulationFolder+folder+'/'+str(j)+"subsamples/repeat\"$i\"_"+str(j)+"samples_"+folder
							#pathToTreeFile=args.pathToSimulationFolder+folder+'/'+str(j)+"subsamples/repeat\"$i\"_"+str(j)+"samples_"+folder

							#remove existing result if already exists - this helps in case of re-running an analysis, and in case the new version fails, to wrongly use the old tree as estimate of the new version.
							fileNamePre=folderNameSimu+methodName+methodOptionsNames[method][option]
							rateVarNames=["","_rateVar"]
							outputNames=["_console_output.txt","_console_error.txt"]
							IQtreeLKoptionNames=["_IQtreeLK","_IQtreeLK_rateVar"]
							mapleLKoptionNames=["_MapleLK","_MapleLK_rateVar","_MapleLK_unrest","_MapleLK_unrest_rateVar"]
							for rateVarName in rateVarNames:
								for outputName in outputNames:
									file.write("rm -f "+fileNamePre+rateVarName+outputName+" || true\n")
									for analysesName in IQtreeLKoptionNames:
										fileIQtreeEvaluation.write("rm -f "+fileNamePre+rateVarName+analysesName+outputName+" || true\n")
										fileIQtreeEvaluation.write("rm -f "+fileNamePre+rateVarName+analysesName+"_tree.tree || true\n")
									for analysesName in mapleLKoptionNames:
										fileMapleEvaluation.write("rm -f "+fileNamePre+rateVarName+analysesName+outputName+" || true\n")
									fileRF.write("rm -f "+fileNamePre+rateVarName+"_RF"+outputName+" || true\n")
									fileParsimony.write("rm -f "+fileNamePre+rateVarName+"_parsimony"+outputName+" || true\n")
							if method==0:
								treeFiles=[folderNameSimu+"UShER_output/final-tree.nh"]
							elif method==1:
								treeFiles=[folderNameSimu+"matOptimize_output.nh"]
							elif method==2:
								treeFiles=[fileNamePre+".treefile"]
								if scenario==0 or scenario==2 or scenario==3:
									treeFiles.append(fileNamePre+"_rateVar.treefile")
							elif method==3:
								treeFiles=[fileNamePre+"_console_output.txt"]
								if scenario==0 or scenario==2 or scenario==3:
									treeFiles.append(fileNamePre+"_rateVar_console_output.txt")
							elif method==4:
								treeFiles=[folderNameSimu+"RAxML_result."+methodName+methodOptionsNames[method][option]+"_output"]
								if scenario==0 or scenario==2 or scenario==3:
									treeFiles.append(folderNameSimu+"RAxML_result."+methodName+methodOptionsNames[method][option]+"_rateVar_output")
									file.write("rm -f "+folderNameSimu+"RAxML_log."+methodName+methodOptionsNames[method][option]+"_rateVar_output || true\n")
									file.write("rm -f "+folderNameSimu+"RAxML_info."+methodName+methodOptionsNames[method][option]+"_rateVar_output || true\n")
									file.write("rm -f "+folderNameSimu+"RAxML_bestTree."+methodName+methodOptionsNames[method][option]+"_rateVar_output || true\n")
									file.write("rm -f "+folderNameSimu+"RAxML_parsimonyTree."+methodName+methodOptionsNames[method][option]+"_rateVar_output || true\n")
								file.write("rm -f "+folderNameSimu+"RAxML_log."+methodName+methodOptionsNames[method][option]+"_output || true\n")
								file.write("rm -f "+folderNameSimu+"RAxML_info."+methodName+methodOptionsNames[method][option]+"_output || true\n")
								file.write("rm -f "+folderNameSimu+"RAxML_bestTree."+methodName+methodOptionsNames[method][option]+"_output || true\n")
								file.write("rm -f "+folderNameSimu+"RAxML_parsimonyTree."+methodName+methodOptionsNames[method][option]+"_output || true\n")	
							elif method==5:
								treeFiles=[fileNamePre+"_output.raxml.bestTree"]
								if scenario==0 or scenario==2 or scenario==3:
									treeFiles.append(fileNamePre+"_rateVar_output.raxml.bestTree")
							elif method==6:
								treeFiles=[fileNamePre+"_tree.tree"]
							for treeFile in treeFiles:
								file.write("rm -f "+treeFile+" || true\n")
										
							#run method
							#UShER
							if method==0:
								file.write("bsub -g /"+methodName+str(numSamples[0])+" -M "+str(int(400+j/50))+" -o "+fileNamePre+"_console_output.txt -e "+fileNamePre+"_console_error.txt "
								+"usher -v "+pathToFile+".vcf -d "+folderNameSimu+"UShER_output/ -T 1 -t "+pathToFile+"_initialTree.tree \n")
							
							#matOptimize
							elif method==1:
								file.write("bsub -g /"+methodName+str(numSamples[0])+" -M "+str(int(400+j/5))+" -o "+fileNamePre+"_console_output.txt -e "+fileNamePre+"_console_error.txt "
								+"matOptimize -n -v "+pathToFile+".vcf -o "+folderNameSimu+"matOptimize_output.pb -T 1 -t "+folderNameSimu+"UShER_output/final-tree.nh \n")
								fileMatConv.write("bsub -g /matConv"+str(numSamples[0])+" -M "+str(int(400+j/5))+" -o "+fileNamePre+"_conv_console_output.txt -e "+fileNamePre+"_conv_console_error.txt "
								+"matUtils extract -i "+folderNameSimu+"matOptimize_output.pb -T 1 -d "+folderNameSimu+" -t matOptimize_output.nh \n")
							
							#IQtree
							elif method==2:
								file.write("cd "+folderNameSimu+"\n\t"+"bsub -g /"+methodName+str(numSamples[0])+" -M "+str(int(j*1.7))+" -o "+fileNamePre+"_console_output.txt -e "+fileNamePre+"_console_error.txt "
								+""+iqTreePath+" -s "+pathToFile+".phy -st DNA -pre "+methodName+methodOptionsNames[method][option]+" -m GTR -quiet -redo -nt 1 "+methodOptions[method][option]+" \n")
								if scenario==0 or scenario==2 or scenario==3:
									file.write("cd "+folderNameSimu+"\n\t"+"bsub -g /"+methodName+str(numSamples[0])+" -M "+str(int(j*3.0))+" -o "+fileNamePre+"_rateVar_console_output.txt -e "+fileNamePre+"_rateVar_console_error.txt "
									+""+iqTreePath+" -s "+pathToFile+".phy -st DNA -pre "+methodName+methodOptionsNames[method][option]+"_rateVar -m GTR+G -quiet -redo -nt 1 "+methodOptions[method][option]+" \n")
							
							#FastTree
							elif method==3:
								file.write("cd "+folderNameSimu+"\n\t"+"bsub -g /"+methodName+str(numSamples[0])+" -M "+str(int(j*1.7+500))+" -o "+fileNamePre+"_console_output.txt -e "+fileNamePre+"_console_error.txt "
								+""+fastTreePath+" -quiet -nosupport -nt -gtr -nocat "+methodOptions[method][option]+" "+pathToFile+".fa \n")
								if scenario==0 or scenario==3:
									file.write("cd "+folderNameSimu+"\n\t"+"bsub -g /"+methodName+str(numSamples[0])+" -M "+str(int(j*1.7+500))+" -o "+fileNamePre+"_rateVar_console_output.txt -e "+fileNamePre+"_rateVar_console_error.txt "
									+""+fastTreePath+" -quiet -nosupport -nt -gtr "+methodOptions[method][option]+" "+pathToFile+".fa \n")
								if scenario==2:
									file.write("cd "+folderNameSimu+"\n\t"+"bsub -g /"+methodName+str(numSamples[0])+" -M "+str(int(j*1.7+500))+" -o "+fileNamePre+"_rateVar_console_output.txt -e "+fileNamePre+"_rateVar_console_error.txt "
									+""+fastTreePath+" -quiet -nosupport -nt -gtr -cat 4 "+methodOptions[method][option]+" "+pathToFile+".fa  \n")
							
							#RAxML
							elif method==4:
								file.write("bsub -g /"+methodName+str(numSamples[0])+" -M "+str(int(j*3+500))+" -o "+fileNamePre+"_console_output.txt -e "+fileNamePre+"_console_error.txt "
								+"raxmlHPC "+methodOptions[method][option]+" -s "+pathToFile+".phy -p 1 -c 1 -m GTRCAT -V -n "+methodName+methodOptionsNames[method][option]+"_output -w "+folderNameSimu+" \n")
								if scenario==0 or scenario==2 or scenario==3:
									file.write("bsub -g /"+methodName+str(numSamples[0])+" -M "+str(int(j*3+500))+" -o "+fileNamePre+"_rateVar_console_output.txt -e "+fileNamePre+"_rateVar_console_error.txt "
									+"raxmlHPC "+methodOptions[method][option]+" -s "+pathToFile+".phy -p 1 -c 4 -m GTRCAT -V -n "+methodName+methodOptionsNames[method][option]+"_rateVar_output -w "+folderNameSimu+" \n")
							
							#RAxML-NG
							elif method==5:
								file.write("bsub -g /"+methodName+str(numSamples[0])+" -M "+str(int(j+500))+" -o "+fileNamePre+"_console_output.txt -e "+fileNamePre+"_console_error.txt "
								+"raxml-ng --search "+methodOptions[method][option]+" --msa "+pathToFile+".fa --msa-format FASTA --data-type DNA --redo --prefix "+fileNamePre+"_output --model GTR --threads 1 \n")
								if scenario==0 or scenario==2 or scenario==3:
									file.write("bsub -g /"+methodName+str(numSamples[0])+" -M "+str(int(j+500))+" -o "+fileNamePre+"_rateVar_console_output.txt -e "+fileNamePre+"_rateVar_console_error.txt "
									+"raxml-ng --search "+methodOptions[method][option]+" --msa "+pathToFile+".fa --msa-format FASTA --data-type DNA --redo --prefix "+fileNamePre+"_rateVar_output --model GTR+G --threads 1 \n")
							
							#MAPLE
							elif method==6:
								file.write("bsub -g /"+methodName+str(numSamples[0])+" -M "+str(int(400+j/20))+" -o "+fileNamePre+"_console_output.txt -e "+fileNamePre+"_console_error.txt "
								+""+pypy3Path+" "+folderCluster+"/fastLK/code/MAPLEv"+version+".py --reference "+refFile+methodOptions[method][option]+" --input "
								+pathToFile+".txt --overwrite --output "+fileNamePre+"\n")

							#run IQtree LK evaluation for all methods
							fileIQtreeEvaluation.write("cd "+folderNameSimu+"\n\t"+"bsub -M "+str(int(500+j*1.7))+" -o "+fileNamePre+"_IQtreeLK_console_output.txt -e "+fileNamePre+"_IQtreeLK_console_error.txt "
							+""+iqTreePath+" -s "+pathToFile+".phy -st DNA -te "+treeFiles[0]+" -pre "+methodName+methodOptionsNames[method][option]+"_IQtreeLK -m GTR -quiet -nt 1 -keep-ident -redo -blmin 0.000000005 \n")
							if len(treeFiles)>1:
								fileIQtreeEvaluation.write("cd "+folderNameSimu+"\n\t"+"bsub -M "+str(int(500+j*1.7))+" -o "+fileNamePre+"_rateVar_IQtreeLK_console_output.txt -e "+fileNamePre+"_rateVar_IQtreeLK_console_error.txt "
								+""+iqTreePath+" -s "+pathToFile+".phy -st DNA -te "+treeFiles[1]+" -pre "+methodName+methodOptionsNames[method][option]+"_rateVar_IQtreeLK -m GTR -quiet -nt 1 -keep-ident -redo -blmin 0.000000005 \n")
							if scenario==0 or scenario==2 or scenario==3:
								fileIQtreeEvaluation.write("cd "+folderNameSimu+"\n\t"+"bsub -M "+str(int(500+j*3.0))+" -o "+fileNamePre+"_IQtreeLK_rateVar_console_output.txt -e "+fileNamePre+"_IQtreeLK_rateVar_console_error.txt "
								+""+iqTreePath+" -s "+pathToFile+".phy -st DNA -te "+treeFiles[0]+" -pre "+methodName+methodOptionsNames[method][option]+"_IQtreeLK_rateVar -m GTR+G -quiet -nt 1 -keep-ident -redo -blmin 0.000000005 \n")
								if len(treeFiles)>1:
									fileIQtreeEvaluation.write("cd "+folderNameSimu+"\n\t"+"bsub -M "+str(int(500+j*3.0))+" -o "+fileNamePre+"_rateVar_IQtreeLK_rateVar_console_output.txt -e "+fileNamePre+"_rateVar_IQtreeLK_rateVar_console_error.txt "
									+""+iqTreePath+" -s "+pathToFile+".phy -st DNA -te "+treeFiles[1]+" -pre "+methodName+methodOptionsNames[method][option]+"_rateVar_IQtreeLK_rateVar -m GTR+G -quiet -nt 1 -keep-ident -redo -blmin 0.000000005 \n")
						
							#Maple LK evaluation for all methods
							mapleLKoptionNames=["_MapleLK","_MapleLK_unrest"]
							mapleLKoptions=[""," --model UNREST"]
							bLenScaleOption=""
							if method<2:
								bLenScaleOption=" --normalizeInputBLen 0.000033 "
							for mapleLKoptionIndex in range(len(mapleLKoptionNames)):
								fileMapleEvaluation.write("bsub -M "+str(int(4000+j/20))+" -o "+fileNamePre+mapleLKoptionNames[mapleLKoptionIndex]+"_console_output.txt -e "+fileNamePre+mapleLKoptionNames[mapleLKoptionIndex]+"_console_error.txt "
								+""+pypy3Path+" "+folderCluster+"/fastLK/code/MAPLEv"+versionForLK+".py --reference "+refFile+mapleLKoptions[mapleLKoptionIndex]+bLenScaleOption+" --inputTree "+treeFiles[0]+" --input "
								+pathToFile+".txt --calculateLKfinalTree --numTopologyImprovements 0 --noFastTopologyInitialSearch  --overwrite --output "+fileNamePre+mapleLKoptionNames[mapleLKoptionIndex]+"\n")
								if len(treeFiles)>1:
									fileMapleEvaluation.write("bsub -M "+str(int(4000+j/20))+" -o "+fileNamePre+"_rateVar"+mapleLKoptionNames[mapleLKoptionIndex]+"_console_output.txt -e "+fileNamePre+"_rateVar"+mapleLKoptionNames[mapleLKoptionIndex]+"_rateVar_console_error.txt "
									+""+pypy3Path+" "+folderCluster+"/fastLK/code/MAPLEv"+versionForLK+".py --reference "+refFile+mapleLKoptions[mapleLKoptionIndex]+bLenScaleOption+" --inputTree "+treeFiles[1]+" --input "
									+pathToFile+".txt --calculateLKfinalTree --numTopologyImprovements 0 --noFastTopologyInitialSearch  --overwrite --output "+fileNamePre+"_rateVar"+mapleLKoptionNames[mapleLKoptionIndex]+"\n")
								if scenario==0 or scenario==2 or scenario==3:
									fileMapleEvaluation.write("bsub -M "+str(int(8000+j/20))+" -o "+fileNamePre+mapleLKoptionNames[mapleLKoptionIndex]+"_rateVar_console_output.txt -e "+fileNamePre+mapleLKoptionNames[mapleLKoptionIndex]+"_rateVar_console_error.txt "
									+""+pypy3Path+" "+folderCluster+"/fastLK/code/MAPLEv"+versionForLK+".py --reference "+refFile+mapleLKoptions[mapleLKoptionIndex]+bLenScaleOption+" --inputTree "+treeFiles[0]+" --input "
									+pathToFile+".txt --calculateLKfinalTree --numTopologyImprovements 0 --noFastTopologyInitialSearch  --overwrite --rateVariation --output "+fileNamePre+mapleLKoptionNames[mapleLKoptionIndex]+"_rateVar\n")
									if len(treeFiles)>1:
										fileMapleEvaluation.write("bsub -M "+str(int(8000+j/20))+" -o "+fileNamePre+"_rateVar"+mapleLKoptionNames[mapleLKoptionIndex]+"_rateVar_console_output.txt -e "+fileNamePre+"_rateVar"+mapleLKoptionNames[mapleLKoptionIndex]+"_rateVar_console_error.txt "
										+""+pypy3Path+" "+folderCluster+"/fastLK/code/MAPLEv"+versionForLK+".py --reference "+refFile+mapleLKoptions[mapleLKoptionIndex]+bLenScaleOption+" --inputTree "+treeFiles[1]+" --input "
										+pathToFile+".txt --calculateLKfinalTree --numTopologyImprovements 0 --noFastTopologyInitialSearch  --overwrite --rateVariation --output "+fileNamePre+"_rateVar"+mapleLKoptionNames[mapleLKoptionIndex]+"_rateVar\n")

							#RF calculations
							fileRF.write("bsub -M "+str(int(400+j/20))+" -o "+fileNamePre+"_RF_console_output.txt -e "+fileNamePre+"_RF_console_error.txt "
							+""+pypy3Path+" "+folderCluster+"/fastLK/code/MAPLEv"+versionForRF+".py --inputRFtrees "+treeFiles[0]+" --inputTree "+pathToFile+"_realized.nw "
							+" --overwrite --output "+fileNamePre+"_RF\n")
							if len(treeFiles)>1:
								fileRF.write("bsub -M "+str(int(400+j/20))+" -o "+fileNamePre+"_rateVar_RF_console_output.txt -e "+fileNamePre+"_rateVar_RF_console_error.txt "
								+""+pypy3Path+" "+folderCluster+"/fastLK/code/MAPLEv"+versionForRF+".py --inputRFtrees "+treeFiles[1]+" --inputTree "+pathToFile+"_realized.nw "
								+" --overwrite --output "+fileNamePre+"_rateVar_RF\n")

							#parsimony calculations
							fileParsimony.write("bsub -M "+str(int(400+j/50))+" -o "+fileNamePre+"_parsimony_console_output.txt -e "+fileNamePre+"_parsimony_console_error.txt "
							+"usher -v "+pathToFile+".vcf -t "+treeFiles[0]+" -o "+fileNamePre+"_parsimony.mat \n")
							if len(treeFiles)>1:
								fileParsimony.write("bsub -M "+str(int(400+j/50))+" -o "+fileNamePre+"_rateVar_parsimony_console_output.txt -e "+fileNamePre+"_rateVar_parsimony_console_error.txt "
								+"usher -v "+pathToFile+".vcf -t "+treeFiles[0]+" -o "+fileNamePre+"_rateVar_parsimony.mat \n")

						file.write("done\n\n")
						fileIQtreeEvaluation.write("done\n\n")
						fileMapleEvaluation.write("done\n\n")
						fileRF.write("done\n\n")
						fileParsimony.write("done\n\n")
						if method==1:
							fileMatConv.write("done\n\n")
			file.close()
			fileIQtreeEvaluation.close()
			fileMapleEvaluation.close()
			fileRF.close()
			fileParsimony.close()
			print("Created "+methodName+" bash script "+args.pathToSimulationFolder+"submit"+methodName+".sh")
			print("Created "+methodName+" bash script "+args.pathToSimulationFolder+"submitIQtreeLK_"+methodName+".sh")
			print("Created "+methodName+" bash script "+args.pathToSimulationFolder+"submitMapleLK_"+methodName+".sh")
			print("Created "+methodName+" bash script "+args.pathToSimulationFolder+"submitRF_"+methodName+".sh")
			print("Created "+methodName+" bash script "+args.pathToSimulationFolder+"submitParsimony_"+methodName+".sh")
		fileMatConv.close()
		print("Created matConvert bash script "+args.pathToSimulationFolder+"submitMatOptimizeConversion.sh")





def printMinMeanMax(valueList):
	print(valueList)
	print(str(min(valueList))+" \t"+str(sum(valueList)/len(valueList))+" \t"+str(max(valueList)))








if args.collectResults:
	if args.runExtraAnalysesDivergence:
		ref = collectReference(args.inputSimulationReference)
		scenario=args.scenario
		subsampleTreeInference=args.subSampleNum
		divergences=["01X","02X","05X","1X","2X","5X","10X","20X","50X","100X","200X","500X","1000X"]
		folder="divergence"
		seed=args.repeat
		repeat=seed

		phylos=[args.pathToSimulationFolder+"phastSim_genomes_"+divergences[scenario]+"_realizedTree.tree"]
		pathToRepeat=args.pathToSimulationFolder+folder+'/'+divergences[scenario]+"/repeat"+str(repeat)+"_"+divergences[scenario]+"_"+str(subsampleTreeInference)+"samples"
		
		if args.numSamples==0:
			numSamples=[2000, 5000]
			numSamplesIQtreeLK=[2000, 5000]
		else:
			numSamples=[args.numSamples]
			numSamplesIQtreeLK=[args.numSamples]

		#MAPLEversions=["0.1.4","0.1.5","0.1.6","0.1.7","0.1.9"]
		MAPLEversions=["0.2.0"]
		#folders=["realDataSubsamples","simulationsSubsamples","simulations4catSubsamples","simulationsAlphaSubsamples","simulationsNsSubsamples"]

		MAPLEoptions=[""]
		MAPLEoptionsNames=[""]
		testOptions=False
		
		parameterValues=[""]
		parameterNames=[""]
		testParameters=False

		#collect results for different methods
		valueNames=["time","memory","maxMemory","LK","IQLK","RF","parsimony"]
		results={}
		for res in valueNames:
			results[res]={}
		methods=["UShER","matOptimize","IQtree","FastTree","RAxML-NG",'Maple']
		if args.methodToFocusOn==-1:
			methodsToFocusOn=methods
		else:
			methodsToFocusOn=[methods[args.methodToFocusOn]]
		methodOptions=[[""],[""],[" -fast"],[" -fastest"],[" --blmin 0.000005 --tree pars\{1\}"],[""]]
		methodOptionsNames=[[""],[""],["_fast"],["_fastest"],["_fast"],[""]]
		#methodOptions=[[""],[""],[" -pers 0.1 -nstop 500 -blmin 0.000000005"," -blmin 0.000000005"," -fast"],[" -spr 4 -mlacc 2 -slownni"," -spr 4",""," -fastest"],[""," -D"],[" --blmin 0.000000005 --tree pars\{3\}"," --blmin 0.000005 --tree pars\{1\}"],[""," --fast"," --rateVariation"," --model UNREST"," --model UNREST --rateVariation"]]
		#methodOptionsNames=[[""],[""],["_slow","_medium","_fast"],["_slow","_medium","_fast","_fastest"],["_slow","_fast"],["_slow","_fast"],["","_fast","_rateVar","_unrest","_unrest_rateVar"]]
		#maxNumSamples=[[200000],[200000],[5000,10000,20000],[20000,20000,20000,20000],[5000,5000],[5000,5000],[500000,500000,500000,500000,500000]]
		for method in range(len(methods)):
			methodName=methods[method]
			if methodName in methodsToFocusOn:
				for res in valueNames:
					results[res][methodName]={}
				for option in range(len(methodOptions[method])):
					for res in valueNames:
						results[res][methodName][option]={}
					for j in numSamples:
						#if j<=maxNumSamples[method][option]:
							for res in valueNames:
								results[res][methodName][option][j]={}
							for scenario in range(len(divergences)):
								for res in valueNames:
									results[res][methodName][option][j][scenario]=[[],[]]
								#folder=folders[scenario]
								#pathToFile=args.pathToSimulationFolder+folder+'/'+str(j)+"subsamples/repeat\"$i\"_"+str(j)+"samples_"+folder
								#pathToTreeFile=args.pathToSimulationFolder+folder+'/'+str(j)+"subsamples/repeat\"$i\"_"+str(j)+"samples_"+folder
								pathToFile=args.pathToSimulationFolder+folder+'/'+divergences[scenario]+"/repeat\"$i\"_"+divergences[scenario]+"_"+str(j)+"samples"
								#pathToFile=args.pathToSimulationFolder+folder+'/'+str(j)+"subsamples/repeat\"$i\"_"+str(j)+"samples_"+folder
								rateVarStr=[""]
								outputNames=["_console_output.txt","_console_error.txt"]
								IQtreeLKoptionNames=["_IQtreeLK"]
								mapleLKoptionNames=["_MapleLK"]

								refFile=args.inputSimulationReference
								
								notFound=[["",""],["",""],["",""],["",""],["",""],["",""],["",""],["",""],["",""],["",""],["",""],["",""],["",""]]
								notFoundCount=[[0,0],[0,0],[0,0],[0,0],[0,0],[0,0],[0,0],[0,0],[0,0],[0,0],[0,0],[0,0],[0,0],[0,0]]
								outputFiles=[]
								for res in valueNames:
									outputFiles.append([open(args.pathToSimulationFolder+folder+'/'+divergences[scenario]+"/"+methodName+methodOptionsNames[method][option]+"_"+str(j)+"samples_"+res+"_summary.txt","w")])

								for i in range(10):
									#folderNameSimu=args.pathToSimulationFolder+folder+'/'+str(j)+"subsamples/output_repl"+str(i+1)+"/"
									folderNameSimu=args.pathToSimulationFolder+folder+'/'+divergences[scenario]+"/output_repl"+str(i+1)+"/"
									#fileNamePre=folderNameSimu+methodName+methodOptionsNames[method][option]
									fileNamePre=folderNameSimu+methodName+methodOptionsNames[method][option]+"_"+str(j)+"samples"

									if method==0:
										treeFiles=[folderNameSimu+str(j)+"samples_UShER_output/final-tree.nh"]
									elif method==1:
										treeFiles=[folderNameSimu+str(j)+"samples_matOptimize_output.nh"]
									elif method==2:
										treeFiles=[fileNamePre+".treefile"]
									elif method==3:
										treeFiles=[fileNamePre+"_console_output.txt"]
									elif method==4:
										treeFiles=[fileNamePre+"_output.raxml.bestTree"]
									elif method==5:
										treeFiles=[fileNamePre+"_tree.tree"]
									for treeFileNum in range(len(treeFiles)):
										treeFile=treeFiles[treeFileNum]

										if not os.path.isfile(treeFile):
											notFound[0][treeFileNum]+="-"+str(i+1)
											notFoundCount[0][treeFileNum]+=1
											for res in valueNames:
												results[res][methodName][option][j][scenario][treeFileNum].append(float("nan"))
										else:
											#collect runtime and memory demands
											file=open(fileNamePre+"_console_output.txt")
											line=file.readline()
											while line!="Resource usage summary:\n" and line!="":
												line=file.readline()
											if line!="":
												line=file.readline()
												line=file.readline()
												timeRun=float(line.split()[3])
												results["time"][methodName][option][j][scenario][treeFileNum].append(float(line.split()[3]))
												line=file.readline()
												if line.split()[3]=="-":
													results["memory"][methodName][option][j][scenario][treeFileNum].append(float("nan"))
													results["maxMemory"][methodName][option][j][scenario][treeFileNum].append(float("nan"))
													#maxMemoryUsed=float("NaN")
													#aveMemoryUsed=float("NaN")
												else:
													results["maxMemory"][methodName][option][j][scenario][treeFileNum].append(float(line.split()[3]))
													#maxMemoryUsed=float(line.split()[3])
													line=file.readline()
													if line.split()[3]=="-":
														results["memory"][methodName][option][j][scenario][treeFileNum].append(float("nan"))
													else:
														results["memory"][methodName][option][j][scenario][treeFileNum].append(float(line.split()[3]))
											else:
												results["time"][methodName][option][j][scenario][treeFileNum].append(float("nan"))
												results["memory"][methodName][option][j][scenario][treeFileNum].append(float("nan"))
												results["maxMemory"][methodName][option][j][scenario][treeFileNum].append(float("nan"))
												notFound[1][treeFileNum]+="-"+str(i+1)
												notFoundCount[1][treeFileNum]+=1
											file.close()
											
											#collect LKs
											mapleLKoptionNames=["_MapleLK"]
											rateVarStrLK=[""]
											for mapleLKoptionIndex in range(len(mapleLKoptionNames)):
												for rateVarStrNum in range(len(rateVarStrLK)):
													if not os.path.isfile(fileNamePre+rateVarStr[treeFileNum]+mapleLKoptionNames[mapleLKoptionIndex]+rateVarStr[rateVarStrNum]+"_LK.txt"):
														notFound[2+mapleLKoptionIndex+2*rateVarStrNum][treeFileNum]+="-"+str(i+1)
														notFoundCount[2+mapleLKoptionIndex+2*rateVarStrNum][treeFileNum]+=1
														results[valueNames[3+mapleLKoptionIndex+2*rateVarStrNum]][methodName][option][j][scenario][treeFileNum].append(float('nan'))
													else:
														file=open(fileNamePre+rateVarStr[treeFileNum]+mapleLKoptionNames[mapleLKoptionIndex]+rateVarStr[rateVarStrNum]+"_LK.txt")
														line=file.readline()
														results[valueNames[3+mapleLKoptionIndex+2*rateVarStrNum]][methodName][option][j][scenario][treeFileNum].append(float(line.split()[0]))
														#lk=float(line.split()[0])
														#lks[versNum][j][scenario][option][paramNum].append(lk)
														file.close()

											#collect IQtree LKs
											for rateVarStrNum in range(len(rateVarStrLK)):
												if not os.path.isfile(fileNamePre+rateVarStr[treeFileNum]+"_IQtreeLK"+rateVarStr[rateVarStrNum]+".log"):
													notFound[3+rateVarStrNum][treeFileNum]+="-"+str(i+1)
													notFoundCount[3+rateVarStrNum][treeFileNum]+=1
													results[valueNames[4+rateVarStrNum]][methodName][option][j][scenario][treeFileNum].append(float('nan'))
													#notFoundIQLK+="-"+str(i+1)
													#notFoundIQLKCount+=1
												else:
													file=open(fileNamePre+rateVarStr[treeFileNum]+"_IQtreeLK"+rateVarStr[rateVarStrNum]+".log")
													line=file.readline()
													while not ("BEST SCORE FOUND :" in line) and line!="":
														line=file.readline()
													if line=="":
														notFound[3+rateVarStrNum][treeFileNum]+="-"+str(i+1)
														notFoundCount[3+rateVarStrNum][treeFileNum]+=1
														results[valueNames[4+rateVarStrNum]][methodName][option][j][scenario][treeFileNum].append(float('nan'))
														#notFoundIQLK+="-"+str(i+1)
														#notFoundIQLKCount+=1
													else:
														#IQlk=float(line.split()[4])
														#IQlks[versNum][j][scenario][option][paramNum].append(IQlk)
														results[valueNames[4+rateVarStrNum]][methodName][option][j][scenario][treeFileNum].append(float(line.split()[4]))
													file.close()

											#collect RF distances
											if not os.path.isfile(fileNamePre+rateVarStr[treeFileNum]+"_RF_RFdistances.txt"):
												#notFoundRF+="-"+str(i+1)
												#notFoundRFCount+=1
												notFound[4][treeFileNum]+="-"+str(i+1)
												notFoundCount[4][treeFileNum]+=1
												results[valueNames[5]][methodName][option][j][scenario][treeFileNum].append(float('nan'))
											else:
												file=open(fileNamePre+rateVarStr[treeFileNum]+"_RF_RFdistances.txt")
												line=file.readline()
												line=file.readline()
												if len(line.split())<2 or line.split()[1]=='None':
													print(fileNamePre+rateVarStr[treeFileNum]+"_RF_RFdistances.txt")
													print(line)
													#notFoundRF+="-"+str(i+1)
													#notFoundRFCount+=1
													notFound[4][treeFileNum]+="-"+str(i+1)
													notFoundCount[4][treeFileNum]+=1
													results[valueNames[5]][methodName][option][j][scenario][treeFileNum].append(float('nan'))
												else:
													results[valueNames[5]][methodName][option][j][scenario][treeFileNum].append(float(line.split()[1]))
													#RFdistance=float(line.split()[1])
													#RFs[versNum][j][scenario][option][paramNum].append(RFdistance)
												file.close()

											#valueNames=["time","memory","maxMemory","LK","LKunrest","LKsiteVar","LKunrestSiteVar","IQLK","IQLKsiteVar","RF","parsimony"]
											if not os.path.isfile(fileNamePre+rateVarStr[treeFileNum]+"_parsimony_console_error.txt"):
												notFound[5][treeFileNum]+="-"+str(i+1)
												notFoundCount[5][treeFileNum]+=1
												results[valueNames[6]][methodName][option][j][scenario][treeFileNum].append(float('nan'))
											else:
												file=open(fileNamePre+rateVarStr[treeFileNum]+"_parsimony_console_error.txt")
												line=file.readline()
												while (not ("The parsimony score for this tree is:" in line)) and line!="":
													line=file.readline()
												if line=="":
													notFound[5][treeFileNum]+="-"+str(i+1)
													notFoundCount[5][treeFileNum]+=1
													results[valueNames[6]][methodName][option][j][scenario][treeFileNum].append(float('nan'))
												else:
													results[valueNames[6]][methodName][option][j][scenario][treeFileNum].append(float(line.split()[7]))
												file.close()

								fileTypes=["tree","time-memory","LK","IQLK","RF","parsimony"]
								for treeFileNum in range(1):
									if notFoundCount[0][treeFileNum]<10:
										for fileTypeNum in range(len(fileTypes)):
											if notFoundCount[fileTypeNum][treeFileNum]>0:
												print("Scenario "+folder+" "+divergences[scenario]+" "+methodName+methodOptionsNames[method][option]+rateVarStr[treeFileNum]+", "+str(j)+" samples, "+fileTypes[fileTypeNum]+" file not found for replicates "+notFound[fileTypeNum][treeFileNum])
									
								#for treeFileNum in range(2):
										if len(results[valueNames[0]][methodName][option][j][scenario][treeFileNum])>0:
											print("\n Scenario "+folder+" "+divergences[scenario]+" "+methodName+methodOptionsNames[method][option]+rateVarStr[treeFileNum]+", "+str(j)+" samples : times, memory, max memory, LKs, IQlk, RF, parsimony")
											for valueNum in range(len(valueNames)):
												print(valueNames[valueNum])
												try:
													printMinMeanMax(results[valueNames[valueNum]][methodName][option][j][scenario][treeFileNum])
												except:
													print(results[valueNames[valueNum]][methodName][option][j][scenario][treeFileNum])
											print("\n")
									else:
										print("no results files found for scenario "+folder+" "+divergences[scenario]+" "+methodName+methodOptionsNames[method][option]+rateVarStr[treeFileNum]+", "+str(j)+" samples \n")

								for treeFileNum in range(1):
									for valueNum in range(len(valueNames)):
										for val in results[valueNames[valueNum]][methodName][option][j][scenario][treeFileNum]:
											outputFiles[valueNum][treeFileNum].write(str(val)+"\t")
										outputFiles[valueNum][treeFileNum].close()
									#outputFiles.append([open(args.pathToSimulationFolder+folder+'/'+str(j)+"subsamples/"+methodName+methodOptionsNames[method][option]+"_"+res+"_summary.txt","w"),open(args.pathToSimulationFolder+folder+'/'+str(j)+"subsamples/"+methodName+methodOptionsNames[method][option]+"_rateVar_"+res+"_summary.txt","w")])

						




	else:
		#printOptions=True
		#MAPLEversions=["0.1.4","0.1.5","0.1.6","0.1.7","0.1.8"]
		MAPLEversions=["0.1.9"]
		if args.numSamples==0:
			numSamples=[1000, 2000, 5000, 10000, 20000, 50000, 100000, 200000, 500000]
		else:
			numSamples=[args.numSamples]

		folders=["realDataSubsamples","simulationsSubsamples","simulations4catSubsamples","simulationsAlphaSubsamples","simulationsNsSubsamples"]

		testOptions=False
		MAPLEoptionsNames=["","_fast"]
		if testOptions:
			MAPLEoptionsNames=["","_fast","_rateVar","_fast_rateVar","_unrest","_fast_unrest","_unrest_rateVar","_fast_unrest_rateVar"]
		
		testParameters=False
		parameterNames=[""]
		if testParameters:
			#parameterValues=[" --minBLenSensitivity 0.1"," --minBLenSensitivity 0.01"," --minBLenSensitivity 0.001"," --minBLenSensitivity 0.0001"," --minBLenSensitivity 0.00001"," --minBLenSensitivity 0.000001"," --factorOptimizePlacementLKvsSearchLK 0.004"," --factorOptimizePlacementLKvsSearchLK 0.01"," --factorOptimizePlacementLKvsSearchLK 0.04"," --factorOptimizePlacementLKvsSearchLK 0.1"," --factorOptimizePlacementLKvsSearchLK 0.2"," --factorOptimizePlacementLKvsSearchLK 0.4"]
			parameterNames=["_minBLen01","_minBLen001","_minBLen0001","_LK18","_LK12","_LK6","_LK1","_LK01","_LKtopology18","_LKtopology12","_LKtopology6","_LKtopology1","_LKtopology01"]
			folders=["realDataSubsamples"]

		if testOptions or testParameters:
			#collecting the results from the benchmarking
			times={}
			memories={}
			maxMemories={}
			lks={}
			lksUnrest={}
			lksSiteVar={}
			lksSiteVarUnrest={}
			IQlks={}
			RFs={}
			for versNum in range(len(MAPLEversions)):
				version=MAPLEversions[versNum]
				times[versNum]={}
				memories[versNum]={}
				maxMemories[versNum]={}
				lks[versNum]={}
				lksUnrest[versNum]={}
				lksSiteVar[versNum]={}
				lksSiteVarUnrest[versNum]={}
				IQlks[versNum]={}
				RFs[versNum]={}
				for j in numSamples:
					times[versNum][j]={}
					memories[versNum][j]={}
					maxMemories[versNum][j]={}
					lks[versNum][j]={}
					lksUnrest[versNum][j]={}
					lksSiteVar[versNum][j]={}
					lksSiteVarUnrest[versNum][j]={}
					IQlks[versNum][j]={}
					RFs[versNum][j]={}
					for scenario in range(len(folders)):
						times[versNum][j][scenario]={}
						memories[versNum][j][scenario]={}
						maxMemories[versNum][j][scenario]={}
						lks[versNum][j][scenario]={}
						lksUnrest[versNum][j][scenario]={}
						lksSiteVar[versNum][j][scenario]={}
						lksSiteVarUnrest[versNum][j][scenario]={}
						IQlks[versNum][j][scenario]={}
						RFs[versNum][j][scenario]={}
						folder=folders[scenario]
						for option in range(len(MAPLEoptionsNames)):
							times[versNum][j][scenario][option]={}
							memories[versNum][j][scenario][option]={}
							maxMemories[versNum][j][scenario][option]={}
							lks[versNum][j][scenario][option]={}
							lksUnrest[versNum][j][scenario][option]={}
							lksSiteVar[versNum][j][scenario][option]={}
							lksSiteVarUnrest[versNum][j][scenario][option]={}
							IQlks[versNum][j][scenario][option]={}
							RFs[versNum][j][scenario][option]={}
							for paramNum in range(len(parameterNames)):
								times[versNum][j][scenario][option][paramNum]=[]
								memories[versNum][j][scenario][option][paramNum]=[]
								maxMemories[versNum][j][scenario][option][paramNum]=[]
								lks[versNum][j][scenario][option][paramNum]=[]
								lksUnrest[versNum][j][scenario][option][paramNum]=[]
								lksSiteVar[versNum][j][scenario][option][paramNum]=[]
								lksSiteVarUnrest[versNum][j][scenario][option][paramNum]=[]
								IQlks[versNum][j][scenario][option][paramNum]=[]
								RFs[versNum][j][scenario][option][paramNum]=[]
								notFound=""
								notFoundCount=0
								notFoundTime=""
								notFoundTimeCount=0
								notFoundRF=""
								notFoundRFCount=0
								notFoundLK=""
								notFoundLKCount=0
								notFoundLKunrest=""
								notFoundLKunrestCount=0
								notFoundLKsiteVarUnrest=""
								notFoundLKsiteVarUnrestCount=0
								notFoundLKsiteVar=""
								notFoundLKsiteVarCount=0
								notFoundIQLK=""
								notFoundIQLKCount=0
								timeFile=open(args.pathToSimulationFolder+folder+'/'+str(j)+"subsamples/MAPLE"+version+MAPLEoptionsNames[option]+parameterNames[paramNum]+"_times.txt","w")
								memoryFile=open(args.pathToSimulationFolder+folder+'/'+str(j)+"subsamples/MAPLE"+version+MAPLEoptionsNames[option]+parameterNames[paramNum]+"_memory.txt","w")
								maxMemoryFile=open(args.pathToSimulationFolder+folder+'/'+str(j)+"subsamples/MAPLE"+version+MAPLEoptionsNames[option]+parameterNames[paramNum]+"_maxMemory.txt","w")
								lkFile=open(args.pathToSimulationFolder+folder+'/'+str(j)+"subsamples/MAPLE"+version+MAPLEoptionsNames[option]+parameterNames[paramNum]+"_lk.txt","w")
								lkUnrestFile=open(args.pathToSimulationFolder+folder+'/'+str(j)+"subsamples/MAPLE"+version+MAPLEoptionsNames[option]+parameterNames[paramNum]+"_lkUnrest.txt","w")
								if scenario==0:
									lkSiteVarUnrestFile=open(args.pathToSimulationFolder+folder+'/'+str(j)+"subsamples/MAPLE"+version+MAPLEoptionsNames[option]+parameterNames[paramNum]+"_lkSiteVarUnrest.txt","w")
									lkSiteVarFile=open(args.pathToSimulationFolder+folder+'/'+str(j)+"subsamples/MAPLE"+version+MAPLEoptionsNames[option]+parameterNames[paramNum]+"_lkSiteVar.txt","w")
								iqtreeLKFile=open(args.pathToSimulationFolder+folder+'/'+str(j)+"subsamples/MAPLE"+version+MAPLEoptionsNames[option]+parameterNames[paramNum]+"_IQtreeLK.txt","w")
								if scenario>0:
									rfFile=open(args.pathToSimulationFolder+folder+'/'+str(j)+"subsamples/MAPLE"+version+MAPLEoptionsNames[option]+parameterNames[paramNum]+"_RF.txt","w")
								for i in range(10):
									folderNameSimu=args.pathToSimulationFolder+folder+'/'+str(j)+"subsamples/output_repl"+str(i+1)+"/"
									if not os.path.isfile(folderNameSimu+"MAPLE"+version+MAPLEoptionsNames[option]+parameterNames[paramNum]+"_tree.tree"):
										notFound+="-"+str(i+1)
										notFoundCount+=1
									else:

										#collect runtime and memory demands
										file=open(folderNameSimu+"MAPLE"+version+MAPLEoptionsNames[option]+parameterNames[paramNum]+"_console_output.txt")
										line=file.readline()
										timeRun=-1.0
										while line!="":
											if "Final Substitution matrix:" in line:
												#bestLK=float(line.split()[1].replace(" ","").replace("\n",""))
												while line!="Resource usage summary:\n":
													line=file.readline()
												line=file.readline()
												line=file.readline()
												timeRun=float(line.split()[3])
												line=file.readline()
												#print("FastLK tree file with "+str(j)+" samples, replicate "+str(i+1)+" speed "+s+" ")
												#print(line)
												if line.split()[3]=="-":
													maxMemoryUsed=float("NaN")
													aveMemoryUsed=float("NaN")
													continue
												maxMemoryUsed=float(line.split()[3])
												line=file.readline()
												if line.split()[3]=="-":
													aveMemoryUsed=float("NaN")
													continue
												aveMemoryUsed=float(line.split()[3])
											line=file.readline()
										if timeRun<0.0:
											notFoundTime+="-"+str(i+1)
											notFoundTimeCount+=1
											#print("MAPLE"+versions[s]+" with "+str(j)+" samples "+simulationsMAPLE[simu]+", replicate "+str(i+1)+" - running time not found.")
										else:
											times[versNum][j][scenario][option][paramNum].append(timeRun)
											memories[versNum][j][scenario][option][paramNum].append(aveMemoryUsed)
											maxMemories[versNum][j][scenario][option][paramNum].append(maxMemoryUsed)
											#lks[versNum][j][scenario][option]=[]
											#IQlks[versNum][j][scenario][option]=[]
											#lks[j][s].append(bestLK)
										file.close()

										#collect LKs
										if not os.path.isfile(folderNameSimu+"MAPLE"+version+MAPLEoptionsNames[option]+parameterNames[paramNum]+"_LKestimation_LK.txt"):
											notFoundLK+="-"+str(i+1)
											notFoundLKCount+=1
										else:
											file=open(folderNameSimu+"MAPLE"+version+MAPLEoptionsNames[option]+parameterNames[paramNum]+"_LKestimation_LK.txt")
											line=file.readline()
											lk=float(line.split()[0])
											lks[versNum][j][scenario][option][paramNum].append(lk)
											file.close()

										#collect LKs estimated under UNREST model
										if not os.path.isfile(folderNameSimu+"MAPLE"+version+MAPLEoptionsNames[option]+parameterNames[paramNum]+"_LKestimation_unrest_LK.txt"):
											notFoundLKunrest+="-"+str(i+1)
											notFoundLKunrestCount+=1
										else:
											file=open(folderNameSimu+"MAPLE"+version+MAPLEoptionsNames[option]+parameterNames[paramNum]+"_LKestimation_unrest_LK.txt")
											line=file.readline()
											lk=float(line.split()[0])
											lksUnrest[versNum][j][scenario][option][paramNum].append(lk)
											file.close()

										if scenario==0:
											#collect also likelihoods estimated under rate variation for the real data
											if not os.path.isfile(folderNameSimu+"MAPLE"+version+MAPLEoptionsNames[option]+parameterNames[paramNum]+"_LKestimation_unrest_rateVar_LK.txt"):
												notFoundLKsiteVarUnrest+="-"+str(i+1)
												notFoundLKsiteVarUnrestCount+=1
											else:
												file=open(folderNameSimu+"MAPLE"+version+MAPLEoptionsNames[option]+parameterNames[paramNum]+"_LKestimation_unrest_rateVar_LK.txt")
												line=file.readline()
												lk=float(line.split()[0])
												lksSiteVarUnrest[versNum][j][scenario][option][paramNum].append(lk)
												file.close()

											if not os.path.isfile(folderNameSimu+"MAPLE"+version+MAPLEoptionsNames[option]+parameterNames[paramNum]+"_LKestimation_rateVar_LK.txt"):
												notFoundLKsiteVar+="-"+str(i+1)
												notFoundLKsiteVarCount+=1
											else:
												file=open(folderNameSimu+"MAPLE"+version+MAPLEoptionsNames[option]+parameterNames[paramNum]+"_LKestimation_rateVar_LK.txt")
												line=file.readline()
												lk=float(line.split()[0])
												lksSiteVar[versNum][j][scenario][option][paramNum].append(lk)
												file.close()

										#collect IQtree LKs
										if not os.path.isfile(folderNameSimu+"MAPLE"+version+MAPLEoptionsNames[option]+parameterNames[paramNum]+"_IQtreeLKestimation.log"):
											notFoundIQLK+="-"+str(i+1)
											notFoundIQLKCount+=1
										else:
											file=open(folderNameSimu+"MAPLE"+version+MAPLEoptionsNames[option]+parameterNames[paramNum]+"_IQtreeLKestimation.log")
											line=file.readline()
											while not ("BEST SCORE FOUND :" in line) and line!="":
												line=file.readline()
											if line=="":
												notFoundIQLK+="-"+str(i+1)
												notFoundIQLKCount+=1
											else:
												IQlk=float(line.split()[4])
												IQlks[versNum][j][scenario][option][paramNum].append(IQlk)
											file.close()

										#collect RF distances
										if scenario>0:
											if not os.path.isfile(folderNameSimu+"MAPLE"+version+MAPLEoptionsNames[option]+parameterNames[paramNum]+"_RFdistances.txt"):
												notFoundRF+="-"+str(i+1)
												notFoundRFCount+=1
											else:
												file=open(folderNameSimu+"MAPLE"+version+MAPLEoptionsNames[option]+parameterNames[paramNum]+"_RFdistances.txt")
												line=file.readline()
												line=file.readline()
												if line.split()[1]=='None':
													print(folderNameSimu+"MAPLE"+version+MAPLEoptionsNames[option]+parameterNames[paramNum]+"_RFdistances.txt")
													print(line)
													notFoundRF+="-"+str(i+1)
													notFoundRFCount+=1
												else:
													RFdistance=float(line.split()[1])
													RFs[versNum][j][scenario][option][paramNum].append(RFdistance)
												file.close()

								if notFoundCount>0:
									print("Scenario "+folder+" MAPLE"+version+", "+str(j)+" samples, options "+MAPLEoptionsNames[option]+parameterNames[paramNum]+", tree file not found for replicates "+notFound)
								if notFoundTimeCount>0:
									print("Scenario "+folder+" MAPLE"+version+", "+str(j)+" samples, options "+MAPLEoptionsNames[option]+parameterNames[paramNum]+", runtime not found for replicates "+notFoundTime)
								if notFoundLKCount>0:
									print("Scenario "+folder+" MAPLE"+version+", "+str(j)+" samples, options "+MAPLEoptionsNames[option]+parameterNames[paramNum]+", LK file not found for replicates "+notFoundLK)
								if notFoundLKunrestCount>0:
									print("Scenario "+folder+" MAPLE"+version+", "+str(j)+" samples, options "+MAPLEoptionsNames[option]+parameterNames[paramNum]+", LK under unrest file not found for replicates "+notFoundLKunrest)
								if notFoundLKsiteVarUnrestCount>0:
									print("Scenario "+folder+" MAPLE"+version+", "+str(j)+" samples, options "+MAPLEoptionsNames[option]+parameterNames[paramNum]+", LK under siteVar-unrest file not found for replicates "+notFoundLKsiteVarUnrest)
								if notFoundLKsiteVarCount>0:
									print("Scenario "+folder+" MAPLE"+version+", "+str(j)+" samples, options "+MAPLEoptionsNames[option]+parameterNames[paramNum]+", LK under siteVar file not found for replicates "+notFoundLKsiteVar)
								if notFoundIQLKCount>0:
									print("Scenario "+folder+" MAPLE"+version+", "+str(j)+" samples, options "+MAPLEoptionsNames[option]+parameterNames[paramNum]+", IQtree LK file not found for replicates "+notFoundIQLK)
								if scenario>0 and notFoundRFCount>0:
									print("Scenario "+folder+" MAPLE"+version+", "+str(j)+" samples, options "+MAPLEoptionsNames[option]+parameterNames[paramNum]+", RF file not found for replicates "+notFoundRF)
								if len(times[versNum][j][scenario][option][paramNum])>0:
									print("Scenario "+folder+" MAPLE"+version+", "+str(j)+" samples, options "+MAPLEoptionsNames[option]+parameterNames[paramNum]+" : times, memory, max memory, LKs, unrest  LKs, IQtree LKs, RFdistances")
									printMinMeanMax(times[versNum][j][scenario][option][paramNum])
									printMinMeanMax(memories[versNum][j][scenario][option][paramNum])
									printMinMeanMax(maxMemories[versNum][j][scenario][option][paramNum])
								if len(lks[versNum][j][scenario][option][paramNum])>0:
									printMinMeanMax(lks[versNum][j][scenario][option][paramNum])
								if len(lksUnrest[versNum][j][scenario][option][paramNum])>0:
									printMinMeanMax(lksUnrest[versNum][j][scenario][option][paramNum])
								if scenario==0:
									if len(lksSiteVarUnrest[versNum][j][scenario][option][paramNum])>0:
										printMinMeanMax(lksSiteVarUnrest[versNum][j][scenario][option][paramNum])
									if len(lksSiteVar[versNum][j][scenario][option][paramNum])>0:
										printMinMeanMax(lksSiteVar[versNum][j][scenario][option][paramNum])
								if len(IQlks[versNum][j][scenario][option][paramNum])>0:
									printMinMeanMax(IQlks[versNum][j][scenario][option][paramNum])
								if scenario>0 and len(RFs[versNum][j][scenario][option][paramNum])>0:
									printMinMeanMax(RFs[versNum][j][scenario][option][paramNum])
								print("\n")


								for time in times[versNum][j][scenario][option][paramNum]:
									timeFile.write(str(time)+"\t")
								timeFile.close()
								for time in memories[versNum][j][scenario][option][paramNum]:
									memoryFile.write(str(time)+"\t")
								memoryFile.close()
								for time in maxMemories[versNum][j][scenario][option][paramNum]:
									maxMemoryFile.write(str(time)+"\t")
								maxMemoryFile.close()
								for time in lks[versNum][j][scenario][option][paramNum]:
									lkFile.write(str(time)+"\t")
								lkFile.close()
								for time in lksUnrest[versNum][j][scenario][option][paramNum]:
									lkUnrestFile.write(str(time)+"\t")
								lkUnrestFile.close()
								if scenario==0:
									for time in lksSiteVarUnrest[versNum][j][scenario][option][paramNum]:
										lkSiteVarUnrestFile.write(str(time)+"\t")
									lkSiteVarUnrestFile.close()
									for time in lksSiteVar[versNum][j][scenario][option][paramNum]:
										lkSiteVarFile.write(str(time)+"\t")
									lkSiteVarFile.close()
								for time in IQlks[versNum][j][scenario][option][paramNum]:
									iqtreeLKFile.write(str(time)+"\t")
								iqtreeLKFile.close()
								if scenario>0:
									for time in RFs[versNum][j][scenario][option][paramNum]:
										rfFile.write(str(time)+"\t")
									rfFile.close()
								


		#collect results for different methods
		valueNames=["time","memory","maxMemory","LK","LKunrest","LKsiteVar","LKunrestSiteVar","IQLK","IQLKsiteVar","RF","parsimony"]
		results={}
		for res in valueNames:
			results[res]={}
		methods=["UShER","matOptimize","IQtree","FastTree","RAxML","RAxML-NG",'Maple']
		if args.methodToFocusOn==-1:
			methodsToFocusOn=methods
		else:
			methodsToFocusOn=[methods[args.methodToFocusOn]]
		methodOptions=[[""],[""],[" -pers 0.1 -nstop 500 -blmin 0.000000005"," -blmin 0.000000005"," -fast"],[" -spr 4 -mlacc 2 -slownni"," -spr 4",""," -fastest"],[""," -D"],[" --blmin 0.000000005 --tree pars\{3\}"," --blmin 0.000005 --tree pars\{1\}"],[""," --fast"," --rateVariation"," --model UNREST"," --model UNREST --rateVariation"]]
		methodOptionsNames=[[""],[""],["_slow","_medium","_fast"],["_slow","_medium","_fast","_fastest"],["_slow","_fast"],["_slow","_fast"],["","_fast","_rateVar","_unrest","_unrest_rateVar"]]
		maxNumSamples=[[200000],[200000],[5000,10000,20000],[20000,20000,20000,20000],[5000,5000],[5000,5000],[500000,500000,500000,500000,500000]]
		for method in range(len(methods)):
			methodName=methods[method]
			if methodName in methodsToFocusOn:
				for res in valueNames:
					results[res][methodName]={}
				for option in range(len(methodOptions[method])):
					for res in valueNames:
						results[res][methodName][option]={}
					for j in numSamples:
						if j<=maxNumSamples[method][option]:
							for res in valueNames:
								results[res][methodName][option][j]={}
							for scenario in range(len(folders)):
								for res in valueNames:
									results[res][methodName][option][j][scenario]=[[],[]]
								folder=folders[scenario]
								pathToFile=args.pathToSimulationFolder+folder+'/'+str(j)+"subsamples/repeat\"$i\"_"+str(j)+"samples_"+folder
								rateVarStr=["","_rateVar"]
								outputNames=["_console_output.txt","_console_error.txt"]
								IQtreeLKoptionNames=["_IQtreeLK","_IQtreeLK_rateVar"]
								mapleLKoptionNames=["_MapleLK","_MapleLK_rateVar","_MapleLK_unrest","_MapleLK_unrest_rateVar"]
								
								notFound=[["",""],["",""],["",""],["",""],["",""],["",""],["",""],["",""],["",""],["",""],["",""],["",""],["",""]]
								notFoundCount=[[0,0],[0,0],[0,0],[0,0],[0,0],[0,0],[0,0],[0,0],[0,0],[0,0],[0,0],[0,0],[0,0],[0,0]]
								outputFiles=[]
								for res in valueNames:
									outputFiles.append([open(args.pathToSimulationFolder+folder+'/'+str(j)+"subsamples/"+methodName+methodOptionsNames[method][option]+"_"+res+"_summary.txt","w"),open(args.pathToSimulationFolder+folder+'/'+str(j)+"subsamples/"+methodName+methodOptionsNames[method][option]+"_rateVar_"+res+"_summary.txt","w")])

								for i in range(10):
									folderNameSimu=args.pathToSimulationFolder+folder+'/'+str(j)+"subsamples/output_repl"+str(i+1)+"/"
									fileNamePre=folderNameSimu+methodName+methodOptionsNames[method][option]
									if method==0:
										treeFiles=[folderNameSimu+"UShER_output/final-tree.nh"]
									elif method==1:
										treeFiles=[folderNameSimu+"matOptimize_output.nh"]
									elif method==2:
										treeFiles=[fileNamePre+".treefile"]
										if scenario==0 or scenario==2 or scenario==3:
											treeFiles.append(fileNamePre+"_rateVar.treefile")
									elif method==3:
										treeFiles=[fileNamePre+"_console_output.txt"]
										if scenario==0 or scenario==2 or scenario==3:
											treeFiles.append(fileNamePre+"_rateVar_console_output.txt")
									elif method==4:
										treeFiles=[folderNameSimu+"RAxML_result."+methodName+methodOptionsNames[method][option]+"_output"]
										if scenario==0 or scenario==2 or scenario==3:
											treeFiles.append(folderNameSimu+"RAxML_result."+methodName+methodOptionsNames[method][option]+"_rateVar_output")
									elif method==5:
										treeFiles=[fileNamePre+"_output.raxml.bestTree"]
										if scenario==0 or scenario==2 or scenario==3:
											treeFiles.append(fileNamePre+"_rateVar_output.raxml.bestTree")
									elif method==6:
										treeFiles=[fileNamePre+"_tree.tree"]
									for treeFileNum in range(len(treeFiles)):
										treeFile=treeFiles[treeFileNum]

										if not os.path.isfile(treeFile):
											notFound[0][treeFileNum]+="-"+str(i+1)
											notFoundCount[0][treeFileNum]+=1
											for res in valueNames:
												results[res][methodName][option][j][scenario][treeFileNum].append(float("nan"))
										else:
											#collect runtime and memory demands
											if treeFileNum==0:
												file=open(fileNamePre+"_console_output.txt")
											else:
												file=open(fileNamePre+"_rateVar_console_output.txt")
											line=file.readline()
											while line!="Resource usage summary:\n" and line!="":
												line=file.readline()
											if line!="":
												line=file.readline()
												line=file.readline()
												timeRun=float(line.split()[3])
												results["time"][methodName][option][j][scenario][treeFileNum].append(float(line.split()[3]))
												line=file.readline()
												if line.split()[3]=="-":
													results["memory"][methodName][option][j][scenario][treeFileNum].append(float("nan"))
													results["maxMemory"][methodName][option][j][scenario][treeFileNum].append(float("nan"))
													#maxMemoryUsed=float("NaN")
													#aveMemoryUsed=float("NaN")
												else:
													results["maxMemory"][methodName][option][j][scenario][treeFileNum].append(float(line.split()[3]))
													#maxMemoryUsed=float(line.split()[3])
													line=file.readline()
													if line.split()[3]=="-":
														results["memory"][methodName][option][j][scenario][treeFileNum].append(float("nan"))
													else:
														results["memory"][methodName][option][j][scenario][treeFileNum].append(float(line.split()[3]))
											else:
												results["time"][methodName][option][j][scenario][treeFileNum].append(float("nan"))
												results["memory"][methodName][option][j][scenario][treeFileNum].append(float("nan"))
												results["maxMemory"][methodName][option][j][scenario][treeFileNum].append(float("nan"))
												notFound[1][treeFileNum]+="-"+str(i+1)
												notFoundCount[1][treeFileNum]+=1
											file.close()
											
											#collect LKs
											mapleLKoptionNames=["_MapleLK","_MapleLK_unrest"]
											if scenario==0 or scenario==2 or scenario==3:
												rateVarStrLK=["","_rateVar"]
											else:
												rateVarStrLK=[""]
											for mapleLKoptionIndex in range(len(mapleLKoptionNames)):
												for rateVarStrNum in range(len(rateVarStrLK)):
													if not os.path.isfile(fileNamePre+rateVarStr[treeFileNum]+mapleLKoptionNames[mapleLKoptionIndex]+rateVarStr[rateVarStrNum]+"_LK.txt"):
														notFound[2+mapleLKoptionIndex+2*rateVarStrNum][treeFileNum]+="-"+str(i+1)
														notFoundCount[2+mapleLKoptionIndex+2*rateVarStrNum][treeFileNum]+=1
														results[valueNames[3+mapleLKoptionIndex+2*rateVarStrNum]][methodName][option][j][scenario][treeFileNum].append(float('nan'))
													else:
														file=open(fileNamePre+rateVarStr[treeFileNum]+mapleLKoptionNames[mapleLKoptionIndex]+rateVarStr[rateVarStrNum]+"_LK.txt")
														line=file.readline()
														results[valueNames[3+mapleLKoptionIndex+2*rateVarStrNum]][methodName][option][j][scenario][treeFileNum].append(float(line.split()[0]))
														#lk=float(line.split()[0])
														#lks[versNum][j][scenario][option][paramNum].append(lk)
														file.close()

											#collect IQtree LKs
											for rateVarStrNum in range(len(rateVarStrLK)):
												if not os.path.isfile(fileNamePre+rateVarStr[treeFileNum]+"_IQtreeLK"+rateVarStr[rateVarStrNum]+".log"):
													notFound[6+rateVarStrNum][treeFileNum]+="-"+str(i+1)
													notFoundCount[6+rateVarStrNum][treeFileNum]+=1
													results[valueNames[7+rateVarStrNum]][methodName][option][j][scenario][treeFileNum].append(float('nan'))
													#notFoundIQLK+="-"+str(i+1)
													#notFoundIQLKCount+=1
												else:
													file=open(fileNamePre+rateVarStr[treeFileNum]+"_IQtreeLK"+rateVarStr[rateVarStrNum]+".log")
													line=file.readline()
													while not ("BEST SCORE FOUND :" in line) and line!="":
														line=file.readline()
													if line=="":
														notFound[6+rateVarStrNum][treeFileNum]+="-"+str(i+1)
														notFoundCount[6+rateVarStrNum][treeFileNum]+=1
														results[valueNames[7+rateVarStrNum]][methodName][option][j][scenario][treeFileNum].append(float('nan'))
														#notFoundIQLK+="-"+str(i+1)
														#notFoundIQLKCount+=1
													else:
														#IQlk=float(line.split()[4])
														#IQlks[versNum][j][scenario][option][paramNum].append(IQlk)
														results[valueNames[7+rateVarStrNum]][methodName][option][j][scenario][treeFileNum].append(float(line.split()[4]))
													file.close()

											#collect RF distances
											if scenario>0:
												if not os.path.isfile(fileNamePre+rateVarStr[treeFileNum]+"_RF_RFdistances.txt"):
													#notFoundRF+="-"+str(i+1)
													#notFoundRFCount+=1
													notFound[8][treeFileNum]+="-"+str(i+1)
													notFoundCount[8][treeFileNum]+=1
													results[valueNames[9]][methodName][option][j][scenario][treeFileNum].append(float('nan'))
												else:
													file=open(fileNamePre+rateVarStr[treeFileNum]+"_RF_RFdistances.txt")
													line=file.readline()
													line=file.readline()
													if len(line.split())<2 or line.split()[1]=='None':
														print(fileNamePre+rateVarStr[treeFileNum]+"_RF_RFdistances.txt")
														print(line)
														#notFoundRF+="-"+str(i+1)
														#notFoundRFCount+=1
														notFound[8][treeFileNum]+="-"+str(i+1)
														notFoundCount[8][treeFileNum]+=1
														results[valueNames[9]][methodName][option][j][scenario][treeFileNum].append(float('nan'))
													else:
														results[valueNames[9]][methodName][option][j][scenario][treeFileNum].append(float(line.split()[1]))
														#RFdistance=float(line.split()[1])
														#RFs[versNum][j][scenario][option][paramNum].append(RFdistance)
													file.close()

											#valueNames=["time","memory","maxMemory","LK","LKunrest","LKsiteVar","LKunrestSiteVar","IQLK","IQLKsiteVar","RF","parsimony"]
											if not os.path.isfile(fileNamePre+rateVarStr[treeFileNum]+"_parsimony_console_error.txt"):
												notFound[9][treeFileNum]+="-"+str(i+1)
												notFoundCount[9][treeFileNum]+=1
												results[valueNames[10]][methodName][option][j][scenario][treeFileNum].append(float('nan'))
											else:
												file=open(fileNamePre+rateVarStr[treeFileNum]+"_parsimony_console_error.txt")
												line=file.readline()
												while not ("The parsimony score for this tree is:" in line) and line!="":
													line=file.readline()
												if line=="":
													notFound[9][treeFileNum]+="-"+str(i+1)
													notFoundCount[9][treeFileNum]+=1
													results[valueNames[10]][methodName][option][j][scenario][treeFileNum].append(float('nan'))
												else:
													results[valueNames[10]][methodName][option][j][scenario][treeFileNum].append(float(line.split()[7]))
												file.close()

								fileTypes=["tree","time-memory","LK","LKunrest","LKsiteVar","LKunrestSiteVar","IQLK","IQLKsiteVar","RF","parsimony"]
								for treeFileNum in range(2):
									if notFoundCount[0][treeFileNum]<10:
										for fileTypeNum in range(len(fileTypes)):
											if notFoundCount[fileTypeNum][treeFileNum]>0:
												print("Scenario "+folder+" "+methodName+methodOptionsNames[method][option]+rateVarStr[treeFileNum]+", "+str(j)+" samples, "+fileTypes[fileTypeNum]+" file not found for replicates "+notFound[fileTypeNum][treeFileNum])
									
								#for treeFileNum in range(2):
										if len(results[valueNames[0]][methodName][option][j][scenario][treeFileNum])>0:
											print("\n Scenario "+folder+" "+methodName+methodOptionsNames[method][option]+rateVarStr[treeFileNum]+", "+str(j)+" samples : times, memory, max memory, LKs, unrestLK, siteVarLK, siteVar-unrestLK, IQlk, siteVarIQlk  LKs, RF, parsimony")
											for valueNum in range(len(valueNames)):
												print(valueNames[valueNum])
												try:
													printMinMeanMax(results[valueNames[valueNum]][methodName][option][j][scenario][treeFileNum])
												except:
													print(results[valueNames[valueNum]][methodName][option][j][scenario][treeFileNum])
											print("\n")
									else:
										print("no results files found for scenario "+folder+" "+methodName+methodOptionsNames[method][option]+rateVarStr[treeFileNum]+", "+str(j)+" samples \n")

								for treeFileNum in range(2):
									for valueNum in range(len(valueNames)):
										for val in results[valueNames[valueNum]][methodName][option][j][scenario][treeFileNum]:
											outputFiles[valueNum][treeFileNum].write(str(val)+"\t")
										outputFiles[valueNum][treeFileNum].close()
									#outputFiles.append([open(args.pathToSimulationFolder+folder+'/'+str(j)+"subsamples/"+methodName+methodOptionsNames[method][option]+"_"+res+"_summary.txt","w"),open(args.pathToSimulationFolder+folder+'/'+str(j)+"subsamples/"+methodName+methodOptionsNames[method][option]+"_rateVar_"+res+"_summary.txt","w")])

						






#collect data for generating the figures
if args.createFigures:
	if args.runExtraAnalysesDivergence:
		ref = collectReference(args.inputSimulationReference)
		scenario=args.scenario
		subsampleTreeInference=args.subSampleNum
		divergences=["01X","02X","05X","1X","2X","5X","10X","20X","50X","100X","200X","500X","1000X"]
		folder="divergence"
		seed=args.repeat
		repeat=seed

		#collect results for different methods
		valueNames=["time","memory","maxMemory","LK","IQLK","RF","parsimony"]
		results={}
		for res in valueNames:
			results[res]={}
		methods=["UShER","matOptimize","IQtree","FastTree","RAxML-NG",'Maple']
		if args.methodToFocusOn==-1:
			methodsToFocusOn=methods
		else:
			methodsToFocusOn=[methods[args.methodToFocusOn]]
		methodOptions=[[""],[""],[" -fast"],[" -fastest"],[" --blmin 0.000005 --tree pars\{1\}"],[""]]
		methodOptionsNames=[[""],[""],["_fast"],["_fastest"],["_fast"],[""]]

		folders=[folder]
		leavesRange=[2000]
		nLeaves=leavesRange
		figureFolder=args.pathToSimulationFolder+"Figures/"
		if not os.path.isdir(figureFolder):
			os.mkdir(figureFolder)

		#create data files for all other figures
		for value in range(len(valueNames)):
			res=valueNames[value]
			fileResults=open(figureFolder+folder+res+"_results.txt","w")
			for i in range(len(divergences)):
				for method in range(len(methods)):
					methodName=methods[method]
					for option in range(len(methodOptionsNames[method])):
						fileResults.write(folder+" "+str(nLeaves[0])+" samples, divergence "+divergences[i]+" "+methodName+" "+methodOptionsNames[method][option]+" : "+res+" \n")
						if os.path.isfile(args.pathToSimulationFolder+folder+'/'+divergences[i]+"/"+methodName+methodOptionsNames[method][option]+"_"+str(leavesRange[0])+"samples_"+res+"_summary.txt"):
							#file=open(args.pathToSimulationFolder+folder+'/'+str(nLeaves[i])+"subsamples/"+methodName+methodOptionsNames[method][option]+"_"+res+"_summary.txt")
							file=open(args.pathToSimulationFolder+folder+'/'+divergences[i]+"/"+methodName+methodOptionsNames[method][option]+"_"+str(leavesRange[0])+"samples_"+res+"_summary.txt")
							line=file.readline()
							file.close()
							fileResults.write(line+"\n")
						else:
							fileResults.write((("nan"+"\t")*10)+"\n")
			fileResults.close()



	else:
		folders=["realDataSubsamples","simulationsSubsamples","simulations4catSubsamples","simulationsAlphaSubsamples","simulationsNsSubsamples"]
		leavesRange=[1000, 2000, 5000, 10000, 20000, 50000, 100000, 200000, 500000]
		nLeaves=leavesRange
		methodOptions=[[""],[""],[" -pers 0.1 -nstop 500 -blmin 0.000000005"," -blmin 0.000000005"," -fast"],[" -spr 4 -mlacc 2 -slownni"," -spr 4",""," -fastest"],[""," -D"],[" --blmin 0.000000005 --tree pars\{3\}"," --blmin 0.000005 --tree pars\{1\}"],[""," --fast"," --rateVariation"," --model UNREST"," --model UNREST --rateVariation"]]
		methodOptionsNames=[[""],[""],["_slow","_medium","_fast"],["_slow","_medium","_fast","_fastest"],["_slow","_fast"],["_slow","_fast"],["","_fast","_rateVar","_unrest","_unrest_rateVar"]]
		maxNumSamples=[[200000],[200000],[5000,10000,20000],[20000,20000,20000,20000],[5000,5000],[5000,5000],[500000,500000,500000,500000,500000]]
		methods=["UShER","matOptimize","IQtree","FastTree","RAxML","RAxML-NG",'Maple']
		valueNames=["time","memory","maxMemory","LK","LKunrest","LKsiteVar","LKunrestSiteVar","IQLK","IQLKsiteVar","RF","parsimony"]
		rateVars=[[""],[""],["","_rateVar"],["","_rateVar"],["","_rateVar"],["","_rateVar"],[""]]

		figureFolder=args.pathToSimulationFolder+"Figures/"
		if not os.path.isdir(figureFolder):
			os.mkdir(figureFolder)

		#create data file for the plot size figure
		createPlotFileSizes=False
		if createPlotFileSizes:
			#compare file sizes
			folder="realDataSubsamples"
			formats=["Fasta","VCF","MAPLE"]
			sizesDiff={}
			sizesVCF={}
			sizesFasta={}
			referenceSize=os.stat(args.inputRealDataReference).st_size
			n=0
			for j in leavesRange:
				sizesDiff[j]=[]
				sizesVCF[j]=[]
				sizesFasta[j]=[]
				for i in range(10):
					pathToFile=args.pathToSimulationFolder+folder+'/'+str(j)+"subsamples/repeat"+str(i+1)+"_"+str(j)+"samples_"+folder
					sizesVCF[j].append(os.stat(pathToFile+".vcf").st_size + referenceSize)
					sizesFasta[j].append(os.stat(pathToFile+".fa").st_size)
					sizesDiff[j].append(os.stat(pathToFile+".txt").st_size + referenceSize)
				n+=1

			fileResults=open(figureFolder+"comparisonResuts_sizes.txt","w")
			for j in leavesRange:
				fileResults.write("fasta, "+str(j)+" leaves\n")
				for k in range(len(sizesFasta[j])):
					fileResults.write(str(sizesFasta[j][k])+"\t")
				fileResults.write("\n")
			for j in leavesRange:
				fileResults.write("VCF, "+str(j)+" leaves\n")
				for k in range(len(sizesVCF[j])):
					fileResults.write(str(sizesVCF[j][k])+"\t")
				fileResults.write("\n")
			for j in leavesRange:
				fileResults.write("diffFile, "+str(j)+" leaves\n")
				for k in range(len(sizesDiff[j])):
					fileResults.write(str(sizesDiff[j][k])+"\t")
				fileResults.write("\n")
			fileResults.close()

		#create data files for all other figures
		for scenario in range(len(folders)):
			folder=folders[scenario]
			for value in range(len(valueNames)):
				res=valueNames[value]
				fileResults=open(figureFolder+folder+res+"_results.txt","w")
				for i in range(len(nLeaves)):
					for method in range(len(methods)):
						methodName=methods[method]
						rateVarStr=rateVars[method]
						for option in range(len(methodOptionsNames[method])):
							if nLeaves[i]<= maxNumSamples[method][option]:
								for rateVar in range(len(rateVarStr)):
									fileResults.write(folder+" "+str(nLeaves[i])+" samples, "+methodName+" "+methodOptionsNames[method][option]+rateVarStr[rateVar]+" : "+res+" \n")
									if os.path.isfile(args.pathToSimulationFolder+folder+'/'+str(nLeaves[i])+"subsamples/"+methodName+methodOptionsNames[method][option]+rateVarStr[rateVar]+"_"+res+"_summary.txt"):
										file=open(args.pathToSimulationFolder+folder+'/'+str(nLeaves[i])+"subsamples/"+methodName+methodOptionsNames[method][option]+rateVarStr[rateVar]+"_"+res+"_summary.txt")
										line=file.readline()
										file.close()
										fileResults.write(line+"\n")
									else:
										fileResults.write((("nan"+"\t")*10)+"\n")
				fileResults.close()





#Actually create the figures.
#this has to be executed somewhere where there is matpolotlib installed - so not on the cluster
if args.runFigureGeneration:
	if args.runExtraAnalysesDivergence:
		scenario=args.scenario
		subsampleTreeInference=args.subSampleNum
		divergences=["01X","02X","05X","1X","2X","5X","10X","20X","50X","100X","200X","500X","1000X"]
		divergencesText=["0.1","0.2","0.5","1","2","5","10","20","50","100","200","500","1000"]
		divergencesNumbers=[0.1,0.2,0.5,1.0,2.0,5.0,10.0,20.0,50.0,100.0,200.0,500.0,1000.0]
		folder="divergence"
		seed=args.repeat
		repeat=seed

		#collect results for different methods
		valueNames=["time","memory","maxMemory","LK","IQLK","RF","parsimony"]
		results={}
		for res in valueNames:
			results[res]={}
		methods=["UShER","matOptimize","IQtree 2","FastTree 2","RAxML-NG",'Maple']
		if args.methodToFocusOn==-1:
			methodsToFocusOn=methods
		else:
			methodsToFocusOn=[methods[args.methodToFocusOn]]
		methodOptions=[[""],[""],[" -fast"],[" -fastest"],[" --blmin 0.000005 --tree pars\{1\}"],[""]]
		#methodOptionsNames=[[""],[""],["_fast"],["_fastest"],["_fast"],[""]]
		methodOptionsNames=[[""],[""],[""],[""],[""],[""]]

		folders=[folder]
		leavesRange=[2000]
		nLeaves=leavesRange
		figureFolder=args.pathToSimulationFolder+"Figures/"

		
		minV=-0.001
		minY=0.01
		maxY=8000
		logYlist=[True,True,True,False,False,False,False]
		compareToMaxList=[False,False,False,True,True,False,True]

		#TODO hardcoded for now
		figureFolder=folderName+"/FiguresNew/"
		summaryfolder=figureFolder+"resultsSummaries/"

		import numpy as np
		import matplotlib.pyplot as plt
		import matplotlib.patches as mpatches
		import matplotlib.ticker as ticker

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

		def errplot(times,labels,plotFileName,n_leaves,colors,topPlot,linestyles,title="Comparison of Tree Estimation Running Times",yAxisLabel="Time (seconds)",logY=True,ymin=None,ymax=None,violin=True,legendSize=14,tickSize=12,labelSize=18,drawLegend=True,legendLoc='upper left',valuesXticks=None,baseLocator=None,linewidth=3):
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
					#print([float(t) for t in t2])
					#print(np.mean([float(t) for t in t2]))
					#print(np.std([float(t) for t in t2]))
					#print("\n")

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

			#ax1.yaxis.set_major_formatter(ticker.FormatStrFormatter('%0.0001f'))
			
			#reordered_numbers = [0, 1, 2, 4, 6, 5, 3]
			for i in range(len(labels)):
				#ax1.errorbar(x, mean_times[i, :], yerr=errors[i, :], marker='x', c=colors[i], label=labels[i], capsize=2)
				if violin:
					#print(values[i])
					#print(x[:len(values[i])])
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
							violin_parts = ax1.violinplot(valuesPlot, positions=positionsPlot, vert=True, widths=widsPlot, showmeans=False, showextrema=False, showmedians=False, quantiles=None, points=100, bw_method=None, data=None)
						
							for pc in violin_parts['bodies']:
								pc.set_facecolor(colors[i])
								pc.set_edgecolor(colors[i])
							#ax1.errorbar(x, mean_times[i, :], yerr=errors[i, :], marker='.', c=colors[i], label=labels[i], capsize=2, linestyle=linestyles[i])
							if linestyles!=None:
								ax1.plot(x, mean_times[i, :], marker='.', c=colors[i], label=labels[i], linestyle=linestyles[i],linewidth=linewidth)
							else:
								ax1.plot(x, mean_times[i, :], marker='.', c=colors[i], label=labels[i],linewidth=linewidth)
				else:
					if linestyles!=None:
						ax1.errorbar(x, mean_times[i, :], yerr=errors[i, :], marker='x', c=colors[i], label=labels[i], capsize=2, linestyle=linestyles[i],linewidth=linewidth)
					else:
						ax1.errorbar(x, mean_times[i, :], yerr=errors[i, :], marker='x', c=colors[i], label=labels[i], capsize=2,linewidth=linewidth)
			#if violin:
			#	plt.legend(patches, labels, loc='upper left')
			#else:
			if drawLegend:
				# reordering the labels
				handles, labels = plt.gca().get_legend_handles_labels()
				handles.reverse()
				labels.reverse()
				leg = plt.legend(handles, labels,loc=legendLoc, prop={'size': legendSize})
				for i in range(len(leg.get_lines())):
					leg.get_lines()[i].set_linewidth(linewidth)
		
			if valuesXticks!=None:
				plt.xticks(x,valuesXticks[:len(x)],fontsize=tickSize)
			fig.savefig(plotFileName,bbox_inches='tight')
			plt.close('all')

		#now make all other plots
		plotTypes=valueNames
		titles=valueNames
		yAxisLabels=["Time (seconds)","Average memory (Mb)","Maximum memory (Mb)","Topology likelihood difference (Maple)","Topology LK difference (IQtree)","Robinson-Foulds distance from truth","Negative parsimony score difference"]
		nans=[]
		for i in range(10):
			nans.append(float("NaN"))
		speeds=methodOptionsNames
		colors=["black","black","red","green","purple","blue"]
		linestyles=[   (0, ())   ,  (0, (1, 1))   ,   (0, (1, 10))    ,    (5, (10, 3))   ,   (0, (5, 10))    ,   (0, (5, 5))    ,    (0, (5, 1))    ,   (0, (3, 10, 1, 10))    ,    (0, (3, 5, 1, 5, 1, 5))  ]
		
		#create lists of method names etc for the plotting
		namesAll=[]
		colorsAll=[]
		linestylesAll=[]
		for method in range(len(methods)):
			for option in range(len(methodOptionsNames[method])):
				namesAll.append(methods[method]+methodOptionsNames[method][option])
				colorsAll.append(colors[method])
				lineStyle=linestyles[method]
				linestylesAll.append(lineStyle)

		#collect all the data values
		for value in range(len(valueNames)):
			res=valueNames[value]
			fileResults=open(figureFolder+"resultsSummaries/"+folder+res+"_results.txt")
			#first collect the maximum value for each replicate - this is useful for plots where one compares values to the maximum achieved in that replicate
			maxs=[]
			for i in range(len(divergences)):
				maxs.append([])
				for rep in range(10):
					maxs[i].append(float("-inf"))
				for method in range(len(methods)):
					for option in range(len(methodOptionsNames[method])):
						line=fileResults.readline()
						line=fileResults.readline()
						linelist=line.split()
						for l in range(len(linelist)):
							if value==6:
								linelist[l]=-float(linelist[l])
							else:
								linelist[l]=float(linelist[l])
							if linelist[l]>maxs[i][l]:
								maxs[i][l]=linelist[l]
			#create lists of values to plot
			fileResults.close()
			fileResults=open(figureFolder+"resultsSummaries/"+folder+res+"_results.txt")
			plotValues=[]
			for i in range(len(divergences)):
				plotValues.append([])
				for method in range(len(methods)):
					for option in range(len(methodOptionsNames[method])):
						line=fileResults.readline()
						line=fileResults.readline()
						linelist=line.split()
						valueList=[]
						for l in range(len(linelist)):
							if linelist[l]!="nan":
								if method!=1 or i<=10:
									if value==6:
										linelist[l]="-"+linelist[l]
									if compareToMaxList[value]:
										valueList.append(-(maxs[i][l]-float(linelist[l])))
									else:
										if value==0 and method==1:
											if len(plotValues[i][0])>0:
												pos=len(valueList)
												if pos>=len(plotValues[i][0]):
													pos=-1
												valueList.append(float(linelist[l])+plotValues[i][0][pos])
										else:
											valueList.append(float(linelist[l]))
						plotValues[i].append(valueList)
			fileResults.close()

			#create plot with all methods/options crammed together
			#valueNames=["time","memory","maxMemory","LK","LKunrest","LKsiteVar","LKunrestSiteVar","IQLK","IQLKsiteVar","RF","parsimony"]
			#valueNames=["time","memory","maxMemory","LK","IQLK","RF","parsimony"]
			legendLoc='lower right'
			if value==5:
				legendLoc='upper right'
			elif value>2 and value<=6:
				legendLoc='lower left'
			legendSize=10
			tickSize=22
			labelSize=27
			if value>0:
				#TODO for the manuscript might want to turn this off
				#drawLegend=False
				drawLegend=False
			else:
				drawLegend=True
			errplot(plotValues,namesAll,figureFolder+folder+plotTypes[value]+"_all.pdf",divergencesNumbers,colorsAll,'Divergence level',linestylesAll,title='',yAxisLabel=yAxisLabels[value],legendSize=27,tickSize=tickSize,labelSize=labelSize,logY=logYlist[value],drawLegend=drawLegend,legendLoc=legendLoc,valuesXticks=divergencesText,linewidth=2)


	else:
		minV=-0.001
		minY=0.01
		maxY=8000

		leavesRange=[1000, 2000, 5000, 10000, 20000, 50000, 100000, 200000, 500000]
		nLeaves=leavesRange
		nLeavesText=["1000","2000","5000","10000","20000","50000","100000","200000","500000"]
		nLeavesTextTicks=['1000','2000','5000',"10⁴","2x10⁴","5x10⁴","10⁵","2x10⁵","5x10⁵"]

		numSamplesIQtreeLK=[1000, 2000, 5000, 10000, 20000]

		folders=["realDataSubsamples","simulationsSubsamples","simulations4catSubsamples","simulationsAlphaSubsamples","simulationsNsSubsamples"]
		methodOptions=[[""],[""],[" -pers 0.1 -nstop 500 -blmin 0.000000005"," -blmin 0.000000005"," -fast"],[" -spr 4 -mlacc 2 -slownni"," -spr 4",""," -fastest"],[""," -D"],[" --blmin 0.000000005 --tree pars\{3\}"," --blmin 0.000005 --tree pars\{1\}"],[""," --fast"," --rateVariation"," --model UNREST"," --model UNREST --rateVariation"]]
		methodOptionsNames=[[""],[""],["_slow","_medium","_fast"],["_slow","_medium","_fast","_fastest"],["_slow","_fast"],["_slow","_fast"],["","_fast","_rateVar","_unrest","_unrest_rateVar"]]
		maxNumSamples=[[200000],[200000],[5000,10000,20000],[20000,20000,20000,20000],[5000,5000],[5000,5000],[500000,500000,500000,500000,500000]]
		methods=["UShER","matOptimize","IQtree","FastTree","RAxML","RAxML-NG",'Maple']
		valueNames=["time","memory","maxMemory","LK","LKunrest","LKsiteVar","LKunrestSiteVar","IQLK","IQLKsiteVar","RF","parsimony"]
		logYlist=[True,True,True,False,False,False,False,False,False,False,False]
		compareToMaxList=[False,False,False,True,True,True,True,True,True,False,True]
		rateVars=[[""],[""],["","_rateVar"],["","_rateVar"],["","_rateVar"],["","_rateVar"],[""]]

		figureFolder=folderName+"/FiguresNew/"
		summaryfolder=figureFolder+"resultsSummaries/"

		import numpy as np
		import matplotlib.pyplot as plt
		import matplotlib.patches as mpatches
		import matplotlib.ticker as ticker

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

		def errplot(times,labels,plotFileName,n_leaves,colors,topPlot,linestyles,title="Comparison of Tree Estimation Running Times",yAxisLabel="Time (seconds)",logY=True,ymin=None,ymax=None,violin=True,legendSize=14,tickSize=12,labelSize=18,drawLegend=True,legendLoc='upper left',valuesXticks=None,baseLocator=None,linewidth=3):
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
					#print([float(t) for t in t2])
					#print(np.mean([float(t) for t in t2]))
					#print(np.std([float(t) for t in t2]))
					#print("\n")

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

			#ax1.yaxis.set_major_formatter(ticker.FormatStrFormatter('%0.0001f'))
			
			#reordered_numbers = [0, 1, 2, 4, 6, 5, 3]
			for i in range(len(labels)):
				#ax1.errorbar(x, mean_times[i, :], yerr=errors[i, :], marker='x', c=colors[i], label=labels[i], capsize=2)
				if violin:
					#print(values[i])
					#print(x[:len(values[i])])
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
							violin_parts = ax1.violinplot(valuesPlot, positions=positionsPlot, vert=True, widths=widsPlot, showmeans=False, showextrema=False, showmedians=False, quantiles=None, points=100, bw_method=None, data=None)
						
							for pc in violin_parts['bodies']:
								pc.set_facecolor(colors[i])
								pc.set_edgecolor(colors[i])
							#ax1.errorbar(x, mean_times[i, :], yerr=errors[i, :], marker='.', c=colors[i], label=labels[i], capsize=2, linestyle=linestyles[i])
							if linestyles!=None:
								ax1.plot(x, mean_times[i, :], marker='.', c=colors[i], label=labels[i], linestyle=linestyles[i],linewidth=linewidth)
							else:
								ax1.plot(x, mean_times[i, :], marker='.', c=colors[i], label=labels[i],linewidth=linewidth)
				else:
					if linestyles!=None:
						ax1.errorbar(x, mean_times[i, :], yerr=errors[i, :], marker='x', c=colors[i], label=labels[i], capsize=2, linestyle=linestyles[i],linewidth=linewidth)
					else:
						ax1.errorbar(x, mean_times[i, :], yerr=errors[i, :], marker='x', c=colors[i], label=labels[i], capsize=2,linewidth=linewidth)
			#if violin:
			#	plt.legend(patches, labels, loc='upper left')
			#else:
			if drawLegend:
				# reordering the labels
				handles, labels = plt.gca().get_legend_handles_labels()
				handles.reverse()
				labels.reverse()
				leg = plt.legend(handles, labels,loc=legendLoc, prop={'size': legendSize})
				for i in range(len(leg.get_lines())):
					leg.get_lines()[i].set_linewidth(linewidth)
		
			if valuesXticks!=None:
				plt.xticks(x,valuesXticks[:len(x)],fontsize=tickSize)
			fig.savefig(plotFileName,bbox_inches='tight')
			plt.close('all')

		createSizePlot=True
		if createSizePlot:
			#make plot of file sizes
			folder="realDataSubsamples"
			formats=["Fasta","VCF","MAPLE"]
			fileName=summaryfolder+"comparisonResuts_sizes.txt"
			title="Comparison of file sizes for different alignment formats"
			yAxisLabel="Size of file in bytes"
			colors=["green","black","blue"]
			file=open(fileName)
			values=[]
			nLeaves=leavesRange
			for m in range(len(formats)):
				values.append([])
				leavesRange=range(len(nLeaves))
				for n in leavesRange:
					values[m].append([])
					#values[m][n].append([])
					line=file.readline()
					line=file.readline()
					linelist=line.split()
					#if len(linelist)>4:
					for v in linelist:
						if not np.isnan(float(v)):
							values[m][n].append(float(v))
			file.close()
			plotValues=[]
			for n in range(len(nLeaves)):
				plotValues.append([])
				for m in range(len(formats)):
					plotValues[n].append(values[m][n])

			errplot(plotValues[1:],formats,figureFolder+"errplot_sizes.pdf",nLeaves[1:],colors,'Number of tips',None,title=title,yAxisLabel=yAxisLabel,legendSize=27,tickSize=22,labelSize=27,legendLoc='upper left',valuesXticks=['2000','5000',"10⁴","2x10⁴","5x10⁴","10⁵","2x10⁵","5x10⁵"])

		#now make all other plots
		plotTypes=valueNames
		titles=valueNames
		yAxisLabels=["Time in seconds","Average memory (Mb)","Maximum memory (Mb)","Topology likelihood difference (Maple)","Topology LK diff (Maple, UNREST)","Topology LK diff (MAPLE, site var)"," Topology LK diff (MAPLE, UNREST, site var)","Topology LK difference (IQtree)","Topology LK difference (IQtree, site variation)","Robinson-Foulds distance from truth","Negative parsimony score difference"]
		nans=[]
		for i in range(10):
			nans.append(float("NaN"))
		speeds=methodOptionsNames
		colors=["black","black","red","green","purple","purple","blue"]
		linestyles=[   (0, ())   ,  (0, (1, 1))   ,   (0, (1, 10))    ,    (5, (10, 3))   ,   (0, (5, 10))    ,   (0, (5, 5))    ,    (0, (5, 1))    ,   (0, (3, 10, 1, 10))    ,    (0, (3, 5, 1, 5))    ,    (0, (3, 5, 1, 5, 1, 5))  ]
		indecesVeryShort=[[[False]],[[True]],[[False,False],[False,False],[True,False]],[[False,False],[False,False],[False,False],[True,False]],[[False,False],[False,False]],[[False,False],[True,False]],[[True],[False],[False],[False],[False]]]
		indecesRateVar=[[[False]],[[True]],[[False,False],[False,False],[True,True]],[[False,False],[False,False],[False,False],[True,True]],[[False,False],[False,False]],[[False,False],[True,True]],[[True],[False],[True],[True],[True]]]
		
		#create lists of method names etc for the plotting
		namesAll=[]
		colorsAll=[]
		linestylesAll=[]
		namesVeryShort=[]
		colorsVeryShort=[]
		linestylesVeryShort=[]
		namesRateVar=[]
		colorsRateVar=[]
		linestylesRateVar=[]
		for method in range(len(methods)):
			for option in range(len(methodOptionsNames[method])):
				for rateVar in range(len(rateVars[method])):
					namesAll.append(methods[method]+methodOptionsNames[method][option]+rateVars[method][rateVar])
					colorsAll.append(colors[method])
					if method==1:
						lineStyle=linestyles[1+option+rateVar*len(methodOptionsNames[method])]
					elif method==5:
						lineStyle=linestyles[4+option+rateVar*len(methodOptionsNames[method])]
					else:
						lineStyle=linestyles[option+rateVar*len(methodOptionsNames[method])]
					linestylesAll.append(lineStyle)
					if indecesVeryShort[method][option][rateVar]:
						namesVeryShort.append(methods[method]+methodOptionsNames[method][option]+rateVars[method][rateVar])
						colorsVeryShort.append(colors[method])
						linestylesVeryShort.append(lineStyle)
					if indecesRateVar[method][option][rateVar]:
						namesRateVar.append(methods[method]+methodOptionsNames[method][option]+rateVars[method][rateVar])
						colorsRateVar.append(colors[method])
						linestylesRateVar.append(lineStyle)

		#collect all the data values
		for scenario in range(len(folders)):
			folder=folders[scenario]
			for value in range(len(valueNames)):
				res=valueNames[value]
				fileResults=open(figureFolder+"resultsSummaries/"+folder+res+"_results.txt")
				#first collect the maximum value for each replicate - this is useful for plots where one compares values to the maximum achieved in that replicate
				maxs=[]
				maxsVeryShort=[]
				maxsRateVar=[]
				for i in range(len(nLeaves)):
					maxs.append([])
					maxsVeryShort.append([])
					maxsRateVar.append([])
					for rep in range(10):
						maxs[i].append(float("-inf"))
						maxsVeryShort[i].append(float("-inf"))
						maxsRateVar[i].append(float("-inf"))
					for method in range(len(methods)):
						rateVarStr=rateVars[method]
						for option in range(len(methodOptionsNames[method])):
							if nLeaves[i]<= maxNumSamples[method][option]:
								for rateVar in range(len(rateVarStr)):
									line=fileResults.readline()
									line=fileResults.readline()
									linelist=line.split()
									for l in range(len(linelist)):
										if value==10:
											linelist[l]=-float(linelist[l])
										else:
											linelist[l]=float(linelist[l])
										if linelist[l]>maxs[i][l]:
											maxs[i][l]=linelist[l]
										if indecesVeryShort[method][option][rateVar] and linelist[l]>maxsVeryShort[i][l]:
											maxsVeryShort[i][l]=linelist[l]
										if indecesRateVar[method][option][rateVar] and linelist[l]>maxsRateVar[i][l]:
											maxsRateVar[i][l]=linelist[l]
				#create lists of values to plot
				fileResults.close()
				fileResults=open(figureFolder+"resultsSummaries/"+folder+res+"_results.txt")
				plotValues=[]
				plotValuesVeryShort=[]
				plotValuesRateVar=[]
				for i in range(len(nLeaves)):
					plotValues.append([])
					plotValuesVeryShort.append([])
					plotValuesRateVar.append([])
					for method in range(len(methods)):
						rateVarStr=rateVars[method]
						for option in range(len(methodOptionsNames[method])):
							if nLeaves[i]<= maxNumSamples[method][option]:
								for rateVar in range(len(rateVarStr)):
									line=fileResults.readline()
									line=fileResults.readline()
									linelist=line.split()
									valueList=[]
									valueListVeryShort=[]
									valueListRateVar=[]
									for l in range(len(linelist)):
										if linelist[l]!="nan":
											if value==10:
												linelist[l]="-"+linelist[l]
											if compareToMaxList[value]:
												valueList.append(-(maxs[i][l]-float(linelist[l])))
												if indecesVeryShort[method][option][rateVar]:
													valueListVeryShort.append(-(maxsVeryShort[i][l]-float(linelist[l])))
												if indecesRateVar[method][option][rateVar]:
													valueListRateVar.append(-(maxsRateVar[i][l]-float(linelist[l])))
											else:
												if value==0 and method==1:
													pos=len(valueList)
													if pos>=len(plotValues[i][0]):
														pos=-1
													valueList.append(float(linelist[l])+plotValues[i][0][pos])
													if indecesVeryShort[method][option][rateVar]:
														valueListVeryShort.append(float(linelist[l])+plotValues[i][0][pos])
													if indecesRateVar[method][option][rateVar]:
														valueListRateVar.append(float(linelist[l])+plotValues[i][0][pos])
												else:
													valueList.append(float(linelist[l]))
													if indecesVeryShort[method][option][rateVar]:
														valueListVeryShort.append(float(linelist[l]))
													if indecesRateVar[method][option][rateVar]:
														valueListRateVar.append(float(linelist[l]))
									plotValues[i].append(valueList)
									if indecesVeryShort[method][option][rateVar]:
										plotValuesVeryShort[i].append(valueListVeryShort)
									if indecesRateVar[method][option][rateVar]:
										plotValuesRateVar[i].append(valueListRateVar)
							else:
								for rateVar in range(len(rateVarStr)):
									plotValues[i].append([])
									if indecesVeryShort[method][option][rateVar]:
										plotValuesVeryShort[i].append([])
									if indecesRateVar[method][option][rateVar]:
										plotValuesRateVar[i].append([])	
				fileResults.close()

				#create plot with all methods/options crammed together
				#valueNames=["time","memory","maxMemory","LK","LKunrest","LKsiteVar","LKunrestSiteVar","IQLK","IQLKsiteVar","RF","parsimony"]
				legendLoc='lower right'
				if value==9:
					legendLoc='upper right'
				elif value>2 and value<=10:
					legendLoc='lower left'
				if value==7 or value==8:
					maxLeaves=5
				else:
					maxLeaves=100
				legendSize=10
				tickSize=22
				labelSize=27
				errplot(plotValues[1:maxLeaves],namesAll,figureFolder+folder+plotTypes[value]+"_all.pdf",nLeaves[1:maxLeaves],colorsAll,'',linestylesAll,title='',yAxisLabel=yAxisLabels[value],legendSize=legendSize,tickSize=tickSize,labelSize=labelSize,logY=logYlist[value],drawLegend=True,legendLoc=legendLoc,valuesXticks=nLeavesTextTicks[1:maxLeaves],linewidth=2)
				
				#create also plot with only a few methods/options
				doVeryShort=True
				if doVeryShort:
					if value>0:
						#TODO for the manuscript might want to turn this off
						#drawLegend=False
						drawLegend=True
					else:
						drawLegend=True
					errplot(plotValuesVeryShort[1:maxLeaves],namesVeryShort,figureFolder+folder+plotTypes[value]+"_veryFew.pdf",nLeaves[1:maxLeaves],colorsVeryShort,'',linestylesVeryShort,title='',yAxisLabel=yAxisLabels[value],legendSize=27,tickSize=tickSize,labelSize=labelSize,logY=logYlist[value],drawLegend=drawLegend,legendLoc=legendLoc,valuesXticks=nLeavesTextTicks[1:maxLeaves])
			
				#create also plots specifically for unrest and sitevar models
				doRateVar=True
				if doRateVar:
					if value>0:
						#TODO for the manuscript might want to turn this off
						#drawLegend=False
						drawLegend=True
					else:
						drawLegend=True
					errplot(plotValuesRateVar[1:maxLeaves],namesRateVar,figureFolder+folder+plotTypes[value]+"_rateVar.pdf",nLeaves[1:maxLeaves],colorsRateVar,'',linestylesRateVar,title='',yAxisLabel=yAxisLabels[value],legendSize=20,tickSize=tickSize,labelSize=labelSize,logY=logYlist[value],drawLegend=drawLegend,legendLoc=legendLoc,valuesXticks=nLeavesTextTicks[1:maxLeaves])
			












exit()


















