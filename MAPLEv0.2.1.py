import sys
from math import log
import argparse
from time import time
import os.path
from tkinter import E

#©EMBL-European Bioinformatics Institute, 2021

# MAPLE code to estimate a tree by maximum likelihood from a MAPLE format input.

#things in progress
#TODO C++ implementation - parallelization
#TODO Model sequencing errors

#longer-term plans
#TODO codon models, selection
#TODO indels, alignment
#TODO Timed trees
#TODO Bayesian implementation in BEAST
#TODO Recombination
#TODO phylogeography
#TODO phylodynamics/selection

#things one could MAYBE also try in the future:
#TODO This would change the structure of a lot of the code - however one could try to represent internal genomes/likelihoods a mutation-annotated tree instead of representing all nucleotides different from the reference at each internal node.
#TODO branch lenth optimization could maybe be made more efficient by representing the monomer (in t) factors in the derivative all jointly as a polymonial, so to reduce numbers of divisions to be performed (at the cost of increasing numbers of products and sums).
#TODO use slots (https://www.geeksforgeeks.org/slots-in-python/) to reduce memory cost?
#TODO develop more nuanced exploration of topology space: reattempt to re-place only subtrees whose root have been affected in the last round of topology search (which is now already done), 
# and their neighbours (which is not implemented yet)?
#TODO for increased accuracy, one could calculate probVectTot genome lists in a way that would not depend by the choice of which two genome lists to merge first. 
#For now implemented consistent calculation of tot likelihoods so to reduce inconsinstencies met.
#TODO would it be convenient to remove the genome position element from entries of type ACGTO ? Sounds like a lot of changes, and might make things slower and less readable?
#TODO if the model is reversible, no need to consider complicated case of root position, so the genome vector lists can be simpler and I can avoid additional calculations?
#TODO try to replace appendProb() with appendProbNode()? it might be slightly slower, but it would be one fewer function in the code.
#TODO do multiple rounds of SPR search with increasing thresholds? Now performing a quick short-range SPR tree traverse before the proper one.
#TODO create an initial stage for the placement with very low thresholds; then, a second placement stage would use the default thresholds and work 
# like the SPR search but starting the placement search from the node found by the first stage; finally, the last stage would be as usual.
# I tried this but it comes at little computational advantage and at large accuracy cost.


parser = argparse.ArgumentParser(description='Estimate a tree from a diff format and using iterative approximate maximum likelihood sample placement.')
parser.add_argument('--input',default="/Users/demaio/Desktop/GISAID-hCoV-19-phylogeny-2021-03-12/phylogenetic_inference/2021-03-31_unmasked_differences_reduced.txt_consensus-based.txt", help='input MAPLE file name; should contain first the reference genome and then the difference of all samples with respet to the reference.')
parser.add_argument('--reference',default="", help='optional input reference file name. By default it assumes instead that the reference is part of the MAPLE format input.')
parser.add_argument('--output',default="/Users/demaio/Desktop/GISAID-hCoV-19-phylogeny-2021-03-12/phylogenetic_inference/MAPLE", help='output path and identifier to be used for newick output file.')
parser.add_argument('--inputTree',default="", help='input newick tree file name; this is optional, and is used for online inference (or for Robinson-Foulds distance calculation if option --inputRFtrees is also used).')
parser.add_argument("--largeUpdate", help="When using option --inputTree to do online inference, by default the tree is only updated locally where the new sequences are inserted. Use this option instead to perform a thorough update of the phylogeny.", action="store_true")
parser.add_argument('--inputRFtrees',default="", help='file name with input newick trees to be compared with the input tree in option --inputTree to calculate Robinson-Foulds distances; this option will turn off the normal MAPLE estimation - only RF distances will be calculated. Each newick tree contained in this file will be compared to the single tree specified with option --inputTree .')
parser.add_argument("--onlyNambiguities", help="Treat all ambiguities as N (total missing information).", action="store_true")
parser.add_argument("--thresholdProb",help="relative probability threshold used to ignore possible states with very low probabilities.",  type=float, default=0.00000001)
parser.add_argument("--thresholdLogLK",help="logLK difference threshold to consider a logLk close to optimal.",  type=float, default=200.0)
parser.add_argument("--thresholdLogLKtopology",help="logLK difference threshold to consider a logLk close to optimal when looking for topology improvements.",  type=float, default=160.0)
parser.add_argument("--allowedFails",help="Number of times one can go down the tree without increasing placement likelihood before the tree traversal is stopped (only applies to non-0 branch lengths).",  type=int, default=5)
parser.add_argument("--allowedFailsTopology",help="Number of times one can crawl along the tree decreasing placement likelihood before the tree traversal is stopped during topology search (only applies to non-0 branch lengths).",  type=int, default=4)
parser.add_argument("--verbose", help="Print to screen a lot of stuff.", action="store_true")
parser.add_argument("--debugging", help="Test that likelihoods are calculated and updated as expected - time consuming and only meant for small trees for debugging purposes.", action="store_true")
parser.add_argument("--model", help="Which substitution model should be used. Allowed models so far are JC, GTR (default) or UNREST.", default="GTR")
parser.add_argument("--overwrite", help="Overwrite previous results if already present.", action="store_true")
parser.add_argument("--nonBinaryTree", help="Write output tree with multifurcations - by default the tree is written as binary so to avoid problems reading the tree in other software.", action="store_true")
parser.add_argument("--numTopologyImprovements",help="Number of times we traverse the tree looking for topological improvements.",  type=int, default=1)
parser.add_argument("--thresholdTopologyPlacement",help="Don't try to re-place nodes that have current appending logLK cost above this threshold.",  type=float, default=-0.01)
parser.add_argument("--updateSubstMatrixEveryThisSamples",help="How many new samples to place before each update of the substitution rate matrix.",  type=int, default=25)
parser.add_argument("--nonStrictInitialStopRules", help="If specified, then during the initial placement stage, a slower, non-strict rule for stopping the placement search is applied: the search is stopped if enough many consencutive LK worsening are observed, AND if LK is below the considered threshold.", action="store_true")
parser.add_argument("--strictTopologyStopRules", help="If specified, then during the topological improvement stage, a slower, non-strict rule for stopping the SPR search is applied: the search is stopped if enough many consencutive LK worsening are observed, AND if LK is below the considered threshold.", action="store_true")
parser.add_argument("--thresholdDiffForUpdate",help="Consider the probability of a new partial changed if the difference between old and new is above this threshold.",  type=float, default=0.0000001)
parser.add_argument("--thresholdFoldChangeUpdate",help="Consider the probability of a new partial changed, if the fold difference between old and new if above this threshold.",  type=float, default=1.001)
parser.add_argument("--thresholdLogLKconsecutivePlacement",help="logLK difference threshold to consider something as a significant decrease in log-LK when considering consecutive likelihood decreases.",  type=float, default=0.01)
parser.add_argument("--thresholdLogLKwholeTopologyImprovement",help="logLK difference threshold to consider something as a significant decrease in log-LK when considering whole rounds of SPR moves.",  type=float, default=1.0)
parser.add_argument("--calculateLKfinalTree", help="Calculate the log-LK of the final inferred tree.", action="store_true")
parser.add_argument("--fast", help="Set parameters so to run tree inference faster; this will be less accurate in cases of high complexity, for example with recombination, sequencing errors, etc. It will overrule user choices for options --thresholdLogLK , --thresholdLogLKtopology , --allowedFails , --allowedFailsTopology .", action="store_true")
parser.add_argument("--noFastTopologyInitialSearch", help="Don't run a fast short-range topology search before the extensive one.", action="store_true")
parser.add_argument("--noOptimizeBranchLengths", help="Don't run a final round of detailed branch length optimization.", action="store_true")
parser.add_argument("--rateVariation", help="Estimate and use rate variation: the model assumes one rate per site, and the rates are assumed independently (no rate categories). This might cause overfitting if the dataset is not large enough, but in any case one would probably only use MAPLE for large enough datasets.", action="store_true")
parser.add_argument("--minBLenSensitivity",help="Fraction of a mutation to be considered as a precision for branch length estimation (default 0.0001, which means branch lengths estimated up to a 10000th of a mutation precision).",  type=float, default=0.0001)
parser.add_argument("--thresholdLogLKoptimization",help="logLK difference threshold to consider a logLk close to optimal when deciding for which possible placements to perform branch length optimization.",  type=float, default=12.0)
parser.add_argument("--thresholdLogLKoptimizationTopology",help="logLK difference threshold to consider a logLk close to optimal when deciding for which possible placements to perform branch length optimization during the SPR stage.",  type=float, default=12.0)
parser.add_argument('--assignmentFileCSV',default="", help='give path and name to a .csv file containing reference sample names and their lineage assignment. When using this option, option --inputTree should also be used. Then MAPLE will assign all samples in the tree to a lineage following the reference lineages file. Each sample is assigned a lineage same as its closest reference parent.')
parser.add_argument('--assignmentFile',default="", help='Like --assignmentFileCSV but expects an alignment file in any format (as long as sequence names follow a > character as in Fasta and Maple formats).')
parser.add_argument("--defaultBLen",help="Default length of branches, for example when the input tree has no branch length information.",  type=float, default=0.000033)
parser.add_argument("--normalizeInputBLen",help="For the case the input tree has branch lengths expressed not in a likelihood fashion (that is, expected number of substitutions per site), then multiply the input branch lengths by this factor. Particularly useful when using parsimony-based input trees.",  type=float, default=1.0)
parser.add_argument("--multipleInputRFTrees", help="Use this option if the file specified by option --inputRFtrees contains multiple trees - otherwise only read the first tree from that file.", action="store_true")
parser.add_argument("--maxReplacements",help="Maximum number of replacements attempts per node per SPR round (prevents loops).",  type=int, default=10)
parser.add_argument("--useFixedThresholdLogLKoptimizationTopology", help="Use this option if you want to specify the value in --thresholdLogLKoptimizationTopology instead of estimating it from the data.", action="store_true")
parser.add_argument("--writeTreesToFileEveryTheseSteps", help="By default, don't write intermediate trees to file. If however a positive integer is specified with this option, intermediate trees will be written to file every this many topological changes.",  type=int, default=0)
parser.add_argument("--writeLKsToFileEveryTheseSteps", help="By default, don't write likelihoods of intermediate trees to file. If however a positive integer is specified with this option, likelihoods of intermediate trees will be written to file every this many topological changes.",  type=int, default=0)
args = parser.parse_args()

onlyNambiguities=args.onlyNambiguities
thresholdProb=args.thresholdProb
verbose=args.verbose
debugging=args.debugging
inputFile=args.input
outputFile=args.output
refFile=args.reference
allowedFails=args.allowedFails
allowedFailsTopology=args.allowedFailsTopology
model=args.model
thresholdLogLK=args.thresholdLogLK
thresholdLogLKtopology=args.thresholdLogLKtopology
overwrite=args.overwrite
binaryTree=(not args.nonBinaryTree)
numTopologyImprovements=args.numTopologyImprovements
thresholdTopologyPlacement=args.thresholdTopologyPlacement
updateSubstMatrixEveryThisSamples=args.updateSubstMatrixEveryThisSamples
strictInitialStopRules=(not args.nonStrictInitialStopRules)
strictTopologyStopRules=args.strictTopologyStopRules
thresholdDiffForUpdate=args.thresholdDiffForUpdate
thresholdFoldChangeUpdate=args.thresholdFoldChangeUpdate
thresholdLogLKconsecutivePlacement=args.thresholdLogLKconsecutivePlacement
thresholdLogLKwholeTopologyImprovement=args.thresholdLogLKwholeTopologyImprovement
calculateLKfinalTree=args.calculateLKfinalTree

thresholdLogLKoptimization=args.thresholdLogLKoptimization
thresholdLogLKoptimizationTopology=args.thresholdLogLKoptimizationTopology
minBLenSensitivity=args.minBLenSensitivity

fastTopologyInitialSearch=(not args.noFastTopologyInitialSearch)
optimizeBranchLengths=(not args.noOptimizeBranchLengths)
example=False
runFluTest=False
runFast=args.fast
if runFast:
	thresholdLogLK=160.0
	allowedFails=4
	allowedFailsTopology=2
	thresholdLogLKtopology=80.0
	thresholdTopologyPlacement=-1.0
	minBLenSensitivity=0.001

strictTopologyStopRulesInitial=True
allowedFailsTopologyInitial=1
thresholdLogLKtopologyInitial=40.0
thresholdTopologyPlacementInitial=-1.0

defaultBLen=args.defaultBLen
normalizeInputBLen=args.normalizeInputBLen
multipleInputRFTrees=args.multipleInputRFTrees

inputTree=args.inputTree
inputRFtrees=args.inputRFtrees
largeUpdate=args.largeUpdate
rateVariation=args.rateVariation
assignmentFile=args.assignmentFile
assignmentFileCSV=args.assignmentFileCSV
maxReplacements=args.maxReplacements
writeTreesToFileEveryTheseSteps=args.writeTreesToFileEveryTheseSteps
writeLKsToFileEveryTheseSteps=args.writeLKsToFileEveryTheseSteps

warnedBLen=[False]
warnedTotDiv=[False]
sumChildLKs=[0.0]
numChildLKs=[0]
useFixedThresholdLogLKoptimizationTopology=args.useFixedThresholdLogLKoptimizationTopology
totDivFromRef=[0.0]



#Class defininf nodes of the tree
class Tree(object):
	def __init__(self, name='', children=None, dist=defaultBLen):
		if name!='':
			self.name = name
		self.dist = dist
		self.replacements=0
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

		branchLengthsIn=(":" in nwString)
		index=0
		node=Tree()
		#node.replacements=0
		node.dirty=dirtiness
		name=""
		distStr=""
		finished=False
		while index<len(nwString):
			if nwString[index]=="(":
				newNode=Tree()
				#newNode.replacements=0
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
					node.dist=float(distStr)*normalizeInputBLen
					distStr=""
				else:
					node.dist=defaultBLen
				newNode=Tree()
				#newNode.replacements=0
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
					node.dist=float(distStr)*normalizeInputBLen
					distStr=""
				else:
					node.dist=defaultBLen
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
			



#Robinson-Foulds distance (1981) using a simplification of the algorithm from Day 1985.
#this function prepare the data to compare trees to a reference one t1.
#I split in two functions (the second one is RobinsonFouldsWithDay1985() below) so that I don't have to repeat these steps for the reference tree if I compare to the same reference tree multiple times.
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
# Example usage: leafNameDict, nodeTable, leafCount, numBranches = prepareTreeComparison(phyloTrue,rooted=False)
# numDiffs, normalisedRF, leafCount, foundBranches, missedBranches, notFoundBranches = RobinsonFouldsWithDay1985(phyloEstimated,leafNameDict,nodeTable,leafCount,numBranches,rooted=False)
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
					#######
					#del leafNameDict[node.name]
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
		print("There are leaves in the reference that have not been found in this new tree - leafCount "+str(leafCount)+" visitedLeaves "+str(visitedLeaves))
		#print(leafNameDict)
		return None, None, None, None, None, None
	#first value is number of differences, second value is max number of differences just in case one wants the normalized values; 
	#the other values are there just in case on wants more detail.
	numDiffs=((numBranches-foundBranches)+missedBranches)
	return numDiffs, float(numDiffs)/(2*(leafCount-3)), leafCount, foundBranches, missedBranches, (numBranches-foundBranches)

topologyChanges=[0]
if writeTreesToFileEveryTheseSteps>0:
	lineSplit=outputFile.split("/")
	lineSplit[-1]=""
	outFolder="/".join(lineSplit)
	if not os.path.isdir(outFolder):
		print("Path to output file "+outFolder+" does not exist, quitting MAPLE. Use option --output to specify a valid output file path and output file name.")
		exit()
	if os.path.isfile(outputFile+"_intermediateTrees.tree")  and (not overwrite):
		print("File "+outputFile+"_intermediateTrees.tree already exists, quitting MAPLE. Use option --overwrite if you want to overwrite previous inference.")
		exit()
	intermediateTreesFile=open(outputFile+"_intermediateTrees.tree","w")
if writeLKsToFileEveryTheseSteps>0:
	lineSplit=outputFile.split("/")
	lineSplit[-1]=""
	outFolder="/".join(lineSplit)
	if not os.path.isdir(outFolder):
		print("Path to output file "+outFolder+" does not exist, quitting MAPLE. Use option --output to specify a valid output file path and output file name.")
		exit()
	if os.path.isfile(outputFile+"_intermediateLKs.txt")  and (not overwrite):
		print("File "+outputFile+"_intermediateLKs.txt already exists, quitting MAPLE. Use option --overwrite if you want to overwrite previous inference.")
		exit()
	intermediateLKsFile=open(outputFile+"_intermediateLKs.txt","w")


#run Robinson-Foulds distance calculations
if inputRFtrees!="":
	if not os.path.isfile(inputTree):
		print("Input tree in newick format "+inputTree+" not found, quitting MAPLE RF distance calculation. Use option --inputTree to specify a valid input newick tree file.")
		exit()
	if not os.path.isfile(inputRFtrees):
		print("Input trees in newick format "+inputRFtrees+" not found, quitting MAPLE RF distance calculation. Use option --inputRFtrees to specify valid file with input newick trees.")
		exit()
	lineSplit=outputFile.split("/")
	lineSplit[-1]=""
	outFolder="/".join(lineSplit)
	if not os.path.isdir(outFolder):
		print("Path to output file "+outFolder+" does not exist, quitting MAPLE RF calculation. Use option --output to specify a valid output file path and output file name.")
		exit()
	if os.path.isfile(outputFile+"_RFdistances.txt")  and (not overwrite):
		print("File "+outputFile+"_RFdistances.txt already exists, quitting MAPLE RF calculation. Use option --overwrite if you want to overwrite previous inference.")
		exit()
	tree1=readNewick(inputTree)[0]
	print("Read input newick tree")
	leafNameDict, nodeTable, leafCount, numBranches = prepareTreeComparison(tree1,rooted=False)
	otherTrees=readNewick(inputRFtrees,multipleTrees=multipleInputRFTrees)
	print("Read other input newick trees to be compared to the first one")
	file=open(outputFile+"_RFdistances.txt","w")
	file.write("RF\t"+"normalisedRF\t"+"leaves\t"+"foundBranches\t"+"missedBranches\t"+"notFoundBranches\n")
	for tree in otherTrees:
		numDiffs, normalisedRF, leafCount, foundBranches, missedBranches, notFoundBranches = RobinsonFouldsWithDay1985(tree,leafNameDict, nodeTable, leafCount, numBranches,rooted=False)
		file.write(str(numDiffs)+"\t"+str(normalisedRF)+"\t"+str(leafCount)+"\t"+str(foundBranches)+"\t"+str(missedBranches)+"\t"+str(notFoundBranches)+"\n")
	print("Comparison ended")
	file.close()
	exit()

#run lineage assignment of samples in input tree following the input references
if assignmentFile!="" or assignmentFileCSV!="":
	if assignmentFile!="" and assignmentFileCSV!="":
		print("Please only use one between options --assignmentFile and --assignmentFileCSV .")
		exit()
	if not os.path.isfile(inputTree):
		print("Input tree in newick format "+inputTree+" not found, quitting MAPLE lineage assignment. Use option --inputTree to specify a valid input newick tree file.")
		exit()
	if assignmentFile!="" and (not os.path.isfile(assignmentFile)):
		print("Input references assignments file "+assignmentFile+" not found, quitting MAPLE lineage assignment. Use option --assignmentFile to specify a valid csv file in reference sequences and their lineages.")
		exit()
	if assignmentFileCSV!="" and (not os.path.isfile(assignmentFileCSV)):
		print("Input references assignments file "+assignmentFileCSV+" not found, quitting MAPLE lineage assignment. Use option --assignmentFile to specify a valid csv file in reference sequences and their lineages.")
		exit()
	lineSplit=outputFile.split("/")
	lineSplit[-1]=""
	outFolder="/".join(lineSplit)
	if not os.path.isdir(outFolder):
		print("Path to output file "+outFolder+" does not exist, quitting MAPLE lineage assignment. Use option --output to specify a valid output file path and output file name.")
		exit()
	if os.path.isfile(outputFile+"_lineageAssignments.csv")  and (not overwrite):
		print("File "+outputFile+"_lineageAssignments.csv already exists, quitting MAPLE lineage assignment. Use option --overwrite if you want to overwirte previous inference.")
		exit()
	tree1=readNewick(inputTree)[0]
	print("Read input newick tree")
	if assignmentFileCSV!="":
		referencesFile=open(assignmentFileCSV)
		line=referencesFile.readline()
		references={}
		while line!="":
			linelist=line.split(",")
			if len(linelist)==2:
				references[linelist[0]]=linelist[1].replace("\n","")
			line=referencesFile.readline()
	else:
		referencesFile=open(assignmentFile)
		line=referencesFile.readline()
		references={}
		while line!="":
			if line[0]==">":
				name=line.replace("\n","").replace(">","")
				references[name]=name
			line=referencesFile.readline()
	file=open(outputFile+"_lineageAssignments.csv","w")

	node=tree1
	direction=0 #0=from parent, 1+ from children
	lineage=""
	#for each internal node with positive bLen, check all descendants with a distance of 0 from it; if any of them is a reference,
	#then assign that lineage to the node and all its descendants unless the assignment is changed downstream.
	#If no reference at distance 0 is found, then use parent assignment.
	while node!=None:
		#case of an internal node
		if node.children:
			#coming from parent node
			if direction==0:
				if node.dist:
					#first find 0-distance nodes and collect any reference lineage among them
					nodesToVisit=list(node.children)
					while nodesToVisit:
						nextNode=nodesToVisit.pop()
						if not nextNode.dist:
							if nextNode.children:
								for child in nextNode.children:
									nodesToVisit.append(child)
							elif nextNode.name in references:
								lineage=references[nextNode.name]
				node.lineage=lineage
				node=node.children[0]
			else:
				if direction==len(node.children):
					if node.up!=None:
						direction=-1
						for childNum in range(len(node.up.children)):
							if node==node.up.children[childNum]:
								direction=childNum+1
						if direction==-1:
							print("Error, not found child node in parent node list!")
							exit()
					node=node.up
				else:
					lineage=node.lineage
					node=node.children[direction]
					direction=0
		
		#case of a sample
		else:
			if node.name in references:
				file.write(node.name+","+references[node.name]+"\n")
			else:
				file.write(node.name+","+lineage+"\n")
			if node.up!=None:
				direction=-1
				for childNum in range(len(node.up.children)):
					if node==node.up.children[childNum]:
						direction=childNum+1
				if direction==-1:
					print("Error2, not found child node in parent node list!")
					exit()
			node=node.up

	print("Lineage assignment completed")
	file.close()
	exit()

	



minimumCarryOver=sys.float_info.min*(1e50)

if os.path.isfile(outputFile+"_tree.tree")  and (not overwrite):
	print("File "+outputFile+"_tree.tree already exists, quitting MAPLE tree inference. Use option --overwrite if you want to overwirte previous inference.")
	exit()
if not os.path.isfile(inputFile):
	print("Input file in Maple format "+inputFile+" not found, quitting MAPLE tree inference. Use option --input to specify a valid input file.")
	exit()
if refFile!="" and (not os.path.isfile(refFile)):
	print("Input reference fasta file "+refFile+" not found, quitting MAPLE tree inference. Use option --reference to specify a valid input reference file.")
	exit()
lineSplit=outputFile.split("/")
lineSplit[-1]=""
outFolder="/".join(lineSplit)
if not os.path.isdir(outFolder):
	print("Path to output file "+outFolder+" does not exist, quitting MAPLE tree inference. Use option --output to specify a valid output file path and output file name.")
	exit()
if inputTree!="":
	if not os.path.isfile(inputTree):
		print("Input tree in newick format "+inputTree+" not found, quitting MAPLE. Use option --inputTree to specify a valid input newick tree file.")
		exit()
	tree1=readNewick(inputTree,dirtiness=largeUpdate)[0]
	print("Read input newick tree")
	makeTreeBinary(tree1)


alleles={"A":0,"C":1,"G":2,"T":3}
allelesList=["A","C","G","T"]
allelesLow={"a":0,"c":1,"g":2,"t":3}
allelesUpOrLow={"a":0,"c":1,"g":2,"t":3,"A":0,"C":1,"G":2,"T":3}
allelesListLow=["a","c","g","t"]
ambiguities={"y":[0.0,0.5,0.0,0.5],"r":[0.5,0.0,0.5,0.0],"w":[0.5,0.0,0.0,0.5],"s":[0.0,0.5,0.5,0.0],"k":[0.0,0.0,0.5,0.5],"m":[0.5,0.5,0.0,0.0],"d":[1.0/3,0.0,1.0/3,1.0/3],"v":[1.0/3,1.0/3,1.0/3,0.0],"h":[1.0/3,1.0/3,0.0,1.0/3],"b":[0.0,1.0/3,1.0/3,1.0/3]}



#collect reference and calculate vector of cumulative numbers of nucleotides and calculate root frequencies.
def collectReference(fileName):
	file=open(fileName)
	line=file.readline()
	ref=""
	while line!="":
		line=file.readline()
		ref+=line.replace("\n","")
	ref=ref.lower()
	file.close()
	return ref



#Read input file
def readConciseAlignment(fileName,extractReference=True,ref="",extractNames=False):
	fileI=open(fileName)
	line=fileI.readline()
	if extractReference:
		line=fileI.readline()
		ref=""
		while line!="" and line[0]!=">":
			ref+=line.replace("\n","")
			line=fileI.readline()
		ref=ref.lower()
	nSeqs=0
	if extractNames:
		data={}
	else:
		data=[]
	while line!="" and line!="\n":
		seqList=[]
		if extractNames:
			name=line.replace(">","").replace("\n","")
		line=fileI.readline()
		pos=0
		while line!="" and line!="\n" and line[0]!=">":
			linelist=line.split()
			if len(linelist)>2:
				entry=(linelist[0].lower(),int(linelist[1]),int(linelist[2]))
			elif len(linelist)<2:
				print("In input file "+fileName+" found line with only one column: \n"+line+"ERROR Please check for errors in the alignment format; if the reference is included at the top of the alignment, then please don't use option --reference.")
				exit()
			else:
				entry=(linelist[0].lower(),int(linelist[1]))
			if ref[entry[1]-1]==entry[0]:
				print("Mutation observed into reference nucleotide at position "+str(entry[1])+" , nucleotide "+entry[0]+". Wrong reference and/or diff file?")
				exit()
			if entry[1]<=pos:
				print("WARNING, at sample number "+str(nSeqs+1)+" found entry")
				print(line.replace("\n",""))
				print("which is inconsistent since the position is already represented by another entry:")
				print(seqList[-1])
				#print(" ignoring the entry in the file and proceeding with the analysis.")
				exit()
			else:
				seqList.append(entry)
				if len(entry)==2:
					pos=entry[1]
				else:
					pos=entry[1]+entry[2]-1
			line=fileI.readline()
		if extractNames:
			data[name]=seqList
		else:
			data.append(seqList)
		nSeqs+=1
	fileI.close()
	print(str(nSeqs)+" sequences in diff file.")
	if extractReference:
		return ref, data
	else:
		return data

runOnlyExample=False
if not runOnlyExample:	
	#read sequence data from file
	extractNamesFlag=False
	if (inputTree!="") or (writeTreesToFileEveryTheseSteps>0):
		extractNamesFlag=True
	if refFile=="":
		
		ref, data=readConciseAlignment(inputFile,extractNames=extractNamesFlag)
	else:
		ref=collectReference(refFile)
		data=readConciseAlignment(inputFile, extractReference=False, ref=ref,extractNames=extractNamesFlag)
		
lRef=len(ref)
print("Length of reference genome: "+str(lRef))
#print(str(len(data))+" Samples in the input file.")
#vector to count how many bases of each type are cumulatively in the reference genome up to a certain position
cumulativeBases=[[0,0,0,0]]
for i in range(lRef):
	cumulativeBases.append(list(cumulativeBases[i]))
	cumulativeBases[i+1][allelesUpOrLow[ref[i]]]+=1
rootFreqs=[0.0,0.0,0.0,0.0]
rootFreqsLog=[0.0,0.0,0.0,0.0]
for i in range(4):
	rootFreqs[i]=cumulativeBases[-1][i]/float(lRef)
	rootFreqsLog[i]=log(rootFreqs[i])
refIndeces=[]
for i in range(lRef):
	refIndeces.append(allelesLow[ref[i]])
if model=="JC":
	rootFreqs=[0.25,0.25,0.25,0.25]
	rootFreqsLog=[log(0.25),log(0.25),log(0.25),log(0.25)]
oneMutBLen=1.0/lRef
minBLenSensitivity*=oneMutBLen

range4=range(4)
#stricter thresholds than the input one for probabilities.
thresholdProb2=thresholdProb*thresholdProb
thresholdProb4=thresholdProb2*thresholdProb2

#IMPORTANT definitions of genome list entries structure.
#The first element of a genome list entry represnts its type: 0="A", 1="C", 2="G", 3="T", 4="R", 5="N", 6="O"
#the second element represents the last position of the stretch of genome that they represent; 
# a value of 1 refers to the first position of the genome.
#For entries of type "O"/6 the last element is always a normalized vector of likelihoods (they always have to sum to 1.0).
#If present, another element after the position one represents the evolutionary distance since the considered type was observed (this is always missing with type "N"/5);
#If the element is not not present, this distance is assumed to be 0.0.
#If the distance (branch length) value is present, another distance value can also be present for entries of type <5 ; 
#this is to account for the fact that the observation might have occurred on the other side of the phylogeny with respect to the root; 
#when this additional branch length is also present, for example if the entry is (1, 234, 0.0001, 0.0002) then this means that a "C" was observed at genome position 234, and the observation
#is separated from the root by a distance of 0.0001, while a distance of 0.0002 separates the root from the current node (or position along a branch) considered.



#if probability mass is concentrated in one nucleotide, simplify the entry from "O" to another type.
def simplfy(vec,refA):
	maxP=0.0
	maxI=0
	numA=0
	for i in range4:
		if vec[i]>maxP:
			maxP=vec[i]
			maxI=i
		if vec[i]>thresholdProb:
			numA+=1
	if maxP<thresholdProb4:
		print("Inside simplify(), all values in vector are too small - something wrong?")
		print(vec)
		exit()
	if numA==1:
		if maxI==refA:
			return 4
		else:
			return maxI
	else:
		return 6



#Shorten genome list by merging together R entries that are mergeable
def shorten(vec):
	entryOld=vec[0]
	index=0
	while index<len(vec)-1:
		newVec=vec[index+1]
		if newVec[0]==4 and entryOld[0]==4 and len(newVec)==len(entryOld):
			if len(newVec)==2:
				vec.pop(index)
			elif abs(newVec[2]-entryOld[2])>thresholdProb:
				index+=1
				entryOld=vec[index]
			elif len(newVec)==3:
				vec.pop(index)
			elif abs(newVec[3]-entryOld[3])<thresholdProb:
				vec.pop(index)
			else:
				index+=1
				entryOld=vec[index]
		else:
			index+=1
			entryOld=vec[index]
	return


#Like shorten(), but specific to lower likelihoods that don't need to consider root frequencies. 
#UNNECESSARY, one could always use shorten() instead of it.
def shortenLower(vec):
	entryOld=vec[0]
	index=0
	while index<len(vec)-1:
		newVec=vec[index+1]
		if newVec[0]==4 and entryOld[0]==4 and len(newVec)==len(entryOld):
			if len(newVec)==2:
				vec.pop(index)
			elif abs(newVec[2]-entryOld[2])>thresholdProb:
				index+=1
				entryOld=vec[index]
			else:
				vec.pop(index)
		else:
			index+=1
			entryOld=vec[index]
	return



#define partial likelihood vector for a sample given its data
def probVectTerminalNode(diffs):
	if diffs is None:
		probVect=[(5,lRef)]
		return probVect
	pos=1
	probVect=[]
	for m in diffs:
			currPos=m[1]
			if currPos>pos:
				#region where the node with branch length bLen is identical to the ref.
				probVect.append((4,currPos-1))
				pos=currPos
			if m[0]=="n" or m[0]=="-":
				if len(m)>2:
					length=m[2]
				else:
					length=1
				#region with no info, store last position and length.
				probVect.append((5,currPos+length-1))
				pos=currPos+length
			elif m[0] in allelesLow:
				#position at which node allele is sure but is different from the reference.
				probVect.append((allelesLow[m[0]],currPos))
				pos=currPos+1
			else:
				# non-"n" ambiguity character; for now interpret this as ambiguity instead of as a polymorphism.
				if onlyNambiguities:
					# if user asks to, to make things easier, interpret any ambiguity as an "n".
					probVect.append((5,currPos))
				else:
					#otherwise, store as an "other" scenario, where each nucleotide has its own partial likelihood.
					probVect.append((6,currPos,ambiguities[m[0]]))
				pos=currPos+1
	if pos<=lRef:
		probVect.append((4,lRef))
	return probVect




#Update and normalize the mutation rate matrix, given new mutation counts
def updateSubMatrix(pseudoMutCounts,model,oldMutMatrix): #,mutMatrixBLen,mutMatrixBLenAdjusted
	mutMatrix=[[0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0]]
	if model=="UNREST":
		for i in range4:
			tot=0.0
			for j in range4:
				if j!=i:
					mutMatrix[i][j]=pseudoMutCounts[i][j]/rootFreqs[i]
					tot+=mutMatrix[i][j]
			mutMatrix[i][i]=-tot
	elif model=="GTR":
		for i in range4:
			tot=0.0
			for j in range4:
				if j!=i:
					mutMatrix[i][j]=(pseudoMutCounts[i][j]+pseudoMutCounts[j][i])/rootFreqs[i]
					tot+=mutMatrix[i][j]
			mutMatrix[i][i]=-tot
	else:
		print("Substitution model not recognised! Exiting")
		exit()
	totRate=-(rootFreqs[0]*mutMatrix[0][0]+ rootFreqs[1]*mutMatrix[1][1]+ rootFreqs[2]*mutMatrix[2][2]+ rootFreqs[3]*mutMatrix[3][3] )
	for i in range4:
		for j in range4:
			mutMatrix[i][j]=mutMatrix[i][j]/totRate
	matChange=0.0
	for i in range4:
		for j in range4:
			if j!=i:
				matChange+=abs(mutMatrix[i][j]-oldMutMatrix[i][j])
	if matChange > 0.001: #the matrix has changed significantly, update it
		for i in range4:
			for j in range4:
				oldMutMatrix[i][j]=mutMatrix[i][j]
		#the matrix has changed significantly, return True so that cumulative rates can be updated as well.
		return True
	else:
		return False







#merge two partial likelihood vectors, one from above, probVect1, and one from below, probVect2
#unlike appendProb(), this function is not used on a large part of the tree at each placement, but only in a small neighbourhood;
def mergeVectorsUpDown(probVect1,bLenUp,probVect2,bLenDown,mutMatrix,useRateVariation=False,mutMatrices=None):
	indexEntry1, indexEntry2, pos = 0, 0, 0
	probVect=[]
	entry1=probVect1[indexEntry1]
	entry2=probVect2[indexEntry2]
	while True:
		if entry1[0]==5:
			if entry2[0]==5:
				pos=min(entry1[1],entry2[1])
				probVect.append((5,pos))
				
			elif entry2[0]<5:
				pos=min(entry1[1],entry2[1])
				if len(entry2)==3:
					if bLenDown:
						probVect.append((entry2[0],pos,entry2[2]+bLenDown,0.0))
					else:
						probVect.append((entry2[0],pos,entry2[2],0.0))
				else:
					if bLenDown:
						probVect.append((entry2[0],pos,bLenDown,0.0))
					else:
						probVect.append((entry2[0],pos))
			else: # entry2 case "O", entry 1 is "N"
				if useRateVariation:
					mutMatrix=mutMatrices[pos]
				pos+=1
				if len(entry2)==4:
					totBLen=entry2[2]
					if bLenDown:
						totBLen+=bLenDown
				else:
					totBLen=bLenDown
				newVect=[]
				if totBLen:
					for i in range4:
						tot=0.0
						for j in range4:
							tot+=mutMatrix[i][j]*entry2[-1][j]
						tot*=totBLen
						tot+=entry2[-1][i]
						newVect.append(tot*rootFreqs[i])
					totSum=sum(newVect)
					for i in range4:
						newVect[i]/=totSum
					probVect.append((6,pos,newVect))
				else:
					for i in range4:
						newVect.append(entry2[-1][i]*rootFreqs[i])
					totSum=sum(newVect)
					for i in range4:
						newVect[i]/=totSum
					probVect.append((6,pos,newVect))
		elif entry2[0]==5: #entry2 is N
			if entry1[0]<5:
				pos=min(entry1[1],entry2[1])
				if len(entry1)==2:
					if bLenUp:
						probVect.append((entry1[0],pos,bLenUp))
					else:
						probVect.append((entry1[0],pos))
				elif len(entry1)==3:
					if bLenUp:
						probVect.append((entry1[0],pos,entry1[2]+bLenUp))
					else:
						probVect.append((entry1[0],pos,entry1[2]))
				else:
					if bLenUp:
						probVect.append((entry1[0],pos,entry1[2],entry1[3]+bLenUp))
					else:
						probVect.append((entry1[0],pos,entry1[2],entry1[3]))

			else: #entry1 is "O", entry 2 is "N"
				if useRateVariation:
					mutMatrix=mutMatrices[pos]
				pos+=1
				if len(entry1)==4:
					totBLen=entry1[2]
					if bLenUp:
						totBLen+=bLenUp
				elif bLenUp:
					totBLen=bLenUp
				else:
					totBLen=False
				newVect=[]
				if totBLen:
					for i in range4:
						tot=0.0
						for j in range4:
							tot+=entry1[-1][j]*mutMatrix[j][i]
						tot*=totBLen
						tot+=entry1[-1][i]
						newVect.append(tot)
					totSum=sum(newVect)
					for i in range4:
						newVect[i]/=totSum
					probVect.append((6,pos,newVect))
				else:
					probVect.append((6,pos,entry1[-1]))
		elif entry2[0]==entry1[0] and entry1[0]<5:
			pos=min(entry1[1],entry2[1])
			probVect.append((entry2[0],pos))
		else: #cases where the new genome list entry will likely be of type "O"
			if entry1[0]<5:
				if len(entry1)==2:
					totLen1=bLenUp
				else:
					totLen1=entry1[2]
					if bLenUp:
						totLen1+=bLenUp
					if len(entry1)==4:
						totLen1+=entry1[3]
			else:
				if len(entry1)==3:
					totLen1=bLenUp
				else:
					totLen1=entry1[2]
					if bLenUp:
						totLen1+=bLenUp
			
			if entry2[0]<5:
				if len(entry2)==2:
					totLen2=bLenDown
				else:
					totLen2=entry2[2]
					if bLenDown:
						totLen2+=bLenDown
			else:
				if len(entry2)==3:
					totLen2=bLenDown
				else:
					totLen2=entry2[2]
					if bLenDown:
						totLen2+=bLenDown
			if entry2[0]<5 and (not totLen2): #due to 0 distance, the entry will be of same type as entry2
				if (not totLen1) and entry1[0]<5:
					return None
				pos=min(entry1[1],entry2[1])
				probVect.append((entry2[0],pos))
			elif entry1[0]<5 and (not totLen1): #due to 0 distance, the entry will be of same type as entry1
				pos=min(entry1[1],entry2[1])
				probVect.append((entry1[0],pos))
			elif entry1[0]<5:
				if useRateVariation:
					mutMatrix=mutMatrices[pos]
				if entry1[0]==4:
					i1=refIndeces[pos]
				else:
					i1=entry1[0]
				newVec=[]
				if len(entry1)==4:
					rootVec=list(rootFreqs)
					for i in range4:
						if i==i1:
							rootVec[i]*=(1.0+mutMatrix[i1][i1]*(entry1[2]))
						else:
							rootVec[i]*=mutMatrix[i][i1]*(entry1[2])
					if bLenUp:
						lenToRoot=entry1[3]+bLenUp
					else:
						lenToRoot=entry1[3]
					for j in range4:
						tot=0.0
						for i in range4:
							tot+=mutMatrix[i][j]*rootVec[i]
						tot*=lenToRoot
						tot+=rootVec[j]
						newVec.append(tot)
				else:
					if totLen1:
						for i in range4:
							if i==i1:
								newVec.append(1.0+mutMatrix[i][i]*totLen1)
							else:
								newVec.append(mutMatrix[i1][i]*totLen1)
					else:
						for i in range4:
							if i==i1:
								newVec.append(1.0)
							else:
								newVec.append(0.0)
				if entry2[0]==6: #entry 2 is "O" and entry1 is a nucleotide
					for j in range4:
						tot=0.0
						if totLen2:
							for i in range4:
								tot+=mutMatrix[j][i]*entry2[-1][i]
							tot*=totLen2
						tot+=entry2[-1][j]
						newVec[j]*=tot
					sumV=sum(newVec)
					for i in range4:
						newVec[i]=newVec[i]/sumV
					state =simplfy(newVec,refIndeces[pos])
					pos+=1
					if state==6:
						probVect.append((6,pos,newVec))
					else:
						probVect.append((state,pos))
				else: #entry 1 and entry 2 are different nucleotides (possibly reference)
					if entry2[0]==4:
						i2=refIndeces[pos]
					else:
						i2=entry2[0]
					if totLen2:
						for i in range4:
							if i==i2:
								newVec[i]*=1.0+mutMatrix[i][i]*totLen2
							else:
								newVec[i]*=mutMatrix[i][i2]*totLen2
					else:
						for i in range4:
							if i!=i2:
								newVec[i]=0
					sumV=sum(newVec)
					if not sumV:
						print("situation")
						print(entry1)
						print(entry2)
						print(bLenUp)
						print(bLenDown)
						print(totLen1)
						print(totLen2)
						print(newVec)
					for i in range4:
						newVec[i]=newVec[i]/sumV
					pos+=1
					probVect.append((6,pos,newVec))
			else: #entry1[0]==6:
				if useRateVariation:
					mutMatrix=mutMatrices[pos]
				if totLen1:
					newVec=[]
					for i in range4:
						tot=0.0
						for j in range4:
							tot+=mutMatrix[j][i]*entry1[-1][j]
						tot*=totLen1
						tot+=entry1[-1][i]
						newVec.append(tot)
				else:
					newVec=list(entry1[-1])

				if entry2[0]==6:
					if totLen2:
						for i in range4:
							tot=0.0
							for j in range4:
								tot+=mutMatrix[i][j]*entry2[-1][j]
							tot*=totLen2
							tot+=entry2[-1][i]
							newVec[i]*=tot
					else:
						for i in range4:
							newVec[i]*=entry2[-1][i]
				else: #entry1 is "O" and entry2 is a nucleotide
					if entry2[0]==4:
						i2=refIndeces[pos]
					else:
						i2=entry2[0]
					if totLen2:
						for i in range4:
							if i==i2:
								newVec[i]*=(1.0+mutMatrix[i][i]*totLen2)
							else:
								newVec[i]*=mutMatrix[i][i2]*totLen2
					else:
						for i in range4:
							if i!=i2:
								newVec[i]=0.0
				sumV=sum(newVec)
				if not sumV:
					return None
				for i in range4:
					newVec[i]=newVec[i]/sumV
				state =simplfy(newVec,refIndeces[pos])
				pos+=1
				if state==6:
					probVect.append((6,pos,newVec))
				else:
					probVect.append((state,pos))

		if pos==lRef:
			break
		if pos==entry1[1]:
			indexEntry1+=1
			entry1=probVect1[indexEntry1]
		if pos==entry2[1]:
			indexEntry2+=1
			entry2=probVect2[indexEntry2]

	if verbose:
		print("Merged up-down ")
		print(probVect)
	#check if the final  probVect can be simplified by merging consecutive entries
	shorten(probVect)
	if verbose:
		print("Shortened up-down merging ")
		print(probVect)
	return probVect







#merge two lower child partial likelihood vectors to create a new one 
#(and also calculate the logLk of the merging if necessary, which is currently only useful for the root but could also be useful to calculate the overall total likelihood of the tree).
def mergeVectors(probVect1,bLen1,probVect2,bLen2,mutMatrix,returnLK=False,useRateVariation=False,mutMatrices=None):
	indexEntry1, indexEntry2, pos = 0, 0, 0
	probVect=[]
	cumulPartLk=0.0
	entry1=probVect1[indexEntry1]
	entry2=probVect2[indexEntry2]
	end=min(entry1[1],entry2[1])
	while True:
		if entry1[0]==5:
			if entry2[0]==5:
				pos=min(entry1[1],entry2[1])
				probVect.append((5,pos))
			elif entry2[0]<5:
				pos=min(entry1[1],entry2[1])
				if len(entry2)==2:
					if bLen2:
						probVect.append((entry2[0],pos,bLen2))
					else:
						probVect.append((entry2[0],pos))
				else:
					if bLen2:
						probVect.append((entry2[0],pos,entry2[2]+bLen2))
					else:
						probVect.append((entry2[0],pos,entry2[2]))		
			else: # case entry2 is "O" and entry1 is "N"
				pos+=1
				if len(entry2)==3:
					if bLen2:
						probVect.append((6,pos,bLen2,entry2[-1]))
					else:
						probVect.append((6,pos,entry2[-1]))
				else:
					if bLen2:
						probVect.append((6,pos,entry2[2]+bLen2,entry2[-1]))
					else:
						probVect.append((6,pos,entry2[2],entry2[-1]))
		elif entry2[0]==5: #entry2 is N
			if entry1[0]<5:
				pos=min(entry1[1],entry2[1])
				if len(entry1)==2:
					if bLen1:
						probVect.append((entry1[0],pos,bLen1))
					else:
						probVect.append((entry1[0],pos))
				else:
					if bLen1:
						probVect.append((entry1[0],pos,entry1[2]+bLen1))
					else:
						probVect.append((entry1[0],pos,entry1[2]))
			else: #entry1 is "O" and entry2 is "N"
				pos+=1
				if len(entry1)==3:
					if bLen1:
						probVect.append((6,pos,bLen1,entry1[-1]))
					else:
						probVect.append((6,pos,entry1[-1]))
				else:
					if bLen1:
						probVect.append((6,pos,entry1[2]+bLen1,entry1[-1]))
					else:
						probVect.append((6,pos,entry1[2],entry1[-1]))
					
		else: #entry1 and entry2 are not "N"
			if entry1[0]<5:
				if len(entry1)==2:
					totLen1=bLen1
				else:
					totLen1=entry1[2]
					if bLen1:
						totLen1+=bLen1
			else:
				if len(entry1)==3:
					totLen1=bLen1
				else:
					totLen1=entry1[2]
					if bLen1:
						totLen1+=bLen1
			
			if entry2[0]<5:
				if len(entry2)==2:
					totLen2=bLen2
				else:
					totLen2=entry2[2]
					if bLen2:
						totLen2+=bLen2
			else:
				if len(entry2)==3:
					totLen2=bLen2
				else:
					totLen2=entry2[2]
					if bLen2:
						totLen2+=bLen2

			if entry2[0]==entry1[0] and entry2[0]<5: #entry1 and entry2 are two identical nucleotides
				end=min(entry1[1],entry2[1])
				probVect.append((entry2[0],end))
				if returnLK:
					if entry2[0]==4:
						cumulPartLk+=(totLen1+totLen2)*(cumulativeRate[end]-cumulativeRate[pos])
					else:
						if useRateVariation:
							cumulPartLk+=mutMatrices[pos][entry1[0]][entry1[0]]*(totLen1+totLen2)
						else:
							cumulPartLk+=nonMutRates[entry1[0]]*(totLen1+totLen2)
				pos=end
			elif (not totLen1) and (not totLen2) and entry1[0]<5 and entry2[0]<5: #0 distance between different nucleotides: merge is not possible
				if returnLK:
					print("mergeVectors() returning None 1")
					return None, float("-inf")
				else:
					return None
			elif entry1[0]<5: #entry1 is a nucleotide
				if useRateVariation:
					mutMatrix=mutMatrices[pos]
				if entry1[0]==4:
					i1=refIndeces[pos]
				else:
					i1=entry1[0]
				newVec=[]
				if totLen1:
					for i in range4:
						if i==i1:
							newVec.append(1.0+mutMatrix[i][i]*totLen1)
						else:
							newVec.append(mutMatrix[i][i1]*totLen1)
				else:
					newVec=[0.0,0.0,0.0,0.0]
					newVec[i1]=1.0

				if entry2[0]==6: #entry1 is a nucleotide and entry2 is "O"
					if totLen2:
						for j in range4:
							tot=0.0
							for i in range4:
								tot+=mutMatrix[j][i]*entry2[-1][i]
							tot*=totLen2
							tot+=entry2[-1][j]
							newVec[j]*=tot
					else:
						for j in range4:
							newVec[j]*=entry2[-1][j]
					sumV=sum(newVec)
					if not sumV:
						if returnLK:
							print("mergeVectors() returning None 2")
							return None, float("-inf")
						else:
							return None
					for i in range4:
						newVec[i]=newVec[i]/sumV
					state =simplfy(newVec,refIndeces[pos])
					pos+=1
					if state==6:
						probVect.append((6,pos,newVec))
					else:
						probVect.append((state,pos))
					if returnLK:
						cumulPartLk+=log(sumV)
				else: #entry1 and entry2 are nucleotides
					if entry2[0]==4:
						i2=refIndeces[pos]
					else:
						i2=entry2[0]
					if totLen2:
						for i in range4:
							if i==i2:
								newVec[i]*=1.0+mutMatrix[i][i]*totLen2
							else:
								newVec[i]*=mutMatrix[i][i2]*totLen2
						sumV=sum(newVec)
						for i in range4:
							newVec[i]=newVec[i]/sumV
						state =simplfy(newVec,refIndeces[pos])
						pos+=1
						if state==6:
							probVect.append((6,pos,newVec))
						else:
							probVect.append((state,pos))
						#probVect.append((6,pos,newVec))
						if returnLK:
							cumulPartLk+=log(sumV)
					else:
						pos+=1
						probVect.append((entry2[0],pos))
						if returnLK:
							cumulPartLk+=log(newVec[i2])
				
			else: #entry1[0]==6:
				if useRateVariation:
					mutMatrix=mutMatrices[pos]
				if totLen1:
					newVec=[]
					for i in range4:
						tot=0.0
						for j in range4:
							tot+=mutMatrix[i][j]*entry1[-1][j]
						tot*=totLen1
						tot+=entry1[-1][i]
						newVec.append(tot)
				else:
					newVec=list(entry1[-1])

				if entry2[0]==6:
					if totLen2:
						for i in range4:
							tot=0.0
							for j in range4:
								tot+=mutMatrix[i][j]*entry2[-1][j]
							tot*=totLen2
							tot+=entry2[-1][i]
							newVec[i]*=tot
					else:
						for i in range4:
							newVec[i]*=entry2[-1][i]
					sumV=sum(newVec)
					if not sumV:
						if returnLK:
							print("mergeVectors() returning None 3")
							return None, float("-inf")
						else:
							return None
					for i in range4:
						newVec[i]=newVec[i]/sumV
					state =simplfy(newVec,refIndeces[pos])
					pos+=1
					if state==6:
						probVect.append((6,pos,newVec))
					else:
						probVect.append((state,pos))
					if returnLK:
						cumulPartLk+=log(sumV)
				else: #entry2 is a nucleotide and entry1 is "O"
					if entry2[0]==4:
						i2=refIndeces[pos]
					else:
						i2=entry2[0]
					if totLen2:
						for i in range4:
							if i==i2:
								newVec[i]*=(1.0+mutMatrix[i][i]*totLen2)
							else:
								newVec[i]*=mutMatrix[i][i2]*totLen2
						sumV=sum(newVec)
						for i in range4:
							newVec[i]=newVec[i]/sumV
						state =simplfy(newVec,refIndeces[pos])
						pos+=1
						if state==6:
							probVect.append((6,pos,newVec))
						else:
							probVect.append((state,pos))
						if returnLK:
							cumulPartLk+=log(sumV)
					else:
						if not newVec[i2]:
							if returnLK:
								print("mergeVectors() returning None 4")
								return None, float("-inf")
							else:
								return None
						pos+=1
						probVect.append((entry2[0],pos))
						if returnLK:
							cumulPartLk+=log(newVec[i2])

		if pos==lRef:
			break
		if pos==entry1[1]:
			indexEntry1+=1
			entry1=probVect1[indexEntry1]
		if pos==entry2[1]:
			indexEntry2+=1
			entry2=probVect2[indexEntry2]

	if verbose:
		print("Merged vector at root ")
		print(probVect)
	#check if the final  probVect can be simplified by merging consecutive entries
	shorten(probVect)
	if verbose:
		print("Shortened root merging ")
		print(probVect)
	if returnLK:
		return probVect, cumulPartLk
	else:
		return probVect











#calculate the probability that results from combining a lower likelihood genome list of the root with root frequencies.
#Could maybe be combined with function rootVector() ?
def findProbRoot(probVect):
	logLK=0.0
	logFactor=1.0
	pos=0
	for entry in probVect:
		if entry[0]==4:
			for i in range4:
				logLK+=rootFreqsLog[i]*(cumulativeBases[entry[1]][i]-cumulativeBases[pos][i])
		elif entry[0]<4:
			logLK+=rootFreqsLog[entry[0]]
		elif entry[0]==6:
			tot=0.0
			for i in range4:
				tot+=rootFreqs[i]*entry[-1][i]
			logFactor*=tot
			if logFactor<=minimumCarryOver:
				try:
					if logFactor<sys.float_info.min:
						return float("-inf")
				except:
					print("In findProbRoot() value error")
					print(logFactor)
					return float("-inf")
				logLK+=log(logFactor)
				logFactor=1.0
		pos=entry[1]
	logLK+=log(logFactor)
	return logLK









#for the root, take lower likelihood genome list probVect, and create an overall likelihood (or upper right or upper left) genome list by multiplying likelihoods by root frequencies.
def rootVector(probVect,bLen,mutMatrix,useRateVariation=False,mutMatrices=None):
	newProbVect=[]
	for entry in probVect:
		if entry[0]==5:
			newProbVect.append(entry)
		elif entry[0]==6:
			if len(entry)==4:
				totBLen=entry[2]
				if bLen:
					totBLen+=bLen
			else:
				totBLen=bLen
			newVect=[]
			if totBLen:
				if useRateVariation:
					mutMatrix=mutMatrices[entry[1]-1]
				for i in range4:
					tot=0.0
					for j in range4:
						tot+=mutMatrix[i][j]*entry[-1][j]
					tot*=totBLen
					tot+=entry[-1][i]
					newVect.append(tot*rootFreqs[i])
				totSum=sum(newVect)
				for i in range4:
					newVect[i]/=totSum
				newProbVect.append((6,entry[1],newVect))
			else:
				for i in range4:
					newVect.append(entry[-1][i]*rootFreqs[i])
				totSum=sum(newVect)
				for i in range4:
					newVect[i]/=totSum
				newProbVect.append((6,entry[1],newVect))
		else:
			if len(entry)==3:
				if bLen:
					newProbVect.append((entry[0],entry[1],entry[2]+bLen,0.0))
				else:
					newProbVect.append((entry[0],entry[1],entry[2],0.0))
			else:
				if bLen:
					newProbVect.append((entry[0],entry[1],bLen,0.0))
				else:
					newProbVect.append((entry[0],entry[1]))
	return newProbVect







#function to add new mutation events in new sample to the pre-exisitng pseudocounts so to improve the estimate of the substitution rates.
# probVect1 is the genome list at the node where the appending happens; probVect2 is the genome list for the new sample.
def updatePesudoCounts(probVect1,probVect2,pseudoMutCounts):
	if model!="JC":
		indexEntry1, indexEntry2, pos = 0, 0, 0
		entry1=probVect1[indexEntry1]
		entry2=probVect2[indexEntry2]
		while True:
			if entry1[0]!=entry2[0] and entry1[0]<5 and entry2[0]<5:
				if entry1[0]==4:
					pseudoMutCounts[refIndeces[pos]][entry2[0]]+=1		
				elif entry2[0]==4:
					pseudoMutCounts[entry1[0]][refIndeces[pos]]+=1
				else:
					pseudoMutCounts[entry1[0]][entry2[0]]+=1
				pos+=1
			else:
				pos=min(entry1[1],entry2[1])

			if pos==lRef:
				break
			if pos==entry1[1]:
				indexEntry1+=1
				entry1=probVect1[indexEntry1]
			if pos==entry2[1]:
				indexEntry2+=1
				entry2=probVect2[indexEntry2]










#function to optimize branch lengths.
#calculate features of the derivative of the likelihood cost function wrt the branch length, then finds branch length that minimizes likelihood cost.
def estimateBranchLengthWithDerivativeOld(probVectP,probVectC,mutMatrix,useRateVariation=False,mutMatrices=None):
	c1=0.0
	ais=[]
	indexEntry1, indexEntry2, pos = 0, 0, 0
	entry1=probVectP[indexEntry1]
	entry2=probVectC[indexEntry2]
	end=min(entry1[1],entry2[1])
	while True:
		if entry2[0]==5: # case entry1 is N
			pos=min(entry1[1],entry2[1])
			pass
		elif entry1[0]==5: # case entry2 is N
			#if parent node is type "N", in theory we might have to calculate the contribution of root nucleotides; 
			# however, if this node is "N" then every other node in the current tree is "N", so we can ignore this since this contribution cancels out in relative terms.
			pos=min(entry1[1],entry2[1])
			pass
		else:
			#contribLength will be here the total length from the root or from the upper node, down to the down node.
			if entry1[0]<5:
				if len(entry1)==2:
					contribLength=False
				elif len(entry1)==3:
					contribLength=entry1[2]
				else:
					contribLength=entry1[3]
			else:
				if len(entry1)==3:
					contribLength=False
				else:
					contribLength=entry1[2]
			if entry2[0]<5:
				if len(entry2)==3:
					if contribLength:
						contribLength+=entry2[2]
					else:
						contribLength=entry2[2]
			else:
				if len(entry2)==4:
					if contribLength:
						contribLength+=entry2[2]
					else:
						contribLength=entry2[2]

			#below, when necessary, we represent the likelihood as coeff0*l +coeff1, where l is the branch length to be optimized.
			if entry1[0]==4: # case entry1 is R	
				if entry2[0]==4:
					end=min(entry1[1],entry2[1])
					c1+=(cumulativeRate[end]-cumulativeRate[pos])
					pos=end

				#entry1 is reference and entry2 is of type "O"
				elif entry2[0]==6:
					if useRateVariation:
						mutMatrix=mutMatrices[pos]
						#print(mutMatrix)
						#print(pos)
					i1=refIndeces[pos]
					if len(entry1)==4:
						coeff0=rootFreqs[i1]*entry2[-1][i1] 
						coeff1=0.0
						for i in range4:
							coeff0+=rootFreqs[i]*mutMatrix[i][i1]*entry1[2]*entry2[-1][i]
							coeff1+=mutMatrix[i1][i]*entry2[-1][i]
						coeff1*=rootFreqs[i1]
						if contribLength:
							coeff0+=coeff1*contribLength
					else:
						coeff0=entry2[-1][i1]
						coeff1=0.0
						for j in range4:
							coeff1+=mutMatrix[i1][j]*entry2[-1][j]
						if contribLength:
							coeff0+=coeff1*contribLength
					if coeff1<0.0:
						c1+=coeff1/coeff0
					elif coeff1:
						coeff0=coeff0/coeff1
						ais.append(coeff0)
					pos+=1

				else: #entry1 is R and entry2 is a different but single nucleotide
					if len(entry1)==4:
						if useRateVariation:
							mutMatrix=mutMatrices[pos]
						i1=refIndeces[pos]
						i2=entry2[0]
						if contribLength:
							coeff0=rootFreqs[i1]*mutMatrix[i1][i2]*contribLength+rootFreqs[i2]*mutMatrix[i2][i1]*entry1[2]
							coeff1=rootFreqs[i1]*mutMatrix[i1][i2]
						else:
							coeff0=rootFreqs[i2]*mutMatrix[i2][i1]*entry1[2]
							coeff1=rootFreqs[i1]*mutMatrix[i1][i2]
						coeff0=coeff0/coeff1
					else:
						if contribLength:
							coeff0=contribLength
						else:
							coeff0=0.0
					ais.append(coeff0)
					pos+=1

			# entry1 is of type "O"
			elif entry1[0]==6:
				if useRateVariation:
					mutMatrix=mutMatrices[pos]
				if entry2[0]==6:
					coeff0=entry1[-1][0]*entry2[-1][0]+entry1[-1][1]*entry2[-1][1]+entry1[-1][2]*entry2[-1][2]+entry1[-1][3]*entry2[-1][3]
					coeff1=0.0
					for i in range4:
						for j in range4:
							coeff1+=entry1[-1][i]*entry2[-1][j]*mutMatrix[i][j]
					if contribLength:
						coeff0+=coeff1*contribLength
				else: #entry1 is "O" and entry2 is a nucleotide
					if entry2[0]==4:
						i2=refIndeces[pos]
					else:
						i2=entry2[0]
					coeff0=entry1[-1][i2]
					coeff1=0.0
					for i in range4:
						coeff1+=entry1[-1][i]*mutMatrix[i][i2]
					if contribLength:
						coeff0+=coeff1*contribLength
				if coeff1<0.0:
					c1+=coeff1/coeff0
				elif coeff1:
					coeff0=coeff0/coeff1
					ais.append(coeff0)
				pos+=1

			else: #entry1 is a non-ref nuc
				if useRateVariation:
					mutMatrix=mutMatrices[pos]
				if entry2[0]==entry1[0]:
					c1+=mutMatrix[entry1[0]][entry1[0]]
				else: #entry1 is a nucleotide and entry2 is not the same as entry1
					i1=entry1[0]
					if entry2[0]<5: #entry2 is a nucleotide
						if entry2[0]==4:
							i2=refIndeces[pos]
						else:
							i2=entry2[0]

						if len(entry1)==4:
							if contribLength:
								coeff0=rootFreqs[i1]*mutMatrix[i1][i2]*contribLength+rootFreqs[i2]*mutMatrix[i2][i1]*entry1[2]
							else:
								coeff0=rootFreqs[i2]*mutMatrix[i2][i1]*entry1[2]
							coeff1=rootFreqs[i1]*mutMatrix[i1][i2]
							coeff0=coeff0/coeff1
						else:
							if contribLength:
								coeff0=contribLength
							else:
								coeff0=0.0
						ais.append(coeff0)

					else: #entry1 is a nucleotide and entry2 is of type "O"
						if len(entry1)==4:
							coeff0=rootFreqs[i1]*entry2[-1][i1] 
							coeff1=0.0
							for i in range4:
								coeff0+=rootFreqs[i]*mutMatrix[i][i1]*entry1[2]*entry2[-1][i]
								coeff1+=mutMatrix[i1][i]*entry2[-1][i]
							coeff1*=rootFreqs[i1]
							if contribLength:
								coeff0+=coeff1*contribLength
						else:
							coeff0=entry2[-1][i1]
							coeff1=0.0
							for j in range4:
								coeff1+=mutMatrix[i1][j]*entry2[-1][j]
							if contribLength:
								coeff0+=coeff1*contribLength
						if coeff1<0.0:
							c1+=coeff1/coeff0
						elif coeff1:
							coeff0=coeff0/coeff1
							ais.append(coeff0)
				pos+=1

		if pos==lRef:
			break
		if pos==entry1[1]:
			indexEntry1+=1
			entry1=probVectP[indexEntry1]
		if pos==entry2[1]:
			indexEntry2+=1
			entry2=probVectC[indexEntry2]

	#now optimized branch length based on coefficients
	c1=-c1
	n=len(ais)
	if n==0:
		return False
	else:
		if min(ais)<0.0:
			return 0.1
		tDown=min(0.1,n/c1-min(ais))
		if tDown<=0.0:
			return False
		vDown=0.0
		for ai in ais:
			vDown+=1.0/(ai+tDown)
		tUp=min(0.1,n/c1-max(ais))
		if tUp>=0.1:
			return 0.1
		#TODO Myrthe, changed here
		if tUp<=minBLenSensitivity:
			if min(ais):
				tUp=0.0
			else:
				tUp=minBLenSensitivity
		vUp=0.0
		for ai in ais:
			vUp+=1.0/(ai+tUp)
	if vDown>c1+minBLenSensitivity or vUp<c1-minBLenSensitivity:
		if vUp<c1-minBLenSensitivity and (not tUp):
			return False
		if (vDown>c1+minBLenSensitivity) and tDown>=0.1:
			return 0.1
		print(c1,ais)
		print("Initial values")
		print(tDown,vDown,tUp,vUp)
		print("Initial border parameters don't fit expectations:")
		#print(probVectP)
		#print(probVectC)
		#print(mutMatrix)
	
	while tDown-tUp>minBLenSensitivity:
		tMiddle=(tUp+tDown)/2
		vMiddle=0.0
		for ai in ais:
			vMiddle+=1.0/(ai+tMiddle)
		if vMiddle>c1:
			tUp=tMiddle
		else:
			tDown=tMiddle
	return tUp




#function to optimize branch lengths.
#calculate features of the derivative of the likelihood cost function wrt the branch length, then finds branch length that minimizes likelihood cost.
def estimateBranchLengthWithDerivative(probVectP,probVectC,mutMatrix,useRateVariation=False,mutMatrices=None):
	c1=0.0
	ais=[]
	indexEntry1, indexEntry2, pos = 0, 0, 0
	entry1=probVectP[indexEntry1]
	entry2=probVectC[indexEntry2]
	end=min(entry1[1],entry2[1])
	nZeros=0
	while True:
		if entry2[0]==5: # case entry1 is N
			pos=min(entry1[1],entry2[1])
			pass
		elif entry1[0]==5: # case entry2 is N
			#if parent node is type "N", in theory we might have to calculate the contribution of root nucleotides; 
			# however, if this node is "N" then every other node in the current tree is "N", so we can ignore this since this contribution cancels out in relative terms.
			pos=min(entry1[1],entry2[1])
			pass
		else:
			#contribLength will be here the total length from the root or from the upper node, down to the down node.
			if entry1[0]<5:
				if len(entry1)==2:
					contribLength=False
				elif len(entry1)==3:
					contribLength=entry1[2]
				else:
					contribLength=entry1[3]
			else:
				if len(entry1)==3:
					contribLength=False
				else:
					contribLength=entry1[2]
			if entry2[0]<5:
				if len(entry2)==3:
					if contribLength:
						contribLength+=entry2[2]
					else:
						contribLength=entry2[2]
			else:
				if len(entry2)==4:
					if contribLength:
						contribLength+=entry2[2]
					else:
						contribLength=entry2[2]

			#below, when necessary, we represent the likelihood as coeff0*l +coeff1, where l is the branch length to be optimized.
			if entry1[0]==4: # case entry1 is R	
				if entry2[0]==4:
					end=min(entry1[1],entry2[1])
					c1+=(cumulativeRate[end]-cumulativeRate[pos])
					pos=end

				#entry1 is reference and entry2 is of type "O"
				elif entry2[0]==6:
					if useRateVariation:
						mutMatrix=mutMatrices[pos]
						#print(mutMatrix)
						#print(pos)
					i1=refIndeces[pos]
					if len(entry1)==4:
						coeff0=rootFreqs[i1]*entry2[-1][i1] 
						coeff1=0.0
						for i in range4:
							coeff0+=rootFreqs[i]*mutMatrix[i][i1]*entry1[2]*entry2[-1][i]
							coeff1+=mutMatrix[i1][i]*entry2[-1][i]
						coeff1*=rootFreqs[i1]
						if contribLength:
							coeff0+=coeff1*contribLength
					else:
						coeff0=entry2[-1][i1]
						coeff1=0.0
						for j in range4:
							coeff1+=mutMatrix[i1][j]*entry2[-1][j]
						if contribLength:
							coeff0+=coeff1*contribLength
					if coeff1<0.0:
						c1+=coeff1/coeff0
					elif coeff1:
						coeff0=coeff0/coeff1
						ais.append(coeff0)
					pos+=1

				else: #entry1 is R and entry2 is a different but single nucleotide
					if len(entry1)==4:
						if useRateVariation:
							mutMatrix=mutMatrices[pos]
						i1=refIndeces[pos]
						i2=entry2[0]
						if contribLength:
							coeff0=rootFreqs[i1]*mutMatrix[i1][i2]*contribLength+rootFreqs[i2]*mutMatrix[i2][i1]*entry1[2]
							coeff1=rootFreqs[i1]*mutMatrix[i1][i2]
						else:
							coeff0=rootFreqs[i2]*mutMatrix[i2][i1]*entry1[2]
							coeff1=rootFreqs[i1]*mutMatrix[i1][i2]
						coeff0=coeff0/coeff1
					else:
						if contribLength:
							coeff0=contribLength
						else:
							coeff0=0.0
					if coeff0:
						ais.append(coeff0)
					else:
						nZeros+=1
					pos+=1

			# entry1 is of type "O"
			elif entry1[0]==6:
				if useRateVariation:
					mutMatrix=mutMatrices[pos]
				if entry2[0]==6:
					coeff0=entry1[-1][0]*entry2[-1][0]+entry1[-1][1]*entry2[-1][1]+entry1[-1][2]*entry2[-1][2]+entry1[-1][3]*entry2[-1][3]
					coeff1=0.0
					for i in range4:
						for j in range4:
							coeff1+=entry1[-1][i]*entry2[-1][j]*mutMatrix[i][j]
					if contribLength:
						coeff0+=coeff1*contribLength
				else: #entry1 is "O" and entry2 is a nucleotide
					if entry2[0]==4:
						i2=refIndeces[pos]
					else:
						i2=entry2[0]
					coeff0=entry1[-1][i2]
					coeff1=0.0
					for i in range4:
						coeff1+=entry1[-1][i]*mutMatrix[i][i2]
					if contribLength:
						coeff0+=coeff1*contribLength
				if coeff1<0.0:
					c1+=coeff1/coeff0
				elif coeff1:
					coeff0=coeff0/coeff1
					ais.append(coeff0)
				pos+=1

			else: #entry1 is a non-ref nuc
				if useRateVariation:
					mutMatrix=mutMatrices[pos]
				if entry2[0]==entry1[0]:
					c1+=mutMatrix[entry1[0]][entry1[0]]
				else: #entry1 is a nucleotide and entry2 is not the same as entry1
					i1=entry1[0]
					if entry2[0]<5: #entry2 is a nucleotide
						if entry2[0]==4:
							i2=refIndeces[pos]
						else:
							i2=entry2[0]

						if len(entry1)==4:
							if contribLength:
								coeff0=rootFreqs[i1]*mutMatrix[i1][i2]*contribLength+rootFreqs[i2]*mutMatrix[i2][i1]*entry1[2]
							else:
								coeff0=rootFreqs[i2]*mutMatrix[i2][i1]*entry1[2]
							coeff1=rootFreqs[i1]*mutMatrix[i1][i2]
							coeff0=coeff0/coeff1
						else:
							if contribLength:
								coeff0=contribLength
							else:
								coeff0=0.0
						if coeff0:
							ais.append(coeff0)
						else:
							nZeros+=1
						#ais.append(coeff0)

					else: #entry1 is a nucleotide and entry2 is of type "O"
						if len(entry1)==4:
							coeff0=rootFreqs[i1]*entry2[-1][i1] 
							coeff1=0.0
							for i in range4:
								coeff0+=rootFreqs[i]*mutMatrix[i][i1]*entry1[2]*entry2[-1][i]
								coeff1+=mutMatrix[i1][i]*entry2[-1][i]
							coeff1*=rootFreqs[i1]
							if contribLength:
								coeff0+=coeff1*contribLength
						else:
							coeff0=entry2[-1][i1]
							coeff1=0.0
							for j in range4:
								coeff1+=mutMatrix[i1][j]*entry2[-1][j]
							if contribLength:
								coeff0+=coeff1*contribLength
						if coeff1<0.0:
							c1+=coeff1/coeff0
						elif coeff1:
							coeff0=coeff0/coeff1
							ais.append(coeff0)
				pos+=1

		if pos==lRef:
			break
		if pos==entry1[1]:
			indexEntry1+=1
			entry1=probVectP[indexEntry1]
		if pos==entry2[1]:
			indexEntry2+=1
			entry2=probVectC[indexEntry2]

	#now optimized branch length based on coefficients
	c1=-c1
	n=len(ais)+nZeros
	if n==0:
		return False
	else:
		if len(ais):
			minAis=min(ais)
		else:
			minAis=0.0
		if nZeros:
			minAis=min(0.0,minAis)
		if minAis<0.0:
			return 0.1
		tDown=min(0.1,n/c1-minAis)
		if tDown<=0.0:
			return False
		if nZeros:
			vDown=nZeros/tDown
		else:
			vDown=0.0
		for ai in ais:
			vDown+=1.0/(ai+tDown)
		if len(ais):
			maxAis=max(ais)
		else:
			maxAis=0.0
		tUp=min(0.1,n/c1-maxAis)
		if tUp>=0.1:
			return 0.1
		#TODO Myrthe, changed here
		if tUp<=minBLenSensitivity:
			if minAis:
				tUp=0.0
			else:
				tUp=minBLenSensitivity
		if nZeros:
			vUp=nZeros/tUp
		else:
			vUp=0.0
		for ai in ais:
			vUp+=1.0/(ai+tUp)
	if vDown>c1+minBLenSensitivity or vUp<c1-minBLenSensitivity:
		if vUp<c1-minBLenSensitivity and (not tUp):
			return False
		if (vDown>c1+minBLenSensitivity) and tDown>=0.1:
			return 0.1
		print(c1,ais,nZeros)
		print("Initial values")
		print(tDown,vDown,tUp,vUp)
		print("Initial border parameters don't fit expectations:")
		#print(probVectP)
		#print(probVectC)
		#print(mutMatrix)
	
	while tDown-tUp>minBLenSensitivity:
		tMiddle=(tUp+tDown)/2
		if nZeros:
			vMiddle=nZeros/tMiddle
		else:
			vMiddle=0.0
		for ai in ais:
			vMiddle+=1.0/(ai+tMiddle)
		if vMiddle>c1:
			tUp=tMiddle
		else:
			tDown=tMiddle
	return tUp









#if updating genome lists in updatePartials() creates an inconsistency, this function can increase the length of a 0-length branch to resolve the inconsistency.
#In doing so, it updates, the input list of nodes to visit and update.
def updateBLen(node,mutMatrix,useRateVariation=False,mutMatrices=None,addToList=False,nodeList=None):
	cNode=node
	node=node.up
	if cNode==node.children[0]:
		vectUp=node.probVectUpRight
		cNum=0
	else:
		vectUp=node.probVectUpLeft
		cNum=1
	bestLength=estimateBranchLengthWithDerivative(vectUp,cNode.probVect,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices)
	#if not bestLength:
	#	print("Estimated 0 branch length in updateBLen()")
	cNode.dist=bestLength
	node.dirty=True
	cNode.dirty=True
	if addToList:
		nodeList.append((cNode,2))
		nodeList.append((node,cNum))




# update the partials iteratively starting from the nodes in nodeList
#each entry in nodeList contains the node it refers to, and the direction where the update comes from (0 is left child, 1 is right child, 2 is parent)
def updatePartials(nodeList,mutMatrix,useRateVariation=False,mutMatrices=None):
	while nodeList:
		updatedBLen=False # if there has been an inconsistency, function updateBLen() has been called, and so there is no point continuing with some updates.
		madeChange=False # some change has been made, so continue traversal
		node, direction = nodeList.pop()
		node.dirty=True
		if node.up != None:
			try:
				if node==node.up.children[0]:
					childNumUp=0
					vectUpUp=node.up.probVectUpRight
				else:
					childNumUp=1
					vectUpUp=node.up.probVectUpLeft
			except AttributeError:
				vectUpUp=None
		#change in likelihoods is coming from parent node
		if direction==2:
			if node.dist : #if necessary, update the total probabilities at the mid node.
				newTot=mergeVectorsUpDown(vectUpUp,node.dist/2,node.probVect,node.dist/2,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices)
				if newTot==None:
					if node.dist>1e-100:
						print("inside updatePartials(), from parent: should not have happened since node.dist>0")
					updateBLen(node,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices)
					nodeList.append((node.up,childNumUp))
					newTot=mergeVectorsUpDown(vectUpUp,node.dist/2,node.probVect,node.dist/2,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices)
					madeChange=True
				#else:
				node.probVectTotUp=newTot
			if len(node.children)>0:# and (not updatedBLen): #at valid internal node, update upLeft and upRight, and if necessary add children to nodeList.
				child0Vect=node.children[0].probVect
				child1Vect=node.children[1].probVect
				dist0=node.children[0].dist
				dist1=node.children[1].dist
				newUpRight=mergeVectorsUpDown(vectUpUp,node.dist,child1Vect,dist1,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices)
				if newUpRight==None:
					if (not node.dist) and (not dist1):
						updateBLen(node,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices)
						if not node.dist:
							updateBLen(node.children[1],mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices,addToList=True,nodeList=nodeList)
							updatedBLen=True
						else:
							node.probVectTotUp=mergeVectorsUpDown(vectUpUp,node.dist/2,node.probVect,node.dist/2,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices)
							newUpRight=mergeVectorsUpDown(vectUpUp,node.dist,child1Vect,dist1,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices)
							nodeList.append((node.up,childNumUp))
							madeChange=True
					else:
						print("Strange: None vector from non-zero distances in updatePartials() from parent direction.")
						exit()
				if not updatedBLen:
					newUpLeft=mergeVectorsUpDown(vectUpUp,node.dist,child0Vect,dist0,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices)
					if newUpLeft==None:
						if (not node.dist) and (not dist0) :
							updateBLen(node,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices)
							if not node.dist:
								updateBLen(node.children[0],mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices,addToList=True,nodeList=nodeList)
								updatedBLen=True
							else:
								node.probVectTotUp=mergeVectorsUpDown(vectUpUp,node.dist/2,node.probVect,node.dist/2,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices)
								newUpRight=mergeVectorsUpDown(vectUpUp,node.dist,child1Vect,dist1,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices)
								newUpLeft=mergeVectorsUpDown(vectUpUp,node.dist,child0Vect,dist0,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices)
								nodeList.append((node.up,childNumUp))
								madeChange=True
						else:
							print("Strange: None vector from non-zero distances in updatePartials() from parent direction, child0.")
							exit()
				if not updatedBLen:
					if madeChange or areVectorsDifferent(node.probVectUpRight,newUpRight):
						node.probVectUpRight=newUpRight
						nodeList.append((node.children[0],2))
					if madeChange or areVectorsDifferent(node.probVectUpLeft,newUpLeft):
						node.probVectUpLeft=newUpLeft
						nodeList.append((node.children[1],2))


		else: #change in likelihoods is coming from child number "direction".
			childNum=direction
			otherChildNum=1-childNum
			childDist=node.children[childNum].dist
			otherChildDist=node.children[otherChildNum].dist
			otherChildVect=node.children[otherChildNum].probVect
			probVectDown=node.children[childNum].probVect
			try:
				if childNum:
					otherVectUp=node.probVectUpRight
				else:
					otherVectUp=node.probVectUpLeft
			except AttributeError:
				otherVectUp=None

			#update lower likelihoods
			newVect=mergeVectors(otherChildVect,otherChildDist,probVectDown,childDist,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices)
			if newVect==None:
				if (not childDist) and (not otherChildDist):
					updateBLen(node.children[childNum],mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices)
					if not node.children[childNum].dist:
						updateBLen(node.children[otherChildNum],mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices,addToList=True,nodeList=nodeList)
						updatedBLen=True
					else:
						childDist=node.children[childNum].dist
						node.probVect=mergeVectors(otherChildVect,otherChildDist,probVectDown,childDist,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices)
						nodeList.append((node.children[childNum],2))
						madeChange=True
				else:
					print("Strange: None vector from non-zero distances in updatePartials() from child direction.")
					print(childDist)
					print(probVectDown)
					print(otherChildDist)
					print(otherChildVect)
					exit()
			else:
				try:
					oldProbVect=node.probVect
				except AttributeError:
					oldProbVect=None
				node.probVect=newVect

			#update total mid-branches likelihood
			if (not updatedBLen) and node.dist and (node.up != None) and (vectUpUp!=None):
				newTot=mergeVectorsUpDown(vectUpUp,node.dist/2,node.probVect,node.dist/2,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices)
				if newTot==None:
					updateBLen(node,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices)
					node.probVect=mergeVectors(otherChildVect,otherChildDist,probVectDown,childDist,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices)
					nodeList.append((node.children[childNum],2))
					node.probVectTotUp=mergeVectorsUpDown(vectUpUp,node.dist/2,node.probVect,node.dist/2,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices)
					madeChange=True
					print("inside updatePartials(), from child: should not have happened since node.dist>0")
				else:
					node.probVectTotUp=newTot
					
			if (not updatedBLen) and (otherVectUp!=None):
				#update likelihoods at sibling node
				if node.up != None:
					if verbose:
						print("Trying to move update to sibling while updating partials")
					newUpVect=mergeVectorsUpDown(vectUpUp,node.dist,probVectDown,childDist,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices)
				else:
					if verbose:
						print("reached root while moving up updating likelihoods and trying to move to sibling")
					newUpVect=rootVector(probVectDown,childDist,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices)
				if newUpVect==None:
					if (not node.dist) and (not childDist):
						updateBLen(node,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices)
						if not node.dist:
							updateBLen(node.children[childNum],mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices,addToList=True,nodeList=nodeList)
							updatedBLen=True
						else:
							node.probVectTotUp=mergeVectorsUpDown(vectUpUp,node.dist/2,node.probVect,node.dist/2,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices)
							nodeList.append((node.children[childNum],2))
							madeChange=True
							newUpVect=mergeVectorsUpDown(vectUpUp,node.dist,probVectDown,childDist,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices)
					else:
						print("Strange: None vector from non-zero distances in updatePartials() from child direction, newUpVect.")
						exit()
			if (not updatedBLen) and (otherVectUp!=None):
				if madeChange or areVectorsDifferent(otherVectUp,newUpVect):
					if verbose:
						print("Moving update to sibling; newUpVect and otherVectUp")
						print(newUpVect)
						print(otherVectUp)
					if childNum:
						node.probVectUpRight=newUpVect
					else:
						node.probVectUpLeft=newUpVect
					nodeList.append((node.children[otherChildNum],2))

			if (not updatedBLen):
				#update likelihoods at parent node
				if madeChange or areVectorsDifferent(node.probVect,oldProbVect):
					if node.up != None:
						nodeList.append((node.up,childNumUp))
					if verbose:
						print("Old upRight and UpLeft:")
						print(node.probVectUpRight)
						print(node.probVectUpLeft)
		if verbose:
			print("New down, upRight, UpLeft and tot:")
			print(node.probVect)
			print(node.probVectUpRight)
			print(node.probVectUpLeft)












#Check if two genome lists represent the same partial likelihoods or not.
#This is used to traverse the tree and only update genome lists if they changed after a change in the tree.
def areVectorsDifferent(probVect1,probVect2):
	if probVect2==None:
		return True
	indexEntry1, indexEntry2, pos = 0, 0, 0
	entry1=probVect1[indexEntry1]
	entry2=probVect2[indexEntry2]
	while True:
		if entry1[0]!=entry2[0]:	
			return True
		if len(entry1)!=len(entry2):
			return True
		if entry1[0]<5:
			if len(entry1)>2:
				if abs(entry1[2] - entry2[2])>thresholdProb:
					return True
				if len(entry1)==4:
					if abs(entry1[3] - entry2[3])>thresholdProb:
						return True
		if entry1[0]==6:
			if len(entry1)==4:
				if abs(entry1[2] - entry2[2])>thresholdProb:
					return True
			for i in range4:
				diffVal=abs(entry1[-1][i] - entry2[-1][i])
				if diffVal:
					if (not entry1[-1][i]) or (not entry2[-1][i]):
						return True
					if diffVal>thresholdDiffForUpdate or (diffVal>thresholdProb and ( (diffVal/entry1[-1][i]>thresholdFoldChangeUpdate)  or  (diffVal/entry2[-1][i]>thresholdFoldChangeUpdate))):
						return True
		#update pos, end, etc
		pos=min(entry1[1],entry2[1])
		if pos==lRef:
			break
		if pos==entry1[1]:
			indexEntry1+=1
			entry1=probVect1[indexEntry1]
		if pos==entry2[1]:
			indexEntry2+=1
			entry2=probVect2[indexEntry2]
	return False



#Check if two genome lists represent the same partial likelihoods or not.
#This is a less strict version used for debugging purposes.
def areVectorsDifferentDebugging(probVect1,probVect2,threshold=0.00001):
	if probVect2==None:
		print('none case')
		return True
	indexEntry1, indexEntry2, pos = 0, 0, 0
	entry1=probVect1[indexEntry1]
	entry2=probVect2[indexEntry2]
	while True:
		if entry1[0]!=6 and entry2[0]!=6:
			if entry1[0]!=entry2[0]:
				print('diff type case')
				return True
			if len(entry1)!=len(entry2):
				print('diff len case')
				return True
			if entry1[0]<5:
				if len(entry1)>2:
					if abs(entry1[2] - entry2[2])>thresholdProb:
						print('diff bLen case')
						return True
					if len(entry1)==4:
						if abs(entry1[3] - entry2[3])>thresholdProb:
							print('diff bLen2 case')
							return True
		elif entry1[0]==6 and entry2[0]==6 :
			if len(entry1)==4 and len(entry2)==4 :
				if abs(entry1[2] - entry2[2])>thresholdProb:
					print("66 bLen case")
					return True
			elif len(entry1)!=len(entry2):
				print("66 entry len")
				return True
			for i in range4:
				diffVal=abs(entry1[-1][i] - entry2[-1][i])
				if diffVal:
					if (not entry1[-1][i]) or (not entry2[-1][i]):
						print("66 0 like case")
						return True
					if diffVal>0.01 or (diffVal>threshold and ( (diffVal/entry1[-1][i]>thresholdFoldChangeUpdate)  or  (diffVal/entry2[-1][i]>thresholdFoldChangeUpdate))):
						print("66 like diff case")
						print(entry1)
						print(entry2)
						return True
		else:
			if not (entry1[0]==5 and entry2[0]==5):
				if (entry1[0]==5 or entry2[0]==5):
					print("N case")
					return True
				if entry1[0]<5:
					if entry1[0]==4:
						i1=refIndeces[pos]
					else:
						i1=entry1[0]
					if entry2[-1][i1]+threshold<1.0:
						print("i1 case")
						print(entry1)
						print(entry2)
						return True
				elif entry2[0]<5:
					if entry2[0]==4:
						i2=refIndeces[pos]
					else:
						i2=entry2[0]
					if entry1[-1][i2]+threshold<1.0:
						print("i2 case")
						print(entry1)
						print(entry2)
						return True
		#update pos, end, etc
		pos=min(entry1[1],entry2[1])
		if pos==lRef:
			break
		if pos==entry1[1]:
			indexEntry1+=1
			entry1=probVect1[indexEntry1]
		if pos==entry2[1]:
			indexEntry2+=1
			entry2=probVect2[indexEntry2]
	return False












numNodes=[0,0,0,0,0]

#Given a tree, and a final substitution rate matrix, re-calculate all genome lists within the tree according to this matrix.
# this is useful once the matrix estimation has finished, to make sure all genome lists replect this matrix. 
def reCalculateAllGenomeLists(root,mutMatrix, checkExistingAreCorrect=False,countNodes=False,countPseudoCounts=False,pseudoMutCounts=None,data=None,useRateVariation=False,mutMatrices=None,firstSetUp=False):
	#first pass to update all lower likelihoods.
	node=root
	#direction 0 means from parent, direction 1 means from a child
	lastNode=None
	direction=0
	dataNamesConverted=False
	while node!=None:
		if direction==0:
			if node.children:
				node=node.children[0]
			else:
				if data!=None:
					if node.name in data:
						node.probVect=probVectTerminalNode(data[node.name])
						del data[node.name]
					elif not dataNamesConverted:
						print("Sample name "+node.name+" not found in the input sequence data - all samples in the input tree need to have a corresponding sequence entry. Converting sequence names by replacing ? and & characters with _ and trying again.")
						names=list(data.keys())
						for name in names:
							newName=name.replace("?","_").replace("&","_")
							if newName!=name:
								data[newName]=data[name]
						dataNamesConverted=True
						names=None
						if node.name in data:
							node.probVect=probVectTerminalNode(data[node.name])
							del data[node.name]
						else:
							print("Error: sample name "+node.name+" not found in the input sequence data - all samples in the input tree need to have a corresponding sequence entry.")
							exit()
					else:
						print("Error: sample name "+node.name+" not found in the input sequence data - all samples in the input tree need to have a corresponding sequence entry.")
						exit()
				if countNodes:
					numNodes[0]+=1
					for entry in node.probVect:
						if entry[0]<4:
							numNodes[1]+=1
						elif entry[0]==4:
							numNodes[2]+=1
						elif entry[0]==5:
							numNodes[3]+=1
						else:
							numNodes[4]+=1
				lastNode=node
				node=node.up
				direction=1
		else :
			if lastNode==node.children[0]:
				node=node.children[1]
				direction=0
			else:
				newLower=mergeVectors(node.children[0].probVect,node.children[0].dist,node.children[1].probVect,node.children[1].dist,mutMatrix,returnLK=False,useRateVariation=useRateVariation,mutMatrices=mutMatrices)
				if checkExistingAreCorrect:
					if areVectorsDifferentDebugging(newLower,node.probVect):
						print("Inside reCalculateAllGenomeLists(), new lower at node is different from the old one, and it shouldn't be.")
						print(newLower)
						print(node.probVect)
						exit()
				if newLower==None:
					if (not node.children[0].dist) and (not node.children[1].dist):
						if firstSetUp:
							node.children[0].dist=oneMutBLen/2
							node.children[1].dist=oneMutBLen/2
						else:
							updateBLen(node.children[0],mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices)
							if (not node.children[0].dist):
								updateBLen(node.children[1],mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices)
						node.probVect=mergeVectors(node.children[0].probVect,node.children[0].dist,node.children[1].probVect,node.children[1].dist,mutMatrix,returnLK=False,useRateVariation=useRateVariation,mutMatrices=mutMatrices)
						if node.probVect==None:
							node.children[0].dist=oneMutBLen/2
							node.children[1].dist=oneMutBLen/2
							node.probVect=mergeVectors(node.children[0].probVect,node.children[0].dist,node.children[1].probVect,node.children[1].dist,mutMatrix,returnLK=False,useRateVariation=useRateVariation,mutMatrices=mutMatrices)
							if node.probVect==None:
								print(node.children[0].probVect)
								print(node.children[1].probVect)
								print("None vector when merging these two vectors in reCalculateAllGenomeLists(), despite updating branch lengths.")
								exit()
					else:
						print("Strange, distances>0 but inconsistent lower genome list creation in reCalculateAllGenomeLists()")
						exit()
				else:
					node.probVect=newLower
				if countNodes:
					numNodes[0]+=1
					for entry in node.probVect:
						if entry[0]<4:
							numNodes[1]+=1
						elif entry[0]==4:
							numNodes[2]+=1
						elif entry[0]==5:
							numNodes[3]+=1
						else:
							numNodes[4]+=1

				lastNode=node
				node=node.up
				direction=1

	#now update the other genome lists for the root
	node=root
	if node.children:
		newVect=rootVector(node.children[1].probVect,node.children[1].dist,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices)
		if checkExistingAreCorrect:
			if areVectorsDifferentDebugging(newVect,node.probVectUpRight):
				print("new probVectUpRight at root is different from the old one, and it shouldn't be.")
				print(newVect)
				print(node.probVectUpRight)
				exit()
		node.probVectUpRight=newVect
		newVect=rootVector(node.children[0].probVect,node.children[0].dist,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices)
		if checkExistingAreCorrect:
			if areVectorsDifferentDebugging(newVect,node.probVectUpLeft):
				print("new probVectUpLeft at root is different from the old one, and it shouldn't be; distance of child node "+str(node.children[0].dist))
				print(node.dist)
				print(node.children)
				print(len(node.children))
				print(node.children[0].dist)
				print(node.children[1].dist)
				print(node.children[0].probVect)
				print(createBinaryNewick(node))
				print(newVect)
				print(node.probVectUpLeft)
				exit()
		node.probVectUpLeft=newVect
		
		#now traverse the tree downward and update the non-lower genome lists for all other nodes of the tree.
		totNodeList=[]
		lastNode=None
		node=node.children[0]
		direction=0
		while node!=None:
			if direction==0:
				if node==node.up.children[0]:
					vectUp=node.up.probVectUpRight
					nodeChildNum=0
				else:
					vectUp=node.up.probVectUpLeft
					nodeChildNum=1
				if node.dist:
					if countPseudoCounts:
						updatePesudoCounts(vectUp,node.probVect,pseudoMutCounts)
					newVect=mergeVectorsUpDown(vectUp,node.dist/2,node.probVect,node.dist/2,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices)
					if checkExistingAreCorrect:
						if areVectorsDifferentDebugging(newVect,node.probVectTotUp):
							print("new probVectTotUp at node is different from the old one, and it shouldn't be.")
							print(newVect)
							print(node.probVectTotUp)
							exit()
					node.probVectTotUp=newVect
				if node.children:
					newUpRight=mergeVectorsUpDown(vectUp,node.dist,node.children[1].probVect,node.children[1].dist,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices)
					if newUpRight==None:
						if (not node.children[1].dist) and (not node.dist):
							#nodeList=[]
							updateBLen(node,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices)
							if not node.dist:
								if firstSetUp:
									node.probVectUpLeft=mergeVectorsUpDown(vectUp,node.dist,node.children[0].probVect,node.children[0].dist,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices)
								updateBLen(node.children[1],mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices)
								totNodeList.append((node,1))
							else:
								node.probVectTotUp=mergeVectorsUpDown(vectUp,node.dist/2,node.probVect,node.dist/2,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices)
								totNodeList.append((node.up,nodeChildNum))
							node.probVectUpRight=mergeVectorsUpDown(vectUp,node.dist,node.children[1].probVect,node.children[1].dist,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices)
						else:
							print("Strange, distances>0 but inconsistent upRight genome list creation in reCalculateAllGenomeLists()")
							exit()
					else:
						if checkExistingAreCorrect:
							if areVectorsDifferentDebugging(newUpRight,node.probVectUpRight):
								print("new probVectUpRight at node is different from the old one, and it shouldn't be.")
								print(node.children)
								print(node.dist)
								print(node.children[1].dist)
								print()
								print(newUpRight)
								print()
								print(node.probVectUpRight)
								print()
								print(vectUp)
								print()
								print(node.children[1].probVect)
								exit()
						node.probVectUpRight=newUpRight
					newUpLeft=mergeVectorsUpDown(vectUp,node.dist,node.children[0].probVect,node.children[0].dist,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices)
					if newUpLeft==None:
						if(not node.children[0].dist) and (not node.dist):
							updateBLen(node.children[0],mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices)
							if not node.children[0].dist:
								updateBLen(node,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices)
								totNodeList.append((node.up,nodeChildNum))
								node.probVectTotUp=mergeVectorsUpDown(vectUp,node.dist/2,node.probVect,node.dist/2,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices)
								node.probVectUpRight=mergeVectorsUpDown(vectUp,node.dist,node.children[1].probVect,node.children[1].dist,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices)
							else:
								totNodeList.append((node,0))
							node.probVectUpLeft=mergeVectorsUpDown(vectUp,node.dist,node.children[0].probVect,node.children[0].dist,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices)
						else:
							print("Strange, distances>0 but inconsistent upLeft genome list creation in reCalculateAllGenomeLists()")
							exit()
					else:
						if checkExistingAreCorrect:
							if areVectorsDifferentDebugging(newUpLeft,node.probVectUpLeft):
								print("new probVectUpLeft at node is different from the old one, and it shouldn't be.")
								print(node.children)
								print(node.dist)
								print(node.children[0].dist)
								print(newUpLeft)
								print()
								print(node.probVectUpLeft)
								print()
								print(vectUp)
								print()
								print(node.children[0].probVect)
								exit()
						node.probVectUpLeft=newUpLeft
					node=node.children[0]
				else:
					lastNode=node
					node=node.up
					direction=1
			else:
				if lastNode==node.children[0]:
					node=node.children[1]
					direction=0
				else:
					lastNode=node
					node=node.up
					direction=1

		updatePartials(totNodeList,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices)





#Initialize the mutation rate matrix
#If a fixed rate matrix is needed for SARS-CoV-2, this is the nucleotide mutation rate matrix from De Maio et al 2021
#mutMatrix=[[0.0,0.039,0.310,0.123],[0.140,0.0,0.022,3.028],[0.747,0.113,0.0,2.953],[0.056,0.261,0.036,0.0]]
pseudoMutCounts=[[0.0,1.0,5.0,2.0],[2.0,0.0,1.0,40.0],[5.0,2.0,0.0,20.0],[2.0,3.0,1.0,0.0]]
if model=="JC":
	mutMatrix=[[-1.0,1.0/3,1.0/3,1.0/3],[1.0/3,-1.0,1.0/3,1.0/3],[1.0/3,1.0/3,-1.0,1.0/3],[1.0/3,1.0/3,1.0/3,-1.0]]
else:
	mutMatrix=[[0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0]]
	updateSubMatrix(pseudoMutCounts,model,mutMatrix)

nonMutRates=[0,0,0,0]
for i in range4:
	nonMutRates[i]=mutMatrix[i][i]

if runFluTest:
	mutMatrix  =[[-0.8397839380152881, 0.1478749673910128, 0.6118225025955456, 0.08008646802872968], [0.2585686827322552, -1.435769845796703, 0.021540539923035597, 1.1556606231414122], [0.8457231374881425, 0.01702856651052903, -0.9138505913154702, 0.05109888731679858], [0.09943753682215804, 0.8206165504318554, 0.04589864812565988, -0.9659527353796733]]
	nodesToVisit=[tree1]
	while nodesToVisit:
		node=nodesToVisit.pop()
		if node.children:
			nodesToVisit.append(node.children[0])
			nodesToVisit.append(node.children[1])
		else:
			if not node.name in data.keys():
				data[node.name]=[]




cumulativeRate=[0.0]
for i in range(lRef):
	ind=refIndeces[i]
	cumulativeRate.append(cumulativeRate[-1]+nonMutRates[ind])



#In case an input tree is given, calculate all genome lists, pseudocounts and rates, and then recalculate genome lists again according to the new rates.
if inputTree!="" and (not runFluTest):
	reCalculateAllGenomeLists(tree1,mutMatrix,countPseudoCounts=True,pseudoMutCounts=pseudoMutCounts,data=data,firstSetUp=True)
	if model!="JC" and updateSubMatrix(pseudoMutCounts,model,mutMatrix):
		for i in range(lRef):
			cumulativeRate[i+1]=cumulativeRate[i]+nonMutRates[refIndeces[i]]
		for i in range4:
			nonMutRates[i]=mutMatrix[i][i]
	reCalculateAllGenomeLists(tree1,mutMatrix)
	print("Genome list for initial tree and initial pseudocounts calculated.")
if runFluTest:
	reCalculateAllGenomeLists(tree1,mutMatrix,data=data)


#Sort samples based on distance from reference, but punishing more isolated N's and ambiguity characters.
#more ambiguous sequences are placed last this way - this is useful since ambiguous sequences are harder to place and are more likely to be less informative (and so be removed from the analysis altogether)
#def distancesFromRefPunishNs(data,samples):
def distancesFromRefPunishNs(data,samples=None):
	sampleDistances=[]
	if samples==None:
		rangeInd=range(len(data))
	else:
		rangeInd=samples
	#for sample in samples:
	for diffIndex in rangeInd:
		diffs=data[diffIndex]
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
				if len(m)>2:
					pos=currPos+m[2]
				else:
					pos=currPos+1
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
		sampleDistances.append((diffNum*1000+lRef-comparisons,diffIndex))
		totDivFromRef[0]+=float(diffNum)/comparisons
		if ((float(diffNum)/comparisons)>0.1) and (not warnedTotDiv[0]):
			warnedTotDiv[0]=True
			print("\n WARNING!!!!!!!!!!!\n\n Found sample with divergence "+str(float(diffNum)/comparisons)+" from the reference ; at high divergence, MAPLE will struggle both in terms of accuracy and computational demand, so it might make sense to use a traditional phylogenetic approach instead.\n\n End of warning.\n")

	from operator import itemgetter
	print("Now doing sorting")
	sampleDistances.sort(reverse=True,key=itemgetter(0))
	return sampleDistances



#function to check is one sequence is less informative than another;
#returns 0 when the 2 sequences are not comparable, otherwise returns 1 if the first is more informative or if they are identical, and 2 otherwise.
def isMinorSequence(probVect1,probVect2):
	indexEntry1, indexEntry2, pos = 0, 0, 0
	entry1=probVect1[indexEntry1]
	entry2=probVect2[indexEntry2]
	found1bigger=False
	found2bigger=False
	while True:
		if entry1[0]!=entry2[0]:
			if entry1[0]==5:
				pos=min(entry1[1],entry2[1])
				found2bigger=True
			elif entry2[0]==5:
				pos=min(entry1[1],entry2[1])
				found1bigger=True
			elif entry1[0]==6:
				if entry2[0]==4:
					i2=refIndeces[pos]
				else:
					i2=entry2[0]
				if entry1[-1][i2]>0.1:
					found2bigger=True
				else:
					return 0
				pos+=1
			elif entry2[0]==6:
				if entry1[0]==4:
					i1=refIndeces[pos]
				else:
					i1=entry1[0]
				if entry2[-1][i1]>0.1:
					found1bigger=True
				else:
					return 0
				pos+=1
			else:
				return 0
		elif entry1[0]==6:
			for j in range4:
				if entry2[-1][j]>0.1 and entry1[-1][j]<0.1:
					found1bigger=True
				elif entry1[-1][j]>0.1 and entry2[-1][j]<0.1:
					found2bigger=True
			pos+=1
		else:
			pos=min(entry1[1],entry2[1])
		if found1bigger and found2bigger:
			return 0
		if pos==lRef:
			break
		if pos==entry1[1]:
			indexEntry1+=1
			entry1=probVect1[indexEntry1]
		if pos==entry2[1]:
			indexEntry2+=1
			entry2=probVect2[indexEntry2]

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






#function to calculate likelihood cost of appending sample to parent node
#UNNECESSARY this can be considered an optimized version of function appendProbNode(), although appendProb() is likely marginally faster.
#as such, when rewriting the code, one may get rid of this function altogether and use just appendProbNode() instead to reduce amount of code required.
def appendProb(probVectP,probVectC,bLen,mutMatrix):
	if not bLen:
		bLen=0.0
	Lkcost, indexEntry1, indexEntry2, totalFactor, pos = 0.0, 0, 0, 1.0, 0
	entry1=probVectP[indexEntry1]
	entry2=probVectC[indexEntry2]
	end=min(entry1[1],entry2[1])
	contribLength=bLen
	while True:
		if entry2[0]==5: # case N
			pos=min(entry1[1],entry2[1])
		elif entry1[0]==5: # case N
			#if parent node is type "N", in theory we might have to calculate the contribution of root nucleotides; 
			# however, if this node is "N" then every other node in the current tree is "N", so we can ignore this since this contribution cancels out in relative terms.
			pos=min(entry1[1],entry2[1])
		elif entry1[0]==4: # case entry1 is R
			if entry2[0]==4:
				end=min(entry1[1],entry2[1])
				if len(entry1)==2:
					Lkcost+=bLen*(cumulativeRate[end]-cumulativeRate[pos])
				else:
					contribLength=bLen+entry1[2]
					if len(entry1)==3:
						Lkcost+=contribLength*(cumulativeRate[end]-cumulativeRate[pos])
					else:
						#here contribution from root frequency gets added and subtracted so it's ignored
						Lkcost+=(contribLength+entry1[3])*(cumulativeRate[end]-cumulativeRate[pos])
				pos=end
			elif entry2[0]==6:
				i1=refIndeces[pos]
				if len(entry1)==4:
					contribLength=bLen+entry1[3]
					if entry2[2][i1]>0.1:
						contribLength+=entry1[2]
						#here contribution from root frequency can also be also ignored
						Lkcost+=nonMutRates[i1]*contribLength 
					else:
						tot=0.0
						for i in range4:
							if i1==i:
								tot2=rootFreqs[i]*(1.0+nonMutRates[i]*entry1[2])
							else:
								tot2=rootFreqs[i]*mutMatrix[i][i1]*entry1[2]
							tot3=0.0
							for j in range4:
								if entry2[2][j]>0.1:
									tot3+=mutMatrix[i][j]
							tot3*=contribLength
							if entry2[2][i]>0.1:
								tot3+=1.0
							tot+=tot2*tot3
						totalFactor*=(tot/rootFreqs[i1])
				else:
					if entry2[2][i1]>0.1:
						if len(entry1)==3:
							Lkcost+=nonMutRates[i1]*(bLen+entry1[2])
						else:
							Lkcost+=nonMutRates[i1]*bLen
					else:
						tot=0.0
						for j in range4:
							if entry2[2][j]>0.1:
								tot+=mutMatrix[i1][j]
						if len(entry1)==3:
							totalFactor*=tot*(bLen+entry1[2])
						else:
							totalFactor*=tot*bLen
				pos+=1

			else: #entry1 is reference and entry2 is a different but single nucleotide
				if len(entry1)==2:
					totalFactor*=mutMatrix[refIndeces[pos]][entry2[0]]*bLen
				elif len(entry1)==3:
					totalFactor*=mutMatrix[refIndeces[pos]][entry2[0]]*(bLen+entry1[2])
				else:
					i1=refIndeces[pos]
					i2=entry2[0]
					totalFactor*=((rootFreqs[i1]*mutMatrix[i1][i2]*(bLen+entry1[3])*(1.0+nonMutRates[i1]*entry1[2])+rootFreqs[i2]*mutMatrix[i2][i1]*entry1[2]*(1.0+nonMutRates[i2]*(bLen+entry1[3])))/rootFreqs[i1])
				pos+=1

		# entry1 is of type "O"
		elif entry1[0]==6:
			if len(entry1)==3:
				bLen13=bLen
			else:
				bLen13=bLen+entry1[2]
			if entry2[0]==6:
				tot=0.0
				for j in range4:
					tot2=0.0
					for j2 in range4:
						if entry2[2][j2]>0.1:
							tot2+=mutMatrix[j][j2]
					tot2*=bLen13
					if entry2[2][j]>0.1:
						tot2+=1.0
					tot+=tot2*entry1[-1][j]
				totalFactor*=tot
			else:
				if entry2[0]==4:
					i2=refIndeces[pos]
				else:
					i2=entry2[0]
				totalFactor*=(entry1[-1][i2]+bLen13*(entry1[-1][0]*mutMatrix[0][i2]+entry1[-1][1]*mutMatrix[1][i2]+entry1[-1][2]*mutMatrix[2][i2]+entry1[-1][3]*mutMatrix[3][i2]))
			pos+=1

		else: #entry1 is a non-ref nuc
			i1=entry1[0]
			if entry2[0]==i1:
				if len(entry1)==2:
					Lkcost+=nonMutRates[i1]*bLen
				elif len(entry1)==3:
					Lkcost+=nonMutRates[i1]*(bLen+entry1[2])
				else:
					Lkcost+=nonMutRates[i1]*(bLen+entry1[2]+entry1[3])
			else: #entry1 and entry2 are of different types
				if entry2[0]==6:
					if len(entry1)==4:
						bLen15=bLen+entry1[3]
						if entry2[2][i1]>0.1:
							Lkcost+=nonMutRates[i1]*(bLen15+entry1[2])
						else:
							tot=0.0
							for i in range4:
								if i1==i:
									tot2=rootFreqs[i]*(1.0+nonMutRates[i1]*entry1[2])
								else:
									tot2=rootFreqs[i]*mutMatrix[i][i1]*entry1[2]
								tot3=0.0
								for j in range4:
									if entry2[2][j]>0.1:
										tot3+=mutMatrix[i][j]
								if entry2[2][i]>0.1:
									tot+=tot2*(1.0+bLen15*tot3)
								else:
									tot+=tot2*bLen15*tot3
							totalFactor*=(tot/rootFreqs[i1])
					else:
						if entry2[2][i1]>0.1:
							if len(entry1)==2:
								Lkcost+=nonMutRates[i1]*bLen
							else:
								Lkcost+=nonMutRates[i1]*(bLen+entry1[2])
						else:
							tot=0.0
							for j in range4:
								if entry2[2][j]>0.1:
									tot+=mutMatrix[i1][j]
							if len(entry1)==2:
								totalFactor*=tot*bLen
							else:
								totalFactor*=tot*(bLen+entry1[2])
				#entry2 is a nucleotide type (like entry1)		
				else:
					if entry2[0]==4:
						i2=refIndeces[pos]
					else:
						i2=entry2[0]
					if len(entry1)==2:
						totalFactor*=mutMatrix[i1][i2]*bLen
					elif len(entry1)==3:
						totalFactor*=mutMatrix[i1][i2]*(bLen+entry1[2])
					else:
						#here we ignore contribution of non-parsimonious mutational histories
						totalFactor*=((rootFreqs[i1]*mutMatrix[i1][i2]*(bLen+entry1[3])*(1.0+nonMutRates[i1]*entry1[2]) + rootFreqs[i2]*mutMatrix[i2][i1]*entry1[2]*(1.0+nonMutRates[i2]*(bLen+entry1[3])) )/rootFreqs[i1])						
			pos+=1

		if totalFactor<=minimumCarryOver:
			try:
				if totalFactor<sys.float_info.min:
					return float("-inf")
			except:
				print("In appendProb() value error")
				print(totalFactor)
				return float("-inf")
			Lkcost+=log(totalFactor)
			totalFactor=1.0

		if pos==lRef:
			break
		if pos==entry1[1]:
			indexEntry1+=1
			entry1=probVectP[indexEntry1]
		if pos==entry2[1]:
			indexEntry2+=1
			entry2=probVectC[indexEntry2]
	return Lkcost+log(totalFactor)





#number of samples that could have been placed as major of another sample but weren't due to sample placement order
totalMissedMinors=[0]








#function traversing the tree to find the best node in the tree where to re-append the given subtree (rooted at node.children[child]) to improve the topology of the current tree.
# bestLKdiff is the best likelihood cost found for the current placement (after optimizing the branch length).
# removedBLen is such branch length that optimizes the current placement - it will be used to place the subtree attached at other nodes of the tree.
def findBestParentTopology(node,child,bestLKdiff,removedBLen,mutMatrix,strictTopologyStopRules=strictTopologyStopRules,allowedFailsTopology=allowedFailsTopology,thresholdLogLKtopology=thresholdLogLKtopology,useRateVariation=False,mutMatrices=None):
	bestNode=node.children[1-child]
	bestNodes=[]
	nodesToVisit=[]
	removedPartials=node.children[child].probVect
	originalLK=bestLKdiff
	originalPlacement=bestNode
	if debugging:
		print("findBestParentTopology() with parameters thresholdLogLKtopology "+str(thresholdLogLKtopology))
	if node.up!=None:
		if node.up.children[0]==node:
			childUp=1
			vectUpUp=node.up.probVectUpRight
		else:
			childUp=2
			vectUpUp=node.up.probVectUpLeft
		#the list nodesToVisit keeps trak of the nodes of the tree we wil need to traverse, (first element of each entry),
		# the direction from where we are visitng them (0=from parent, and 1,2=from one of the children),
		# an updated genome list from the direction where we come from (taking into account the removal of the given subtree),
		# a branch length value separating the node from this updated genome list (useful for the fact that the removal of the subtree changes the branch length at the removal node),
		# a flag that says if the updated genome list passed needs still updating, or if it has become identical to the pre-existing genome list in the tree (which usually happens after a while),
		# the likelihood cost of appending at the last node encountered in this direction,
		# a number of consecutively failed traversal steps since the last improvement found (if this number goes beyond a threshold, traversal in the considered direction might be halted).
		nodesToVisit.append((node.up,childUp,node.children[1-child].probVect,node.children[1-child].dist+node.dist,True,bestLKdiff,0))
		nodesToVisit.append((node.children[1-child],0,vectUpUp,node.children[1-child].dist+node.dist,True,bestLKdiff,0))
		originalBLens=(node.dist,node.children[1-child].dist,removedBLen)
	else:
		# case node is root
		if node.children[1-child].children: # case there is only one sample outside of the subtree doesn't need to be considered
			child1=node.children[1-child].children[0]
			child2=node.children[1-child].children[1]
			vectUp1=rootVector(child2.probVect,child2.dist,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices)
			nodesToVisit.append((child1,0,vectUp1,child1.dist,True,bestLKdiff,0))
			vectUp2=rootVector(child1.probVect,child1.dist,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices)
			nodesToVisit.append((child2,0,vectUp2,child2.dist,True,bestLKdiff,0))
			originalPlacement=node.children[1-child].children[0]
			originalBLens=(0.0,node.children[1-child].children[0].dist,removedBLen)
		else:
			originalBLens=(0.0,node.children[1-child].dist,removedBLen)


	while nodesToVisit:
		t1,direction,passedPartials,distance,needsUpdating,lastLK,failedPasses=nodesToVisit.pop()
		if direction==0:
			#consider the case we are moving from a parent to a child
			if t1.dist and (not (t1.up==node or t1.up==None)):
				if needsUpdating:
					midTot=mergeVectorsUpDown(passedPartials,distance/2,t1.probVect,distance/2,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices)
					if not areVectorsDifferent(midTot,t1.probVectTotUp):
						needsUpdating=False
				else:
					midTot=t1.probVectTotUp
				if midTot==None:
					continue
				midProb=appendProbNode(midTot,removedPartials,removedBLen,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices)
				if midProb>bestLKdiff-thresholdLogLKoptimizationTopology:
					#if needsUpdating, then add to the tuple also the information on the up and down genome lists to use to recalculate intermediate genome lists at varying branch lengths
					if needsUpdating:
						bestNodes.append((t1,midProb,passedPartials,t1.probVect,distance,midTot))
					else:
						bestNodes.append((t1,midProb))
				if midProb>bestLKdiff:
					bestLKdiff=midProb
					bestNode=t1
					failedPasses=0
				if midProb<(lastLK-thresholdLogLKconsecutivePlacement): #placement at current node is considered failed if placement likelihood is not improved by a certain margin compared to best placement so far for the nodes above it.
					failedPasses+=1
			else:
				midProb=lastLK
			
			#keep crawling down into children nodes unless the stop criteria for the traversal are satisfied.
			traverseChildren=False
			if strictTopologyStopRules:
				if failedPasses<=allowedFailsTopology and midProb>(bestLKdiff-thresholdLogLKtopology) and t1.children:
					traverseChildren=True
			else:
				if failedPasses<=allowedFailsTopology or midProb>(bestLKdiff-thresholdLogLKtopology):
					if t1.children:
						traverseChildren=True
			if traverseChildren:
				child=t1.children[0]
				otherChild=t1.children[1]
				if needsUpdating:
					vectUpRight=mergeVectorsUpDown(passedPartials,distance,otherChild.probVect,otherChild.dist,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices)
				else:
					vectUpRight=t1.probVectUpRight
				if vectUpRight!=None:
					nodesToVisit.append((child,0,vectUpRight,child.dist,needsUpdating,midProb,failedPasses))
				child=t1.children[1]
				otherChild=t1.children[0]
				if needsUpdating:
					vectUpLeft=mergeVectorsUpDown(passedPartials,distance,otherChild.probVect,otherChild.dist,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices)
				else:
					vectUpLeft=t1.probVectUpLeft
				if vectUpLeft!=None:
					nodesToVisit.append((child,0,vectUpLeft,child.dist,needsUpdating,midProb,failedPasses))

		else: #case when crawling up from child to parent
			otherChild=t1.children[2-direction]
			midBottom=None
			if t1.dist and t1.up!=None: #try appending mid-branch
				if needsUpdating:
					midBottom=mergeVectors(otherChild.probVect,otherChild.dist,passedPartials,distance,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices)
					if midBottom==None:
						continue
					if t1==t1.up.children[0]:
						vectUp=t1.up.probVectUpRight
					else:
						vectUp=t1.up.probVectUpLeft
					midTot=mergeVectorsUpDown(vectUp,t1.dist/2,midBottom,t1.dist/2,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices)
					if not areVectorsDifferent(midTot,t1.probVectTotUp):
						needsUpdating=False
				else:
					midTot=t1.probVectTotUp
				if midTot==None:
					continue
				midProb=appendProbNode(midTot,removedPartials,removedBLen,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices)
				if midProb>bestLKdiff:
					bestLKdiff=midProb
					bestNode=t1
					failedPasses=0
				if midProb>=(bestLKdiff-thresholdLogLKoptimizationTopology):
					#if needsUpdating, then add to the tuple also the information on the up and down genome lists to use to recalculate intermediate genome lists at varying branch lengths
					if needsUpdating:
						bestNodes.append((t1,midProb,vectUp,midBottom,t1.dist,midTot))
					else:
						bestNodes.append((t1,midProb))
				if midProb<(lastLK-thresholdLogLKconsecutivePlacement): #placement at current node is considered failed if placement likelihood is not improved by a certain margin compared to best placement so far for the nodes above it.
					failedPasses+=1
			else:
				midProb=lastLK
			
			#testing for the stop rule of the traversal process
			keepTraversing=False
			if strictTopologyStopRules:
				if failedPasses<=allowedFailsTopology and midProb>(bestLKdiff-thresholdLogLKtopology):
					keepTraversing=True
			else:
				if failedPasses<=allowedFailsTopology or midProb>(bestLKdiff-thresholdLogLKtopology):
					keepTraversing=True
			if keepTraversing:
				# keep crawling up into parent and sibling node
				if t1.up!=None: #case the node is not the root
					#first pass the crawling down the other child
					if t1==t1.up.children[0]:
						upChild=0
						if needsUpdating:
							vectUpUp=t1.up.probVectUpRight
					else:
						upChild=1
						if needsUpdating:
							vectUpUp=t1.up.probVectUpLeft
					if needsUpdating:
						vectUp=mergeVectorsUpDown(vectUpUp,t1.dist,passedPartials,distance,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices)
					else:
						if direction==1:
							vectUp=t1.probVectUpLeft
						else:
							vectUp=t1.probVectUpRight

					if vectUp==None:
						continue
					else:
						nodesToVisit.append((otherChild,0,vectUp,otherChild.dist,needsUpdating,midProb,failedPasses))
					#now pass the crawling up to the parent node
					if needsUpdating:
						if midBottom==None:
							midBottom=mergeVectors(otherChild.probVect,otherChild.dist,passedPartials,distance,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices)
							if midBottom==None:
								continue
					else:
						midBottom=t1.probVect
					nodesToVisit.append((t1.up,upChild+1,midBottom,t1.dist,needsUpdating,midProb,failedPasses))
				#now consider case of root node
				else:
					if needsUpdating:
						vectUp=rootVector(passedPartials,distance,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices)
					else:
						if direction==1:
							vectUp=t1.probVectUpLeft
						else:
							vectUp=t1.probVectUpRight
					nodesToVisit.append((otherChild,0,vectUp,otherChild.dist,needsUpdating,midProb,failedPasses))

	if debugging:
		print("Within findBestParentTopology(), bestNodes found dunring the first pass (no branch length optimization): ")
		print(bestNodes)

	#Initial exploration is finished.
	#Now, for each branch within threshold likelihood distance from the best found, optimize branch lengths.
	#Use optimized scores to select final best branch
	#TODO Myrthe, changed this
	bestBranchLengths=originalBLens
	bestScore=bestLKdiff
	compensanteForBranchLengthChange=True
	secondBranchLengthOptimizationRound=False
	if not bestNodes:
		return originalPlacement, originalLK ,originalBLens
	BLenHaveBeenOptimized=False
	for nodePair in bestNodes:
		score=nodePair[1]
		if score>=bestLKdiff-thresholdLogLKoptimizationTopology:
			t1=nodePair[0]
			if len(nodePair)==2:
				#optimize branch lengths of appendage
				if t1==t1.up.children[0]:
					upVect=t1.up.probVectUpRight
				else:
					upVect=t1.up.probVectUpLeft
				downVect=t1.probVect
				distance=t1.dist
				midTot=t1.probVectTotUp
			else:
				upVect=nodePair[2]
				downVect=nodePair[3]
				distance=nodePair[4]
				midTot=nodePair[5]
			bestAppendingLength=estimateBranchLengthWithDerivative(midTot,removedPartials,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices)
			#now optimize appending location
			midLowerVector=mergeVectors(downVect,distance/2,removedPartials,bestAppendingLength,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices)
			bestTopLength=estimateBranchLengthWithDerivative(upVect,midLowerVector,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices)
			midTopVector=mergeVectorsUpDown(upVect,bestTopLength,removedPartials,bestAppendingLength,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices)
			bestBottomLength=estimateBranchLengthWithDerivative(midTopVector,downVect,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices)
			newMidVector=mergeVectorsUpDown(upVect,bestTopLength,downVect,bestBottomLength,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices)
			if secondBranchLengthOptimizationRound: #if wanted, do a second round of branch length optimization
				bestAppendingLength=estimateBranchLengthWithDerivative(newMidVector,removedPartials,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices)
				midLowerVector=mergeVectors(downVect,bestBottomLength,removedPartials,bestAppendingLength,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices)
				bestTopLength=estimateBranchLengthWithDerivative(upVect,midLowerVector,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices)
				midTopVector=mergeVectorsUpDown(upVect,bestTopLength,removedPartials,bestAppendingLength,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices)
				bestBottomLength=estimateBranchLengthWithDerivative(midTopVector,downVect,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices)
				newMidVector=mergeVectorsUpDown(upVect,bestTopLength,downVect,bestBottomLength,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices)
			appendingCost=appendProbNode(newMidVector,removedPartials,bestAppendingLength,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices)
			if debugging:
				print("Within findBestParentTopology(), second pass (branch length optimization), found appendingCost "+str(appendingCost)+" and branch lengths "+str(bestTopLength)+" "+str(bestBottomLength)+" "+str(bestAppendingLength))
			if compensanteForBranchLengthChange: #if wanted, do a more thorough examination of the appending cost, taking into account a possible change in branch length of the branch on which to be appended.
				initialCost=appendProbNode(upVect,downVect,distance,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices)
				newPartialCost=appendProbNode(upVect,downVect,bestBottomLength+bestTopLength,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices)
				optimizedScore=appendingCost+newPartialCost-initialCost
				if debugging:
					print("Within findBestParentTopology(), initialCost "+str(initialCost)+" for branch length "+str(distance)+", newPartialCost "+str(newPartialCost)+" for distance "+str(bestBottomLength+bestTopLength)+", leads to optimizedScore "+str(optimizedScore))
			else:
				optimizedScore=appendingCost
			if optimizedScore>=bestScore:
				BLenHaveBeenOptimized=True
				bestNode=t1
				bestScore=optimizedScore
				bestBranchLengths=(bestTopLength,bestBottomLength,bestAppendingLength)
				testNonsensicalBLens=True
				if testNonsensicalBLens:
					if (not bestBottomLength) and (not bestAppendingLength):
						midLowerVector=mergeVectors(downVect,bestBottomLength,removedPartials,bestAppendingLength,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices)
						if midLowerVector==None:
							print("Problem in findBestParentTopology(), midLowerVector is None")
							print(downVect)
							print(bestBottomLength)
							print(removedPartials)
							print(bestAppendingLength)
							print(upVect)
							print(bestTopLength)
							print("AppendingCost "+str(appendingCost)+" initialCost "+str(initialCost)+" newPartialCost "+str(newPartialCost)+" optimizedScore "+str(optimizedScore))
							exit()
					if (not bestBottomLength) and (not bestTopLength):
						newMidVector=mergeVectorsUpDown(upVect,bestTopLength,downVect,bestBottomLength,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices)
						if newMidVector==None:
							print("Problem in findBestParentTopology(), newMidVector is None")
							print(downVect)
							print(bestBottomLength)
							print(removedPartials)
							print(bestAppendingLength)
							print(upVect)
							print(bestTopLength)
							print("AppendingCost "+str(appendingCost)+" initialCost "+str(initialCost)+" newPartialCost "+str(newPartialCost)+" optimizedScore "+str(optimizedScore))
							exit()
					if (not bestAppendingLength) and (not bestTopLength):
						midTopVector=mergeVectorsUpDown(upVect,bestTopLength,removedPartials,bestAppendingLength,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices)
						if midTopVector==None:
							print("Problem in findBestParentTopology(), midTopVector is None")
							print(downVect)
							print(bestBottomLength)
							print(removedPartials)
							print(bestAppendingLength)
							print(upVect)
							print(bestTopLength)
							print("AppendingCost "+str(appendingCost)+" initialCost "+str(initialCost)+" newPartialCost "+str(newPartialCost)+" optimizedScore "+str(optimizedScore))
							exit()
	#TODO Myrthe, added this
	if not BLenHaveBeenOptimized:
		bestBranchLengths=(bestNode.dist/2,bestNode.dist/2,removedBLen)
	return bestNode, bestScore, bestBranchLengths












#function to find the best node in the tree where to append the new sample; traverses the tree and tries to append the sample at each node and mid-branch nodes, 
# but stops traversing when certain criteria are met.
def findBestParentForNewSample(root,diffs,sample,mutMatrix):
	bestNodes=[]
	bestNode=root
	bestBranchLengths=(False,False,oneMutBLen)
	if not root.children: #check if the new leaf is strictly less informative than already placed leaf
		comparison=isMinorSequence(root.probVect,diffs)
		if comparison==1:
			root.minorSequences.append(sample)
			return root, 1.0, None
		elif comparison==2:
			totalMissedMinors[0]+=1
	rootVect=rootVector(root.probVect,False,mutMatrix)
	bestLKdiff=appendProb(rootVect,diffs,oneMutBLen,mutMatrix)
	nodesToVisit=[]
	for child in root.children:
		nodesToVisit.append((child,bestLKdiff,0))
	while nodesToVisit:
		t1,parentLK,failedPasses=nodesToVisit.pop()
		if not t1.children: #check if the new leaf is strictly less informative than already placed leaf
			comparison=isMinorSequence(t1.probVect,diffs)
			if comparison==1:
				t1.minorSequences.append(sample)
				return t1, 1.0, None
			elif comparison==2:
				totalMissedMinors[0]+=1

		if t1.dist and t1.up!=None: # try first placing as a descendant of the mid-branch point of the branch above the current node.
			LKdiff=appendProb(t1.probVectTotUp,diffs,oneMutBLen,mutMatrix)
			if LKdiff>=bestLKdiff:
				bestLKdiff=LKdiff
				bestNode=t1
				failedPasses=0
				bestNodes.append((t1,LKdiff))
			elif LKdiff>bestLKdiff-thresholdLogLKoptimization:
				bestNodes.append((t1,LKdiff))
			if LKdiff<(parentLK-thresholdLogLKconsecutivePlacement): #placement at current node is considered failed if placement likelihood is not improved by a certain margin compared to best placement so far for the nodes above it.
				failedPasses+=1
		else:
			LKdiff=parentLK
		#keep trying to place at children nodes, unless placement has failed too many times already or the likelihood cost suboptimality is above a certain threshold
		if strictInitialStopRules:
			if failedPasses<=allowedFails and LKdiff>(bestLKdiff-thresholdLogLK): 
				for c in t1.children:
					nodesToVisit.append((c,LKdiff,failedPasses))
		else:
			if failedPasses<=allowedFails or LKdiff>(bestLKdiff-thresholdLogLK): 
				for c in t1.children:
					nodesToVisit.append((c,LKdiff,failedPasses))
	#Initial exploration is finished.
	#Now, for each branch within threshold likelihood distance from the best found, optimize branch lengths.
	#Use optimized scores to select final best branch
	if bestNode!=root:
		bestBranchLengths=(bestNode.dist/2,bestNode.dist/2,oneMutBLen)
	bestScore=bestLKdiff
	compensanteForBranchLengthChange=True
	secondBranchLengthOptimizationRound=False
	for nodePair in bestNodes:
		score=nodePair[1]
		if score>=bestLKdiff-thresholdLogLKoptimization:
			node=nodePair[0]
			#optimize branch lengths of appendage
			if node==node.up.children[0]:
				upVect=node.up.probVectUpRight
			else:
				upVect=node.up.probVectUpLeft
			bestAppendingLength=estimateBranchLengthWithDerivative(node.probVectTotUp,diffs,mutMatrix)
			#now optimize appending location
			midLowerVector=mergeVectors(node.probVect,node.dist/2,diffs,bestAppendingLength,mutMatrix)
			bestTopLength=estimateBranchLengthWithDerivative(upVect,midLowerVector,mutMatrix)
			midTopVector=mergeVectorsUpDown(upVect,bestTopLength,diffs,bestAppendingLength,mutMatrix)
			bestBottomLength=estimateBranchLengthWithDerivative(midTopVector,node.probVect,mutMatrix)
			newMidVector=mergeVectorsUpDown(upVect,bestTopLength,node.probVect,bestBottomLength,mutMatrix)
			if secondBranchLengthOptimizationRound: #if wanted, do a second round of branch length optimization
				bestAppendingLength=estimateBranchLengthWithDerivative(newMidVector,diffs,mutMatrix)
				midLowerVector=mergeVectors(node.probVect,bestBottomLength,diffs,bestAppendingLength,mutMatrix)
				bestTopLength=estimateBranchLengthWithDerivative(upVect,midLowerVector,mutMatrix)
				midTopVector=mergeVectorsUpDown(upVect,bestTopLength,diffs,bestAppendingLength,mutMatrix)
				bestBottomLength=estimateBranchLengthWithDerivative(midTopVector,node.probVect,mutMatrix)
				newMidVector=mergeVectorsUpDown(upVect,bestTopLength,node.probVect,bestBottomLength,mutMatrix)
			appendingCost=appendProb(newMidVector,diffs,bestAppendingLength,mutMatrix)
			if compensanteForBranchLengthChange: #if wante, do a more thorough examination of the appending cost, taking into account a possible change in branch length of the branch on which to be appended.
				initialCost=appendProbNode(upVect,node.probVect,node.dist,mutMatrix)
				newPartialCost=appendProbNode(upVect,node.probVect,bestBottomLength+bestTopLength,mutMatrix)
				optimizedScore=appendingCost+newPartialCost-initialCost
			else:
				optimizedScore=appendingCost
			if optimizedScore>=bestScore:
				bestNode=node
				bestScore=optimizedScore
				bestBranchLengths=(bestTopLength,bestBottomLength,bestAppendingLength)

	return bestNode, bestScore, bestBranchLengths







#we know that sample "sample", with partials "newPartials", is best placed near a node resulting in logLK contribution of newChildLK
# now explore exactly which position near node is best for placement (direct child, sibling or sibling of existing child), and which are the best branch lengths,
#then add the sample at that position of the tree, and update all the internal probability vectors.
#UNNECESSARY? Could probably just be replaced by the more general placeSubtreeOnTree().
def placeSampleOnTree(node,newPartials,sample,newChildLK, bestUpLength, bestDownLength, bestAppendingLength,mutMatrix,pseudoMutCounts):
	tryNewRoot=False
	if newChildLK<-0.01:
		sumChildLKs[0]+=newChildLK
		numChildLKs[0]+=1
	if node.up==None:
		tryNewRoot=True
		totRoot=rootVector(node.probVect,False,mutMatrix)
		bestAppendingLength=estimateBranchLengthWithDerivative(totRoot,newPartials,mutMatrix)
		root=node
		newChildLK=appendProb(totRoot,newPartials,bestAppendingLength,mutMatrix)
	else:
		if node.up.children[0]==node:
			child=0
			vectUp=node.up.probVectUpRight
		else:
			child=1
			vectUp=node.up.probVectUpLeft
		if not bestUpLength:
			pNode=node.up
			while (not pNode.dist) and (pNode.up!=None):
				pNode=pNode.up
			if pNode.up==None:
				root=pNode
				tryNewRoot=True
				#TODO Myrthe, added this for improvement, similar to the SPR case
				if (not bestDownLength) or (bestDownLength>1.01*node.dist) or (bestDownLength<0.99*node.dist):
					node.dist=bestDownLength
					nodeList=[(node,2),(node.up,child)]
					updatePartials(nodeList,mutMatrix)
	#in case of best placement as a descendant appended exactly at the root node, attempt also to create new root
	if tryNewRoot:
		node=root
		probOldRoot = findProbRoot(node.probVect)
		rootUpLeft=rootVector(node.probVect,bestAppendingLength/2,mutMatrix)
		bestRightLength=estimateBranchLengthWithDerivative(rootUpLeft,newPartials,mutMatrix)
		rootUpRight=rootVector(newPartials,bestRightLength,mutMatrix)
		bestLeftLength=estimateBranchLengthWithDerivative(rootUpRight,node.probVect,mutMatrix)
		secondBranchLengthOptimizationRound=True
		if secondBranchLengthOptimizationRound: #if wanted, do a second round of branch length optimization
			rootUpLeft=rootVector(node.probVect,bestLeftLength,mutMatrix)
			bestRightLength=estimateBranchLengthWithDerivative(rootUpLeft,newPartials,mutMatrix)
			rootUpRight=rootVector(newPartials,bestRightLength,mutMatrix)
			bestLeftLength=estimateBranchLengthWithDerivative(rootUpRight,node.probVect,mutMatrix)
		probVectRoot,probRoot = mergeVectors(node.probVect,bestLeftLength,newPartials,bestRightLength,mutMatrix,returnLK=True)
		probRoot+= findProbRoot(probVectRoot)
		parentLKdiff=probRoot-probOldRoot
		if parentLKdiff<=newChildLK: #best is just placing as descendant of the root
			bestRightLength=bestAppendingLength
			bestLeftLength=False
			probVectRoot=mergeVectors(node.probVect,bestLeftLength,newPartials,bestRightLength,mutMatrix)
			rootUpRight=rootVector(newPartials,bestRightLength,mutMatrix)
		#now add new root to the tree
		newRoot=Tree()
		if probVectRoot==None:
			print("Issue with new root probVect inside placeSampleOnTree()")
			exit()
		#newRoot.replacements=0
		newRoot.probVect=probVectRoot
		newRoot.probVectUpRight=rootUpRight
		newRoot.probVectUpLeft=rootVector(node.probVect,bestLeftLength,mutMatrix)
		node.up=newRoot
		node.dist=bestLeftLength
		newRoot.add_child(node)
		newNode=Tree(name=sample,dist=bestRightLength)
		if bestRightLength>0.01 and (not warnedBLen[0]):
			warnedBLen[0]=True
			print("\n WARNING!!!!!!!!!!!\n\n Found branch of length "+str(bestRightLength)+" ; at high divergence, MAPLE will struggle both in terms of accuracy and computational demand, so it might make sense to use a traditional phylogenetic approach.\n\n End of warning.\n")
		#newNode.replacements=0
		newNode.minorSequences=[]
		newNode.up=newRoot
		newRoot.add_child(newNode)
		newNode.probVect=newPartials
		if bestRightLength:
			newNode.probVectTotUp=mergeVectorsUpDown(newRoot.probVectUpLeft,bestRightLength/2,newPartials,bestRightLength/2,mutMatrix)
		if verbose:
			print("new root added to tree")
			print(newRoot.probVect)
			print(newRoot.children[0].probVect)
			print(newNode.probVect)
		nodeList=[(node,2)]
		updatePartials(nodeList,mutMatrix)
		return newRoot

	#print("adding internal node")
	#in all other cases (not attempting to add a new root) create a new internal node in the tree and add sample as a descendant.
	newInternalNode=Tree()
	#newInternalNode.replacements=0
	node.up.children[child]=newInternalNode
	newInternalNode.up=node.up
	newInternalNode.add_child(node)
	node.up=newInternalNode
	node.dist=bestDownLength
	newNode=Tree(name=sample,dist=bestAppendingLength)
	if bestAppendingLength>0.01 and (not warnedBLen[0]):
		warnedBLen[0]=True
		print("\n WARNING!!!!!!!!!!!\n\n Found branch of length "+str(bestAppendingLength)+" ; at high divergence, MAPLE will struggle both in terms of accuracy and computational demand, so it might make sense to use a more traditional phylogenetic approach (e.g. FastTree, IQtree or RAxML).\n\n End of warning.\n")
	#newNode.replacements=0
	newNode.minorSequences=[]
	newNode.up=newInternalNode
	newInternalNode.add_child(newNode)
	newInternalNode.dist=bestUpLength
	newNode.probVect=newPartials
	newInternalNode.probVect=mergeVectors(node.probVect,bestDownLength,newPartials,bestAppendingLength,mutMatrix)
	newInternalNode.probVectUpRight=mergeVectorsUpDown(vectUp,bestUpLength,newPartials,bestAppendingLength,mutMatrix)
	newInternalNode.probVectUpLeft=mergeVectorsUpDown(vectUp,bestUpLength,node.probVect,bestDownLength,mutMatrix)
	if newInternalNode.probVect==None:
		print("Problem in placeSampleOnTree(), probVect is None")
		print(node.probVect)
		print(bestDownLength)
		print(newPartials)
		print(bestAppendingLength)
		exit()
	if newInternalNode.probVectUpRight==None:
		print("Problem in placeSampleOnTree(), probVectUpRight is None")
		print(vectUp)
		print(bestUpLength)
		print(newPartials)
		print(bestAppendingLength)
		exit()
	if newInternalNode.probVectUpLeft==None:
		print("Problem in placeSampleOnTree(), probVectUpLeft is None")
		print(vectUp)
		print(bestUpLength)
		print(node.probVect)
		print(bestDownLength)
		exit()
	if bestUpLength:
		newInternalNode.probVectTotUp=mergeVectorsUpDown(vectUp,bestUpLength/2,newInternalNode.probVect,bestUpLength/2,mutMatrix)
	if bestAppendingLength:
		newNode.probVectTotUp=mergeVectorsUpDown(newInternalNode.probVectUpLeft,bestAppendingLength/2,newPartials,bestAppendingLength/2,mutMatrix)
		updatePesudoCounts(newInternalNode.probVectUpLeft,newPartials,pseudoMutCounts)
	if not bestDownLength:
		node.probVectTotUp=None
	if verbose:
		print("new internal node added to tree")
		print(newInternalNode.probVect)
	nodeList=[(node,2),(newInternalNode.up,child)]
	updatePartials(nodeList,mutMatrix)

	return None











#set all descendant nodes to dirty.
#So far thses flags are used to prevent traversing the same part of the tree multiple times.
def setAllDirty(node):
	nextLeaves=[node]
	while nextLeaves:
		nextNode=nextLeaves.pop()
		nextNode.dirty=True
		#number of times this node has been replaced through SPR during this round.
		nextNode.replacements=0
		for c in nextNode.children:
			nextLeaves.append(c)







#function to calculate likelihood cost of appending node to parent node 
#differently from appendProb, this allows the bottom node to be internal, not just a sample.
def appendProbNode(probVectP,probVectC,bLen,mutMatrix,useRateVariation=False,mutMatrices=None):
	Lkcost, indexEntry1, indexEntry2, totalFactor, pos = 0.0, 0, 0, 1.0, 0
	entry1=probVectP[indexEntry1]
	entry2=probVectC[indexEntry2]
	end=min(entry1[1],entry2[1])
	contribLength=bLen
	while True:
		if entry2[0]==5: # case entry1 is N
			pos=min(entry1[1],entry2[1])
			pass
		elif entry1[0]==5: # case entry2 is N
			#if parent node is type "N", in theory we might have to calculate the contribution of root nucleotides; 
			# however, if this node is "N" then every other node in the current tree is "N", so we can ignore this since this contribution cancels out in relative terms.
			pos=min(entry1[1],entry2[1])
			pass
		else:
			#contribLength will be here the total length from the root or from the upper node, down to the down node.
			if entry1[0]<5:
				if len(entry1)==2:
					contribLength=bLen
				elif len(entry1)==3:
					contribLength=entry1[2]
					if bLen:
						contribLength+=bLen
				else:
					contribLength=entry1[3]
					if bLen:
						contribLength+=bLen
			else:
				if len(entry1)==3:
					contribLength=bLen
				else:
					contribLength=entry1[2]
					if bLen:
						contribLength+=bLen
			if entry2[0]<5:
				if len(entry2)==3:
					if contribLength:
						contribLength+=entry2[2]
					else:
						contribLength=entry2[2]
			else:
				if len(entry2)==4:
					if contribLength:
						contribLength+=entry2[2]
					else:
						contribLength=entry2[2]

			if entry1[0]==4: # case entry1 is R	
				if entry2[0]==4:
					if len(entry1)==4:
						end=min(entry1[1],entry2[1])
						contribLength+=entry1[2]
						Lkcost+=contribLength*(cumulativeRate[end]-cumulativeRate[pos])
						pos=end
					else:
						if contribLength:
							end=min(entry1[1],entry2[1])
							Lkcost+=contribLength*(cumulativeRate[end]-cumulativeRate[pos])
							pos=end
						else:
							pos=min(entry1[1],entry2[1])

				#entry1 is reference and entry2 is of type "O"
				elif entry2[0]==6:
					if useRateVariation:
						mutMatrix=mutMatrices[pos]
					i1=refIndeces[pos]
					if len(entry1)==4:
						tot=0.0
						for i in range4:
							if i1==i:
								tot2=rootFreqs[i]*(1.0+mutMatrix[i][i]*entry1[2])
							else:
								tot2=rootFreqs[i]*mutMatrix[i][i1]*entry1[2]
							if contribLength:
								tot3=0.0
								for j in range4:
									tot3+=mutMatrix[i][j]*entry2[-1][j]
								tot+=tot2*(entry2[-1][i] + contribLength*tot3)
							else:
								tot+=tot2*entry2[-1][i]
						tot/=rootFreqs[i1]
					else:
						if contribLength:
							tot=0.0
							for j in range4:
								tot+=mutMatrix[i1][j]*entry2[-1][j]
							tot*=contribLength
							tot+=entry2[-1][i1]
						else:
							tot=entry2[-1][i1]
					totalFactor*=tot
					pos+=1

				else: #entry1 is R and entry2 is a different but single nucleotide
					if useRateVariation:
						mutMatrix=mutMatrices[pos]
					if len(entry1)==4:
						i1=refIndeces[pos]
						i2=entry2[0]
						if contribLength:
							totalFactor*=((rootFreqs[i1]*mutMatrix[i1][i2]*contribLength*(1.0+mutMatrix[i1][i1]*entry1[2])+rootFreqs[i2]*mutMatrix[i2][i1]*entry1[2]*(1.0+mutMatrix[i2][i2]*contribLength))/rootFreqs[i1])
						else:
							totalFactor*=((rootFreqs[i2]*mutMatrix[i2][i1]*entry1[2])/rootFreqs[i1])
					else:
						if contribLength:
							totalFactor*=mutMatrix[refIndeces[pos]][entry2[0]]*contribLength
						else:
							return float("-inf")
					pos+=1

			# entry1 is of type "O"
			elif entry1[0]==6:
				if useRateVariation:
					mutMatrix=mutMatrices[pos]
				if entry2[0]==6:
					if contribLength:
						tot=0.0
						for j in range4:
							tot+=entry1[-1][j]*(entry2[-1][j] + contribLength*(mutMatrix[j][0]*entry2[-1][0]+mutMatrix[j][1]*entry2[-1][1]+mutMatrix[j][2]*entry2[-1][2]+mutMatrix[j][3]*entry2[-1][3]) )
						totalFactor*=tot
					else:
						tot=0.0
						for j in range4:
							tot+=entry1[-1][j]*entry2[-1][j]
						totalFactor*=tot
				else: #entry1 is "O" and entry2 is a nucleotide
					if entry2[0]==4:
						i2=refIndeces[pos]
					else:
						i2=entry2[0]
					if contribLength:
						totalFactor*=(entry1[-1][i2]+contribLength*(entry1[-1][0]*mutMatrix[0][i2]+entry1[-1][1]*mutMatrix[1][i2]+entry1[-1][2]*mutMatrix[2][i2]+entry1[-1][3]*mutMatrix[3][i2]))
					else:
						totalFactor*=entry1[-1][i2]
				pos+=1

			else: #entry1 is a non-ref nuc
				if useRateVariation:
					mutMatrix=mutMatrices[pos]
				if entry2[0]==entry1[0]:
					if len(entry1)==4:
						contribLength+=entry1[2]
					if contribLength:
						Lkcost+=mutMatrix[entry1[0]][entry1[0]]*contribLength
				else: #entry1 is a nucleotide and entry2 is not the same as entry1
					i1=entry1[0]
					if entry2[0]<5: #entry2 is a nucleotide
						if entry2[0]==4:
							i2=refIndeces[pos]
						else:
							i2=entry2[0]
						if len(entry1)==4:
							if contribLength:
								totalFactor*=(( rootFreqs[i1]*mutMatrix[i1][i2]*contribLength*(1.0+mutMatrix[i1][i1]*entry1[2]) + rootFreqs[i2]*mutMatrix[i2][i1]*entry1[2]*(1.0+mutMatrix[i2][i2]*contribLength) )/rootFreqs[i1])
							else:
								totalFactor*=((rootFreqs[i2]*mutMatrix[i2][i1]*entry1[2] )/rootFreqs[i1])
						else:
							if contribLength:
								totalFactor*=mutMatrix[i1][i2]*contribLength
							else:
								return float("-inf")
					else: #entry1 is a nucleotide and entry2 is of type "O"
						if len(entry1)==4:
							tot=0.0
							for i in range4:
								if i1==i:
									tot2=rootFreqs[i]*(1.0+mutMatrix[i][i]*entry1[2])
								else:
									tot2=rootFreqs[i]*mutMatrix[i][i1]*entry1[2]
								tot3=0.0
								for j in range4:
									tot3+=mutMatrix[i][j]*entry2[-1][j]
								tot+=tot2*(entry2[-1][i] + contribLength*tot3)
							totalFactor*=(tot/rootFreqs[i1])
						else:
							tot=0.0
							for j in range4:
								tot+=mutMatrix[i1][j]*entry2[-1][j]
							tot*=contribLength
							tot+=entry2[-1][i1]
							totalFactor*=tot
				pos+=1

		if totalFactor<=minimumCarryOver:
			try:
				if totalFactor<sys.float_info.min:
					return float("-inf")
			except:
				print("In appendProbNode() value error")
				print(totalFactor)
				return float("-inf")
			Lkcost+=log(totalFactor)
			totalFactor=1.0
		if pos==lRef:
			break
		if pos==entry1[1]:
			indexEntry1+=1
			entry1=probVectP[indexEntry1]
		if pos==entry2[1]:
			indexEntry2+=1
			entry2=probVectC[indexEntry2]
	return Lkcost+log(totalFactor)









#traverse the tree to optimize (only) the length of the branches using the derivative approach.
def traverseTreeToOptimizeBranchLengths(root,mutMatrix,testing=False,useRateVariation=False,mutMatrices=None,fastPass=False):
	totalLKimprovementBL=0.0
	updates=0
	nodesTraversed=0
	if root.children:
		nodesToTraverse=[root.children[0],root.children[1]]
	else:
		return 0
	while nodesToTraverse:
		nodesTraversed+=1
		node=nodesToTraverse.pop()
		if node==node.up.children[0]:
			upVect=node.up.probVectUpRight
			child=0
		else:
			upVect=node.up.probVectUpLeft
			child=1
		if node.dirty:
			bestLength=estimateBranchLengthWithDerivative(upVect,node.probVect,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices)
			if bestLength or node.dist:
				if testing:
					currentCost=appendProbNode(upVect,node.probVect,node.dist,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices)
					newCost=appendProbNode(upVect,node.probVect,bestLength,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices)
					totalLKimprovementBL+=newCost-currentCost
				if (not bestLength) or (not node.dist) or node.dist/bestLength>1.01 or node.dist/bestLength<0.99:
					node.dist=bestLength
					updates+=1
					if not fastPass:
						nodeList=[(node,2),(node.up,child)]
						updatePartials(nodeList,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices)
				else:
					node.dirty=False
			else:
				node.dirty=False
		for child in node.children:
			nodesToTraverse.append(child)
	if testing:
		return totalLKimprovementBL
	else:
		return updates
		










#we know that subtree "appendedNode", with partials "newPartials", is best placed as child of "node" resulting in logLK contribution of newChildLK
# now explore exactly which position near node is best for placement (direct child, sibling or sibling of existing child), and which are the best branch lengths,
#then add the subtree at that position of the tree, and update all the internal probability vectors.
def placeSubtreeOnTree(node,newPartials,appendedNode,newChildLK,bestBranchLengths,mutMatrix,useRateVariation=False,mutMatrices=None,parentNodeReplacements=1):
	bestAppendingLength=bestBranchLengths[2]
	bestUpLength=bestBranchLengths[0]
	bestDownLength=bestBranchLengths[1]
	tryNewRoot=False
	if node.up.children[0]==node:
		child=0
		vectUp=node.up.probVectUpRight
	else:
		child=1
		vectUp=node.up.probVectUpLeft
	
	#test if we should also attempt placing the subtree as child of a new root
	if not bestUpLength:
		pNode=node.up
		while (not pNode.dist) and (pNode.up!=None):
			pNode=pNode.up
		if pNode.up==None:
			root=pNode
			tryNewRoot=True
			#TODO Myrthe, this was the bit missing that was causing the incongruence in score
			if (not bestDownLength) or (bestDownLength>1.01*node.dist) or (bestDownLength<0.99*node.dist):
				node.dist=bestDownLength
				nodeList=[(node,2),(node.up,child)]
				updatePartials(nodeList,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices)
			
	#in case of best placement as a descendant appended exactly at the root node, attempt also to create new root
	if tryNewRoot:
		node=root
		probOldRoot = findProbRoot(node.probVect)
		rootUpLeft=rootVector(node.probVect,bestAppendingLength/2,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices)
		bestRightLength=estimateBranchLengthWithDerivative(rootUpLeft,newPartials,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices)
		rootUpRight=rootVector(newPartials,bestRightLength,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices)
		bestLeftLength=estimateBranchLengthWithDerivative(rootUpRight,node.probVect,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices)
		secondBranchLengthOptimizationRound=True
		if secondBranchLengthOptimizationRound: #if wanted, do a second round of branch length optimization
			rootUpLeft=rootVector(node.probVect,bestLeftLength,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices)
			bestRightLength=estimateBranchLengthWithDerivative(rootUpLeft,newPartials,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices)
			rootUpRight=rootVector(newPartials,bestRightLength,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices)
			bestLeftLength=estimateBranchLengthWithDerivative(rootUpRight,node.probVect,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices)
		probVectRoot,probRoot = mergeVectors(node.probVect,bestLeftLength,newPartials,bestRightLength,mutMatrix,returnLK=True,useRateVariation=useRateVariation,mutMatrices=mutMatrices)
		probRoot+= findProbRoot(probVectRoot)
		parentLKdiff=probRoot-probOldRoot
		if parentLKdiff<=newChildLK: #best is just placing as descendant of the root
			bestRightLength=bestAppendingLength
			bestLeftLength=False
			probVectRoot=mergeVectors(node.probVect,bestLeftLength,newPartials,bestRightLength,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices)
			rootUpRight=rootVector(newPartials,bestRightLength,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices)
		#now add new root to the tree
		newRoot=Tree()
		newRoot.dirty=True
		newRoot.replacements=parentNodeReplacements+1
		if probVectRoot==None:
			print("Issue with new root probVect inside placeSubtreeOnTree()")
			exit()
		newRoot.probVect=probVectRoot
		newRoot.probVectUpRight=rootUpRight
		newRoot.probVectUpLeft=rootVector(node.probVect,bestLeftLength,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices)
		node.up=newRoot
		node.dist=bestLeftLength
		newRoot.add_child(node)
		appendedNode.up=newRoot
		newRoot.add_child(appendedNode)
		appendedNode.dist=bestRightLength
		appendedNode.replacements+=1
		if verbose:
			print("new root added to tree")
			print(newRoot.probVect)
			print(newRoot.children[0].probVect)
			print(appendedNode.probVect)
		nodeList=[(node,2),(appendedNode,2)]
		updatePartials(nodeList,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices)
		return newRoot

	#print("adding internal node")
	#in all other cases (not attempting to add a new root) create a new internal node in the tree and add sample as a descendant.
	newInternalNode=Tree()
	newInternalNode.dirty=True
	newInternalNode.replacements=parentNodeReplacements+1
	node.up.children[child]=newInternalNode
	newInternalNode.up=node.up
	newInternalNode.add_child(node)
	node.up=newInternalNode
	appendedNode.up=newInternalNode
	appendedNode.replacements+=1
	newInternalNode.add_child(appendedNode)
	newInternalNode.probVect=mergeVectors(node.probVect,bestDownLength,newPartials,bestAppendingLength,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices)
	if newInternalNode.probVect==None:
		newInternalNode.probVectUpLeft=mergeVectorsUpDown(vectUp,bestUpLength,node.probVect,bestDownLength,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices)
		bestAppendingLength=estimateBranchLengthWithDerivative(newInternalNode.probVectUpLeft,newPartials,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices)
		newInternalNode.probVectUpRight=mergeVectorsUpDown(vectUp,bestUpLength,newPartials,bestAppendingLength,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices)
		bestDownLength=estimateBranchLengthWithDerivative(newInternalNode.probVectUpRight,node.probVect,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices)
		newInternalNode.probVect=mergeVectors(node.probVect,bestDownLength,newPartials,bestAppendingLength,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices)
		if newInternalNode.probVect==None:
			print("newInternalNode.probVect is None after updating the optimal branch lengths")
			bestAppendingLength=oneMutBLen/5
			bestDownLength=oneMutBLen/5
			newInternalNode.probVect=mergeVectors(node.probVect,bestDownLength,newPartials,bestAppendingLength,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices)
	newInternalNode.probVectUpRight=mergeVectorsUpDown(vectUp,bestUpLength,newPartials,bestAppendingLength,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices)
	if newInternalNode.probVectUpRight==None:
		bestUpLength=estimateBranchLengthWithDerivative(vectUp,newInternalNode.probVect,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices)
		newInternalNode.probVectUpLeft=mergeVectorsUpDown(vectUp,bestUpLength,node.probVect,bestDownLength,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices)
		bestAppendingLength=estimateBranchLengthWithDerivative(newInternalNode.probVectUpLeft,newPartials,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices)
		newInternalNode.probVectUpRight=mergeVectorsUpDown(vectUp,bestUpLength,newPartials,bestAppendingLength,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices)
		if newInternalNode.probVectUpRight==None:
			print("newInternalNode.probVectUpRight is None after updating the optimal branch lengths")
			bestUpLength=oneMutBLen/5
			bestAppendingLength=oneMutBLen/5
			newInternalNode.probVectUpRight=mergeVectorsUpDown(vectUp,bestUpLength,newPartials,bestAppendingLength,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices)
		newInternalNode.probVect=mergeVectors(node.probVect,bestDownLength,newPartials,bestAppendingLength,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices)
	newInternalNode.probVectUpLeft=mergeVectorsUpDown(vectUp,bestUpLength,node.probVect,bestDownLength,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices)
	if newInternalNode.probVectUpLeft==None:
		bestUpLength=estimateBranchLengthWithDerivative(vectUp,newInternalNode.probVect,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices)
		bestDownLength=estimateBranchLengthWithDerivative(newInternalNode.probVectUpRight,node.probVect,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices)
		newInternalNode.probVectUpLeft=mergeVectorsUpDown(vectUp,bestUpLength,node.probVect,bestDownLength,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices)
		if newInternalNode.probVectUpLeft==None:
			print("newInternalNode.probVectUpRight is None after updating the optimal branch lengths")
			bestUpLength=oneMutBLen/5
			bestDownLength=oneMutBLen/5
			newInternalNode.probVectUpLeft=mergeVectorsUpDown(vectUp,bestUpLength,node.probVect,bestDownLength,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices)
		newInternalNode.probVect=mergeVectors(node.probVect,bestDownLength,newPartials,bestAppendingLength,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices)
		newInternalNode.probVectUpRight=mergeVectorsUpDown(vectUp,bestUpLength,newPartials,bestAppendingLength,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices)
	appendedNode.dist=bestAppendingLength
	newInternalNode.dist=bestUpLength
	node.dist=bestDownLength
	if not bestAppendingLength:
		appendedNode.probVectTotUp=None
	if bestUpLength:
		newInternalNode.probVectTotUp=mergeVectorsUpDown(vectUp,bestUpLength/2,newInternalNode.probVect,bestUpLength/2,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices)
	if not bestDownLength:
		node.probVectTotUp=None
	if verbose:
		print("new internal node added to tree")
		print(newInternalNode.probVect)
	nodeList=[(node,2),(newInternalNode.up,child),(appendedNode,2)]
	updatePartials(nodeList,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices)

	return None





#remove node from the current position in the tree and re-attach it at a new given place new bestNode.
# First remove node from the tree, then update the genome lists;
# then find the exact best reattachment of node and update the genome lists again using function placeSubtreeOnTree().
def cutAndPasteNode(node,bestNode,bestBranchLengths,bestLK,mutMatrix,useRateVariation=False,mutMatrices=None):
	#remove node from the tree
	if verbose or debugging:
		print("In cutAndPasteNode() removing subtree from the tree, subtree root partials: ")
	parentNode=node.up
	parentNodeReplacements=parentNode.replacements
	if node==parentNode.children[0]:
		sibling=parentNode.children[1]
	else:
		sibling=parentNode.children[0]
	if parentNode.up!=None:
		if parentNode==parentNode.up.children[0]:
			childP=0
		else:
			childP=1
		parentNode.up.children[childP]=sibling
	sibling.up=parentNode.up
	if sibling.dist:
		if parentNode.dist:
			sibling.dist+=parentNode.dist
	else:
		sibling.dist=parentNode.dist
	#update likelihood lists after node removal
	if sibling.up==None:
		sibling.dist=1.0
		if verbose:
			print("cutAndPasteNode(), sibling node is root")
		if sibling.children:
			sibling.probVectUpRight=rootVector(sibling.children[1].probVect,sibling.children[1].dist,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices)
			sibling.probVectUpLeft=rootVector(sibling.children[0].probVect,sibling.children[0].dist,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices)
			nodeList=[(sibling.children[0],2),(sibling.children[1],2)]
			updatePartials(nodeList,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices)
	else:
		nodeList=[(sibling,2),(sibling.up,childP)]
		updatePartials(nodeList,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices)
		if verbose:
			if sibling.dist:
				print("cutAndPasteNode(), sibling node is not root. New overall likelihoods at previously sibling node: ")
				print(sibling.probVectTot)
	#re-place the node and re-update the vector lists
	newRoot = placeSubtreeOnTree(bestNode,node.probVect,node,bestLK,bestBranchLengths, mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices,parentNodeReplacements=parentNodeReplacements)
	
	topologyChanges[0]+=1
	if (writeTreesToFileEveryTheseSteps>0) and (topologyChanges[0]%writeTreesToFileEveryTheseSteps)==0:
		currentRoot=sibling
		while currentRoot.up!=None:
			currentRoot=currentRoot.up
		intermediateTreesFile.write("Topology "+str(topologyChanges[0])+"\n")
		if binaryTree:
			intermediateTreesFile.write(createBinaryNewick(currentRoot)+"\n")
		else:
			intermediateTreesFile.write(createNewick(currentRoot)+"\n")
		
	if (writeLKsToFileEveryTheseSteps>0) and (topologyChanges[0]%writeLKsToFileEveryTheseSteps)==0:
		currentRoot=sibling
		while currentRoot.up!=None:
			currentRoot=currentRoot.up
		totalLK=calculateTreeLikelihood(currentRoot,mutMatrix,useRateVariation=rateVariation,mutMatrices=mutMatrices)
		intermediateLKsFile.write("Topology "+str(topologyChanges[0])+", LK: "+str(totalLK)+"\n")


	#if the root of the tree has changed, return the new root
	if sibling.up==None:
		if newRoot!=None:
			return newRoot
		return sibling
	return newRoot


# try to find a re-placement of a dirty node of the tree to improve the topology.
# Cut out the subtree at this node, and look for somewhere else in the tree where to attach it (an SPR move).
# To find the best location of the new re-attachment, we traverse the tree starting at the current attachment node, and for each node visited we evaluate the reattachment,
# as done by the findBestParentTopology() function.
# To avoid traversing the whole tree for each SPR move, we use stopping conditions similar to those used in the initial sample placement process.
# After we find the best SPR move for the given node, we execute it using the cutAndPasteNode() function.
def traverseTreeForTopologyUpdate(node,mutMatrix,strictTopologyStopRules=strictTopologyStopRules,allowedFailsTopology=allowedFailsTopology,thresholdLogLKtopology=thresholdLogLKtopology,thresholdTopologyPlacement=thresholdTopologyPlacement,useRateVariation=False,mutMatrices=None):
	#track if the root has changed, so that the new root node can be returned.
	newRoot=None
	# has the branch length been updated so far? If we change the branch length above the current node, even if we don't perform any SPR, we still need to update the genome lists in the tree.
	bLenChanged=False
	totalImprovement=0.0
	# we avoid the root node since it cannot be re-placed with SPR moves
	if node.up!=None:
		#evaluate current placement
		parentNode=node.up
		if parentNode.children[0]==node:
			child=0
			vectUp=parentNode.probVectUpRight
		else:
			child=1
			vectUp=parentNode.probVectUpLeft
		#score of current tree
		bestCurrenBLen=node.dist
		originalLK=appendProbNode(vectUp,node.probVect,bestCurrenBLen,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices)
		if debugging:
			print("Original LK "+str(originalLK)+" from length "+str(bestCurrenBLen)+" for lists ")
			print(vectUp)
			print(node.probVect)
		bestCurrentLK=originalLK
		if bestCurrentLK<thresholdTopologyPlacement:
			bestCurrenBLen=estimateBranchLengthWithDerivative(vectUp,node.probVect,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices)
			if bestCurrenBLen or node.dist:
				if (not bestCurrenBLen) or (not node.dist) or node.dist/bestCurrenBLen>1.01 or node.dist/bestCurrenBLen<0.99:
					bLenChanged=True
				bestCurrentLK=appendProbNode(vectUp,node.probVect,bestCurrenBLen,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices)
				if bestCurrentLK<originalLK:
					#print("Inside traverseTreeForTopologyUpdate(), BL optimization brought to a decrease in LK.")
					bestCurrenBLen=node.dist
					bestCurrentLK=originalLK
					bLenChanged=False
				if debugging:
					print("New LK after modifing bLen: "+str(bestCurrentLK)+" from length "+str(bestCurrenBLen))
					print(mutMatrix)
				if bestCurrentLK==float("-inf"):
					print("Found an infinite cost of bestCurrentLK "+str(bestCurrentLK)+" using appendProbNode()")
					print(vectUp)
					print(node.probVect)
					print(bestCurrenBLen)
					exit()

		topologyUpdated=False
		if bestCurrentLK<thresholdTopologyPlacement:
			#now find the best place on the tree where to re-attach the subtree rooted ar "node"
			#but to do that we need to consider new vector probabilities after removing the node that we want to replace
			# this is done using findBestParentTopology().
			if debugging:
				if bLenChanged:
					print("Better branch length found at node "+str(node)+", from "+str(node.dist)+" (LK "+str(originalLK)+") to "+str(bestCurrenBLen)+" (LK "+str(bestCurrentLK)+"). Children")
					print(node.children)

			bestNodeSoFar , bestLKdiff , bestBranchLengths = findBestParentTopology(parentNode,child,bestCurrentLK,bestCurrenBLen,mutMatrix,strictTopologyStopRules=strictTopologyStopRules,allowedFailsTopology=allowedFailsTopology,thresholdLogLKtopology=thresholdLogLKtopology,useRateVariation=useRateVariation,mutMatrices=mutMatrices)
			if bestLKdiff==float("inf"):
				print("Found an infinite improvement of "+str(bestLKdiff)+" using findBestParentTopology()")
				print(bestNodeSoFar.probVectTotUp)
				print(node.probVect)
				print(bestBranchLengths)
				exit()
			if bestLKdiff<-1e50:
				print("Error: found likelihood cost is very heavy, this might mean that the reference used is not the same used to generate the input diff file")
				exit()
			if debugging:
				print("Best SPR move found for this node has cost "+str(bestLKdiff)+" affecting node with children")
				print(bestNodeSoFar.children)
				print("With appending distance "+str(bestCurrenBLen))
			if bestLKdiff+thresholdTopologyPlacement>bestCurrentLK:
				topologyUpdated=True
				topNode=node.up
				if bestNodeSoFar==topNode:
					topologyUpdated=False
				while (not topNode.dist) and (topNode.up!=None):
					topNode=topNode.up
				if bestNodeSoFar==topNode and (not bestBranchLengths[1]):
					topologyUpdated=False
				parentNode=node.up
				if node==parentNode.children[0]:
					sibling=parentNode.children[1]
				else:
					sibling=parentNode.children[0]
				if bestNodeSoFar==sibling:
					topologyUpdated=False
				#TODO Myrthe, modified here and added this
				if bestNodeSoFar.up==sibling and (not bestBranchLengths[0]):
					topologyUpdated=False

				if topologyUpdated:
					totalImprovement=(bestLKdiff-originalLK)
					if originalLK==float("-inf"):
						totalImprovement=(bestLKdiff-bestCurrentLK)
					if totalImprovement==float("inf"):
						print("Found an infinite topology improvement of "+str(totalImprovement)+" from originalLK "+str(originalLK)+", bestCurrentLK "+str(bestCurrentLK)+" and bestLKdiff "+str(bestLKdiff))
						exit()
					if verbose:
						print("\n\n In traverseTreeForTopologyUpdate() found SPR move with improvement "+str(totalImprovement))
					if debugging:
						print("Performing SPR move, detaching node "+str(node)+" with children")
						print(node.children)
						print("and reattaching it around node "+str(bestNodeSoFar)+" with children")
						print(bestNodeSoFar.children)
						print("And branch lengths")
						print(bestBranchLengths)
						print("for supposed improvement "+str(totalImprovement))
					newRoot = cutAndPasteNode(node,bestNodeSoFar,bestBranchLengths,bestLKdiff,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices)
					bLenChanged=False
		if (not topologyUpdated) and bLenChanged:
			if debugging:
				print("Changing branch length (pos 3) from "+str(node.dist)+" to "+str(bestCurrenBLen)+" at node "+str(node)+" with children ")
				print(node.children)
			node.dist=bestCurrenBLen
			nodeList=[(node,2),(node.up,child)]
			updatePartials(nodeList,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices)
			totalImprovement=(bestCurrentLK-originalLK)
			if originalLK==float("-inf"):
				totalImprovement=0
			if totalImprovement==float("inf"):
				print("Found an infinite branch length improvement of "+str(totalImprovement)+" from originalLK "+str(originalLK)+" and bestCurrentLK "+str(bestCurrentLK))
				exit()

	return newRoot,totalImprovement





#traverse the tree (here the input "node" will usually be the root), and for each dirty node ancountered, call traverseTreeForTopologyUpdate() 
# to attempt an SPR move by cutting the subtree rooted at this dirty node and trying to re-append it elsewhere.
def startTopologyUpdates(node,mutMatrix,checkEachSPR=False,strictTopologyStopRules=strictTopologyStopRules,allowedFailsTopology=allowedFailsTopology,thresholdLogLKtopology=thresholdLogLKtopology,thresholdTopologyPlacement=thresholdTopologyPlacement,useRateVariation=False,mutMatrices=None):
	nodesToVisit=[node]
	totalImprovement=0.0
	newRoot=None
	numNodes=0
	while nodesToVisit:
		newNode=nodesToVisit.pop()
		for c in newNode.children:
			nodesToVisit.append(c)
		if newNode.dirty and newNode.replacements<=maxReplacements:
			newNode.dirty=False
			if checkEachSPR:
				root=newNode
				while root.up!=None:
					root=root.up
				oldTreeLK=calculateTreeLikelihood(root,mutMatrix,checkCorrectness=True,useRateVariation=useRateVariation,mutMatrices=mutMatrices)
				reCalculateAllGenomeLists(root,mutMatrix, checkExistingAreCorrect=True,useRateVariation=useRateVariation,mutMatrices=mutMatrices)
			newRoot2,improvement=traverseTreeForTopologyUpdate(newNode,mutMatrix,strictTopologyStopRules=strictTopologyStopRules,allowedFailsTopology=allowedFailsTopology,thresholdLogLKtopology=thresholdLogLKtopology,thresholdTopologyPlacement=thresholdTopologyPlacement,useRateVariation=useRateVariation,mutMatrices=mutMatrices)
			if checkEachSPR:
				#print(" apparent improvement "+str(improvement))
				root=newNode
				while root.up!=None:
					root=root.up
				#print("Post-SPR tree: "+createBinaryNewick(root))
				newTreeLK=calculateTreeLikelihood(root,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices)
				reCalculateAllGenomeLists(root,mutMatrix, checkExistingAreCorrect=True,useRateVariation=useRateVariation,mutMatrices=mutMatrices)
				#print("Post-SPR tree likelihood: "+str(newTreeLK))
				print("In startTopologyUpdates, LK score of improvement "+str(newTreeLK)+" - "+str(oldTreeLK)+" = "+str(newTreeLK-oldTreeLK)+", was supposed to be "+str(improvement))
				if newTreeLK-oldTreeLK < improvement-1.0:
					print("In startTopologyUpdates, LK score of improvement "+str(newTreeLK)+" - "+str(oldTreeLK)+" = "+str(newTreeLK-oldTreeLK)+" is less than what is supposd to be "+str(improvement))
					exit()
			totalImprovement+=improvement
			if newRoot2!=None:
				newRoot=newRoot2
			numNodes+=1
			if (numNodes%1000)==0:
				print("Processed topology for "+str(numNodes)+" nodes.")
	return newRoot,totalImprovement




#create newick string of a given tree (input node is assumed to be the root).
#This function can create multifurcations, which are not welcome by some software - use createBinaryNewick() if you want only binary trees
def createNewick(node):
	nextNode=node
	stringList=[]
	direction=0
	while nextNode!=None:
		if nextNode.children:
			if direction==0:
				stringList.append("(")
				nextNode=nextNode.children[0]
			elif direction==1:
				stringList.append(",")
				nextNode=nextNode.children[1]
				direction=0
			else:
				if nextNode.dist:
					stringList.append("):"+str(nextNode.dist))
				else:
					stringList.append("):0.0")
				if nextNode.up!=None:
					if nextNode.up.children[0]==nextNode:
						direction=1
					else:
						direction=2
				nextNode=nextNode.up
		else:
			if len(nextNode.minorSequences)>0:
				stringList.append("("+nextNode.name+":0.0")
				for s2 in nextNode.minorSequences:
					stringList.append(","+s2+":0.0")
				if nextNode.dist:
					stringList.append("):"+str(nextNode.dist))
				else:
					stringList.append("):0.0")
			else:
				if nextNode.dist:
					stringList.append(nextNode.name+":"+str(nextNode.dist))
				else:
					stringList.append(nextNode.name+":0.0")
			if nextNode.up!=None:
				if nextNode.up.children[0]==nextNode:
					direction=1
				else:
					direction=2
			nextNode=nextNode.up
	stringList.append(";")
	return "".join(stringList)


def createBinaryNewick(node):
	nextNode=node
	stringList=[]
	direction=0
	numLeaves=0
	while nextNode!=None:
		if nextNode.children:
			if direction==0:
				stringList.append("(")
				nextNode=nextNode.children[0]
			elif direction==1:
				stringList.append(",")
				nextNode=nextNode.children[1]
				direction=0
			else:
				if nextNode.dist:
					stringList.append("):"+str(nextNode.dist))
				else:
					stringList.append("):"+str(0.0))
				if nextNode.up!=None:
					if nextNode.up.children[0]==nextNode:
						direction=1
					else:
						direction=2
				nextNode=nextNode.up
		else:
			numLeaves+=(1+len(nextNode.minorSequences))
			if len(nextNode.minorSequences)>0:
				for i in nextNode.minorSequences:
					stringList.append("(")
				stringList.append(nextNode.name+":")
				for s2 in nextNode.minorSequences:
					stringList.append(str(0.0)+","+s2+":"+str(0.0)+"):")
				if nextNode.dist:
					stringList.append(str(nextNode.dist))
				else:
					stringList.append(str(0.0))
			else:
				if nextNode.dist:
					stringList.append(nextNode.name+":"+str(nextNode.dist))
				else:
					stringList.append(nextNode.name+":"+str(0.0))
			if nextNode.up!=None:
				if nextNode.up.children[0]==nextNode:
					direction=1
				else:
					direction=2
			nextNode=nextNode.up
	stringList.append(";")
	print("created newick string for tree with "+str(numLeaves)+" leaves.")
	return "".join(stringList)

def countLeavesInTree(root):
	nodesToVisit=[root]
	numLeaves=0
	while nodesToVisit:
		node=nodesToVisit.pop()
		if node.children:
			for child in node.children:
				nodesToVisit.append(child)
		else:
			numLeaves+=1
			numLeaves+=len(node.minorSequences)
	return numLeaves





#Given a tree, and a final substitution rate matrix, calculate the likelihood of the tree
def calculateTreeLikelihood(root,mutMatrix,checkCorrectness=False,useRateVariation=False,mutMatrices=None):
	node=root
	#direction 0 means from parent, direction 1 means from a child
	lastNode=None
	direction=0
	totalLK=0.0
	while node!=None:
		if direction==0:
			if node.children:
				node=node.children[0]
			else:
				lastNode=node
				node=node.up
				direction=1
				#we can ignore the likelihood contribution from the normalization of the likelihoods of the ambiguity characters at terminal nodes:
				# these are anyway the same for all trees and substitution models!
		else :
			if lastNode==node.children[0]:
				node=node.children[1]
				direction=0
			else:
				newLower, Lkcontribution=mergeVectors(node.children[0].probVect,node.children[0].dist,node.children[1].probVect,node.children[1].dist,mutMatrix,returnLK=True,useRateVariation=useRateVariation,mutMatrices=mutMatrices)
				totalLK+=Lkcontribution
				if newLower==None:
					print("Strange, inconsistent lower genome list creation in calculateTreeLikelihood(); old list, and children lists")
					print(node.probVect)
					print(node.children[0].probVect)
					print(node.children[1].probVect)
					exit()
				elif numTopologyImprovements and checkCorrectness and areVectorsDifferentDebugging(node.probVect,newLower):
					print("Strange, while calculating tree likelihood encountered non-updated lower likelihood genome list at node "+str(node)+" with children ")
					print(node.children)
					print("new list")
					print(newLower)
					print("old list")
					print(node.probVect)
					print("child 0 list")
					print(node.children[0].probVect)
					print("child 1 list")
					print(node.children[1].probVect)
					exit()
				lastNode=node
				node=node.up
				direction=1
	#now add contribution from the root
	totalLK+=findProbRoot(root.probVect)
	return totalLK





#Given a tree and its genome lists, calculate mutation counts and waiting times for expectation maximization estimation of substitution rates
def expectationMaximizationCalculationRates(root,useRateVariation=False):
	node=root
	#direction 0 means from parent, direction 1 means from a child
	lastNode=None
	direction=0
	counts=[[0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0]]
	waitingTimes=[0.0,0.0,0.0,0.0]
	if useRateVariation:
		totTreeLength=0.0
		waitingTimesSites=[]
		countsSites=[]
		#list used to track missed waiting times across the genome due to Ns
		trackingNs=[]
		for i in range(lRef):
			waitingTimesSites.append([0.0,0.0,0.0,0.0])
			countsSites.append(0.0)
			trackingNs.append(0.0)
		trackingNs.append(0.0)

	while node!=None:
		if direction==0:
			if node.dist and node.up!=None:
				if useRateVariation:
					totTreeLength+=node.dist
				#update counts and waiting times
				if node==node.up.children[0]:
					probVectP=node.up.probVectUpRight
				else:
					probVectP=node.up.probVectUpLeft
				probVectC=node.probVect
				indexEntry1, indexEntry2, pos = 0, 0, 0
				entry1=probVectP[indexEntry1]
				entry2=probVectC[indexEntry2]
				end=min(entry1[1],entry2[1])
				contribLength=node.dist

				while True:
					if entry2[0]==5: # case N
						if useRateVariation:
							trackingNs[pos]-=node.dist
						pos=min(entry1[1],entry2[1])
						if useRateVariation:
							trackingNs[pos]+=node.dist
					elif entry1[0]==5: # case N
						if useRateVariation:
							trackingNs[pos]-=node.dist
						pos=min(entry1[1],entry2[1])
						if useRateVariation:
							trackingNs[pos]+=node.dist
					else:
						if entry1[0]<5:
							if len(entry1)==2:
								totLen1=node.dist
							elif len(entry1)==3:
								totLen1=entry1[2]
								if node.dist:
									totLen1+=node.dist
							else:
								totLen1=entry1[3]
								if node.dist:
									totLen1+=node.dist
								#consider contribution across the root twice, each time only considering the bit on the own side of the root,
						else:
							if len(entry1)==3:
								totLen1=node.dist
							else:
								totLen1=entry1[2]
								if node.dist:
									totLen1+=node.dist
						
						if entry2[0]<5:
							if len(entry2)==2:
								totLen2=False
							else:
								totLen2=entry2[2]
						else:
							if len(entry2)==3:
								totLen2=False
							else:
								totLen2=entry2[2]
						if entry1[0]==4: # case entry1 is R
							if entry2[0]==4:
								end=min(entry1[1],entry2[1])
								if not totLen2:
									for i in range4:
										waitingTimes[i]+=totLen1*(cumulativeBases[end][i]-cumulativeBases[pos][i])
								pos=end
							
							elif entry2[0]==6:
								if not totLen2:
									i1=refIndeces[pos]
									normalization=0.0
									if len(entry1)==4:
										contribLength=node.dist+entry1[3]
										if useRateVariation:
											waitingTimesSites[pos][i1]-=contribLength
										for i in range4:
											if i1==i:
												prob=rootFreqs[i]*(1.0+nonMutRates[i]*entry1[2])
												tot3=0.0
												for j in range4:
													tot3+=mutMatrix[i][j]*entry2[-1][j]
												tot3*=contribLength
												tot3+=entry2[-1][i]
												normalization+=prob*tot3
											else:
												prob=rootFreqs[i]*mutMatrix[i][i1]*entry1[2]*(1.0+nonMutRates[i]*contribLength)*entry2[-1][i]
												normalization+=prob
										for i in range4:
											if i1==i:
												prob=rootFreqs[i]*(1.0+nonMutRates[i]*entry1[2])
												for j in range4:
													if j==i:
														tot3=prob*(1.0+nonMutRates[i]*contribLength)*entry2[-1][j]
														waitingTimes[i]+=contribLength*tot3/normalization
														if useRateVariation:
															waitingTimesSites[pos][i]+=contribLength*tot3/normalization
													else:
														tot3=prob*mutMatrix[i][j]*contribLength*entry2[-1][j]/normalization
														waitingTimes[i]+=(contribLength/2)*tot3
														waitingTimes[j]+=(contribLength/2)*tot3
														counts[i][j]+=tot3
														if useRateVariation:
															waitingTimesSites[pos][i]+=(contribLength/2)*tot3
															waitingTimesSites[pos][j]+=(contribLength/2)*tot3
															countsSites[pos]+=tot3
											else:
												prob=rootFreqs[i]*mutMatrix[i][i1]*entry1[2]*(1.0+nonMutRates[i]*contribLength)*entry2[-1][i]
												waitingTimes[i]+=contribLength*prob/normalization
												if useRateVariation:
													waitingTimesSites[pos][i]+=contribLength*prob/normalization

									else:
										if useRateVariation:
											waitingTimesSites[pos][i1]-=totLen1
										for i in range4:
											if i1==i:
												normalization+=(1.0+nonMutRates[i]*totLen1)*entry2[-1][i]
											else:
												normalization+=mutMatrix[i1][i]*totLen1*entry2[-1][i]
										for i in range4:
											if i1==i:
												prob=(1.0+nonMutRates[i]*totLen1)*entry2[-1][i]
												waitingTimes[i]+=totLen1*prob/normalization
												if useRateVariation:
													waitingTimesSites[pos][i]+=totLen1*prob/normalization
											else:
												prob=mutMatrix[i1][i]*totLen1*entry2[-1][i]/normalization
												waitingTimes[i1]+=(totLen1/2)*prob
												waitingTimes[i]+=(totLen1/2)*prob
												counts[i1][i]+=prob
												if useRateVariation:
													waitingTimesSites[pos][i1]+=(totLen1/2)*prob
													waitingTimesSites[pos][i]+=(totLen1/2)*prob
													countsSites[pos]+=prob
								pos+=1
							
							else: #entry1 is reference and entry2 is a different but single nucleotide
								if not totLen2:
									i1=refIndeces[pos]
									i2=entry2[0]
									if len(entry1)<4:
										if useRateVariation:
											waitingTimesSites[pos][i1]-=totLen1/2
											waitingTimesSites[pos][i2]+=totLen1/2
											countsSites[pos]+=1
										waitingTimes[i1]+=(totLen1/2)
										waitingTimes[i2]+=(totLen1/2)
										counts[i1][i2]+=1
									else:
										contribLength=node.dist+entry1[3]
										prob1=rootFreqs[i1]*mutMatrix[i1][i2]*contribLength*(1.0+nonMutRates[i1]*entry1[2])
										prob2=rootFreqs[i2]*mutMatrix[i2][i1]*entry1[2]*(1.0+nonMutRates[i2]*contribLength)
										normalization=prob1+prob2
										prob1=prob1/normalization
										prob2=prob2/normalization
										waitingTimes[i1]+=(contribLength/2)*prob1
										waitingTimes[i2]+=(contribLength/2)*prob1
										counts[i1][i2]+=prob1
										waitingTimes[i2]+=contribLength*prob2
										if useRateVariation:
											waitingTimesSites[pos][i1]-=contribLength
											waitingTimesSites[pos][i1]+=(contribLength/2)*prob1
											waitingTimesSites[pos][i2]+=(contribLength/2)*prob1
											waitingTimesSites[pos][i2]+=contribLength*prob2
											countsSites[pos]+=prob1
								pos+=1

						# entry1 is of type "O"
						elif entry1[0]==6:
							if not totLen2:
								normalization=0.0
								if useRateVariation:
									waitingTimesSites[pos][refIndeces[pos]]-=totLen1
								if entry2[0]==6:
									for i in range4:
										for j in range4:
											if i==j:
												normalization+=entry1[-1][i]*(1.0+nonMutRates[i]*totLen1)*entry2[-1][j]
											else:
												normalization+=entry1[-1][i]*mutMatrix[i][j]*totLen1*entry2[-1][j]
									for i in range4:
										for j in range4:
											if i==j:
												prob=entry1[-1][i]*(1.0+nonMutRates[i]*totLen1)*entry2[-1][j]
												waitingTimes[i]+=totLen1*prob/normalization
												if useRateVariation:
													waitingTimesSites[pos][i]+=totLen1*prob/normalization
											else:
												prob=entry1[-1][i]*mutMatrix[i][j]*totLen1*entry2[-1][j]/normalization
												waitingTimes[i]+=(totLen1/2)*prob
												waitingTimes[j]+=(totLen1/2)*prob
												counts[i][j]+=prob
												if useRateVariation:
													waitingTimesSites[pos][i]+=(totLen1/2)*prob
													waitingTimesSites[pos][j]+=(totLen1/2)*prob
													countsSites[pos]+=prob
								else:
									if entry2[0]==4:
										i2=refIndeces[pos]
									else:
										i2=entry2[0]
									for i in range4:
										if i==i2:
											normalization+=entry1[-1][i]*(1.0+nonMutRates[i]*totLen1)
										else:
											normalization+=entry1[-1][i]*mutMatrix[i][i2]*totLen1
									for i in range4:
										if i==i2:
											prob=entry1[-1][i]*(1.0+nonMutRates[i]*totLen1)
											waitingTimes[i]+=totLen1*prob/normalization
											if useRateVariation:
												waitingTimesSites[pos][i]+=totLen1*prob/normalization
										else:
											prob=entry1[-1][i]*mutMatrix[i][i2]*totLen1/normalization
											waitingTimes[i]+=(totLen1/2)*prob
											waitingTimes[i2]+=(totLen1/2)*prob
											counts[i][i2]+=prob
											if useRateVariation:
												waitingTimesSites[pos][i]+=(totLen1/2)*prob
												waitingTimesSites[pos][i2]+=(totLen1/2)*prob
												countsSites[pos]+=prob
							pos+=1

						else: #entry1 is a non-ref nuc
							i1=entry1[0]
							if entry2[0]==i1:
								if not totLen2:
									waitingTimes[i1]+=totLen1
									if useRateVariation:
										waitingTimesSites[pos][i1]+=totLen1
										waitingTimesSites[pos][refIndeces[pos]]-=totLen1
							else: #entry1 and entry2 are of different types
								if entry2[0]==6:
									if not totLen2:
										normalization=0.0
										if len(entry1)==4:
											contribLength=node.dist+entry1[3]
											if useRateVariation:
												waitingTimesSites[pos][refIndeces[pos]]-=contribLength
											for i in range4:
												if i1==i:
													prob=rootFreqs[i]*(1.0+nonMutRates[i]*entry1[2])
													tot3=0.0
													for j in range4:
														tot3+=mutMatrix[i][j]*entry2[-1][j]
													tot3*=contribLength
													tot3+=entry2[-1][i]
													normalization+=prob*tot3
												else:
													prob=rootFreqs[i]*mutMatrix[i][i1]*entry1[2]*(1.0+nonMutRates[i]*contribLength)*entry2[-1][i]
													normalization+=prob
											for i in range4:
												if i1==i:
													prob=rootFreqs[i]*(1.0+nonMutRates[i]*entry1[2])
													for j in range4:
														if j==i:
															tot3=prob*(1.0+nonMutRates[i]*contribLength)*entry2[-1][j]
															waitingTimes[i]+=contribLength*tot3/normalization
															if useRateVariation:
																waitingTimesSites[pos][i]+=contribLength*tot3/normalization
														else:
															tot3=prob*mutMatrix[i][j]*contribLength*entry2[-1][j]/normalization
															waitingTimes[i]+=(contribLength/2)*tot3
															waitingTimes[j]+=(contribLength/2)*tot3
															counts[i][j]+=tot3
															if useRateVariation:
																waitingTimesSites[pos][i]+=(contribLength/2)*tot3
																waitingTimesSites[pos][j]+=(contribLength/2)*tot3
																countsSites[pos]+=tot3
												else:
													prob=rootFreqs[i]*mutMatrix[i][i1]*entry1[2]*(1.0+nonMutRates[i]*contribLength)*entry2[-1][i]
													waitingTimes[i]+=contribLength*prob/normalization
													if useRateVariation:
														waitingTimesSites[pos][i]+=contribLength*prob/normalization
										else:
											if useRateVariation:
												waitingTimesSites[pos][refIndeces[pos]]-=totLen1
											for i in range4:
												if i1==i:
													normalization+=(1.0+nonMutRates[i]*totLen1)*entry2[-1][i]
												else:
													normalization+=mutMatrix[i1][i]*totLen1*entry2[-1][i]
											for i in range4:
												if i1==i:
													prob=(1.0+nonMutRates[i]*totLen1)*entry2[-1][i]
													waitingTimes[i]+=totLen1*prob/normalization
													if useRateVariation:
														waitingTimesSites[pos][i]+=totLen1*prob/normalization
												else:
													prob=mutMatrix[i1][i]*totLen1*entry2[-1][i]/normalization
													waitingTimes[i1]+=(totLen1/2)*prob
													waitingTimes[i]+=(totLen1/2)*prob
													counts[i1][i]+=prob
													if useRateVariation:
														waitingTimesSites[pos][i1]+=(totLen1/2)*prob
														waitingTimesSites[pos][i]+=(totLen1/2)*prob
														countsSites[pos]+=prob
								#entry2 is a nucleotide type (like entry1)		
								else:
									if not totLen2:
										if entry2[0]==4:
											i2=refIndeces[pos]
										else:
											i2=entry2[0]
										if len(entry1)<4:
											if useRateVariation:
												waitingTimesSites[pos][refIndeces[pos]]-=totLen1
												waitingTimesSites[pos][i1]+=(totLen1/2)
												waitingTimesSites[pos][i2]+=(totLen1/2)
												countsSites[pos]+=1
											waitingTimes[i1]+=(totLen1/2)
											waitingTimes[i2]+=(totLen1/2)
											counts[i1][i2]+=1
										else:
											contribLength=node.dist+entry1[3]
											prob1=rootFreqs[i1]*mutMatrix[i1][i2]*contribLength*(1.0+nonMutRates[i1]*entry1[2])
											prob2=rootFreqs[i2]*mutMatrix[i2][i1]*entry1[2]*(1.0+nonMutRates[i2]*contribLength)
											normalization=prob1+prob2
											prob1=prob1/normalization
											prob2=prob2/normalization
											waitingTimes[i1]+=(contribLength/2)*prob1
											waitingTimes[i2]+=(contribLength/2)*prob1
											counts[i1][i2]+=prob1
											waitingTimes[i2]+=contribLength*prob2
											if useRateVariation:
												waitingTimesSites[pos][refIndeces[pos]]-=contribLength
												waitingTimesSites[pos][i1]+=(contribLength/2)*prob1
												waitingTimesSites[pos][i2]+=(contribLength/2)*prob1
												countsSites[pos]+=prob1
												waitingTimesSites[pos][i2]+=contribLength*prob2

							pos+=1

					if pos==lRef:
						break
					if pos==entry1[1]:
						indexEntry1+=1
						entry1=probVectP[indexEntry1]
					if pos==entry2[1]:
						indexEntry2+=1
						entry2=probVectC[indexEntry2]

			if node.children:
				node=node.children[0]
			else:
				lastNode=node
				node=node.up
				direction=1
			
		else :
			if lastNode==node.children[0]:
				node=node.children[1]
				direction=0
			else:
				lastNode=node
				node=node.up
				direction=1
			
	if model=="UNREST":
		for i in range4:
			if not waitingTimes[i]:
				for j in range4:
					counts[i][j]=0.0
			else:
				for j in range4:
					if i!=j:
						counts[i][j]/=waitingTimes[i]
				counts[i][i]=-sum(counts[i])
	elif model=="GTR":
		newRates=[[0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0]]
		for i in range4:
			if not waitingTimes[i]:
				for j in range4:
					newRates[i][j]=0.0
			else:
				for j in range4:
					if i!=j:
						newRates[i][j]=(counts[i][j]+counts[j][i])/waitingTimes[i]
				newRates[i][i]=-sum(newRates[i])
		counts=newRates
	else:
		print("Expectation Maximization for given model "+model+" not implemented yet")
		exit()
	totRate=-(rootFreqs[0]*counts[0][0]+ rootFreqs[1]*counts[1][1]+ rootFreqs[2]*counts[2][2]+ rootFreqs[3]*counts[3][3] )
	if totRate:
		for i in range4:
			for j in range4:
				counts[i][j]=counts[i][j]/totRate

	if useRateVariation:
		siteRates=[]
		totRate=0.0
		for i in range(lRef):
			waitingTimesSites[i][refIndeces[i]]+=totTreeLength+trackingNs[i]
			totExpected=0.0
			for j in range4:
				totExpected-=waitingTimesSites[i][j]*counts[j][j]
			if not totExpected:
				siteRates.append(1.0)
			else:
				siteRates.append(countsSites[i]/totExpected)
			totRate+=siteRates[-1]
		totRate=totRate/lRef
		for i in range(lRef):
			siteRates[i]=siteRates[i]/totRate
		return counts, siteRates
	else:
		return counts, None





if runFluTest:
	nodesToVisit=[tree1]
	while nodesToVisit:
		node=nodesToVisit.pop()
		if node.children:
			nodesToVisit.append(node.children[0])
			nodesToVisit.append(node.children[1])
		else:
			if node.name=="S104" or node.name=="S22":
				nodeToBeReplaced=node.up
				break
	debugging=True
	checkEachSPR=True

	if checkEachSPR:
		root=nodeToBeReplaced
		while root.up!=None:
			root=root.up
		oldTreeLK=calculateTreeLikelihood(root,mutMatrix,checkCorrectness=True)
		print("Pre-SPR tree likelihood: "+str(oldTreeLK))
	newRoot2,improvement=traverseTreeForTopologyUpdate(nodeToBeReplaced,mutMatrix)
	if checkEachSPR:
		print(" apparent improvement "+str(improvement))
		root=nodeToBeReplaced
		while root.up!=None:
			root=root.up
		newTreeLK=calculateTreeLikelihood(root,mutMatrix)
		print("Post-SPR tree likelihood: "+str(newTreeLK))
		reCalculateAllGenomeLists(root,mutMatrix, checkExistingAreCorrect=True)
		print("In startTopologyUpdates, LK score of improvement "+str(newTreeLK)+" - "+str(oldTreeLK)+" = "+str(newTreeLK-oldTreeLK)+", was supposed to be "+str(improvement))
		exit()


if runOnlyExample:
	data={}
	data["refSeq"]=[]
	data["Seq1"]=[('t',1),('-',10,1)]
	data["Seq1bis"]=[('t',1),('-',11,1)]
	data["Seq2"]=[('t',1),('a',3)]
	data["Seq3"]=[('t',1),('a',3),('t',5)]
	data["Seq4"]=[('a',2)]
	data["Seq5"]=[('a',2),('t',4)]
	data["Seq6"]=[('a',2),('t',6)]
	print(data)
	print(ref[:10])
	samples=data.keys()


if inputTree=="" and (writeTreesToFileEveryTheseSteps==0):
	print("Calculating distances with no sample names")
	distances=distancesFromRefPunishNs(data)
	print("Average divergence from refererence: "+str(totDivFromRef[0]/len(distances)))
else:
	print("Calculating distances with sample names")
	distances=distancesFromRefPunishNs(data,samples=data.keys())
print("Distances from the reference calculated")

if inputTree=="": #initialize tree to just the initial root sample
	#extract root genome among those closest to the reference but not empty
	firstSample=distances.pop()
	t1=Tree(name=firstSample[1])
	t1.probVect=probVectTerminalNode(data[firstSample[1]])
	data[firstSample[1]]=None
	t1.probVectTot=rootVector(t1.probVect,False,mutMatrix)
	t1.minorSequences=[]
else:
	t1=tree1

timeFinding=0.0
timePlacing=0.0
numSamples=0
while distances:
	d=distances.pop()
	numSamples+=1
	sample=d[1]
	newPartials=probVectTerminalNode(data[sample])
	data[sample]=None
	if (numSamples%updateSubstMatrixEveryThisSamples)==0:
		if model!="JC":
			if updateSubMatrix(pseudoMutCounts,model,mutMatrix):
				for i in range4:
					nonMutRates[i]=mutMatrix[i][i]
				for i in range(lRef):
					cumulativeRate[i+1]=cumulativeRate[i]+nonMutRates[refIndeces[i]]
	if (numSamples%1000)==0:
		print("Sample num "+str(numSamples))
	start=time()
	bestNode , bestScore, bestBranchLengths = findBestParentForNewSample(t1,newPartials,sample,mutMatrix)
	timeFinding+=(time()-start)
	#TODO Myrthe, this assumption that the cost is negative is not necessarily met anymore due to possible improvements to existing branch lengths in the tree.
	#if bestScore<0.5:
	if bestBranchLengths!=None:
		start=time()
		newRoot=placeSampleOnTree(bestNode,newPartials,sample,bestScore, bestBranchLengths[0], bestBranchLengths[1], bestBranchLengths[2],mutMatrix,pseudoMutCounts)
		if newRoot!=None:
			t1=newRoot
		timePlacing+=(time()-start)
	

if (numChildLKs[0]>0) and (not useFixedThresholdLogLKoptimizationTopology):
	aveChildLK=sumChildLKs[0]/numChildLKs[0]
	thresholdLogLKoptimizationTopology=max(thresholdLogLKoptimizationTopology,-0.2*aveChildLK)
	print("thresholdLogLKoptimizationTopology set to "+str(thresholdLogLKoptimizationTopology))
	



if runOnlyExample:
	print("Tree after initial placement:")
	newickString=createNewick(t1)
	print(newickString)
	print("Now making change to the tree to create imperfection. Moving")
	nodeToReplace=t1.children[1].children[0].children[0].children[1]
	destination=t1.children[0].children[1].children[1]
	newickString=createNewick(nodeToReplace)
	print(newickString)
	print("to")
	newickString=createNewick(destination)
	print(newickString)
	cutAndPasteNode(nodeToReplace,destination,False,0.0001,-100.0,mutMatrix)
	newickString=createNewick(t1)
	print(newickString)

data.clear()


#put sample names in the tree
if debugging and inputTree=="":
	nextLeaves=[t1]
	while nextLeaves:
		node=nextLeaves.pop()
		if not node.children:
			node.name="S"+str(node.name)
			print(node.name)
			print(node.probVect)
			for m in range(len(node.minorSequences)):
				node.minorSequences[m]="S"+str(node.minorSequences[m])
		else:
			for c in node.children:
				nextLeaves.append(c)


#recalculate all genome lists according to the final substitution model
if inputTree=="" or largeUpdate or rateVariation:
	start=time()
	reCalculateAllGenomeLists(t1,mutMatrix,countNodes=True)
	if model!="JC" or rateVariation:
		mutMatrix , siteRates = expectationMaximizationCalculationRates(t1,useRateVariation=rateVariation)
		if rateVariation:
			for i in range4:
				nonMutRates[i]=mutMatrix[i][i]
			for i in range(lRef):
				cumulativeRate[i+1]=cumulativeRate[i]+nonMutRates[refIndeces[i]]*siteRates[i]
			mutMatrices=[]
			for i in range(lRef):
				mutMatrices.append([])
				for j in range4:
					mutMatrices[i].append(list(mutMatrix[j]))
					for k in range4:
						mutMatrices[i][j][k]*=siteRates[i]
		else:
			for i in range4:
				nonMutRates[i]=mutMatrix[i][i]
			for i in range(lRef):
				cumulativeRate[i+1]=cumulativeRate[i]+nonMutRates[refIndeces[i]]
			mutMatrices=None
		reCalculateAllGenomeLists(t1,mutMatrix,countNodes=False,useRateVariation=rateVariation,mutMatrices=mutMatrices)
	else:
		mutMatrices=None
	timeRecalculation=time()-start
	print("Time to recalculate all genome lists: "+str(timeRecalculation))

	print("Number of nodes: "+str(numNodes[0]))
	print("Os per node: "+str(float(numNodes[4])/numNodes[0]))
	print("Nucs per node: "+str(float(numNodes[1])/numNodes[0]))
	print("Ns per node: "+str(float(numNodes[3])/numNodes[0]))
	#print("R per node: "+str(float(numNodes[2])/numNodes[0]))
	#print("Non-O per node: "+str(float(numNodes[1]+numNodes[2]+numNodes[3])/numNodes[0]))

	start=time()
	setAllDirty(t1)
	print("Post-model optimization branch length optimization, initial fast pass")
	improvement=traverseTreeToOptimizeBranchLengths(t1,mutMatrix,useRateVariation=rateVariation,mutMatrices=mutMatrices, fastPass=True)
	reCalculateAllGenomeLists(t1,mutMatrix,countNodes=False,useRateVariation=rateVariation,mutMatrices=mutMatrices)
	print("Now proper branch length optimization, first pass")
	setAllDirty(t1)
	improvement=traverseTreeToOptimizeBranchLengths(t1,mutMatrix,useRateVariation=rateVariation,mutMatrices=mutMatrices)
	print("Branch length optimization round 1, number of changes: "+str(improvement))
	subRound=0
	while subRound<20:
		if (not improvement):
			break
		subRound+=1
		improvement=traverseTreeToOptimizeBranchLengths(t1,mutMatrix,useRateVariation=rateVariation,mutMatrices=mutMatrices)
		print("branch length finalization subround "+str(subRound+1)+", number of changes "+str(improvement))
	timeForBranchOptimization=(time()-start)
	print("Time for updating branch lengths: "+str(timeForBranchOptimization))
else:
	mutMatrices=None

if (writeTreesToFileEveryTheseSteps>0):
	currentRoot=t1
	intermediateTreesFile.write("Topology 0\n")
	if binaryTree:
		intermediateTreesFile.write(createBinaryNewick(currentRoot)+"\n")
	else:
		intermediateTreesFile.write(createNewick(currentRoot)+"\n")
if (writeLKsToFileEveryTheseSteps>0):
	currentRoot=t1
	totalLK=calculateTreeLikelihood(currentRoot,mutMatrix,useRateVariation=rateVariation,mutMatrices=mutMatrices)
	intermediateLKsFile.write("Topology 0, LK: "+str(totalLK)+"\n")

#if asked, first run a round of short range topology search
timeTopology=0.0
if fastTopologyInitialSearch and (inputTree=="" or largeUpdate):
	#first update just branch lengths
	start=time()
	setAllDirty(t1)
	print("Preliminary branch length optimization")
	improvement=traverseTreeToOptimizeBranchLengths(t1,mutMatrix,useRateVariation=rateVariation,mutMatrices=mutMatrices)
	print("Branch length optimization round 1, number of changes: "+str(improvement))
	subRound=0
	while subRound<20:
		if (not improvement):
			break
		subRound+=1
		improvement=traverseTreeToOptimizeBranchLengths(t1,mutMatrix,useRateVariation=rateVariation,mutMatrices=mutMatrices)
		print("branch length finalization subround "+str(subRound+1)+" number of changes "+str(improvement))
	timeForBranchOptimization=(time()-start)
	print("Time for updating branch lengths: "+str(timeForBranchOptimization))

	#now update the topology with a shallow search
	start=time()
	setAllDirty(t1)
	newRoot,improvement=startTopologyUpdates(t1,mutMatrix,checkEachSPR=debugging,strictTopologyStopRules=strictTopologyStopRulesInitial,allowedFailsTopology=allowedFailsTopologyInitial,thresholdLogLKtopology=thresholdLogLKtopologyInitial,thresholdTopologyPlacement=thresholdTopologyPlacementInitial,useRateVariation=rateVariation,mutMatrices=mutMatrices)
	if newRoot!=None:
		t1=newRoot
	timeForUpdatingTopology=(time()-start)
	print("Time for initial traversal of the tree for only short range topology changes or branch length: "+str(timeForUpdatingTopology))
	timeTopology+=timeForUpdatingTopology
	print("LK improvement apparently brought: "+str(improvement))

	#run improvements only on the nodes that have been affected by some changes in the last round, and so on
	start=time()
	subRound=0
	while subRound<20:
		print("Topological subround "+str(subRound+1))
		newRoot,improvement=startTopologyUpdates(t1,mutMatrix,checkEachSPR=debugging,strictTopologyStopRules=strictTopologyStopRulesInitial,allowedFailsTopology=allowedFailsTopologyInitial,thresholdLogLKtopology=thresholdLogLKtopologyInitial,thresholdTopologyPlacement=thresholdTopologyPlacementInitial,useRateVariation=rateVariation,mutMatrices=mutMatrices)
		if newRoot!=None:
			t1=newRoot
		print("LK improvement apparently brought: "+str(improvement))
		if improvement<thresholdLogLKwholeTopologyImprovement:
			break
		subRound+=1
	timeForUpdatingTopology=(time()-start)
	print("Time for the subrounds of this traversal of the tree: "+str(timeForUpdatingTopology))
	timeTopology+=timeForUpdatingTopology

### testing
#print(str(countLeavesInTree(t1))+" leaves in tree after preliminary SPR.")

#now run deep topological improvements
debugging=False
for i in range(numTopologyImprovements):
	print("Starting topological impromevement attempt traversing number "+str(i+1))
	start=time()
	if inputTree=="" or largeUpdate:
		setAllDirty(t1)
	newRoot,improvement=startTopologyUpdates(t1,mutMatrix,checkEachSPR=debugging,useRateVariation=rateVariation,mutMatrices=mutMatrices)
	if newRoot!=None:
		t1=newRoot
	timeForUpdatingTopology=(time()-start)
	print("Time for this traversal of the tree to update the topology or branch length: "+str(timeForUpdatingTopology))
	timeTopology+=timeForUpdatingTopology
	print("LK improvement apparently brought: "+str(improvement))

	if runOnlyExample:
		print("Tree after improvement "+str(i))
		newickString=createNewick(t1)
		print(newickString)

	if improvement<thresholdLogLKwholeTopologyImprovement:
		print("Small improvement, stopping topological search.")
		break
	
	#run improvements only on the nodes that have been affected by some changes in the last round, and so on
	start=time()
	subRound=0
	while subRound<20:
		print("Topological subround "+str(subRound+1))
		newRoot,improvement=startTopologyUpdates(t1,mutMatrix,checkEachSPR=debugging,useRateVariation=rateVariation,mutMatrices=mutMatrices)
		if newRoot!=None:
			t1=newRoot
		print("LK improvement apparently brought: "+str(improvement))
		if improvement<thresholdLogLKwholeTopologyImprovement:
			break
		subRound+=1
	timeForUpdatingTopology=(time()-start)
	print("Time for the subrounds of this traversal of the tree: "+str(timeForUpdatingTopology))
	timeTopology+=timeForUpdatingTopology
	if not (inputTree=="" or largeUpdate):
		break



#optimize branch lengths of the final tree
if optimizeBranchLengths:
	start=time()
	setAllDirty(t1)
	print("Final branch length optimization")
	improvement=traverseTreeToOptimizeBranchLengths(t1,mutMatrix,useRateVariation=rateVariation,mutMatrices=mutMatrices)
	print("Branch length optimization round 1, number of changes: "+str(improvement))
	subRound=0
	while subRound<20:
		if (not improvement):
			break
		subRound+=1
		improvement=traverseTreeToOptimizeBranchLengths(t1,mutMatrix,useRateVariation=rateVariation,mutMatrices=mutMatrices)
		print("branch length finalization subround "+str(subRound+1)+", number of changes "+str(improvement))
	timeForBranchOptimization=(time()-start)
	print("Time for updating branch lengths: "+str(timeForBranchOptimization))


#calculate total likelihood
if calculateLKfinalTree:
	totalLK=calculateTreeLikelihood(t1,mutMatrix,useRateVariation=rateVariation,mutMatrices=mutMatrices)
	print("totalLK: "+str(totalLK))
	file=open(outputFile+"_LK.txt","w")
	file.write(str(totalLK)+"\n")
	file.close()

#Free space by deleting the genome lists in the tree.
nextLeaves=[t1]
while nextLeaves:
	node=nextLeaves.pop()
	node.probVect=None
	if node.dist:
		node.probVectTot=None
		if node.up!=None:
			node.probVectTotUp=None
			node.furtherMidNodes=None
	if node.children:
		node.probVectUpRight=None
		node.probVectUpLeft=None
		for c in node.children:
			nextLeaves.append(c)
print("Deleted genome lists.")

#put sample names in the tree
if (not debugging) and (writeTreesToFileEveryTheseSteps==0) and inputTree=="":
	#Read input file to collect all sample names
	fileI=open(inputFile)
	line=fileI.readline()
	if refFile=="":
		line=fileI.readline()
		while line[0]!=">":
			line=fileI.readline()
	sampleNames=[]
	while line!="" and line!="\n":
		name=line.replace(">","").replace("\n","")
		line=fileI.readline()
		while line!="" and line!="\n" and line[0]!=">":
			line=fileI.readline()
		sampleNames.append(name)
	fileI.close()
	print("Sample names read.")

	#put sample names in the tree
	nextLeaves=[t1]
	name=""
	while nextLeaves:
		node=nextLeaves.pop()
		if not node.children:
			name=sampleNames[node.name]
			sampleNames[node.name]=None
			node.name=name
			for m in range(len(node.minorSequences)):
				name=sampleNames[node.minorSequences[m]]
				sampleNames[node.minorSequences[m]]=None
				node.minorSequences[m]=name
		else:
			for c in node.children:
				nextLeaves.append(c)
	print("Sample names assigned.")


if binaryTree:
	newickString=createBinaryNewick(t1)
else:
	newickString=createNewick(t1)

if runOnlyExample:
	print("Final tree:")
	print(newickString)
	exit()

file=open(outputFile+"_tree.tree","w")
file.write(newickString)
file.close()
file=open(outputFile+"_subs.txt","w")
for i in range4:
	for j in range4:
		file.write(str(mutMatrix[i][j])+"\t")
	file.write("\n")
if rateVariation:
	file.write("\n\n"+"Site rates:")
	for i in range(lRef):
		file.write(str(i+1)+"\n"+str(siteRates[i])+"\n")
file.close()
print("Missed minor samples: "+str(totalMissedMinors[0]))
print("Final Substitution matrix:")
print(mutMatrix)
print("Time spent finding placement nodes: "+str(timeFinding))
print("Time spent placing samples on the tree: "+str(timePlacing))
print("Time spent in total updating the topology and branch lengths: "+str(timeTopology))
exit()





