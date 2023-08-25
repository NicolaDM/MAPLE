import sys
from math import log
import argparse
from time import time
import os.path

#©EMBL-European Bioinformatics Institute, 2021-2023
#Developd by Nicola De Maio, with contributions from Myrthe Willemsen.

# MAPLE code to estimate a tree by maximum likelihood from a MAPLE format input.

#things in progress
#TODO C++ implementation 
#TODO parallelization
#TODO branch support measures

#longer-term plans
#TODO MAT reference
#TODO codon models, selection
#TODO indels, alignment
#TODO Timed trees
#TODO Bayesian implementation in BEAST
#TODO Recombination
#TODO phylogeography
#TODO phylodynamics/selection
#TODO relativistic phylogenetics for arbitrary large trees

#things one could MAYBE also try in the future:
#TODO This would change the structure of a lot of the code - however one could try to represent internal genomes/likelihoods a mutation-annotated tree instead of representing all nucleotides different from the reference at each internal node.
#TODO branch lenth optimization could maybe be made more efficient by representing the monomer (in t) factors in the derivative all jointly as a polymonial, so to reduce numbers of divisions to be performed (at the cost of increasing numbers of products and sums).
#TODO use slots (https://www.geeksforgeeks.org/slots-in-python/) to reduce memory cost?
#TODO develop more nuanced exploration of topology space: reattempt to re-place only subtrees whose root have been affected in the last round of topology search (which is now already done), 
# and their neighbours (which is not implemented yet)?
#TODO for increased accuracy, one could calculate probVectTot genome lists in a way that would not depend on the choice of which two genome lists to merge first. 
#For now implemented consistent calculation of tot likelihoods so to reduce inconsinstencies met.
#TODO would it be convenient to remove the genome position element from entries of type ACGTO ? Sounds like a lot of changes, and might make things slower and less readable.
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
parser.add_argument("--thresholdLogLKwholeTopologyImprovement",help="logLK difference threshold to consider something as a significant decrease in log-LK when considering whole rounds of SPR moves.",  type=float, default=3.0)
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
parser.add_argument('--inputNexusTree',default="", help='input nexus tree file name; this is optional, and is used for lineage assignment. The nexus tree is supposed to be the output of MAPLE, so that it has alternativePlacements annotation to represent topological uncertainty.')

parser.add_argument("--defaultBLen",help="Default length of branches, for example when the input tree has no branch length information.",  type=float, default=0.000033)
parser.add_argument("--normalizeInputBLen",help="For the case the input tree has branch lengths expressed not in a likelihood fashion (that is, expected number of substitutions per site), then multiply the input branch lengths by this factor. Particularly useful when using parsimony-based input trees.",  type=float, default=1.0)
parser.add_argument("--multipleInputRFTrees", help="Use this option if the file specified by option --inputRFtrees contains multiple trees - otherwise only read the first tree from that file.", action="store_true")
parser.add_argument("--maxReplacements",help="Maximum number of replacements attempts per node per SPR round (prevents loops).",  type=int, default=10)
parser.add_argument("--useFixedThresholdLogLKoptimizationTopology", help="Use this option if you want to specify the value in --thresholdLogLKoptimizationTopology instead of estimating it from the data.", action="store_true")

parser.add_argument("--writeTreesToFileEveryTheseSteps", help="By default, don't write intermediate trees to file. If however a positive integer is specified with this option, intermediate trees will be written to file every this many topological changes.",  type=int, default=0)
parser.add_argument("--writeLKsToFileEveryTheseSteps", help="By default, don't write likelihoods of intermediate trees to file. If however a positive integer is specified with this option, likelihoods of intermediate trees will be written to file every this many topological changes.",  type=int, default=0)

parser.add_argument("--estimateErrorRate", help="Estimate a single error rate for the whole genome. Input value is used as starting value", action="store_true")
parser.add_argument("--estimateSiteSpecificErrorRate", help="Estimate a separate error rate for each genome genome. Input value is used as starting value", action="store_true")
parser.add_argument("--errorRateInitial", help="Initial value used for estimating the error rate. The default is the inverse of the reference genome length (one error expected per genome).", type=float, default=0.0)
parser.add_argument("--errorRateFixed", help="Fix the error rate to a given input value", type=float, default=0.0)
parser.add_argument("--errorRateSiteSpecificFile", help="provide a file directory to a file that contains the siteSpecific error rates", type=str, default=None)
parser.add_argument("--estimateErrors", help="Estimate erroneous positions in the input sequences. This option is only allowed if option --estimateSiteSpecificErrorRate or errorRateSiteSpecificFile is also used.", action="store_true")
parser.add_argument("--minErrorProb", help="Minimum error probability to be reported when using option --estimateErrors", type=float, default=0.01)

parser.add_argument("--SPRTA", help="Calculate branch support values with a modification of the aBayes approach of Anisimova et al 2011 (De Maio et al 2023, in prep.).", action="store_true")
parser.add_argument("--aBayesPlus", help="Synonim for option --SPRTA.", action="store_true")
parser.add_argument("--networkOutput", help="Include in the output tree the alternative branches with their support (like in a weighted phylogenetic network).", action="store_true")
parser.add_argument("--minBranchSupport", help="Minimum branch support to be considered when using option --networkOutput", type=float, default=0.01)
parser.add_argument("--supportFor0Branches", help="When calculating branch support, also consider supports of branches of length 0, such as samples less informative than other samples which might have multiple placements on the tree.", action="store_true")
parser.add_argument("--supportForIdenticalSequences", help="Also write support measures multiple times for identical sequences (useful when one wants supports written on every branch).", action="store_true")
parser.add_argument("--estimateMAT", help="Estimate mutation events on the final tree, and write them on the nexus tree.", action="store_true")
parser.add_argument("--minMutProb", help="Minimum mutation probability to be written to output when using option --estimateMAT", type=float, default=0.01)
parser.add_argument("--doNotImproveTopology", help="Do not perform SPR moves, despite searching for them; this is useful if one wants to analyse a tree and calculate branch supports without changing the input tree.", action="store_true")
parser.add_argument("--keepInputIQtreeSupports", help="Assumes the input tree is from IQTREE, reads the support values on the tree branches, and prints the same values in the output nexus tree.", action="store_true")
args = parser.parse_args()

onlyNambiguities=args.onlyNambiguities
thresholdProb=args.thresholdProb
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

assignmentFile=args.assignmentFile
assignmentFileCSV=args.assignmentFileCSV
inputNexusTree=args.inputNexusTree
maxReplacements=args.maxReplacements
writeTreesToFileEveryTheseSteps=args.writeTreesToFileEveryTheseSteps
writeLKsToFileEveryTheseSteps=args.writeLKsToFileEveryTheseSteps

#errorTesting = args.errorTesting
mutMatrixGlobal=None
estimateErrorRate=args.estimateErrorRate
estimateSiteSpecificErrorRate = args.estimateSiteSpecificErrorRate
errorRateInitial = args.errorRateInitial
errorRateFixed = args.errorRateFixed
errorRateSiteSpecificFile = args.errorRateSiteSpecificFile
estimateErrors=args.estimateErrors
if estimateErrors and (not estimateSiteSpecificErrorRate) and (not errorRateSiteSpecificFile):
	print("Because option --estimateErrors has been selected, I am switching on --estimateSiteSpecificErrorRate so to allow estimation of site-specific error probabilities.")
	estimateSiteSpecificErrorRate=True
minErrorProb=args.minErrorProb
errorRateGlobal=0.0
usingErrorRate=False
errorRateSiteSpecific=False
rateVariation=args.rateVariation
useRateVariation=False
mutMatrices=None

aBayesPlus=args.aBayesPlus
sprta=args.SPRTA
aBayesPlus=aBayesPlus or sprta
networkOutput=args.networkOutput
minBranchSupport=args.minBranchSupport
supportFor0Branches=args.supportFor0Branches
supportForIdenticalSequences=args.supportForIdenticalSequences
estimateMAT=args.estimateMAT
minMutProb=args.minMutProb
doNotImproveTopology=args.doNotImproveTopology
keepInputIQtreeSupports=args.keepInputIQtreeSupports
aBayesPlusOn=False
if aBayesPlus:
	import math

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


# now also reads IQTREE branch supports if keepInputIQtreeSupports==True
# different types of branch support in IQTREE:
# )63:
# )0.650000:
# )/0.125:
# )75.4:
# )/0.999:
# )75.4/67.3:

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
					if node.dist<0.0:
						print("Warning: negative branch length in the input tree: "+distStr+" ; converting it to positive.")
						node.dist=abs(node.dist)
					distStr=""
				else:
					node.dist=defaultBLen
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
					node.dist=float(distStr)*normalizeInputBLen
					distStr=""
				else:
					node.dist=defaultBLen
				if keepInputIQtreeSupports:
					suppStr=""
					index+=1
					while nwString[index]!=":" and nwString[index]!=")" and nwString[index]!=";":
						suppStr+=nwString[index]
						index+=1
					if suppStr!="":
						suppVal=float(suppStr.split("/")[-1])
						if suppVal>1:
							suppVal=suppVal/100.0
						node.up.IQsupport=suppVal
				else:
					index+=1
				node=node.up
			else:
				name+=nwString[index]
				index+=1
		if not finished:
			print("Error, final character ; not found in newick string in file "+nwFile+".")
			raise Exception("exit")

		if not multipleTrees:
			break
		line=phyloFile.readline()

	phyloFile.close()
	return trees

def assignNodeBLen(node,distStr):
	if distStr!="":
		node.dist=float(distStr)*normalizeInputBLen
		if node.dist<0.0:
			print("Warning: negative branch length in the input tree: "+distStr+" ; converting it to positive.")
			node.dist=abs(node.dist)
		distStr=""
	else:
		node.dist=defaultBLen

def assignNodeFeatures(node,annotationString):
	#print("Assigning fatures "+annotationString)
	st=annotationString.replace("[","").replace("]","")
	features={}
	index=0
	while index<len(st):
		oldIndex=index
		while st[index]!="=":
			index+=1
		featureName=st[oldIndex:index]
		featureName=featureName.replace("&","")
		index+=1
		if st[index]=="{":
			oldIndex=index
			while st[index]!="}":
				index+=1
			featureStr=st[oldIndex:index].replace("{","").replace("}","")
			featureList=featureStr.split(",")
			featureDict={}
			for featureElement in featureList:
				if featureElement!="":
					featureElementList=featureElement.split(":")
					if len(featureElementList)==2:
						featureDict[featureElementList[0]]=float(featureElementList[1])
					else:
						featureDict[featureElementList[0]]=None
			features[featureName]=featureDict
			index+=1
		else:
			oldIndex=index
			while st[index]!="}" and st[index]!=",":
				index+=1
			featureStr=st[oldIndex:index]
			try:
				featureValue=float(featureStr)
				features[featureName]=featureValue
			except ValueError:
				features[featureName]=featureStr
		if index<len(st) and st[index]==",":
			index+=1
	node.features=features
	#if (not node.children) and ("alternativePlacements" in features) and (features["alternativePlacements"]):
	#	print(node.name)
	#	print(features["alternativePlacements"])
	#	print()
	#print("Assigned fatures ")
	#print(features)
	#print()




#function to read input newick string
def readNexus(nxFile,dirtiness=True):
	print("Reading Nexus file "+nxFile)
	phyloFile=open(nxFile)
	line=phyloFile.readline()
	while line!="begin trees;\n":
		line=phyloFile.readline()
		if line=="":
			print("Error, no tree found in nexus file "+nxFile+".")
			raise Exception("exit")
	line=phyloFile.readline()
	nwString=line.replace("\n","").split()[4]
	index=0
	node=Tree()
	node.dirty=dirtiness
	name=""
	distStr=""
	annotationString=""
	finished=False
	madeUpNameCount=0
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
			if name!="":
				node.name=name
				name=""
			else:
				madeUpNameCount+=1
				node.name="madeUpNodeName"+str(madeUpNameCount)
			assignNodeBLen(node,distStr)
			distStr=""
			assignNodeFeatures(node,annotationString)
			annotationString=""
			tree=node
			finished=True
			break
		elif nwString[index]=="[":
			firstInd=index
			while nwString[index]!="]":
				index+=1
			lastInd=index
			annotationString=nwString[firstInd:lastInd+1]
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
			else:
				madeUpNameCount+=1
				node.name="madeUpNodeName"+str(madeUpNameCount)
			assignNodeBLen(node,distStr)
			distStr=""
			assignNodeFeatures(node,annotationString)
			annotationString=""
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
			else:
				madeUpNameCount+=1
				node.name="madeUpNodeName"+str(madeUpNameCount)
			assignNodeBLen(node,distStr)
			distStr=""
			assignNodeFeatures(node,annotationString)
			annotationString=""
			index+=1
			node=node.up
		else:
			name+=nwString[index]
			index+=1
	if not finished:
		print("Error, final character ; not found in newick string in file "+nxFile+".")
		raise Exception("exit")

	phyloFile.close()
	print("Finished reading Nexus file.")
	return tree

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
	"""
	Myrthe:
	I extended this to calculate the RFL distance (see the Robinson Foulds paper) where they describe this.
	Basically, for every branch that is shared among the two trees, add the absolute difference of the two branches to the RFL.
	(This is first stored seperately as the KL value first).
	For every branch that is only in either of the two trees, add the length to the RFL distance.
	For the tree that is 'prepared for comparison', usually the true tree, I will store two hash tables (python dictionaries) together containing all branch lengths.
	One contains all leaf-branches. These are indexed by the leaf's name, usually a number (node.name).
	The other dictionary contains all inner branches, indexed by the <Left, Right> values of the subtree that is attached by the branch of interest.
	(One could improve the implementation by writing one's own hash maps and functions. Given that the indices have certain limits to what values they can become, it would be easy to obtain a perfect hash function)
	Furthermore, I calculate the sum of all the inner branch lengths. Whenever in the second function of the Day's algorithm (the comparison with the estimated tree),
	a branch is found in both trees, I add the absolute difference between the two branches, and substract the branch length of the true tree.
	I don't do this for leaf branches, as we can be sure that a leaf branch exists in both trees (if higher than the minimum length).
	(I do this because it is less straight forward to add branch lengths of the true tree when they are not found during the comparison with the estimated tree)

	Similar to the normal RF distance function, I do not take into account branches with length smaller than the minimum branch length.
	However, cases may occur that a branch exists in both trees, but in one of them it is lower than the minimum branch length, making the RFL less accurate.

	:param t1: input tree
	:param rooted:
	:param minimumBLen:
	:param addRootRFL: Should the distance from the root node up be taken into account? I think this has usually no significance, and is always set to 1.0.
	:return: Various metrics, among which the RFL.
	"""
	#dictionary of values given to sequence names
	leafNameDict={}
	#list of sequence names sorted according to value
	leafNameDictReverse=[]
	#table containing clusters in the tree
	nodeTable=[]
	#dictionaries with branch lengths
	branchLengthDict, leafDistDict = {}, {}
	sumBranchLengths = 0 # sum of all branch lengths, except from the leaves
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
				leafDistDict[node.name] = node.dist
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
				if node!=t1:
					sumBranchLengths += node.dist
				if node==t1:
					nodeTable[lastR][0]=lastL
					nodeTable[lastR][1]=lastR
				else:
					if (not rooted) and (node.up==t1) and (len(t1.children)==2):
						if node==t1.children[1]:
							currentBL=node.dist+t1.children[0].dist
							addBranch=True
						else:
							addBranch=False
					else:
						currentBL=node.dist
						addBranch=True
					if addBranch and currentBL>minimumBLen:
						numBranches+=1
						if rooted or lastL>0:
							if node==node.up.children[-1]:
								if nodeTable[lastL][0]==0 and nodeTable[lastL][1]==0:
									nodeTable[lastL][0]=lastL
									nodeTable[lastL][1]=lastR
								else:
									nodeTable[lastR][0]=lastL
									nodeTable[lastR][1]=lastR
							else:
								nodeTable[lastR][0]=lastL
								nodeTable[lastR][1]=lastR
							branchLengthDict[(lastL, lastR)] = currentBL
						else: # re-root at leaf 0, so flip the values for the current branch if it contains leaf 0.
								flippedL=lastR+1
								flippedR=nLeaves-1
								nodeTable[flippedL][0]=flippedL
								nodeTable[flippedL][1]=flippedR
								branchLengthDict[(flippedL, flippedR)] = currentBL
			else:
				nextNode=node.children[node.exploredChildren]
				movingFrom=0
		node=nextNode
	return leafNameDict, nodeTable, leafCount, numBranches, leafDistDict, branchLengthDict, sumBranchLengths

#Robinson-Foulds distance (1981) using a simplification of the algorithm from Day 1985.
#this function compares the current tree t2 to a previous one for which prepareTreeComparison() was run.
# Example usage: leafNameDict, nodeTable, leafCount, numBranches = prepareTreeComparison(phyloTrue,rooted=False)
# numDiffs, normalisedRF, leafCount, foundBranches, missedBranches, notFoundBranches = RobinsonFouldsWithDay1985(phyloEstimated,leafNameDict,nodeTable,leafCount,numBranches,rooted=False)
def RobinsonFouldsWithDay1985(t2,leafNameDict,nodeTable,leafCount,numBranches,leafDistDict, branchLengthDict, sumBranchLengths,rooted=False,minimumBLen=0.000006):
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
	RFL = sumBranchLengths
	KF = 0
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
				trueDist = leafDistDict[node.name]
				KF += abs(trueDist- node.dist)   #As described, I have not added branch lengths of leaf nodes to sumBranchLengths, so no need for: RFL=- trueDist
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
					if (not rooted) and (node.up==t2) and (len(t2.children)==2):
						if node==t2.children[1]:
							currentBL=node.dist+t2.children[0].dist
							searchBranch=True
						else:
							searchBranch=False
					else:
						currentBL=node.dist
						searchBranch=True
					if searchBranch and currentBL>minimumBLen:
						if (lastR+1-lastL)==lastDesc:
							if rooted or lastL>0:						
								if nodeTable[lastL][0]==lastL and nodeTable[lastL][1]==lastR:
									foundBranches+=1
									trueDist = branchLengthDict[(lastL, lastR)]
									KF += abs(trueDist - currentBL)
									RFL-= trueDist
								elif nodeTable[lastR][0]==lastL and nodeTable[lastR][1]==lastR:
									foundBranches+=1
									trueDist = branchLengthDict[(lastL, lastR)]
									KF += abs(trueDist - currentBL)
									RFL-= trueDist
								else:
									missedBranches+=1
									RFL += currentBL
							else: # re-root at leaf 0, so flip the values for the current branch if it contains leaf 0.
								flippedL=lastR+1
								flippedR=leafCount-1
								if nodeTable[flippedL][0]==flippedL and nodeTable[flippedL][1]==flippedR:
									foundBranches+=1
									trueDist = branchLengthDict[(flippedL, flippedR)]
									KF += abs(trueDist - currentBL)
									RFL -= trueDist
								elif nodeTable[flippedR][0]==flippedL and nodeTable[flippedR][1]==flippedR:
									foundBranches+=1
									trueDist = branchLengthDict[(flippedL, flippedR)]
									KF += abs(trueDist - currentBL)
									RFL -= trueDist
								else:
									missedBranches+=1
									RFL += currentBL
						else:
							missedBranches+=1
							RFL += currentBL
			else:
				nextNode=node.children[node.exploredChildren]
				movingFrom=0
		node=nextNode
	if visitedLeaves<leafCount:
		print("There are leaves in the reference that have not been found in this new tree - leafCount "+str(leafCount)+" visitedLeaves "+str(visitedLeaves))
		return None, None, None, None, None, None
	#first value is number of differences, second value is max number of differences just in case one wants the normalized values; 
	#the other values are there just in case one wants more detail.
	numDiffs=((numBranches-foundBranches)+missedBranches)
	RFL += KF
	if rooted:
		normalization=numBranches+leafCount-2
	else:
		normalization=numBranches+leafCount-3
	return numDiffs, float(numDiffs)/(normalization), leafCount, foundBranches, missedBranches, (numBranches-foundBranches), RFL


performLineageAssignment=False
if assignmentFile!="" or assignmentFileCSV!="":
	performLineageAssignment=True

#generate the string corresponding to a node, taking into account support measures.
def stringForNode(nextNode,name,dist):
	#print("start stringForNode")
	#print(name)
	#print(dist)
	stringList=[]
	if networkOutput or (not nextNode.children):
		stringList.append(name)
	printIQtreeSupportForNode=False
	if keepInputIQtreeSupports:
		if hasattr(nextNode, 'IQsupport'):
			printIQtreeSupportForNode=True
	if aBayesPlusOn or estimateMAT or printIQtreeSupportForNode or printIQtreeSupportForNode:
		#print("first if")
		if nextNode.up!=None and (dist or supportFor0Branches or usingErrorRate):
			stringList.append("[&")
			if aBayesPlusOn and (dist or supportFor0Branches):
				stringList.append("support="+str(nextNode.support))
				if networkOutput:
					stringList.append(",alternativePlacements={")
					for iNode in range(len(nextNode.alternativePlacements)):
						stringList.append(nextNode.alternativePlacements[iNode][0].name+":"+str(nextNode.alternativePlacements[iNode][1]))
						if iNode<(len(nextNode.alternativePlacements)-1):
							stringList.append(",")
					stringList.append("}")
				if estimateMAT and (dist or usingErrorRate):
					stringList.append(",")
			if estimateMAT and (dist or usingErrorRate):
				stringList.append("mutations={")
				for iNode in range(len(nextNode.mutations)):
					mutation=nextNode.mutations[iNode]
					stringList.append(allelesListExt[mutation[0]]+str(mutation[1])+allelesListExt[mutation[2]]+":"+str(mutation[3]))
					if iNode<(len(nextNode.mutations)-1):
						stringList.append(",")
				stringList.append("},Ns={")
				for iNode in range(len(nextNode.Ns)):
					mutation=nextNode.Ns[iNode]
					if type(mutation)==int:
						stringList.append(str(mutation))
					else:
						stringList.append(str(mutation[0])+"-"+str(mutation[1]))
					if iNode<(len(nextNode.Ns)-1):
						stringList.append(",")
				stringList.append("}")
				if usingErrorRate and (not nextNode.children) and node.errors:
					stringList.append(",errors={")
					for iNode in range(len(nextNode.errors)):
						mutation=nextNode.errors[iNode]
						stringList.append(allelesListExt[mutation[0]]+str(mutation[1])+allelesListExt[mutation[2]]+":"+str(mutation[3]))
						if iNode<(len(nextNode.errors)-1):
							stringList.append(",")
					stringList.append("}")
			if printIQtreeSupportForNode:
				stringList.append(",")
			else:
				stringList.append("]")
		elif nextNode.up==None and estimateMAT:
			stringList.append("[&rootState={")
			currentPos=0
			firstDoneEntry=False
			rootVect=rootVector(nextNode.probVect,False)
			for entry in rootVect:
				if entry[0]!=4:
					if not firstDoneEntry:
						firstDoneEntry=True
					else:
						stringList.append(",")
				if entry[0]==5:
					stringList.append("N"+str(currentPos+1)+"-"+str(entry[1]))
				elif entry[0]==6:
					vect=entry[-1]
					firstDone=False
					for i in range4:
						if vect[i]>minMutProb:
							if not firstDone:
								firstDone=True
							else:
								stringList.append(",")
							stringList.append(allelesList[i]+str(entry[1])+":"+str(vect[i]))
				elif entry[0]<4:
					stringList.append(allelesList[entry[0]]+str(entry[1])+":1.0")
				currentPos=entry[1]
			stringList.append("}")
			if printIQtreeSupportForNode:
				stringList.append(",")
			else:
				stringList.append("]")
		elif printIQtreeSupportForNode:
			stringList.append("[&")
		if printIQtreeSupportForNode:
			stringList.append("IQsupport="+str(nextNode.IQsupport))
			stringList.append("]")
	elif performLineageAssignment and (hasattr(nextNode,"lineage") or hasattr(nextNode,"lineages")):
		#print("creating lineage string for node")
		stringList.append("[&")
		if hasattr(nextNode,"lineage"):
			stringList.append("lineage="+nextNode.lineage)
			if hasattr(nextNode,"lineages"):
				stringList.append(",")
		if hasattr(nextNode,"lineages"):
			stringList.append("lineages={")
			for lineageName in nextNode.lineages.keys():
				stringList.append(lineageName+":"+str(nextNode.lineages[lineageName]))
				stringList.append(",")
			stringList.pop()
			stringList.append("}")
		stringList.append("]")
		#print("".join(stringList))
	return "".join(stringList)




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
				stringList.append(")"+nextNode.name)
				if aBayesPlusOn or estimateMAT or performLineageAssignment:
					name=""
					#if networkOutput:
					#	name=nextNode.name
					#else:
					#	name=""
					stringList.append(stringForNode(nextNode,name,nextNode.dist))
				#else:
				#	stringList.append(nextNode.name)
				if nextNode.dist:
					stringList.append(":"+str(nextNode.dist))
				else:
					stringList.append(":0.0")
				if nextNode.up!=None:
					if nextNode.up.children[0]==nextNode:
						direction=1
					else:
						direction=2
				nextNode=nextNode.up
		else:
			if len(nextNode.minorSequences)>0:
				stringList.append("(")
				if supportForIdenticalSequences or performLineageAssignment:
					stringList.append(stringForNode(nextNode,nextNode.name,0.0))
				else:
					stringList.append(nextNode.name)
				stringList.append(":0.0")
				for s2 in nextNode.minorSequences:
					stringList.append(",")
					if supportForIdenticalSequences or performLineageAssignment:
						stringList.append(stringForNode(nextNode,s2,0.0))
					else:
						stringList.append(s2)
					stringList.append(":0.0")
				stringList.append(")"+nextNode.name+"_MinorSeqsClade")
			else:
				stringList.append(nextNode.name)
			if aBayesPlusOn or estimateMAT or performLineageAssignment:
				#name=nextNode.name+"_MinorSeqsClade"
				#if nextNode.minorSequences and networkOutput:
				#	name=nextNode.name+"_MinorSeqsClade"
				#else:
				#	name=""
				stringList.append(stringForNode(nextNode,"",nextNode.dist))
			#else:
			#	stringList.append(nextNode.name+"_MinorSeqsClade")
			if nextNode.dist:
				stringList.append(":"+str(nextNode.dist))
			else:
				stringList.append(":0.0")
			if nextNode.up!=None:
				if nextNode.up.children[0]==nextNode:
					direction=1
				else:
					direction=2
			nextNode=nextNode.up
	stringList.append(";")
	return "".join(stringList)



#create newick string of a given tree (input node is assumed to be the root) - this time the generated tree is binary (polytomies are represented with branches of length 0).
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
				stringList.append(")"+nextNode.name)
				if aBayesPlusOn or estimateMAT or performLineageAssignment:
					name=""
					#if networkOutput:
					#	name=nextNode.name
					#else:
					#	name=""
					stringList.append(stringForNode(nextNode,name,nextNode.dist))
				#else:
				#	stringList.append(nextNode.name)
				if nextNode.dist:
					stringList.append(":"+str(nextNode.dist))
				else:
					stringList.append(":"+str(0.0))
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
				if supportForIdenticalSequences or performLineageAssignment:
					stringList.append(stringForNode(nextNode,nextNode.name,0.0))
				else:
					stringList.append(nextNode.name)
				stringList.append(":")
				for s2 in nextNode.minorSequences[:-1]:
					stringList.append("0.0,")
					if supportForIdenticalSequences or performLineageAssignment:
						stringList.append(stringForNode(nextNode,s2,0.0))
					else:
						stringList.append(s2)
					stringList.append(":0.0):")
				stringList.append("0.0,")
				if supportForIdenticalSequences or performLineageAssignment:
					stringList.append(stringForNode(nextNode,nextNode.minorSequences[-1],0.0))
				else:
					stringList.append(nextNode.minorSequences[-1])
				stringList.append(":0.0)"+nextNode.name+"_MinorSeqsClade")
			else:
				stringList.append(nextNode.name)
			if aBayesPlusOn or estimateMAT or performLineageAssignment:
				name=""
				#if nextNode.minorSequences and networkOutput:
				#	name=nextNode.name+"_MinorSeqsClade"
				#else:
				#	name=""
				stringList.append(stringForNode(nextNode,name,nextNode.dist))
			#else:
			#	stringList.append(nextNode.name+"_MinorSeqsClade")
			if nextNode.dist:
				stringList.append(":"+str(nextNode.dist))
			else:
				stringList.append(":"+str(0.0))
			if nextNode.up!=None:
				if nextNode.up.children[0]==nextNode:
					direction=1
				else:
					direction=2
			nextNode=nextNode.up
	stringList.append(";")
	print("created newick string for tree with "+str(numLeaves)+" leaves.")
	return "".join(stringList)

#count how many samples are in the (sub)tree rooted at "node".
def countTips(node):
	nextNode=node
	direction=0
	numSamples=0
	while nextNode!=None:
		if nextNode.children:
			if direction==0:
				nextNode=nextNode.children[0]
			elif direction==1:
				nextNode=nextNode.children[1]
				direction=0
			else:
				if nextNode.up!=None:
					if nextNode.up.children[0]==nextNode:
						direction=1
					else:
						direction=2
				nextNode=nextNode.up
		else:
			numSamples+=1+len(nextNode.minorSequences)
			if nextNode.up!=None:
				if nextNode.up.children[0]==nextNode:
					direction=1
				else:
					direction=2
			nextNode=nextNode.up
	return numSamples

#count how many samples are in the (sub)tree rooted at "node".
def writeTaxaNames(file,node):
	nextNode=node
	direction=0
	numSamples=0
	while nextNode!=None:
		if nextNode.children:
			if direction==0:
				nextNode=nextNode.children[0]
			elif direction==1:
				nextNode=nextNode.children[1]
				direction=0
			else:
				if nextNode.up!=None:
					if nextNode.up.children[0]==nextNode:
						direction=1
					else:
						direction=2
				nextNode=nextNode.up
		else:
			file.write("	"+nextNode.name+"\n")
			for samName in nextNode.minorSequences:
				file.write("	"+samName+"\n")
			if nextNode.up!=None:
				if nextNode.up.children[0]==nextNode:
					direction=1
				else:
					direction=2
			nextNode=nextNode.up
	return numSamples


#set all descendant nodes to dirty.
#So far these flags are used to prevent traversing the same part of the tree multiple times.
def giveInternalNodeNames(node):
	#print("giveInternalNodeNames")
	counter=1
	nextLeaves=[node]
	while nextLeaves:
		nextNode=nextLeaves.pop()
		#list of node names where possible alternative placements of the current node are
		nextNode.alternativePlacements=[]
		if nextNode.children:
			#if ((not hasattr(nextNode, "name")) or nextNode.name==""):
			nextNode.name="in"+str(counter)
				#print("New name given")
				#print(nextNode.name)
			counter+=1
			#else:
			#	print("Name not changed")
			#	print(nextNode.name)
		for c in nextNode.children:
			nextLeaves.append(c)
	#print("Counter")
	#print(counter)



topologyChanges=[0]
if writeTreesToFileEveryTheseSteps>0:
	lineSplit=outputFile.split("/")
	lineSplit[-1]=""
	outFolder="/".join(lineSplit)
	if not os.path.isdir(outFolder):
		print("Path to output file "+outFolder+" does not exist, quitting MAPLE. Use option --output to specify a valid output file path and output file name.")
		raise Exception("exit")
	if os.path.isfile(outputFile+"_intermediateTrees.tree")  and (not overwrite):
		print("File "+outputFile+"_intermediateTrees.tree already exists, quitting MAPLE. Use option --overwrite if you want to overwrite previous inference.")
		raise Exception("exit")
	intermediateTreesFile=open(outputFile+"_intermediateTrees.tree","w")
if writeLKsToFileEveryTheseSteps>0:
	lineSplit=outputFile.split("/")
	lineSplit[-1]=""
	outFolder="/".join(lineSplit)
	if not os.path.isdir(outFolder):
		print("Path to output file "+outFolder+" does not exist, quitting MAPLE. Use option --output to specify a valid output file path and output file name.")
		raise Exception("exit")
	if os.path.isfile(outputFile+"_intermediateLKs.txt")  and (not overwrite):
		print("File "+outputFile+"_intermediateLKs.txt already exists, quitting MAPLE. Use option --overwrite if you want to overwrite previous inference.")
		raise Exception("exit")
	intermediateLKsFile=open(outputFile+"_intermediateLKs.txt","w")


#run Robinson-Foulds distance calculations
if inputRFtrees!="":
	if not os.path.isfile(inputTree):
		print("Input tree in newick format "+inputTree+" not found, quitting MAPLE RF distance calculation. Use option --inputTree to specify a valid input newick tree file.")
		raise Exception("exit")
	if not os.path.isfile(inputRFtrees):
		print("Input trees in newick format "+inputRFtrees+" not found, quitting MAPLE RF distance calculation. Use option --inputRFtrees to specify valid file with input newick trees.")
		raise Exception("exit")
	lineSplit=outputFile.split("/")
	lineSplit[-1]=""
	outFolder="/".join(lineSplit)
	if not os.path.isdir(outFolder):
		print("Path to output file "+outFolder+" does not exist, quitting MAPLE RF calculation. Use option --output to specify a valid output file path and output file name.")
		raise Exception("exit")
	if os.path.isfile(outputFile+"_RFdistances.txt")  and (not overwrite):
		print("File "+outputFile+"_RFdistances.txt already exists, quitting MAPLE RF calculation. Use option --overwrite if you want to overwrite previous inference.")
		raise Exception("exit")
	tree1=readNewick(inputTree)[0]
	print("Read input newick tree")
	leafNameDict, nodeTable, leafCount, numBranches, leafDistDict, branchLengthDict, sumBranchLengths = prepareTreeComparison(tree1,rooted=False)
	otherTrees=readNewick(inputRFtrees,multipleTrees=multipleInputRFTrees)
	print("Read other input newick trees to be compared to the first one")
	file=open(outputFile+"_RFdistances.txt","w")
	file.write("RF\t"+"normalisedRF\t"+"leaves\t"+"foundBranches\t"+"missedBranches\t"+"notFoundBranches\t"+"RFL\n")
	for tree in otherTrees:
		numDiffs, normalisedRF, leafCount, foundBranches, missedBranches, notFoundBranches, RFL = RobinsonFouldsWithDay1985(tree,leafNameDict, nodeTable, leafCount, numBranches,leafDistDict, branchLengthDict, sumBranchLengths,rooted=False)
		file.write(str(numDiffs)+"\t"+str(normalisedRF)+"\t"+str(leafCount)+"\t"+str(foundBranches)+"\t"+str(missedBranches)+"\t"+str(notFoundBranches)+"\t"+str(RFL)+"\n")
	print("Comparison ended")
	file.close()
	raise Exception("exit")

#run lineage assignment of samples in input tree following the input references
if performLineageAssignment:
	if assignmentFile!="" and assignmentFileCSV!="":
		print("Please only use one between options --assignmentFile and --assignmentFileCSV .")
		raise Exception("exit")
	if (not os.path.isfile(inputNexusTree)) and (not os.path.isfile(inputTree)):
		print("Input tree in newick format "+inputTree+" or nexus format "+inputNexusTree+" not found, quitting MAPLE lineage assignment. Use option --inputTree or --inputNexusTree to specify a valid input tree file.")
		raise Exception("exit")
	if assignmentFile!="" and (not os.path.isfile(assignmentFile)):
		print("Input references assignments file "+assignmentFile+" not found, quitting MAPLE lineage assignment. Use option --assignmentFile to specify a valid file of reference sequences and their lineages.")
		raise Exception("exit")
	if assignmentFileCSV!="" and (not os.path.isfile(assignmentFileCSV)):
		print("Input references assignments file "+assignmentFileCSV+" not found, quitting MAPLE lineage assignment. Use option --assignmentFile to specify a valid csv file of reference sequences and their lineages.")
		raise Exception("exit")
	lineSplit=outputFile.split("/")
	lineSplit[-1]=""
	outFolder="/".join(lineSplit)
	if not os.path.isdir(outFolder):
		print("Path to output file "+outFolder+" does not exist, quitting MAPLE lineage assignment. Use option --output to specify a valid output file path and output file name.")
		raise Exception("exit")
	if os.path.isfile(outputFile+"_lineageAssignments.csv")  and (not overwrite):
		print("File "+outputFile+"_lineageAssignments.csv already exists, quitting MAPLE lineage assignment. Use option --overwrite if you want to overwirte previous inference.")
		raise Exception("exit")
	print("Reading input tree")
	if os.path.isfile(inputNexusTree):
		tree1=readNexus(inputNexusTree)
	else:
		tree1=readNewick(inputTree)[0]
	print("Input tree read")
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
	print("Input lineage definition read")
	file=open(outputFile+"_lineageAssignments.csv","w")

	node=tree1
	direction=0 #0=from parent, 1+ from children
	lineage=""
	mostAncestralLineages=[]
	allLineages=[]
	if os.path.isfile(inputNexusTree):
		uncertaintyFlag=True
		nodeDict={}
	else:
		giveInternalNodeNames(tree1)
		uncertaintyFlag=False
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
					mostAncestralLineages2=[]
					allLineages2=[]
					nodesToVisit=list(node.children)
					while nodesToVisit:
						nextNode=nodesToVisit.pop()
						if not nextNode.dist:
							if nextNode.children:
								for child in nextNode.children:
									nodesToVisit.append(child)
							elif nextNode.name in references:
								lineage=references[nextNode.name]
								allLineages2.append(lineage)
								indMost=0
								foundAncestor=False
								while indMost<len(mostAncestralLineages2):
									if mostAncestralLineages2[indMost] in lineage:
										foundAncestor=True
										break
									elif lineage in mostAncestralLineages2[indMost]:
										del mostAncestralLineages2[indMost]
									else:
										indMost+=1
								if not foundAncestor:
									mostAncestralLineages2.append(lineage)
					if len(mostAncestralLineages2)>0:
						lineage=mostAncestralLineages2[0]
						mostAncestralLineages=mostAncestralLineages2
						allLineages=allLineages2
				if lineage=="" or len(mostAncestralLineages)==0 or len(allLineages)==0:
					print(lineage)
					print(mostAncestralLineages)
					print(allLineages)
					print()
				node.lineage=lineage
				node.mostAncestralLineages=mostAncestralLineages
				node.allLineages=allLineages
				if uncertaintyFlag:
					nodeDict[node.name]=node
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
							raise Exception("exit")
					node=node.up
				else:
					lineage=node.lineage
					mostAncestralLineages=node.mostAncestralLineages
					allLineages=node.allLineages
					node=node.children[direction]
					direction=0
		
		#case of a sample
		else:
			if uncertaintyFlag:
				nodeDict[node.name]=node
				if node.name in references:
					node.lineage=references[node.name]
					if node.dist:
						node.mostAncestralLineages=[node.lineage]
						node.allLineages=[node.lineage]
					else:
						node.mostAncestralLineages=mostAncestralLineages
						node.allLineages=allLineages
				else:
					node.lineage=lineage
					node.mostAncestralLineages=mostAncestralLineages
					node.allLineages=allLineages
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
					raise Exception("exit")
			node=node.up

	print("Finished tree pass for lineage assignment")
	#If nexus tree is input, do one more round using alternative placements to inform lineage probabilities.
	#if os.path.isfile(inputNexusTree):
	if uncertaintyFlag:
		#for indTraverse in range(5):
		node=tree1
		numRefsTest=0
		numNoSupportTest=0
		numSupportTest=0
		lineage=""
		direction=0 #0=from parent, 1+ from children
		#for each internal node with positive bLen, check all descendants with a distance of 0 from it; if any of them is a reference,
		#then assign that lineage to the node and all its descendants unless the assignment is changed downstream.
		#If no reference at distance 0 is found, then use parent assignment.
		while node!=None:
			#case of an internal node
			if node.children:
				#coming from parent node
				if direction==0:
					#TODO modify so to account for overlapping references
					lineages={}
					if hasattr(node,"features") and ("support" in node.features):
					#if hasattr(node,"support"):
						for lin in node.allLineages:
							lineages[lin]=node.features["support"]/len(node.allLineages)
						#lineages[node.lineage]=node.features["support"]
						if "alternativePlacements" in node.features:
						#if hasattr(node,"alternativePlacements"):
							for alterPl in node.features["alternativePlacements"].keys():
								alterNode=nodeDict[alterPl]
								#alterLin=alterNode.lineage
								alterLins=alterNode.allLineages
								alterProb=node.features["alternativePlacements"][alterPl]/len(alterLins)
								for alterLin in alterLins:
									if alterLin in lineages:
										lineages[alterLin]+=alterProb
									else:
										lineages[alterLin]=alterProb
					else:
						for lin in node.allLineages:
							lineages[lin]=1.0/len(node.allLineages)
						#lineages[node.lineage]=1.0
					node.lineages=lineages
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
								raise Exception("exit")
						node=node.up
					else:
						lineage=node.lineage
						node=node.children[direction]
						direction=0
			
			#case of a sample
			else:
				lineages={}
				if node.name in references:
					file.write(node.name+","+references[node.name]+":1.0\n")
					lineages[references[node.name]]=1.0
					numRefsTest+=1
				else:
					if hasattr(node,"features") and ("support" in node.features):
					#if hasattr(node,"support"):
						for lin in node.allLineages:
							lineages[lin]=node.features["support"]/len(node.allLineages)
						numSupportTest+=1
						#lineages[node.lineage]=node.features["support"]
						if "alternativePlacements" in node.features:
						#if hasattr(node,"alternativePlacements"):
							for alterPl in node.features["alternativePlacements"].keys():
								alterNode=nodeDict[alterPl]
								#alterLin=alterNode.lineage
								alterLins=alterNode.allLineages
								alterProb=node.features["alternativePlacements"][alterPl]/len(alterNode.allLineages)
								for alterLin in alterLins:
									if alterLin in lineages:
										lineages[alterLin]+=alterProb
									else:
										lineages[alterLin]=alterProb
						file.write(node.name)
						for alterPl in lineages.keys():
							file.write(","+alterPl+":"+str(lineages[alterPl]))
						file.write("\n")
					else:
						numNoSupportTest+=1
						for lin in node.allLineages:
							lineages[lin]=1.0/len(node.allLineages)
						file.write(node.name)
						for alterPl in lineages.keys():
							file.write(","+alterPl+":"+str(lineages[alterPl]))
						file.write("\n")
						#file.write(node.name+","+lineage+":1.0\n")
						#lineages[node.lineage]=1.0
				node.lineages=lineages
				if node.up!=None:
					direction=-1
					for childNum in range(len(node.up.children)):
						if node==node.up.children[childNum]:
							direction=childNum+1
					if direction==-1:
						print("Error2, not found child node in parent node list!")
						raise Exception("exit")
				node=node.up
		print("Finished second tree pass for lineage assignment with uncertainty")
		#print(numRefsTest)
		#print(numNoSupportTest)
		#print(numSupportTest)
	print("Lineage assignment completed")
	file.close()
	if binaryTree:
		newickString=createBinaryNewick(tree1)
	else:
		newickString=createNewick(tree1)
	file=open(outputFile+"_nexusTree.tree","w")
	file.write("#NEXUS\nbegin taxa;\n	dimensions ntax="+str(countTips(tree1))+";\n	taxlabels\n")
	writeTaxaNames(file,tree1)
	file.write(";\nend;\n\nbegin trees;\n	tree TREE1 = [&R] ")
	file.write(newickString)
	file.write("\nend;\n")
	file.close()
	print("Output nexus tree with lineage assignments created.")
	exit()

	



minimumCarryOver=sys.float_info.min*(1e50)

if os.path.isfile(outputFile+"_tree.tree")  and (not overwrite):
	print("File "+outputFile+"_tree.tree already exists, quitting MAPLE tree inference. Use option --overwrite if you want to overwirte previous inference.")
	raise Exception("exit")
if not os.path.isfile(inputFile):
	print("Input file in Maple format "+inputFile+" not found, quitting MAPLE tree inference. Use option --input to specify a valid input file.")
	raise Exception("exit")
if refFile!="" and (not os.path.isfile(refFile)):
	print("Input reference fasta file "+refFile+" not found, quitting MAPLE tree inference. Use option --reference to specify a valid input reference file.")
	raise Exception("exit")
lineSplit=outputFile.split("/")
lineSplit[-1]=""
outFolder="/".join(lineSplit)
if not os.path.isdir(outFolder):
	print("Path to output file "+outFolder+" does not exist, quitting MAPLE tree inference. Use option --output to specify a valid output file path and output file name.")
	raise Exception("exit")
if inputTree!="":
	if not os.path.isfile(inputTree):
		print("Input tree in newick format "+inputTree+" not found, quitting MAPLE. Use option --inputTree to specify a valid input newick tree file.")
		raise Exception("exit")
	tree1=readNewick(inputTree,dirtiness=largeUpdate)[0]
	print("Read input newick tree")
	makeTreeBinary(tree1)

alleles={"A":0,"C":1,"G":2,"T":3}
allelesList=["A","C","G","T"]
allelesListExt=["A","C","G","T","?"]
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
				raise Exception("exit")
			else:
				entry=(linelist[0].lower(),int(linelist[1]))
			if ref[entry[1]-1]==entry[0] and entry[0]!="n" and entry[0]!="-":
				print("Mutation observed into reference nucleotide at position "+str(entry[1])+" , nucleotide "+entry[0]+". Wrong reference and/or diff file?")
				raise Exception("exit")
			if entry[1]<=pos:
				print("WARNING, at sample number "+str(nSeqs+1)+" found entry")
				print(line.replace("\n",""))
				print("which is inconsistent since the position is already represented by another entry:")
				print(seqList[-1])
				raise Exception("exit")
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

#vector to count how many bases of each type are cumulatively in the reference genome up to a certain position
cumulativeBases=[[0,0,0,0]]
for i in range(lRef):
	cumulativeBases.append(list(cumulativeBases[i]))
	#if (not (ref[i] in allelesList)) and (not (ref[i] in allelesListLow)):
	#	print("ERROR: ambiguity character in the reference genome. Exiting MAPLE.")
	#	raise Exception("exit")
	if (ref[i] in allelesList) or (ref[i] in allelesListLow):
		cumulativeBases[i+1][allelesUpOrLow[ref[i]]]+=1
rootFreqs=[0.0,0.0,0.0,0.0]
rootFreqsLog=[0.0,0.0,0.0,0.0]
for i in range(4):
	rootFreqs[i]=cumulativeBases[-1][i]/float(lRef)
	rootFreqsLog[i]=log(rootFreqs[i])
refIndeces=[]
for i in range(lRef):
	if (ref[i] in allelesList) or (ref[i] in allelesListLow):
		refIndeces.append(allelesLow[ref[i]])
	else:
		refIndeces.append(0)
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
def simplify(vec,refA):
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
		raise Exception("exit")
	if numA==1:
		if maxI==refA:
			return 4
		else:
			return maxI
	else:
		return 6



#Shorten genome list by merging together R entries that are mergeable - now extended also for the scenario with errors.
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
			elif abs(newVec[3]-entryOld[3])>thresholdProb:
				index+=1
				entryOld=vec[index]
			elif len(newVec)==4 or (newVec[4]==entryOld[4]):
				vec.pop(index)
			else:
				index+=1
				entryOld=vec[index]
		else:
			index+=1
			entryOld=vec[index]
	return






#define partial likelihood vector for a sample given its data - now extended to account also for sequence errors
def probVectTerminalNode(diffs,numMinorSeqs=0):
	if not errorRateSiteSpecific: 
		errorRate=errorRateGlobal
	if diffs is None:
		probVect=[(5,lRef)]
		return probVect
	pos=1
	probVect=[]
	for m in diffs:
			currPos=m[1]
			if currPos>pos: #region where the node with branch length bLen is identical to the ref.
				if usingErrorRate and (numMinorSeqs==0):
					probVect.append((4,currPos-1,True))
				else:
					probVect.append((4,currPos-1))
				pos=currPos
			if m[0]=="n" or m[0]=="-": #region with no info, store last position and length.
				if len(m)>2:
					length=m[2]
				else:
					length=1
				probVect.append((5,currPos+length-1))
				pos=currPos+length
			elif m[0] in allelesLow:
				#position at which node allele is sure but is different from the reference.
				if usingErrorRate and (numMinorSeqs==0):
					probVect.append((allelesLow[m[0]],currPos,True))
				else:
					probVect.append((allelesLow[m[0]],currPos))
				pos=currPos+1
			else:
				# non-"n" ambiguity character; for now interpret this as ambiguity instead of as a polymorphism.
				if onlyNambiguities:
					# if user asks to, to make things easier, interpret any ambiguity as an "n".
					probVect.append((5,currPos))
				else:
					#otherwise, store as an "other" scenario, where each nucleotide has its own partial likelihood.
					if usingErrorRate and (numMinorSeqs==0):
						ambigVect=list(ambiguities[m[0]])
						sumUnnormalizedVector =  sum([bool(item) for item in ambigVect])
						if errorRateSiteSpecific: errorRate = errorRates[currPos-1]
						if sumUnnormalizedVector ==2:
							for i in range4:             # M, instead of [0.5, 0.5, 0, 0] we will now get [0.5- ⅓ε, 0.5- ⅓ε,  ⅓ε,  ⅓ε]
								if ambigVect[i]==0:
									ambigVect[i] = errorRate*0.33333
								else: #if entry[-1][i]==0.5:
									ambigVect[i] -= errorRate*0.33333
						elif sumUnnormalizedVector == 3:
							for i in range4:             # for V instead of [⅓, ⅓, ⅓, 0] we get, [ ⅓ - ε/9,  ⅓ -ε/9,  ⅓ -ε/9,  ⅓ε]
								if ambigVect[i] == 0:
									ambigVect[i] = errorRate*0.33333
								else:  # if entry[-1][i]==⅓:
									ambigVect[i] -= errorRate/9 # ⅓ - ε/9
						probVect.append((6,currPos,ambigVect))
					else:
						probVect.append((6,currPos,ambiguities[m[0]]))
				pos=currPos+1
	if pos<=lRef:
		if usingErrorRate and (numMinorSeqs==0):
			probVect.append((4,lRef,True))
		else:
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
		raise Exception("exit")
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



# get the likelihoods as they change by moving along a branch
def getPartialVec(i12, totLen, mutMatrix, errorRate, vect=None, upNode=False, flag=False):
	"""
	When merging a CATGR-entry with an O entry, or for the case that T1=/=T2, it is useful to
	first create the partial vectors of the entries seperately, instead of working with one 'newVec'.
	With this function we create the partial likelihood vector for the considered entry:
	(if needed) we account for the error rates, and (if needed) we propagate the vector along its branch.
	:param i12: i1 or i2;  the nucleotide that was stored in the original entry
	:param flag: whether we should take into account the error rate
	:param totLen: branch length along which the vector should be 'mutated'
	:param upNode: should be True if the nucleotide was observed above the branch.
	:return: the new likelihood vector for this entry.
	"""
	if i12==6:
		if (not totLen):
			return list(vect)
		newVect=[]
		if upNode:
			for i in range4:
				tot=0.0
				for j in range4:
					tot+=mutMatrix[j][i]*vect[j]
				tot*=totLen
				tot+=vect[i]
				if tot<0:
					return [0.25,0.25,0.25,0.25]
				newVect.append(tot)
		else:
			for i in range4:
				tot=0.0
				for j in range4:
					tot+=mutMatrix[i][j]*vect[j]
				tot*=totLen
				tot+=vect[i]
				if tot<0:
					return [0.25,0.25,0.25,0.25]
				newVect.append(tot)
		return newVect
	elif usingErrorRate and flag:
		newVect = [errorRate*0.33333] * 4  # without error rate [0.0, 0.0, 0.0, 0.0]
		newVect[i12] = 1.0 - errorRate  # without error rate: 1.0
		if (not totLen):
			return newVect
		mutatedPartialVec = []
		for j in range4:
			tot = 0.0
			for i in range4:
				tot += mutMatrix[j][i] * newVect[i]
			tot *= totLen
			tot += newVect[j]
			if tot<0:
				return [0.25,0.25,0.25,0.25]
			mutatedPartialVec.append(tot)
		return mutatedPartialVec
	else:  # no flag (no error rate)
		if (not totLen):
			newVect=[0.0,0.0,0.0,0.0]
			newVect[i12]+=1.0
			return newVect
		newVect = []
		if upNode:
			for i in range4:
				newVect.append(mutMatrix[i12][i] * totLen)
		else:
			for i in range4:
				newVect.append(mutMatrix[i][i12] * totLen)
		newVect[i12]+=1.0
		if newVect[i12]<0:
			return [0.25,0.25,0.25,0.25]
		return newVect





#merge two partial likelihood vectors, one from above, probVect1, and one from below, probVect2
#unlike appendProb(), this function is not used on a large part of the tree at each placement, but only in a small neighbourhood;
def mergeVectorsUpDown(probVect1,bLenUp,probVect2,bLenDown):
	if not errorRateSiteSpecific: 
		errorRate=errorRateGlobal
	if not useRateVariation:
		mutMatrix=mutMatrixGlobal
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
				if usingErrorRate:
					if len(entry2)==3:
						probVect.append((entry2[0],pos,bLenDown,0.0,entry2[2]))
					elif len(entry2)==4:
						probVect.append((entry2[0],pos,entry2[2]+bLenDown,0.0,entry2[3]))
					else:
						if bLenDown:
							probVect.append((entry2[0],pos,bLenDown,0.0,False))
						else:
							probVect.append((entry2[0],pos))
				else:
					if len(entry2)>2:
						probVect.append((entry2[0],pos,entry2[2]+bLenDown,0.0))
					else:
						if bLenDown:
							probVect.append((entry2[0],pos,bLenDown,0.0))
						else:
							probVect.append((entry2[0],pos))
			else: # entry2 case "O", entry 1 is "N"
				if useRateVariation:
					mutMatrix=mutMatrices[pos]
				pos+=1
				totBLen=bLenDown
				if len(entry2)>3:
					totBLen+=entry2[2]
				newVect=[]
				if totBLen:
					newVect=getPartialVec(6, totBLen, mutMatrix, 0, vect=entry2[-1])
				else:
					newVect=list(entry2[-1])
				for i in range4:
					newVect[i]*=rootFreqs[i]
				totSum=sum(newVect)
				for i in range4:
					newVect[i]/=totSum
				probVect.append((6,pos,newVect))
		elif entry2[0]==5: #entry2 is N
			if entry1[0]<5:
				pos=min(entry1[1],entry2[1])
				if usingErrorRate:
					if len(entry1)==2:
						if bLenUp:
							probVect.append((entry1[0],pos,bLenUp,False))
						else:
							probVect.append((entry1[0],pos))
					elif len(entry1)==3:
						if bLenUp:
							probVect.append((entry1[0],pos,bLenUp,entry1[2]))
						else:
							probVect.append((entry1[0],pos,entry1[2]))
					elif len(entry1)==4:
						probVect.append((entry1[0],pos,entry1[2]+bLenUp, entry1[3]))
					else:
						probVect.append((entry1[0],pos,entry1[2],entry1[3]+bLenUp, entry1[4]))
				else:
					if len(entry1)==2:
						if bLenUp:
							probVect.append((entry1[0],pos,bLenUp))
						else:
							probVect.append((entry1[0],pos))
					elif len(entry1)==3:
						probVect.append((entry1[0],pos,entry1[2]+bLenUp))
					else:
						probVect.append((entry1[0],pos,entry1[2],entry1[3]+bLenUp))

			else: #entry1 is "O", entry 2 is "N"
				if useRateVariation:
					mutMatrix=mutMatrices[pos]
				pos+=1
				totBLen=bLenUp
				if len(entry1)==4:
					totBLen+=entry1[2]
				if totBLen:
					newVect=getPartialVec(6, totBLen, mutMatrix, 0, vect=entry1[-1],upNode=True)
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
			flag1=(usingErrorRate and (entry1[0]!=6) and (len(entry1)>2) and entry1[-1])
			flag2=(usingErrorRate and (entry2[0]!=6) and (len(entry2)>2) and entry2[-1])
			if errorRateSiteSpecific: errorRate = errorRates[pos]
			totLen1=bLenUp
			if entry1[0]>5:
				if len(entry1)>3:
					totLen1+=entry1[2]
			else:
				if len(entry1)>(2+usingErrorRate):
					totLen1+=entry1[2]
					if len(entry1)>3+usingErrorRate:
						totLen1+=entry1[3]
			totLen2=bLenDown
			if len(entry2)>(2+(usingErrorRate or entry2[0]==6)):
				totLen2+=entry2[2]

			if entry2[0]<5 and (not totLen2) and (not flag2): #due to 0 distance, the entry will be of same type as entry2
				if (not totLen1) and entry1[0]<5 and (not flag1):
					return None
				pos=min(entry1[1],entry2[1])
				probVect.append((entry2[0],pos))
			elif entry1[0]<5 and (not totLen1) and (not flag1): #due to 0 distance, the entry will be of same type as entry1
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
				if len(entry1)==4+usingErrorRate:
					newVec = getPartialVec(i1, entry1[2], mutMatrix, errorRate, flag=flag1)
					for i in range4:
						newVec[i]*=rootFreqs[i]
					if entry1[3]+bLenUp:
						newVec = getPartialVec(6, entry1[3]+bLenUp, mutMatrix, 0, vect=newVec, upNode=True)
				else:
					newVec = getPartialVec(i1, totLen1, mutMatrix, 0, upNode=True)

				if entry2[0]==6: #entry 2 is "O" and entry1 is a nucleotide
					if totLen2:
						newVec2 = getPartialVec(6, totLen2, mutMatrix, 0, vect=entry2[-1])
					else:
						newVec2 = entry2[-1]
					for j in range4:
						newVec[j]*=newVec2[j]
					sumV=sum(newVec)
					for i in range4:
						newVec[i]=newVec[i]/sumV
					state =simplify(newVec,refIndeces[pos])
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
					newVec2=getPartialVec(i2, totLen2, mutMatrix, errorRate, flag=flag2)
					for i in range4:
						newVec[i]*=newVec2[i]
					sumV=sum(newVec)
					if not sumV:
						print("0 sum in mergeVectorsUpDown")
						return None
					for i in range4:
						newVec[i]=newVec[i]/sumV
					pos+=1
					probVect.append((6,pos,newVec))
			else: #entry1[0]==6:
				if useRateVariation:
					mutMatrix=mutMatrices[pos]
				if totLen1:
					newVec = getPartialVec(6, totLen1, mutMatrix, 0, vect=entry1[-1],upNode=True)
				else:
					newVec=list(entry1[-1])

				if entry2[0]==6:
					if totLen2:
						newVec2 = getPartialVec(6, totLen2, mutMatrix, 0, vect=entry2[-1])
					else:
						newVec2=list(entry2[-1])
					for i in range4:
						newVec[i]*=newVec2[i]
				else: #entry1 is "O" and entry2 is a nucleotide
					if entry2[0]==4:
						i2=refIndeces[pos]
					else:
						i2=entry2[0]
					newVec2=getPartialVec(i2, totLen2, mutMatrix, errorRate, flag=flag2)
					for i in range4:
						newVec[i]*=newVec2[i]
				sumV=sum(newVec)
				if not sumV:
					return None
				for i in range4:
					newVec[i]=newVec[i]/sumV
				state =simplify(newVec,refIndeces[pos])
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

	#check if the final  probVect can be simplified by merging consecutive entries
	shorten(probVect)
	return probVect







#merge two lower child partial likelihood vectors to create a new one 
#(and also calculate the logLk of the merging if necessary, which is currently only useful for the root but could also be useful to calculate the overall total likelihood of the tree).
def mergeVectors(probVect1,bLen1,probVect2,bLen2,returnLK=False):
	if not errorRateSiteSpecific: 
		errorRate=errorRateGlobal
	if not useRateVariation:
		mutMatrix=mutMatrixGlobal
	indexEntry1, indexEntry2, pos, totalFactor = 0, 0, 0, 1.0
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
				if usingErrorRate:
					if len(entry2) == 2: # entry 2 is a child node or has bLength element of 0 with its predecessor.
						if bLen2: #To account for cases with 0 branch length elements, when the flag needs to be set to True.
							probVect.append((entry2[0], pos, bLen2, False))
						else: #0-branch length: no need to set the flag
							probVect.append((entry2[0], pos))
					elif len(entry2) == 3:
						if bLen2:
							probVect.append((entry2[0], pos, bLen2, entry2[2]))
						else:
							probVect.append((entry2[0], pos, entry2[2]))
					else:
						probVect.append((entry2[0], pos, entry2[2] + bLen2, entry2[3]))
				else:
					if len(entry2) == 2: # entry 2 is a child node or has bLength element of 0 with its predecessor.
						if bLen2: #To account for cases with 0 branch length elements, when the flag needs to be set to True.
							probVect.append((entry2[0], pos, bLen2)) #flag = flag2 = node2isleaf in this case. If node 2 was not a leaf, it means it is a node with a 0 branch length element.
						else: #0-branch length: no need to set the flag
							probVect.append((entry2[0], pos))
					else:
						probVect.append((entry2[0], pos, entry2[2] + bLen2))
			else: # case entry2 is "O" and entry1 is "N"
				pos+=1
				if len(entry2)==3:
					if bLen2:
						probVect.append((6,pos,bLen2,entry2[-1]))
					else:
						probVect.append((6,pos,entry2[-1]))
				else:
					probVect.append((6,pos,entry2[2]+bLen2,entry2[-1]))
		elif entry2[0]==5: #entry2 is N
			if entry1[0]<5:
				pos=min(entry1[1],entry2[1])
				if usingErrorRate:
					if len(entry1)==2:
						if bLen1:
							probVect.append((entry1[0],pos,bLen1,False))
						else:
							probVect.append((entry1[0],pos))
					elif len(entry1)==3:
						if bLen1:
							probVect.append((entry1[0],pos,bLen1,entry1[2]))
						else:
							probVect.append((entry1[0],pos,entry1[2]))
					else:
						probVect.append((entry1[0],pos,entry1[2]+bLen1,entry1[3]))
				else:
					if len(entry1)==2:
						if bLen1:
							probVect.append((entry1[0],pos,bLen1))
						else:
							probVect.append((entry1[0],pos))
					else:
						probVect.append((entry1[0],pos,entry1[2]+bLen1))
			else: #entry1 is "O" and entry2 is "N"
				pos+=1
				if len(entry1)==3:
					if bLen1:
						probVect.append((6,pos,bLen1,entry1[-1]))
					else:
						probVect.append((6,pos,entry1[-1]))
				else:
					probVect.append((6,pos,entry1[2]+bLen1,entry1[-1]))
				
		else: #entry1 and entry2 are not "N"
			totLen1=bLen1
			if len(entry1)>(2+(usingErrorRate or entry1[0]==6)):
				totLen1+=entry1[2]
			totLen2=bLen2
			if len(entry2)>(2+(usingErrorRate or entry2[0]==6)):
				totLen2+=entry2[2]

			flag1=(usingErrorRate and (entry1[0]!=6) and (len(entry1)>2) and entry1[-1])
			flag2=(usingErrorRate and (entry2[0]!=6) and (len(entry2)>2) and entry2[-1])

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
					if flag1 or flag2:
						if errorRateSiteSpecific: cumErrorRate = cumulativeErrorRate[pos] - cumulativeErrorRate[end] #both this and the next one should end up being negative
						else: cumErrorRate = errorRate * (pos-end) 
						cumulPartLk += cumErrorRate * (flag1 + flag2)
				pos=end
                
			elif (not totLen1) and (not totLen2) and entry1[0]<5 and entry2[0]<5 and (not flag1) and (not flag2): #0 distance between different nucleotides: merge is not possible
				if returnLK:
					print("mergeVectors() returning None 1")
					raise Exception("exit")
				else:
					return None
			elif entry1[0]<5: #entry1 is a nucleotide
				if errorRateSiteSpecific: errorRate = errorRates[pos]
				if useRateVariation:
					mutMatrix=mutMatrices[pos]
				if entry1[0]==4:
					i1=refIndeces[pos]
				else:
					i1=entry1[0]
				if totLen1 or flag1:
					newVec=getPartialVec(i1, totLen1, mutMatrix, errorRate, flag=flag1)
				else:
					newVec=[0.0,0.0,0.0,0.0]
					newVec[i1]=1.0

				if entry2[0]==6: #entry1 is a nucleotide and entry2 is "O"
					if totLen2:
						newVec2=getPartialVec(6, totLen2, mutMatrix, 0, vect=entry2[-1])
					else:
						newVec2=entry2[-1]
					for j in range4:
						newVec[j]*=newVec2[j]
					sumV=sum(newVec)
					if not sumV:
						if returnLK:
							print("mergeVectors() returning None 2")
							raise Exception("exit")
						else:
							return None
					for i in range4:
						newVec[i]=newVec[i]/sumV
					state =simplify(newVec,refIndeces[pos])
					pos+=1
					if state==6:
						probVect.append((6,pos,newVec))
					else:
						probVect.append((state,pos))
					if returnLK:
						totalFactor*=sumV
				else: #entry1 and entry2 are nucleotides
					if entry2[0]==4:
						i2=refIndeces[pos]
					else:
						i2=entry2[0]
					if totLen2 or flag2:
						partialVec2 = getPartialVec(i2, totLen2, mutMatrix, errorRate, flag=flag2)
						for i in range4:
							newVec[i]*=partialVec2[i]
						sumV=sum(newVec)
						for i in range4:
							newVec[i]=newVec[i]/sumV
						state =simplify(newVec,refIndeces[pos])
						pos+=1
						if state==6:
							probVect.append((6,pos,newVec))
						else:
							probVect.append((state,pos))
						if returnLK:
							if sumV<=0.0:
								print("mergeVectors() different nucleotides, 0 or negative LK?")
								raise Exception("exit")
							totalFactor*=sumV
					else:
						pos+=1
						probVect.append((entry2[0],pos))
						if returnLK:
							if newVec[i2]<=0:
								print("mergeVectors() different nucleotides, 0 bLen, 0 LK?")
								raise Exception("exit")
							totalFactor*=newVec[i2]
				
			else: #entry1[0]==6:
				if useRateVariation:
					mutMatrix=mutMatrices[pos]
				if totLen1:
					newVec=getPartialVec(6, totLen1, mutMatrix, 0, vect=entry1[-1])
				else:
					newVec=list(entry1[-1])

				if entry2[0]==6:
					if totLen2:
						newVec2=getPartialVec(6, totLen2, mutMatrix, 0, vect=entry2[-1])
					else:
						newVec2=entry2[-1]
					for i in range4:
						newVec[i]*=newVec2[i]
					sumV=sum(newVec)
					if not sumV:
						if returnLK:
							print("mergeVectors() returning None 3")
							raise Exception("exit")
						else:
							return None
					for i in range4:
						newVec[i]=newVec[i]/sumV
					state =simplify(newVec,refIndeces[pos])
					pos+=1
					if state==6:
						probVect.append((6,pos,newVec))
					else:
						probVect.append((state,pos))
					if returnLK:
						totalFactor*=sumV
				else: #entry2 is a nucleotide and entry1 is "O"
					if entry2[0]==4:
						i2=refIndeces[pos]
					else:
						i2=entry2[0]
					if totLen2 or flag2:
						if errorRateSiteSpecific: errorRate = errorRates[pos]
						partialVec2 = getPartialVec(i2, totLen2, mutMatrix, errorRate, flag=flag2)
						for i in range4:
							newVec[i]*=partialVec2[i]
						sumV=sum(newVec)
						for i in range4:
							newVec[i]=newVec[i]/sumV
						state =simplify(newVec,refIndeces[pos])
						pos+=1
						if state==6:
							probVect.append((6,pos,newVec))
						else:
							probVect.append((state,pos))
						if returnLK:
							totalFactor*=sumV
					else:
						if not newVec[i2]:
							if returnLK:
								print("mergeVectors() returning None 4")
								raise Exception("exit")
							else:
								return None
						pos+=1
						probVect.append((entry2[0],pos))
						if returnLK:
							totalFactor*=newVec[i2]

		if returnLK and totalFactor<=minimumCarryOver:
			try:
				if totalFactor<sys.float_info.min:
					print("In mergeVectors() too small LK")
					raise Exception("exit")
			except:
				print("In mergeVectors() value error")
				raise Exception("exit")
			cumulPartLk+=log(totalFactor)
			totalFactor=1.0

		if pos==lRef:
			break
		if pos==entry1[1]:
			indexEntry1+=1
			entry1=probVect1[indexEntry1]
		if pos==entry2[1]:
			indexEntry2+=1
			entry2=probVect2[indexEntry2]

	#check if the final  probVect can be simplified by merging consecutive entries
	shorten(probVect)

	if returnLK:
		return probVect, cumulPartLk+log(totalFactor)
	else:
		return probVect











#calculate the probability that results from combining a lower likelihood genome list of the root with root frequencies.
#Could maybe be combined with function rootVector() ?
def findProbRoot(probVect):
	if not errorRateSiteSpecific: 
		errorRate=errorRateGlobal
	logLK=0.0
	logFactor=1.0
	pos=0
	for entry in probVect:
		if usingErrorRate and (entry[0]<5) and (len(entry)>2) and entry[-1]:
			if entry[0]==4:
				logLK+=rootFreqsLogErrorCumulative[entry[1]]-rootFreqsLogErrorCumulative[pos]
			else:
				if errorRateSiteSpecific: errorRate = errorRates[pos]
				logFactor*=(rootFreqs[entry[0]]*(1.0-1.33333*errorRate)+0.33333*errorRate)
		else:
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
def rootVector(probVect,bLen):
	newProbVect=[]
	for entry in probVect:
		if entry[0]==5:
			newProbVect.append(entry)
		elif entry[0]==6:
			totBLen=bLen
			if len(entry)>3:
				totBLen+=entry[2]
			if totBLen:
				if useRateVariation:
					newVect=getPartialVec(6, totBLen, mutMatrices[entry[1]-1], 0, vect=entry[-1])
				else:
					newVect=getPartialVec(6, totBLen, mutMatrixGlobal, 0, vect=entry[-1])
				for i in range4:
					newVect[i]*=rootFreqs[i]
			else:
				newVect=[]
				for i in range4:
					newVect.append(entry[-1][i]*rootFreqs[i])
			totSum=sum(newVect)
			for i in range4:
				newVect[i]/=totSum
			newProbVect.append((6,entry[1],newVect))
		else:
			if usingErrorRate:
				flag1 = ((len(entry)>2) and entry[-1])
				if len(entry)>3:
					newProbVect.append((entry[0],entry[1],entry[2]+bLen,0.0,flag1))
				else:
					if bLen or flag1:
						newProbVect.append((entry[0],entry[1],bLen,0.0,flag1))
					else:
						newProbVect.append((entry[0],entry[1]))
			else:
				if len(entry)==3:
					newProbVect.append((entry[0],entry[1],entry[2]+bLen,0.0))
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
def estimateBranchLengthWithDerivative(probVectP,probVectC):
	if not errorRateSiteSpecific: 
		errorRate=errorRateGlobal
	if not useRateVariation:
		mutMatrix=mutMatrixGlobal
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
		elif entry1[0]==5: # case entry2 is N
			#if parent node is type "N", in theory we might have to calculate the contribution of root nucleotides; 
			# however, if this node is "N" then every other node in the current tree is "N", so we can ignore this since this contribution cancels out in relative terms.
			pos=min(entry1[1],entry2[1])
		else:
			#below, when necessary, we represent the likelihood as coeff0*l +coeff1, where l is the branch length to be optimized.
			if entry1[0]==4 and entry2[0]==4: # case entry1 is R	
				end=min(entry1[1],entry2[1])
				c1+=(cumulativeRate[end]-cumulativeRate[pos])
				pos=end
			else:
				if useRateVariation:
					mutMatrix=mutMatrices[pos]
				flag1 = (usingErrorRate and (entry1[0]!=6) and (len(entry1)>2) and entry1[-1])
				flag2 = (usingErrorRate and (entry2[0]!=6) and (len(entry2)>2) and entry2[-1])
				if errorRateSiteSpecific: errorRate = errorRates[pos]

				#contribLength will be here the total length from the root or from the upper node, down to the down node.
				contribLength=False
				if entry1[0]<5:
					if len(entry1)==3+usingErrorRate:
						contribLength=entry1[2]
					elif len(entry1)==4+usingErrorRate:
						contribLength=entry1[3]
				else:
					if len(entry1)>3:
						contribLength=entry1[2]
				if entry2[0]<5:
					if len(entry2)>2+usingErrorRate:
						contribLength+=entry2[2]
				else:
					if len(entry2)>3:
						contribLength+=entry2[2]

				if entry1[0]==4:
					#entry1 is reference and entry2 is of type "O"
					if entry2[0]==6:
						i1=refIndeces[pos]
						if len(entry1)==(4+usingErrorRate):
							coeff0=rootFreqs[i1]*entry2[-1][i1] 
							coeff1=0.0
							for i in range4:
								coeff0+=rootFreqs[i]*mutMatrix[i][i1]*entry1[2]*entry2[-1][i]
								coeff1+=mutMatrix[i1][i]*entry2[-1][i]
							coeff1*=rootFreqs[i1]
							if contribLength:
								coeff0+=coeff1*contribLength
							if flag1:
								coeff0-=1.33333*errorRate*rootFreqs[i1]*entry2[-1][i1]
								for i in range4:
									coeff0+=rootFreqs[i]*entry2[-1][i]*0.33333*errorRate
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
						if len(entry1)==4+usingErrorRate:
							i1=refIndeces[pos]
							i2=entry2[0]
							coeff0=rootFreqs[i2]*mutMatrix[i2][i1]*entry1[2]
							if contribLength:
								coeff0+=rootFreqs[i1]*mutMatrix[i1][i2]*contribLength
							if flag2:
								coeff0+=rootFreqs[i1]*0.33333*errorRate
							if flag1:
								coeff0+=rootFreqs[i2]*0.33333*errorRate
							coeff1=rootFreqs[i1]*mutMatrix[i1][i2]
							if coeff1:
								coeff0=coeff0/coeff1
							else:
								coeff0=None
						else:
							coeff0=contribLength
							if flag2:
								if mutMatrix[refIndeces[pos]][entry2[0]]:
									coeff0+=errorRate*0.33333/mutMatrix[refIndeces[pos]][entry2[0]]
								else:
									coeff0=None
						if coeff0!=None:	
							if coeff0:
								ais.append(coeff0)
							else:
								nZeros+=1
						pos+=1

				# entry1 is of type "O"
				elif entry1[0]==6:
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
						if flag2:
							coeff0+=errorRate*0.33333
					if coeff1<0.0:
						c1+=coeff1/coeff0
					elif coeff1:
						coeff0=coeff0/coeff1
						ais.append(coeff0)
					pos+=1

				else: #entry1 is a non-ref nuc
					if entry2[0]==entry1[0]:
						c1+=mutMatrix[entry1[0]][entry1[0]]
					else: #entry1 is a nucleotide and entry2 is not the same as entry1
						i1=entry1[0]
						if entry2[0]<5: #entry2 is a nucleotide
							if entry2[0]==4:
								i2=refIndeces[pos]
							else:
								i2=entry2[0]

							if len(entry1)==4+usingErrorRate:
								coeff0=rootFreqs[i2]*mutMatrix[i2][i1]*entry1[2]
								if contribLength:
									coeff0+=rootFreqs[i1]*mutMatrix[i1][i2]*contribLength
								if flag2:
									coeff0+=rootFreqs[i1]*0.33333*errorRate
								if flag1:
									coeff0+=rootFreqs[i2]*0.33333*errorRate
								coeff1=rootFreqs[i1]*mutMatrix[i1][i2]
								if coeff1:
									coeff0=coeff0/coeff1
								else:
									coeff0=None
							else:
								coeff0=contribLength
								if flag2:
									coeff0+=errorRate*0.33333/mutMatrix[i1][i2]
							if coeff0!=None:
								if coeff0:
									ais.append(coeff0)
								else:
									nZeros+=1

						else: #entry1 is a nucleotide and entry2 is of type "O"
							if len(entry1)==4+usingErrorRate:
								coeff0=rootFreqs[i1]*entry2[-1][i1] 
								coeff1=0.0
								for i in range4:
									coeff0+=rootFreqs[i]*mutMatrix[i][i1]*entry1[2]*entry2[-1][i]
									coeff1+=mutMatrix[i1][i]*entry2[-1][i]
								coeff1*=rootFreqs[i1]
								if contribLength:
									coeff0+=coeff1*contribLength
								if flag1:
									coeff0-=1.33333*errorRate*rootFreqs[i1]*entry2[-1][i1]
									for i in range4:
										coeff0+=rootFreqs[i]*entry2[-1][i]*0.33333*errorRate
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
		print("Initial border parameters don't fit expectations")
	
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
def updateBLen(node,addToList=False,nodeList=None):
	cNode=node
	node=node.up
	if cNode==node.children[0]:
		vectUp=node.probVectUpRight
		cNum=0
	else:
		vectUp=node.probVectUpLeft
		cNum=1
	bestLength=estimateBranchLengthWithDerivative(vectUp,cNode.probVect)
	cNode.dist=bestLength
	node.dirty=True
	cNode.dirty=True
	if addToList:
		nodeList.append((cNode,2))
		nodeList.append((node,cNum))




# update the partials iteratively starting from the nodes in nodeList
#each entry in nodeList contains the node it refers to, and the direction where the update comes from (0 is left child, 1 is right child, 2 is parent)
def updatePartials(nodeList):
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
				newTot=mergeVectorsUpDown(vectUpUp,node.dist/2,node.probVect,node.dist/2)
				if newTot==None:
					if node.dist>1e-100:
						print("inside updatePartials(), from parent: should not have happened since node.dist>0")
					updateBLen(node)
					nodeList.append((node.up,childNumUp))
					newTot=mergeVectorsUpDown(vectUpUp,node.dist/2,node.probVect,node.dist/2)
					madeChange=True
				node.probVectTotUp=newTot
			else:
				node.probVectTotUp=None
			if len(node.children)>0:# and (not updatedBLen): #at valid internal node, update upLeft and upRight, and if necessary add children to nodeList.
				child0Vect=node.children[0].probVect
				child1Vect=node.children[1].probVect
				dist0=node.children[0].dist
				dist1=node.children[1].dist
				newUpRight=mergeVectorsUpDown(vectUpUp,node.dist,child1Vect,dist1)
				if newUpRight==None:
					if (not node.dist) and (not dist1):
						updateBLen(node)
						if not node.dist:
							updateBLen(node.children[1],addToList=True,nodeList=nodeList)
							updatedBLen=True
						else:
							node.probVectTotUp=mergeVectorsUpDown(vectUpUp,node.dist/2,node.probVect,node.dist/2)
							newUpRight=mergeVectorsUpDown(vectUpUp,node.dist,child1Vect,dist1)
							nodeList.append((node.up,childNumUp))
							madeChange=True
					else:
						print("Strange: None vector from non-zero distances in updatePartials() from parent direction.")
						raise Exception("exit")
				if not updatedBLen:
					newUpLeft=mergeVectorsUpDown(vectUpUp,node.dist,child0Vect,dist0)
					if newUpLeft==None:
						if (not node.dist) and (not dist0) :
							updateBLen(node)
							if not node.dist:
								updateBLen(node.children[0],addToList=True,nodeList=nodeList)
								updatedBLen=True
							else:
								node.probVectTotUp=mergeVectorsUpDown(vectUpUp,node.dist/2,node.probVect,node.dist/2)
								newUpRight=mergeVectorsUpDown(vectUpUp,node.dist,child1Vect,dist1)
								newUpLeft=mergeVectorsUpDown(vectUpUp,node.dist,child0Vect,dist0)
								nodeList.append((node.up,childNumUp))
								madeChange=True
						else:
							print("Strange: None vector from non-zero distances in updatePartials() from parent direction, child0.")
							raise Exception("exit")
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
			newVect=mergeVectors(otherChildVect,otherChildDist,probVectDown,childDist)
			if newVect==None:
				if (not childDist) and (not otherChildDist):
					updateBLen(node.children[childNum])
					if not node.children[childNum].dist:
						updateBLen(node.children[otherChildNum],addToList=True,nodeList=nodeList)
						updatedBLen=True
					else:
						childDist=node.children[childNum].dist
						node.probVect=mergeVectors(otherChildVect,otherChildDist,probVectDown,childDist)
						nodeList.append((node.children[childNum],2))
						madeChange=True
				else:
					print("Strange: None vector from non-zero distances in updatePartials() from child direction.")
					raise Exception("exit")
			else:
				try:
					oldProbVect=node.probVect
				except AttributeError:
					oldProbVect=None
				node.probVect=newVect

			#update total mid-branches likelihood
			if (not updatedBLen) and node.dist and (node.up != None) and (vectUpUp!=None):
				newTot=mergeVectorsUpDown(vectUpUp,node.dist/2,node.probVect,node.dist/2)
				if newTot==None:
					updateBLen(node)
					node.probVect=mergeVectors(otherChildVect,otherChildDist,probVectDown,childDist)
					nodeList.append((node.children[childNum],2))
					node.probVectTotUp=mergeVectorsUpDown(vectUpUp,node.dist/2,node.probVect,node.dist/2)
					madeChange=True
					print("inside updatePartials(), from child: should not have happened since node.dist>0")
				else:
					node.probVectTotUp=newTot
					
			if (not updatedBLen) and (otherVectUp!=None):
				#update likelihoods at sibling node
				if node.up != None:
					newUpVect=mergeVectorsUpDown(vectUpUp,node.dist,probVectDown,childDist)
				else:
					newUpVect=rootVector(probVectDown,childDist)
				if newUpVect==None:
					if (not node.dist) and (not childDist):
						updateBLen(node)
						if not node.dist:
							updateBLen(node.children[childNum],addToList=True,nodeList=nodeList)
							updatedBLen=True
						else:
							node.probVectTotUp=mergeVectorsUpDown(vectUpUp,node.dist/2,node.probVect,node.dist/2)
							nodeList.append((node.children[childNum],2))
							madeChange=True
							newUpVect=mergeVectorsUpDown(vectUpUp,node.dist,probVectDown,childDist)
					else:
						print("Strange: None vector from non-zero distances in updatePartials() from child direction, newUpVect.")
						raise Exception("exit")
			if (not updatedBLen) and (otherVectUp!=None):
				if madeChange or areVectorsDifferent(otherVectUp,newUpVect):
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
				if len(entry1)>3:
					if abs(entry1[3] - entry2[3])>thresholdProb:
						return True
					if len(entry1)>4:
						if abs(entry1[4] - entry2[4])>thresholdProb:
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
					if len(entry1)>3:
						if abs(entry1[3] - entry2[3])>thresholdProb:
							print('diff bLen2 or flag case')
							return True
						if len(entry1)>4:
							if abs(entry1[4] - entry2[4])>thresholdProb:
								print('diff flag case')
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
# this is useful once the matrix estimation has finished, to make sure all genome lists reflect this matrix. 
def reCalculateAllGenomeLists(root, checkExistingAreCorrect=False,countNodes=False,countPseudoCounts=False,pseudoMutCounts=None,data=None,firstSetUp=False,checkSamplesIntree=False):
	#first pass to update all lower likelihoods.
	node=root
	#direction 0 means from parent, direction 1 means from a child
	lastNode=None
	direction=0
	dataNamesConverted=False
	if checkSamplesIntree:
		samplesInTree=set()
	while node!=None:
		if direction==0:
			if node.children:
				node=node.children[0]
			else:
				if data!=None:
					if (not extractNamesFlag) or (extractNamesFlag and (node.name in data)):
						if usingErrorRate:
							node.probVect=probVectTerminalNode(data[node.name],numMinorSeqs=len(node.minorSequences))
						else:
							node.probVect=probVectTerminalNode(data[node.name])
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
							if usingErrorRate:
								node.probVect=probVectTerminalNode(data[node.name],numMinorSeqs=len(node.minorSequences))
							else:
								node.probVect=probVectTerminalNode(data[node.name])
						else:
							print("Error: sample name "+node.name+" not found in the input sequence data - all samples in the input tree need to have a corresponding sequence entry.")
							raise Exception("exit")
					else:
						print("Error: sample name "+node.name+" not found in the input sequence data - all samples in the input tree need to have a corresponding sequence entry.")
						raise Exception("exit")
					if checkSamplesIntree:
						samplesInTree.add(node.name)

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
				newLower=mergeVectors(node.children[0].probVect,node.children[0].dist,node.children[1].probVect,node.children[1].dist)
				if checkExistingAreCorrect:
					if areVectorsDifferentDebugging(newLower,node.probVect):
						print("Inside reCalculateAllGenomeLists(), new lower at node is different from the old one, and it shouldn't be.")
						raise Exception("exit")
				if newLower==None:
					if (not node.children[0].dist) and (not node.children[1].dist):
						if firstSetUp:
							node.children[0].dist=oneMutBLen/2
							node.children[1].dist=oneMutBLen/2
						else:
							updateBLen(node.children[0])
							if (not node.children[0].dist):
								updateBLen(node.children[1])
						node.probVect=mergeVectors(node.children[0].probVect,node.children[0].dist,node.children[1].probVect,node.children[1].dist)
						if node.probVect==None:
							node.children[0].dist=oneMutBLen/2
							node.children[1].dist=oneMutBLen/2
							node.probVect=mergeVectors(node.children[0].probVect,node.children[0].dist,node.children[1].probVect,node.children[1].dist)
							if node.probVect==None:
								print("None vector when merging two vectors in reCalculateAllGenomeLists(), despite updating branch lengths.")
								raise Exception("exit")
					else:
						print("Strange, distances>0 but inconsistent lower genome list creation in reCalculateAllGenomeLists()")
						raise Exception("exit")
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
		newVect=rootVector(node.children[1].probVect,node.children[1].dist)
		if checkExistingAreCorrect:
			if areVectorsDifferentDebugging(newVect,node.probVectUpRight):
				print("new probVectUpRight at root is different from the old one, and it shouldn't be.")
				raise Exception("exit")
		node.probVectUpRight=newVect
		newVect=rootVector(node.children[0].probVect,node.children[0].dist)
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
				raise Exception("exit")
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
					newVect=mergeVectorsUpDown(vectUp,node.dist/2,node.probVect,node.dist/2)
					if checkExistingAreCorrect:
						if areVectorsDifferentDebugging(newVect,node.probVectTotUp):
							print("new probVectTotUp at node is different from the old one, and it shouldn't be.")
							print(newVect)
							print(node.probVectTotUp)
							raise Exception("exit")
					node.probVectTotUp=newVect
				if node.children:
					newUpRight=mergeVectorsUpDown(vectUp,node.dist,node.children[1].probVect,node.children[1].dist)
					if newUpRight==None:
						if (not node.children[1].dist) and (not node.dist):
							#nodeList=[]
							updateBLen(node)
							if not node.dist:
								if firstSetUp:
									node.probVectUpLeft=mergeVectorsUpDown(vectUp,node.dist,node.children[0].probVect,node.children[0].dist)
								updateBLen(node.children[1])
								totNodeList.append((node,1))
							else:
								node.probVectTotUp=mergeVectorsUpDown(vectUp,node.dist/2,node.probVect,node.dist/2)
								totNodeList.append((node.up,nodeChildNum))
							node.probVectUpRight=mergeVectorsUpDown(vectUp,node.dist,node.children[1].probVect,node.children[1].dist)
						else:
							print("Strange, distances>0 but inconsistent upRight genome list creation in reCalculateAllGenomeLists()")
							raise Exception("exit")
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
								raise Exception("exit")
						node.probVectUpRight=newUpRight
					newUpLeft=mergeVectorsUpDown(vectUp,node.dist,node.children[0].probVect,node.children[0].dist)
					if newUpLeft==None:
						if(not node.children[0].dist) and (not node.dist):
							updateBLen(node.children[0])
							if not node.children[0].dist:
								updateBLen(node)
								totNodeList.append((node.up,nodeChildNum))
								node.probVectTotUp=mergeVectorsUpDown(vectUp,node.dist/2,node.probVect,node.dist/2)
								node.probVectUpRight=mergeVectorsUpDown(vectUp,node.dist,node.children[1].probVect,node.children[1].dist)
							else:
								totNodeList.append((node,0))
							node.probVectUpLeft=mergeVectorsUpDown(vectUp,node.dist,node.children[0].probVect,node.children[0].dist)
						else:
							print("Strange, distances>0 but inconsistent upLeft genome list creation in reCalculateAllGenomeLists()")
							raise Exception("exit")
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
								raise Exception("exit")
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

		updatePartials(totNodeList)

	if checkSamplesIntree:
		return samplesInTree





#Initialize the mutation rate matrix
#If a fixed rate matrix is needed for SARS-CoV-2, this is the nucleotide mutation rate matrix from De Maio et al 2021 (now obsolete)
#mutMatrix=[[0.0,0.039,0.310,0.123],[0.140,0.0,0.022,3.028],[0.747,0.113,0.0,2.953],[0.056,0.261,0.036,0.0]]
pseudoMutCounts=[[0.0,1.0,5.0,2.0],[2.0,0.0,1.0,40.0],[5.0,2.0,0.0,20.0],[2.0,3.0,1.0,0.0]]
if model=="JC":
	mutMatrixGlobal=[[-1.0,1.0/3,1.0/3,1.0/3],[1.0/3,-1.0,1.0/3,1.0/3],[1.0/3,1.0/3,-1.0,1.0/3],[1.0/3,1.0/3,1.0/3,-1.0]]
else:
	mutMatrixGlobal=[[0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0]]
	updateSubMatrix(pseudoMutCounts,model,mutMatrixGlobal)

nonMutRates=[0,0,0,0]
for i in range4:
	nonMutRates[i]=mutMatrixGlobal[i][i]

cumulativeRate=[0.0]
for i in range(lRef):
	ind=refIndeces[i]
	cumulativeRate.append(cumulativeRate[-1]+nonMutRates[ind])



#In case an input tree is given, calculate all genome lists, pseudocounts and rates, and then recalculate genome lists again according to the new rates.
if inputTree!="":
	samplesAlreadyInTree=reCalculateAllGenomeLists(tree1,countPseudoCounts=True,pseudoMutCounts=pseudoMutCounts,data=data,firstSetUp=True,checkSamplesIntree=True)
	if model!="JC" and updateSubMatrix(pseudoMutCounts,model,mutMatrixGlobal):
		for i in range(lRef):
			cumulativeRate[i+1]=cumulativeRate[i]+nonMutRates[refIndeces[i]]
		for i in range4:
			nonMutRates[i]=mutMatrixGlobal[i][i]
	reCalculateAllGenomeLists(tree1)
	print("Genome list for initial tree and initial pseudocounts calculated.")
else:
	samplesAlreadyInTree=set()


#Sort samples based on distance from reference, but punishing more isolated N's and ambiguity characters.
#more ambiguous sequences are placed last this way - this is useful since ambiguous sequences are harder to place and are more likely to be less informative (and so be removed from the analysis altogether)
#def distancesFromRefPunishNs(data,samples):
def distancesFromRefPunishNs(data,samples=None,samplesInInitialTree=set()):
	sampleDistances=[]
	if samples==None:
		rangeInd=range(len(data))
	else:
		rangeInd=samples
	#for sample in samples:
	for diffIndex in rangeInd:
		if (samples==None) or (not (diffIndex in samplesInInitialTree)):
			diffs=data[diffIndex]
			pos=1
			comparisons=0
			diffNum=0
			for m in diffs:
				currPos=m[1]
				if currPos>pos: #region identical to the reference
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
			if not comparisons:
				print("\n WARNING!!!!!!!!!!!\n\n Sample number "+str(diffIndex+1)+" appears to be completely non-informative. It will still be included in the MAPLE analysis, but please consider removing it from the alignment.")
				print(diffs)
			else:
				totDivFromRef[0]+=float(diffNum)/comparisons
				if ((float(diffNum)/comparisons)>0.1) and (not warnedTotDiv[0]):
					warnedTotDiv[0]=True
					print("\n WARNING!!!!!!!!!!!\n\n Found sample with divergence "+str(float(diffNum)/comparisons)+" from the reference ; at high divergence, MAPLE will struggle both in terms of accuracy and computational demand, so it might make sense to use a traditional phylogenetic approach instead.\n\n End of warning.\n")
					print(diffs)
					print(diffIndex)

	from operator import itemgetter
	print("Now doing sorting")
	sampleDistances.sort(reverse=True,key=itemgetter(0))
	return sampleDistances



#function to check is one sequence is less informative than another;
#returns 0 when the 2 sequences are not comparable, otherwise returns 1 if the first is more informative or if they are identical, and 2 otherwise.
#if onlyFindIdentical (the case sequencing errors are considered) returns 1 if sequences are identical, and 0 otherwise.
def isMinorSequence(probVect1,probVect2,onlyFindIdentical=False):
	indexEntry1, indexEntry2, pos = 0, 0, 0
	entry1=probVect1[indexEntry1]
	entry2=probVect2[indexEntry2]
	found1bigger=False
	found2bigger=False
	while True:
		if entry1[0]!=entry2[0]:
			if onlyFindIdentical:
				return 0
			elif entry1[0]==5:
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
				if onlyFindIdentical:
					if entry2[-1][j]!=entry1[-1][j]:
						return 0
				elif entry2[-1][j]>0.1 and entry1[-1][j]<0.1:
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
#TODO (not necessary) this does not implement rate variation or sequence error models yet - the reason being that the initial placement is always made so far under a simple substitution model. This could however be extended in the future.
#TODO (not necessary) including errors and rate variation could be made simpler by using getPartialVec substantially.
def appendProb(probVectP,probVectC,bLen):
	mutMatrix=mutMatrixGlobal
	Lkcost, indexEntry1, indexEntry2, totalFactor, pos = 0.0, 0, 0, 1.0, 0
	entry1=probVectP[indexEntry1]
	entry2=probVectC[indexEntry2]
	end=min(entry1[1],entry2[1])
	contribLength=bLen
	while True:
		if entry2[0]==5 or entry1[0]==5: # case N
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
def findBestParentTopology(node,child,bestLKdiff,removedBLen,strictTopologyStopRules=strictTopologyStopRules,allowedFailsTopology=allowedFailsTopology,thresholdLogLKtopology=thresholdLogLKtopology):
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
			vectUp1=rootVector(child2.probVect,child2.dist)
			nodesToVisit.append((child1,0,vectUp1,child1.dist,True,bestLKdiff,0))
			vectUp2=rootVector(child1.probVect,child1.dist)
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
					midTot=mergeVectorsUpDown(passedPartials,distance/2,t1.probVect,distance/2)
					if not areVectorsDifferent(midTot,t1.probVectTotUp):
						needsUpdating=False
				else:
					midTot=t1.probVectTotUp
				if midTot==None:
					continue
				midProb=appendProbNode(midTot,removedPartials,removedBLen)
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
				child1=t1.children[0]
				otherChild=t1.children[1]
				if needsUpdating:
					vectUpRight=mergeVectorsUpDown(passedPartials,distance,otherChild.probVect,otherChild.dist)
				else:
					vectUpRight=t1.probVectUpRight
				if vectUpRight!=None:
					nodesToVisit.append((child1,0,vectUpRight,child1.dist,needsUpdating,midProb,failedPasses))
				child1=t1.children[1]
				otherChild=t1.children[0]
				if needsUpdating:
					vectUpLeft=mergeVectorsUpDown(passedPartials,distance,otherChild.probVect,otherChild.dist)
				else:
					vectUpLeft=t1.probVectUpLeft
				if vectUpLeft!=None:
					nodesToVisit.append((child1,0,vectUpLeft,child1.dist,needsUpdating,midProb,failedPasses))

		else: #case when crawling up from child to parent
			otherChild=t1.children[2-direction]
			midBottom=None
			if t1.dist and t1.up!=None: #try appending mid-branch
				if needsUpdating:
					midBottom=mergeVectors(otherChild.probVect,otherChild.dist,passedPartials,distance)
					if midBottom==None:
						continue
					if t1==t1.up.children[0]:
						vectUp=t1.up.probVectUpRight
					else:
						vectUp=t1.up.probVectUpLeft
					midTot=mergeVectorsUpDown(vectUp,t1.dist/2,midBottom,t1.dist/2)
					if not hasattr(t1, 'probVectTotUp'):
						print("Node has no probVectTotUp "+str(t1.name)+" with dist "+str(t1.dist)+" calculating new one")
						print(t1.probVect)
						print(vectUp)
						t1.probVectTotUp=mergeVectorsUpDown(vectUp,t1.dist/2,t1.probVect,t1.dist/2)
						print(t1.probVectTotUp)
					if not areVectorsDifferent(midTot,t1.probVectTotUp):
						needsUpdating=False
				else:
					midTot=t1.probVectTotUp
				if midTot==None:
					continue
				midProb=appendProbNode(midTot,removedPartials,removedBLen)
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
						vectUp=mergeVectorsUpDown(vectUpUp,t1.dist,passedPartials,distance)
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
							midBottom=mergeVectors(otherChild.probVect,otherChild.dist,passedPartials,distance)
							if midBottom==None:
								continue
					else:
						midBottom=t1.probVect
					nodesToVisit.append((t1.up,upChild+1,midBottom,t1.dist,needsUpdating,midProb,failedPasses))
				#now consider case of root node
				else:
					if needsUpdating:
						vectUp=rootVector(passedPartials,distance)
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
	bestBranchLengths=originalBLens
	bestScore=bestLKdiff
	compensanteForBranchLengthChange=True
	secondBranchLengthOptimizationRound=False
	if not bestNodes:
		return originalPlacement, originalLK ,originalBLens, [], 1.0
	BLenHaveBeenOptimized=False
	if aBayesPlusOn:
		if networkOutput:
			listOfProbableNodes=[]
		listofLKcosts=[]
		rootAlreadyConsidered=False
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
			bestAppendingLength=estimateBranchLengthWithDerivative(midTot,removedPartials)
			#now optimize appending location
			midLowerVector=mergeVectors(downVect,distance/2,removedPartials,bestAppendingLength)
			bestTopLength=estimateBranchLengthWithDerivative(upVect,midLowerVector)
			midTopVector=mergeVectorsUpDown(upVect,bestTopLength,removedPartials,bestAppendingLength)
			bestBottomLength=estimateBranchLengthWithDerivative(midTopVector,downVect)
			newMidVector=mergeVectorsUpDown(upVect,bestTopLength,downVect,bestBottomLength)
			if secondBranchLengthOptimizationRound: #if wanted, do a second round of branch length optimization
				bestAppendingLength=estimateBranchLengthWithDerivative(newMidVector,removedPartials)
				midLowerVector=mergeVectors(downVect,bestBottomLength,removedPartials,bestAppendingLength)
				bestTopLength=estimateBranchLengthWithDerivative(upVect,midLowerVector)
				midTopVector=mergeVectorsUpDown(upVect,bestTopLength,removedPartials,bestAppendingLength)
				bestBottomLength=estimateBranchLengthWithDerivative(midTopVector,downVect)
				newMidVector=mergeVectorsUpDown(upVect,bestTopLength,downVect,bestBottomLength)
			appendingCost=appendProbNode(newMidVector,removedPartials,bestAppendingLength)
			if debugging:
				print("Within findBestParentTopology(), second pass (branch length optimization), found appendingCost "+str(appendingCost)+" and branch lengths "+str(bestTopLength)+" "+str(bestBottomLength)+" "+str(bestAppendingLength))
			if compensanteForBranchLengthChange: #if wanted, do a more thorough examination of the appending cost, taking into account a possible change in branch length of the branch on which to be appended.
				initialCost=appendProbNode(upVect,downVect,distance)
				newPartialCost=appendProbNode(upVect,downVect,bestBottomLength+bestTopLength)
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
				testNonsensicalBLens=False
				if testNonsensicalBLens:
					if (not bestBottomLength) and (not bestAppendingLength):
						midLowerVector=mergeVectors(downVect,bestBottomLength,removedPartials,bestAppendingLength)
						if midLowerVector==None:
							print("Problem in findBestParentTopology(), midLowerVector is None")
							raise Exception("exit")
					if (not bestBottomLength) and (not bestTopLength):
						newMidVector=mergeVectorsUpDown(upVect,bestTopLength,downVect,bestBottomLength)
						if newMidVector==None:
							print("Problem in findBestParentTopology(), newMidVector is None")
							raise Exception("exit")
					if (not bestAppendingLength) and (not bestTopLength):
						midTopVector=mergeVectorsUpDown(upVect,bestTopLength,removedPartials,bestAppendingLength)
						if midTopVector==None:
							print("Problem in findBestParentTopology(), midTopVector is None")
							raise Exception("exit")
			if aBayesPlusOn:
				#check that the placement locations is effectively different from the original node
				differentNode=True
				topNode=node
				if t1==topNode:
					differentNode=False
				if (not bestBottomLength):
					while (not topNode.dist) and (topNode.up!=None):
						topNode=topNode.up
					if t1==topNode:
						differentNode=False
				if t1==node.children[1-child]:
					differentNode=False
				#check that placement is not redundant
				if (not bestTopLength):
					differentNode=False
				# check if this is a root placement
				if (not rootAlreadyConsidered) and (not bestTopLength):
					topNode=t1.up
					while (not topNode.dist) and (topNode.up!=None):
						topNode=topNode.up
					if topNode.up==None:
						rootAlreadyConsidered=True
						listofLKcosts.append(optimizedScore)
						if networkOutput:
							listOfProbableNodes.append(topNode)
				elif differentNode: #add placement to the list of legit ones
					listofLKcosts.append(optimizedScore)
					if networkOutput :
						listOfProbableNodes.append(t1)
						# if optimizedScore>(bestScore-2.0) and optimizedScore<(bestScore-1.0) and node.children[child].name=="EPI_ISL_914595":
						# 	print("low-prob Alternative placement")
						# 	print(node.children[child].name)
						# 	print(t1.name)
						# 	print(optimizedScore)
						# 	print(bestScore)
						# 	print(bestTopLength)
						# 	print(bestBottomLength)
						# 	print(bestAppendingLength)
						# 	print(removedPartials)
						# 	print("")
						# 	print(newMidVector)
						# 	print("")
						# 	if child==0:
						# 		print(node.probVectUpRight)
						# 	else:
						# 		print(node.probVectUpLeft)


	if not BLenHaveBeenOptimized:
		bestBranchLengths=(bestNode.dist/2,bestNode.dist/2,removedBLen)

	if aBayesPlusOn:
		# calculate support(s) and possibly add nodes to the list of alternative placements
		support=math.exp(originalLK)
		totSupport=support
		for i in range(len(listofLKcosts)):
			listofLKcosts[i]=math.exp(listofLKcosts[i])
			totSupport+=listofLKcosts[i]
		support=support/totSupport
		finalListOfNodes=[]
		if networkOutput:
			for i in range(len(listofLKcosts)):
				listofLKcosts[i]=listofLKcosts[i]/totSupport
			for i in range(len(listofLKcosts)):
				if listofLKcosts[i]>=minBranchSupport:
					finalListOfNodes.append((listOfProbableNodes[i],listofLKcosts[i]))
		return bestNode, bestScore, bestBranchLengths, finalListOfNodes, support

	else:
		return bestNode, bestScore, bestBranchLengths, [], 1.0












#function to find the best node in the tree where to append the new sample; traverses the tree and tries to append the sample at each node and mid-branch nodes, 
# but stops traversing when certain criteria are met.
#TODO (not necessary) not included rate variation or sequence errors yet since the initial placement so far is only done under a simple model.
def findBestParentForNewSample(root,diffs,sample):
	bestNodes=[]
	bestNode=root
	bestBranchLengths=(False,False,oneMutBLen)
	if not root.children: #check if the new leaf is strictly less informative than already placed leaf
		if errorRateSiteSpecificFile or errorRateFixed or estimateErrorRate or estimateSiteSpecificErrorRate or supportFor0Branches:
			comparison=isMinorSequence(root.probVect,diffs,onlyFindIdentical=True)
		else:
			comparison=isMinorSequence(root.probVect,diffs)
		if comparison==1:
			root.minorSequences.append(sample)
			return root, 1.0, None
		elif comparison==2:
			totalMissedMinors[0]+=1
	rootVect=rootVector(root.probVect,False)
	bestLKdiff=appendProb(rootVect,diffs,oneMutBLen)
	nodesToVisit=[]
	for child in root.children:
		nodesToVisit.append((child,bestLKdiff,0))
	while nodesToVisit:
		t1,parentLK,failedPasses=nodesToVisit.pop()
		if not t1.children: #check if the new leaf is strictly less informative than already placed leaf
			if errorRateSiteSpecificFile or errorRateFixed or estimateErrorRate or estimateSiteSpecificErrorRate or supportFor0Branches:
				comparison=isMinorSequence(t1.probVect,diffs,onlyFindIdentical=True)
			else:
				comparison=isMinorSequence(t1.probVect,diffs)
			if comparison==1:
				t1.minorSequences.append(sample)
				return t1, 1.0, None
			elif comparison==2:
				totalMissedMinors[0]+=1

		if t1.dist and t1.up!=None: # try first placing as a descendant of the mid-branch point of the branch above the current node.
			LKdiff=appendProb(t1.probVectTotUp,diffs,oneMutBLen)
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
			bestAppendingLength=estimateBranchLengthWithDerivative(node.probVectTotUp,diffs)
			#now optimize appending location
			midLowerVector=mergeVectors(node.probVect,node.dist/2,diffs,bestAppendingLength)
			bestTopLength=estimateBranchLengthWithDerivative(upVect,midLowerVector)
			midTopVector=mergeVectorsUpDown(upVect,bestTopLength,diffs,bestAppendingLength)
			bestBottomLength=estimateBranchLengthWithDerivative(midTopVector,node.probVect)
			newMidVector=mergeVectorsUpDown(upVect,bestTopLength,node.probVect,bestBottomLength)
			if secondBranchLengthOptimizationRound: #if wanted, do a second round of branch length optimization
				bestAppendingLength=estimateBranchLengthWithDerivative(newMidVector,diffs)
				midLowerVector=mergeVectors(node.probVect,bestBottomLength,diffs,bestAppendingLength)
				bestTopLength=estimateBranchLengthWithDerivative(upVect,midLowerVector)
				midTopVector=mergeVectorsUpDown(upVect,bestTopLength,diffs,bestAppendingLength)
				bestBottomLength=estimateBranchLengthWithDerivative(midTopVector,node.probVect)
				newMidVector=mergeVectorsUpDown(upVect,bestTopLength,node.probVect,bestBottomLength)
			appendingCost=appendProb(newMidVector,diffs,bestAppendingLength)
			if compensanteForBranchLengthChange: #if wanted, do a more thorough examination of the appending cost, taking into account a possible change in branch length of the branch on which to be appended.
				initialCost=appendProbNode(upVect,node.probVect,node.dist)
				newPartialCost=appendProbNode(upVect,node.probVect,bestBottomLength+bestTopLength)
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
#TODO (not necessary) not included rate variation or sequence errors yet since the initial placement so far is only done under a simple model.
def placeSampleOnTree(node,newPartials,sample,newChildLK, bestUpLength, bestDownLength, bestAppendingLength,pseudoMutCounts):
	tryNewRoot=False
	if newChildLK<-0.01:
		sumChildLKs[0]+=newChildLK
		numChildLKs[0]+=1
	if node.up==None:
		tryNewRoot=True
		totRoot=rootVector(node.probVect,False)
		bestAppendingLength=estimateBranchLengthWithDerivative(totRoot,newPartials)
		root=node
		newChildLK=appendProb(totRoot,newPartials,bestAppendingLength)
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
				if (not bestDownLength) or (bestDownLength>1.01*node.dist) or (bestDownLength<0.99*node.dist):
					node.dist=bestDownLength
					nodeList=[(node,2),(node.up,child)]
					updatePartials(nodeList)
	#in case of best placement as a descendant appended exactly at the root node, attempt also to create new root
	if tryNewRoot:
		node=root
		probOldRoot = findProbRoot(node.probVect)
		rootUpLeft=rootVector(node.probVect,bestAppendingLength/2)
		bestRightLength=estimateBranchLengthWithDerivative(rootUpLeft,newPartials)
		rootUpRight=rootVector(newPartials,bestRightLength)
		bestLeftLength=estimateBranchLengthWithDerivative(rootUpRight,node.probVect)
		secondBranchLengthOptimizationRound=True
		if secondBranchLengthOptimizationRound: #if wanted, do a second round of branch length optimization
			rootUpLeft=rootVector(node.probVect,bestLeftLength)
			bestRightLength=estimateBranchLengthWithDerivative(rootUpLeft,newPartials)
			rootUpRight=rootVector(newPartials,bestRightLength)
			bestLeftLength=estimateBranchLengthWithDerivative(rootUpRight,node.probVect)
		probVectRoot,probRoot = mergeVectors(node.probVect,bestLeftLength,newPartials,bestRightLength,returnLK=True)
		probRoot+= findProbRoot(probVectRoot)
		parentLKdiff=probRoot-probOldRoot
		if parentLKdiff<=newChildLK: #best is just placing as descendant of the root
			bestRightLength=bestAppendingLength
			bestLeftLength=False
			probVectRoot=mergeVectors(node.probVect,bestLeftLength,newPartials,bestRightLength)
			rootUpRight=rootVector(newPartials,bestRightLength)
		#now add new root to the tree
		newRoot=Tree()
		if probVectRoot==None:
			print("Issue with new root probVect inside placeSampleOnTree()")
			raise Exception("exit")
		newRoot.probVect=probVectRoot
		newRoot.probVectUpRight=rootUpRight
		newRoot.probVectUpLeft=rootVector(node.probVect,bestLeftLength)
		node.up=newRoot
		node.dist=bestLeftLength
		newRoot.add_child(node)
		newNode=Tree(name=sample,dist=bestRightLength)
		if bestRightLength>0.01 and (not warnedBLen[0]):
			warnedBLen[0]=True
			print("\n WARNING!!!!!!!!!!!\n\n Found branch of length "+str(bestRightLength)+" ; at high divergence, MAPLE will struggle both in terms of accuracy and computational demand, so it might make sense to use a traditional phylogenetic approach.\n\n End of warning.\n")
		newNode.minorSequences=[]
		newNode.up=newRoot
		newRoot.add_child(newNode)
		newNode.probVect=newPartials
		if bestRightLength:
			newNode.probVectTotUp=mergeVectorsUpDown(newRoot.probVectUpLeft,bestRightLength/2,newPartials,bestRightLength/2)
		nodeList=[(node,2)]
		updatePartials(nodeList)
		return newRoot

	#print("adding internal node")
	#in all other cases (not attempting to add a new root) create a new internal node in the tree and add sample as a descendant.
	newInternalNode=Tree()
	node.up.children[child]=newInternalNode
	newInternalNode.up=node.up
	newInternalNode.add_child(node)
	node.up=newInternalNode
	node.dist=bestDownLength
	newNode=Tree(name=sample,dist=bestAppendingLength)
	if bestAppendingLength>0.01 and (not warnedBLen[0]):
		warnedBLen[0]=True
		print("\n WARNING!!!!!!!!!!!\n\n Found branch of length "+str(bestAppendingLength)+" ; at high divergence, MAPLE will struggle both in terms of accuracy and computational demand, so it might make sense to use a more traditional phylogenetic approach (e.g. FastTree, IQtree or RAxML).\n\n End of warning.\n")
	newNode.minorSequences=[]
	newNode.up=newInternalNode
	newInternalNode.add_child(newNode)
	newInternalNode.dist=bestUpLength
	newNode.probVect=newPartials
	newInternalNode.probVect=mergeVectors(node.probVect,bestDownLength,newPartials,bestAppendingLength)
	newInternalNode.probVectUpRight=mergeVectorsUpDown(vectUp,bestUpLength,newPartials,bestAppendingLength)
	newInternalNode.probVectUpLeft=mergeVectorsUpDown(vectUp,bestUpLength,node.probVect,bestDownLength)
	if newInternalNode.probVect==None:
		print("Problem in placeSampleOnTree(), probVect is None")
		raise Exception("exit")
	if newInternalNode.probVectUpRight==None:
		print("Problem in placeSampleOnTree(), probVectUpRight is None")
		raise Exception("exit")
	if newInternalNode.probVectUpLeft==None:
		print("Problem in placeSampleOnTree(), probVectUpLeft is None")
		raise Exception("exit")
	if bestUpLength:
		newInternalNode.probVectTotUp=mergeVectorsUpDown(vectUp,bestUpLength/2,newInternalNode.probVect,bestUpLength/2)
	else:
		newInternalNode.probVectTotUp=None
	if bestAppendingLength:
		newNode.probVectTotUp=mergeVectorsUpDown(newInternalNode.probVectUpLeft,bestAppendingLength/2,newPartials,bestAppendingLength/2)
		updatePesudoCounts(newInternalNode.probVectUpLeft,newPartials,pseudoMutCounts)
	else:
		newNode.probVectTotUp=None
	if not bestDownLength:
		node.probVectTotUp=None
	nodeList=[(node,2),(newInternalNode.up,child)]
	updatePartials(nodeList)

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
def appendProbNode(probVectP,probVectC,bLen):
	if not errorRateSiteSpecific: 
		errorRate=errorRateGlobal
	if not useRateVariation:
		mutMatrix=mutMatrixGlobal
	Lkcost, indexEntry1, indexEntry2, totalFactor, pos = 0.0, 0, 0, 1.0, 0
	entry1=probVectP[indexEntry1] #parent
	entry2=probVectC[indexEntry2] #child
	end=min(entry1[1],entry2[1])
	contribLength=bLen
	while True:
		if entry2[0]==5 or entry1[0]==5: # case entry1 is N
			#if parent node is type "N", in theory we might have to calculate the contribution of root nucleotides; 
			# however, if this node is "N" then every other node in the current tree is "N", so we can ignore this since this contribution cancels out in relative terms.
			pos=min(entry1[1],entry2[1])
			pass
		else:
			#contribLength will be here the total length from the root or from the upper node, down to the down node.
			contribLength=bLen
			if entry1[0]<5:
				if len(entry1)==3+usingErrorRate:
					contribLength+=entry1[2]
				elif len(entry1)==4+usingErrorRate:
					contribLength+=entry1[3]
			else:
				if len(entry1)==4:
					contribLength+=entry1[2]
			if entry2[0]<5:
				if len(entry2)==3+usingErrorRate:
					contribLength+=entry2[2]
			else:
				if len(entry2)==4:
					contribLength+=entry2[2]

			if entry1[0]==4: # case entry1 is R	
				flag1 = (usingErrorRate and (len(entry1)>2) and entry1[-1] )
				if entry2[0]==4:
					flag2 = (usingErrorRate and (len(entry2)>2) and entry2[-1] )
					if len(entry1)==4+usingErrorRate:
						contribLength+=entry1[2]
					end=min(entry1[1],entry2[1])
					if contribLength:
						Lkcost+=contribLength*(cumulativeRate[end]-cumulativeRate[pos])
					if flag1 or flag2:
						if errorRateSiteSpecific: cumErrorRate = cumulativeErrorRate[pos] - cumulativeErrorRate[end] #both this and the next one should end up being negative
						else: cumErrorRate = errorRate * (pos-end) 
						Lkcost += cumErrorRate * (flag1 + flag2)
					pos=end

				#entry1 is reference and entry2 is of type "O"
				elif entry2[0]==6:
					if useRateVariation:
						mutMatrix=mutMatrices[pos]
					i1=refIndeces[pos]
					if len(entry1)==4+usingErrorRate:
						tot=0.0
						if errorRateSiteSpecific: errorRate = errorRates[pos]
						tot3=getPartialVec(6, contribLength, mutMatrix, None, vect=entry2[-1])
						tot2=getPartialVec(i1, entry1[2], mutMatrix, errorRate, flag=flag1)
						for i in range4:
							tot+=tot3[i]*tot2[i]*rootFreqs[i]
						tot/=rootFreqs[i1]
					else:
						if contribLength:
							tot3=getPartialVec(6, contribLength, mutMatrix, None, vect=entry2[-1])
							tot=tot3[i1]
						else:
							tot=entry2[-1][i1]
					totalFactor*=tot
					pos+=1

				else: #entry1 is R and entry2 is a different but single nucleotide
					flag2 = (usingErrorRate and (len(entry2)>2) and entry2[-1] )
					if useRateVariation:
						mutMatrix=mutMatrices[pos]
					if len(entry1)==4+usingErrorRate:
						i1=refIndeces[pos]
						i2=entry2[0]
						if errorRateSiteSpecific: errorRate = errorRates[pos]
						tot3=getPartialVec(i2, contribLength, mutMatrix, errorRate, flag=flag2)
						tot2=getPartialVec(i1, entry1[2], mutMatrix, errorRate, flag=flag1)
						tot=0.0
						for i in range4:
							tot+=tot3[i]*tot2[i]*rootFreqs[i]
						totalFactor*=tot/rootFreqs[i1]
					else:
						if flag2:
							if errorRateSiteSpecific: errorRate = errorRates[pos]
							totalFactor*=min(0.25,mutMatrix[refIndeces[pos]][entry2[0]]*contribLength)+errorRate*0.33333
						else:
							if contribLength:
								totalFactor*=min(0.25,mutMatrix[refIndeces[pos]][entry2[0]]*contribLength)
							else:
								return float("-inf")
					pos+=1

			# entry1 is of type "O"
			elif entry1[0]==6:
				if useRateVariation:
					mutMatrix=mutMatrices[pos]
				if entry2[0]==6:
					tot=0.0
					if contribLength:
						tot3=getPartialVec(6, contribLength, mutMatrix, None, vect=entry2[-1])
						for j in range4:
							tot+=entry1[-1][j]*tot3[j]
					else:
						for j in range4:
							tot+=entry1[-1][j]*entry2[-1][j]
					totalFactor*=tot
				else: #entry1 is "O" and entry2 is a nucleotide
					if entry2[0]==4:
						i2=refIndeces[pos]
					else:
						i2=entry2[0]
					if (usingErrorRate and (len(entry2)>2) and entry2[-1] ):
						if errorRateSiteSpecific: errorRate = errorRates[pos]
						tot3=getPartialVec(i2, contribLength, mutMatrix, errorRate, flag=True)
					else:
						tot3=getPartialVec(i2, contribLength, mutMatrix, None, flag=False)
					tot=0.0
					for j in range4:
						tot+=entry1[-1][j]*tot3[j]
					totalFactor*=tot
				pos+=1

			else: #entry1 is a non-ref nuc
				flag1 = (usingErrorRate and (len(entry1)>2) and entry1[-1] )
				if useRateVariation:
					mutMatrix=mutMatrices[pos]
				if entry2[0]==entry1[0]:
					flag2 = (usingErrorRate and (len(entry2)>2) and entry2[-1] )
					if len(entry1)==4+usingErrorRate:
						contribLength+=entry1[2]
					if contribLength:
						Lkcost+=mutMatrix[entry1[0]][entry1[0]]*contribLength
					if flag1 or flag2:
						if errorRateSiteSpecific: errorRate = errorRates[pos]
						Lkcost -= errorRate * (flag1 + flag2)
				else: #entry1 is a nucleotide and entry2 is not the same as entry1
					i1=entry1[0]
					if entry2[0]<5: #entry2 is a nucleotide
						if entry2[0]==4:
							i2=refIndeces[pos]
						else:
							i2=entry2[0]
						flag2 = (usingErrorRate and (len(entry2)>2) and entry2[-1] )
						if len(entry1)==4+usingErrorRate:
							if errorRateSiteSpecific: errorRate = errorRates[pos]
							tot3=getPartialVec(i2, contribLength, mutMatrix, errorRate, flag=flag2)
							tot2=getPartialVec(i1, entry1[2], mutMatrix, errorRate, flag=flag1)
							tot=0.0
							for j in range4:
								tot+=rootFreqs[j]*tot3[j]*tot2[j]
							totalFactor*=tot/rootFreqs[i1]
						else:
							if flag1 or flag2:
								if errorRateSiteSpecific: errorRate = errorRates[pos]
								totalFactor*=(min(0.25,mutMatrix[i1][i2]*contribLength)+(flag1 + flag2)*0.33333*errorRate)
							else:
								if contribLength:
									totalFactor*=min(0.25,mutMatrix[i1][i2]*contribLength)
								else:
									return float("-inf")
					else: #entry1 is a nucleotide and entry2 is of type "O"
						if errorRateSiteSpecific: errorRate = errorRates[pos]
						if len(entry1)==4+usingErrorRate:
							tot2=getPartialVec(i1, entry1[2], mutMatrix, errorRate, flag=flag1)
							tot3=getPartialVec(6, contribLength, mutMatrix, errorRate, vect=entry2[-1])
							tot=0.0
							for i in range4:
								tot+=tot2[i]*tot3[i]*rootFreqs[i]
							totalFactor*=(tot/rootFreqs[i1])
						else:
							if contribLength:
								tot3=getPartialVec(6, contribLength, mutMatrix, None, vect=entry2[-1])
								totalFactor*=tot3[i1]
							else:
								totalFactor*=entry2[-1][i1]
				pos+=1

		if totalFactor<=minimumCarryOver:
			try:
				if totalFactor<sys.float_info.min:
					if bLen:
						print("In appendProbNode() strangely small probability")
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
def traverseTreeToOptimizeBranchLengths(root,testing=False,fastPass=False):
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
			bestLength=estimateBranchLengthWithDerivative(upVect,node.probVect)
			if bestLength or node.dist:
				if testing:
					currentCost=appendProbNode(upVect,node.probVect,node.dist)
					newCost=appendProbNode(upVect,node.probVect,bestLength)
					totalLKimprovementBL+=newCost-currentCost
				if (not bestLength) or (not node.dist) or node.dist/bestLength>1.01 or node.dist/bestLength<0.99:
					node.dist=bestLength
					updates+=1
					if not fastPass:
						nodeList=[(node,2),(node.up,child)]
						updatePartials(nodeList)
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
def placeSubtreeOnTree(node,newPartials,appendedNode,newChildLK,bestBranchLengths): #,parentNode,parentNodeReplacements=1
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
			if (not bestDownLength) or (bestDownLength>1.01*node.dist) or (bestDownLength<0.99*node.dist):
				node.dist=bestDownLength
				nodeList=[(node,2),(node.up,child)]
				updatePartials(nodeList)
			
	#in case of best placement as a descendant appended exactly at the root node, attempt also to create new root
	if tryNewRoot:
		node=root
		probOldRoot = findProbRoot(node.probVect)
		rootUpLeft=rootVector(node.probVect,bestAppendingLength/2)
		bestRightLength=estimateBranchLengthWithDerivative(rootUpLeft,newPartials)
		rootUpRight=rootVector(newPartials,bestRightLength)
		bestLeftLength=estimateBranchLengthWithDerivative(rootUpRight,node.probVect)
		secondBranchLengthOptimizationRound=True
		if secondBranchLengthOptimizationRound: #if wanted, do a second round of branch length optimization
			rootUpLeft=rootVector(node.probVect,bestLeftLength)
			bestRightLength=estimateBranchLengthWithDerivative(rootUpLeft,newPartials)
			rootUpRight=rootVector(newPartials,bestRightLength)
			bestLeftLength=estimateBranchLengthWithDerivative(rootUpRight,node.probVect)
		probVectRoot,probRoot = mergeVectors(node.probVect,bestLeftLength,newPartials,bestRightLength,returnLK=True)
		probRoot+= findProbRoot(probVectRoot)
		parentLKdiff=probRoot-probOldRoot
		if parentLKdiff<=newChildLK: #best is just placing as descendant of the root
			bestRightLength=bestAppendingLength
			bestLeftLength=False
			probVectRoot=mergeVectors(node.probVect,bestLeftLength,newPartials,bestRightLength)
			rootUpRight=rootVector(newPartials,bestRightLength)
		#now add new root to the tree
		#newRoot=Tree()
		newRoot=appendedNode.up
		newRoot.up=None
		newRoot.dirty=True
		newRoot.dist=defaultBLen
		#newRoot.replacements=parentNodeReplacements+1
		newRoot.replacements+=1
		if probVectRoot==None:
			print("Issue with new root probVect inside placeSubtreeOnTree()")
			raise Exception("exit")
		newRoot.probVect=probVectRoot
		newRoot.probVectUpRight=rootUpRight
		newRoot.probVectUpLeft=rootVector(node.probVect,bestLeftLength)
		node.up=newRoot
		node.dist=bestLeftLength
		newRoot.children[0]=node
		#newRoot.add_child(node)
		#appendedNode.up=newRoot
		#newRoot.add_child(appendedNode)
		newRoot.children[1]=appendedNode
		appendedNode.dist=bestRightLength
		appendedNode.replacements+=1
		nodeList=[(node,2),(appendedNode,2)]
		updatePartials(nodeList)
		return newRoot

	#print("adding internal node")
	#in all other cases (not attempting to add a new root) create a new internal node in the tree and add sample as a descendant.
	#newInternalNode=Tree()
	newInternalNode=appendedNode.up
	newInternalNode.dirty=True
	#newInternalNode.replacements=parentNodeReplacements+1
	newInternalNode.replacements+=1
	node.up.children[child]=newInternalNode
	newInternalNode.up=node.up
	newInternalNode.children[0]=node
	#newInternalNode.add_child(node)
	node.up=newInternalNode
	#appendedNode.up=newInternalNode
	appendedNode.replacements+=1
	#newInternalNode.add_child(appendedNode)
	newInternalNode.children[1]=appendedNode
	newInternalNode.probVect=mergeVectors(node.probVect,bestDownLength,newPartials,bestAppendingLength)
	if newInternalNode.probVect==None:
		newInternalNode.probVectUpLeft=mergeVectorsUpDown(vectUp,bestUpLength,node.probVect,bestDownLength)
		bestAppendingLength=estimateBranchLengthWithDerivative(newInternalNode.probVectUpLeft,newPartials)
		newInternalNode.probVectUpRight=mergeVectorsUpDown(vectUp,bestUpLength,newPartials,bestAppendingLength)
		bestDownLength=estimateBranchLengthWithDerivative(newInternalNode.probVectUpRight,node.probVect)
		newInternalNode.probVect=mergeVectors(node.probVect,bestDownLength,newPartials,bestAppendingLength)
		if newInternalNode.probVect==None:
			print("newInternalNode.probVect is None after updating the optimal branch lengths")
			bestAppendingLength=oneMutBLen/5
			bestDownLength=oneMutBLen/5
			newInternalNode.probVect=mergeVectors(node.probVect,bestDownLength,newPartials,bestAppendingLength)
	newInternalNode.probVectUpRight=mergeVectorsUpDown(vectUp,bestUpLength,newPartials,bestAppendingLength)
	if newInternalNode.probVectUpRight==None:
		bestUpLength=estimateBranchLengthWithDerivative(vectUp,newInternalNode.probVect)
		newInternalNode.probVectUpLeft=mergeVectorsUpDown(vectUp,bestUpLength,node.probVect,bestDownLength)
		bestAppendingLength=estimateBranchLengthWithDerivative(newInternalNode.probVectUpLeft,newPartials)
		newInternalNode.probVectUpRight=mergeVectorsUpDown(vectUp,bestUpLength,newPartials,bestAppendingLength)
		if newInternalNode.probVectUpRight==None:
			print("newInternalNode.probVectUpRight is None after updating the optimal branch lengths")
			bestUpLength=oneMutBLen/5
			bestAppendingLength=oneMutBLen/5
			newInternalNode.probVectUpRight=mergeVectorsUpDown(vectUp,bestUpLength,newPartials,bestAppendingLength)
		newInternalNode.probVect=mergeVectors(node.probVect,bestDownLength,newPartials,bestAppendingLength)
	newInternalNode.probVectUpLeft=mergeVectorsUpDown(vectUp,bestUpLength,node.probVect,bestDownLength)
	if newInternalNode.probVectUpLeft==None:
		bestUpLength=estimateBranchLengthWithDerivative(vectUp,newInternalNode.probVect)
		bestDownLength=estimateBranchLengthWithDerivative(newInternalNode.probVectUpRight,node.probVect)
		newInternalNode.probVectUpLeft=mergeVectorsUpDown(vectUp,bestUpLength,node.probVect,bestDownLength)
		if newInternalNode.probVectUpLeft==None:
			print("newInternalNode.probVectUpRight is None after updating the optimal branch lengths")
			bestUpLength=oneMutBLen/5
			bestDownLength=oneMutBLen/5
			newInternalNode.probVectUpLeft=mergeVectorsUpDown(vectUp,bestUpLength,node.probVect,bestDownLength)
		newInternalNode.probVect=mergeVectors(node.probVect,bestDownLength,newPartials,bestAppendingLength)
		newInternalNode.probVectUpRight=mergeVectorsUpDown(vectUp,bestUpLength,newPartials,bestAppendingLength)
	appendedNode.dist=bestAppendingLength
	newInternalNode.dist=bestUpLength
	node.dist=bestDownLength
	if not bestAppendingLength:
		appendedNode.probVectTotUp=None
	if bestUpLength:
		newInternalNode.probVectTotUp=mergeVectorsUpDown(vectUp,bestUpLength/2,newInternalNode.probVect,bestUpLength/2)
	if not bestDownLength:
		node.probVectTotUp=None
	nodeList=[(node,2),(newInternalNode.up,child),(appendedNode,2)]
	updatePartials(nodeList)

	return None





#remove node from the current position in the tree and re-attach it at a new given place new bestNode.
# First remove node from the tree, then update the genome lists;
# then find the exact best reattachment of node and update the genome lists again using function placeSubtreeOnTree().
def cutAndPasteNode(node,bestNode,bestBranchLengths,bestLK):
	#remove node from the tree
	if debugging:
		print("In cutAndPasteNode() removing subtree from the tree, subtree root partials: ")
	parentNode=node.up
	#parentNodeReplacements=parentNode.replacements
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
		if sibling.children:
			sibling.probVectUpRight=rootVector(sibling.children[1].probVect,sibling.children[1].dist)
			sibling.probVectUpLeft=rootVector(sibling.children[0].probVect,sibling.children[0].dist)
			nodeList=[(sibling.children[0],2),(sibling.children[1],2)]
			updatePartials(nodeList)
	else:
		nodeList=[(sibling,2),(sibling.up,childP)]
		updatePartials(nodeList)
	#re-place the node and re-update the vector lists
	newRoot = placeSubtreeOnTree(bestNode,node.probVect,node,bestLK,bestBranchLengths)#,parentNode,parentNodeReplacements=parentNodeReplacements
	
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
		totalLK=calculateTreeLikelihood(currentRoot)
		intermediateLKsFile.write("Topology "+str(topologyChanges[0])+", LK: "+str(totalLK)+"\n")

	#if the root of the tree has changed, return the new root
	if sibling.up==None:
		if newRoot!=None:
			return newRoot
		return sibling
	return newRoot




# try to find a re-placement of a dirty node of the tree to improve the topology.
# Pretend to cut out the subtree at this node (so to keep track of the effect of subtree removal on the likelihoods), and look for somewhere else in the tree where to attach it (an SPR move).
# To find the best location of the new re-attachment, we traverse the tree starting at the current attachment node, and for each node visited we evaluate the reattachment,
# as done by the findBestParentTopology() function.
# To avoid traversing the whole tree for each SPR move, we use stopping conditions similar to those used in the initial sample placement process.
# After we find the best SPR move for the given node, we execute it using the cutAndPasteNode() function.
def traverseTreeForTopologyUpdate(node,strictTopologyStopRules=strictTopologyStopRules,allowedFailsTopology=allowedFailsTopology,thresholdLogLKtopology=thresholdLogLKtopology,thresholdTopologyPlacement=thresholdTopologyPlacement):
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
		originalLK=appendProbNode(vectUp,node.probVect,bestCurrenBLen)
		bestCurrentLK=originalLK
		if (bestCurrentLK<thresholdTopologyPlacement) or (supportFor0Branches and aBayesPlusOn):
			bestCurrenBLen=estimateBranchLengthWithDerivative(vectUp,node.probVect)
			if bestCurrenBLen or node.dist:
				if (not bestCurrenBLen) or (not node.dist) or node.dist/bestCurrenBLen>1.01 or node.dist/bestCurrenBLen<0.99:
					bLenChanged=True
				bestCurrentLK=appendProbNode(vectUp,node.probVect,bestCurrenBLen)
				if bestCurrentLK<originalLK:
					#print("Inside traverseTreeForTopologyUpdate(), BL optimization brought to a decrease in LK.")
					bestCurrenBLen=node.dist
					bestCurrentLK=originalLK
					bLenChanged=False
				if debugging:
					print("New LK after modifing bLen: "+str(bestCurrentLK)+" from length "+str(bestCurrenBLen))
					print(mutMatrixGlobal)
				if bestCurrentLK==float("-inf"):
					print("Found an infinite cost of bestCurrentLK "+str(bestCurrentLK)+" using appendProbNode()")
					raise Exception("exit")

		topologyUpdated=False
		if (bestCurrentLK<thresholdTopologyPlacement) or (supportFor0Branches and aBayesPlusOn):
			#now find the best place on the tree where to re-attach the subtree rooted at "node"
			#but to do that we need to consider new vector probabilities after removing the node that we want to replace
			# this is done using findBestParentTopology().
			if debugging:
				if bLenChanged:
					print("Better branch length found at node "+str(node)+", from "+str(node.dist)+" (LK "+str(originalLK)+") to "+str(bestCurrenBLen)+" (LK "+str(bestCurrentLK)+"). Children")
					print(node.children)

			bestNodeSoFar , bestLKdiff , bestBranchLengths, listOfBestPlacements, branchSupport = findBestParentTopology(parentNode,child,bestCurrentLK,bestCurrenBLen,strictTopologyStopRules=strictTopologyStopRules,allowedFailsTopology=allowedFailsTopology,thresholdLogLKtopology=thresholdLogLKtopology)
			if bestLKdiff==float("inf"):
				print("Found an infinite improvement of "+str(bestLKdiff)+" using findBestParentTopology()")
				raise Exception("exit")
			if bestLKdiff<-1e50:
				print("Error: found likelihood cost is very heavy, this might mean that the reference used is not the same used to generate the input diff file")
				raise Exception("exit")
			if (bestLKdiff+thresholdTopologyPlacement>bestCurrentLK) and (not doNotImproveTopology):
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
				if bestNodeSoFar.up==sibling and (not bestBranchLengths[0]):
					topologyUpdated=False

				if topologyUpdated:
					totalImprovement=(bestLKdiff-originalLK)
					if originalLK==float("-inf"):
						totalImprovement=(bestLKdiff-bestCurrentLK)
					if totalImprovement==float("inf"):
						print("Found an infinite topology improvement of "+str(totalImprovement)+" from originalLK "+str(originalLK)+", bestCurrentLK "+str(bestCurrentLK)+" and bestLKdiff "+str(bestLKdiff))
						raise Exception("exit")
					if debugging:
						print("Performing SPR move, detaching node "+str(node)+" with children")
						print(node.children)
						print("and reattaching it around node "+str(bestNodeSoFar)+" with children")
						print(bestNodeSoFar.children)
						print("And branch lengths")
						print(bestBranchLengths)
						print("for supposed improvement "+str(totalImprovement))
					newRoot = cutAndPasteNode(node,bestNodeSoFar,bestBranchLengths,bestLKdiff)
					bLenChanged=False
			if (not topologyUpdated) and aBayesPlusOn:
				if networkOutput:
					node.alternativePlacements=listOfBestPlacements
				node.support=branchSupport
					
		if (not topologyUpdated) and bLenChanged:
			if debugging:
				print("Changing branch length (pos 3) from "+str(node.dist)+" to "+str(bestCurrenBLen)+" at node "+str(node)+" with children ")
				print(node.children)
			node.dist=bestCurrenBLen
			nodeList=[(node,2),(node.up,child)]
			updatePartials(nodeList)
			totalImprovement=(bestCurrentLK-originalLK)
			if originalLK==float("-inf"):
				totalImprovement=0
			if totalImprovement==float("inf"):
				print("Found an infinite branch length improvement of "+str(totalImprovement)+" from originalLK "+str(originalLK)+" and bestCurrentLK "+str(bestCurrentLK))
				raise Exception("exit")

	return newRoot,totalImprovement









#traverse the tree (here the input "node" will usually be the root), and for each dirty node ancountered, call traverseTreeForTopologyUpdate() 
# to attempt an SPR move by cutting the subtree rooted at this dirty node and trying to re-append it elsewhere.
def startTopologyUpdates(node,checkEachSPR=False,strictTopologyStopRules=strictTopologyStopRules,allowedFailsTopology=allowedFailsTopology,thresholdLogLKtopology=thresholdLogLKtopology,thresholdTopologyPlacement=thresholdTopologyPlacement,printEvery=1000):
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
				oldTreeLK=calculateTreeLikelihood(root,checkCorrectness=True)
				reCalculateAllGenomeLists(root, checkExistingAreCorrect=True)
			if aBayesPlusOn and networkOutput:
				newNode.alternativePlacements=[]
			newRoot2,improvement=traverseTreeForTopologyUpdate(newNode,strictTopologyStopRules=strictTopologyStopRules,allowedFailsTopology=allowedFailsTopology,thresholdLogLKtopology=thresholdLogLKtopology,thresholdTopologyPlacement=thresholdTopologyPlacement)
			if checkEachSPR:
				root=newNode
				while root.up!=None:
					root=root.up
				newTreeLK=calculateTreeLikelihood(root)
				reCalculateAllGenomeLists(root, checkExistingAreCorrect=True)
				print("In startTopologyUpdates, LK score of improvement "+str(newTreeLK)+" - "+str(oldTreeLK)+" = "+str(newTreeLK-oldTreeLK)+", was supposed to be "+str(improvement))
				if newTreeLK-oldTreeLK < improvement-1.0:
					print("In startTopologyUpdates, LK score of improvement "+str(newTreeLK)+" - "+str(oldTreeLK)+" = "+str(newTreeLK-oldTreeLK)+" is less than what is supposed to be "+str(improvement))
					raise Exception("exit")
			totalImprovement+=improvement
			if newRoot2!=None:
				newRoot=newRoot2
			numNodes+=1
			#if (numNodes%1000)==0:
			if (numNodes%printEvery)==0:
				print("Processed topology for "+str(numNodes)+" nodes.")
	return newRoot,totalImprovement





#Given a tree, and a final substitution rate matrix, calculate the likelihood of the tree
def calculateTreeLikelihood(root,checkCorrectness=False):
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
				newLower, Lkcontribution=mergeVectors(node.children[0].probVect,node.children[0].dist,node.children[1].probVect,node.children[1].dist,returnLK=True)
				totalLK+=Lkcontribution
				if newLower==None:
					print("Strange, inconsistent lower genome list creation in calculateTreeLikelihood(); old list, and children lists")
					raise Exception("exit")
				elif numTopologyImprovements and checkCorrectness and areVectorsDifferentDebugging(node.probVect,newLower):
					print("Strange, while calculating tree likelihood encountered non-updated lower likelihood genome list at node "+str(node)+" with children ")
					raise Exception("exit")
				lastNode=node
				node=node.up
				direction=1
	#now add contribution from the root
	totalLK+=findProbRoot(root.probVect)
	return totalLK




# find sequence positions with high probability of being errors, and write them to file.
def calculateErrorProbabilities(root,file,minErrorProb):
	if not useRateVariation:
		mutMatrix=mutMatrixGlobal
	node=root
	#direction 0 means from parent, direction 1 means from a child
	lastNode=None
	direction=0
	while node!=None:
		if direction==0:
			if len(node.children)==0:
				file.write(">"+node.name+"\n")
				if len(node.minorSequences)>0:
					for idNode in node.minorSequences:
						file.write(">"+idNode+"\n")
				else:
					if node==node.up.children[0]:
						probVectP=node.up.probVectUpRight
					else:
						probVectP=node.up.probVectUpLeft
					probVectC=node.probVect
					indexEntry1, indexEntry2, pos = 0, 0, 0
					entry1=probVectP[indexEntry1]
					entry2=probVectC[indexEntry2]
					contribLength=node.dist

					while True:
						if entry2[0]==5 or entry1[0]==5: # case N
							pos=min(entry1[1],entry2[1])
						else:
							totLen1=node.dist
							if entry1[0]<5:
								if len(entry1)==3+usingErrorRate:
									totLen1+=entry1[2]
								elif len(entry1)==4+usingErrorRate:
									totLen1+=entry1[3]
									#consider contribution across the root twice, each time only considering the bit on the own side of the root,
							else:
								if len(entry1)>3:
									totLen1+=entry1[2]
							
							totLen2=0.0
							if entry2[0]<5:
								if len(entry2)>2+usingErrorRate:
									totLen2+=entry2[2]
							else:
								if len(entry2)>3:
									totLen2+=entry2[2]
							if entry1[0]==4: # case entry1 is R
								if entry2[0]==4:
									pos=min(entry1[1],entry2[1])
								
								elif entry2[0]==6:
									i1=refIndeces[pos]
									if entry2[-1][i1]<0.1:
										if useRateVariation:
											mutMatrix=mutMatrices[pos]
										errorRate = errorRates[pos]
										numAltNucs=0
										for i in range4:
											if entry2[-1][i]>0.1:
												numAltNucs+=1
										if len(entry1)==4+usingErrorRate:
											errProb=rootFreqs[i1]*(1.0+mutMatrix[i1][i1]*(totLen1+entry1[2]))*errorRate*0.33333*numAltNucs
											mutProb=0.0
											i1rootProb=rootFreqs[i1]*(1.0+mutMatrix[i1][i1]*entry1[2])
											for i in range4:
												if entry2[-1][i]>0.1:
													mutProb+=i1rootProb*mutMatrix[i1][i]*totLen1
													mutProb+=rootFreqs[i]*(1.0+mutMatrix[i][i]*totLen1)*mutMatrix[i][i1]*entry1[2]
											errProb=errProb/(errProb+mutProb)
										else:
											errProb=(1.0+mutMatrix[i1][i1]*totLen1)*errorRate*0.33333*numAltNucs
											mutProb=0.0
											for i in range4:
												if entry2[-1][i]>0.1:
													mutProb+=mutMatrix[i1][i]*totLen1
											errProb=errProb/(errProb+mutProb)
										if errProb>=minErrorProb:
											file.write(str(entry2[1])+"\t"+"X"+"\t"+str(errProb)+"\n")
									pos+=1
								
								else: #entry1 is reference and entry2 is a different but single nucleotide
									i1=refIndeces[pos]
									i2=entry2[0]
									if useRateVariation:
										mutMatrix=mutMatrices[pos]
									if errorRateSiteSpecific: errorRate = errorRates[pos]
									if len(entry1)<4+usingErrorRate:
										errorProb=errorRate*0.33333
										mutProb=mutMatrix[i1][i2]*totLen1
										errorProb=errorProb/(errorProb+mutProb)
									else:
										contribLength=node.dist+entry1[3]
										mutprob1=rootFreqs[i1]*mutMatrix[i1][i2]*contribLength
										mutprob2=rootFreqs[i2]*mutMatrix[i2][i1]*entry1[2]
										errorProb=rootFreqs[i1]*errorRate*0.33333
										normalization=mutprob1+mutprob2+errorProb
										errorProb=errorProb/normalization
									if errorProb>=minErrorProb:
										file.write(str(entry2[1])+"\t"+allelesList[i2]+"\t"+str(errorProb)+"\n")
									pos+=1

							# entry1 is of type "O"
							elif entry1[0]==6:
									normalization=0.0
									if useRateVariation:
										mutMatrix=mutMatrices[pos]
									if entry2[0]==6:
										if errorRateSiteSpecific: errorRate = errorRates[pos]
										noMutProb=0.0
										mutProb=0.0
										errorProb=0.0
										for j in range4:
											if entry2[-1][j]>0.1:
												noMutProb+=entry1[-1][j]
												errorProb+=(1.0-entry1[-1][j])*errorRate*0.33333
												for i in range4:
													if j!=i:
														mutProb+=entry1[-1][i]*mutMatrix[i][j]*totLen1
										normalization=errorProb+noMutProb+mutProb
										errorProb=errorProb/normalization
										if errorProb>=minErrorProb:
											file.write(str(entry2[1])+"\t"+"X"+"\t"+str(errorProb)+"\n")
									else:
										if entry2[0]==4:
											i2=refIndeces[pos]
										else:
											i2=entry2[0]
										if errorRateSiteSpecific: errorRate = errorRates[pos]
										errorProb=(1.0-entry1[-1][i2])*errorRate*0.33333
										noMutProb=entry1[-1][i2]
										mutProb=0.0
										for i in range4:
											if i!=i2:
												mutProb+=entry1[-1][i]*mutMatrix[i][i2]*totLen1
										normalization=errorProb+noMutProb+mutProb
										errorProb=errorProb/normalization
										if errorProb>=minErrorProb:
											file.write(str(entry1[1])+"\t"+allelesList[i2]+"\t"+str(errorProb)+"\n")
									pos+=1

							else: #entry1 is a non-ref nuc
								i1=entry1[0]
								if entry2[0]==i1:
									pass
								else: #entry1 and entry2 are of different types
									if useRateVariation:
										mutMatrix=mutMatrices[pos]
									if entry2[0]==6:
										if entry2[-1][i1]>0.1:#in case the reference allele is possible, we ignore the other possibilities since they are unlikely
											pass
										else:
											if errorRateSiteSpecific: errorRate = errorRates[pos]
											numAltNucs=0
											for i in range4:
												if entry2[-1][i]>0.1:
													numAltNucs+=1

											if len(entry1)==4+usingErrorRate:
												errProb=rootFreqs[i1]*(1.0+mutMatrix[i1][i1]*(totLen1+entry1[2]))*errorRate*0.33333*numAltNucs
												mutProb=0.0
												i1rootProb=rootFreqs[i1]*(1.0+mutMatrix[i1][i1]*entry1[2])
												for i in range4:
													if entry2[-1][i]>0.1:
														mutProb+=i1rootProb*mutMatrix[i1][i]*totLen1
														mutProb+=rootFreqs[i]*(1.0+mutMatrix[i][i]*totLen1)*mutMatrix[i][i1]*entry1[2]
												errProb=errProb/(errProb+mutProb)

											else:
												errProb=(1.0+mutMatrix[i1][i1]*totLen1)*errorRate*0.33333*numAltNucs
												mutProb=0.0
												for i in range4:
													if entry2[-1][i]>0.1:
														mutProb+=mutMatrix[i1][i]*totLen1
												errProb=errProb/(errProb+mutProb)
											if errProb>=minErrorProb:
												file.write(str(entry2[1])+"\t"+"X"+"\t"+str(errProb)+"\n")

									#entry2 is a nucleotide type (like entry1)		
									else:
										if entry2[0]==4:
											i2=refIndeces[pos]
										else:
											i2=entry2[0]
										if errorRateSiteSpecific: errorRate = errorRates[pos]
										if len(entry1)<4+usingErrorRate:
											errorProb=errorRate*0.33333
											mutProb=mutMatrix[i1][i2]*totLen1
											errorProb=errorProb/(errorProb+mutProb)
										else:
											contribLength=node.dist+entry1[3]
											mutprob1=rootFreqs[i1]*mutMatrix[i1][i2]*contribLength
											mutprob2=rootFreqs[i2]*mutMatrix[i2][i1]*entry1[2]
											errorProb=rootFreqs[i1]*errorRate*0.33333
											normalization=mutprob1+mutprob2+errorProb
											errorProb=errorProb/normalization
										if errorProb>=minErrorProb:
											file.write(str(entry1[1])+"\t"+allelesList[i2]+"\t"+str(errorProb)+"\n")
											
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




#Given a tree and its genome lists, calculate mutation counts and waiting times for expectation maximization estimation of substitution rates
def expectationMaximizationCalculationRates(root,trackMutations=False):
	if not useRateVariation:
		mutMatrix=mutMatrixGlobal
	if usingErrorRate and (not errorRateSiteSpecific):
		errorRate=errorRateGlobal
	node=root
	#direction 0 means from parent, direction 1 means from a child
	lastNode=None
	direction=0
	counts=[[0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0]]
	waitingTimes=[0.0,0.0,0.0,0.0]
	numTips=0
	if usingErrorRate:
		errorCount=0.0
		observedTotNucsSites=0
		#if errorRateSiteSpecific: errorRate = errorRates[pos]
		if errorRateSiteSpecific:
			observedNucsSites=[]
			errorCountSites=[]
			for i in range(lRef):
				observedNucsSites.append(0.0)
				errorCountSites.append(0.0)
			observedNucsSites.append(0.0)
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
			if len(node.children)==0:
				nodeIsLeaf=True
				numTips+=1
			else:
				nodeIsLeaf=False

			if (node.dist or usingErrorRate) and node.up!=None:
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
				if trackMutations:
					node.Ns=[]
					node.mutations=[]
					if usingErrorRate and nodeIsLeaf:
						node.errors=[]

				while True:
					if entry2[0]==5: # case N
						if usingErrorRate and nodeIsLeaf:
							if errorRateSiteSpecific:
								observedNucsSites[pos]-=1
							else:
								observedTotNucsSites-=(min(entry1[1],entry2[1])-pos)
						if useRateVariation:
							trackingNs[pos]-=node.dist
						if trackMutations:
							if (not node.Ns) or (type(node.Ns[-1])==int or node.Ns[-1][1]!=entry2[1]):
								node.Ns.append((pos+1,entry2[1]))
						pos=min(entry1[1],entry2[1])
						if useRateVariation:
							trackingNs[pos]+=node.dist
						if errorRateSiteSpecific and nodeIsLeaf:
							observedNucsSites[pos]+=1

					elif entry1[0]==5: # case N
						if useRateVariation:
							trackingNs[pos]-=node.dist
						if usingErrorRate and nodeIsLeaf: #we treat leaf nodes as if they were not observed, since there is no information regarding errors in this case.
							if errorRateSiteSpecific:
								observedNucsSites[pos]-=1
							else:
								observedTotNucsSites-=(min(entry1[1],entry2[1])-pos)
						pos=min(entry1[1],entry2[1])
						if errorRateSiteSpecific and nodeIsLeaf:
							observedNucsSites[pos]+=1
						if useRateVariation:
							trackingNs[pos]+=node.dist
					else:
						totLen1=node.dist
						if entry1[0]<5:
							if len(entry1)==3+usingErrorRate:
								totLen1+=entry1[2]
							elif len(entry1)==4+usingErrorRate:
								totLen1+=entry1[3]
								#consider contribution across the root twice, each time only considering the bit on the own side of the root,
						else:
							if len(entry1)>3:
								totLen1+=entry1[2]
						
						totLen2=0.0
						if entry2[0]<5:
							if len(entry2)>2+usingErrorRate:
								totLen2+=entry2[2]
						else:
							if len(entry2)>3:
								totLen2+=entry2[2]
						if entry1[0]==4: # case entry1 is R
							if entry2[0]==4:
								end=min(entry1[1],entry2[1])
								if (not totLen2) and node.dist:
									for i in range4:
										waitingTimes[i]+=totLen1*(cumulativeBases[end][i]-cumulativeBases[pos][i])
								pos=end
							
							elif entry2[0]==6:
								i1=refIndeces[pos]
								if trackMutations and nodeIsLeaf:
									node.Ns.append(entry2[1])
								if entry2[-1][i1]>0.1:#in case the reference allele is possible, we ignore the other possibilities since they are unlikely
									waitingTimes[i1]+=totLen1
								elif nodeIsLeaf and usingErrorRate:
									if useRateVariation:
										mutMatrix=mutMatrices[pos]
									if errorRateSiteSpecific: errorRate = errorRates[pos]
									numAltNucs=0
									for i in range4:
										if entry2[-1][i]>0.1:
											numAltNucs+=1

									if len(entry1)==4+usingErrorRate:
										stayProb1=1.0+mutMatrix[i1][i1]*totLen1
										if stayProb1<0:
											approxFailed1=True
											stayProb1=0.25
										else:
											approxFailed1=False
										stayProb2=1.0+mutMatrix[i1][i1]*entry1[2]
										if stayProb2<0:
											approxFailed2=True
											stayProb2=0.25
										else:
											approxFailed2=False
										errProb=rootFreqs[i1]*stayProb1*stayProb2*errorRate*0.33333*numAltNucs
										mutProb=0.0
										i1rootProb=rootFreqs[i1]*stayProb2
										for i in range4:
											if entry2[-1][i]>0.1:
												stayProb1=1.0+mutMatrix[i][i]*totLen1
												if stayProb1<0:
													approxFailed1=True
													stayProb1=0.25
												else:
													approxFailed1=False
												if approxFailed1:
													mutProb+=i1rootProb*0.25
												else:
													mutProb+=i1rootProb*mutMatrix[i1][i]*totLen1
												if approxFailed2:
													mutProb+=rootFreqs[i]*stayProb1*0.25
												else:
													mutProb+=rootFreqs[i]*stayProb1*mutMatrix[i][i1]*entry1[2]
										normalization=errProb+mutProb
										errProb=errProb/normalization
										mutProb=mutProb/normalization
										contribLength=node.dist+entry1[3]
										if useRateVariation:
											waitingTimesSites[pos][i1]-=contribLength*(1.0-errProb)
										waitingTimes[i1]+=contribLength*errProb
										errorCount+=errProb
										if errorRateSiteSpecific:
											errorCountSites[pos]+=errProb
										for i in range4:#i is the root nucleotide
											if entry2[-1][i]>0.1:
												stayProb1=1.0+mutMatrix[i][i]*totLen1
												if stayProb1<0:
													approxFailed1=True
													stayProb1=0.25
												else:
													approxFailed1=False
												if approxFailed1:
													prob1=i1rootProb*0.25/normalization
												else:
													prob1=i1rootProb*mutMatrix[i1][i]*totLen1/normalization
												if approxFailed2:
													probi=rootFreqs[i]*stayProb1*0.25/normalization
												else:
													probi=rootFreqs[i]*stayProb1*mutMatrix[i][i1]*entry1[2]/normalization
												waitingTimes[i]+=contribLength*(probi+prob1/2)
												waitingTimes[i1]+=contribLength*prob1/2
												counts[i1][i]+=prob1
												if useRateVariation:
													waitingTimesSites[pos][i]+=contribLength*(probi+prob1/2)
													waitingTimesSites[pos][i1]+=contribLength*prob1/2
													countsSites[pos]+=prob1
													if prob1<0.0:
														raise Exception("exit")

									else:
										stayProb=1.0+mutMatrix[i1][i1]*totLen1
										if stayProb<0:
											approxFailed=True
											stayProb=0.25
										else:
											approxFailed=False
										errProb=stayProb*errorRate*0.33333*numAltNucs
										mutProb=0.0
										for i in range4:
											if entry2[-1][i]>0.1:
												if approxFailed:
													mutProb+=0.25
												else:
													mutProb+=mutMatrix[i1][i]*totLen1
										normalization=errProb+mutProb
										errProb=errProb/normalization
										mutProb=mutProb/normalization
										if useRateVariation:
											waitingTimesSites[pos][i1]-=totLen1*(1.0-errProb)
										waitingTimes[i1]+=totLen1*errProb
										errorCount+=errProb
										if errorRateSiteSpecific:
											errorCountSites[pos]+=errProb
										for i in range4:
											if entry2[-1][i]>0.1:
												prob=mutMatrix[i1][i]*totLen1/normalization
												waitingTimes[i1]+=(totLen1/2)*prob
												waitingTimes[i]+=(totLen1/2)*prob
												counts[i1][i]+=prob
												if useRateVariation:
													waitingTimesSites[pos][i1]+=(totLen1/2)*prob
													waitingTimesSites[pos][i]+=(totLen1/2)*prob
													countsSites[pos]+=prob
													if prob<0.0:
														raise Exception("exit")

								elif not totLen2:
									if useRateVariation:
										mutMatrix=mutMatrices[pos]
									normalization=0.0
									if len(entry1)==4+usingErrorRate:
										contribLength=node.dist+entry1[3]
										if useRateVariation:
											waitingTimesSites[pos][i1]-=contribLength
										stayProb1=1.0+mutMatrix[i1][i1]*entry1[2]
										if stayProb1<0:
											approxFailed1=True
											stayProb1=0.25
										else:
											approxFailed1=False
										for i in range4:
											stayProb2=1.0+mutMatrix[i][i]*contribLength
											if stayProb2<0:
												approxFailed2=True
												stayProb2=0.25
											else:
												approxFailed2=False
											if i1==i:
												prob=rootFreqs[i]*stayProb1
												if approxFailed2:
													tot3=0.25
												else:
													tot3=0.0
													for j in range4:
														tot3+=mutMatrix[i][j]*entry2[-1][j]
													tot3*=contribLength
													tot3+=entry2[-1][i]
												normalization+=prob*tot3
											else:
												if approxFailed1:
													prob=rootFreqs[i]*0.25*stayProb2*entry2[-1][i]
												else:
													prob=rootFreqs[i]*mutMatrix[i][i1]*entry1[2]*stayProb2*entry2[-1][i]
												normalization+=prob
										for i in range4:
											stayProb2=1.0+mutMatrix[i][i]*contribLength
											if stayProb2<0:
												approxFailed2=True
												stayProb2=0.25
											else:
												approxFailed2=False
											if i1==i:
												prob=rootFreqs[i]*stayProb1
												for j in range4:
													if j==i:
														tot3=prob*stayProb2*entry2[-1][j]/normalization
														waitingTimes[i]+=contribLength*tot3
														if useRateVariation:
															waitingTimesSites[pos][i]+=contribLength*tot3
													else:
														if approxFailed2:
															tot3=prob*0.25*entry2[-1][j]/normalization
														else:
															tot3=prob*mutMatrix[i][j]*contribLength*entry2[-1][j]/normalization
														waitingTimes[i]+=(contribLength/2)*tot3
														waitingTimes[j]+=(contribLength/2)*tot3
														counts[i][j]+=tot3
														if trackMutations and (not nodeIsLeaf) and (tot3>minMutProb):
															node.mutations.append((i1,pos+1,j,tot3))
														if useRateVariation:
															waitingTimesSites[pos][i]+=(contribLength/2)*tot3
															waitingTimesSites[pos][j]+=(contribLength/2)*tot3
															countsSites[pos]+=tot3
															if tot3<0.0:
																raise Exception("exit")
											else:
												if approxFailed1:
													prob=rootFreqs[i]*0.25*stayProb2*entry2[-1][i]/normalization
												else:
													prob=rootFreqs[i]*mutMatrix[i][i1]*entry1[2]*stayProb2*entry2[-1][i]/normalization
												waitingTimes[i]+=contribLength*prob
												if useRateVariation:
													waitingTimesSites[pos][i]+=contribLength*prob

									else:
										if useRateVariation:
											waitingTimesSites[pos][i1]-=totLen1
										stayProb=1.0+mutMatrix[i1][i1]*totLen1
										if stayProb<0:
											normalization=0.25
											approxFailed=True
										else:
											approxFailed=False
											for i in range4:
												if i1==i:
													normalization+=stayProb*entry2[-1][i]
												else:
													normalization+=mutMatrix[i1][i]*totLen1*entry2[-1][i]
										for i in range4:
											if i1==i:
												if approxFailed:
													prob=entry2[-1][i]
												else:
													prob=(1.0+mutMatrix[i][i]*totLen1)*entry2[-1][i]/normalization
												waitingTimes[i]+=totLen1*prob
												if useRateVariation:
													waitingTimesSites[pos][i]+=totLen1*prob
											else:
												if approxFailed:
													prob=entry2[-1][i]
												else:
													prob=mutMatrix[i1][i]*totLen1*entry2[-1][i]/normalization
												waitingTimes[i1]+=(totLen1/2)*prob
												waitingTimes[i]+=(totLen1/2)*prob
												counts[i1][i]+=prob
												if trackMutations and (not nodeIsLeaf) and (prob>minMutProb):
													node.mutations.append((i1,pos+1,i,prob))
												if useRateVariation:
													waitingTimesSites[pos][i1]+=(totLen1/2)*prob
													waitingTimesSites[pos][i]+=(totLen1/2)*prob
													countsSites[pos]+=prob
													if prob<0.0:
														raise Exception("exit")
								pos+=1
							
							else: #entry1 is reference and entry2 is a different but single nucleotide
								i1=refIndeces[pos]
								i2=entry2[0]
								if useRateVariation:
									mutMatrix=mutMatrices[pos]
									
								if nodeIsLeaf and usingErrorRate:
									if errorRateSiteSpecific: errorRate = errorRates[pos]
									if len(entry1)<4+usingErrorRate:
										errorProb=errorRate*0.33333
										mutProb=mutMatrix[i1][i2]*totLen1
										normalization=errorProb+mutProb
										errorProb=errorProb/normalization
										mutProb=mutProb/normalization
										if useRateVariation:
											waitingTimesSites[pos][i1]+=totLen1*(mutProb/2-1.0)
											waitingTimesSites[pos][i2]+=totLen1*(errorProb+mutProb/2)
											countsSites[pos]+=mutProb
											if mutProb<0.0:
												raise Exception("exit")
										waitingTimes[i1]+=totLen1*(errorProb+mutProb/2)
										waitingTimes[i2]+=(totLen1*mutProb/2)
										counts[i1][i2]+=mutProb
										if trackMutations:
											if (mutProb>minMutProb):
												node.mutations.append((i1,pos+1,i2,mutProb))
											if (errorProb>minMutProb):
												node.errors.append((i1,pos+1,i2,errorProb))
										errorCount+=errorProb
										if errorRateSiteSpecific:
											errorCountSites[pos]+=errorProb
									else:
										contribLength=node.dist+entry1[3]
										mutprob1=rootFreqs[i1]*mutMatrix[i1][i2]*contribLength
										mutprob2=rootFreqs[i2]*mutMatrix[i2][i1]*entry1[2]
										errorProb=rootFreqs[i1]*errorRate*0.33333
										normalization=mutprob1+mutprob2+errorProb
										mutprob1=mutprob1/normalization
										mutprob2=mutprob2/normalization
										errorProb=errorProb/normalization
										waitingTimes[i1]+=(contribLength)*(mutprob1/2+errorProb)
										waitingTimes[i2]+=(contribLength)*(mutprob2+mutprob1/2)
										counts[i1][i2]+=mutprob1
										if trackMutations:
											if (mutprob1>minMutProb):
												node.mutations.append((i1,pos+1,i2,mutprob1))
											if (errorProb>minMutProb):
												node.errors.append((i1,pos+1,i2,errorProb))
										errorCount+=errorProb
										if errorRateSiteSpecific:
											errorCountSites[pos]+=errorProb
										if useRateVariation:
											waitingTimesSites[pos][i1]+=(contribLength)*(mutprob1/2+errorProb-1.0)
											waitingTimesSites[pos][i2]+=(contribLength)*(mutprob2+mutprob1/2)
											countsSites[pos]+=mutprob1
											if mutprob1<0.0:
												raise Exception("exit")

								elif not totLen2:
									if len(entry1)<4+usingErrorRate:
										if useRateVariation:
											waitingTimesSites[pos][i1]-=totLen1/2
											waitingTimesSites[pos][i2]+=totLen1/2
											countsSites[pos]+=1
										waitingTimes[i1]+=(totLen1/2)
										waitingTimes[i2]+=(totLen1/2)
										counts[i1][i2]+=1
										if trackMutations:
											node.mutations.append((i1,pos+1,i2,1.0))
									else:
										contribLength=node.dist+entry1[3]
										noMutProb1=1.0+mutMatrix[i1][i1]*entry1[2]
										if noMutProb1<0:
											noMutProb1=0.25
										noMutProb2=1.0+mutMatrix[i2][i2]*contribLength
										if noMutProb2<0:
											noMutProb2=0.25
										prob1=rootFreqs[i1]*mutMatrix[i1][i2]*contribLength*noMutProb1
										prob2=rootFreqs[i2]*mutMatrix[i2][i1]*entry1[2]*noMutProb2
										normalization=prob1+prob2
										prob1=prob1/normalization
										prob2=prob2/normalization
										waitingTimes[i1]+=(contribLength/2)*prob1
										waitingTimes[i2]+=(contribLength/2)*prob1
										counts[i1][i2]+=prob1
										if trackMutations:
											if (prob1>minMutProb):
												node.mutations.append((i1,pos+1,i2,prob1))
										waitingTimes[i2]+=contribLength*prob2
										if useRateVariation:
											waitingTimesSites[pos][i1]-=contribLength
											waitingTimesSites[pos][i1]+=(contribLength/2)*prob1
											waitingTimesSites[pos][i2]+=(contribLength/2)*prob1
											waitingTimesSites[pos][i2]+=contribLength*prob2
											countsSites[pos]+=prob1
											if prob1<0.0:
												raise Exception("exit")
								pos+=1

						# entry1 is of type "O"
						elif entry1[0]==6:
							if not totLen2:
								normalization=0.0
								if useRateVariation:
									mutMatrix=mutMatrices[pos]
									waitingTimesSites[pos][refIndeces[pos]]-=totLen1
								if entry2[0]==6:
									if trackMutations and nodeIsLeaf:
										node.Ns.append(pos+1)
									if nodeIsLeaf and usingErrorRate:
										if errorRateSiteSpecific: errorRate = errorRates[pos]
										noMutProb=0.0
										mutProb=0.0
										errorProb=0.0
										for j in range4:
											if entry2[-1][j]>0.1:
												noMutProb+=entry1[-1][j]
												errorProb+=(1.0-entry1[-1][j])*errorRate*0.33333
												for i in range4:
													if j!=i:
														mutProb+=entry1[-1][i]*mutMatrix[i][j]*totLen1
										normalization=errorProb+noMutProb+mutProb
										errorProb=errorProb/normalization
										noMutProb=noMutProb/normalization
										mutProb=mutProb/normalization
										errorCount+=errorProb
										if errorRateSiteSpecific:
											errorCountSites[pos]+=errorProb
										for j in range4:
											if entry2[-1][j]>0.1:
												waitingTimes[j]+=totLen1*entry1[-1][j]/normalization
												if useRateVariation:
													waitingTimesSites[pos][j]+=totLen1*entry1[-1][j]/normalization
												for i in range4:
													if j!=i:
														mutProbij=entry1[-1][i]*mutMatrix[i][j]*totLen1/normalization
														waitingTimes[j]+=totLen1*mutProbij/2
														waitingTimes[i]+=totLen1*mutProbij/2
														counts[i][j]+=mutProbij
														if useRateVariation:
															waitingTimesSites[pos][j]+=totLen1*mutProbij/2
															waitingTimesSites[pos][i]+=totLen1*mutProbij/2
															countsSites[pos]+=mutProbij
															if mutProbij<0.0:
																raise Exception("exit")
									else:
										approxFailed=[False,False,False,False]
										for i in range4:
											stayProb=1.0+mutMatrix[i][i]*totLen1
											if stayProb<0:
												for j in range4:
													normalization+=entry1[-1][i]*0.25*entry2[-1][j]
												approxFailed[i]=True
											else:
												for j in range4:
													if i==j:
														normalization+=entry1[-1][i]*stayProb*entry2[-1][j]
													else:
														normalization+=entry1[-1][i]*mutMatrix[i][j]*totLen1*entry2[-1][j]
										for i in range4:
											for j in range4:
												if i==j:
													if approxFailed[i]:
														prob=entry1[-1][i]*0.25*entry2[-1][j]/normalization
													else:
														prob=entry1[-1][i]*(1.0+mutMatrix[i][i]*totLen1)*entry2[-1][j]/normalization
													waitingTimes[i]+=totLen1*prob
													if useRateVariation:
														waitingTimesSites[pos][i]+=totLen1*prob
												else:
													if approxFailed[i]:
														prob=entry1[-1][i]*0.25*entry2[-1][j]/normalization
													else:
														prob=entry1[-1][i]*mutMatrix[i][j]*totLen1*entry2[-1][j]/normalization
													waitingTimes[i]+=(totLen1/2)*prob
													waitingTimes[j]+=(totLen1/2)*prob
													counts[i][j]+=prob
													if trackMutations:
														if (prob>minMutProb):
															node.mutations.append((i,pos+1,j,prob))
													if useRateVariation:
														waitingTimesSites[pos][i]+=(totLen1/2)*prob
														waitingTimesSites[pos][j]+=(totLen1/2)*prob
														countsSites[pos]+=prob
														if prob<0.0:
															raise Exception("exit")
								else:
									if entry2[0]==4:
										i2=refIndeces[pos]
									else:
										i2=entry2[0]
									if nodeIsLeaf and usingErrorRate:
										if errorRateSiteSpecific: errorRate = errorRates[pos]
										errorProb=(1.0-entry1[-1][i2])*errorRate*0.33333
										noMutProb=entry1[-1][i2]
										mutProb=0.0
										for i in range4:
											if i!=i2:
												mutProb+=entry1[-1][i]*mutMatrix[i][i2]*totLen1
										normalization=errorProb+noMutProb+mutProb
										errorProb=errorProb/normalization
										noMutProb=noMutProb/normalization
										mutProb=mutProb/normalization
										errorCount+=errorProb
										if trackMutations:
											if (errorProb>minMutProb):
												node.errors.append((5,pos+1,i2,errorProb))
										if errorRateSiteSpecific:
											errorCountSites[pos]+=errorProb
										waitingTimes[i2]+=totLen1*noMutProb
										waitingTimes[i2]+=(totLen1/2)*mutProb
										if useRateVariation:
											waitingTimesSites[pos][i2]+=totLen1*noMutProb
											waitingTimesSites[pos][i2]+=totLen1*mutProb/2
											countsSites[pos]+=mutProb
											if mutProb<0.0:
												raise Exception("exit")
										for i in range4:
											if i!=i2:
												prob=entry1[-1][i]*mutMatrix[i][i2]*totLen1/normalization
												probErr=entry1[-1][i]*errorRate*0.33333/normalization
												waitingTimes[i]+=totLen1*(probErr+prob/2)
												counts[i][i2]+=prob
												if trackMutations:
													if (prob>minMutProb):
														node.mutations.append((i,pos+1,i2,prob))
												if useRateVariation:
													waitingTimesSites[pos][i]+=totLen1*(probErr+prob/2)
									else:
										stayProb=1.0+mutMatrix[i2][i2]*totLen1
										if stayProb<0:
											normalization=0.25
											approxFailed=True
										else:
											approxFailed=False
											for i in range4:
												if i==i2:
													normalization+=entry1[-1][i]*stayProb
												else:
													normalization+=entry1[-1][i]*mutMatrix[i][i2]*totLen1
										for i in range4:
											if i==i2:
												if approxFailed:
													prob=entry1[-1][i]
												else:
													prob=entry1[-1][i]*(1.0+mutMatrix[i][i]*totLen1)/normalization
												waitingTimes[i]+=totLen1*prob
												if useRateVariation:
													waitingTimesSites[pos][i]+=totLen1*prob
											else:
												if approxFailed:
													prob=entry1[-1][i]
												else:
													prob=entry1[-1][i]*mutMatrix[i][i2]*totLen1/normalization
												waitingTimes[i]+=(totLen1/2)*prob
												waitingTimes[i2]+=(totLen1/2)*prob
												counts[i][i2]+=prob
												if trackMutations:
													if (prob>minMutProb):
														node.mutations.append((i,pos+1,i2,prob))
												if useRateVariation:
													waitingTimesSites[pos][i]+=(totLen1/2)*prob
													waitingTimesSites[pos][i2]+=(totLen1/2)*prob
													countsSites[pos]+=prob
													if prob<0.0:
														print("Negative value, numerical issue")
														print(node.dist)
														print(node)
														for child in node.children:
															print(child)
														print(entry1[-1][i])
														print(mutMatrix[i][i2])
														print(totLen1)
														print(normalization)
														print(mutMatrix)
														print(entry1[-1])
														print(i2)
														raise Exception("exit")
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
								if useRateVariation:
									mutMatrix=mutMatrices[pos]
								if entry2[0]==6:
									if trackMutations and nodeIsLeaf:
										node.Ns.append(pos+1)
									if entry2[-1][i1]>0.1:#in case the reference allele is possible, we ignore the other possibilities since they are unlikely
										waitingTimes[i1]+=totLen1
										if useRateVariation:
											waitingTimesSites[pos][i1]+=totLen1
											waitingTimesSites[pos][refIndeces[pos]]-=totLen1
									elif nodeIsLeaf and usingErrorRate:
										if errorRateSiteSpecific: errorRate = errorRates[pos]
										numAltNucs=0
										for i in range4:
											if entry2[-1][i]>0.1:
												numAltNucs+=1

										if len(entry1)==4+usingErrorRate:
											stayProb1=1.0+mutMatrix[i1][i1]*totLen1
											if stayProb1<0:
												approxFailed1=True
												stayProb1=0.25
											else:
												approxFailed1=False
											stayProb2=1.0+mutMatrix[i1][i1]*entry1[2]
											if stayProb2<0:
												approxFailed2=True
												stayProb2=0.25
											else:
												approxFailed2=False
											errProb=rootFreqs[i1]*stayProb1*stayProb2*errorRate*0.33333*numAltNucs
											mutProb=0.0
											i1rootProb=rootFreqs[i1]*stayProb2
											for i in range4:
												if entry2[-1][i]>0.1:
													stayProb1=1.0+mutMatrix[i][i]*totLen1
													if stayProb1<0:
														approxFailed1=True
														stayProb1=0.25
													else:
														approxFailed1=False
													if approxFailed1:
														mutProb+=i1rootProb*0.25
													else:
														mutProb+=i1rootProb*mutMatrix[i1][i]*totLen1
													if approxFailed2:
														mutProb+=rootFreqs[i]*stayProb1*0.25
													else:
														mutProb+=rootFreqs[i]*stayProb1*mutMatrix[i][i1]*entry1[2]
											normalization=errProb+mutProb
											errProb=errProb/normalization
											mutProb=mutProb/normalization
											contribLength=node.dist+entry1[3]
											if useRateVariation:
												waitingTimesSites[pos][i1]-=contribLength*(1.0-errProb)
											waitingTimes[i1]+=contribLength*errProb
											errorCount+=errProb
											if errorRateSiteSpecific:
												errorCountSites[pos]+=errProb
											for i in range4:#i is the root nucleotide
												if entry2[-1][i]>0.1:
													stayProb1=1.0+mutMatrix[i][i]*totLen1
													if stayProb1<0:
														approxFailed1=True
														stayProb1=0.25
													else:
														approxFailed1=False
													if approxFailed1:
														prob1=i1rootProb*0.25/normalization
													else:
														prob1=i1rootProb*mutMatrix[i1][i]*totLen1/normalization
													if approxFailed2:
														probi=rootFreqs[i]*stayProb1*0.25/normalization
													else:
														probi=rootFreqs[i]*stayProb1*mutMatrix[i][i1]*entry1[2]/normalization
													waitingTimes[i]+=contribLength*(probi+prob1/2)
													waitingTimes[i1]+=contribLength*prob1/2
													counts[i1][i]+=prob1
													if useRateVariation:
														waitingTimesSites[pos][i]+=contribLength*(probi+prob1/2)
														waitingTimesSites[pos][i1]+=contribLength*prob1/2
														countsSites[pos]+=prob1
														if prob1<0.0:
															raise Exception("exit")

										else:
											stayProb=1.0+mutMatrix[i1][i1]*totLen1
											if stayProb<0:
												approxFailed=True
												stayProb=0.25
											else:
												approxFailed=False
											errProb=stayProb*errorRate*0.33333*numAltNucs
											mutProb=0.0
											for i in range4:
												if entry2[-1][i]>0.1:
													if approxFailed:
														mutProb+=0.25
													else:
														mutProb+=mutMatrix[i1][i]*totLen1
											normalization=errProb+mutProb
											errProb=errProb/normalization
											mutProb=mutProb/normalization
											if useRateVariation:
												waitingTimesSites[pos][i1]-=totLen1*(1.0-errProb)
											waitingTimes[i1]+=totLen1*errProb
											errorCount+=errProb
											if errorRateSiteSpecific:
												errorCountSites[pos]+=errProb
											for i in range4:
												if entry2[-1][i]>0.1:
													if approxFailed:
														prob=0.25/normalization
													else:
														prob=mutMatrix[i1][i]*totLen1/normalization
													waitingTimes[i1]+=(totLen1/2)*prob
													waitingTimes[i]+=(totLen1/2)*prob
													counts[i1][i]+=prob
													if useRateVariation:
														waitingTimesSites[pos][i1]+=(totLen1/2)*prob
														waitingTimesSites[pos][i]+=(totLen1/2)*prob
														countsSites[pos]+=prob
														if prob<0.0:
															raise Exception("exit")

									elif not totLen2:
										normalization=0.0
										if len(entry1)==4+usingErrorRate:
											contribLength=node.dist+entry1[3]
											if useRateVariation:
												waitingTimesSites[pos][refIndeces[pos]]-=contribLength
											stayProb2=1.0+mutMatrix[i1][i1]*entry1[2]
											if stayProb2<0:
												approxFailed2=True
												stayProb2=0.25
											else:
												approxFailed2=False
											for i in range4:
												stayProb1=1.0+mutMatrix[i][i]*contribLength
												if stayProb1<0:
													approxFailed1=True
													stayProb1=0.25
												else:
													approxFailed1=False
												if i1==i:
													prob=rootFreqs[i]*stayProb2
													if approxFailed1:
														tot3=0.25
													else:
														tot3=0.0
														for j in range4:
															tot3+=mutMatrix[i][j]*entry2[-1][j]
														tot3*=contribLength
														tot3+=entry2[-1][i]
													normalization+=prob*tot3
												else:
													prob=rootFreqs[i]*mutMatrix[i][i1]*entry1[2]*stayProb1*entry2[-1][i]
													normalization+=prob
											for i in range4:
												stayProb1=1.0+mutMatrix[i][i]*contribLength
												if stayProb1<0:
													approxFailed1=True
													stayProb1=0.25
												else:
													approxFailed1=False
												if i1==i:
													prob=rootFreqs[i]*stayProb2
													for j in range4:
														if j==i:
															tot3=prob*stayProb1*entry2[-1][j]/normalization
															waitingTimes[i]+=contribLength*tot3
															if useRateVariation:
																waitingTimesSites[pos][i]+=contribLength*tot3
														else:
															if approxFailed1:
																tot3=prob*0.25*entry2[-1][j]/normalization
															else:
																tot3=prob*mutMatrix[i][j]*contribLength*entry2[-1][j]/normalization
															waitingTimes[i]+=(contribLength/2)*tot3
															waitingTimes[j]+=(contribLength/2)*tot3
															counts[i][j]+=tot3
															if trackMutations:
																if (tot3>minMutProb):
																	node.mutations.append((i,pos+1,j,tot3))
															if useRateVariation:
																waitingTimesSites[pos][i]+=(contribLength/2)*tot3
																waitingTimesSites[pos][j]+=(contribLength/2)*tot3
																countsSites[pos]+=tot3
																if tot3<0.0:
																	raise Exception("exit")
												else:
													if approxFailed2:
														prob=rootFreqs[i]*0.25*stayProb1*entry2[-1][i]
													else:
														prob=rootFreqs[i]*mutMatrix[i][i1]*entry1[2]*stayProb1*entry2[-1][i]
													waitingTimes[i]+=contribLength*prob/normalization
													if useRateVariation:
														waitingTimesSites[pos][i]+=contribLength*prob/normalization
										else:
											if useRateVariation:
												waitingTimesSites[pos][refIndeces[pos]]-=totLen1
											stayProb=1.0+mutMatrix[i1][i1]*totLen1
											if stayProb<0:
												approxFailed=True
												stayProb=0.25
											else:
												approxFailed=False
											for i in range4:
												if i1==i:
													normalization+=stayProb*entry2[-1][i]
												else:
													if approxFailed:
														normalization+=0.25*entry2[-1][i]
													else:
														normalization+=mutMatrix[i1][i]*totLen1*entry2[-1][i]
											for i in range4:
												if i1==i:
													prob=stayProb*entry2[-1][i]/normalization
													waitingTimes[i]+=totLen1*prob
													if useRateVariation:
														waitingTimesSites[pos][i]+=totLen1*prob
												else:
													if approxFailed:
														prob=0.25*entry2[-1][i]/normalization
													else:
														prob=mutMatrix[i1][i]*totLen1*entry2[-1][i]/normalization
													waitingTimes[i1]+=(totLen1/2)*prob
													waitingTimes[i]+=(totLen1/2)*prob
													counts[i1][i]+=prob
													if trackMutations:
														if (prob>minMutProb):
															node.mutations.append((i1,pos+1,i,prob))
													if useRateVariation:
														waitingTimesSites[pos][i1]+=(totLen1/2)*prob
														waitingTimesSites[pos][i]+=(totLen1/2)*prob
														countsSites[pos]+=prob
														if prob<0.0:
															raise Exception("exit")

								#entry2 is a nucleotide type (like entry1)		
								else:
									if not totLen2:
										if entry2[0]==4:
											i2=refIndeces[pos]
										else:
											i2=entry2[0]
										if usingErrorRate and nodeIsLeaf:
											if errorRateSiteSpecific: errorRate = errorRates[pos]
											if len(entry1)<4+usingErrorRate:
												errorProb=errorRate*0.33333
												mutProb=mutMatrix[i1][i2]*totLen1
												normalization=errorProb+mutProb
												errorProb=errorProb/normalization
												mutProb=mutProb/normalization
												if useRateVariation:
													waitingTimesSites[pos][i1]+=totLen1*(-1.0+errorProb+mutProb/2)
													waitingTimesSites[pos][i2]+=totLen1*(mutProb/2)
													countsSites[pos]+=mutProb
													if mutProb<0.0:
														raise Exception("exit")
												waitingTimes[i1]+=totLen1*(errorProb+mutProb/2)
												waitingTimes[i2]+=(totLen1*mutProb/2)
												counts[i1][i2]+=mutProb
												if trackMutations:
													if (mutProb>minMutProb):
														node.mutations.append((i1,pos+1,i2,mutProb))
													if (errorProb>minMutProb):
														node.errors.append((i1,pos+1,i2,errorProb))
												errorCount+=errorProb
												if errorRateSiteSpecific:
													errorCountSites[pos]+=errorProb
											else:
												contribLength=node.dist+entry1[3]
												mutprob1=rootFreqs[i1]*mutMatrix[i1][i2]*contribLength
												mutprob2=rootFreqs[i2]*mutMatrix[i2][i1]*entry1[2]
												errorProb=rootFreqs[i1]*errorRate*0.33333
												normalization=mutprob1+mutprob2+errorProb
												mutprob1=mutprob1/normalization
												mutprob2=mutprob2/normalization
												errorProb=errorProb/normalization
												waitingTimes[i1]+=(contribLength)*(mutprob1/2+errorProb)
												waitingTimes[i2]+=(contribLength)*(mutprob2+mutprob1/2)
												counts[i1][i2]+=mutprob1
												if trackMutations:
													if (mutprob1>minMutProb):
														node.mutations.append((i1,pos+1,i2,mutprob1))
													if (errorProb>minMutProb):
														node.errors.append((i1,pos+1,i2,errorProb))
												errorCount+=errorProb
												if errorRateSiteSpecific:
													errorCountSites[pos]+=errorProb
												if useRateVariation:
													waitingTimesSites[pos][i1]+=(contribLength)*(mutprob1/2+errorProb-1.0)
													waitingTimesSites[pos][i2]+=(contribLength)*(mutprob2+mutprob1/2)
													countsSites[pos]+=mutprob1
													if mutprob1<0.0:
														raise Exception("exit")

										else:
											if len(entry1)<4+usingErrorRate:
												if useRateVariation:
													waitingTimesSites[pos][refIndeces[pos]]-=totLen1
													waitingTimesSites[pos][i1]+=(totLen1/2)
													waitingTimesSites[pos][i2]+=(totLen1/2)
													countsSites[pos]+=1
												waitingTimes[i1]+=(totLen1/2)
												waitingTimes[i2]+=(totLen1/2)
												counts[i1][i2]+=1
												if trackMutations:
													node.mutations.append((i1,pos+1,i2,1.0))
											else:
												contribLength=node.dist+entry1[3]
												prob1=rootFreqs[i1]*mutMatrix[i1][i2]*contribLength
												prob2=rootFreqs[i2]*mutMatrix[i2][i1]*entry1[2]
												if usingErrorRate and entry1[4]:
													if errorRateSiteSpecific: errorRate = errorRates[pos]
													errorProb=rootFreqs[i2]*errorRate*0.33333
													normalization=prob1+prob2+errorProb
												else:
													normalization=prob1+prob2
												if not normalization:
													print("0 normalization issues in EM")
													raise Exception("exit")
												prob1=prob1/normalization
												prob2=prob2/normalization
												if usingErrorRate and entry1[4]:
													errorProb=errorProb/normalization
													waitingTimes[i2]+=contribLength*errorProb
												waitingTimes[i1]+=(contribLength/2)*prob1
												waitingTimes[i2]+=(contribLength/2)*prob1
												counts[i1][i2]+=prob1
												if trackMutations:
													node.mutations.append((i1,pos+1,i2,prob1))
												waitingTimes[i2]+=contribLength*prob2
												if useRateVariation:
													waitingTimesSites[pos][refIndeces[pos]]-=contribLength
													waitingTimesSites[pos][i1]+=(contribLength/2)*prob1
													waitingTimesSites[pos][i2]+=(contribLength/2)*prob1
													waitingTimesSites[pos][i2]+=contribLength*prob2
													countsSites[pos]+=prob1
													if prob1<0.0:
														raise Exception("exit")
													if usingErrorRate and entry1[4]:
														waitingTimesSites[pos][i2]+=contribLength*errorProb

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

	if usingErrorRate:
		observedTotNucsSites+=lRef*numTips
			
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
	elif (not trackMutations) and (not usingErrorRate):
		print("Expectation Maximization for given model "+model+" not implemented yet")
		raise Exception("exit")
	totRate=-(rootFreqs[0]*counts[0][0]+ rootFreqs[1]*counts[1][1]+ rootFreqs[2]*counts[2][2]+ rootFreqs[3]*counts[3][3] )
	if totRate:
		for i in range4:
			for j in range4:
				counts[i][j]=counts[i][j]/totRate
	
	if usingErrorRate:
		errorRateEstimate=errorCount/observedTotNucsSites
		if errorRateSiteSpecific:
			siteErrRates=[]
			observedNuc=numTips
			for i in range(lRef):
				observedNuc+=observedNucsSites[i]
				if observedNuc>0:
					siteErrRates.append(errorCountSites[i]/observedNuc)
				else:
					siteErrRates.append(0.0)
		else:
			siteErrRates=None
	else:
		errorRateEstimate=None
		siteErrRates=None

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
			siteRates[i]=min(100.0,max(0.001,siteRates[i]/totRate))
	else:
		siteRates=None
	return counts, siteRates, errorRateEstimate, siteErrRates





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

if errorRateSiteSpecificFile or errorRateFixed or estimateErrorRate or estimateSiteSpecificErrorRate:
	forgetData=False
else:
	forgetData=True

if inputTree=="" and (writeTreesToFileEveryTheseSteps==0):
	print("Calculating distances with no sample names")
	distances=distancesFromRefPunishNs(data)
	print("Average divergence from refererence: "+str(totDivFromRef[0]/len(distances)))
else:
	print("Calculating distances with sample names")
	#don't samples already in the tree back to the tree again
	distances=distancesFromRefPunishNs(data,samples=data.keys(),samplesInInitialTree=samplesAlreadyInTree)
print("Distances from the reference calculated")



if inputTree=="": #initialize tree to just the initial root sample
	#extract root genome among those closest to the reference but not empty
	firstSample=distances.pop()
	t1=Tree(name=firstSample[1])
	t1.probVect=probVectTerminalNode(data[firstSample[1]])
	if forgetData:
		data[firstSample[1]]=None
	t1.probVectTot=rootVector(t1.probVect,False)
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
	if forgetData:
		data[sample]=None
	if (numSamples%updateSubstMatrixEveryThisSamples)==0:
		if model!="JC":
			if updateSubMatrix(pseudoMutCounts,model,mutMatrixGlobal):
				for i in range4:
					nonMutRates[i]=mutMatrixGlobal[i][i]
				for i in range(lRef):
					cumulativeRate[i+1]=cumulativeRate[i]+nonMutRates[refIndeces[i]]
	if (numSamples%1000)==0:
		print("Sample num "+str(numSamples))
	start=time()
	bestNode , bestScore, bestBranchLengths = findBestParentForNewSample(t1,newPartials,sample)
	timeFinding+=(time()-start)
	if bestBranchLengths!=None:
		start=time()
		newRoot=placeSampleOnTree(bestNode,newPartials,sample,bestScore, bestBranchLengths[0], bestBranchLengths[1], bestBranchLengths[2],pseudoMutCounts)
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
	cutAndPasteNode(nodeToReplace,destination,False,0.0001,-100.0)
	newickString=createNewick(t1)
	print(newickString)

if forgetData:
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


def updateErrorRates(errorRate,errorRates=None):
	global rootFreqsLogErrorCumulative
	if errorRates!=None:
		global cumulativeErrorRate
		cumulativeErrorRate = [0]*(lRef+1) #creating cumulative error rates.
		cumulativeErrorRate[0] = 0.0
		for i in range(0,len(errorRates)):
			cumulativeErrorRate[i+1] = cumulativeErrorRate[i] + errorRates[i] # CodeImprovement: to make estimations more accurate, one could scale these log(1-errorRate) here, which would be of little extra computational demand.
		rootFreqsLogErrorCumulative=[0]*(lRef+1)
		rootFreqsLogErrorCumulative[0]=0.0
		for i in range(0,len(errorRates)):
			rootFreqsLogErrorCumulative[i+1] = rootFreqsLogErrorCumulative[i] + log(rootFreqs[refIndeces[i]]*(1.0-1.33333*errorRates[i])+0.333333*errorRates[i])
	else:
		rootFreqsLogErrorCumulative=[0]*(lRef+1)
		rootFreqsLogErrorCumulative[0]=0.0
		for i in range(0,len(ref)):
			rootFreqsLogErrorCumulative[i+1] = rootFreqsLogErrorCumulative[i] + log(rootFreqs[refIndeces[i]]*(1.0-1.33333*errorRate)+0.333333*errorRate)





if errorRateSiteSpecificFile:
	#check if file exists
	if not os.path.isfile(errorRateSiteSpecificFile):
		print("errorRateSiteSpecific file " + errorRateSiteSpecificFile + " not found, quitting MAPLE tree inference. Use option --errorRateSiteSpecificFile to specify a valid input file.")
		raise Exception("exit")
	fileI = open(errorRateSiteSpecificFile, "r")
	line = fileI.readline()
	fileI.close()
	errorRates = [float(errR) for errR in line.split(", ")]
	if len(errorRates)!=lRef:
		print("Number of error rates in errorRateSiteSpecific file " + errorRateSiteSpecificFile + " is different from the length of the reference - exiting MAPLE.")
		raise Exception("exit")
	errorRateGlobal=sum(errorRates)/lRef
	updateErrorRates(errorRateGlobal,errorRates=errorRates)
	usingErrorRate=True
	errorRateSiteSpecific=True
elif errorRateFixed:
	errorRateGlobal=errorRateFixed
	updateErrorRates(errorRateGlobal)
	usingErrorRate=True
elif estimateErrorRate:
	if errorRateInitial:
		errorRateGlobal=errorRateInitial
	else:
		errorRateGlobal=1.0/lRef
	updateErrorRates(errorRateGlobal)
	usingErrorRate=True
elif estimateSiteSpecificErrorRate:
	if errorRateInitial:
		errorRateGlobal=errorRateInitial
	else:
		errorRateGlobal=1.0/lRef
	errorRates=[errorRateGlobal]*lRef
	updateErrorRates(errorRateGlobal,errorRates=errorRates)
	usingErrorRate=True
	errorRateSiteSpecific=True




def updateMutMatrices(mutMatrix,siteRates=None):
	global nonMutRates
	global cumulativeRate
	for i in range4:
		nonMutRates[i]=mutMatrix[i][i]
	if siteRates!=None:
		global mutMatrices
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
		for i in range(lRef):
			cumulativeRate[i+1]=cumulativeRate[i]+nonMutRates[refIndeces[i]]


def variance(numList):
	# calculate mean
	m = sum(numList) / len(numList)
	# calculate variance using a list comprehension
	var_res = sum((xi - m) ** 2 for xi in numList) / len(numList)
	return m, var_res



#recalculate all genome lists according to the final substitution model
if inputTree=="" or largeUpdate or rateVariation or usingErrorRate:
	start=time()
	if usingErrorRate:
		reCalculateAllGenomeLists(t1,countNodes=True, data=data)
		if not (estimateErrorRate or estimateSiteSpecificErrorRate):
			data.clear()
	else:
		reCalculateAllGenomeLists(t1,countNodes=True)
	
	#if using error rates, run EM iteratively until convergence in likelihood, otherwise only run one EM step
	if model!="JC" or rateVariation or estimateErrorRate or estimateSiteSpecificErrorRate:
		useRateVariation=rateVariation
		if useRateVariation:
			siteRates=[1.0]*lRef
			mutMatrices=[]
			for i in range(lRef):
				mutMatrices.append([])
				for j in range4:
					mutMatrices[i].append(list(mutMatrixGlobal[j]))
		mutMatrixGlobal, siteRates, errorRateGlobal, errorRates = expectationMaximizationCalculationRates(t1)
		print("Initial EM terminated.  ")
		if estimateErrorRate:
			print("Error rate: "+str(errorRateGlobal))
		if rateVariation:
			meanRate, varianceRate = variance(siteRates)
			for i in range(lRef):
				if siteRates[i]<0.0:
					print("negative rate "+str(siteRates[i])+" "+str(i))
			print("Rate variation variance: "+str(varianceRate))
		if errorRateSiteSpecific:
			meanErrorRate, varianceErrorRate = variance(errorRates)
			for i in range(lRef):
				if errorRates[i]<0.0:
					print("negative error rate "+str(errorRates[i])+" "+str(i))
			print("Error rate variation, mean: "+str(meanErrorRate)+" , variance: "+str(varianceErrorRate))
		updateMutMatrices(mutMatrixGlobal,siteRates=siteRates)
		if usingErrorRate:
			updateErrorRates(errorRateGlobal,errorRates=errorRates)
		setAllDirty(t1)
		reCalculateAllGenomeLists(t1)
		improvement=traverseTreeToOptimizeBranchLengths(t1, fastPass=True)
		reCalculateAllGenomeLists(t1)
		if estimateErrorRate or estimateSiteSpecificErrorRate:
			oldLK=float("-inf")
			newLk=calculateTreeLikelihood(t1)
			print("EM estimation of error rates. Initial LK after first pass: "+str(newLk))
			numEMsteps=0
			while (newLk-oldLK>1.0) and numEMsteps<20:
				setAllDirty(t1)
				improvement=traverseTreeToOptimizeBranchLengths(t1)
				newLkBranch=calculateTreeLikelihood(t1)
				print("Updated "+str(improvement)+" branch lengths leading to LK "+str(newLkBranch))

				mutMatrixGlobal, siteRates, errorRateGlobal, errorRates = expectationMaximizationCalculationRates(t1)
				updateMutMatrices(mutMatrixGlobal,siteRates=siteRates)
				updateErrorRates(errorRateGlobal,errorRates=errorRates)
				reCalculateAllGenomeLists(t1)
				oldLK=newLk
				newLk=calculateTreeLikelihood(t1)
				print("New LK step "+str(numEMsteps)+": "+str(newLk))
				if rateVariation:
					meanRate, varianceRate = variance(siteRates)
					print("Rate variation variance: "+str(varianceRate))
				if errorRateSiteSpecific:
					meanErrorRate, varianceErrorRate = variance(errorRates)
					print("Error rate variation, mean: "+str(meanErrorRate)+" , variance: "+str(varianceErrorRate))
				print("Error rate: "+str(errorRateGlobal))
				numEMsteps+=1

	timeRecalculation=time()-start
	print("Time to run EM pre-topology estimation: "+str(timeRecalculation))

	print("Number of nodes: "+str(numNodes[0]))
	print("Os per node: "+str(float(numNodes[4])/numNodes[0]))
	print("Nucs per node: "+str(float(numNodes[1])/numNodes[0]))
	print("Ns per node: "+str(float(numNodes[3])/numNodes[0]))

	start=time()
	setAllDirty(t1)
	print("Post-model optimization branch length optimization, initial fast pass")
	improvement=traverseTreeToOptimizeBranchLengths(t1, fastPass=True)
	reCalculateAllGenomeLists(t1)
	print("Now proper branch length optimization, first pass")
	setAllDirty(t1)
	improvement=traverseTreeToOptimizeBranchLengths(t1)
	print("Branch length optimization round 1, number of changes: "+str(improvement))
	subRound=0
	while subRound<20:
		if (not improvement):
			break
		subRound+=1
		improvement=traverseTreeToOptimizeBranchLengths(t1)
		print("branch length finalization subround "+str(subRound+1)+", number of changes "+str(improvement))
	timeForBranchOptimization=(time()-start)
	print("Time for updating branch lengths: "+str(timeForBranchOptimization))

if (writeTreesToFileEveryTheseSteps>0):
	currentRoot=t1
	intermediateTreesFile.write("Topology 0\n")
	if binaryTree:
		intermediateTreesFile.write(createBinaryNewick(currentRoot)+"\n")
	else:
		intermediateTreesFile.write(createNewick(currentRoot)+"\n")
if (writeLKsToFileEveryTheseSteps>0):
	currentRoot=t1
	totalLK=calculateTreeLikelihood(currentRoot)
	intermediateLKsFile.write("Topology 0, LK: "+str(totalLK)+"\n")








internalNodeNamesGiven=False
giveInternalNodeNames(t1)
internalNodeNamesGiven=True

#if asked, first run a round of short range topology search
timeTopology=0.0
if fastTopologyInitialSearch and (inputTree=="" or largeUpdate or aBayesPlus):
	if aBayesPlus:
		aBayesPlusOn=True
		#if networkOutput and (not internalNodeNamesGiven):
		#	giveInternalNodeNames(t1)
		#	internalNodeNamesGiven=True
	#first update just branch lengths
	start=time()
	setAllDirty(t1)
	print("Preliminary branch length optimization")
	improvement=traverseTreeToOptimizeBranchLengths(t1)
	print("Branch length optimization round 1, number of changes: "+str(improvement))
	subRound=0
	while subRound<20:
		if (not improvement):
			break
		subRound+=1
		improvement=traverseTreeToOptimizeBranchLengths(t1)
		print("branch length finalization subround "+str(subRound+1)+" number of changes "+str(improvement))
	timeForBranchOptimization=(time()-start)
	print("Time for updating branch lengths: "+str(timeForBranchOptimization))

	#now update the topology with a shallow search
	start=time()
	setAllDirty(t1)
	newRoot,improvement=startTopologyUpdates(t1,checkEachSPR=debugging,strictTopologyStopRules=strictTopologyStopRulesInitial,allowedFailsTopology=allowedFailsTopologyInitial,thresholdLogLKtopology=thresholdLogLKtopologyInitial,thresholdTopologyPlacement=thresholdTopologyPlacementInitial)
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
		newRoot,improvement=startTopologyUpdates(t1,checkEachSPR=debugging,strictTopologyStopRules=strictTopologyStopRulesInitial,allowedFailsTopology=allowedFailsTopologyInitial,thresholdLogLKtopology=thresholdLogLKtopologyInitial,thresholdTopologyPlacement=thresholdTopologyPlacementInitial)
		if newRoot!=None:
			t1=newRoot
		print("LK improvement apparently brought: "+str(improvement))
		if improvement<thresholdLogLKwholeTopologyImprovement:
			break
		subRound+=1
	timeForUpdatingTopology=(time()-start)
	print("Time for the subrounds of this traversal of the tree: "+str(timeForUpdatingTopology))
	timeTopology+=timeForUpdatingTopology


	#if estimating error rates, repeating EM after shallow topological search
	if estimateErrorRate or estimateSiteSpecificErrorRate:
		oldLK=float("-inf")
		newLk=calculateTreeLikelihood(t1)
		print("New EM estimation of error rates after shallow topological search. Initial LK: "+str(newLk))
		numEMsteps=0
		while (newLk-oldLK>1.0) and numEMsteps<20:
			setAllDirty(t1)
			improvement=traverseTreeToOptimizeBranchLengths(t1)
			newLkBranch=calculateTreeLikelihood(t1)
			print("Updated "+str(improvement)+" branch lengths leading to LK "+str(newLkBranch))

			mutMatrixGlobal, siteRates, errorRateGlobal, errorRates = expectationMaximizationCalculationRates(t1)
			updateMutMatrices(mutMatrixGlobal,siteRates=siteRates)
			updateErrorRates(errorRateGlobal,errorRates=errorRates)
			reCalculateAllGenomeLists(t1)
			oldLK=newLk
			newLk=calculateTreeLikelihood(t1)
			print("New LK step "+str(numEMsteps)+": "+str(newLk))
			if rateVariation:
				meanRate, varianceRate = variance(siteRates)
				print("Rate variation variance: "+str(varianceRate))
			if errorRateSiteSpecific:
				meanErrorRate, varianceErrorRate = variance(errorRates)
				print("Error rate variation, mean: "+str(meanErrorRate)+" , variance: "+str(varianceErrorRate))
			print("Error rate: "+str(errorRateGlobal))
			numEMsteps+=1

#writing to output the substitution model (possibly with rate variation and error rates)
file=open(outputFile+"_subs.txt","w")
for i in range4:
	for j in range4:
		file.write(str(mutMatrixGlobal[i][j])+"\t")
	file.write("\n")
if rateVariation:
	file.write("\n\n"+"Site rates:\n")
	for i in range(lRef):
		file.write(str(i+1)+"\t"+str(siteRates[i])+"\n")
if estimateSiteSpecificErrorRate:
	file.write("\n\n"+"Site error rates:\n")
	for i in range(lRef):
		file.write(str(i+1)+"\t"+str(errorRates[i])+"\n")
elif estimateErrorRate:
	file.write("\n\n"+"Error rate: "+str(errorRateGlobal)+"\n")
file.close()
print("Missed minor samples: "+str(totalMissedMinors[0]))
print("Final Substitution matrix:")
print(mutMatrixGlobal)


#now run deep topological improvements
debugging=False
for i in range(numTopologyImprovements):
	if aBayesPlus:
		aBayesPlusOn=True
		if networkOutput and (not internalNodeNamesGiven):
			giveInternalNodeNames(t1)
			internalNodeNamesGiven=True
	print("Starting topological impromevement attempt traversing number "+str(i+1))
	start=time()
	if inputTree=="" or largeUpdate or aBayesPlus:
		setAllDirty(t1)
	newRoot,improvement=startTopologyUpdates(t1,checkEachSPR=debugging)
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
		newRoot,improvement=startTopologyUpdates(t1,checkEachSPR=debugging)
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
	improvement=traverseTreeToOptimizeBranchLengths(t1)
	print("Branch length optimization round 1, number of changes: "+str(improvement))
	subRound=0
	while subRound<20:
		if (not improvement):
			break
		subRound+=1
		improvement=traverseTreeToOptimizeBranchLengths(t1)
		print("branch length finalization subround "+str(subRound+1)+", number of changes "+str(improvement))
	timeForBranchOptimization=(time()-start)
	print("Time for updating branch lengths: "+str(timeForBranchOptimization))


#calculate total likelihood
if calculateLKfinalTree:
	totalLK=calculateTreeLikelihood(t1)
	print("totalLK: "+str(totalLK))
	file=open(outputFile+"_LK.txt","w")
	file.write(str(totalLK)+"\n")
	file.close()

#put sample names in the tree
def addSampleNamesToTree(inputAlignmentFile,treeRoot):
	#Read input file to collect all sample names
	fileI=open(inputAlignmentFile)
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
	nextLeaves=[treeRoot]
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

if estimateErrors:
	if (not debugging) and (writeTreesToFileEveryTheseSteps==0) and inputTree=="":
		addSampleNamesToTree(inputFile,t1)
	print("Estimating errors in input alignment")
	file=open(outputFile+"_estimatedErrors.txt","w")
	calculateErrorProbabilities(t1,file,minErrorProb)
	file.close()
	print("Errors estimated")

if estimateMAT:
	expectationMaximizationCalculationRates(t1,trackMutations=True)



#Free space by deleting the genome lists in the tree.
nextLeaves=[t1]
while nextLeaves:
	node=nextLeaves.pop()
	if node.up!=None:
		node.probVect=None
	if node.dist:
		if node.up!=None:
			node.probVectTot=None
			node.probVectTotUp=None
			node.furtherMidNodes=None
	if node.children:
		node.probVectUpRight=None
		node.probVectUpLeft=None
		for c in node.children:
			nextLeaves.append(c)
print("Deleted genome lists.")

#put sample names in the tree
if (not debugging) and (writeTreesToFileEveryTheseSteps==0) and inputTree=="" and (not estimateErrors):
	addSampleNamesToTree(inputFile,t1)




#generate the string corresponding to a line of the tsv file for use in Taxonium.
def tsvForNode(node,name,featureList,identicalTo=""):
	stringList=[name+"\t"]
	if identicalTo!="":
		stringList.append(identicalTo)
	stringList.append("\t")
	#print("tsvForNode "+name)
	for feat in featureList:
		if node!=None:
			if hasattr(node, feat):
				#print(feat)
				if feat=="support" or feat=="IQsupport":
					stringList.append(str(getattr(node, feat)))
				elif feat=="alternativePlacements":
					#print(node.alternativePlacements)
					for iNode in range(len(node.alternativePlacements)):
						stringList.append(node.alternativePlacements[iNode][0].name+":"+str(node.alternativePlacements[iNode][1]))
						if iNode<(len(node.alternativePlacements)-1):
							stringList.append(",")
					#print(stringList)
				elif feat=="mutations":
					for iNode in range(len(node.mutations)):
						mutation=node.mutations[iNode]
						stringList.append(allelesListExt[mutation[0]]+str(mutation[1])+allelesListExt[mutation[2]]+":"+str(mutation[3]))
						if iNode<(len(node.mutations)-1):
							stringList.append(",")
				elif feat=="Ns":
					for iNode in range(len(node.Ns)):
						mutation=node.Ns[iNode]
						if type(mutation)==int:
							stringList.append(str(mutation))
						else:
							stringList.append(str(mutation[0])+"-"+str(mutation[1]))
						if iNode<(len(node.Ns)-1):
							stringList.append(",")
				elif feat=="errors":
					for iNode in range(len(node.errors)):
						mutation=node.errors[iNode]
						stringList.append(allelesListExt[mutation[0]]+str(mutation[1])+allelesListExt[mutation[2]]+":"+str(mutation[3]))
						if iNode<(len(node.errors)-1):
							stringList.append(",")
				elif feat=="lineage":
					stringList.append(node.lineage)
				elif feat=="lineages":
					for lineageName in node.lineages.keys():
						stringList.append(lineageName+":"+str(node.lineages[lineageName]))
						stringList.append(",")
					stringList.pop()
		stringList.append("\t")
	stringList.append("\n")
	return "".join(stringList)

#write tsv file with metadata/annotations for tree nodes and taxa
def writeTSVfile(node,file):
	featureNames={}
	if keepInputIQtreeSupports:
		featureNames['IQsupport']='IQsupport'
	if aBayesPlusOn:
		featureNames['support']='support'
		if networkOutput:
			featureNames['alternativePlacements']='uncertainty'
	if estimateMAT:
		featureNames['mutations']='mutations'
		featureNames['Ns']='Ns'
	if usingErrorRate:
		featureNames['errors']='errors'
	if performLineageAssignment:
		featureNames['lineage']='lineage'
		featureNames['lineages']='lineages'
	featureList=list(featureNames.keys())
	file.write("node"+"\t"+"collapsedTo")
	for feat in featureList:
		file.write("\t"+featureNames[feat])
	file.write("\n")
	#now write to file the features for each node of the tree.
	nextNode=node
	direction=0
	numLeaves=0
	while nextNode!=None:
		if nextNode.children:
			if direction==0:
				nextNode=nextNode.children[0]
			elif direction==1:
				nextNode=nextNode.children[1]
				direction=0
			else:
				if aBayesPlusOn or estimateMAT or performLineageAssignment:
					file.write(tsvForNode(nextNode,nextNode.name,featureList))
				if nextNode.up!=None:
					if nextNode.up.children[0]==nextNode:
						direction=1
					else:
						direction=2
				nextNode=nextNode.up
		else:
			numLeaves+=(1+len(nextNode.minorSequences))
			file.write(tsvForNode(nextNode,nextNode.name,featureList))
			if len(nextNode.minorSequences)>0:
				for s2 in nextNode.minorSequences:
					if supportForIdenticalSequences or performLineageAssignment:
						file.write(tsvForNode(nextNode,s2,featureList))
					else:
						file.write(tsvForNode(None,s2,featureList,identicalTo=nextNode.name))
				if aBayesPlusOn or estimateMAT or performLineageAssignment:
					if supportForIdenticalSequences or performLineageAssignment:
						file.write(tsvForNode(nextNode,nextNode.name+"_MinorSeqsClade",featureList))
					else:
						file.write(tsvForNode(None,nextNode.name+"_MinorSeqsClade",featureList,identicalTo=nextNode.name))
			if nextNode.up!=None:
				if nextNode.up.children[0]==nextNode:
					direction=1
				else:
					direction=2
			nextNode=nextNode.up
	file.close()




if binaryTree:
	newickString=createBinaryNewick(t1)
else:
	newickString=createNewick(t1)

if aBayesPlus or estimateMAT:
	file=open(outputFile+"_nexusTree.tree","w")
	file.write("#NEXUS\nbegin taxa;\n	dimensions ntax="+str(countTips(t1))+";\n	taxlabels\n")
	writeTaxaNames(file,t1)
	file.write(";\nend;\n\nbegin trees;\n	tree TREE1 = [&R] ")
	file.write(newickString)
	file.write("\nend;\n")
	file.close()
	#TODO test and then copy this into new version
	file=open(outputFile+"_metaData.tsv","w")
	writeTSVfile(t1,file)
	file.close()
	aBayesPlusOn=False
	networkOutput=False
	estimateMAT=False
	if binaryTree:
		newickString=createBinaryNewick(t1)
	else:
		newickString=createNewick(t1)

if runOnlyExample:
	print("Final tree:")
	print(newickString)
	raise Exception("exit")

file=open(outputFile+"_tree.tree","w")
file.write(newickString)
file.close()
print("Time spent finding placement nodes: "+str(timeFinding))
print("Time spent placing samples on the tree: "+str(timePlacing))
print("Time spent in total updating the topology and branch lengths: "+str(timeTopology))
exit()





