import sys
from math import log
import argparse
from time import time
import os.path
from operator import itemgetter

#©EMBL-European Bioinformatics Institute, 2021-2023
#Developd by Nicola De Maio, with contributions from Myrthe Willemsen.

# MAPLE code to estimate a tree by maximum likelihood from a MAPLE format input.

#TODO test HnZ modifier to favour placement onto successful lineages (those with many samples and direct descendants)

#longer-term plans
#TODO parallelization of initial tree search, probably only feasible in C++
#TODO Timed trees
#TODO further detailed substitution models that has a specific rate for each position, from each nucleotide, to each nucleotide (12 times more parameters than position-specific model).
#TODO phylogeography
#TODO codon models, selection
#TODO indels, alignment via alignment-free pair-HMM
#TODO Bayesian implementation in BEAST
#TODO Recombination via pair-HMM based SPR search of second parent
#TODO phylodynamics/selection
#TODO relativistic phylogenetics for arbitrary large trees

#things one could MAYBE also try in the future:
#TODO implement parsimony?
#TODO create an initial stage for sample placement with very low thresholds; then, a second placement stage would use the default thresholds and work 
# like the SPR search but starting the placement search from the node found by the first stage; finally, the last stage would be as usual.
# I tried this but it comes at little computational advantage and at large accuracy cost.

parser = argparse.ArgumentParser(description='Estimate a tree from a diff format and using iterative approximate maximum likelihood sample placement.')
#important options
parser.add_argument('--input',default="/Users/demaio/Desktop/GISAID-hCoV-19-phylogeny-2021-03-12/phylogenetic_inference/2021-03-31_unmasked_differences_reduced.txt_consensus-based.txt", help='Input MAPLE file name; should contain first the reference genome and then the difference of all samples with respet to the reference.')
parser.add_argument('--reference',default="", help='Optional input reference file name. By default it assumes instead that the reference is part of the MAPLE format input.')
parser.add_argument("--model", help="Which substitution model should be used. Allowed models so far are JC, GTR (default) or UNREST.", default="GTR")
parser.add_argument('--output',default="/Users/demaio/Desktop/GISAID-hCoV-19-phylogeny-2021-03-12/phylogenetic_inference/MAPLE", help='Output path and identifier to be used for newick output file.')
parser.add_argument('--inputTree',default="", help='Input newick tree file name; this is optional, and is used for online inference (or for Robinson-Foulds distance calculation if option --inputRFtrees is also used).')
parser.add_argument('--inputRates',default="", help='Name of input containing pre-estimated substitution rates and possibly other model parameters; this is optional, and is only used for online inference. The same format as the MAPLE output is expected, and information about all parameters of the selected substitution model is expected.')
parser.add_argument("--largeUpdate", help="When using option --inputTree to do online inference, by default the tree is only updated locally where the new sequences are inserted. Use this option instead to perform a thorough update of the phylogeny.", action="store_true")
parser.add_argument('--inputRFtrees',default="", help='file name with input newick trees to be compared with the input tree in option --inputTree to calculate Robinson-Foulds distances; this option will turn off the normal MAPLE estimation - only RF distances will be calculated. Each newick tree contained in this file will be compared to the single tree specified with option --inputTree .')
parser.add_argument("--overwrite", help="Overwrite previous results if already present.", action="store_true")
parser.add_argument("--fast", help="Set parameters so to run tree inference faster; this will be less accurate in cases of high complexity, for example with recombination, sequencing errors, etc. It will overrule user choices for options --thresholdLogLK , --thresholdLogLKtopology , --allowedFails , --allowedFailsTopology .", action="store_true")
parser.add_argument("--rateVariation", help="Estimate and use rate variation: the model assumes one rate per site, and the rates are assumed independently (no rate categories). This might cause overfitting if the dataset is not large enough, but in any case one would probably only use MAPLE for large enough datasets.", action="store_true")
parser.add_argument("--estimateMAT", help="Estimate mutation events on the final tree, and write them on the nexus tree.", action="store_true")
parser.add_argument("--doNotImproveTopology", help="Do not perform SPR moves, despite searching for them; this is useful if one wants to analyse a tree and calculate branch supports without changing the input tree.", action="store_true")
parser.add_argument("--saveInitialTreeEvery",help="Every these many samples placed (default 50,000) save the current tree to file. This way if there is any problem, the building of the initial tree can be restarted from the last saved tree.",  type=int, default=50000)
parser.add_argument("--doNotPlaceNewSamples", help="Given an input tree, skip the placement of samples on the tree (in case the input alignment contained more samples than on the tree), so keep only the samples originally already on the tree.", action="store_true")
parser.add_argument("--doNotReroot", help="Skip rooting optimization.", action="store_true")
parser.add_argument("--noLocalRef", help="Do not use local references (this will usually take longer to run).", action="store_true")
#parallelization options
parser.add_argument("--numCores",help="Maximum number of cores to use if parallelizing. By default (0) use as many cores as available.",  type=int, default=1)
#threshold options
parser.add_argument("--minNumNon4",help="Minimum number of mutations threshold to define a subreference in the MAT - smaller values will cause denser references in the MAT.",  type=int, default=1)
parser.add_argument("--maxNumDescendantsForMATClade",help="Number of positive-branch-length descendants allowed before triggering mutation list (i.e. local reference) creation in the tree.",  type=int, default=50)
parser.add_argument("--noFastTopologyInitialSearch", help="Don't run a first fast short-range topology search (before the extensive one in case the latter is performed).", action="store_true")
parser.add_argument("--thresholdProb",help="relative probability threshold used to ignore possible states with very low probabilities.",  type=float, default=0.00000001)
parser.add_argument("--thresholdLogLK",help="logLK difference threshold (in number of mutations) to consider a logLk close to optimal.",  type=float, default=18.0)
parser.add_argument("--thresholdLogLKtopology",help="Maximum logLK difference threshold (in number of mutations) to consider a logLk close to optimal when looking for topology improvements.",  type=float, default=14.0)
parser.add_argument("--allowedFails",help="Number of times one can go down the tree without increasing placement likelihood before the tree traversal is stopped (only applies to non-0 branch lengths).",  type=int, default=5)
parser.add_argument("--allowedFailsTopology",help="Maximum number of times one can crawl along the tree without increasing placement likelihood before the tree traversal is stopped during topology search (only applies to non-0 branch lengths).",  type=int, default=4)
parser.add_argument("--numTopologyImprovements",help="Number of times we traverse the tree looking for deep (and slow) topological improvements. Default is 1, select 0 to skip deep topological search. values >1 are not recommended.",  type=int, default=1)
parser.add_argument("--thresholdTopologyPlacement",help="Don't try to re-place nodes that have current appending logLK cost above this threshold.",  type=float, default=-0.1)
parser.add_argument("--updateSubstMatrixEveryThisSamples",help="How many new samples to place before each update of the substitution rate matrix.",  type=int, default=25)
parser.add_argument("--nonStrictStopRules", help="If specified, then during the initial placement stage, a slower, non-strict rule for stopping the placement search is applied: the search is stopped if enough many consencutive LK worsening are observed, AND if LK is below the considered threshold.", action="store_true")
parser.add_argument("--strictTopologyStopRules", help="If specified, then during the topological improvement stage, a faster, strict rule for stopping the SPR search is applied: the search is stopped if enough many consencutive LK worsening are observed, OR if LK is below the considered threshold.", action="store_true")
parser.add_argument("--thresholdDiffForUpdate",help="Consider the probability of a new partial changed if the difference between old and new is above this threshold.",  type=float, default=0.00001)
parser.add_argument("--thresholdFoldChangeUpdate",help="Consider the probability of a new partial changed, if the fold difference between old and new if above this threshold.",  type=float, default=1.01)
parser.add_argument("--thresholdLogLKconsecutivePlacement",help="logLK difference threshold to consider something as a significant decrease in log-LK when considering consecutive likelihood decreases.",  type=float, default=0.01)
parser.add_argument("--thresholdLogLKTopologySubRoundImprovement",help="logLK difference threshold to consider something as a significant decrease in log-LK when considering sub-rounds of SPR moves.",  type=float, default=3.0)
parser.add_argument("--minBLenSensitivity",help="Fraction of a mutation to be considered as a precision for branch length estimation (default 0.001, which means branch lengths estimated up to a 1000th of a mutation precision).",  type=float, default=0.001)
parser.add_argument("--thresholdLogLKoptimization",help="logLK difference threshold (in number of mutations) to consider a logLk close to optimal when deciding for which possible placements to perform branch length optimization.",  type=float, default=1.0)
parser.add_argument("--thresholdLogLKoptimizationTopology",help="logLK difference threshold (in number of mutations) to consider a logLk close to optimal when deciding for which possible placements to perform branch length optimization during the SPR stage.",  type=float, default=1.0)
parser.add_argument("--maxReplacements",help="Maximum number of replacements attempts per node per SPR round (prevents loops).",  type=int, default=10)
parser.add_argument("--useFixedThresholdLogLKoptimizationTopology", help="Use this option if you want to specify the value in --thresholdLogLKoptimizationTopology instead of estimating it from the data.", action="store_true")
parser.add_argument("--minNumSamplesForRateVar",help="When creating the initial tree, start using the rate variation model only after these many samples have been added (better not to decrease this too much to avoid overfitting).",  type=int, default=510000)
parser.add_argument("--minNumSamplesForErrorModel",help="When creating the initial tree, start using the error model only after these many samples have been added (better not to decrease this too much to avoid overfitting).",  type=int, default=510000)
#lineage assignment options
parser.add_argument('--assignmentFileCSV',default="", help='give path and name to a .csv file containing reference sample names and their lineage assignment. When using this option, option --inputTree should also be used. Then MAPLE will assign all samples in the tree to a lineage following the reference lineages file. Each sample is assigned a lineage same as its closest reference parent.')
parser.add_argument('--assignmentFile',default="", help='Like --assignmentFileCSV but expects an alignment file in any format (as long as sequence names follow a > character as in Fasta and Maple formats).')
parser.add_argument('--inputNexusTree',default="", help='input nexus tree file name; this is optional, and is used for lineage assignment. The nexus tree is supposed to be the output of MAPLE, so that it has alternativePlacements annotation to represent topological uncertainty.')
parser.add_argument('--reRoot',default="", help='Re-root the input newick tree so that the specified sample/lineage is root. By default, no specified lineage/sample, so no rerooting.')
#rarer options
parser.add_argument("--defaultBLen",help="Default length of branches, for example when the input tree has no branch length information.",  type=float, default=0.000033)
parser.add_argument("--normalizeInputBLen",help="For the case the input tree has branch lengths expressed not in a likelihood fashion (that is, expected number of substitutions per site), then multiply the input branch lengths by this factor. Particularly useful when using parsimony-based input trees.",  type=float, default=1.0)
parser.add_argument("--multipleInputRFTrees", help="Use this option if the file specified by option --inputRFtrees contains multiple trees - otherwise only read the first tree from that file.", action="store_true")
parser.add_argument("--debugging", help="Test that likelihoods are calculated and updated as expected - time consuming and only meant for small trees for debugging purposes.", action="store_true")
parser.add_argument("--onlyNambiguities", help="Treat all ambiguities as N (total missing information).", action="store_true")
parser.add_argument("--nonBinaryTree", help="Write output tree with multifurcations - by default the tree is written as binary so to avoid problems reading the tree in other software.", action="store_true")
parser.add_argument("--writeTreesToFileEveryTheseSteps", help="By default, don't write intermediate trees to file. If however a positive integer is specified with this option, intermediate trees will be written to file every this many topological changes.",  type=int, default=0)
parser.add_argument("--writeLKsToFileEveryTheseSteps", help="By default, don't write likelihoods of intermediate trees to file. If however a positive integer is specified with this option, likelihoods of intermediate trees will be written to file every this many topological changes.",  type=int, default=0)
parser.add_argument("--noSubroundTrees", help="Do not write to file subround treees.", action="store_true")
#error model options
parser.add_argument("--estimateErrorRate", help="Estimate a single error rate for the whole genome. Input value is used as starting value", action="store_true")
parser.add_argument("--estimateSiteSpecificErrorRate", help="Estimate a separate error rate for each genome genome. Input value is used as starting value", action="store_true")
parser.add_argument("--errorRateInitial", help="Initial value used for estimating the error rate. The default is the inverse of the reference genome length (one error expected per genome).", type=float, default=0.0)
parser.add_argument("--errorRateFixed", help="Fix the error rate to a given input value", type=float, default=0.0)
parser.add_argument("--errorRateSiteSpecificFile", help="provide a file directory to a file that contains the siteSpecific error rates", type=str, default=None)
parser.add_argument("--estimateErrors", help="Estimate erroneous positions in the input sequences. This option is only allowed if option --estimateSiteSpecificErrorRate or errorRateSiteSpecificFile is also used.", action="store_true")
parser.add_argument("--minErrorProb", help="Minimum error probability to be reported when using option --estimateErrors", type=float, default=0.01)
#SPRTA options
parser.add_argument("--SPRTA", help="Calculate branch support values with a modification of the aBayes approach of Anisimova et al 2011 (De Maio et al 2024, in prep.).", action="store_true")
parser.add_argument("--aBayesPlus", help="Synonim for option --SPRTA.", action="store_true")
parser.add_argument("--networkOutput", help="Include in the output tree the alternative branches with their support (like in a weighted phylogenetic network).", action="store_true")
parser.add_argument("--minBranchSupport", help="Minimum branch support to be considered when using option --networkOutput", type=float, default=0.01)
parser.add_argument("--supportFor0Branches", help="When calculating branch support, also consider supports of branches of length 0, such as samples less informative than other samples which might have multiple placements on the tree.", action="store_true")
#parser.add_argument("--supportForIdenticalSequences", help="Also write support measures multiple times for identical sequences (useful when one wants supports written on every branch).", action="store_true")
parser.add_argument("--minMutProb", help="Minimum mutation probability to be written to output when using option --estimateMAT", type=float, default=0.01)
parser.add_argument("--keepInputIQtreeSupports", help="Assumes the input tree is from IQTREE, reads the support values on the tree branches, and prints the same values in the output nexus tree.", action="store_true")
#TODO HorseNotZebra modifiers option
parser.add_argument("--HnZ", help="By default (0), don't use modifiers. If 1, use the topological HorseNotZebra modifier; if 2, use the abundance HnZ modifier.",  type=int, default=0)
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
strictStopRules=(not args.nonStrictStopRules)
strictTopologyStopRules=args.strictTopologyStopRules
thresholdDiffForUpdate=args.thresholdDiffForUpdate
thresholdFoldChangeUpdate=args.thresholdFoldChangeUpdate
thresholdLogLKconsecutivePlacement=args.thresholdLogLKconsecutivePlacement
thresholdLogLKTopologySubRoundImprovement=args.thresholdLogLKTopologySubRoundImprovement
calculateLKfinalTree=True
maxNumDescendantsForMATClade=args.maxNumDescendantsForMATClade
minNumNon4=args.minNumNon4
saveInitialTreeEvery=args.saveInitialTreeEvery
doNotPlaceNewSamples=args.doNotPlaceNewSamples
doNotReroot=args.doNotReroot
noSubroundTrees=args.noSubroundTrees
useLocalReference= (not args.noLocalRef)

numCores=args.numCores
parallelize=False
if numCores>1:
	parallelize=True
	from multiprocessing import Pool
	from os import cpu_count
	numAvailCores=cpu_count()
	print('Number of CPUs available: '+str(numAvailCores))
	if numCores>numAvailCores:
		numCores=numAvailCores
	print('Number of CPUs that will be used: '+str(numCores))

thresholdLogLKoptimization=args.thresholdLogLKoptimization
thresholdLogLKoptimizationTopology=args.thresholdLogLKoptimizationTopology
minBLenSensitivity=args.minBLenSensitivity

example=False
runFast=args.fast
if runFast:
	thresholdLogLK=14.0
	allowedFails=4
	allowedFailsTopology=3
	thresholdLogLKtopology=7.0
	thresholdTopologyPlacement=-1.0
	minBLenSensitivity=0.001

fastTopologyInitialSearch=(not args.noFastTopologyInitialSearch)
strictTopologyStopRulesInitial=True
allowedFailsTopologyInitial=2
thresholdLogLKtopologyInitial=6.0
thresholdTopologyPlacementInitial=-0.1
minNumSamplesForRateVar=args.minNumSamplesForRateVar
minNumSamplesForErrorModel=args.minNumSamplesForErrorModel

defaultBLen=args.defaultBLen
normalizeInputBLen=args.normalizeInputBLen
multipleInputRFTrees=args.multipleInputRFTrees

inputTree=args.inputTree
inputRates=args.inputRates
inputRFtrees=args.inputRFtrees
largeUpdate=args.largeUpdate

assignmentFile=args.assignmentFile
assignmentFileCSV=args.assignmentFileCSV
inputNexusTree=args.inputNexusTree
maxReplacements=args.maxReplacements
writeTreesToFileEveryTheseSteps=args.writeTreesToFileEveryTheseSteps
writeLKsToFileEveryTheseSteps=args.writeLKsToFileEveryTheseSteps
reRoot=args.reRoot

mutMatrixGlobal=None
estimateErrorRate=args.estimateErrorRate
estimateSiteSpecificErrorRate = args.estimateSiteSpecificErrorRate
errorRateInitial = args.errorRateInitial
errorRateFixed = args.errorRateFixed
errorRateSiteSpecificFile = args.errorRateSiteSpecificFile
estimateErrors=args.estimateErrors
if estimateErrors and (not estimateSiteSpecificErrorRate) and (not errorRateSiteSpecificFile):
	print("Since option --estimateErrors has been selected, I am switching on --estimateSiteSpecificErrorRate so to allow estimation of site-specific error probabilities.")
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
if supportFor0Branches:
	supportForIdenticalSequences=True
else:
	supportForIdenticalSequences=False
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

#TODO HnZ check and HnZvector definition and update 
HnZ=args.HnZ
if HnZ>2 or HnZ<0:
	print("Option --HnZ only allows values 0 (no HnZ), 1 (topological HnZ) or 2 (abundance HnZ).")
	exit()
if HnZ==1:
	HnZvector=[0,0,0]
elif HnZ==2:
	HnZvector=[0,0,2*log(2)]
def updateHnZvector(n):
	currentN=len(HnZvector)
	while  currentN<=n:
		if HnZ==1:
			HnZvector.append(HnZvector[-1]+log( (2*currentN) - 3))
		elif HnZ==2:
			HnZvector.append(currentN*log(currentN))
		currentN+=1
def getHnZ(n):
	if n>=len(HnZvector):
		updateHnZvector(n)
	if n<=0:
		print("Requested HnZ score for non-positive nDesc0? "+str(n))
		raise Exception("exit")
	return HnZvector[n]

	

class Tree(object):
	def __init__(self):
		self.dist = []
		self.replacements=[]
		self.children = []
		self.mutations=[]
		self.up=[]
		self.dirty=[]
		self.name=[]
		self.minorSequences=[]
		self.probVect=[]
		self.probVectUpRight=[]
		self.probVectUpLeft=[]
		self.probVectTotUp=[]
		self.nDesc=[]
		#TODO number of branches descending from node after collapsing 0-length branches
		self.nDesc0=[]
	def __repr__(self):
		return "Tree object"
	def addNode(self,dirtiness=True):
		self.up.append(None)
		self.children.append([])
		self.dirty.append(dirtiness)
		self.name.append("")
		self.minorSequences.append([])
		self.mutations.append([])
		self.replacements.append(0)
		self.dist.append(0.0)
		self.probVect.append(None)
		self.probVectUpRight.append(None)
		self.probVectUpLeft.append(None)
		self.probVectTotUp.append(None)
		self.nDesc.append(0)
		#TODO make sure this si updated correctly throughout
		if HnZ:
			self.nDesc0.append(1)

#IMPORTANT definitions of genome list entries structure.
#The first element of a genome list entry represnts its type: 0="A", 1="C", 2="G", 3="T", 4="R", 5="N", 6="O"
# for entries of type 4="R" or 5="N" the second element represents the last position of the stretch of genome that they represent; 
# a value of 1 refers to the first position of the genome. In all other cases the second element represents the mutation-annotated tree reference nucleotide 
# at the considered position at the considered node.
#For entries of type 6="O" the last element is always a normalized vector of likelihoods (they always have to sum to 1.0).
# For elements of type <5, when considering error probabilities, a last additional flag element informs if the observation comes from a tip 
# without minor sequences (in which case we need to account for probabilities of errors), or not.
#If present, another element in third position represents the evolutionary distance since the considered type was observed (this is always missing with type "N"=5);
#If the element is not not present, this distance is assumed to be 0.0.
#If the distance (branch length) value is present, another distance value can also be present for entries of type <5 ; 
#this is to account for the fact that the observation might have occurred on the other side of the phylogeny with respect to the root; 
#when this additional branch length is also present, for example if the entry is (1, 234, 0.0001, 0.0002) then this means that a "C" was observed at genome position 234, and the observation
#is separated from the root by a distance of 0.0001, while a distance of 0.0002 separates the root from the current node (or position along a branch) considered.

# readNewick() now also reads IQTREE branch supports if keepInputIQtreeSupports==True
# different types of branch support in IQTREE:
# )63:
# )0.650000:
# )/0.125:
# )75.4:
# )/0.999:
# )75.4/67.3:

#function to read input newick string
def readNewick(nwFile,multipleTrees=False,dirtiness=True,createDict=False,inputDictNames=None,keepNames=False):
	phyloFile=open(nwFile)
	trees=[]
	line=phyloFile.readline()
	if inputDictNames==None and (not keepNames):
		namesInTree=[]
	if createDict:
		namesInTreeDict={}
	sampleNum=0
	while line!="":
		while line=="\n":
			line=phyloFile.readline()
		if line=="":
			break
		tree=Tree()
		tree.addNode(dirtiness=dirtiness)
		if keepInputIQtreeSupports:
			tree.IQsupport=[0.0]
		nwString=line.replace("\n","")
		index=0
		nodeIndex=len(tree.name)-1
		name=""
		distStr=""
		finished=False
		while index<len(nwString):
			if nwString[index]=="(":
				tree.children[nodeIndex].append(len(tree.up))
				tree.addNode(dirtiness=dirtiness)
				if keepInputIQtreeSupports:
					tree.IQsupport.append(None)
				tree.up[-1]=nodeIndex
				nodeIndex=len(tree.up)-1
				index+=1
			elif nwString[index]==";":
				trees.append((tree,nodeIndex))
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
					if keepNames:
						tree.name[nodeIndex]=name
					elif inputDictNames==None:
						tree.name[nodeIndex]=sampleNum
						if createDict:
							namesInTreeDict[name]=sampleNum
						sampleNum+=1
						namesInTree.append(name)
					else:
						tree.name[nodeIndex]=inputDictNames[name]
					name=""
				if distStr!="":
					tree.dist[nodeIndex]=float(distStr)*normalizeInputBLen
					if tree.dist[nodeIndex]<0.0:
						print("Warning: negative branch length in the input tree: "+distStr+" ; converting it to positive.")
						tree.dist[nodeIndex]=abs(tree.dist[nodeIndex])
					distStr=""
				else:
					tree.dist[nodeIndex]=defaultBLen
				nodeIndex=tree.up[nodeIndex]
				tree.children[nodeIndex].append(len(tree.up))
				tree.addNode(dirtiness=dirtiness)
				if keepInputIQtreeSupports:
					tree.IQsupport.append(None)
				tree.up[-1]=nodeIndex
				nodeIndex=len(tree.up)-1
				index+=1
			elif nwString[index]==")":
				if name!="":
					if keepNames:
						tree.name[nodeIndex]=name
					elif inputDictNames==None:
						tree.name[nodeIndex]=sampleNum
						if createDict:
							namesInTreeDict[name]=sampleNum
						sampleNum+=1
						namesInTree.append(name)
					else:
						tree.name[nodeIndex]=inputDictNames[name]
					name=""
				if distStr!="":
					tree.dist[nodeIndex]=float(distStr)*normalizeInputBLen
					distStr=""
				else:
					tree.dist[nodeIndex]=defaultBLen
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
						tree.IQsupport[tree.up[nodeIndex]]=suppVal
				else:
					index+=1
				nodeIndex=tree.up[nodeIndex]
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
	if keepNames:
		return trees
	elif createDict:
		return trees, namesInTree, namesInTreeDict
	elif inputDictNames==None:
		return trees, namesInTree
	else:
		return trees


#assign branch length to a node from a substring of a newick file
def assignNodeBLen(dist,nodeIndex,distStr):
	if distStr!="":
		dist[nodeIndex]=float(distStr)*normalizeInputBLen
		if dist[nodeIndex]<0.0:
			print("Warning: negative branch length in the input tree: "+distStr+" ; converting it to positive.")
			dist[nodeIndex]=abs(dist[nodeIndex])
		distStr=""
	else:
		dist[nodeIndex]=defaultBLen


#assign node features from a substring of a nexus file
def assignNodeFeatures(featureDicts,nodeIndex,annotationString):
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
	featureDicts[nodeIndex]=features


#function to read input nexus file
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
	tree=Tree()
	featureDicts=[]
	tree.addNode(dirtiness=dirtiness)
	featureDicts.append(None)
	nodeIndex=len(tree.name)-1
	name=""
	distStr=""
	annotationString=""
	finished=False
	madeUpNameCount=0
	while index<len(nwString):
		if nwString[index]=="(":
			tree.children[nodeIndex].append(len(tree.up))
			tree.addNode(dirtiness=dirtiness)
			featureDicts.append(None)
			tree.up[-1]=nodeIndex
			nodeIndex=len(tree.up)-1
			index+=1
		elif nwString[index]==";":
			if name!="":
				tree.name[nodeIndex]=name
				name=""
			else:
				madeUpNameCount+=1
				tree.name[nodeIndex]="madeUpNodeName"+str(madeUpNameCount)
			assignNodeBLen(tree.dist,nodeIndex,distStr)
			distStr=""
			assignNodeFeatures(featureDicts,nodeIndex,annotationString)
			annotationString=""
			root=nodeIndex
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
				tree.name[nodeIndex]=name
				name=""
			else:
				madeUpNameCount+=1
				tree.name[nodeIndex]="madeUpNodeName"+str(madeUpNameCount)
			assignNodeBLen(tree.dist,nodeIndex,distStr)
			distStr=""
			assignNodeFeatures(featureDicts,nodeIndex,annotationString)
			annotationString=""

			nodeIndex=tree.up[nodeIndex]
			tree.children[nodeIndex].append(len(tree.up))
			tree.addNode(dirtiness=dirtiness)
			featureDicts.append(None)
			tree.up[-1]=nodeIndex
			nodeIndex=len(tree.up)-1
			index+=1
		elif nwString[index]==")":
			if name!="":
				tree.name[nodeIndex]=name
				name=""
			else:
				madeUpNameCount+=1
				tree.name[nodeIndex]="madeUpNodeName"+str(madeUpNameCount)
			assignNodeBLen(tree.dist,nodeIndex,distStr)
			distStr=""
			assignNodeFeatures(featureDicts,nodeIndex,annotationString)
			annotationString=""
			index+=1
			nodeIndex=tree.up[nodeIndex]
		else:
			name+=nwString[index]
			index+=1
	if not finished:
		print("Error, final character ; not found in newick string in file "+nxFile+".")
		raise Exception("exit")

	phyloFile.close()
	tree.featureDicts=featureDicts
	print("Finished reading Nexus file.")
	return tree,root


#function that changes multifurcating tree structure into a binary tree by adding 0-length branches/nodes
def makeTreeBinary(tree,root):
	nodesToVisit=[root]
	while nodesToVisit:
		node=nodesToVisit.pop()
		if tree.children[node]:
			while len(tree.children[node])>2:
				child2=tree.children[node].pop()
				child1=tree.children[node].pop()
				tree.up[child1]=len(tree.up)
				tree.up[child2]=len(tree.up)
				tree.addNode()
				tree.children[-1].append(child1)
				tree.children[-1].append(child2)
				tree.up[-1]=node
				tree.children[node].append(len(tree.up)-1)
			nodesToVisit.append(tree.children[node][0])
			nodesToVisit.append(tree.children[node][1])


#given a list of mutation events, create a reverted list to be used for re-rooting (when the direction of branches changes).
def flipMutations(mutationList):
	newMutationList=[]
	for mut in mutationList:
		newMutationList.append((mut[0],mut[2],mut[1]))
	return newMutationList


#check inconsistent references between two vectors
def checkInconsistency(probVectP,probVectC):
	indexEntry1, indexEntry2, pos = 0, 0, 0
	entry1=probVectP[indexEntry1]
	entry2=probVectC[indexEntry2] 
	while True:
		if entry1[0]==4 or entry1[0]==5 : # case entry1 is R or N
			if entry2[0]==4 or entry2[0]==5:
				pos=min(entry1[1],entry2[1])
				if pos==lRef:
					break
				if entry2[1]==pos:
					indexEntry2+=1
					entry2=probVectC[indexEntry2]
			else: 
				pos+=1
				if pos==lRef:
					break
				indexEntry2+=1
				entry2=probVectC[indexEntry2]
			if entry1[1]==pos:
				indexEntry1+=1
				entry1=probVectP[indexEntry1]

		else: #entry1 is a non-ref nuc or O
			if entry2[0]!=4 and entry2[0]!=5:
				if entry1[1]!=entry2[1]:
					return True

			pos+=1
			if pos==lRef:
				break
			indexEntry1+=1
			entry1=probVectP[indexEntry1]
			if (entry2[0]!=4 and entry2[0]!=5) or entry2[1]==pos:
				indexEntry2+=1
				entry2=probVectC[indexEntry2]

	return False


#function to merge 2 MAT mutation lists (for example for when removing an internal node from the tree). First list is the top one.
#downward=True is for merging across different sides of the MRCA.
def mergeMutationLists(mutations1,mutations2,downward=False):
	ind1, ind2, pos1=0, 0, 0
	newMutations=[]
	while True:
		if ind1<len(mutations1):
			pos1=mutations1[ind1][0]
			if ind2<len(mutations2):
				pos2=mutations2[ind2][0]
				if pos1<pos2:
					if downward:
						newMutations.append((pos1,mutations1[ind1][2],mutations1[ind1][1]))
					else:
						newMutations.append(mutations1[ind1])
					ind1+=1
				elif pos2<pos1:
					newMutations.append(mutations2[ind2])
					ind2+=1
				else:
					if downward:
						sourceNuc=mutations1[ind1][2]
						endNuc=mutations1[ind1][1]
					else:
						sourceNuc=mutations1[ind1][1]
						endNuc=mutations1[ind1][2]
					if endNuc!=mutations2[ind2][1]:
						print("WARNING, inconsistent mutations in the MAT? "+str(ind1)+" "+str(ind2))
						print(mutations1)
						print(mutations2)
						print()
					if sourceNuc!=mutations2[ind2][2]:
						newMutations.append((pos2,sourceNuc,mutations2[ind2][2]))
					
					ind2+=1
					ind1+=1
			else:
				if downward:
					newMutations.append((pos1,mutations1[ind1][2],mutations1[ind1][1]))
				else:
					newMutations.append(mutations1[ind1])
				ind1+=1
		else:
			if ind2<len(mutations2):
				newMutations.append(mutations2[ind2])
				ind2+=1
			else:
				break
	return newMutations


#re-root a tree based on a given sample name, making the sample the root.
# if reRootAtInternalNode==True , flip also mutations in the local references. In this case sample represents the internal node whose parent branch to use as new root.
#TODO update also nDesc0 while re-rooting
def reRootTree(tree,root,sample,reRootAtInternalNode=False):
	sampleNode=None
	up=tree.up
	children=tree.children
	dist=tree.dist
	#TODO
	nDesc0=tree.nDesc0
	minorSequences=tree.minorSequences
	if reRootAtInternalNode:
		sampleNode=sample
		mutations=tree.mutations
		#mutations at the new root, representing differences wrt the reference
		rootMuts=mutations[root]
		childrenList=[up[sampleNode]]
		while up[childrenList[-1]]!=root:
			childrenList.append(up[childrenList[-1]])
		while childrenList:
			newNode=childrenList.pop()
			if mutations[newNode]:
				rootMuts=mergeMutationLists(rootMuts,mutations[newNode])
	else:
		#find the node associated to the sample
		nodesToVisit=[root]
		while nodesToVisit:
			node=nodesToVisit.pop()
			try:
				if tree.name[node]==sample:
					sampleNode=node
					break
			except:
				pass
			if children[node]:
				nodesToVisit.append(children[node][0])
				nodesToVisit.append(children[node][1])
	if sampleNode==None:
		print("Input lineage/sample for rerooting not found, please specify a lineage/sample contained in the tree with option --reRoot.")
		return root
	#now create new root
	if up[sampleNode]!=None:
		if up[up[sampleNode]]==None:
			if sampleNode==children[up[sampleNode]][0]:
				sibling=children[up[sampleNode]][1]
			else:
				sibling=children[up[sampleNode]][0]
			dist[sibling]+=dist[sampleNode]
			dist[sampleNode]=False
			#TODO
			if HnZ:
				nDesc0[up[sampleNode]]=nDesc0[sampleNode]
				if dist[sibling]>effectivelyNon0BLen:
					nDesc0[up[sampleNode]]+=1
				else:
					nDesc0[up[sampleNode]]+=nDesc0[sibling]
			return up[sampleNode]
		tree.addNode()
		children[-1].append(sampleNode)
		children[-1].append(up[sampleNode])
		newRoot=len(dist)-1
		oldDist=dist[sampleNode]
		oldDistUp=dist[up[sampleNode]]
		oldUp=up[sampleNode]
		oldUpUp=up[up[sampleNode]]
		dist[-1]=0.00000001
		if reRootAtInternalNode:
			dist[oldUp]=dist[sampleNode]/2
			dist[sampleNode]=dist[sampleNode]/2
		else:
			dist[sampleNode]=0.0
			dist[oldUp]=oldDist
		up[sampleNode]=newRoot
		up[oldUp]=newRoot
		
		currentNode=oldUpUp
		currentBLen=oldDistUp
		currentChild=oldUp
		currentChildChild=sampleNode
		if reRootAtInternalNode:
			oldMutations=mutations[currentChild]
			mutations[currentChild]=[]
		#now reiteratively flip the direction of branches until one reaches the old root
		while up[currentNode]!=None:
			if currentChildChild==children[currentChild][0]:
				numChildChild=0
			else:
				numChildChild=1
			children[currentChild][numChildChild]=currentNode
			if reRootAtInternalNode:
				newMutations=flipMutations(oldMutations)
				oldMutations=mutations[currentNode]
				mutations[currentNode]=newMutations
			oldBLen=dist[currentNode]
			oldPNode=up[currentNode]
			dist[currentNode]=currentBLen
			up[currentNode]=currentChild

			currentChildChild=currentChild
			currentChild=currentNode
			currentNode=oldPNode
			currentBLen=oldBLen
		#Now finalize by removing the old root
		if currentChildChild==children[currentChild][0]:
			numChildChild=0
		else:
			numChildChild=1
		if currentChild==children[currentNode][0]:
			numChild=0
		else:
			numChild=1
		if reRootAtInternalNode:
			#deal with mutations at the children of the old root
			newMutations=flipMutations(oldMutations)
			mutations[children[currentNode][1-numChild]]=mergeMutationLists(newMutations,mutations[children[currentNode][1-numChild]])
			mutations[newRoot]=rootMuts
		children[currentChild][numChildChild]=children[currentNode][1-numChild]
		up[children[currentNode][1-numChild]]=currentChild
		dist[children[currentNode][1-numChild]]+=currentBLen
		#TODO update nDesc0 of nodes afected by the re-rooting
		if HnZ:
			node0=currentChild
			while node0!=None:
				if children[node0]:
					if dist[children[node0][0]]>effectivelyNon0BLen:
						nDesc0[node0]=1
					else:
						nDesc0[node0]=nDesc0[children[node0][0]]
					if dist[children[node0][1]]>effectivelyNon0BLen:
						nDesc0[node0]+=1
					else:
						nDesc0[node0]+=nDesc0[children[node0][1]]
				else:
					nDesc0[node0]=1+len(minorSequences[node0])
				node0=up[node0]
					
		return newRoot
	else:
		return sampleNode


#Robinson-Foulds distance (1981) using a simplification of the algorithm from Day 1985.
#this function prepare the data to compare trees to a reference one t1.
#I split in two functions (the second one is RobinsonFouldsWithDay1985() below) so that I don't have to repeat these steps for the reference tree if I compare to the same reference tree multiple times.
def prepareTreeComparison(tree,t1,rooted=False,minimumBLen=0.000006):
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
	:param tree: input tree
	:param t1: input root
	:param rooted:
	:param minimumBLen:
	:param addRootRFL: Should the distance from the root node up be taken into account? I think this has usually no significance, and is always set to 1.0.
	:return: Various metrics, among which the RFL.
	"""
	children=tree.children
	up=tree.up
	dist=tree.dist
	name=tree.name
	exploredChildren = [0] * len(up)
	maxSoFar=[float("-inf")] * len(up)
	minSoFar=[float("inf")] * len(up) 
	nDescendants = [0] * len(up)
	tree.exploredChildren=exploredChildren
	tree.maxSoFar=maxSoFar
	tree.minSoFar=minSoFar
	tree.nDescendants=nDescendants
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
		while node!=up[t1]:
			if movingFrom==0: #0 means reaching node from parent, 1 means coming back from a child
				if len(children[node])==0:
					nLeaves+=1
					nextNode=up[node]
					movingFrom=1
					nodeTable.append([0,0])
				else:
					nextNode=children[node][0]
					movingFrom=0
			else:
				nChildren=len(children[node])
				exploredChildren[node]+=1
				if exploredChildren[node]==nChildren:
					nextNode=up[node]
					movingFrom=1
				else:
					nextNode=children[node][exploredChildren[node]]
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
	while node!=up[t1]:
		if movingFrom==0: #0 means reaching node from parent, 1 means coming back from a child
			if len(children[node])==0:
				name[node]=(name[node]).replace("?","_").replace("&","_")
				leafNameDict[name[node]]=leafCount
				leafNameDictReverse.append(name[node])
				if rooted:
					nodeTable.append([0,0])
				lastL=leafCount
				lastR=leafCount
				lastDesc=1
				leafCount+=1
				nextNode=up[node]
				movingFrom=1
				leafDistDict[name[node]] = dist[node]
			else:
				exploredChildren[node]=0
				nextNode=children[node][0]
				movingFrom=0
		else:
			nChildren=len(children[node])
			exploredChildren[node]+=1
			if lastL<minSoFar[node]:
				minSoFar[node]=lastL
			if lastR>maxSoFar[node]:
				maxSoFar[node]=lastR
			nDescendants[node]+=lastDesc
			if exploredChildren[node]==nChildren:
				nextNode=up[node]
				movingFrom=1
				lastL=minSoFar[node]
				lastR=maxSoFar[node]
				lastDesc=nDescendants[node]
				if node!=t1:
					sumBranchLengths += dist[node]
				if node==t1:
					nodeTable[lastR][0]=lastL
					nodeTable[lastR][1]=lastR
				else:
					if (not rooted) and (up[node]==t1) and (len(children[t1])==2):
						if node==children[t1][1]:
							currentBL=dist[node]+dist[children[t1][0]]
							addBranch=True
						else:
							addBranch=False
					else:
						currentBL=dist[node]
						addBranch=True
					if addBranch and currentBL>minimumBLen:
						numBranches+=1
						if rooted or lastL>0:
							if node==children[up[node]][-1]:
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
				nextNode=children[node][exploredChildren[node]]
				movingFrom=0
		node=nextNode
	return leafNameDict, nodeTable, leafCount, numBranches, leafDistDict, branchLengthDict, sumBranchLengths


#Robinson-Foulds distance (1981) using a simplification of the algorithm from Day 1985.
#this function compares the current tree t2 to a previous one for which prepareTreeComparison() was run.
# Example usage: leafNameDict, nodeTable, leafCount, numBranches = prepareTreeComparison(phyloTrue,rooted=False)
# numDiffs, normalisedRF, leafCount, foundBranches, missedBranches, notFoundBranches = RobinsonFouldsWithDay1985(phyloEstimated,leafNameDict,nodeTable,leafCount,numBranches,rooted=False)
def RobinsonFouldsWithDay1985(tree,t2,leafNameDict,nodeTable,leafCount,numBranches,leafDistDict, branchLengthDict, sumBranchLengths,rooted=False,minimumBLen=0.000006):
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
	children=tree.children
	up=tree.up
	dist=tree.dist
	name=tree.name
	exploredChildren = [0] * len(up)
	maxSoFar=[float("-inf")] * len(up)
	minSoFar=[float("inf")] * len(up) 
	nDescendants = [0] * len(up)
	tree.exploredChildren=exploredChildren
	tree.maxSoFar=maxSoFar
	tree.minSoFar=minSoFar
	tree.nDescendants=nDescendants
	RFL = sumBranchLengths
	KF = 0
	while node!=up[t2]:
		if movingFrom==0: #0 means reaching node from parent, 1 means coming back from a child
			if len(children[node])==0:
				name[node]=(name[node]).replace("?","_").replace("&","_")
				if name[node] in leafNameDict:
					leafNum=leafNameDict[name[node]]
				else:
					print(name[node]+" not in reference tree - aborting RF distance")
					return None, None, None, None, None, None
				lastL=leafNum
				lastR=leafNum
				lastDesc=1
				nextNode=up[node]
				movingFrom=1
				visitedLeaves+=1
				trueDist = leafDistDict[name[node]]
				KF += abs(trueDist- dist[node])   #As described, I have not added branch lengths of leaf nodes to sumBranchLengths, so no need for: RFL=- trueDist
			else:
				nextNode=children[node][0]
				movingFrom=0
		else:
			nChildren=len(children[node])
			exploredChildren[node]+=1
			if lastL<minSoFar[node]:
				minSoFar[node]=lastL
			if lastR>maxSoFar[node]:
				maxSoFar[node]=lastR
			nDescendants[node]+=lastDesc
			if exploredChildren[node]==nChildren:
				nextNode=up[node]
				movingFrom=1
				lastL=minSoFar[node]
				lastR=maxSoFar[node]
				lastDesc=nDescendants[node]
				if node!=t2:
					if (not rooted) and (up[node]==t2) and (len(children[t2])==2):
						if node==children[t2][1]:
							currentBL=dist[node]+dist[children[t2][0]]
							searchBranch=True
						else:
							searchBranch=False
					else:
						currentBL=dist[node]
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
				nextNode=children[node][exploredChildren[node]]
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


#generate the string corresponding to a node, taking into account support measures and other possible node features.
def stringForNode(tree,nextNode,nameNode,distB,estimateMAT=estimateMAT,networkOutput=networkOutput,aBayesPlusOn=aBayesPlusOn,namesInTree=None):
	children=tree.children
	up=tree.up
	name=tree.name
	aBayesPlusActive=False
	if aBayesPlusOn and hasattr(tree, 'alternativePlacements') and hasattr(tree, 'support'):
		aBayesPlusActive=True
		alternativePlacements=tree.alternativePlacements
		support=tree.support
		if hasattr(tree, 'rootSupport'):
			rootSupport=tree.rootSupport
		else:
			rootSupport=None
	estimateMATon=False
	errorsOn=False
	if estimateMAT and hasattr(tree, 'mutationsInf') and hasattr(tree, 'Ns'):
		estimateMATon=True
		mutationsInf=tree.mutationsInf
		Ns=tree.Ns
		if usingErrorRate and hasattr(tree, 'errors'):
			errorsOn=True
			errors=tree.errors
	if performLineageAssignment:
		lineage=tree.lineage
		lineages=tree.lineages
	printIQtreeSupportForNode=False
	if keepInputIQtreeSupports:
		if hasattr(tree, 'IQsupport'):
			printIQtreeSupportForNode=True
			IQsupport=tree.IQsupport
	strings=[]
	if aBayesPlusActive or estimateMATon or printIQtreeSupportForNode or printIQtreeSupportForNode:
		if (up[nextNode]!=None and ( (distB>effectivelyNon0BLen) or supportFor0Branches or errorsOn)) :
			#stringList.append("[&")
			if aBayesPlusActive and rootSupport!=None and rootSupport[nextNode]!=None:
				strings.append("rootSupport="+str(rootSupport[nextNode]))
				#stringList.append("rootSupport="+str(rootSupport[nextNode])+",")
			if aBayesPlusActive and ( (distB>effectivelyNon0BLen) or supportFor0Branches) and support[nextNode]!=None:
				#stringList.append("support="+str(support[nextNode]))
				strings.append("support="+str(support[nextNode]))
				if networkOutput and alternativePlacements[nextNode]:
					newString="alternativePlacements={"
					#stringList.append(",alternativePlacements={")
					for iNode in range(len(alternativePlacements[nextNode])):
						newString+=namesInTree[name[alternativePlacements[nextNode][iNode][0]]]+":"+str(alternativePlacements[nextNode][iNode][1])
						#stringList.append(namesInTree[name[alternativePlacements[nextNode][iNode][0]]]+":"+str(alternativePlacements[nextNode][iNode][1]))
						if iNode<(len(alternativePlacements[nextNode])-1):
							#stringList.append(",")
							newString+=","
					#stringList.append("}")
					newString+="}"
					strings.append(newString)
				#if estimateMATon and (distB or errorsOn):
				#	stringList.append(",")
			if estimateMATon and ( (distB) or errorsOn or (not children[nextNode])):
				if mutationsInf[nextNode]:
					newString="mutationsInf={"
					#stringList.append("mutationsInf={")
					for iNode in range(len(mutationsInf[nextNode])):
						mutation=mutationsInf[nextNode][iNode]
						#stringList.append(allelesListExt[mutation[0]]+str(mutation[1])+allelesListExt[mutation[2]]+":"+str(mutation[3]))
						newString+=allelesListExt[mutation[0]]+str(mutation[1])+allelesListExt[mutation[2]]+":"+str(mutation[3])
						if iNode<(len(mutationsInf[nextNode])-1):
							#stringList.append(",")
							newString+=","
					#stringList.append("},Ns={")
					newString+="}"
					strings.append(newString)
				if Ns[nextNode]:
					newString="Ns={"
					for iNode in range(len(Ns[nextNode])):
						mutation=Ns[nextNode][iNode]
						if type(mutation)==int:
							#stringList.append(str(mutation))
							newString+=str(mutation)
						else:
							#stringList.append(str(mutation[0])+"-"+str(mutation[1]))
							newString+=str(mutation[0])+"-"+str(mutation[1])
						if iNode<(len(Ns[nextNode])-1):
							#stringList.append(",")
							newString+=","
					#stringList.append("}")
					newString+="}"
					strings.append(newString)
				if errorsOn and (not children[nextNode]) and errors[nextNode]:
					#stringList.append(",errors={")
					newString="errors={"
					for iNode in range(len(errors[nextNode])):
						mutation=errors[nextNode][iNode]
						#stringList.append(allelesListExt[mutation[0]]+str(mutation[1])+allelesListExt[mutation[2]]+":"+str(mutation[3]))
						newString+=allelesListExt[mutation[0]]+str(mutation[1])+allelesListExt[mutation[2]]+":"+str(mutation[3])
						if iNode<(len(errors[nextNode])-1):
							#stringList.append(",")
							newString+=","
					#stringList.append("}")
					newString+="}"
					strings.append(newString)
			#if printIQtreeSupportForNode:
			#	stringList.append(",")
			#else:
			#	stringList.append("]")
		elif up[nextNode]==None and estimateMATon:
			#stringList.append("[&rootState={")
			newString="rootState={"
			currentPos=0
			firstDoneEntry=False
			rootVect=rootVector(tree.probVect[nextNode],False,(len(children[nextNode])==0 and len(tree.minorSequences[nextNode])==0),tree,nextNode )
			for entry in rootVect:
				if entry[0]!=4:
					if not firstDoneEntry:
						firstDoneEntry=True
					else:
						#stringList.append(",")
						newString+=","
				if entry[0]==5:
					#stringList.append("N"+str(currentPos+1)+"-"+str(entry[1]))
					newString+="N"+str(currentPos+1)+"-"+str(entry[1])
					currentPos=entry[1]
				elif entry[0]==6:
					vect=entry[-1]
					firstDone=False
					for i in range4:
						if vect[i]>minMutProb:
							if not firstDone:
								firstDone=True
							else:
								#stringList.append(",")
								newString+=","
							#stringList.append(allelesList[i]+str(currentPos+1)+":"+str(vect[i]))
							newString+=allelesList[i]+str(currentPos+1)+":"+str(vect[i])
					currentPos+=1
				elif entry[0]<4:
					#stringList.append(allelesList[entry[0]]+str(currentPos+1)+":1.0")
					newString+=allelesList[entry[0]]+str(currentPos+1)+":1.0"
					currentPos+=1
				else:
					currentPos=entry[1]
			#stringList.append("}")
			newString+="}"
			strings.append(newString)
			if aBayesPlusActive and rootSupport!=None and rootSupport[nextNode]!=None:
				#stringList.append(",rootSupport="+str(rootSupport[nextNode]))
				strings.append("rootSupport="+str(rootSupport[nextNode]))
			#if printIQtreeSupportForNode:
			#	stringList.append(",")
			#else:
			#	stringList.append("]")
		elif up[nextNode]==None and aBayesPlusActive and rootSupport!=None and rootSupport[nextNode]!=None:
			#stringList.append("[&rootSupport="+str(rootSupport[nextNode]))
			strings.append("rootSupport="+str(rootSupport[nextNode]))
			#if printIQtreeSupportForNode:
			#	stringList.append(",")
			#else:
			#	stringList.append("]")
		#elif printIQtreeSupportForNode:
		#	stringList.append("[&")
		if printIQtreeSupportForNode:
			#stringList.append("IQsupport="+str(IQsupport[nextNode]))
			strings.append("IQsupport="+str(IQsupport[nextNode]))
			#stringList.append("]")
	elif performLineageAssignment and (lineage[nextNode]!=None or lineages[nextNode]!=None):
		#stringList.append("[&")
		if lineage[nextNode]!=None:
			#stringList.append("lineage="+lineage[nextNode])
			strings.append("lineage="+lineage[nextNode])
			#if lineages[nextNode]!=None:
			#	stringList.append(",")
		if lineages[nextNode]!=None and lineages:
			#stringList.append("lineages={")
			newString="lineages={"
			for lineageName in lineages[nextNode].keys():
				#stringList.append(lineageName+":"+str(lineages[nextNode][lineageName]))
				newString+=lineageName+":"+str(lineages[nextNode][lineageName])
				#stringList.append(",")
				newString+=","
			#stringList.pop()
			newString=newString[:-1]
			#stringList.append("}")
			newString+="}"
			strings.append(newString)
		#stringList.append("]")
	#finalString="".join(stringList)
	#stringList=[]
	finalString=""
	if networkOutput or (not children[nextNode]):
		#stringList.append(nameNode)
		finalString=nameNode
	if strings:
		finalString+="[&"
		for s in range(len(strings)):
			finalString+=strings[s]
			if s<(len(strings)-1):
				finalString+=","
		finalString+="]"
	return finalString


#create newick string of a given tree (input node is assumed to be the root) - with option "binary" the generated tree is binary (polytomies are represented with branches of length 0).
def createNewick(tree,node,binary=True,namesInTree=None,includeMinorSeqs=True,estimateMAT=estimateMAT,networkOutput=networkOutput,aBayesPlusOn=aBayesPlusOn):
	nextNode=node
	stringList=[]
	direction=0
	numLeaves=0
	up=tree.up
	children=tree.children
	dist=tree.dist
	name=tree.name
	minorSequences=tree.minorSequences
	while nextNode!=None:
		if children[nextNode]:
			if direction==0:
				if dist[nextNode] or binary or up[nextNode]==None:
					stringList.append("(")
				nextNode=children[nextNode][0]
			elif direction==1:
				stringList.append(",")
				nextNode=children[nextNode][1]
				direction=0
			else:
				if dist[nextNode] or binary or up[nextNode]==None:
					if namesInTree==None:
						stringList.append(")"+name[nextNode])
					else:
						if name[nextNode]=="":
							stringList.append(")")
						else:
							stringList.append(")"+namesInTree[name[nextNode]])
					if aBayesPlusOn or estimateMAT or performLineageAssignment:
						stringList.append(stringForNode(tree,nextNode,"",dist[nextNode],estimateMAT=estimateMAT,networkOutput=networkOutput,aBayesPlusOn=aBayesPlusOn,namesInTree=namesInTree))
					if dist[nextNode]:
						stringList.append(":"+str(dist[nextNode]))
					else:
						stringList.append(":"+str(0.0))
				if up[nextNode]!=None:
					if children[up[nextNode]][0]==nextNode:
						direction=1
					else:
						direction=2
				nextNode=up[nextNode]
		else:
			numLeaves+=(1+len(minorSequences[nextNode]))
			if len(minorSequences[nextNode])>0 and includeMinorSeqs:
				if binary:
					for i in minorSequences[nextNode]:
						stringList.append("(")
					if supportForIdenticalSequences or performLineageAssignment:
						if namesInTree==None:
							stringList.append(stringForNode(tree,nextNode,name[nextNode],0.0,estimateMAT=estimateMAT,networkOutput=networkOutput,aBayesPlusOn=aBayesPlusOn,namesInTree=namesInTree))
						else:
							stringList.append(stringForNode(tree,nextNode,namesInTree[name[nextNode]],0.0,estimateMAT=estimateMAT,networkOutput=networkOutput,aBayesPlusOn=aBayesPlusOn,namesInTree=namesInTree))
					else:
						if namesInTree==None:
							stringList.append(name[nextNode])
						else:
							if name[nextNode]!="":
								stringList.append(namesInTree[name[nextNode]])
					stringList.append(":")
					for s2 in minorSequences[nextNode][:-1]:
						stringList.append("0.0,")
						if supportForIdenticalSequences or performLineageAssignment:
							if namesInTree==None:
								stringList.append(stringForNode(tree,nextNode,s2,0.0,estimateMAT=estimateMAT,networkOutput=networkOutput,aBayesPlusOn=aBayesPlusOn,namesInTree=namesInTree))
							else:
								stringList.append(stringForNode(tree,nextNode,namesInTree[s2],0.0,estimateMAT=estimateMAT,networkOutput=networkOutput,aBayesPlusOn=aBayesPlusOn,namesInTree=namesInTree))
						else:
							if namesInTree==None:
								stringList.append(s2)
							else:
								stringList.append(namesInTree[s2])
						stringList.append(":0.0):")
					stringList.append("0.0,")
					if supportForIdenticalSequences or performLineageAssignment:
						if namesInTree==None:
							stringList.append(stringForNode(tree,nextNode,minorSequences[nextNode][-1],0.0,estimateMAT=estimateMAT,networkOutput=networkOutput,aBayesPlusOn=aBayesPlusOn,namesInTree=namesInTree))
						else:
							stringList.append(stringForNode(tree,nextNode,namesInTree[minorSequences[nextNode][-1]],0.0,estimateMAT=estimateMAT,networkOutput=networkOutput,aBayesPlusOn=aBayesPlusOn,namesInTree=namesInTree))
					else:
						if namesInTree==None:
							stringList.append(minorSequences[nextNode][-1])
						else:
							stringList.append(namesInTree[minorSequences[nextNode][-1]])
					if namesInTree==None:
						stringList.append(":0.0)"+name[nextNode]+"_MinorSeqsClade")
					else:
						stringList.append(":0.0)"+namesInTree[name[nextNode]]+"_MinorSeqsClade")
				else:
					if dist[nextNode] or up[nextNode]==None:
						stringList.append("(")
					if supportForIdenticalSequences or performLineageAssignment:
						if namesInTree==None:
							stringList.append(stringForNode(tree,nextNode,name[nextNode],0.0,estimateMAT=estimateMAT,networkOutput=networkOutput,aBayesPlusOn=aBayesPlusOn,namesInTree=namesInTree))
						else:
							stringList.append(stringForNode(tree,nextNode,namesInTree[name[nextNode]],0.0,estimateMAT=estimateMAT,networkOutput=networkOutput,aBayesPlusOn=aBayesPlusOn,namesInTree=namesInTree))
					else:
						if namesInTree==None:
							stringList.append(name[nextNode])
						else:
							if name[nextNode]!="":
								stringList.append(namesInTree[name[nextNode]])
					stringList.append(":0.0")
					for s2 in minorSequences[nextNode]:
						stringList.append(",")
						if supportForIdenticalSequences or performLineageAssignment:
							if namesInTree==None:
								stringList.append(stringForNode(tree,nextNode,s2,0.0,estimateMAT=estimateMAT,networkOutput=networkOutput,aBayesPlusOn=aBayesPlusOn,namesInTree=namesInTree))
							else:
								stringList.append(stringForNode(tree,nextNode,namesInTree[s2],0.0,estimateMAT=estimateMAT,networkOutput=networkOutput,aBayesPlusOn=aBayesPlusOn,namesInTree=namesInTree))
						else:
							if namesInTree==None:
								stringList.append(s2)
							else:
								stringList.append(namesInTree[s2])
						stringList.append(":0.0")
					if dist[nextNode] or up[nextNode]==None:
						if namesInTree==None:
							stringList.append(")"+name[nextNode]+"_MinorSeqsClade")
						else:
							stringList.append(")"+namesInTree[name[nextNode]]+"_MinorSeqsClade")
			else:
				if namesInTree==None:
					stringList.append(name[nextNode])
				else:
					if name[nextNode]!="":
						stringList.append(namesInTree[name[nextNode]])
			if aBayesPlusOn or estimateMAT or performLineageAssignment:
				stringList.append(stringForNode(tree,nextNode,"",dist[nextNode],estimateMAT=estimateMAT,networkOutput=networkOutput,aBayesPlusOn=aBayesPlusOn,namesInTree=namesInTree))
			if dist[nextNode]:
				stringList.append(":"+str(dist[nextNode]))
			else:
				stringList.append(":"+str(0.0))
			if up[nextNode]!=None:
				if children[up[nextNode]][0]==nextNode:
					direction=1
				else:
					direction=2
			nextNode=up[nextNode]
	stringList.append(";")
	print("created newick string for tree with "+str(numLeaves)+" leaves.")
	return "".join(stringList)


#count how many samples are in the (sub)tree rooted at "node".
def countTips(tree,node):
	nextNode=node
	direction=0
	numSamples=0
	children=tree.children
	up=tree.up
	minorSequences=tree.minorSequences
	while nextNode!=None:
		if children[nextNode]:
			if direction==0:
				nextNode=children[nextNode][0]
			elif direction==1:
				nextNode=children[nextNode][1]
				direction=0
			else:
				if up[nextNode]!=None:
					if children[up[nextNode]][0]==nextNode:
						direction=1
					else:
						direction=2
				nextNode=up[nextNode]
		else:
			numSamples+=1+len(minorSequences[nextNode])
			if up[nextNode]!=None:
				if children[up[nextNode]][0]==nextNode:
					direction=1
				else:
					direction=2
			nextNode=up[nextNode]
	return numSamples


#write all taxa names to file - used for writing nexus files in BEAST fashion.
def writeTaxaNames(file,tree,node):
	nextNode=node
	direction=0
	numSamples=0
	children=tree.children
	up=tree.up
	name=tree.name
	minorSequences=tree.minorSequences
	while nextNode!=None:
		if children[nextNode]:
			if direction==0:
				nextNode=children[nextNode][0]
			elif direction==1:
				nextNode=children[nextNode][1]
				direction=0
			else:
				if up[nextNode]!=None:
					if children[up[nextNode]][0]==nextNode:
						direction=1
					else:
						direction=2
				nextNode=up[nextNode]
		else:
			file.write("	"+name[nextNode]+"\n")
			for samName in minorSequences[nextNode]:
				file.write("	"+samName+"\n")
			if up[nextNode]!=None:
				if children[up[nextNode]][0]==nextNode:
					direction=1
				else:
					direction=2
			nextNode=up[nextNode]
	return numSamples


#give names to internal nodes of the tree so to identify them in case one uses SPRTA etc.
def giveInternalNodeNames(tree,node,replaceNames=True,namesInTree=None):
	counter=1
	nextLeaves=[node]
	children=tree.children
	#alternativePlacements=[[]]*len(children)
	alternativePlacements=[]
	for i in range(len(children)):
		alternativePlacements.append([])
	tree.alternativePlacements=alternativePlacements
	name=tree.name
	while nextLeaves:
		nextNode=nextLeaves.pop()
		#list of node names where possible alternative placements of the current node are
		if children[nextNode]:
			if namesInTree!=None:
				if  (not replaceNames) and isinstance(name[nextNode], int):
					pass
				elif name[nextNode]!="" and (not replaceNames):
					namesInTree.append(name[nextNode])
					name[nextNode]=len(namesInTree)-1
				else:
					name[nextNode]=len(namesInTree)
					namesInTree.append("in"+str(len(namesInTree)))
			else:
				if replaceNames or (name[nextNode]==""):
					name[nextNode]="in"+str(counter)
				counter+=1
		for c in children[nextNode]:
			nextLeaves.append(c)


#TODO function to calculate number of branches under each multifurcation - to be used for HnZ modifier
#calculate and assign nDesc0 (number of effective branches under a multifurcation) for all nodes of a given tree
def calculateNDesc0(tree,root,checkExisting=False):
	children=tree.children
	minorSequences=tree.minorSequences
	up=tree.up
	dist=tree.dist
	nDesc0=tree.nDesc0
	node=root
	#direction 0 means from parent, direction 1 means from a child
	lastNode=None
	numProblems=0
	direction=0
	while node!=None:
		if direction==0:
			if children[node]:
				node=children[node][0]
			else:
				if minorSequences[node]:
					if checkExisting:
						if nDesc0[node]!=(1+len(minorSequences[node])):
							print("tip with minor sequences problem ")
							raise Exception("exit")
					nDesc0[node]=1+len(minorSequences[node])
				else:
					if checkExisting:
						if nDesc0[node]!=1:
							print("simple tip problem ")
							raise Exception("exit")
					nDesc0[node]=1
				lastNode=node
				node=up[node]
				direction=1
		else :
			if lastNode==children[node][0]:
				node=children[node][1]
				direction=0
			else:
				if checkExisting:
					oldNDesc0=nDesc0[node]
				nDesc0[node]=0
				for i in range(2):
					if dist[children[node][i]]:
						nDesc0[node]+=1
					else:
						nDesc0[node]+=(nDesc0[children[node][i]])
				if nDesc0[node]<=0:
					print("?")
					raise Exception("exit")
				if checkExisting:
					if oldNDesc0!=nDesc0[node]:
						numProblems+=1
						print("nDesc0 internal node problem ")
						print(nDesc0[children[node][0]])
						print(nDesc0[children[node][1]])
						print(dist[children[node][0]])
						print(dist[children[node][1]])
						print(oldNDesc0)
						print(nDesc0[node])
						print()
						raise Exception("exit")
				lastNode=node
				node=up[node]
				direction=1
	if checkExisting and numProblems:
		print("number of internal node problems: "+str(numProblems))

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
	trees, namesInTree, namesInTreeDict=readNewick(inputTree,createDict=True)
	tree1, rootIndex1=trees[0]
	print("Read input newick tree")
	leafNameDict, nodeTable, leafCount, numBranches, leafDistDict, branchLengthDict, sumBranchLengths = prepareTreeComparison(tree1,rootIndex1,rooted=False)
	otherTrees=readNewick(inputRFtrees,multipleTrees=multipleInputRFTrees,inputDictNames=namesInTreeDict)
	print("Read other input newick trees to be compared to the first one")
	file=open(outputFile+"_RFdistances.txt","w")
	file.write("RF\t"+"normalisedRF\t"+"leaves\t"+"foundBranches\t"+"missedBranches\t"+"notFoundBranches\t"+"RFL\n")
	for treePair in otherTrees:
		tree, rootIndex=treePair
		numDiffs, normalisedRF, leafCount, foundBranches, missedBranches, notFoundBranches, RFL = RobinsonFouldsWithDay1985(tree,rootIndex,leafNameDict, nodeTable, leafCount, numBranches,leafDistDict, branchLengthDict, sumBranchLengths,rooted=False)
		file.write(str(numDiffs)+"\t"+str(normalisedRF)+"\t"+str(leafCount)+"\t"+str(foundBranches)+"\t"+str(missedBranches)+"\t"+str(notFoundBranches)+"\t"+str(RFL)+"\n")
	print("Comparison ended")
	file.close()
	raise Exception("exit")


#run lineage assignment of samples in input tree following the input references
#beacuse this is run on an input tree, we don't need to consider minorSequences.
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
		tree1,rootIndex1=readNexus(inputNexusTree)
	else:
		tree1,rootIndex1=readNewick(inputTree,keepNames=True)[0]
		print("Input tree read")
		if reRoot!="":
			rootIndex1=reRootTree(tree1,rootIndex1,reRoot)

	giveInternalNodeNames(tree1,rootIndex1,replaceNames=False)
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

	node=rootIndex1
	direction=0 #0=from parent, 1+ from children
	lineage=""
	mostAncestralLineages=[]
	allLineages=[]
	if os.path.isfile(inputNexusTree):
		uncertaintyFlag=True
		nodeDict={}
	else:
		giveInternalNodeNames(tree1,rootIndex1)
		uncertaintyFlag=False
	#for each internal node with positive bLen, check all descendants with a distance of 0 from it; if any of them is a reference,
	#then assign that lineage to the node and all its descendants unless the assignment is changed downstream.
	#If no reference at distance 0 is found, then use parent assignment.
	children=tree1.children
	dist=tree1.dist
	up=tree1.up
	name=tree1.name
	tree1.lineage=[None]*len(up)
	lineageList=tree1.lineage
	tree1.mostAncestralLineages=[None]*len(up)
	mostAncestralLineagesList=tree1.mostAncestralLineages
	tree1.allLineages=[None]*len(up)
	allLineagesList=tree1.allLineages
	tree1.lineages=[None]*len(up)
	lineagesList=tree1.lineages
	while node!=None:
		#case of an internal node
		if children[node]:
			#coming from parent node
			if direction==0:
				if dist[node]:
					#first find 0-distance nodes and collect any reference lineage among them
					mostAncestralLineages2=[]
					allLineages2=[]
					nodesToVisit=list(children[node])
					while nodesToVisit:
						nextNode=nodesToVisit.pop()
						if not dist[nextNode]:
							if children[nextNode]:
								for child in children[nextNode]:
									nodesToVisit.append(child)
							elif name[nextNode] in references:
								lineage=references[name[nextNode]]
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
				lineageList[node]=lineage
				mostAncestralLineagesList[node]=mostAncestralLineages
				allLineagesList[node]=allLineages
				if uncertaintyFlag:
					nodeDict[name[node]]=node
				node=children[node][0]
			else:
				if direction==len(children[node]):
					if up[node]!=None:
						direction=-1
						for childNum in range(len(children[up[node]])):
							if node==children[up[node]][childNum]:
								direction=childNum+1
						if direction==-1:
							print("Error, not found child node in parent node list!")
							raise Exception("exit")
					node=up[node]
				else:
					lineage=lineageList[node]
					mostAncestralLineages=mostAncestralLineagesList[node]
					allLineages=allLineagesList[node]
					node=children[node][direction]
					direction=0
		
		#case of a sample
		else:
			if uncertaintyFlag:
				nodeDict[name[node]]=node
				if name[node] in references:
					lineageList[node]=references[name[node]]
					if dist[node]:
						mostAncestralLineagesList[node]=[lineageList[node]]
						allLineagesList[node]=[lineageList[node]]
					else:
						mostAncestralLineagesList[node]=mostAncestralLineages
						allLineagesList[node]=allLineages
				else:
					lineageList[node]=lineage
					mostAncestralLineagesList[node]=mostAncestralLineages
					allLineagesList[node]=allLineages
			else:
				if name[node] in references:
					file.write(name[node]+","+references[name[node]]+"\n")
				else:
					file.write(name[node]+","+lineage+"\n")
			if up[node]!=None:
				direction=-1
				for childNum in range(len(children[up[node]])):
					if node==children[up[node]][childNum]:
						direction=childNum+1
				if direction==-1:
					print("Error2, not found child node in parent node list!")
					raise Exception("exit")
			node=up[node]

	print("Finished tree pass for lineage assignment")
	#If nexus tree is input, do one more round using alternative placements to inform lineage probabilities.
	if uncertaintyFlag:
		if hasattr(tree1,"features") and ("support" in tree1.features):
			features=tree1.features
			support=features["support"]
			if "alternativePlacements" in features:
				alternativePlacements=features["alternativePlacements"]
			else:
				alternativePlacements=False
		else:
			features=False
			support=False
			alternativePlacements=False
		node=rootIndex1
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
			if children[node]:
				#coming from parent node
				if direction==0:
					lineages={}
					if support:
						for lin in allLineagesList[node]:
							lineages[lin]=support[node]/len(allLineagesList[node])
						if alternativePlacements:
							for alterPl in alternativePlacements[node].keys():
								alterNode=nodeDict[alterPl]
								alterLins=allLineagesList[alterNode]
								alterProb=alternativePlacements[node][alterPl]/len(alterLins)
								for alterLin in alterLins:
									if alterLin in lineages:
										lineages[alterLin]+=alterProb
									else:
										lineages[alterLin]=alterProb
					else:
						for lin in allLineagesList[node]:
							lineages[lin]=1.0/len(allLineagesList[node])
					lineagesList[node]=lineages
					node=children[node][0]
				else:
					if direction==len(children[node]):
						if up[node]!=None:
							direction=-1
							for childNum in range(len(children[up[node]])):
								if node==children[up[node]][childNum]:
									direction=childNum+1
							if direction==-1:
								print("Error, not found child node in parent node list!")
								raise Exception("exit")
						node=up[node]
					else:
						lineage=lineageList[node]
						node=children[node][direction]
						direction=0
			
			#case of a sample
			else:
				lineages={}
				if name[node] in references:
					file.write(name[node]+","+references[name[node]]+":1.0\n")
					lineages[references[name[node]]]=1.0
					numRefsTest+=1
				else:
					if support:
						for lin in allLineagesList[node]:
							lineages[lin]=support[node]/len(allLineagesList[node])
						numSupportTest+=1
						if alternativePlacements:
							for alterPl in alternativePlacements[node].keys():
								alterNode=nodeDict[alterPl]
								alterLins=allLineagesList[alterNode]
								alterProb=alternativePlacements[node][alterPl]/len(allLineagesList[alterNode])
								for alterLin in alterLins:
									if alterLin in lineages:
										lineages[alterLin]+=alterProb
									else:
										lineages[alterLin]=alterProb
						file.write(name[node])
						for alterPl in lineages.keys():
							file.write(","+alterPl+":"+str(lineages[alterPl]))
						file.write("\n")
					else:
						numNoSupportTest+=1
						for lin in allLineagesList[node]:
							lineages[lin]=1.0/len(allLineagesList[node])
						file.write(name[node])
						for alterPl in lineages.keys():
							file.write(","+alterPl+":"+str(lineages[alterPl]))
						file.write("\n")
				lineagesList[node]=lineages
				if up[node]!=None:
					direction=-1
					for childNum in range(len(children[up[node]])):
						if node==children[up[node]][childNum]:
							direction=childNum+1
					if direction==-1:
						print("Error2, not found child node in parent node list!")
						raise Exception("exit")
				node=up[node]
		print("Finished second tree pass for lineage assignment with uncertainty")
	print("Lineage assignment completed")
	file.close()
	newickString=createNewick(tree1,rootIndex1,binary=binaryTree,namesInTree=None)
	file=open(outputFile+"_nexusTree.tree","w")
	file.write("#NEXUS\nbegin taxa;\n	dimensions ntax="+str(countTips(tree1,rootIndex1))+";\n	taxlabels\n")
	writeTaxaNames(file,tree1,rootIndex1)
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
	trees, namesInTree, namesInTreeDict=readNewick(inputTree,dirtiness=largeUpdate,createDict=True)
	tree1,rootIndex1=trees[0]
	print("Read input newick tree")
	makeTreeBinary(tree1,rootIndex1)
	#TODO
	if HnZ:
		calculateNDesc0(tree1,rootIndex1)
else:
	namesInTree=[]
	namesInTreeDict={}

alleles={"A":0,"C":1,"G":2,"T":3}
allelesList=["A","C","G","T"]
allelesListExt=["A","C","G","T","?"]
allelesLow={"a":0,"c":1,"g":2,"t":3}
allelesUpOrLow={"a":0,"c":1,"g":2,"t":3,"A":0,"C":1,"G":2,"T":3}
allelesListLow=["a","c","g","t"]
ambiguities={"y":[0.0,0.5,0.0,0.5],"r":[0.5,0.0,0.5,0.0],"w":[0.5,0.0,0.0,0.5],"s":[0.0,0.5,0.5,0.0],"k":[0.0,0.0,0.5,0.5],"m":[0.5,0.5,0.0,0.0],"d":[1.0/3,0.0,1.0/3,1.0/3],"v":[1.0/3,1.0/3,1.0/3,0.0],"h":[1.0/3,1.0/3,0.0,1.0/3],"b":[0.0,1.0/3,1.0/3,1.0/3]}


#collect reference
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
def readConciseAlignment(fileName,extractReference=True,ref=""): #extractNames=False
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
	data={}
	while line!="" and line!="\n":
		seqList=[]
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
		data[name]=seqList
		nSeqs+=1
	fileI.close()
	print(str(nSeqs)+" sequences in diff file.")
	if extractReference:
		return ref, data
	else:
		return data


if refFile=="":
	ref, data=readConciseAlignment(inputFile) #extractNames=extractNamesFlag
else:
	ref=collectReference(refFile)
	data=readConciseAlignment(inputFile, extractReference=False, ref=ref) #,extractNames=extractNamesFlag
		
lRef=len(ref)
globalTotRate=-float(lRef)
print("Length of reference genome: "+str(lRef))
logLRef=log(lRef)
thresholdLogLKoptimizationTopology*=logLRef
thresholdLogLKoptimization*=logLRef
thresholdLogLKtopology*=logLRef
thresholdLogLK*=logLRef
thresholdLogLKtopologyInitial*=logLRef
effectivelyNon0BLen=1.0/(10*lRef)

#vector to count how many bases of each type are cumulatively in the reference genome up to a certain position
cumulativeBases=[[0,0,0,0]]
for i in range(lRef):
	cumulativeBases.append(list(cumulativeBases[i]))
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


#Shorten genome list by merging together R entries that are mergeable (that have the same details and are neighbours).
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


# Create a modified genome list, or modify an existing one, by taking into consideration the direction of the move, and the mutations on the MAT on a specific branch to traverse
def passGenomeListThroughBranch(probVect, mutations, dirIsUp=False): #, modifyCurrentList=False
	lMut=len(mutations)
	indexMAT, indexEntry, lastPos = 0, 0, 0
	#if not modifyCurrentList:
	newProbVect=[]
	entry=probVect[indexEntry]
	newEl, nucToPass = -1, -1
	while True:
		if entry[0]==5 :
			#if not modifyCurrentList:
			newProbVect.append(entry)
			lastPos=entry[1]
			if lastPos==lRef:
				break
			while indexMAT<lMut and mutations[indexMAT][0]<=lastPos:
				indexMAT+=1

			indexEntry+=1
			entry=probVect[indexEntry]

		elif entry[0]<4:
			lastPos+=1
			if indexMAT<lMut and mutations[indexMAT][0]<=lastPos:
				if dirIsUp:
					if entry[0]==mutations[indexMAT][1]:
						nucToPass=4
						newEl=lastPos
					else:
						nucToPass=entry[0]
						newEl=mutations[indexMAT][1]

				else:
					if entry[0]==mutations[indexMAT][2]:
						nucToPass=4
						newEl=lastPos
					else:
						nucToPass=entry[0]
						newEl=mutations[indexMAT][2]
				
				indexMAT+=1

				if len(entry)==2:
					newEntry=(nucToPass,newEl)
				elif len(entry)==3:
					newEntry=(nucToPass,newEl,entry[2])
				elif len(entry)==4:
					newEntry=(nucToPass,newEl,entry[2],entry[3])
				else:
					newEntry=(nucToPass,newEl,entry[2],entry[3],entry[4])

				newProbVect.append(newEntry)
					
			else:
				newProbVect.append(entry)

			if lastPos==lRef:
				break
			indexEntry+=1
			entry=probVect[indexEntry]

		elif entry[0]==4:
			while indexMAT<lMut and mutations[indexMAT][0]<=entry[1]:
				if mutations[indexMAT][0]>lastPos+1:
					lastPos=mutations[indexMAT][0]-1

					if len(entry)==2:
						newEntry=(4,lastPos)
					elif len(entry)==3:
						newEntry=(4,lastPos,entry[2])
					elif len(entry)==4:
						newEntry=(4,lastPos,entry[2],entry[3])
					else:
						newEntry=(4,lastPos,entry[2],entry[3],entry[4])

					newProbVect.append(newEntry)

				lastPos+=1
				if dirIsUp:
					nucToPass=mutations[indexMAT][2]
					newEl=mutations[indexMAT][1]
				else:
					nucToPass=mutations[indexMAT][1]
					newEl=mutations[indexMAT][2]
				
				indexMAT+=1

				if len(entry)==2:
					newEntry=(nucToPass,newEl)
				elif len(entry)==3:
					newEntry=(nucToPass,newEl,entry[2])
				elif len(entry)==4:
					newEntry=(nucToPass,newEl,entry[2],entry[3])
				else:
					newEntry=(nucToPass,newEl,entry[2],entry[3],entry[4])

				newProbVect.append(newEntry)
					
			if lastPos<entry[1]:
				lastPos=entry[1]
				newProbVect.append(entry)
			if lastPos==lRef:
				break
			indexEntry+=1
			entry=probVect[indexEntry]
				
		else: #entry is "O"
			lastPos+=1
			if indexMAT<lMut and mutations[indexMAT][0]<=lastPos:
				if dirIsUp:
					newEl=mutations[indexMAT][1]
				else:
					newEl=mutations[indexMAT][2]
				indexMAT+=1

				if len(entry)==3:
					newEntry=(6,newEl,entry[2])
				elif len(entry)==4:
					newEntry=(6,newEl,entry[2],entry[3])

				newProbVect.append(newEntry)
			else:
				newProbVect.append(entry)

			if lastPos==lRef:
				break
			indexEntry+=1
			entry=probVect[indexEntry]

	return newProbVect


#define partial likelihood vector for a sample given its data - now extended to account also for sequence errors
def probVectTerminalNode(diffs,tree,node):
	if node==None:
		numMinSeqs=0
	else:
		numMinSeqs=len(tree.minorSequences[node])
	if not errorRateSiteSpecific: 
		errorRate=errorRateGlobal
	if diffs is None:
		return [(5,lRef)]
	pos=1
	probVect=[]
	for m in diffs:
		currPos=m[1]
		if currPos>pos: #region where the node with branch length bLen is identical to the ref.
			probVect.append((4,currPos-1))
			pos=currPos
		if m[0]=="n" or m[0]=="-": #region with no info, store last position and length.
			if len(m)>2:
				length=m[2]
			else:
				length=1
			entry=(5,currPos+length-1)
			pos=currPos+length
		elif m[0] in allelesLow:
			#position at which node allele is sure but is different from the reference.
			if allelesLow[m[0]]==refIndeces[currPos-1]:
				print("Warning: alternative nucleotide entry coincides with reference: is the wrong reference being used?")
				print(m)
				print()
				entry=(4,currPos)
			else:
				entry=(allelesLow[m[0]],refIndeces[currPos-1])
			pos=currPos+1
		else:
			# non-"n" ambiguity character; for now interpret this as ambiguity instead of as a polymorphism.
			if onlyNambiguities:
				# if user asks to, to make things easier, interpret any ambiguity as an "n".
				entry=(5,currPos)
			else:
				#otherwise, store as an "other" scenario, where each nucleotide has its own partial likelihood.
				if usingErrorRate and (numMinSeqs==0):
					ambigVect=list(ambiguities[m[0]])
					sumUnnormalizedVector =  sum([bool(item) for item in ambigVect])
					if errorRateSiteSpecific: errorRate = errorRates[currPos-1]
					if sumUnnormalizedVector ==2:
						for i in range4:             # M, instead of [0.5, 0.5, 0, 0] we will now get [0.5- ⅓ε, 0.5- ⅓ε,  ⅓ε,  ⅓ε]
							if ambigVect[i]==0:
								ambigVect[i] = errorRate*0.33333
							else: #if entry.probs[i]==0.5:
								ambigVect[i] -= errorRate*0.33333
					elif sumUnnormalizedVector == 3:
						for i in range4:             # for V instead of [⅓, ⅓, ⅓, 0] we get, [ ⅓ - ε/9,  ⅓ -ε/9,  ⅓ -ε/9,  ⅓ε]
							if ambigVect[i] == 0:
								ambigVect[i] = errorRate*0.33333
							else:  # if entry.probs[i]==⅓:
								ambigVect[i] -= errorRate/9 # ⅓ - ε/9
					entry=(6,refIndeces[currPos-1],ambigVect)
				else:
					entry=(6,refIndeces[currPos-1],ambiguities[m[0]])
			pos=currPos+1
		probVect.append(entry)
	if pos<=lRef:
		probVect.append((4,lRef))

	if node!=None:
		#making sure that the genome list is relative to the current position of the node in the MAT
		up=tree.up
		mutations=tree.mutations
		listNodes=[node]
		nextNode=node
		while up[nextNode]!=None:
			nextNode=up[nextNode]
			listNodes.append(nextNode)
		while listNodes:
			nextNode=listNodes.pop()
			if mutations[nextNode]:
				probVect=passGenomeListThroughBranch(probVect,mutations[nextNode])#,modifyCurrentList=True

		shorten(probVect)

	return probVect


#update the genome list of terminal nodes (useful in case one wants to model errors and the error rate has changed)
def updateProbVectTerminalNode(probVect,numMinSeqs):
	if not errorRateSiteSpecific: 
		errorRate=errorRateGlobal
	if probVect is None:
		return None
	currentPos=0
	for m in probVect:
		if m[0]==6 :
			probs=m[2]
			sumUnnormalizedVector=0
			for i in range4:
				if probs[i]>0.2:
					sumUnnormalizedVector+=1
			if errorRateSiteSpecific: errorRate = errorRates[currentPos]
			if sumUnnormalizedVector ==2:
				for i in range4:             # M, instead of [0.5, 0.5, 0, 0] we will now get [0.5- ⅓ε, 0.5- ⅓ε,  ⅓ε,  ⅓ε]
					if probs[i]<0.2:
						if numMinSeqs:
							probs[i]=0.0
						else:
							probs[i] = errorRate*0.33333
					else: #if entry.probs[i]==0.5:
						if numMinSeqs:
							probs[i] = 0.5
						else:
							probs[i] = 0.5 - errorRate*0.33333
			elif sumUnnormalizedVector == 3:
				for i in range4:             # for V instead of [⅓, ⅓, ⅓, 0] we get, [ ⅓ - ε/9,  ⅓ -ε/9,  ⅓ -ε/9,  ⅓ε]
					if probs[i]<0.2:
						if numMinSeqs:
							probs[i] = 0.0
						else:
							probs[i] = errorRate*0.33333
					else:  # if entry.probs[i]==⅓:
						if numMinSeqs:
							probs[i] = (1.0/3)
						else:
							probs[i] = (1.0/3) - errorRate/9 # ⅓ - ε/9
			currentPos+=1
		elif m[0]<4:
			currentPos+=1
		else:
			currentPos=m[1]


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


numRefs=[0]
# this function is called only for input trees, when probVect lists are in the tree already but the other genome lists are not there yet.
#Initialize the MAT references and re-define the probVect's accordingly.
def setUpMAT(tree,node):
	lastNode=None
	direction=0
	mutationsAdded=[]
	probVect=tree.probVect
	isRef=tree.isRef
	mutations=tree.mutations
	children=tree.children
	up=tree.up
	while node!=None:
		if direction==0:
			newProbVect=[]
			if isRef[node]:
				newMutationsAdded=[]
				numRefs[0]+=1
			indProb=0
			indMut=0
			lastPos=0
			probVect1=probVect[node]
			mutations1=mutationsAdded
			entry=probVect1[indProb]
			posEntry=1
			if entry[0]==4 or entry[0]==5:
				posEntry=entry[1]
			if mutations1:
				mut=mutations1[0]
				posMut=mut[0]
			else:
				mut=None
				posMut=lRef+1
			while True:
				if posEntry<posMut:
					if entry[0]<4 and isRef[node]:
						#modify the mutation list to add a new mutation
						newMutationsAdded.append((posEntry,entry[0]))
						mutations[node].append((posEntry,entry[1],entry[0]))
						if len(entry)==2:
							newEntry=(4,posEntry)
						elif len(entry)==3:
							newEntry=(4,posEntry,entry[2])
						elif len(entry)==4:
							newEntry=(4,posEntry,entry[2],entry[3])
						newProbVect.append(newEntry)
					else:
						newProbVect.append(entry)

					if posEntry==lRef:
						break
					
					lastPos=posEntry
					indProb+=1
					entry=probVect1[indProb]
					if entry[0]==4 or entry[0]==5:
						posEntry=entry[1]
					else:
						posEntry+=1
				
				elif posEntry>posMut:
					if entry[0]==4 and isRef[node]:
						mutations[node].append((posMut,mut[1],refIndeces[posMut-1]))
					elif entry[0]==4 :
						if (posMut-1)>lastPos:
							if len(entry)==2:
								newEntry=(4,posMut-1)
							elif len(entry)==3:
								newEntry=(4,posMut-1,entry[2])
							elif len(entry)==4:
								newEntry=(4,posMut-1,entry[2],entry[3])
							newProbVect.append(newEntry)
						if len(entry)==2:
							newEntry=(refIndeces[posMut-1],mut[1])
						elif len(entry)==3:
							newEntry=(refIndeces[posMut-1],mut[1],entry[2])
						elif len(entry)==4:
							newEntry=(refIndeces[posMut-1],mut[1],entry[2],entry[3])
						newProbVect.append(newEntry)
						lastPos=posMut
					elif isRef[node]:
						newMutationsAdded.append(mut)

					indMut+=1
					if indMut<len(mutations1):
						mut=mutations1[indMut]
						posMut=mut[0]
					else:
						mut=None
						posMut=lRef+1

				else: #posEntry==posMut
					if entry[0]==6:
						if len(entry)==3:
							newEntry=(6,mut[1],entry[2])
						else:
							newEntry=(6,mut[1],entry[2],entry[3])
						newProbVect.append(newEntry)
						if isRef[node]:
							newMutationsAdded.append(mut)
					elif entry[0]==5:
						newProbVect.append(entry)
						if isRef[node]:
							newMutationsAdded.append(mut)
					elif entry[0]==mut[1]:
						if len(entry)==2:
							newEntry=(4,posEntry)
						elif len(entry)==3:
							newEntry=(4,posEntry,entry[2])
						elif len(entry)==4:
							newEntry=(4,posEntry,entry[2],entry[3])
						newProbVect.append(newEntry)
						if isRef[node]:
							newMutationsAdded.append(mut)
					else:
						if entry[0]==4 and isRef[node]:
							newProbVect.append(entry)
							mutations[node].append((posMut,mut[1],refIndeces[posMut-1]))
						elif entry[0]==4:
							if (posMut-1)>lastPos:
								if len(entry)==2:
									newEntry=(4,posMut-1)
								elif len(entry)==3:
									newEntry=(4,posMut-1,entry[2])
								elif len(entry)==4:
									newEntry=(4,posMut-1,entry[2],entry[3])
								newProbVect.append(newEntry)
							
							if len(entry)==2:
								newEntry=(refIndeces[posMut-1],mut[1])
							elif len(entry)==3:
								newEntry=(refIndeces[posMut-1],mut[1],entry[2])
							elif len(entry)==4:
								newEntry=(refIndeces[posMut-1],mut[1],entry[2],entry[3])
							newProbVect.append(newEntry)
						else:
							if isRef[node]:
								if len(entry)==2:
									newEntry=(4,posMut)
								elif len(entry)==3:
									newEntry=(4,posMut,entry[2])
								elif len(entry)==4:
									newEntry=(4,posMut,entry[2],entry[3])
								newProbVect.append(newEntry)
								newMutationsAdded.append((posMut,entry[0]))
								mutations[node].append((posMut,mut[1],entry[0]))
							else:
								if len(entry)==2:
									newEntry=(entry[0],mut[1])
								elif len(entry)==3:
									newEntry=(entry[0],mut[1],entry[2])
								elif len(entry)==4:
									newEntry=(entry[0],mut[1],entry[2],entry[3])
								newProbVect.append(newEntry)

					indMut+=1
					lastPos=posMut
					if indMut<len(mutations1):
						mut=mutations1[indMut]
						posMut=mut[0]
					else:
						mut=None
						posMut=lRef+1
					if posEntry==lRef:
						break
					indProb+=1
					entry=probVect1[indProb]
					if entry[0]==4 or entry[0]==5:
						posEntry=entry[1]
					else:
						posEntry+=1
			
			shorten(newProbVect)
			probVect[node]=newProbVect

			if children[node]:
				if isRef[node]:
					mutationsAdded=newMutationsAdded
				node=children[node][0]
			else:
				lastNode=node
				node=up[node]
				direction=1
		else:
			if lastNode==children[node][0]:
				node=children[node][1]
				direction=0
			else:
				#if reference node, remove mutations from mutationsAdded
				if isRef[node]:
					newMutationsAdded=[]
					indexMut=0
					indexAdded=0
					if mutations[node]:
						mut=mutations[node][0]
						posMut=mut[0]
					else:
						mut=None
						posMut=lRef+1
					if mutationsAdded:
						added=mutationsAdded[0]
						posAdded=added[0]
					else:
						added=None
						posAdded=lRef+1
					while posAdded<=lRef or posMut<=lRef:
						if posMut<posAdded:
							newMutationsAdded.append((posMut,mut[1]))
							indexMut+=1
							if indexMut<len(mutations[node]):
								mut=mutations[node][indexMut]
								posMut=mut[0]
							else:
								mut=None
								posMut=lRef+1
						elif posMut>posAdded:
							newMutationsAdded.append(added)
							indexAdded+=1
							if indexAdded<len(mutationsAdded):
								added=mutationsAdded[indexAdded]
								posAdded=added[0]
							else:
								added=None
								posAdded=lRef+1
						else:
							if mut[1]!=refIndeces[posMut-1]:
								newMutationsAdded.append((posMut,mut[1]))
							indexMut+=1
							if indexMut<len(mutations[node]):
								mut=mutations[node][indexMut]
								posMut=mut[0]
							else:
								mut=None
								posMut=lRef+1
							indexAdded+=1
							if indexAdded<len(mutationsAdded):
								added=mutationsAdded[indexAdded]
								posAdded=added[0]
							else:
								added=None
								posAdded=lRef+1

					mutationsAdded=newMutationsAdded

				lastNode=node
				node=up[node]
				direction=1


# Update .mutations for a severed and re-appended node, by traversing all the way up to the MRCA and then back again,
# all the way updating the mutation list.
def traverseTreeToUpdateMutationList(tree,appendedNode,node):
	#First find MRCA of the two nodes
	distFromRootAppended=0
	up=tree.up
	mutations=tree.mutations
	pNodeAppended=up[appendedNode]
	while pNodeAppended!=None:
		pNodeAppended=up[pNodeAppended]
		distFromRootAppended+=1
	distFromRoot=0
	pNode=up[node]
	while pNode!=None:
		pNode=up[pNode]
		distFromRoot+=1

	nodeList=[node]
	pNode=node
	pNodeAppended=appendedNode
	while distFromRootAppended>distFromRoot:
		pNodeAppended=up[pNodeAppended]
		distFromRootAppended-=1
	while distFromRootAppended<distFromRoot:
		pNode=up[pNode]
		nodeList.append(pNode)
		distFromRoot-=1
	
	while pNodeAppended!=pNode:
		pNode=up[pNode]
		nodeList.append(pNode)
		pNodeAppended=up[pNodeAppended]
	MRCAnode=pNodeAppended
	nodeList.pop()

	#now move mutation list up to the MRCA
	pNodeAppended=up[appendedNode]
	while pNodeAppended!=MRCAnode:
		if mutations[pNodeAppended]:
			mutations[appendedNode] = mergeMutationLists(mutations[pNodeAppended],mutations[appendedNode])
		pNodeAppended=up[pNodeAppended]
	while nodeList:
		newNode=nodeList.pop()
		if mutations[newNode]:
			mutations[appendedNode] = mergeMutationLists(mutations[newNode],mutations[appendedNode],downward=True)
	return


#merge two likelihood vectors to create a new one 
# if isUpDown, then probVect1 is assumed to be an upper vector and probVect2 a lower vector.
#(and also calculate the logLk of the merging if necessary, which is currently only useful for the root but could also be useful to calculate the overall total likelihood of the tree).
# fromTip1 and fromTip1 tell us if the two vectors come from tip nodes.
def mergeVectors(probVect1,bLen1,fromTip1,probVect2,bLen2,fromTip2,returnLK=False,isUpDown=False,numMinor1=0,numMinor2=0):
	indexEntry1, indexEntry2, pos, totalFactor = 0, 0, 0, 1.0
	probVect=[]
	entry1=probVect1[indexEntry1]
	entry2=probVect2[indexEntry2]
	totSum, cumErrorRate, refNucToPass, flag1, flag2, i, j ,i1, i2, newPos = 0.0, 0.0, -1, False, False, 0, 0, 0, 0, 0
	newVec, newVec2 = [], []

	if not (usingErrorRate and errorRateSiteSpecific): 
		errorRate=errorRateGlobal
	if not useRateVariation:
		mutMatrix=mutMatrixGlobal

	#contribution to non-mutation and non-error probabilities for the whole genome
	if returnLK:
		cumulPartLk=(bLen1+bLen2)*globalTotRate
		if usingErrorRate:
			# account for minor sequences not being erroneous (it helps not over-estimating error rates when there are many identical sequences)
			if fromTip1 or numMinor1:
				cumulPartLk+=totError*(1+numMinor1)
			if fromTip2 or numMinor2:
				cumulPartLk+=totError*(1+numMinor2)

	while True:
		if entry1[0]==5:
			if entry2[0]==5:
				newPos=min(entry1[1],entry2[1])
				probVect.append((5,newPos))
			elif entry2[0]<5:
				if entry2[0]<4:
					newPos=pos+1
					newEl=entry2[1]
				else:
					newPos=min(entry1[1],entry2[1])
					newEl=newPos
				if isUpDown:
					if usingErrorRate:
						if len(entry2)==2:
							if bLen2 or fromTip2:
								probVect.append((entry2[0],newEl,bLen2,0.0,fromTip2))
							else:
								probVect.append((entry2[0],newEl))
						elif len(entry2)==3:
							probVect.append((entry2[0],newEl,bLen2,0.0,entry2[3]))
						else:
							probVect.append((entry2[0],newEl,entry2[2]+bLen2,0.0,entry2[3]))
					else:
						if len(entry2)>2:
							probVect.append((entry2[0],newEl,entry2[2]+bLen2,0.0))
						else:
							if bLen2:
								probVect.append((entry2[0],newEl,bLen2,0.0))
							else:
								probVect.append((entry2[0],newEl))
				else:
					if usingErrorRate:
						if len(entry2)==2:
							if bLen2 or fromTip2:
								probVect.append((entry2[0],newEl,bLen2,fromTip2))
							else:
								probVect.append((entry2[0],newEl))
						elif len(entry2)==3:
							if bLen2:
								probVect.append((entry2[0],newEl,bLen2,entry2[3]))
							else:
								probVect.append((entry2[0],newEl,entry2[3]))
						else:
							probVect.append((entry2[0],newEl,entry2[2]+bLen2,entry2[3]))
					else:
						if len(entry2)>2:
							probVect.append((entry2[0],newEl,entry2[2]+bLen2))
						else:
							if bLen2:
								probVect.append((entry2[0],newEl,bLen2))
							else:
								probVect.append((entry2[0],newEl))

			else: # case entry2 is "O" and entry1 is "N"
				newPos=pos+1
				if isUpDown:
					if useRateVariation:
						mutMatrix=mutMatrices[pos]
					totBLen=bLen2
					if len(entry2)>3:
						totBLen+=entry2[2]
					if totBLen:
						newVec=getPartialVec(6, totBLen, mutMatrix, 0, vect=entry2[-1])
					else:
						newVec=list(entry2[-1])
					for i in range4:
						newVec[i]*=rootFreqs[i]
					totSum=sum(newVec)
					for i in range4:
						newVec[i]/=totSum
					probVect.append((6,entry2[1],newVec))

				else:
					if len(entry2)>3:
						probVect.append((6,entry2[1],entry2[2]+bLen2,entry2[3]))
					else:
						if bLen2:
							probVect.append((6,entry2[1],bLen2,entry2[2]))
						else:
							probVect.append((6,entry2[1],entry2[2]))

			if returnLK:
				cumulPartLk+=(bLen1+bLen2)*(cumulativeRate[pos]-cumulativeRate[newPos])
				if usingErrorRate:
					if fromTip1 or fromTip2:
						if errorRateSiteSpecific: cumErrorRate =  cumulativeErrorRate[newPos] - cumulativeErrorRate[pos] #both this and the next one should end up being positive
						else: cumErrorRate = errorRate * (newPos-pos)
					if fromTip1: #here we do not remove the contribution in case numMinor1>0 since this would anyway have to be re-added at some point
						cumulPartLk+=cumErrorRate
					if fromTip2:
						cumulPartLk+=cumErrorRate
			pos=newPos

		elif entry2[0]==5: #entry2 is N
			if entry1[0]<5:
				if entry1[0]<4:
					newPos=pos+1
					newEl=entry1[1]
				else:
					newPos=min(entry1[1],entry2[1])
					newEl=newPos
				if isUpDown:
					if usingErrorRate:
						if len(entry1)==2:
							if bLen1:
								probVect.append((entry1[0],newEl,bLen1,False))
							else:
								probVect.append((entry1[0],newEl))
						elif len(entry1)==3:
							probVect.append((entry1[0],newEl,bLen1,entry1[2]))
						elif len(entry1)==4:
							probVect.append((entry1[0],newEl,entry1[2]+bLen1, entry1[3]))
						else:
							probVect.append((entry1[0],newEl,entry1[2],entry1[3]+bLen1, entry1[4]))
					else:
						if len(entry1)==2:
							if bLen1:
								probVect.append((entry1[0],newEl,bLen1))
							else:
								probVect.append((entry1[0],newEl))
						elif len(entry1)==3:
							probVect.append((entry1[0],newEl,entry1[2]+bLen1))
						else:
							probVect.append((entry1[0],newEl,entry1[2],entry1[3]+bLen1))
							
				else:
					if usingErrorRate:
						if len(entry1)==2:
							if bLen1 or fromTip1:
								probVect.append((entry1[0],newEl,bLen1,fromTip1))
							else:
								probVect.append((entry1[0],newEl))
						elif len(entry1)==3:
							if bLen1:
								probVect.append((entry1[0],newEl,bLen1,entry1[3]))
							else:
								probVect.append((entry1[0],newEl,entry1[3]))
						else:
							probVect.append((entry1[0],newEl,entry1[2]+bLen1,entry1[3]))
					else:
						if len(entry1)>2:
							probVect.append((entry1[0],newEl,entry1[2]+bLen1))
						else:
							if bLen1:
								probVect.append((entry1[0],newEl,bLen1))
							else:
								probVect.append((entry1[0],newEl))
			else: #entry1 is "O" and entry2 is "N"
				newPos=pos+1
				refNucToPass=entry1[1]
				if isUpDown and ((len(entry1)==4 and entry1[2]>0) or bLen1):
					if useRateVariation:
						mutMatrix=mutMatrices[pos]
					totBLen=bLen1
					if len(entry1)>3:
						totBLen+=entry1[2]
					if totBLen:
						newVec=getPartialVec(6, totBLen, mutMatrix, 0, vect=entry1[-1],upNode=True)
					else:
						newVec=list(entry1[-1])
					totSum=sum(newVec)
					for i in range4:
						newVec[i]/=totSum
					probVect.append((6,entry1[1],newVec))
				else:
					if len(entry1)>3:
						probVect.append((6,entry1[1],entry1[2]+bLen1,entry1[3]))
					else:
						if bLen1:
							probVect.append((6,entry1[1],bLen1,entry1[2]))
						else:
							probVect.append((6,entry1[1],entry1[2]))

			if returnLK:
				cumulPartLk+=(bLen1+bLen2)*(cumulativeRate[pos]-cumulativeRate[newPos])
				if usingErrorRate:
					if fromTip1 or fromTip2:
						if errorRateSiteSpecific: cumErrorRate = cumulativeErrorRate[newPos] - cumulativeErrorRate[pos] #both this and the next one should end up being positive
						else: cumErrorRate = errorRate * (newPos-pos)
					if fromTip1:
						cumulPartLk+=cumErrorRate
					if fromTip2:
						cumulPartLk+=cumErrorRate
			pos=newPos
				
		else: #entry1 and entry2 are not "N"
			totLen1=bLen1
			if entry1[0]==6:
				if len(entry1)>3:
					totLen1+=entry1[2]
			elif len(entry1)>(2+usingErrorRate):
				totLen1+=entry1[2]
				if len(entry1)>(3+usingErrorRate):
					totLen1+=entry1[3]
			totLen2=bLen2
			if len(entry2)>(2+(usingErrorRate or entry2[0]==6)):
				totLen2+=entry2[2]

			flag1=(usingErrorRate and (entry1[0]!=6) and ((len(entry1)>2 and entry1[-1]) or fromTip1))
			flag2=(usingErrorRate and (entry2[0]!=6) and ((len(entry2)>2 and entry2[-1]) or fromTip2))

			if entry1[0]==4 and entry2[0]==4:
				newPos=min(entry1[1],entry2[1])
			else:
				newPos=pos+1

			if returnLK:
				if entry1[0]==4 and entry2[0]==4:
					if totLen2>bLen2 or totLen1>bLen1:
						cumulPartLk+=(totLen2-bLen2+totLen1-bLen1)*(cumulativeRate[newPos]-cumulativeRate[pos])
						if usingErrorRate:
							if ((not fromTip1) and flag1) or ((not fromTip2) and flag2):
								if errorRateSiteSpecific: cumErrorRate = cumulativeErrorRate[pos] - cumulativeErrorRate[newPos] #both this and the next one should end up being negative
								else: cumErrorRate = errorRate * (pos-newPos)
								if ((not fromTip1) and flag1):
									cumulPartLk+=cumErrorRate
								if ((not fromTip2) and flag2):
									cumulPartLk+=cumErrorRate
				else: # removing pre-calculated overall contributions to the likelihood from this position
					if entry1[0]!=4:
						refNucToPass=entry1[1]
					else:
						refNucToPass=entry2[1]
					if useRateVariation:
						cumulPartLk-=mutMatrices[pos][refNucToPass][refNucToPass]*(bLen2+bLen1)
					else:
						cumulPartLk-=mutMatrix[refNucToPass][refNucToPass]*(bLen2+bLen1)
					if usingErrorRate and ((entry1[0] != entry2[0]) or entry1[0]==6) and (fromTip1 or fromTip2):
						if errorRateSiteSpecific: cumErrorRate = errorRates[pos]
						else: cumErrorRate = errorRate 
						if fromTip1:
							cumulPartLk+=cumErrorRate
						if fromTip2:
							cumulPartLk+=cumErrorRate

			if entry2[0]==entry1[0] and entry2[0]<5: #entry1 and entry2 are two identical nucleotides
				
				if entry1[0]==4:
					probVect.append((4,newPos))
				else:
					probVect.append((entry1[0],entry1[1]))
					
					if returnLK :
						if useRateVariation:
							cumulPartLk+=mutMatrices[pos][entry1[0]][entry1[0]]*(totLen1+totLen2)
						else:
							cumulPartLk+=mutMatrix[entry1[0]][entry1[0]]*(totLen1+totLen2)
						if usingErrorRate:
							if ((not fromTip1) and flag1) or ((not fromTip2) and flag2):
								if errorRateSiteSpecific: cumErrorRate = errorRates[pos]
								else: cumErrorRate = errorRate
								if ((not fromTip1) and flag1):
									cumulPartLk-=cumErrorRate
								if ((not fromTip2) and flag2):
									cumulPartLk-=cumErrorRate
                
			elif (not totLen1) and (not totLen2) and entry1[0]<5 and entry2[0]<5 and (not flag1) and (not flag2): #0 distance between different nucleotides: merge is not possible
				if returnLK:
					print("mergeVectors() returning None 1")
					raise Exception("exit")
				else:
					return None
			else:
				if usingErrorRate and errorRateSiteSpecific: errorRate = errorRates[pos]
				if useRateVariation:	mutMatrix=mutMatrices[pos]

				if entry1[0]==4:
					refNucToPass=entry2[1]
					i1=refNucToPass
				else:
					refNucToPass=entry1[1]
					i1=entry1[0]
				if i1<=4:
					if totLen1 or flag1:
						if isUpDown and len(entry1)>3+usingErrorRate:
							newVec = getPartialVec(i1, entry1[2], mutMatrix, errorRate, flag=flag1)
							for i in range4:
								newVec[i]*=rootFreqs[i]
							if entry1[3]+bLen1:
								newVec = getPartialVec(6, entry1[3]+bLen1, mutMatrix, 0, vect=newVec, upNode=True)
						else:
							newVec=getPartialVec(i1, totLen1, mutMatrix, errorRate, flag=flag1, upNode=isUpDown)
					else:
						newVec=[0.0,0.0,0.0,0.0]
						newVec[i1]=1.0
				else:
					if totLen1:
						newVec=getPartialVec(6, totLen1, mutMatrix, 0, vect=entry1[-1],upNode=isUpDown)
					else:
						newVec=list(entry1[-1])

				if entry2[0]==4:
					i2=refNucToPass
				else:
					i2=entry2[0]
				if i2==6: #entry1 is a nucleotide and entry2 is "O"
					if totLen2:
						newVec2=getPartialVec(6, totLen2, mutMatrix, 0, vect=entry2[-1])
					else:
						newVec2=entry2[-1]
				else:
					if totLen2 or flag2:
						newVec2 = getPartialVec(i2, totLen2, mutMatrix, errorRate, flag=flag2)
					else:
						newVec2 = [0.0,0.0,0.0,0.0]
						newVec2[i2]=1.0

				for j in range4:
					newVec[j]*=newVec2[j]
				totSum=sum(newVec)
				if not totSum:
					if returnLK:
						print("mergeVectors() returning None 2")
						raise Exception("exit")
					else:
						return None
				for i in range4:
					newVec[i]/=totSum

				state =simplify(newVec,refNucToPass)
				if state==6:
					probVect.append((6,refNucToPass,newVec))
				else:
					if state==4:
						probVect.append((4,newPos))
					else:
						probVect.append((state,refNucToPass))

				if returnLK:
					totalFactor*=totSum

			pos=newPos

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
		if entry1[0]<4 or entry1[0]==6:
			indexEntry1+=1
			entry1=probVect1[indexEntry1]
		elif pos==entry1[1]:
			indexEntry1+=1
			entry1=probVect1[indexEntry1]
		if entry2[0]<4 or entry2[0]==6:
			indexEntry2+=1
			entry2=probVect2[indexEntry2]
		elif pos==entry2[1]:
			indexEntry2+=1
			entry2=probVect2[indexEntry2]

	if returnLK:
		return probVect, cumulPartLk+log(totalFactor)
	else:
		return probVect


#calculate the probability that results from combining a lower likelihood genome list of the root with root frequencies.
#assumes that when this is used, R regions are wrt the reference genome, not other parts of the MAT.
# this should be achievable since this function is always called close to the root, and the root should keep track of its differences wrt the reference genome.
def findProbRoot(probVect,node=None,mutations=None,up=None):
	while node!=None:
		if mutations[node]:
			probVect=passGenomeListThroughBranch(probVect,mutations[node],dirIsUp=True)
		node=up[node]
	if not (usingErrorRate and errorRateSiteSpecific): 
		errorRate=errorRateGlobal
	logLK=0.0
	logFactor=1.0
	pos=0
	for entry in probVect:
		if usingErrorRate and (entry[0]<5) and len(entry)>2 and entry[-1]:
			if entry[0]==4:
				logLK+=rootFreqsLogErrorCumulative[entry[1]]-rootFreqsLogErrorCumulative[pos]
				pos=entry[1]
			else:
				if errorRateSiteSpecific: errorRate = errorRates[pos]
				logFactor*=(rootFreqs[entry[0]]*(1.0-1.33333*errorRate)+0.33333*errorRate)
				pos+=1
		else:
			if entry[0]==4:
				for i in range4:
					logLK+=rootFreqsLog[i]*(cumulativeBases[entry[1]][i]-cumulativeBases[pos][i])
				pos=entry[1]
			elif entry[0]<4:
				logLK+=rootFreqsLog[entry[0]]
				pos+=1
			elif entry[0]==6:
				tot=0.0
				for i in range4:
					tot+=rootFreqs[i]*entry[-1][i]
				logFactor*=tot
				pos+=1
			else:
				pos=entry[1]
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
	logLK+=log(logFactor)

	return logLK


#for the root, take lower likelihood genome list probVect, and create an overall likelihood (or upper right or upper left) genome list by multiplying likelihoods by root frequencies.
def rootVector(probVect,bLen,isFromTip,tree,node):
	#first representing the genome list wrt the reference genome - this is so to that mutationLists can be accurate.
	nodeList=[]
	mutations=tree.mutations
	up=tree.up
	if mutations[node]:
		probVect=passGenomeListThroughBranch(probVect, mutations[node], dirIsUp=True) #, modifyCurrentList=False
	nodeList.append(node)
	node=up[node]
	while node!=None:
		nodeList.append(node)
		if mutations[node]:
			probVect=passGenomeListThroughBranch(probVect, mutations[node], dirIsUp=True) #, modifyCurrentList=True
		node=up[node]
	newProbVect=[]
	newPos=0
	for entry in probVect:
		if entry[0]==5:
			newProbVect.append(entry)
			newPos=entry[1]
		elif entry[0]==6:
			totBLen=bLen
			if len(entry)>3:
				totBLen+=entry[2]
			if totBLen:
				if useRateVariation:
					newVect=getPartialVec(6, totBLen, mutMatrices[newPos], 0, vect=entry[-1])
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
			newPos+=1
		else:
			if usingErrorRate:
				flag1 = ((len(entry)>2) and entry[-1]) or isFromTip
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
			if entry[0]<4:
				newPos+=1
			else:
				newPos=entry[1]

	while nodeList:
		node=nodeList.pop()
		if mutations[node]:
			newProbVect=passGenomeListThroughBranch(newProbVect, mutations[node], dirIsUp=False)
	shorten(newProbVect)
					
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
					pseudoMutCounts[entry2[1]][entry2[0]]+=1		
				elif entry2[0]==4:
					pseudoMutCounts[entry1[0]][entry1[1]]+=1
				else:
					pseudoMutCounts[entry1[0]][entry2[0]]+=1
				pos+=1
			else:
				if (entry1[0]==4 or entry1[0]==5) and (entry2[0]==4 or entry2[0]==5):
					pos=min(entry1[1],entry2[1])
				else:
					pos+=1

			if pos==lRef:
				break
			if entry1[0]<4 or entry1[0]==6:
				indexEntry1+=1
				entry1=probVect1[indexEntry1]
			elif pos==entry1[1]:
				indexEntry1+=1
				entry1=probVect1[indexEntry1]
			if entry2[0]<4 or entry2[0]==6:
				indexEntry2+=1
				entry2=probVect2[indexEntry2]
			elif pos==entry2[1]:
				indexEntry2+=1
				entry2=probVect2[indexEntry2]


#function to optimize branch lengths.
#calculate features of the derivative of the likelihood cost function wrt the branch length, then finds branch length that minimizes likelihood cost.
def estimateBranchLengthWithDerivative(probVectP,probVectC,fromTipC=False):
	if not (usingErrorRate and errorRateSiteSpecific): 
		errorRate=errorRateGlobal
	if not useRateVariation:
		mutMatrix=mutMatrixGlobal
	c1=globalTotRate
	ais=[]
	indexEntry1, indexEntry2, pos, contribLength, nZeros = 0, 0, 0, 0.0, 0
	entry1=probVectP[indexEntry1]
	entry2=probVectC[indexEntry2]
	while True:
		if entry2[0]==5:
			if entry1[0]==4 or entry1[0]==5:
				end=min(entry1[1],entry2[1])
			else:
				end=pos+1
			c1+=(cumulativeRate[pos]-cumulativeRate[end])
			pos=end
		elif entry1[0]==5: # case entry1 or entry2 is N
			#if parent node is type "N", in theory we might have to calculate the contribution of root nucleotides; 
			# however, if this node is "N" then every other node in the current tree is "N", so we can ignore this since this contribution cancels out in relative terms.
			if entry2[0]==4:
				end=min(entry1[1],entry2[1])
			else:
				end=pos+1
			c1+=(cumulativeRate[pos]-cumulativeRate[end])
			pos=end
		else:
			#below, when necessary, we represent the likelihood as coeff0*l +coeff1, where l is the branch length to be optimized.
			if entry1[0]==4 and entry2[0]==4: # case entry1 and entry2 are R	
				pos=min(entry1[1],entry2[1])
			else:
				if useRateVariation:
					mutMatrix=mutMatrices[pos]
				
				if entry1[0]==4:
					c1-=mutMatrix[entry2[1]][entry2[1]]
				else:
					c1-=mutMatrix[entry1[1]][entry1[1]]
				flag1 = (usingErrorRate and (entry1[0]!=6) and len(entry1)>2 and entry1[-1])
				flag2 = (usingErrorRate and (entry2[0]!=6) and (fromTipC or (len(entry2)>2 and entry2[-1])))
				if usingErrorRate and errorRateSiteSpecific: errorRate = errorRates[pos]

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
						i1=entry2[1]
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
							i1=entry2[1]
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
								if mutMatrix[entry2[1]][entry2[0]]:
									coeff0+=errorRate*0.33333/mutMatrix[entry2[1]][entry2[0]]
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
							i2=entry1[1]
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
								i2=entry1[1]
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
		if entry1[0]<4 or entry1[0]==6:
			indexEntry1+=1
			entry1=probVectP[indexEntry1]
		elif pos==entry1[1]:
			indexEntry1+=1
			entry1=probVectP[indexEntry1]
		if entry2[0]<4 or entry2[0]==6:
			indexEntry2+=1
			entry2=probVectC[indexEntry2]
		elif pos==entry2[1]:
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
#TODO when chaging branch lengths, update nDesc0
def updateBLen(tree,cNode,addToList=False,nodeList=None):
	up=tree.up
	probVectUpRight=tree.probVectUpRight
	probVectUpLeft=tree.probVectUpLeft
	probVect=tree.probVect
	mutations=tree.mutations
	children=tree.children
	dirty=tree.dirty
	dist=tree.dist
	minorSequences=tree.minorSequences
	#TODO
	nDesc0=tree.nDesc0
	node=up[cNode]
	if cNode==children[node][0]:
		vectUp=probVectUpRight[node]
		cNum=0
	else:
		vectUp=probVectUpLeft[node]
		cNum=1
	if mutations[cNode]:
		vectUp=passGenomeListThroughBranch(vectUp, mutations[cNode])
	fromTipC=(len(children[cNode])==0 and len(minorSequences[cNode])==0)
	bestLength=estimateBranchLengthWithDerivative(vectUp,probVect[cNode],fromTipC=fromTipC)
	#TODO updating in case of branch length change
	if HnZ:
		#print("run updateBLen, initial len "+str(dist[cNode])+" final: "+str(bestLength))
		if dist[cNode]>effectivelyNon0BLen and (bestLength<=effectivelyNon0BLen):
			addendum0=(nDesc0[cNode]-1)
		elif (dist[cNode]<=effectivelyNon0BLen) and bestLength>effectivelyNon0BLen:
			addendum0=(1-nDesc0[cNode])
		else:
			addendum0=0
		if addendum0:
			nDesc0[up[cNode]]+=addendum0
			parent0=up[cNode]
			while up[parent0]!=None and (dist[parent0]<=effectivelyNon0BLen):
				parent0=up[parent0]
				nDesc0[parent0]+=addendum0

	dist[cNode]=bestLength
	dirty[node]=True
	dirty[cNode]=True
	if addToList:
		nodeList.append((cNode,2))
		nodeList.append((node,cNum))


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
			if entry1[0]<4:
				pos+=1
			else:
				pos=min(entry1[1],entry2[1])
		elif entry1[0]==6:
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
			pos+=1
		else:
			pos=min(entry1[1],entry2[1])
		if pos==lRef:
			break
		if entry1[0]<4 or entry1[0]==6:
			indexEntry1+=1
			entry1=probVect1[indexEntry1]
		elif pos==entry1[1]:
			indexEntry1+=1
			entry1=probVect1[indexEntry1]
		if entry2[0]<4 or entry2[0]==6:
			indexEntry2+=1
			entry2=probVect2[indexEntry2]
		elif pos==entry2[1]:
			indexEntry2+=1
			entry2=probVect2[indexEntry2]
	return False


# update the partials iteratively starting from the nodes in nodeList
#each entry in nodeList contains the node it refers to, and the direction where the update comes from (0 is left child, 1 is right child, 2 is parent)
def updatePartials(tree,nodeList):
	dirty=tree.dirty
	up=tree.up
	children=tree.children
	probVectUpRight=tree.probVectUpRight
	probVectUpLeft=tree.probVectUpLeft
	mutations=tree.mutations
	minorSequences=tree.minorSequences
	dist=tree.dist
	probVect=tree.probVect
	probVectTotUp=tree.probVectTotUp
	while nodeList:
		updatedBLen=False # if there has been an inconsistency, function updateBLen() has been called, and so there is no point continuing with some updates.
		madeChange=False # some change has been made, so continue traversal
		node, direction = nodeList.pop()
		dirty[node]=True
		if up[node] != None:
			try:
				if node==children[up[node]][0]:
					childNumUp=0
					vectUpUp=probVectUpRight[up[node]]
				else:
					childNumUp=1
					vectUpUp=probVectUpLeft[up[node]]
				if mutations[node]:
					vectUpUp=passGenomeListThroughBranch(vectUpUp,mutations[node])
			except AttributeError:
				vectUpUp=None
		isTip=(len(children[node])==0 and len(minorSequences[node])==0)
		#change in likelihoods is coming from parent node
		if direction==2:
			if dist[node] : #if necessary, update the total probabilities at the mid node.
				newTot=mergeVectors(vectUpUp,dist[node]/2,False,probVect[node],dist[node]/2,isTip,isUpDown=True)
				if newTot==None:
					if dist[node]>1e-100:
						print("inside updatePartials(), from parent: should not have happened since node.dist>0")
					updateBLen(tree,node)
					nodeList.append((up[node],childNumUp))
					newTot=mergeVectors(vectUpUp,dist[node]/2,False,probVect[node],dist[node]/2,isTip,isUpDown=True)
					madeChange=True
				probVectTotUp[node]=newTot
				shorten(probVectTotUp[node])
			else:
				probVectTotUp[node]=None
			if children[node]:# and (not updatedBLen): #at valid internal node, update upLeft and upRight, and if necessary add children to nodeList.
				child0Vect=probVect[children[node][0]]
				if mutations[children[node][0]]:
					child0Vect=passGenomeListThroughBranch(child0Vect,mutations[children[node][0]],dirIsUp=True)
				child1Vect=probVect[children[node][1]]
				if mutations[children[node][1]]:
					child1Vect=passGenomeListThroughBranch(child1Vect,mutations[children[node][1]],dirIsUp=True)
				dist0=dist[children[node][0]]
				dist1=dist[children[node][1]]
				isTip0=(len(children[children[node][0]])==0 and len(minorSequences[children[node][0]])==0)
				isTip1=(len(children[children[node][1]])==0 and len(minorSequences[children[node][1]])==0)
				newUpRight=mergeVectors(vectUpUp,dist[node],False,child1Vect,dist1,isTip1,isUpDown=True)
				if newUpRight==None:
					if (not dist[node]) and (not dist1):
						updateBLen(tree,node)
						if not dist[node]:
							updateBLen(tree,children[node][1],addToList=True,nodeList=nodeList)
							updatedBLen=True
						else:
							probVectTotUp[node]=mergeVectors(vectUpUp,dist[node]/2,False,probVect[node],dist[node]/2,isTip,isUpDown=True)
							newUpRight=mergeVectors(vectUpUp,dist[node],False,child1Vect,dist1,isTip1,isUpDown=True)
							nodeList.append((up[node],childNumUp))
							madeChange=True
					else:
						print("Strange: None vector from non-zero distances in updatePartials() from parent direction.")
						raise Exception("exit")
				if not updatedBLen:
					newUpLeft=mergeVectors(vectUpUp,dist[node],False,child0Vect,dist0,isTip0,isUpDown=True)
					if newUpLeft==None:
						if (not dist[node]) and (not dist0) :
							updateBLen(tree,node)
							if not dist[node]:
								updateBLen(tree,children[node][0],addToList=True,nodeList=nodeList)
								updatedBLen=True
							else:
								probVectTotUp[node]=mergeVectors(vectUpUp,dist[node]/2,False,probVect[node],dist[node]/2,isTip,isUpDown=True)
								newUpRight=mergeVectors(vectUpUp,dist[node],False,child1Vect,dist1,isTip1,isUpDown=True)
								newUpLeft=mergeVectors(vectUpUp,dist[node],False,child0Vect,dist0,isTip0,isUpDown=True)
								nodeList.append((up[node],childNumUp))
								madeChange=True
						else:
							print("Strange: None vector from non-zero distances in updatePartials() from parent direction, child0.")
							raise Exception("exit")
				if not updatedBLen:
					if madeChange or areVectorsDifferent(probVectUpRight[node],newUpRight):
						probVectUpRight[node]=newUpRight
						shorten(probVectUpRight[node])
						nodeList.append((children[node][0],2))
					if madeChange or areVectorsDifferent(probVectUpLeft[node],newUpLeft):
						probVectUpLeft[node]=newUpLeft
						shorten(probVectUpLeft[node])
						nodeList.append((children[node][1],2))

		else: #change in likelihoods is coming from child number "direction".
			childNum=direction
			otherChildNum=1-childNum
			childDist=dist[children[node][childNum]]
			otherChildDist=dist[children[node][otherChildNum]]
			otherChildVect=probVect[children[node][otherChildNum]]
			if mutations[children[node][otherChildNum]]:
				otherChildVect=passGenomeListThroughBranch(otherChildVect,mutations[children[node][otherChildNum]],dirIsUp=True)
			probVectDown=probVect[children[node][childNum]]
			if mutations[children[node][childNum]]:
				probVectDown=passGenomeListThroughBranch(probVectDown,mutations[children[node][childNum]],dirIsUp=True)
			isTip=(len(children[children[node][childNum]])==0 and len(minorSequences[children[node][childNum]])==0)
			otherIsTip=(len(children[children[node][otherChildNum]])==0 and len(minorSequences[children[node][otherChildNum]])==0)
			try:
				if childNum:
					otherVectUp=probVectUpRight[node]
				else:
					otherVectUp=probVectUpLeft[node]
			except AttributeError:
				otherVectUp=None

			#update lower likelihoods
			newVect=mergeVectors(otherChildVect,otherChildDist,otherIsTip,probVectDown,childDist,isTip)
			if newVect==None:
				if (not childDist) and (not otherChildDist):
					updateBLen(tree,children[node][childNum])
					if not dist[children[node][childNum]]:
						updateBLen(tree,children[node][otherChildNum],addToList=True,nodeList=nodeList)
						updatedBLen=True
					else:
						childDist=dist[children[node][childNum]]
						probVect[node]=mergeVectors(otherChildVect,otherChildDist,otherIsTip,probVectDown,childDist,isTip)
						nodeList.append((children[node][childNum],2))
						madeChange=True
				else:
					print("Strange: None vector from non-zero distances in updatePartials() from child direction.")
					raise Exception("exit")
			else:
				try:
					oldProbVect=probVect[node]
				except AttributeError:
					oldProbVect=None
				probVect[node]=newVect
				shorten(probVect[node])

			#update total mid-branches likelihood
			if (not updatedBLen) and dist[node] and (up[node] != None) and (vectUpUp!=None):
				newTot=mergeVectors(vectUpUp,dist[node]/2,False,probVect[node],dist[node]/2,False,isUpDown=True)
				if newTot==None:
					updateBLen(tree,node)
					probVect[node]=mergeVectors(otherChildVect,otherChildDist,otherIsTip,probVectDown,childDist,isTip)
					nodeList.append((children[node][childNum],2))
					probVectTotUp[node]=mergeVectors(vectUpUp,dist[node]/2,False,probVect[node],dist[node]/2,False,isUpDown=True)
					madeChange=True
					print("inside updatePartials(), from child: should not have happened since node.dist>0")
				else:
					probVectTotUp[node]=newTot
					shorten(probVectTotUp[node])
					
			if (not updatedBLen) and (otherVectUp!=None):
				#update likelihoods at sibling node
				if up[node] != None:
					newUpVect=mergeVectors(vectUpUp,dist[node],False,probVectDown,childDist,isTip,isUpDown=True)
				else:
					newUpVect=rootVector(probVectDown,childDist,isTip,tree,node)
				if newUpVect==None:
					if (not dist[node]) and (not childDist):
						updateBLen(tree,node)
						if not dist[node]:
							updateBLen(tree,children[node][childNum],addToList=True,nodeList=nodeList)
							updatedBLen=True
						else:
							probVectTotUp[node]=mergeVectors(vectUpUp,dist[node]/2,False,probVect[node],dist[node]/2,False,isUpDown=True)
							nodeList.append((children[node][childNum],2))
							madeChange=True
							newUpVect=mergeVectors(vectUpUp,dist[node],False,probVectDown,childDist,isTip,isUpDown=True)
					else:
						print("Strange: None vector from non-zero distances in updatePartials() from child direction, newUpVect.")
						raise Exception("exit")
			if (not updatedBLen) and (otherVectUp!=None):
				if madeChange or areVectorsDifferent(otherVectUp,newUpVect):
					if childNum:
						probVectUpRight[node]=newUpVect
						shorten(probVectUpRight[node])
					else:
						probVectUpLeft[node]=newUpVect
						shorten(probVectUpLeft[node])
					nodeList.append((children[node][otherChildNum],2))

			if (not updatedBLen):
				#update likelihoods at parent node
				if madeChange or areVectorsDifferent(probVect[node],oldProbVect):
					if up[node] != None:
						nodeList.append((up[node],childNumUp))


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
						print('diff bLen or flag case')
						return True
					if len(entry1)>3:
						if abs(entry1[3] - entry2[3])>thresholdProb:
							print('diff bLen2 or flag case')
							return True
						if len(entry1)>4:
							if abs(entry1[4] - entry2[4])>thresholdProb:
								print('diff flag case')
								return True
			if entry1[0]<4 or entry2[0]<4:
				pos+=1
			else:
				pos=min(entry1[1],entry2[1])
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
			pos+=1
		else:
			if not (entry1[0]==5 and entry2[0]==5):
				if (entry1[0]==5 or entry2[0]==5):
					print("N case")
					return True
				if entry1[0]<5:
					if entry1[0]==4:
						i1=entry2[1]
					else:
						i1=entry1[0]
					if entry2[-1][i1]+threshold<1.0:
						print("i1 case")
						print(entry1)
						print(entry2)
						return True
				elif entry2[0]<5:
					if entry2[0]==4:
						i2=entry1[1]
					else:
						i2=entry2[0]
					if entry1[-1][i2]+threshold<1.0:
						print("i2 case")
						print(entry1)
						print(entry2)
						return True
			pos+=1
		#update pos, end, etc
		if pos==lRef:
			break
		if entry1[0]<4 or entry1[0]==6:
			indexEntry1+=1
			entry1=probVect1[indexEntry1]
		elif pos==entry1[1]:
			indexEntry1+=1
			entry1=probVect1[indexEntry1]
		if entry2[0]<4 or entry2[0]==6:
			indexEntry2+=1
			entry2=probVect2[indexEntry2]
		elif pos==entry2[1]:
			indexEntry2+=1
			entry2=probVect2[indexEntry2]
	return False


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
				if entry2[0]==4:
					pos=min(entry1[1],entry2[1])
				else:
					pos+=1
				found2bigger=True
			elif entry2[0]==5:
				if entry1[0]==4:
					pos=min(entry1[1],entry2[1])
				else:
					pos+=1
				found1bigger=True
			elif entry1[0]==6:
				if entry2[0]==4:
					i2=entry1[1]
				else:
					i2=entry2[0]
				if entry1[-1][i2]>0.1:
					found2bigger=True
				else:
					return 0
				pos+=1
			elif entry2[0]==6:
				if entry1[0]==4:
					i1=entry2[1]
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
			if entry1[0]<4:
				pos+=1
			else:
				pos=min(entry1[1],entry2[1])
		if found1bigger and found2bigger:
			return 0
		if pos==lRef:
			break
		if entry1[0]<4 or entry1[0]==6:
			indexEntry1+=1
			entry1=probVect1[indexEntry1]
		elif pos==entry1[1]:
			indexEntry1+=1
			entry1=probVect1[indexEntry1]
		if entry2[0]<4 or entry2[0]==6:
			indexEntry2+=1
			entry2=probVect2[indexEntry2]
		elif pos==entry2[1]:
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


numMinorsRemoved=[0]
numNodes=[0,0,0,0,0,0,0,0]
#Given a tree, and a substitution rate matrix, re-calculate all genome lists within the tree according to this matrix.
# this is useful once the matrix estimation has finished, to make sure all genome lists reflect this matrix. 
def reCalculateAllGenomeLists(tree,root, checkExistingAreCorrect=False,countNodes=False,countPseudoCounts=False,pseudoMutCounts=None,data=None,names=None,firstSetUp=False): #checkSamplesIntree=False
	up=tree.up
	children=tree.children
	probVectUpRight=tree.probVectUpRight
	probVectUpLeft=tree.probVectUpLeft
	mutations=tree.mutations
	minorSequences=tree.minorSequences
	dist=tree.dist
	probVect=tree.probVect
	probVectTotUp=tree.probVectTotUp
	name=tree.name
	if firstSetUp:
		isRef=[False]*len(up)
		tree.isRef=isRef
		nDesc=tree.nDesc
	#first pass to update all lower likelihoods.
	node=root
	#direction 0 means from parent, direction 1 means from a child
	lastNode=None
	direction=0
	dataNamesConverted=False
	while node!=None:
		if direction==0:
			if children[node]:
				node=children[node][0]
			else:
				if firstSetUp:
					if data==None:
						print("Problem, initializing genome lists but there is no data.")
						exit()
					if (names[name[node]] in data):
						probVect[node]=probVectTerminalNode(data[names[name[node]]],tree,node)
						shorten(probVect[node])
					elif not dataNamesConverted:
						print("Sample name "+name[node]+" not found in the input sequence data - all samples in the input tree need to have a corresponding sequence entry. Converting sequence names by replacing ? and & characters with _ and trying again.")
						oldNames=list(data.keys())
						for name in oldNames:
							newName=name.replace("?","_").replace("&","_")
							if newName!=name:
								data[newName]=data[name]
						dataNamesConverted=True
						oldNames=None
						if names[name[node]] in data:
							probVect[node]=probVectTerminalNode(data[names[name[node]]],tree,node)
							shorten(probVect[node])
						else:
							print("Error: sample name "+names[name[node]]+" not found in the input sequence data - all samples in the input tree need to have a corresponding sequence entry.")
							raise Exception("exit")
					else:
						print("Error: sample name "+names[name[node]]+" not found in the input sequence data - all samples in the input tree need to have a corresponding sequence entry.")
						raise Exception("exit")
					
					# removed minor sequences from the tree
					tryRemovingMinor=False
					if children[up[node]][1]==node and (not dist[node]):
						sibling=children[up[node]][0]
						if (not dist[sibling]) and (not children[sibling]):
							tryRemovingMinor=True
					while tryRemovingMinor:
						#test if one of the samples is minor of the other
						#TODO also accounting for HnZ since minor sequences might have different modifiers for different placements
						if errorRateSiteSpecificFile or errorRateFixed or estimateErrorRate or estimateSiteSpecificErrorRate or supportFor0Branches or HnZ:
							comparison=isMinorSequence(probVect[node],probVect[sibling],onlyFindIdentical=True)
						else:
							comparison=isMinorSequence(probVect[node],probVect[sibling])
						if comparison==1:
							majorNode=node
							minorNode=sibling
						elif comparison==2:
							majorNode=sibling
							minorNode=node
						else:
							break
						numMinorsRemoved[0]+=1
						#remove the minor sample and its parent
						minorSequences[majorNode].append(name[minorNode])
						for mino in minorSequences[minorNode]:
							minorSequences[majorNode].append(mino)
						probVect[minorNode]=None
						parentNode=up[majorNode]
						up[majorNode]=up[parentNode]
						dist[majorNode]=dist[parentNode]
						if up[majorNode]!=None:
							if children[up[majorNode]][0]==parentNode:
								children[up[majorNode]][0]=majorNode
							else:
								children[up[majorNode]][1]=majorNode
						children[parentNode]=None
							
						#update tryRemovingMinor to check if continuing removal
						tryRemovingMinor=False
						node=majorNode
						if up[node]!=None:
							if children[up[node]][1]==node and (not dist[node]):
								sibling=children[up[node]][0]
								if (not dist[sibling]) and (not children[sibling]):
									tryRemovingMinor=True
					

				if (not onlyNambiguities) and usingErrorRate:
					updateProbVectTerminalNode(probVect[node],minorSequences[node])

				if countNodes:
					numNodes[0]+=1
					for entry in probVect[node]:
						if entry[0]<4:
							numNodes[1]+=1
						elif entry[0]==4:
							numNodes[2]+=1
						elif entry[0]==5:
							numNodes[3]+=1
						else:
							numNodes[4]+=1
					numNodes[5]+=len(mutations[node])
				lastNode=node
				node=up[node]
				direction=1
		else :
			if lastNode==children[node][0]:
				node=children[node][1]
				direction=0
			else:
				if firstSetUp:
					if children[children[node][0]] and (not isRef[children[node][0]]):
						nDesc[node]+=nDesc[children[node][0]]
					if children[children[node][1]] and (not isRef[children[node][1]]):
						nDesc[node]+=nDesc[children[node][1]]
					if dist[children[node][0]]:
						nDesc[node]+=1
					if dist[children[node][0]]:
						nDesc[node]+=1
					if nDesc[node]>=maxNumDescendantsForMATClade and dist[node]:
						nDesc[node]=0
						isRef[node]=True

				isTip0=(len(children[children[node][0]])==0 ) and (len(minorSequences[children[node][0]])==0 )
				isTip1=(len(children[children[node][1]])==0 ) and (len(minorSequences[children[node][1]])==0 )
				probVect0=probVect[children[node][0]]
				if mutations[children[node][0]]:
					probVect0=passGenomeListThroughBranch(probVect0,mutations[children[node][0]],dirIsUp=True)
				probVect1=probVect[children[node][1]]
				if mutations[children[node][1]]:
					probVect1=passGenomeListThroughBranch(probVect1,mutations[children[node][1]],dirIsUp=True)
				newLower=mergeVectors(probVect0,dist[children[node][0]],isTip0,probVect1,dist[children[node][1]],isTip1)
				if checkExistingAreCorrect:
					if areVectorsDifferentDebugging(newLower,probVect[node]):
						print("Inside reCalculateAllGenomeLists(), new lower at node is different from the old one, and it shouldn't be.")
						raise Exception("exit")
				if newLower==None:
					if (not dist[children[node][0]]) and (not dist[children[node][1]]):
						if firstSetUp:
							dist[children[node][0]]=oneMutBLen/2
							dist[children[node][1]]=oneMutBLen/2
						else:
							updateBLen(tree,children[node][0])
							if (not dist[children[node][0]]):
								updateBLen(tree,children[node][1])
						probVect[node]=mergeVectors(probVect0,dist[children[node][0]],isTip0,probVect1,dist[children[node][1]],isTip1)
						if probVect[node]==None:
							dist[children[node][0]]=oneMutBLen/2
							dist[children[node][1]]=oneMutBLen/2
							probVect[node]=mergeVectors(probVect0,dist[children[node][0]],isTip0,probVect1,dist[children[node][1]],isTip1)
							if probVect[node]==None:
								print("None vector when merging two vectors in reCalculateAllGenomeLists(), despite updating branch lengths.")
								raise Exception("exit")
					else:
						print("Strange, distances>0 but inconsistent lower genome list creation in reCalculateAllGenomeLists()")
						raise Exception("exit")
				else:
					probVect[node]=newLower
					shorten(probVect[node])
				if countNodes:
					numNodes[0]+=1
					for entry in probVect[node]:
						if entry[0]<4:
							numNodes[1]+=1
						elif entry[0]==4:
							numNodes[2]+=1
						elif entry[0]==5:
							numNodes[3]+=1
						else:
							numNodes[4]+=1
					numNodes[5]+=len(mutations[node])
				lastNode=node
				node=up[node]
				direction=1

	#if reading an input tree, set up the MAT
	if firstSetUp and useLocalReference:
		setUpMAT(tree,root)	

	#now update the other genome lists for the root
	node=root
	if children[node]:
		probVect1=probVect[children[node][1]]
		if mutations[children[node][1]]:
			probVect1=passGenomeListThroughBranch(probVect1,mutations[children[node][1]],dirIsUp=True)
		newVect=rootVector(probVect1,dist[children[node][1]],(len(children[children[node][1]])==0 and len(minorSequences[children[node][1]])==0),tree,node)
		if checkExistingAreCorrect:
			if areVectorsDifferentDebugging(newVect,probVectUpRight[node]):
				print("new probVectUpRight at root is different from the old one, and it shouldn't be.")
				raise Exception("exit")
		probVectUpRight[node]=newVect
		probVect0=probVect[children[node][0]]
		if mutations[children[node][0]]:
			probVect0=passGenomeListThroughBranch(probVect0,mutations[children[node][0]],dirIsUp=True)
		newVect=rootVector(probVect0,dist[children[node][0]],(len(children[children[node][0]])==0 and len(minorSequences[children[node][0]])==0),tree,node)
		if checkExistingAreCorrect:
			if areVectorsDifferentDebugging(newVect,probVectUpLeft[node]):
				print("new probVectUpLeft at root is different from the old one, and it shouldn't be; distance of child node "+str(dist[children[node][0]]))
				raise Exception("exit")
		probVectUpLeft[node]=newVect
		
		#now traverse the tree downward and update the non-lower genome lists for all other nodes of the tree.
		totNodeList=[]
		lastNode=None
		node=children[node][0]
		direction=0
		while node!=None:
			if direction==0:
				if node==children[up[node]][0]:
					vectUp=probVectUpRight[up[node]]
					nodeChildNum=0
				else:
					vectUp=probVectUpLeft[up[node]]
					nodeChildNum=1
				if mutations[node]:
					vectUp=passGenomeListThroughBranch(vectUp,mutations[node])
				if dist[node]:
					isTip=(len(children[node])==0 ) and (len(minorSequences[node])==0 )
					if countPseudoCounts:
						updatePesudoCounts(vectUp,probVect[node],pseudoMutCounts)
					newVect=mergeVectors(vectUp,dist[node]/2,False,probVect[node],dist[node]/2,isTip,isUpDown=True)
					if checkExistingAreCorrect:
						if areVectorsDifferentDebugging(newVect,probVectTotUp[node]):
							print("new probVectTotUp at node is different from the old one, and it shouldn't be.")
							print(newVect)
							print(probVectTotUp[node])
							raise Exception("exit")
					shorten(newVect)
					probVectTotUp[node]=newVect
				if children[node]:
					isTip0=(len(children[children[node][0]])==0 ) and (len(minorSequences[children[node][0]])==0 )
					isTip1=(len(children[children[node][1]])==0 ) and (len(minorSequences[children[node][1]])==0 )
					probVect0=probVect[children[node][0]]
					if mutations[children[node][0]]:
						probVect0=passGenomeListThroughBranch(probVect0,mutations[children[node][0]],dirIsUp=True)
					probVect1=probVect[children[node][1]]
					if mutations[children[node][1]]:
						probVect1=passGenomeListThroughBranch(probVect1,mutations[children[node][1]],dirIsUp=True)
					newUpRight=mergeVectors(vectUp,dist[node],False,probVect1,dist[children[node][1]],isTip1,isUpDown=True)
					if newUpRight==None:
						if (not dist[children[node][1]]) and (not dist[node]):
							updateBLen(tree,node)
							if not dist[node]:
								if firstSetUp:
									probVectUpLeft[node]=mergeVectors(vectUp,dist[node],False,probVect0,dist[children[node][0]],isTip0,isUpDown=True)
								updateBLen(tree,children[node][1])
								totNodeList.append((node,1))
							else:
								probVectTotUp[node]=mergeVectors(vectUp,dist[node]/2,False,probVect[node],dist[node]/2,False,isUpDown=True)
								totNodeList.append((up[node],nodeChildNum))
							probVectUpRight[node]=mergeVectors(vectUp,dist[node],False,probVect1,dist[children[node][1]],isTip1,isUpDown=True)
						else:
							print("Strange, distances>0 but inconsistent upRight genome list creation in reCalculateAllGenomeLists()")
							raise Exception("exit")
					else:
						if checkExistingAreCorrect:
							if areVectorsDifferentDebugging(newUpRight,probVectUpRight[node]):
								print("new probVectUpRight at node is different from the old one, and it shouldn't be.")
								raise Exception("exit")
						shorten(newUpRight)
						probVectUpRight[node]=newUpRight
					newUpLeft=mergeVectors(vectUp,dist[node],False,probVect0,dist[children[node][0]],isTip0,isUpDown=True)
					if newUpLeft==None:
						if(not dist[children[node][0]]) and (not dist[node]):
							updateBLen(tree,children[node][0])
							if not dist[children[node][0]]:
								updateBLen(tree,node)
								totNodeList.append((up[node],nodeChildNum))
								probVectTotUp[node]=mergeVectors(vectUp,dist[node]/2,False,probVect[node],dist[node]/2,isTip,isUpDown=True)
								probVectUpRight[node]=mergeVectors(vectUp,dist[node],False,probVect1,dist[children[node][1]],isTip1,isUpDown=True)
							else:
								totNodeList.append((node,0))
							probVectUpLeft[node]=mergeVectors(vectUp,dist[node],False,probVect0,dist[children[node][0]],isTip0,isUpDown=True)
						else:
							print("Strange, distances>0 but inconsistent upLeft genome list creation in reCalculateAllGenomeLists()")
							raise Exception("exit")
					else:
						if checkExistingAreCorrect:
							if areVectorsDifferentDebugging(newUpLeft,probVectUpLeft[node]):
								print("new probVectUpLeft at node is different from the old one, and it shouldn't be.")
								raise Exception("exit")
						shorten(newUpLeft)
						probVectUpLeft[node]=newUpLeft
					node=children[node][0]
				else:
					lastNode=node
					node=up[node]
					direction=1
			else:
				if lastNode==children[node][0]:
					node=children[node][1]
					direction=0
				else:
					lastNode=node
					node=up[node]
					direction=1

		updatePartials(tree,totNodeList)


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


def updateErrorRates(errorRate,errorRates=None):
	global rootFreqsLogErrorCumulative
	global totError
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
		totError=-cumulativeErrorRate[-1]
	else:
		rootFreqsLogErrorCumulative=[0]*(lRef+1)
		rootFreqsLogErrorCumulative[0]=0.0
		for i in range(0,len(ref)):
			rootFreqsLogErrorCumulative[i+1] = rootFreqsLogErrorCumulative[i] + log(rootFreqs[refIndeces[i]]*(1.0-1.33333*errorRate)+0.333333*errorRate)
		totError=-errorRate*lRef


#reading input subst model to be used as initial values for parameter estimation and tree building.
if inputRates!="":
	if not os.path.isfile(inputRates):
		print("Input rates file "+inputRates+" not found, quitting MAPLE. Use option --inputRates to specify a valid input rates file, identical in format to MAPLE's output file with model parameters.")
		raise Exception("exit")
	fileR=open(inputRates)
	for i in range4:
		line=fileR.readline()
		linelist=line.split()
		for j in range4:
			mutMatrixGlobal[i][j]=float(linelist[j])
	if rateVariation:
		siteRates=[]
		while line!="Site rates:\n":
			line=fileR.readline()
		for i in range(lRef):
			line=fileR.readline()
			linelist=line.split()
			siteRates.append(float(linelist[1]))
		useRateVariation=True
		mutMatrices=[]
		updateMutMatrices(mutMatrixGlobal,siteRates=siteRates)
	if estimateSiteSpecificErrorRate:
		usingErrorRate=True
		errorRates=[]
		while line!="Site error rates:\n":
			line=fileR.readline()
		for i in range(lRef):
			line=fileR.readline()
			linelist=line.split()
			errorRates.append(float(linelist[1]))
		errorRateGlobal=sum(errorRates)/lRef
		updateErrorRates(errorRateGlobal,errorRates=errorRates)
	fileR.close()
	print("Read input rates")


#In case an input tree is given, calculate all genome lists, pseudocounts and rates, and then recalculate genome lists again according to the new rates.
if inputTree!="":
	#samplesAlreadyInTree=reCalculateAllGenomeLists(tree1,countPseudoCounts=True,pseudoMutCounts=pseudoMutCounts,data=data,firstSetUp=True,checkSamplesIntree=True)
	if not inputRates:
		reCalculateAllGenomeLists(tree1,rootIndex1,countPseudoCounts=True,pseudoMutCounts=pseudoMutCounts,data=data,names=namesInTree,firstSetUp=True)
		if (model!="JC" and updateSubMatrix(pseudoMutCounts,model,mutMatrixGlobal)):
			for i in range4:
				nonMutRates[i]=mutMatrixGlobal[i][i]
			for i in range(lRef):
				cumulativeRate[i+1]=cumulativeRate[i]+nonMutRates[refIndeces[i]]
		reCalculateAllGenomeLists(tree1,rootIndex1)
	else:
		reCalculateAllGenomeLists(tree1,rootIndex1,data=data,names=namesInTree,firstSetUp=True)
	print("Genome list for initial tree and initial pseudocounts calculated.")


#Sort samples based on distance from reference, but punishing more isolated N's and ambiguity characters.
#more ambiguous sequences are placed last this way - this is useful since ambiguous sequences are harder to place and are more likely to be less informative (and so be removed from the analysis altogether)
#def distancesFromRefPunishNs(data,samples):
def distancesFromRefPunishNs(data,samples=None,samplesInInitialTree=set(),forgetData=False):
	sampleDistances=[]
	if samples==None:
		rangeInd=range(len(data))
	else:
		rangeInd=samples
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
		elif forgetData:
			data[diffIndex]=None

	print("Now doing sorting")
	sampleDistances.sort(reverse=True,key=itemgetter(0))
	return sampleDistances


#function to calculate likelihood cost of appending node to parent node 
# this now mostly ignores the differences in contributions of non-mutating nucleotides.
def appendProbNode(probVectP,probVectC,isTipC,bLen):
	if not (usingErrorRate and errorRateSiteSpecific): 
		errorRate=errorRateGlobal
	if not useRateVariation:
		mutMatrix=mutMatrixGlobal
	indexEntry1, indexEntry2, totalFactor, pos = 0, 0, 1.0, 0
	entry1=probVectP[indexEntry1] #parent
	entry2=probVectC[indexEntry2] #child
	contribLength=bLen
	Lkcost=bLen*globalTotRate
	if usingErrorRate and isTipC:
		Lkcost+=totError
	while True:
		if entry2[0]==5:
			if entry1[0]==4 or entry1[0]==5:
				pos=min(entry1[1],entry2[1])
				if pos==lRef:
					break
				if entry1[1]==pos:
					indexEntry1+=1
					entry1=probVectP[indexEntry1]
			else:
				pos+=1
				if pos==lRef:
					break
				indexEntry1+=1
				entry1=probVectP[indexEntry1]
			if entry2[1]==pos:
				indexEntry2+=1
				entry2=probVectC[indexEntry2]

		elif entry1[0]==5: # case entry1 is N
			#if parent node is type "N", in theory we might have to calculate the contribution of root nucleotides; 
			# however, if this node is "N" then every other node in the current tree is "N", so we can ignore this since this contribution cancels out in relative terms.
			#Here we are now ignoring the effect of N's on the probabilities of non-mutations: these are almost the same independent of the placement!
			if entry2[0]==4:
				pos=min(entry1[1],entry2[1])
				if pos==lRef:
					break
				if entry2[1]==pos:
					indexEntry2+=1
					entry2=probVectC[indexEntry2]
			else:
				pos+=1
				if pos==lRef:
					break
				indexEntry2+=1
				entry2=probVectC[indexEntry2]
			if entry1[1]==pos:
				indexEntry1+=1
				entry1=probVectP[indexEntry1]

		else:
			#contribLength will be here the total length from the root or from the upper node, down to the down node.
			if entry1[0]!=entry2[0] or entry1[0]==6:
				contribLength=bLen
				if entry1[0]<5:
					if len(entry1)==3+usingErrorRate:
						contribLength+=entry1[2]
					elif len(entry1)==4+usingErrorRate:
						contribLength+=entry1[3]
				elif len(entry1)==4:
					contribLength+=entry1[2]
				if entry2[0]<5:
					if len(entry2)==3+usingErrorRate:
						contribLength+=entry2[2]
				elif len(entry2)==4:
					contribLength+=entry2[2]

			if entry1[0]==4: # case entry1 is R	
				if entry2[0]==4:
					pos=min(entry1[1],entry2[1])
					if pos==lRef:
						break
					if entry2[1]==pos:
						indexEntry2+=1
						entry2=probVectC[indexEntry2]

				#entry1 is reference and entry2 is of type "O"
				elif entry2[0]==6:
					if useRateVariation:
						mutMatrix=mutMatrices[pos]
					i1=entry2[1]
					if entry2[-1][i1]>0.02:
						totalFactor*=entry2[-1][i1]
					else:
						if len(entry1)==4+usingErrorRate:
							flag1 = (usingErrorRate and (len(entry1)>2) and entry1[-1] )
							tot=0.0
							if usingErrorRate and errorRateSiteSpecific: errorRate = errorRates[pos]
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
					if pos==lRef:
						break
					indexEntry2+=1
					entry2=probVectC[indexEntry2]

				else: #entry1 is R and entry2 is a different but single nucleotide
					flag2 = (usingErrorRate and (isTipC or (len(entry2)>2) and entry2[-1] ))
					if useRateVariation:
						mutMatrix=mutMatrices[pos]
					if len(entry1)==4+usingErrorRate:
						flag1 = (usingErrorRate and (len(entry1)>2) and entry1[-1] )
						i1=entry2[1]
						i2=entry2[0]
						if usingErrorRate and errorRateSiteSpecific: errorRate = errorRates[pos]
						tot3=getPartialVec(i2, contribLength, mutMatrix, errorRate, flag=flag2)
						tot2=getPartialVec(i1, entry1[2], mutMatrix, errorRate, flag=flag1)
						tot=0.0
						for i in range4:
							tot+=tot3[i]*tot2[i]*rootFreqs[i]
						totalFactor*=tot/rootFreqs[i1]
					else:
						if flag2:
							if usingErrorRate and errorRateSiteSpecific: errorRate = errorRates[pos]
							totalFactor*=min(0.25,mutMatrix[entry2[1]][entry2[0]]*contribLength)+errorRate*0.33333
						else:
							if contribLength:
								totalFactor*=min(0.25,mutMatrix[entry2[1]][entry2[0]]*contribLength)
							else:
								return float("-inf")
					pos+=1
					if pos==lRef:
						break
					indexEntry2+=1
					entry2=probVectC[indexEntry2]
				if entry1[1]==pos:
					indexEntry1+=1
					entry1=probVectP[indexEntry1]

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
						i2=entry1[1]
					else:
						i2=entry2[0]
					if entry1[-1][i2]>0.02:
						totalFactor*=entry1[-1][i2]
					else:
						if (usingErrorRate and (isTipC or (len(entry2)>2) and entry2[-1] )):
							if errorRateSiteSpecific: errorRate = errorRates[pos]
							tot3=getPartialVec(i2, contribLength, mutMatrix, errorRate, flag=True)
						else:
							tot3=getPartialVec(i2, contribLength, mutMatrix, None, flag=False)
						tot=0.0
						for j in range4:
							tot+=entry1[-1][j]*tot3[j]
						totalFactor*=tot
				pos+=1
				if pos==lRef:
					break
				indexEntry1+=1
				entry1=probVectP[indexEntry1]
				if entry2[0]!=4 or entry2[1]==pos:
					indexEntry2+=1
					entry2=probVectC[indexEntry2]

			else: #entry1 is a non-ref nuc
				if entry2[0]!=entry1[0]:
					flag1 = (usingErrorRate and (len(entry1)>2) and entry1[-1] )
					if useRateVariation:
						mutMatrix=mutMatrices[pos]

					i1=entry1[0]
					if entry2[0]<5: #entry2 is a nucleotide
						if entry2[0]==4:
							i2=entry1[1]
						else:
							i2=entry2[0]
						flag2 = (usingErrorRate and (isTipC or (len(entry2)>2) and entry2[-1] ))
						if len(entry1)==4+usingErrorRate:
							if usingErrorRate and errorRateSiteSpecific: errorRate = errorRates[pos]
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
						if usingErrorRate and errorRateSiteSpecific: errorRate = errorRates[pos]
						if entry2[-1][i1]>0.02:
							totalFactor*=entry2[-1][i1]
						else:
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
				if pos==lRef:
					break
				indexEntry1+=1
				entry1=probVectP[indexEntry1]
				if entry2[0]!=4 or entry2[1]==pos:
					indexEntry2+=1
					entry2=probVectC[indexEntry2]

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

	return Lkcost+log(totalFactor)


#number of samples that could have been placed as major of another sample but weren't due to sample placement order
totalMissedMinors=[0]
totalTimeFindingParent=[0.0]
#function traversing the tree to find the best node in the tree where to re-append the given subtree (rooted at node.children[child]) to improve the topology of the current tree.
# bestLKdiff is the best likelihood cost found for the current placement (after optimizing the branch length).
# removedBLen is such branch length that optimizes the current placement - it will be used to place the subtree attached at other nodes of the tree.
#TODO account for HnZ modifiers in SPR search
def findBestParentTopology(tree,node,child,bestLKdiff,removedBLen,strictTopologyStopRules=strictTopologyStopRules,allowedFailsTopology=allowedFailsTopology,thresholdLogLKtopology=thresholdLogLKtopology,sequentialSearch=True):
	timeStartParentTopology=time()
	up=tree.up
	children=tree.children
	probVectUpRight=tree.probVectUpRight
	probVectUpLeft=tree.probVectUpLeft
	mutations=tree.mutations
	minorSequences=tree.minorSequences
	dist=tree.dist
	probVect=tree.probVect
	probVectTotUp=tree.probVectTotUp
	name=tree.name
	#TODO
	nDesc0=tree.nDesc0
	if HnZ:
		originalParent0=node
		while (dist[originalParent0]<=effectivelyNon0BLen) and up[originalParent0]!=None:
			originalParent0=up[originalParent0]
	bestNode=children[node][1-child]
	bestNodes=[]
	nodesToVisit=[]
	removedPartialsRelative=probVect[children[node][child]]
	lastLK, distance, direction, needsUpdating, failedPasses = 0.0, 0.0, 0, True, 0
	if mutations[children[node][child]]:
		removedPartialsRelative=passGenomeListThroughBranch(removedPartialsRelative,mutations[children[node][child]],dirIsUp=True)
	bestRemovedPartials=removedPartialsRelative
	if mutations[bestNode]:
		bestRemovedPartials=passGenomeListThroughBranch(bestRemovedPartials,mutations[bestNode],dirIsUp=False)
	isRemovedTip=(len(children[children[node][child]])==0) and (len(minorSequences[children[node][child]])==0)
	originalLK=bestLKdiff
	originalPlacement=bestNode
	originalRemovedPartialsRelative=bestRemovedPartials

	if up[node]!=None:
		if children[up[node]][0]==node:
			childUp=1
			vectUpUp=probVectUpRight[up[node]]
		else:
			childUp=2
			vectUpUp=probVectUpLeft[up[node]]
		#the list nodesToVisit keeps track of the nodes of the tree we wil need to traverse, (first element of each entry),
		# the direction from where we are visitng them (0=from parent, and 1,2=from one of the children),
		# an updated genome list from the direction where we come from (taking into account the removal of the given subtree),
		# a branch length value separating the node from this updated genome list (useful for the fact that the removal of the subtree changes the branch length at the removal node),
		# a flag that says if the updated genome list passed needs still updating, or if it has become identical to the pre-existing genome list in the tree (which usually happens after a while),
		# the likelihood cost of appending at the last node encountered in this direction,
		# a number of consecutively failed traversal steps since the last improvement found (if this number goes beyond a threshold, traversal in the considered direction might be halted).
		probVect1=probVect[bestNode]
		if mutations[bestNode]:
			probVect1=passGenomeListThroughBranch(probVect1,mutations[bestNode],dirIsUp=True)
		removedPartialsRelative1=removedPartialsRelative
		if mutations[node]:
			probVect1=passGenomeListThroughBranch(probVect1,mutations[node],dirIsUp=True)
			removedPartialsRelative1=passGenomeListThroughBranch(removedPartialsRelative,mutations[node],dirIsUp=True)
		nodesToVisit.append((up[node],childUp,probVect1,dist[bestNode]+dist[node],bestLKdiff,0,removedPartialsRelative1))
		if mutations[node]:
			vectUpUp=passGenomeListThroughBranch(vectUpUp,mutations[node])
		removedPartialsRelative1=removedPartialsRelative
		if mutations[bestNode]:
			vectUpUp=passGenomeListThroughBranch(vectUpUp,mutations[bestNode])
			removedPartialsRelative1=passGenomeListThroughBranch(removedPartialsRelative,mutations[bestNode])
		nodesToVisit.append((bestNode,0,vectUpUp,dist[bestNode]+dist[node],bestLKdiff,0,removedPartialsRelative1))
		originalBLens=(dist[node],dist[bestNode],removedBLen)
	else:
		# case node is root
		if children[bestNode]: # case there is only one sample outside of the subtree doesn't need to be considered
			child1=children[bestNode][0]
			child2=children[bestNode][1]
			vectUp1=probVect[child2]
			if mutations[child2]:
				vectUp1=passGenomeListThroughBranch(vectUp1,mutations[child2],dirIsUp=True)
			vectUp1=rootVector(vectUp1,dist[child2],(len(children[child2])==0 and len(minorSequences[child2])==0),tree,node)
			if mutations[child1]:
				removedPartialsRelative1=passGenomeListThroughBranch(bestRemovedPartials,mutations[child1])
				vectUp1=passGenomeListThroughBranch(vectUp1,mutations[child1],dirIsUp=False)
			else:
				removedPartialsRelative1=bestRemovedPartials
			nodesToVisit.append((child1,0,vectUp1,dist[child1],bestLKdiff,0,removedPartialsRelative1))
			vectUp2=probVect[child1]
			if mutations[child1]:
				vectUp2=passGenomeListThroughBranch(vectUp2,mutations[child1],dirIsUp=True)
			vectUp2=rootVector(vectUp2,dist[child1],(len(children[child1])==0 and len(minorSequences[child1])==0),tree,node)
			if mutations[child2]:
				removedPartialsRelative2=passGenomeListThroughBranch(bestRemovedPartials,mutations[child2])
				vectUp2=passGenomeListThroughBranch(vectUp2,mutations[child2],dirIsUp=False)
			else:
				removedPartialsRelative2=bestRemovedPartials
			nodesToVisit.append((child2,0,vectUp2,dist[child2],bestLKdiff,0,removedPartialsRelative2))
		originalBLens=(0.0,dist[bestNode],removedBLen)
	bestBranchLengths=originalBLens

	while nodesToVisit:
		nodeToVisitInfo=nodesToVisit.pop()
		if len(nodeToVisitInfo)>5:
			needsUpdating=True
			t1,direction,passedPartials,distance,lastLK,failedPasses,removedPartialsRelative=nodeToVisitInfo #modifiable
		else:
			needsUpdating=False
			t1,direction,lastLK,failedPasses,removedPartialsRelative=nodeToVisitInfo #modifiable

		if direction==0:
			#consider the case we are moving from a parent to a child
			if dist[t1] and (not (up[t1]==node or up[t1]==None)):
				if needsUpdating:
					isTip=(len(children[t1])==0) and (len(minorSequences[t1])==0)
					midTot=mergeVectors(passedPartials,distance/2,False,probVect[t1],distance/2,isTip,isUpDown=True)
					if not areVectorsDifferent(midTot,probVectTotUp[t1]):
						needsUpdating=False
				else:
					midTot=probVectTotUp[t1]
					distance=dist[t1]
				if midTot==None:
					continue
				midProb=appendProbNode(midTot,removedPartialsRelative,isRemovedTip,removedBLen)
				#TODO correct new placement likelihood by the change in HnZ
				if HnZ:
					if removedBLen>effectivelyNon0BLen:
						midProb+=getHnZ(2) - getHnZ(1)
					else:
						midProb+=getHnZ(nDesc0[children[node][child]]+1) - getHnZ(nDesc0[children[node][child]])
				if midProb>bestLKdiff-thresholdLogLKoptimizationTopology:
					#if needsUpdating, then add to the tuple also the information on the up and down genome lists to use to recalculate intermediate genome lists at varying branch lengths
					if needsUpdating:
						bestNodes.append((t1,midProb,passedPartials,probVect[t1],distance,midTot,removedPartialsRelative))
					else:
						bestNodes.append((t1,midProb,removedPartialsRelative))
				if midProb>bestLKdiff:
					bestLKdiff=midProb
					bestNode=t1
					failedPasses=0
					shorten(removedPartialsRelative)
					bestRemovedPartials=removedPartialsRelative
					bestBranchLengths=(distance/2,distance/2,removedBLen)
				elif midProb<(lastLK-thresholdLogLKconsecutivePlacement): #placement at current node is considered failed if placement likelihood is not improved by a certain margin compared to best placement so far for the nodes above it.
					failedPasses+=1
			else:
				midProb=lastLK
			
			#keep crawling down into children nodes unless the stop criteria for the traversal are satisfied.
			traverseChildren=False
			if strictTopologyStopRules:
				if failedPasses<=allowedFailsTopology and midProb>(bestLKdiff-thresholdLogLKtopology) and children[t1]:
					traverseChildren=True
			else:
				if failedPasses<=allowedFailsTopology or midProb>(bestLKdiff-thresholdLogLKtopology):
					if children[t1]:
						traverseChildren=True
			if traverseChildren:
				child1=children[t1][0]
				if needsUpdating:
					otherChild=children[t1][1]
					isTip=(len(children[otherChild])==0) and (len(minorSequences[otherChild])==0)
					otherProbVect=probVect[otherChild]
					if mutations[otherChild]:
						otherProbVect=passGenomeListThroughBranch(otherProbVect,mutations[otherChild],dirIsUp=True)
					vectUpRight=mergeVectors(passedPartials,distance,False,otherProbVect,dist[otherChild],isTip,isUpDown=True)
				else:
					vectUpRight=probVectUpRight[t1]
				if vectUpRight!=None:
					removedPartialsRelative1=removedPartialsRelative
					if mutations[child1]:
						removedPartialsRelative1=passGenomeListThroughBranch(removedPartialsRelative,mutations[child1])
					if needsUpdating:			
						if mutations[child1]:
							vectUpRight=passGenomeListThroughBranch(vectUpRight,mutations[child1])
						nodesToVisit.append((child1,0,vectUpRight,dist[child1],midProb,failedPasses,removedPartialsRelative1))
					else:
						nodesToVisit.append((child1,0,midProb,failedPasses,removedPartialsRelative1))

				child1=children[t1][1]
				if needsUpdating:
					otherChild=children[t1][0]
					isTip=(len(children[otherChild])==0) and (len(minorSequences[otherChild])==0)
					otherProbVect=probVect[otherChild]
					if mutations[otherChild]:
						otherProbVect=passGenomeListThroughBranch(otherProbVect,mutations[otherChild],dirIsUp=True)
					vectUpLeft=mergeVectors(passedPartials,distance,False,otherProbVect,dist[otherChild],isTip,isUpDown=True)
				else:
					vectUpLeft=probVectUpLeft[t1]
				if vectUpLeft!=None:
					removedPartialsRelative1=removedPartialsRelative
					if mutations[child1]:
						removedPartialsRelative1=passGenomeListThroughBranch(removedPartialsRelative,mutations[child1])
					if needsUpdating:
						if mutations[child1]:
							vectUpLeft=passGenomeListThroughBranch(vectUpLeft,mutations[child1])
						nodesToVisit.append((child1,0,vectUpLeft,dist[child1],midProb,failedPasses,removedPartialsRelative1))
					else:
						nodesToVisit.append((child1,0,midProb,failedPasses,removedPartialsRelative1))

		else: #case when crawling up from child to parent
			otherChild=children[t1][2-direction]
			midBottom=None
			if dist[t1] and up[t1]!=None: #try appending mid-branch
				if needsUpdating:
					otherProbVect=probVect[otherChild]
					if mutations[otherChild]:
						otherProbVect=passGenomeListThroughBranch(otherProbVect,mutations[otherChild],dirIsUp=True)
					isTip=(len(children[otherChild])==0) and (len(minorSequences[otherChild])==0)
					midBottom=mergeVectors(passedPartials,distance,False,otherProbVect,dist[otherChild],isTip)
					if midBottom==None:
						continue
					if t1==children[up[t1]][0]:
						vectUp=probVectUpRight[up[t1]]
					else:
						vectUp=probVectUpLeft[up[t1]]
					if mutations[t1]:
						vectUp=passGenomeListThroughBranch(vectUp,mutations[t1])
					midTot=mergeVectors(vectUp,dist[t1]/2,False,midBottom,dist[t1]/2,False,isUpDown=True)
					if not probVectTotUp[t1]:
						print("Node has no probVectTotUp "+str(name[t1])+" with dist "+str(dist[t1])+" calculating new one")
						probVectTotUp[t1]=mergeVectors(vectUp,dist[t1]/2,False,probVect[t1],dist[t1]/2,False,isUpDown=True)
					if not areVectorsDifferent(midTot,probVectTotUp[t1]):
						needsUpdating=False
				else:
					midTot=probVectTotUp[t1]
				if midTot==None:
					continue
				midProb=appendProbNode(midTot,removedPartialsRelative,isRemovedTip,removedBLen)
				#TODO correct new placement likelihood by the change in HnZ
				if HnZ:
					if removedBLen>effectivelyNon0BLen:
						midProb+=getHnZ(2) - getHnZ(1)
					else:
						midProb+=getHnZ(nDesc0[children[node][child]]+1) - getHnZ(nDesc0[children[node][child]])
				if midProb>=(bestLKdiff-thresholdLogLKoptimizationTopology):
					#if needsUpdating, then add to the tuple also the information on the up and down genome lists to use to recalculate intermediate genome lists at varying branch lengths
					if needsUpdating:
						bestNodes.append((t1,midProb,vectUp,midBottom,dist[t1],midTot,removedPartialsRelative))
					else:
						bestNodes.append((t1,midProb,removedPartialsRelative))
				if midProb>bestLKdiff:
					bestLKdiff=midProb
					bestNode=t1
					failedPasses=0
					bestRemovedPartials=removedPartialsRelative
					bestBranchLengths=(dist[t1]/2,dist[t1]/2,removedBLen)
				elif midProb<(lastLK-thresholdLogLKconsecutivePlacement): #placement at current node is considered failed if placement likelihood is not improved by a certain margin compared to best placement so far for the nodes above it.
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
				if up[t1]!=None: #case the node is not the root
					#first pass the crawling down the other child
					if t1==children[up[t1]][0]:
						upChild=0
						if needsUpdating:
							vectUpUp=probVectUpRight[up[t1]]
					else:
						upChild=1
						if needsUpdating:
							vectUpUp=probVectUpLeft[up[t1]]
					if needsUpdating:
						if mutations[t1]:
							vectUpUp=passGenomeListThroughBranch(vectUpUp,mutations[t1])
						vectUp=mergeVectors(vectUpUp,dist[t1],False,passedPartials,distance,False,isUpDown=True)
					else:
						if direction==1:
							vectUp=probVectUpLeft[t1]
						else:
							vectUp=probVectUpRight[t1]

					if vectUp==None:
						continue
					else:
						removedPartialsRelative1=removedPartialsRelative
						if mutations[otherChild]:
							removedPartialsRelative1=passGenomeListThroughBranch(removedPartialsRelative,mutations[otherChild])
						if needsUpdating:
							if mutations[otherChild]:
								vectUp=passGenomeListThroughBranch(vectUp,mutations[otherChild])
							nodesToVisit.append((otherChild,0,vectUp,dist[otherChild],midProb,failedPasses,removedPartialsRelative1))
						else:
							nodesToVisit.append((otherChild,0,midProb,failedPasses,removedPartialsRelative1))
					#now pass the crawling up to the parent node
					if needsUpdating:
						if midBottom==None:
							otherProbVect=probVect[otherChild]
							if mutations[otherChild]:
								otherProbVect=passGenomeListThroughBranch(otherProbVect,mutations[otherChild],dirIsUp=True)
							isTip=(len(children[otherChild])==0) and (len(minorSequences[otherChild])==0)
							midBottom=mergeVectors(passedPartials,distance,False,otherProbVect,dist[otherChild],isTip)
							if midBottom==None:
								continue
					removedPartialsRelative1=removedPartialsRelative
					if mutations[t1]:
						removedPartialsRelative1=passGenomeListThroughBranch(removedPartialsRelative,mutations[t1],dirIsUp=True)
					if needsUpdating:
						if mutations[t1]:
							midBottom=passGenomeListThroughBranch(midBottom,mutations[t1],dirIsUp=True)
						nodesToVisit.append((up[t1],upChild+1,midBottom,dist[t1],midProb,failedPasses,removedPartialsRelative1))
					else:
						nodesToVisit.append((up[t1],upChild+1,midProb,failedPasses,removedPartialsRelative1))
				#now consider case of root node
				else:
					if needsUpdating:
						vectUp=rootVector(passedPartials,distance,False,tree,t1)
						if mutations[otherChild]:
							vectUp=passGenomeListThroughBranch(vectUp,mutations[otherChild])
					removedPartialsRelative1=removedPartialsRelative
					if mutations[otherChild]:
						removedPartialsRelative1=passGenomeListThroughBranch(removedPartialsRelative,mutations[otherChild])
					if needsUpdating:
						nodesToVisit.append((otherChild,0,vectUp,dist[otherChild],midProb,failedPasses,removedPartialsRelative1))
					else:
						nodesToVisit.append((otherChild,0,midProb,failedPasses,removedPartialsRelative1))

	#Initial exploration is finished.
	#Now, for each branch within threshold likelihood distance from the best found, optimize branch lengths.
	#Use optimized scores to select final best branch
	bestScore=bestLKdiff
	compensanteForBranchLengthChange=True
	if not bestNodes:
		totalTimeFindingParent[0]+=(time()-timeStartParentTopology)
		return originalPlacement, originalLK ,originalBLens, [], 1.0, originalRemovedPartialsRelative
	if aBayesPlusOn and sequentialSearch:
		if networkOutput:
			listOfProbableNodes=[]
		listofLKcosts=[]
		rootAlreadyConsidered=False
	for nodePair in bestNodes:
		score=nodePair[1]
		if score>=bestLKdiff-thresholdLogLKoptimizationTopology:
			t1=nodePair[0]
			if len(nodePair)==3:
				#optimize branch lengths of appendage
				if t1==children[up[t1]][0]:
					upVect=probVectUpRight[up[t1]]
				else:
					upVect=probVectUpLeft[up[t1]]
				if mutations[t1]:
					upVect=passGenomeListThroughBranch(upVect,mutations[t1])
				downVect=probVect[t1]
				distance=dist[t1]
				midTot=probVectTotUp[t1]
			else:
				upVect=nodePair[2]
				downVect=nodePair[3]
				distance=nodePair[4]
				midTot=nodePair[5]
			removedPartials=nodePair[-1]
			bestAppendingLength=estimateBranchLengthWithDerivative(midTot,removedPartials,fromTipC=isRemovedTip)
			#now optimize appending location
			fromTip1=(len(children[t1])==0) and (len(minorSequences[t1])==0)
			midLowerVector=mergeVectors(downVect,distance/2,fromTip1,removedPartials,bestAppendingLength,isRemovedTip)
			bestTopLength=estimateBranchLengthWithDerivative(upVect,midLowerVector,fromTipC=False)
			midTopVector=mergeVectors(upVect,bestTopLength,False,removedPartials,bestAppendingLength,isRemovedTip,isUpDown=True)
			if midTopVector==None:
				print("Inconsistent branch length found")
				print(upVect)
				bestTopLength=defaultBLen*0.1
				midTopVector=mergeVectors(upVect,bestTopLength,False,removedPartials,bestAppendingLength,isRemovedTip,isUpDown=True)
			bestBottomLength=estimateBranchLengthWithDerivative(midTopVector,downVect,fromTipC=fromTip1)
			newMidVector=mergeVectors(upVect,bestTopLength,False,downVect,bestBottomLength,fromTip1,isUpDown=True)
			appendingCost=appendProbNode(newMidVector,removedPartials,isRemovedTip,bestAppendingLength)
			if compensanteForBranchLengthChange: #if wanted, do a more thorough examination of the appending cost, taking into account a possible change in branch length of the branch on which to be appended.
				initialCost=appendProbNode(upVect,downVect,fromTip1,distance)
				newPartialCost=appendProbNode(upVect,downVect,fromTip1,bestBottomLength+bestTopLength)
				optimizedScore=appendingCost+newPartialCost-initialCost
			else:
				optimizedScore=appendingCost
			#TODO correct new placement likelihood by the change in HnZ
			#branch lengths are bestTopLength, bestBottomLength, bestAppendingLength, 
			if HnZ:
				parentNode0=up[t1]
				if bestTopLength>effectivelyNon0BLen and bestBottomLength>effectivelyNon0BLen:
					if bestAppendingLength>effectivelyNon0BLen:
						addendum0=getHnZ(2) - getHnZ(1)
					else:
						addendum0=getHnZ(nDesc0[children[node][child]]+1) - getHnZ(nDesc0[children[node][child]])
				elif bestBottomLength>effectivelyNon0BLen:
					while (dist[parentNode0]<=effectivelyNon0BLen) and up[parentNode0]!=None:
						parentNode0=up[parentNode0]
					if parentNode0==originalParent0:
						addendum0=float("-inf")
					else:
						if bestAppendingLength>effectivelyNon0BLen:
							addendum0=getHnZ(nDesc0[parentNode0]+1) - getHnZ(nDesc0[parentNode0])
						else:
							addendum0=getHnZ(nDesc0[parentNode0]+nDesc0[children[node][child]]) - ( getHnZ(nDesc0[children[node][child]]) + getHnZ(nDesc0[parentNode0]) )
				elif bestTopLength>effectivelyNon0BLen:
					if t1==originalParent0:
						addendum0=float("-inf")
					else:
						if bestAppendingLength>effectivelyNon0BLen:
							addendum0=getHnZ(nDesc0[t1]+1) - getHnZ(nDesc0[t1])
						else:
							addendum0=getHnZ(nDesc0[t1]+nDesc0[children[node][child]]) - ( getHnZ(nDesc0[children[node][child]]) + getHnZ(nDesc0[t1]) )
				else:
					while (dist[parentNode0]<=effectivelyNon0BLen) and up[parentNode0]!=None:
						parentNode0=up[parentNode0]
					if parentNode0==originalParent0 or t1==originalParent0:
						addendum0=float("-inf")
					else:
						if bestAppendingLength>effectivelyNon0BLen:
							addendum0=getHnZ(nDesc0[parentNode0]+nDesc0[t1]+1) - (getHnZ(nDesc0[parentNode0]) + getHnZ(nDesc0[t1]))
						else:
							addendum0=getHnZ(nDesc0[parentNode0]+nDesc0[t1]+nDesc0[children[node][child]]) - ( getHnZ(nDesc0[children[node][child]]) + getHnZ(nDesc0[parentNode0]) + getHnZ(nDesc0[t1]) )
				#if (optimizedScore+addendum0)>=bestScore:
				#	bestVect=[addendum0,optimizedScore,bestBottomLength,bestTopLength,bestAppendingLength,dist[t1],nDesc0[t1],nDesc0[parentNode0],nDesc0[children[node][child]]]
				optimizedScore+=addendum0

			if optimizedScore>=bestScore:
				bestNode=t1
				bestScore=optimizedScore
				bestBranchLengths=(bestTopLength,bestBottomLength,bestAppendingLength)
				bestRemovedPartials=removedPartials

			if aBayesPlusOn and sequentialSearch:
				#check that the placement location is effectively different from the original node
				differentNode=True
				topNode=node
				if t1==topNode:
					differentNode=False
				if (not bestBottomLength):
					while (dist[topNode]<=effectivelyNon0BLen) and (up[topNode]!=None):
						topNode=up[topNode]
					if t1==topNode:
						differentNode=False
				if t1==children[node][1-child]:
					differentNode=False
				#check that placement is not redundant
				if (not bestTopLength):
					differentNode=False
				if dist[t1]<=effectivelyNon0BLen:
					differentNode=False
				# check if this is a root placement
				if (not rootAlreadyConsidered) and (not bestTopLength):
					topNode=up[t1]
					while (dist[topNode]<=effectivelyNon0BLen) and (up[topNode]!=None):
						topNode=up[topNode]
					if up[topNode]==None:
						rootAlreadyConsidered=True
						listofLKcosts.append(optimizedScore)
						if networkOutput:
							listOfProbableNodes.append(topNode)
				elif differentNode: #add placement to the list of legit ones
					listofLKcosts.append(optimizedScore)
					if networkOutput :
						listOfProbableNodes.append(t1)

	if aBayesPlusOn and sequentialSearch:
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
		totalTimeFindingParent[0]+=(time()-timeStartParentTopology)
		return bestNode, bestScore, bestBranchLengths, finalListOfNodes, support, bestRemovedPartials

	else:
		totalTimeFindingParent[0]+=(time()-timeStartParentTopology)
		#TODO just for debugging
		#if HnZ and bestScore>originalLK:
		#	print("Details about best found node:")
		#	print(bestVect)
		return bestNode, bestScore, bestBranchLengths, [], None, bestRemovedPartials


#function traversing the tree to find the best new root.
def findBestRoot(tree,root,strictTopologyStopRules=strictTopologyStopRules,allowedFailsTopology=allowedFailsTopology,thresholdLogLKtopology=thresholdLogLKtopology,aBayesPlusOn=aBayesPlusOn):
	up=tree.up
	children=tree.children
	mutations=tree.mutations
	minorSequences=tree.minorSequences
	dist=tree.dist
	probVect=tree.probVect
	bestNode=root
	nodesToVisit=[]
	lastLK, distance, failedPasses = 0.0, 0.0, 0
	#best relative lk cost with respect to original root
	bestLKdiff=0.0
	bestNodes={}
	bestNodes[root]=0.0
	#bestNodes=[(root,0.0)]

	#initialize the tree search starting with the children of the root
	if children[root]: # case there is only one sample outside of the subtree doesn't need to be considered
		child1=children[root][0]
		child2=children[root][1]
		vectUp1=probVect[child2]
		if mutations[child2]:
			vectUp1=passGenomeListThroughBranch(vectUp1,mutations[child2],dirIsUp=True)
		vectUp2=probVect[child1]
		if mutations[child1]:
			vectUp2=passGenomeListThroughBranch(vectUp2,mutations[child1],dirIsUp=True)
		originalLKcost=findProbRoot(probVect[root],node=root,mutations=mutations,up=up)
		#originalLKcost=findProbRoot(probVect[root])
		isTip2=(len(children[child2])==0 and len(minorSequences[child2])==0)
		isTip1=(len(children[child1])==0 and len(minorSequences[child1])==0)
		#print(vectUp1)
		#print()
		#print(vectUp2)
		#print()
		#print(probVect[root])
		#print()
		try:
			newLower, Lkcontribution=mergeVectors(vectUp1,dist[child2],isTip2,vectUp2,dist[child1],isTip1,returnLK=True,numMinor1=len(minorSequences[child2]),numMinor2=len(minorSequences[child1]))
		except:
				print("Error in merging 0")
				print(vectUp1)
				exit()
		originalLKcost+=Lkcontribution
		#print("initial node "+str(root)+" originalLKcost "+str(originalLKcost)+" Lkcontribution "+str(Lkcontribution)+" originalLKcost-Lkcontribution "+str(originalLKcost-Lkcontribution)+" mutations:")
		#print(mutations[root])
		#print("children: ")
		#print(children[root])
		#exit()
		if mutations[child1]:
			vectUp1=passGenomeListThroughBranch(vectUp1,mutations[child1],dirIsUp=False)
		# the first LK is the likelihood cost of merges associated with the original root, these will need to be subtracted
		# the second likelihood (initialized here at 0) is the relative LK benefit of the new root - this is passed on so that one can assess if there has been an inprovement with the next branch over the previous one.
		if children[child1]:
			nodesToVisit.append((child1,vectUp1,dist[child1]+dist[child2],isTip2,len(minorSequences[child2]),originalLKcost,bestLKdiff,0))
		if mutations[child2]:
			vectUp2=passGenomeListThroughBranch(vectUp2,mutations[child2],dirIsUp=False)
		if children[child2]:
			nodesToVisit.append((child2,vectUp2,dist[child2]+dist[child1],isTip1,len(minorSequences[child1]),originalLKcost,bestLKdiff,0))

	#now visit nodes, and for each node visited, attempt to re-root at the two branches below it.
	nodesVisitedRoot=1
	while nodesToVisit:
		nodesVisitedRoot+=1
		t1,passedPartials,distance,isTip,numMinor,LKtoRemove,lastLK,failedPasses=nodesToVisit.pop()
		childs=[children[t1][0],children[t1][1]]
		probVects=[probVect[childs[0]],probVect[childs[1]]]
		dists=[dist[childs[0]],dist[childs[1]]]
		isTips=[]
		numMinors=[len(minorSequences[childs[0]]),len(minorSequences[childs[1]])]
		for i in range(2):
			if mutations[childs[i]]:
				probVects[i]=passGenomeListThroughBranch(probVects[i],mutations[childs[i]],dirIsUp=True)
			isTips.append((len(children[childs[i]])==0 and len(minorSequences[childs[i]])==0))
		newLKtoRemove=LKtoRemove
		try:
			newLower, Lkcontribution=mergeVectors(probVects[0],dists[0],isTips[0],probVects[1],dists[1],isTips[1],returnLK=True,numMinor1=numMinors[0],numMinor2=numMinors[1])
		except:
			print("Error in merging 1")
			print(probVects)
			exit()
		newLKtoRemove+=Lkcontribution
		#print("node "+str(t1)+" newLKtoRemove "+str(newLKtoRemove))
		#for each child node, assess how good the rooting is
		for i in range(2):
			try:
				upVect, lk=mergeVectors(probVects[1-i],dists[1-i],isTips[1-i],passedPartials,distance,isTip,returnLK=True,numMinor1=numMinors[1-i],numMinor2=numMinor)
			except:
				print("Error in merging 2")
				print(probVects)
				exit()
			newLKtoRemoveToPass=newLKtoRemove-lk
			try:
				newRootProbVect, lkRoot=mergeVectors(upVect,dists[i]/2,False,probVects[i],dists[i]/2,isTips[i],returnLK=True,numMinor1=0,numMinor2=numMinors[i])
			except:
				print("Error in merging 3")
				print(probVects)
				exit()
			#rootProbLK=findProbRoot(newRootProbVect)
			rootProbLK=findProbRoot(newRootProbVect,node=t1,mutations=mutations,up=up)
			score=rootProbLK+lkRoot+lk-newLKtoRemove
			#print("child "+str(i)+" is node "+str(childs[i])+" score "+str(score)+" mutations ")
			#print(mutations[childs[i]])
			failedPassesNew=failedPasses
			if score>bestLKdiff:
				shorten(upVect)
				bestLKdiff=score
				bestNode=childs[i]
				failedPassesNew=0
			elif score<(lastLK-thresholdLogLKconsecutivePlacement):
				failedPassesNew+=1
			if score>=(bestLKdiff-thresholdLogLKoptimizationTopology):
				#bestNodes.append((childs[i],score))
				bestNodes[childs[i]]=score

			#assess if to pass the search onto the children nodes
			traverseChildren=False
			if children[childs[i]]:
				if strictTopologyStopRules:
					if failedPassesNew<=allowedFailsTopology and score>(bestLKdiff-thresholdLogLKtopology):
						traverseChildren=True
				else:
					if failedPassesNew<=allowedFailsTopology or score>(bestLKdiff-thresholdLogLKtopology):
						traverseChildren=True
			if traverseChildren:
				if mutations[childs[i]]:
					vectToPass=passGenomeListThroughBranch(upVect,mutations[childs[i]],dirIsUp=False)
					shorten(vectToPass)
				else:
					vectToPass=upVect
				nodesToVisit.append((childs[i],vectToPass,dists[i],False,0,newLKtoRemoveToPass,score,failedPassesNew))

			#if childs[i]==6:
			#	exit()

	#print("Nodes visited looking for the best rooting: "+str(nodesVisitedRoot))
	#exit()

	#print(bestNode)
	#re-root the tree.
	if bestNode!=root:
		#print(bestNodes)
		#reassign the lkcost of 0 for the old root in the list bestNodes to one of its children (the one that is below it in the newly rooted tree).
		rootChild=bestNode
		nodesToInvert=[]
		while up[rootChild]!=root:
			rootChild=up[rootChild]
			#print(rootChild)
			if up[rootChild]!=root:
				nodesToInvert.append(rootChild)
		if rootChild==children[root][0]:
			sibling=children[root][1]
		else:
			sibling=children[root][0]
		#bestNodes[0]=(sibling,0.0)
		bestNodes[sibling]=bestNodes.pop(root)

		#re-assign scores for nodes whose orientation will change with the new rooting
		currentNode=up[bestNode]
		while nodesToInvert:
			currentNode=nodesToInvert.pop()
			if currentNode in bestNodes:
				bestNodes[up[currentNode]]=bestNodes.pop(currentNode)

		#re-root
		newRoot=reRootTree(tree,root,bestNode,reRootAtInternalNode=True)
	
		#assign the new best score in bestNodes to the new root, instead of its child.
		#for pairI in range(len(bestNodes)):
		#	if bestNodes[pairI][0]==bestNode:
		#		bestScore=bestNodes[pairI][1]
		#		bestNodes[pairI]=(newRoot,bestScore)
		#		break
		bestNodes[newRoot]=bestNodes.pop(bestNode)

		#print(bestNode)
		
		#recalculate the genome lists following the new rooting
		reCalculateAllGenomeLists(tree,newRoot)

	else:
		newRoot=root
	
	#calculate support scores for alternative roots
	if aBayesPlusOn:
		#listofLKcosts=[]
		totSupport=0.0
		tree.rootSupport=[]
		for i in range(len(up)):
			tree.rootSupport.append(None)
		for currentNode in bestNodes:
			bestNodes[currentNode]=math.exp(bestNodes[currentNode])
			totSupport+=bestNodes[currentNode]
		#for i in range(len(bestNodes)):
		#	listofLKcosts.append(math.exp(bestNodes[i][1]))
		#	totSupport+=listofLKcosts[i]
		for currentNode in bestNodes:
			bestNodes[currentNode]=bestNodes[currentNode]/totSupport
			if bestNodes[currentNode]>=minBranchSupport:
				tree.rootSupport[currentNode]=bestNodes[currentNode]
		#for i in range(len(bestNodes)):
		#	listofLKcosts[i]=listofLKcosts[i]/totSupport
		#	if listofLKcosts[i]>=minBranchSupport:
		#		tree.rootSupport[bestNodes[i][0]]=listofLKcosts[i]
		#print(bestNodes)

	return newRoot


numMinorsFound=[0]
#function to find the best node in the tree where to append the new sample; traverses the tree and tries to append the sample at each node and mid-branch nodes, 
# but stops traversing when certain criteria are met.
#TODO account for HnZ modifiers in placement search
def findBestParentForNewSample(tree,root,diffs,sample):
	up=tree.up
	children=tree.children
	probVectUpRight=tree.probVectUpRight
	probVectUpLeft=tree.probVectUpLeft
	mutations=tree.mutations
	minorSequences=tree.minorSequences
	dist=tree.dist
	probVect=tree.probVect
	probVectTotUp=tree.probVectTotUp
	#TODO
	nDesc0=tree.nDesc0
	bestNodes=[]
	bestNode=root
	bestBranchLengths=(False,False,oneMutBLen)
	if mutations[root]:
		diffs=passGenomeListThroughBranch(diffs,mutations[root])
	bestDiffs=diffs
	if not children[root]: #check if the new leaf is strictly less informative than already placed leaf
		#TODO do not collapse less informative sequences also in case of HnZ
		if errorRateSiteSpecificFile or errorRateFixed or estimateErrorRate or estimateSiteSpecificErrorRate or supportFor0Branches or HnZ:
			comparison=isMinorSequence(probVect[root],diffs,onlyFindIdentical=True)
		else:
			comparison=isMinorSequence(probVect[root],diffs)
		if comparison==1:
			minorSequences[root].append(sample)
			#TODO
			nDesc0[root]+=1
			numMinorsFound[0]+=1
			return root, 1.0, None, diffs
		elif comparison==2:
			totalMissedMinors[0]+=1
	rootVect=rootVector(probVect[root],False,False,tree,root)
	bestLKdiff=appendProbNode(rootVect,diffs,True,oneMutBLen)
	#TODO
	if HnZ:
		bestLKdiff+=getHnZ(nDesc0[root]+1) - getHnZ(nDesc0[root])
	nodesToVisit=[]
	for child in children[root]:
		diffsChild=diffs
		if mutations[child]:
			diffsChild=passGenomeListThroughBranch(diffs,mutations[child])
		nodesToVisit.append((child,bestLKdiff,0,diffsChild))
	while nodesToVisit:
		t1,parentLK,failedPasses,diffs=nodesToVisit.pop()
		if not children[t1]: #check if the new leaf is strictly less informative than already placed leaf
			#TODO do not collapse less informative sequences also in case of HnZ
			if errorRateSiteSpecificFile or errorRateFixed or estimateErrorRate or estimateSiteSpecificErrorRate or supportFor0Branches or HnZ:
				comparison=isMinorSequence(probVect[t1],diffs,onlyFindIdentical=True)
			else:
				comparison=isMinorSequence(probVect[t1],diffs)
			if comparison==1:
				minorSequences[t1].append(sample)
				#TODO
				if HnZ:
					nDesc0[t1]+=1
					if (dist[t1]<=effectivelyNon0BLen):
						parent0=t1
						while (dist[parent0]<=effectivelyNon0BLen) and up[parent0]!=None:
							parent0=up[parent0]
							nDesc0[parent0]+=1
				numMinorsFound[0]+=1
				if (not onlyNambiguities) and usingErrorRate:
					updateProbVectTerminalNode(probVect[t1],minorSequences[t1])
				if usingErrorRate:
					nodeList=[(t1,2)]
					if up[t1]!=None:
						if t1==children[up[t1]][0]:
							nodeList.append((up[t1],0))
						else:
							nodeList.append((up[t1],1))
					updatePartials(tree,nodeList)
				return t1, 1.0, None, diffs
			elif comparison==2:
				totalMissedMinors[0]+=1

		if dist[t1] and up[t1]!=None: # try first placing as a descendant of the mid-branch point of the branch above the current node.
			LKdiff=appendProbNode(probVectTotUp[t1],diffs,True,oneMutBLen)
			#TODO
			if HnZ:
				LKdiff+=getHnZ(2) - getHnZ(1)
			if LKdiff>=bestLKdiff:
				shorten(diffs)
				bestLKdiff=LKdiff
				bestNode=t1
				failedPasses=0
				bestNodes.append((t1,LKdiff,diffs))
				bestDiffs=diffs
			elif LKdiff>bestLKdiff-thresholdLogLKoptimization:
				bestNodes.append((t1,LKdiff,diffs))
			if LKdiff<(parentLK-thresholdLogLKconsecutivePlacement): #placement at current node is considered failed if placement likelihood is not improved by a certain margin compared to best placement so far for the nodes above it.
				failedPasses+=1
		else:
			LKdiff=parentLK
		#keep trying to place at children nodes, unless placement has failed too many times already or the likelihood cost suboptimality is above a certain threshold
		if strictStopRules:
			if failedPasses<=allowedFails and LKdiff>(bestLKdiff-thresholdLogLK): 
				for c in children[t1]:
					diffsChild=diffs
					if mutations[c]:
						diffsChild=passGenomeListThroughBranch(diffs,mutations[c])
					nodesToVisit.append((c,LKdiff,failedPasses,diffsChild))
		else:
			if failedPasses<=allowedFails or LKdiff>(bestLKdiff-thresholdLogLK):
				for c in children[t1]:
					diffsChild=diffs
					if mutations[c]:
						diffsChild=passGenomeListThroughBranch(diffs,mutations[c])
					nodesToVisit.append((c,LKdiff,failedPasses,diffsChild))
	#Initial exploration is finished.
	#Now, for each branch within threshold likelihood distance from the best found, optimize branch lengths.
	#Use optimized scores to select final best branch
	if bestNode!=root:
		bestBranchLengths=(dist[bestNode]/2,dist[bestNode]/2,oneMutBLen)
	bestScore=bestLKdiff
	compensanteForBranchLengthChange=True
	for nodePair in bestNodes:
		score=nodePair[1]
		if score>=bestLKdiff-thresholdLogLKoptimization:
			node=nodePair[0]
			#optimize branch lengths of appendage
			if node==children[up[node]][0]:
				upVect=probVectUpRight[up[node]]
			else:
				upVect=probVectUpLeft[up[node]]
			if mutations[node]:
				upVect=passGenomeListThroughBranch(upVect,mutations[node])
			diffs=nodePair[-1]
			isTip=(len(children[node])==0) and (len(minorSequences[node])==0)
			bestAppendingLength=estimateBranchLengthWithDerivative(probVectTotUp[node],diffs,fromTipC=True)
			#now optimize appending location
			midLowerVector=mergeVectors(probVect[node],dist[node]/2,isTip,diffs,bestAppendingLength,True)
			bestTopLength=estimateBranchLengthWithDerivative(upVect,midLowerVector,fromTipC=False)
			midTopVector=mergeVectors(upVect,bestTopLength,False,diffs,bestAppendingLength,True,isUpDown=True)
			try:
				bestBottomLength=estimateBranchLengthWithDerivative(midTopVector,probVect[node],fromTipC=isTip)
			except:
				print("Problem with estimateBranchLengthWithDerivative()")
				exit()
			newMidVector=mergeVectors(upVect,bestTopLength,False,probVect[node],bestBottomLength,isTip,isUpDown=True)
			appendingCost=appendProbNode(newMidVector,diffs,True,bestAppendingLength)
			if compensanteForBranchLengthChange: #if wanted, do a more thorough examination of the appending cost, taking into account a possible change in branch length of the branch on which to be appended.
				initialCost=appendProbNode(upVect,probVect[node],isTip,dist[node])
				newPartialCost=appendProbNode(upVect,probVect[node],isTip,bestBottomLength+bestTopLength)
				optimizedScore=appendingCost+newPartialCost-initialCost
			else:
				optimizedScore=appendingCost

			#TODO
			if HnZ:
				bestLKdiff+=getHnZ(nDesc0[root]+1) - getHnZ(nDesc0[root])

				if bestTopLength>effectivelyNon0BLen and bestBottomLength>effectivelyNon0BLen:
					optimizedScore+=getHnZ(2) - getHnZ(1)
				elif bestBottomLength>effectivelyNon0BLen:
					parentNode0=up[node]
					while (dist[parentNode0]<=effectivelyNon0BLen) and up[parentNode0]!=None:
						parentNode0=up[parentNode0]
						optimizedScore+=getHnZ(nDesc0[parentNode0]+1) - getHnZ(nDesc0[parentNode0])
				elif bestTopLength>effectivelyNon0BLen:
					optimizedScore+=getHnZ(nDesc0[node]+1) - getHnZ(nDesc0[node])

			if optimizedScore>=bestScore:
				bestNode=node
				bestScore=optimizedScore
				bestBranchLengths=(bestTopLength,bestBottomLength,bestAppendingLength)
				bestDiffs=diffs

	return bestNode, bestScore, bestBranchLengths, bestDiffs


# Change the genome lists of the node and edescendants after making the node reference. Also change nDesc of the ancestor nodes.
def makeNodeReference(tree,node,oldValue=0):
	up=tree.up
	children=tree.children
	probVectUpRight=tree.probVectUpRight
	probVectUpLeft=tree.probVectUpLeft
	mutations=tree.mutations
	dist=tree.dist
	probVect=tree.probVect
	probVectTotUp=tree.probVectTotUp
	nDesc=tree.nDesc
	numRefs[0]+=1
	#reduce nDesc of parent nodes according to oldValue
	if oldValue:
		pNode=up[node]
		while pNode!=None:
			nDesc[pNode]-=oldValue
			if mutations[pNode]:
				break
			pNode=up[pNode]
	#define node.mutations
	pos=0
	for entry in probVect[node]:
		if entry[0]<4:
			pos+=1
			mutations[node].append((pos,entry[1],entry[0]))
		elif entry[0]==6:
			pos+=1
		else:
			pos=entry[1]
	#update node's genome lists accordingly
	probVect[node]=passGenomeListThroughBranch(probVect[node],mutations[node])
	shorten(probVect[node])
	if dist[node] and up[node]!=None:
		probVectTotUp[node]=passGenomeListThroughBranch(probVectTotUp[node],mutations[node])
		shorten(probVectTotUp[node])
	probVectUpRight[node]=passGenomeListThroughBranch(probVectUpRight[node],mutations[node])
	shorten(probVectUpRight[node])
	probVectUpLeft[node]=passGenomeListThroughBranch(probVectUpLeft[node],mutations[node])
	shorten(probVectUpLeft[node])
	#now traverse descendant nodes to update their genome lists and mutation lists
	nodesToVisit=[children[node][0],children[node][1]]
	while nodesToVisit:
		newNode=nodesToVisit.pop()
		if mutations[newNode]:
			mutations[newNode]=mergeMutationLists(mutations[node],mutations[newNode],downward=True)
		else:
			probVect[newNode]=passGenomeListThroughBranch(probVect[newNode],mutations[node])
			shorten(probVect[newNode])
			if dist[newNode]:
				probVectTotUp[newNode]=passGenomeListThroughBranch(probVectTotUp[newNode],mutations[node])
				shorten(probVectTotUp[newNode])
			if children[newNode]:
				probVectUpRight[newNode]=passGenomeListThroughBranch(probVectUpRight[newNode],mutations[node])
				shorten(probVectUpRight[newNode])
				probVectUpLeft[newNode]=passGenomeListThroughBranch(probVectUpLeft[newNode],mutations[node])
				shorten(probVectUpLeft[newNode])
				nodesToVisit.append(children[newNode][0])
				nodesToVisit.append(children[newNode][1])


# number of non-reference nucleotides in the genome list.
def numNon4(probVect):
	num=0
	for entry in probVect:
		if entry[0]<4:
			num+=1
	return num


#we know that sample "sample", with partials "newPartials", is best placed near a node resulting in logLK contribution of newChildLK
# now explore exactly which position near node is best for placement (direct child, sibling or sibling of existing child), and which are the best branch lengths,
#then add the sample at that position of the tree, and update all the internal probability vectors.
#UNNECESSARY? Could probably just be replaced by the more general placeSubtreeOnTree().
#TODO update nDesc0 as more samples are added
def placeSampleOnTree(tree,node,newPartials,sample,newChildLK, bestUpLength, bestDownLength, bestAppendingLength,pseudoMutCounts):
	up=tree.up
	children=tree.children
	probVectUpRight=tree.probVectUpRight
	probVectUpLeft=tree.probVectUpLeft
	mutations=tree.mutations
	dist=tree.dist
	probVect=tree.probVect
	probVectTotUp=tree.probVectTotUp
	nDesc=tree.nDesc
	minorSequences=tree.minorSequences
	name=tree.name
	#TODO
	nDesc0=tree.nDesc0
	tryNewRoot=False
	if newChildLK<-0.01:
		sumChildLKs[0]+=newChildLK
		numChildLKs[0]+=1
	if up[node]==None:
		tryNewRoot=True
		rootNewPartials=newPartials
		totRoot=rootVector(probVect[node],False,False,tree,node)
		bestAppendingLength=estimateBranchLengthWithDerivative(totRoot,newPartials,fromTipC=True)
		root=node
		newChildLK=appendProbNode(totRoot,newPartials,True,bestAppendingLength)
	else:
		if children[up[node]][0]==node:
			child=0
			vectUp=probVectUpRight[up[node]]
		else:
			child=1
			vectUp=probVectUpLeft[up[node]]
		if mutations[node]:
			vectUp=passGenomeListThroughBranch(vectUp,mutations[node],dirIsUp=False)
		if not bestUpLength:
			pNode=up[node]
			while (not dist[pNode]) and (up[pNode]!=None):
				pNode=up[pNode]
			if up[pNode]==None:
				root=pNode
				tryNewRoot=True
				if (not bestDownLength) or (bestDownLength>1.01*dist[node]) or (bestDownLength<0.99*dist[node]):
					dist[node]=bestDownLength
					nodeList=[(node,2),(up[node],child)]
					updatePartials(tree,nodeList)
			if tryNewRoot:
				pNode=up[node]
				rootNewPartials=newPartials
				if mutations[node]:
					rootNewPartials=passGenomeListThroughBranch(newPartials,mutations[node],dirIsUp=True)
				while (not dist[pNode]) and (up[pNode]!=None):
					if mutations[pNode]:
						rootNewPartials=passGenomeListThroughBranch(rootNewPartials,mutations[pNode],dirIsUp=True)
					pNode=up[pNode]
	isTip=(len(children[node])==0) and (len(minorSequences[node])==0)
	#in case of best placement as a descendant appended exactly at the root node, attempt also to create new root
	if tryNewRoot:
		node=root
		probVectRoot=probVect[node]
		if mutations[node]:
			probVectRoot=passGenomeListThroughBranch(probVectRoot,mutations[node],dirIsUp=True)
		probOldRoot = findProbRoot(probVectRoot)
		rootUpLeft=rootVector(probVect[node],bestAppendingLength/2,isTip,tree,node)
		bestRightLength=estimateBranchLengthWithDerivative(rootUpLeft,rootNewPartials,fromTipC=True)
		rootUpRight=rootVector(rootNewPartials,bestRightLength,True,tree,node)
		bestLeftLength=estimateBranchLengthWithDerivative(rootUpRight,probVect[node],fromTipC=isTip)
		secondBranchLengthOptimizationRound=True
		if secondBranchLengthOptimizationRound: #if wanted, do a second round of branch length optimization
			rootUpLeft=rootVector(probVect[node],bestLeftLength,isTip,tree,node)
			bestRightLength=estimateBranchLengthWithDerivative(rootUpLeft,rootNewPartials,fromTipC=True)
			rootUpRight=rootVector(rootNewPartials,bestRightLength,True,tree,node)
			bestLeftLength=estimateBranchLengthWithDerivative(rootUpRight,probVect[node],fromTipC=isTip)
		probVectRoot,probRoot = mergeVectors(probVect[node],bestLeftLength,isTip,rootNewPartials,bestRightLength,True,returnLK=True,numMinor1=len(minorSequences[node]),numMinor2=0)
		probVectRootUp=probVectRoot
		if mutations[node]:
			probVectRootUp=passGenomeListThroughBranch(probVectRoot,mutations[node],dirIsUp=True)
		probRoot+= findProbRoot(probVectRootUp)
		parentLKdiff=probRoot-probOldRoot
		if parentLKdiff<=newChildLK: #best is just placing as descendant of the root
			bestRightLength=bestAppendingLength
			bestLeftLength=False
			probVectRoot=mergeVectors(probVect[node],bestLeftLength,isTip,rootNewPartials,bestRightLength,True)
			rootUpRight=rootVector(rootNewPartials,bestRightLength,True,tree,node)
		#now add new root to the tree
		newRoot=len(up)
		tree.addNode()
		if probVectRoot==None:
			print("Issue with new root probVect inside placeSampleOnTree()")
			raise Exception("exit")
		shorten(probVectRoot)
		probVect[newRoot]=probVectRoot
		shorten(rootUpRight)
		probVectUpRight[newRoot]=rootUpRight
		probVectUpLeft[newRoot]=rootVector(probVect[node],bestLeftLength,isTip,tree,node)
		shorten(probVectUpLeft[newRoot])
		mutations[newRoot]=mutations[node]
		mutations[node]=[]
		up[node]=newRoot
		dist[node]=bestLeftLength
		#TODO 
		if HnZ:
			if bestLeftLength>effectivelyNon0BLen:
				nDesc0[newRoot]=2
			else:
				nDesc0[newRoot]=nDesc0[node]+1
		children[newRoot].append(node)
		if children[node]:
			nDesc[newRoot]+=nDesc[node]
		if bestLeftLength:
			nDesc[newRoot]+=1
		if bestRightLength:
			nDesc[newRoot]+=1
		newNode=len(up)
		tree.addNode()
		name[-1]=sample
		dist[-1]=bestRightLength
		if bestRightLength>0.01 and (not warnedBLen[0]):
			warnedBLen[0]=True
			print("\n WARNING!!!!!!!!!!!\n\n Found branch of length "+str(bestRightLength)+" ; at high divergence, MAPLE will struggle both in terms of accuracy and computational demand, so it might make sense to use a traditional phylogenetic approach.\n\n End of warning.\n")
		minorSequences[newNode]=[]
		up[newNode]=newRoot
		children[newRoot].append(newNode)
		shorten(rootNewPartials)
		probVect[newNode]=rootNewPartials
		mutations[newNode]=[]
		if bestRightLength:
			probVectTotUp[newNode]=mergeVectors(probVectUpLeft[newRoot],bestRightLength/2,False,rootNewPartials,bestRightLength/2,True,isUpDown=True)
			shorten(probVectTotUp[newNode])
		nodeList=[(node,2)]
		updatePartials(tree,nodeList)
		if (not mutations[newRoot]) and nDesc[newRoot]>=maxNumDescendantsForMATClade and numNon4(probVect[newRoot])>minNumNon4:
			makeNodeReference(tree,newRoot)
		return newRoot

	#in all other cases (not attempting to add a new root) create a new internal node in the tree and add sample as a descendant.
	newInternalNode=len(up)
	tree.addNode()
	children[up[node]][child]=newInternalNode
	up[newInternalNode]=up[node]
	children[newInternalNode].append(node)
	up[node]=newInternalNode
	dist[node]=bestDownLength
	#TODO 
	if HnZ:
		if bestDownLength>effectivelyNon0BLen:
			nDesc0[newInternalNode]=2
		else:
			nDesc0[newInternalNode]=nDesc0[node]+1
	passUpMutations=False
	#use descendantsToPass to update nDesc of parent nodes (up to the next reference met), and decide if and which internal nodes will become references
	if mutations[node] and (not bestDownLength):
		mutations[newInternalNode]=mutations[node]
		nDesc[newInternalNode]=nDesc[node]
		if bestAppendingLength:
			nDesc[newInternalNode]+=1
		mutations[node]=[]
		descendantsToPass=0
	else:
		if mutations[node]:
			#use this flag to remember to pass genome lists through node.mutations
			passUpMutations=True
			nDesc[newInternalNode]=1
			descendantsToPass=1
		else:
			if children[node]:
				nDesc[newInternalNode]=nDesc[node]
			else:
				nDesc[newInternalNode]=0
			descendantsToPass=0
			if bestDownLength:
				descendantsToPass+=1
				nDesc[newInternalNode]+=1
		mutations[newInternalNode]=[]
		if bestAppendingLength:
			nDesc[newInternalNode]+=1
			descendantsToPass+=1
		if bestDownLength and (not bestUpLength):
			descendantsToPass-=1
			
	#newNode=Tree(name=sample,dist=bestAppendingLength)
	newNode=len(up)
	tree.addNode()
	name[-1]=sample
	dist[-1]=bestAppendingLength
	if bestAppendingLength>0.01 and (not warnedBLen[0]):
		warnedBLen[0]=True
		print("\n WARNING!!!!!!!!!!!\n\n Found branch of length "+str(bestAppendingLength)+" ; at high divergence, MAPLE will struggle both in terms of accuracy and computational demand, so it might make sense to use a more traditional phylogenetic approach (e.g. FastTree, IQtree or RAxML).\n\n End of warning.\n")
	minorSequences[newNode]=[]
	up[newNode]=newInternalNode
	children[newInternalNode].append(newNode)
	dist[newInternalNode]=bestUpLength
	#TODO 
	if HnZ:
		parent0=newInternalNode
		while up[parent0]!=None and (dist[parent0]<=effectivelyNon0BLen):
			parent0=up[parent0]
			nDesc0[parent0]+=1

	probVect[newNode]=newPartials
	if passUpMutations:
		probVect[newNode]=passGenomeListThroughBranch(probVect[newNode],mutations[node],dirIsUp=True)
	shorten(probVect[newNode])
	mutations[newNode]=[]
	probVect[newInternalNode]=mergeVectors(probVect[node],bestDownLength,isTip,newPartials,bestAppendingLength,True,isUpDown=False)
	if passUpMutations:
		probVect[newInternalNode]=passGenomeListThroughBranch(probVect[newInternalNode],mutations[node],dirIsUp=True)
	shorten(probVect[newInternalNode])
	probVectUpRight[newInternalNode]=mergeVectors(vectUp,bestUpLength,False,newPartials,bestAppendingLength,True,isUpDown=True)
	if passUpMutations:
		probVectUpRight[newInternalNode]=passGenomeListThroughBranch(probVectUpRight[newInternalNode],mutations[node],dirIsUp=True)
	shorten(probVectUpRight[newInternalNode])
	probVectUpLeft[newInternalNode]=mergeVectors(vectUp,bestUpLength,False,probVect[node],bestDownLength,isTip,isUpDown=True)
	if passUpMutations:
		probVectUpLeft[newInternalNode]=passGenomeListThroughBranch(probVectUpLeft[newInternalNode],mutations[node],dirIsUp=True)
	shorten(probVectUpLeft[newInternalNode])
	if probVect[newInternalNode]==None:
		print("Problem in placeSampleOnTree(), probVect is None")
		raise Exception("exit")
	if probVectUpRight[newInternalNode]==None:
		print("Problem in placeSampleOnTree(), probVectUpRight is None")
		raise Exception("exit")
	if probVectUpLeft[newInternalNode]==None:
		print("Problem in placeSampleOnTree(), probVectUpLeft is None")
		raise Exception("exit")
	if bestUpLength:
		probVectTotUp[newInternalNode]=mergeVectors(vectUp,bestUpLength/2,False,probVect[newInternalNode],bestUpLength/2,False,isUpDown=True)
		if passUpMutations:
			probVectTotUp[newInternalNode]=passGenomeListThroughBranch(probVectTotUp[newInternalNode],mutations[node],dirIsUp=True)
		shorten(probVectTotUp[newInternalNode])
	else:
		probVectTotUp[newInternalNode]=None
	if bestAppendingLength:
		probVectTotUp[newNode]=mergeVectors(probVectUpLeft[newInternalNode],bestAppendingLength/2,False,newPartials,bestAppendingLength/2,True,isUpDown=True)
		if passUpMutations:
			probVectTotUp[newNode]=passGenomeListThroughBranch(probVectTotUp[newNode],mutations[node],dirIsUp=True)
		shorten(probVectTotUp[newNode])
		updatePesudoCounts(probVectUpLeft[newInternalNode],newPartials,pseudoMutCounts)
	else:
		probVectTotUp[newNode]=None
	if not bestDownLength:
		probVectTotUp[node]=None
	
	if descendantsToPass:
		pNode=up[newInternalNode]
		nDesc[pNode]+=descendantsToPass
		while (not mutations[pNode]):
			if nDesc[pNode]>=maxNumDescendantsForMATClade and numNon4(probVect[pNode])>minNumNon4:
				makeNodeReference(tree,pNode,oldValue=(nDesc[pNode]-descendantsToPass))
				break
			pNode=up[pNode]
			if pNode==None:
				break
			nDesc[pNode]+=descendantsToPass
	nodeList=[(node,2),(up[newInternalNode],child)]
	updatePartials(tree,nodeList)

	return None


#set all descendant nodes to dirty.
#So far thses flags are used to prevent traversing the same part of the tree multiple times.
def setAllDirty(tree,node,dirtiness=True):
	dirty=tree.dirty
	replacements=tree.replacements
	children=tree.children
	nextLeaves=[node]
	while nextLeaves:
		nextNode=nextLeaves.pop()
		dirty[nextNode]=dirtiness
		#number of times this node has been replaced through SPR during this round.
		replacements[nextNode]=0
		for c in children[nextNode]:
			nextLeaves.append(c)


#traverse the tree to optimize (only) the length of the branches using the derivative approach.
#TODO update nDesc0 when branch lengths become 0 or are not 0 anymore?
def traverseTreeToOptimizeBranchLengths(tree,root,testing=False,fastPass=False):
	up=tree.up
	children=tree.children
	probVectUpRight=tree.probVectUpRight
	probVectUpLeft=tree.probVectUpLeft
	mutations=tree.mutations
	dist=tree.dist
	probVect=tree.probVect
	minorSequences=tree.minorSequences
	dirty=tree.dirty
	nDesc0=tree.nDesc0
	totalLKimprovementBL=0.0
	updates=0
	nodesTraversed=0
	if children[root]:
		nodesToTraverse=[children[root][0],children[root][1]]
	else:
		return 0
	while nodesToTraverse:
		nodesTraversed+=1
		node=nodesToTraverse.pop()
		if dirty[node]:
			if node==children[up[node]][0]:
				upVect=probVectUpRight[up[node]]
				child=0
			else:
				upVect=probVectUpLeft[up[node]]
				child=1
			if mutations[node]:
				upVect=passGenomeListThroughBranch(upVect,mutations[node],dirIsUp=False)
			isTip=(len(children[node])==0) and (len(minorSequences[node])==0)
			bestLength=estimateBranchLengthWithDerivative(upVect,probVect[node],fromTipC=isTip)
			if bestLength or dist[node]:
				if testing:
					currentCost=appendProbNode(upVect,probVect[node],isTip,dist[node])
					newCost=appendProbNode(upVect,probVect[node],isTip,bestLength)
					totalLKimprovementBL+=newCost-currentCost
				if (not bestLength) or (not dist[node]) or dist[node]/bestLength>1.01 or dist[node]/bestLength<0.99:
					#TODO updating nDesc0 in case branch lengths move to 0 or from 0
					if HnZ:
						if bestLength>effectivelyNon0BLen and (dist[node]<=effectivelyNon0BLen):
							parent0=node
							#print("initial nDesc0: "+str(nDesc0[parent0]))
							while (dist[parent0]<=effectivelyNon0BLen) and up[parent0]!=None:
								parent0=up[parent0]
								#print("Old nDesc0: "+str(nDesc0[parent0]))
								nDesc0[parent0]-=(nDesc0[node]-1)
								#print("new nDesc0: "+str(nDesc0[parent0])+" after removing "+str(nDesc0[node]))
								if nDesc0[parent0]<=0:
									raise Exception("exit")
						elif dist[node]>effectivelyNon0BLen and (bestLength<=effectivelyNon0BLen):
							parent0=up[node]
							nDesc0[parent0]+=(nDesc0[node]-1)
							while (dist[parent0]<=effectivelyNon0BLen) and up[parent0]!=None:
								parent0=up[parent0]
								nDesc0[parent0]+=(nDesc0[node]-1)
					dist[node]=bestLength
					updates+=1
					if not fastPass:
						nodeList=[(node,2),(up[node],child)]
						updatePartials(tree,nodeList)
				else:
					#if inputTree=="" or largeUpdate or aBayesPlus:
					dirty[node]=False
			else:
				#if inputTree=="" or largeUpdate or aBayesPlus:
				dirty[node]=False
		for child in children[node]:
			nodesToTraverse.append(child)
	if testing:
		return totalLKimprovementBL
	else:
		return updates
		

#we know that subtree "appendedNode", with partials "newPartials", is best placed as child of "node" resulting in logLK contribution of newChildLK
# now explore exactly which position near node is best for placement (direct child, sibling or sibling of existing child), and which are the best branch lengths,
#then add the subtree at that position of the tree, and update all the internal probability vectors.
#TODO update nDesc0 as a subtrees is added to the tree.
def placeSubtreeOnTree(tree,node,newPartials,appendedNode,newChildLK,bestBranchLengths): #,parentNode,parentNodeReplacements=1
	up=tree.up
	children=tree.children
	probVectUpRight=tree.probVectUpRight
	probVectUpLeft=tree.probVectUpLeft
	mutations=tree.mutations
	dist=tree.dist
	probVect=tree.probVect
	minorSequences=tree.minorSequences
	dirty=tree.dirty
	replacements=tree.replacements
	probVectTotUp=tree.probVectTotUp
	#TODO
	nDesc0=tree.nDesc0
	bestAppendingLength=bestBranchLengths[2]
	bestUpLength=bestBranchLengths[0]
	bestDownLength=bestBranchLengths[1]
	tryNewRoot=False
	if children[up[node]][0]==node:
		child=0
		vectUp=probVectUpRight[up[node]]
	else:
		child=1
		vectUp=probVectUpLeft[up[node]]
	
	#test if we should also attempt placing the subtree as child of a new root
	if not bestUpLength:
		pNode=up[node]
		while (not dist[pNode]) and (up[pNode]!=None):
			pNode=up[pNode]
		if up[pNode]==None:
			root=pNode
			tryNewRoot=True
			if (not bestDownLength) or (bestDownLength>1.01*dist[node]) or (bestDownLength<0.99*dist[node]):
				dist[node]=bestDownLength
				nodeList=[(node,2),(up[node],child)]
				updatePartials(tree,nodeList)
		if tryNewRoot:
			pNode=up[node]
			rootNewPartials=newPartials
			if mutations[node]:
				rootNewPartials=passGenomeListThroughBranch(newPartials,mutations[node],dirIsUp=True)
			while (not dist[pNode]) and (up[pNode]!=None):
				if mutations[pNode]:
					rootNewPartials=passGenomeListThroughBranch(rootNewPartials,mutations[pNode],dirIsUp=True)
				pNode=up[pNode]
	appendedIsTip=(len(children[appendedNode])==0) and (len(minorSequences[appendedNode])==0)
			
	#in case of best placement as a descendant appended exactly at the root node, attempt also to create new root
	if tryNewRoot:
		node=root
		isTip=(len(children[node])==0) and (len(minorSequences[node])==0)
		probVectRootUp=probVect[node]
		if mutations[node]:
			probVectRootUp=passGenomeListThroughBranch(probVect[node],mutations[node],dirIsUp=True)
		probOldRoot = findProbRoot(probVectRootUp)
		rootUpLeft=rootVector(probVect[node],bestAppendingLength/2,isTip,tree,node)
		bestRightLength=estimateBranchLengthWithDerivative(rootUpLeft,rootNewPartials,fromTipC=appendedIsTip)
		rootUpRight=rootVector(rootNewPartials,bestRightLength,appendedIsTip,tree,node)
		bestLeftLength=estimateBranchLengthWithDerivative(rootUpRight,probVect[node],fromTipC=isTip)
		secondBranchLengthOptimizationRound=True
		if secondBranchLengthOptimizationRound: #if wanted, do a second round of branch length optimization
			rootUpLeft=rootVector(probVect[node],bestLeftLength,isTip,tree,node)
			bestRightLength=estimateBranchLengthWithDerivative(rootUpLeft,rootNewPartials,fromTipC=appendedIsTip)
			rootUpRight=rootVector(rootNewPartials,bestRightLength,appendedIsTip,tree,node)
			bestLeftLength=estimateBranchLengthWithDerivative(rootUpRight,probVect[node],fromTipC=isTip)
		probVectRoot,probRoot = mergeVectors(probVect[node],bestLeftLength,isTip,rootNewPartials,bestRightLength,appendedIsTip,returnLK=True,isUpDown=False,numMinor1=len(minorSequences[node]),numMinor2=len(minorSequences[appendedNode]))
		probVectRootUp=probVectRoot
		if mutations[node]:
			probVectRootUp=passGenomeListThroughBranch(probVectRoot,mutations[node],dirIsUp=True)
		probRoot+= findProbRoot(probVectRootUp)
		parentLKdiff=probRoot-probOldRoot
		if parentLKdiff<=newChildLK: #best is just placing as descendant of the root
			bestRightLength=bestAppendingLength
			bestLeftLength=False
			probVectRoot=mergeVectors(probVect[node],bestLeftLength,isTip,rootNewPartials,bestRightLength,appendedIsTip)
			rootUpRight=rootVector(rootNewPartials,bestRightLength,appendedIsTip,tree,node)
		#now add new root to the tree
		#calculate new .mutations for appendedNode by traversing up and down until MRCA is found (in this case MRCA is just the root)
		if mutations[appendedNode]:
			numRefs[0]-=1
		traverseTreeToUpdateMutationList(tree,appendedNode,node)
		if mutations[appendedNode]:
			numRefs[0]+=1
		newRoot=up[appendedNode]
		up[newRoot]=None
		dirty[newRoot]=True
		dist[newRoot]=defaultBLen
		replacements[newRoot]+=1
		if probVectRoot==None:
			print("Issue with new root probVect inside placeSubtreeOnTree()")
			raise Exception("exit")
		shorten(probVectRoot)
		probVect[newRoot]=probVectRoot
		shorten(rootUpRight)
		probVectUpRight[newRoot]=rootUpRight
		probVectUpLeft[newRoot]=rootVector(probVect[node],bestLeftLength,isTip,tree,node)
		shorten(probVectUpLeft[newRoot])
		mutations[newRoot]=mutations[node]
		mutations[node]=[]
		up[node]=newRoot
		dist[node]=bestLeftLength
		children[newRoot][0]=node
		children[newRoot][1]=appendedNode
		dist[appendedNode]=bestRightLength
		replacements[appendedNode]+=1
		nodeList=[(node,2),(appendedNode,2)]
		updatePartials(tree,nodeList)
		#TODO
		if HnZ:
			#print("Placing a new root")
			if dist[node]>effectivelyNon0BLen:
				nDesc0[newRoot]=1
			else:
				nDesc0[newRoot]=nDesc0[node]
			if dist[appendedNode]>effectivelyNon0BLen:
				nDesc0[newRoot]+=1
			else:
				nDesc0[newRoot]+=nDesc0[appendedNode]
		return newRoot

	#in all other cases (not attempting to add a new root) create a new internal node in the tree and add subtree as a descendant.
	if mutations[node]:
		vectUp=passGenomeListThroughBranch(vectUp,mutations[node])
	isTip=(len(children[node])==0) and (len(minorSequences[node])==0)
	#calculate new .mutations for appendedNode by traversing up and down until MRCA is found
	if mutations[appendedNode]:
		numRefs[0]-=1
	traverseTreeToUpdateMutationList(tree,appendedNode,node)
	if mutations[appendedNode]:
		numRefs[0]+=1
	newInternalNode=up[appendedNode]
	mutations[newInternalNode]=mutations[node]
	mutations[node]=[]
	dirty[newInternalNode]=True
	replacements[newInternalNode]+=1
	children[up[node]][child]=newInternalNode
	up[newInternalNode]=up[node]
	children[newInternalNode][0]=node
	up[node]=newInternalNode
	replacements[appendedNode]+=1
	children[newInternalNode][1]=appendedNode
	probVect[newInternalNode]=mergeVectors(probVect[node],bestDownLength,isTip,newPartials,bestAppendingLength,appendedIsTip)
	if probVect[newInternalNode]==None:
		probVectUpLeft[newInternalNode]=mergeVectors(vectUp,bestUpLength,False,probVect[node],bestDownLength,isTip,isUpDown=True)
		if probVectUpLeft[newInternalNode]==None:
			probVectUpRight[newInternalNode]=mergeVectors(vectUp,bestUpLength,False,newPartials,bestAppendingLength,appendedIsTip,isUpDown=True)
			bestDownLength=estimateBranchLengthWithDerivative(probVectUpRight[newInternalNode],probVect[node],fromTipC=isTip)
			probVectUpLeft[newInternalNode]=mergeVectors(vectUp,bestUpLength,False,probVect[node],bestDownLength,isTip,isUpDown=True)
			bestAppendingLength=estimateBranchLengthWithDerivative(probVectUpLeft[newInternalNode],newPartials,fromTipC=appendedIsTip)
		else:
			bestAppendingLength=estimateBranchLengthWithDerivative(probVectUpLeft[newInternalNode],newPartials,fromTipC=appendedIsTip)
			probVectUpRight[newInternalNode]=mergeVectors(vectUp,bestUpLength,False,newPartials,bestAppendingLength,appendedIsTip,isUpDown=True)
			bestDownLength=estimateBranchLengthWithDerivative(probVectUpRight[newInternalNode],probVect[node],fromTipC=isTip)
		probVect[newInternalNode]=mergeVectors(probVect[node],bestDownLength,isTip,newPartials,bestAppendingLength,appendedIsTip)
		if probVect[newInternalNode]==None:
			print("newInternalNode.probVect is None after updating the optimal branch lengths")
			bestAppendingLength=oneMutBLen/5
			bestDownLength=oneMutBLen/5
			probVect[newInternalNode]=mergeVectors(probVect[node],bestDownLength,isTip,newPartials,bestAppendingLength,appendedIsTip)
	shorten(probVect[newInternalNode])
	probVectUpRight[newInternalNode]=mergeVectors(vectUp,bestUpLength,False,newPartials,bestAppendingLength,appendedIsTip,isUpDown=True)
	if probVectUpRight[newInternalNode]==None:
		bestUpLength=estimateBranchLengthWithDerivative(vectUp,probVect[newInternalNode],fromTipC=False)
		probVectUpLeft[newInternalNode]=mergeVectors(vectUp,bestUpLength,False,probVect[node],bestDownLength,isTip,isUpDown=True)
		bestAppendingLength=estimateBranchLengthWithDerivative(probVectUpLeft[newInternalNode],newPartials,fromTipC=appendedIsTip)
		probVectUpRight[newInternalNode]=mergeVectors(vectUp,bestUpLength,False,newPartials,bestAppendingLength,appendedIsTip,isUpDown=True)
		if probVectUpRight[newInternalNode]==None:
			print("newInternalNode.probVectUpRight is None after updating the optimal branch lengths")
			bestUpLength=oneMutBLen/5
			bestAppendingLength=oneMutBLen/5
			probVectUpRight[newInternalNode]=mergeVectors(vectUp,bestUpLength,False,newPartials,bestAppendingLength,appendedIsTip,isUpDown=True)
		probVect[newInternalNode]=mergeVectors(probVect[node],bestDownLength,isTip,newPartials,bestAppendingLength,appendedIsTip)
	shorten(probVectUpRight[newInternalNode])
	probVectUpLeft[newInternalNode]=mergeVectors(vectUp,bestUpLength,False,probVect[node],bestDownLength,isTip,isUpDown=True)
	if probVectUpLeft[newInternalNode]==None:
		bestUpLength=estimateBranchLengthWithDerivative(vectUp,probVect[newInternalNode],fromTipC=False)
		bestDownLength=estimateBranchLengthWithDerivative(probVectUpRight[newInternalNode],probVect[node],fromTipC=isTip)
		probVectUpLeft[newInternalNode]=mergeVectors(vectUp,bestUpLength,False,probVect[node],bestDownLength,isTip,isUpDown=True)
		if probVectUpLeft[newInternalNode]==None:
			print("newInternalNode.probVectUpRight is None after updating the optimal branch lengths")
			bestUpLength=oneMutBLen/5
			bestDownLength=oneMutBLen/5
			probVectUpLeft[newInternalNode]=mergeVectors(vectUp,bestUpLength,False,probVect[node],bestDownLength,isTip,isUpDown=True)
		probVect[newInternalNode]=mergeVectors(probVect[node],bestDownLength,isTip,newPartials,bestAppendingLength,appendedIsTip)
		probVectUpRight[newInternalNode]=mergeVectors(vectUp,bestUpLength,False,newPartials,bestAppendingLength,appendedIsTip,isUpDown=True)
	shorten(probVectUpLeft[newInternalNode])
	dist[appendedNode]=bestAppendingLength
	dist[newInternalNode]=bestUpLength
	dist[node]=bestDownLength
	#TODO
	if HnZ:
		if dist[node]<=effectivelyNon0BLen:
			nDesc0[newInternalNode]=nDesc0[node]
		else:
			nDesc0[newInternalNode]=1
		if dist[appendedNode]>effectivelyNon0BLen:
			nDesc0[newInternalNode]+=1
		else:
			nDesc0[newInternalNode]+=nDesc0[appendedNode]
		if dist[newInternalNode]<=effectivelyNon0BLen:
			parent0=newInternalNode
			while (dist[parent0]<=effectivelyNon0BLen) and up[parent0]!=None:
				parent0=up[parent0]
				#print("going up "+str(nDesc0[parent0])+" "+str(nDesc0[children[parent0][0]])+" "+str(nDesc0[children[parent0][1]])+" "+str(dist[children[parent0][0]])+" "+str(dist[children[parent0][1]])+" ")
				nDesc0[parent0]+=nDesc0[newInternalNode]-1
				#TODO remove after debugging
				# supposedValue=1
				# if not dist[children[parent0][0]]:
				# 	supposedValue=nDesc0[children[parent0][0]]
				# if dist[children[parent0][1]]:
				# 	supposedValue+=1
				# else:
				# 	supposedValue+=nDesc0[children[parent0][1]]
				# if supposedValue!=nDesc0[parent0]:
				# 	print("node nDesc0 wrongly updated")
				# 	print(dist[node])
				# 	print(dist[appendedNode])
				# 	print(dist[newInternalNode])
				# 	print(nDesc0[node])
				# 	print(nDesc0[appendedNode])
				# 	print(nDesc0[newInternalNode])
				# 	print(nDesc0[parent0])
				# 	print(supposedValue)
				# 	raise Exception("exit")

	if not bestAppendingLength:
		probVectTotUp[appendedNode]=None
	if bestUpLength:
		probVectTotUp[newInternalNode]=mergeVectors(vectUp,bestUpLength/2,False,probVect[newInternalNode],bestUpLength/2,False,isUpDown=True)
		shorten(probVectTotUp[newInternalNode])
	if not bestDownLength:
		probVectTotUp[node]=None
	nodeList=[(node,2),(up[newInternalNode],child),(appendedNode,2)]
	updatePartials(tree,nodeList)
	return None


#remove node from the current position in the tree and re-attach it at a new given place new bestNode.
# First remove node from the tree, then update the genome lists;
# then find the exact best reattachment of node and update the genome lists again using function placeSubtreeOnTree().
#TODO update nDesc0 as a subtree is removed
def cutAndPasteNode(tree,node,bestNode,bestBranchLengths,bestLK,passedProbVect):
	up=tree.up
	children=tree.children
	probVectUpRight=tree.probVectUpRight
	probVectUpLeft=tree.probVectUpLeft
	mutations=tree.mutations
	dist=tree.dist
	probVect=tree.probVect
	minorSequences=tree.minorSequences
	#TODO
	nDesc0=tree.nDesc0
	#TODO - remove, only for debugging
	#if HnZ:
	#	print("Re-calculating the nDesc0 before removing subtree from the tree")
	#	calculateNDesc0(tree,t1,checkExisting=True)
	#remove node from the tree
	parentNode=up[node]
	if node==children[parentNode][0]:
		sibling=children[parentNode][1]
	else:
		sibling=children[parentNode][0]
	if up[parentNode]!=None:
		if parentNode==children[up[parentNode]][0]:
			childP=0
		else:
			childP=1
		children[up[parentNode]][childP]=sibling
		#TODO
		if HnZ:
			if dist[parentNode]<=effectivelyNon0BLen:
				if dist[node]>effectivelyNon0BLen:
					toBeRemoved=-1
				else:
					toBeRemoved=-nDesc0[node]
				parent0=parentNode
				while (dist[parent0]<=effectivelyNon0BLen) and up[parent0]!=None:
					nDesc0[up[parent0]]+=toBeRemoved
					parent0=up[parent0]
					#print("Removing subtree, nDesc0 after: "+str(nDesc0[parent0]))
					if nDesc0[parent0]<=0:
						print("problem removing subtree")
						raise Exception("exit") 
						
	up[sibling]=up[parentNode]
	dist[sibling]=dist[sibling]+dist[parentNode]
	#TODO - remove, only for debugging
	#if HnZ:
	#	print("Re-calculating the nDesc0 after removing subtree from the tree")
	#	calculateNDesc0(tree,t1,checkExisting=True)
	if mutations[parentNode]:
		mutations[sibling]=mergeMutationLists(mutations[parentNode],mutations[sibling])
	#update likelihood lists after node removal
	if up[sibling]==None:
		dist[sibling]=1.0
		if children[sibling]:
			probVect1=probVect[children[sibling][1]]
			if mutations[children[sibling][1]]:
				probVect1=passGenomeListThroughBranch(probVect1,mutations[children[sibling][1]],dirIsUp=True)
			probVectUpRight[sibling]=rootVector(probVect1,dist[children[sibling][1]],(len(children[children[sibling][1]])==0 and len(minorSequences[children[sibling][1]])==0),tree,sibling)
			probVect1=probVect[children[sibling][0]]
			if mutations[children[sibling][0]]:
				probVect1=passGenomeListThroughBranch(probVect1,mutations[children[sibling][0]],dirIsUp=True)
			probVectUpLeft[sibling]=rootVector(probVect1,dist[children[sibling][0]],(len(children[children[sibling][0]])==0 and len(minorSequences[children[sibling][0]])==0),tree,sibling)
			nodeList=[(children[sibling][0],2),(children[sibling][1],2)]
			updatePartials(tree,nodeList)
	else:
		nodeList=[(sibling,2),(up[sibling],childP)]
		updatePartials(tree,nodeList)
	#re-place the node and re-update the vector lists
	newRoot = placeSubtreeOnTree(tree,bestNode,passedProbVect,node,bestLK,bestBranchLengths)

	#TODO - remove, only for debugging
	#if HnZ:
	#	print("Re-calculating the nDesc0 after one SPR move")
	#	calculateNDesc0(tree,t1,checkExisting=True)
	
	topologyChanges[0]+=1
	if (writeTreesToFileEveryTheseSteps>0) and (topologyChanges[0]%writeTreesToFileEveryTheseSteps)==0:
		currentRoot=sibling
		while up[currentRoot]!=None:
			currentRoot=up[currentRoot]
		intermediateTreesFile.write("Topology "+str(topologyChanges[0])+"\n")
		intermediateTreesFile.write(createNewick(tree,currentRoot,binary=binaryTree,namesInTree=namesInTree)+"\n")
		
	if (writeLKsToFileEveryTheseSteps>0) and (topologyChanges[0]%writeLKsToFileEveryTheseSteps)==0:
		currentRoot=sibling
		while up[currentRoot]!=None:
			currentRoot=up[currentRoot]
		totalLK=calculateTreeLikelihood(tree,currentRoot)
		intermediateLKsFile.write("Topology "+str(topologyChanges[0])+", LK: "+str(totalLK)+"\n")

	#if the root of the tree has changed, return the new root
	if up[sibling]==None:
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
#TODO accounting for branch length changes affecting nDesc0
topologyUpdates=[0]
bLenUpdates=[0]
def traverseTreeForTopologyUpdate(tree,node,strictTopologyStopRules=strictTopologyStopRules,allowedFailsTopology=allowedFailsTopology,thresholdLogLKtopology=thresholdLogLKtopology,thresholdTopologyPlacement=thresholdTopologyPlacement):
	up=tree.up
	children=tree.children
	probVectUpRight=tree.probVectUpRight
	probVectUpLeft=tree.probVectUpLeft
	mutations=tree.mutations
	dist=tree.dist
	probVect=tree.probVect
	minorSequences=tree.minorSequences
	#TODO
	nDesc0=tree.nDesc0
	if aBayesPlusOn:
		if networkOutput:
			alternativePlacements=tree.alternativePlacements
		support=tree.support
	#track if the root has changed, so that the new root node can be returned.
	newRoot=None
	# has the branch length been updated so far? If we change the branch length above the current node, even if we don't perform any SPR, we still need to update the genome lists in the tree.
	bLenChanged=False
	totalImprovement=0.0
	# we avoid the root node since it cannot be re-placed with SPR moves
	if up[node]!=None:
		#TODO - remove, only for debugging
		#if HnZ:
		#	print("Re-calculating the nDesc0 at the start of searching SPR move")
		#	calculateNDesc0(tree,t1,checkExisting=True)
		#evaluate current placement
		parentNode=up[node]
		if children[parentNode][0]==node:
			child=0
			vectUp=probVectUpRight[parentNode]
		else:
			child=1
			vectUp=probVectUpLeft[parentNode]
		if mutations[node]:
			vectUp=passGenomeListThroughBranch(vectUp,mutations[node],dirIsUp=False)
		#score of current tree
		bestCurrenBLen=dist[node]
		isTip=(len(children[node])==0) and (len(minorSequences[node])==0)
		originalLK=appendProbNode(vectUp,probVect[node],isTip,bestCurrenBLen)
		bestCurrentLK=originalLK
		if (bestCurrentLK<thresholdTopologyPlacement) or (supportFor0Branches and aBayesPlusOn):
			bestCurrenBLen=estimateBranchLengthWithDerivative(vectUp,probVect[node],fromTipC=isTip)
			if bestCurrenBLen or dist[node]:
				if (not bestCurrenBLen) or (not dist[node]) or dist[node]/bestCurrenBLen>1.01 or dist[node]/bestCurrenBLen<0.99:
					bLenChanged=True
				bestCurrentLK=appendProbNode(vectUp,probVect[node],isTip,bestCurrenBLen)
				if bestCurrentLK<originalLK:
					bestCurrenBLen=dist[node]
					bestCurrentLK=originalLK
					bLenChanged=False
				if bestCurrentLK==float("-inf"):
					print("Found an infinite cost of bestCurrentLK "+str(bestCurrentLK)+" using appendProbNode()")
					raise Exception("exit")
				
		#TODO correct initial placement likelihood by the change in HnZ
		if HnZ:
			parentNode0=up[node]
			while (dist[parentNode0]<=effectivelyNon0BLen) and up[parentNode0]!=None:
				parentNode0=up[parentNode0]
			if dist[node]>effectivelyNon0BLen:
				modifier0= getHnZ(nDesc0[parentNode0]) - getHnZ(nDesc0[parentNode0]-1)
			else:
				modifier0= getHnZ(nDesc0[parentNode0]) - ( getHnZ(nDesc0[parentNode0]-nDesc0[node]) + getHnZ(nDesc0[node]) )
			#print("Initial details before SPR search: number of replacements "+str(tree.replacements[node])+"   nDesc0[node] "+str(nDesc0[node])+"   nDesc0[parentNode0] "+str(nDesc0[parentNode0])+"   dist[node] "+str(dist[node])+"   modifier0 "+str(modifier0)+"   bestCurrentLK "+str(bestCurrentLK)+"   originalLK "+str(originalLK)+"   getHnZ(nDesc0[node]) "+str(getHnZ(nDesc0[node]))+"   getHnZ(nDesc0[parentNode0]) "+str(getHnZ(nDesc0[parentNode0])))
			bestCurrentLK+=modifier0
			originalLK+=modifier0


		topologyUpdated=False
		#TODO in case of HnZ, also try to replace 0-cost subtrees, as they might have better modifers somewhere else
		if ((bestCurrentLK<thresholdTopologyPlacement) and (not doNotImproveTopology)) or dist[node] or (supportFor0Branches and aBayesPlusOn) or HnZ:
			#now find the best place on the tree where to re-attach the subtree rooted at "node"
			#but to do that we need to consider new vector probabilities after removing the node that we want to replace
			# this is done using findBestParentTopology().
			#TODO - remove, only for debugging
			#if HnZ:
			#	print("Re-calculating the nDesc0 before searching SPR move")
			#	calculateNDesc0(tree,t1,checkExisting=True)
			bestNodeSoFar , bestLKdiff , bestBranchLengths, listOfBestPlacements, branchSupport, passedProbVect = findBestParentTopology(tree,parentNode,child,bestCurrentLK,bestCurrenBLen,strictTopologyStopRules=strictTopologyStopRules,allowedFailsTopology=allowedFailsTopology,thresholdLogLKtopology=thresholdLogLKtopology)
			#TODO - remove, only for debugging
			#if HnZ:
			#	print("Re-calculating the nDesc0 after searching SPR move")
			#	calculateNDesc0(tree,t1,checkExisting=True)
			if bestLKdiff==float("inf"):
				print("Found an infinite improvement of "+str(bestLKdiff)+" using findBestParentTopology()")
				raise Exception("exit")
			if bestLKdiff<-1e50:
				print("Error: found likelihood cost is very heavy, this might mean that the reference used is not the same used to generate the input diff file")
				raise Exception("exit")
			if (bestLKdiff+thresholdTopologyPlacement>bestCurrentLK) and (not doNotImproveTopology):
				topologyUpdated=True
				topNode=up[node]
				if bestNodeSoFar==topNode:
					topologyUpdated=False
				while (not dist[topNode]) and (up[topNode]!=None):
					topNode=up[topNode]
				if bestNodeSoFar==topNode and (not bestBranchLengths[1]):
					topologyUpdated=False
				parentNode=up[node]
				if node==children[parentNode][0]:
					sibling=children[parentNode][1]
				else:
					sibling=children[parentNode][0]
				if bestNodeSoFar==sibling:
					topologyUpdated=False
				if up[bestNodeSoFar]==sibling and (not bestBranchLengths[0]):
					topologyUpdated=False

				if topologyUpdated:
					#TODO
					topologyUpdates[0]+=1
					totalImprovement=(bestLKdiff-originalLK)
					if originalLK==float("-inf"):
						totalImprovement=(bestLKdiff-bestCurrentLK)
					if totalImprovement==float("inf"):
						print("Found an infinite topology improvement of "+str(totalImprovement)+" from originalLK "+str(originalLK)+", bestCurrentLK "+str(bestCurrentLK)+" and bestLKdiff "+str(bestLKdiff))
						raise Exception("exit")
					newRoot = cutAndPasteNode(tree,node,bestNodeSoFar,bestBranchLengths,bestLKdiff,passedProbVect)
					bLenChanged=False
			if (not topologyUpdated) and aBayesPlusOn:
				if networkOutput:
					alternativePlacements[node]=listOfBestPlacements
				support[node]=branchSupport
					
		if (not topologyUpdated) and bLenChanged:
			#TODO updating nDesc0
			bLenUpdates[0]+=1
			if HnZ:
				if dist[node]>effectivelyNon0BLen and (bestCurrenBLen<=effectivelyNon0BLen):
					addendum0=(nDesc0[node]-1)
				elif (dist[node]<=effectivelyNon0BLen) and bestCurrenBLen>effectivelyNon0BLen:
					addendum0=(1-nDesc0[node])
				else:
					addendum0=0
				if addendum0:
					nDesc0[up[node]]+=addendum0
					parent0=up[node]
					while up[parent0]!=None and (dist[parent0]<=effectivelyNon0BLen):
						parent0=up[parent0]
						nDesc0[parent0]+=addendum0

			dist[node]=bestCurrenBLen
			nodeList=[(node,2),(up[node],child)]
			updatePartials(tree,nodeList)
			totalImprovement=(bestCurrentLK-originalLK)
			if originalLK==float("-inf"):
				totalImprovement=0
			if totalImprovement==float("inf"):
				print("Found an infinite branch length improvement of "+str(totalImprovement)+" from originalLK "+str(originalLK)+" and bestCurrentLK "+str(bestCurrentLK))
				raise Exception("exit")

	return newRoot,totalImprovement


#Function to apply the SPR moves identified by the parallelized SPR search (startTopologyUpdatesParallel).
#Applies changes starting from the last ones in the input list (those with highest LK improvement).
#Input list contains for each entry 1)pruned node 2) regraft node 3) LK improvement.
def applySPRMovesParallel(tree,results,strictTopologyStopRules=strictTopologyStopRules,allowedFailsTopology=allowedFailsTopology,thresholdLogLKtopology=thresholdLogLKtopology,thresholdTopologyPlacement=thresholdTopologyPlacement):
	if aBayesPlusOn and networkOutput:
		alternativePlacements=tree.alternativePlacements
	cumulativeImprovement=0.0
	#track if the root has changed, so that the new root node can be returned.
	newRoot=None
	while results:
		node, placement, improvementOld=results.pop()
		if aBayesPlusOn and networkOutput:
			alternativePlacements[node]=[]
		newRoot2,improvement=traverseTreeForTopologyUpdate(tree,node,strictTopologyStopRules=strictTopologyStopRules,allowedFailsTopology=allowedFailsTopology,thresholdLogLKtopology=thresholdLogLKtopology,thresholdTopologyPlacement=thresholdTopologyPlacement)
		cumulativeImprovement+=improvement
		if newRoot2!=None:
			newRoot=newRoot2
	return newRoot, cumulativeImprovement


#traverse the tree (here the input "node" will usually be the root), and for each dirty node ancountered, call traverseTreeForTopologyUpdate() 
# to attempt an SPR move by cutting the subtree rooted at this dirty node and trying to re-append it elsewhere.
#TODO count how many nodes are re-placed, how many are investigated, etc.
def startTopologyUpdates(tree,node,checkEachSPR=False,strictTopologyStopRules=strictTopologyStopRules,allowedFailsTopology=allowedFailsTopology,thresholdLogLKtopology=thresholdLogLKtopology,thresholdTopologyPlacement=thresholdTopologyPlacement,printEvery=10000):
	up=tree.up
	children=tree.children
	dirty=tree.dirty
	replacements=tree.replacements
	if aBayesPlusOn and networkOutput:
		alternativePlacements=tree.alternativePlacements
	nodesToVisit=[node]
	totalImprovement=0.0
	newRoot=None
	numNodes=0
	#TODO remove - just for debugging
	#notDirtyNodes=0
	#numMaxRepNodes=0
	#positiveImprov=0
	#higherThan1=0
	#maxImprov=0.0
	topologyUpdates[0]=0
	bLenUpdates[0]=0
	while nodesToVisit:
		newNode=nodesToVisit.pop()
		for c in children[newNode]:
			nodesToVisit.append(c)
		if dirty[newNode] and replacements[newNode]<=maxReplacements:
			dirty[newNode]=False
			if checkEachSPR:
				root=newNode
				while up[root]!=None:
					root=up[root]
				oldTreeLK=calculateTreeLikelihood(tree,root,checkCorrectness=True)
				reCalculateAllGenomeLists(tree,root, checkExistingAreCorrect=True)
			if aBayesPlusOn and networkOutput:
				alternativePlacements[newNode]=[]
			newRoot2,improvement=traverseTreeForTopologyUpdate(tree,newNode,strictTopologyStopRules=strictTopologyStopRules,allowedFailsTopology=allowedFailsTopology,thresholdLogLKtopology=thresholdLogLKtopology,thresholdTopologyPlacement=thresholdTopologyPlacement)
			# if improvement:
			# 	positiveImprov+=1
			# 	if improvement>1.0:
			# 		higherThan1+=1
			# 		if improvement>maxImprov:
			# 			maxImprov=improvement
			if checkEachSPR:
				root=newNode
				while up[root]!=None:
					root=up[root]
				newTreeLK=calculateTreeLikelihood(tree,root)
				reCalculateAllGenomeLists(tree,root, checkExistingAreCorrect=True)
				print("In startTopologyUpdates, LK score of improvement "+str(newTreeLK)+" - "+str(oldTreeLK)+" = "+str(newTreeLK-oldTreeLK)+", was supposed to be "+str(improvement))
				if newTreeLK-oldTreeLK < improvement-1.0:
					print("In startTopologyUpdates, LK score of improvement "+str(newTreeLK)+" - "+str(oldTreeLK)+" = "+str(newTreeLK-oldTreeLK)+" is less than what is supposed to be "+str(improvement))
					raise Exception("exit")
			totalImprovement+=improvement
			if newRoot2!=None:
				newRoot=newRoot2
			numNodes+=1
			if (numNodes%printEvery)==0:
				print("Processed topology for "+str(numNodes)+" nodes.")
		#elif not dirty[newNode]:
		#	notDirtyNodes+=1
		#elif replacements[newNode]>maxReplacements:
		#	numMaxRepNodes+=1
	#print("Attempted re-placement for "+str(numNodes)+" nodes. "+str(notDirtyNodes)+" nodes not dirty. "+str(numMaxRepNodes)+" nodes reached max re-placements.")
	#print("NUmber of positive improvements to the LK: "+str(positiveImprov)+" , o which >1: "+str(higherThan1)+" , and with maximum "+str(maxImprov))
	print("Topology updates "+str(topologyUpdates[0])+" ; bLen updates "+str(bLenUpdates[0]))
	return newRoot,totalImprovement


#parallelized version of startTopologyUpdates(): takes a number "corNum" in input, and only performs search on the nodes that are assinged to that core.
#traverse the tree (here the input "node" will usually be the root), and for each appropriate dirty node encountered, search the best SPR move re-placing that node.
# attempt an SPR move by cutting the subtree rooted at this dirty node and trying to re-append it elsewhere.
def startTopologyUpdatesParallel(inputTuple):
	tree, startingNode,corNum,strictTopologyStopRules,allowedFailsTopology,thresholdLogLKtopology,thresholdTopologyPlacement= inputTuple
	up=tree.up
	children=tree.children
	dirty=tree.dirty
	replacements=tree.replacements
	coreNum=tree.coreNum
	probVectUpRight=tree.probVectUpRight
	probVectUpLeft=tree.probVectUpLeft
	minorSequences=tree.minorSequences
	probVect=tree.probVect
	mutations=tree.mutations
	dist=tree.dist
	nodesToVisit=[startingNode]
	proposedMoves=[]
	nodesSearched=0
	print("Starting SPR search within core "+str(corNum))
	while nodesToVisit:
		node=nodesToVisit.pop()
		for c in children[node]:
			nodesToVisit.append(c)
		if dirty[node] and replacements[node]<=maxReplacements and coreNum[node]==corNum:
			#placement,improvement=traverseTreeForTopologyUpdateParallel(newNode,strictTopologyStopRules=strictTopologyStopRules,allowedFailsTopology=allowedFailsTopology,thresholdLogLKtopology=thresholdLogLKtopology,thresholdTopologyPlacement=thresholdTopologyPlacement)
			placement=None
			improvement=0
			# we avoid the root node since it cannot be re-placed with SPR moves
			if up[node]!=None:
				nodesSearched+=1
				#evaluate current placement
				parentNode=up[node]
				if children[parentNode][0]==node:
					child=0
					vectUp=probVectUpRight[parentNode]
				else:
					child=1
					vectUp=probVectUpLeft[parentNode]
				if mutations[node]:
					vectUp=passGenomeListThroughBranch(vectUp,mutations[node],dirIsUp=False)
				#score of current tree
				bestCurrenBLen=dist[node]
				isTip=(len(children[node])==0) and (len(minorSequences[node])==0)
				bestCurrentLK=appendProbNode(vectUp,probVect[node],isTip,bestCurrenBLen)

				topologyUpdated=False
				#TODO in case of HnZ, also try to replace 0-cost subtrees, as they might have better modifers somewhere else
				if (bestCurrentLK<thresholdTopologyPlacement) or HnZ:# or (supportFor0Branches and aBayesPlusOn):
					#now find the best place on the tree where to re-attach the subtree rooted at "node"
					#but to do that we need to consider new vector probabilities after removing the node that we want to replace
					# this is done using findBestParentTopology().
					bestNodeSoFar , bestLKdiff , bestBranchLengths, listOfBestPlacements, branchSupport, passedProbVect = findBestParentTopology(tree,parentNode,child,bestCurrentLK,bestCurrenBLen,strictTopologyStopRules=strictTopologyStopRules,allowedFailsTopology=allowedFailsTopology,thresholdLogLKtopology=thresholdLogLKtopology,sequentialSearch=False)
					if (bestLKdiff+thresholdTopologyPlacement>bestCurrentLK) and (not doNotImproveTopology):
						topologyUpdated=True
						topNode=up[node]
						if bestNodeSoFar==topNode:
							topologyUpdated=False
						while (not dist[topNode]) and (up[topNode]!=None):
							topNode=up[topNode]
						if bestNodeSoFar==topNode and (not bestBranchLengths[1]):
							topologyUpdated=False
						parentNode=up[node]
						if node==children[parentNode][0]:
							sibling=children[parentNode][1]
						else:
							sibling=children[parentNode][0]
						if bestNodeSoFar==sibling:
							topologyUpdated=False
						if up[bestNodeSoFar]==sibling and (not bestBranchLengths[0]):
							topologyUpdated=False

						if topologyUpdated:
							improvement=(bestLKdiff-bestCurrentLK)
							placement=bestNodeSoFar
					#TODO for now SPRTA not allowed with parallelization 
					# if aBayesPlusOn:
					# 	if networkOutput:
					# 		alternativePlacements[node]=listOfBestPlacements
					# 	support[node]=branchSupport
				
				if placement!=None and (not doNotImproveTopology):
					proposedMoves.append((node,placement,improvement))
	print("Searched "+str(nodesSearched)+" nodes within core "+str(corNum)+" and found "+str(len(proposedMoves))+" proposed SPR moves")
	return proposedMoves


#Given a tree, and a final substitution rate matrix, calculate the likelihood of the tree
#TODO include also for HnZ modifier if required
def calculateTreeLikelihood(tree,root,checkCorrectness=False):
	up=tree.up
	children=tree.children
	minorSequences=tree.minorSequences
	probVect=tree.probVect
	mutations=tree.mutations
	dist=tree.dist
	#TODO
	nDesc0=tree.nDesc0
	node=root
	#direction 0 means from parent, direction 1 means from a child
	lastNode=None
	direction=0
	totalLK=0.0
	while node!=None:
		if direction==0:
			if children[node]:
				node=children[node][0]
			else:
				lastNode=node
				node=up[node]
				direction=1
				#we can ignore the likelihood contribution from the normalization of the likelihoods of the ambiguity characters at terminal nodes:
				# these are anyway the same for all trees and substitution models!
		else :
			if lastNode==children[node][0]:
				node=children[node][1]
				direction=0
			else:
				probVect0=probVect[children[node][0]]
				if mutations[children[node][0]]:
					probVect0=passGenomeListThroughBranch(probVect0,mutations[children[node][0]],dirIsUp=True)
				probVect1=probVect[children[node][1]]
				if mutations[children[node][1]]:
					probVect1=passGenomeListThroughBranch(probVect1,mutations[children[node][1]],dirIsUp=True)
				newLower, Lkcontribution=mergeVectors(probVect0,dist[children[node][0]],(len(children[children[node][0]])==0 and len(minorSequences[children[node][0]])==0),probVect1,dist[children[node][1]],(len(children[children[node][1]])==0 and len(minorSequences[children[node][1]])==0),returnLK=True,numMinor1=len(minorSequences[children[node][0]]),numMinor2=len(minorSequences[children[node][1]]))
				totalLK+=Lkcontribution
				#TODO adding HnZ contribution
				if HnZ and dist[node]>effectivelyNon0BLen:
					totalLK+=getHnZ(nDesc0[node])
				if newLower==None:
					print("Strange, inconsistent lower genome list creation in calculateTreeLikelihood(); old list, and children lists")
					raise Exception("exit")
				elif checkCorrectness and areVectorsDifferentDebugging(probVect[node],newLower):
					print("Strange, while calculating tree likelihood encountered non-updated lower likelihood genome list at node "+str(node)+" with children ")
					raise Exception("exit")
				lastNode=node
				node=up[node]
				direction=1
	#now add contribution from the root
	probVectRootUp=probVect[root]
	if mutations[root]:
		probVectRootUp=passGenomeListThroughBranch(probVectRootUp,mutations[root],dirIsUp=True)
	totalLK+=findProbRoot(probVectRootUp)
	return totalLK


# find sequence positions with high probability of being errors, and write them to file.
def calculateErrorProbabilities(tree,root,file,minErrorProb,namesInTree):
	up=tree.up
	children=tree.children
	probVectUpRight=tree.probVectUpRight
	probVectUpLeft=tree.probVectUpLeft
	minorSequences=tree.minorSequences
	probVect=tree.probVect
	name=tree.name
	mutations=tree.mutations
	dist=tree.dist
	if not useRateVariation:
		mutMatrix=mutMatrixGlobal
	if not (usingErrorRate and errorRateSiteSpecific):
		errorRate=errorRateGlobal
	node=root
	#direction 0 means from parent, direction 1 means from a child
	lastNode=None
	direction=0
	while node!=None:
		if direction==0:
			if len(children[node])==0:
				file.write(">"+namesInTree[name[node]]+"\n")
				if len(minorSequences[node])>0:
					for idNode in minorSequences[node]:
						file.write(">"+namesInTree[idNode]+"\n")
				else:
					if node==children[up[node]][0]:
						probVectP=probVectUpRight[up[node]]
					else:
						probVectP=probVectUpLeft[up[node]]
					if mutations[node]:
						probVectP=passGenomeListThroughBranch(probVectP,mutations[node],dirIsUp=False)
					probVectC=probVect[node]
					indexEntry1, indexEntry2, pos = 0, 0, 0
					entry1=probVectP[indexEntry1]
					entry2=probVectC[indexEntry2]

					while True:
						if entry2[0]==5:
							if entry1[0]==4 or entry1[0]==5:
								pos=min(entry1[1],entry2[1])
							else:
								pos+=1
						elif entry1[0]==5: # case entry1 is N
							if entry2[0]==4:
								pos=min(entry1[1],entry2[1])
							else:
								pos+=1
						else:
							totLen1=dist[node]
							if entry1[0]<5:
								if len(entry1)==3+usingErrorRate:
									totLen1+=entry1[2]
								elif len(entry1)==4+usingErrorRate:
									totLen1+=entry1[3]
							else:
								if len(entry1)>3:
									totLen1+=entry1[2]

							if entry1[0]==4: # case entry1 is R
								if entry2[0]==4:
									pos=min(entry1[1],entry2[1])

								elif entry2[0]==6:
									i1=entry2[1]
									if entry2[-1][i1]<0.1:
										if useRateVariation:
											mutMatrix=mutMatrices[pos]
										if usingErrorRate and errorRateSiteSpecific: errorRate = errorRates[pos]
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
											file.write(str(pos+1)+"\t"+"X"+"\t"+str(errProb)+"\n")
									pos+=1
								
								else: #entry1 is reference and entry2 is a different but single nucleotide
									i1=entry2[1]
									i2=entry2[0]
									if useRateVariation:
										mutMatrix=mutMatrices[pos]
									if usingErrorRate and errorRateSiteSpecific: errorRate = errorRates[pos]
									if len(entry1)<4+usingErrorRate:
										errorProb=errorRate*0.33333
										mutProb=mutMatrix[i1][i2]*totLen1
										errorProb=errorProb/(errorProb+mutProb)
									else:
										mutprob1=rootFreqs[i1]*mutMatrix[i1][i2]*totLen1
										mutprob2=rootFreqs[i2]*mutMatrix[i2][i1]*entry1[2]
										errorProb=rootFreqs[i1]*errorRate*0.33333
										normalization=mutprob1+mutprob2+errorProb
										errorProb=errorProb/normalization
									if errorProb>=minErrorProb:
										file.write(str(pos+1)+"\t"+allelesList[i2]+"\t"+str(errorProb)+"\n")
									pos+=1

							# entry1 is of type "O"
							elif entry1[0]==6:
								normalization=0.0
								if useRateVariation:
									mutMatrix=mutMatrices[pos]
								if entry2[0]==6:
									if usingErrorRate and errorRateSiteSpecific: errorRate = errorRates[pos]
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
										file.write(str(pos+1)+"\t"+"X"+"\t"+str(errorProb)+"\n")
								else:
									if entry2[0]==4:
										i2=entry1[1]
									else:
										i2=entry2[0]
									if usingErrorRate and errorRateSiteSpecific: errorRate = errorRates[pos]
									errorProb=(1.0-entry1[-1][i2])*errorRate*0.33333
									noMutProb=entry1[-1][i2]
									mutProb=0.0
									for i in range4:
										if i!=i2:
											mutProb+=entry1[-1][i]*mutMatrix[i][i2]*totLen1
									normalization=errorProb+noMutProb+mutProb
									errorProb=errorProb/normalization
									if errorProb>=minErrorProb:
										file.write(str(pos+1)+"\t"+allelesList[i2]+"\t"+str(errorProb)+"\n")
								pos+=1

							else: #entry1 is a non-ref nuc
								i1=entry1[0]
								if entry2[0]!=i1: #entry1 and entry2 are of different types
									if useRateVariation:
										mutMatrix=mutMatrices[pos]
									if entry2[0]==6:
										if entry2[-1][i1]<0.1:#in case the reference allele is possible, we ignore the other possibilities since they are unlikely
											if usingErrorRate and errorRateSiteSpecific: errorRate = errorRates[pos]
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
												file.write(str(pos+1)+"\t"+"X"+"\t"+str(errProb)+"\n")

									#entry2 is a nucleotide type (like entry1)		
									else:
										if entry2[0]==4:
											i2=entry1[1]
										else:
											i2=entry2[0]
										if usingErrorRate and errorRateSiteSpecific: errorRate = errorRates[pos]
										if len(entry1)<4+usingErrorRate:
											errorProb=errorRate*0.33333
											mutProb=mutMatrix[i1][i2]*totLen1
											errorProb=errorProb/(errorProb+mutProb)
										else:
											mutprob1=rootFreqs[i1]*mutMatrix[i1][i2]*totLen1
											mutprob2=rootFreqs[i2]*mutMatrix[i2][i1]*entry1[2]
											errorProb=rootFreqs[i1]*errorRate*0.33333
											normalization=mutprob1+mutprob2+errorProb
											errorProb=errorProb/normalization
										if errorProb>=minErrorProb:
											file.write(str(pos+1)+"\t"+allelesList[i2]+"\t"+str(errorProb)+"\n")
											
								pos+=1

						if pos==lRef:
							break

						if entry1[0]<4 or entry1[0]==6:
							indexEntry1+=1
							entry1=probVectP[indexEntry1]
						elif pos==entry1[1]:
							indexEntry1+=1
							entry1=probVectP[indexEntry1]
						if entry2[0]<4 or entry2[0]==6:
							indexEntry2+=1
							entry2=probVectC[indexEntry2]
						elif pos==entry2[1]:
							indexEntry2+=1
							entry2=probVectC[indexEntry2]

			if children[node]:
				node=children[node][0]
			else:
				lastNode=node
				node=up[node]
				direction=1
			
		else :
			if lastNode==children[node][0]:
				node=children[node][1]
				direction=0
			else:
				lastNode=node
				node=up[node]
				direction=1


#moving a mutation list along a branch of a MAT; the first (mutations1) is assumed to be an entry.mutationlist and the second (mutations2) would usually be a node.mutations from the MAT.
# this function is also used to extract the sub-list of mutations within a genome position range (between minPos and maxPos, extremes included).
# when mutations2!=[], mutations in the returned list will be with respect to the reference, not the current MAT location, as in the case of a "N" entry.
#Differently from passGenomeListThroughBranch(), this function uses a different type of mutation list (differences with respect to the reference).
def passMutationListThroughBranch(mutations1,mutations2,dirIsUp=False):
	finalMutations=[]
	ind1, ind2, pos1=0, 0, 0
	while True:
		if ind1<len(mutations1):
			pos1=mutations1[ind1][0]
			if ind2<len(mutations2):
				pos2=mutations2[ind2][0]
				if pos1<pos2:
					finalMutations.append(mutations1[ind1])
					ind1+=1
				else:
					if dirIsUp:
						endNuc=mutations2[ind2][1]
					else:
						endNuc=mutations2[ind2][2]
					if endNuc!=refIndeces[pos2-1]:
						finalMutations.append((pos2,endNuc))

					ind2+=1
					if pos1==pos2:
						ind1+=1
			else:
				finalMutations.append(mutations1[ind1])
				ind1+=1
		else:
			if ind2<len(mutations2):
				pos2=mutations2[ind2][0]
				if dirIsUp:
					endNuc=mutations2[ind2][1]
				else:
					endNuc=mutations2[ind2][2]
				if endNuc!=refIndeces[pos2-1]:
					finalMutations.append((pos2,endNuc))
				ind2+=1
			else:
				break

	return finalMutations


#Given a tree and its genome lists, calculate mutation counts and waiting times for expectation maximization estimation of substitution rates
def expectationMaximizationCalculationRates(tree,root,trackMutations=False):
	up=tree.up
	children=tree.children
	probVectUpRight=tree.probVectUpRight
	probVectUpLeft=tree.probVectUpLeft
	minorSequences=tree.minorSequences
	probVect=tree.probVect
	mutations=tree.mutations
	dist=tree.dist
	minErrorPrb=0.0000000001
	if trackMutations:
		Ns=[]
		mutationsInf=[]
		for i in range(len(up)):
			Ns.append([])
			mutationsInf.append([])
		tree.mutationsInf=mutationsInf
		tree.Ns=Ns
		#Ns=[[]] * len(up)
		#mutationsInf=[[]] * len(up)
		if usingErrorRate:
			errors=[]
			for i in range(len(up)):
				errors.append([])
			#errors=[[]] * len(up)
			tree.errors=errors
	if not useRateVariation:
		mutMatrix=mutMatrixGlobal
	if usingErrorRate and (not errorRateSiteSpecific):
		errorRate=errorRateGlobal
	node=root
	mutationsList=[]
	for mut in mutations[root]:
		mutationsList.append((mut[0],mut[2]))
	#direction 0 means from parent, direction 1 means from a child
	lastNode=None
	direction=0
	counts=[[0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0]]
	waitingTimes=[0.0,0.0,0.0,0.0]
	numTips=0
	if usingErrorRate:
		errorCount=0.0
		observedTotNucsSites=0
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
			if len(children[node])==0:
				nodeIsLeaf=True
				numTips+=1+len(minorSequences[node])
			else:
				nodeIsLeaf=False

			if (dist[node] or (usingErrorRate and nodeIsLeaf)) and up[node]!=None:
				if useRateVariation:
					totTreeLength+=dist[node]
				#update counts and waiting times
				if node==children[up[node]][0]:
					probVectP=probVectUpRight[up[node]]
				else:
					probVectP=probVectUpLeft[up[node]]
				if mutations[node]:
					probVectP=passGenomeListThroughBranch(probVectP,mutations[node],dirIsUp=False)
				probVectC=probVect[node]
				indexEntry1, indexEntry2, pos, indexMutationsList = 0, 0, 0, 0
				entry1=probVectP[indexEntry1]
				entry2=probVectC[indexEntry2]
				if (entry1[0]==4 or entry1[0]==5) and (entry2[0]==4 or entry2[0]==5):
					end=min(entry1[1],entry2[1])
				else:
					end=1

				while True:
					while (indexMutationsList<len(mutationsList)) and (mutationsList[indexMutationsList][0]<pos):
						indexMutationsList+=1
					if entry2[0]==5: # case N
						if (entry1[0]==4 or entry1[0]==5):
							end=min(entry1[1],entry2[1])
						else:
							end=pos+1
						if usingErrorRate and nodeIsLeaf:
							if errorRateSiteSpecific:
								observedNucsSites[pos]-=(1+len(minorSequences[node]))
							else:
								observedTotNucsSites-=(end-pos)*(1+len(minorSequences[node]))
						if useRateVariation:
							trackingNs[pos]-=dist[node]
						if trackMutations:
							if (not Ns[node]) or (type(Ns[node][-1])==int or Ns[node][-1][1]!=entry2[1]):
								Ns[node].append((pos+1,entry2[1]))
						pos=end
						if useRateVariation:
							trackingNs[pos]+=dist[node]
						if usingErrorRate and errorRateSiteSpecific and nodeIsLeaf:
							observedNucsSites[pos]+=(1+len(minorSequences[node]))

					elif entry1[0]==5: # case N
						if entry2[0]==4:
							end=min(entry1[1],entry2[1])
						else:
							end=pos+1
						if useRateVariation:
							trackingNs[pos]-=dist[node]
						#if usingErrorRate and nodeIsLeaf: #we treat leaf nodes as if they were not observed, since there is no information regarding errors in this case - unless there are minor sequences.
							#if not minorSequences[node]:
							#	if errorRateSiteSpecific:
							#		observedNucsSites[pos]-=1
							#	else:
							#		observedTotNucsSites-=(end-pos)
						pos=end
						#if usingErrorRate and errorRateSiteSpecific and nodeIsLeaf and (not minorSequences[node]):
						#	observedNucsSites[pos]+=1
						if useRateVariation:
							trackingNs[pos]+=dist[node]
					else:
						totLen1=dist[node]
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

						if entry1[0]==4 and entry2[0]==4: # case entry1 is R
							end=min(entry1[1],entry2[1])
							if (not totLen2) and dist[node]:
								for i in range4:
									waitingTimes[i]+=totLen1*(cumulativeBases[end][i]-cumulativeBases[pos][i])

								while (indexMutationsList<len(mutationsList)) and (mutationsList[indexMutationsList][0]<end):
									altNuc=mutationsList[indexMutationsList][1]
									altPos=mutationsList[indexMutationsList][0]
									refNuc=refIndeces[altPos]
									waitingTimes[refNuc]-=totLen1
									waitingTimes[altNuc]+=totLen1
									indexMutationsList+=1
									if useRateVariation:
										waitingTimesSites[altPos][altNuc]+=totLen1
										waitingTimesSites[altPos][refNuc]-=totLen1
							pos=end
						else:
							
							# entry1 is of type "O"
							if entry1[0]==6:
								if not totLen2:
									normalization=0.0
									if useRateVariation:
										mutMatrix=mutMatrices[pos]
										waitingTimesSites[pos][refIndeces[pos]]-=totLen1

									if entry2[0]==6:
										if trackMutations and nodeIsLeaf:
											Ns[node].append(pos+1)
										if nodeIsLeaf and usingErrorRate :
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
																mutationsInf[node].append((i,pos+1,j,prob))
														if useRateVariation:
															waitingTimesSites[pos][i]+=(totLen1/2)*prob
															waitingTimesSites[pos][j]+=(totLen1/2)*prob
															countsSites[pos]+=prob
															if prob<0.0:
																raise Exception("exit")

									else: #entry1 is O and entry2 is a nucleotide
										if entry2[0]==4:
											i2=entry1[1]
										else:
											i2=entry2[0]
										if nodeIsLeaf and usingErrorRate and (not minorSequences[node]):
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
													errors[node].append((4,pos+1,i2,errorProb))
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
															mutationsInf[node].append((i,pos+1,i2,prob))
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
															mutationsInf[node].append((i,pos+1,i2,prob))
													if useRateVariation:
														waitingTimesSites[pos][i]+=(totLen1/2)*prob
														waitingTimesSites[pos][i2]+=(totLen1/2)*prob
														countsSites[pos]+=prob
														if prob<0.0:
															print("Negative value, numerical issue")
															raise Exception("exit")

							else: #entry1 is a nucleotide
								if entry1[0]==4:
									i1=entry2[1]
								else:
									i1=entry1[0]

								if entry2[0]==6: #entry1 is a nucleotide and entry2 is O
									if trackMutations and nodeIsLeaf: # and (not minorSequences[node]):
										Ns[node].append(pos+1)
									if entry2[-1][i1]>0.1:#in case the reference allele is possible, we ignore the other possibilities since they are unlikely
										waitingTimes[i1]+=totLen1
										if useRateVariation:
											waitingTimesSites[pos][refIndeces[pos]]-=totLen1
											waitingTimesSites[pos][i1]+=totLen1
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
											if useRateVariation:
												waitingTimesSites[pos][refIndeces[pos]]-=totLen1
												waitingTimesSites[pos][i1]+=totLen1*errProb
											waitingTimes[i1]+=totLen1*errProb
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
													waitingTimes[i]+=totLen1*(probi+prob1/2)
													waitingTimes[i1]+=totLen1*prob1/2
													counts[i1][i]+=prob1
													if useRateVariation:
														waitingTimesSites[pos][i]+=totLen1*(probi+prob1/2)
														waitingTimesSites[pos][i1]+=totLen1*prob1/2
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
												waitingTimesSites[pos][refIndeces[pos]]-=totLen1
												waitingTimesSites[pos][i1]+=totLen1*errProb
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
											if useRateVariation:
												waitingTimesSites[pos][refIndeces[pos]]-=totLen1
											stayProb1=1.0+mutMatrix[i1][i1]*entry1[2]
											if stayProb1<0:
												approxFailed1=True
												stayProb1=0.25
											else:
												approxFailed1=False

											for i in range4:
												stayProb2=1.0+mutMatrix[i][i]*totLen1
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
														tot3*=totLen1
														tot3+=entry2[-1][i]
													normalization+=prob*tot3
												else:
													if approxFailed1:
														prob=rootFreqs[i]*0.25*stayProb2*entry2[-1][i]
													else:
														prob=rootFreqs[i]*mutMatrix[i][i1]*entry1[2]*stayProb2*entry2[-1][i]
													normalization+=prob

											for i in range4:
												stayProb2=1.0+mutMatrix[i][i]*totLen1
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
															waitingTimes[i]+=totLen1*tot3
															if useRateVariation:
																waitingTimesSites[pos][i]+=totLen1*tot3
														else:
															if approxFailed2:
																tot3=prob*0.25*entry2[-1][j]/normalization
															else:
																tot3=prob*mutMatrix[i][j]*totLen1*entry2[-1][j]/normalization
															waitingTimes[i]+=(totLen1/2)*tot3
															waitingTimes[j]+=(totLen1/2)*tot3
															counts[i][j]+=tot3
															if trackMutations and (not nodeIsLeaf) and (tot3>minMutProb):
																mutationsInf[node].append((i1,pos+1,j,tot3))
															if useRateVariation:
																waitingTimesSites[pos][i]+=(totLen1/2)*tot3
																waitingTimesSites[pos][j]+=(totLen1/2)*tot3
																countsSites[pos]+=tot3
																if tot3<0.0:
																	raise Exception("exit")
																
												else:
													if approxFailed1:
														prob=rootFreqs[i]*0.25*stayProb2*entry2[-1][i]/normalization
													else:
														prob=rootFreqs[i]*mutMatrix[i][i1]*entry1[2]*stayProb2*entry2[-1][i]/normalization
													waitingTimes[i]+=totLen1*prob
													if useRateVariation:
														waitingTimesSites[pos][i]+=totLen1*prob

										else:
											if useRateVariation:
												waitingTimesSites[pos][refIndeces[pos]]-=totLen1
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
														mutationsInf[node].append((i1,pos+1,i,prob))
													if useRateVariation:
														waitingTimesSites[pos][i1]+=(totLen1/2)*prob
														waitingTimesSites[pos][i]+=(totLen1/2)*prob
														countsSites[pos]+=prob
														if prob<0.0:
															raise Exception("exit")

								else: #entry1 and entry2 are nucleotides
									if entry2[0]<4:
										i2=entry2[0]
									else:
										i2=entry1[1]
									if useRateVariation:
										mutMatrix=mutMatrices[pos]

									if i2==i1:
										if not totLen2:
											waitingTimes[i1]+=totLen1
											if useRateVariation:
												waitingTimesSites[pos][i1]+=totLen1
												waitingTimesSites[pos][refIndeces[pos]]-=totLen1

									else: #entry1 and entry2 are of different nucleotides
										
										if nodeIsLeaf and usingErrorRate and (not minorSequences[node]):
											if errorRateSiteSpecific: errorRate = errorRates[pos]
											if len(entry1)<4+usingErrorRate:
												errorProb=errorRate*0.33333
												mutProb=mutMatrix[i1][i2]*totLen1
												normalization=errorProb+mutProb
												errorProb=errorProb/normalization
												mutProb=mutProb/normalization
												if useRateVariation:
													waitingTimesSites[pos][refIndeces[pos]]-=totLen1
													waitingTimesSites[pos][i1]+=totLen1*(mutProb/2)
													waitingTimesSites[pos][i2]+=totLen1*(errorProb+mutProb/2)
													countsSites[pos]+=mutProb
													if mutProb<0.0:
														raise Exception("exit")
												waitingTimes[i1]+=totLen1*(errorProb+mutProb/2)
												waitingTimes[i2]+=(totLen1*mutProb/2)
												counts[i1][i2]+=mutProb
												if trackMutations:
													if (mutProb>minMutProb):
														mutationsInf[node].append((i1,pos+1,i2,mutProb))
													if (errorProb>minMutProb):
														errors[node].append((i1,pos+1,i2,errorProb))
												errorCount+=errorProb
												if errorRateSiteSpecific:
													errorCountSites[pos]+=errorProb
									
											else:
												mutprob1=rootFreqs[i1]*mutMatrix[i1][i2]*totLen1
												mutprob2=rootFreqs[i2]*mutMatrix[i2][i1]*entry1[2]
												errorProb=rootFreqs[i1]*errorRate*0.33333
												normalization=mutprob1+mutprob2+errorProb
												mutprob1=mutprob1/normalization
												mutprob2=mutprob2/normalization
												errorProb=errorProb/normalization
												waitingTimes[i1]+=totLen1*(mutprob1/2+errorProb)
												waitingTimes[i2]+=totLen1*(mutprob2+mutprob1/2)
												counts[i1][i2]+=mutprob1
												if trackMutations:
													if (mutprob1>minMutProb):
														mutationsInf[node].append((i1,pos+1,i2,mutprob1))
													if (errorProb>minMutProb):
														errors[node].append((i1,pos+1,i2,errorProb))
												errorCount+=errorProb
												if errorRateSiteSpecific:
													errorCountSites[pos]+=errorProb
												if useRateVariation:
													waitingTimesSites[pos][refIndeces[pos]]-=totLen1
													waitingTimesSites[pos][i1]+=(totLen1)*(mutprob1/2+errorProb)
													waitingTimesSites[pos][i2]+=(totLen1)*(mutprob2+mutprob1/2)
													countsSites[pos]+=mutprob1
													if mutprob1<0.0:
														raise Exception("exit")

										elif not totLen2:
											if len(entry1)<4+usingErrorRate:
												if useRateVariation:
													waitingTimesSites[pos][refIndeces[pos]]-=totLen1/2
													waitingTimesSites[pos][i1]+=totLen1/2
													waitingTimesSites[pos][i2]+=totLen1/2
													countsSites[pos]+=1
												waitingTimes[i1]+=(totLen1/2)
												waitingTimes[i2]+=(totLen1/2)
												counts[i1][i2]+=1
												if trackMutations:
													mutationsInf[node].append((i1,pos+1,i2,1.0))

											else:
												noMutProb1=1.0+mutMatrix[i1][i1]*entry1[2]
												if noMutProb1<0:
													noMutProb1=0.25
												noMutProb2=1.0+mutMatrix[i2][i2]*totLen1
												if noMutProb2<0:
													noMutProb2=0.25
												prob1=rootFreqs[i1]*mutMatrix[i1][i2]*totLen1*noMutProb1
												prob2=rootFreqs[i2]*mutMatrix[i2][i1]*entry1[2]*noMutProb2
												normalization=prob1+prob2
												prob1=prob1/normalization
												prob2=prob2/normalization
												waitingTimes[i1]+=(totLen1/2)*prob1
												waitingTimes[i2]+=(totLen1/2)*prob1
												counts[i1][i2]+=prob1
												if trackMutations:
													if (prob1>minMutProb):
														mutationsInf[node].append((i1,pos+1,i2,prob1))
												waitingTimes[i2]+=totLen1*prob2
												if useRateVariation:
													waitingTimesSites[pos][refIndeces[pos]]-=totLen1
													waitingTimesSites[pos][i1]+=(totLen1/2)*prob1
													waitingTimesSites[pos][i2]+=(totLen1/2)*prob1
													waitingTimesSites[pos][i2]+=totLen1*prob2
													countsSites[pos]+=prob1
													if prob1<0.0:
														raise Exception("exit")

							pos+=1

					if pos==lRef:
						break

					if entry1[0]<4 or entry1[0]==6:
						indexEntry1+=1
						entry1=probVectP[indexEntry1]
					elif pos==entry1[1]:
						indexEntry1+=1
						entry1=probVectP[indexEntry1]
					if entry2[0]<4 or entry2[0]==6:
						indexEntry2+=1
						entry2=probVectC[indexEntry2]
					elif pos==entry2[1]:
						indexEntry2+=1
						entry2=probVectC[indexEntry2]

			else: #building Ns vectors for nodes with branch length 0
				if trackMutations:
					pos = 0
					for entry2 in probVect[node]:
						if entry2[0]==5: # case N
							if trackMutations:
								if entry2[1]>(pos+1):
									Ns[node].append((pos+1,entry2[1]))
								else:
									Ns[node].append(pos+1)
							pos=entry2[1]
						else:
							if entry2[0]==4: # case entry2 is R
								pos=entry2[1]
							else:
								if entry2[0]==6: #entry2 is a nucleotide and entry2 is O
									if trackMutations and nodeIsLeaf:
										Ns[node].append(pos+1)
								pos+=1

			if children[node]:
				node=children[node][0]
				if mutations[node]:
					mutationsList = passMutationListThroughBranch(mutationsList,mutations[node],dirIsUp=False)
			else:
				lastNode=node
				if mutations[node]:
					mutationsList = passMutationListThroughBranch(mutationsList,mutations[node],dirIsUp=True)
				node=up[node]
				direction=1
			
		else :
			if lastNode==children[node][0]:
				node=children[node][1]
				if mutations[node]:
					mutationsList = passMutationListThroughBranch(mutationsList,mutations[node],dirIsUp=False)
				direction=0
			else:
				lastNode=node
				if mutations[node]:
					mutationsList = passMutationListThroughBranch(mutationsList,mutations[node],dirIsUp=True)
				node=up[node]
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
					siteErrRates.append(max(minErrorPrb,errorCountSites[i]/observedNuc))
				else:
					siteErrRates.append(minErrorPrb)
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


def variance(numList):
	# calculate mean
	m = sum(numList) / len(numList)
	# calculate variance using a list comprehension
	var_res = sum((xi - m) ** 2 for xi in numList) / len(numList)
	return m, var_res


print("Calculating distances with sample names", flush=True)
#don't put samples already in the tree back to the tree again
distances=distancesFromRefPunishNs(data,samples=data.keys(),samplesInInitialTree=namesInTreeDict,forgetData=True)
if len(distances)>0:
	print("Distances from the reference calculated. Average divergence from refererence: "+str(totDivFromRef[0]/len(distances)), flush=True)
namesInTreeDict.clear()

if rateVariation and (not inputRates):
	useRateVariation=True
	siteRates=[1.0]*lRef
	mutMatrices=[]
	for i in range(lRef):
		mutMatrices.append([])
		for j in range4:
			mutMatrices[i].append(list(mutMatrixGlobal[j]))
if not rateVariation:
	siteRates=None

if inputTree=="": #initialize tree to just the initial root sample
	#extract root genome among those closest to the reference but not empty
	firstSample=distances.pop()
	namesInTree.append(firstSample[1])
	tree1=Tree()
	tree1.addNode()
	tree1.name[-1]=0
	t1=0
	tree1.probVect[-1]=probVectTerminalNode(data[firstSample[1]],tree1,t1)
	data[firstSample[1]]=None
	numSamples=1
else:
	t1=rootIndex1
	numSamples=len(namesInTree)
	if numSamples>minNumSamplesForErrorModel or (not largeUpdate):
		if errorRateSiteSpecificFile or errorRateFixed or estimateErrorRate or estimateSiteSpecificErrorRate:
			usingErrorRate=True
	print(str(len(namesInTree))+" named samples in the initial tree")
tree=tree1

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
	errorRateSiteSpecific=True
elif errorRateFixed:
	errorRateGlobal=errorRateFixed
	updateErrorRates(errorRateGlobal)
elif estimateErrorRate:
	if errorRateInitial:
		errorRateGlobal=errorRateInitial
	else:
		errorRateGlobal=1.0/lRef
	updateErrorRates(errorRateGlobal)
elif estimateSiteSpecificErrorRate:
	if not inputRates:
		if errorRateInitial:
			errorRateGlobal=errorRateInitial
		else:
			errorRateGlobal=1.0/lRef
		errorRates=[errorRateGlobal]*lRef
		updateErrorRates(errorRateGlobal,errorRates=errorRates)
	errorRateSiteSpecific=True

# initial EM round to estimate rate variation etc on the initial tree
if (numSamples>1) and (model!="JC" or ((numSamples>=minNumSamplesForRateVar) and useRateVariation ) or ((numSamples>=minNumSamplesForErrorModel) and usingErrorRate)):
	start=time()
	mutMatrixGlobal, siteRatesEst, errorRateGlobal, errorRatesEst = expectationMaximizationCalculationRates(tree,t1)
	print("EM from initial tree terminated, using rate variation "+str(useRateVariation)+", using error rates "+str(usingErrorRate)+".  ")
	updateMutMatrices(mutMatrixGlobal,siteRates=siteRatesEst)
	if usingErrorRate:
		errorRates=errorRatesEst
		updateErrorRates(errorRateGlobal,errorRates=errorRates)
	reCalculateAllGenomeLists(tree,t1)
	newLk=calculateTreeLikelihood(tree,t1)
	print("LK after first EM: "+str(newLk))
	if usingErrorRate and (estimateErrorRate or estimateSiteSpecificErrorRate):
		oldLK=float("-inf")
		numEMsteps=0
		while (newLk-oldLK>1.0) and numEMsteps<20:
			setAllDirty(tree,t1)
			improvement=traverseTreeToOptimizeBranchLengths(tree,t1)
			reCalculateAllGenomeLists(tree,t1)
			newLkBranch=calculateTreeLikelihood(tree,t1)
			print("Updated "+str(improvement)+" branch lengths leading to LK "+str(newLkBranch))
			mutMatrixGlobal, siteRates, errorRateGlobal, errorRatesEst = expectationMaximizationCalculationRates(tree,t1)
			updateMutMatrices(mutMatrixGlobal,siteRates=siteRates)
			if usingErrorRate:
				errorRates=errorRatesEst
				updateErrorRates(errorRateGlobal,errorRates=errorRates)
			reCalculateAllGenomeLists(tree,t1)
			oldLK=newLk
			newLk=calculateTreeLikelihood(tree,t1)
			print("New LK step "+str(numEMsteps)+": "+str(newLk))
			numEMsteps+=1

	timeRecalculation=time()-start
	print("Time to run initial tree EM estimation: "+str(timeRecalculation))


#Place input samples to create an initial tree (or extend the input tree).
timeFinding=0.0
timePlacing=0.0
lastUpdateNumSamples=numSamples
if not doNotPlaceNewSamples:
	while distances:
		d=distances.pop()
		sample=d[1]
		namesInTree.append(sample)
		newPartials=probVectTerminalNode(data[sample],None,None)
		data[sample]=None
		if (numSamples<minNumSamplesForRateVar or (not useRateVariation)) and ((numSamples%updateSubstMatrixEveryThisSamples)==0):
			if model!="JC":
				if updateSubMatrix(pseudoMutCounts,model,mutMatrixGlobal):
					updateMutMatrices(mutMatrixGlobal,siteRates=siteRates)
		if (numSamples%50000)==0:
			print("Sample num "+str(numSamples), flush=True)
		if (useRateVariation and (numSamples>minNumSamplesForRateVar) and (numSamples>2*lastUpdateNumSamples)):
			start=time()
			lastUpdateNumSamples=numSamples
			if numSamples>minNumSamplesForErrorModel:
				if errorRateSiteSpecificFile or errorRateFixed or estimateErrorRate or estimateSiteSpecificErrorRate:
					usingErrorRate=True
			reCalculateAllGenomeLists(tree,t1)
			mutMatrixGlobal, siteRates, errorRateGlobal, errorRatesEst = expectationMaximizationCalculationRates(tree,t1)
			updateMutMatrices(mutMatrixGlobal,siteRates=siteRates)
			if usingErrorRate:
				errorRates=errorRatesEst
				updateErrorRates(errorRateGlobal,errorRates=errorRates)
			reCalculateAllGenomeLists(tree,t1)
			improvement=traverseTreeToOptimizeBranchLengths(tree,t1)
			reCalculateAllGenomeLists(tree,t1)
			print(" EM to update parameters during initial placement terminated, time taken: "+str(time()-start))

		start=time()
		bestNode , bestScore, bestBranchLengths, bestPassedVect = findBestParentForNewSample(tree,t1,newPartials,numSamples)
		timeFinding+=(time()-start)
		if bestBranchLengths!=None:
			start=time()
			newRoot=placeSampleOnTree(tree,bestNode,bestPassedVect,numSamples,bestScore, bestBranchLengths[0], bestBranchLengths[1], bestBranchLengths[2],pseudoMutCounts)
			if newRoot!=None:
				t1=newRoot
			timePlacing+=(time()-start)
		numSamples+=1

		if (numSamples%saveInitialTreeEvery)==0:
			newickString=createNewick(tree,t1,binary=binaryTree,namesInTree=namesInTree)
			fileName=outputFile+"_initialTree_"+str(numSamples)+"samples.tree"
			file=open(fileName,"w")
			file.write(newickString)
			file.close()
			print("Partial initial tree written to file "+fileName)
print("Sample placement completed", flush=True)
print("Number of minor samples removed from the starting tree: "+str(numMinorsRemoved[0]))
print("Number of placed samples that have become minor sequences: "+str(numMinorsFound[0]))

if HnZ:
	#print("Re-calculating the nDesc0")
	calculateNDesc0(tree,t1) # ,checkExisting=True

if (numChildLKs[0]>0) and (not useFixedThresholdLogLKoptimizationTopology):
	aveChildLK=sumChildLKs[0]/numChildLKs[0]
	thresholdLogLKoptimizationTopology=max(thresholdLogLKoptimizationTopology,-0.2*aveChildLK)
	print("thresholdLogLKoptimizationTopology set to "+str(thresholdLogLKoptimizationTopology))

print("Number of references in the MAT: "+str(numRefs[0]), flush=True)

#After placing all the samples, re-estimate the substitution model and recalculate the likelihoods.
reCalculateAllGenomeLists(tree,t1,countNodes=True)
if errorRateSiteSpecificFile or errorRateFixed or estimateErrorRate or estimateSiteSpecificErrorRate:
	newLk=calculateTreeLikelihood(tree,t1)
	print("Tree LK before error rates EM: "+str(newLk))
	usingErrorRate=True
	mutMatrixGlobal, siteRates, errorRateGlobal, errorRates = expectationMaximizationCalculationRates(tree,t1)
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
	reCalculateAllGenomeLists(tree,t1)
	newLk=calculateTreeLikelihood(tree,t1)
	print("Tree LK after first errors EM: "+str(newLk))
	improvement=traverseTreeToOptimizeBranchLengths(tree,t1)
	reCalculateAllGenomeLists(tree,t1)
	newLk=calculateTreeLikelihood(tree,t1)
	print("Tree LK after branch length optimization: "+str(newLk))


data.clear()

#put sample names in the tree
if debugging and inputTree=="":
	nextLeaves=[t1]
	while nextLeaves:
		node=nextLeaves.pop()
		if not tree.children[node]:
			tree.name[node]="S"+str(tree.name[node])
			print(tree.name[node])
			print(tree.probVect[node])
			for m in range(len(tree.minorSequences[node])):
				tree.minorSequences[node][m]="S"+str(tree.minorSequences[node][m])
		else:
			for c in tree.children[node]:
				nextLeaves.append(c)


#recalculate all genome lists according to the final substitution model
if inputTree=="" or largeUpdate or rateVariation or usingErrorRate:
	start=time()
	reCalculateAllGenomeLists(tree,t1)
	
	#if using error rates, run EM iteratively until convergence in likelihood, otherwise only run one EM step
	if model!="JC" or rateVariation or estimateErrorRate or estimateSiteSpecificErrorRate:
		newLk=calculateTreeLikelihood(tree,t1)
		print("Tree LK before EM: "+str(newLk))
		mutMatrixGlobal, siteRates, errorRateGlobal, errorRates = expectationMaximizationCalculationRates(tree,t1)
		updateMutMatrices(mutMatrixGlobal,siteRates=siteRates)
		if usingErrorRate:
			updateErrorRates(errorRateGlobal,errorRates=errorRates)
		reCalculateAllGenomeLists(tree,t1)
		newLk=calculateTreeLikelihood(tree,t1)
		print("Tree LK after EM: "+str(newLk))
		setAllDirty(tree,t1)
		improvement=traverseTreeToOptimizeBranchLengths(tree,t1)
		reCalculateAllGenomeLists(tree,t1)
		newLk=calculateTreeLikelihood(tree,t1)
		print("Tree LK after branch length optimization: "+str(newLk))
		if estimateErrorRate or estimateSiteSpecificErrorRate:
			oldLK=float("-inf")
			#newLk=calculateTreeLikelihood(tree,t1)
			#print("EM estimation of error rates. Initial LK after first pass: "+str(newLk))
			numEMsteps=0
			while (newLk-oldLK>1.0) and numEMsteps<20:
				setAllDirty(tree,t1)
				improvement=traverseTreeToOptimizeBranchLengths(tree,t1)
				reCalculateAllGenomeLists(tree,t1)
				newLkBranch=calculateTreeLikelihood(tree,t1)
				print("Updated "+str(improvement)+" branch lengths leading to LK "+str(newLkBranch))

				mutMatrixGlobal, siteRates, errorRateGlobal, errorRates = expectationMaximizationCalculationRates(tree,t1)
				updateMutMatrices(mutMatrixGlobal,siteRates=siteRates)
				updateErrorRates(errorRateGlobal,errorRates=errorRates)
				reCalculateAllGenomeLists(tree,t1)
				oldLK=newLk
				newLk=calculateTreeLikelihood(tree,t1)
				print("New LK step "+str(numEMsteps)+": "+str(newLk))
				numEMsteps+=1
			if rateVariation:
				meanRate, varianceRate = variance(siteRates)
				print("Rate variation variance: "+str(varianceRate))
			if errorRateSiteSpecific:
				meanErrorRate, varianceErrorRate = variance(errorRates)
				print("Error rate variation, mean: "+str(meanErrorRate)+" , variance: "+str(varianceErrorRate))
			print("Error rate: "+str(errorRateGlobal))

	timeRecalculation=time()-start
	print("Time to run EM pre-topology estimation: "+str(timeRecalculation))

	print("Number of nodes: "+str(numNodes[0]))
	print("Rs per node: "+str(float(numNodes[2])/numNodes[0]))
	print("Os per node: "+str(float(numNodes[4])/numNodes[0]))
	print("Nucs per node: "+str(float(numNodes[1])/numNodes[0]))
	print("Ns per node: "+str(float(numNodes[3])/numNodes[0]))
	print("MAT mutations per node: "+str(float(numNodes[5])/numNodes[0]))

	start=time()
	newLk=calculateTreeLikelihood(tree,t1)
	print("Now proper branch length optimization, LK before: "+str(newLk))
	setAllDirty(tree,t1)
	improvement=traverseTreeToOptimizeBranchLengths(tree,t1)
	subRound=0
	while subRound<20:
		if (not improvement):
			break
		subRound+=1
		improvement=traverseTreeToOptimizeBranchLengths(tree,t1)
	newLk=calculateTreeLikelihood(tree,t1)
	print("branch length finalization subround "+str(subRound+1)+", number of changes "+str(improvement)+" final LK: "+str(newLk))
	timeForBranchOptimization=(time()-start)
	print("Time for updating branch lengths: "+str(timeForBranchOptimization))

#TODO
if HnZ:
	#print("Re-calculating the nDesc0")
	calculateNDesc0(tree,t1,checkExisting=True)

#find best tree root
if not doNotReroot:
	print("Looking for possible better root", flush=True)
	newT1=findBestRoot(tree,t1,strictTopologyStopRules=strictTopologyStopRules,allowedFailsTopology=allowedFailsTopology,thresholdLogLKtopology=thresholdLogLKtopology,aBayesPlusOn=aBayesPlus)
	if newT1!=t1:
		print("Better root found")
		t1=newT1

		newLk=calculateTreeLikelihood(tree,t1)
		print("Tree LK before EM: "+str(newLk))
		mutMatrixGlobal, siteRates, errorRateGlobal, errorRates = expectationMaximizationCalculationRates(tree,t1)
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
		reCalculateAllGenomeLists(tree,t1)
		newLk=calculateTreeLikelihood(tree,t1)
		print("Tree LK after first errors EM: "+str(newLk))
		improvement=traverseTreeToOptimizeBranchLengths(tree,t1)
		reCalculateAllGenomeLists(tree,t1)
		newLk=calculateTreeLikelihood(tree,t1)
		print("Tree LK after branch length optimization: "+str(newLk))

		print("Looking a second time for possible better root", flush=True)
		newT1=findBestRoot(tree,t1,strictTopologyStopRules=strictTopologyStopRules,allowedFailsTopology=allowedFailsTopology,thresholdLogLKtopology=thresholdLogLKtopology,aBayesPlusOn=aBayesPlus)
		if newT1!=t1:
			print("Better root found again")
			t1=newT1



if (writeTreesToFileEveryTheseSteps>0):
	currentRoot=t1
	intermediateTreesFile.write("Topology 0\n")
	intermediateTreesFile.write(createNewick(tree,currentRoot,binary=binaryTree,namesInTree=namesInTree)+"\n")
if (writeLKsToFileEveryTheseSteps>0):
	currentRoot=t1
	totalLK=calculateTreeLikelihood(tree,currentRoot)
	intermediateLKsFile.write("Topology 0, LK: "+str(totalLK)+"\n")


print(str(len(namesInTree))+" named samples in the tree", flush=True)
internalNodeNamesGiven=False
giveInternalNodeNames(tree,t1,namesInTree=namesInTree,replaceNames=False)
internalNodeNamesGiven=True
print(str(len(namesInTree))+" named nodes in the tree after assigning internal node names", flush=True)

#generate the string corresponding to a line of the tsv file for use in Taxonium.
# TODO should I switch off sample collapsing altogether if calculating support also for 0-dist samples?
# TODO if representative node is 0-dist, then its support is also the support of all represented nodes.
# TODO if not, can we assume that the support is of the represented nodes is 1? Yes.
# TODO pass information about the dist, support and Ns of representative nodes below. Include support of represented nodes only if supportFor0Branches is true, otherwise smpty string.
def tsvForNode(tree,node,name,featureList,namesInTree,identicalTo=""):
	dist=tree.dist
	stringList=[name+"\t"]
	if identicalTo!="":
		stringList.append(identicalTo)
	stringList.append("\t")
	for feat in featureList:
		if node!=None:
			if hasattr(tree, feat):
				feature=getattr(tree, feat)
				if (feat=="support" or feat=="IQsupport") :
					if feature[node]!=None:
					#if feat=="support" : #dist[node]<=effectivelyNon0BLen:
						if feat=="support":
							if identicalTo!="":
								if supportForIdenticalSequences:
									if dist[node]<=effectivelyNon0BLen:
										stringList.append(str(feature[node]))
									else:
										stringList.append("1.0")
							else:
								stringList.append(str(feature[node]))
						else:
							stringList.append(str(feature[node]))
				#use this column to highlight which nodes could be placed (with probability above threshold) on the branch above the current node - used to highlight alternative placements of a given node on the tree.
				elif feat=="supportTo" and identicalTo=="":
					for iNode in range(len(feature[node])):
						stringList.append(namesInTree[tree.name[feature[node][iNode][0]]]+":"+str(feature[node][iNode][1]))
						if iNode<(len(feature[node])-1):
							stringList.append(",")
				# elif feat=="alternativePlacements":
				# 	for iNode in range(len(feature[node])):
				# 		stringList.append(namesInTree[tree.name[feature[node][iNode][0]]]+":"+str(feature[node][iNode][1]))
				# 		if iNode<(len(feature[node])-1):
				# 			stringList.append(",")
				elif feat=="mutationsInf" and identicalTo=="":
					for iNode in range(len(feature[node])):
						mutation=feature[node][iNode]
						stringList.append(allelesListExt[mutation[0]]+str(mutation[1])+allelesListExt[mutation[2]]+":"+str(mutation[3]))
						if iNode<(len(feature[node])-1):
							stringList.append(",")
				elif feat=="Ns":
					if identicalTo=="" or supportFor0Branches:
						for iNode in range(len(feature[node])):
							mutation=feature[node][iNode]
							if type(mutation)==int:
								stringList.append(str(mutation))
							else:
								stringList.append(str(mutation[0])+"-"+str(mutation[1]))
							if iNode<(len(feature[node])-1):
								stringList.append(",")
				elif feat=="errors":
					for iNode in range(len(feature[node])):
						mutation=feature[node][iNode]
						stringList.append(allelesListExt[mutation[0]]+str(mutation[1])+allelesListExt[mutation[2]]+":"+str(mutation[3]))
						if iNode<(len(feature[node])-1):
							stringList.append(",")
				elif feat=="lineage":
					stringList.append(feature[node])
				elif feat=="lineages":
					for lineageName in feature[node].keys():
						stringList.append(lineageName+":"+str(feature[node][lineageName]))
						stringList.append(",")
					stringList.pop()
				elif feat=="rootSupport" and identicalTo=="":
					if feature[node]!=None:
						stringList.append(str(feature[node]))
			# this column highlights nodes with support below threshold and with number of descendants above threshold
			elif feat=="supportGroup":
				if tree.support[node]!=None:
					if tree.support[node]<0.9:
						nDescString="nDesc<11_"
						if identicalTo=="":
							if tree.nDesc[node]>100000:
								nDescString="nDesc>100000_"
							elif tree.nDesc[node]>10000:
								nDescString="nDesc>10000_"
							elif tree.nDesc[node]>1000:
								nDescString="nDesc>1000_"
							elif tree.nDesc[node]>100:
								nDescString="nDesc>100_"
							elif tree.nDesc[node]>10:
								nDescString="nDesc>10_"
						if tree.support[node]<0.5:
							nDescString+="support<0.5"
						else:
							nDescString+="support<0.9"
					else:
						nDescString=""
					stringList.append(nDescString)
		stringList.append("\t")
	stringList[-1]="\n"
	#stringList.append("\n")
	return "".join(stringList)


#calculate number of descendants for each node.
def calculateNDesc(tree,node):
	children=tree.children
	nDesc=tree.nDesc
	minorSequences=tree.minorSequences
	nextLeaves=[node]
	for i in range(len(nDesc)):
		nDesc[i]=0
	while nextLeaves:
		nextNode=nextLeaves.pop()
		if children[nextNode]:
			if nDesc[children[nextNode][0]]:
				for c in children[nextNode]:
					nDesc[nextNode]+=nDesc[c]
			else:
				nextLeaves.append(nextNode)
				for c in children[nextNode]:
					nextLeaves.append(c)
		else:
			nDesc[nextNode]=1+len(minorSequences[nextNode])

#"Inverts" the alternativePlacements feature, recording, for each node, which other nodes might be placed on the branch above it. This is for better graphical representation in Taxonium.
def defineSupportedNodes(tree,node):
	children=tree.children
	alternativePlacements=tree.alternativePlacements
	tree.supportTo=[]
	supportTo=tree.supportTo
	for i in range(len(alternativePlacements)):
		supportTo.append([])
	nextLeaves=[node]
	while nextLeaves:
		nextNode=nextLeaves.pop()
		#number of times this node has been replaced through SPR during this round.
		if children[nextNode]:
			for c in children[nextNode]:
				nextLeaves.append(c)
		if alternativePlacements[nextNode]:
			for nodePair in alternativePlacements[nextNode]:
				supportTo[nodePair[0]].append((nextNode,nodePair[1]))


#write tsv file with metadata/annotations for tree nodes and taxa
def writeTSVfile(tree,node,file,namesInTree):
	children=tree.children
	up=tree.up
	name=tree.name
	minorSequences=tree.minorSequences
	featureNames={}
	if keepInputIQtreeSupports:
		featureNames['IQsupport']='IQsupport'
	if aBayesPlusOn:
		featureNames['support']='support'
		featureNames['rootSupport']='rootSupport'
		if networkOutput:
			# this one now turned off - still available in the nexus tree for graph representation. Now using "supportTo" so that one can more easily highlight the locations where a given node might be placed.
			#featureNames['alternativePlacements']='uncertainty'
			calculateNDesc(tree,node)
			# "invert" alternativePlacements so that each node has the list of other nodes that it supports as a placement.
			defineSupportedNodes(tree,node)
			#use this column to highlight nodes with support below threshold and with number of descendants above threshold
			featureNames['supportGroup']='supportGroup'
			# use this column to highlight which nodes could be placed (with probability above threshold) on the branch above the current node - used to highlight alternative placements of a given node on the tree.
			featureNames['supportTo']='supportTo'
	if estimateMAT:
		featureNames['mutationsInf']='mutationsInf'
		featureNames['Ns']='Ns'
	if usingErrorRate:
		featureNames['errors']='errors'
	if performLineageAssignment:
		featureNames['lineage']='lineage'
		featureNames['lineages']='lineages'
	featureList=list(featureNames.keys())
	file.write("strain"+"\t"+"collapsedTo")
	for feat in featureList:
		file.write("\t"+featureNames[feat])
	file.write("\n")
	#now write to file the features for each node of the tree.
	nextNode=node
	direction=0
	numLeaves=0
	while nextNode!=None:
		if children[nextNode]:
			if direction==0:
				nextNode=children[nextNode][0]
			elif direction==1:
				nextNode=children[nextNode][1]
				direction=0
			else:
				if aBayesPlusOn or estimateMAT or performLineageAssignment:
					file.write(tsvForNode(tree,nextNode,namesInTree[name[nextNode]],featureList,namesInTree))
				if up[nextNode]!=None:
					if children[up[nextNode]][0]==nextNode:
						direction=1
					else:
						direction=2
				nextNode=up[nextNode]
		else:
			numLeaves+=(1+len(minorSequences[nextNode]))
			if len(minorSequences[nextNode])>0:
				if supportForIdenticalSequences or performLineageAssignment:
					file.write(tsvForNode(tree,nextNode,namesInTree[name[nextNode]],featureList,namesInTree,identicalTo=namesInTree[name[nextNode]]+"_MinorSeqsClade"))
				else:
					file.write(tsvForNode(tree,None,namesInTree[name[nextNode]],featureList,namesInTree,identicalTo=namesInTree[name[nextNode]]+"_MinorSeqsClade"))
				for s2 in minorSequences[nextNode]:
					if supportForIdenticalSequences or performLineageAssignment:
						file.write(tsvForNode(tree,nextNode,namesInTree[s2],featureList,namesInTree,identicalTo=namesInTree[name[nextNode]]+"_MinorSeqsClade"))
					else:
						file.write(tsvForNode(tree,None,namesInTree[s2],featureList,namesInTree,identicalTo=namesInTree[name[nextNode]]+"_MinorSeqsClade"))
				if aBayesPlusOn or estimateMAT or performLineageAssignment:
					#if supportForIdenticalSequences or performLineageAssignment:
					file.write(tsvForNode(tree,nextNode,namesInTree[name[nextNode]]+"_MinorSeqsClade",featureList,namesInTree))
					#else:
					#	file.write(tsvForNode(tree,None,namesInTree[name[nextNode]]+"_MinorSeqsClade",featureList,namesInTree,identicalTo=namesInTree[name[nextNode]]))
			else:
				file.write(tsvForNode(tree,nextNode,namesInTree[name[nextNode]],featureList,namesInTree))
			if up[nextNode]!=None:
				if children[up[nextNode]][0]==nextNode:
					direction=1
				else:
					direction=2
			nextNode=up[nextNode]
	file.close()


#run topology search
timeTopology=0.0
threshValues=[]
threshNums=[]
threshStricts=[]
threshPlaces=[]
if fastTopologyInitialSearch:
	threshValues.append(thresholdLogLKtopologyInitial)
	threshNums.append(allowedFailsTopologyInitial)
	threshStricts.append(strictTopologyStopRulesInitial)
	threshPlaces.append(thresholdTopologyPlacementInitial)
if (inputTree=="" or largeUpdate or aBayesPlus):
	for i in range(numTopologyImprovements):
		threshValues.append(thresholdLogLKtopology)
		threshNums.append(allowedFailsTopology)
		threshStricts.append(strictTopologyStopRules)
		threshPlaces.append(thresholdTopologyPlacement)
nRounds=len(threshNums)


#Assign a core to each node of the tree, so that parallelization will run search for different nodes in different cores.
def assignCoreNumbers(tree,root,numCores):
	children=tree.children
	up=tree.up
	coreNum=[None]*len(up)
	tree.coreNum=coreNum
	node=root
	#direction 0 means from parent, direction 1 means from a child
	lastNode=None
	direction=0
	currentCore=0
	numNodes=0
	while node!=None:
		if direction==0:
			numNodes+=1
			coreNum[node]=currentCore
			currentCore+=1
			currentCore=currentCore%numCores
			if children[node]:
				node=children[node][0]
			else:
				lastNode=node
				node=up[node]
				direction=1
		else :
			if lastNode==children[node][0]:
				node=children[node][1]
				direction=0
			else:
				lastNode=node
				node=up[node]
				direction=1
	print("Assigned "+str(numCores)+" cores for "+str(numNodes)+" nodes.")


#count number of dirty nodes to decide if to parallelize or not.
def countDirtyNodes(tree,root):
	children=tree.children
	up=tree.up
	dirty=tree.dirty
	#dist=tree.dist
	node=root
	#direction 0 means from parent, direction 1 means from a child
	lastNode=None
	direction=0
	numNodes=0
	numDirty=0
	while node!=None:
		if direction==0:
			#if dist[node]:
			numNodes+=1
			if dirty[node]:
				numDirty+=1
			if children[node]:
				node=children[node][0]
			else:
				lastNode=node
				node=up[node]
				direction=1
		else :
			if lastNode==children[node][0]:
				node=children[node][1]
				direction=0
			else:
				lastNode=node
				node=up[node]
				direction=1
	print(str(numDirty)+" dirty nodes out of "+str(numNodes)+".")
	return numDirty,numNodes


if aBayesPlus:
	tree.support=[None]*len(tree.up)
	if networkOutput:
		tree.alternativePlacements=[]
		for i in range(len(tree.up)):
			tree.alternativePlacements.append([])


#Run rounds of SPR topological improvements
for nRound in range(nRounds):
	if aBayesPlus:
		aBayesPlusOn=True
	print("Starting topological impromevement traversal number "+str(nRound+1), flush=True)
	start=time()
	setAllDirty(tree,t1)
	
	reCalculateAllGenomeLists(tree,t1)
	newLk=calculateTreeLikelihood(tree,t1)
	print("Preliminary branch length optimization from LK: "+str(newLk))
	#TODO
	if HnZ:
		#print("Re-calculating the nDesc0")
		calculateNDesc0(tree,t1,checkExisting=True)
	improvement=traverseTreeToOptimizeBranchLengths(tree,t1)
	print("Branch length optimization, number of changes: "+str(improvement))
	subRound=0
	while subRound<20:
		if (not improvement):
			break
		subRound+=1
		improvement=traverseTreeToOptimizeBranchLengths(tree,t1)
	newLk=calculateTreeLikelihood(tree,t1)
	print("branch length finalization subround "+str(subRound+1)+" number of changes "+str(improvement)+" final LK: "+str(newLk), flush=True)
	timeForBranchOptimization=(time()-start)
	print("Time for updating branch lengths: "+str(timeForBranchOptimization))

	#TODO
	if HnZ:
		#print("Re-calculating the nDesc0")
		calculateNDesc0(tree,t1,checkExisting=True)

	#now update the topology
	start=time()
	setAllDirty(tree,t1)
	print("Starting tolopogical improvement round "+str(nRound+1)+" with allowed fails "+str(threshNums[nRound])+" and LK threshold "+str(threshValues[nRound]), flush=True)
	reCalculateAllGenomeLists(tree,t1)
	preLK=calculateTreeLikelihood(tree,t1)
	print("Likelihood before SPR moves: "+str(preLK), flush=True)

	if parallelize:
	#if parallelRounds[nRound]:
		#assign numbers to nodes so that each core will only investigate re-placement of its own nodes
		if nRound==0:
			assignCoreNumbers(tree,t1,numCores)
		paralleleInputs=[]
		for i in range(numCores):
			paralleleInputs.append((tree,t1,i,threshStricts[nRound],threshNums[nRound],threshValues[nRound],threshPlaces[nRound]))
		if __name__ == "__main__":
			with Pool() as pool:
				results = pool.map(startTopologyUpdatesParallel, paralleleInputs)
		for i in range(numCores-1):
			results[0].extend(results[i+1])
		results[0].sort(reverse=False,key=itemgetter(2))
		totalTimeFindingParent[0]+=time()-start
		print("Found proposed SPR moves, merged, and sorted.")
		setAllDirty(tree,t1,dirtiness=False)
		newRoot, improvement = applySPRMovesParallel(tree,results[0],strictTopologyStopRules=threshStricts[nRound],allowedFailsTopology=threshNums[nRound],thresholdLogLKtopology=threshValues[nRound],thresholdTopologyPlacement=threshPlaces[nRound])
	else:
		newRoot,improvement=startTopologyUpdates(tree,t1,checkEachSPR=debugging,strictTopologyStopRules=threshStricts[nRound],allowedFailsTopology=threshNums[nRound],thresholdLogLKtopology=threshValues[nRound],thresholdTopologyPlacement=threshPlaces[nRound])
	if newRoot!=None:
		t1=newRoot
	timeForUpdatingTopology=(time()-start)
	print("Time for round "+str(nRound+1)+" traversal of the tree for topology changes: "+str(timeForUpdatingTopology), flush=True)
	timeTopology+=timeForUpdatingTopology
	print("LK improvement apparently brought: "+str(improvement))
	reCalculateAllGenomeLists(tree,t1)
	postLK=calculateTreeLikelihood(tree,t1)
	print("Likelihood after SPR moves: "+str(postLK))
	
	print("Writing preliminary tree to file: "+outputFile+"_round"+str(nRound+1)+"_preliminary_tree.tree", flush=True)
	newickString=createNewick(tree,t1,binary=binaryTree,namesInTree=namesInTree,estimateMAT=False,networkOutput=False,aBayesPlusOn=False)
	file=open(outputFile+"_round"+str(nRound+1)+"_preliminary_tree.tree","w")
	file.write(newickString)
	file.close()
	

	#run improvements only on the nodes that have been affected by some changes in the last round, and so on
	start=time()
	subRound=0
	while subRound<20:
		print("Topological subround "+str(subRound+1), flush=True)
		numDirty,numNodes=countDirtyNodes(tree,t1)
		#if HnZ:
			#print("Re-calculating the nDesc0")
			#calculateNDesc0(tree,t1,checkExisting=True)
		if parallelize and (numDirty>0.1*numNodes):
			paralleleInputs=[]
			for i in range(numCores):
				paralleleInputs.append((tree,t1,i,threshStricts[nRound],threshNums[nRound],threshValues[nRound],threshPlaces[nRound]))
			if __name__ == "__main__":
				with Pool() as pool:
					results = pool.map(startTopologyUpdatesParallel, paralleleInputs)
			for i in range(numCores-1):
				results[0].extend(results[i+1])
			results[0].sort(reverse=False,key=itemgetter(2))
			totalTimeFindingParent[0]+=time()-start
			print("Found proposed SPR moves, merged, and sorted.")
			setAllDirty(tree,t1,dirtiness=False)
			newRoot, improvement = applySPRMovesParallel(tree,results[0],strictTopologyStopRules=threshStricts[nRound],allowedFailsTopology=threshNums[nRound],thresholdLogLKtopology=threshValues[nRound],thresholdTopologyPlacement=threshPlaces[nRound])
		else:
			newRoot,improvement=startTopologyUpdates(tree,t1,checkEachSPR=debugging,strictTopologyStopRules=threshStricts[nRound],allowedFailsTopology=threshNums[nRound],thresholdLogLKtopology=threshValues[nRound],thresholdTopologyPlacement=threshPlaces[nRound])
		if newRoot!=None:
			t1=newRoot
		print("LK improvement apparently brought: "+str(improvement), flush=True)
		if not noSubroundTrees:
			print("Writing preliminary tree to file: "+outputFile+"_round"+str(nRound+1)+"_subround"+str(subRound+1)+"_preliminary_tree.tree", flush=True)
			newickString=createNewick(tree,t1,binary=binaryTree,namesInTree=namesInTree,estimateMAT=False,networkOutput=False,aBayesPlusOn=False)
			file=open(outputFile+"_round"+str(nRound+1)+"_subround"+str(subRound+1)+"_preliminary_tree.tree","w")
			file.write(newickString)
			file.close()
		if improvement<thresholdLogLKTopologySubRoundImprovement:
			break
		subRound+=1
	timeForUpdatingTopology=(time()-start)
	reCalculateAllGenomeLists(tree,t1)
	postLK=calculateTreeLikelihood(tree,t1)
	print("Likelihood after SPR subrounds: "+str(postLK), flush=True)
	print("Time for the subrounds of this traversal of the tree: "+str(timeForUpdatingTopology), flush=True)
	timeTopology+=timeForUpdatingTopology

	#if estimating error rates, repeating EM after shallow topological search
	if estimateErrorRate or estimateSiteSpecificErrorRate:
		oldLK=float("-inf")
		newLk=calculateTreeLikelihood(tree,t1)
		print("Initial LK before error rates EM: "+str(newLk), flush=True)
		mutMatrixGlobal, siteRates, errorRateGlobal, errorRates = expectationMaximizationCalculationRates(tree,t1)
		updateMutMatrices(mutMatrixGlobal,siteRates=siteRates)
		updateErrorRates(errorRateGlobal,errorRates=errorRates)
		reCalculateAllGenomeLists(tree,t1)
		newLk=calculateTreeLikelihood(tree,t1)
		print("Initial LK after first error rates EM: "+str(newLk))
		numEMsteps=0
		while (newLk-oldLK>1.0) and numEMsteps<20:
			setAllDirty(tree,t1)
			improvement=traverseTreeToOptimizeBranchLengths(tree,t1)
			reCalculateAllGenomeLists(tree,t1)
			newLkBranch=calculateTreeLikelihood(tree,t1)
			print("Updated "+str(improvement)+" branch lengths leading to LK "+str(newLkBranch), flush=True)

			mutMatrixGlobal, siteRates, errorRateGlobal, errorRates = expectationMaximizationCalculationRates(tree,t1)
			updateMutMatrices(mutMatrixGlobal,siteRates=siteRates)
			updateErrorRates(errorRateGlobal,errorRates=errorRates)
			reCalculateAllGenomeLists(tree,t1)
			oldLK=newLk
			newLk=calculateTreeLikelihood(tree,t1)
			print("New LK step "+str(numEMsteps)+": "+str(newLk))
			numEMsteps+=1
		if rateVariation:
			meanRate, varianceRate = variance(siteRates)
			print("Rate variation variance: "+str(varianceRate))
		if errorRateSiteSpecific:
			meanErrorRate, varianceErrorRate = variance(errorRates)
			print("Error rate variation, mean: "+str(meanErrorRate)+" , variance: "+str(varianceErrorRate))
		print("Error rate: "+str(errorRateGlobal), flush=True)


	# update just branch lengths
	start=time()
	reCalculateAllGenomeLists(tree,t1)
	newLk=calculateTreeLikelihood(tree,t1)
	setAllDirty(tree,t1)
	print(" branch length optimization starting from LK "+str(newLk), flush=True)
	improvement=traverseTreeToOptimizeBranchLengths(tree,t1)
	print("Branch length optimization round 1, number of changes: "+str(improvement))
	subRound=0
	while subRound<20:
		if (not improvement):
			break
		subRound+=1
		improvement=traverseTreeToOptimizeBranchLengths(tree,t1)
	reCalculateAllGenomeLists(tree,t1)
	newLk=calculateTreeLikelihood(tree,t1)
	print("branch length finalization subround "+str(subRound+1)+" number of changes "+str(improvement)+" final LK: "+str(newLk))
	timeForBranchOptimization=(time()-start)
	print("Time for updating branch lengths: "+str(timeForBranchOptimization), flush=True)

	#writing to output the substitution model (possibly with rate variation and error rates)
	if nRound<(nRounds-1):
		fileNameAdd="_round"+str(nRound+1)
	else:
		fileNameAdd=""
	fileName=outputFile+fileNameAdd+"_subs.txt"
	file=open(fileName,"w")
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
	print("Current Substitution matrix:")
	print(mutMatrixGlobal)
	print("All substitution rates and error probabilities written to file "+fileName)

	#calculate total likelihood
	totalLK=calculateTreeLikelihood(tree,t1)
	print("totalLK round "+str(nRound+1)+": "+str(totalLK), flush=True)
	fileName=outputFile+fileNameAdd+"_LK.txt"
	file=open(fileName,"w")
	print("totalLK written to file "+fileName)
	file.write(str(totalLK)+"\n")
	file.close()

	if estimateErrors:
		print("Estimating errors in input alignment round "+str(nRound+1), flush=True)
		fileName=outputFile+fileNameAdd+"_estimatedErrors.txt"
		file=open(fileName,"w")
		calculateErrorProbabilities(tree,t1,file,minErrorProb,namesInTree)
		file.close()
		print("Errors estimated, written to file "+fileName)

	if estimateMAT:
		expectationMaximizationCalculationRates(tree,t1,trackMutations=True)

	newickString=createNewick(tree,t1,binary=binaryTree,namesInTree=namesInTree,estimateMAT=estimateMAT,networkOutput=networkOutput,aBayesPlusOn=aBayesPlusOn)

	if aBayesPlus or estimateMAT:
		fileName=outputFile+fileNameAdd+"_nexusTree.tree"
		file=open(fileName,"w")
		file.write("#NEXUS\nbegin taxa;\n	dimensions ntax="+str(len(namesInTree))+";\n	taxlabels\n")
		for name in namesInTree:
			file.write("	"+name+"\n")
		file.write(";\nend;\n\nbegin trees;\n	tree TREE1 = [&R] ")
		file.write(newickString)
		file.write("\nend;\n")
		file.close()
		file=open(outputFile+fileNameAdd+"_metaData.tsv","w")
		writeTSVfile(tree,t1,file,namesInTree)
		file.close()
		print("Nexus tree written to file "+fileName, flush=True)
		newickString=createNewick(tree,t1,binary=binaryTree,namesInTree=namesInTree,estimateMAT=False,networkOutput=False,aBayesPlusOn=False)

	fileName=outputFile+fileNameAdd+"_tree.tree"
	file=open(fileName,"w")
	file.write(newickString)
	file.close()
	
	print("Tree written to file "+fileName, flush=True)


#If no rounds were run, still output estimates
if nRounds==0:
	#writing to output the substitution model (possibly with rate variation and error rates)
	fileNameAdd=""
	fileName=outputFile+fileNameAdd+"_subs.txt"
	file=open(fileName,"w")
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
	print("Substitution matrix:")
	print(mutMatrixGlobal)
	print("All substitution rates and error probabilities written to file "+fileName)

	#calculate total likelihood
	totalLK=calculateTreeLikelihood(tree,t1)
	print("totalLK: "+str(totalLK))
	fileName=outputFile+fileNameAdd+"_LK.txt"
	file=open(fileName,"w")
	print("totalLK written to file "+fileName)
	file.write(str(totalLK)+"\n")
	file.close()

	if estimateErrors:
		print("Estimating errors in input alignment")
		fileName=outputFile+fileNameAdd+"_estimatedErrors.txt"
		file=open(fileName,"w")
		calculateErrorProbabilities(tree,t1,file,minErrorProb,namesInTree)
		file.close()
		print("Errors estimated, written to file "+fileName)


	if estimateMAT:
		expectationMaximizationCalculationRates(tree,t1,trackMutations=True)

	newickString=createNewick(tree,t1,binary=binaryTree,namesInTree=namesInTree)

	if aBayesPlus or estimateMAT:
		fileName=outputFile+fileNameAdd+"_nexusTree.tree"
		file=open(fileName,"w")
		file.write("#NEXUS\nbegin taxa;\n	dimensions ntax="+str(len(namesInTree))+";\n	taxlabels\n")
		for name in namesInTree:
			file.write("	"+name+"\n")
		file.write(";\nend;\n\nbegin trees;\n	tree TREE1 = [&R] ")
		file.write(newickString)
		file.write("\nend;\n")
		file.close()
		file=open(outputFile+fileNameAdd+"_metaData.tsv","w")
		writeTSVfile(tree,t1,file,namesInTree)
		file.close()
		print("Nexus tree written to file "+fileName)
		aBayesPlusOn=False
		networkOutput=False
		estimateMAT=False
		newickString=createNewick(tree,t1,binary=binaryTree,namesInTree=namesInTree)

	fileName=outputFile+fileNameAdd+"_tree.tree"
	file=open(fileName,"w")
	file.write(newickString)
	file.close()
	print("Tree written to file "+fileName)



print("Number of final references in the MAT: "+str(numRefs[0]), flush=True)
print("Time spent finding placement nodes: "+str(timeFinding))
print("Time spent placing samples on the tree: "+str(timePlacing))
print("Time spent in total updating the topology and branch lengths: "+str(timeTopology))
print("Of which looking for placements for better topologies: "+str(totalTimeFindingParent[0]))

exit()














