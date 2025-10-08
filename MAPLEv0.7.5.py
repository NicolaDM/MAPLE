import sys
from math import log
import argparse
from time import time
import os.path
from operator import itemgetter

#©EMBL-European Bioinformatics Institute, 2021-2025
#Developd by Nicola De Maio, with contributions from Myrthe Willemsen.

# MAPLE code to estimate a tree by maximum likelihood from a MAPLE format input.

# TODO currently time trees ar being implemented
#TODO TODO TODO infer posterior node time ranges and write time tree to output file.
#TODO TODO TODO investigate if any input time is plausible wrong and mask it.

#long-term plans
#TODO phylogeography
#TODO indels, alignment via alignment-free pair-HMM
#TODO Recombination via pair-HMM based SPR search of second parent
#TODO relativistic phylogenetics for arbitrary trees

parser = argparse.ArgumentParser(description='Estimate a tree from a diff format and using iterative approximate maximum likelihood sample placement.')
#important options
parser.add_argument('--input',default="MAPLE_input.txt", help='Input MAPLE file name; should contain first the reference genome and then the difference of all samples with respect to the reference.')
parser.add_argument('--reference',default="", help='Optional input reference file name. By default it assumes instead that the reference is part of the MAPLE format input.')
parser.add_argument("--model", help="Which substitution model should be used. Allowed models so far are JC, GTR (default) or UNREST.", default="GTR")
parser.add_argument('--output',default="MAPLE_output", help='Output path and identifier to be used for newick output file.')
parser.add_argument('--inputTree',default="", help='Input newick tree file name; this is optional, and is used for online inference (or for Robinson-Foulds distance calculation if option --inputRFtrees is also used).')
parser.add_argument('--inputRates',default="", help='Name of input containing pre-estimated substitution rates and possibly other model parameters; this is optional, and is only used for online inference. The same format as the MAPLE output is expected, and information about all parameters of the selected substitution model is expected.')
parser.add_argument("--largeUpdate", help="When using option --inputTree to do online inference, by default the tree is only updated locally where the new sequences are inserted. Use this option instead to perform a thorough update of the phylogeny.", action="store_true")
parser.add_argument('--inputRFtrees',default="", help='file name with input newick trees to be compared with the input tree in option --inputTree to calculate Robinson-Foulds distances; this option will turn off the normal MAPLE estimation - only RF distances will be calculated. Each newick tree contained in this file will be compared to the single tree specified with option --inputTree .')
parser.add_argument("--overwrite", help="Overwrite previous results if already present.", action="store_true")
parser.add_argument("--fast", help="Set parameters so to run tree inference faster; this will be less accurate in cases of high complexity, for example with recombination, sequencing errors, etc. It will overrule user choices for options --thresholdLogLK , --thresholdLogLKtopology , --allowedFails , --allowedFailsTopology .", action="store_true")
parser.add_argument("--rateVariation", help="Estimate and use rate variation: the model assumes one rate per site, and the rates are assumed independently (no rate categories). This might cause overfitting if the dataset is not large enough, but in any case one would probably only use MAPLE for large enough datasets.", action="store_true")
parser.add_argument("--estimateMAT", help="Estimate mutation events on the final tree, and write them on the nexus tree.", action="store_true")
parser.add_argument("--doNotImproveTopology", help="Do not perform SPR moves, despite searching for them; this is useful if one wants to analyse a tree and calculate branch supports without changing the input tree.", action="store_true")
parser.add_argument("--saveInitialTreeEvery",help="Every this many samples placed (default 50,000) save the current tree to file. This way if there is any problem, the building of the initial tree can be restarted from the last saved tree.",  type=int, default=50000)
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
parser.add_argument("--nonStrictStopRules", help="If specified, then during the initial placement stage, a slower, non-strict rule for stopping the placement search is applied: the search is stopped if enough consecutive LK worsening are observed, AND if LK is below the considered threshold.", action="store_true")
parser.add_argument("--strictTopologyStopRules", help="If specified, then during the topological improvement stage, a faster, strict rule for stopping the SPR search is applied: the search is stopped if enough consecutive LK worsening are observed, OR if LK is below the considered threshold.", action="store_true")
parser.add_argument("--thresholdDiffForUpdate",help="Consider the probability of a new partial changed if the difference between old and new is above this threshold.",  type=float, default=0.00001)
parser.add_argument("--thresholdFoldChangeUpdate",help="Consider the probability of a new partial changed, if the fold difference between old and new if above this threshold.",  type=float, default=1.01)
parser.add_argument("--thresholdLogLKconsecutivePlacement",help="logLK difference threshold to consider something as a significant decrease in log-LK when considering consecutive likelihood decreases.",  type=float, default=1.0)#0.01
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
#lineage assignment from reference genomes which are not in the tree yet - NHAN
parser.add_argument('--lineageRefs',default="", help='give path and name to an alignment file (in MAPLE format) containing reference genomes, each represents one lineage. When using this option, option --inputTree should also be used. Then MAPLE will find the best placement for each lineage reference (note that these lineage references do not exist in the input tree). Each sample is assigned a lineage same as its closest reference parent.')
parser.add_argument('--lineageRefsThresh',default=0.2, help='The threshold (in terms of #mutation) to check whether a reference lineage genome could be considered as the parent of a subtree. Default: 0.2 mutation', type = float)
parser.add_argument('--lineageRefsSupportThresh',default=0.95, help='A lineage will be assigned to a subtree only if the SPRTA support for that lineage placement exceeds this threshold. Default: 0.95', type = float)
parser.add_argument('--allowMultiLineagesPerNode', help='When a node is selected as the best placements for multiple lineages, whether we allow assigning all of these lineages (or only the closest lineage) to the subtree. Default: assigning the closest lineage', action="store_true")
# find placements (in an input tree) for new samples (without changing the input tree ~ find placements only) - NHAN
parser.add_argument('--findSamplePlacements', help='Find placements (in an input tree) for new samples (without changing the input tree ~ find placements only)', action="store_true")
parser.add_argument('--threshMutation',default=0.01, help='Threshold for detecting a mutation from an entry O to a new nucleotide X: if the probability of X in O < threshMutation, the change is counted as a mutation. Default: 0.01', type = float)


#rarer options
parser.add_argument("--defaultBLen",help="Default length of branches, for example when the input tree has no branch length information.",  type=float, default=0.000033)
parser.add_argument("--normalizeInputBLen",help="For the case the input tree has branch lengths expressed not in a likelihood fashion (that is, expected number of substitutions per site), then multiply the input branch lengths by this factor. Particularly useful when using parsimony-based input trees.",  type=float, default=1.0)
parser.add_argument("--multipleInputRFTrees", help="Use this option if the file specified by option --inputRFtrees contains multiple trees - otherwise only read the first tree from that file.", action="store_true")
parser.add_argument("--debugging", help="Test that likelihoods are calculated and updated as expected - time consuming and only meant for small trees for debugging purposes.", action="store_true")
parser.add_argument("--onlyNambiguities", help="Treat all ambiguities as N (total missing information).", action="store_true")
parser.add_argument("--nonBinaryTree", help="Write output tree with multifurcations - by default the tree is written as binary so to avoid problems reading the tree in other software.", action="store_true")
parser.add_argument("--writeTreesToFileEveryTheseSteps", help="By default, don't write intermediate trees to file. If however a positive integer is specified with this option, intermediate trees will be written to file every this many topological changes.",  type=int, default=0)
parser.add_argument("--writeLKsToFileEveryTheseSteps", help="By default, don't write likelihoods of intermediate trees to file. If however a positive integer is specified with this option, likelihoods of intermediate trees will be written to file every this many topological changes.",  type=int, default=0)
parser.add_argument("--noSubroundTrees", help="Do not write to file subround trees.", action="store_true")
parser.add_argument("--doNotOptimiseBLengths", help="Do not optimise the branch lengths of the tree (useful if the input tree is already optimal or doesn't need changing).", action="store_true")
parser.add_argument("--forgetInputTreeInternalNodeNames", help="Do keep the names of the internal tree nodes in the input tree.", action="store_true")
#error model options
parser.add_argument("--estimateErrorRate", help="Estimate a single error rate for the whole genome. Input value is used as starting value", action="store_true")
parser.add_argument("--estimateSiteSpecificErrorRate", help="Estimate a separate error rate for each genome. Input value is used as starting value", action="store_true")
parser.add_argument("--errorRateInitial", help="Initial value used for estimating the error rate. The default is the inverse of the reference genome length (one error expected per genome).", type=float, default=0.0)
parser.add_argument("--errorRateFixed", help="Fix the error rate to a given input value", type=float, default=0.0)
parser.add_argument("--errorRateSiteSpecificFile", help="provide a file path to a file that contains the siteSpecific error rates", type=str, default=None)
parser.add_argument("--estimateErrors", help="Estimate erroneous positions in the input sequences. This option is only allowed if option --estimateSiteSpecificErrorRate or errorRateSiteSpecificFile is also used.", action="store_true")
parser.add_argument("--minErrorProb", help="Minimum error probability to be reported when using option --estimateErrors", type=float, default=0.01)
#SPRTA options
parser.add_argument("--SPRTA", help="Calculate branch support values with a modification of the aBayes approach of Anisimova et al 2011 (De Maio et al 2024, https://doi.org/10.1101/2024.10.21.619398 ).", action="store_true")
parser.add_argument("--aBayesPlus", help="Synonym for option --SPRTA.", action="store_true")
parser.add_argument("--networkOutput", help="Include in the output tree the alternative branches with their support (like in a weighted phylogenetic network).", action="store_true")
parser.add_argument("--minBranchSupport", help="Minimum branch support to be considered when using option --networkOutput", type=float, default=0.01)
parser.add_argument("--supportFor0Branches", help="When calculating branch support, also consider supports of branches of length 0, such as samples less informative than other samples which might have multiple placements on the tree.", action="store_true")
parser.add_argument("--minMutProb", help="Minimum mutation probability to be written to output when using option --estimateMAT", type=float, default=0.01)
parser.add_argument("--keepInputIQtreeSupports", help="Assumes the input tree is from IQTREE, reads the support values on the tree branches, and prints the same values in the output nexus tree.", action="store_true")
#HorseNotZebra modifiers option
parser.add_argument("--HnZ", help="By default (0), don't use HnZ modifiers. If 1, use the topological HorseNotZebra modifier; if 2, use the sampling likelihood HnZ modifier.",  type=int, default=0)
#Time tree options
parser.add_argument('--datesFile',default=None, help='Input dates file name; should contain sample names in the first column and dates in the second column.')
parser.add_argument("--intervalLength", help="Number of days forming a time interval (default 7). smaller numbers (minimum 1) correspond to finer grained and more computationally intensive inference. Only used with option --datesFile.",  type=int, default=7)
parser.add_argument('--strainName',default="strain", help='Column name for sample names in the input time metadata (for option --datesFile). Default is \"strain\".')
parser.add_argument('--dateName',default="date", help='Column name for date values in the input time metadata (for option --datesFile). Default is \"date\".')
parser.add_argument("--minSamplingYear", help="Earliest possible sampling year - earlier dates in the metadata will be considered errors. Only used with option --datesFile.",  type=int, default=None)
parser.add_argument("--maxSamplingYear", help="Latest possible sampling year - later dates in the metadata will be considered errors. Only used with option --datesFile.",  type=int, default=None)
parser.add_argument("--mutRate", help="Initial value of the mutation rate (per genome per day). Default is 0.09, close to SARS-CoV-2. This is also used to inform a prior on the mutation rate. Only used with option --datesFile.", type=float, default=0.09)
parser.add_argument("--minMutRate", help="Minimum value of the mutation rate (per genome per day). Default is 0.03, used to prevent too low estimates and too long runtime. Only used with option --datesFile.", type=float, default=0.03)
parser.add_argument("--timeProbThreshold", help="Minimum, threshold probability for ignoring unlikely dates (default is 0.01). Only used with option --datesFile. Smaller values might result in more accurate and computationally intensive inference.", type=float, default=0.0001)
parser.add_argument("--minNumSamplesForMutRate",help="When creating the initial tree, start re-estimating the time-scaled mutation rate only after these many samples have been added (better not to decrease this too much to avoid overfitting).",  type=int, default=1000)
args = parser.parse_args()

onlyNambiguities=args.onlyNambiguities
thresholdProb=args.thresholdProb
debugging=args.debugging
inputFile=args.input
outputFile=args.output
refFile=args.reference
lineageRefs=args.lineageRefs
lineageRefsThresh = args.lineageRefsThresh
lineageRefsSupportThresh = args.lineageRefsSupportThresh
allowMultiLineagesPerNode = args.allowMultiLineagesPerNode
performLineageAssignmentByRefPlacement = (lineageRefs != "")
findSamplePlacements = args.findSamplePlacements
threshMutation = args.threshMutation
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
datesFile=args.datesFile
intervalLength=args.intervalLength
strainName=args.strainName
dateName=args.dateName
minSamplingYear=args.minSamplingYear
maxSamplingYear=args.maxSamplingYear
timeProbThreshold=args.timeProbThreshold
minNumSamplesForMutRate=args.minNumSamplesForMutRate
forgetInputTreeInternalNodeNames=args.forgetInputTreeInternalNodeNames
if datesFile!=None:
	doTimeTree=True
else:
	doTimeTree=False
if __name__ == "__main__":
	if doTimeTree:
		mutRate=args.mutRate
		mutRate*=intervalLength
		minMutRate=args.minMutRate
		minMutRate*=intervalLength
	else:
		mutRate=None
		minMutRate=None
else:
	mutRate=None
	minMutRate=None

numCores=args.numCores
parallelize=False
if numCores>1:
	parallelize=True
	from multiprocessing import Pool
	if __name__ == "__main__":
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
doNotOptimiseBLengths=args.doNotOptimiseBLengths

assignmentFile=args.assignmentFile
assignmentFileCSV=args.assignmentFileCSV
inputNexusTree=args.inputNexusTree
maxReplacements=args.maxReplacements
writeTreesToFileEveryTheseSteps=args.writeTreesToFileEveryTheseSteps
writeLKsToFileEveryTheseSteps=args.writeLKsToFileEveryTheseSteps
reRoot=args.reRoot

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
rateVariation=args.rateVariation
mutMatrixGlobal=None
rootFreqsLogErrorCumulative=None
totError=None
cumulativeErrorRate=None
if __name__ == "__main__":
	errorRateGlobal=0.0
	usingErrorRate=False
	errorRateSiteSpecific=False
	useRateVariation=False
else:
	usingErrorRate=None
	errorRateGlobal=None
	errorRateSiteSpecific=None
	useRateVariation=None

aBayesPlus=args.aBayesPlus or args.SPRTA
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
if __name__ == "__main__":
	aBayesPlusOn=False
else:
	aBayesPlusOn=aBayesPlus
if aBayesPlus or doTimeTree or performLineageAssignmentByRefPlacement or findSamplePlacements:
	from math import exp

warnedBLen=[False]
warnedTotDiv=[False]
sumChildLKs=[0.0]
numChildLKs=[0]
useFixedThresholdLogLKoptimizationTopology=args.useFixedThresholdLogLKoptimizationTopology
totDivFromRef=[0.0]

nonMutRates=None
cumulativeRate=None
siteRates=None
errorRates=None
mutMatrices=None

#HnZ check and HnZvector definition and update 
HnZ=args.HnZ
if HnZ>2 or HnZ<0:
	print("Option --HnZ only allows values 0 (no HnZ), 1 (topological HnZ) or 2 (abundance HnZ).")
	raise Exception("exit")
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
		#number of branches descending from node after collapsing 0-length branches
		self.nDesc0=[]
		self.probVectTime=[]
		self.probVectUpRightTime=[]
		self.probVectUpLeftTime=[]
		self.probVectTotUpTime=[]
		self.dateData=[]
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
		if HnZ:
			self.nDesc0.append(1)
		if doTimeTree:
			self.probVectTime.append(None)
			self.probVectUpRightTime.append(None)
			self.probVectUpLeftTime.append(None)
			self.probVectTotUpTime.append(None)
			self.dateData.append(False)

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



dates=None
if doTimeTree:
	if __name__ == "__main__":
		from calendar import isleap
	from math import floor, ceil


	poissonCoeff=[[1.0]]
	# TODO get cached Poisson probabilities, and add more cached values if needed. 
	# TODO TODO TODO these values depend on mutRate, so if this changes, the probabilities also change. TODO TODO
	def getPoissonCoeff(b,t,mutRate):
		if t==0:
			if b>0:
				return 0.0
			else:
				return 1.0
		if t<0:#TODO remove after testing
			print("Negative time requested for Poisson probability?")
			raise Exception("exit")
		if b<0:#TODO remove after testing
			print("Negative mutations in Poisson probability?")
			raise Exception("exit")
		try:
			return poissonCoeff[t][b]
		except:
			try:
				poiT=poissonCoeff[t]
			except:
				lenPoiT=len(poissonCoeff)
				for i in range(t+1-lenPoiT):
					poissonCoeff.append([exp(-(lenPoiT+i)*mutRate)])
				poiT=poissonCoeff[t]
			lenPoiB=len(poiT)
			for i in range(b+1-lenPoiB):
				poiT.append(poiT[-1]*t*mutRate/float(lenPoiB+i))
			return poissonCoeff[t][b]


	# TODO get the time partial likelihoods as they change by moving along a branch
	# TODO TODO TODO make sure the resulting vector is normalised before being stored in the tree.
	def getPartialVecTime(probVectTime, mutRate, totLen, upNode=False, otherMinT=None, otherMaxT=None, tryMin=None):
		"""
		moving the time probability vector probVectTime along a branch of genetic distance totLen.
		:param probVectTime: the time probability vector to be moved.
		:param mutRate: the mutation rate per time interval per site.
		:param totLen: branch length in terms of genetic distance
		:param upNode: should be True if moving downward from above the branch.
		:param otherMaxT: if there are time constraints from another branch, impose them with this parameter.
		:param otherMinT: if there are time constraints from another branch, impose them with this parameter.
		:return: the new time probability vector.
		"""
		if probVectTime==None:
			return None
		elif len(probVectTime)==1:
			minT=probVectTime[0]
			maxT=probVectTime[0]
			probV=None
		elif len(probVectTime)==2:
			minT=probVectTime[1]
			maxT=probVectTime[0]
			probV=None
		else:
			minT=probVectTime[1]
			maxT=probVectTime[0]
			probV=probVectTime[2]
		if minT>maxT: #TODO remove after testing?
			print("Improper formatting of input time vector in getPartialVecTime")
			print(probVectTime)
			raise Exception("exit")

		newProbV=[]
		if upNode:
			if abs(round(totLen)-totLen)>timeProbThreshold:
				totLens=[floor(totLen),ceil(totLen)]
				totLenProbs=[ceil(totLen)-totLen,totLen-floor(totLen)]
				if totLens[0]>0:
					newMinT=minT+1
				else:
					newMinT=minT
			else:
				totLens=None
				totLen=round(totLen)
				if totLen:
					newMinT=minT+1
				else:
					newMinT=minT
			#if otherMinT!=None and otherMinT>newMinT:
			#	newMinT=otherMinT
			if otherMaxT!=None and otherMaxT<newMinT:
				print("Merging incompatible time vectors in getPartialVecTime")
				return otherMaxT
			highestProb=0.0
			lastProb=1.0
			currentT=newMinT
			while (otherMaxT!=None and currentT<=otherMaxT) or (otherMaxT==None and lastProb>=highestProb*timeProbThreshold):
				lastProb=0.0
				if probV:
					if totLens:
						for timeUp in range(minT,min(currentT,maxT)+1):					
							lastProb+=( getPoissonCoeff(totLens[0],currentT-timeUp,mutRate)*totLenProbs[0] + getPoissonCoeff(totLens[1],currentT-timeUp,mutRate)*totLenProbs[1] ) * probV[maxT-timeUp]
					else:
						for timeUp in range(minT,min(currentT,maxT)+1):	
							lastProb+=getPoissonCoeff(totLen,currentT-timeUp,mutRate)*probV[maxT-timeUp]
				else:
					if totLens:
						for timeUp in range(minT,min(currentT,maxT)+1):
							lastProb+=( getPoissonCoeff(totLens[0],currentT-timeUp,mutRate)*totLenProbs[0] + getPoissonCoeff(totLens[1],currentT-timeUp,mutRate)*totLenProbs[1] ) 
					else:
						for timeUp in range(minT,min(currentT,maxT)+1):
							lastProb+=getPoissonCoeff(totLen,currentT-timeUp,mutRate)
				if lastProb>highestProb:
					highestProb=lastProb
				newProbV.append(lastProb)
				currentT+=1
			newProbV.reverse()
			newMaxT=currentT-1
		else:
			if otherMinT==None:
				otherMinT=float("-inf")
			elif tryMin!=None and otherMinT>(tryMin-1):
				otherMinT=tryMin-1
			if abs(round(totLen)-totLen)>timeProbThreshold:
				totLens=[floor(totLen),ceil(totLen)]
				totLenProbs=[ceil(totLen)-totLen,totLen-floor(totLen)]
				if totLens[0]>0:
					newMaxT=maxT-1
				else:
					newMaxT=maxT
			else:
				totLens=None
				totLen=round(totLen)
				if totLen:
					newMaxT=maxT-1
				else:
					newMaxT=maxT
			if otherMaxT!=None and otherMaxT<newMaxT:
				newMaxT=otherMaxT
			if otherMinT>newMaxT:
				print("Merging incompatible time vectors in getPartialVecTime 2")
				print(probVectTime)
				#print(otherMinT)
				#print(otherMaxT)
				return newMaxT
			highestProb=0.0
			lastProb=1.0
			currentT=newMaxT
			while (otherMinT!=float("-inf") and currentT>=otherMinT) or ((tryMin!=None) and (currentT>=(tryMin-1))) or (otherMinT==float("-inf") and lastProb>=highestProb*timeProbThreshold):
				lastProb=0.0
				if 	probV:
					if totLens:
						for timeDown in range(max(currentT,minT),maxT+1):
							lastProb+=( getPoissonCoeff(totLens[0],timeDown-currentT,mutRate)*totLenProbs[0] + getPoissonCoeff(totLens[1],timeDown-currentT,mutRate)*totLenProbs[1] ) * probV[maxT-timeDown]
					else:
						for timeDown in range(max(currentT,minT),maxT+1):
							lastProb+=getPoissonCoeff(totLen,timeDown-currentT,mutRate)*probV[maxT-timeDown]
				else:
					if totLens:
						for timeDown in range(max(currentT,minT),maxT+1):
							lastProb+=( getPoissonCoeff(totLens[0],timeDown-currentT,mutRate)*totLenProbs[0] + getPoissonCoeff(totLens[1],timeDown-currentT,mutRate)*totLenProbs[1] ) 
					else:
						for timeDown in range(max(currentT,minT),maxT+1):
							lastProb+=getPoissonCoeff(totLen,timeDown-currentT,mutRate)
				if lastProb>highestProb:
					highestProb=lastProb
				newProbV.append(lastProb)
				currentT-=1
			newMinT=currentT+1
		return (newMaxT,newMinT,newProbV)


	# TODO TODO TODO
	# traverse the tree upward to recalculate upRight or upLeft time likelihoods above the starting node to extend the vectors so to be 
	# compatible with an exceptionally early down time vector (with minimum value newMin)
	def resolveTimeInconsistency(tree,node,newMin,mutRate):
		print("Resolving Time inconsistency, newMin "+str(newMin)+" node "+str(node))
		#newMin=-20
		#print("Actually let's do newMin = -20")
		probVectTime=tree.probVectTime
		probVectUpRightTime=tree.probVectUpRightTime
		probVectUpLeftTime=tree.probVectUpLeftTime
		dist=tree.dist
		up=tree.up
		children=tree.children
		nodesToUpdate=[node]
		while nodesToUpdate:
			currentNode=nodesToUpdate[-1]
			print(currentNode)
			if dist[currentNode]:
				newMin-=1
			pNode=up[currentNode]
			if currentNode==children[pNode][0]:
				siblingVector=probVectTime[children[pNode][1]]
				siblingDist=dist[children[pNode][1]]
			else:
				siblingVector=probVectTime[children[pNode][0]]
				siblingDist=dist[children[pNode][0]]
			#print("\n node "+str(currentNode)+" parent "+str(pNode)+" sibling vector:")
			#print(siblingVector)
			if up[pNode]!=None:
				if pNode==children[up[pNode]][0]:
					upVect=probVectUpRightTime[up[pNode]]
				else:
					upVect=probVectUpLeftTime[up[pNode]]
				#print("upVect:")
				#print(upVect)
				if upVect!=None and upVect[1]>(newMin-1):
					#print("consider also parent")
					nodesToUpdate.append(pNode)
				else:
					newVectUpTime=mergeVectorsTime(upVect,dist[pNode],siblingVector,siblingDist,mutRate,isUpDown=True,tryMin=newMin)
					#print("Vector extension successful:")
					#print(newVectUpTime)
					if currentNode==children[pNode][0]:
						probVectUpRightTime[pNode]=newVectUpTime
					else:
						probVectUpLeftTime[pNode]=newVectUpTime
					nodesToUpdate.pop()
					break
			else:
				newVectUpTime=rootVectorTime(siblingVector,siblingDist,mutRate,tryMin=newMin)
				#print("Reached root, new vector:")
				#print(newVectUpTime)
				if currentNode==children[pNode][0]:
					probVectUpRightTime[pNode]=newVectUpTime
				else:
					probVectUpLeftTime[pNode]=newVectUpTime
				nodesToUpdate.pop()
				break
		while nodesToUpdate:
			currentNode=nodesToUpdate.pop()
			pNode=up[currentNode]
			#print("Updating node "+str(currentNode))
			if currentNode==children[pNode][0]:
				siblingVector=probVectTime[children[pNode][1]]
				siblingDist=dist[children[pNode][1]]
			else:
				siblingVector=probVectTime[children[pNode][0]]
				siblingDist=dist[children[pNode][0]]
			if pNode==children[up[pNode]][0]:
				upVect=probVectUpRightTime[up[pNode]]
			else:
				upVect=probVectUpLeftTime[up[pNode]]
			newVectUpTime=mergeVectorsTime(upVect,dist[pNode],siblingVector,siblingDist,mutRate,isUpDown=True,tryMin=newMin)
			#print("upVect:")
			#print(upVect)
			#print("siblingVector:")
			#print(siblingVector)
			#print("newVectUpTime")
			#print(newVectUpTime)
			if currentNode==children[pNode][0]:
				probVectUpRightTime[pNode]=newVectUpTime
			else:
				probVectUpLeftTime[pNode]=newVectUpTime
		return
	

	#TODO merge two time likelihood vectors to create a new one 
	# if isUpDown, then probVect1 is assumed to be an upper vector and probVect2 a lower vector.
	#(and also calculate the logLk of the merging if necessary).
	def mergeVectorsTime(probVect1,bLen1,probVect2,bLen2,mutRate,returnLK=False,isUpDown=False,tryMin=None,verbose=False):
		bLen1*=lRef
		bLen2*=lRef
		if probVect1==None:
			if probVect2==None:
				if returnLK:
					return None,0.0
				else:
					return None
			passedVect2=getPartialVecTime(probVect2, mutRate, bLen2, otherMinT=None, otherMaxT=None, upNode=False, tryMin=tryMin)
			probVect=passedVect2[2]
			totSum=sum(probVect)
			for pos in range(len(probVect)):
				probVect[pos]=probVect[pos]/totSum
			if returnLK:
				return (passedVect2[0],passedVect2[1],probVect),log(totSum)
			else:
				return (passedVect2[0],passedVect2[1],probVect)
		elif probVect2==None:
			passedVect1=getPartialVecTime(probVect1, mutRate, bLen1, otherMinT=None, otherMaxT=None, upNode=isUpDown, tryMin=tryMin)
			probVect=passedVect1[2]
			totSum=sum(probVect)
			for pos in range(len(probVect)):
				probVect[pos]=probVect[pos]/totSum
			if returnLK:
				return (passedVect1[0],passedVect1[1],probVect),log(totSum)
			else:
				return (passedVect1[0],passedVect1[1],probVect)
		
		if bLen2>=1:
			maxT2=probVect2[0]-1
		else:
			maxT2=probVect2[0]
		if verbose:
			print("maxT2 "+str(maxT2))
		if isUpDown:
			if len(probVect1)==1:
				minT1=probVect1[0]
			else:
				minT1=probVect1[1]
			if bLen1>=1:
				minT1+=1
			#print("\n mergeVectorsTime")
			#print(probVect2)
			#print(probVect1)
			#print("")
			passedVect2=getPartialVecTime(probVect2, mutRate, bLen2, otherMinT=minT1, otherMaxT=None, upNode=False, tryMin=tryMin)
			if isinstance(passedVect2, int):
				print("mergeVectorsTime returning "+str(passedVect2))
				if returnLK:
					return float("-inf")
				else:
					return passedVect2
			passedVect1=getPartialVecTime(probVect1, mutRate, bLen1, otherMinT=None, otherMaxT=maxT2, upNode=True, tryMin=tryMin)
			#print(passedVect2)
			#print(passedVect1)
			#print("")
		else:
			if bLen1>=1:
				maxT1=probVect1[0]-1
			else:
				maxT1=probVect1[0]
			passedVect2=getPartialVecTime(probVect2, mutRate, bLen2, otherMinT=None, otherMaxT=maxT1, upNode=False, tryMin=tryMin)
			passedVect1=getPartialVecTime(probVect1, mutRate, bLen1, otherMinT=None, otherMaxT=maxT2, upNode=False, tryMin=tryMin)
		if verbose:
			print("initial vectors:")
			print(probVect1)
			print(probVect2)
			print("passedVects:")
			print(passedVect1)
			print(passedVect2)
		minT1=passedVect1[1]
		maxT1=passedVect1[0]
		minT2=passedVect2[1]
		maxT2=passedVect2[0]
		minT=max(minT1,minT2)
		maxT=min(maxT1,maxT2)
		#print(mutRate)
		#print(passedVect1)
		#print(passedVect2)
		#print(passedVect1)
		#print(passedVect2)
		#print(maxT)
		#print(minT)
		probVect=[]
		for pos in range(maxT,minT-1,-1):
			probVect.append(passedVect1[2][maxT1-pos]*passedVect2[2][maxT2-pos])
		maxValue=max(probVect)
		if tryMin==None:
			while (probVect[-1]<maxValue*timeProbThreshold):
				probVect.pop()
				minT+=1
		if probVect[0]<maxValue*timeProbThreshold:
			newProbVect=[]
			toReduce=True
			for i in range(len(probVect)):
				if toReduce and (probVect[i]<maxValue*timeProbThreshold):
					maxT-=1
				else:
					toReduce=False
					newProbVect.append(probVect[i])
			probVect=newProbVect
	
		totSum=sum(probVect)
		for pos in range(len(probVect)):
			probVect[pos]=probVect[pos]/totSum
		if returnLK:
			return (maxT,minT,probVect), log(totSum)
		else:
			return (maxT,minT,probVect)


	#TODO calculate the probability that results from combining a lower time likelihood vector at the root with root frequencies.
	# assumes root frequencies are very broad and uniform (uninformative flat prior).
	# The root probabilities should therefore be equal to 0<delta<<1 , but in reality the contribution of delta to the total likelihood is ignored for simplicity.
	# This is why this function is useless - but leaving here for future reference, in case one might want to use other root time priors.
	def findProbRootTime(probVect):
		return 0.0
	

	#TODO for the root, take lower time likelihoods, and create an overall likelihood (or upper right or upper left) by multiplying likelihoods by root frequencies.
	# As before, for now we assume an improper flat prior on root frequencies, hence the simplicity of this function.
	def rootVectorTime(probVect,bLen,mutRate,tryMin=None,returnLK=False):
		bLen*=lRef
		if probVect==None:
			return None
		passedVect=getPartialVecTime(probVect, mutRate, bLen, upNode=False,tryMin=tryMin)
		probVect2=passedVect[2]
		maxT=passedVect[0]
		maxValue=max(probVect2)
		if probVect2[0]<maxValue*timeProbThreshold:
			newProbVect=[]
			toReduce=True
			for i in range(len(probVect2)):
				if toReduce and (probVect2[i]<maxValue*timeProbThreshold):
					maxT-=1
				else:
					toReduce=False
					newProbVect.append(probVect2[i])
			probVect2=newProbVect

		totSum=sum(probVect2)
		for pos in range(len(probVect2)):
			probVect2[pos]=probVect2[pos]/totSum
		#print("rootVectorTime")
		#print(probVect)
		#print(bLen)
		#print(tryMin)
		#print((maxT,passedVect[1],probVect2))
		if returnLK:
			return (maxT,passedVect[1],probVect2),log(totSum)
		else:
			return (maxT,passedVect[1],probVect2)


	timeProbThreshold2=timeProbThreshold*timeProbThreshold
	#timeProbThreshold2=timeProbThreshold
	#TODO Check if two time likelihood vectors are effectively the same.
	#This is used to traverse the tree and only update vectors if they changed after a change in the tree.
	def areVectorsDifferentTime(probVect1,probVect2):
		if probVect2==None :
			if probVect1==None :
				return False
			else:
				return True
		elif probVect1==None:
			return True
		if len(probVect1)!=len(probVect2):
			return True
		if len(probVect1)==1:
			if probVect1[0]!=probVect2[0]:
				return True
			else:
				return False
		elif len(probVect1)==2:
			if probVect1[0]!=probVect2[0] or probVect1[1]!=probVect2[1]:
				return True
			else:
				return False
		for i in range(max(probVect1[0],probVect2[0]),min(probVect1[1],probVect2[1])-1,-1):
			if i>=probVect1[1] and i<=probVect1[0]:
				value1=probVect1[2][probVect1[0]-i]
			else:
				value1=None
			if i>=probVect2[1] and i<=probVect2[0]:
				value2=probVect2[2][probVect2[0]-i]
			else:
				value2=None
			if value1==None:
				if value2>=timeProbThreshold2:
					return True
			elif value2==None:
				if value1>=timeProbThreshold2:
					return True
			else:
				if abs(value1 - value2)>=timeProbThreshold2:
					return True
		return False

	if __name__ == "__main__":
		#TODO function to check is one date is less informative than another;
		#returns 0 when the 2 dates are not comparable, otherwise returns 1 if the first is more informative or if they are identical, and 2 otherwise.
		#if onlyFindIdentical, returns 1 if dates are identical, and 0 otherwise.
		def isMinorDate(date1,date2,onlyFindIdentical=False):
			found1bigger=False
			found2bigger=False
			if date2==None:
				if onlyFindIdentical:
					if date1==None:
						return 1
					else:
						return 0
				else:
					return 1
			elif date1==None:
				if onlyFindIdentical:
					return 0
				else:
					return 2
			max1=date1[0]
			max2=date2[0]
			if len(date1)>1:
				min1=date1[1]
			else:
				min1=max1
			if len(date2)>1:
				min2=date2[1]
			else:
				min2=max2
			if min1!=min2:
				if onlyFindIdentical:
					return 0
				if min1<min2:
					found2bigger=True
				else:
					found1bigger=True
			if max1!=max2:
				if onlyFindIdentical:
					return 0
				if max1>max2:
					found2bigger=True
				else:
					found1bigger=True

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
				

		# Check if two time likelihoods are the same or not.
		#This is a less strict version used for debugging purposes.
		# def areVectorsDifferentDebuggingTime(probVect1,probVect2,threshold=timeProbThreshold):
		# 	if len(probVect1)!=len(probVect2):
		# 		print("Different length")
		# 		return True
		# 	if probVect1[0]!=probVect2[0]:
		# 		print("Different max")
		# 		return True
		# 	if len(probVect1)>1 and probVect1[1]!=probVect2[1]:
		# 		print("Different min")
		# 		return True
		# 	if len(probVect1)>2:
		# 		for i in range(len(probVect1[2])):
		# 			if abs(probVect1[2][i]-probVect2[2][i])>threshold:
		# 				print("Different prob")
		# 				return True
		# 	return False
		

		#TODO Sort samples based on their latest possible sampling dates - those sampled later will be added later to the initial tree.
		def sortSamplesByDate(dates,data,samples=None,samplesInInitialTree=set(),forgetData=False):
			latestDates=[]
			if samples==None:
				rangeInd=range(len(data))
			else:
				rangeInd=samples
			for diffIndex in rangeInd:
				if (samples==None) or (not (diffIndex in samplesInInitialTree)):
					if (diffIndex in dates) and (dates[diffIndex]!=None):
						latestDates.append((dates[diffIndex][0],diffIndex))
					else:
						latestDates.append((float("inf"),diffIndex))
				elif forgetData:
					data[diffIndex]=None

			print("Now doing sorting")
			latestDates.sort(reverse=True,key=itemgetter(0))
			return latestDates
		

	# TODO function to calculate time likelihood cost of appending node to parent node 
	def appendProbNodeTime(probVectP,probVectC,mutRate,totLen,verbose=False):
		totLen*=lRef
		if probVectP==None or probVectC==None:
			return 0.0
		if len(probVectC)==1:
			minTC=probVectC[0]
		else:
			minTC=probVectC[1]
		minTP=probVectP[1]
		maxTP=probVectP[0]
		probV=probVectP[2]
		maxTC=probVectC[0]

		if abs(round(totLen)-totLen)>timeProbThreshold:
			totLens=[floor(totLen),ceil(totLen)]
			totLenProbs=[ceil(totLen)-totLen,totLen-floor(totLen)]
			if totLens[0]>0:
				newMinT=max(minTP+1,minTC)
			else:
				newMinT=max(minTP,minTC)
		else:
			totLens=None
			totLen=round(totLen)
			if totLen:
				newMinT=max(minTP+1,minTC)
			else:
				newMinT=max(minTP,minTC)
		if maxTC<newMinT:
			print("appendProbNodeTime() incompatible time vectors.")
			return float("-inf")
		if verbose:
			print("totLens:")
			print(totLens)
			print(totLen)
			print("Vectors:")
			print(probVectP)
			print(probVectC)
		currentT=newMinT
		totSum=0.0
		if totLens:
			if len(probVectC)>2:
				while currentT<=maxTC:
					lastProb=0.0
					for timeUp in range(minTP,min(currentT,maxTP)+1):
						lastProb+=( getPoissonCoeff(totLens[0],currentT-timeUp,mutRate)*totLenProbs[0] + getPoissonCoeff(totLens[1],currentT-timeUp,mutRate)*totLenProbs[1] ) * probV[maxTP-timeUp]
					totSum+=probVectC[2][maxTC-currentT]*lastProb
					currentT+=1
			else:
				while currentT<=maxTC:
					for timeUp in range(minTP,min(currentT,maxTP)+1):
						totSum+=( getPoissonCoeff(totLens[0],currentT-timeUp,mutRate)*totLenProbs[0] + getPoissonCoeff(totLens[1],currentT-timeUp,mutRate)*totLenProbs[1] ) * probV[maxTP-timeUp]
					currentT+=1
		else:
			if len(probVectC)>2:
				while currentT<=maxTC:
					lastProb=0.0
					for timeUp in range(minTP,min(currentT,maxTP)+1):
						lastProb+=getPoissonCoeff(totLen,currentT-timeUp,mutRate)*probV[maxTP-timeUp]
					totSum+=probVectC[2][maxTC-currentT]*lastProb
					currentT+=1
			else:
				while currentT<=maxTC:
					for timeUp in range(minTP,min(currentT,maxTP)+1):
						totSum+=getPoissonCoeff(totLen,currentT-timeUp,mutRate)*probV[maxTP-timeUp]
					currentT+=1

		return log(totSum)


	if __name__ == "__main__":

		# TODO adding more minor sequences can sligtly affect the time likelihood - this function updates the time likelihood as more identical samples are placed
		def updateProbVectTerminalNodeTime(tree,node,sampleTimeLK,numMinSeqs,mutRate,onlyAddOne=False):
			if onlyAddOne:
				tree.probVectTime[node]=mergeVectorsTime(tree.probVectTime[node],0.0,sampleTimeLK,0.0,mutRate,returnLK=False,isUpDown=False)
			else:
				if sampleTimeLK is None:
					tree.probVectTime[node]=None
				else:
					newProbVect=sampleTimeLK
					for i in range(numMinSeqs):
						newProbVect=mergeVectorsTime(newProbVect,0.0,sampleTimeLK,0.0,mutRate,returnLK=False,isUpDown=False)
					tree.probVectTime[node]=newProbVect
				#print("\n updateProbVectTerminalNodeTime "+str(node))
				#print(sampleTimeLK)
				#print(tree.probVectTime[node])


		# TODO Given a tree and its time vectors, calculate mutation counts and waiting times for expectation maximization estimation of the mutation rate
		def expectationMaximizationCalculationRatesTime(tree,root,mutRate):
			up=tree.up
			children=tree.children
			probVectUpRightTime=tree.probVectUpRightTime
			probVectUpLeftTime=tree.probVectUpLeftTime
			probVectTime=tree.probVectTime
			dist=tree.dist
			node=root
			#direction 0 means from parent, direction 1 means from a child
			lastNode=None
			direction=0
			#Initializing with pseudo-counts based on initial value of mutation rate
			waitingTimes=20.0
			counts=args.mutRate*intervalLength*waitingTimes
			while node!=None:
				if direction==0:
					if up[node]!=None:
						#update counts and waiting times
						if node==children[up[node]][0]:
							probVectP=probVectUpRightTime[up[node]]
						else:
							probVectP=probVectUpLeftTime[up[node]]
						probVectC=probVectTime[node]

						if probVectP!=None and probVectC!=None:
							totLen=dist[node]*lRef
							if len(probVectC)==1:
								minTC=probVectC[0]
							else:
								minTC=probVectC[1]
							minTP=probVectP[1]
							maxTP=probVectP[0]
							probV=probVectP[2]
							maxTC=probVectC[0]

							if abs(round(totLen)-totLen)>timeProbThreshold:
								totLens=[floor(totLen),ceil(totLen)]
								totLenProbs=[ceil(totLen)-totLen,totLen-floor(totLen)]
								if totLens[0]>0:
									newMinT=max(minTP+1,minTC)
								else:
									newMinT=max(minTP,minTC)
							else:
								totLens=None
								totLen=round(totLen)
								if totLen:
									newMinT=max(minTP+1,minTC)
								else:
									newMinT=max(minTP,minTC)
							if maxTC<newMinT:
								print("expectationMaximizationCalculationRatesTime() incompatible time vectors.")
								raise Exception("exit")
							#first calculate the normalization factor
							totSum=0.0
							currentT=newMinT
							if totLens:
								if len(probVectC)>2:
									while currentT<=maxTC:
										lastProb=0.0
										for timeUp in range(minTP,min(currentT,maxTP)+1):
											lastProb+=( getPoissonCoeff(totLens[0],currentT-timeUp,mutRate)*totLenProbs[0] + getPoissonCoeff(totLens[1],currentT-timeUp,mutRate)*totLenProbs[1] ) * probV[maxTP-timeUp]
										totSum+=probVectC[2][maxTC-currentT]*lastProb
										currentT+=1
								else:
									while currentT<=maxTC:
										for timeUp in range(minTP,min(currentT,maxTP)+1):
											totSum+=( getPoissonCoeff(totLens[0],currentT-timeUp,mutRate)*totLenProbs[0] + getPoissonCoeff(totLens[1],currentT-timeUp,mutRate)*totLenProbs[1] ) * probV[maxTP-timeUp]
										currentT+=1
							else:
								if len(probVectC)>2:
									while currentT<=maxTC:
										lastProb=0.0
										for timeUp in range(minTP,min(currentT,maxTP)+1):
											lastProb+=getPoissonCoeff(totLen,currentT-timeUp,mutRate)*probV[maxTP-timeUp]
										totSum+=probVectC[2][maxTC-currentT]*lastProb
										currentT+=1
								else:
									while currentT<=maxTC:
										for timeUp in range(minTP,min(currentT,maxTP)+1):
											totSum+=getPoissonCoeff(totLen,currentT-timeUp,mutRate)*probV[maxTP-timeUp]
										currentT+=1
							#now add the contributions to mutation counts and waiting times
							currentT=newMinT
							if totLens:
								if len(probVectC)>2:
									while currentT<=maxTC:
										for timeUp in range(minTP,min(currentT,maxTP)+1):
											lastProb= getPoissonCoeff(totLens[0],currentT-timeUp,mutRate)*totLenProbs[0] * probV[maxTP-timeUp] * probVectC[2][maxTC-currentT] / totSum
											waitingTimes+=lastProb*(currentT-timeUp)
											counts+=lastProb*totLens[0]
											lastProb= getPoissonCoeff(totLens[1],currentT-timeUp,mutRate)*totLenProbs[1] * probV[maxTP-timeUp] * probVectC[2][maxTC-currentT] / totSum
											waitingTimes+=lastProb*(currentT-timeUp)
											counts+=lastProb*totLens[1]
										currentT+=1
								else:
									while currentT<=maxTC:
										for timeUp in range(minTP,min(currentT,maxTP)+1):
											lastProb=getPoissonCoeff(totLens[0],currentT-timeUp,mutRate)*totLenProbs[0] * probV[maxTP-timeUp] / totSum
											waitingTimes+=lastProb*(currentT-timeUp)
											counts+=lastProb*totLens[0]
											lastProb=getPoissonCoeff(totLens[1],currentT-timeUp,mutRate)*totLenProbs[1] * probV[maxTP-timeUp] / totSum
											waitingTimes+=lastProb*(currentT-timeUp)
											counts+=lastProb*totLens[1]
										currentT+=1
							else:
								if len(probVectC)>2:
									while currentT<=maxTC:
										for timeUp in range(minTP,min(currentT,maxTP)+1):
											lastProb=getPoissonCoeff(totLen,currentT-timeUp,mutRate) * probV[maxTP-timeUp] * probVectC[2][maxTC-currentT] / totSum
											waitingTimes+=lastProb*(currentT-timeUp)
											counts+=lastProb*totLen
										currentT+=1
								else:
									while currentT<=maxTC:
										for timeUp in range(minTP,min(currentT,maxTP)+1):
											lastProb=getPoissonCoeff(totLen,currentT-timeUp,mutRate)*probV[maxTP-timeUp] / totSum
											waitingTimes+=lastProb*(currentT-timeUp)
											counts+=lastProb*totLen
										currentT+=1

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
			if counts/waitingTimes<minMutRate:
				print("WARNING Mutation rate estimate reached the input minimum allowed "+str(minMutRate))
				return counts, waitingTimes, minMutRate
			return counts, waitingTimes, counts/waitingTimes


		# TODO Given a tree, and a mutation rate, calculate the time likelihood of the tree
		def calculateTreeLikelihoodTime(tree,root,mutRate,checkCorrectness=False):
			up=tree.up
			children=tree.children
			probVectTime=tree.probVectTime
			dist=tree.dist
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
						probVect0=probVectTime[children[node][0]]
						probVect1=probVectTime[children[node][1]]
						newLowerTime, LkcontributionTime=mergeVectorsTime(probVect0,dist[children[node][0]],probVect1,dist[children[node][1]],mutRate,returnLK=True)
						totalLK+=LkcontributionTime
						#print(probVect0)
						#print(dist[children[node][0]])
						#print(probVect1)
						#print(dist[children[node][1]])
						#print("TotalLk: "+str(totalLK))
						#print(" ")
						#if newLowerTime==None:
						#	print("Strange, inconsistent lower time LK list creation in calculateTreeLikelihoodTime().")
						#	raise Exception("exit")
						if checkCorrectness and areVectorsDifferentTime(probVectTime[node],newLowerTime): #areVectorsDifferentDebuggingTime(probVectTime[node],newLowerTime):
							print("Strange, while calculating tree time likelihood encountered non-updated lower time likelihood list at node "+str(node))
							raise Exception("exit")
						lastNode=node
						node=up[node]
						direction=1
			#now add contribution from the root
			#print("Before root contribution: "+str(totalLK))
			totalLK+=findProbRootTime(probVectTime[root])
			#print("After root contribution: "+str(totalLK))
			return totalLK
		

		# TODO we need at least 15 very informative tips/dates to defeat the impact of the improper root prior
		# TODO which pushes down the inferred mutation rate. 
		# TODO below I test this using a caterpillar tree with perfectly informative branch lengths and times.
		# lRef=10000
		# bLen=0.0
		# bLen2=0.0001
		# mutRates=[0.1,1.0]
		# for mutRate in mutRates:
		# 	poissonCoeff=[[1.0]]
		# 	newVect=None
		# 	totalLk=0.0
		# 	for i in range(15):
		# 		#print(i)
		# 		#print(newVect)
		# 		#print(totalLk)
		# 		newVect,Lkcost=mergeVectorsTime(newVect,0.0001,(10-i,),0.0,mutRate,returnLK=True)
		# 		totalLk+=Lkcost
		# 	print(mutRate)
		# 	print(totalLk)
		#raise Exception("exit")

		#used to make sure that likelihoods calculated from different nodes are consistent
		def calculateTreeLikelihoodTimeDebugging(tree,root,mutRate,checkCorrectness=False):
			up=tree.up
			children=tree.children
			probVectTime=tree.probVectTime
			probVectUpRightTime=tree.probVectUpRightTime
			probVectUpLeftTime=tree.probVectUpLeftTime
			dist=tree.dist
			node=root
			#direction 0 means from parent, direction 1 means from a child
			lastNode=None
			direction=0
			totalLK=0.0
			totalLKdown=[0.0]*len(up)
			totalLKupRight=[0.0]*len(up)
			totalLKupLeft=[0.0]*len(up)
			totalLKnode=[0.0]*len(up)
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
						probVect0=probVectTime[children[node][0]]
						probVect1=probVectTime[children[node][1]]
						newLowerTime, LkcontributionTime=mergeVectorsTime(probVect0,dist[children[node][0]],probVect1,dist[children[node][1]],mutRate,returnLK=True)
						totalLKdown[node]=totalLKdown[children[node][0]]+totalLKdown[children[node][1]]+LkcontributionTime
						totalLK+=LkcontributionTime
						if newLowerTime==None:
							print("Strange, inconsistent lower time LK list creation in calculateTreeLikelihoodTime().")
							raise Exception("exit")
						elif checkCorrectness and areVectorsDifferentTime(probVectTime[node],newLowerTime): #areVectorsDifferentDebuggingTime(probVectTime[node],newLowerTime):
							print("Strange, while calculating tree time likelihood encountered non-updated lower time likelihood list at node "+str(node))
							raise Exception("exit")
						lastNode=node
						node=up[node]
						direction=1
			totalLKnode[root]=totalLKdown[root]
			print("TotalLKdebugging root node "+str(root)+" LK "+str(totalLKnode[root]))

			#now update the other LKs for the root
			node=root
			if children[node]:
				probVect1=probVectTime[children[node][1]]
				#print("Recalcualting upRightTime for node "+str(node))
				newVect,LK=rootVectorTime(probVect1,dist[children[node][1]],mutRate,returnLK=True)
				totalLKupRight[node]=totalLKdown[children[node][1]]+LK
				probVect0=probVectTime[children[node][0]]
				#print("Recalcualting upLeftTime for node "+str(node))
				newVect,LK=rootVectorTime(probVect0,dist[children[node][0]],mutRate,returnLK=True)
				totalLKupLeft[node]=totalLKdown[children[node][0]]+LK
				
				#now traverse the tree downward and update the non-lower genome lists for all other nodes of the tree.
				lastNode=None
				node=children[node][0]
				direction=0
				while node!=None:
					if direction==0:
						if node==children[up[node]][0]:
							vectUpTime=probVectUpRightTime[up[node]]
							LK=totalLKupRight[up[node]]
						else:
							vectUpTime=probVectUpLeftTime[up[node]]
							LK=totalLKupLeft[up[node]]
						#print("merging")
						#print(vectUpTime)
						#print(probVectTime[node])
						#print(node)
						#print(dist[node])
						#if node==6:
						#	verbose=True
						#else:
						verbose=False
						LkcontributionTime=appendProbNodeTime(vectUpTime,probVectTime[node],mutRate,dist[node],verbose=verbose)
						totalLKnode[node]=totalLKdown[node]+LK+LkcontributionTime
						print("TotalLKdebugging node "+str(node)+" totalLK "+str(totalLKnode[node]))
						#newVect,LkcontributionTime2=mergeVectorsTime(vectUpTime,dist[node]/2,probVectTime[node],dist[node]/2,mutRate,isUpDown=True,returnLK=True,verbose=verbose)
						#newVect,LkcontributionTime3=mergeVectorsTime(vectUpTime,0.0,probVectTime[node],dist[node],mutRate,isUpDown=True,returnLK=True)
						#newVect,LkcontributionTime4=mergeVectorsTime(vectUpTime,dist[node],probVectTime[node],0.0,mutRate,isUpDown=True,returnLK=True)
						#print(str(LkcontributionTime)+" "+str(LkcontributionTime2)+" "+str(LkcontributionTime3)+" "+str(LkcontributionTime4)+" ")

						if children[node]:
							probVect0=probVectTime[children[node][0]]
							probVect1=probVectTime[children[node][1]]
							newUpRight,LkcontributionTime=mergeVectorsTime(vectUpTime,dist[node],probVect1,dist[children[node][1]],mutRate,isUpDown=True,returnLK=True)
							totalLKupRight[node]=totalLKdown[children[node][1]]+LK+LkcontributionTime
							newUpLeft,LkcontributionTime=mergeVectorsTime(vectUpTime,dist[node],probVect0,dist[children[node][0]],mutRate,isUpDown=True,returnLK=True)
							totalLKupLeft[node]=totalLKdown[children[node][0]]+LK+LkcontributionTime
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

			#now add contribution from the root
			totalLK+=findProbRootTime(probVectTime[root])
			return totalLK
		

		#TODO Given a tree, and a mutation rate, re-calculate all time LK lists within the tree.
		def reCalculateAllGenomeListsTime(tree,root,mutRate, checkExistingAreCorrect=False): #checkSamplesIntree=False
			#print("reCalculateAllGenomeListsTime")
			up=tree.up
			children=tree.children
			minorSequences=tree.minorSequences
			dist=tree.dist
			probVectTime=tree.probVectTime
			probVectTotUpTime=tree.probVectTotUpTime
			probVectUpRightTime=tree.probVectUpRightTime
			probVectUpLeftTime=tree.probVectUpLeftTime
			dateData=tree.dateData
			#first pass to update all lower likelihoods.
			node=root
			#direction 0 means from parent, direction 1 means from a child
			lastNode=None
			direction=0
			while node!=None:
				if direction==0:
					if children[node]:
						node=children[node][0]
					else:
						updateProbVectTerminalNodeTime(tree,node,dateData[node],len(minorSequences[node]),mutRate,onlyAddOne=False)
						lastNode=node
						node=up[node]
						direction=1
				else :
					if lastNode==children[node][0]:
						node=children[node][1]
						direction=0
					else:
						probVect0=probVectTime[children[node][0]]
						probVect1=probVectTime[children[node][1]]
						newLower=mergeVectorsTime(probVect0,dist[children[node][0]],probVect1,dist[children[node][1]],mutRate)
						if checkExistingAreCorrect:
							#if areVectorsDifferentDebuggingTime(newLower,probVectTime[node]):
							if areVectorsDifferentTime(newLower,probVectTime[node]):
								print("Inside reCalculateAllGenomeLists(), new lower at node is different from the old one, and it shouldn't be.")
								raise Exception("exit")
						probVectTime[node]=newLower
						lastNode=node
						node=up[node]
						direction=1

			#now update the other LKs for the root
			node=root
			if children[node]:
				probVect1=probVectTime[children[node][1]]
				#print("Recalcualting upRightTime for node "+str(node))
				newVect=rootVectorTime(probVect1,dist[children[node][1]],mutRate)
				if checkExistingAreCorrect:
					if areVectorsDifferentTime(newVect,probVectUpRightTime[node]):
					#if areVectorsDifferentDebuggingTime(newVect,probVectUpRightTime[node]):
						print("new probVectUpRightTime at root is different from the old one, and it shouldn't be.")
						raise Exception("exit")
				probVectUpRightTime[node]=newVect
				probVect0=probVectTime[children[node][0]]
				#print("Recalcualting upLeftTime for node "+str(node))
				newVect=rootVectorTime(probVect0,dist[children[node][0]],mutRate)
				if checkExistingAreCorrect:
					if areVectorsDifferentTime(newVect,probVectUpLeftTime[node]):
					#if areVectorsDifferentDebuggingTime(newVect,probVectUpLeftTime[node]):
						print("new probVectUpLeftTime at root is different from the old one, and it shouldn't be; distance of child node "+str(dist[children[node][0]]))
						raise Exception("exit")
				probVectUpLeftTime[node]=newVect
				
				#now traverse the tree downward and update the non-lower genome lists for all other nodes of the tree.
				lastNode=None
				node=children[node][0]
				direction=0
				while node!=None:
					if direction==0:
						if node==children[up[node]][0]:
							vectUpTime=probVectUpRightTime[up[node]]
						else:
							vectUpTime=probVectUpLeftTime[up[node]]
						#print("merging")
						#print(vectUpTime)
						#print(probVectTime[node])
						#print(node)
						#print(dist[node])
						newVect,newVectProb=mergeVectorsTime(vectUpTime,dist[node]/2,probVectTime[node],dist[node]/2,mutRate,isUpDown=True,returnLK=True)
						if isinstance(newVect, int):
							resolveTimeInconsistency(tree,node,newVect,mutRate)
							if node==children[up[node]][0]:
								vectUpTime=probVectUpRightTime[up[node]]
							else:
								vectUpTime=probVectUpLeftTime[up[node]]
							newVect,newVectProb=mergeVectorsTime(vectUpTime,dist[node]/2,probVectTime[node],dist[node]/2,mutRate,isUpDown=True,returnLK=True)
						if checkExistingAreCorrect:
							#if areVectorsDifferentDebuggingTime(newVect,probVectTotUpTime[node]):
							if areVectorsDifferentTime(newVect,probVectTotUpTime[node][0]):
								print("new probVectTotUpTime at node is different from the old one, and it shouldn't be.")
								raise Exception("exit")
						newVectProb-=appendProbNodeTime(vectUpTime,probVectTime[node],mutRate,dist[node])
						probVectTotUpTime[node]=(newVect,newVectProb)
						if children[node]:
							probVect0=probVectTime[children[node][0]]
							probVect1=probVectTime[children[node][1]]
							newUpRight=mergeVectorsTime(vectUpTime,dist[node],probVect1,dist[children[node][1]],mutRate,isUpDown=True)
							if isinstance(newUpRight, int):
								resolveTimeInconsistency(tree,node,newUpRight,mutRate)
								if node==children[up[node]][0]:
									vectUpTime=probVectUpRightTime[up[node]]
								else:
									vectUpTime=probVectUpLeftTime[up[node]]
								newUpRight=mergeVectorsTime(vectUpTime,dist[node],probVect1,dist[children[node][1]],mutRate,isUpDown=True)
							if checkExistingAreCorrect:
								#if areVectorsDifferentDebuggingTime(newUpRight,probVectUpRightTime[node]):
								if areVectorsDifferentTime(newUpRight,probVectUpRightTime[node]):
									print("new probVectUpRightTime at node is different from the old one, and it shouldn't be.")
									print(vectUpTime)
									print(probVect1)
									print(dist[node])
									print(dist[children[node][1]])
									print(newUpRight)
									print(probVectUpRightTime[node])
									print(node)
									raise Exception("exit")
							probVectUpRightTime[node]=newUpRight
							newUpLeft=mergeVectorsTime(vectUpTime,dist[node],probVect0,dist[children[node][0]],mutRate,isUpDown=True)
							if isinstance(newUpLeft, int):
								resolveTimeInconsistency(tree,node,newUpLeft,mutRate)
								if node==children[up[node]][0]:
									vectUpTime=probVectUpRightTime[up[node]]
								else:
									vectUpTime=probVectUpLeftTime[up[node]]
								newUpLeft=mergeVectorsTime(vectUpTime,dist[node],probVect0,dist[children[node][0]],mutRate,isUpDown=True)
							if checkExistingAreCorrect:
								#if areVectorsDifferentDebuggingTime(newUpLeft,probVectUpLeftTime[node]):
								if areVectorsDifferentTime(newUpLeft,probVectUpLeftTime[node]):
									print("new probVectUpLeftTime at node is different from the old one, and it shouldn't be.")
									print(newUpLeft)
									print(probVectUpLeftTime[node])
									print(node)
									raise Exception("exit")
							probVectUpLeftTime[node]=newUpLeft

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
			#print("finished reCalculateAllGenomeListsTime")
			return
		

		# TODO process time data. At the end, dates is a dictionary of the date data of all samples, represented as interval numbers.
		#  minYear and minDay is the year and day number (in that year) for the earliest possible sampling time observed.
		# MAPLE represents this earliest sampling time internally as interval 0, but a final translation into real dates needs to be performed before creating the output files. 
		print("Reading input dates")
		dates = {}
		if not os.path.isfile(datesFile):
			print("\n\tERROR: file %s does not exist"%datesFile)
			raise Exception("exit")
		else:
			full_sep = ',' if datesFile.endswith('.csv') else '\t'

			file=open(datesFile)
			line=file.readline()
			line=line.replace("\n","")
			linelist=line.split(full_sep)
			print(line)
			print(linelist)
			if len(linelist)<2:
				print("ERROR: at least two columns are needed in the time metadata file: one for sample names, and one for dates.")
				print(line)
				print('attempted separator is '+full_sep)
				raise Exception("exit")
			# Find columns corresponding to sample names and times
			index_date=None
			index_name=None
			potential_index_date=None
			potential_index_name=None
			date_col=dateName.lower()
			name_col=strainName.lower()
			columns=[]
			for entry in linelist:
				columns.append(entry.lower())
			if date_col:
				for entry in range(len(columns)):
					if columns[entry]==date_col:
						index_date=entry
						break
					elif date_col in columns[entry]:
						potential_index_date=entry
				if index_date==None:
					if potential_index_date!=None:
						index_date=potential_index_date
				if index_date==None:
					print("ERROR: specified column for dates does not exist. Will look for other suitable columns. The first row of the time metadata file is \n "+line+" You specified '%s'"%date_col)
			if index_date==None:
				timeStrList=['date','time']
				for entry in range(len(columns)):
					if columns[entry] in timeStrList:
						index_date=entry
						break
				if index_date!=None:
					print("Chosen as date column "+linelist[index_date])
			if name_col:
				for entry in range(len(columns)):
					if columns[entry]==name_col:
						index_name=entry
						break
					elif name_col in columns[entry]:
						potential_index_name=entry
				if index_name==None:
					if potential_index_name!=None:
						index_name=potential_index_name
				if index_name==None:
					print("ERROR: specified column for sample names does not exist. Will look for other suitable columns. The first row of the time metadata file is \n "+line+" You specified '%s'"%name_col)
			if index_name==None:
				nameStrList=['name','strain','accession','id','sample','names','strains','ids','accessions','samples']
				for entry in range(len(columns)):
					if columns[entry] in nameStrList:
						index_name=entry
						break
				if index_name==None:
					for entry in range(len(columns)):
						for name in nameStrList:
							if name in columns[entry]:
								index_name=entry
								break	
				if index_name!=None:
					print("Chosen as name column "+linelist[index_name])
			if index_name==None:
				if index_date==None:
					print("Suitable column names not found in time metadata. I will assume that the first column contains sample names, the second dates, and that column names are missing")
					index_name=0
					index_date=1
				else:
					line=file.readline()
					line=line.replace("\n","")
					if index_date==0:
						print("Suitable sample names column not found in time metadata. I will assume that the first column contains dates, the second sample names.")
						index_name=1
					else:
						print("Suitable sample names column not found in time metadata. I will assume that the first column contains sample names.")
						index_name=0
			else:
				line=file.readline()
				line=line.replace("\n","")
				if index_date==None:
					if index_name==0:
						print("Suitable date column not found in time metadata. I will assume that the first column contains sample names, the second dates.")
						index_date=1
					else:
						print("Suitable date column not found in time metadata. I will assume that the first column contains dates.")
						index_date=0
			# Read time data and translate it into floats
			dates = {}
			minLength=1+max(index_date,index_name)
			minDate=float("inf")
			maxDate=0
			months={"01":1,"02":2,"03":3,"04":4,"05":5,"06":6,"07":7,"08":8,"09":9,"10":10,"11":11,"12":12
				,"1":1,"2":2,"3":3,"4":4,"5":5,"6":6,"7":7,"8":8,"9":9
				,"january":1,"february":2,"march":3,"april":4,"may":5,"june":6,"july":7,"august":8,"september":9,"october":10,"november":11,"december":12
				,"jan":1,"feb":2,"mar":3,"apr":4,"may":5,"jun":6,"jul":7,"aug":8,"sep":9,"oct":10,"nov":11,"dec":12}
			days={"01":1,"02":2,"03":3,"04":4,"05":5,"06":6,"07":7,"08":8,"09":9,"10":10,"11":11,"12":12,"13":13,"14":14,"15":15,"16":16,"17":17,"18":18,"19":9,"20":20
				,"21":21,"22":22,"23":23,"24":24,"25":25,"26":26,"27":27,"28":28,"29":29,"30":30,"31":31
				,"1":1,"2":2,"3":3,"4":4,"5":5,"6":6,"7":7,"8":8,"9":9}
			monthdays=[[0,31,59,90,120,151,181,212,243,273,304,334,365],[0,31,60,91,121,152,182,213,244,274,305,335,366]]
			while line!="" and line!="\n":
				linelist=line.split(full_sep)
				if len(linelist)<minLength:
					print("Found fewer entries in a date metadata file row than needed. Terminating reading the metadata file. Row: \n"+line)
					print(line)
					print(linelist)
					print(index_name)
					print(index_date)
					break
				date_str=linelist[index_date]
				name=linelist[index_name]
				if date_str and date_str!="." and date_str.lower()!="unknown"  and date_str.lower()!="not applicable" and date_str.lower()!="not provided"  and date_str.lower()!="not collected" and date_str.lower()!="missing":
					try:
						date=int(date_str)
						date=(float(date),float(date+1))
						if ((minSamplingYear!=None) and (int(date_str)<minSamplingYear)) or ((maxSamplingYear!=None) and (int(date_str)>maxSamplingYear)):
							print("Sampling date outside of sampling range, treated as an error and ignored:")
							print(line)
							date=None
					except:
						try:
							date = (float(date_str),)
							if ((minSamplingYear!=None) and (float(date_str)<minSamplingYear)) or ((maxSamplingYear!=None) and ((float(date_str)-1)>maxSamplingYear)):
								print("Sampling date outside of sampling range, treated as an error and ignored:")
								print(line)
								date=None
						except:
							if "-" in date_str:
								sep="-"
							elif " " in date_str:
								sep=" "
							elif "\t" in date_str:
								sep="\t"
							elif "\\" in date_str:
								sep="\\"
							elif "/" in date_str:
								sep="/"
							elif "." in date_str:
								sep="."
							datelist=date_str.split(sep)
							if len(datelist)>3:
								print("0 Unrecognized date format: "+date_str+" ignoring this date for sample "+name+" and treating it as unknown date.")
								date=None
							elif len(datelist)==2:
								year=None
								month=None
								try:
									year=int(datelist[0])
									if year>100:
										month=months[datelist[1].lower()]
									else:
										year=int(datelist[1])
										month=months[datelist[0].lower()]
									if isleap(year):
										date=(year+monthdays[1][month-1]/366.0,year+(monthdays[1][month]-1)/366.0)
									else:
										date=(year+monthdays[0][month-1]/365.0,year+(monthdays[0][month]-1)/365.0)

									if ((minSamplingYear!=None) and (year<minSamplingYear)) or ((maxSamplingYear!=None) and (year>maxSamplingYear)):
										print("Sampling date outside of sampling range, treated as an error and ignored:")
										print(line)
										date=None
									
								except:
									print("1 Unrecognized date format: "+date_str+" ignoring this date for sample "+name+" and treating it as unknown date.")
									date=None
							elif len(datelist)==3:
								year=None
								month=None
								day=None
								try:
									year=int(datelist[0])
									if year>100:
										month=months[datelist[1].lower()]
										day=days[datelist[2]]
									else:
										year=int(datelist[2])
										month=months[datelist[1].lower()]
										day=days[datelist[0]]
									if isleap(year):
										date=(year+(monthdays[1][month-1]+day-1)/366.0,)
									else:
										date=(year+(monthdays[0][month-1]+day-1)/365.0,)

									if ((minSamplingYear!=None) and (year<minSamplingYear)) or ((maxSamplingYear!=None) and (year>maxSamplingYear)):
										print("Sampling date outside of sampling range, treated as an error and ignored:")
										print(line)
										date=None

								except:
									print("2 Unrecognized date format: "+date_str+" ignoring this date for sample "+name+" and treating it as unknown date.")
									date=None
							else:
								print("3 Unrecognized date format: "+date_str+" ignoring this date for sample "+name+" and treating it as unknown date.")
								date=None
				else:
					date = None

				if date:
					if len(date)>1:
						if date[0]<minDate:
							minDate=date[0]
						if date[1]>maxDate:
							maxDate=date[1]
					else:
						if date[0]<minDate:
							minDate=date[0]
						if date[0]>maxDate:
							maxDate=date[0]
				dates[name]=date
				line=file.readline()
				line=line.replace("\n","")
			minYear=floor(minDate)
			if isleap(minYear):
				minDay=round((minDate-minYear)*366)
			else:
				minDay=round((minDate-minYear)*365)
			carryOverDays={minYear:0}
			maxYear=floor(maxDate)
			carryOver=0
			for i in range(maxYear-minYear):
				if isleap(minYear+i):
					carryOver+=366
				else:
					carryOver+=365
				carryOverDays[minYear+i+1]=carryOver
			for name in dates:
				if dates[name]:
					if len(dates[name])==2:
						newDates=[]
						for i in range(2):
							date=dates[name][i]
							year=floor(date)
							if isleap(minYear):
								days=round((date-year)*366)
							else:
								days=round((date-year)*365)
							days+=carryOverDays[year]
							days-=minDay
							days=floor(days/intervalLength)
							newDates.append(days)
						dates[name]=(newDates[1],newDates[0])
					else:
						date=dates[name][0]
						year=floor(date)
						if isleap(minYear):
							days=round((date-year)*366)
						else:
							days=round((date-year)*365)
						days+=carryOverDays[year]
						days-=minDay
						dates[name]=(floor(days/intervalLength),)

			print("Finished reading time metadata. Earliest day is "+str(minDay)+" of year "+str(minYear))







if __name__ == "__main__":
	#function to read input newick string
	def readNewick(nwFile,multipleTrees=False,dirtiness=True,createDict=False,inputDictNames=None,keepNames=False,onlyTerminalNodeName=False):
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
			isInternalName=False
			while index<len(nwString):
				if nwString[index]=="(":
					tree.children[nodeIndex].append(len(tree.up))
					tree.addNode(dirtiness=dirtiness)
					if keepInputIQtreeSupports:
						tree.IQsupport.append(None)
					tree.up[-1]=nodeIndex
					nodeIndex=len(tree.up)-1
					index+=1
					isInternalName=False
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
							if (not onlyTerminalNodeName) or (not isInternalName):
								tree.name[nodeIndex]=sampleNum
								if createDict:
									namesInTreeDict[name]=sampleNum
								sampleNum+=1
								namesInTree.append(name)
						else:
							if (not onlyTerminalNodeName) or (not isInternalName):
								name=name.replace("?","_").replace("&","_")
								try:
									tree.name[nodeIndex]=inputDictNames[name]
								except:
									print("Error: sample "+name+" not found in the original tree. Exiting.")
									raise Exception("exit")
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
					isInternalName=False
				elif nwString[index]==")":
					if name!="":
						if keepNames:
							tree.name[nodeIndex]=name
						elif inputDictNames==None:
							if (not onlyTerminalNodeName) or (not isInternalName):
								tree.name[nodeIndex]=sampleNum
								if createDict:
									namesInTreeDict[name]=sampleNum
								sampleNum+=1
								namesInTree.append(name)
						else:
							if (not onlyTerminalNodeName) or (not isInternalName):
								name=name.replace("?","_").replace("&","_")
								try:
									tree.name[nodeIndex]=inputDictNames[name]
								except:
									print("Error: sample "+name+" not found in the original tree. Exiting.")
									raise Exception("exit")
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
					isInternalName=True
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
	def reRootTree(tree,root,sample,reRootAtInternalNode=False):
		sampleNode=None
		up=tree.up
		children=tree.children
		dist=tree.dist
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
			#update nDesc0 of nodes afected by the re-rooting
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
	def prepareTreeComparison(tree,t1,namesInTree,namesInTreeDict,rooted=False,minimumBLen=0.000006):
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
		#leafNameDictReverse=[]
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
					newname=(namesInTree[name[node]]).replace("?","_").replace("&","_")
					if newname!=namesInTree[name[node]]:
						namesInTreeDict[newname]=namesInTreeDict[namesInTree[name[node]]]
						namesInTree[name[node]]=newname
					#name[node]=(name[node]).replace("?","_").replace("&","_")
					leafNameDict[name[node]]=leafCount
					#leafNameDictReverse.append(name[node])
					if rooted:
						nodeTable.append([0,0])
					lastL=leafCount
					lastR=leafCount
					lastDesc=1
					leafCount+=1
					nextNode=up[node]
					movingFrom=1
					#leafDistDict[namesInTree[name[node]]] = dist[node]
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
					#name[node]=(name[node]).replace("?","_").replace("&","_")
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

if __name__ == "__main__":
	#TODO write also inferred dates or date intervals in case of flag doTimeTree
	#TODO
	#TODO
	#generate the string corresponding to a node, taking into account support measures and other possible node features.
	def stringForNode(tree,nextNode,nameNode,distB,estimateMAT=estimateMAT,networkOutput=networkOutput,aBayesPlusOn=aBayesPlusOn,namesInTree=None, performLineageAssignmentByRefPlacement = performLineageAssignmentByRefPlacement):
		children=tree.children
		up=tree.up
		name=tree.name
		aBayesPlusActive=False
		writeLineageAssignment = performLineageAssignment or performLineageAssignmentByRefPlacement
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
		if writeLineageAssignment:
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
				if aBayesPlusActive and rootSupport!=None and rootSupport[nextNode]!=None:
					strings.append("rootSupport="+str(rootSupport[nextNode]))
				if aBayesPlusActive and ( (distB>effectivelyNon0BLen) or supportFor0Branches) and support[nextNode]!=None:
					strings.append("support="+str(support[nextNode]))
					if networkOutput and alternativePlacements[nextNode]:
						newString="alternativePlacements={"
						for iNode in range(len(alternativePlacements[nextNode])):
							newString+=namesInTree[name[alternativePlacements[nextNode][iNode][0]]]+":"+str(alternativePlacements[nextNode][iNode][1])
							if iNode<(len(alternativePlacements[nextNode])-1):
								newString+=","
						newString+="}"
						strings.append(newString)
				if estimateMATon and ( (distB) or errorsOn or (not children[nextNode])):
					if mutationsInf[nextNode]:
						newString="mutationsInf={"
						for iNode in range(len(mutationsInf[nextNode])):
							mutation=mutationsInf[nextNode][iNode]
							newString+=allelesListExt[mutation[0]]+str(mutation[1])+allelesListExt[mutation[2]]+":"+str(mutation[3])
							if iNode<(len(mutationsInf[nextNode])-1):
								newString+=","
						newString+="}"
						strings.append(newString)
					if Ns[nextNode]:
						newString="Ns={"
						for iNode in range(len(Ns[nextNode])):
							mutation=Ns[nextNode][iNode]
							if type(mutation)==int:
								newString+=str(mutation)
							else:
								newString+=str(mutation[0])+"-"+str(mutation[1])
							if iNode<(len(Ns[nextNode])-1):
								newString+=","
						newString+="}"
						strings.append(newString)
					if errorsOn and (not children[nextNode]) and errors[nextNode]:
						newString="errors={"
						for iNode in range(len(errors[nextNode])):
							mutation=errors[nextNode][iNode]
							newString+=allelesListExt[mutation[0]]+str(mutation[1])+allelesListExt[mutation[2]]+":"+str(mutation[3])
							if iNode<(len(errors[nextNode])-1):
								newString+=","
						newString+="}"
						strings.append(newString)
			elif up[nextNode]==None and estimateMATon:
				newString="rootState={"
				currentPos=0
				firstDoneEntry=False
				rootVect=rootVector(tree.probVect[nextNode],False,(len(children[nextNode])==0 and len(tree.minorSequences[nextNode])==0),tree,nextNode )
				for entry in rootVect:
					if entry[0]!=4:
						if not firstDoneEntry:
							firstDoneEntry=True
						else:
							newString+=","
					if entry[0]==5:
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
									newString+=","
								newString+=allelesList[i]+str(currentPos+1)+":"+str(vect[i])
						currentPos+=1
					elif entry[0]<4:
						newString+=allelesList[entry[0]]+str(currentPos+1)+":1.0"
						currentPos+=1
					else:
						currentPos=entry[1]
				newString+="}"
				strings.append(newString)
				if aBayesPlusActive and rootSupport!=None and rootSupport[nextNode]!=None:
					strings.append("rootSupport="+str(rootSupport[nextNode]))
			elif up[nextNode]==None and aBayesPlusActive and rootSupport!=None and rootSupport[nextNode]!=None:
				strings.append("rootSupport="+str(rootSupport[nextNode]))
			if printIQtreeSupportForNode:
				strings.append("IQsupport="+str(IQsupport[nextNode]))
		elif writeLineageAssignment and (lineage[nextNode]!=None or lineages[nextNode]!=None):
			if lineage[nextNode]!=None:
				strings.append("lineage="+lineage[nextNode])
			if lineages[nextNode]!=None and lineages:
				newString="lineages={"
				for lineageName in lineages[nextNode].keys():
					newString+=lineageName+":"+str(lineages[nextNode][lineageName])
					newString+=","
				newString=newString[:-1]
				newString+="}"
				strings.append(newString)
		finalString=""
		if networkOutput or (not children[nextNode]):
			finalString=nameNode
		if strings:
			finalString+="[&"
			for s in range(len(strings)):
				finalString+=strings[s]
				if s<(len(strings)-1):
					finalString+=","
			finalString+="]"
		return finalString


	#TODO write also function to write a timed tree, with branch lengths representing time.
	#TODO
	#TODO
	#create newick string of a given tree (input node is assumed to be the root) - with option "binary" the generated tree is binary (polytomies are represented with branches of length 0).
	def createNewick(tree,node,binary=True,namesInTree=None,includeMinorSeqs=True,estimateMAT=estimateMAT,networkOutput=networkOutput,aBayesPlusOn=aBayesPlusOn, performLineageAssignmentByRefPlacement= performLineageAssignmentByRefPlacement):
		nextNode=node
		stringList=[]
		direction=0
		numLeaves=0
		up=tree.up
		children=tree.children
		dist=tree.dist
		name=tree.name
		minorSequences=tree.minorSequences
		writeLineageAssignment = performLineageAssignment or performLineageAssignmentByRefPlacement
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
						if aBayesPlusOn or estimateMAT or writeLineageAssignment:
							stringList.append(stringForNode(tree,nextNode,"",dist[nextNode],estimateMAT=estimateMAT,networkOutput=networkOutput,aBayesPlusOn=aBayesPlusOn,namesInTree=namesInTree, performLineageAssignmentByRefPlacement = performLineageAssignmentByRefPlacement))
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
						if supportForIdenticalSequences or writeLineageAssignment:
							if namesInTree==None:
								stringList.append(stringForNode(tree,nextNode,name[nextNode],0.0,estimateMAT=estimateMAT,networkOutput=networkOutput,aBayesPlusOn=aBayesPlusOn,namesInTree=namesInTree, performLineageAssignmentByRefPlacement = performLineageAssignmentByRefPlacement))
							else:
								stringList.append(stringForNode(tree,nextNode,namesInTree[name[nextNode]],0.0,estimateMAT=estimateMAT,networkOutput=networkOutput,aBayesPlusOn=aBayesPlusOn,namesInTree=namesInTree, performLineageAssignmentByRefPlacement = performLineageAssignmentByRefPlacement))
						else:
							if namesInTree==None:
								stringList.append(name[nextNode])
							else:
								if name[nextNode]!="":
									stringList.append(namesInTree[name[nextNode]])
						stringList.append(":")
						for s2 in minorSequences[nextNode][:-1]:
							stringList.append("0.0,")
							if supportForIdenticalSequences or writeLineageAssignment:
								if namesInTree==None:
									stringList.append(stringForNode(tree,nextNode,s2,0.0,estimateMAT=estimateMAT,networkOutput=networkOutput,aBayesPlusOn=aBayesPlusOn,namesInTree=namesInTree, performLineageAssignmentByRefPlacement = performLineageAssignmentByRefPlacement))
								else:
									stringList.append(stringForNode(tree,nextNode,namesInTree[s2],0.0,estimateMAT=estimateMAT,networkOutput=networkOutput,aBayesPlusOn=aBayesPlusOn,namesInTree=namesInTree, performLineageAssignmentByRefPlacement = performLineageAssignmentByRefPlacement))
							else:
								if namesInTree==None:
									stringList.append(s2)
								else:
									stringList.append(namesInTree[s2])
							stringList.append(":0.0):")
						stringList.append("0.0,")
						if supportForIdenticalSequences or writeLineageAssignment:
							if namesInTree==None:
								stringList.append(stringForNode(tree,nextNode,minorSequences[nextNode][-1],0.0,estimateMAT=estimateMAT,networkOutput=networkOutput,aBayesPlusOn=aBayesPlusOn,namesInTree=namesInTree, performLineageAssignmentByRefPlacement = performLineageAssignmentByRefPlacement))
							else:
								stringList.append(stringForNode(tree,nextNode,namesInTree[minorSequences[nextNode][-1]],0.0,estimateMAT=estimateMAT,networkOutput=networkOutput,aBayesPlusOn=aBayesPlusOn,namesInTree=namesInTree, performLineageAssignmentByRefPlacement = performLineageAssignmentByRefPlacement))
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
						if supportForIdenticalSequences or writeLineageAssignment:
							if namesInTree==None:
								stringList.append(stringForNode(tree,nextNode,name[nextNode],0.0,estimateMAT=estimateMAT,networkOutput=networkOutput,aBayesPlusOn=aBayesPlusOn,namesInTree=namesInTree, performLineageAssignmentByRefPlacement = performLineageAssignmentByRefPlacement))
							else:
								stringList.append(stringForNode(tree,nextNode,namesInTree[name[nextNode]],0.0,estimateMAT=estimateMAT,networkOutput=networkOutput,aBayesPlusOn=aBayesPlusOn,namesInTree=namesInTree, performLineageAssignmentByRefPlacement = performLineageAssignmentByRefPlacement))
						else:
							if namesInTree==None:
								stringList.append(name[nextNode])
							else:
								if name[nextNode]!="":
									stringList.append(namesInTree[name[nextNode]])
						stringList.append(":0.0")
						for s2 in minorSequences[nextNode]:
							stringList.append(",")
							if supportForIdenticalSequences or writeLineageAssignment:
								if namesInTree==None:
									stringList.append(stringForNode(tree,nextNode,s2,0.0,estimateMAT=estimateMAT,networkOutput=networkOutput,aBayesPlusOn=aBayesPlusOn,namesInTree=namesInTree, performLineageAssignmentByRefPlacement = performLineageAssignmentByRefPlacement))
								else:
									stringList.append(stringForNode(tree,nextNode,namesInTree[s2],0.0,estimateMAT=estimateMAT,networkOutput=networkOutput,aBayesPlusOn=aBayesPlusOn,namesInTree=namesInTree, performLineageAssignmentByRefPlacement = performLineageAssignmentByRefPlacement))
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
				if aBayesPlusOn or estimateMAT or writeLineageAssignment:
					stringList.append(stringForNode(tree,nextNode,"",dist[nextNode],estimateMAT=estimateMAT,networkOutput=networkOutput,aBayesPlusOn=aBayesPlusOn,namesInTree=namesInTree, performLineageAssignmentByRefPlacement = performLineageAssignmentByRefPlacement))
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


	#Function to calculate number of branches under each multifurcation - to be used for HnZ modifier
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
						if dist[children[node][i]]>effectivelyNon0BLen:
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
if __name__ == "__main__":
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
		trees, namesInTree, namesInTreeDict=readNewick(inputTree,createDict=True,onlyTerminalNodeName=True)
		tree1, rootIndex1=trees[0]
		print("Read input newick tree")
		leafNameDict, nodeTable, leafCount, numBranches, leafDistDict, branchLengthDict, sumBranchLengths = prepareTreeComparison(tree1,rootIndex1,namesInTree,namesInTreeDict,rooted=False)
		otherTrees=readNewick(inputRFtrees,multipleTrees=multipleInputRFTrees,inputDictNames=namesInTreeDict,onlyTerminalNodeName=True)
		print("Read other input newick trees to be compared to the first one")
		file=open(outputFile+"_RFdistances.txt","w")
		file.write("RF\t"+"normalisedRF\t"+"leaves\t"+"foundBranches\t"+"missedBranches\t"+"notFoundBranches\t"+"RFL\n")
		for treePair in otherTrees:
			tree, rootIndex=treePair
			numDiffs, normalisedRF, leafCount, foundBranches, missedBranches, notFoundBranches, RFL = RobinsonFouldsWithDay1985(tree,rootIndex,leafNameDict, nodeTable, leafCount, numBranches,leafDistDict, branchLengthDict, sumBranchLengths,rooted=False)
			file.write(str(numDiffs)+"\t"+str(normalisedRF)+"\t"+str(leafCount)+"\t"+str(foundBranches)+"\t"+str(missedBranches)+"\t"+str(notFoundBranches)+"\t"+str(RFL)+"\n")
		print("Comparison ended and results written to "+outputFile+"_RFdistances.txt")
		file.close()
		exit()


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
def readConciseAlignment(fileName,extractReference=True,ref="",onlyRef=False): #extractNames=False
	if fileName.endswith(".gz"):
		import gzip
		fileI=gzip.open(fileName, 'rt')
	else:
		fileI=open(fileName)
	line=fileI.readline()
	if extractReference:
		line=fileI.readline()
		ref=""
		while line!="" and line[0]!=">":
			ref+=line.replace("\n","")
			line=fileI.readline()
		ref=ref.lower()
	if onlyRef:
		return ref
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


#read reference genome
if __name__ == "__main__":
	if refFile=="":
		ref, data=readConciseAlignment(inputFile) #extractNames=extractNamesFlag
	else:
		ref=collectReference(refFile)
		data=readConciseAlignment(inputFile, extractReference=False, ref=ref) #,extractNames=extractNamesFlag

	# read lineage references
	#NHAN
	if performLineageAssignmentByRefPlacement:
		from joblib import Parallel, delayed

		# don't allow running two lineage assignment methods at the same time
		if assignmentFile != "" and assignmentFileCSV != "":
			print("Please only use one among these options: --assignmentFile or --assignmentFileCSV or --lineageRefs.")
			raise Exception("exit")

		# make sure users specify a tree
		if (not os.path.isfile(inputTree)):
			print("Input tree in newick format " + inputTree + " not found, quitting MAPLE lineage assignment. Use option --inputTree to specify a valid input tree file.")
			raise Exception("exit")

		# check if file exits
		if not os.path.isfile(lineageRefs):
			print("Lineage reference file in Maple format " + lineageRefs + " not found.")
			raise Exception("exit")

		# don't allow rerooting the tree -> the program terminates immediately after lineage assignment
		# doNotReroot = True

		# read the lineage reference genomes
		if refFile == "":
			ref2, lineageRefData = readConciseAlignment(lineageRefs)
		else:
			ref2 = collectReference(refFile)
			lineageRefData = readConciseAlignment(lineageRefs, extractReference=False, ref=ref2)

		# make sure lineageRefs uses the same ref genome with the input alignment
		if ref2 != ref:
			refSrcFile = refFile
			if refSrcFile == "":
				refSrcFile = inputFile
			print("Reference genome in ", lineageRefs, " is different from that of ", refSrcFile)
			raise Exception("exit")
else:
	if refFile=="":
		ref=readConciseAlignment(inputFile,onlyRef=True) #extractNames=extractNamesFlag
	else:
		ref=collectReference(refFile)
lRef=len(ref)
print("Length of reference genome: "+str(lRef))
globalTotRate=-float(lRef)
logLRef=log(lRef)
thresholdLogLKoptimizationTopology*=logLRef
thresholdLogLKoptimization*=logLRef
thresholdLogLKtopology*=logLRef
thresholdLogLK*=logLRef
thresholdLogLKtopologyInitial*=logLRef
effectivelyNon0BLen=1.0/(10*lRef)
lineageRefsThresh /= lRef
oneMutBLen=1.0/lRef
minBLenSensitivity*=oneMutBLen
if errorRateInitial:
	errorRateGlobal=errorRateInitial
else:
	errorRateGlobal=1.0/lRef
minimumCarryOver=sys.float_info.min*(1e50)


#Read input tree
if __name__ == "__main__":
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
		trees, namesInTree, namesInTreeDict=readNewick(inputTree,dirtiness=largeUpdate,createDict=True,onlyTerminalNodeName=forgetInputTreeInternalNodeNames)
		tree1,rootIndex1=trees[0]
		print("Read input newick tree")
		makeTreeBinary(tree1,rootIndex1)
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
#ambiguities={"y":[0.0,0.5,0.0,0.5],"r":[0.5,0.0,0.5,0.0],"w":[0.5,0.0,0.0,0.5],"s":[0.0,0.5,0.5,0.0],"k":[0.0,0.0,0.5,0.5],"m":[0.5,0.5,0.0,0.0],"d":[1.0/3,0.0,1.0/3,1.0/3],"v":[1.0/3,1.0/3,1.0/3,0.0],"h":[1.0/3,1.0/3,0.0,1.0/3],"b":[0.0,1.0/3,1.0/3,1.0/3]}
ambiguities={"y":[0.0,1.0,0.0,1.0],"r":[1.0,0.0,1.0,0.0],"w":[1.0,0.0,0.0,1.0],"s":[0.0,1.0,1.0,0.0],"k":[0.0,0.0,1.0,1.0],"m":[1.0,1.0,0.0,0.0],"d":[1.0,0.0,1.0,1.0],"v":[1.0,1.0,1.0,0.0],"h":[1.0,1.0,0.0,1.0],"b":[0.0,1.0,1.0,1.0]}


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


if __name__ == "__main__":
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
			if currPos>pos: #region where the node is identical to the ref.
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
if __name__ == "__main__":
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
#(and also calculate the logLk of the merging if necessary).
# fromTip1 and fromTip1 tell us if the two vectors come from tip nodes.
def mergeVectors(probVect1,bLen1,fromTip1,probVect2,bLen2,fromTip2,returnLK=False,isUpDown=False,numMinor1=0,numMinor2=0,errorRateGlobalPassed=None,mutMatrixGlobalPassed=None,errorRatesGlobal=None,mutMatricesGlobal=None,cumulativeRateGlobal=None,cumulativeErrorRateGlobal=None):
	indexEntry1, indexEntry2, pos, totalFactor = 0, 0, 0, 1.0
	probVect=[]
	entry1=probVect1[indexEntry1]
	entry2=probVect2[indexEntry2]
	totSum, cumErrorRate, refNucToPass, flag1, flag2, i, j ,i1, i2, newPos = 0.0, 0.0, -1, False, False, 0, 0, 0, 0, 0
	newVec, newVec2 = [], []
	if useRateVariation:
		if mutMatricesGlobal!=None:
			mutMatricesUsed=mutMatricesGlobal
		else:
			mutMatricesUsed=mutMatrices
	else:
		if mutMatrixGlobalPassed!=None:
			mutMatrix=mutMatrixGlobalPassed
		else:
			mutMatrix=mutMatrixGlobal

	if usingErrorRate and errorRateSiteSpecific:
		if errorRatesGlobal!=None:
			errorRatesUsed=errorRatesGlobal
		else:
			errorRatesUsed=errorRates
	else:
		if errorRateGlobalPassed==None:
			errorRate=errorRateGlobal
		else:
			errorRate=errorRateGlobalPassed

	if returnLK:
		if cumulativeRateGlobal!=None:
			cumulativeRateUsed=cumulativeRateGlobal
		else:
			cumulativeRateUsed=cumulativeRate
		if usingErrorRate and errorRateSiteSpecific:
			if cumulativeErrorRateGlobal!=None:
				cumulativeErrorRateUsed=cumulativeErrorRateGlobal
			else:
				cumulativeErrorRateUsed=cumulativeErrorRate

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
						mutMatrix=mutMatricesUsed[pos]
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
				cumulPartLk+=(bLen1+bLen2)*(cumulativeRateUsed[pos]-cumulativeRateUsed[newPos])
				if usingErrorRate:
					if fromTip1 or fromTip2:
						if errorRateSiteSpecific: cumErrorRate =  cumulativeErrorRateUsed[newPos] - cumulativeErrorRateUsed[pos] #both this and the next one should end up being positive
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
						mutMatrix=mutMatricesUsed[pos]
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
				cumulPartLk+=(bLen1+bLen2)*(cumulativeRateUsed[pos]-cumulativeRateUsed[newPos])
				if usingErrorRate:
					if fromTip1 or fromTip2:
						if errorRateSiteSpecific: cumErrorRate = cumulativeErrorRateUsed[newPos] - cumulativeErrorRateUsed[pos] #both this and the next one should end up being positive
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
						cumulPartLk+=(totLen2-bLen2+totLen1-bLen1)*(cumulativeRateUsed[newPos]-cumulativeRateUsed[pos])
						if usingErrorRate:
							if ((not fromTip1) and flag1) or ((not fromTip2) and flag2):
								if errorRateSiteSpecific: cumErrorRate = cumulativeErrorRateUsed[pos] - cumulativeErrorRateUsed[newPos] #both this and the next one should end up being negative
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
						cumulPartLk-=mutMatricesUsed[pos][refNucToPass][refNucToPass]*(bLen2+bLen1)
					else:
						cumulPartLk-=mutMatrix[refNucToPass][refNucToPass]*(bLen2+bLen1)
					if usingErrorRate and ((entry1[0] != entry2[0]) or entry1[0]==6) and (fromTip1 or fromTip2):
						if errorRateSiteSpecific: cumErrorRate = errorRatesUsed[pos]
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
							cumulPartLk+=mutMatricesUsed[pos][entry1[0]][entry1[0]]*(totLen1+totLen2)
						else:
							cumulPartLk+=mutMatrix[entry1[0]][entry1[0]]*(totLen1+totLen2)
						if usingErrorRate:
							if ((not fromTip1) and flag1) or ((not fromTip2) and flag2):
								if errorRateSiteSpecific: cumErrorRate = errorRatesUsed[pos]
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
				if usingErrorRate and errorRateSiteSpecific: errorRate = errorRatesUsed[pos]
				if useRateVariation:	mutMatrix=mutMatricesUsed[pos]

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
def rootVector(probVect,bLen,isFromTip,tree,node,mutMatrixGlobalPassed=None,mutMatricesGlobal=None):
	if useRateVariation:
		if mutMatricesGlobal!=None:
			mutMatricesUsed=mutMatricesGlobal
		else:
			mutMatricesUsed=mutMatrices
	else:
		if mutMatrixGlobalPassed!=None:
			mutMatrix=mutMatrixGlobalPassed
		else:
			mutMatrix=mutMatrixGlobal

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
					newVect=getPartialVec(6, totBLen, mutMatricesUsed[newPos], 0, vect=entry[-1])
				else:
					newVect=getPartialVec(6, totBLen, mutMatrix, 0, vect=entry[-1])
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


if __name__ == "__main__":
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
def estimateBranchLengthWithDerivative(probVectP,probVectC,fromTipC=False,errorRateGlobalPassed=None,mutMatrixGlobalPassed=None,errorRatesGlobal=None,mutMatricesGlobal=None,cumulativeRateGlobal=None):
	if useRateVariation:
		if mutMatricesGlobal!=None:
			mutMatricesUsed=mutMatricesGlobal
		else:
			mutMatricesUsed=mutMatrices
	else:
		if mutMatrixGlobalPassed!=None:
			mutMatrix=mutMatrixGlobalPassed
		else:
			mutMatrix=mutMatrixGlobal

	if usingErrorRate and errorRateSiteSpecific:
		if errorRatesGlobal!=None:
			errorRatesUsed=errorRatesGlobal
		else:
			errorRatesUsed=errorRates
	else:
		if errorRateGlobalPassed==None:
			errorRate=errorRateGlobal
		else:
			errorRate=errorRateGlobalPassed
	if cumulativeRateGlobal!=None:
		cumulativeRateUsed=cumulativeRateGlobal
	else:
		cumulativeRateUsed=cumulativeRate
	
	# if not (usingErrorRate and errorRateSiteSpecific): 
	# 	errorRate=errorRateGlobal
	# if not useRateVariation:
	# 	mutMatrix=mutMatrixGlobal
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
			c1+=(cumulativeRateUsed[pos]-cumulativeRateUsed[end])
			pos=end
		elif entry1[0]==5: # case entry1 or entry2 is N
			#if parent node is type "N", in theory we might have to calculate the contribution of root nucleotides; 
			# however, if this node is "N" then every other node in the current tree is "N", so we can ignore this since this contribution cancels out in relative terms.
			if entry2[0]==4:
				end=min(entry1[1],entry2[1])
			else:
				end=pos+1
			c1+=(cumulativeRateUsed[pos]-cumulativeRateUsed[end])
			pos=end
		else:
			#below, when necessary, we represent the likelihood as coeff0*l +coeff1, where l is the branch length to be optimized.
			if entry1[0]==4 and entry2[0]==4: # case entry1 and entry2 are R	
				pos=min(entry1[1],entry2[1])
			else:
				if useRateVariation:
					mutMatrix=mutMatricesUsed[pos]
				
				if entry1[0]==4:
					c1-=mutMatrix[entry2[1]][entry2[1]]
				else:
					c1-=mutMatrix[entry1[1]][entry1[1]]
				flag1 = (usingErrorRate and (entry1[0]!=6) and len(entry1)>2 and entry1[-1])
				flag2 = (usingErrorRate and (entry2[0]!=6) and (fromTipC or (len(entry2)>2 and entry2[-1])))
				if usingErrorRate and errorRateSiteSpecific: errorRate = errorRatesUsed[pos]

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

if __name__ == "__main__":
	#function to update nDesc when changing a branch length in the tree
	def updateNDesc0whenChangingDist(tree,node,newDist):
		dist=tree.dist
		nDesc0=tree.nDesc0
		up=tree.up
		if dist[node]>effectivelyNon0BLen and (newDist<=effectivelyNon0BLen):
			addendum0=(nDesc0[node]-1)
		elif (dist[node]<=effectivelyNon0BLen) and newDist>effectivelyNon0BLen:
			addendum0=(1-nDesc0[node])
		else:
			addendum0=0
		if addendum0:
			parent0=up[node]
			nDesc0[parent0]+=addendum0
			while up[parent0]!=None and (dist[parent0]<=effectivelyNon0BLen):
				parent0=up[parent0]
				nDesc0[parent0]+=addendum0
				if nDesc0[parent0]<=0:
					print("negative nDesc0 in updateNDesc0whenChangingDist")
					raise Exception("exit")


	#if updating genome lists in updatePartials() creates an inconsistency, this function can increase the length of a 0-length branch to resolve the inconsistency.
	#In doing so, it updates, the input list of nodes to visit and update.
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
		#updating nDesc0 in case of branch length change
		if HnZ:
			updateNDesc0whenChangingDist(tree,cNode,bestLength)
		dist[cNode]=bestLength
		dirty[node]=True
		dirty[cNode]=True
		if addToList:
			nodeList.append((cNode,2,True,doTimeTree))
			nodeList.append((node,cNum,True,doTimeTree))


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

if __name__ == "__main__":
	#TODO now updates also time likelihoods along the way, if needed.
	#TODO TODO TODO test changes
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
		if doTimeTree:
			probVectUpRightTime=tree.probVectUpRightTime
			probVectUpLeftTime=tree.probVectUpLeftTime
			probVectTime=tree.probVectTime
			probVectTotUpTime=tree.probVectTotUpTime
		while nodeList:
			updatedBLen=False # if there has been an inconsistency, function updateBLen() has been called, and so there is no point continuing with some updates.
			madeChange=False # some change has been made, so continue traversal
			#TODO adding flags for the time tree case, so to update only genetic or time likelihoods as needed. Needs to be implemented througout, also at the time of use of the function. 
			node, direction, Lkdirty, timeLKdirty = nodeList.pop()
			dirty[node]=True
			if up[node] != None:
				try:
					if node==children[up[node]][0]:
						childNumUp=0
						vectUpUp=probVectUpRight[up[node]]
						if doTimeTree:
							vectUpUpTime=probVectUpRightTime[up[node]]
					else:
						childNumUp=1
						vectUpUp=probVectUpLeft[up[node]]
						if doTimeTree:
							vectUpUpTime=probVectUpLeftTime[up[node]]
					if mutations[node] and Lkdirty:
						vectUpUp=passGenomeListThroughBranch(vectUpUp,mutations[node])
				except AttributeError:
					vectUpUp=None
			isTip=(len(children[node])==0 and len(minorSequences[node])==0)
			#change in likelihoods is coming from parent node
			if direction==2:
				if dist[node] or doTimeTree: #if necessary, update the total probabilities at the mid node.
					if Lkdirty:
						newTot=mergeVectors(vectUpUp,dist[node]/2,False,probVect[node],dist[node]/2,isTip,isUpDown=True)
						if newTot==None:
							if dist[node]>1e-100:
								print("inside updatePartials(), from parent: should not have happened since node.dist>0")
							updateBLen(tree,node)
							nodeList.append((up[node],childNumUp,True,doTimeTree))
							newTot=mergeVectors(vectUpUp,dist[node]/2,False,probVect[node],dist[node]/2,isTip,isUpDown=True)
							madeChange=True
							if doTimeTree:
								newVectTime,newVectTimeProb=mergeVectorsTime(vectUpUpTime,dist[node]/2,probVectTime[node],dist[node]/2,mutRate,isUpDown=True,returnLK=True)
								if isinstance(newVectTime, int):
									resolveTimeInconsistency(tree,node,newVectTime,mutRate)
									if node==children[up[node]][0]:
										vectUpUpTime=probVectUpRightTime[up[node]]
									else:
										vectUpUpTime=probVectUpLeftTime[up[node]]
									newVectTime,newVectTimeProb==mergeVectorsTime(vectUpUpTime,dist[node]/2,probVectTime[node],dist[node]/2,mutRate,isUpDown=True,returnLK=True)
								newVectTimeProb-=appendProbNodeTime(vectUpUpTime,probVectTime[node],mutRate,dist[node])
								probVectTotUpTime[node]=(newVectTime,newVectTimeProb)
						probVectTotUp[node]=newTot
						shorten(probVectTotUp[node])
					if timeLKdirty:
						newVectTime,newVectTimeProb=mergeVectorsTime(vectUpUpTime,dist[node]/2,probVectTime[node],dist[node]/2,mutRate,isUpDown=True,returnLK=True)
						if isinstance(newVectTime, int):
							resolveTimeInconsistency(tree,node,newVectTime,mutRate)
							if node==children[up[node]][0]:
								vectUpUpTime=probVectUpRightTime[up[node]]
							else:
								vectUpUpTime=probVectUpLeftTime[up[node]]
							newVectTime,newVectTimeProb=mergeVectorsTime(vectUpUpTime,dist[node]/2,probVectTime[node],dist[node]/2,mutRate,isUpDown=True,returnLK=True)
						newVectTimeProb-=appendProbNodeTime(vectUpUpTime,probVectTime[node],mutRate,dist[node])
						probVectTotUpTime[node]=(newVectTime,newVectTimeProb)
				else:
					probVectTotUp[node]=None

				if children[node]:# at valid internal node, update upLeft and upRight, and if necessary add children to nodeList.
					dist0=dist[children[node][0]]
					dist1=dist[children[node][1]]
					if Lkdirty:
						child0Vect=probVect[children[node][0]]
						if mutations[children[node][0]]:
							child0Vect=passGenomeListThroughBranch(child0Vect,mutations[children[node][0]],dirIsUp=True)
						child1Vect=probVect[children[node][1]]
						if mutations[children[node][1]]:
							child1Vect=passGenomeListThroughBranch(child1Vect,mutations[children[node][1]],dirIsUp=True)
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
									nodeList.append((up[node],childNumUp,True,doTimeTree))
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
										nodeList.append((up[node],childNumUp,True,doTimeTree))
										madeChange=True
								else:
									print("Strange: None vector from non-zero distances in updatePartials() from parent direction, child0.")
									raise Exception("exit")

					if not updatedBLen:
						upRightChangedTime=False
						upLeftChangedTime=False
						if doTimeTree:
							if madeChange:
								newVectTime,newVectTimeProb=mergeVectorsTime(vectUpUpTime,dist[node]/2,probVectTime[node],dist[node]/2,mutRate,isUpDown=True,returnLK=True)
								if isinstance(newVectTime, int):
									resolveTimeInconsistency(tree,node,newVectTime,mutRate)
									if node==children[up[node]][0]:
										vectUpUpTime=probVectUpRightTime[up[node]]
									else:
										vectUpUpTime=probVectUpLeftTime[up[node]]
									newVectTime,newVectTimeProb=mergeVectorsTime(vectUpUpTime,dist[node]/2,probVectTime[node],dist[node]/2,mutRate,isUpDown=True,returnLK=True)
								newVectTimeProb-=appendProbNodeTime(vectUpUpTime,probVectTime[node],mutRate,dist[node])
								probVectTotUpTime[node]=(newVectTime,newVectTimeProb)
							if timeLKdirty or madeChange:
								child0Vect=probVectTime[children[node][0]]
								child1Vect=probVectTime[children[node][1]]
								newUpRightTime=mergeVectorsTime(vectUpUpTime,dist[node],child1Vect,dist1,mutRate,isUpDown=True)
								if isinstance(newUpRightTime, int):
									resolveTimeInconsistency(tree,node,newUpRightTime,mutRate)
									if node==children[up[node]][0]:
										vectUpUpTime=probVectUpRightTime[up[node]]
									else:
										vectUpUpTime=probVectUpLeftTime[up[node]]
									newUpRightTime=mergeVectorsTime(vectUpUpTime,dist[node],child1Vect,dist1,mutRate,isUpDown=True)
								newUpLeftTime=mergeVectorsTime(vectUpUpTime,dist[node],child0Vect,dist0,mutRate,isUpDown=True)
								if isinstance(newUpLeftTime, int):
									resolveTimeInconsistency(tree,node,newUpLeftTime,mutRate)
									if node==children[up[node]][0]:
										vectUpUpTime=probVectUpRightTime[up[node]]
									else:
										vectUpUpTime=probVectUpLeftTime[up[node]]
									newUpLeftTime=mergeVectorsTime(vectUpUpTime,dist[node],child0Vect,dist0,mutRate,isUpDown=True)
								if areVectorsDifferentTime(probVectUpRightTime[node],newUpRightTime):
									upRightChangedTime=True
									probVectUpRightTime[node]=newUpRightTime
								if areVectorsDifferentTime(probVectUpLeftTime[node],newUpLeftTime):
									upLeftChangedTime=True
									probVectUpLeftTime[node]=newUpLeftTime

						upRightChanged=False
						upLeftChanged=False
						if Lkdirty:
							if madeChange or areVectorsDifferent(probVectUpRight[node],newUpRight):
								probVectUpRight[node]=newUpRight
								shorten(probVectUpRight[node])
								upRightChanged=True
							if madeChange or areVectorsDifferent(probVectUpLeft[node],newUpLeft):
								probVectUpLeft[node]=newUpLeft
								shorten(probVectUpLeft[node])
								upLeftChanged=True
						
						if upRightChanged or upRightChangedTime:
							nodeList.append((children[node][0],2,upRightChanged,upRightChangedTime))
						if upLeftChanged or upLeftChangedTime:
							nodeList.append((children[node][1],2,upLeftChanged,upLeftChangedTime))
					
			else: #change in likelihoods is coming from child number "direction".
				childNum=direction
				otherChildNum=1-childNum
				childDist=dist[children[node][childNum]]
				otherChildDist=dist[children[node][otherChildNum]]
				if Lkdirty:
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
								nodeList.append((children[node][childNum],2,True,doTimeTree))
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
					if (not updatedBLen) and (dist[node] or doTimeTree) and (up[node] != None) and (vectUpUp!=None):
						newTot=mergeVectors(vectUpUp,dist[node]/2,False,probVect[node],dist[node]/2,False,isUpDown=True)
						if newTot==None:
							updateBLen(tree,node)
							probVect[node]=mergeVectors(otherChildVect,otherChildDist,otherIsTip,probVectDown,childDist,isTip)
							nodeList.append((children[node][childNum],2,True,doTimeTree))
							probVectTotUp[node]=mergeVectors(vectUpUp,dist[node]/2,False,probVect[node],dist[node]/2,False,isUpDown=True)
							madeChange=True
							print("inside updatePartials(), from child: should not have happened since node.dist>0")
						else:
							probVectTotUp[node]=newTot
							shorten(probVectTotUp[node])
					elif not dist[node]:
						probVectTotUp[node]=None
							
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
									nodeList.append((children[node][childNum],2,True,doTimeTree))
									madeChange=True
									newUpVect=mergeVectors(vectUpUp,dist[node],False,probVectDown,childDist,isTip,isUpDown=True)
							else:
								print("Strange: None vector from non-zero distances in updatePartials() from child direction, newUpVect.")
								raise Exception("exit")
							
				if not updatedBLen:
					upChangedTime=False
					downChangedTime=False
					if doTimeTree:
						if timeLKdirty or madeChange:
							otherChildVectTime=probVectTime[children[node][otherChildNum]]
							probVectDownTime=probVectTime[children[node][childNum]]
							try:
								if childNum:
									otherVectUpTime=probVectUpRightTime[node]
								else:
									otherVectUpTime=probVectUpLeftTime[node]
							except AttributeError:
								otherVectUpTime=None
							#update lower likelihoods
							newVect=mergeVectorsTime(otherChildVectTime,otherChildDist,probVectDownTime,childDist,mutRate)
							try:
								oldProbVectTime=probVectTime[node]
							except AttributeError:
								oldProbVectTime=None
							probVectTime[node]=newVect
							if (up[node] != None):
								newTot,newTotProb=mergeVectorsTime(vectUpUpTime,dist[node]/2,probVectTime[node],dist[node]/2,mutRate,isUpDown=True,returnLK=True)
								if isinstance(newTot, int):
									resolveTimeInconsistency(tree,node,newTot,mutRate)
									if node==children[up[node]][0]:
										vectUpUpTime=probVectUpRightTime[up[node]]
									else:
										vectUpUpTime=probVectUpLeftTime[up[node]]
									newTot,newTotProb=mergeVectorsTime(vectUpUpTime,dist[node]/2,probVectTime[node],dist[node]/2,mutRate,isUpDown=True,returnLK=True)
								newTotProb-=appendProbNodeTime(vectUpUpTime,probVectTime[node],mutRate,dist[node])
								probVectTotUpTime[node]=(newTot,newTotProb)
								newUpVectTime=mergeVectorsTime(vectUpUpTime,dist[node],probVectDownTime,childDist,mutRate,isUpDown=True)
								if isinstance(newUpVectTime, int):
									resolveTimeInconsistency(tree,node,newUpVectTime,mutRate)
									if node==children[up[node]][0]:
										vectUpUpTime=probVectUpRightTime[up[node]]
									else:
										vectUpUpTime=probVectUpLeftTime[up[node]]
									newUpVectTime=mergeVectorsTime(vectUpUpTime,dist[node],probVectDownTime,childDist,mutRate,isUpDown=True)
							else:
								newUpVectTime=rootVectorTime(probVectDownTime,childDist,mutRate)

							if areVectorsDifferentTime(otherVectUpTime,newUpVectTime):
								upChangedTime=True
							if areVectorsDifferentTime(probVectTime[node],oldProbVectTime):
								downChangedTime=True
							if childNum:
								probVectUpRightTime[node]=newUpVectTime
							else:
								probVectUpLeftTime[node]=newUpVectTime

					upChanged=False
					downChanged=False
					if Lkdirty:
						if otherVectUp!=None:
							if madeChange or areVectorsDifferent(otherVectUp,newUpVect):
								upChanged=True
								if childNum:
									probVectUpRight[node]=newUpVect
									shorten(probVectUpRight[node])
								else:
									probVectUpLeft[node]=newUpVect
									shorten(probVectUpLeft[node])
						#update likelihoods at parent node
						if madeChange or areVectorsDifferent(probVect[node],oldProbVect):
							downChanged=True
					if up[node] != None:
						if downChanged or downChangedTime:
							nodeList.append((up[node],childNumUp,downChanged,downChangedTime))
					if upChanged or upChangedTime:
						nodeList.append((children[node][otherChildNum],2,upChanged,upChangedTime))
		return


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
if __name__ == "__main__":
	#Given a tree, and a substitution rate matrix, re-calculate all genome lists within the tree according to this matrix.
	# this is useful once the matrix estimation has finished, to make sure all genome lists reflect this matrix. 
	#TODO TODO TODO check if additions made for time probabilities work
	def reCalculateAllGenomeLists(tree,root, checkExistingAreCorrect=False,countNodes=False,countPseudoCounts=False,pseudoMutCounts=None,data=None,dates=None,names=None,firstSetUp=False):
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
			if doTimeTree:
				tree.dateData=[False]*len(up)
				dateData=tree.dateData
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
							raise Exception("exit")
						if doTimeTree and dates==None:
							print("Problem, initializing dates but there is no date data.")
							raise Exception("exit")
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
						if doTimeTree:
							if (names[name[node]] in dates):
								dateData[node]=dates[names[name[node]]]
							else:
								print("No date for sample "+names[name[node]]+", treating it as an unknown date.")
								dateData[node]=None
						
						# removed minor sequences from the tree
						tryRemovingMinor=False
						if children[up[node]][1]==node and (not dist[node]):
							sibling=children[up[node]][0]
							if (not dist[sibling]) and (not children[sibling]):
								tryRemovingMinor=True
						while tryRemovingMinor:
							#test if one of the samples is minor of the other
							#also accounting for HnZ since minor sequences might have different modifiers for different placements
							comparison2=0
							if errorRateSiteSpecificFile or errorRateFixed or estimateErrorRate or estimateSiteSpecificErrorRate or supportFor0Branches or HnZ:
								comparison=isMinorSequence(probVect[node],probVect[sibling],onlyFindIdentical=True)
								if doTimeTree:
									comparison2=isMinorDate(dateData[node],dateData[sibling],onlyFindIdentical=True)
							else:
								comparison=isMinorSequence(probVect[node],probVect[sibling])
								if doTimeTree:
									comparison2=isMinorDate(dateData[node],dateData[sibling])

							if comparison==1 and ((not doTimeTree) or (comparison2==1)):
								majorNode=node
								minorNode=sibling
							elif comparison==2 and ((not doTimeTree) or (comparison2==2)):
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
		#re-set up nDesc0 since some sequences might have become minor sequences
		if firstSetUp and HnZ:
			calculateNDesc0(tree,root)

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
					if dist[node] or doTimeTree:
						isTip=(len(children[node])==0 ) and (len(minorSequences[node])==0 )
						if dist[node] and countPseudoCounts:
							updatePesudoCounts(vectUp,probVect[node],pseudoMutCounts)
						newVect=mergeVectors(vectUp,dist[node]/2,False,probVect[node],dist[node]/2,isTip,isUpDown=True)
						shorten(newVect)
						if checkExistingAreCorrect:
							if areVectorsDifferentDebugging(newVect,probVectTotUp[node]):
								print("new probVectTotUp at node is different from the old one, and it shouldn't be.")
								print(newVect)
								print(probVectTotUp[node])
								print(node)
								raise Exception("exit")
						probVectTotUp[node]=newVect
					else:
						probVectTotUp[node]=None
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
									totNodeList.append((node,1,True,doTimeTree))
								else:
									probVectTotUp[node]=mergeVectors(vectUp,dist[node]/2,False,probVect[node],dist[node]/2,False,isUpDown=True)
									totNodeList.append((up[node],nodeChildNum,True,doTimeTree))
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
									totNodeList.append((up[node],nodeChildNum,True,doTimeTree))
									probVectTotUp[node]=mergeVectors(vectUp,dist[node]/2,False,probVect[node],dist[node]/2,isTip,isUpDown=True)
									probVectUpRight[node]=mergeVectors(vectUp,dist[node],False,probVect1,dist[children[node][1]],isTip1,isUpDown=True)
								else:
									totNodeList.append((node,0,True,doTimeTree))
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
		nonMutRates=[0.0]*4
		for i in range4:
			nonMutRates[i]=mutMatrix[i][i]
		cumulativeRate=[0.0]*(lRef+1)
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
			rootFreqsLogErrorCumulative=[0.0]*(lRef+1)
			for i in range(0,len(errorRates)):
				rootFreqsLogErrorCumulative[i+1] = rootFreqsLogErrorCumulative[i] + log(rootFreqs[refIndeces[i]]*(1.0-1.33333*errorRates[i])+0.333333*errorRates[i])
			totError=-cumulativeErrorRate[-1]
		else:
			rootFreqsLogErrorCumulative=[0.0]*(lRef+1)
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
		if not inputRates:
			reCalculateAllGenomeLists(tree1,rootIndex1,countPseudoCounts=True,pseudoMutCounts=pseudoMutCounts,data=data,dates=dates,names=namesInTree,firstSetUp=True)
			if (model!="JC" and updateSubMatrix(pseudoMutCounts,model,mutMatrixGlobal)):
				for i in range4:
					nonMutRates[i]=mutMatrixGlobal[i][i]
				for i in range(lRef):
					cumulativeRate[i+1]=cumulativeRate[i]+nonMutRates[refIndeces[i]]
			reCalculateAllGenomeLists(tree1,rootIndex1)
		else:
			reCalculateAllGenomeLists(tree1,rootIndex1,data=data,dates=dates,names=namesInTree,firstSetUp=True)
		print("Genome list for initial tree and initial pseudocounts calculated.")
		if doTimeTree:
			reCalculateAllGenomeListsTime(tree1,rootIndex1,mutRate)
			print("Time LKs for initial tree calculated.")


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
# errorRates errorRateGlobal mutMatrixGlobal mutMatrices
def appendProbNode(probVectP,probVectC,isTipC,bLen, errorRateGlobalPassed=None,mutMatrixGlobalPassed=None,errorRatesGlobal=None,mutMatricesGlobal=None,totErrorPassed=None):
	if useRateVariation:
		if mutMatricesGlobal!=None:
			mutMatricesUsed=mutMatricesGlobal
		else:
			mutMatricesUsed=mutMatrices
	else:
		if mutMatrixGlobalPassed!=None:
			mutMatrix=mutMatrixGlobalPassed
		else:
			mutMatrix=mutMatrixGlobal

	if usingErrorRate and errorRateSiteSpecific:
		if errorRatesGlobal!=None:
			errorRatesUsed=errorRatesGlobal
		else:
			errorRatesUsed=errorRates
	else:
		if errorRateGlobalPassed==None:
			errorRate=errorRateGlobal
		else:
			errorRate=errorRateGlobalPassed
	if usingErrorRate:
		if totErrorPassed!=None:
			totErrorUsed=totErrorPassed
		else:
			totErrorUsed=totError
	
	# if not (usingErrorRate and errorRateSiteSpecific): 
	# 	errorRate=errorRateGlobal
	# if not useRateVariation:
	# 	mutMatrix=mutMatrixGlobal
	indexEntry1, indexEntry2, totalFactor, pos = 0, 0, 1.0, 0
	entry1=probVectP[indexEntry1] #parent
	entry2=probVectC[indexEntry2] #child
	contribLength=bLen
	Lkcost=bLen*globalTotRate
	if usingErrorRate and isTipC:
		Lkcost+=totErrorUsed
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
						mutMatrix=mutMatricesUsed[pos]
					i1=entry2[1]
					if entry2[-1][i1]>0.02:
						totalFactor*=entry2[-1][i1]
					else:
						if len(entry1)==4+usingErrorRate:
							flag1 = (usingErrorRate and (len(entry1)>2) and entry1[-1] )
							tot=0.0
							if usingErrorRate and errorRateSiteSpecific: errorRate = errorRatesUsed[pos]
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
						mutMatrix=mutMatricesUsed[pos]
					if len(entry1)==4+usingErrorRate:
						flag1 = (usingErrorRate and (len(entry1)>2) and entry1[-1] )
						i1=entry2[1]
						i2=entry2[0]
						if usingErrorRate and errorRateSiteSpecific: errorRate = errorRatesUsed[pos]
						tot3=getPartialVec(i2, contribLength, mutMatrix, errorRate, flag=flag2)
						tot2=getPartialVec(i1, entry1[2], mutMatrix, errorRate, flag=flag1)
						tot=0.0
						for i in range4:
							tot+=tot3[i]*tot2[i]*rootFreqs[i]
						totalFactor*=tot/rootFreqs[i1]
					else:
						if flag2:
							if usingErrorRate and errorRateSiteSpecific: errorRate = errorRatesUsed[pos]
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
					mutMatrix=mutMatricesUsed[pos]
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
							if errorRateSiteSpecific: errorRate = errorRatesUsed[pos]
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
						mutMatrix=mutMatricesUsed[pos]

					i1=entry1[0]
					if entry2[0]<5: #entry2 is a nucleotide
						if entry2[0]==4:
							i2=entry1[1]
						else:
							i2=entry2[0]
						flag2 = (usingErrorRate and (isTipC or (len(entry2)>2) and entry2[-1] ))
						if len(entry1)==4+usingErrorRate:
							if usingErrorRate and errorRateSiteSpecific: errorRate = errorRatesUsed[pos]
							tot3=getPartialVec(i2, contribLength, mutMatrix, errorRate, flag=flag2)
							tot2=getPartialVec(i1, entry1[2], mutMatrix, errorRate, flag=flag1)
							tot=0.0
							for j in range4:
								tot+=rootFreqs[j]*tot3[j]*tot2[j]
							totalFactor*=tot/rootFreqs[i1]
						else:
							if flag1 or flag2:
								if errorRateSiteSpecific: errorRate = errorRatesUsed[pos]
								totalFactor*=(min(0.25,mutMatrix[i1][i2]*contribLength)+(flag1 + flag2)*0.33333*errorRate)
							else:
								if contribLength:
									totalFactor*=min(0.25,mutMatrix[i1][i2]*contribLength)
								else:
									return float("-inf")

					else: #entry1 is a nucleotide and entry2 is of type "O"
						if usingErrorRate and errorRateSiteSpecific: errorRate = errorRatesUsed[pos]
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
# TODO function traversing the tree to find the best node in the tree where to re-append the given subtree (rooted at node.children[child]) to improve the topology of the current tree.
# bestLKdiff is the best likelihood cost found for the current placement (after optimizing the branch length).
# removedBLen is such branch length that optimizes the current placement - it will be used to place the subtree attached at other nodes of the tree.
#TODO TODO TODO test changes
#TODO TODO TODO in case doTimeTree search within politomies but collapse these placements when doing SPRTA
# TODO TODO TODO in the rest of the code, in case doTimeTree store probVectTotUp also for branches of length 0 - different nodes in the same multifurcation might have different time likelihoods.
def findBestParentTopology(tree,node,child,bestLKdiff,removedBLen,mutRate=mutRate,strictTopologyStopRules=strictTopologyStopRules,allowedFailsTopology=allowedFailsTopology,thresholdLogLKtopology=thresholdLogLKtopology,errorRateGlobalPassed=None,mutMatrixGlobalPassed=None,errorRatesGlobal=None,mutMatricesGlobal=None,cumulativeRateGlobal=None,cumulativeErrorRateGlobal=None,totErrorPassed=None): #,sequentialSearch=True   
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
	nDesc0=tree.nDesc0
	if doTimeTree:
		probVectTime=tree.probVectTime
		probVectTotUpTime=tree.probVectTotUpTime
		probVectUpRightTime=tree.probVectUpRightTime
		probVectUpLeftTime=tree.probVectUpLeftTime
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
	if doTimeTree:
		removedPartialsRelativeTime=probVectTime[children[node][child]]

	if up[node]!=None:
		vectUpUpTime=None
		if children[up[node]][0]==node:
			childUp=1
			vectUpUp=probVectUpRight[up[node]]
			if doTimeTree:
				vectUpUpTime=probVectUpRightTime[up[node]]
		else:
			childUp=2
			vectUpUp=probVectUpLeft[up[node]]
			if doTimeTree:
				vectUpUpTime=probVectUpLeftTime[up[node]]
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
		if doTimeTree:
			probVect1Time=probVectTime[bestNode]
		else:
			probVect1Time=None
		nDesc0ToAdd=0
		if HnZ:
			if dist[node]<effectivelyNon0BLen:
				if dist[children[node][child]]>=effectivelyNon0BLen:
					nDesc0ToAdd=-1
				else:
					nDesc0ToAdd=-nDesc0[children[node][child]]
		if doTimeTree:
			nodesToVisit.append((up[node],childUp,probVect1,probVect1Time,dist[bestNode]+dist[node],bestLKdiff,0,removedPartialsRelative1,nDesc0ToAdd))
		else:
			nodesToVisit.append((up[node],childUp,probVect1,dist[bestNode]+dist[node],bestLKdiff,0,removedPartialsRelative1,nDesc0ToAdd))
		if mutations[node]:
			vectUpUp=passGenomeListThroughBranch(vectUpUp,mutations[node])
		removedPartialsRelative1=removedPartialsRelative
		if mutations[bestNode]:
			vectUpUp=passGenomeListThroughBranch(vectUpUp,mutations[bestNode])
			removedPartialsRelative1=passGenomeListThroughBranch(removedPartialsRelative,mutations[bestNode])
		nDesc0ToAdd=0
		if HnZ:
			if dist[bestNode]<effectivelyNon0BLen:
				if dist[children[node][child]]>=effectivelyNon0BLen:
					nDesc0ToAdd=-1
				else:
					nDesc0ToAdd=-nDesc0[children[node][child]]
		if doTimeTree:
			nodesToVisit.append((bestNode,0,vectUpUp,vectUpUpTime,dist[bestNode]+dist[node],bestLKdiff,0,removedPartialsRelative1,nDesc0ToAdd))
		else:
			nodesToVisit.append((bestNode,0,vectUpUp,dist[bestNode]+dist[node],bestLKdiff,0,removedPartialsRelative1,nDesc0ToAdd))
		originalBLens=(dist[node],dist[bestNode],removedBLen)
	else:
		# case node is root
		if children[bestNode]: # case there is only one sample outside of the subtree doesn't need to be considered
			child1=children[bestNode][0]
			child2=children[bestNode][1]
			vectUp1=probVect[child2]
			if mutations[child2]:
				vectUp1=passGenomeListThroughBranch(vectUp1,mutations[child2],dirIsUp=True)
			vectUp1=rootVector(vectUp1,dist[child2],(len(children[child2])==0 and len(minorSequences[child2])==0),tree,node,mutMatrixGlobalPassed=mutMatrixGlobalPassed,mutMatricesGlobal=mutMatricesGlobal)
			if mutations[child1]:
				removedPartialsRelative1=passGenomeListThroughBranch(bestRemovedPartials,mutations[child1])
				vectUp1=passGenomeListThroughBranch(vectUp1,mutations[child1],dirIsUp=False)
			else:
				removedPartialsRelative1=bestRemovedPartials
			nDesc0ToAdd=0
			if HnZ:
				if dist[child1]<effectivelyNon0BLen and dist[bestNode]<effectivelyNon0BLen:
					if dist[children[node][child]]>=effectivelyNon0BLen:
						nDesc0ToAdd=-1
					else:
						nDesc0ToAdd=-nDesc0[children[node][child]]
			if doTimeTree:
				vectUp1Time=rootVectorTime(probVectTime[child2],dist[child2],mutRate)
				nodesToVisit.append((child1,0,vectUp1,vectUp1Time,dist[child1],bestLKdiff,0,removedPartialsRelative1,nDesc0ToAdd))
			else:
				nodesToVisit.append((child1,0,vectUp1,dist[child1],bestLKdiff,0,removedPartialsRelative1,nDesc0ToAdd))
			vectUp2=probVect[child1]
			if mutations[child1]:
				vectUp2=passGenomeListThroughBranch(vectUp2,mutations[child1],dirIsUp=True)
			vectUp2=rootVector(vectUp2,dist[child1],(len(children[child1])==0 and len(minorSequences[child1])==0),tree,node,mutMatrixGlobalPassed=mutMatrixGlobalPassed,mutMatricesGlobal=mutMatricesGlobal)
			if mutations[child2]:
				removedPartialsRelative2=passGenomeListThroughBranch(bestRemovedPartials,mutations[child2])
				vectUp2=passGenomeListThroughBranch(vectUp2,mutations[child2],dirIsUp=False)
			else:
				removedPartialsRelative2=bestRemovedPartials
			nDesc0ToAdd=0
			if HnZ:
				if dist[child2]<effectivelyNon0BLen and dist[bestNode]<effectivelyNon0BLen:
					if dist[children[node][child]]>=effectivelyNon0BLen:
						nDesc0ToAdd=-1
					else:
						nDesc0ToAdd=-nDesc0[children[node][child]]
			if doTimeTree:
				vectUp2Time=rootVectorTime(probVectTime[child1],dist[child1],mutRate)
				nodesToVisit.append((child2,0,vectUp2,vectUp2Time,dist[child2],bestLKdiff,0,removedPartialsRelative2,nDesc0ToAdd))
			else:
				nodesToVisit.append((child2,0,vectUp2,dist[child2],bestLKdiff,0,removedPartialsRelative2,nDesc0ToAdd))
		originalBLens=(0.0,dist[bestNode],removedBLen)
	bestBranchLengths=originalBLens

	while nodesToVisit:
		nodeToVisitInfo=nodesToVisit.pop()
		if len(nodeToVisitInfo)==9:
			t1,direction,passedPartials,passedPartialsTime,distance,lastLK,failedPasses,removedPartialsRelative,nDesc0ToAdd=nodeToVisitInfo #modifiable
			if passedPartials!=None:
				needsUpdating=True
			else:
				needsUpdating=False
			needsUpdatingTime=True
		elif len(nodeToVisitInfo)==8:
			t1,direction,passedPartials,distance,lastLK,failedPasses,removedPartialsRelative,nDesc0ToAdd=nodeToVisitInfo #modifiable
			needsUpdating=True
			needsUpdatingTime=False
		else:
			needsUpdating=False
			needsUpdatingTime=False
			t1,direction,lastLK,failedPasses,removedPartialsRelative,nDesc0ToAdd=nodeToVisitInfo #modifiable

		if direction==0:
			#consider the case we are moving from a parent to a child
			if (not (up[t1]==node or up[t1]==None)) and (dist[t1]>effectivelyNon0BLen or doTimeTree or (up[up[t1]]==None)):
				if needsUpdating:
					isTip=(len(children[t1])==0) and (len(minorSequences[t1])==0)
					midTot=mergeVectors(passedPartials,distance/2,False,probVect[t1],distance/2,isTip,isUpDown=True   ,errorRateGlobalPassed=errorRateGlobalPassed,mutMatrixGlobalPassed=mutMatrixGlobalPassed,errorRatesGlobal=errorRatesGlobal,mutMatricesGlobal=mutMatricesGlobal,cumulativeRateGlobal=cumulativeRateGlobal,cumulativeErrorRateGlobal=cumulativeErrorRateGlobal)
					if midTot==None:
						continue
					if not areVectorsDifferent(midTot,probVectTotUp[t1]):
						needsUpdating=False
				else:
					midTot=probVectTotUp[t1]
					distance=dist[t1]
				if midTot==None:
					continue
				midProb=appendProbNode(midTot,removedPartialsRelative,isRemovedTip,removedBLen,   errorRateGlobalPassed=errorRateGlobalPassed,mutMatrixGlobalPassed=mutMatrixGlobalPassed,errorRatesGlobal=errorRatesGlobal,mutMatricesGlobal=mutMatricesGlobal,totErrorPassed=totErrorPassed)
				if doTimeTree:
					if needsUpdatingTime:
						midTotTime,midTotTimeCost=mergeVectorsTime(passedPartialsTime,distance/2,probVectTime[t1],distance/2,mutRate,isUpDown=True,returnLK=True)
						if isinstance(midTotTime, int):
							continue
						if not areVectorsDifferentTime(midTotTime,probVectTotUpTime[t1][0]):
							needsUpdatingTime=False
					else:
						midTotTime=probVectTotUpTime[t1][0]
						midTotTimeCost=probVectTotUpTime[t1][1]
					midProb+=midTotTimeCost
					midProb+=appendProbNodeTime(midTotTime,removedPartialsRelativeTime,mutRate,removedBLen)
				#correct new placement likelihood by the change in HnZ
				if HnZ:
					if (doTimeTree or up[up[t1]]==None) and distance<=effectivelyNon0BLen:
						parentNode0=t1
						while dist[parentNode0]<=effectivelyNon0BLen and up[parentNode0]!=None:
							parentNode0=up[parentNode0]
						if removedBLen>effectivelyNon0BLen:
							midProb+=getHnZ(nDesc0[parentNode0]+nDesc0ToAdd+1) - getHnZ(nDesc0[parentNode0]+nDesc0ToAdd)
						else:
							midProb+=getHnZ(nDesc0[children[node][child]]+nDesc0[parentNode0]+nDesc0ToAdd) - (getHnZ(nDesc0[children[node][child]])+getHnZ(nDesc0[parentNode0]+nDesc0ToAdd))
					else:
						if removedBLen>effectivelyNon0BLen:
							midProb+=getHnZ(2) - getHnZ(1)
						else:
							midProb+=getHnZ(nDesc0[children[node][child]]+1) - getHnZ(nDesc0[children[node][child]])
				if midProb>bestLKdiff-thresholdLogLKoptimizationTopology:
					#if needsUpdating, then add to the tuple also the information on the up and down genome lists to use to recalculate intermediate genome lists at varying branch lengths
					if needsUpdating:
						if needsUpdatingTime:
							bestNodes.append((t1,midProb,passedPartials,passedPartialsTime,probVect[t1],probVectTime[t1],distance,midTot,removedPartialsRelative))
						else:
							bestNodes.append((t1,midProb,passedPartials,probVect[t1],distance,midTot,removedPartialsRelative))
					else:
						if needsUpdatingTime:
							bestNodes.append((t1,midProb,None,passedPartialsTime,None,probVectTime[t1],distance,None,removedPartialsRelative))
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
					vectUpRight=mergeVectors(passedPartials,distance,False,otherProbVect,dist[otherChild],isTip,isUpDown=True   ,errorRateGlobalPassed=errorRateGlobalPassed,mutMatrixGlobalPassed=mutMatrixGlobalPassed,errorRatesGlobal=errorRatesGlobal,mutMatricesGlobal=mutMatricesGlobal,cumulativeRateGlobal=cumulativeRateGlobal,cumulativeErrorRateGlobal=cumulativeErrorRateGlobal)
				else:
					vectUpRight=probVectUpRight[t1]
				if vectUpRight!=None:
					removedPartialsRelative1=removedPartialsRelative
					if mutations[child1]:
						removedPartialsRelative1=passGenomeListThroughBranch(removedPartialsRelative,mutations[child1])
					if needsUpdatingTime:
						otherChild=children[t1][1]
						otherProbVectTime=probVectTime[otherChild]
						vectUpRightTime=mergeVectorsTime(passedPartialsTime,distance,otherProbVectTime,dist[otherChild],mutRate,isUpDown=True)
						if isinstance(vectUpRightTime, int):
							continue
					nDesc0ToAddToPass=0
					if nDesc0ToAdd:
						if dist[child1]<effectivelyNon0BLen:
							nDesc0ToAddToPass=nDesc0ToAdd
					if needsUpdating:			
						if mutations[child1]:
							vectUpRight=passGenomeListThroughBranch(vectUpRight,mutations[child1])
						if needsUpdatingTime:
							nodesToVisit.append((child1,0,vectUpRight,vectUpRightTime,dist[child1],midProb,failedPasses,removedPartialsRelative1,nDesc0ToAddToPass))
						else:
							nodesToVisit.append((child1,0,vectUpRight,dist[child1],midProb,failedPasses,removedPartialsRelative1,nDesc0ToAddToPass))
					else:
						if needsUpdatingTime:
							nodesToVisit.append((child1,0,None,vectUpRightTime,dist[child1],midProb,failedPasses,removedPartialsRelative1,nDesc0ToAddToPass))
						else:
							nodesToVisit.append((child1,0,midProb,failedPasses,removedPartialsRelative1,nDesc0ToAddToPass))

				child1=children[t1][1]
				if needsUpdating:
					otherChild=children[t1][0]
					isTip=(len(children[otherChild])==0) and (len(minorSequences[otherChild])==0)
					otherProbVect=probVect[otherChild]
					if mutations[otherChild]:
						otherProbVect=passGenomeListThroughBranch(otherProbVect,mutations[otherChild],dirIsUp=True)
					vectUpLeft=mergeVectors(passedPartials,distance,False,otherProbVect,dist[otherChild],isTip,isUpDown=True   ,errorRateGlobalPassed=errorRateGlobalPassed,mutMatrixGlobalPassed=mutMatrixGlobalPassed,errorRatesGlobal=errorRatesGlobal,mutMatricesGlobal=mutMatricesGlobal,cumulativeRateGlobal=cumulativeRateGlobal,cumulativeErrorRateGlobal=cumulativeErrorRateGlobal)
				else:
					vectUpLeft=probVectUpLeft[t1]
				if vectUpLeft!=None:
					removedPartialsRelative1=removedPartialsRelative
					if mutations[child1]:
						removedPartialsRelative1=passGenomeListThroughBranch(removedPartialsRelative,mutations[child1])
					if needsUpdatingTime:
						otherChild=children[t1][0]
						otherProbVectTime=probVectTime[otherChild]
						vectUpLeftTime=mergeVectorsTime(passedPartialsTime,distance,otherProbVectTime,dist[otherChild],mutRate,isUpDown=True)
						if isinstance(vectUpLeftTime, int):
							continue
					nDesc0ToAddToPass=0
					if nDesc0ToAdd:
						if dist[child1]<effectivelyNon0BLen:
							nDesc0ToAddToPass=nDesc0ToAdd
					if needsUpdating:
						if mutations[child1]:
							vectUpLeft=passGenomeListThroughBranch(vectUpLeft,mutations[child1])
						if needsUpdatingTime:
							nodesToVisit.append((child1,0,vectUpLeft,vectUpLeftTime,dist[child1],midProb,failedPasses,removedPartialsRelative1,nDesc0ToAddToPass))
						else:
							nodesToVisit.append((child1,0,vectUpLeft,dist[child1],midProb,failedPasses,removedPartialsRelative1,nDesc0ToAddToPass))
					else:
						if needsUpdatingTime:
							nodesToVisit.append((child1,0,None,vectUpLeftTime,dist[child1],midProb,failedPasses,removedPartialsRelative1,nDesc0ToAddToPass))
						else:
							nodesToVisit.append((child1,0,midProb,failedPasses,removedPartialsRelative1,nDesc0ToAddToPass))

		else: #case when crawling up from child to parent
			otherChild=children[t1][2-direction]
			midBottom=None
			if up[t1]!=None and (dist[t1]>effectivelyNon0BLen or doTimeTree or up[up[t1]]==None): #try appending mid-branch
				if needsUpdating:
					otherProbVect=probVect[otherChild]
					if mutations[otherChild]:
						otherProbVect=passGenomeListThroughBranch(otherProbVect,mutations[otherChild],dirIsUp=True)
					isTip=(len(children[otherChild])==0) and (len(minorSequences[otherChild])==0)
					midBottom=mergeVectors(passedPartials,distance,False,otherProbVect,dist[otherChild],isTip   ,errorRateGlobalPassed=errorRateGlobalPassed,mutMatrixGlobalPassed=mutMatrixGlobalPassed,errorRatesGlobal=errorRatesGlobal,mutMatricesGlobal=mutMatricesGlobal,cumulativeRateGlobal=cumulativeRateGlobal,cumulativeErrorRateGlobal=cumulativeErrorRateGlobal)
					if midBottom==None:
						continue
					if t1==children[up[t1]][0]:
						vectUp=probVectUpRight[up[t1]]
					else:
						vectUp=probVectUpLeft[up[t1]]
					if mutations[t1]:
						vectUp=passGenomeListThroughBranch(vectUp,mutations[t1])
					midTot=mergeVectors(vectUp,dist[t1]/2,False,midBottom,dist[t1]/2,False,isUpDown=True   ,errorRateGlobalPassed=errorRateGlobalPassed,mutMatrixGlobalPassed=mutMatrixGlobalPassed,errorRatesGlobal=errorRatesGlobal,mutMatricesGlobal=mutMatricesGlobal,cumulativeRateGlobal=cumulativeRateGlobal,cumulativeErrorRateGlobal=cumulativeErrorRateGlobal)
					if not probVectTotUp[t1]:
						print("Node has no probVectTotUp "+str(name[t1])+" with dist "+str(dist[t1])+" calculating new one")
						probVectTotUp[t1]=mergeVectors(vectUp,dist[t1]/2,False,probVect[t1],dist[t1]/2,False,isUpDown=True   ,errorRateGlobalPassed=errorRateGlobalPassed,mutMatrixGlobalPassed=mutMatrixGlobalPassed,errorRatesGlobal=errorRatesGlobal,mutMatricesGlobal=mutMatricesGlobal,cumulativeRateGlobal=cumulativeRateGlobal,cumulativeErrorRateGlobal=cumulativeErrorRateGlobal)
					if midTot==None:
						continue
					if not areVectorsDifferent(midTot,probVectTotUp[t1]):
						needsUpdating=False
				else:
					midTot=probVectTotUp[t1]
				if midTot==None:
					continue
				midProb=appendProbNode(midTot,removedPartialsRelative,isRemovedTip,removedBLen,   errorRateGlobalPassed=errorRateGlobalPassed,mutMatrixGlobalPassed=mutMatrixGlobalPassed,errorRatesGlobal=errorRatesGlobal,mutMatricesGlobal=mutMatricesGlobal,totErrorPassed=totErrorPassed)
				if doTimeTree:
					if needsUpdatingTime:
						otherProbVectTime=probVectTime[otherChild]
						midBottomTime=mergeVectorsTime(passedPartialsTime,distance,otherProbVectTime,dist[otherChild],mutRate)
						if t1==children[up[t1]][0]:
							vectUpTime=probVectUpRightTime[up[t1]]
						else:
							vectUpTime=probVectUpLeftTime[up[t1]]
						midTotTime,midTotTimeCost=mergeVectorsTime(vectUpTime,dist[t1]/2,midBottomTime,dist[t1]/2,mutRate,isUpDown=True,returnLK=True)
						if isinstance(midTotTime, int):
							midProb=float("-inf")
						if not areVectorsDifferentTime(midTotTime,probVectTotUpTime[t1][0]):
							needsUpdatingTime=False
					else:
						midTotTime=probVectTotUpTime[t1][0]
						midTotTimeCost=probVectTotUpTime[t1][1]
					if not isinstance(midTotTime, int):
						midProb+=appendProbNodeTime(midTotTime,removedPartialsRelativeTime,mutRate,removedBLen)
						midProb+=midTotTimeCost
				#correct new placement likelihood by the change in HnZ
				if HnZ:
					if (doTimeTree or up[up[t1]]==None) and distance<=effectivelyNon0BLen:
						parentNode0=t1
						while dist[parentNode0]<=effectivelyNon0BLen and up[parentNode0]!=None:
							parentNode0=up[parentNode0]
						if removedBLen>effectivelyNon0BLen:
							midProb+=getHnZ(nDesc0[parentNode0]+nDesc0ToAdd+1) - getHnZ(nDesc0[parentNode0]+nDesc0ToAdd)
						else:
							midProb+=getHnZ(nDesc0[children[node][child]]+nDesc0[parentNode0]+nDesc0ToAdd) - (getHnZ(nDesc0[children[node][child]])+getHnZ(nDesc0[parentNode0]+nDesc0ToAdd))
					else:
						if removedBLen>effectivelyNon0BLen:
							midProb+=getHnZ(2) - getHnZ(1)
						else:
							midProb+=getHnZ(nDesc0[children[node][child]]+1) - getHnZ(nDesc0[children[node][child]])
				if midProb>=(bestLKdiff-thresholdLogLKoptimizationTopology):
					#if needsUpdating, then add to the tuple also the information on the up and down genome lists to use to recalculate intermediate genome lists at varying branch lengths
					if needsUpdating:
						if needsUpdatingTime:
							bestNodes.append((t1,midProb,vectUp,vectUpTime,midBottom,midBottomTime,dist[t1],midTot,removedPartialsRelative))
						else:
							bestNodes.append((t1,midProb,vectUp,midBottom,dist[t1],midTot,removedPartialsRelative))
					else:
						if needsUpdatingTime:
							bestNodes.append((t1,midProb,None,vectUpTime,None,midBottomTime,dist[t1],None,removedPartialsRelative))
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
						if needsUpdatingTime:
							vectUpUpTime=probVectUpRightTime[up[t1]]
					else:
						upChild=1
						if needsUpdating:
							vectUpUp=probVectUpLeft[up[t1]]
						if needsUpdatingTime:
							vectUpUpTime=probVectUpLeftTime[up[t1]]
					if needsUpdating:
						if mutations[t1]:
							vectUpUp=passGenomeListThroughBranch(vectUpUp,mutations[t1])
						vectUp=mergeVectors(vectUpUp,dist[t1],False,passedPartials,distance,False,isUpDown=True   ,errorRateGlobalPassed=errorRateGlobalPassed,mutMatrixGlobalPassed=mutMatrixGlobalPassed,errorRatesGlobal=errorRatesGlobal,mutMatricesGlobal=mutMatricesGlobal,cumulativeRateGlobal=cumulativeRateGlobal,cumulativeErrorRateGlobal=cumulativeErrorRateGlobal)
					else:
						if direction==1:
							vectUp=probVectUpLeft[t1]
						else:
							vectUp=probVectUpRight[t1]
					if needsUpdatingTime:
						vectUpTime=mergeVectorsTime(vectUpUpTime,dist[t1],passedPartialsTime,distance,mutRate,isUpDown=True)
						if isinstance(vectUpTime, int):
							continue

					if vectUp==None:
						continue
					else:
						removedPartialsRelative1=removedPartialsRelative
						if mutations[otherChild]:
							removedPartialsRelative1=passGenomeListThroughBranch(removedPartialsRelative,mutations[otherChild])
						nDesc0ToAddToPass=0
						if nDesc0ToAdd:
							if dist[otherChild]<effectivelyNon0BLen:
								nDesc0ToAddToPass=nDesc0ToAdd
						if needsUpdating:
							if mutations[otherChild]:
								vectUp=passGenomeListThroughBranch(vectUp,mutations[otherChild])
							if needsUpdatingTime:
								nodesToVisit.append((otherChild,0,vectUp,vectUpTime,dist[otherChild],midProb,failedPasses,removedPartialsRelative1,nDesc0ToAddToPass))
							else:
								nodesToVisit.append((otherChild,0,vectUp,dist[otherChild],midProb,failedPasses,removedPartialsRelative1,nDesc0ToAddToPass))
						else:
							if needsUpdatingTime:
								nodesToVisit.append((otherChild,0,None,vectUpTime,dist[otherChild],midProb,failedPasses,removedPartialsRelative1,nDesc0ToAddToPass))
							else:
								nodesToVisit.append((otherChild,0,midProb,failedPasses,removedPartialsRelative1,nDesc0ToAddToPass))
					#now pass the crawling up to the parent node
					if needsUpdating:
						if midBottom==None:
							otherProbVect=probVect[otherChild]
							if mutations[otherChild]:
								otherProbVect=passGenomeListThroughBranch(otherProbVect,mutations[otherChild],dirIsUp=True)
							isTip=(len(children[otherChild])==0) and (len(minorSequences[otherChild])==0)
							midBottom=mergeVectors(passedPartials,distance,False,otherProbVect,dist[otherChild],isTip   ,errorRateGlobalPassed=errorRateGlobalPassed,mutMatrixGlobalPassed=mutMatrixGlobalPassed,errorRatesGlobal=errorRatesGlobal,mutMatricesGlobal=mutMatricesGlobal,cumulativeRateGlobal=cumulativeRateGlobal,cumulativeErrorRateGlobal=cumulativeErrorRateGlobal)
							if midBottom==None:
								continue
					removedPartialsRelative1=removedPartialsRelative
					if mutations[t1]:
						removedPartialsRelative1=passGenomeListThroughBranch(removedPartialsRelative,mutations[t1],dirIsUp=True)
					nDesc0ToAddToPass=0
					if nDesc0ToAdd:
						if dist[t1]<effectivelyNon0BLen:
							nDesc0ToAddToPass=nDesc0ToAdd
					if needsUpdating:
						if mutations[t1]:
							midBottom=passGenomeListThroughBranch(midBottom,mutations[t1],dirIsUp=True)
						if needsUpdatingTime:
							nodesToVisit.append((up[t1],upChild+1,midBottom,midBottomTime,dist[t1],midProb,failedPasses,removedPartialsRelative1,nDesc0ToAddToPass))
						else:
							nodesToVisit.append((up[t1],upChild+1,midBottom,dist[t1],midProb,failedPasses,removedPartialsRelative1,nDesc0ToAddToPass))
					else:
						if needsUpdatingTime:
							nodesToVisit.append((up[t1],upChild+1,None,midBottomTime,dist[t1],midProb,failedPasses,removedPartialsRelative1,nDesc0ToAddToPass))
						else:
							nodesToVisit.append((up[t1],upChild+1,midProb,failedPasses,removedPartialsRelative1,nDesc0ToAddToPass))
				#now consider case of root node
				else:
					if needsUpdating:
						vectUp=rootVector(passedPartials,distance,False,tree,t1,mutMatrixGlobalPassed=mutMatrixGlobalPassed,mutMatricesGlobal=mutMatricesGlobal)
						if mutations[otherChild]:
							vectUp=passGenomeListThroughBranch(vectUp,mutations[otherChild])
					removedPartialsRelative1=removedPartialsRelative
					if mutations[otherChild]:
						removedPartialsRelative1=passGenomeListThroughBranch(removedPartialsRelative,mutations[otherChild])
					nDesc0ToAddToPass=0
					if nDesc0ToAdd:
						if dist[otherChild]<effectivelyNon0BLen:
							nDesc0ToAddToPass=nDesc0ToAdd
					if needsUpdating:
						if needsUpdatingTime:
							vectUpTime=rootVectorTime(passedPartialsTime,distance,mutRate)
							nodesToVisit.append((otherChild,0,vectUp,vectUpTime,dist[otherChild],midProb,failedPasses,removedPartialsRelative1,nDesc0ToAddToPass))
						else:
							nodesToVisit.append((otherChild,0,vectUp,dist[otherChild],midProb,failedPasses,removedPartialsRelative1,nDesc0ToAddToPass))
					else:
						if needsUpdatingTime:
							vectUpTime=rootVectorTime(passedPartialsTime,distance,mutRate)
							nodesToVisit.append((otherChild,0,None,vectUpTime,dist[otherChild],midProb,failedPasses,removedPartialsRelative1,nDesc0ToAddToPass))
						else:
							nodesToVisit.append((otherChild,0,midProb,failedPasses,removedPartialsRelative1,nDesc0ToAddToPass))

	#Initial exploration is finished.
	#Now, for each branch within threshold likelihood distance from the best found, optimize branch lengths.
	#Use optimized scores to select final best branch
	bestScore=bestLKdiff
	compensanteForBranchLengthChange=True
	if not bestNodes:
		totalTimeFindingParent[0]+=(time()-timeStartParentTopology)
		return originalPlacement, originalLK ,originalBLens, [], 1.0, originalRemovedPartialsRelative
	#if aBayesPlusOn and sequentialSearch:
	if aBayesPlusOn:
		if networkOutput:
			listOfProbableNodes=[]
		listofLKcosts=[]
		rootAlreadyConsidered=False
		if up[originalParent0]==None:
			rootAlreadyConsidered=True
		if up[node]==None or (up[up[node]]==None and dist[children[node][1-child]]>effectivelyNon0BLen):
			rootAlreadyConsidered=True
	if doTimeTree:
		topNodes={}
		originalNode=children[node][1-child]
		if dist[children[node][1-child]]<=effectivelyNon0BLen:
			originalNode=node
		if dist[node]<=effectivelyNon0BLen:
			originalNode=originalParent0
		if up[node]!=None and up[up[node]]==None and dist[children[node][1-child]]>effectivelyNon0BLen:
			originalNode=up[node]
		topNodes[originalNode]=originalLK
	for nodePair in bestNodes:
		score=nodePair[1]
		if score>=bestLKdiff-thresholdLogLKoptimizationTopology:
			t1=nodePair[0]
			if len(nodePair)==3 or nodePair[2]==None:
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
			elif len(nodePair)==7:
				upVect=nodePair[2]
				downVect=nodePair[3]
				distance=nodePair[4]
				midTot=nodePair[5]
			else:
				upVect=nodePair[2]
				downVect=nodePair[4]
				distance=nodePair[6]
				midTot=nodePair[7]
			if doTimeTree:
				if len(nodePair)<9:
					if t1==children[up[t1]][0]:
						upVectTime=probVectUpRightTime[up[t1]]
					else:
						upVectTime=probVectUpLeftTime[up[t1]]
					downVectTime=probVectTime[t1]
				else:
					upVectTime=nodePair[3]
					downVectTime=nodePair[5]

			removedPartials=nodePair[-1]
			bestAppendingLength=estimateBranchLengthWithDerivative(midTot,removedPartials,fromTipC=isRemovedTip   ,errorRateGlobalPassed=errorRateGlobalPassed,mutMatrixGlobalPassed=mutMatrixGlobalPassed,errorRatesGlobal=errorRatesGlobal,mutMatricesGlobal=mutMatricesGlobal,cumulativeRateGlobal=cumulativeRateGlobal)
			#now optimize appending location
			fromTip1=(len(children[t1])==0) and (len(minorSequences[t1])==0)
			midLowerVector=mergeVectors(downVect,distance/2,fromTip1,removedPartials,bestAppendingLength,isRemovedTip   ,errorRateGlobalPassed=errorRateGlobalPassed,mutMatrixGlobalPassed=mutMatrixGlobalPassed,errorRatesGlobal=errorRatesGlobal,mutMatricesGlobal=mutMatricesGlobal,cumulativeRateGlobal=cumulativeRateGlobal,cumulativeErrorRateGlobal=cumulativeErrorRateGlobal)
			if midLowerVector==None:
				print("Issue with merging vectors in findBestParentTopology, midLowerVector is None")
				print(t1)
				print(up[t1])
				print(up[up[t1]])
				print("UpVects")
				print(upVect)
				if t1==children[up[t1]][0]:
					print(probVectUpRight[up[t1]])
				else:
					print(probVectUpLeft[up[t1]])
				print(distance)
				print("removedPartials")
				print(removedPartials)
				print(bestAppendingLength)
				print("downVects")
				print(downVect)
				print(probVect[t1])
				print("probVectTotUp")
				print(probVectTotUp[t1])
			bestTopLength=estimateBranchLengthWithDerivative(upVect,midLowerVector,fromTipC=False   ,errorRateGlobalPassed=errorRateGlobalPassed,mutMatrixGlobalPassed=mutMatrixGlobalPassed,errorRatesGlobal=errorRatesGlobal,mutMatricesGlobal=mutMatricesGlobal,cumulativeRateGlobal=cumulativeRateGlobal)
			midTopVector=mergeVectors(upVect,bestTopLength,False,removedPartials,bestAppendingLength,isRemovedTip,isUpDown=True   ,errorRateGlobalPassed=errorRateGlobalPassed,mutMatrixGlobalPassed=mutMatrixGlobalPassed,errorRatesGlobal=errorRatesGlobal,mutMatricesGlobal=mutMatricesGlobal,cumulativeRateGlobal=cumulativeRateGlobal,cumulativeErrorRateGlobal=cumulativeErrorRateGlobal)
			if midTopVector==None:
				print("Inconsistent branch length found")
				print(upVect)
				bestTopLength=defaultBLen*0.1
				midTopVector=mergeVectors(upVect,bestTopLength,False,removedPartials,bestAppendingLength,isRemovedTip,isUpDown=True   ,errorRateGlobalPassed=errorRateGlobalPassed,mutMatrixGlobalPassed=mutMatrixGlobalPassed,errorRatesGlobal=errorRatesGlobal,mutMatricesGlobal=mutMatricesGlobal,cumulativeRateGlobal=cumulativeRateGlobal,cumulativeErrorRateGlobal=cumulativeErrorRateGlobal)
			bestBottomLength=estimateBranchLengthWithDerivative(midTopVector,downVect,fromTipC=fromTip1   ,errorRateGlobalPassed=errorRateGlobalPassed,mutMatrixGlobalPassed=mutMatrixGlobalPassed,errorRatesGlobal=errorRatesGlobal,mutMatricesGlobal=mutMatricesGlobal,cumulativeRateGlobal=cumulativeRateGlobal)
			newMidVector=mergeVectors(upVect,bestTopLength,False,downVect,bestBottomLength,fromTip1,isUpDown=True   ,errorRateGlobalPassed=errorRateGlobalPassed,mutMatrixGlobalPassed=mutMatrixGlobalPassed,errorRatesGlobal=errorRatesGlobal,mutMatricesGlobal=mutMatricesGlobal,cumulativeRateGlobal=cumulativeRateGlobal,cumulativeErrorRateGlobal=cumulativeErrorRateGlobal)
			appendingCost=appendProbNode(newMidVector,removedPartials,isRemovedTip,bestAppendingLength,   errorRateGlobalPassed=errorRateGlobalPassed,mutMatrixGlobalPassed=mutMatrixGlobalPassed,errorRatesGlobal=errorRatesGlobal,mutMatricesGlobal=mutMatricesGlobal,totErrorPassed=totErrorPassed)
			if doTimeTree:
				newMidVectorTime,appendingCostTime=mergeVectorsTime(upVectTime,bestTopLength,downVectTime,bestBottomLength,mutRate,isUpDown=True,returnLK=True)
				appendingCostTime-=appendProbNodeTime(upVectTime,downVectTime,mutRate,distance)
				if isinstance(newMidVectorTime, int):
					appendingCost=float("-inf")
				else:
					appendingCostTime+=appendProbNodeTime(newMidVectorTime,removedPartialsRelativeTime,mutRate,bestAppendingLength)
					appendingCost+=appendingCostTime
			if compensanteForBranchLengthChange: #if wanted, do a more thorough examination of the appending cost, taking into account a possible change in branch length of the branch on which to be appended.
				initialCost=appendProbNode(upVect,downVect,fromTip1,distance,   errorRateGlobalPassed=errorRateGlobalPassed,mutMatrixGlobalPassed=mutMatrixGlobalPassed,errorRatesGlobal=errorRatesGlobal,mutMatricesGlobal=mutMatricesGlobal,totErrorPassed=totErrorPassed)
				newPartialCost=appendProbNode(upVect,downVect,fromTip1,bestBottomLength+bestTopLength,   errorRateGlobalPassed=errorRateGlobalPassed,mutMatrixGlobalPassed=mutMatrixGlobalPassed,errorRatesGlobal=errorRatesGlobal,mutMatricesGlobal=mutMatricesGlobal,totErrorPassed=totErrorPassed)
				optimizedScore=appendingCost+newPartialCost-initialCost
			else:
				optimizedScore=appendingCost
			
			#correct new placement likelihood by the change in HnZ
			#branch lengths are bestTopLength, bestBottomLength, bestAppendingLength, 
			if HnZ:
				belowT1=False
				originalParentNode0=node
				if originalParentNode0==t1:
					belowT1=True
				while (dist[originalParentNode0]<=effectivelyNon0BLen) and up[originalParentNode0]!=None:
					originalParentNode0=up[originalParentNode0]
					if originalParentNode0==t1:
						belowT1=True
				parentNode0=up[t1]
				while (dist[parentNode0]<=effectivelyNon0BLen) and up[parentNode0]!=None:
					parentNode0=up[parentNode0]
				toBeAddedToCompensate=0
				if parentNode0==originalParentNode0:
					if dist[children[node][child]]:
						toBeAddedToCompensate=-1
					else:
						toBeAddedToCompensate=-nDesc0[children[node][child]]
				toBeAddedToCompensateT1=0
				if belowT1:
					if dist[children[node][child]]:
						toBeAddedToCompensateT1=-1
					else:
						toBeAddedToCompensateT1=-nDesc0[children[node][child]]
				if bestTopLength>effectivelyNon0BLen and bestBottomLength>effectivelyNon0BLen:
					if bestAppendingLength>effectivelyNon0BLen:
						addendum0=getHnZ(2) - getHnZ(1)
					else:
						addendum0=getHnZ(nDesc0[children[node][child]]+1) - getHnZ(nDesc0[children[node][child]])
					if dist[t1]<=effectivelyNon0BLen:
						addendum0+=getHnZ(nDesc0[parentNode0]+1-toBeAddedToCompensateT1+toBeAddedToCompensate-nDesc0[t1])+getHnZ(nDesc0[t1]+toBeAddedToCompensateT1)-getHnZ(nDesc0[parentNode0]+toBeAddedToCompensate)
				elif bestBottomLength>effectivelyNon0BLen :
					if parentNode0==originalParent0 and (not doTimeTree):
						addendum0=float("-inf")
					else:
						if bestAppendingLength>effectivelyNon0BLen:
							if dist[t1]<=effectivelyNon0BLen:
								addendum0=getHnZ(nDesc0[parentNode0]+toBeAddedToCompensate+2-toBeAddedToCompensateT1-nDesc0[t1])+getHnZ(nDesc0[t1]+toBeAddedToCompensateT1) - getHnZ(nDesc0[parentNode0]+toBeAddedToCompensate)
							else:
								addendum0=getHnZ(nDesc0[parentNode0]+toBeAddedToCompensate+1) - getHnZ(nDesc0[parentNode0]+toBeAddedToCompensate)
						else:
							if dist[t1]<=effectivelyNon0BLen:
								addendum0=getHnZ(nDesc0[parentNode0]+toBeAddedToCompensate+1-toBeAddedToCompensateT1+nDesc0[children[node][child]]-nDesc0[t1])+getHnZ(nDesc0[t1]+toBeAddedToCompensateT1) - (getHnZ(nDesc0[children[node][child]]) + getHnZ(nDesc0[parentNode0]+toBeAddedToCompensate))
							else:
								addendum0=getHnZ(nDesc0[parentNode0]+toBeAddedToCompensate+nDesc0[children[node][child]]) - ( getHnZ(nDesc0[children[node][child]]) + getHnZ(nDesc0[parentNode0]+toBeAddedToCompensate) )

				elif bestTopLength>effectivelyNon0BLen:
					if t1==originalParent0 and (not doTimeTree):
						addendum0=float("-inf")
					else:
						if dist[t1]<=effectivelyNon0BLen:
							if bestAppendingLength>effectivelyNon0BLen:
								addendum0=getHnZ(nDesc0[t1]+toBeAddedToCompensateT1+1) + getHnZ(nDesc0[parentNode0]+1+toBeAddedToCompensate-toBeAddedToCompensateT1-nDesc0[t1])-getHnZ(nDesc0[parentNode0]+toBeAddedToCompensate)
							else:
								addendum0=getHnZ(nDesc0[t1]+toBeAddedToCompensateT1+nDesc0[children[node][child]]) + getHnZ(nDesc0[parentNode0]+1+toBeAddedToCompensate-toBeAddedToCompensateT1-nDesc0[t1]) - ( getHnZ(nDesc0[children[node][child]]) + getHnZ(nDesc0[parentNode0]+toBeAddedToCompensate) )
						else:
							if bestAppendingLength>effectivelyNon0BLen:
								addendum0=getHnZ(nDesc0[t1]+toBeAddedToCompensateT1+1) - getHnZ(nDesc0[t1]+toBeAddedToCompensateT1)
							else:
								addendum0=getHnZ(nDesc0[t1]+toBeAddedToCompensateT1+nDesc0[children[node][child]]) - ( getHnZ(nDesc0[children[node][child]]) + getHnZ(nDesc0[t1]+toBeAddedToCompensateT1) )
				else:
					if (parentNode0==originalParent0 or t1==originalParent0)  and (not doTimeTree):
						addendum0=float("-inf")
					else:
						if dist[t1]<=effectivelyNon0BLen:
							if bestAppendingLength>effectivelyNon0BLen:
								addendum0=getHnZ(nDesc0[parentNode0]+toBeAddedToCompensate+1) - getHnZ(nDesc0[parentNode0]+toBeAddedToCompensate)
							else:
								addendum0=getHnZ(nDesc0[parentNode0]+toBeAddedToCompensate+nDesc0[children[node][child]]) - ( getHnZ(nDesc0[children[node][child]]) + getHnZ(nDesc0[parentNode0]+toBeAddedToCompensate) )
						else:
							if bestAppendingLength>effectivelyNon0BLen:
								addendum0=getHnZ(nDesc0[parentNode0]+toBeAddedToCompensate+nDesc0[t1]+toBeAddedToCompensateT1+1) - (getHnZ(nDesc0[parentNode0]+toBeAddedToCompensate) + getHnZ(nDesc0[t1]+toBeAddedToCompensateT1))
							else:
								addendum0=getHnZ(nDesc0[parentNode0]+toBeAddedToCompensate+nDesc0[t1]+toBeAddedToCompensateT1+nDesc0[children[node][child]]) - ( getHnZ(nDesc0[children[node][child]]) + getHnZ(nDesc0[parentNode0]+toBeAddedToCompensate) + getHnZ(nDesc0[t1]+toBeAddedToCompensateT1) )
				optimizedScore+=addendum0

				#try also placing with a 0 bottom length in case it increase HnZ score.
				if bestBottomLength>effectivelyNon0BLen and dist[t1]>effectivelyNon0BLen:
					altNewMidVector=mergeVectors(upVect,bestTopLength+bestBottomLength,False,downVect,0.0,fromTip1,isUpDown=True   ,errorRateGlobalPassed=errorRateGlobalPassed,mutMatrixGlobalPassed=mutMatrixGlobalPassed,errorRatesGlobal=errorRatesGlobal,mutMatricesGlobal=mutMatricesGlobal,cumulativeRateGlobal=cumulativeRateGlobal,cumulativeErrorRateGlobal=cumulativeErrorRateGlobal)
					altAppendingCost=appendProbNode(altNewMidVector,removedPartials,isRemovedTip,bestAppendingLength,   errorRateGlobalPassed=errorRateGlobalPassed,mutMatrixGlobalPassed=mutMatrixGlobalPassed,errorRatesGlobal=errorRatesGlobal,mutMatricesGlobal=mutMatricesGlobal,totErrorPassed=totErrorPassed)
					if doTimeTree:
						altNewMidVectorTime,altAppendingCostTime=mergeVectorsTime(upVectTime,bestTopLength+bestBottomLength,downVectTime,0.0,mutRate,isUpDown=True,returnLK=True)
						altAppendingCostTime-=appendProbNodeTime(upVectTime,downVectTime,mutRate,distance)
						if isinstance(altNewMidVectorTime, int):
							altAppendingCost=float("-inf")
						else:
							altAppendingCostTime+=appendProbNodeTime(altNewMidVectorTime,removedPartialsRelativeTime,mutRate,bestAppendingLength)
							altAppendingCost+=altAppendingCostTime
					if compensanteForBranchLengthChange: #if wanted, do a more thorough examination of the appending cost, taking into account a possible change in branch length of the branch on which to be appended.
						initialCost=appendProbNode(upVect,downVect,fromTip1,distance,   errorRateGlobalPassed=errorRateGlobalPassed,mutMatrixGlobalPassed=mutMatrixGlobalPassed,errorRatesGlobal=errorRatesGlobal,mutMatricesGlobal=mutMatricesGlobal,totErrorPassed=totErrorPassed)
						newPartialCost=appendProbNode(upVect,downVect,fromTip1,bestBottomLength+bestTopLength,   errorRateGlobalPassed=errorRateGlobalPassed,mutMatrixGlobalPassed=mutMatrixGlobalPassed,errorRatesGlobal=errorRatesGlobal,mutMatricesGlobal=mutMatricesGlobal,totErrorPassed=totErrorPassed)
						altOptimizedScore=altAppendingCost+newPartialCost-initialCost
					else:
						altOptimizedScore=altAppendingCost
					if (bestTopLength+bestBottomLength)>effectivelyNon0BLen:
						if t1==originalParent0 and (not doTimeTree):
							addendum0=float("-inf")
						else:
							if bestAppendingLength>effectivelyNon0BLen:
								addendum0=getHnZ(nDesc0[t1]+toBeAddedToCompensateT1+1) - getHnZ(nDesc0[t1]+toBeAddedToCompensateT1)
							else:
								addendum0=getHnZ(nDesc0[t1]+toBeAddedToCompensateT1+nDesc0[children[node][child]]) - ( getHnZ(nDesc0[children[node][child]]) + getHnZ(nDesc0[t1]+toBeAddedToCompensateT1) )
					else:
						if (parentNode0==originalParent0 or t1==originalParent0)  and (not doTimeTree):
							addendum0=float("-inf")
						else:
							if bestAppendingLength>effectivelyNon0BLen:
								addendum0=getHnZ(nDesc0[parentNode0]+toBeAddedToCompensate+nDesc0[t1]+toBeAddedToCompensateT1+1) - (getHnZ(nDesc0[parentNode0]+toBeAddedToCompensate) + getHnZ(nDesc0[t1]+toBeAddedToCompensateT1))
							else:
								addendum0=getHnZ(nDesc0[parentNode0]+toBeAddedToCompensate+nDesc0[t1]+toBeAddedToCompensateT1+nDesc0[children[node][child]]) - ( getHnZ(nDesc0[children[node][child]]) + getHnZ(nDesc0[parentNode0]+toBeAddedToCompensate) + getHnZ(nDesc0[t1]+toBeAddedToCompensateT1) )
					altOptimizedScore+=addendum0

					if altOptimizedScore>optimizedScore:
						optimizedScore=altOptimizedScore
						bestTopLength=bestTopLength+bestBottomLength
						bestBottomLength=0.0

			if optimizedScore>=bestScore:
				bestNode=t1
				bestScore=optimizedScore
				bestBranchLengths=(bestTopLength,bestBottomLength,bestAppendingLength)
				bestRemovedPartials=removedPartials

			#if aBayesPlusOn and sequentialSearch:
			if aBayesPlusOn:
				#check that the placement location is effectively different from the original node
				if doTimeTree:
					if bestTopLength<=effectivelyNon0BLen:
						topNode=up[t1]
						while (dist[topNode]<=effectivelyNon0BLen) and (up[topNode]!=None):
							topNode=up[topNode]
					else:
						topNode=t1
					if (up[up[t1]]==None and bestBottomLength>effectivelyNon0BLen):
						topNode=up[t1]
					if (up[node]==None and up[topNode]==node):
						topNode=node
					if (topNode in topNodes):
						if (optimizedScore>topNodes[topNode]):
							topNodes[topNode]=optimizedScore
					else:
						topNodes[topNode]=optimizedScore

				else:
					differentNode=True
					if t1==node:
						differentNode=False
					elif t1==children[node][1-child]:
						if dist[node]>=effectivelyNon0BLen or bestTopLength<=effectivelyNon0BLen:
							differentNode=False
					if (bestBottomLength<=effectivelyNon0BLen):
						if t1==originalParent0:
							differentNode=False
					#check that placement is not redundant
					if (bestTopLength<=effectivelyNon0BLen):
						differentNode=False
					if dist[t1]<=effectivelyNon0BLen and up[up[t1]]!=None:
						differentNode=False
					if (not rootAlreadyConsidered) and up[up[t1]]==None and (bestBottomLength>=effectivelyNon0BLen or bestTopLength <= effectivelyNon0BLen):
						rootAlreadyConsidered=True
						listofLKcosts.append(optimizedScore)
						if networkOutput:
							listOfProbableNodes.append(t1)
					elif differentNode: #add placement to the list of legit ones
						listofLKcosts.append(optimizedScore)
						if networkOutput :
							listOfProbableNodes.append(t1)

	#print("Best placement (node "+str(bestNode)+") scores found: genetic "+str(bestScoreGenetic)+" , HnZ "+str(bestScore0)+" , Time "+str(bestScoreTime)+" , sum "+str(bestScoreGenetic+bestScore0+bestScoreTime)+" , total "+str(bestScore))
	#if aBayesPlusOn and sequentialSearch:
	if aBayesPlusOn:
		# calculate support(s) and possibly add nodes to the list of alternative placements
		finalListOfNodes=[]
		if doTimeTree:
			totSupport=0.0
			for i in topNodes:
				topNodes[i]=exp(topNodes[i])
				totSupport+=topNodes[i]
			support=topNodes[originalNode]/totSupport
			if networkOutput:
				for i in topNodes:
					topNodes[i]=topNodes[i]/totSupport
					if (i!=originalNode) and (topNodes[i]>=minBranchSupport):
						finalListOfNodes.append((i,topNodes[i]))
		else:
			support=exp(originalLK)
			totSupport=support
			for i in range(len(listofLKcosts)):
				listofLKcosts[i]=exp(listofLKcosts[i])
				totSupport+=listofLKcosts[i]
			if not totSupport:
				support=1.0
				print("Found no support?")
				print(probVect[children[node][child]])
				print(mutations[children[node][child]])
				print(dist[children[node][child]])
				if child==0:
					print(probVectUpRight[node])
				else:
					print(probVectUpLeft[node])
			else:
				support=support/totSupport
				if networkOutput:
					for i in range(len(listofLKcosts)):
						listofLKcosts[i]=listofLKcosts[i]/totSupport
						if listofLKcosts[i]>=minBranchSupport:
							finalListOfNodes.append((listOfProbableNodes[i],listofLKcosts[i]))
		totalTimeFindingParent[0]+=(time()-timeStartParentTopology)
		return bestNode, bestScore, bestBranchLengths, finalListOfNodes, support, bestRemovedPartials

	else:
		totalTimeFindingParent[0]+=(time()-timeStartParentTopology)
		return bestNode, bestScore, bestBranchLengths, [], None, bestRemovedPartials


if __name__ == "__main__":
	#TODO TODO TODO here and in the rest of the code, if doTimeTree try placing also inside multifurcations - they might not be time-multifurcations.
	#TODO function traversing the tree to find the best new root.
	def findBestRoot(tree,root,strictTopologyStopRules=strictTopologyStopRules,allowedFailsTopology=allowedFailsTopology,thresholdLogLKtopology=thresholdLogLKtopology,aBayesPlusOn=aBayesPlusOn):
		up=tree.up
		children=tree.children
		mutations=tree.mutations
		minorSequences=tree.minorSequences
		dist=tree.dist
		probVect=tree.probVect
		if doTimeTree:
			probVectTime=tree.probVectTime
		bestNode=root
		nodesToVisit=[]
		lastLK, distance, failedPasses = 0.0, 0.0, 0
		#best relative lk cost with respect to original root
		bestLKdiff=0.0
		bestNodes={}
		bestNodes[root]=0.0

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
			isTip2=(len(children[child2])==0 and len(minorSequences[child2])==0)
			isTip1=(len(children[child1])==0 and len(minorSequences[child1])==0)
			newLower, Lkcontribution=mergeVectors(vectUp1,dist[child2],isTip2,vectUp2,dist[child1],isTip1,returnLK=True,numMinor1=len(minorSequences[child2]),numMinor2=len(minorSequences[child1]))
			originalLKcost+=Lkcontribution
			vectUp1Time=None
			vectUp2Time=None
			if doTimeTree:
				vectUp1Time=probVectTime[child2]
				vectUp2Time=probVectTime[child1]
				originalLKcost+=findProbRootTime(probVectTime[root])
				originalLKcost+=mergeVectorsTime(vectUp1Time,dist[child2],vectUp2Time,dist[child1],mutRate,returnLK=True)[1]
			if mutations[child1]:
				vectUp1=passGenomeListThroughBranch(vectUp1,mutations[child1],dirIsUp=False)
			# the first LK is the likelihood cost of merges associated with the original root, these will need to be subtracted
			# the second likelihood (initialized here at 0) is the relative LK benefit of the new root - this is passed on so that one can assess if there has been an inprovement with the next branch over the previous one.
			if children[child1]:
				nodesToVisit.append((child1,vectUp1,vectUp1Time,dist[child1]+dist[child2],isTip2,len(minorSequences[child2]),originalLKcost,bestLKdiff,0))
			if mutations[child2]:
				vectUp2=passGenomeListThroughBranch(vectUp2,mutations[child2],dirIsUp=False)
			if children[child2]:
				nodesToVisit.append((child2,vectUp2,vectUp2Time,dist[child2]+dist[child1],isTip1,len(minorSequences[child1]),originalLKcost,bestLKdiff,0))

		#now visit nodes, and for each node visited, attempt to re-root at the two branches below it.
		nodesVisitedRoot=1
		while nodesToVisit:
			nodesVisitedRoot+=1
			t1,passedPartials,passedPartialsTime,distance,isTip,numMinor,LKtoRemove,lastLK,failedPasses=nodesToVisit.pop()
			childs=[children[t1][0],children[t1][1]]
			probVects=[probVect[childs[0]],probVect[childs[1]]]
			if doTimeTree:
				probVectsTime=[probVectTime[childs[0]],probVectTime[childs[1]]]
			dists=[dist[childs[0]],dist[childs[1]]]
			isTips=[]
			numMinors=[len(minorSequences[childs[0]]),len(minorSequences[childs[1]])]
			for i in range(2):
				if mutations[childs[i]]:
					probVects[i]=passGenomeListThroughBranch(probVects[i],mutations[childs[i]],dirIsUp=True)
				isTips.append((len(children[childs[i]])==0 and len(minorSequences[childs[i]])==0))
			newLKtoRemove=LKtoRemove
			newLower, Lkcontribution=mergeVectors(probVects[0],dists[0],isTips[0],probVects[1],dists[1],isTips[1],returnLK=True,numMinor1=numMinors[0],numMinor2=numMinors[1])
			newLKtoRemove+=Lkcontribution
			if doTimeTree:
				newLKtoRemove+=mergeVectorsTime(probVectsTime[0],dists[0],probVectsTime[1],dists[1],mutRate,returnLK=True)[1]
			#for each child node, assess how good the rooting is
			for i in range(2):
				try:
					upVect, lk=mergeVectors(probVects[1-i],dists[1-i],isTips[1-i],passedPartials,distance,isTip,returnLK=True,numMinor1=numMinors[1-i],numMinor2=numMinor)
					newLKtoRemoveToPass=newLKtoRemove-lk
					newRootProbVect, lkRoot=mergeVectors(upVect,dists[i]/2,False,probVects[i],dists[i]/2,isTips[i],returnLK=True,numMinor1=0,numMinor2=numMinors[i])
					rootProbLK=findProbRoot(newRootProbVect,node=t1,mutations=mutations,up=up)
					score=rootProbLK+lkRoot+lk-newLKtoRemove
					if doTimeTree:
						upVectTime, lk=mergeVectorsTime(probVectsTime[1-i],dists[1-i],passedPartialsTime,distance,mutRate,returnLK=True)
						newLKtoRemoveToPass-=lk
						newRootProbVectTime, lkRoot=mergeVectorsTime(upVectTime,dists[i]/2,probVectsTime[i],dists[i]/2,mutRate,returnLK=True)
						rootProbLK=findProbRootTime(newRootProbVectTime)
						score+=rootProbLK+lkRoot+lk
					else:
						upVectTime=None
					failedPassesNew=failedPasses
					if score>bestLKdiff:
						shorten(upVect)
						bestLKdiff=score
						bestNode=childs[i]
						failedPassesNew=0
					elif score<(lastLK-thresholdLogLKconsecutivePlacement):
						failedPassesNew+=1
					if score>=(bestLKdiff-thresholdLogLKoptimizationTopology):
						bestNodes[childs[i]]=score
						#print(str(score)+" "+str(rootProbLK)+" "+str(lkRoot)+" "+str(lk)+" "+str(newLKtoRemove)+" ")

					#assess if to pass the search onto the children nodes
					traverseChildren=False
					if children[childs[i]]:
						if strictTopologyStopRules:
							if failedPassesNew<=allowedFailsTopology and score>(bestLKdiff-thresholdLogLKtopology):
								traverseChildren=True
						else:
							if failedPassesNew<=allowedFailsTopology or score>(bestLKdiff-thresholdLogLKtopology):
								traverseChildren=True
				except:
					print("Stopping root search at node "+str(t1)+" due to error")
					traverseChildren=False

				if traverseChildren:
					if mutations[childs[i]]:
						vectToPass=passGenomeListThroughBranch(upVect,mutations[childs[i]],dirIsUp=False)
						shorten(vectToPass)
					else:
						vectToPass=upVect
					nodesToVisit.append((childs[i],vectToPass,upVectTime,dists[i],False,0,newLKtoRemoveToPass,score,failedPassesNew))

		#re-root the tree.
		if bestNode!=root:
			#reassign the lkcost of 0 for the old root in the list bestNodes to one of its children (the one that is below it in the newly rooted tree).
			rootChild=bestNode
			nodesToInvert=[]
			while up[rootChild]!=root:
				rootChild=up[rootChild]
				if up[rootChild]!=root:
					nodesToInvert.append(rootChild)
			if rootChild==children[root][0]:
				sibling=children[root][1]
			else:
				sibling=children[root][0]
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
			bestNodes[newRoot]=bestNodes.pop(bestNode)
			
			#recalculate the genome lists following the new rooting
			reCalculateAllGenomeLists(tree,newRoot)
			if doTimeTree:
				reCalculateAllGenomeListsTime(tree,newRoot,mutRate)

		else:
			newRoot=root
		
		#calculate support scores for alternative roots
		if aBayesPlusOn:
			totSupport=0.0
			tree.rootSupport=[]
			for i in range(len(up)):
				tree.rootSupport.append(None)
			normalizationV=bestNodes[newRoot]
			for currentNode in bestNodes:
				#print(bestNodes[currentNode])
				bestNodes[currentNode]=exp(bestNodes[currentNode]-normalizationV)
				totSupport+=bestNodes[currentNode]
			for currentNode in bestNodes:
				bestNodes[currentNode]=bestNodes[currentNode]/totSupport
				if bestNodes[currentNode]>=minBranchSupport:
					tree.rootSupport[currentNode]=bestNodes[currentNode]

		return newRoot


numMinorsFound=[0]
topologyUpdates=[0]
bLenUpdates=[0]
if __name__ == "__main__":
	#TODO function to find the best node in the tree where to append the new sample; traverses the tree and tries to append the sample at each node and mid-branch nodes, 
	# but stops traversing when certain criteria are met
	#TODO TODO TODO test changes
	def findBestParentForNewSample(tree,root,diffs,sample,computePlacementSupportOnly, diffsTime=None):
		up=tree.up
		children=tree.children
		probVectUpRight=tree.probVectUpRight
		probVectUpLeft=tree.probVectUpLeft
		mutations=tree.mutations
		minorSequences=tree.minorSequences
		dist=tree.dist
		probVect=tree.probVect
		probVectTotUp=tree.probVectTotUp
		nDesc0=tree.nDesc0
		if doTimeTree:
			probVectTime=tree.probVectTime
			probVectTotUpTime=tree.probVectTotUpTime
			probVectUpRightTime=tree.probVectUpRightTime
			probVectUpLeftTime=tree.probVectUpLeftTime
		bestNodes=[]
		bestNode=root
		bestBranchLengths=(False,False,oneMutBLen)
		if mutations[root]:
			diffs=passGenomeListThroughBranch(diffs,mutations[root])
		bestDiffs=diffs
		if not children[root]: #check if the new leaf is strictly less informative than already placed leaf
			# do not collapse less informative sequences also in case of HnZ
			comparison2=0
			if errorRateSiteSpecificFile or errorRateFixed or estimateErrorRate or estimateSiteSpecificErrorRate or supportFor0Branches or HnZ:
				comparison=isMinorSequence(probVect[root],diffs,onlyFindIdentical=True)
				if doTimeTree:
					comparison2=isMinorDate(probVectTime[root],diffsTime,onlyFindIdentical=True)
			else:
				comparison=isMinorSequence(probVect[root],diffs)
				if doTimeTree:
					comparison2=isMinorDate(probVectTime[root],diffsTime)
			# if we need to compute the placement support only, we have to process further and can't add the new sample as a minor sequence
			if comparison==1 and ((not doTimeTree) or (comparison2==1)) and (not computePlacementSupportOnly):
				minorSequences[root].append(sample)
				if HnZ:
					nDesc0[root]+=1
				numMinorsFound[0]+=1
				if (not onlyNambiguities) and usingErrorRate:
					updateProbVectTerminalNode(probVect[root],minorSequences[root])
				if doTimeTree:
					updateProbVectTerminalNodeTime(tree,root,diffsTime,len(minorSequences[root]),mutRate,onlyAddOne=True)
				return root, 1.0, None, diffs
			elif comparison==2 and ((not doTimeTree) or (comparison2==2)):
				totalMissedMinors[0]+=1
		rootVect=rootVector(probVect[root],False,False,tree,root)
		bestLKdiff=appendProbNode(rootVect,diffs,True,oneMutBLen)
		if doTimeTree:
			rootVectTime,rootVectTimeCost=mergeVectorsTime(probVectTime[root],0.0,diffsTime,oneMutBLen,mutRate,returnLK=True)
			bestLKdiff+=rootVectTimeCost
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
				# do not collapse less informative sequences also in case of HnZ
				comparison2=0
				if errorRateSiteSpecificFile or errorRateFixed or estimateErrorRate or estimateSiteSpecificErrorRate or supportFor0Branches or HnZ:
					comparison=isMinorSequence(probVect[t1],diffs,onlyFindIdentical=True)
					if doTimeTree:
						comparison2=isMinorDate(probVectTime[t1],diffsTime,onlyFindIdentical=True)
				else:
					comparison=isMinorSequence(probVect[t1],diffs)
					if doTimeTree:
						comparison2=isMinorDate(probVectTime[t1],diffsTime)
				# if we need to compute the placement support only, we have to process further and can't add the new sample as a minor sequence
				if comparison==1 and ((not doTimeTree) or (comparison2==1)) and (not computePlacementSupportOnly):
					minorSequences[t1].append(sample)
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
					if doTimeTree:
						updateProbVectTerminalNodeTime(tree,t1,diffsTime,len(minorSequences[t1]),mutRate,onlyAddOne=True)
					if usingErrorRate or doTimeTree:
						nodeList=[(t1,2,True,False)]
						if up[t1]!=None:
							if t1==children[up[t1]][0]:
								nodeList.append((up[t1],0,True,False))
							else:
								nodeList.append((up[t1],1,True,False))
						updatePartials(tree,nodeList)
					return t1, 1.0, None, diffs
				elif comparison==2 and ((not doTimeTree) or (comparison2==2)):
					totalMissedMinors[0]+=1

			if(dist[t1]>effectivelyNon0BLen or doTimeTree) and up[t1]!=None: # try first placing as a descendant of the mid-branch point of the branch above the current node.
				LKdiff=appendProbNode(probVectTotUp[t1],diffs,True,oneMutBLen)
				if doTimeTree:
					LKdiff+=probVectTotUpTime[t1][1]
					LKdiff+=appendProbNodeTime(probVectTotUpTime[t1][0],diffsTime,mutRate,oneMutBLen)
				if HnZ:
					if dist[t1]<=effectivelyNon0BLen:
						pNodeHnZ=up[t1]
						while dist[pNodeHnZ]<=effectivelyNon0BLen and up[pNodeHnZ]!=None:
							pNodeHnZ=up[pNodeHnZ]
						LKdiff+=getHnZ(nDesc0[pNodeHnZ]+1) - getHnZ(nDesc0[pNodeHnZ])
					else:
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
		if computePlacementSupportOnly:
			listOfProbableNodes = []
			listofLKcosts = []
			rootAlreadyConsidered = False
			listOfOptBlengths = []
			listOfPlacementTotalLhs = []
			placementAtRoot = None
			bestMidVector = []
		for nodePair in bestNodes:
			score=nodePair[1]
			if (score>=bestLKdiff-thresholdLogLKoptimization) or (computePlacementSupportOnly and score>=bestLKdiff-thresholdLogLKoptimizationTopology):
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
				bestBottomLength=estimateBranchLengthWithDerivative(midTopVector,probVect[node],fromTipC=isTip)
				newMidVector=mergeVectors(upVect,bestTopLength,False,probVect[node],bestBottomLength,isTip,isUpDown=True)
				appendingCost=appendProbNode(newMidVector,diffs,True,bestAppendingLength)
				if doTimeTree:
					if node==children[up[node]][0]:
						upVectTime=probVectUpRightTime[up[node]]
					else:
						upVectTime=probVectUpLeftTime[up[node]]
					newMidVectorTime,newMidVectorTimeCost=mergeVectorsTime(upVectTime,bestTopLength,probVectTime[node],bestBottomLength,mutRate,isUpDown=True,returnLK=True)
					if isinstance(newMidVectorTime, int):
						appendingCost+=float("-inf")
					appendingCost+=appendProbNodeTime(newMidVectorTime,diffsTime,mutRate,bestAppendingLength)
					appendingCost+=newMidVectorTimeCost
					appendingCost-=appendProbNodeTime(upVectTime,probVectTime[node],mutRate,dist[node])
				if compensanteForBranchLengthChange: #if wanted, do a more thorough examination of the appending cost, taking into account a possible change in branch length of the branch on which to be appended.
					initialCost=appendProbNode(upVect,probVect[node],isTip,dist[node])
					newPartialCost=appendProbNode(upVect,probVect[node],isTip,bestBottomLength+bestTopLength)
					optimizedScore=appendingCost+newPartialCost-initialCost
				else:
					optimizedScore=appendingCost

				if HnZ:
					if bestTopLength>effectivelyNon0BLen and bestBottomLength>effectivelyNon0BLen:
						optimizedScore+=getHnZ(2) - getHnZ(1)
					elif bestTopLength>effectivelyNon0BLen:
						optimizedScore+=getHnZ(nDesc0[node]+1) - getHnZ(nDesc0[node])
					else:
						parentNode0=up[node]
						while (dist[parentNode0]<=effectivelyNon0BLen) and up[parentNode0]!=None:
							parentNode0=up[parentNode0]
						optimizedScore+=getHnZ(nDesc0[parentNode0]+1) - getHnZ(nDesc0[parentNode0])

					#try also placing with a 0 bottom length in case it increase HnZ score.
					if bestBottomLength>effectivelyNon0BLen and dist[node]>effectivelyNon0BLen:
						altNewMidVector=mergeVectors(upVect,bestTopLength+bestBottomLength,False,probVect[node],0.0,isTip,isUpDown=True)
						altAppendingCost=appendProbNode(altNewMidVector,diffs,True,bestAppendingLength)
						if doTimeTree:
							altNewMidVectorTime,altNewMidVectorTimeCost=mergeVectorsTime(upVectTime,bestTopLength+bestBottomLength,probVectTime[node],0.0,mutRate,isUpDown=True,returnLK=True)
							if isinstance(altNewMidVectorTime, int):
								altAppendingCost+=float("-inf")
							altAppendingCost+=appendProbNodeTime(altNewMidVectorTime,diffsTime,mutRate,bestAppendingLength)
							altAppendingCost+=altNewMidVectorTimeCost
							altAppendingCost-=appendProbNodeTime(upVectTime,probVectTime[node],mutRate,dist[node])
						if compensanteForBranchLengthChange: #if wanted, do a more thorough examination of the appending cost, taking into account a possible change in branch length of the branch on which to be appended.
							initialCost=appendProbNode(upVect,probVect[node],isTip,dist[node])
							newPartialCost=appendProbNode(upVect,probVect[node],isTip,bestBottomLength+bestTopLength)
							altOptimizedScore=altAppendingCost+newPartialCost-initialCost
						altOptimizedScore+=getHnZ(nDesc0[node]+1) - getHnZ(nDesc0[node])
						if altOptimizedScore>optimizedScore:
							optimizedScore=altOptimizedScore
							bestTopLength=bestTopLength+bestBottomLength
							bestBottomLength=0.0

				if optimizedScore>=bestScore:
					bestNode=node
					bestScore=optimizedScore
					bestBranchLengths=(bestTopLength,bestBottomLength,bestAppendingLength)
					bestDiffs=diffs
					bestMidVector = newMidVector
				if computePlacementSupportOnly:
					t1 = node
					# check that the placement location is effectively different from the original node
					differentNode = True
					# topNode=node
					# if t1 == topNode:
					#	differentNode = False
					# if (not bestBottomLength):
					#	while (dist[topNode] <= effectivelyNon0BLen) and (up[topNode] != None):
					#		topNode = up[topNode]
					#	if t1 == topNode:
					#		differentNode = False
					# if t1 == children[node][1 - child]:
					#	differentNode = False
					# check that placement is not redundant
					# a placement at node X with bestTopLength = 0 could be presented by another placement:
					# at the parent of X with bestBottomLength = 0
					# or at the sibling of X with bestTopLength = 0
					if bestTopLength <= effectivelyNon0BLen:
						differentNode = False
					if dist[t1] <= effectivelyNon0BLen and up[up[t1]] != None:
						differentNode = False
					# check if this is a root placement
					if (not rootAlreadyConsidered) and (bestTopLength <= effectivelyNon0BLen):
						topNode = up[t1]
						while (dist[topNode] <= effectivelyNon0BLen) and (up[topNode] != None):
							topNode = up[topNode]
						if up[topNode] == None:
							rootAlreadyConsidered = True
							# listofLKcosts.append(optimizedScore)
							# listOfProbableNodes.append(topNode)
							# listOfOptBlengths.append((bestTopLength,bestBottomLength,bestAppendingLength))
							# record the placement at root
							placementAtRoot = (
							topNode, optimizedScore, (bestTopLength, bestBottomLength, bestAppendingLength), newMidVector)
					elif differentNode:  # add placement to the list of legit ones
						listofLKcosts.append(optimizedScore)
						listOfProbableNodes.append(t1)
						listOfOptBlengths.append((bestTopLength, bestBottomLength, bestAppendingLength))
						listOfPlacementTotalLhs.append(newMidVector)

		if computePlacementSupportOnly:
			# add the placement at root if no placements at any of its children has been recored
			if placementAtRoot:
				addPlacementAtRoot = True
				if children[root]:
					rootChild1 = children[root][0]
					rootChild2 = children[root][1]
					for placement in listOfProbableNodes:
						if placement == rootChild1 or placement == rootChild2:
							addPlacementAtRoot = False
							break
				# if no placements at any of its children has been recored
				# add the placement at root
				if addPlacementAtRoot:
					t1, optimizedScore, bestBlengths, placementTotalLh = placementAtRoot
					listofLKcosts.append(optimizedScore)
					listOfProbableNodes.append(t1)
					listOfOptBlengths.append(bestBlengths)
					listOfPlacementTotalLhs.append(placementTotalLh)

			# make sure at least one placement was found
			# because there are cases where all placements were considered as redundant
			if len(listOfProbableNodes) == 0:
				listofLKcosts.append(bestScore)
				listOfProbableNodes.append(bestNode)
				listOfOptBlengths.append(bestBranchLengths)
				listOfPlacementTotalLhs.append(bestMidVector)

			# Loop over all placements, if topBlength <= effectivelyNon0BLen, record the parent node instead of the original one
			for i in range(len(listOfOptBlengths)):
				topBlength, bottomBlength, appendingBlength = listOfOptBlengths[i]
				if topBlength <= effectivelyNon0BLen:
					topNode = listOfProbableNodes[i]
					# go to the top of the polytomy
					while (dist[topNode] <= effectivelyNon0BLen) and (up[topNode] != None):
						topNode = up[topNode]
					# go to the top of the polytomy of the parent node
					if up[topNode] != None:
						topNode = up[topNode]
						while (dist[topNode] <= effectivelyNon0BLen) and (up[topNode] != None):
							topNode = up[topNode]
						# update the placement
						listOfProbableNodes[i] = topNode
						listOfOptBlengths[i] = dist[topNode], topBlength, appendingBlength

			# calculate support(s) of possible placements
			totSupport = 0
			for i in range(len(listofLKcosts)):
				listofLKcosts[i] = exp(listofLKcosts[i])
				totSupport += listofLKcosts[i]

			possiblePlacements = []
			for i in range(len(listofLKcosts)):
				if totSupport:
					listofLKcosts[i] = listofLKcosts[i] / totSupport
				else:
					listofLKcosts[i] = 0.0

			# extract bestPlacementTotalLh
			bestPlacementTotalLh = []
			highest_support = 0
			for i in range(len(listofLKcosts)):
				if listofLKcosts[i] >= minBranchSupport:
					possiblePlacements.append((listOfProbableNodes[i], listofLKcosts[i], listOfOptBlengths[i]))

				# extract bestPlacementTotalLh
				if listofLKcosts[i] > highest_support:
					highest_support = listofLKcosts[i]
					bestPlacementTotalLh = listOfPlacementTotalLhs[i]

			# return possiblePlacements
			return possiblePlacements, bestPlacementTotalLh
		else:
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


	#TODO we know that sample "sample", with partials "newPartials", is best placed near a node resulting in logLK contribution of newChildLK
	# now explore exactly which position near node is best for placement (direct child, sibling or sibling of existing child), and which are the best branch lengths,
	#then add the sample at that position of the tree, and update all the internal probability vectors.
	#UNNECESSARY? Could probably just be replaced by the more general placeSubtreeOnTree().
	#TODO TODO TODO test
	def placeSampleOnTree(tree,node,newPartials,sample,newChildLK, bestUpLength, bestDownLength, bestAppendingLength,pseudoMutCounts,newPartialsTime=None):
		up=tree.up
		children=tree.children
		probVectUpRight=tree.probVectUpRight
		probVectUpLeft=tree.probVectUpLeft
		mutations=tree.mutations
		dist=tree.dist
		probVect=tree.probVect
		probVectTotUp=tree.probVectTotUp
		if doTimeTree:
			probVectUpRightTime=tree.probVectUpRightTime
			probVectUpLeftTime=tree.probVectUpLeftTime
			probVectTime=tree.probVectTime
			probVectTotUpTime=tree.probVectTotUpTime
			dateData=tree.dateData
		nDesc=tree.nDesc
		minorSequences=tree.minorSequences
		name=tree.name
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
			if doTimeTree:
				totRootTime,totRootTimeCost=mergeVectorsTime(probVectTime[node],0.0,newPartialsTime,bestAppendingLength,mutRate,returnLK=True)
				newChildLK+=totRootTimeCost
		else:
			if children[up[node]][0]==node:
				child=0
				vectUp=probVectUpRight[up[node]]
				if doTimeTree:
					vectUpTime=probVectUpRightTime[up[node]]
			else:
				child=1
				vectUp=probVectUpLeft[up[node]]
				if doTimeTree:
					vectUpTime=probVectUpLeftTime[up[node]]
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
						#updating nDesc0 in case of branch length change
						if HnZ:
							updateNDesc0whenChangingDist(tree,node,bestDownLength)
						dist[node]=bestDownLength
						nodeList=[(node,2,True,doTimeTree),(up[node],child,True,doTimeTree)]
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
			if doTimeTree:
				probVectRootTime=probVectTime[node]
				probOldRoot += findProbRootTime(probVectRootTime)
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
			if doTimeTree:
				probVectRootTime,probRootTime=mergeVectorsTime(probVectTime[node],bestLeftLength,newPartialsTime,bestRightLength,mutRate,returnLK=True)
				probRoot+=probRootTime
				probRoot+= findProbRootTime(probVectRootTime)

				rootUpRightTime=rootVectorTime(newPartialsTime,bestRightLength,mutRate)
			if HnZ:
				probRoot+=getHnZ(2)-getHnZ(1)
			parentLKdiff=probRoot-probOldRoot
			if parentLKdiff<=newChildLK: #best is just placing as descendant of the root
				bestRightLength=bestAppendingLength
				bestLeftLength=False
				probVectRoot=mergeVectors(probVect[node],bestLeftLength,isTip,rootNewPartials,bestRightLength,True)
				rootUpRight=rootVector(rootNewPartials,bestRightLength,True,tree,node)
				if doTimeTree:
					probVectRootTime=mergeVectorsTime(probVectTime[node],bestLeftLength,newPartialsTime,bestRightLength,mutRate)
					rootUpRightTime=rootVectorTime(newPartialsTime,bestRightLength,mutRate)
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
			if doTimeTree:
				probVectTime[newRoot]=probVectRootTime
				probVectUpRightTime[newRoot]=rootUpRightTime
				probVectUpLeftTime[newRoot]=rootVectorTime(probVectTime[node],bestLeftLength,mutRate)
			mutations[newRoot]=mutations[node]
			mutations[node]=[]
			up[node]=newRoot
			dist[node]=bestLeftLength
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
			if bestRightLength or doTimeTree:
				probVectTotUp[newNode]=mergeVectors(probVectUpLeft[newRoot],bestRightLength/2,False,rootNewPartials,bestRightLength/2,True,isUpDown=True)
				shorten(probVectTotUp[newNode])
			if doTimeTree:
				probVectTime[newNode]=newPartialsTime
				dateData[newNode]=newPartialsTime
				newTot,newTotProb=mergeVectorsTime(probVectUpLeftTime[newRoot],bestRightLength/2,newPartialsTime,bestRightLength/2,mutRate,isUpDown=True,returnLK=True)
				if isinstance(newTot, int):
					resolveTimeInconsistency(tree,newNode,newTot,mutRate)
					newTot,newTotProb=mergeVectorsTime(probVectUpLeftTime[newRoot],bestRightLength/2,newPartialsTime,bestRightLength/2,mutRate,isUpDown=True,returnLK=True)
				newTotProb-=appendProbNodeTime(probVectUpLeftTime[newRoot],newPartialsTime,mutRate,bestRightLength)
				probVectTotUpTime[newNode]=(newTot,newTotProb)
			nodeList=[(node,2,True,doTimeTree)]
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
		oldLen=dist[node]
		dist[node]=bestDownLength
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
		if HnZ and bestUpLength<=effectivelyNon0BLen:
			parent0=newInternalNode
			addendum=1
			if bestDownLength<=effectivelyNon0BLen and oldLen>effectivelyNon0BLen:
				addendum=nDesc0[node]
			while up[parent0]!=None and (dist[parent0]<=effectivelyNon0BLen):
				parent0=up[parent0]
				nDesc0[parent0]+=addendum

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
		if doTimeTree:
			probVectTime[newNode]=newPartialsTime
			dateData[newNode]=newPartialsTime
			probVectTime[newInternalNode]=mergeVectorsTime(probVectTime[node],bestDownLength,newPartialsTime,bestAppendingLength,mutRate,isUpDown=False)
			probVectUpRightTime[newInternalNode]=mergeVectorsTime(vectUpTime,bestUpLength,newPartialsTime,bestAppendingLength,mutRate,isUpDown=True)
			if isinstance(probVectUpRightTime[newInternalNode], int):
				resolveTimeInconsistency(tree,newInternalNode,probVectUpRightTime[newInternalNode],mutRate)
				if children[up[node]][0]==node:
					vectUpTime=probVectUpRightTime[up[node]]
				else:
					vectUpTime=probVectUpLeftTime[up[node]]
				probVectUpRightTime[newInternalNode]=mergeVectorsTime(vectUpTime,bestUpLength,newPartialsTime,bestAppendingLength,mutRate,isUpDown=True)
			probVectUpLeftTime[newInternalNode]=mergeVectorsTime(vectUpTime,bestUpLength,probVectTime[node],bestDownLength,mutRate,isUpDown=True)
			if isinstance(probVectUpLeftTime[newInternalNode], int):
				resolveTimeInconsistency(tree,newInternalNode,probVectUpLeftTime[newInternalNode],mutRate)
				if children[up[node]][0]==node:
					vectUpTime=probVectUpRightTime[up[node]]
				else:
					vectUpTime=probVectUpLeftTime[up[node]]
				probVectUpLeftTime[newInternalNode]=mergeVectorsTime(vectUpTime,bestUpLength,probVectTime[node],bestDownLength,mutRate,isUpDown=True)
			newTot,newTotProb=mergeVectorsTime(vectUpTime,bestUpLength/2,probVectTime[newInternalNode],bestUpLength/2,mutRate,isUpDown=True,returnLK=True)
			if isinstance(newTot, int):
				resolveTimeInconsistency(tree,newInternalNode,newTot,mutRate)
				if children[up[node]][0]==node:
					vectUpTime=probVectUpRightTime[up[node]]
				else:
					vectUpTime=probVectUpLeftTime[up[node]]
				newTot,newTotProb=mergeVectorsTime(vectUpTime,bestUpLength/2,probVectTime[newInternalNode],bestUpLength/2,mutRate,isUpDown=True,returnLK=True)
			newTotProb-=appendProbNodeTime(vectUpTime,probVectTime[newInternalNode],mutRate,bestUpLength)
			probVectTotUpTime[newInternalNode]=(newTot,newTotProb)
			newTot,newTotProb=mergeVectorsTime(probVectUpLeftTime[newInternalNode],bestAppendingLength/2,newPartialsTime,bestAppendingLength/2,mutRate,isUpDown=True,returnLK=True)
			if isinstance(newTot, int):
				resolveTimeInconsistency(tree,newNode,newTot,mutRate)
				newTot,newTotProb=mergeVectorsTime(probVectUpLeftTime[newInternalNode],bestAppendingLength/2,newPartialsTime,bestAppendingLength/2,mutRate,isUpDown=True,returnLK=True)
			newTotProb-=appendProbNodeTime(probVectUpLeftTime[newInternalNode],newPartialsTime,mutRate,bestAppendingLength)
			probVectTotUpTime[newNode]=(newTot,newTotProb)
		if probVect[newInternalNode]==None:
			print("Problem in placeSampleOnTree(), probVect is None")
			raise Exception("exit")
		if probVectUpRight[newInternalNode]==None:
			print("Problem in placeSampleOnTree(), probVectUpRight is None")
			raise Exception("exit")
		if probVectUpLeft[newInternalNode]==None:
			print("Problem in placeSampleOnTree(), probVectUpLeft is None")
			raise Exception("exit")
		if bestUpLength or doTimeTree:
			probVectTotUp[newInternalNode]=mergeVectors(vectUp,bestUpLength/2,False,probVect[newInternalNode],bestUpLength/2,False,isUpDown=True)
			if passUpMutations:
				probVectTotUp[newInternalNode]=passGenomeListThroughBranch(probVectTotUp[newInternalNode],mutations[node],dirIsUp=True)
			shorten(probVectTotUp[newInternalNode])
		else:
			probVectTotUp[newInternalNode]=None
		if bestAppendingLength or doTimeTree:
			probVectTotUp[newNode]=mergeVectors(probVectUpLeft[newInternalNode],bestAppendingLength/2,False,newPartials,bestAppendingLength/2,True,isUpDown=True)
			if passUpMutations:
				probVectTotUp[newNode]=passGenomeListThroughBranch(probVectTotUp[newNode],mutations[node],dirIsUp=True)
			shorten(probVectTotUp[newNode])
			if bestAppendingLength:
				updatePesudoCounts(probVectUpLeft[newInternalNode],newPartials,pseudoMutCounts)
		else:
			probVectTotUp[newNode]=None
		if (not bestDownLength) and (not doTimeTree):
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
		nodeList=[(node,2,True,doTimeTree),(up[newInternalNode],child,True,doTimeTree)]
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


	# TODO traverse the tree to optimize (only) the length of the branches using the derivative approach.
	#TODO TODO TODO test
	def traverseTreeToOptimizeBranchLengths(tree,root,testing=False,fastPass=False):
		up=tree.up
		children=tree.children
		probVectUpRight=tree.probVectUpRight
		probVectUpLeft=tree.probVectUpLeft
		mutations=tree.mutations
		dist=tree.dist
		probVect=tree.probVect
		if doTimeTree:
			probVectTime=tree.probVectTime
		minorSequences=tree.minorSequences
		dirty=tree.dirty
		nDesc0=tree.nDesc0
		totalLKimprovementBL=0.0
		updates=0
		nodesTraversed=0
		if children[root]:
			child1=children[root][0]
			child2=children[root][1]
			if (dist[child1]>effectivelyNon0BLen) or (dist[child2]>effectivelyNon0BLen):
				totDist=(dist[child1]+dist[child2])*lRef
				fromTip1=False
				if len(children[child1])==0 and len(minorSequences[child1])==0:
					fromTip1=True
				fromTip2=False
				if len(children[child2])==0 and len(minorSequences[child2])==0:
					fromTip2=True
				probVect1=probVect[child1]
				if mutations[child1]:
					probVect1=passGenomeListThroughBranch(probVect1,mutations[child1],dirIsUp=True)
				probVect2=probVect[child2]
				if mutations[child2]:
					probVect2=passGenomeListThroughBranch(probVect2,mutations[child2],dirIsUp=True)
				bestCost=float("-inf")
				bestBL1=None
				for i in range(max(1,round(totDist))*2+1):
					bLen1=min(totDist,float(i)/2)
					bLen2=max(totDist-bLen1,0.0)
					bLen1=bLen1/lRef
					bLen2=bLen2/lRef
					rootVector,cost=mergeVectors(probVect1,bLen1,fromTip1,probVect2,bLen2,fromTip2,returnLK=True)
					if mutations[root]:
						rootVector=passGenomeListThroughBranch(rootVector,mutations[root],dirIsUp=True)
					cost+=findProbRoot(rootVector)
					if HnZ:
						if bLen1<effectivelyNon0BLen:
							cost+=getHnZ(nDesc0[child1]+1)-getHnZ(nDesc0[child1])
						if bLen2<effectivelyNon0BLen:
							cost+=getHnZ(nDesc0[child2]+1)-getHnZ(nDesc0[child2])
					if doTimeTree:
						rootVectorTime,costTime=mergeVectorsTime(probVectTime[child1],bLen1,probVectTime[child2],bLen2,mutRate,returnLK=True)
						costTime+=findProbRootTime(rootVectorTime)
						cost+=costTime
					if cost>bestCost:
						bestCost=cost
						bestBL1=bLen1

				bestBL2=max(dist[child1]+dist[child2]-bestBL1,0.0)
				try:
					if HnZ:
						updateNDesc0whenChangingDist(tree,child1,bestBL1)
					dist[child1]=bestBL1
					if not fastPass:
						nodeList=[(child1,2,True,doTimeTree),(root,0,True,doTimeTree)]
						updatePartials(tree,nodeList)
					if HnZ:
						updateNDesc0whenChangingDist(tree,child2,bestBL2)
					dist[child2]=bestBL2
					if not fastPass:
						nodeList=[(child2,2,True,doTimeTree),(root,0,True,doTimeTree)]
						updatePartials(tree,nodeList)
				except:
					if HnZ:
						updateNDesc0whenChangingDist(tree,child2,bestBL2)
					dist[child2]=bestBL2
					if not fastPass:
						nodeList=[(child2,2,True,doTimeTree),(root,1,True,doTimeTree)]
						updatePartials(tree,nodeList)
					if HnZ:
						updateNDesc0whenChangingDist(tree,child1,bestBL1)
					dist[child1]=bestBL1
					if not fastPass:
						nodeList=[(child1,2,True,doTimeTree),(root,0,True,doTimeTree)]
						updatePartials(tree,nodeList)

			if children[children[root][0]]:
				nodesToTraverse=[children[children[root][0]][0],children[children[root][0]][1]]
			else:
				nodesToTraverse=[]
			if children[children[root][1]]:
				nodesToTraverse.append(children[children[root][1]][0])
				nodesToTraverse.append(children[children[root][1]][1])
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
					if testing or doTimeTree or HnZ:
						currentCost=appendProbNode(upVect,probVect[node],isTip,dist[node])
						newCost=appendProbNode(upVect,probVect[node],isTip,bestLength)
						if HnZ:
							parentNode0=up[node]
							while (dist[parentNode0]<=effectivelyNon0BLen) and up[parentNode0]!=None:
								parentNode0=up[parentNode0]
							if dist[node]>effectivelyNon0BLen:
								currentCost+=getHnZ(nDesc0[parentNode0])+getHnZ(nDesc0[node])
								if bestLength>effectivelyNon0BLen:
									newCost+=getHnZ(nDesc0[parentNode0])+getHnZ(nDesc0[node])
								else:
									newCost+=getHnZ(nDesc0[parentNode0]+nDesc0[node]-1)
							else:
								currentCost+=getHnZ(nDesc0[parentNode0])
								if bestLength>effectivelyNon0BLen:
									newCost+=getHnZ(nDesc0[parentNode0]+1-nDesc0[node])+getHnZ(nDesc0[node])
								else:
									newCost+=getHnZ(nDesc0[parentNode0])
						if testing:
							totalLKimprovementBL+=newCost-currentCost
					# test also length 0 
					if HnZ and dist[node]>effectivelyNon0BLen and bestLength>effectivelyNon0BLen:
						cost0=appendProbNode(upVect,probVect[node],isTip,0.0)
						if cost0>-1000000:
							cost0+=getHnZ(nDesc0[parentNode0]+nDesc0[node]-1)
							if cost0 > newCost:
								bestLength=0.0
								newCost=cost0
					if (doTimeTree or HnZ) and currentCost>newCost:
						bestLength=dist[node]
						newCost=currentCost
					
					if bestLength or dist[node]:
						if (not bestLength) or (not dist[node]) or dist[node]/bestLength>1.01 or dist[node]/bestLength<0.99:
							if HnZ: #updating nDesc0 in case branch lengths move to 0 or from 0
								updateNDesc0whenChangingDist(tree,node,bestLength)
							dist[node]=bestLength
							updates+=1
							if not fastPass:
								nodeList=[(node,2,True,doTimeTree),(up[node],child,True,doTimeTree)]
								updatePartials(tree,nodeList)
						else:
							dirty[node]=False
					else:
						dirty[node]=False
				else:
					dirty[node]=False
			for child in children[node]:
				nodesToTraverse.append(child)
		if testing:
			return totalLKimprovementBL
		else:
			return updates
			

	# TODO we know that subtree "appendedNode", with partials "newPartials", is best placed as child of "node" resulting in logLK contribution of newChildLK
	# now explore exactly which position near node is best for placement (direct child, sibling or sibling of existing child), and which are the best branch lengths,
	#then add the subtree at that position of the tree, and update all the internal probability vectors.
	#TODO TODO TODO test
	def placeSubtreeOnTree(tree,node,newPartials,appendedNode,newChildLK,bestBranchLengths,newPartialsTime=None): #,parentNode,parentNodeReplacements=1
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
		if doTimeTree:
			probVectUpRightTime=tree.probVectUpRightTime
			probVectUpLeftTime=tree.probVectUpLeftTime
			probVectTime=tree.probVectTime
			probVectTotUpTime=tree.probVectTotUpTime
		nDesc0=tree.nDesc0
		bestAppendingLength=bestBranchLengths[2]
		bestUpLength=bestBranchLengths[0]
		bestDownLength=bestBranchLengths[1]
		tryNewRoot=False
		if children[up[node]][0]==node:
			child=0
			vectUp=probVectUpRight[up[node]]
			if doTimeTree:
				vectUpTime=probVectUpRightTime[up[node]]
		else:
			child=1
			vectUp=probVectUpLeft[up[node]]
			if doTimeTree:
				vectUpTime=probVectUpLeftTime[up[node]]
		
		#test if we should also attempt placing the subtree as child of a new root
		if not bestUpLength:
			pNode=up[node]
			if not doTimeTree:
				while (not dist[pNode]) and (up[pNode]!=None):
					pNode=up[pNode]
			if up[pNode]==None:
				root=pNode
				tryNewRoot=True
				if (not bestDownLength) or (bestDownLength>1.01*dist[node]) or (bestDownLength<0.99*dist[node]):
					if HnZ:
						#updating nDesc0 in case branch lengths move to 0 or from 0
						updateNDesc0whenChangingDist(tree,node,bestDownLength)
					dist[node]=bestDownLength
					nodeList=[(node,2,True,doTimeTree),(up[node],child,True,doTimeTree)]
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
			if doTimeTree:
				probOldRoot += findProbRootTime(probVectTime[node])
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
			if doTimeTree:
				probVectRootTime,probRootTime = mergeVectorsTime(probVectTime[node],bestLeftLength,newPartialsTime,bestRightLength,mutRate,returnLK=True,isUpDown=False)
				probRoot+=probRootTime
				probRoot+= findProbRootTime(probVectRootTime)
				rootUpRightTime=rootVectorTime(newPartialsTime,bestRightLength,mutRate)
			parentLKdiff=probRoot-probOldRoot
			if parentLKdiff<=newChildLK: #best is just placing as descendant of the root
				bestRightLength=bestAppendingLength
				bestLeftLength=False
				probVectRoot=mergeVectors(probVect[node],bestLeftLength,isTip,rootNewPartials,bestRightLength,appendedIsTip)
				rootUpRight=rootVector(rootNewPartials,bestRightLength,appendedIsTip,tree,node)
				if doTimeTree:
					probVectRootTime=mergeVectorsTime(probVectTime[node],bestLeftLength,newPartialsTime,bestRightLength,mutRate)
					rootUpRightTime=rootVectorTime(newPartialsTime,bestRightLength,mutRate)
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
			if doTimeTree:
				probVectTime[newRoot]=probVectRootTime
				probVectUpRightTime[newRoot]=rootUpRightTime
				probVectUpLeftTime[newRoot]=rootVectorTime(probVectTime[node],bestLeftLength,mutRate)
			mutations[newRoot]=mutations[node]
			mutations[node]=[]
			up[node]=newRoot
			dist[node]=bestLeftLength
			children[newRoot][0]=node
			children[newRoot][1]=appendedNode
			dist[appendedNode]=bestRightLength
			replacements[appendedNode]+=1
			nodeList=[(node,2,True,doTimeTree),(appendedNode,2,True,doTimeTree)]
			if HnZ:
				if dist[node]>effectivelyNon0BLen:
					nDesc0[newRoot]=1
				else:
					nDesc0[newRoot]=nDesc0[node]
				if dist[appendedNode]>effectivelyNon0BLen:
					nDesc0[newRoot]+=1
				else:
					nDesc0[newRoot]+=nDesc0[appendedNode]
			updatePartials(tree,nodeList)
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
		oldDist=dist[node]
		dist[appendedNode]=bestAppendingLength
		dist[newInternalNode]=bestUpLength
		dist[node]=bestDownLength
		if doTimeTree:
			probVectTime[newInternalNode]=mergeVectorsTime(probVectTime[node],bestDownLength,newPartialsTime,bestAppendingLength,mutRate)
			probVectUpRightTime[newInternalNode]=mergeVectorsTime(vectUpTime,bestUpLength,newPartialsTime,bestAppendingLength,mutRate,isUpDown=True)
			if isinstance(probVectUpRightTime[newInternalNode], int):
				resolveTimeInconsistency(tree,newInternalNode,probVectUpRightTime[newInternalNode],mutRate)
				if children[up[node]][0]==node:
					vectUpTime=probVectUpRightTime[up[node]]
				else:
					vectUpTime=probVectUpLeftTime[up[node]]
				probVectUpRightTime[newInternalNode]=mergeVectorsTime(vectUpTime,bestUpLength,newPartialsTime,bestAppendingLength,mutRate,isUpDown=True)
			probVectUpLeftTime[newInternalNode]=mergeVectorsTime(vectUpTime,bestUpLength,probVectTime[node],bestDownLength,mutRate,isUpDown=True)
			if isinstance(probVectUpLeftTime[newInternalNode], int):
				resolveTimeInconsistency(tree,newInternalNode,probVectUpLeftTime[newInternalNode],mutRate)
				if children[up[node]][0]==node:
					vectUpTime=probVectUpRightTime[up[node]]
				else:
					vectUpTime=probVectUpLeftTime[up[node]]
				probVectUpLeftTime[newInternalNode]=mergeVectorsTime(vectUpTime,bestUpLength,probVectTime[node],bestDownLength,mutRate,isUpDown=True)
		if HnZ:
			if dist[node]<=effectivelyNon0BLen:
				nDesc0[newInternalNode]=nDesc0[node]
			else:
				nDesc0[newInternalNode]=1
			nDesc0ToBeAdded=0
			if dist[appendedNode]>effectivelyNon0BLen:
				nDesc0ToBeAdded=1
			else:
				nDesc0ToBeAdded=nDesc0[appendedNode]
			nDesc0[newInternalNode]+=nDesc0ToBeAdded
			#TODO TODO TODO alter toBeAddd to consider the case oldDist was 0 and now is not, and vice-versa
			nDesc0ToBeAdded=0
			if (oldDist>effectivelyNon0BLen) and (dist[newInternalNode]<=effectivelyNon0BLen):
				nDesc0ToBeAdded=nDesc0[newInternalNode]-1
			elif (oldDist<=effectivelyNon0BLen) and (dist[newInternalNode]>effectivelyNon0BLen):
				nDesc0ToBeAdded=1-nDesc0[node]
			elif (oldDist<=effectivelyNon0BLen) and (dist[newInternalNode]<=effectivelyNon0BLen):
				nDesc0ToBeAdded=nDesc0[newInternalNode]-nDesc0[node]
			if nDesc0ToBeAdded!=0:
				parent0=up[newInternalNode]
				while True:
					nDesc0[parent0]+=nDesc0ToBeAdded
					if dist[parent0]>effectivelyNon0BLen:
						break
					parent0=up[parent0]
					if parent0==None:
						break

		if (not bestAppendingLength) and (not doTimeTree):
			probVectTotUp[appendedNode]=None
		if bestUpLength or doTimeTree:
			probVectTotUp[newInternalNode]=mergeVectors(vectUp,bestUpLength/2,False,probVect[newInternalNode],bestUpLength/2,False,isUpDown=True)
			shorten(probVectTotUp[newInternalNode])
			if doTimeTree:
				newTot,newTotProb=mergeVectorsTime(vectUpTime,bestUpLength/2,probVectTime[newInternalNode],bestUpLength/2,mutRate,isUpDown=True,returnLK=True)
				if isinstance(newTot, int):
					resolveTimeInconsistency(tree,newInternalNode,newTot,mutRate)
					if children[up[node]][0]==node:
						vectUpTime=probVectUpRightTime[up[node]]
					else:
						vectUpTime=probVectUpLeftTime[up[node]]
					newTot,newTotProb=mergeVectorsTime(vectUpTime,bestUpLength/2,probVectTime[newInternalNode],bestUpLength/2,mutRate,isUpDown=True,returnLK=True)
				newTotProb-=appendProbNodeTime(vectUpTime,probVectTime[newInternalNode],mutRate,bestUpLength)
				probVectTotUpTime[newInternalNode]=(newTot,newTotProb)
				
		if (not bestDownLength) and (not doTimeTree):
			probVectTotUp[node]=None
		nodeList=[(node,2,True,doTimeTree),(up[newInternalNode],child,True,doTimeTree),(appendedNode,2,True,doTimeTree)]
		updatePartials(tree,nodeList)		
		return None


	# TODO remove node from the current position in the tree and re-attach it at a new given place new bestNode.
	# First remove node from the tree, then update the genome lists;
	# then find the exact best reattachment of node and update the genome lists again using function placeSubtreeOnTree().
	#TODO TODO TODO test
	def cutAndPasteNode(tree,node,bestNode,bestBranchLengths,bestLK,passedProbVect,passedProbVectTime=None):
		up=tree.up
		children=tree.children
		probVectUpRight=tree.probVectUpRight
		probVectUpLeft=tree.probVectUpLeft
		mutations=tree.mutations
		dist=tree.dist
		probVect=tree.probVect
		if doTimeTree:
			probVectUpRightTime=tree.probVectUpRightTime
			probVectUpLeftTime=tree.probVectUpLeftTime
			probVectTime=tree.probVectTime
		minorSequences=tree.minorSequences
		nDesc0=tree.nDesc0
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
			if HnZ:
				if dist[parentNode]<=effectivelyNon0BLen:
					if dist[node]>effectivelyNon0BLen:
						toBeRemoved=-1
					else:
						toBeRemoved=-nDesc0[node]
					if dist[sibling]<=effectivelyNon0BLen and (dist[sibling]+dist[parentNode])>effectivelyNon0BLen:
						toBeRemoved+=(1-nDesc0[sibling])
					parent0=parentNode
					while (dist[parent0]<=effectivelyNon0BLen) and up[parent0]!=None:
						parent0=up[parent0]
						nDesc0[parent0]+=toBeRemoved
						if nDesc0[parent0]<=0:
							print("problem removing subtree")
							raise Exception("exit")
							
		up[sibling]=up[parentNode]
		dist[sibling]=dist[sibling]+dist[parentNode]
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
				if doTimeTree:
					probVectUpRightTime[sibling]=rootVectorTime(probVectTime[children[sibling][1]],dist[children[sibling][1]],mutRate)
					probVectUpLeftTime[sibling]=rootVectorTime(probVectTime[children[sibling][0]],dist[children[sibling][0]],mutRate)
				nodeList=[(children[sibling][0],2,True,doTimeTree),(children[sibling][1],2,True,doTimeTree)]
				updatePartials(tree,nodeList)
		else:
			nodeList=[(sibling,2,True,doTimeTree),(up[sibling],childP,True,doTimeTree)]
			updatePartials(tree,nodeList)
		newRoot = placeSubtreeOnTree(tree,bestNode,passedProbVect,node,bestLK,bestBranchLengths,newPartialsTime=passedProbVectTime)
		
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
			if doTimeTree:
				totalLK+=calculateTreeLikelihoodTime(tree,currentRoot,mutRate)
			intermediateLKsFile.write("Topology "+str(topologyChanges[0])+", LK: "+str(totalLK)+"\n")

		#if the root of the tree has changed, return the new root
		if up[sibling]==None:
			if newRoot!=None:
				return newRoot
			return sibling
		return newRoot


	# TODO try to find a re-placement of a dirty node of the tree to improve the topology.
	# Pretend to cut out the subtree at this node (so to keep track of the effect of subtree removal on the likelihoods), and look for somewhere else in the tree where to attach it (an SPR move).
	# To find the best location of the new re-attachment, we traverse the tree starting at the current attachment node, and for each node visited we evaluate the reattachment,
	# as done by the findBestParentTopology() function.
	# To avoid traversing the whole tree for each SPR move, we use stopping conditions similar to those used in the initial sample placement process.
	# After we find the best SPR move for the given node, we execute it using the cutAndPasteNode() function.
	# TODO TODO TODO test
	def traverseTreeForTopologyUpdate(tree,node,strictTopologyStopRules=strictTopologyStopRules,allowedFailsTopology=allowedFailsTopology,thresholdLogLKtopology=thresholdLogLKtopology,thresholdTopologyPlacement=thresholdTopologyPlacement):
		up=tree.up
		children=tree.children
		probVectUpRight=tree.probVectUpRight
		probVectUpLeft=tree.probVectUpLeft
		mutations=tree.mutations
		dist=tree.dist
		probVect=tree.probVect
		if doTimeTree:
			probVectUpRightTime=tree.probVectUpRightTime
			probVectUpLeftTime=tree.probVectUpLeftTime
			probVectTime=tree.probVectTime
		minorSequences=tree.minorSequences
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
			#evaluate current placement
			parentNode=up[node]
			if children[parentNode][0]==node:
				child=0
				vectUp=probVectUpRight[parentNode]
				if doTimeTree:
					vectUpTime=probVectUpRightTime[parentNode]
			else:
				child=1
				vectUp=probVectUpLeft[parentNode]
				if doTimeTree:
					vectUpTime=probVectUpLeftTime[parentNode]
			sibling=children[parentNode][1-child]
			if mutations[node]:
				vectUp=passGenomeListThroughBranch(vectUp,mutations[node],dirIsUp=False)
			#score of current tree
			bestCurrenBLen=dist[node]
			isTip=(len(children[node])==0) and (len(minorSequences[node])==0)
			originalLK=appendProbNode(vectUp,probVect[node],isTip,bestCurrenBLen)
			geneticLK=originalLK
			originalLKTime=0
			if doTimeTree:
				if up[parentNode]==None:
					originalLKTime=mergeVectorsTime(probVectTime[node],dist[node],probVectTime[sibling],dist[sibling],mutRate,returnLK=True,isUpDown=False)[1]
					originalLKTime+=findProbRootTime(probVectTime[parentNode])
				else:
					originalLKTime=appendProbNodeTime(vectUpTime,probVectTime[node],mutRate,bestCurrenBLen)
					if children[up[parentNode]][0]==parentNode:
						vectUpUpTime=probVectUpRightTime[up[parentNode]]
					else:
						vectUpUpTime=probVectUpLeftTime[up[parentNode]]
					originalVectLKTime,originalVectLKTimeCost = mergeVectorsTime(vectUpUpTime,dist[parentNode],probVectTime[sibling],dist[sibling],mutRate,returnLK=True,isUpDown=True)
					originalLKTime+=originalVectLKTimeCost
					originalLKTime-=appendProbNodeTime(vectUpUpTime,probVectTime[sibling],mutRate,dist[sibling]+dist[parentNode])
				originalLK+=originalLKTime
			originalLK0=0
			if HnZ:
				parentNode0=up[node]
				while (dist[parentNode0]<=effectivelyNon0BLen) and up[parentNode0]!=None:
					parentNode0=up[parentNode0]
				if dist[node]>effectivelyNon0BLen:
					originalLK0=getHnZ(nDesc0[parentNode0]) - getHnZ(nDesc0[parentNode0]-1)
				else:
					originalLK0=getHnZ(nDesc0[parentNode0]) - ( getHnZ(nDesc0[parentNode0]-nDesc0[node]) + getHnZ(nDesc0[node]) )
				originalLK+= originalLK0
			bestCurrentLK=originalLK
			if ((geneticLK<thresholdTopologyPlacement) or (supportFor0Branches and aBayesPlusOn)) and up[up[node]]!=None:
				bestCurrenBLen=estimateBranchLengthWithDerivative(vectUp,probVect[node],fromTipC=isTip)
				if bestCurrenBLen or dist[node]:
					if (not bestCurrenBLen) or (not dist[node]) or dist[node]/bestCurrenBLen>1.01 or dist[node]/bestCurrenBLen<0.99:
						bLenChanged=True
					bestCurrentLK=appendProbNode(vectUp,probVect[node],isTip,bestCurrenBLen)
					bestCurrentLKGen=bestCurrentLK
					if doTimeTree:
						bestCurrentLKTime=appendProbNodeTime(vectUpTime,probVectTime[node],mutRate,bestCurrenBLen)
						bestCurrentLK+=bestCurrentLKTime
					if HnZ:
						if bestCurrenBLen>effectivelyNon0BLen:
							if dist[node]>effectivelyNon0BLen:
								bestCurrentLK0= getHnZ(nDesc0[parentNode0]) - getHnZ(nDesc0[parentNode0]-1)
							else:
								bestCurrentLK0= getHnZ(nDesc0[parentNode0]+1-nDesc0[node]) - getHnZ(nDesc0[parentNode0]-nDesc0[node])
						else:
							if dist[node]>effectivelyNon0BLen:
								bestCurrentLK0= getHnZ(nDesc0[parentNode0] + nDesc0[node] -1) - (getHnZ(nDesc0[parentNode0]) + getHnZ(nDesc0[node]))
							else:
								bestCurrentLK0= getHnZ(nDesc0[parentNode0]) - ( getHnZ(nDesc0[parentNode0]-nDesc0[node]) + getHnZ(nDesc0[node]) )
						bestCurrentLK+=bestCurrentLK0
					if bestCurrentLK<originalLK:
						bestCurrenBLen=dist[node]
						bestCurrentLK=originalLK
						bLenChanged=False
					else:
						if doTimeTree:
							originalLKTime=bestCurrentLKTime
						if HnZ:
							originalLK0=bestCurrentLK0
						geneticLK=bestCurrentLKGen
					if bestCurrentLK==float("-inf"):
						print("Found an infinite cost of bestCurrentLK "+str(bestCurrentLK)+" using appendProbNode()")
						raise Exception("exit")

			topologyUpdated=False
			# in case of HnZ, also try to replace 0-cost subtrees, as they might have better modifers somewhere else
			if ((bestCurrentLK<thresholdTopologyPlacement or dist[node] or HnZ or doTimeTree) and (not doNotImproveTopology)) or ((dist[node] or supportFor0Branches) and aBayesPlusOn) :
			#if ((bestCurrentLK<thresholdTopologyPlacement) and (not doNotImproveTopology)) or dist[node] or (supportFor0Branches and aBayesPlusOn) or HnZ or doTimeTree:
				#now find the best place on the tree where to re-attach the subtree rooted at "node"
				#but to do that we need to consider new vector probabilities after removing the node that we want to replace
				# this is done using findBestParentTopology().
				bestNodeSoFar , bestLKdiff , bestBranchLengths, listOfBestPlacements, branchSupport, passedProbVect = findBestParentTopology(tree,parentNode,child,bestCurrentLK,bestCurrenBLen,strictTopologyStopRules=strictTopologyStopRules,allowedFailsTopology=allowedFailsTopology,thresholdLogLKtopology=thresholdLogLKtopology)
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
					if bestNodeSoFar==topNode and (not bestBranchLengths[1]) and (not doTimeTree):
						topologyUpdated=False
					parentNode=up[node]
					if node==children[parentNode][0]:
						sibling=children[parentNode][1]
					else:
						sibling=children[parentNode][0]
					if bestNodeSoFar==sibling:
						topologyUpdated=False
					if up[bestNodeSoFar]==sibling and (not bestBranchLengths[0]) and (not doTimeTree):
						topologyUpdated=False

					if topologyUpdated:
						topologyUpdates[0]+=1
						totalImprovement=(bestLKdiff-originalLK)
						if originalLK==float("-inf"):
							totalImprovement=(bestLKdiff-bestCurrentLK)
						if totalImprovement==float("inf"):
							print("Found an infinite topology improvement of "+str(totalImprovement)+" from originalLK "+str(originalLK)+", bestCurrentLK "+str(bestCurrentLK)+" and bestLKdiff "+str(bestLKdiff))
							raise Exception("exit")
						if doTimeTree:
							passedProbVectTime=probVectTime[node]
						else:
							passedProbVectTime=None
						newRoot = cutAndPasteNode(tree,node,bestNodeSoFar,bestBranchLengths,bestLKdiff,passedProbVect,passedProbVectTime=passedProbVectTime)
						bLenChanged=False
				if (not topologyUpdated) and aBayesPlusOn:
					if networkOutput:
						alternativePlacements[node]=listOfBestPlacements
					support[node]=branchSupport
						
			if (not topologyUpdated) and bLenChanged:
				bLenUpdates[0]+=1
				if HnZ: # updating nDesc0
					updateNDesc0whenChangingDist(tree,node,bestCurrenBLen)
				dist[node]=bestCurrenBLen
				nodeList=[(node,2,True,doTimeTree),(up[node],child,True,doTimeTree)]
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
					oldTreeLK,oldTreeLK0=calculateTreeLikelihood(tree,root,checkCorrectness=True,separate=True)
					oldTreeLKTime=0
					if doTimeTree:
						oldTreeLKTime=calculateTreeLikelihoodTime(tree,root,mutRate,checkCorrectness=True)
					reCalculateAllGenomeLists(tree,root, checkExistingAreCorrect=True)
					if doTimeTree:
						#print("Before")
						#calculateTreeLikelihoodTimeDebugging(tree,root)
						reCalculateAllGenomeListsTime(tree,root,mutRate, checkExistingAreCorrect=True)
				if aBayesPlusOn and networkOutput:
					alternativePlacements[newNode]=[]
				newRoot2,improvement=traverseTreeForTopologyUpdate(tree,newNode,strictTopologyStopRules=strictTopologyStopRules,allowedFailsTopology=allowedFailsTopology,thresholdLogLKtopology=thresholdLogLKtopology,thresholdTopologyPlacement=thresholdTopologyPlacement)
				if checkEachSPR:
					root=newNode
					while up[root]!=None:
						root=up[root]
					newTreeLK,newTreeLK0=calculateTreeLikelihood(tree,root,separate=True)
					newTreeLKTime=0
					if doTimeTree:
						newTreeLKTime=calculateTreeLikelihoodTime(tree,root,mutRate)
					reCalculateAllGenomeLists(tree,root, checkExistingAreCorrect=True)
					if doTimeTree:
						#print("after:")
						#calculateTreeLikelihoodTimeDebugging(tree,root)
						reCalculateAllGenomeListsTime(tree,root,mutRate, checkExistingAreCorrect=True)
					print("In startTopologyUpdates for node "+str(newNode)+", new LK scores "+str(newTreeLK)+" "+str(newTreeLK0)+" "+str(newTreeLKTime)+" , old LK scores "+str(oldTreeLK)+" "+str(oldTreeLK0)+" "+str(oldTreeLKTime)+" , with improvement "+str(newTreeLK+newTreeLK0+newTreeLKTime)+" - "+str(oldTreeLK+oldTreeLK0+oldTreeLKTime)+" = "+str(newTreeLK+newTreeLK0+newTreeLKTime-(oldTreeLK+oldTreeLK0+oldTreeLKTime))+", was supposed to be "+str(improvement))
					#if (newTreeLK+newTreeLK0+newTreeLKTime-(oldTreeLK+oldTreeLK0+oldTreeLKTime)) < improvement-1.0:
					#print("nDesc0")
					#print(tree.nDesc0)
					if ( (newTreeLK+newTreeLK0+newTreeLKTime-(oldTreeLK+oldTreeLK0+oldTreeLKTime)) < improvement-0.5) or ( (newTreeLK+newTreeLK0+newTreeLKTime-(oldTreeLK+oldTreeLK0+oldTreeLKTime)) > improvement+0.5):
						print("In startTopologyUpdates, LK score of improvement "+str(newTreeLK+newTreeLK0+newTreeLKTime)+" - "+str(oldTreeLK+oldTreeLK0+oldTreeLKTime)+" = "+str(newTreeLK+newTreeLK0+newTreeLKTime-(oldTreeLK+oldTreeLK0+oldTreeLKTime))+" is less than what is supposed to be "+str(improvement))
						print("names")
						print(tree.name)
						print("nDesc0")
						print(tree.nDesc0)
						print("up")
						print(tree.up)
						print("children")
						print(tree.children)
						print("probVectTime")
						print(tree.probVectTime)
						print("dateData")
						print(tree.dateData)
						print("probVectUpRightTime")
						print(tree.probVectUpRightTime)
						print("probVectUpLeftTime")
						print(tree.probVectUpLeftTime)
						print("probVectTotUpTime")
						print(tree.probVectTotUpTime)
						print("dist")
						print(tree.dist)
						print("mutRate")
						print(mutRate)
						raise Exception("exit")
				totalImprovement+=improvement
				if newRoot2!=None:
					newRoot=newRoot2
				numNodes+=1
				if (numNodes%printEvery)==0:
					print("Processed topology for "+str(numNodes)+" nodes.", flush=True)
		print("Topology updates "+str(topologyUpdates[0])+" ; bLen updates "+str(bLenUpdates[0]))
		return newRoot,totalImprovement


# TODO parallelized version of startTopologyUpdates(): takes a number "corNum" in input, and only performs search on the nodes that are assinged to that core.
#traverse the tree (here the input "node" will usually be the root), and for each appropriate dirty node encountered, search the best SPR move re-placing that node.
# attempt an SPR move by cutting the subtree rooted at this dirty node and trying to re-append it elsewhere.
# TODO TODO TODO test
def startTopologyUpdatesParallel(inputTuple):
	tree, startingNode,corNum,strictTopologyStopRules,allowedFailsTopology,thresholdLogLKtopology,thresholdTopologyPlacement,mutRate,errorRateGlobal,mutMatrixGlobal,errorRates,mutMatrices,cumulativeRate,cumulativeErrorRate= inputTuple
	up=tree.up
	children=tree.children
	dirty=tree.dirty
	replacements=tree.replacements
	coreNum=tree.coreNum
	probVectUpRight=tree.probVectUpRight
	probVectUpLeft=tree.probVectUpLeft
	minorSequences=tree.minorSequences
	probVect=tree.probVect
	if doTimeTree:
		probVectUpRightTime=tree.probVectUpRightTime
		probVectUpLeftTime=tree.probVectUpLeftTime
		probVectTime=tree.probVectTime
	mutations=tree.mutations
	dist=tree.dist
	nDesc0=tree.nDesc0
	nodesToVisit=[startingNode]
	proposedMoves=[]
	aBayesPlusOn=False
	if aBayesPlus:
		aBayesPlusOn=True
		SPRTAreporting=[]
	nodesSearched=0

	if usingErrorRate:
		if errorRates!=None:
			totErrorPassed=-cumulativeErrorRate[-1]
		else:
			totErrorPassed=-errorRateGlobal*lRef
	else:
		totErrorPassed=None

	print("Starting SPR search within core "+str(corNum))
	while nodesToVisit:
		node=nodesToVisit.pop()
		for c in children[node]:
			nodesToVisit.append(c)
		if dirty[node] and replacements[node]<=maxReplacements and coreNum[node]==corNum:
			#placement,improvement=traverseTreeForTopologyUpdateParallel(newNode,strictTopologyStopRules=strictTopologyStopRules,allowedFailsTopology=allowedFailsTopology,thresholdLogLKtopology=thresholdLogLKtopology,thresholdTopologyPlacement=thresholdTopologyPlacement)
			placement=None
			branchSupport=None
			listOfBestPlacements=[]
			improvement=0
			# we avoid the root node since it cannot be re-placed with SPR moves
			if up[node]!=None:
				nodesSearched+=1
				#evaluate current placement
				parentNode=up[node]
				if children[parentNode][0]==node:
					child=0
					vectUp=probVectUpRight[parentNode]
					if doTimeTree:
						vectUpTime=probVectUpRightTime[parentNode]
				else:
					child=1
					vectUp=probVectUpLeft[parentNode]
					if doTimeTree:
						vectUpTime=probVectUpLeftTime[parentNode]
				sibling=children[parentNode][1-child]
				if mutations[node]:
					vectUp=passGenomeListThroughBranch(vectUp,mutations[node],dirIsUp=False)
				#score of current tree
				bestCurrenBLen=dist[node]
				isTip=(len(children[node])==0) and (len(minorSequences[node])==0)
				bestCurrentLK=appendProbNode(vectUp,probVect[node],isTip,bestCurrenBLen,errorRateGlobalPassed=errorRateGlobal,mutMatrixGlobalPassed=mutMatrixGlobal,errorRatesGlobal=errorRates,mutMatricesGlobal=mutMatrices,totErrorPassed=totErrorPassed)
				if doTimeTree:
					if up[parentNode]==None:
						bestCurrentLK+=mergeVectorsTime(probVectTime[node],dist[node],probVectTime[sibling],dist[sibling],mutRate,returnLK=True,isUpDown=False)[1]
						bestCurrentLK+=findProbRootTime(probVectTime[parentNode])
					else:
						bestCurrentLK+=appendProbNodeTime(vectUpTime,probVectTime[node],mutRate,bestCurrenBLen)
						if children[up[parentNode]][0]==parentNode:
							vectUpUpTime=probVectUpRightTime[up[parentNode]]
						else:
							vectUpUpTime=probVectUpLeftTime[up[parentNode]]
						originalVectLKTime,originalVectLKTimeCost = mergeVectorsTime(vectUpUpTime,dist[parentNode],probVectTime[sibling],dist[sibling],mutRate,returnLK=True,isUpDown=True)
						bestCurrentLK+=originalVectLKTimeCost
						bestCurrentLK-=appendProbNodeTime(vectUpUpTime,probVectTime[sibling],mutRate,dist[sibling]+dist[parentNode])
				#correct initial placement likelihood by the change in HnZ
				if HnZ:
					parentNode0=up[node]
					while (dist[parentNode0]<=effectivelyNon0BLen) and up[parentNode0]!=None:
						parentNode0=up[parentNode0]
					if dist[node]>effectivelyNon0BLen:
						bestCurrentLK+= getHnZ(nDesc0[parentNode0]) - getHnZ(nDesc0[parentNode0]-1)
					else:
						bestCurrentLK+= getHnZ(nDesc0[parentNode0]) - ( getHnZ(nDesc0[parentNode0]-nDesc0[node]) + getHnZ(nDesc0[node]) )

				topologyUpdated=False
				#in case of HnZ, also try to replace 0-cost subtrees, as they might have better modifers somewhere else
				if ((bestCurrentLK<thresholdTopologyPlacement or dist[node] or HnZ or doTimeTree) and (not doNotImproveTopology)) or ((dist[node] or supportFor0Branches) and aBayesPlusOn) :
				#if (bestCurrentLK<thresholdTopologyPlacement) or HnZ:# or (supportFor0Branches and aBayesPlusOn):
					#now find the best place on the tree where to re-attach the subtree rooted at "node"
					#but to do that we need to consider new vector probabilities after removing the node that we want to replace
					# this is done using findBestParentTopology().
					try:
						bestNodeSoFar , bestLKdiff , bestBranchLengths, listOfBestPlacements, branchSupport, passedProbVect = findBestParentTopology(tree,parentNode,child,bestCurrentLK,bestCurrenBLen,mutRate=mutRate,strictTopologyStopRules=strictTopologyStopRules,allowedFailsTopology=allowedFailsTopology,thresholdLogLKtopology=thresholdLogLKtopology,errorRateGlobalPassed=errorRateGlobal,mutMatrixGlobalPassed=mutMatrixGlobal,errorRatesGlobal=errorRates,mutMatricesGlobal=mutMatrices,cumulativeRateGlobal=cumulativeRate,cumulativeErrorRateGlobal=cumulativeErrorRate,totErrorPassed=totErrorPassed) # ,sequentialSearch=False
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
					except:
						placement=None
				if branchSupport!=None and aBayesPlusOn:
					if networkOutput:
						SPRTAreporting.append((node,branchSupport,listOfBestPlacements))
						#print("added ")
						#print((node,branchSupport,listOfBestPlacements))
						#print("",flush=True)
						#alternativePlacements[node]=listOfBestPlacements
					else:
						SPRTAreporting.append((node,branchSupport,None))
					#support[node]=branchSupport
				if placement!=None and (not doNotImproveTopology):
					proposedMoves.append((node,placement,improvement))
	print("Searched "+str(nodesSearched)+" nodes within core "+str(corNum)+" and found "+str(len(proposedMoves))+" proposed SPR moves")
	if aBayesPlusOn:
		return (proposedMoves,SPRTAreporting)
	else:
		return proposedMoves


if __name__ == "__main__":
	# Given a tree, and a final substitution rate matrix, calculate the likelihood of the tree
	def calculateTreeLikelihood(tree,root,checkCorrectness=False,separate=False):
		up=tree.up
		children=tree.children
		minorSequences=tree.minorSequences
		probVect=tree.probVect
		mutations=tree.mutations
		dist=tree.dist
		nDesc0=tree.nDesc0
		node=root
		#direction 0 means from parent, direction 1 means from a child
		lastNode=None
		direction=0
		totalLK=0.0
		totalLK0=0.0
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
					#adding HnZ contribution
					if HnZ and ((dist[node]>effectivelyNon0BLen) or up[node]==None):
						totalLK0+=getHnZ(nDesc0[node])
						#print("totalLK0 "+str(totalLK0))
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
		if separate:
			return totalLK, totalLK0
		else:
			return totalLK+totalLK0


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




if __name__ == "__main__":
	#print("expectationMaximizationCalculationRates pre")
	#print(errorRateGlobal)
	#Given a tree and its genome lists, calculate mutation counts and waiting times for expectation maximization estimation of substitution rates
	def expectationMaximizationCalculationRates(tree,root,trackMutations=False):
		#print("expectationMaximizationCalculationRates")
		#print(errorRateGlobal)
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
			if usingErrorRate:
				errors=[]
				for i in range(len(up)):
					errors.append([])
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


	if doTimeTree:
		print("Sorting samples based on dates", flush=True)
		distances=sortSamplesByDate(dates,data,samples=data.keys(),samplesInInitialTree=namesInTreeDict,forgetData=False)
	else:
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
		if doTimeTree:
			tree1.probVectTime[-1]=dates[firstSample[1]]
			tree1.dateData[-1]=dates[firstSample[1]]
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
	print("errorRateGlobal now "+str(errorRateGlobal))

	# initial EM round to estimate rate variation etc on the initial tree
	if (numSamples>1) and (model!="JC" or ((numSamples>=minNumSamplesForRateVar) and useRateVariation ) or ((numSamples>=minNumSamplesForErrorModel) and usingErrorRate)):
		start=time()
		mutMatrixGlobal, siteRatesEst, errorRateGlobalEst, errorRatesEst = expectationMaximizationCalculationRates(tree,t1)
		print("EM from initial tree terminated, using rate variation "+str(useRateVariation)+", using error rates "+str(usingErrorRate)+".  ")
		updateMutMatrices(mutMatrixGlobal,siteRates=siteRatesEst)
		if usingErrorRate:
			errorRates=errorRatesEst
			errorRateGlobal=errorRateGlobalEst
			updateErrorRates(errorRateGlobal,errorRates=errorRates)
		#if useRateVariation:
		#	print(mutMatrices[0])
		reCalculateAllGenomeLists(tree,t1)
		newLk=calculateTreeLikelihood(tree,t1)
		print("LK after first EM: "+str(newLk))
		if usingErrorRate and (estimateErrorRate or estimateSiteSpecificErrorRate):
			oldLK=float("-inf")
			numEMsteps=0
			while (newLk-oldLK>1.0) and numEMsteps<20:
				if not doNotOptimiseBLengths:
					setAllDirty(tree,t1)
					improvement=traverseTreeToOptimizeBranchLengths(tree,t1)
					if doTimeTree:
						reCalculateAllGenomeListsTime(tree,t1,mutRate)
						print("Time LK: "+str(calculateTreeLikelihoodTime(tree,t1,mutRate)))
				reCalculateAllGenomeLists(tree,t1)
				newLkBranch=calculateTreeLikelihood(tree,t1)
				print("Updated "+str(improvement)+" branch lengths leading to LK "+str(newLkBranch))
				mutMatrixGlobal, siteRates, errorRateGlobalEst, errorRatesEst = expectationMaximizationCalculationRates(tree,t1)
				updateMutMatrices(mutMatrixGlobal,siteRates=siteRates)
				if usingErrorRate:
					errorRates=errorRatesEst
					errorRateGlobal=errorRateGlobalEst
					updateErrorRates(errorRateGlobal,errorRates=errorRates)
				reCalculateAllGenomeLists(tree,t1)
				oldLK=newLk
				newLk=calculateTreeLikelihood(tree,t1)
				print("New LK step "+str(numEMsteps)+": "+str(newLk))
				numEMsteps+=1

		timeRecalculation=time()-start
		print("Time to run initial tree EM estimation: "+str(timeRecalculation))

	# This function was moved upward to be used in Lineage Assignments by Lineage reference genomes
	#generate the string corresponding to a line of the tsv file for use in Taxonium.
	# If representative node is 0-dist, then its support is also the support of all represented nodes. If not, can we assume that the support of the represented nodes is 1.
	# Includes support of represented nodes only if supportFor0Branches is true, otherwise smpty string.
	def tsvForNode(tree,node,name,featureList,namesInTree,identicalTo=""):
		dist=tree.dist
		minorSequences=tree.minorSequences
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
							#TODO TODO TODO
							#TODO TODO TODO I added this "if" here to make sure that supportTo refers to the clade containg identical genomes, rather than to one of those genomes itself.
							if len(minorSequences[feature[node][iNode][0]])>0:
								stringList.append(namesInTree[tree.name[feature[node][iNode][0]]]+"_MinorSeqsClade:"+str(feature[node][iNode][1]))
							else:
								stringList.append(namesInTree[tree.name[feature[node][iNode][0]]]+":"+str(feature[node][iNode][1]))
							if iNode<(len(feature[node])-1):
								stringList.append(",")
					# use this column to highlight which lineages could be placed (with probability above threshold) on the branch above the current node - used to highlight alternative placements of a given lineage on the tree.
					elif feat == "supportToLineages" and identicalTo == "":
						for iNode in range(len(feature[node])):
							stringList.append(feature[node][iNode][0] + ":" + str(
								feature[node][iNode][1]))
							if iNode < (len(feature[node]) - 1):
								stringList.append(";")
					elif feat=="lineageParent":
						stringList.append(feature[node])
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
		return "".join(stringList)

	# seek placements for a chunk that contains a subset of lineage reference genomes
	#NHAN
	def process_chunk(job_id, start, end, lineageRefNames, tree, t1, lineageRefData):
		chunk_output = []
		numSamples = 0
		totalSamples = end - start
		# Process a chunk of tasks and return results
		for lineageRefName in lineageRefNames[start:end]:
			# extract the partial of the lineage reference genome
			newPartials = probVectTerminalNode(lineageRefData[lineageRefName], None, None)

			# find the best placement for the lineage reference genome
			possiblePlacements, bestPlacementTotalLh = findBestParentForNewSample(tree, t1, newPartials, numSamples, computePlacementSupportOnly=True)

			# finetune the placement position
			if len(possiblePlacements):
				# sort possiblePlacements by placement's support descending
				sortedPlacements = sorted(possiblePlacements, key=lambda x: x[1], reverse=True)

			else:
				print(f"PossiblePlacements for {lineageRefName} is empty - probably many very low probability placements were found",flush=True)
				#raise Exception("exit")
				sortedPlacements = []
				bestPlacementTotalLh = []

			chunk_output.append((lineageRefName, sortedPlacements, bestPlacementTotalLh))

			# update progress
			numSamples += 1
			if numSamples % 50 == 0 or numSamples == totalSamples:
				print(f"JobID {job_id} processed {numSamples}/{totalSamples}",flush=True)

		return chunk_output


	# extract list of mutations from two genomes
	def extractMutations(probVect1, probVect2):
		mutations_list = []
		indexEntry1, indexEntry2, pos = 0, 0, 0
		entry1 = probVect1[indexEntry1]
		entry2 = probVect2[indexEntry2]
		while True:
			if entry1[0] != entry2[0] and entry1[0] < 5 and entry2[0] < 5:
				if entry1[0] == 4:
					mutations_list.append((entry2[1], entry2[0], pos + 1, None))
				elif entry2[0] == 4:
					mutations_list.append((entry1[0], entry1[1], pos + 1, None))
				else:
					mutations_list.append((entry1[0], entry2[0], pos + 1, None))
				pos += 1
			else:
				# handle special case when entry1 is O
				if entry1[0] != entry2[0] and entry1[0] == 6 and entry2[0] < 5:
					# extract the nucleotide of entry2
					entry2_nuc = entry2[0]
					if entry2[0] == 4:
						entry2_nuc = entry1[1]
					
					# only consider cases where the probability of the nucleotide in the sample is lower than threshMutation
					if entry1[-1][entry2_nuc] < threshMutation:
						mutations_list.append((entry1[0], entry2_nuc, pos + 1, entry1[-1]))

				if (entry1[0] == 4 or entry1[0] == 5) and (entry2[0] == 4 or entry2[0] == 5):
					pos = min(entry1[1], entry2[1])
				else:
					pos += 1

			if pos == lRef:
				break
			if entry1[0] < 4 or entry1[0] == 6:
				indexEntry1 += 1
				entry1 = probVect1[indexEntry1]
			elif pos == entry1[1]:
				indexEntry1 += 1
				entry1 = probVect1[indexEntry1]
			if entry2[0] < 4 or entry2[0] == 6:
				indexEntry2 += 1
				entry2 = probVect2[indexEntry2]
			elif pos == entry2[1]:
				indexEntry2 += 1
				entry2 = probVect2[indexEntry2]
		return mutations_list

	# seek placements for lineage reference genomes
	#NHAN
	def seekPlacementOfLineageRefs(tree, t1, lineageRefData, numCores, findPlacementOnly):
		#dist = tree.dist
		#up = tree.up
		probVect = tree.probVect
		# create a map from a lineage to its possible placements
		tree.lineagePlacements = {}
		lineageRefNames = list(lineageRefData.keys())
		chunk_size = (len(lineageRefNames) + numCores - 1) // numCores  # Ensures rounding up
		chunks = [(i, min(i + chunk_size, len(lineageRefNames))) for i in range(0, len(lineageRefNames), chunk_size)]

		# Parallelize over chunks
		results = Parallel(n_jobs=numCores)(
			delayed(process_chunk)(job_id, start, end, lineageRefNames, tree, t1, lineageRefData) for
			job_id, (start, end) in enumerate(chunks)
		)

		for chunk_output in results:
			for lineageRefName, sortedPlacements, bestPlacementTotalLh in chunk_output:
				# extract the best placement (with the highest support)
				if sortedPlacements:
					selectedPlacement = sortedPlacements[0][0]
				# conduct lineage assignment if needed
				if not findPlacementOnly:
					lineageRootPosition = None
					# extract support and optimized blengths of the best placement
					if sortedPlacements:
						selectedPlacementSupport = sortedPlacements[0][1]
						topBlength, bottomBlength, appendingBlength = sortedPlacements[0][2]

						# append the lineage assignment into the selected node
						if appendingBlength <= lineageRefsThresh and selectedPlacementSupport >= lineageRefsSupportThresh:
							# if topBlength == 0, we already record the parent instead of the original placement, so no further processing needed
							# if not topBlength and up[selectedPlacement]:
							#	selectedPlacement = up[selectedPlacement]
							# traverse upward to the top of the polytomy
							#	while (dist[selectedPlacement] <= effectivelyNon0BLen) and (up[selectedPlacement] != None):
							#		selectedPlacement = up[selectedPlacement]
							tree.lineageAssignments[selectedPlacement].append([lineageRefName, bottomBlength])
							lineageRootPosition = selectedPlacement

					# update the list of possible placements for this lineage
					tree.lineagePlacements[lineageRefName] = (sortedPlacements, lineageRootPosition)

				# otherwise, extract the list of mutations that separate the sample from the placement
				else:
					if sortedPlacements:
						# extract sample genome
						samplePartials = probVectTerminalNode(lineageRefData[lineageRefName], None, None)

						# extract list of mutations
						mutations_list = extractMutations(bestPlacementTotalLh, samplePartials)

						# update the list of possible placements for this sample
						tree.lineagePlacements[lineageRefName] = (sortedPlacements, mutations_list)
					else:
						tree.lineagePlacements[lineageRefName] = (sortedPlacements, None)

				# delete lineage genome that is already processed
				lineageRefData[lineageRefName] = None

		if not findPlacementOnly:
			# a node may be assigned multiple lineages
			for node in range(len(tree.lineageAssignments)):
				lineageAssignments = tree.lineageAssignments[node]

				# assign the list of lineages to this node
				if len(lineageAssignments) > 0:
					# when a node is descendant of multiple references,
					# assign all these lineages to this node
					if allowMultiLineagesPerNode:
						tree.lineage[node] = "/".join(lineageRefName for lineageRefName, _ in lineageAssignments)
					# otherwise, we assign the phylogenetically closest one to it
					else:
						closestLineage = lineageAssignments[0][0]
						closetDistance = lineageAssignments[0][1]
						for i in range(1, len(lineageAssignments)):
							if lineageAssignments[i][1] < closetDistance:
								closestLineage = lineageAssignments[i][0]
								closetDistance = lineageAssignments[i][1]
						tree.lineage[node] = closestLineage

		return tree


	# Annotate nodes by their lineage assignments
	#NHAN
	def annotateLineageAssignments(tree, root):
		children = tree.children
		lineages = tree.lineage
		lineagesParent = tree.lineageParent

		# traverse the tree and annotate nodes with their lineage assignments
		# start from the root
		nodesToVisit = []
		# If root is not assigned a lineage
		# set lineages[root] = "-" instead of None
		if not lineages[root]:
			lineages[root] = "-"
		lineagesParent[root] = "-"
		for child in children[root]:
			nodesToVisit.append((child, lineages[root]))
		while nodesToVisit:
			node, lineage = nodesToVisit.pop()

			# record the lineage of the parent node
			lineagesParent[node] = lineage

			# if this node has NOT been already assigned any lineage,
			# inherit the assignment from its parent node
			if not lineages[node]:
				lineages[node] = lineage

			# traverse downward to children nodes
			for child in children[node]:
				nodesToVisit.append((child, lineages[node]))

		# synchronize lineages to tree
		tree.lineage = lineages
		tree.lineageParent = lineagesParent

		# return the updated tree
		return tree

	#NHAN
	def defineSupportedToLineages(tree):
		# init supportToLineages
		numNodes = len(tree.up)
		tree.supportToLineages = [[] for _ in range(numNodes)]
		lineagePlacements = tree.lineagePlacements
		for key in lineagePlacements:
			plausiblePlacements, lineageRootPosition = lineagePlacements[key]
			for placement, support, optimizedBlengths in plausiblePlacements:
				topBlength, bottomBlength, appendingBlength = optimizedBlengths
				if appendingBlength <= lineageRefsThresh:
					tree.supportToLineages[placement].append([key, support])
		return tree


	# Write lineage assignments to output file
	#NHAN
	def outputLineageAssignments(outputFile, tree, root):
		tree = defineSupportedToLineages(tree)
		# ------------ write TSV file ------------------
		giveInternalNodeNames(tree, t1, namesInTree=namesInTree, replaceNames=False)
		file = open(outputFile + "_metaData_lineageAssignment.tsv", "w")

		children = tree.children
		up = tree.up
		name = tree.name
		minorSequences = tree.minorSequences
		featureNames = {}
		featureNames['lineage'] = 'lineage'
		featureNames['supportToLineages'] = 'supportToLineages'
		featureNames['lineageParent'] = 'lineageParent'
		featureList = list(featureNames.keys())
		file.write("strain" + "\t" + "collapsedTo")
		for feat in featureList:
			file.write("\t" + featureNames[feat])
		file.write("\n")
		# now write to file the features for each node of the tree.
		nextNode = root
		direction = 0
		numLeaves = 0
		while nextNode != None:
			if children[nextNode]:
				if direction == 0:
					nextNode = children[nextNode][0]
				elif direction == 1:
					nextNode = children[nextNode][1]
					direction = 0
				else:
					file.write(tsvForNode(tree, nextNode, namesInTree[name[nextNode]], featureList, namesInTree))
					if up[nextNode] != None:
						if children[up[nextNode]][0] == nextNode:
							direction = 1
						else:
							direction = 2
					nextNode = up[nextNode]
			else:
				numLeaves += (1 + len(minorSequences[nextNode]))
				if len(minorSequences[nextNode]) > 0:
					file.write(tsvForNode(tree, nextNode, namesInTree[name[nextNode]], featureList, namesInTree,
										  identicalTo=namesInTree[name[nextNode]] + "_MinorSeqsClade"))

					for s2 in minorSequences[nextNode]:
						file.write(tsvForNode(tree, nextNode, namesInTree[s2], featureList, namesInTree,
											  identicalTo=namesInTree[name[nextNode]] + "_MinorSeqsClade"))

					file.write(
						tsvForNode(tree, nextNode, namesInTree[name[nextNode]] + "_MinorSeqsClade", featureList,
								   namesInTree))
				else:
					file.write(tsvForNode(tree, nextNode, namesInTree[name[nextNode]], featureList, namesInTree))
				if up[nextNode] != None:
					if children[up[nextNode]][0] == nextNode:
						direction = 1
					else:
						direction = 2
				nextNode = up[nextNode]

		# close the output file
		file.close()

		print(f"Output lineage assignments at {outputFile}_metaData_lineageAssignment.tsv.",flush=True)

		# write TSV mapping from lineage to its possible placements
		file = open(outputFile + "_metaData_lineagePlacements.tsv", "w")
		lineagePlacements = tree.lineagePlacements
		file.write("lineage\tplacements\toptimizedBlengths\tlineageRootPosition\n")
		for key in lineagePlacements:
			placementStrVec = []
			placementBlengthsVec = []
			plausiblePlacements, lineageRootPosition = lineagePlacements[key]
			for placement, support, optimizedBlengths in plausiblePlacements:
				placementStrVec.append(f"{namesInTree[name[placement]]}:{str(support)}")
				blengthsVec = []
				for blength in optimizedBlengths:
					if blength:
						blengthsVec.append(str(blength))
					else:
						blengthsVec.append("0")
				blengthsStr = "/".join(blengthsVec)
				placementBlengthsVec.append(f"{namesInTree[name[placement]]}:({blengthsStr})")

			# extract the root position of the lineage (if the current lineage is considered as an ancestral of a subtree)
			lineageRootPositionStr = "-"
			if lineageRootPosition != None:
				lineageRootPositionStr = namesInTree[name[lineageRootPosition]]

			placementStr = ";".join(placementStrVec)
			placementBlengthsStr = ";".join(placementBlengthsVec)
			file.write(
				key + "\t" + placementStr + "\t" + placementBlengthsStr + "\t" + lineageRootPositionStr + "\n")

		# close the output file
		file.close()

		print(f"Output a map from lineages to their placements at {outputFile}_metaData_lineagePlacements.tsv.",flush=True)
		# ------------ end of write TSV file ------------------

		# ------------ write Nexus treefile ------------------
		newickString = createNewick(tree, root, binary=binaryTree, namesInTree=namesInTree)
		file = open(outputFile + "_lineageAssignment.tree", "w")
		file.write("#NEXUS\nbegin taxa;\n	dimensions ntax=" + str(len(namesInTree)) + ";\n	taxlabels\n")
		for name in namesInTree:
			file.write("	" + name + "\n")
		file.write(";\nend;\n\nbegin trees;\n	tree TREE1 = [&R] ")
		file.write(newickString)
		file.write("\nend;\n")
		file.close()
		print(f"Output Nexus tree with lineage assignments at {outputFile}_lineageAssignment.tree.",flush=True)
		# ------------ end of write Nexus treefile ------------------

		# ------------ write Newick treefile ------------------
		newickString = createNewick(tree, root, binary=binaryTree, namesInTree=namesInTree, estimateMAT=False,
									networkOutput=False, aBayesPlusOn=False,
									performLineageAssignmentByRefPlacement=False)
		file = open(outputFile + "_updatedBlengths.tree", "w")
		file.write(newickString)
		file.close()
		print(f"Output Newick tree with updated branch lengths at {outputFile}_updatedBlengths.tree.",flush=True)
		# ------------ end of write Newick treefile ------------------
		# return success
		return tree

	# Write sample placements to output file
	# NHAN
	def outputSamplePlacements(outputFile, tree, root):
		nucletides = "ACGTRNO"
		# write TSV mapping from lineage to its possible placements
		giveInternalNodeNames(tree, t1, namesInTree=namesInTree, replaceNames=False)
		name = tree.name
		file = open(outputFile + "_metaData_samplePlacements.tsv", "w")
		lineagePlacements = tree.lineagePlacements
		# sample: Names of samples in S
		# possiblePlacements: a set of possible placements for each sample. Each placement is presented as <placementNode>:<support>.
		# optimizedBlengths: the corresponding set of optimized branch lengths of possible placements. Each placement contains a set of three branch lengths and is presented as: <placementNode>:(<topBlength>/<bottomBlength>/<sampleBlength>).
		# mutations: a list of mutations that separates the sample from the most similar genome (i.e., placement's partials). Each mutation is present as <stateAtPlacement><Position><stateAtSample>
		file.write("sample\tplacements\toptimizedBlengths\tmutations\n")
		for key in lineagePlacements:
			placementStrVec = []
			placementBlengthsVec = []
			plausiblePlacements, mutations_list = lineagePlacements[key]
			for placement, support, optimizedBlengths in plausiblePlacements:
				placementStrVec.append(f"{namesInTree[name[placement]]}:{str(support)}")
				blengthsVec = []
				for blength in optimizedBlengths:
					if blength:
						blengthsVec.append(str(blength))
					else:
						blengthsVec.append("0")
				blengthsStr = "/".join(blengthsVec)
				placementBlengthsVec.append(f"{namesInTree[name[placement]]}:({blengthsStr})")

			placementStr = ";".join(placementStrVec)
			placementBlengthsStr = ";".join(placementBlengthsVec)

			# convert list of mutations into a string
			mutationStrVec = []
			if mutations_list!=None:
				for from_state, to_state, position, prob_nuc in mutations_list:
					# case when from_state = O => return O(prob_A/prob_C/prob_G/prob_T)<pos><to_state>
					if from_state == 6:
						# normalize the probabilities over all nucleotides
						prob_nuc_vec = []
						total_prob = sum(prob_nuc)
						for i in range(len(prob_nuc)):
							prob_nuc_vec.append(f"{prob_nuc[i]/total_prob:.6f}")
						prob_nuc_str = "/".join(prob_nuc_vec)
						mutationStrVec.append(f"{nucletides[from_state]}({prob_nuc_str}){position}{nucletides[to_state]}")
					# otherwise, => return <from_state><pos><to_state>
					else:
						mutationStrVec.append(f"{nucletides[from_state]}{position}{nucletides[to_state]}")
			mutationStr = ";".join(mutationStrVec)

			file.write(key + "\t" + placementStr + "\t" + placementBlengthsStr + "\t" + mutationStr + "\n")

		# close the output file
		file.close()

		print(f"Output a map from lineages to their placements at {outputFile}_metaData_metaData_samplePlacements.tsv.",flush=True)

		# ------------ write Newick treefile ------------------
		newickString = createNewick(tree, root, binary=binaryTree, namesInTree=namesInTree, estimateMAT=False,
									networkOutput=False, aBayesPlusOn=False,
									performLineageAssignmentByRefPlacement=False)
		file = open(outputFile + "_updatedBlengths.tree", "w")
		file.write(newickString)
		file.close()
		print(f"Output Newick tree with updated branch lengths at {outputFile}_updatedBlengths.tree.",flush=True)
		# ------------ end of write Newick treefile ------------------
		# return success
		return tree


	# Process linage assignments by reference genomes
	# Input: tree and lineage reference genomes
	# Output: assignments of nodes (tip and internal nodes) to lineages
	# 1. find a placement for each lineage reference
	# 2. locate the subtree rooted at the placement of each lineage reference; assign all children of that subtree to that lineage
	#NHAN
	def assignLineageByReferencePlacement(tree, t1, lineageRefData, numCores):
		numNodes = len(tree.up)
		tree.lineageAssignments = [[] for _ in range(numNodes)]
		tree.lineage = [None] * numNodes
		tree.lineageParent = [None] * numNodes # the lineage of the parent node of this node
		tree.lineages = [None] * numNodes  # don't use but need to add to reuse other functions

		# 1. Find a placement for each lineage reference
		tree = seekPlacementOfLineageRefs(tree, t1, lineageRefData, numCores, findPlacementOnly = False)

		# 2. Annotate nodes by their lineage assignments
		tree = annotateLineageAssignments(tree, t1)

		# Write lineage assignments to output file
		outputLineageAssignments(outputFile, tree, t1)

		# terminate the program
		exit(0)

	# Find placements for new samples
	# Input: tree and new sample genomes
	# Output: possible placements of each new sample
	# NHAN
	def findPlacementsForSamples(tree, t1, distances, numCores):
		# extract sampleGenomes from distances
		sampleGenomes = {}
		while distances:
			sample=distances.pop()[1]
			sampleGenomes[sample] = data[sample]

		# 1. Find a placement for each new sample
		tree = seekPlacementOfLineageRefs(tree, t1, sampleGenomes, numCores, findPlacementOnly = True)

		# Write sample placements to output file
		outputSamplePlacements(outputFile, tree, t1)

		# terminate the program
		exit(0)

	# Process linage assignments by reference genomes
	#NHAN
	if performLineageAssignmentByRefPlacement:
		assignLineageByReferencePlacement(tree, t1, lineageRefData, numCores)

	# If users only want to find placements for new samples
	if findSamplePlacements:
		from joblib import Parallel, delayed
		findPlacementsForSamples(tree, t1, distances, numCores)

	# initial EM round to estimate the time-scaled mutation rate
	if doTimeTree and (numSamples>=minNumSamplesForMutRate):
		start=time()
		oldLK=calculateTreeLikelihoodTime(tree,t1,mutRate)
		print("pre-EM mutation rate "+str(mutRate)+" time LK before first EM: "+str(oldLK))
		counts, waitingTimes, mutRate = expectationMaximizationCalculationRatesTime(tree,t1,mutRate)
		poissonCoeff=[[1.0]]
		reCalculateAllGenomeListsTime(tree,t1,mutRate)
		newLk=calculateTreeLikelihoodTime(tree,t1,mutRate)
		print("EM from initial tree terminated, using mutation rate "+str(mutRate)+" time LK: "+str(newLk))
		numEMsteps=0
		while (newLk-oldLK>0.1) and numEMsteps<20:
			counts, waitingTimes, mutRate = expectationMaximizationCalculationRatesTime(tree,t1,mutRate)
			poissonCoeff=[[1.0]]
			reCalculateAllGenomeListsTime(tree,t1,mutRate)
			oldLK=newLk
			newLk=calculateTreeLikelihoodTime(tree,t1,mutRate)
			numEMsteps+=1
		timeRecalculation=time()-start
		print("New time LK step "+str(numEMsteps)+" mutRate "+str(mutRate)+": "+str(newLk))
		print("Time to run initial tree time EM estimation: "+str(timeRecalculation))

	#Place input samples to create an initial tree (or extend the input tree).
	timeFinding=0.0
	timePlacing=0.0
	lastUpdateNumSamples=numSamples
	lastUpdateNumSamplesTime=numSamples
	missingDateWarned=False
	if not doNotPlaceNewSamples:
		while distances:
			d=distances.pop()
			sample=d[1]
			namesInTree.append(sample)
			newPartials=probVectTerminalNode(data[sample],None,None)
			if doTimeTree:
				if sample in dates:
					newPartialsTime=dates[sample]
				else:
					if not missingDateWarned:
						print("WARNING Some samples have no date data (e.g. "+sample+"), they will be considered as having no date information.")
						missingDateWarned=True
					newPartialsTime=None
			else:
				newPartialsTime=None
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
				mutMatrixGlobal, siteRates, errorRateGlobalEst, errorRatesEst = expectationMaximizationCalculationRates(tree,t1)
				updateMutMatrices(mutMatrixGlobal,siteRates=siteRates)
				if usingErrorRate:
					errorRates=errorRatesEst
					errorRateGlobal=errorRateGlobalEst
					updateErrorRates(errorRateGlobal,errorRates=errorRates)
				reCalculateAllGenomeLists(tree,t1)
				improvement=traverseTreeToOptimizeBranchLengths(tree,t1)
				reCalculateAllGenomeLists(tree,t1)
				if doTimeTree:
					reCalculateAllGenomeListsTime(tree,t1,mutRate)
					print("Time LK: "+str(calculateTreeLikelihoodTime(tree,t1,mutRate)))
				print(" EM to update parameters during initial placement terminated, time taken: "+str(time()-start))

			if (doTimeTree and (numSamples>minNumSamplesForMutRate) and (numSamples>2*lastUpdateNumSamplesTime)):
				start=time()
				lastUpdateNumSamplesTime=numSamples
				reCalculateAllGenomeListsTime(tree,t1,mutRate)
				counts, waitingTimes, mutRate = expectationMaximizationCalculationRatesTime(tree,t1,mutRate)
				reCalculateAllGenomeListsTime(tree,t1,mutRate)
				print(" EM to update mutRate during initial placement terminated, new mutRate "+str(mutRate)+" time taken: "+str(time()-start))

			start=time()
			bestNode , bestScore, bestBranchLengths, bestPassedVect = findBestParentForNewSample(tree,t1,newPartials,numSamples,computePlacementSupportOnly=False, diffsTime=newPartialsTime)
			timeFinding+=(time()-start)
			if bestBranchLengths!=None:
				start=time()
				newRoot=placeSampleOnTree(tree,bestNode,bestPassedVect,numSamples,bestScore, bestBranchLengths[0], bestBranchLengths[1], bestBranchLengths[2],pseudoMutCounts,newPartialsTime=newPartialsTime)
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
		print("Re-calculating the nDesc0")
		calculateNDesc0(tree,t1,checkExisting=False)

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
		mutMatrixGlobal, siteRates, errorRateGlobalEst, errorRates = expectationMaximizationCalculationRates(tree,t1)
		if estimateErrorRate:
			errorRateGlobal=errorRateGlobalEst
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
		if not doNotOptimiseBLengths:
			improvement=traverseTreeToOptimizeBranchLengths(tree,t1)
			reCalculateAllGenomeLists(tree,t1)
			newLk=calculateTreeLikelihood(tree,t1)
			print("Tree LK after branch length optimization: "+str(newLk))

			if doTimeTree:
				reCalculateAllGenomeListsTime(tree,t1,mutRate)
				print("Time LK: "+str(calculateTreeLikelihoodTime(tree,t1,mutRate)))


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
			mutMatrixGlobal, siteRates, errorRateGlobalEst, errorRates = expectationMaximizationCalculationRates(tree,t1)
			updateMutMatrices(mutMatrixGlobal,siteRates=siteRates)
			if usingErrorRate:
				errorRateGlobal=errorRateGlobalEst
				updateErrorRates(errorRateGlobal,errorRates=errorRates)
			reCalculateAllGenomeLists(tree,t1)
			newLk=calculateTreeLikelihood(tree,t1)
			print("Tree LK after EM: "+str(newLk))
			if not doNotOptimiseBLengths:
				setAllDirty(tree,t1)
				improvement=traverseTreeToOptimizeBranchLengths(tree,t1)
				reCalculateAllGenomeLists(tree,t1)
				newLk=calculateTreeLikelihood(tree,t1)
				print("Tree LK after branch length optimization: "+str(newLk))
			if estimateErrorRate or estimateSiteSpecificErrorRate:
				oldLK=float("-inf")
				numEMsteps=0
				while (newLk-oldLK>1.0) and numEMsteps<20:
					if not doNotOptimiseBLengths:
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

		if not doNotOptimiseBLengths:
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
			#reCalculateAllGenomeLists(tree,t1)
			newLk=calculateTreeLikelihood(tree,t1)
			print("branch length finalization subround "+str(subRound+1)+", number of changes "+str(improvement)+" final LK: "+str(newLk))
			timeForBranchOptimization=(time()-start)
			print("Time for updating branch lengths: "+str(timeForBranchOptimization))

	if HnZ:
		print("Re-calculating the nDesc0")
		calculateNDesc0(tree,t1,checkExisting=False)


	# EM round to estimate the time-scaled mutation rate
	if doTimeTree:
		start=time()
		reCalculateAllGenomeListsTime(tree,t1,mutRate)
		oldLK=calculateTreeLikelihoodTime(tree,t1,mutRate)
		print("pre-EM mutation rate "+str(mutRate)+" time LK before post-initial-tree EM: "+str(oldLK))
		counts, waitingTimes, mutRate = expectationMaximizationCalculationRatesTime(tree,t1,mutRate)
		poissonCoeff=[[1.0]]
		reCalculateAllGenomeListsTime(tree,t1,mutRate)
		newLk=calculateTreeLikelihoodTime(tree,t1,mutRate)
		print("EM post-initial-tree terminated, using mutation rate "+str(mutRate)+" time LK: "+str(newLk))
		numEMsteps=0
		while (newLk-oldLK>0.1) and numEMsteps<20:
			counts, waitingTimes, mutRate = expectationMaximizationCalculationRatesTime(tree,t1,mutRate)
			poissonCoeff=[[1.0]]
			reCalculateAllGenomeListsTime(tree,t1,mutRate)
			oldLK=newLk
			newLk=calculateTreeLikelihoodTime(tree,t1,mutRate)
			numEMsteps+=1
		timeRecalculation=time()-start
		print("New time LK step "+str(numEMsteps)+" mutRate "+str(mutRate)+": "+str(newLk))
		print("Time to run post-initial-tree time EM estimation: "+str(timeRecalculation))


	#find best tree root
	if not doNotReroot:
		print("Looking for possible better root", flush=True)
		print("LK before looking for root: "+str(calculateTreeLikelihood(tree,t1)))
		if doTimeTree:
			print("Time LK before looking for root: "+str(calculateTreeLikelihoodTime(tree,t1,mutRate)))
		newT1=findBestRoot(tree,t1,strictTopologyStopRules=strictTopologyStopRules,allowedFailsTopology=allowedFailsTopology,thresholdLogLKtopology=thresholdLogLKtopology,aBayesPlusOn=aBayesPlus)
		if newT1!=t1:
			print("Better root found")
			t1=newT1

			if model!="JC" or rateVariation or estimateErrorRate or estimateSiteSpecificErrorRate:
				newLk=calculateTreeLikelihood(tree,t1)
				print("Tree LK before EM: "+str(newLk))
				mutMatrixGlobal, siteRates, errorRateGlobalEst, errorRates = expectationMaximizationCalculationRates(tree,t1)
				if estimateErrorRate:
					errorRateGlobal=errorRateGlobalEst
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
				print("Tree LK after EM: "+str(newLk))
			if not doNotOptimiseBLengths:
				improvement=traverseTreeToOptimizeBranchLengths(tree,t1)
				reCalculateAllGenomeLists(tree,t1)
				newLk=calculateTreeLikelihood(tree,t1)
				print("Tree LK after branch length optimization: "+str(newLk))

			print("Looking a second time for possible better root", flush=True)
			newT1=findBestRoot(tree,t1,strictTopologyStopRules=strictTopologyStopRules,allowedFailsTopology=allowedFailsTopology,thresholdLogLKtopology=thresholdLogLKtopology,aBayesPlusOn=aBayesPlus)
			if newT1!=t1:
				print("Better root found again")
				t1=newT1
				reCalculateAllGenomeLists(tree,t1)
				newLk=calculateTreeLikelihood(tree,t1)
				print("Tree LK after further new root: "+str(newLk))
			
			if doTimeTree:
				reCalculateAllGenomeListsTime(tree,t1,mutRate)
				oldLK=calculateTreeLikelihoodTime(tree,t1,mutRate)
				print("pre-EM mutation rate "+str(mutRate)+" time LK before post-initial-tree EM: "+str(oldLK))
				counts, waitingTimes, mutRate = expectationMaximizationCalculationRatesTime(tree,t1,mutRate)
				poissonCoeff=[[1.0]]
				reCalculateAllGenomeListsTime(tree,t1,mutRate)
				newLk=calculateTreeLikelihoodTime(tree,t1,mutRate)
				print("EM terminated, using mutation rate "+str(mutRate)+" time LK: "+str(newLk))


	if (writeTreesToFileEveryTheseSteps>0):
		currentRoot=t1
		intermediateTreesFile.write("Topology 0\n")
		intermediateTreesFile.write(createNewick(tree,currentRoot,binary=binaryTree,namesInTree=namesInTree)+"\n")
	if (writeLKsToFileEveryTheseSteps>0):
		currentRoot=t1
		totalLK=calculateTreeLikelihood(tree,currentRoot)
		if doTimeTree:
			totalLK+=calculateTreeLikelihoodTime(tree,currentRoot,mutRate)
		intermediateLKsFile.write("Topology 0, LK: "+str(totalLK)+"\n")


	print(str(len(namesInTree))+" named samples in the tree", flush=True)
	internalNodeNamesGiven=False
	giveInternalNodeNames(tree,t1,namesInTree=namesInTree,replaceNames=False)
	internalNodeNamesGiven=True
	print(str(len(namesInTree))+" named nodes in the tree after assigning internal node names", flush=True)

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
						file.write(tsvForNode(tree,nextNode,namesInTree[name[nextNode]]+"_MinorSeqsClade",featureList,namesInTree))
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
		node=root
		#direction 0 means from parent, direction 1 means from a child
		lastNode=None
		direction=0
		numNodes=0
		numDirty=0
		while node!=None:
			if direction==0:
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
		if HnZ:
			print("Re-calculating the nDesc0")
			calculateNDesc0(tree,t1,checkExisting=True)
		if not doNotOptimiseBLengths:
			newLk=calculateTreeLikelihood(tree,t1)
			print("Preliminary branch length optimization from LK: "+str(newLk))
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
		if HnZ:
			print("Re-calculating the nDesc0")
			calculateNDesc0(tree,t1,checkExisting=True)

		#now update the topology
		start=time()
		setAllDirty(tree,t1)
		print("Starting tolopogical improvement round "+str(nRound+1)+" with allowed fails "+str(threshNums[nRound])+" and LK threshold "+str(threshValues[nRound]), flush=True)
		reCalculateAllGenomeLists(tree,t1)
		preLK=calculateTreeLikelihood(tree,t1)
		print("Likelihood before SPR moves: "+str(preLK), flush=True)
		if doTimeTree:
			reCalculateAllGenomeListsTime(tree,t1,mutRate)
			print("Time LK: "+str(calculateTreeLikelihoodTime(tree,t1,mutRate)))

		if parallelize:
			#assign numbers to nodes so that each core will only investigate re-placement of its own nodes
			if nRound==0:
				assignCoreNumbers(tree,t1,numCores)
			paralleleInputs=[]
			for i in range(numCores):
				paralleleInputs.append((tree,t1,i,threshStricts[nRound],threshNums[nRound],threshValues[nRound],threshPlaces[nRound],mutRate,   errorRateGlobal,mutMatrixGlobal,errorRates,mutMatrices,cumulativeRate,cumulativeErrorRate))
			if __name__ == "__main__":
				#with Pool(initializer=init_worker, initargs=(ref,mutMatrixGlobal,siteRates,errorRates,errorRateGlobal)) as pool:
				with Pool() as pool:
					results = pool.map(startTopologyUpdatesParallel, paralleleInputs)
			if aBayesPlusOn:
				improvementsFound=[]
				for i in range(numCores):
					improvementsFound.extend(results[i][0])
					for altPlac in results[i][1]:
						tree.support[altPlac[0]]=altPlac[1]
						if networkOutput:
							tree.alternativePlacements[altPlac[0]]=altPlac[2]
						#print("Number of SPRTA results found in parallel search")
						#print(len(results[i][1]))
						#for exa in range(10):
						#	print(results[i][1][exa])
			else:
				#print("aBayesPlusOn not on outside of parallelized bit")
				for i in range(numCores-1):
					results[0].extend(results[i+1])
				improvementsFound=results[0]

			improvementsFound.sort(reverse=False,key=itemgetter(2))
			totalTimeFindingParent[0]+=time()-start
			print("Found proposed SPR moves, merged, and sorted.")
			setAllDirty(tree,t1,dirtiness=False)
			newRoot, improvement = applySPRMovesParallel(tree,improvementsFound,strictTopologyStopRules=threshStricts[nRound],allowedFailsTopology=threshNums[nRound],thresholdLogLKtopology=threshValues[nRound],thresholdTopologyPlacement=threshPlaces[nRound])
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
		if doTimeTree:
			reCalculateAllGenomeListsTime(tree,t1,mutRate)
			print("Time LK: "+str(calculateTreeLikelihoodTime(tree,t1,mutRate)))
		
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
			if HnZ:
				print("Re-calculating the nDesc0")
				calculateNDesc0(tree,t1,checkExisting=True)
			if parallelize and (numDirty>0.1*numNodes):
				paralleleInputs=[]
				for i in range(numCores):
					paralleleInputs.append((tree,t1,i,threshStricts[nRound],threshNums[nRound],threshValues[nRound],threshPlaces[nRound],mutRate,   errorRateGlobal,mutMatrixGlobal,errorRates,mutMatrices,cumulativeRate,cumulativeErrorRate))
				if __name__ == "__main__":
					#with Pool(initializer=init_worker, initargs=(ref,mutMatrixGlobal,siteRates,errorRates,errorRateGlobal)) as pool:
					with Pool() as pool:
						results = pool.map(startTopologyUpdatesParallel, paralleleInputs)
				if aBayesPlusOn:
					improvementsFound=[]
					for i in range(numCores):
						improvementsFound.extend(results[i][0])
						for altPlac in results[i][1]:
							tree.support[altPlac[0]]=altPlac[1]
							if networkOutput:
								tree.alternativePlacements[altPlac[0]]=altPlac[2]
				else:
					for i in range(numCores-1):
						results[0].extend(results[i+1])
					improvementsFound=results[0]
				improvementsFound.sort(reverse=False,key=itemgetter(2))
				totalTimeFindingParent[0]+=time()-start
				print("Found proposed SPR moves, merged, and sorted.")
				setAllDirty(tree,t1,dirtiness=False)
				newRoot, improvement = applySPRMovesParallel(tree,improvementsFound,strictTopologyStopRules=threshStricts[nRound],allowedFailsTopology=threshNums[nRound],thresholdLogLKtopology=threshValues[nRound],thresholdTopologyPlacement=threshPlaces[nRound])
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
		if doTimeTree:
			reCalculateAllGenomeListsTime(tree,t1,mutRate)
			print("Time LK: "+str(calculateTreeLikelihoodTime(tree,t1,mutRate)))
		print("Time for the subrounds of this traversal of the tree: "+str(timeForUpdatingTopology), flush=True)
		timeTopology+=timeForUpdatingTopology

		#if estimating error rates, repeating EM after shallow topological search
		oldLK=float("-inf")
		newLk=calculateTreeLikelihood(tree,t1)
		print("Initial LK before EM: "+str(newLk), flush=True)
		mutMatrixGlobal, siteRates, errorRateGlobal, errorRates = expectationMaximizationCalculationRates(tree,t1)
		updateMutMatrices(mutMatrixGlobal,siteRates=siteRates)
		if estimateErrorRate or estimateSiteSpecificErrorRate:
			updateErrorRates(errorRateGlobal,errorRates=errorRates)
		reCalculateAllGenomeLists(tree,t1)
		newLk=calculateTreeLikelihood(tree,t1)
		print("LK after one round of EM: "+str(newLk))
		if estimateErrorRate or estimateSiteSpecificErrorRate:
			numEMsteps=0
			while (newLk-oldLK>1.0) and numEMsteps<20:
				if not doNotOptimiseBLengths:
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
		if not doNotOptimiseBLengths:
			start=time()
			reCalculateAllGenomeLists(tree,t1)
			newLk=calculateTreeLikelihood(tree,t1)
			setAllDirty(tree,t1)
			print(" branch length optimization starting from LK "+str(newLk), flush=True)
			improvement=traverseTreeToOptimizeBranchLengths(tree,t1)
			print("Branch length optimization round 1, number of changes: "+str(improvement), flush=True)
			subRound=0
			while subRound<20:
				if (not improvement):
					break
				subRound+=1
				improvement=traverseTreeToOptimizeBranchLengths(tree,t1)
			reCalculateAllGenomeLists(tree,t1)
			newLk=calculateTreeLikelihood(tree,t1)
			print("branch length finalization subround "+str(subRound+1)+" number of changes "+str(improvement)+" final LK: "+str(newLk))
			if doTimeTree:
				reCalculateAllGenomeListsTime(tree,t1,mutRate)
				print("Time LK: "+str(calculateTreeLikelihoodTime(tree,t1,mutRate)))
			timeForBranchOptimization=(time()-start)
			print("Time for updating branch lengths: "+str(timeForBranchOptimization), flush=True)
		if HnZ:
			print("Re-calculating the nDesc0")
			calculateNDesc0(tree,t1,checkExisting=True)

		# EM round to estimate the time-scaled mutation rate
		if doTimeTree:
			reCalculateAllGenomeListsTime(tree,t1,mutRate)
			oldLK=calculateTreeLikelihoodTime(tree,t1,mutRate)
			print("After SPR round "+str(nRound+1)+", pre-EM mutation rate "+str(mutRate)+" with time LK: "+str(oldLK))
			counts, waitingTimes, mutRate = expectationMaximizationCalculationRatesTime(tree,t1,mutRate)
			poissonCoeff=[[1.0]]
			reCalculateAllGenomeListsTime(tree,t1,mutRate)
			newLk=calculateTreeLikelihoodTime(tree,t1,mutRate)
			print("After SPR round "+str(nRound+1)+", estimated mutation rate "+str(mutRate)+" with time LK: "+str(newLk))
			numEMsteps=0
			while (newLk-oldLK>0.1) and numEMsteps<20:
				counts, waitingTimes, mutRate = expectationMaximizationCalculationRatesTime(tree,t1,mutRate)
				poissonCoeff=[[1.0]]
				reCalculateAllGenomeListsTime(tree,t1,mutRate)
				oldLK=newLk
				newLk=calculateTreeLikelihoodTime(tree,t1,mutRate)
				numEMsteps+=1
			print("New time LK step "+str(numEMsteps)+" mutRate "+str(mutRate)+": "+str(newLk))

		#writing to output the substitution model (possibly with rate variation and error rates)
		if nRound<(nRounds-1):
			fileNameAdd="_round"+str(nRound+1)
		else:
			fileNameAdd=""
		#print("\n nRound "+str(nRound)+" nRounds "+str(nRounds)+"\n")
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
		if doTimeTree:
			timeLK=calculateTreeLikelihoodTime(tree,t1,mutRate)
			print("Time LK: "+str(timeLK))
			totalLK+=timeLK
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

if rateVariation:
	useRateVariation=True
else:
	useRateVariation=False
if errorRateSiteSpecificFile or estimateSiteSpecificErrorRate:
	errorRateSiteSpecific=True
else:
	errorRateSiteSpecific=False
if errorRateSiteSpecificFile or errorRateFixed or estimateErrorRate or estimateSiteSpecificErrorRate:
	usingErrorRate=True
else:
	usingErrorRate=False
if aBayesPlus:
	aBayesPlusOn=True










