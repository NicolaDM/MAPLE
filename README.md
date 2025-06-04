# MAPLE

MAPLE is a new approach for maximum likelihood phylogenetic estimation for genomic epidemiology, or otherwise for datasets with short divergence.
For a description of the methods used, see De Maio et al. 2023 Nature Genetics [https://doi.org/10.1101/2022.03.22.485312](https://www.nature.com/articles/s41588-023-01368-0).
For a description of new features including rate variation, a model of sequence errors, parallelization, and an algorithmic improvement, see De Maio et al 2024a [https://doi.org/10.1101/2024.07.12.603240](https://doi.org/10.1101/2024.07.12.603240). 
For a description of SPRTA, our new approach for assessing phylogenetic uncertainty, see De Maio et al. 2024b [https://doi.org/10.1101/2024.10.21.619398](https://doi.org/10.1101/2024.10.21.619398).
CMAPLE, the C++ implementation of MAPLE within IQ-TREE, can be found at [https://github.com/iqtree/cmaple](https://github.com/iqtree/cmaple) see Ly-Trong et al. 2024 MBE [https://doi.org/10.1093/molbev/msae134](https://doi.org/10.1093/molbev/msae134).

For a full(er) documentation, check the "MAPLE documentation.pdf" file in the repository.

For a tutorial, check the "MAPLE tutorial" file.


### Installation

The code doesn't need installation: just download the latest python script (MAPLEv...py) and use it to run MAPLE.
It is strongly recommended to execute the code using pypy3 to achieve best performance; running MAPLE without pypy3 might result in about 10 times longer runtime.
To install pypy3, see https://www.pypy.org/ . It is recommended that you use at least pypy v3.10 if you want to use parallelization in python.
If pypy3 cannot be installed, it is also possible to execute it with python3, but be aware that it will be substantially slower.


### Basic usage

The main script ("MAPLEv...py") takes in input an alignment file in MAPLE format (see e.g. file "MAPLE_alignment_example.txt"), containing a reference genome sequence, followed by all the considered genome sequences represented in terms of differences with respect to the reference. MAPLE's basic usage is:

    pypy3 MAPLEv0.6.8.py --input inputMapleFile.txt --output MAPLE_outputFilePrefix

The --output option is used by MAPLE to name the output files: in this case the final output tree will be named "MAPLE_outputFilePrefix_tree.tree", the file containing the estimated model parameters will be "MAPLE_outputFilePrefix_subs.txt", and so on. You can us option --overwrite to overwrite existing files with those names.

We provide in this repository a very small example input file "MAPLE_alignment_example.txt" which MAPLE should be able to analyse in just a few seconds to infer a phylogenetic tree and substitution rates (and possibly other outputs if other options are used).


### Creating an input MAPLE alignment file

The python script "createMapleFile.py", included in this repository, can be used to translate a fasta alignment file into a MAPLE format alignment file.

We recommend using MAPLE only with closely related genomes. When analysing non-closely related genomes (e.g. branch lengths >0.01) the software will be both slower and less accurate.
A multiple sequence alignment can be obtained by aligning every considered genome to the same reference, and removing inserted material, for example using MAFFT with options --auto --keeplength --addfragment.

We also recommend masking unreliable genome positions; MAPLE includes a substitution model that can be used to infer these positions, see below the section "Substitution models".
Finally, we noticed that it can be useful to mask deletions in the unput alignment, that is, replacing gap "-" characters with reference nucleotides. This is because errors (either alignment or consensus calling) at positions with common deletions can cause high ancestral sequence ancertainty and errors.


### Online tree inference (adding sequences to existing tree)

Given a tree previously estimated, and given an alignment containing the sequences of the samples in the tree, plus possibly some additional sequences, one can use MAPLE to add these additional samples to the given tree, and/or to improve the topology of the input tree. To do this, run:

    pypy3 MAPLEv0.6.8.py --inputTree inputTreeFile.tree --input inputMapleFile.txt --output MAPLE_outputFile

By default, MAPLE will only update the topology of the parts of the tree affected by the addition of the new sequences - this will typically be much faster than running a new inference anew, unless many sequences are added to the tree. In the case one wants to not only add sequences to the tree, but also perform a full topological update, then option --largeUpdate can be used to force an extensive topological search over the whole tree.


### Substitution models

So far we have implemented only nucleotide substitution models JC69, GTR, and UNREST. The model can be specified by the user with option --model.

We have now also developed in MAPLE a model of rate variation, see De Maio et al 2024a [https://doi.org/10.1101/2024.07.12.603240](https://doi.org/10.1101/2024.07.12.603240). This model assigns a free rate parameter to each genome positions; since this model is parameter-rich, its use is recommended only for larger datasets. To use this rate variation model, use option --rateVariation . We only recommend using this model with the UNREST substitution model (option --model UNREST).

Additionally, MAPLE includes a model of heterogeneous recurrent sequence errors, that can account for and estimate recurrent sequence errors (option --estimateSiteSpecificErrorRate), see De Maio et al 2024a [https://doi.org/10.1101/2024.07.12.603240](https://doi.org/10.1101/2024.07.12.603240). This model also assigns a free parameter (the error probability) to each genome positions; its use is recommended only for larger datasets and in conjunction with the rate variation model (option --rateVariation) and the UNREST substitution model (option --model UNREST).

In summary, to run the most advanced model in MAPLE, you can use options

    pypy3 MAPLEv0.6.8.py --input inputMapleFile.txt --output MAPLE_outputFile --model UNREST --rateVariation --estimateSiteSpecificErrorRate

Further, when using the sequence error model in MAPLE, it is possible to estimate individual sequence errors in the input alignment with option --estimateErrors . An output file will then contain estimated sequence errors, each with its posterior probability of being an error.

Note however that only part of the highly recurrent errors might be identified by MAPLE if the error rates are too high or correlated with one another. In these cases it might be better to first identify higly recurrent errors with MAPLE, then mask these columns from the alignment, then re-run MAPLE.
It is also important to start inference from a high-quality alignment, and in particular we recommend using Viridian genomes as they prevent calling many wrong reversions to the reference, see Hunt et al. 2024 [https://doi.org/10.1101/2024.04.29.591666](https://doi.org/10.1101/2024.04.29.591666).


### Parallelization

The most time-demanding part of MAPLE is the SPR search to improve the topology of the initial tree.
SPR search can now be run in parallel in MAPLE using multiple cores using option --numCores , see De Maio et al 2024a [https://doi.org/10.1101/2024.07.12.603240](https://doi.org/10.1101/2024.07.12.603240) . For example:

    pypy3 MAPLEv0.6.8.py --input inputMapleFile.txt --output MAPLE_outputFile --numCores 10

will parallelize the SPR search over 10 cores. It is not recommended to try to parallelize over an excessive number of cores (e.g. >20) since this can currently deteriorate the method's performance.
No matter the number of cores used, the initial stepwise addition will still be run sequentially on 1 core.


### Inferring mutation events and a mutation-annotated tree

MAPLE can infer mutation events on the final tree (that is, estimate a mutation-annotated tree) using option --estimateMAT .
Inferred mutations will be annotated with posterior probabilities, so that the same mutation might be annotated on multiple branches in case of non-negligeable uncertainty in its inference.


### Branch support

MAPLE can estimate branch support with a new pandemic-scale approach (SPRTA, see De Maio et al. 2024b [https://doi.org/10.1101/2024.10.21.619398](https://doi.org/10.1101/2024.10.21.619398)) using option --SPRTA.
MAPLE will then approximate the posterior probabilities of branches in the tree with positive length.
Note that these support scores are not intended as the posterior brobabilities of clades intended as sets of taxa, a typical for branch support methods like Felsenstein's bootstrap.
Instead, SPRTA scores are intended as a measure of confidence in the placements of ancestral genomes, or similarly in the inferred genome evolution history.
SPRTA scores are also assigned to terminal branches of the tree - these can be interpreted as placement probabilities for the corresponding genome sequences, similarly to pplacer [https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-11-538](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-11-538).

An advantage of the SPRTA approach in MAPLE is that its output can be interpreted as a phylogenetic network: for each tree branch we not only define a support score, but also a list of alternative placements (other branches or nodes in the tree the considered ancestral genome might have evolved from), each with its own estimated support probability. This network-like output can be obtained with option --networkOutput , for example:

    pypy3 MAPLEv0.6.11.py --input inputMapleFile.txt --output MAPLE_outputFile --SPRTA --networkOutput

MAPLE will create a metadata .tsv output file containing SPRTA support scores (in the "support" column) as well as, for each branch A, the list of branches that A is a plausible placement locations of ("supportTo" column). This format is particularly useful for visualizing alternative placements in Taxonium [https://taxonium.org/](https://taxonium.org/).
To do this, first create a json file to be used in Taxonium combining the MAPLE tree with the SPRTA metadata:

     newick_to_taxonium -i MAPLE_outputFile_tree.tree -m MAPLE_outputFile_metaData.tsv -o Taxonium_inputFile.jsonl -c support,supportGroup,supportTo

(other metadata columns like rootSupport and mutationsInf can also be included in the json file, as well as further metadata, if available, after being combined to the SPRTA metadata file).
Then open the resulting .jsonl file in Taxonium. For any branch, one can then highlight alternative placements of that branch by searching for the branch name in the "supportTo" search field in Taxonium.

Option --supportFor0Branches will make SPRTA also evaluate the support of 0-length branches, and therefore in particular the placement of all samples. Expect slightly longer runtime when using this option.


### Accounting for lineage abundance

As of version 0.6.12 MAPLE can account for lineage abundance during phylogenetic inference. The idea is that, for example when placing a genome on a phylogeny, it is more likely that the genome is a copy (or a descendant) of an abundant lineage, than of a rare one.
This is normally not accounted for in maximum likelihood phylogenetics. To achieve this, we have implemented two diferrent approaches collectively called HnZ (Horse not Zebra). 
The first one (HnZ1, corresponding to option --HnZ 1 ) accounts for the number of bifurcating tree topologies embedded within a multifurcating tree, and multiplies the multifurcating tree likelihood by this number.
The second approach (HnZ2, corresponding to option --HnZ 2 ) multiplies the phylogenetic tree likelihood by a sampling likelihood, representing the probability that genomes were sampled given their abundance (as inferred by the tree itself).
In our simulations, these HnZ methods (and in particular HnZ1) improve substantially the accuracy of MAPLE on simulated SARS-CoV-2 data, but they also require longer runtime.
A manuscript describing HnZ1 and HnZ2 is in preparation.



### Robinson-Foulds distance calculation

MAPLE can also be used to perform fast Robinson-Foulds distance calculation (using the algorithm from Day 1985) instead of performing tree inference.
This can be done by running:

    pypy3 MAPLEv0.3.6.py --inputTree inputTreeFile.tree --inputRFtrees otherInputTreeFiles.tree

The tree contained in the file specified with option --inputTree will be compared to all the trees in the file specified with option --inputRFtrees.
Having multiple trees in this second file is faster than running the script many times with 2 trees at the time, but you will need to specify option --multipleInputRFTrees to prevent MAPLE from only reading the first tree in the file.


<!--
%### Benchmarking

We also include in this repository a script that we used to run benchmarking analyses of MAPLE: MAPLE_benchmarking.py .
Comments at the top of this script will give instructions on how to use it within the context of benchmark analyses.
In short:
1) First obtain a tree, for example from http://hgdownload.soe.ucsc.edu/goldenPath/wuhCor1/UShER_SARS-CoV-2/ , and simulate genome evolution along the tree, for example using phastSim https://github.com/NicolaDM/phastSim .
2) If needed, add ambiguities to the simulated alignment and extract the tree where mutation-less branches are collapsed with

        pypy3 MAPLE_benchmarking.py --createTotalData
        
3) Create files containing subsamples of the global alignments using

        pypy3 MAPLE_benchmarking.py --createBashScript
        sh createSubsampleInputFiles.sh
        
(the latter command parallelizes file creating with bsub on a computational cluster).
4) For any subsample size (here 100000 as an example), first create the corresponding bash scripts with

        pypy3 MAPLE_benchmarking.py --createBashScript --numSamples 100000
        
and then submit execution of phylogenetic inference for all methods using

        sh submitUShER.sh ; sh submitFastTree.sh ; sh submitIQtree.sh ; sh submitRAxML.sh ; sh submitRAxML-NG.sh ; sh submitMaple.sh
        
when execution of UShER is finished, run matOptimize:

        sh submitmatOptimize.sh
        
and when this is finished, convert the output to newick:

        sh submitMatOptimizeConversion.sh
        
5) When the phylogenetic ineference methods are finished, run data collection on the results (measure execution time/memory, RF distances, etc):

        sh submitIQtreeLK_UShER.sh ; sh submitMapleLK_UShER.sh ; sh submitRF_UShER.sh ; sh submitParsimony_UShER.sh ; sh submitIQtreeLK_matOptimize.sh ; sh submitMapleLK_matOptimize.sh ; sh submitRF_matOptimize.sh ; sh submitParsimony_matOptimize.sh ; sh submitIQtreeLK_IQtree.sh ; sh submitMapleLK_IQtree.sh ; sh submitRF_IQtree.sh ; sh submitParsimony_IQtree.sh ; sh submitIQtreeLK_FastTree.sh ; sh submitMapleLK_FastTree.sh ; sh submitRF_FastTree.sh ; sh submitParsimony_FastTree.sh ; sh submitIQtreeLK_RAxML.sh ; sh submitMapleLK_RAxML.sh ; sh submitRF_RAxML.sh ; sh submitParsimony_RAxML.sh ; sh submitIQtreeLK_RAxML-NG.sh ; sh submitMapleLK_RAxML-NG.sh ; sh submitRF_RAxML-NG.sh ; sh submitParsimony_RAxML-NG.sh ; sh submitIQtreeLK_Maple.sh ; sh submitMapleLK_Maple.sh ; sh submitRF_Maple.sh ; sh submitParsimony_Maple.sh
        
6) Then to collect all results run:

        pypy3 MAPLE_benchmarking.py --collectResults
        
Then to prepare the input files for figure generation:

        pypy3 MAPLE_benchmarking.py --createFigures
        
And finally to generate the figures (will require matplotlib):

        python3 MAPLE_benchmarking.py --runFigureGeneration
    -->


<br />
<br />
<br />
©EMBL-European Bioinformatics Institute, 2023
