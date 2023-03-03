# MAPLE

MAPLE is a new approach for maximum likelihood phylogenetic estimation for genomic epidemiology, or otherwise for datasets with short divergence.
For a description of the methods used, see https://doi.org/10.1101/2022.03.22.485312


### Installation

The code doesn't need installation: just download the latest python script (MAPLEv...py) and use it to run MAPLE.
It is strongly recommended to execute the code using pypy3 to achieve best performance; running MAPLE without pypy3 might result in about 10 times longer runtime.
To install pypy3, see https://www.pypy.org/ .
If pypy3 cannot be installed, it is also possible to execute it with python3, but be aware that it will be substantially slower.


### Basic usage

The main script ("MAPLEv...py") takes in input a reference genome, and an alignment file in MAPLE format (see e.g. file "MAPLE_input_example.txt"), and estimates a maximum likelihood phylogeny. For example, if both the reference and the alignment are contained in the same input file, with the reference first, as in the example file MAPLE_input_example.txt contained in this repository, then the standard usage to perform pylogenetic tree estimation would be:

    pypy3 MAPLEv0.2.1.py --input inputMapleFile.txt --output MAPLE_outputFile

If the reference is not contained in the alignment MAPLE file, it can be specified with option --reference:

    pypy3 MAPLEv0.2.1.py --input inputMapleFile.txt --reference referenceFile.fa --output MAPLE_outputFile

The specified --output option is used by MAPLE to name the output files: the output tree will be named "MAPLE_outputFile_tree.tree" and the output containing the estimated model parameters will be "MAPLE_outputFile_subs.txt". You can us option --overwrite to overwrite existing files with those names.


### Creating an input MAPLE alignment file

The python script "createMapleFile.py" included in the repository can translate a fasta file into a MAPLE format file.


### Robinson-Foulds distance calculation

MAPLE can also be used to perform fast Robinson-Foulds distance calculation (using the algorithm from Day 1985) instead of performing tree inference.
This can be done by running:

    pypy3 MAPLEv0.2.1.py --inputTree inputTreeFile.tree --inputRFtrees otherInputTreeFiles.tree

The tree contained in the file specified with option --inputTree will be compared to all the trees in the file specified with option --inputRFtrees.
Having multiple trees in this second file is faster than performing one comparison at the time by running MAPLE on only 2 trees at the time.


### Online tree update (adding sequences to existing tree)

Given a tree previously estimated, and given an alignment containing the sequences of the samples in the tree, plus possibly some additional sequences (the same option can be used to run analyses on a given complete starting tree, if no input alignment samples are missing from the input tree), it's now possible to use MAPLE to add these additional samples to the given tree. To do this, run:

    pypy3 MAPLEv0.2.1.py --inputTree inputTreeFile.tree --input inputMapleFile.txt --output MAPLE_outputFile

By default, MAPLE will only update the topology of the parts of the tree affected by the addition of the new sequences - this will typically be much faster than running a new inference anew, unless many sequences are added to the tree. In the case one wants to not only add sequences to the tree, but also perform a full topological update, then option --largeUpdate can be used to force an extensive topological search over the whole tree.


### Substitution models

Thus far the substitution models JC69, GTR (default) and UNREST have been implemented. The model can be specified by the user with option --model.

As of version 0.2.1, MAPLE also includes a model of rate variation. This model assigns a free rate parameter to each genome positions; since this model is parameter-rich, its use is recommended only for larger datasets. To use this rate variation model, use option --rateVariation .


### Final tree likelihood

The total likelihood of the final tree can be calculated using option --calculateLKfinalTree, and will be printed to screen and written to file with name ending "_LK.txt".     


### Error probabilities estimation

As of version 0.3.1, MAPLE can model sequence errors by assigning an error probability free parameter to each genome position. Using option --estimateSiteSpecificErrorRate MAPLE will estimate these error probabilities and perform phylogenetic inference while accounting for these.



### Lineage assignment

MAPLE can now assign samples to given lineages. This is done by first estimating a tree (for example with MAPLE) relating all the input query sequences and the reference genomes (one reference genome to define each lineage). Then, this tree can be used as input for MAPLE, together with a file specifying which names correspond to the reference, which can be given in input to MAPLE using option --assignmentFile . MAPLE will then assign each sample to the lineage whose reference is its closest direct ancestor in the input tree.
To perform this analysis, run for example:

    pypy3 MAPLEv0.2.1.py --inputTree inputTreeFile.tree --assignmentFile pango_consensus_sequences.maple --output MAPLE_outputFile 

It is necessary that the names in the file used for --assignmentFile are contained in the input tree as sample names.



### Benchmarking




<br />
<br />
<br />
©EMBL-European Bioinformatics Institues, 2021
