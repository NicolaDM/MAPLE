# MAPLE

The code in this repository implements the MAPLE approach for phylogenetic maximum likelihood estimation for genomic epidemiology datasets.
For a description of the methods used, see https://doi.org/10.1101/2022.03.22.485312


### Installation

The code doesn't need installation, just download the python scripts (in particular the MAPLEv...py which runs MAPLE) to use them.
However, it is highly recommended to execute the code using pypy3 to achieve best performance; running MAPLE without pypy3 might result in about 10 times higher runtime.
To install pypy3, see https://www.pypy.org/ .
If pypy3 cannot be installed, it is also possible to execute it with python3, but be aware that it will be substantially slower.


### Basic usage

The main script ("MAPLEv...py") takes in input a reference genome, and an alignment file in MAPLE format (see e.g. file "MAPLE_input_example.txt"), and estimates a maximum likelihood phylogeny. For example, if both the reference and the alignment are contained in the same input file, typical usage might be:

    pypy3 MAPLEv0.1.5.py --input inputMapleFile.txt --output MAPLE_outputFile

This will perform a thorough phylogenetic inference.
If the reference is not contained in the alignment MAPLE file, then it can be specified with option --reference:

    pypy3 MAPLEv0.1.5.py --input inputMapleFile.txt --reference referenceFile.fa --output MAPLE_outputFile

Given the specified --output option, the output tree will be named "MAPLE_outputFile_tree.tree". Using option --overwrite the previous tree of the same name will be overwritten.


### Creating an input MAPLE alignment file

The python script "createMapleFile.py" included in the repository can translate a fasta file into a MAPLE format file.


### Faster runs

MAPLE has many parameters that can be changed to make the inference deeper and slower, or shallower and faster.
For simplicity, option --fast sets parameters so to have substantially faster inference, sacrificing only a small accuracy (under our simulations, MAPLE with option --fast is still more accurate on SARS-CoV-2 data than the other methods considered in our manuscript).


### Robinson-Foulds distance calculation

MAPLE can also be used to perform fast Robinson-Foulds distance calculation (using the algorithm from Day 1985) instead of performing tree inference.
This can be done by running:

    pypy3 MAPLEv0.1.5.py --inputTree inputTreeFile.tree --inputRFtrees otherInputTreeFiles.tree

The tree contained in the file specified with option --inputTree will be compared to all the trees in the file specified with option --inputRFtrees.
Having multiple trees in this second file is faster than performing one comparison at the time by running MAPLE on only 2 trees at the time.


### Online tree update (adding sequences to existing tree)

Given a tree previously estimated, and given an alignment containing the sequences of the samples in the tree, plus some additional sequences, it's now possible to use MAPLE to add these additional samples to the given tree. To do this, run:

    pypy3 MAPLEv0.1.5.py --inputTree inputTreeFile.tree --input inputMapleFile.txt --output MAPLE_outputFile

By default, MAPLE will only update the topology of the parts of the tree affected by the addition of the new sequences - this will typically be much faster than running a new inference anew, unless many sequences are added to the tree. In the case one wants to add sequencing to the tree, but also perform a full topological update, then option --largeUpdate can be used to force an extensive topological search over the whole tree.


### Substitution models

Thus far the substitution models JC69, GTR (default) and UNREST have been implemented. The model can be specified by the user with option --model.
THe substitution rates inferred by MAPLE will be written to file using the file name specified by option --output, for example when running

    pypy3 MAPLEv0.1.5.py --input inputMapleFile.txt --output MAPLE_outputFile
  
The substitution rates will be written in file MAPLE_outputFile_subs.txt


### Final tree likelihood

The total likelihood of the final tree can be calculated using option --calculateLKfinalTree, and will be printed to screen.     




<br />
<br />
<br />
©EMBL-European Bioinformatics Institues, 2021
