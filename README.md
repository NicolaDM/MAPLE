# MAPLE

The code in this repository implements the MAPLE approach for phylogenetic maximum likelihood estimation for genomic epidemiology datasets.
For a description of the methods used, see https://doi.org/10.1101/2022.03.22.485312

The main script ("MAPLEv...py") takes in input a reference genome, and an alignment file in MAPLE format (see file "MAPLE_exampleInput_100samples.txt" for an example), and estimates a maximum likelihood phylogeny. For example, typical usage might be:

pypy3 MAPLEv0.0.5.py --reference MN908947.3.fasta --input inputDiffFile.txt --output MAPLE_outputFile

So that the output will be preceded by the "MAPLE_outputFile" string.
The code doesn't need installation, just download the python scripts to use them.
However, it is highly recommended to execute the code using pypy3 to achieve best performance.
To install pypy3, see https://www.pypy.org/ .
If pypy3 cannot be installed, it is also possible to execute it with python3, but this will usually be slower.

The python script "createDiffFile.py" included in the repository can translate a fasta file into a MAPLE format file.

©EMBL-European Bioinformatics Institues, 2021
