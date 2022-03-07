# MAPLE

Code for fast likelihood calculation in the case of short-branched trees.

The main script is estimatePhylogenyIterativeFastLK.py which estimates a phylogenetic tree starting from a concise alignment ("MAPLE") file.
Typical use is:

pypy3 estimatePhylogenyIterativeFastLK.py --reference MN908947.3.fasta --input inputDiffFile.txt --model UNREST --numTopologyImprovements 3 --overwrite  --output fastLK_outputFile

The code doesn't need installation, just download the python scripts to use them.
However, it is highly recommended to execute the code using pypy3 to achieve best performance.
To install pypy3, see https://www.pypy.org/ .
If pypy3 cannot be installed, it is also possible to execute it with python3, but this will usually be slower.

The package includes an example MAPLE file (diffFile_repeat1_100samples_Ns.txt) of simulate genomes.

The package also includes other scripts, for example createDiffFile.py can translate a fasta file into a concise alignment ("MAPLE") file.

calculateFastLK.py contains code for calculating the likelihood (but not estimating a phylogeny) on short-branched trees.

extraBits.py contains code for comparing samples, for example to remove one sample if it is less informative than another one in the same dataset.

©EMBL-European Bioinformatics Institues, 2021
