import sys
from math import log
import argparse
from time import time
import os.path
from operator import itemgetter

#©EMBL-European Bioinformatics Institute, 2024
#Developd by Nicola De Maio

#execution:
# pypy3 analyseMNMs.py
# pypy3 analyseMNMs.py --recombinationFile RIVET_final_recombinants_2024-06-23.txt
# python3 analyseMNMs.py --createFigures
# pypy3.10 analyseMNMs.py --cherries
# pypy3.10 analyseMNMs.py --createCherryAlignments

parser = argparse.ArgumentParser(description='Assess recurring multinucleotide substitution on SARS-COV-2 trees.')
parser.add_argument('--inputTree',default="Viridian_2M_noShortDel_deeper_SPRTA_tree.tree", help='Input newick tree.')
parser.add_argument('--inputTSV',default="Viridian_analysis/Viridian_2M_noShortDel_deeper_SPRTA_metaData.tsv", help='Input tsv metadata file.')
parser.add_argument('--inputAl',default="Viridian_2Malignment_filtered_masked_noShortDeletions.maple", help='Input MAPLE alignment.')
parser.add_argument("--createFigures", help="Treat all ambiguities as N (total missing information).", action="store_true")
parser.add_argument("--minNumDescendants",help="Minimum number of descendats for a mutation event to be considered.",  type=int, default=1)
parser.add_argument("--thresholdProb",help="Only consider substitution events with probability above this thrshold.",  type=float, default=0.95)
parser.add_argument("--cherries", help="Investigate cherries with mutations of interest inside, to see which ones happen around indels, low coverage, and heterozygosity.", action="store_true")
parser.add_argument("--recombinationFile", help="RIVET file to check how many of RIVET's recombinations might be due to multi-nuc recurrent mutations.", default="")
parser.add_argument("--createCherryAlignments", help="Create alignments showing the mutations of interest within the cherries.", action="store_true")
parser.add_argument("--stats", help="Run statistical tests.", action="store_true")
parser.add_argument("--hypergeom", help="Run the hypergeometric test (not allowed with pypy).", action="store_true")
parser.add_argument("--testLocation", help="Run the statistical test for significance of association with the consensus TRS motif and of enrichment in nonsynonymous mutations.", action="store_true")
args = parser.parse_args()

inputTree=args.inputTree
inputTSV=args.inputTSV
minNumDescendants=args.minNumDescendants
thresholdProb=args.thresholdProb
createFigures=args.createFigures
inputAl=args.inputAl
cherries=args.cherries
recombinationFile=args.recombinationFile
createCherryAlignments=args.createCherryAlignments
stats=args.stats
hypergeom=args.hypergeom
testLocation=args.testLocation

nucleotidesDict={}
nucleotidesDict["a"]=True
nucleotidesDict["c"]=True
nucleotidesDict["g"]=True
nucleotidesDict["t"]=True


if hypergeom:
	#import numpy as np
	from scipy.stats import hypergeom

	dictSimu={'C21302T-C21304A-G21305A': (0, 0.0, 0), 'C21304A-G21305A': (0, 1.0, 5), 'A28877T-G28878C': (0, 0.0, 2), 'G27382C-A27383T-T27384C': (0, 0.0, 0), 'T26491C-A26492T-T26497C': (0, 0.0, 0), 'G27758A-T27760A': (0, 0.0, 1), 'T27875C-C27881T-G27882C-C27883T': (0, 0.0, 0), 'C27881T-G27882C-C27883T': (0, 0.0, 0), 'C25162A-C25163A': (0, 0.0, 1), 'T21294A-G21295A-G21296A': (0, 0.0, 0), 'A27038T-T27039A-C27040A': (0, 0.0, 0), 'A507T-T508C-G509A': (0, 0.0, 0), 'T28881A-G28882A-G28883C': (0, 0.0, 0), 'A21550C-A21551T': (0, 0.0, 0), 'T21994C-T21995C': (0, 0.0, 1), 'C13423A-C13424A': (0, 0.0, 1), 'A4576T-T4579A': (0, 0.0, 1), 'A20284T-T20285C': (0, 0.0, 0), 'G11083T-C21575T': (52, 66.0, 88), 'G910A-T911A-C912A': (0, 0.0, 0), 'T27381C-G27382T-A27383G': (0, 0.0, 0), 'A5703T-G5704A-T5705A': (0, 0.0, 0), 'A3684T-G3685A-C3686A': (0, 0.0, 0), 'A27400T-A27403C-C27406G': (0, 0.0, 0), 'G28280C-A28281T-T28282A': (0, 0.0, 0), 'A4420G-C4421A-G4422A': (0, 0.0, 0), 'A28131G-C28132A-G28134A': (0, 0.0, 0), 'G29757A-T29758C-G29759C': (0, 0.0, 0), 'T21042G-A21043T-G21044C': (0, 0.0, 0), 'C25572T-A25573C-A25574T': (0, 0.0, 0), 'C21302T-T21304A-G21305A': (0, 0.0, 0), 'G27572A-C27573T-T27576C-C27577A': (0, 0.0, 0), 'A28855T-A28856C-G28857T': (0, 0.0, 0), 'C27481T-T27482C-T27484C': (0, 0.0, 0), 'C21302T-C21304A': (0, 0.0, 3), 'G27382C-A27383T': (0, 0.0, 1), 'C22716A-T22717C': (0, 0.0, 0), 'C11074T-G11083T': (4, 13.0, 22), 'T27672A-C27673A': (0, 0.0, 1), 'G6975T-G6977A': (0, 0.0, 1), 'G6513A-T6515A': (0, 0.0, 1), 'C19977A-C19979A': (0, 0.0, 0), 'T22207G-T22209C': (0, 0.0, 0), 'A27865T-T27866A': (0, 0.0, 0), 'A21892G-G21893A': (0, 0.0, 1), 'A10323G-C21575T': (8, 16.0, 26), 'T26485C-T26486A': (0, 0.0, 0), 'G11071C-C11074T': (0, 0.0, 1), 'C17734T-T17735C': (0, 0.0, 1), 'T27299C-A27300G': (0, 0.0, 1), 'G11083T-C29095T': (8, 16.0, 30), 'G11083T-C16887T': (21, 31.0, 52), 'G11083T-C29614T': (6, 17.0, 26), 'C11950T-C28472T': (0, 0.0, 1), 'G11083T-T27384C': (9, 19.0, 28), 'G11083T-C26681T': (9, 19.0, 28), 'G1820A-G11083T': (9, 19.0, 32), 'G23957T-T23959G': (0, 0.0, 0), 'A25562T-G25563A': (0, 0.0, 0), 'T27760A-T27761C': (0, 0.0, 0), 'T3370G-G3371A': (0, 0.0, 0), 'G11083T-A21137G': (13, 25.0, 37), 'G11083T-C11750T': (7, 15.0, 28), 'A28131G-C28132A': (0, 0.0, 0), 'C7303T-C9430T': (0, 0.0, 2), 'C9165T-G11083T': (2, 8.0, 16), 'A22194G-T22196G': (0, 0.0, 0), 'T27383A-C27384T': (0, 0.0, 1), 'G11083T-C21614T': (4, 9.0, 18), 'C9430T-C25521T': (0, 2.0, 7), 'C9430T-G11083T': (5, 12.0, 20), 'C3096T-C25521T': (0, 2.0, 6), 'C16887T-C21575T': (10, 19.0, 27), 'T21304A-G21305A': (0, 0.0, 1), 'C5869T-G11083T': (3, 9.0, 17), 'C25380T-A25381C': (0, 0.0, 1), 'C635T-C27577T': (0, 1.0, 4), 'A10323G-G11083T': (17, 26.0, 39), 'G11083T-C28657T': (3, 11.0, 21), 'G29688C-G29692C': (0, 0.0, 1), 'G21082C-A21083T': (0, 0.0, 0), 'A6024C-C6027T': (0, 0.0, 1), 'A26709G-T26767C': (0, 0.0, 0), 'G11083T-C23625T': (3, 8.0, 15), 'C635T-G11083T': (4, 11.5, 21), 'G11083T-C29750T': (2, 8.0, 17), 'T967C-G970A': (0, 0.0, 1), 'T16884G-T16885A': (0, 0.0, 0), 'T26491C-A26492T': (0, 0.0, 1), 'T6148A-G6149A': (0, 0.0, 0), 'C203T-G11083T': (5, 12.0, 25), 'C203T-C25521T': (0, 2.0, 5), 'A21137G-C21575T': (6, 15.0, 24), 'G11083T-C28603T': (5, 12.0, 26), 'C635T-C25521T': (0, 1.0, 4), 'G11083T-C25521T': (4, 13.0, 21), 'C3096T-G11083T': (7, 15.0, 26), 'G11083T-C23638T': (1, 7.0, 15), 'T28241A-T28243A': (0, 0.0, 0), 'C8140T-G11083T': (1, 6.0, 12), 'T27257A-G27258A': (0, 0.0, 0), 'T28877A-C28878G': (0, 0.0, 0), 'C683T-C25521T': (0, 1.0, 6), 'G11083T-C21855T': (1, 7.0, 17), 'G11083T-C14724T': (5, 11.0, 21), 'C9474T-G11083T': (1, 5.0, 12), 'C8964T-G11083T': (1, 8.0, 14), 'C1912T-C21575T': (2, 9.0, 17), 'G26466T-G26467C': (0, 0.0, 0), 'G11083T-C24023T': (2, 8.0, 19), 'G11083T-C25413T': (3, 9.0, 16), 'G11083T-C28087T': (4, 10.0, 18), 'C21575T-C25521T': (3, 8.0, 19), 'G11083T-C27476T': (3, 11.0, 21), 'C3096T-C9430T': (0, 1.0, 6), 'G11083T-C17550T': (1, 8.0, 17), 'C4331T-G11083T': (0, 5.0, 13), 'C19524T-C21575T': (1, 6.0, 13), 'C835T-C841A': (0, 0.0, 1), 'G11083T-C29733T': (3, 8.0, 15), 'C5144T-G11083T': (0, 8.0, 17), 'C7528T-G11083T': (1, 7.0, 15), 'C21575T-C24023T': (1, 5.0, 14), 'C1684T-G11083T': (3, 8.0, 16), 'C9438T-G11083T': (1, 6.0, 15), 'C9430T-C21614T': (0, 1.0, 6), 'C8140T-T8141C': (0, 0.0, 1), 'C1684T-C21575T': (1, 6.0, 16), 'C25797T-A25798C': (0, 0.0, 1), 'G11083T-G24933T': (2, 9.0, 16), 'G11083T-C29200T': (2, 7.0, 14), 'G11083T-C19983T': (5, 13.0, 24), 'G11083T-C20178T': (4, 9.0, 18), 'G11083T-C19524T': (3, 9.0, 18), 'T25704C-C25708T': (0, 0.0, 0), 'C21762T-C21846T': (0, 0.0, 2), 'C3096T-C21575T': (2, 9.0, 15), 'G11083T-C17676T': (2, 8.0, 16), 'G11083T-C29253T': (3, 9.0, 15), 'C6781T-T6782C': (0, 0.0, 1), 'C21575T-C26111T': (0, 3.0, 12), 'T12612C-G12614A': (0, 0.0, 0), 'C21575T-C27389T': (1, 9.5, 16), 'G11083T-C13665T': (2, 6.0, 17), 'G11083T-C18086T': (1, 6.0, 13), 'A11782G-C21575T': (0, 4.0, 11), 'C337T-G11083T': (3, 7.5, 14), 'G11083T-G25217T': (0, 3.0, 9), 'C1594T-G11083T': (2, 8.0, 16), 'C11750T-C21575T': (3, 10.0, 19), 'A10323G-C11750T': (0, 4.0, 13), 'T8141C-G8144A': (0, 0.0, 0), 'C27881T-G27882C': (0, 0.0, 1), 'G11083T-C14599T': (0, 7.0, 16), 'G28882A-G28883C': (0, 0.0, 0), 'C16887T-C25904T': (0, 2.0, 8), 'T27752C-G27758A': (0, 0.0, 1), 'T5077C-A5078T': (0, 0.0, 0), 'C21575T-C29095T': (4, 9.0, 17), 'C683T-G11083T': (4, 10.0, 19), 'C9491T-G11083T': (1, 7.0, 15), 'C4795T-G11083T': (1, 5.0, 10), 'C27879T-A27880C': (0, 0.0, 0), 'T25629C-G25630C': (0, 0.0, 0), 'C16887T-A21137G': (2, 7.0, 15), 'C635T-C9430T': (0, 1.0, 7), 'G11083T-C27741T': (2, 7.0, 15), 'C16887T-C27389T': (0, 4.0, 10), 'C21575T-C28093T': (3, 8.0, 14), 'C27739A-C27741T': (0, 0.0, 1), 'C203T-C635T': (0, 1.0, 4), 'C18744T-C21575T': (1, 5.0, 13), 'G24219A-G24220T': (0, 0.0, 0), 'G22335T-G22336T': (0, 0.0, 1), 'C9803T-G11083T': (1, 5.0, 18), 'G11083T-C12459T': (1, 5.5, 12), 'G11083T-C23277T': (0, 6.0, 12), 'C6285T-G11083T': (1, 4.0, 10), 'C21575T-G28079T': (0, 5.0, 11), 'G21295A-G21296A': (0, 0.0, 1), 'C21575T-T27384C': (4, 11.0, 20), 'C21575T-C26681T': (4, 12.0, 21), 'G11083T-C25916T': (2, 8.0, 15), 'T27039A-C27040A': (0, 0.0, 1), 'G4583A-T4585A': (0, 0.0, 0)}


# statistical tests of combinations of mutations
focusedMultiNucMuts={}
focusedMultiNucMuts["C21302T-C21304A-G21305A"]=0
focusedMultiNucMuts["C21304A-G21305A"]=0
focusedMultiNucMuts["A28877T-G28878C"]=0
focusedMultiNucMuts["G27382C-A27383T-T27384C"]=0
focusedMultiNucMuts["T26491C-A26492T-T26497C"]=0
focusedMultiNucMuts["G27758A-T27760A"]=0
focusedMultiNucMuts["T27875C-C27881T-G27882C-C27883T"]=0
focusedMultiNucMuts["C27881T-G27882C-C27883T"]=0
focusedMultiNucMuts["C25162A-C25163A"]=0
focusedMultiNucMuts["T21294A-G21295A-G21296A"]=0
focusedMultiNucMuts["A27038T-T27039A-C27040A"]=0
focusedMultiNucMuts["A507T-T508C-G509A"]=0#problematic
focusedMultiNucMuts["T28881A-G28882A-G28883C"]=0#problematic
focusedMultiNucMuts["A21550C-A21551T"]=0
focusedMultiNucMuts["T21994C-T21995C"]=0#problematic
focusedMultiNucMuts["C13423A-C13424A"]=0
focusedMultiNucMuts["A4576T-T4579A"]=0
focusedMultiNucMuts["A20284T-T20285C"]=0
focusedMultiNucMuts["G11083T-C21575T"]=0 #control
#rarer ones
focusedMultiNucMuts["G910A-T911A-C912A"]=0
focusedMultiNucMuts["T27381C-G27382T-A27383G"]=0
focusedMultiNucMuts["A5703T-G5704A-T5705A"]=0
focusedMultiNucMuts["A3684T-G3685A-C3686A"]=0
focusedMultiNucMuts["A27400T-A27403C-C27406G"]=0
focusedMultiNucMuts["G28280C-A28281T-T28282A"]=0
focusedMultiNucMuts["A4420G-C4421A-G4422A"]=0
focusedMultiNucMuts["A28131G-C28132A-G28134A"]=0
focusedMultiNucMuts["G29757A-T29758C-G29759C"]=0
focusedMultiNucMuts["T21042G-A21043T-G21044C"]=0
focusedMultiNucMuts["C25572T-A25573C-A25574T"]=0
focusedMultiNucMuts["C21302T-T21304A-G21305A"]=0
focusedMultiNucMuts["G27572A-C27573T-T27576C-C27577A"]=0
focusedMultiNucMuts["A28855T-A28856C-G28857T"]=0
focusedMultiNucMuts["C27481T-T27482C-T27484C"]=0
mutCombs=focusedMultiNucMuts.keys()
recordedCombs=list(mutCombs)

# G27382C-A27383T-T27384C	253	135	75	1724
# G27382C-A27383T	46	135	75

# C21302T-C21304A-G21305A	284	219	369	290
# C21304A-G21305A	452	369	290
# C21302T-C21304A	52	219	369

# T26491C-A26492T-T26497C	160	21	36	137
# T26491C-A26492T	14	21	36

# C22716A-T22717C	45	48	52
# T27672A-C27673A	37	47	179
# G6975T-G6977A	33	46	62
# G6513A-T6515A	33	46	34
# C19977A-C19979A	32	40	57
# T22207G-T22209C	31	53	43
# A27865T-T27866A	29	41	33
# A21892G-G21893A	28	44	72
# T26485C-T26486A	27	45	31
# G11071C-C11074T	26	79	1355
# C17734T-T17735C	25	94	26
# T27299C-A27300G	24	140	27
# C11950T-C28472T	21	172	232 #only one made of non-contiguous rare mutations
# G23957T-T23959G	20	22	23
# A25562T-G25563A	20	43	51
# T27760A-T27761C	20	60	61
# T3370G-G3371A	19	25	36
# A22194G-T22196G	18	88	24
# C25380T-A25381C	15	116	27
# G29688C-G29692C	15	81	48
# G21082C-A21083T	15	18	15
# A6024C-C6027T	15	20	570
# A26709G-T26767C	15	25	55
# T967C-G970A	14	70	83
# T16884G-T16885A	14	18	14
# T6148A-G6149A	14	34	30
# T28241A-T28243A	13	13	14

#run statistical tests
if stats:

	if testLocation:

		#statistically test enrichment in nonsynonymous substitutions among recurrent MNMs.
		file=open(inputAl)
		line=file.readline()
		line=file.readline()
		ref=""
		while line!="" and line!="\n" and line[0]!=">":
			ref+=line.replace("\n","")
			line=file.readline()
		ref=ref.upper()
		file.close()
		ORF1a=ref[265:13468]
		print(ORF1a[:3])
		print(len(ORF1a))
		print(float(len(ORF1a))/3)
		#print(ORF1a[-30:])
		#266..13468

		codonTranslation={
		"TTT":"F","TTC":"F","TTA":"L","TTG":"L","TCT":"S","TCA":"S","TCC":"S","TCG":"S","TAT":"Y","TAC":"Y","TAA":"*","TAG":"*","TGT":"C","TGC":"C","TGA":"*","TGG":"W",
		"CTT":"L","CTC":"L","CTA":"L","CTG":"L","CCT":"P","CCC":"P","CCA":"P","CCG":"P","CAT":"H","CAC":"H","CAA":"Q","CAG":"Q","CGT":"R","CGC":"R","CGA":"R","CGG":"R",
		"ATT":"I","ATC":"I","ATA":"I","ATG":"M","ACT":"T","ACC":"T","ACA":"T","ACG":"T","AAT":"N","AAC":"N","AAA":"K","AAG":"K","AGT":"S","AGC":"S","AGA":"R","AGG":"R",
		"GTT":"V","GTC":"V","GTA":"V","GTG":"V","GCT":"A","GCC":"A","GCA":"A","GCG":"A","GAT":"D","GAC":"D","GAA":"E","GAG":"E","GGT":"G","GGC":"G","GGA":"G","GGG":"G"}
		nucs=["A","C","G","T"]
		
		import random
		import math
		random.seed(1)
		gaps=[1,2,3]
		probs=[5.0/7.0,1.0/7.0,1.0/7.0]
		iGap=0
		averageNonsynRate=0.0
		averageNonsynRateNuc=0.0
		for gap in gaps:
			stops=0
			nonsyn=0
			syn=0
			synNuc=0
			nonsynNuc=0
			for i in range(10000000):
				i1=random.randint(3, 13197)
				#print(i1)
				codon1=math.floor(i1/3)
				#print(codon1)
				seq1=ORF1a[codon1*3:codon1*3+3]
				#print(seq1)
				# if codonTranslation[seq1]=="*":
				# 	print(codonTranslation[seq1])
				# 	print("ERROR")
				# 	exit()
				codPos1=i1%3
				oldNuc=seq1[codPos1]
				#print(oldNuc)
				newNuc=oldNuc
				while newNuc==oldNuc:
					iNuc=random.randint(0, 3)
					newNuc=nucs[iNuc]
				#print(newNuc)
				seq1new=list(seq1)
				seq1new[codPos1]=newNuc
				seq1new="".join(seq1new)
				nonsynNucFlag=False
				synNucFlag=False
				#print(seq1new)
				if codonTranslation[seq1new]=="*":
					stops+=1
				elif codonTranslation[seq1new]==codonTranslation[seq1]:
					syn+=1
					synNuc+=1
					synNucFlag=True
				else:
					nonsynNuc+=1
					nonsynNucFlag=True
				if synNucFlag or nonsynNucFlag:
					#now sample second mutation
					i2=i1+gap
					#print(i2)
					codon2=math.floor(i2/3)
					#print(codon2)
					if codon2==codon1:
						seq2=seq1new
					else:
						seq2=ORF1a[codon2*3:codon2*3+3]
					#print(seq2)
					if codonTranslation[seq2]=="*":
						print(codonTranslation[seq2])
						print("ERROR")
						exit()
					codPos2=i2%3
					oldNuc=seq2[codPos2]
					#print(oldNuc)
					newNuc=oldNuc
					while newNuc==oldNuc:
						iNuc=random.randint(0, 3)
						newNuc=nucs[iNuc]
					#print(newNuc)
					seq2new=list(seq2)
					seq2new[codPos2]=newNuc
					seq2new="".join(seq2new)
					#print(seq2new)
					if codonTranslation[seq2new]=="*":
						stops+=1
					elif codonTranslation[seq2new]!=codonTranslation[seq2]:
						if nonsynNucFlag:
							nonsyn+=1
						nonsynNuc+=1
					else:
						if nonsynNucFlag:
							syn+=1	
						synNuc+=1

				#print()
			print("Gap "+str(gap))
			print(stops)
			print(nonsyn)
			print(syn)
			print("Ratio nonsyn: "+str(nonsyn/float(syn+nonsyn)))
			print(nonsynNuc)
			print(synNuc)
			print("Ratio nonsyn individual substitutions: "+str(nonsynNuc/float(synNuc+nonsynNuc)))
			averageNonsynRate+=probs[iGap]*nonsyn/float(syn+nonsyn)
			averageNonsynRateNuc+=probs[iGap]*nonsynNuc/float(synNuc+nonsynNuc)
			iGap+=1
			print()
		print("Average nonsyn rate "+str(averageNonsynRate))
		from scipy.stats import binomtest
		print(binomtest(4, 7, p=averageNonsynRate, alternative='greater'))

		print("Average nonsyn rate individual substitutions "+str(averageNonsynRateNuc))
		print(binomtest(10, 14, p=averageNonsynRateNuc, alternative='greater'))
		exit()

#            FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
#   Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
#   Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
#   Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG

		consideredMutsLoc={}
		consideredMutsLoc["C21302T-C21304A-G21305A"]=0
		consideredMutsLoc["C21304A-G21305A"]=0
		consideredMutsLoc["A28877T-G28878C"]=0
		consideredMutsLoc["G27382C-A27383T-T27384C"]=0
		consideredMutsLoc["T26491C-A26492T-T26497C"]=0
		consideredMutsLoc["G27758A-T27760A"]=0
		consideredMutsLoc["T27875C-C27881T-G27882C-C27883T"]=0
		consideredMutsLoc["C27881T-G27882C-C27883T"]=0
		consideredMutsLoc["C25162A-C25163A"]=0
		consideredMutsLoc["T21294A-G21295A-G21296A"]=0
		consideredMutsLoc["A27038T-T27039A-C27040A"]=0
		consideredMutsLoc["A21550C-A21551T"]=0
		consideredMutsLoc["C13423A-C13424A"]=0
		consideredMutsLoc["A4576T-T4579A"]=0
		consideredMutsLoc["A20284T-T20285C"]=0
		numMNMs=15
		numFit=6
		GL=29903
		target=10*5
		p=float(target)/GL
		print(p)
		expected=15*p
		print(expected)
		
		print(binomtest(6, 15, p=p, alternative='greater'))

		# MNM1a is CGC-AAC, which is nonsyn requires 2 nonsyn substitutions
		#consideredMutsLoc["C21304A-G21305A"]=0

		# MNM1 contains CGC-AAC like MNM1a, plus CCA-CTA which is also nonsyn
		#consideredMutsLoc["C21302T-C21304A-G21305A"]=0

		# MNM2 involves 2 nonsyn substitutions together creating a syn substitution.
		#consideredMutsLoc["A28877T-G28878C"]=0

		# MNM3 is GAT->CTC
		# GAT-CAT-CTT-CTC (2 nonsyn 1 syn)
		# GAT-CAT-CAC-CTC (2 nonsyn 1 syn)
		# GAT-GTT-CTT-CTC (2 nonsyn 1 syn)
		# GAT-GTT-GTC-CTC (2 nonsyn 1 syn)
		# GAT-GAC-GTC-CTC (2 nonsyn 1 syn)
		# GAT-GAC-CAC-CTC (2 nonsyn 1 syn)
		#consideredMutsLoc["G27382C-A27383T-T27384C"]=0

		# MNM4 is non-coding
		#consideredMutsLoc["T26491C-A26492T-T26497C"]=0

		# MNM5 involves 2 nonsyn on separate codons.
		#consideredMutsLoc["G27758A-T27760A"]=0

		# MNM6 involves 1 syn and 1 nonsyn substitution
		#consideredMutsLoc["C25162A-C25163A"]=0

		# MNM7/7a disrupts ORF7b, and involves 3/4 substitutions and so we ignore them.
		#consideredMutsLoc["T27875C-C27881T-G27882C-C27883T"]=0
		#consideredMutsLoc["C27881T-G27882C-C27883T"]=0

		# MNM8 involves 1 syn and GGC-AAC (2 nonsyn)
		#consideredMutsLoc["T21294A-G21295A-G21296A"]=0

		# MNM9 involves 1 syn substitution and TCA->AAA (2 nonsyn)
		#consideredMutsLoc["A27038T-T27039A-C27040A"]=0

		# MNM10 involves the deletion of 1 aa and 2 out of 2 nonsyn substitutions.
		#consideredMutsLoc["A21550C-A21551T"]=0

		# MNM11 involves 1 syn and 1 nonsyn substitutions
		#consideredMutsLoc["C13423A-C13424A"]=0

		# MNM12 2 syn mutations
		#consideredMutsLoc["A4576T-T4579A"]=0

		# MNM13 2 nonsyn
		#consideredMutsLoc["A20284T-T20285C"]=0



		exit()


	file=open(inputTSV)
	line=file.readline()
	line=file.readline()
	numMutations={}
	numMutationPairs={}
	mutationsList=[]
	branchesList=[]
	while line!="\n" and line!="":
		linelist=line.split("\t")
		mutations=linelist[6]
		nodeName=linelist[0]
		if mutations!="":
			mutationlist=mutations.split(",")
			passedMutations=[]
			for mutation in mutationlist:
				mutationPair=mutation.split(":")
				support=float(mutationPair[1])
				if support>=thresholdProb:
					mutationName=mutationPair[0]
					passedMutations.append(mutationName)
					mutationsList.append(mutationName)
			branchesList.append(len(passedMutations))
			foundComb=False
			if len(passedMutations)>1:
				foundComb=False
				for mutComb in mutCombs:
					mutCombList=mutComb.split("-")
					found=True
					for mut in mutCombList:
						if not (mut in passedMutations):
							found=False
							break
					if found:
						foundComb=True
						focusedMultiNucMuts[mutComb]+=1
						break
				
				#count other mutation pairs that are not under focus
				for i in range(len(passedMutations)):
					if (not foundComb) or ( not (passedMutations[i] in mutCombList)):
						for j in range(len(passedMutations)):
							if i<j:
								if (not foundComb) or ( not (passedMutations[j] in mutCombList)):
									mutationPairName=passedMutations[i]+"-"+passedMutations[j]
									if mutationPairName in numMutationPairs:
										numMutationPairs[mutationPairName]+=1
									else:
										numMutationPairs[mutationPairName]=1

			#count individual single-nucleotide mutations
			for mutationName in passedMutations:
					if mutationName in numMutations:
						numMutations[mutationName]+=1
					else:
						numMutations[mutationName]=1
				
		line=file.readline()
	file.close()


	M=4144221

	for pair in mutCombs:
		occurrencies=focusedMultiNucMuts[pair]
		mutList=pair.split("-")
		stringToPrint=pair+"\t"+str(occurrencies)+"\t"
		min1=100000000
		min2=100000000
		mut1=False
		mut2=False
		for mut in mutList:
			stringToPrint+=mut+"\t"+str(numMutations[mut])+"\t"
			if numMutations[mut] <min1:
				min2=min1
				mut2=mut1
				min1=numMutations[mut]
				mut1=mut
			elif numMutations[mut] <min2:
				min2=numMutations[mut]
				mut2=mut
		n=min2
		N=min1
		k=occurrencies
		if hypergeom:
			rv = hypergeom(M, n, N)
			pValue = rv.sf(k-1)
		else:
			pValue=0.0
		print(stringToPrint+"\t"+str(pValue)+"\t"+str(pValue*29903*29903*3*3)+"\t"+str(dictSimu[pair]))

	print()

	sortedListMutationPairs=sorted(numMutationPairs.items(), key=lambda item: item[1], reverse=True)
	sorted100=sortedListMutationPairs[:1000]
	listSign=[]
	listNonSign=[]
	for pair in sorted100:
		occurrencies=int(pair[1])
		mut1=pair[0].split("-")[0]
		mut2=pair[0].split("-")[1]
		num1=numMutations[mut1]
		num2=numMutations[mut2]
		n=num1
		N=num2
		k=occurrencies
		if hypergeom:
			rv = hypergeom(M, n, N)
			pValue = rv.sf(k-1)
		else:
			pValue=0.0
		if occurrencies>9:
			#if pValue<10**(-25):
			stringToPrint=pair[0]+"\t"+str(pair[1])+"\t"+str(num1)+"\t"+str(num2)+"\t"+str(pValue)+"\t"+str(pValue*29903*29903*3*3)+"\t"+str(dictSimu[pair[0]])
			print(stringToPrint)
			#else:
			#	listNonSign.append(stringToPrint)
			recordedCombs.append(pair[0])

	print()

	if not hypergeom:
		

		#here results from simulations have been hardcoded below, but simulations can stil be run with the code further below.
		alreadySimulated=True

		if alreadySimulated:
			simuDict={'C21302T-C21304A-G21305A': [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
			'C21304A-G21305A': [2, 1, 1, 3, 1, 2, 2, 2, 4, 1, 0, 3, 1, 0, 2, 0, 1, 3, 3, 1, 0, 2, 0, 0, 1, 1, 2, 0, 1, 2, 1, 3, 1, 0, 2, 0, 0, 1, 2, 2, 1, 0, 2, 0, 0, 0, 1, 1, 1, 0, 0, 3, 0, 0, 0, 0, 0, 2, 1, 0, 0, 1, 0, 1, 2, 1, 0, 2, 0, 2, 1, 3, 1, 5, 1, 1, 0, 0, 0, 1, 0, 1, 0, 1, 1, 1, 2, 2, 2, 0, 1, 0, 0, 1, 0, 1, 0, 0, 2, 0], 
			'A28877T-G28878C': [0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 0, 0, 2, 0, 1, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
			'G27382C-A27383T-T27384C': [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
			'T26491C-A26492T-T26497C': [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
			'G27758A-T27760A': [0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0], 
			'T27875C-C27881T-G27882C-C27883T': [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
			'C27881T-G27882C-C27883T': [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
			'C25162A-C25163A': [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
			'T21294A-G21295A-G21296A': [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
			'A27038T-T27039A-C27040A': [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
			'A507T-T508C-G509A': [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
			'T28881A-G28882A-G28883C': [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
			'A21550C-A21551T': [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
			'T21994C-T21995C': [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0], 
			'C13423A-C13424A': [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
			'A4576T-T4579A': [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
			'A20284T-T20285C': [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
			'G11083T-C21575T': [64, 67, 60, 73, 59, 65, 77, 85, 63, 68, 57, 82, 75, 57, 54, 71, 56, 58, 73, 65, 68, 56, 59, 58, 70, 62, 76, 66, 68, 69, 70, 63, 71, 60, 69, 65, 59, 60, 62, 77, 78, 77, 65, 69, 72, 78, 72, 70, 72, 70, 81, 63, 65, 67, 69, 69, 67, 81, 79, 64, 66, 59, 55, 62, 86, 68, 59, 76, 58, 64, 69, 68, 60, 72, 66, 58, 64, 60, 72, 63, 53, 71, 67, 64, 78, 62, 67, 64, 64, 58, 64, 66, 64, 73, 71, 59, 52, 88, 55, 71], 
			'G910A-T911A-C912A': [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
			'T27381C-G27382T-A27383G': [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
			'A5703T-G5704A-T5705A': [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
			'A3684T-G3685A-C3686A': [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
			'A27400T-A27403C-C27406G': [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
			'G28280C-A28281T-T28282A': [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
			'A4420G-C4421A-G4422A': [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
			'A28131G-C28132A-G28134A': [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
			'G29757A-T29758C-G29759C': [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
			'T21042G-A21043T-G21044C': [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
			'C25572T-A25573C-A25574T': [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
			'C21302T-T21304A-G21305A': [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
			'G27572A-C27573T-T27576C-C27577A': [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
			'A28855T-A28856C-G28857T': [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
			'C27481T-T27482C-T27484C': [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
			'C21302T-C21304A': [0, 1, 0, 0, 0, 1, 1, 0, 2, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 1, 1, 0, 0, 0, 1, 0, 2, 1, 1, 0, 1, 0, 0, 1, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 3, 0, 0, 0, 0, 1, 2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 1, 0, 0, 0], 
			'G27382C-A27383T': [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1], 
			'C22716A-T22717C': [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
			'C11074T-G11083T': [13, 7, 10, 17, 15, 11, 16, 12, 12, 17, 10, 12, 14, 10, 13, 11, 19, 14, 16, 10, 14, 12, 7, 14, 16, 8, 11, 15, 18, 17, 15, 8, 10, 11, 4, 10, 14, 10, 16, 4, 10, 18, 15, 17, 16, 11, 12, 13, 6, 13, 20, 6, 22, 15, 16, 13, 18, 22, 12, 8, 14, 14, 12, 19, 9, 11, 9, 18, 13, 10, 10, 11, 12, 15, 10, 12, 10, 11, 12, 17, 11, 18, 15, 16, 16, 11, 21, 6, 14, 13, 10, 19, 10, 15, 21, 17, 10, 20, 7, 13], 
			'T27672A-C27673A': [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
			'G6975T-G6977A': [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
			'G6513A-T6515A': [0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
			'C19977A-C19979A': [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
			'T22207G-T22209C': [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
			'A27865T-T27866A': [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
			'A21892G-G21893A': [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
			'A10323G-C21575T': [18, 14, 14, 11, 16, 14, 19, 12, 19, 16, 20, 16, 16, 15, 12, 13, 18, 25, 14, 11, 18, 18, 9, 15, 15, 19, 14, 19, 14, 10, 17, 10, 12, 24, 18, 19, 11, 13, 9, 25, 17, 10, 19, 10, 11, 17, 13, 10, 17, 11, 16, 13, 12, 11, 16, 19, 16, 20, 24, 14, 20, 16, 16, 23, 12, 17, 22, 19, 21, 10, 20, 19, 18, 25, 22, 21, 23, 15, 12, 19, 19, 21, 11, 21, 19, 14, 12, 16, 19, 18, 19, 19, 16, 8, 15, 21, 16, 8, 26, 10], 
			'T26485C-T26486A': [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
			'G11071C-C11074T': [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0], 
			'C17734T-T17735C': [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
			'T27299C-A27300G': [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
			'G11083T-C29095T': [13, 12, 11, 11, 17, 18, 14, 19, 18, 11, 8, 16, 16, 14, 12, 15, 15, 12, 13, 17, 16, 11, 18, 20, 20, 9, 14, 16, 19, 16, 16, 11, 20, 14, 12, 17, 9, 30, 9, 24, 10, 21, 12, 17, 12, 14, 15, 24, 16, 10, 16, 20, 14, 18, 20, 14, 10, 10, 14, 23, 23, 14, 15, 14, 17, 18, 10, 14, 11, 22, 15, 20, 20, 30, 19, 17, 18, 25, 16, 15, 13, 23, 13, 25, 17, 17, 15, 18, 17, 14, 14, 9, 10, 16, 11, 15, 22, 19, 18, 25], 
			'G11083T-C16887T': [28, 25, 34, 24, 24, 34, 40, 25, 31, 33, 31, 35, 29, 34, 29, 33, 33, 33, 29, 36, 26, 48, 25, 33, 23, 29, 35, 26, 33, 34, 29, 26, 33, 34, 27, 34, 32, 30, 35, 34, 24, 30, 26, 23, 22, 30, 33, 23, 21, 27, 41, 31, 27, 31, 25, 33, 31, 25, 33, 34, 32, 35, 32, 24, 35, 26, 38, 24, 31, 36, 36, 32, 33, 52, 26, 34, 32, 31, 33, 32, 35, 30, 30, 23, 32, 22, 24, 40, 34, 24, 29, 28, 25, 38, 40, 40, 32, 28, 28, 27], 
			'G11083T-C29614T': [14, 17, 12, 20, 11, 22, 18, 15, 22, 14, 12, 12, 14, 15, 15, 17, 24, 24, 18, 18, 6, 22, 17, 25, 20, 15, 19, 18, 15, 20, 22, 21, 18, 24, 16, 13, 16, 22, 14, 17, 16, 18, 11, 16, 17, 20, 18, 18, 14, 16, 10, 18, 21, 16, 22, 14, 18, 8, 9, 16, 11, 11, 12, 17, 12, 17, 23, 11, 10, 19, 20, 18, 9, 23, 18, 15, 14, 21, 25, 18, 18, 20, 16, 26, 22, 16, 15, 7, 15, 23, 19, 19, 23, 13, 24, 18, 15, 18, 17, 15], 
			'C11950T-C28472T': [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0], 
			'G11083T-T27384C': [20, 17, 24, 18, 17, 22, 25, 22, 24, 13, 9, 22, 19, 22, 24, 22, 15, 20, 24, 25, 18, 21, 15, 19, 22, 14, 17, 20, 20, 17, 15, 13, 26, 20, 14, 18, 22, 16, 20, 20, 16, 22, 25, 15, 20, 21, 16, 20, 15, 15, 18, 14, 19, 16, 10, 16, 13, 19, 25, 13, 19, 28, 19, 23, 20, 26, 19, 20, 25, 12, 22, 21, 15, 19, 22, 21, 17, 18, 16, 16, 19, 14, 24, 14, 21, 23, 23, 15, 24, 12, 21, 19, 15, 22, 23, 15, 16, 19, 16, 14], 
			'G11083T-C26681T': [15, 20, 14, 15, 22, 15, 14, 20, 17, 24, 23, 16, 12, 20, 17, 17, 22, 20, 19, 21, 27, 20, 16, 12, 23, 22, 26, 16, 14, 9, 24, 13, 15, 18, 21, 15, 23, 16, 24, 19, 18, 23, 14, 19, 21, 9, 15, 20, 15, 13, 12, 22, 18, 23, 28, 12, 25, 19, 15, 28, 20, 26, 19, 21, 16, 16, 17, 17, 17, 26, 22, 24, 18, 22, 20, 25, 10, 25, 15, 16, 22, 19, 20, 22, 19, 14, 22, 15, 13, 28, 19, 20, 15, 20, 16, 19, 15, 15, 20, 20], 
			'G1820A-G11083T': [18, 31, 20, 15, 15, 15, 17, 19, 19, 17, 12, 21, 17, 24, 19, 22, 23, 18, 24, 16, 15, 17, 21, 23, 18, 20, 24, 13, 16, 28, 13, 22, 16, 18, 22, 13, 20, 21, 16, 24, 29, 18, 17, 22, 26, 16, 19, 26, 16, 18, 22, 23, 24, 14, 22, 20, 24, 20, 15, 17, 16, 21, 18, 17, 21, 15, 18, 19, 19, 18, 15, 17, 17, 19, 26, 21, 22, 22, 21, 11, 14, 13, 15, 19, 15, 23, 22, 32, 30, 22, 20, 18, 21, 27, 26, 14, 9, 23, 16, 13], 
			'G23957T-T23959G': [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
			'A25562T-G25563A': [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
			'T27760A-T27761C': [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
			'T3370G-G3371A': [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
			'G11083T-A21137G': [21, 27, 25, 24, 37, 25, 25, 20, 26, 22, 20, 30, 23, 23, 31, 32, 20, 20, 25, 24, 22, 20, 34, 31, 30, 25, 32, 27, 27, 23, 25, 22, 26, 26, 20, 32, 25, 18, 30, 22, 31, 29, 20, 19, 15, 22, 27, 26, 20, 23, 28, 22, 24, 30, 30, 22, 30, 29, 26, 29, 25, 24, 14, 24, 31, 28, 29, 22, 23, 23, 24, 16, 23, 24, 24, 23, 26, 13, 23, 35, 22, 26, 30, 22, 28, 17, 21, 26, 27, 30, 27, 21, 28, 26, 28, 27, 19, 23, 23, 16], 
			'G11083T-C11750T': [15, 12, 17, 7, 14, 14, 16, 15, 14, 11, 15, 12, 18, 15, 13, 14, 13, 14, 12, 22, 21, 17, 14, 19, 18, 10, 18, 16, 12, 21, 20, 9, 22, 13, 15, 13, 8, 22, 14, 12, 10, 20, 17, 14, 13, 16, 13, 12, 13, 22, 17, 16, 19, 16, 14, 21, 16, 17, 15, 11, 19, 15, 18, 20, 21, 13, 12, 20, 15, 22, 28, 18, 17, 15, 15, 14, 12, 16, 14, 9, 16, 13, 15, 15, 17, 14, 15, 12, 17, 14, 13, 15, 21, 13, 13, 20, 16, 19, 27, 15], 
			'A28131G-C28132A': [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
			'C7303T-C9430T': [1, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 1, 0, 0, 1, 2, 0, 1, 1, 0, 0, 0, 0, 0, 2, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 2, 2, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 2, 1, 0, 1, 1, 1, 0, 0, 0], 
			'C9165T-G11083T': [6, 9, 9, 11, 9, 11, 12, 7, 6, 7, 7, 6, 7, 13, 11, 7, 6, 10, 2, 10, 9, 10, 9, 8, 4, 7, 5, 10, 6, 10, 8, 7, 11, 8, 10, 12, 8, 6, 10, 9, 7, 8, 10, 5, 13, 7, 16, 5, 11, 10, 8, 6, 2, 10, 7, 8, 13, 11, 10, 9, 10, 8, 9, 4, 8, 4, 9, 10, 13, 6, 9, 7, 4, 10, 10, 6, 10, 7, 4, 8, 6, 10, 9, 8, 8, 3, 9, 7, 4, 8, 4, 5, 8, 4, 7, 5, 10, 13, 8, 14], 
			'A22194G-T22196G': [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
			'T27383A-C27384T': [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
			'G11083T-C21614T': [5, 12, 6, 14, 6, 10, 12, 4, 7, 9, 17, 8, 11, 10, 9, 9, 9, 11, 7, 15, 10, 5, 7, 11, 14, 10, 9, 9, 7, 10, 7, 7, 8, 11, 13, 13, 11, 9, 9, 8, 9, 12, 10, 11, 10, 9, 6, 18, 9, 9, 10, 10, 7, 7, 6, 9, 4, 5, 8, 6, 13, 8, 7, 11, 10, 9, 10, 6, 12, 11, 7, 5, 6, 15, 8, 15, 8, 11, 9, 6, 16, 7, 10, 12, 7, 5, 13, 12, 6, 7, 6, 6, 8, 9, 8, 11, 13, 11, 7, 7], 
			'C9430T-C25521T': [3, 3, 3, 0, 4, 2, 3, 6, 2, 1, 2, 4, 1, 1, 4, 1, 3, 0, 2, 3, 2, 0, 3, 1, 2, 1, 2, 0, 4, 2, 1, 1, 3, 5, 0, 3, 2, 2, 1, 7, 1, 1, 2, 0, 1, 0, 0, 1, 1, 3, 1, 2, 1, 3, 1, 1, 3, 3, 3, 2, 0, 2, 1, 0, 3, 0, 1, 1, 3, 1, 1, 2, 2, 1, 3, 4, 4, 2, 1, 3, 2, 4, 0, 1, 2, 2, 2, 1, 1, 0, 3, 0, 2, 2, 0, 6, 5, 1, 1, 0], 
			'C9430T-G11083T': [13, 18, 9, 11, 11, 10, 15, 14, 14, 7, 19, 14, 13, 11, 15, 13, 13, 11, 12, 14, 12, 11, 8, 6, 17, 10, 12, 12, 15, 8, 12, 8, 19, 17, 15, 15, 16, 11, 18, 6, 16, 8, 15, 12, 9, 10, 18, 7, 17, 13, 15, 9, 12, 12, 11, 8, 19, 11, 10, 19, 9, 11, 13, 17, 10, 13, 13, 20, 10, 10, 12, 13, 8, 5, 17, 14, 8, 12, 18, 11, 8, 12, 11, 9, 8, 10, 11, 9, 11, 12, 20, 15, 9, 11, 8, 7, 12, 15, 16, 15], 
			'C3096T-C25521T': [3, 1, 3, 3, 4, 3, 4, 5, 2, 4, 2, 1, 4, 2, 1, 2, 3, 3, 2, 2, 2, 2, 3, 1, 1, 2, 4, 0, 4, 1, 0, 2, 4, 0, 2, 4, 2, 5, 0, 5, 3, 5, 0, 2, 3, 4, 1, 3, 3, 2, 2, 2, 1, 1, 3, 1, 4, 2, 2, 5, 0, 2, 2, 2, 0, 3, 2, 3, 3, 2, 2, 4, 2, 2, 2, 0, 4, 3, 2, 3, 5, 2, 4, 6, 0, 2, 3, 1, 4, 1, 0, 1, 2, 1, 2, 3, 1, 3, 2, 3], 
			'C16887T-C21575T': [13, 21, 17, 16, 22, 19, 20, 18, 20, 27, 15, 11, 17, 13, 23, 21, 24, 20, 13, 18, 13, 19, 20, 18, 22, 18, 19, 20, 17, 20, 16, 21, 16, 22, 18, 13, 12, 22, 21, 25, 25, 19, 27, 25, 22, 22, 21, 19, 22, 19, 23, 15, 23, 25, 24, 16, 15, 21, 27, 21, 20, 17, 16, 19, 27, 17, 15, 27, 15, 13, 18, 22, 16, 17, 16, 19, 19, 16, 25, 13, 13, 25, 10, 17, 15, 21, 17, 21, 21, 22, 13, 14, 23, 18, 25, 19, 24, 21, 16, 17], 
			'T21304A-G21305A': [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
			'C5869T-G11083T': [8, 8, 6, 11, 16, 17, 4, 6, 8, 13, 9, 9, 8, 8, 9, 9, 11, 8, 14, 10, 8, 11, 12, 13, 5, 11, 11, 6, 12, 8, 7, 8, 4, 8, 11, 4, 8, 7, 8, 6, 9, 9, 9, 11, 6, 6, 12, 11, 5, 13, 3, 7, 6, 5, 6, 13, 7, 10, 11, 6, 5, 11, 10, 8, 11, 10, 3, 10, 10, 9, 11, 8, 5, 7, 10, 13, 8, 10, 6, 8, 6, 9, 7, 7, 7, 9, 10, 9, 4, 9, 10, 11, 9, 12, 9, 7, 11, 11, 10, 11], 
			'C25380T-A25381C': [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
			'C635T-C27577T': [0, 1, 1, 1, 0, 3, 1, 3, 2, 0, 4, 3, 0, 1, 1, 2, 1, 0, 0, 0, 0, 2, 1, 0, 1, 2, 3, 3, 2, 1, 2, 0, 1, 0, 1, 1, 1, 3, 0, 4, 0, 1, 0, 4, 0, 2, 3, 0, 2, 1, 0, 0, 2, 0, 1, 1, 0, 0, 0, 2, 1, 1, 1, 2, 0, 1, 2, 1, 0, 2, 2, 2, 2, 1, 0, 1, 0, 1, 2, 2, 0, 2, 2, 1, 0, 1, 0, 2, 1, 1, 2, 1, 2, 3, 1, 2, 0, 2, 2, 1], 
			'A10323G-G11083T': [25, 28, 29, 29, 29, 28, 35, 19, 30, 25, 26, 23, 38, 34, 26, 34, 32, 38, 30, 31, 24, 24, 18, 37, 32, 31, 27, 31, 27, 29, 39, 22, 27, 27, 29, 20, 19, 33, 29, 31, 23, 23, 18, 24, 26, 19, 27, 30, 25, 24, 27, 32, 32, 22, 39, 22, 34, 27, 26, 24, 25, 22, 22, 29, 19, 20, 22, 19, 28, 19, 22, 20, 23, 26, 27, 25, 17, 21, 28, 20, 22, 20, 21, 24, 17, 17, 30, 27, 25, 26, 25, 34, 26, 30, 30, 24, 34, 22, 27, 25], 
			'G11083T-C28657T': [14, 10, 21, 12, 11, 14, 6, 10, 15, 10, 11, 9, 13, 8, 15, 11, 10, 10, 13, 7, 11, 12, 17, 20, 16, 13, 10, 11, 8, 11, 10, 8, 15, 14, 9, 9, 14, 12, 14, 14, 14, 13, 8, 15, 9, 18, 10, 11, 11, 13, 12, 7, 10, 6, 14, 9, 12, 13, 12, 13, 11, 10, 12, 17, 9, 9, 17, 18, 13, 11, 10, 8, 9, 12, 21, 11, 11, 12, 13, 10, 10, 11, 12, 5, 19, 8, 10, 5, 12, 8, 13, 3, 13, 10, 11, 15, 10, 20, 14, 10], 
			'G29688C-G29692C': [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
			'G21082C-A21083T': [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
			'A6024C-C6027T': [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
			'A26709G-T26767C': [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
			'G11083T-C23625T': [3, 7, 3, 9, 12, 7, 10, 9, 9, 6, 5, 14, 11, 6, 9, 12, 6, 5, 7, 9, 11, 5, 9, 4, 6, 8, 8, 12, 11, 14, 9, 5, 6, 9, 9, 9, 8, 6, 7, 5, 6, 6, 11, 7, 3, 7, 9, 6, 7, 13, 12, 7, 10, 11, 8, 8, 5, 10, 7, 7, 5, 9, 4, 4, 3, 9, 7, 7, 9, 10, 15, 6, 9, 8, 9, 6, 3, 11, 8, 11, 5, 5, 11, 3, 10, 7, 7, 10, 5, 5, 8, 10, 9, 11, 6, 6, 11, 9, 8, 6], 
			'C635T-G11083T': [11, 10, 7, 12, 11, 12, 8, 12, 10, 19, 11, 14, 15, 17, 13, 10, 21, 8, 11, 10, 14, 8, 8, 7, 13, 7, 19, 13, 12, 15, 6, 10, 8, 10, 13, 12, 6, 17, 5, 18, 17, 7, 8, 17, 16, 5, 9, 10, 16, 10, 8, 7, 17, 6, 8, 9, 13, 12, 13, 14, 11, 12, 11, 14, 7, 13, 12, 9, 8, 18, 4, 8, 9, 15, 9, 16, 11, 9, 14, 12, 12, 15, 13, 11, 13, 9, 13, 8, 15, 14, 10, 13, 17, 11, 10, 11, 15, 14, 13, 18], 
			'G11083T-C29750T': [6, 17, 9, 7, 4, 15, 8, 7, 6, 9, 6, 4, 7, 7, 9, 7, 9, 6, 7, 8, 7, 6, 11, 14, 5, 12, 6, 8, 7, 8, 9, 9, 5, 8, 12, 10, 12, 7, 5, 17, 11, 7, 6, 5, 2, 2, 10, 11, 10, 11, 11, 5, 9, 5, 8, 8, 6, 6, 2, 7, 7, 3, 8, 6, 5, 4, 10, 14, 12, 7, 5, 9, 9, 5, 8, 7, 8, 13, 12, 7, 9, 12, 8, 8, 6, 11, 9, 13, 5, 12, 7, 5, 6, 6, 5, 12, 10, 9, 13, 13], 
			'T967C-G970A': [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
			'T16884G-T16885A': [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
			'T26491C-A26492T': [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
			'T6148A-G6149A': [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
			'C203T-G11083T': [16, 11, 11, 18, 11, 8, 9, 21, 15, 8, 12, 19, 13, 14, 7, 10, 12, 20, 17, 12, 13, 15, 15, 14, 12, 11, 12, 12, 6, 9, 9, 12, 7, 10, 11, 11, 12, 12, 10, 14, 17, 14, 12, 9, 17, 12, 19, 9, 21, 13, 9, 16, 14, 12, 11, 22, 12, 12, 19, 8, 8, 17, 19, 15, 15, 15, 8, 19, 18, 18, 21, 7, 14, 11, 5, 25, 11, 10, 17, 11, 11, 11, 13, 13, 10, 13, 14, 13, 12, 11, 8, 9, 11, 13, 7, 9, 15, 11, 12, 10], 
			'C203T-C25521T': [2, 3, 1, 1, 3, 0, 3, 4, 3, 0, 3, 4, 2, 2, 2, 1, 1, 0, 1, 0, 3, 2, 3, 0, 2, 4, 1, 4, 1, 3, 4, 4, 3, 0, 4, 2, 2, 0, 1, 4, 2, 1, 1, 1, 2, 0, 1, 2, 2, 2, 0, 0, 1, 2, 1, 2, 2, 3, 2, 3, 5, 3, 2, 1, 1, 1, 2, 2, 5, 0, 2, 2, 0, 2, 1, 0, 1, 2, 3, 1, 2, 2, 2, 4, 1, 0, 1, 4, 4, 2, 1, 4, 0, 1, 1, 1, 2, 1, 2, 0], 
			'A21137G-C21575T': [15, 13, 13, 16, 16, 12, 12, 9, 15, 10, 16, 18, 13, 21, 22, 19, 21, 10, 14, 15, 12, 9, 23, 21, 15, 12, 19, 12, 21, 10, 14, 14, 14, 17, 17, 7, 16, 12, 10, 15, 13, 12, 12, 15, 20, 17, 13, 15, 15, 16, 18, 11, 15, 18, 12, 17, 21, 13, 16, 14, 11, 18, 10, 14, 11, 11, 23, 17, 18, 16, 11, 14, 15, 13, 14, 6, 23, 10, 20, 13, 15, 16, 14, 18, 22, 22, 9, 23, 24, 18, 13, 11, 18, 17, 12, 9, 19, 18, 15, 17], 
			'G11083T-C28603T': [7, 10, 14, 14, 10, 18, 10, 12, 9, 10, 12, 13, 12, 16, 13, 15, 5, 6, 17, 19, 17, 12, 8, 6, 8, 11, 10, 12, 8, 12, 9, 12, 9, 8, 13, 12, 9, 12, 14, 14, 21, 9, 19, 9, 7, 16, 13, 18, 9, 14, 11, 11, 5, 9, 9, 12, 18, 13, 19, 12, 12, 6, 12, 11, 13, 12, 11, 11, 8, 14, 13, 10, 12, 9, 8, 17, 19, 12, 8, 11, 15, 12, 13, 12, 10, 26, 10, 5, 16, 13, 7, 16, 12, 13, 9, 13, 17, 11, 9, 13], 
			'C635T-C25521T': [2, 2, 1, 1, 2, 0, 2, 1, 2, 2, 1, 0, 2, 1, 0, 2, 3, 1, 3, 1, 2, 0, 1, 2, 3, 1, 3, 2, 1, 2, 2, 2, 2, 1, 0, 1, 1, 0, 1, 3, 1, 2, 1, 1, 0, 1, 1, 2, 1, 2, 1, 1, 3, 1, 2, 1, 3, 0, 0, 1, 1, 1, 3, 3, 1, 2, 2, 2, 3, 2, 3, 2, 0, 2, 2, 1, 0, 0, 2, 2, 1, 0, 3, 3, 1, 4, 3, 0, 1, 2, 1, 1, 3, 1, 1, 0, 0, 3, 0, 3], 
			'G11083T-C25521T': [13, 15, 10, 15, 12, 16, 16, 10, 15, 14, 14, 13, 13, 15, 12, 11, 12, 9, 19, 12, 9, 14, 14, 11, 12, 8, 8, 17, 16, 15, 10, 15, 10, 21, 10, 19, 11, 14, 6, 11, 16, 18, 17, 14, 13, 15, 13, 13, 15, 11, 13, 13, 7, 17, 21, 6, 13, 11, 8, 13, 13, 10, 7, 8, 14, 4, 12, 18, 17, 12, 17, 10, 19, 12, 11, 10, 10, 11, 14, 15, 10, 7, 6, 18, 9, 9, 11, 15, 14, 13, 16, 16, 17, 17, 17, 11, 19, 20, 17, 6], 
			'C3096T-G11083T': [17, 19, 17, 16, 17, 13, 14, 17, 13, 19, 11, 16, 15, 13, 14, 18, 9, 10, 13, 16, 13, 12, 16, 15, 11, 17, 16, 13, 11, 20, 12, 26, 8, 14, 12, 13, 10, 14, 19, 22, 19, 14, 17, 16, 16, 19, 11, 20, 13, 19, 17, 10, 12, 19, 12, 13, 16, 22, 13, 15, 21, 17, 10, 14, 17, 18, 15, 10, 16, 7, 13, 12, 14, 9, 21, 7, 17, 12, 15, 13, 19, 21, 15, 17, 15, 15, 19, 15, 12, 14, 11, 13, 18, 16, 14, 19, 20, 10, 13, 10], 
			'G11083T-C23638T': [6, 6, 9, 7, 4, 11, 3, 3, 11, 6, 6, 4, 8, 12, 7, 7, 10, 5, 6, 5, 4, 1, 15, 7, 6, 9, 6, 7, 11, 6, 5, 10, 9, 7, 9, 9, 7, 6, 6, 10, 8, 10, 5, 6, 7, 6, 5, 7, 9, 4, 8, 7, 7, 9, 3, 5, 4, 6, 7, 8, 10, 8, 10, 2, 7, 6, 6, 6, 7, 14, 9, 8, 8, 5, 2, 10, 4, 7, 9, 7, 11, 2, 10, 5, 8, 6, 4, 5, 8, 7, 6, 6, 8, 6, 6, 6, 5, 4, 7, 1], 
			'T28241A-T28243A': [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
			'C8140T-G11083T': [9, 5, 4, 4, 8, 5, 2, 4, 9, 6, 10, 4, 4, 8, 4, 4, 7, 6, 4, 5, 7, 4, 5, 9, 9, 5, 7, 3, 5, 7, 6, 7, 5, 3, 9, 7, 8, 4, 8, 7, 2, 2, 6, 7, 9, 1, 4, 7, 7, 7, 3, 5, 3, 6, 5, 5, 7, 8, 6, 4, 4, 9, 3, 11, 11, 4, 4, 7, 12, 7, 6, 4, 5, 2, 9, 3, 5, 8, 10, 6, 5, 6, 6, 6, 3, 8, 7, 6, 3, 9, 2, 4, 6, 7, 5, 3, 7, 7, 3, 9], 
			'T27257A-G27258A': [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
			'T28877A-C28878G': [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
			'C683T-C25521T': [0, 2, 2, 2, 1, 2, 2, 0, 2, 4, 2, 1, 4, 0, 2, 0, 0, 0, 1, 3, 0, 5, 1, 2, 0, 2, 1, 4, 1, 2, 0, 0, 3, 1, 0, 2, 0, 3, 1, 2, 0, 4, 2, 1, 0, 3, 2, 1, 0, 3, 3, 1, 0, 2, 1, 2, 4, 3, 0, 0, 2, 1, 1, 1, 3, 2, 1, 2, 1, 1, 0, 0, 0, 2, 6, 2, 3, 3, 0, 1, 2, 3, 1, 2, 0, 0, 2, 1, 0, 1, 1, 2, 1, 2, 3, 0, 0, 1, 3, 1], 
			'G11083T-C21855T': [4, 7, 1, 6, 8, 4, 4, 11, 6, 5, 6, 11, 8, 11, 6, 8, 7, 6, 9, 7, 4, 4, 6, 8, 10, 7, 4, 7, 6, 5, 6, 7, 8, 7, 1, 8, 9, 9, 5, 10, 8, 12, 8, 9, 9, 10, 8, 5, 5, 12, 7, 10, 4, 3, 8, 5, 17, 5, 9, 8, 4, 5, 5, 4, 4, 9, 8, 9, 11, 6, 7, 5, 8, 6, 4, 7, 4, 6, 3, 6, 5, 9, 10, 8, 9, 9, 4, 5, 9, 6, 8, 9, 9, 4, 8, 5, 5, 6, 4, 10], 
			'G11083T-C14724T': [14, 9, 12, 14, 15, 7, 11, 13, 12, 11, 17, 11, 6, 11, 15, 7, 8, 13, 9, 7, 8, 15, 11, 9, 8, 9, 15, 10, 12, 11, 7, 10, 12, 9, 12, 11, 13, 9, 12, 10, 8, 12, 13, 12, 12, 8, 12, 10, 14, 15, 10, 12, 8, 12, 8, 11, 7, 10, 9, 5, 14, 10, 8, 16, 13, 15, 9, 12, 7, 13, 11, 7, 12, 9, 13, 11, 10, 12, 13, 11, 8, 7, 9, 8, 7, 8, 12, 8, 9, 9, 9, 11, 8, 12, 13, 5, 21, 10, 7, 12], 
			'C9474T-G11083T': [7, 9, 2, 3, 3, 4, 9, 11, 6, 6, 6, 3, 9, 4, 6, 7, 4, 7, 3, 6, 8, 3, 4, 5, 3, 9, 3, 8, 7, 10, 4, 11, 8, 4, 9, 5, 8, 3, 9, 6, 7, 6, 6, 9, 7, 5, 8, 3, 6, 3, 6, 6, 6, 2, 4, 5, 7, 2, 4, 7, 3, 4, 5, 3, 6, 5, 7, 4, 5, 4, 5, 3, 4, 6, 2, 4, 5, 1, 6, 3, 10, 3, 5, 1, 5, 6, 6, 5, 12, 2, 8, 5, 8, 9, 1, 6, 3, 3, 4, 9], 
			'C8964T-G11083T': [6, 8, 11, 3, 4, 9, 7, 6, 10, 5, 4, 7, 14, 14, 6, 10, 4, 5, 4, 6, 10, 10, 12, 9, 6, 5, 10, 9, 9, 6, 5, 8, 2, 7, 9, 8, 8, 7, 7, 11, 11, 10, 2, 6, 7, 6, 6, 11, 9, 9, 11, 8, 9, 11, 7, 7, 10, 4, 8, 5, 9, 4, 3, 4, 13, 5, 10, 10, 6, 6, 4, 13, 6, 13, 4, 8, 8, 11, 5, 5, 5, 4, 7, 5, 10, 8, 9, 8, 8, 4, 9, 11, 10, 9, 10, 1, 6, 12, 5, 11], 
			'C1912T-C21575T': [5, 4, 8, 9, 8, 10, 12, 6, 8, 11, 8, 10, 7, 13, 7, 12, 10, 11, 7, 3, 13, 4, 17, 11, 7, 9, 12, 13, 5, 6, 9, 8, 5, 7, 8, 9, 13, 14, 5, 7, 10, 11, 6, 6, 8, 10, 6, 6, 6, 11, 9, 8, 12, 12, 5, 9, 12, 2, 11, 11, 7, 7, 12, 9, 7, 2, 8, 9, 9, 11, 15, 11, 8, 11, 10, 3, 9, 7, 8, 9, 7, 10, 12, 6, 10, 13, 13, 7, 5, 9, 10, 11, 15, 9, 14, 11, 6, 10, 7, 12], 
			'G26466T-G26467C': [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
			'G11083T-C24023T': [19, 6, 7, 12, 7, 10, 6, 13, 5, 6, 4, 13, 11, 15, 15, 8, 9, 9, 8, 5, 8, 5, 11, 5, 8, 7, 7, 8, 8, 6, 9, 8, 8, 10, 5, 7, 8, 12, 6, 5, 12, 7, 10, 8, 11, 12, 10, 6, 12, 13, 10, 10, 8, 13, 11, 7, 15, 7, 10, 4, 4, 6, 7, 9, 5, 8, 11, 6, 6, 8, 7, 5, 2, 10, 13, 8, 8, 9, 8, 5, 7, 5, 11, 9, 4, 4, 11, 6, 10, 11, 6, 9, 16, 9, 9, 11, 5, 11, 7, 8], 
			'G11083T-C25413T': [10, 7, 5, 9, 15, 16, 9, 10, 6, 10, 9, 6, 9, 7, 11, 8, 8, 8, 8, 3, 9, 9, 6, 10, 9, 7, 7, 14, 8, 4, 7, 10, 13, 13, 15, 10, 5, 8, 4, 5, 10, 8, 11, 11, 5, 9, 6, 13, 7, 12, 10, 11, 8, 9, 9, 15, 8, 8, 6, 10, 11, 10, 3, 12, 14, 6, 9, 5, 11, 5, 11, 13, 14, 6, 4, 7, 7, 7, 6, 7, 8, 8, 9, 9, 9, 9, 5, 11, 7, 9, 7, 9, 9, 13, 9, 10, 5, 8, 5, 12], 
			'G11083T-C28087T': [8, 7, 8, 7, 9, 11, 11, 14, 10, 4, 7, 4, 14, 6, 16, 8, 11, 9, 10, 7, 12, 8, 11, 5, 11, 12, 8, 12, 15, 9, 5, 8, 10, 8, 10, 12, 13, 10, 9, 11, 10, 11, 8, 8, 8, 10, 8, 15, 11, 7, 7, 11, 6, 13, 10, 12, 7, 11, 11, 12, 13, 10, 9, 12, 13, 10, 12, 10, 6, 12, 7, 11, 6, 11, 12, 12, 10, 7, 11, 11, 5, 7, 8, 9, 8, 14, 14, 7, 10, 13, 10, 18, 14, 10, 16, 4, 12, 4, 7, 10], 
			'C21575T-C25521T': [5, 9, 16, 7, 9, 13, 14, 8, 8, 6, 13, 6, 3, 5, 14, 8, 9, 19, 9, 7, 10, 9, 6, 9, 10, 8, 8, 13, 11, 11, 10, 3, 6, 5, 16, 10, 8, 8, 5, 7, 4, 10, 10, 6, 7, 7, 6, 4, 9, 7, 6, 6, 11, 5, 9, 11, 6, 6, 11, 11, 6, 10, 11, 9, 10, 6, 5, 8, 14, 9, 7, 8, 10, 4, 5, 6, 8, 9, 10, 7, 6, 9, 12, 3, 10, 17, 7, 6, 11, 9, 7, 6, 9, 12, 11, 14, 6, 10, 9, 7], 
			'G11083T-C27476T': [16, 16, 8, 3, 9, 15, 11, 10, 12, 12, 6, 12, 11, 9, 13, 11, 11, 15, 16, 7, 15, 10, 14, 13, 10, 7, 11, 8, 4, 11, 9, 7, 12, 12, 11, 9, 18, 21, 12, 8, 12, 9, 10, 8, 14, 10, 11, 7, 11, 7, 9, 6, 13, 15, 11, 10, 8, 15, 7, 11, 15, 12, 9, 9, 8, 12, 11, 11, 14, 5, 7, 10, 12, 12, 15, 14, 10, 12, 14, 8, 10, 18, 9, 11, 10, 11, 10, 8, 8, 12, 12, 12, 10, 6, 10, 16, 13, 10, 7, 8], 
			'C3096T-C9430T': [0, 0, 2, 4, 0, 2, 1, 2, 2, 0, 0, 0, 0, 1, 2, 1, 2, 2, 2, 0, 1, 1, 4, 1, 3, 2, 4, 3, 2, 1, 1, 2, 2, 1, 2, 1, 1, 0, 1, 3, 3, 1, 2, 1, 3, 3, 1, 3, 0, 2, 2, 2, 4, 1, 0, 1, 1, 3, 2, 2, 0, 2, 6, 3, 0, 3, 1, 2, 1, 1, 1, 2, 3, 2, 1, 1, 2, 2, 1, 0, 1, 1, 2, 0, 3, 2, 2, 2, 0, 1, 1, 0, 3, 1, 1, 1, 1, 3, 2, 1], 
			'G11083T-C17550T': [7, 4, 6, 8, 7, 9, 4, 5, 7, 11, 9, 8, 11, 13, 9, 4, 11, 8, 9, 7, 9, 9, 15, 9, 10, 7, 11, 8, 17, 5, 10, 7, 9, 5, 9, 7, 8, 13, 14, 8, 7, 10, 9, 9, 11, 5, 7, 13, 7, 9, 9, 8, 4, 8, 5, 11, 8, 15, 6, 7, 4, 8, 10, 6, 8, 8, 4, 8, 6, 9, 10, 9, 8, 7, 5, 6, 6, 11, 8, 8, 10, 9, 3, 14, 5, 4, 10, 10, 8, 8, 6, 9, 6, 1, 9, 4, 13, 12, 12, 3], 
			'C4331T-G11083T': [1, 5, 4, 7, 2, 3, 7, 3, 6, 6, 3, 6, 6, 6, 10, 9, 6, 9, 4, 5, 3, 4, 6, 8, 4, 2, 7, 5, 6, 10, 4, 11, 4, 6, 5, 4, 5, 3, 6, 3, 11, 6, 13, 5, 0, 4, 9, 4, 6, 8, 6, 6, 4, 4, 9, 2, 2, 6, 8, 4, 4, 10, 4, 10, 9, 7, 5, 4, 3, 10, 3, 3, 5, 7, 7, 11, 4, 5, 7, 5, 7, 2, 5, 8, 4, 7, 8, 7, 4, 2, 1, 13, 6, 6, 5, 5, 5, 5, 6, 5], 
			'C19524T-C21575T': [4, 8, 8, 7, 7, 6, 3, 6, 9, 8, 5, 6, 9, 10, 6, 10, 5, 4, 6, 4, 5, 3, 8, 6, 4, 12, 5, 7, 4, 4, 6, 10, 6, 8, 9, 8, 8, 10, 11, 10, 7, 3, 9, 7, 7, 3, 3, 3, 2, 9, 5, 7, 4, 7, 6, 5, 4, 4, 10, 1, 4, 3, 5, 6, 2, 3, 5, 8, 5, 4, 2, 6, 4, 4, 4, 7, 6, 4, 10, 8, 10, 8, 3, 4, 13, 4, 7, 7, 8, 5, 5, 5, 6, 7, 10, 6, 12, 9, 7, 5], 
			'C835T-C841A': [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
			'G11083T-C29733T': [6, 8, 13, 5, 9, 8, 6, 12, 15, 8, 9, 5, 13, 14, 7, 6, 14, 8, 6, 9, 5, 11, 8, 6, 11, 6, 11, 11, 8, 7, 12, 13, 9, 8, 9, 6, 8, 5, 6, 10, 8, 7, 8, 7, 7, 6, 10, 4, 9, 10, 11, 11, 10, 9, 10, 5, 12, 11, 10, 7, 7, 7, 11, 10, 4, 8, 10, 10, 13, 3, 8, 8, 8, 4, 9, 4, 8, 9, 6, 10, 5, 10, 5, 6, 9, 7, 10, 8, 8, 8, 14, 8, 9, 4, 7, 10, 4, 10, 11, 7], 
			'C5144T-G11083T': [8, 12, 6, 12, 6, 12, 4, 13, 12, 9, 5, 6, 11, 4, 6, 2, 12, 17, 12, 11, 12, 6, 7, 11, 0, 8, 8, 10, 6, 8, 8, 2, 5, 6, 6, 7, 9, 7, 12, 10, 5, 7, 8, 8, 4, 15, 14, 10, 5, 9, 10, 7, 12, 13, 10, 10, 8, 12, 10, 7, 15, 8, 7, 10, 6, 10, 4, 10, 8, 7, 13, 14, 12, 8, 6, 9, 4, 10, 8, 10, 7, 11, 4, 17, 10, 7, 10, 9, 4, 10, 9, 6, 7, 10, 11, 5, 7, 6, 13, 10], 
			'C7528T-G11083T': [2, 10, 5, 12, 4, 9, 6, 7, 7, 1, 5, 7, 12, 12, 6, 5, 7, 10, 7, 6, 2, 11, 3, 10, 10, 7, 5, 4, 8, 8, 4, 2, 14, 1, 2, 10, 5, 4, 12, 5, 2, 4, 8, 10, 7, 10, 6, 4, 12, 8, 10, 6, 7, 4, 7, 4, 10, 14, 7, 6, 8, 9, 10, 5, 2, 7, 6, 9, 7, 4, 7, 5, 12, 6, 5, 6, 5, 13, 11, 8, 8, 8, 3, 6, 8, 4, 7, 11, 12, 15, 9, 11, 12, 11, 5, 11, 11, 9, 9, 5], 
			'C21575T-C24023T': [4, 11, 2, 2, 8, 11, 5, 7, 8, 5, 10, 8, 4, 2, 2, 2, 4, 5, 4, 9, 10, 6, 4, 6, 7, 2, 10, 2, 3, 3, 3, 4, 6, 10, 8, 4, 2, 8, 6, 5, 5, 5, 6, 6, 4, 4, 2, 3, 2, 5, 9, 2, 3, 5, 1, 14, 5, 7, 7, 7, 6, 4, 4, 6, 4, 3, 4, 7, 4, 3, 6, 5, 8, 4, 9, 7, 6, 7, 3, 4, 6, 7, 8, 7, 6, 5, 5, 9, 2, 7, 8, 6, 6, 2, 4, 6, 9, 2, 6, 9], 
			'C1684T-G11083T': [5, 3, 8, 6, 11, 7, 5, 10, 13, 7, 7, 4, 14, 4, 8, 8, 6, 9, 7, 7, 12, 5, 6, 8, 11, 6, 8, 9, 9, 15, 13, 6, 8, 7, 15, 10, 11, 4, 8, 8, 12, 10, 8, 6, 6, 6, 8, 7, 6, 7, 9, 9, 10, 10, 12, 8, 7, 10, 6, 16, 6, 6, 10, 9, 11, 6, 8, 5, 4, 9, 11, 7, 10, 9, 7, 7, 12, 10, 12, 9, 10, 12, 10, 8, 8, 5, 7, 12, 9, 13, 10, 7, 9, 9, 11, 10, 9, 10, 6, 13], 
			'C9438T-G11083T': [10, 4, 15, 7, 6, 9, 5, 5, 5, 6, 7, 11, 5, 5, 6, 9, 5, 5, 5, 4, 5, 4, 6, 7, 7, 3, 7, 10, 7, 4, 6, 13, 8, 4, 6, 15, 4, 4, 5, 6, 5, 9, 11, 12, 1, 5, 6, 6, 6, 8, 8, 7, 6, 6, 7, 8, 5, 8, 2, 3, 9, 5, 9, 7, 4, 5, 9, 3, 7, 6, 6, 3, 5, 4, 6, 6, 6, 3, 3, 6, 3, 10, 4, 2, 5, 5, 4, 6, 10, 8, 5, 5, 4, 6, 6, 5, 8, 5, 8, 8], 
			'C9430T-C21614T': [0, 0, 2, 1, 4, 3, 0, 3, 0, 1, 3, 2, 1, 1, 1, 1, 2, 0, 1, 0, 1, 4, 0, 0, 0, 0, 1, 1, 0, 0, 2, 1, 1, 0, 2, 0, 0, 1, 0, 1, 1, 1, 0, 0, 1, 5, 3, 0, 1, 1, 6, 2, 0, 0, 0, 0, 1, 1, 4, 1, 3, 1, 1, 2, 4, 1, 2, 1, 2, 0, 1, 1, 2, 1, 0, 0, 0, 1, 0, 0, 1, 1, 1, 0, 3, 3, 1, 0, 1, 1, 1, 1, 1, 2, 3, 0, 0, 2, 2, 2], 
			'C8140T-T8141C': [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
			'C1684T-C21575T': [3, 1, 8, 6, 5, 4, 9, 5, 2, 6, 7, 10, 5, 6, 4, 6, 6, 8, 2, 6, 8, 7, 6, 9, 6, 6, 5, 8, 7, 9, 7, 12, 6, 3, 3, 8, 7, 7, 8, 10, 5, 7, 6, 4, 7, 8, 16, 2, 5, 3, 10, 8, 5, 6, 9, 8, 6, 8, 5, 8, 8, 5, 4, 7, 5, 6, 2, 7, 7, 4, 3, 2, 9, 11, 6, 4, 4, 9, 8, 4, 4, 3, 4, 6, 7, 2, 8, 11, 9, 4, 6, 4, 7, 5, 6, 8, 8, 4, 3, 4], 
			'C25797T-A25798C': [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
			'G11083T-G24933T': [8, 10, 7, 9, 9, 13, 10, 7, 11, 12, 12, 11, 5, 6, 15, 12, 10, 10, 10, 3, 7, 11, 8, 9, 9, 8, 15, 6, 8, 6, 7, 6, 10, 10, 7, 5, 2, 12, 10, 8, 9, 5, 12, 10, 11, 5, 13, 9, 5, 7, 13, 8, 12, 10, 9, 8, 6, 7, 3, 9, 5, 8, 10, 4, 7, 13, 5, 8, 13, 10, 5, 5, 10, 5, 16, 11, 5, 13, 8, 8, 9, 7, 6, 3, 8, 9, 13, 11, 13, 12, 9, 10, 14, 9, 4, 13, 3, 10, 14, 9], 
			'G11083T-C29200T': [5, 7, 12, 11, 14, 11, 3, 7, 5, 7, 5, 3, 3, 6, 7, 6, 11, 7, 7, 2, 7, 5, 7, 3, 5, 9, 4, 6, 3, 5, 12, 5, 4, 9, 8, 10, 7, 3, 6, 7, 10, 4, 5, 10, 5, 10, 5, 13, 6, 6, 6, 5, 11, 9, 7, 3, 6, 10, 4, 7, 7, 8, 7, 10, 7, 6, 5, 5, 4, 8, 11, 12, 8, 8, 13, 7, 9, 9, 6, 6, 10, 7, 9, 9, 7, 3, 8, 4, 6, 6, 5, 10, 5, 5, 9, 7, 5, 6, 6, 3], 
			'G11083T-C19983T': [13, 12, 17, 14, 22, 14, 11, 14, 12, 8, 12, 14, 10, 10, 15, 10, 14, 16, 14, 15, 13, 12, 12, 11, 15, 11, 12, 11, 11, 14, 13, 12, 14, 9, 11, 14, 12, 13, 15, 22, 5, 13, 15, 17, 17, 10, 12, 5, 16, 12, 13, 10, 9, 8, 9, 14, 13, 16, 12, 7, 16, 5, 18, 16, 14, 18, 13, 7, 17, 15, 8, 16, 12, 16, 15, 16, 14, 11, 14, 10, 18, 11, 15, 10, 16, 10, 15, 16, 8, 14, 18, 10, 8, 14, 18, 23, 19, 18, 24, 11], 
			'G11083T-C20178T': [10, 7, 6, 10, 6, 10, 7, 7, 11, 7, 13, 8, 8, 7, 5, 10, 10, 8, 8, 8, 14, 11, 10, 8, 9, 6, 8, 8, 11, 12, 10, 9, 16, 13, 4, 13, 13, 8, 15, 11, 13, 8, 13, 14, 14, 13, 8, 12, 9, 10, 8, 8, 9, 10, 8, 7, 12, 15, 15, 8, 9, 15, 10, 6, 7, 10, 7, 18, 9, 13, 8, 4, 10, 7, 8, 7, 11, 12, 12, 12, 15, 9, 9, 10, 13, 10, 6, 7, 5, 13, 9, 8, 15, 9, 10, 11, 8, 8, 12, 7], 
			'G11083T-C19524T': [11, 13, 7, 8, 4, 9, 9, 12, 7, 13, 11, 9, 6, 7, 7, 12, 6, 12, 10, 6, 10, 6, 13, 11, 9, 13, 10, 7, 9, 8, 10, 8, 10, 11, 12, 6, 14, 9, 10, 9, 12, 4, 12, 8, 12, 10, 13, 10, 14, 7, 8, 9, 11, 10, 18, 4, 12, 6, 6, 12, 4, 7, 7, 6, 7, 10, 5, 11, 11, 14, 13, 12, 10, 7, 8, 12, 11, 10, 9, 10, 9, 8, 9, 12, 9, 14, 11, 5, 7, 9, 3, 6, 5, 6, 13, 12, 11, 9, 8, 16], 
			'T25704C-C25708T': [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
			'C21762T-C21846T': [0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 0, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 2, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 1, 2, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 2, 0, 0, 0, 0, 0], 
			'C3096T-C21575T': [9, 11, 10, 10, 6, 7, 13, 8, 12, 5, 12, 12, 9, 12, 15, 15, 5, 9, 15, 13, 11, 9, 10, 9, 11, 6, 8, 12, 7, 5, 13, 9, 6, 4, 7, 11, 15, 5, 6, 12, 9, 8, 5, 7, 13, 13, 4, 8, 7, 8, 6, 7, 6, 11, 7, 7, 8, 9, 8, 2, 10, 7, 9, 10, 3, 9, 9, 8, 9, 6, 8, 8, 8, 13, 8, 7, 9, 6, 13, 4, 5, 10, 12, 10, 10, 3, 8, 7, 10, 6, 11, 8, 11, 10, 8, 14, 13, 8, 14, 12], 
			'G11083T-C17676T': [5, 10, 8, 9, 10, 10, 9, 6, 13, 6, 11, 10, 7, 12, 8, 10, 3, 16, 4, 8, 5, 10, 6, 9, 9, 7, 12, 11, 15, 4, 7, 9, 5, 6, 8, 9, 7, 7, 14, 12, 2, 6, 3, 9, 4, 4, 4, 15, 10, 9, 5, 16, 13, 11, 8, 7, 6, 11, 4, 8, 9, 5, 12, 7, 9, 10, 9, 9, 9, 10, 8, 5, 8, 10, 8, 5, 11, 6, 9, 10, 11, 11, 7, 7, 15, 8, 9, 6, 7, 5, 8, 11, 9, 7, 5, 7, 8, 8, 12, 8], 
			'G11083T-C29253T': [9, 8, 7, 9, 8, 10, 9, 7, 13, 9, 8, 13, 5, 5, 10, 5, 12, 8, 13, 11, 8, 9, 9, 9, 9, 5, 12, 9, 8, 7, 9, 8, 9, 6, 5, 8, 5, 3, 12, 9, 8, 5, 8, 10, 9, 10, 9, 6, 14, 10, 8, 7, 9, 7, 9, 15, 7, 6, 5, 5, 14, 8, 6, 12, 9, 12, 10, 12, 11, 7, 7, 7, 10, 9, 7, 9, 7, 11, 10, 6, 12, 7, 7, 6, 8, 9, 6, 10, 5, 12, 12, 8, 10, 5, 8, 7, 11, 11, 10, 7], 
			'C6781T-T6782C': [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
			'C21575T-C26111T': [6, 2, 3, 2, 7, 3, 12, 1, 2, 2, 3, 3, 3, 3, 3, 2, 1, 3, 7, 3, 1, 2, 6, 1, 3, 5, 2, 3, 1, 1, 2, 2, 5, 2, 5, 3, 2, 3, 4, 3, 6, 1, 6, 3, 1, 2, 3, 3, 4, 3, 6, 1, 4, 4, 0, 5, 3, 0, 3, 7, 3, 3, 3, 0, 3, 1, 4, 1, 3, 2, 4, 3, 4, 3, 5, 4, 3, 8, 3, 2, 2, 6, 5, 4, 6, 5, 3, 2, 1, 2, 3, 8, 0, 3, 1, 4, 0, 2, 3, 7], 
			'T12612C-G12614A': [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
			'C21575T-C27389T': [11, 5, 4, 7, 6, 15, 5, 10, 7, 7, 7, 11, 6, 6, 9, 10, 10, 10, 1, 13, 10, 12, 8, 5, 6, 7, 5, 12, 10, 7, 10, 8, 8, 16, 12, 11, 9, 13, 14, 10, 14, 5, 7, 8, 10, 13, 12, 2, 5, 8, 4, 6, 13, 3, 12, 11, 14, 6, 7, 8, 10, 8, 15, 8, 13, 7, 10, 4, 14, 5, 12, 12, 11, 11, 6, 10, 9, 8, 5, 9, 14, 14, 11, 5, 8, 5, 14, 12, 5, 15, 13, 14, 14, 4, 11, 16, 10, 7, 3, 12], 
			'G11083T-C13665T': [5, 4, 8, 6, 6, 12, 7, 7, 4, 7, 6, 5, 5, 8, 5, 6, 4, 7, 6, 6, 9, 6, 6, 3, 2, 10, 6, 8, 4, 8, 3, 12, 5, 7, 5, 6, 7, 3, 3, 6, 4, 8, 13, 4, 4, 12, 8, 3, 8, 4, 5, 5, 8, 8, 5, 5, 4, 6, 9, 10, 6, 17, 5, 5, 7, 8, 9, 6, 4, 6, 7, 12, 6, 4, 6, 11, 8, 12, 7, 9, 9, 5, 7, 9, 6, 9, 5, 6, 6, 9, 6, 3, 9, 6, 8, 8, 2, 5, 8, 10], 
			'G11083T-C18086T': [7, 5, 8, 2, 9, 9, 5, 4, 6, 6, 6, 7, 4, 2, 8, 4, 6, 6, 6, 12, 4, 10, 5, 9, 7, 6, 5, 5, 5, 7, 7, 6, 6, 9, 4, 12, 9, 2, 4, 3, 8, 3, 5, 6, 8, 7, 5, 3, 5, 8, 9, 4, 4, 3, 9, 6, 4, 6, 6, 5, 7, 4, 6, 9, 4, 11, 3, 3, 3, 7, 4, 5, 7, 9, 4, 1, 4, 9, 8, 5, 5, 4, 6, 5, 2, 8, 4, 2, 5, 8, 8, 5, 8, 4, 5, 6, 13, 7, 8, 10], 
			'A11782G-C21575T': [3, 4, 5, 5, 3, 8, 3, 3, 3, 5, 4, 5, 6, 8, 4, 2, 6, 2, 5, 6, 2, 7, 7, 3, 4, 5, 8, 9, 6, 5, 4, 3, 1, 7, 4, 7, 3, 3, 5, 0, 5, 4, 6, 4, 1, 4, 4, 6, 2, 4, 4, 5, 6, 1, 5, 0, 3, 2, 4, 3, 9, 5, 11, 4, 5, 7, 6, 4, 4, 5, 7, 5, 1, 5, 4, 4, 2, 5, 3, 1, 5, 5, 3, 3, 4, 2, 3, 7, 5, 4, 2, 4, 5, 4, 4, 8, 5, 3, 5, 6], 
			'C337T-G11083T': [13, 7, 10, 5, 10, 5, 7, 5, 6, 6, 7, 13, 14, 8, 9, 10, 6, 10, 10, 3, 9, 8, 4, 12, 8, 5, 6, 7, 4, 12, 4, 12, 11, 7, 10, 6, 8, 7, 9, 5, 6, 12, 6, 8, 11, 8, 6, 10, 9, 6, 7, 10, 6, 8, 5, 5, 8, 12, 8, 10, 5, 8, 6, 10, 11, 7, 8, 6, 5, 10, 8, 6, 12, 11, 7, 11, 5, 7, 6, 4, 7, 7, 7, 12, 6, 3, 8, 10, 12, 5, 6, 11, 9, 10, 9, 7, 5, 6, 12, 5], 
			'G11083T-G25217T': [3, 8, 2, 4, 2, 0, 6, 3, 5, 5, 3, 0, 4, 4, 1, 2, 3, 8, 4, 2, 7, 3, 6, 3, 9, 3, 3, 4, 5, 4, 3, 2, 1, 5, 5, 5, 2, 3, 4, 1, 1, 2, 5, 5, 5, 4, 4, 2, 4, 6, 1, 3, 1, 6, 4, 2, 2, 5, 5, 3, 3, 4, 4, 4, 3, 3, 4, 3, 3, 7, 5, 3, 4, 3, 2, 2, 3, 6, 0, 2, 5, 1, 3, 1, 2, 6, 3, 3, 3, 3, 7, 4, 2, 2, 3, 4, 4, 0, 7, 2], 
			'C1594T-G11083T': [5, 10, 13, 9, 8, 9, 14, 6, 9, 8, 7, 7, 4, 9, 7, 9, 4, 7, 8, 16, 12, 8, 7, 7, 2, 12, 10, 6, 4, 9, 4, 6, 9, 5, 9, 10, 11, 9, 9, 11, 12, 10, 8, 5, 8, 13, 10, 4, 6, 7, 9, 12, 10, 5, 12, 6, 6, 6, 7, 7, 12, 6, 8, 10, 6, 3, 10, 8, 8, 3, 12, 11, 12, 7, 5, 4, 12, 8, 6, 5, 3, 13, 8, 7, 6, 10, 7, 3, 4, 8, 7, 6, 5, 8, 9, 12, 11, 10, 4, 15], 
			'C11750T-C21575T': [10, 13, 13, 5, 11, 12, 19, 13, 13, 10, 9, 12, 8, 15, 8, 14, 10, 11, 10, 9, 5, 15, 9, 7, 12, 9, 3, 15, 12, 6, 14, 3, 6, 8, 14, 12, 8, 12, 7, 5, 5, 10, 11, 18, 9, 12, 3, 11, 15, 4, 8, 6, 9, 12, 5, 10, 8, 11, 9, 9, 11, 10, 6, 12, 9, 7, 14, 8, 7, 13, 8, 10, 14, 7, 10, 8, 9, 7, 7, 5, 10, 15, 7, 14, 11, 7, 13, 10, 9, 10, 6, 9, 15, 10, 13, 8, 13, 8, 10, 6], 
			'A10323G-C11750T': [5, 5, 6, 0, 8, 5, 6, 4, 4, 3, 4, 5, 2, 4, 5, 6, 1, 6, 13, 2, 3, 4, 6, 3, 3, 1, 4, 3, 4, 2, 4, 5, 4, 3, 4, 3, 2, 3, 4, 6, 4, 3, 4, 3, 3, 2, 5, 5, 5, 6, 3, 4, 6, 1, 2, 5, 5, 4, 2, 3, 3, 3, 5, 3, 3, 1, 7, 5, 6, 1, 3, 4, 1, 6, 5, 3, 4, 4, 3, 3, 3, 7, 2, 5, 7, 1, 5, 2, 3, 4, 6, 6, 6, 4, 4, 4, 4, 2, 6, 4], 
			'T8141C-G8144A': [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
			'C27881T-G27882C': [0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0], 
			'G11083T-C14599T': [4, 8, 6, 7, 9, 10, 8, 10, 11, 3, 13, 7, 8, 10, 12, 10, 5, 6, 3, 10, 6, 8, 3, 4, 13, 8, 3, 16, 5, 7, 12, 8, 7, 1, 7, 10, 8, 5, 8, 13, 6, 5, 5, 8, 10, 10, 8, 6, 7, 3, 8, 7, 8, 3, 7, 6, 8, 4, 2, 5, 5, 11, 5, 7, 4, 7, 11, 10, 6, 4, 6, 7, 9, 6, 6, 10, 3, 5, 13, 8, 10, 6, 6, 5, 5, 7, 11, 6, 8, 9, 11, 3, 3, 4, 5, 5, 0, 7, 5, 10], 
			'G28882A-G28883C': [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
			'C16887T-C25904T': [4, 3, 5, 5, 2, 1, 1, 3, 1, 3, 2, 1, 2, 3, 3, 2, 3, 2, 1, 4, 1, 4, 1, 2, 2, 1, 3, 3, 3, 4, 2, 3, 2, 7, 3, 1, 2, 2, 3, 2, 2, 3, 5, 6, 3, 2, 2, 2, 0, 6, 3, 2, 2, 1, 8, 5, 5, 3, 2, 5, 1, 0, 4, 3, 2, 5, 4, 1, 1, 2, 2, 0, 3, 3, 1, 2, 1, 2, 2, 3, 2, 1, 3, 2, 3, 1, 2, 5, 3, 1, 5, 2, 2, 4, 1, 1, 4, 5, 2, 1], 
			'T27752C-G27758A': [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
			'T5077C-A5078T': [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
			'C21575T-C29095T': [11, 7, 8, 7, 9, 7, 9, 12, 9, 9, 5, 12, 9, 10, 4, 10, 7, 13, 9, 7, 12, 16, 5, 5, 15, 13, 11, 8, 9, 8, 9, 17, 13, 11, 8, 10, 13, 7, 12, 13, 7, 10, 13, 12, 9, 9, 5, 13, 8, 7, 14, 6, 11, 12, 9, 4, 7, 9, 14, 10, 9, 6, 12, 7, 10, 9, 8, 7, 15, 5, 11, 5, 5, 8, 15, 10, 9, 5, 13, 9, 6, 11, 9, 10, 11, 6, 9, 10, 8, 10, 13, 8, 13, 8, 11, 15, 7, 9, 7, 7], 
			'C683T-G11083T': [10, 7, 6, 8, 13, 10, 8, 16, 4, 4, 15, 10, 11, 9, 10, 10, 7, 10, 4, 10, 13, 12, 12, 6, 13, 5, 5, 9, 11, 9, 10, 6, 9, 13, 15, 14, 10, 14, 10, 14, 10, 5, 10, 14, 10, 7, 8, 6, 7, 19, 12, 13, 13, 11, 7, 9, 13, 6, 8, 14, 7, 9, 12, 11, 15, 6, 11, 9, 8, 9, 11, 8, 5, 16, 16, 14, 9, 7, 4, 15, 11, 7, 6, 14, 9, 11, 6, 8, 11, 16, 10, 13, 14, 10, 8, 10, 15, 11, 9, 9], 
			'C9491T-G11083T': [6, 7, 7, 6, 9, 6, 1, 3, 8, 9, 5, 7, 8, 4, 7, 5, 6, 10, 5, 6, 4, 6, 9, 6, 8, 7, 13, 8, 8, 9, 7, 7, 9, 5, 10, 7, 3, 9, 6, 9, 8, 11, 2, 15, 10, 4, 5, 6, 5, 5, 10, 12, 6, 4, 6, 10, 3, 14, 9, 15, 4, 6, 7, 10, 6, 5, 6, 10, 6, 5, 7, 4, 7, 9, 6, 5, 6, 2, 11, 6, 10, 8, 9, 5, 11, 15, 5, 8, 7, 5, 14, 7, 5, 8, 6, 5, 3, 10, 3, 11], 
			'C4795T-G11083T': [9, 1, 2, 3, 4, 10, 8, 2, 5, 7, 6, 8, 10, 4, 6, 3, 6, 7, 6, 7, 6, 3, 9, 4, 10, 2, 5, 8, 5, 1, 9, 3, 4, 1, 6, 7, 6, 6, 4, 9, 4, 4, 5, 6, 9, 2, 4, 6, 2, 4, 7, 5, 9, 1, 7, 1, 3, 5, 4, 1, 5, 9, 5, 6, 6, 5, 6, 5, 5, 4, 4, 3, 7, 8, 4, 7, 7, 2, 7, 7, 6, 6, 4, 7, 3, 5, 8, 2, 6, 5, 7, 7, 4, 1, 6, 6, 1, 4, 5, 5], 
			'C27879T-A27880C': [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
			'T25629C-G25630C': [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
			'C16887T-A21137G': [6, 4, 12, 7, 8, 3, 3, 9, 9, 14, 2, 5, 7, 12, 10, 8, 6, 8, 9, 8, 7, 6, 4, 10, 6, 8, 7, 11, 9, 11, 12, 6, 7, 7, 6, 6, 3, 9, 7, 7, 4, 10, 7, 9, 4, 6, 7, 8, 7, 7, 4, 8, 3, 8, 8, 5, 5, 6, 6, 7, 6, 11, 5, 7, 6, 7, 6, 2, 11, 8, 12, 13, 9, 6, 7, 9, 5, 7, 8, 6, 5, 8, 2, 6, 15, 8, 5, 6, 5, 4, 9, 8, 8, 7, 7, 6, 8, 4, 5, 8], 
			'C635T-C9430T': [1, 3, 3, 0, 0, 3, 2, 3, 0, 1, 0, 1, 3, 1, 2, 3, 1, 0, 1, 1, 0, 0, 0, 2, 0, 1, 1, 1, 3, 1, 2, 1, 0, 2, 4, 2, 4, 3, 0, 1, 1, 1, 0, 0, 1, 2, 2, 1, 2, 1, 1, 7, 3, 2, 1, 2, 2, 1, 1, 1, 3, 1, 0, 1, 0, 0, 0, 0, 1, 0, 4, 1, 2, 1, 1, 2, 0, 2, 2, 2, 3, 0, 2, 1, 1, 1, 1, 2, 1, 1, 1, 2, 2, 1, 2, 3, 2, 4, 1, 3], 
			'G11083T-C27741T': [7, 9, 8, 8, 9, 5, 6, 4, 4, 5, 5, 5, 7, 7, 10, 7, 5, 8, 5, 7, 6, 8, 4, 9, 5, 6, 10, 9, 2, 9, 12, 6, 4, 7, 5, 7, 11, 6, 9, 8, 8, 6, 9, 8, 4, 15, 2, 3, 12, 6, 7, 4, 6, 10, 11, 7, 7, 9, 5, 6, 8, 5, 4, 9, 9, 8, 7, 14, 8, 7, 11, 8, 10, 7, 6, 4, 9, 10, 9, 7, 11, 5, 9, 8, 5, 8, 8, 6, 7, 11, 9, 11, 9, 8, 14, 10, 8, 5, 7, 7], 
			'C16887T-C27389T': [5, 2, 3, 7, 5, 2, 4, 1, 5, 4, 5, 4, 5, 6, 2, 5, 5, 3, 4, 5, 8, 5, 6, 3, 7, 4, 5, 4, 5, 5, 2, 5, 3, 2, 6, 4, 2, 4, 8, 4, 10, 0, 3, 5, 4, 2, 1, 3, 8, 2, 4, 2, 4, 8, 1, 3, 6, 5, 6, 6, 1, 3, 5, 4, 6, 4, 1, 9, 5, 3, 6, 5, 8, 3, 5, 1, 2, 7, 6, 3, 3, 5, 4, 5, 4, 4, 5, 6, 3, 6, 6, 3, 2, 7, 3, 3, 6, 5, 5, 4], 
			'C21575T-C28093T': [7, 6, 9, 9, 11, 6, 4, 9, 9, 8, 5, 4, 9, 7, 7, 8, 11, 8, 9, 6, 10, 6, 5, 10, 7, 7, 9, 8, 7, 8, 4, 8, 11, 8, 7, 5, 9, 4, 7, 10, 7, 5, 13, 11, 12, 4, 7, 8, 6, 10, 10, 10, 10, 6, 9, 8, 9, 4, 11, 8, 7, 5, 7, 11, 7, 8, 7, 7, 7, 9, 7, 6, 9, 8, 6, 11, 4, 6, 8, 4, 9, 12, 8, 7, 9, 9, 6, 7, 10, 10, 11, 14, 3, 8, 7, 12, 5, 7, 8, 4], 
			'C27739A-C27741T': [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
			'C203T-C635T': [3, 0, 1, 1, 1, 1, 1, 1, 2, 1, 2, 0, 3, 4, 1, 3, 2, 1, 0, 2, 3, 2, 1, 1, 2, 1, 0, 0, 0, 2, 1, 0, 0, 0, 1, 1, 1, 1, 1, 2, 1, 1, 2, 4, 0, 0, 4, 1, 1, 1, 1, 1, 2, 1, 1, 0, 0, 1, 1, 3, 1, 3, 1, 3, 1, 1, 2, 2, 1, 0, 0, 1, 1, 0, 1, 0, 1, 3, 2, 0, 1, 0, 0, 2, 1, 3, 0, 1, 3, 1, 1, 0, 1, 1, 4, 0, 2, 1, 0, 1], 
			'C18744T-C21575T': [10, 5, 4, 4, 9, 4, 5, 3, 6, 7, 3, 10, 6, 6, 4, 3, 5, 6, 12, 4, 5, 4, 7, 4, 9, 3, 2, 5, 6, 2, 2, 4, 4, 5, 7, 5, 3, 7, 4, 4, 6, 7, 5, 6, 9, 5, 6, 3, 6, 5, 4, 2, 4, 3, 5, 4, 7, 8, 5, 5, 9, 6, 8, 3, 5, 7, 5, 2, 7, 5, 8, 5, 3, 5, 8, 9, 10, 6, 5, 6, 12, 6, 13, 2, 2, 1, 3, 6, 5, 2, 7, 2, 8, 5, 9, 4, 6, 4, 5, 2], 
			'G24219A-G24220T': [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
			'G22335T-G22336T': [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
			'C9803T-G11083T': [5, 4, 2, 3, 4, 3, 5, 9, 7, 2, 8, 5, 4, 3, 8, 7, 3, 3, 8, 2, 3, 5, 8, 18, 6, 3, 8, 3, 5, 3, 3, 6, 6, 7, 3, 6, 6, 11, 6, 8, 7, 5, 5, 3, 3, 6, 6, 5, 11, 1, 3, 6, 4, 6, 7, 2, 6, 7, 6, 5, 3, 2, 6, 3, 4, 5, 5, 8, 3, 11, 3, 8, 10, 10, 5, 4, 4, 7, 4, 5, 6, 4, 3, 5, 5, 5, 3, 8, 6, 9, 1, 4, 7, 4, 5, 6, 1, 10, 4, 5], 
			'G11083T-C12459T': [8, 8, 11, 3, 6, 5, 4, 7, 5, 10, 4, 7, 2, 12, 5, 6, 5, 7, 4, 9, 8, 2, 5, 10, 3, 3, 10, 4, 7, 5, 8, 5, 8, 3, 7, 4, 8, 3, 5, 6, 4, 3, 5, 8, 9, 7, 8, 7, 6, 5, 4, 5, 3, 5, 3, 9, 7, 3, 3, 6, 7, 7, 11, 3, 11, 5, 5, 12, 8, 7, 6, 4, 8, 5, 5, 6, 5, 5, 9, 10, 1, 7, 2, 3, 6, 8, 5, 4, 1, 7, 5, 3, 5, 4, 7, 5, 10, 11, 8, 8], 
			'G11083T-C23277T': [7, 5, 3, 10, 8, 3, 6, 8, 4, 6, 0, 2, 6, 6, 10, 3, 8, 5, 3, 2, 3, 6, 4, 8, 7, 6, 9, 6, 3, 7, 4, 9, 3, 4, 6, 8, 6, 9, 4, 2, 4, 8, 4, 7, 7, 4, 9, 12, 7, 5, 8, 9, 5, 3, 2, 6, 7, 4, 3, 7, 3, 10, 9, 9, 7, 3, 8, 7, 9, 6, 3, 6, 9, 4, 5, 8, 8, 5, 6, 6, 4, 7, 6, 6, 7, 7, 1, 5, 6, 10, 8, 6, 9, 7, 3, 4, 4, 6, 7, 8], 
			'C6285T-G11083T': [2, 1, 1, 3, 5, 9, 4, 3, 6, 3, 8, 7, 7, 1, 5, 1, 7, 8, 6, 5, 5, 4, 5, 10, 3, 6, 1, 5, 2, 4, 3, 3, 3, 2, 4, 5, 5, 6, 7, 9, 3, 2, 7, 9, 1, 6, 3, 4, 4, 3, 2, 2, 4, 2, 1, 8, 2, 5, 3, 3, 8, 8, 1, 5, 5, 6, 2, 8, 5, 4, 2, 6, 5, 7, 1, 4, 1, 5, 2, 7, 1, 3, 5, 2, 6, 6, 7, 4, 3, 5, 7, 4, 7, 9, 9, 3, 6, 4, 4, 6], 
			'C21575T-G28079T': [6, 7, 2, 5, 2, 6, 9, 6, 5, 5, 8, 9, 3, 4, 7, 6, 6, 9, 3, 2, 6, 6, 4, 2, 2, 3, 3, 6, 2, 1, 4, 7, 6, 3, 3, 9, 2, 7, 0, 6, 3, 3, 0, 4, 5, 6, 4, 4, 9, 4, 3, 8, 6, 3, 4, 2, 4, 5, 7, 6, 3, 4, 5, 4, 5, 9, 5, 5, 2, 7, 6, 5, 2, 5, 8, 6, 8, 3, 7, 4, 3, 11, 4, 7, 6, 7, 5, 8, 5, 7, 7, 6, 3, 3, 7, 6, 4, 3, 3, 1], 
			'G21295A-G21296A': [0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
			'C21575T-T27384C': [6, 8, 10, 13, 10, 7, 10, 13, 10, 13, 15, 9, 12, 12, 9, 9, 11, 9, 9, 14, 5, 7, 14, 8, 10, 12, 6, 4, 11, 8, 9, 12, 10, 13, 16, 12, 16, 5, 14, 10, 6, 12, 15, 14, 17, 10, 13, 9, 8, 15, 16, 13, 14, 6, 13, 6, 12, 6, 20, 5, 12, 11, 10, 9, 17, 16, 13, 13, 12, 14, 16, 9, 9, 10, 10, 15, 11, 14, 12, 10, 11, 10, 8, 13, 10, 4, 12, 13, 15, 14, 7, 10, 16, 14, 10, 14, 11, 9, 5, 11], 
			'C21575T-C26681T': [12, 12, 9, 17, 12, 7, 9, 12, 14, 14, 6, 15, 10, 13, 11, 15, 7, 12, 14, 13, 10, 10, 9, 17, 10, 12, 13, 15, 9, 6, 12, 15, 18, 20, 7, 12, 11, 9, 21, 7, 11, 11, 12, 12, 7, 12, 12, 13, 13, 13, 14, 8, 16, 12, 19, 12, 9, 18, 9, 11, 14, 4, 13, 4, 9, 5, 13, 13, 13, 8, 15, 9, 16, 14, 14, 9, 8, 5, 12, 16, 20, 15, 8, 10, 11, 8, 11, 9, 10, 9, 12, 8, 15, 6, 5, 20, 13, 9, 16, 13], 
			'G11083T-C25916T': [11, 4, 6, 9, 9, 6, 7, 5, 7, 7, 13, 10, 7, 4, 5, 3, 7, 11, 10, 10, 3, 4, 7, 8, 7, 13, 5, 8, 12, 8, 7, 5, 12, 5, 9, 5, 9, 8, 9, 9, 11, 4, 8, 5, 7, 4, 10, 6, 7, 9, 11, 7, 9, 4, 4, 11, 11, 7, 7, 9, 7, 8, 6, 10, 9, 14, 5, 5, 3, 8, 8, 7, 8, 3, 7, 9, 8, 10, 10, 11, 9, 8, 5, 9, 9, 7, 12, 6, 9, 10, 11, 10, 8, 2, 7, 4, 8, 13, 7, 15], 
			'T27039A-C27040A': [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
			'G4583A-T4585A': [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]}

			import statistics
			for mutPair in simuDict.keys():
				minC=min(simuDict[mutPair])
				maxC=max(simuDict[mutPair])
				media=statistics.median(simuDict[mutPair])
				simuDict[mutPair]=(minC,media,maxC)
			print(simuDict)
			exit()
		else: # code for simulations
			import time 
			import random

			random.seed(1)

			recordedNumMutationPairs={}
			for comb in recordedCombs:
				recordedNumMutationPairs[comb]=[]
			numReps=100
			timeStep=time.time()
			mutationsListCopy=mutationsList.copy()
			for rep in range(numReps):
				numMutationPairs={}
				mutCounter=0
				random.shuffle(mutationsListCopy)
				counterBranch=0
				for branch in branchesList:
					mutations=[]
					poss=[]
					for i in range(branch):
						mut=mutationsListCopy[mutCounter]
						mutCounter+=1
						pos=int(mut[1:-1])
						posInsert=0
						while posInsert<len(poss) and poss[posInsert]<pos:
							posInsert+=1
						poss.insert(posInsert, pos)
						mutations.insert(posInsert, mut)
					
					if len(mutations)>1:

						foundComb=False
						for mutComb in recordedCombs:
							mutCombList=mutComb.split("-")
							found=True
							for mut in mutCombList:
								if not (mut in mutations):
									found=False
									break
							if found:
								foundComb=True
								if mutComb in numMutationPairs:
									numMutationPairs[mutComb]+=1
								else:
									numMutationPairs[mutComb]=1
								break
						
						#count other mutation pairs that are not under focus
						for i in range(len(mutations)):
							if (not foundComb) or ( not (mutations[i] in mutCombList)):
								for j in range(len(mutations)):
									if i<j:
										if (not foundComb) or ( not (mutations[j] in mutCombList)):
											mutationPairName=mutations[i]+"-"+mutations[j]
											if mutationPairName in numMutationPairs:
												numMutationPairs[mutationPairName]+=1
											else:
												numMutationPairs[mutationPairName]=1
					counterBranch+=1
					#if (counterBranch%100000)==0:
					#	print(counterBranch)
				#print("\n time taken:")
				#print(time.time()-timeStep)
				for comb in recordedCombs:
					if comb in numMutationPairs:
						recordedNumMutationPairs[comb].append(numMutationPairs[comb])
					else:
						recordedNumMutationPairs[comb].append(0)
				print("\n new round of results")
				for mutPair in numMutationPairs.keys():
					if numMutationPairs[mutPair]>9:
						mut1=mutPair.split("-")[0]
						mut2=mutPair.split("-")[1]
						pos1=int(mut1[1:-1])
						pos2=int(mut2[1:-1])
						if not (mutPair in recordedCombs):
							if abs(pos1-pos2)<10 and abs(pos1-pos2)>0:
								print(mutPair+"\t"+str(numMutationPairs[mutPair]))
							elif numMutationPairs[mutPair]>50 and abs(pos1-pos2)>0:
								print(mutPair+"\t"+str(numMutationPairs[mutPair]))

			print("\n time taken:")
			print(time.time()-timeStep)
			print()
			print(recordedNumMutationPairs)




	exit()









focusedMultiNucMuts={}
focusedMultiNucMuts["C21302T-C21304A-G21305A"]=[]
focusedMultiNucMuts["C21304A-G21305A"]=[]
#focusedMultiNucMuts["T28245G-T28251C-C28253T-A28254C"]=[]#problematic
#focusedMultiNucMuts["T28251C-C28253T-A28254C"]=[]#problematic
#focusedMultiNucMuts["T28251G-C28253T-A28254C"]=[]#problematic
#focusedMultiNucMuts["C28253T-A28254C"]=[]#problematic
#focusedMultiNucMuts["A28249T-C28253A"]=[]#problematic
focusedMultiNucMuts["A28877T-G28878C"]=[]
focusedMultiNucMuts["G27382C-A27383T-T27384C"]=[]
focusedMultiNucMuts["T26491C-A26492T-T26497C"]=[]
focusedMultiNucMuts["G27758A-T27760A"]=[]
focusedMultiNucMuts["T27875C-C27881T-G27882C-C27883T"]=[]
focusedMultiNucMuts["C27881T-G27882C-C27883T"]=[]
focusedMultiNucMuts["C25162A-C25163A"]=[]
focusedMultiNucMuts["T21294A-G21295A-G21296A"]=[]
focusedMultiNucMuts["A27038T-T27039A-C27040A"]=[]
focusedMultiNucMuts["A507T-T508C-G509A"]=[]#problematic
focusedMultiNucMuts["T28881A-G28882A-G28883C"]=[]#problematic
focusedMultiNucMuts["A21550C-A21551T"]=[]
focusedMultiNucMuts["T21994C-T21995C"]=[]#problematic
focusedMultiNucMuts["C13423A-C13424A"]=[]
focusedMultiNucMuts["A4576T-T4579A"]=[]
focusedMultiNucMuts["A20284T-T20285C"]=[]
focusedMultiNucMuts["G11083T-C21575T"]=[] #control
mutCombs=focusedMultiNucMuts.keys()
#rarer ones
rare3={}
rare3["G910A-T911A-C912A"]=[]
rare3["T27381C-G27382T-A27383G"]=[]
rare3["A5703T-G5704A-T5705A"]=[]
rare3["A3684T-G3685A-C3686A"]=[]
rare3["A27400T-A27403C-C27406G"]=[]
rare3["G28280C-A28281T-T28282A"]=[]
rare3["A4420G-C4421A-G4422A"]=[]
rare3["A28131G-C28132A-G28134A"]=[]
rare3["G29757A-T29758C-G29759C"]=[]
rare3["T21042G-A21043T-G21044C"]=[]
rare3["C25572T-A25573C-A25574T"]=[]
rare3["C21302T-T21304A-G21305A"]=[]
mutCombsRare=rare3.keys()

# G27382C-A27383T-T27384C	253	135	75	1724
# G27382C-A27383T	46	135	75

# C21302T-C21304A-G21305A	284	219	369	290
# C21304A-G21305A	452	369	290
# C21302T-C21304A	52	219	369

# T26491C-A26492T-T26497C	160	21	36	137
# T26491C-A26492T	14	21	36

# C22716A-T22717C	45	48	52
# T27672A-C27673A	37	47	179
# G6975T-G6977A	33	46	62
# G6513A-T6515A	33	46	34
# C19977A-C19979A	32	40	57
# T22207G-T22209C	31	53	43
# A27865T-T27866A	29	41	33
# A21892G-G21893A	28	44	72
# T26485C-T26486A	27	45	31
# G11071C-C11074T	26	79	1355
# C17734T-T17735C	25	94	26
# T27299C-A27300G	24	140	27
# C11950T-C28472T	21	172	232 #only one made of non-contiguous rare mutations
# G23957T-T23959G	20	22	23
# A25562T-G25563A	20	43	51
# T27760A-T27761C	20	60	61
# T3370G-G3371A	19	25	36
# A22194G-T22196G	18	88	24
# C25380T-A25381C	15	116	27
# G29688C-G29692C	15	81	48
# G21082C-A21083T	15	18	15
# A6024C-C6027T	15	20	570
# A26709G-T26767C	15	25	55
# T967C-G970A	14	70	83
# T16884G-T16885A	14	18	14
# T6148A-G6149A	14	34	30
# T28241A-T28243A	13	13	14


# investigate cherry pairs on the cluster to check for heterozygosity, indels, coverage
# needs first copying the cherry list to the cluster: 
if cherries:
	import gzip
	cherryFile=open("siblingsWithSingletonCombos.txt")
	results={}
	line=cherryFile.readline()
	numCherries=0
	while line!="" and line!="\n":
		linelist=line.split()
		numCherries+=1
		if (numCherries%100)==0:
			print(numCherries)
		nameMut=linelist[0] 
		if not (nameMut in results):
			results[nameMut]=[0,0,0,0,0,0,0,0]
		samples=[linelist[1],linelist[2]]
		pos=int(nameMut.split("-")[0][1:-1])
		for numSample in range(2):
			name=samples[numSample]
		
			file=None
			if name[0]=="E":
				try:
					file=gzip.open("/nfs/research/zi/mhunt/Viridian_wf_paper/Vdn_all_ena/Reads/"+name[0]+"/"+name[1:6]+"/"+name[6:8]+"/"+name[8:10]+"/vdn.v1.0.0/qc.tsv.gz",'rt')
				except:
					try:
						file=gzip.open("/nfs/research/zi/mhunt/Viridian_wf_paper/Vdn_all_ena/Reads/"+name[0]+"/"+name[1:7]+"/"+name[7:9]+"/"+name[9:11]+"/vdn.v1.0.0/qc.tsv.gz",'rt')
					except:
						print("Sample file could not be opened/found "+name)
			elif name[0]=="S":
				file=gzip.open("/nfs/research/zi/mhunt/Viridian_wf_paper/Vdn_all_ena/Reads/"+name[0]+"/"+name[1:7]+"/"+name[7:9]+"/"+name[9:11]+"/vdn.v1.0.0/qc.tsv.gz",'rt')
			elif name[0]=="D":
				file=gzip.open("/nfs/research/zi/mhunt/Viridian_wf_paper/Vdn_all_ena/Reads/"+name[0]+"/"+name[1:5]+"/"+name[5:7]+"/"+name[7:9]+"/vdn.v1.0.0/qc.tsv.gz",'rt')
			else:
				print("Sample file could not be opened/found "+name)
			if file!=None:
				results[nameMut][numSample*4]+=1
				foundIndel=False
				foundLowCov=False
				foundHet=False
				line=file.readline()
				for po in range(pos-20):
					line=file.readline()
				for po in range(45):
					line=file.readline()
					linelist=line.split()
					if linelist[1]=="-" or linelist[3]=="-" or linelist[4]=="-":
						foundIndel=True
					elif linelist[3]=="N" or linelist[4]=="N" or int(linelist[9])<100:
						foundLowCov=True
					elif ((int(linelist[9])-int(linelist[10]))/float(linelist[9])) >0.05:
						foundHet=True
				if foundIndel:
					results[nameMut][1+numSample*4]+=1
				if foundLowCov:
					results[nameMut][2+numSample*4]+=1
				if foundHet:
					results[nameMut][3+numSample*4]+=1
			else:
				print(line)
				print("Sample file could not be opened/found "+name)
				print("")
			
				

		

		line=cherryFile.readline()
	cherryFile.close()
	print(results)
	print(numCherries)
	exit()


if createCherryAlignments:
	import gzip
	cherryFile=open("siblingsWithSingletonCombos.txt")
	summaryAlignmentFile=open("multi-nucleotide_summaryAlignment.txt","w")
	line=cherryFile.readline()
	numCherries=0
	while line!="" and line!="\n":
		linelist=line.split()
		numCherries+=1
		if (numCherries%100)==0:
			print(numCherries)
		nameMut=linelist[0] 
		samples=[linelist[1],linelist[2]]
		pos=int(nameMut.split("-")[0][1:-1])
		names=[">"+linelist[0]+"_"+samples[0]+"\n",">"+linelist[0]+"_"+samples[1]+"\n"]
		seqs=["",""]
		files=[]
		for numSample in range(2):
			name=samples[numSample]
			if name[0]=="E":
				try:
					file=gzip.open("/nfs/research/zi/mhunt/Viridian_wf_paper/Vdn_all_ena/Reads/"+name[0]+"/"+name[1:6]+"/"+name[6:8]+"/"+name[8:10]+"/vdn.v1.0.0/qc.tsv.gz",'rt')
				except:
					try:
						file=gzip.open("/nfs/research/zi/mhunt/Viridian_wf_paper/Vdn_all_ena/Reads/"+name[0]+"/"+name[1:7]+"/"+name[7:9]+"/"+name[9:11]+"/vdn.v1.0.0/qc.tsv.gz",'rt')
					except:
						print("Sample file could not be opened/found "+name)
			elif name[0]=="S":
				file=gzip.open("/nfs/research/zi/mhunt/Viridian_wf_paper/Vdn_all_ena/Reads/"+name[0]+"/"+name[1:7]+"/"+name[7:9]+"/"+name[9:11]+"/vdn.v1.0.0/qc.tsv.gz",'rt')
			elif name[0]=="D":
				file=gzip.open("/nfs/research/zi/mhunt/Viridian_wf_paper/Vdn_all_ena/Reads/"+name[0]+"/"+name[1:5]+"/"+name[5:7]+"/"+name[7:9]+"/vdn.v1.0.0/qc.tsv.gz",'rt')
			else:
				print("Sample file could not be opened/found "+name)
				break
			files.append(file)
		if len(files)>1:
			line0=files[0].readline()
			line0=files[0].readline()
			line1=files[1].readline()
			line1=files[1].readline()
			while int(line0.split()[0])<(pos-100):
				line0=files[0].readline()
			while int(line1.split()[0])<(pos-100):
				line1=files[1].readline()
			currPos=int(line0.split()[0])
			while currPos<(pos+105):
				linelist0=line0.split()
				linelist1=line1.split()
				currPos=int(line0.split()[0])
				while linelist0[1]==".":
					seqs[0]+=linelist0[4]
					if linelist0[1]==".":
						seqs[1]+=linelist1[4]
						line1=files[1].readline()
						linelist1=line1.split()
					else:
						seqs[1]+="-"
					line0=files[0].readline()
					linelist0=line0.split()
				while linelist1[1]==".":
					seqs[1]+=linelist1[4]
					line1=files[1].readline()
					linelist1=line1.split()
					seqs[0]+="-"
				seqs[0]+=linelist0[4]
				seqs[1]+=linelist1[4]
				line0=files[0].readline()
				line1=files[1].readline()

			summaryAlignmentFile.write(names[0]+seqs[0]+"\n"+names[1]+seqs[1]+"\n")





		else:
			print(line)
			print("Sample file could not be opened/found "+name)
			print("")

		line=cherryFile.readline()
	cherryFile.close()
	summaryAlignmentFile.close()
	exit()



if recombinationFile!="":
	file=open(recombinationFile)
	line=file.readline()
	line=file.readline()
	numClusterFlags=0
	numRecombinations=0
	clustersFound=0
	numClusterFlagsAndClusterFound=0
	while line!="" and line!="\n":
		numRecombinations+=1
		linelist=line.split("\t")
		infoSeq=linelist[17]
		sites=linelist[22].split(",")
		flags=linelist[23].split(",")
		if "cluster" in flags:
			clusterFlag=True
			numClusterFlags+=1
		else:
			clusterFlag=False
		parsimonyImprovement=int(linelist[21])
		sitesA=[]
		sitesB=[]
		numA=0
		numB=0
		for site in range(len(sites)):
			if infoSeq[site]=="A":
				sitesA.append(sites[site])
				numA+=1
			else:
				sitesB.append(sites[site])
				numB+=1
		otherSiteList=False
		if numA<6 and numB<6:
			otherSiteList
		if numA<=numB:
			sitesList=sitesA
			if otherSiteList:
				otherSiteList=sitesB
		else:
			sitesList=sitesB
			if otherSiteList:
				otherSiteList=sitesA
		for mutComb in mutCombs:
			mutCombList=mutComb.split("-")
			found=True
			for mutN in range(len(mutCombList)):
				mut=mutCombList[mutN]
				site=mut[1:-1]
				if not (site in sitesList):
					found=False
					break
			if (not found) and otherSiteList:
				foundOther=True
				for mutN in range(len(mutCombList)):
					mut=mutCombList[mutN]
					site=mut[1:-1]
					if not (site in otherSiteList):
						foundOther=False
						break
			else:
				foundOther=False

			if found :
				focusedMultiNucMuts[mutComb].append((len(sitesList),len(mutCombList),parsimonyImprovement,clusterFlag))
				clustersFound+=1
				if clusterFlag:
					numClusterFlagsAndClusterFound+=1
				if linelist[23]=="PASS":
					print((len(sitesList),len(mutCombList),parsimonyImprovement,clusterFlag))
					print(line)
					print()
				break
			if foundOther :
				print("Found other")
				focusedMultiNucMuts[mutComb].append((len(otherSiteList),len(mutCombList),parsimonyImprovement,clusterFlag))
				clustersFound+=1
				if clusterFlag:
					numClusterFlagsAndClusterFound+=1
				break
		line=file.readline()

	print("Total number of recombinations: "+str(numRecombinations))
	print("Recombinations flagged as clusters: "+str(numClusterFlags))
	print("Recombinations in which clusters contributed: "+str(clustersFound))
	print("Recombinations in which clusters contributed that were also flagged as clusters: "+str(numClusterFlagsAndClusterFound))
	file.close()

	#1st entry: the recomb is made entirely of the cluster
	#2nd entry: the recomb is the cluster plus one more substitution
	#4th entry: the cluster contributed to the recomb
	allClusterRcombs=[0,0,clustersFound]
	multiNucMutsRecombs={}
	for mutComb in mutCombs:
		multiNucMutsRecombs[mutComb]=[0,0,len(focusedMultiNucMuts[mutComb])]
		for recomb in focusedMultiNucMuts[mutComb]:
			if recomb[0]==recomb[1]:
				allClusterRcombs[0]+=1
				multiNucMutsRecombs[mutComb][0]+=1
			elif recomb[0]==(recomb[1]+1):
				allClusterRcombs[1]+=1
				multiNucMutsRecombs[mutComb][1]+=1

	print("Found "+str(allClusterRcombs[0])+" completely contributing clusters and "+str(allClusterRcombs[1])+" almost completely contributing clusters out of total contributing clusters "+str(allClusterRcombs[2]))
	print(multiNucMutsRecombs)

	exit()




file=open(inputAl)
line=file.readline()
line=file.readline()
ref=""
while line!="" and line!="\n" and line[0]!=">":
	ref+=line.replace("\n","")
	line=file.readline()
ref=ref.upper()
file.close()


class Tree(object):
	def __init__(self):
		self.dist = []
		self.children = []
		self.mutations=[]
		self.up=[]
		self.name=[]
		self.nDesc=[]
	def __repr__(self):
		return "Tree object"
	def addNode(self):
		self.up.append(None)
		self.children.append([])
		self.name.append("")
		self.mutations.append([])
		self.dist.append(0.0)
		self.nDesc.append(0)

defaultBLen=0.000033

#function to read input newick string
def readNewick(nwFile):
	phyloFile=open(nwFile)
	line=phyloFile.readline()
	while line!="":
		while line=="\n":
			line=phyloFile.readline()
		if line=="":
			break
		tree=Tree()
		tree.addNode()
		nwString=line.replace("\n","")
		index=0
		nodeIndex=len(tree.name)-1
		name=""
		distStr=""
		finished=False
		while index<len(nwString):
			if nwString[index]=="(":
				tree.children[nodeIndex].append(len(tree.up))
				tree.addNode()
				tree.up[-1]=nodeIndex
				nodeIndex=len(tree.up)-1
				index+=1
			elif nwString[index]==";":
				phyloFile.close()
				return tree, nodeIndex
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
					tree.name[nodeIndex]=name
					name=""
				if distStr!="":
					tree.dist[nodeIndex]=float(distStr)
					if tree.dist[nodeIndex]<0.0:
						print("Warning: negative branch length in the input tree: "+distStr+" ; converting it to positive.")
						tree.dist[nodeIndex]=abs(tree.dist[nodeIndex])
					distStr=""
				else:
					tree.dist[nodeIndex]=defaultBLen
				nodeIndex=tree.up[nodeIndex]
				tree.children[nodeIndex].append(len(tree.up))
				tree.addNode()
				tree.up[-1]=nodeIndex
				nodeIndex=len(tree.up)-1
				index+=1
			elif nwString[index]==")":
				if name!="":
					tree.name[nodeIndex]=name
					name=""
				if distStr!="":
					tree.dist[nodeIndex]=float(distStr)
					distStr=""
				else:
					tree.dist[nodeIndex]=defaultBLen
				index+=1
				nodeIndex=tree.up[nodeIndex]
			else:
				name+=nwString[index]
				index+=1
		if not finished:
			print("Error, final character ; not found in newick string in file "+nwFile+".")
			raise Exception("exit")

tree, root=readNewick(inputTree)

print("Input tree read")

#calculate number of descendants for each node.
def calculateNDesc(tree,node):
	children=tree.children
	nDesc=tree.nDesc
	nextLeaves=[node]
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
			nDesc[nextNode]=1

#find singleton-sibling sample pairs
def findSingletonPairs(tree,node,singletonSamplesForCombos,fileOut):
	children=tree.children
	up=tree.up
	name=tree.name
	nextLeaves=[node]
	countSamples=0
	while nextLeaves:
		nextNode=nextLeaves.pop()
		if children[nextNode]:
			for c in children[nextNode]:
				nextLeaves.append(c)
		else:
			if name[nextNode] in singletonSamplesForCombos:
				if children[up[nextNode]][0]==nextNode:
					childNum=0
				else:
					childNum=1
				if not (children[children[up[nextNode]][1-childNum]]):
					fileOut.write(singletonSamplesForCombos[name[nextNode]]+"\t"+name[nextNode]+"\t"+name[children[up[nextNode]][1-childNum]]+"\n")
					countSamples+=1
	print("Total pairs: "+str(countSamples))

calculateNDesc(tree, root)

print("Number of descendants calculated")

#create dictionary with node names, for each recording the index in the tree lists.
nameDict={}
index=0
for name in tree.name:
	nameDict[name]=index
	index+=1

#'T26491C-A26492T-T26497C': [51, 51, 0, 7, 51, 0, 0, 1], 1-nuc deletion part of the mutation. 
#'C27881T-G27882C-C27883T': [23, 23, 0, 4, 23, 0, 0, 0], 2-nuc deletion as part of the mutation event, but rarely 6-nuc deletion also observed, and also 2-nuc insertion.
#'A28877T-G28878C': [131, 1, 2, 13, 131, 1, 0, 2], sample ERR4597661 has substantial heterozygosity, but seems like an isolated case.
#'T21994C-T21995C': [28, 7, 26, 0, 28, 24, 21, 1], TODO mark as unreliable. Sudden drop in coverage in this area, with some Viridian warnings. In sibling, low coverage and deletion at 21991-21993; it looks like a deletion at these positions in Alpha together with amplicon dropping causes an alignment error and 3 substitutions to replace the indel. 
#'G27382C-A27383T-T27384C': [69, 0, 1, 5, 69, 1, 4, 0],
#'C21304A-G21305A': [97, 0, 4, 6, 97, 0, 11, 4], 
#'G27758A-T27760A': [37, 0, 1, 6, 37, 0, 3, 0], 
#'T21294A-G21295A-G21296A': [21, 0, 1, 0, 21, 0, 1, 1], 
#'C21302T-C21304A-G21305A': [78, 0, 5, 12, 78, 0, 4, 0],
#'G11083T-C21575T': [14, 0, 0, 14, 14, 0, 0, 3], heterozygosity at positions 11074-5? Not important 
#'A20284T-T20285C': [14, 0, 0, 0, 14, 0, 0, 0] 
#'A507T-T508C-G509A': [32, 32, 2, 32, 32, 9, 0, 4],, TODO mark as unreliable. Total mess of low coverage, primers, indels, etc. Nothing in the siblings.
#'A4576T-T4579A': [14, 0, 4, 0, 14, 0, 4, 0], 
#'C13423A-C13424A': [24, 0, 0, 0, 24, 0, 0, 1], 
#'C25162A-C25163A': [29, 0, 0, 3, 29, 0, 1, 2],,
#'A27038T-T27039A-C27040A': [26, 0, 3, 6, 26, 0, 0, 0],
#'A21550C-A21551T': [16, 16, 0, 3, 16, 0, 1, 1], deletion of 3 nuc at the same time as substitutions
#'T28881A-G28882A-G28883C': [16, 12, 0, 3, 16, 0, 0, 0], TODO mark as unreliable. These are mostly short deletions that are masked into the reference allele.
#'T27875C-C27881T-G27882C-C27883T': [9, 9, 0, 0, 9, 0, 0, 0], 2-nuc deletion as 27881-27883, usually in C27874T background, hence enrichment in Delta.


#read input metadata file with inferred mutation events.
file=open(inputTSV)
line=file.readline()
line=file.readline()
thresholdsDesc=[1,2,5,10] #number of descendants required to count the mutation in each bin
inversions={}
inversionsNums={}
for mutComb in mutCombs:
	mutCombList=mutComb.split("-")
	inversion=""
	for mutN in range(len(mutCombList)):
		mut=mutCombList[mutN]
		inversion+=mut[-1]
		inversion+=mut[1:-1]
		inversion+=mut[0]
		if mutN<len(mutCombList)-1:
			inversion+="-"
	inversions[inversion]=mutComb
	inversionsNums[mutComb]=[]
inversions["A28881G-A28882G-C28883G"]="T28881A-G28882A-G28883C"
anyMutNDescDict={}
numMutations={}
numSingletons={}
numMutationPairs={}
singletonSamplesForCombos={}
fileSamplesWithSingletonCombos=open("siblingsWithSingletonCombos.txt","w")
while line!="\n" and line!="":
	linelist=line.split("\t")
	mutations=linelist[6]
	nodeName=linelist[0]
	if mutations!="":
		mutationlist=mutations.split(",")
		passedMutations=[]
		for mutation in mutationlist:
			mutationPair=mutation.split(":")
			support=float(mutationPair[1])
			if support>=thresholdProb and tree.nDesc[nameDict[nodeName]]>=minNumDescendants:
				mutationName=mutationPair[0]
				# if tree.nDesc[nameDict[nodeName]] in anyMutNDescDict:
				# 	anyMutNDescDict[tree.nDesc[nameDict[nodeName]]]+=1
				# else:
				# 	anyMutNDescDict[tree.nDesc[nameDict[nodeName]]]=1
				# if mutationName in numMutations:
				# 	numMutations[mutationName]+=1
				# else:
				# 	numMutations[mutationName]=1
				passedMutations.append(mutationName)
		foundComb=False
		if len(passedMutations)>1:
			foundComb=False
			for mutComb in mutCombs:
				mutCombList=mutComb.split("-")
				found=True
				for mut in mutCombList:
					if not (mut in passedMutations):
						found=False
						break
				if found:
					foundComb=True
					focusedMultiNucMuts[mutComb].append(tree.nDesc[nameDict[nodeName]])
					if tree.nDesc[nameDict[nodeName]]==1:
						#fileSamplesWithSingletonCombos.write(mutComb+"\t"+nodeName+"\n")
						singletonSamplesForCombos[nodeName]=mutComb
						#print(mutComb+"\t"+nodeName)
						#countSamples+=1
					break
			
			#check for reversions of multi-nucleotide mutations under focus
			for mutCombInv in inversions:
				mutCombListInv=mutCombInv.split("-")
				found=True
				for mut in mutCombListInv:
					if not (mut in passedMutations):
						found=False
						break
				if found:
					inversionsNums[inversions[mutCombInv]].append(tree.nDesc[nameDict[nodeName]])
			
			#count other mutation pairs that are not under focus
			for i in range(len(passedMutations)):
				if (not foundComb) or ( not (passedMutations[i] in mutCombList)):
					for j in range(len(passedMutations)):
						if i<j:
							if (not foundComb) or ( not (passedMutations[j] in mutCombList)):
								mutationPairName=passedMutations[i]+"-"+passedMutations[j]
								if mutationPairName in numMutationPairs:
									numMutationPairs[mutationPairName]+=1
								else:
									numMutationPairs[mutationPairName]=1

		#count individual single-nucleotide mutations
		for mutationName in passedMutations:
			if (not foundComb) or ( not (mutationName in mutCombList)):
				if tree.nDesc[nameDict[nodeName]] in anyMutNDescDict:
					anyMutNDescDict[tree.nDesc[nameDict[nodeName]]]+=1
				else:
					anyMutNDescDict[tree.nDesc[nameDict[nodeName]]]=1
				if mutationName in numMutations:
					numMutations[mutationName]+=1
				else:
					numMutations[mutationName]=1
				if tree.nDesc[nameDict[nodeName]]==1:
					if mutationName in numSingletons:
						numSingletons[mutationName]+=1
					else:
						numSingletons[mutationName]=1
			
	line=file.readline()
file.close()

findSingletonPairs(tree,root,singletonSamplesForCombos,fileSamplesWithSingletonCombos)
fileSamplesWithSingletonCombos.close()

#read alignment to find out how prevalent are 1-nuc variants and multi-nuc ones
file=open(inputAl)
line=file.readline()
line=file.readline()
line=file.readline()
numSamplesComb={}
numSamplesAll={}
for mutComb in mutCombs:
	numSamplesComb[mutComb.lower()]=0
mutCombsLower=numSamplesComb.keys()
while line!="" and line!="\n":
	variants={}
	line=file.readline()
	found28881=True
	while line!="" and line!="\n" and line[0]!=">":
		linelist=line.split()
		if linelist[0] in nucleotidesDict:
			variants[linelist[1]+linelist[0]]=True
		if int(linelist[1])>28880 and int(linelist[1])<28884 :
			found28881=False
		line=file.readline()
	found21302=False
	found27875=False
	for mutComb in mutCombsLower:
		mutCombList=mutComb.split("-")
		found=True
		for mut in mutCombList:
			variant=mut[1:]
			if not (variant in variants):
				found=False
				break
		if found:
			if mutComb=="t27875c-c27881t-g27882c-c27883t":
				found27875=True
			if mutComb=="c21302t-c21304a-g21305a":
				found21302=True
			if (mutComb!="c21304a-g21305a" or (not found21302)) and (mutComb!="c27881t-g27882c-c27883t" or (not found27875)):
				numSamplesComb[mutComb]+=1
	if found28881:
		numSamplesComb["t28881a-g28882a-g28883c"]+=1
	for var in variants:
		if var in numSamplesAll:
			numSamplesAll[var]+=1
		else:
			numSamplesAll[var]=1
file.close()

# if createFigures:
# 	import matplotlib.pyplot as plt
# 	legendSize=28
# 	tickSize=26
# 	labelSize=32
# 	#plot numbers of samples carrying variants
# 	fig = plt.figure(figsize=(15, 9))
# 	ax1 = fig.add_subplot(111)
# 	ax1.set_xscale('log')
# 	ax1.set_xlabel("Number of genomes per mutation", fontsize=labelSize)
# 	ax1.set_ylabel("Cumulative proportion", fontsize=labelSize)
# 	plt.yticks(fontsize=tickSize)
# 	plt.xticks(fontsize=tickSize)
# 	ax1.grid(False)
# 	ax1.spines[['right', 'top']].set_visible(False)
# 	for axis in ['top','bottom','left','right']:
# 		ax1.spines[axis].set_linewidth(2.5)
# 	ax1.tick_params(width=2.2)
# 	numSamplesForPlot=[]
# 	for var in numSamplesAll:
# 		numSamplesForPlot.append(float(numSamplesAll[var])/numMutations[var.upper()])
# 	ax1.ecdf(numSamplesForPlot)
# 	for mutComb in mutCombsLower:
# 		if mutComb=="g11083t-c21575t" or mutComb=="a507t-t508c-g509a":
# 			plt.axvline(x = float(numSamplesComb[mutComb])/focusedMultiNucMuts[mutComb.upper()], color = 'r', linewidth = 1, label = mutComb)
# 		else:
# 			plt.axvline(x = float(numSamplesComb[mutComb])/focusedMultiNucMuts[mutComb.upper()], color = 'g', linewidth = 1, label = mutComb)
# 	fig.savefig("cumulativeDistribution_numbersOfSamplesPerVariant.pdf",bbox_inches='tight')
# 	plt.close('all')

#create sorted list with most occurring mutation pairs
sortedListNumSamples=sorted(numSamplesAll.items(), key=lambda item: item[1])
print("Number of variants:")
print(len(sortedListNumSamples))
print("100 most abundant variants:")
print(sortedListNumSamples[-100:])
print("\n\n")
sortedListMutationPairs=sorted(numMutationPairs.items(), key=lambda item: item[1], reverse=True)
print("Number of mutation pairs:")
print(len(sortedListMutationPairs))
print("100 most frequent mutation pairs:")
#print(sortedListMutationPairs[-100:])
print(sortedListMutationPairs[:100])
print("\n\n")
sortedListMutations=sorted(numMutations.items(), key=lambda item: item[1])
print("Number of mutations:")
print(len(numMutations))
print("100 most frequent mutations:")
print(sortedListMutations[-100:])
print("\n\n")
numMutsListRef=[]
numMutsListNonRef=[]
proportionsSingletons=[]
genomesPerMutationList=[]
for mut in numMutations:
	if numMutations[mut]>40:
		if mut in numSingletons:
			proportionsSingletons.append(float(numSingletons[mut])/numMutations[mut])
		else:
			proportionsSingletons.append(0.0)
	pos=int(mut[1:-1])
	if ref[pos-1]==mut[0]:
		numMutsListRef.append(numMutations[mut])
		if numMutations[mut]>40:
			genomesPerMutationList.append(float(numSamplesAll[mut[1:].lower()])/numMutations[mut])
	else:
		numMutsListNonRef.append(numMutations[mut])
print(len(numMutsListRef))
print(len(numMutsListNonRef))
print("Average number of mutations of non-focused 1-nuc mutations from the reference: "+str(float(sum(numMutsListRef))/len(numMutsListRef))+" \t out of "+str(len(numMutsListRef))+" types observed ")
print("Average number of mutations of non-focused 1-nuc mutations from the reference: "+str(float(sum(numMutsListNonRef))/len(numMutsListNonRef))+" \t out of "+str(len(numMutsListNonRef))+" types observed ")
numUnder49=0
numOver49=0
for mut in range(len(numMutsListRef)):
	if numMutsListRef[mut]>49:
		numOver49+=1
	else:
		numUnder49+=1
print("Numbers of 1-nuc substitutions over 49 "+str(numOver49)+" under 49 "+str(numUnder49)+" proportion: "+str(float(numOver49/(numOver49+numUnder49))))


if createFigures:
	import matplotlib.pyplot as plt
	legendSize=28
	tickSize=26
	labelSize=32

	#number of mutations of each type
	fig = plt.figure(figsize=(15, 9))
	ax1 = fig.add_subplot(111)
	ax1.set_xscale('log')
	ax1.set_xlabel("Number of inferred substitutions", fontsize=labelSize)
	ax1.set_ylabel("Cumulative proportion", fontsize=labelSize)
	plt.yticks(fontsize=tickSize)
	plt.xticks(fontsize=tickSize)
	ax1.grid(False)
	ax1.spines[['right', 'top']].set_visible(False)
	for axis in ['top','bottom','left','right']:
		ax1.spines[axis].set_linewidth(2.5)
	ax1.set_ylim([0.0, 1.02])
	ax1.ecdf(numMutsListRef,linewidth = 4, color='black')
	for mutComb in mutCombs:
		if mutComb=="G11083T-C21575T" or mutComb=="A507T-T508C-G509A" or mutComb=="T21994C-T21995C" or mutComb=="T28881A-G28882A-G28883C":
			plt.axvline(x = len(focusedMultiNucMuts[mutComb]), color = 'r', linewidth = 1, label = mutComb)
		else:
			plt.axvline(x = len(focusedMultiNucMuts[mutComb]), color = (0.122,0.471,0.706), linewidth = 1, label = mutComb) #'g'
	ax1.tick_params(width=2.2,length=6.0)
	#ax1.tick_params(length=6.0)
	ax1.tick_params(which='minor',width=1.5,length=4.0)
	#ax1.tick_params(which='minor',length=4.5)
	fig.savefig("cumulativeDistribution_numOf1-nucFromRefMutations.pdf",bbox_inches='tight')
	plt.close('all')

	#number of genomes per mutation event
	fig = plt.figure(figsize=(15, 9))
	ax1 = fig.add_subplot(111)
	ax1.set_xscale('log')
	ax1.set_xlabel("Number of genomes per substitution", fontsize=labelSize)
	ax1.set_ylabel("Cumulative proportion", fontsize=labelSize)
	plt.yticks(fontsize=tickSize)
	plt.xticks(fontsize=tickSize)
	ax1.grid(False)
	ax1.spines[['right', 'top']].set_visible(False)
	for axis in ['top','bottom','left','right']:
		ax1.spines[axis].set_linewidth(2.5)
	#ax1.tick_params(width=2.2)
	ax1.tick_params(width=2.2,length=6.0)
	#ax1.tick_params(length=6.0)
	ax1.tick_params(which='minor',width=1.5,length=4.0)
	#ax1.tick_params(which='minor',length=4.5)
	ax1.set_ylim([0.0, 1.02])
	ax1.ecdf(genomesPerMutationList,linewidth = 4,color='black')
	for mutComb in mutCombs:
		if mutComb=="G11083T-C21575T" or mutComb=="A507T-T508C-G509A" or mutComb=="T21994C-T21995C" or mutComb=="T28881A-G28882A-G28883C":
			plt.axvline(x = float(numSamplesComb[mutComb.lower()])/len(focusedMultiNucMuts[mutComb]), color = 'r', linewidth = 1, label = mutComb)
		else:
			plt.axvline(x = float(numSamplesComb[mutComb.lower()])/len(focusedMultiNucMuts[mutComb]), color = (0.122,0.471,0.706), linewidth = 1, label = mutComb) #'g'
	fig.savefig("cumulativeDistribution_numbersOfSamplesPerVariant.pdf",bbox_inches='tight')
	plt.close('all')

#write statistics for combos and individual mutations within combos
proportionsSingletonsCombs={}
print("\n\n")
for mutComb in mutCombs:
	mutCombList=mutComb.split("-")
	counts=[0,0,0]#1desc, 2-10 desc, >10 desc
	for count in focusedMultiNucMuts[mutComb]:
		if count==1:
			counts[0]+=1
		elif count<11:
			counts[1]+=1
		else:
			counts[2]+=1
	proportionsSingletonsCombs[mutComb]=float(counts[0])/len(focusedMultiNucMuts[mutComb])
	print(mutComb+"\t\t"+str(len(focusedMultiNucMuts[mutComb]))+" mutations,   "+str(float(sum(focusedMultiNucMuts[mutComb]))/len(focusedMultiNucMuts[mutComb]))+" descendants per mut,    singletons: "+str(float(counts[0])/len(focusedMultiNucMuts[mutComb]))+",    2-10 desc: "+str(float(counts[1])/len(focusedMultiNucMuts[mutComb]))+",    >10 desc: "+str(float(counts[2])/len(focusedMultiNucMuts[mutComb]))+",    max desc: "+str(max(focusedMultiNucMuts[mutComb]))+",      num samples "+str(numSamplesComb[mutComb.lower()])+",      num samples per mutation "+str(float(numSamplesComb[mutComb.lower()])/len(focusedMultiNucMuts[mutComb]))+",      num inversions "+str(inversionsNums[mutComb]) )
	for mut in mutCombList:
		if mut in numMutations:
			print(mut+"\t\t"+str(numMutations[mut])+" mutations")
		else:
			print(mut+"\t\t"+"0 mutations")
	print("")

if createFigures:
	#plot proportions of singletons
	fig = plt.figure(figsize=(15, 9))
	ax1 = fig.add_subplot(111)
	#ax1.set_xscale('log')
	ax1.set_xlabel("Proportion of singletons", fontsize=labelSize)
	ax1.set_ylabel("Cumulative proportion", fontsize=labelSize)
	plt.yticks(fontsize=tickSize)
	plt.xticks(fontsize=tickSize)
	ax1.grid(False)
	ax1.spines[['right', 'top']].set_visible(False)
	for axis in ['top','bottom','left','right']:
		ax1.spines[axis].set_linewidth(2.5)
	#ax1.tick_params(width=2.2)
	ax1.set_ylim([0.0, 1.02])
	ax1.tick_params(width=2.2,length=6.0)
	ax1.ecdf(proportionsSingletons,linewidth = 4, color='black')
	for mutComb in mutCombs:
		if mutComb=="G11083T-C21575T" or mutComb=="A507T-T508C-G509A" or mutComb=="T21994C-T21995C" or mutComb=="T28881A-G28882A-G28883C":
			plt.axvline(x = proportionsSingletonsCombs[mutComb], color = 'r', linewidth = 1, label = mutComb)
		else:
			plt.axvline(x = proportionsSingletonsCombs[mutComb], color = (0.122,0.471,0.706), linewidth = 1, label = mutComb) #'g'
	fig.savefig("cumulativeDistribution_proportionOfSingletons.pdf",bbox_inches='tight')
	plt.close('all')

numMNMs=284+452+480+253+160+126+49+61+116+89+87+68+62+60+57
numMNMnuc=284*3+452*2+480*2+253*3+160*3+126*2+49*4+61*3+116*2+89*3+87*3+68*2+62*2+60*2+57*2
print(str(1819371)+" substitutions")
print(str(numMNMs)+" MNM substitutions")
print(str(numMNMnuc)+" nuc subs due to MNMs")
print(str(float(numMNMs)/(1819371+numMNMs))+" proportion of MNM substitutions")
print(str(float(numMNMnuc)/(1819371+numMNMnuc))+" proportion of nuc subs due to MNMs")

exit()





























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


#generate the string corresponding to a line of the tsv file for use in Taxonium.
def tsvForNode(tree,node,name,featureList,namesInTree,identicalTo=""):
	stringList=[name+"\t"]
	if identicalTo!="":
		stringList.append(identicalTo)
	stringList.append("\t")
	for feat in featureList:
		if node!=None:
			if hasattr(tree, feat):
				feature=getattr(tree, feat)
				if feat=="support" or feat=="IQsupport":
					stringList.append(str(feature[node]))
				#use this column to highlight which nodes could be placed (with probability above threshold) on the branch above the current node - used to highlight alternative placements of a given node on the tree.
				elif feat=="supportTo":
					for iNode in range(len(feature[node])):
						stringList.append(namesInTree[tree.name[feature[node][iNode][0]]]+":"+str(feature[node][iNode][1]))
						if iNode<(len(feature[node])-1):
							stringList.append(",")
				# elif feat=="alternativePlacements":
				# 	for iNode in range(len(feature[node])):
				# 		stringList.append(namesInTree[tree.name[feature[node][iNode][0]]]+":"+str(feature[node][iNode][1]))
				# 		if iNode<(len(feature[node])-1):
				# 			stringList.append(",")
				elif feat=="mutationsInf":
					if identicalTo=="":
						for iNode in range(len(feature[node])):
							mutation=feature[node][iNode]
							stringList.append(allelesListExt[mutation[0]]+str(mutation[1])+allelesListExt[mutation[2]]+":"+str(mutation[3]))
							if iNode<(len(feature[node])-1):
								stringList.append(",")
				elif feat=="Ns":
					if identicalTo=="":
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
				elif feat=="rootSupport":
					if feature[node]!=None:
						stringList.append(str(feature[node]))
			#TODO use this column to highlight nodes with support below threshold and with number of descendants above threshold
			elif feat=="supportGroup":
				if tree.support[node]<0.9:
					nDescString="nDesc<11_"
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


print("Number of final references in the MAT: "+str(numRefs[0]), flush=True)
print("Time spent finding placement nodes: "+str(timeFinding))
print("Time spent placing samples on the tree: "+str(timePlacing))
print("Time spent in total updating the topology and branch lengths: "+str(timeTopology))
print("Of which looking for placements for better topologies: "+str(totalTimeFindingParent[0]))

exit()














