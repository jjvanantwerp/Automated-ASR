# Semi-Automation of the Ancestral Sequence Reconstruction Workflow
# This is a program that constructs rough estimates of ancestral sequence reconstruction from an amino acid sequence.
# This program requires instalation of the Bioservieces module
# Written by James VanAntwerp in September 2021 - vanant25@msu.edu
# Written by Pattrick Finneran, Menten AI, Palo Alto, California, United States of America
# Written for the Woldring Lab, Michigan State University in East Lansing, Michigan, USA.

import sys
from Bio.Blast import NCBIWWW
import xml.etree.ElementTree as ET
import os
from string import ascii_letters
import matplotlib.pyplot as plt

Directory_Structure='''The program will make the following directory structure, where ROOT is the directory the script is placed in:
    ROOT/
    |--AutoASR.py                                   This script
    |--ASR/                                         The directory containing all work from the most recent run. The software will overwrite old runs.
    | |--HitsInfo.csv                               Data about the BlastP hits - ID and notes.
    | |--CD-Hit_Inital_Sequences.fasta.clstr        Residual from CD-Hit.
    | |--Sequence_Supplement/                       Files from adding additional sequence data to poorly supported regions of the tree. Not always used.
    | | |--*__Suppliment_HitsInfo.csv               Same as HitInfo.csv, but for each suppliment.
    | | |--*_BlastP_Results.fasta                   Each top 50 hits from sequences that are supplimented.
    | |--Final_Sequences.fasta                      The set of moderns equences used for the ASR and the foundation of all the later calculation
    | |--IQTree_Phylo/                              Directory for IQTree files from the phylogoney.
    | | |--Phylo.*                                  All IQTree files.
    | | |--Supports.*                               Data about the confidence of tree topology - .csv holds all data, .txt is summary.
    | | |--UFB_Confidences.png                      Histogram of UFB node support values.
    | | |--SHaLRT_Confidences.png                   Histogram of SHaLRT node support values
    | |--IQTree_ASR/                                Directory for IQTree files from ASR.
    | | |--ASR.*                                    All IQTree files.
    | | |--Confidences.*                            Data about the confidence of ASR prediction - .csv holds all data, .txt is summary, .png is a histogram.
    | |--IQTree_Binary/                             Directory for IQTree files from gap analysis.
    | | |--Binary_Alignment.fasta                   A binary alinment of modern sequences - gap or not.
    | | |--Binary.*                                 All IQTree files.
    | |--DNA_Library_12.5%_Cutoff.fasta.fasta/      Fasta of DNA for all high-confidence ancestors, with degenerate codons coding for all AAs with more than a 1/8 likelihood.
    | |--Concesus_Ancestors_with_Gaps.fasta         The set of likely ancestors, aligned and with gaps where the binary gap analysis predicts them.
    | |--Final_Tree.treefile                        The final tree which is best to read, as it has the most information in one place.'''

Help_String='''This software sutomates the process of ancestral sequence reconstruction. This software has one required input option and two optional input parameters.
     The first input in the name of the FASTA file containing the input sequence(s) or the input sequence itself.
     The second and third input are optional, but both must be entered for clarity. 
     The second input is the desired number of final sequences in the ASR. It must be a whole number greater than 20. The actual number of final sequences may differ depending on the sequence space. Default is 500 sequences.
     The third input is the sequence suppliment cutoff. This represents the similarity threshold below which singleton sequences will be supplimented to fill out sequence space. It must be a number between 0.4 and 0.9. Default is 0.7.
Options should be entered on the command line in the form
     python AutoASR.py {input.fasta } {Number of Desired Sequences} {Sequence Suppliment Cutoff}
Entering just the option "help" or "h" will display this description. '''

Software_Prerequistes='''This AutoASR script makes use of external software, which must be set up before it can be used. See the GitHub page for this project (github.com/jjvanantwerp/Automated-ASR) for download links.\
     This software requires IQTree 2, CD-Hit, and MAFFT to be downloaded and set up with their standard executable names added to the path. It also requires the python modules BioPython (Bio) and xml to be installed.\
    Optionally, the python module matplotlib can be installed for production of figures from the data.'''

# This dictionary provides the amino acid encoded for by every codon
Codon_to_AA = {
    'ata':'I', 'atc':'I', 'att':'I', 'atg':'M',
    'aca':'T', 'acc':'T', 'acg':'T', 'act':'T',
    'aac':'N', 'aat':'N', 'aaa':'K', 'aag':'K',
    'agc':'S', 'agt':'S', 'aga':'R', 'agg':'R',
    'cta':'L', 'ctc':'L', 'ctg':'L', 'ctt':'L',
    'cca':'P', 'ccc':'P', 'ccg':'P', 'cct':'P',
    'cac':'H', 'cat':'H', 'caa':'Q', 'cag':'Q',
    'cga':'R', 'cgc':'R', 'cgg':'R', 'cgt':'R',
    'gta':'V', 'gtc':'V', 'gtg':'V', 'gtt':'V',
    'gca':'A', 'gcc':'A', 'gcg':'A', 'gct':'A',
    'gac':'D', 'gat':'D', 'gaa':'E', 'gag':'E',
    'gga':'G', 'ggc':'G', 'ggg':'G', 'ggt':'G',
    'tca':'S', 'tcc':'S', 'tcg':'S', 'tct':'S',
    'ttc':'F', 'ttt':'F', 'tta':'L', 'ttg':'L',
    'tac':'Y', 'tat':'Y', 'tgc':'C', 'tgt':'C',
    'taa':' stop ', 'tag':' stop ', 'tga':' stop ', 'tgg':'W'
    }
# This dictionary provides the best codon for every amino acid in E Coli
AA_to_Codon_Ecoli = {
    'A':'gcc','R':'cgt','N':'aac','D':'gat','B':'aac',
    'C':'tgc','E':'gaa','Q':'cag','Z':'cag','G':'ggc',
    'H':'cat','I':'att','L':'ctg','K':'aaa','M':'atg',
    'F':'ttt','P':'ccg','S':'agc','T':'acc','W':'tgg',
    'Y':'tat','V':'gtg',' stop ':'taa'
    }
# This dictionary provides the best codon for every amino acid in humans
AA_to_Codon_Human = {
    'A':'gcc','R':'cgg','N':'aac','D':'gac','B':'aac',
    'C':'tgc','E':'gag','Q':'cag','Z':'cag','G':'ggc',
    'H':'cac','I':'atc','L':'ctg','K':'aag','M':'atg',
    'F':'ttc','P':'ccc','S':'agc','T':'acc','W':'tgg',
    'Y':'tac','V':'gtc',' stop ':'tga'
    }
# This dictionary provides the mixed-base code for every combination of 2-4 nucleic acids
Mixed_Bases_lookup = {
    'a':'a','c':'c','g':'g','t':'t',

    'ag':'r', 'ga':'r',
    'ct':'y', 'tc':'y',
    'ac':'m', 'ca':'m',
    'gt':'k', 'tg':'k',
    'gc':'s', 'cg':'s',
    'at':'w', 'ta':'w',

    'act':'h', 'atc':'h', 'cat':'h', 'cta':'h', 'tac':'h', 'tca':'h',
    'gct':'b', 'gtc':'b', 'ctg':'b', 'cgt':'b', 'tgc':'b', 'tcg':'b',
    'acg':'v', 'agc':'v', 'cag':'v', 'cga':'v', 'gca':'v', 'gac':'v',
    'agt':'d', 'atg':'d', 'gat':'d', 'gta':'d', 'tga':'d', 'tag':'d',

    'acgt':'n','actg':'n','agct':'n','agtc':'n','atcg':'n','atgc':'n',
    'cagt':'n','catg':'n','gact':'n','gatc':'n','tacg':'n','tagc':'n',
    'cgat':'n','ctag':'n','gcat':'n','gtac':'n','tcag':'n','tgac':'n',
    'cgta':'n','ctga':'n','gcta':'n','gtca':'n','tcga':'n','tgca':'n'
    }
# This dictionary provides the best degenerate codon for every combination of two amino acids for Humans
AA_Pair_lookup_Human = { 
    'AC':'ksc', 'AD':'gmc', 'AE':'gma', 'AF':'kyc', 'AG':'gsc', 'AH':'smc', 'AI':'ryc', 'AK':'rma', 'AL':'syc', 'AM':'ryg', 'AN':'rmc', 'AP':'scc', 'AQ':'sma', 'AR':'rsa', 'AS':'kcc', 'AT':'rcc', 'AV':'gyc', 'AW':'ksg', 'AY':'kmc', 
    'CA':'ksc', 'CD':'krc', 'CE':'krs', 'CF':'tkc', 'CG':'kgc', 'CH':'yrc', 'CI':'wkc', 'CK':'wrs', 'CL':'ykc', 'CM':'wks', 'CN':'wrc', 'CP':'ysc', 'CQ':'yrs', 'CR':'ygc', 'CS':'tsc', 'CT':'wsc', 'CV':'kkc', 'CW':'tgs', 'CY':'trc', 
    'DA':'gmc', 'DC':'krc', 'DE':'gas', 'DF':'kwc', 'DG':'grc', 'DH':'sac', 'DI':'rwc', 'DK':'ras', 'DL':'swc', 'DM':'rws', 'DN':'rac', 'DP':'smc', 'DQ':'sas', 'DR':'src', 'DS':'rrc', 'DT':'rmc', 'DV':'gwc', 'DW':'krs', 'DY':'kac', 
    'EA':'gma', 'EC':'krs', 'ED':'gas', 'EF':'kws', 'EG':'grg', 'EH':'sas', 'EI':'rwa', 'EK':'rag', 'EL':'swg', 'EM':'rwg', 'EN':'ras', 'EP':'smg', 'EQ':'sag', 'ER':'rrg', 'ES':'kmg', 'ET':'rmg', 'EV':'gwg', 'EW':'rrg', 'EY':'kas', 
    'FA':'kyc', 'FC':'tkc', 'FD':'kwc', 'FE':'kws', 'FG':'kkc', 'FH':'ywc', 'FI':'wtc', 'FK':'wwm', 'FL':'ytc', 'FM':'wts', 'FN':'wwc', 'FP':'yyc', 'FQ':'yws', 'FR':'ykc', 'FS':'tyc', 'FT':'wyc', 'FV':'ktc', 'FW':'tks', 'FY':'twc', 
    'GA':'gsc', 'GC':'kgc', 'GD':'grc', 'GE':'grg', 'GF':'kkc', 'GH':'src', 'GI':'rkc', 'GK':'rra', 'GL':'skc', 'GM':'rrg', 'GN':'rrc', 'GP':'ssc', 'GQ':'grg', 'GR':'sgg', 'GS':'rgc', 'GT':'rsc', 'GV':'gkc', 'GW':'kgg', 'GY':'krc', 
    'HA':'smc', 'HC':'yrc', 'HD':'sac', 'HE':'sas', 'HF':'ywc', 'HG':'src', 'HI':'mwc', 'HK':'mas', 'HL':'cwc', 'HM':'mws', 'HN':'mac', 'HP':'cmc', 'HQ':'cas', 'HR':'crc', 'HS':'mrc', 'HT':'mmc', 'HV':'swc', 'HW':'yrs', 'HY':'yac', 
    'IA':'ryc', 'IC':'wkc', 'ID':'rwc', 'IE':'rwa', 'IF':'wtc', 'IG':'rkc', 'IH':'mwc', 'IK':'awa', 'IL':'mtc', 'IM':'ats', 'IN':'awc', 'IP':'myc', 'IQ':'mya', 'IR':'aka', 'IS':'akc', 'IT':'ayc', 'IV':'rtc', 'IW':'wks', 'IY':'wwc', 
    'KA':'rma', 'KC':'wrs', 'KD':'ras', 'KE':'rag', 'KF':'wwm', 'KG':'rra', 'KH':'mas', 'KI':'awa', 'KL':'mwa', 'KM':'awg', 'KN':'aas', 'KP':'mma', 'KQ':'maa', 'KR':'arg', 'KS':'ars', 'KT':'ama', 'KV':'rwa', 'KW':'wrg', 'KY':'was', 
    'LA':'syc', 'LC':'ykc', 'LD':'swc', 'LE':'swg', 'LF':'ytc', 'LG':'skc', 'LH':'csc', 'LI':'mtc', 'LK':'mwa', 'LM':'mtg', 'LN':'mwc', 'LP':'cyc', 'LQ':'cyg', 'LR':'ckg', 'LS':'tyr', 'LT':'myc', 'LV':'stg', 'LW':'tkg', 'LY':'ywc', 
    'MA':'ryg', 'MC':'wks', 'MD':'rws', 'ME':'rwg', 'MF':'wts', 'MG':'rrg', 'MH':'mws', 'MI':'ats', 'MK':'awg', 'ML':'mtg', 'MN':'aws', 'MP':'myg', 'MQ':'mwg', 'MR':'akg', 'MS':'aks', 'MT':'ayg', 'MV':'rtg', 'MW':'wkg', 'MY':'wws', 
    'NA':'rmc', 'NC':'wrc', 'ND':'rac', 'NE':'ras', 'NF':'wwc', 'NG':'rrc', 'NH':'mac', 'NI':'awc', 'NK':'aas', 'NL':'mwc', 'NM':'aws', 'NP':'mmc', 'NQ':'mas', 'NR':'ars', 'NS':'arc', 'NT':'amc', 'NV':'rwc', 'NW':'wrs', 'NY':'wac', 
    'PA':'scc', 'PC':'ysc', 'PD':'smc', 'PE':'smg', 'PF':'yyc', 'PG':'ssc', 'PH':'cmc', 'PI':'myc', 'PK':'mma', 'PL':'cyc', 'PM':'myg', 'PN':'mmc', 'PQ':'cmr', 'PR':'csc', 'PS':'ycc', 'PT':'mcc', 'PV':'syc', 'PW':'ysg', 'PY':'ymc', 
    'QA':'sma', 'QC':'yrs', 'QD':'sas', 'QE':'sag', 'QF':'yws', 'QG':'grg', 'QH':'cas', 'QI':'mya', 'QK':'maa', 'QL':'cyg', 'QM':'mwg', 'QN':'mas', 'QP':'cmr', 'QR':'crg', 'QS':'yma', 'QT':'mma', 'QV':'swg', 'QW':'yrg', 'QY':'yas', 
    'RA':'rsa', 'RC':'ygc', 'RD':'src', 'RE':'rrg', 'RF':'ykc', 'RG':'sgg', 'RH':'crc', 'RI':'aka', 'RK':'arg', 'RL':'ckg', 'RM':'akg', 'RN':'ars', 'RP':'csc', 'RQ':'crg', 'RS':'ags', 'RT':'asc', 'RV':'rkc', 'RW':'ygg', 'RY':'yrc', 
    'SA':'kcc', 'SC':'tsc', 'SD':'rrc', 'SE':'kmg', 'SF':'tyc', 'SG':'rgc', 'SH':'mrc', 'SI':'akc', 'SK':'ars', 'SL':'tyr', 'SM':'aks', 'SN':'arc', 'SP':'ycc', 'SQ':'yma', 'SR':'ags', 'ST':'asc', 'SV':'kyc', 'SW':'tsg', 'SY':'tmc', 
    'TA':'rcc', 'TC':'wsc', 'TD':'rmc', 'TE':'rmg', 'TF':'wyc', 'TG':'rsc', 'TH':'mmc', 'TI':'ayc', 'TK':'ama', 'TL':'myc', 'TM':'ayg', 'TN':'amc', 'TP':'mcc', 'TQ':'mma', 'TR':'asc', 'TS':'asc', 'TV':'ryc', 'TW':'wsg', 'TY':'wmc', 
    'VA':'gyc', 'VC':'kkc', 'VD':'gwc', 'VE':'gwg', 'VF':'ktc', 'VG':'gkc', 'VH':'swc', 'VI':'rtc', 'VK':'rwa', 'VL':'stg', 'VM':'rtg', 'VN':'rwc', 'VP':'syc', 'VQ':'swg', 'VR':'rkc', 'VS':'kyc', 'VT':'ryc', 'VW':'kkg', 'VY':'kwc', 
    'WA':'ksg', 'WC':'tgs', 'WD':'krs', 'WE':'rrg', 'WF':'tks', 'WG':'kgg', 'WH':'yrs', 'WI':'wks', 'WK':'wrg', 'WL':'tkg', 'WM':'wkg', 'WN':'wrs', 'WP':'ysg', 'WQ':'yrg', 'WR':'ygg', 'WS':'tsg', 'WT':'wsg', 'WV':'kkg', 'WY':'trs', 
    'YA':'kmc', 'YC':'trc', 'YD':'kac', 'YE':'kas', 'YF':'twc', 'YG':'krc', 'YH':'yac', 'YI':'wwc', 'YK':'was', 'YL':'ywc', 'YM':'wws', 'YN':'wac', 'YP':'ymc', 'YQ':'yas', 'YR':'yrc', 'YS':'tmc', 'YT':'wmc', 'YV':'kwc', 'YW':'trs',
}
# This dictionary provides the best degenerate codon for every combination of two amino acids for EColi
AA_Pair_lookup_EColi = { 
    'AC':'ksc', 'AD':'gmc', 'AE':'gma', 'AF':'kyc', 'AG':'gsc', 'AH':'smc', 'AI':'ryc', 'AK':'rma', 'AL':'syc', 'AM':'ryg', 'AN':'rmc', 'AP':'scc', 'AQ':'sma', 'AR':'rsa', 'AS':'kcc', 'AT':'rcc', 'AV':'gyc', 'AW':'ksg', 'AY':'kmc', 
    'CA':'ksc', 'CD':'krc', 'CE':'krs', 'CF':'tkc', 'CG':'kgc', 'CH':'yrc', 'CI':'wkc', 'CK':'wrs', 'CL':'ykc', 'CM':'wks', 'CN':'wrc', 'CP':'ysc', 'CQ':'yrs', 'CR':'ygc', 'CS':'tsc', 'CT':'wsc', 'CV':'kkc', 'CW':'tgs', 'CY':'trc', 
    'DA':'gmc', 'DC':'krc', 'DE':'gas', 'DF':'kwc', 'DG':'grc', 'DH':'sac', 'DI':'rwc', 'DK':'ras', 'DL':'swc', 'DM':'rws', 'DN':'rac', 'DP':'smc', 'DQ':'sas', 'DR':'src', 'DS':'rrc', 'DT':'rmc', 'DV':'gwc', 'DW':'krs', 'DY':'kac', 
    'EA':'gma', 'EC':'krs', 'ED':'gas', 'EF':'kws', 'EG':'grg', 'EH':'sas', 'EI':'rwa', 'EK':'rag', 'EL':'swg', 'EM':'rwg', 'EN':'ras', 'EP':'smg', 'EQ':'sag', 'ER':'rrg', 'ES':'kmg', 'ET':'rmg', 'EV':'gwg', 'EW':'rrg', 'EY':'kas', 
    'FA':'kyc', 'FC':'tkc', 'FD':'kwc', 'FE':'kws', 'FG':'kkc', 'FH':'ywc', 'FI':'wtc', 'FK':'wwm', 'FL':'ytc', 'FM':'wts', 'FN':'wwc', 'FP':'yyc', 'FQ':'yws', 'FR':'ykc', 'FS':'tyc', 'FT':'wyc', 'FV':'ktc', 'FW':'tks', 'FY':'twc', 
    'GA':'gsc', 'GC':'kgc', 'GD':'grc', 'GE':'grg', 'GF':'kkc', 'GH':'src', 'GI':'rkc', 'GK':'rra', 'GL':'skc', 'GM':'rrg', 'GN':'rrc', 'GP':'ssc', 'GQ':'grg', 'GR':'sgg', 'GS':'rgc', 'GT':'rsc', 'GV':'gkc', 'GW':'kgg', 'GY':'krc', 
    'HA':'smc', 'HC':'yrc', 'HD':'sac', 'HE':'sas', 'HF':'ywc', 'HG':'src', 'HI':'mwc', 'HK':'mas', 'HL':'cwc', 'HM':'mws', 'HN':'mac', 'HP':'cmc', 'HQ':'cas', 'HR':'crc', 'HS':'mrc', 'HT':'mmc', 'HV':'swc', 'HW':'yrs', 'HY':'yac', 
    'IA':'ryc', 'IC':'wkc', 'ID':'rwc', 'IE':'rwa', 'IF':'wtc', 'IG':'rkc', 'IH':'mwc', 'IK':'awa', 'IL':'mtc', 'IM':'ats', 'IN':'awc', 'IP':'myc', 'IQ':'mya', 'IR':'aka', 'IS':'akc', 'IT':'ayc', 'IV':'rtc', 'IW':'wks', 'IY':'wwc', 
    'KA':'rma', 'KC':'wrs', 'KD':'ras', 'KE':'rag', 'KF':'wwm', 'KG':'rra', 'KH':'mas', 'KI':'awa', 'KL':'mwa', 'KM':'awg', 'KN':'aas', 'KP':'mma', 'KQ':'maa', 'KR':'arg', 'KS':'ars', 'KT':'ama', 'KV':'rwa', 'KW':'wrg', 'KY':'was', 
    'LA':'syc', 'LC':'ykc', 'LD':'swc', 'LE':'swg', 'LF':'ytc', 'LG':'skc', 'LH':'csc', 'LI':'mtc', 'LK':'mwa', 'LM':'mtg', 'LN':'mwc', 'LP':'cyc', 'LQ':'cyg', 'LR':'ckg', 'LS':'tyr', 'LT':'myc', 'LV':'stg', 'LW':'tkg', 'LY':'ywc', 
    'MA':'ryg', 'MC':'wks', 'MD':'rws', 'ME':'rwg', 'MF':'wts', 'MG':'rrg', 'MH':'mws', 'MI':'ats', 'MK':'awg', 'ML':'mtg', 'MN':'aws', 'MP':'myg', 'MQ':'mwg', 'MR':'akg', 'MS':'aks', 'MT':'ayg', 'MV':'rtg', 'MW':'wkg', 'MY':'wws', 
    'NA':'rmc', 'NC':'wrc', 'ND':'rac', 'NE':'ras', 'NF':'wwc', 'NG':'rrc', 'NH':'mac', 'NI':'awc', 'NK':'aas', 'NL':'mwc', 'NM':'aws', 'NP':'mmc', 'NQ':'mas', 'NR':'ars', 'NS':'arc', 'NT':'amc', 'NV':'rwc', 'NW':'wrs', 'NY':'wac', 
    'PA':'scc', 'PC':'ysc', 'PD':'smc', 'PE':'smg', 'PF':'yyc', 'PG':'ssc', 'PH':'cmc', 'PI':'myc', 'PK':'mma', 'PL':'cyc', 'PM':'myg', 'PN':'mmc', 'PQ':'cmr', 'PR':'csc', 'PS':'yct', 'PT':'mcc', 'PV':'syc', 'PW':'ysg', 'PY':'ymc', 
    'QA':'sma', 'QC':'yrs', 'QD':'sas', 'QE':'sag', 'QF':'yws', 'QG':'grg', 'QH':'cas', 'QI':'mya', 'QK':'maa', 'QL':'cyg', 'QM':'mwg', 'QN':'mas', 'QP':'cmr', 'QR':'crg', 'QS':'yma', 'QT':'mma', 'QV':'swg', 'QW':'yrg', 'QY':'yas', 
    'RA':'rsa', 'RC':'ygc', 'RD':'src', 'RE':'rrg', 'RF':'ykc', 'RG':'sgg', 'RH':'crc', 'RI':'aka', 'RK':'arg', 'RL':'ckg', 'RM':'akg', 'RN':'ars', 'RP':'csc', 'RQ':'crg', 'RS':'ags', 'RT':'asc', 'RV':'rkc', 'RW':'ygg', 'RY':'yrc', 
    'SA':'kcc', 'SC':'tsc', 'SD':'rrc', 'SE':'kmg', 'SF':'tyc', 'SG':'rgc', 'SH':'mrc', 'SI':'akc', 'SK':'ars', 'SL':'tyr', 'SM':'aks', 'SN':'arc', 'SP':'yct', 'SQ':'yma', 'SR':'ags', 'ST':'asc', 'SV':'kyc', 'SW':'tsg', 'SY':'tmc', 
    'TA':'rcc', 'TC':'wsc', 'TD':'rmc', 'TE':'rmg', 'TF':'wyc', 'TG':'rsc', 'TH':'mmc', 'TI':'ayc', 'TK':'ama', 'TL':'myc', 'TM':'ayg', 'TN':'amc', 'TP':'mcc', 'TQ':'mma', 'TR':'asc', 'TS':'asc', 'TV':'ryc', 'TW':'wsg', 'TY':'wmc', 
    'VA':'gyc', 'VC':'kkc', 'VD':'gwc', 'VE':'gwg', 'VF':'ktc', 'VG':'gkc', 'VH':'swc', 'VI':'rtc', 'VK':'rwa', 'VL':'stg', 'VM':'rtg', 'VN':'rwc', 'VP':'syc', 'VQ':'swg', 'VR':'rkc', 'VS':'kyc', 'VT':'ryc', 'VW':'kkg', 'VY':'kwc', 
    'WA':'ksg', 'WC':'tgs', 'WD':'krs', 'WE':'rrg', 'WF':'tks', 'WG':'kgg', 'WH':'yrs', 'WI':'wks', 'WK':'wrg', 'WL':'tkg', 'WM':'wkg', 'WN':'wrs', 'WP':'ysg', 'WQ':'yrg', 'WR':'ygg', 'WS':'tsg', 'WT':'wsg', 'WV':'kkg', 'WY':'trs', 
    'YA':'kmc', 'YC':'trc', 'YD':'kac', 'YE':'kas', 'YF':'twc', 'YG':'krc', 'YH':'yac', 'YI':'wwc', 'YK':'was', 'YL':'ywc', 'YM':'wws', 'YN':'wac', 'YP':'ymc', 'YQ':'yas', 'YR':'yrc', 'YS':'tmc', 'YT':'wmc', 'YV':'kwc', 'YW':'trs',
}
#An amino-acid key for IQ-Tree *.state files.
AA_key = ['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']

def fasta2dict(fasta_path,return_dict={}): ### Reads in fasta file and renames certain sequences based on forbidden characters in IQ Tree as needed
    # Read in the file and prepare some variables
    with open(fasta_path,'r') as infile:
        fastafile = infile.readlines()
    working_sequence = ''
    key = None
    for line in fastafile: # For every line in the input fasta
        if line[0] == '>': # Check if it's a name
            if working_sequence != '':  # If we already have a working sequence, another name indicates we're done. Otherwise record the name
                if len(working_sequence) >= 0:
                    return_dict[key] = working_sequence
            key = line[1:].rstrip().replace(':','_')
            key = key.replace('|','_')
            key = line[1:].rstrip()
            working_sequence = '' #clear the working sequence
        else: # If the line isn't a name, it's part of the sequence.
            if not all([char for char in working_sequence if (char.isalpha() or char=='-')]):
                raise ValueError(f"The provided file ({fasta_path}) was not a fasta file.")
            working_sequence = working_sequence + line.rstrip()
    if working_sequence != '': # Finally, clean up anything left in the working sequence before returning the dictionary.
        if not all([char for char in working_sequence if (char.isalpha() or char=='-')]):
            raise ValueError(f"The provided file ({fasta_path}) was not a fasta file.")
        else:
            if len(working_sequence) >= 0:
                return_dict[key] = working_sequence
    if len(return_dict)==0:
        raise ValueError(f"The provided file ({fasta_path}) was not a fasta file.")
    return return_dict

def dict2fasta(fasta_d,fpath): ### Saves dictionary to fasta where the dictionary key is the protein name and the value is the sequence
    with open(fpath,'w+') as out_fasta:
        for key,seq in fasta_d.items():
            out_fasta.write(f'>{key}\n{seq}\n')

def Is_Valid_AA(AA): #Is the argurment a valid amino acid or list of amino acids
    not_AAs = ['B','O','J','U','Z']
    if isinstance(AA,str):  # This block lets us evaluate strings of amino acids
        return(((AA in ascii_letters) or( AA == '-')) and (AA not in not_AAs))
    if isinstance(AA,list) and isinstance (AA[0],str): # This block lets us evaluate lists of amino acids
        return all(((i in ascii_letters) or( i == '-')) and (i not in not_AAs) for i in AA)
    if isinstance(AA,list) and isinstance (AA[0],list):  # This block lets us evaluate lists of lists of amino acids
        return all( all(((i in ascii_letters) or( i == '-')) and (i not in not_AAs) for i in lst) for lst in AA)
    else:
        raise ValueError ("A bad type was evaluated as an amino acid list....")

def Is_Valid_Codon(codon): #Is the argurment a valid codon
    dna_leters = ['a','c','t','g','r','y','m','k','s','w','h','b','v','d','n']
    if isinstance(codon,str):
        return (((len(codon)%3)==0) and all(i in dna_leters for i in codon))
    if isinstance(codon,list):
        return all(((len(c)==3) and all(i in dna_leters for i in c)) for c in codon)

def NCBI_to_XML(dirname,sequence,hits=2000,expect_value=0.30, seq_name=''): #Interact with BlastP, and record the XML
    # For the given sequence, we will run a BlastP search and parse the XML return to make a multi-fasta file
    blast_result_handle = NCBIWWW.qblast("blastp","nr", sequence,hitlist_size=hits,expect=expect_value)
    XML_Name=f"BlastP_XML_{seq_name}"
    with open(f"./{dirname}/{XML_Name}","w+") as fout:
        fout.write(blast_result_handle.read())
    blastp_xml = ET.parse(f"./{dirname}/{XML_Name}")
    return(blastp_xml)

def Parse_BlastP_XML(dirname,blastp_xml,sequence,sequence_name=None): #Parse the BlastP XML - record Fasta
    Fasta_Dict={}
    if sequence_name is None: #If no sequence_name is identified (Meaning this is the search with the user sequence)
        if isinstance(sequence,str):
            with open(f"./{dirname}/HitsInfo.csv","a+") as fout:
                fout.write(f"Hit ID,Hit Description,Hit Sequence\n")
                #Parsing the XML object, looking for hits
                for hit in blastp_xml.findall('./BlastOutput_iterations/Iteration/Iteration_hits/Hit'):
                    hitid = (hit.find('Hit_id')).text #I've tried to also add the Hit_accession, but I can't access that form the XML for some reason
                    hitdef=(hit.find('Hit_def')).text
                    hitaccession=(hit.find('Hit_accession')).text.replace(".","_")
                    seq = (hit.find('Hit_hsps/Hsp/Hsp_hseq')).text
                    #If the sequence doesn't have unknowns amino acids (annoying) then record it.
                    #The optional second method also removes exceptionally short or long sequences - be sure to synch with the code ~13 lines below
                    #if (("X" not in seq) and (len(seq)<((1+length_cutoff)*User_Sequence_Length)) and (len(seq)>((1-length_cutoff)*User_Sequence_Length))):
                    if ("X" not in seq):
                        fout.write(f"{hitid},{hitdef},{seq}\n")
                        Fasta_Dict[hitaccession]=seq
            with open(f"{dirname}/BlastP_Results.fasta","a+") as blastp_file:    
                blastp_file.write(f">User_Sequence\n{sequence}\n")
                for key,F_D_Sequence in Fasta_Dict.items():
                    #if (len(Sequence)<((1+length_cutoff)*User_Sequence_Length)) and (len(Sequence)>((1-length_cutoff)*User_Sequence_Length)):
                    blastp_file.write(f">{key}\n")
                    blastp_file.write(f"{F_D_Sequence.replace('-','')}\n")#We remove all gaps, because CD-Hit cannot handle gaps.
            Remove_Duplicate_Sequences_FASTA(dirname,"BlastP_Results.fasta") #Modify the BlastP return to remove duplicate sequences.
            return("BlastP_Results.fasta")
        elif isinstance(sequence,dict): #If we've been given multiple sequences, we have a dict for sequence and a list of xmls for blastp_xml
            for xml in blastp_xml:
                with open(f"./{dirname}/HitsInfo.csv","w+") as fout:
                    fout.write(f"Hit ID,Hit Description,Hit Sequence\n")
                    #Parsing the XML objects, looking for hits
                    for hit in xml.findall('./BlastOutput_iterations/Iteration/Iteration_hits/Hit'):
                        hitid = (hit.find('Hit_id')).text #I've tried to also add the Hit_accession, but I can't access that form the XML for some reason
                        hitdef=(hit.find('Hit_def')).text
                        hitaccession=(hit.find('Hit_accession')).text.replace(".","_")
                        seq = (hit.find('Hit_hsps/Hsp/Hsp_hseq')).text
                        #If the sequence doesn't have unknowns amino acids (annoying) then record it.
                        #The optional second method also removes exceptionally short or long sequences - be sure to synch with the code ~13 lines below
                        #if (("X" not in seq) and (len(seq)<((1+length_cutoff)*User_Sequence_Length)) and (len(seq)>((1-length_cutoff)*User_Sequence_Length))):
                        if ("X" not in seq):
                            fout.write(f"{hitid},{hitdef},{seq}\n")
                            Fasta_Dict[hitaccession]=seq
                with open(f"{dirname}/BlastP_Results.fasta","a+") as blastp_file:    
                    for key,F_D_Sequence in Fasta_Dict.items():
                        #if (len(Sequence)<((1+length_cutoff)*User_Sequence_Length)) and (len(Sequence)>((1-length_cutoff)*User_Sequence_Length)):
                        blastp_file.write(f">{key}\n")
                        blastp_file.write(f"{F_D_Sequence.replace('-','')}\n")#We remove all gaps, because CD-Hit cannot handle gaps.
                Remove_Duplicate_Sequences_FASTA(dirname,f"BlastP_Results.fasta") #Modify the BlastP return to remove duplicate sequences.
            return(f"BlastP_Results.fasta")
    elif isinstance(sequence_name,str): #If a sequence_name has been provided, this means we're doing Supplement searches so our output directory structure should be different.
        with open(f"{dirname}/{sequence_name}_Suppliment.csv","a+") as fout: #DIFFERENT FROM ABOVE
            fout.write(f"Hit ID,Hit Description,Hit Sequence\n")
            #Parsing the XML object, looking for hits
            for hit in blastp_xml.findall('./BlastOutput_iterations/Iteration/Iteration_hits/Hit'):
                hitid = (hit.find('Hit_id')).text #I've tried to also add the Hit_accession, but I can't access that form the XML for some reason
                hitdef=(hit.find('Hit_def')).text
                hitaccession=(hit.find('Hit_accession')).text.replace(".","_")
                seq = (hit.find('Hit_hsps/Hsp/Hsp_hseq')).text
                #If the sequence doesn't have unknowns amino acids (annoying) then record it.
                #The optional second method also removes exceptionally short or long sequences - be sure to synch with the code ~13 lines below
                #if (("X" not in seq) and (len(seq)<((1+length_cutoff)*User_Sequence_Length)) and (len(seq)>((1-length_cutoff)*User_Sequence_Length))):
                if ("X" not in seq):
                    fout.write(f"{hitid},{hitdef},{seq}\n")
                    Fasta_Dict[hitaccession]=seq
        with open(f"{dirname}/{sequence_name}_BlastP_Results.fasta","a+") as blastp_file: #DIFFERENT FROM ABOVE
            blastp_file.write(f">{sequence_name}\n{sequence}\n") #DIFFERENT FROM ABOVE
            for key,Sequence in Fasta_Dict.items():
                #if (len(Sequence)<((1+length_cutoff)*User_Sequence_Length)) and (len(Sequence)>((1-length_cutoff)*User_Sequence_Length)):
                blastp_file.write(f">{key}\n")
                blastp_file.write(f"{Sequence.replace('-','')}\n")#We remove all gaps, because CD-Hit cannot handle gaps.
        Remove_Duplicate_Sequences_FASTA(dirname,f"{sequence_name}_BlastP_Results.fasta") #Modify the BlastP return to remove duplicate sequences. #DIFFERENT FROM ABOVE
        return(f"{sequence_name}_BlastP_Results.fasta")
    else:
        print("sequence_name was not a string type")
        raise ValueError ("sequence_name was not a string type")        

def BlastP(dirname,sequence,hits=2000,expect_value=0.2,sequence_name=None): #This function takes an amino acid sequence, submitts a BlastP search, and records the result in a fasta file
    if not (os.path.isdir(dirname)):
        os.mkdir(dirname)
    if isinstance(sequence,str):
        sequence=sequence.replace('-','')
        sequence=sequence.replace('X','')
        if not all([char for char in sequence if (char.isalpha())]):
            print("Invalid sequence submitted for BlastP search")
            raise ValueError("Invalid sequence submitted for BlastP search")
    elif isinstance(sequence,dict):
        for key, item in sequence.items():
            item=item.replace('-','')
            item=item.replace('X','')
            if not all([char for char in item if (char.isalpha())]):
                print("Invalid sequence submitted for BlastP search ({key})")
                raise ValueError("Invalid sequence submitted for BlastP search ({key})")
    try:
        print("Acessing the NCBI database....")
        if isinstance(sequence,str):
            if not os.path.exists(f"{dirname}/BlastP_XML"):
                blastp_xml=NCBI_to_XML(dirname,sequence,hits,expect_value)
        elif isinstance(sequence,dict):
            xmls=[]
            for key,item in sequence.items():
                if not os.path.exists(f"{dirname}/BlastP_XML_{key}"):
                    xmls.append((NCBI_to_XML(dirname,item,hits,expect_value,key)))
    except:
        print("There was an error fetching the BlastP results")
        raise RuntimeError("There was an error fetching the BlastP results.")
    try:
        # Now, we parse the XML object and make a multi-fasta file. We also write the fasta file which is the result of our BlastP search.
        # The output behavior must be different if this the user sequence or a Supplement search though.
        if sequence_name is None and isinstance(sequence,str):  
            if not os.path.exists(f"{dirname}/BlastP_Results.fasta"):
                return_string = Parse_BlastP_XML(dirname,blastp_xml,sequence)
            return(f"BlastP_Results.fasta")
        elif sequence_name is None and isinstance(sequence,dict):  
            return_string = Parse_BlastP_XML(dirname,xmls,sequence)
            with open(f"{dirname}/{return_string}","a") as fout:
                for key,seq in sequence.items():
                    fout.write(f'>{key}\n{seq}\n')
            return(return_string)
        elif isinstance(sequence_name,str):
            if not os.path.exists(f"{dirname}/{sequence_name}_BlastP_Results.fasta"):
                return_string = Parse_BlastP_XML(dirname,blastp_xml,sequence,sequence_name)
            return(f"{sequence_name}_BlastP_Results.fasta")
        else:
            print("Variable sequence_name was not a string type")
            raise ValueError ("Variable sequence_name was not a string type")
    except:
        print("There was an error recording the BlastP Results")
        raise RuntimeError("There was an error recording the BlastP Results")     

def Sequence_Processing(dirname,finname,sequence):
    Fasta_Dict={}
    try:
        os.system(f"mafft {dirname}/{finname} > {dirname}/Early_Alignment_temp.fasta") #Alignment is needed for the Hamming Distance
    except:
        raise RuntimeError("There was an error running MAFFT.")
    Fasta_Dict=fasta2dict(f"{dirname}/Early_Alignment_temp.fasta",{})
    os.remove(f"{dirname}/Early_Alignment_temp.fasta")
    if isinstance(sequence,dict):
        user_seq_name=list(sequence.keys())[0]
        Hamming_Dict=Fasta_Dict_Hamming(Fasta_Dict,Fasta_Dict[user_seq_name])
    else:
        Hamming_Dict=Fasta_Dict_Hamming(Fasta_Dict,Fasta_Dict["User_Sequence"]) #We have now computed the hamming distance of all sequences.
    for key,item in Hamming_Dict.items():
        if (item/len(sequence))>0.5: #If a given sequence has less than 60% similarity with the user sequence, remove it.
            Fasta_Dict.pop(key)
    for key,item in Fasta_Dict.items():
        Fasta_Dict.update({key:item.replace("-","")}) #We now need to remove all the gaps in all the sequences to use CD-Hit
    Fasta_Dict.update(fasta2dict(f"{dirname}/{Supplement_Sequences(dirname,Fasta_Dict)}")) #This one line suppliments poorly supported areas of the tree, and updates the dictionary with those additional sequences.
    dict2fasta(Fasta_Dict,f"{dirname}/pre-Alignment.fasta") #Save the Supplimented Fasta_Dict
    try:
        os.system(f"mafft {dirname}/pre-Alignment.fasta > {dirname}/Raw_Alignment.fasta") #Align the supplimented sequences - called the Raw alignment
    except:
        raise RuntimeError("There was an error running MAFFT.")
    os.remove(f"{dirname}/pre-Alignment.fasta") #clean up our clutter, a little
    Fasta_Dict_Aligned=fasta2dict(f"{dirname}/Raw_Alignment.fasta") 
    if isinstance(sequence,dict):
        return_stirng = Post_MAFFT_processing(dirname,Fasta_Dict_Aligned,user_seq_name) #The Raw alignment can then be processed into an alignment that's ready for IQTree
    else:
        return_stirng = Post_MAFFT_processing(dirname,Fasta_Dict_Aligned,False) #if only one
    return return_stirng

def CDHit_Cutoff(identity):
    if (identity<=0.4 or identity>1 ):
        raise ValueError("The CD-Hit identity is invalid")
        print("The CD-Hit identity is invalid")
    elif identity>0.7:
        return(5)
    elif identity>0.6:
        return(4)
    elif identity>0.5:
        return(3)
    elif identity>0.4:
        return(2)

def Remove_Duplicate_Sequences_FASTA(dirname,fpath): #This function MODIFIES a fasta file to remove highly-duplicate sequences.
    try:
        os.system(f'cd-hit -i {dirname}/{fpath} -o {dirname}/{fpath[:-6]}_temp.fasta -c 0.99 -n 5') #Run CD-Hit with a 99% cutoff to remove duplicates/near duplicates
    except:
        print("There was an error running CD-Hit.")
        raise RuntimeError("There was an error running CD-Hit.")
    os.remove(f"{dirname}/{fpath}")
    os.remove(f"{dirname}/{fpath[:-6]}_temp.fasta.clstr")
    os.rename(f"{dirname}/{fpath[:-6]}_temp.fasta",f"{dirname}/{fpath}") #Replace the original file with one that has no duplicate sequences.

def Fasta_Dict_Hamming(Fasta_Dict,sequence): #returns the Hamming distance for every sequence in an aligned fasta dictionary.
    Hamming_dict={}
    for key,value in Fasta_Dict.items():
        Hamming_Diff=0
        if len(value)!=len(sequence):
            print("Hamming distance should be computed with properly aligned sequences.")
            raise ValueError("Hamming distance should be computed with properly aligned sequences.")
        for i,char in enumerate(value):
            if sequence[i]!=char:
                Hamming_Diff+=1
        Hamming_dict[key]=Hamming_Diff
    return Hamming_dict

def Supplement_Sequences(dirname,fasta_dict): #Find areas of poor coverage on the current tree and get additional BlastP search results.
    if not (os.path.isdir(f"{dirname}/Sequence_Supplement")): #Make a directory for this sequence Supplementing process
        os.mkdir(f"{dirname}/Sequence_Supplement")
    #Now we can begin to identify sequences which are similar enough to the rest of the dataset to be retained, but dissimilar enough that they would benefit from supplimenting.
    dict2fasta(fasta_dict,f"{dirname}/Sequence_Supplement/Sequences_to_be_Supplimented.fasta") #Be sure there are no gaps in this dictionary, as CD-Hit will reject those.
    sequences_list_for_search=[]   
    try:
        os.system(f'cd-hit -i {dirname}/Sequence_Supplement/Sequences_to_be_Supplimented.fasta -o {dirname}/Sequence_Supplement/SingleClusters_Supplemented.fasta -c {Suppliment_Cutoff} -n {CDHit_Cutoff(Suppliment_Cutoff)}') #Run CD-Hit with a 60% cutoff to identify 'loner' sequences
    except:
        raise RuntimeError("There was an error running CD-Hit.")
    with open(f"{dirname}/Sequence_Supplement/SingleClusters_Supplemented.fasta.clstr", 'r') as fin: #Looking at the CD-Hit clstr file (which holds cluster information)
        clusters=(fin.read()).split('>Cluster ')
    for cluster in clusters: #For each cluster
        seqs=cluster.split('>') #Split off the header from each sequence name
        if len(seqs)==2: #If there is only one sequence in the cluster
            sequences_list_for_search.append((seqs[1])[:-6]) #Seperate and record the name.
    #Now we can submit a BlastP search for all of the sequences that need to be supplimented, and write all sequences together as one file.
    files_list=[f"{dirname}/Sequence_Supplement/Sequences_to_be_Supplimented.fasta"]
    print(f"Supplementing {len(sequences_list_for_search)} sequences with poor neighboorhood coverage.")
    for sequence_name in sequences_list_for_search:
        if not os.path.exists(f"{dirname}/Sequence_Supplement/{sequence_name}_BlastP_Results.fasta"):
            BlastP(f"{dirname}/Sequence_Supplement",fasta_dict[sequence_name],50,0.01,sequence_name)  #Now we submit a BlastP search, but override the default expect and hits to get a more narrow set of sequences.
            files_list.append(f"{dirname}/Sequence_Supplement/{sequence_name}_BlastP_Results.fasta") #We also record a list of all the output fasta files to concatanate them together later.
    with open(f"{dirname}/Sequence_Supplement/Supplemented_BlastP_Sequences.fasta","a+") as fout:   #Now we need to write all of our Supplemented sequence searches together as one fasta file. Just smoosh 'em all together
        for fname in files_list:
            with open (fname,'r') as supfile:
                fout.write(supfile.read())
    Remove_Duplicate_Sequences_FASTA(dirname,"Sequence_Supplement/Supplemented_BlastP_Sequences.fasta")
    return(f"Sequence_Supplement/Supplemented_BlastP_Sequences.fasta")

def Trim_N_C_Termini (fasta_dict,Multiple_User_Sequences): #Trim termini based on User_Sequence
    if Multiple_User_Sequences:# Just get one of them
        User_Sequence= [n for n in fasta_dict.keys() if n[0:13] =="User_Sequence_"][0]
    else:
        User_Sequence=fasta_dict.get("User_Sequence")
    sequences_list = [i for i in fasta_dict.values()]
    len_to_remove_N=0
    len_to_remove_C=0
    for i,char in enumerate(User_Sequence): #Find N-terminus end of User_Sequence
        if char != '-':
            len_to_remove_N=i
            break
    for i in reversed(range(len(User_Sequence))): #Find C-terminus end of User_Sequence
        if User_Sequence[i] != '-':
            len_to_remove_C=i
            break
    for name,seq in fasta_dict.items():#Mutate Dictionary
        fasta_dict.update({name : (seq[len_to_remove_N:len_to_remove_C]) })

def Remove_Insertions (fasta_dict,User_Sequence,deletion_percentage=0.02,termini_length=0.05): #Remove sequnces that cause insertions in the alighment
    num_pos=len(User_Sequence)
    acceptable_num_gaps=round(len(User_Sequence.replace('-',''))*deletion_percentage)
    if acceptable_num_gaps<2:
        acceptable_num_gaps=2
    user_gap_pos=[i for i in range(num_pos) if User_Sequence.startswith('-', i)] #record all positions with a gap in the user sequence
    #The object 'key_sequence_list' contains a list of tuples with the name of the sequence, the sequence, and a list of all positions with an amino acid, in that order.
    key_sequence_gap_list = [(key,seq,[i for i in range(num_pos) if not seq.startswith('-', i)]) for key,seq in fasta_dict.items()] #Sorry for the mess of a line, but it prevents unessecary looping.
    keys_to_pop=[]
    #The ax is already at the root of the trees, and every tree that does not produce good fruit will be cut down and thrown into the fire
    for key,seq,no_gap_list in key_sequence_gap_list: #Now, for all sequences
        count=0
        for i,pos in enumerate(no_gap_list):#Evaluate the positions where it has amino acids
            if (pos in user_gap_pos)and(pos>num_pos*termini_length)and(pos<num_pos*(1-termini_length)): #If it has an amino acid where the User_Sequence has a gap, and isn't in the N- or C- terminus
                count+=1                    #Tally a strike
            else:
                count=0                     #If we come to a place of agreement, reset the count.
            if count > acceptable_num_gaps: #If the number of insertions in a row (count) is higher than what is acceptable
                keys_to_pop.append(key)     #Record this key as one to remove
                break                       #What further testimony do we need? We have heard it ourselves from his own lips.
        else:
            continue
    for key in keys_to_pop:
        fasta_dict.pop(key)

def Remove_Deletions (fasta_dict,User_Sequence,deletion_percentage=0.02,termini_length=0.05): #Remove sequnces that cause insertions in the alighment
    num_pos=len(User_Sequence)
    acceptable_num_gaps=round(len(User_Sequence.replace('-',''))*deletion_percentage)
    if acceptable_num_gaps<2:
        acceptable_num_gaps=2
    user_gap_pos=[i for i in range(num_pos) if User_Sequence.startswith('-', i)] #record all positions with a gap
    #The object 'key_sequence_list' contains a list of tuples with the name of the sequence, the sequence, and a list of all positions with a gap, in that order.
    key_sequence_gap_list = [(key,seq,[i for i in range(num_pos) if seq.startswith('-', i)]) for key,seq in fasta_dict.items()] #Sorry for the mess of a line, but it prevents unessecary looping.
    keys_to_pop=[]
    #The ax is already at the root of the trees, and every tree that does not produce good fruit will be cut down and thrown into the fire
    for key,seq,gap_list in key_sequence_gap_list: #Now, for all sequences
        count=0
        for i,pos in enumerate(gap_list):   #Evaluate the positions where it has gaps
            if (pos not in user_gap_pos)and(pos>num_pos*termini_length)and(pos<num_pos*(1-termini_length)): #If it has gaps where the User_Sequence does not, and isn't in the N- or C- terminus
                count+=1                    #Tally a strike
            else:
                count=0                     #If we come to a place of agreement, reset the count.
            if count > acceptable_num_gaps: #If the number of gaps in a row (count) is higher than what is acceptable
                keys_to_pop.append(key)     #Record this key as one to remove
                break                       #What further testimony do we need? We have heard it ourselves from his own lips.
        else:
            continue
    for key in keys_to_pop:
        fasta_dict.pop(key)    

def Clean_all_gaps (fasta_dict):
    #For each position, go thorugh all sequences. If any sequence has a residue, leave the position in at the end.
    sequences_list = [i for i in fasta_dict.values()]
    pos_to_leave=[]
    for pos in range(len(sequences_list[0])):
        for seq in sequences_list:
            if (seq[pos] != '-'):
                pos_to_leave.append(pos)
                break
        else:
            continue
    #Remove all positions of all gaps from all sequences in the alignment
    for key,sequence in fasta_dict.items(): #Remove all the gaps
        fasta_dict.update({key:''.join([ char for i,char in enumerate(sequence) if i in pos_to_leave ])})

def Post_MAFFT_processing(dirname,fasta_dict,multisequence,dynamic_sequence_reduction=True): #Modifications after the alignment, mostly having to do with gaps.
    # These functions **MODIFY** the fasta_dict by doing what their names say they do.
    old_user_length=0
    if multisequence:# Just get one of them
        User_Sequence=fasta_dict.get(multisequence)
        while len(User_Sequence)!=old_user_length: #Each of these functions mutate the alignment, so it's important to repeat until we've finsihed finding all misalignments.
            old_user_length=len(User_Sequence)#Store length
            Remove_Insertions(fasta_dict,User_Sequence)#Insertions
            Clean_all_gaps(fasta_dict)#Clean
            User_Sequence=fasta_dict.get(multisequence)#Update
            Remove_Deletions(fasta_dict,User_Sequence)#Deletions
            for key, seq in fasta_dict.items():#Now we realign the sequences rather than just removing all the gaps - more relaiable
                fasta_dict.update({key:seq.replace("-",'')})
            dict2fasta(fasta_dict,f"{dirname}/Post_Mafft_Working1.fasta")
            try:
                os.system(f"mafft {dirname}/Post_Mafft_Working1.fasta > {dirname}/Post_Mafft_Working2.fasta") #Alignment is needed for the Hamming Distance
                fasta_dict={}            
            except:
                raise RuntimeError("There was an error running MAFFT.")
            fasta2dict(f"{dirname}/Post_Mafft_Working2.fasta",fasta_dict)
            os.system(f"rm {dirname}/Post_Mafft_Working1.fasta {dirname}/Post_Mafft_Working2.fasta")
            User_Sequence=fasta_dict.get(multisequence)#Update
    else:
        User_Sequence=fasta_dict.get("User_Sequence")
        while len(User_Sequence)!=old_user_length: #Each of these functions mutate the alignment, so it's important to repeat until we've finsihed finding all misalignments.
            old_user_length=len(User_Sequence)#Store length
            Remove_Insertions(fasta_dict,User_Sequence)#Insertions
            Clean_all_gaps(fasta_dict)#Clean
            User_Sequence=fasta_dict.get("User_Sequence")#Update
            Remove_Deletions(fasta_dict,User_Sequence)#Deletions
            for key, seq in fasta_dict.items(): #Now we realign the sequences rather than just removing all the gaps - more relaiable
                fasta_dict.update({key:seq.replace("-",'')})
            dict2fasta(fasta_dict,f"{dirname}/Post_Mafft_Working1.fasta")
            try:
                os.system(f"mafft {dirname}/Post_Mafft_Working1.fasta > {dirname}/Post_Mafft_Working2.fasta") #Alignment is needed for the Hamming Distance
                fasta_dict={}            
            except:
                raise RuntimeError("There was an error running MAFFT.")
            fasta2dict(f"{dirname}/Post_Mafft_Working2.fasta",fasta_dict)
            os.system(f"rm {dirname}/Post_Mafft_Working1.fasta {dirname}/Post_Mafft_Working2.fasta")
            User_Sequence=fasta_dict.get("User_Sequence")#Update
    dict2fasta(fasta_dict,f'{dirname}/Post_MAFFT_Cleaned_Penultimate.fasta')
    if dynamic_sequence_reduction and len(fasta_dict)>=Final_Library_Size : #If on, this code will reduce the sequence alignment down to less than 500 sequences. It's on by default
        keep_going=True
        identity=0.97
        n=5
        for key,item in fasta_dict.items():
            fasta_dict.update({key:item.replace("-","")})
        dict2fasta(fasta_dict,f'{dirname}/Post_MAFFT_No_Gaps.fasta')
        while keep_going: #Otherwise, we'll keep lowering the CD-Hit cutoff until we get below 500 sequences.
            try:
                os.system(f"rm {dirname}/Penultimate_Sequences.fasta")
            except:
                continue
            try:
                os.system(f'cd-hit -i {dirname}/Post_MAFFT_No_Gaps.fasta -o {dirname}/Penultimate_Sequences.fasta -c {round(identity,3)} -n {CDHit_Cutoff(identity)}') #Run CD-Hit
            except:
                print("There was an error running CD-Hit.")
                raise RuntimeError("There was an error running CD-Hit.")
            identity-=0.01
            with open(f"{dirname}/Penultimate_Sequences.fasta","r") as fin:
                lines=fin.readlines()
            if len(lines)<=1000 or identity<0.80:
                print(f"Number Reached. Lines={len(lines)}")
                keep_going=False
        try:
            os.system(f"rm {dirname}/Post_MAFFT_No_Gaps.fasta")
            os.system(f'mafft {dirname}/Penultimate_Sequences.fasta > {dirname}/Final_Sequences.fasta') #align the sequences passed into the function
        except:
            print("There was an error creating the sequence alignemnt.")
            raise RuntimeError("There was an error creating the sequence alignemnt.")
        return("Final_Sequences.fasta")
    else:
        os.system(f'cp {dirname}/Post_MAFFT_Cleaned_Penultimate.fasta {dirname}/Final_Sequences.fasta')
        os.system(f'rm {dirname}/Post_MAFFT_Cleaned_Penultimate.fasta')
        return("Final_Sequences.fasta")

def IQTree_Phylo(dirname,finname): #Construct a phylogenetic tree, and determine a good model.
    if not os.path.isdir(f"{dirname}/IQTree_Phylo"): #Make the directory
        os.mkdir(f"{dirname}/IQTree_Phylo")
    try:
        #IQTree section
        os.system(f'iqtree2 -s {dirname}/{finname} -pre {dirname}/IQTree_Phylo/Phylo -alrt 20000 -bb 20000 -nt AUTO')
    except:
        raise RuntimeError("There was an error building the phylogeny in IQTree.")
        print("There was an error building the phylogeny in IQTree.")

def IQTree_ASR(dirname,finname): #Using the tree/model from the phylogney, do the ASR
    model=''
    try: #this can fail because for some reason HPCC tries to run it before the IQTree_Phylo has finished running.
        with open(f"{dirname}/IQTree_Phylo/Phylo.iqtree","r") as fin: #we need to go to the .iqtree file from the phylogony and find the best model for this data set. This will save a ton of time.
            iqtree_file_lines=fin.readlines()
        if "Best-fit model" in iqtree_file_lines[42]:
            model = iqtree_file_lines[42].split(':')[1].strip()
        else:   #If we can't find the Best-fit model where we expect it to be, then check the whole file.
            for line in iqtree_file_lines:
                if "Best-fit model" in line:
                    model = line.split(':')[1].strip()
    except:
        pass
    if not os.path.isdir(f"{dirname}/IQTree_ASR"): #Make the directory
        os.mkdir(f"{dirname}/IQTree_ASR")
    if model == '': #Preventing errors; if we can't find the model, then we'll just use the default behavior.
        model ='MFP'
    try:
        #IQTree section
        #Ask what model is best for ASR, use that one here.
        os.system(f'iqtree2 -s {dirname}/{finname} -te {dirname}/IQTree_Phylo/Phylo.treefile -pre {dirname}/IQTree_ASR/ASR -asr -m {model} -redo -nt AUTO')
    except:
        raise RuntimeError("There was an error conducting ASR in IQTree.")
        print("There was an error conducting ASR in IQTree.")
    with open(f"{dirname}/IQTree_Phylo/Phylo.treefile","r") as fin:
        Phylo_Tree=fin.read()
    with open(f"{dirname}/IQTree_ASR/ASR.treefile","r") as fin:
        ASR_Tree=fin.read()
    with open(f"{dirname}/Final_Tree.treefile","w+") as fout:
        fout.write(f"{Phylo_Tree}{ASR_Tree}")
    ASR_Statefile_Dict = Statefile_to_Dict(dirname,"IQTree_ASR/ASR.state") #Returns a dictionary out of *.state file
    return(ASR_Statefile_Dict)

def Select_Ancestor_Nodes(dirname): #Select ancestral nodes which have a high enough confidence to resurect ancestors from them.
    Supports={}
    UFB_Supports=[]
    SHALRT_Supports=[]
    #This will evaluate the tree, which is in Newick format, and select nodes of high enough confidence to reconstruct an ancestral library.
    #The Phlyo.treefile has confidence values as (SH-aLRT/UFB), and the ASR.treefile has node names.
    with open(f'{dirname}/IQTree_ASR/ASR.treefile', 'r') as ASRtreefin:
        ASRtreefile = ASRtreefin.read()
    with open(f'{dirname}/IQTree_Phylo/Phylo.treefile', 'r') as Phylotreefin:
        Phylotreefile = Phylotreefin.read()
    confident_nodes=[] #these will be the nodes with good values. This isn't quite the same value between *.contree and *.trefile, but close
    ASRnodes=[n.split(',')[0].split(':') for n in ASRtreefile.split(')')]  #Let's split off just the information for each node, which is stored after every close parenthiesis.
    Phylonodes=[n.split(',')[0].split(':') for n in Phylotreefile.split(')')]
    ASRnodes.pop(0) #the first split is not actually a node.
    Phylonodes.pop(0)#the first split is not actually a node.
    for i in range(len(ASRnodes)-1): #We exclude the last node, because that's the root, which will mess up the below lines.
        UFB_i=float(Phylonodes[i][0].split('/')[1])
        SHALRT_i=float(Phylonodes[i][0].split('/')[0])
        Supports[ASRnodes[i][0]]=((UFB_i,SHALRT_i))
        UFB_Supports.append(UFB_i)
        SHALRT_Supports.append(SHALRT_i)
        if (SHALRT_i > 80) and (UFB_i > 95): #If the SH-aLRT >80% and the ultrafast bootstraping is >95%
            confident_nodes.append(ASRnodes[i][0]) #Record the name of the high-confidence nodes.
    SHALRT_mean = sum(SHALRT_Supports)/len(SHALRT_Supports)
    UFB_mean = sum(UFB_Supports)/len(UFB_Supports)
    with open(f"{dirname}/IQTree_Phylo/Supports.txt",'w+') as fout:
        fout.write(f"The phylogenetic tree has a mean ultra-fast bootstrap support of {round(UFB_mean)}%. For this method, 95% is considered the threshold of quality.\n")
        fout.write(f"{len([n for n in UFB_Supports if n > 95])} out of {len(UFB_Supports)} nodes have an ultra-fast bootstrap support above 95%.\n")
        fout.write(f"The standard deviation of ultra-fast bootstrap support is {round(((1/len(UFB_Supports))*sum([(x-UFB_mean)**2 for x in UFB_Supports]))**0.50,1)}%.\n\n")
        fout.write(f"The phylogenetic tree has a mean SH-aLRT support of {round(SHALRT_mean)}%. For this method, 80% is considered the threshold of quality.\n")
        fout.write(f"{len([n for n in SHALRT_Supports if n > 80])} out of {len(SHALRT_Supports)} nodes have a SH-aLRT support above 80%.\n")
        fout.write(f"The standard deviation of SH-aLRT support is {round(((1/len(SHALRT_Supports))*sum([(x-SHALRT_mean)**2 for x in SHALRT_Supports]))**0.50,1)}%.\n")
    UFB=[]
    SHALRT=[]
    with open(f"{dirname}/IQTree_Phylo/Supports.csv",'w+') as fout:
        fout.write("Node,Ultra-Fast Bootstrap,SH-aLRT\n")
        for i in range(len(ASRnodes)-1):
            fout.write(f"{ASRnodes[i][0]},{Supports[ASRnodes[i][0]][0]},{Supports[ASRnodes[i][0]][1]}\n")
            UFB.append(Supports[ASRnodes[i][0]][0])
            SHALRT.append(Supports[ASRnodes[i][0]][1])
    try:
        # Make the sequence prediction histogram
        plt.rcParams["figure.figsize"] = [7.00, 5.00]
        plt.rcParams["figure.autolayout"] = True
        n_bins=101
        # Plot the histogram
        plt.hist(UFB,n_bins)
        plt.title('Ultra-Fast Bootstrap Node Confidence Values')
        plt.xlabel("Confidence (%)")
        # Save the histogram
        plt.savefig(f"{dirname}/IQTree_Phylo/UFB_Confidences.png")
        plt.clf()
    except:
        pass
    try:
        # Make the sequence prediction histogram
        plt.rcParams["figure.figsize"] = [7.00, 5.00]
        plt.rcParams["figure.autolayout"] = True
        n_bins=101
        # Plot the histogram
        plt.hist(SHALRT,n_bins)
        plt.title('SH-aLRT Node Confidence Values')
        plt.xlabel("Confidence (%)")
        # Save the histogram
        plt.savefig(f"{dirname}/IQTree_Phylo/SHaLRT_Confidences.png")
        plt.clf()
    except:
        pass
    return confident_nodes #return the list of node names whose value is above the cutoff.
    
def Statefile_to_Dict(dirname,fname): #Parse the statefile into a dictionary and record the confidence values.
    statefile_dict={} #This is a dicitonary made out of the statefile - its keys are the node names, 
    #and its values are a list of tuples with the amino acid and list of amino acid distributions at each position.
    node_distribution=[]
    statelines=[]
    #This will parse the .state file for the desired nodes to get their AA distribution
    with open(f'{dirname}/{fname}', 'r') as statefin: #Read in each line, skipping the header.
        for i,line in enumerate(statefin):
            if i>8:
                statelines.append(line)
    #Now let's pull the data from each line
    working_node=statelines[0].split()[0] #prime the working node
    for line in statelines: # For every line in the state file
        line_list = line.split() #Break up the line into columns
        if working_node == line_list[0]: #If we're still working on the same node
            int_line_list = list(map(float, line_list[3:]))
            node_distribution.append(int_line_list) #record the probability distribution and the maximum probability of the position
        else: #If we've come to the end of a node,
            statefile_dict[working_node]=node_distribution #Add a key-value pair to the statefile dictionary that is the node name and the node's list
            working_node = line_list[0] #update the working_node value
            node_distribution=[] #Clear the working node lists
            int_line_list = list(map(float, line_list[3:]))
            node_distribution.append(int_line_list) #add to the node list the probabbilty distribution and the maximum probability at that position
    statefile_dict[working_node]=node_distribution #Be sure to add the last node into the dictionaries too!
    return(statefile_dict)    #Dictionary format looks like {NodeX:[[probs],[probs],[probs],...]}

def Binary_Gap_Analysis(dirname,finname): #Do a binary ASR to determine where gaps in the ancestral sequences reside. Returns the Binary_Statefile_Dict
    if not os.path.isdir(f"{dirname}/IQTree_Binary"): #Make the directory
        os.mkdir(f"{dirname}/IQTree_Binary")
    binary_dict = fasta2dict(f"{dirname}/{finname}")
    #These lines of code are to fix the bizare bug that has been showoing up with the dynamic sequence reduction that I spent all of spring 2022 trying to fix.
    #The bug is that binary_dict somehow imported all of the non-alighed sequences from post_mafft_cleaned.fasta along with the sequences from Final_Sequences.fasta, but not in a redundant way?
    alignment_length=max([len(seq)] for seq in binary_dict.values())
    to_pop=[]
    for key,seq in binary_dict.items():
        if len(seq)!=alignment_length[0]:
            to_pop.append(key)
    for key in to_pop:
        binary_dict.pop(key)
    for key,seq in binary_dict.items(): #Make a binary alignment of the input fasta
        binary_dict.update({key:(''.join(['0' if aa == '-' else '1' for aa in seq]))}) #This line is taken from Ben - how to give proper credit? 
        #Note that a 0 is a gap and a 1 is an AA
    dict2fasta(binary_dict, f"{dirname}/IQTree_Binary/Binary_Alignment.fasta") #Write to file for IQTree
    try:
        os.system(f'iqtree2 -s {dirname}/IQTree_Binary/Binary_Alignment.fasta -te {dirname}/IQTree_Phylo/Phylo.treefile -pre {dirname}/IQTree_Binary/Binary -blfix -asr -m GTR2+FO -redo -nt AUTO')
    except:
        print("There was an error determining gaps in the ancestral sequence")
        raise RuntimeError("There was an error determining gaps in the ancestral sequence")
    try:    
        Binary_Statefile_Dict = Statefile_to_Dict(dirname,"IQTree_Binary/Binary.state")
    except FileNotFoundError:
        print("The binary gap analysis was not successful.")
    ASR_Statefile_Dict = Statefile_to_Dict(dirname,"IQTree_ASR/ASR.state")
    Consensus_Ancestors_with_Gaps={}
    Pos_with_Gaps={} #dictionary of {NodeX:[list of gaps at NodeX]}
    #Find positions that are actually gaps in the ASR
    for node,item in Binary_Statefile_Dict.items():
        gap_pos=[]#list of positions with gaps
        for i,pos in enumerate(item): #each position
            if float(pos[0])>0.5:#If the posotion has majority gap
            #The reason I've written it this was is that when the chances are  close together (0.501 to 0.499) IQTree puts a gap in the binary gap analysis *facepalm*
                gap_pos.append(i)
        Pos_with_Gaps[node]=gap_pos #At each node, record the positions in the ancestral sequence that has majority node.
    #Merge the Sequence ASR with the gap ASR
    for node,cons_list in ASR_Statefile_Dict.items():#node is name, cons_list is list of (list of AA confidence values) for each position
        consensus_seq=''
        for i,pos in enumerate(cons_list):#pos is the position of ancestor at node
            if i in (Pos_with_Gaps[node]): #If this position at this node is likely a gap, add a gap to the consensus sequence
                consensus_seq+='-'
            else:
                consensus_seq+= AA_key[pos.index(max(pos))] #Otherwise, add the amino acid from ASR
        Consensus_Ancestors_with_Gaps[node]=consensus_seq
    Clean_all_gaps(Consensus_Ancestors_with_Gaps)
    dict2fasta(Consensus_Ancestors_with_Gaps,f"{dirname}/Consensus_Ancestors_with_Gaps.fasta")
    Consensus_Ancestors_without_Gaps={}
    for key,item in Consensus_Ancestors_with_Gaps.items():
        Consensus_Ancestors_without_Gaps.update({key:item.replace("-","")})
    dict2fasta(Consensus_Ancestors_without_Gaps,f"{dirname}/Consensus_Ancestors_without_Gaps.fasta")
    return(Binary_Statefile_Dict)

def Degenerate_Nucleotide_Codon(AA_List, source='EColi'):  #This function takes a list of AAs for ONE POSITION and makes a degenerate codon for them.
    # This is an area where the library size could be significantly improved upon - take a more detailed look here - high priority
    if len(AA_List)==0:
        return ''
    while ('X' in AA_List):
        AA_List.remove('X')
        # Unknown amino acids can occur, but we're going to ignore them.
    # This function takes a list of AAs for ONE POSITION and makes a degenerate codon for them.
    for AA in AA_List:
        if not Is_Valid_AA(AA):
            raise ValueError ("That's not an amino acid abreviation....")
    # Codon-based Lookup
    codon_list=[] #A list of all the codons that need to be coded for (eg, ['atc','agc','gta'])
    if (source=='EColi'):
        for AA in AA_List:
            codon_list.append(AA_to_Codon_Ecoli[AA])
    elif (source=='Human'):
        for AA in AA_List:
            codon_list.append(AA_to_Codon_Human[AA])
    else:
        raise NameError("Please Specify EColi or Human as an expression organism.")
    degenerate_codon = '' #The three-charecter degenerate codon made from the codons in codon_list
    for pos in range(3): #For every position in the codon
        bases_at_pos='' #bases_at_pos will be a list of all the nucleotide bases wanted at this position of the codon.
        for codon in codon_list: #For every codon
            if (codon[pos] not in bases_at_pos): #If that codon's base is not in the bases_at_pos already,
                bases_at_pos+=codon[pos] #Add it
        degenerate_codon+=Mixed_Bases_lookup[bases_at_pos] #What degenerate codon codes for the desired bases?
    return degenerate_codon

def Build_DNA_Sequence (Primer_Request, source='EColi'):
    primer_txt = ''
    if source == 'EColi':
        for pos in Primer_Request: # pos is the list of AA at the respective position
            if (len(pos)==1):
                if (pos[0]!='-'): # For a gap in the sequence, we will record nothing.
                    primer_txt+=(AA_to_Codon_Ecoli[pos[0]])
            elif (len(pos)==2): # For a request of two AAs we use a lookup table
                primer_txt+=(AA_Pair_lookup_EColi[ pos[0]+pos[1] ])
            else: #Otherwise we bruteforce the degenerate DNA
                primer_txt+=(Degenerate_Nucleotide_Codon(pos))
    elif source == 'Human': #Same as above
        for pos in Primer_Request: 
            if (len(pos)==1):
                if (pos[0]!='-'): 
                    primer_txt+=(AA_to_Codon_Human[pos[0]])
            elif (len(pos)==2):
                primer_txt+=(AA_Pair_lookup_Human[ pos[0]+pos[1] ])
            else:
                primer_txt+=(Degenerate_Nucleotide_Codon(pos))
    else:
        raise NameError ("The requested organism is not specified - please specify Human or EColi")
    if not Is_Valid_Codon(primer_txt):
        raise ValueError ("The primer generated was not a valid DNA sequence. No clue how that happened. If you're seeing this error, it's proabbly caused by a bug.")
    return (primer_txt)

def Make_Uncertianty_Libraries (dirname,ASR_Statefile_Dict,Binary_Statefile_Dict,Cutoff=0.125): #Make Uncertianty Libraries for all ancestral sequences with a high enough confidence at thier node
    if not (Cutoff>0) and (Cutoff<1):
        raise ValueError("DNA Library cutoff value is invalid. It must be a number between 0 and 1")
    Good_Ancestor_Nodes = Select_Ancestor_Nodes(dirname)
    for node in Good_Ancestor_Nodes: #For every node of sufficently high quality, we're going to make a DNA template with uncertianty cutoffs.
        node_request=[] # A request is a list of (list of amino acids needed) at each position
        Positions = ASR_Statefile_Dict[node]
        for i,pos in enumerate(Positions): # For every position,
            if (Binary_Statefile_Dict[node][i][0]) < 0.5: #If this position isn't a gap as determined by the Binary ASR
                pos_AAs=[]
                for j,prob in enumerate(pos): #For every probability at that position,
                    if (prob) > Cutoff: # If the amino acid is above the threshold
                        pos_AAs.append(AA_key[j]) # Record that AA at that position
                if not bool(pos_AAs): # Empty lists evaluate as false
                    pos_AAs=AA_key[pos.index(max(pos))]# If none of the amino acids have a high enough probability to pass the threshold, we'll just record the most likely one.
                node_request.append(pos_AAs)
        Degenerate_DNA = Build_DNA_Sequence(node_request)#Make a DNA sequence with degnerate bases
        with open(f"{dirname}/DNA_Library_{Cutoff*100}%_Cutoff.fasta",'a+') as fout:
            fout.write(f">{node}\n{Degenerate_DNA}\n")

def Write_Confidences(dirname,ASR_Statefile_Dict,Binary_Statefile_Dict):
    nodes_data={}
    for node,cons_list in ASR_Statefile_Dict.items(): #node is name, cons_list is list of (list of AA confidence values) for each position
        node_confidences=[]
        for i,pos in enumerate(cons_list): #pos is AA distribution the position i of ancestor at node node
            if (Binary_Statefile_Dict[node][i][0])>0.5:  #If this position at this node has a greater than 50% chance of being a gap
                node_confidences.append(Binary_Statefile_Dict[node][i][0]) #append the confidnce that there is a gap.
            else: #if there's an amino acid
                node_confidences.append(max(pos)) #append the confidnce of that amino acid.
        #confidences is now a list of the confidence of the most likely state for each position at the node.
        nodes_data[node]=((sum(node_confidences)/len(node_confidences)), (len([i for i in node_confidences if i < 0.85]))) #Each node has a tuple which is the average confidence and the number of positions where confidence is below 85%
    with open(f"{dirname}/IQTree_ASR/Ancestral_Sequence_Confidences.csv","w+") as fout:
        fout.write("Node,Average Confidence,Positions below 85 percent confidence\n")
        for node,data in nodes_data.items():
            fout.write(f"{node},{round((data[0]*100),3)},{data[1]}\n")
    with open(f"{dirname}/IQTree_ASR/Ancestral_Sequence_Confidences.txt","w+") as fout:
        fout.write(f"For the ASR in {dirname} overall:\n\
            {len(ASR_Statefile_Dict)} nodes have an average confidence of {round((sum([n[0] for n in nodes_data.values()])/len(nodes_data)*100),2)}% \n\
            These sequences have an average of {round(sum([n[1] for n in nodes_data.values()]) / len(nodes_data),1)} positions below 85% confidence out of {len(node_confidences)} total positions.\n\
            The topology of these nodes can be found at {dirname}/IQTree_Phylo/Phylo.contree.\n")
    try:
        # Make the sequence prediction histogram
        plt.rcParams["figure.figsize"] = [7.00, 5.00]
        plt.rcParams["figure.autolayout"] = True
        n_bins=101
        # Plot the histogram
        data=[round((n[0]*100),3) for n in nodes_data.values()]
        plt.hist(data,n_bins)
        plt.title('Average Confidnece of each Ancestral Sequence')
        plt.xlabel("Confidence (%)")
        # Save the histogram
        plt.savefig(f"{dirname}/IQTree_ASR/Confidences.png")
        plt.clf()
    except:
        pass

directory='ASR' #Change this to the name of a directory where you want all your resutls output.

if len(sys.argv)==1:
    print("No input file provided. Use option \"help\" to see how to use this program.")
elif (sys.argv[1] == "help") or (sys.argv[1] == "-h") or (sys.argv[1] == "HELP")or (sys.argv[1] == "h"):
    print(Help_String)
    print("\n\n")
    print(Software_Prerequistes)
    print("\n\n")
    print(Directory_Structure)
    print("\n\n")
else: 
    try:
        if len(sys.argv)==2:
            print("User parameters were not proveded. Proceeding with default parameters:\n   Desired Final Library size: 500 sequences.\n   Suppliment Cutoff: 0.7")
            Final_Library_Size=500
            Suppliment_Cutoff=0.7
        else:
            Final_Library_Size=int(sys.argv[2])
            Suppliment_Cutoff=float(sys.argv[3])
            if not (Final_Library_Size>20) and (0.9>Suppliment_Cutoff>0.4):
                print("User provided parameters were not meaningul values.")
                raise ValueError("User provided parameters were not meaningul values.")
    except:
        print("User parameters could not be used. Proceeding with default parameters:\n   Desired Final Library size: 500 sequences.\n   Suppliment Cutoff: 0.7")
        Final_Library_Size=500
        Suppliment_Cutoff=0.7

    if not '.' in sys.argv[1]:
        sequence = sys.argv[1].upper().replace("-",'')
        if any([True for n in sequence if n not in AA_key]):
            raise ValueError("Sequence option read as raw sequence - The provided sequence could not be used. Please be sure no amino acids are \'X\'")
        Blastp_out_name = BlastP(directory,sequence)
        Final_Name = Sequence_Processing(directory,Blastp_out_name, sequence)
    else: 
        if not (os.path.exists(sys.argv[1])):
            raise ValueError("The specified input fasta file does not exist.")
        try:
            User_sequences={}
            temp_User_sequences=fasta2dict(sys.argv[1])
            for key,item in temp_User_sequences.items():
                if '.' in key:
                    key=(key.split("."))[0]
                User_sequences.update({key:item})
        except:
            raise ValueError("The file could not be read as a fasta file.")
        if len(User_sequences)==1:
            User_sequences=list(User_sequences.values())[0]
        else:
            print("Detected input fasta contained multiple sequences. This is an acceptable input, but is less stable and more prone to errors.")
        Blastp_out_name = BlastP(directory,User_sequences)
        Final_Name = Sequence_Processing(directory,Blastp_out_name,User_sequences)

    IQTree_Phylo(directory, Final_Name)
    ASR_Statefile_Dict = IQTree_ASR(directory, Final_Name)
    Binary_Statefile_Dict = Binary_Gap_Analysis(directory, Final_Name)
    Write_Confidences(directory,ASR_Statefile_Dict,Binary_Statefile_Dict)
    Make_Uncertianty_Libraries(directory,ASR_Statefile_Dict,Binary_Statefile_Dict)
