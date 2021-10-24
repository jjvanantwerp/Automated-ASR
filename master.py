# Semi-Automation of the Ancestral Sequence Reconstruction Workflow
# This is a program that constructs rough estimates of ancestral sequence reconstruction from an amino acid sequence.
# This program requires instalation of the Bioservieces module
# Written by James VanAntwerp in September 2021 - vanant25@msu.edu
# Written for the Woldring Lab, Michigan State University in East Lansing, Michigan, USA.
'''
    The program will make the following directory structure:
    ROOT/
    |--master.py 
    |--AutoASR/
    | |--HexKeytoNames.csv
    | |--CD-Hit.fasta.clstr
    | |--BlastP_XML
    | |--*.fasta
    | |--IQTree_Phylo/
    | | |--Post_MAFFT_Cleaned.fasta
    | | |--Phylo.*
    | |--IQTree_ASR/
    | | |--Final_Sequences.fasta
    | | |--ASR.*
    | |--IQTree_Binary/
    | | |--Binary_Alignment.fasta
    | | |--Binary.*
    | |--Primers/
    | | |--Node*
    '''

from Bio.Blast import NCBIWWW
import xml.etree.ElementTree as ET
import os
from string import ascii_letters

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

def fasta2dict(fasta_path,return_dict={},min_len=0): ### Reads in fasta file and renames certain sequences based on forbidden characters in IQ Tree as needed
    # Read in the file and prepare some variables
    with open(fasta_path) as infile:
        fastafile = infile.readlines()
    working_sequence = ''
    key = None
    for line in fastafile: # For every line in the input fasta
        if line[0] == '>': # Check if it's a name
            if working_sequence != '':  # If we already have a working sequence, another name indicates we're done. Otherwise record the name
                if len(working_sequence.replace('-','')) >= min_len:
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
            if len(working_sequence.replace('-','')) >= min_len:
                return_dict[key] = working_sequence
    if len(return_dict)==0:
        raise ValueError(f"The provided file ({fasta_path}) was not a fasta file.")
    return return_dict

def dict2fasta(fasta_d,fname): ### Saves dictionary to fasta where the dictionary key is the protein name and the value is the sequence
    with open(fname,'w+') as out_fasta:
        for key,seq in fasta_d.items():
            out_fasta.write(f'>{key}\n{seq}\n')

def Is_Valid_AA(AA):
    not_AAs = ['B','O','J','U','Z']
    if isinstance(AA,str):  # This block lets us evaluate strings of amino acids
        return(((AA in ascii_letters) or( AA == '-')) and (AA not in not_AAs))
    if isinstance(AA,list) and isinstance (AA[0],str): # This block lets us evaluate lists of amino acids
        return all(((i in ascii_letters) or( i == '-')) and (i not in not_AAs) for i in AA)
    if isinstance(AA,list) and isinstance (AA[0],list):  # This block lets us evaluate lists of lists of amino acids
        return all( all(((i in ascii_letters) or( i == '-')) and (i not in not_AAs) for i in lst) for lst in AA)
    else:
        raise ValueError ("A bad type was evaluated as an amino acid list....")

def Is_Valid_Codon(codon):
    dna_leters = ['a','c','t','g','r','y','m','k','s','w','h','b','v','d','n']
    if isinstance(codon,str):
        return (((len(codon)%3)==0) and all(i in dna_leters for i in codon))
    if isinstance(codon,list):
        return all(((len(c)==3) and all(i in dna_leters for i in c)) for c in codon)

def rename_keys(fasta_d): ### This funciton will re-write keys (names) to remove unwanted charecters and return a renamed (new keys) dictionary
    names = ['Abl','Src','Yes','LCK','Syk','CSK',
             'MATK','Frk','SRMS','BTK','ITK','ZAP-70',
             'ZAP70','Lyn','BMX','FER','FES','Fyn''FGR',
             'Txk','FLK','HCK',]
    renamed_dict = {}
    char_to_ignore = [ ',' , '[' , ']' , '(' , ')' , ';' ]
    char_to_replace = [ ' ' , '|' ]
    with open('Renamed_Fasta_Keys.txt','w+') as outfile:
        for old_key in fasta_d.keys():
            new_key = ''
            split_key = old_key.split(' ')
            for split in split_key:
                for name in names:
                    if name.lower() in split.lower():
                        new_key += split
            if new_key == '':
                new_key += split_key[0]
            if '[' in old_key:
                new_key += '_' + old_key.split('[')[-1][:-1].replace(' ','_')
            for char in char_to_replace:
                new_key = new_key.replace(char,'_')
            for char in char_to_ignore:
                new_key= new_key.replace(char,'')
            renamed = False
            count = 1
            
            new_key0 = new_key
            
            while not renamed:
                if new_key not in renamed_dict.keys():
                    renamed_dict[new_key] = fasta_d[old_key]
                    outfile.write(old_key)
                    outfile.write('\n Was changed to')
                    outfile.write(new_key)
                    outfile.write('\n\n')
                    renamed = True
                else:
                    new_key = new_key0 + '_' + str(count)
                    count += 1
    return renamed_dict

def seq_name_to_hex(fasta_d): ### This function will give each fasta object in fasta_d a hex key name and return the renamed dictionary
    renamed_dict = {}
    count = 1
    model_organisms = ['homo sapiens', 'mus musculus', 'drosophila melanogaster',
            'caenorhabditis elegans','danio rerio','nicotiana tabacum',
            'Triticum aestivum']
    with open('Hexidecimal_Numbered_Fasta_Keys.txt','w+') as outfile:
        for name in fasta_d.keys():
            hexidecimal_key = ''
            if 'predicted' in name.lower() or 'hypothetical' in name.lower():
                hexidecimal_key += 'P'
            organism = name.split('[')[-1][:-1]
            if organism.lower() in model_organisms:
                hexidecimal_key += 'K'
            hexidecimal_key += str(hex(count).lstrip('0x'))
            count += 1
            renamed_dict[hexidecimal_key] = fasta_d[name]
            outfile.write(name)
            outfile.write('\n')
            outfile.write(hexidecimal_key)
            outfile.write('\n')
            renamed = True
    return renamed_dict

def CD_Hit_from_Fasta(dir_with_fastas='Collected Sequences'): ### This function will take all fasta (or valid txt) files in the given directory and run CD-Hit twice, at 98 and 70 percent.
    fasta_files = [dir_with_fastas+'/'+fin for fin in os.listdir(dir_with_fastas) if (fin.split('.')[-1] == 'txt' or fin.split('.')[-1] == 'fasta')]
    
    fasta_dictionary={}
    
    for fasta in fasta_files: #Turn fasta files into dictionaries (This currently overwrites the dictionary)
        try:
            fasta_dictionary = fasta2dict(fasta,d=fasta_dictionary,not_allow_PDB=False,min_len=400)
        except:
            print(f"ERROR - Unable to load fasta:{fasta}")
            pass
    
    new_fasta_named = rename_keys(fasta_dictionary)
    new_fasta_nums = seq_name_to_hex(fasta_dictionary)
    fname = 'Init_Seq_Results/All_Seqs.fasta'
    dict2fasta(new_fasta_nums,fname)
    #Call CD-Hit at both 0.7 and 0.98 cutoff.
    call(['cd-hit', '-i', fname, '-o', f'{fname[:-6]}_0.70_Cluster.fasta', '-c', '0.7'])
    call(['cd-hit', '-i', fname, '-o', f'{fname[:-6]}_0.98_Cluster.fasta', '-c', '0.98'])    

def find_single_clusters(): ### This function recordes single clusters in 'Init_Seq_Results/SingleClusters.txt' using original names
    cd_fname = 'Init_Seq_Results/All_Seqs_0.70_Cluster.fasta'
    
    with open(cd_fname+'.clstr', 'r') as infile:
        clusterText = infile.readlines()
    
    with open('Numbered_Fasta_Keys.txt') as infile: #Here's a call to the Numbered Fasta Keys if that is ever changed.
        text = infile.readlines()
        names = {text[i*2+1].rstrip(): text[i*2].rstrip() for i in range(int(len(text)/2))}
        
    with open('Init_Seq_Results/SingleClusters.txt','w') as outfile: #This writes the single-clustered sequenes at a 70% cutoff.
        count = 0
        for line in clusterText:
            if line[0] == '>':
                if count == 1:
                    key = prev_line.split('\t')[1].split(',')[1].rstrip().lstrip(' >').split(' ')[0].rstrip('...').lstrip('>')
                    outfile.write(names[key].split(' ')[0])
                    outfile.write('\n')
                count = 0
                
            else:
                count += 1
                prev_line = line
        if count == 1:
            key = prev_line.split('\t')[1].split(',')[1].rstrip().lstrip(' >').split(' ')[0].rstrip('...').lstrip('>')
            outfile.write(names[key].split(' ')[0])

def align_prototypes(): ### This function selects useful or otherwise representative sequences from CD-Hit and MAFFT aligns them.
    fname = 'Init_Seq_Results/All_Seqs.fasta'
    cd_fname = 'Init_Seq_Results/All_Seqs_0.98_Cluster.fasta'
    with open(cd_fname+'.clstr', 'r') as infile:
        clusterText = infile.readlines()
        
    cd_index = [i for i in range(len(clusterText)) if clusterText[i][0] =='>']
    cd_index.append(len(clusterText))
    Clusters = {int(clusterText[cd_index[i]].rstrip().lstrip('>Cluster ')): # Clusters is a dictionary of all 0.98 clusters
               [clusterText[j].split('\t')[1].split(',')[1].rstrip().lstrip(' >').split(' ') for j in range(cd_index[i]+1,cd_index[i+1])]for i in range(len(cd_index)-1)} 
    
    full_fasta = fasta2dict(fname)
    
    reselected_fasta = {}
    for c in Clusters: #Determine the most useful sequence(s) in the cluster, often the representative as chosen by CD-Hit. 
        # These were earmarked back when they were convereted to hexadecimal names.
        choices = []
        scores = []
        reselected = False
        for seq in Clusters[c]:
            name = seq[0].rstrip('...')
            if name[0] == 'K': # Select any with a 'K' (model organism)
                reselected_fasta[name] = full_fasta[name]
                reselected = True
            elif name[0] != 'P': # Predicted sequences have 'P'
                s = seq[-1]
                if s == '*': # this indicates a prototype
                    s = 200
                else:
                    s = float(seq[-1].rstrip('%'))
                choices.append(name)
                scores.append(s)
            elif name[0] == 'P' and seq[-1] == '*': # If we have to, we'll use predicted sequences....
                choices.append(name)
                scores.append(0)
        if not reselected: # If we haven't reselected yet, chose the best one.
            final_selection = choices[scores.index(max(scores))]
            reselected_fasta[final_selection] = full_fasta[final_selection]
    
    dict2fasta(reselected_fasta,'Init_Seq_Results/ReselectedReps.fasta') # Make a fasta with the reselected prototype sequences
    mafft_cline = MafftCommandline('mafft', input='Init_Seq_Results/ReselectedReps.fasta', thread=6) # MAFFT align that fasta
    stdout, stderr = mafft_cline()
    with open('Init_Seq_Results/ReselectedReps_Aligned.fasta', "w") as handle:
        handle.write(stdout) # Write the MAFFT alignment

def NCBI_to_XML(diroutname,sequence,hits=1000,expect_value=0.30): #Interact with BlastP, and record the XML
    # For the given sequence, we will run a BlastP search and parse the XML return to make a multi-fasta file
    print("Acessing the NCBI database....")
    blast_result_handle = NCBIWWW.qblast("blastp","nr", sequence,hitlist_size=hits,expect=expect_value)
    with open(f"./{diroutname}/BlastP_XML","w+") as fout:
        fout.write(blast_result_handle.read())
    blastp_xml = ET.parse(f"./{diroutname}/BlastP_XML")
    #os.remove(f"./{diroutname}/BlastP_XML")
    return(blastp_xml)

def Parse_BlastP_XML(diroutname,blastp_xml,sequence): #Parse the BlastP XML - record Hex keys and Fasta
    Hex_Fasta_Dict={}
    with open(f"./{diroutname}/HexKeytoNames.csv","w+") as fout:
        hexcount=0
        fout.write(f"Hexadecimal Name,BlastP Name\n")
        #Parsing the XML object, looking for hits
        for hit in blastp_xml.findall('./BlastOutput_iterations/Iteration/Iteration_hits/Hit'):
            name = (hit.find('Hit_id')).text #I've tried to also add the Hit_accession, but I can't access that form the XML for some reason
            seq = (hit.find('Hit_hsps/Hsp/Hsp_hseq')).text
            #If the sequence doesn't have unknowns amino acids (annoying) then record it.
            #The optional second method also removes exceptionally short or long sequences - be sure to synch with the code ~13 lines below
            #if (("X" not in seq) and (len(seq)<((1+length_cutoff)*User_Sequence_Length)) and (len(seq)>((1-length_cutoff)*User_Sequence_Length))):
            if ("X" not in seq):
                hexcount+=1
                hexkey='Seq_'
                if 'predicted' in name.lower() or 'hypothetical' in name.lower():
                    hexkey += 'P'
                hexkey += str(hex(hexcount).lstrip('0x'))
                fout.write(f"{hexkey},{name}\n")
                Hex_Fasta_Dict[hexkey]=seq
    with open(f"./{diroutname}/BlastP_Results.fasta","w+") as blastp_file:
        blastp_file.write(f">User_Sequence\n{sequence}\n")
        for Hex,Sequence in Hex_Fasta_Dict.items():
            #if (len(Sequence)<((1+length_cutoff)*User_Sequence_Length)) and (len(Sequence)>((1-length_cutoff)*User_Sequence_Length)):
            blastp_file.write(f">{Hex}\n")
            blastp_file.write(f"{Sequence.replace('-','')}\n")#We remove all gaps, because CD-Hit cannot handle gaps.

def BlastP(diroutname,sequence,hits,expect_value): #This function takes an amino acid sequence, submitts a BlastP search, and records the result in a fasta file
    if not (os.path.isdir(diroutname)):
        os.mkdir(diroutname)
    sequence=sequence.replace('-','')
    sequence=sequence.replace('X','')
    #User_Sequence_Length=len(sequence)
    if all([char for char in sequence if (char.isalpha())]):
        try:
            blastp_xml=NCBI_to_XML(diroutname,sequence,hits,expect_value)
        except:
            raise RuntimeError("There was an error fetching the BlastP results.")
            print("There was an error fetching the BlastP results")
        try:
            # Now, we parse the XML object and make a multi-fasta file. We will assign each sequence a hexadecimal name.
            # We also write the fasta file which is the result of our BlastP search  
            Parse_BlastP_XML(diroutname,blastp_xml,sequence)
        except:
            raise RuntimeError("There was an error recording the BlastP Results")
            print("There was an error recording the BlastP Results")
    else:
        raise ValueError("Invalid sequence submitted for BlastP search")
        print("Invalid sequence submitted for BlastP search")
    return("BlastP_Results.fasta")

def CDHit(dirname,finname,identity=0.95): #Run CD-Hit
    #Determine the right n for the identity
    if (identity<=0.4 or identity>1 ):
        raise ValueError("The CD-Hit identity is invalid")
        print("The CD-Hit identity is invalid")
    elif identity>0.7:
        n=5
    elif identity>0.6:
        n=4
    elif identity>0.5:
        n=3
    elif identity>0.4:
        n=2
    try:
        os.system(f'cd-hit -i {dirname}/{finname} -o {dirname}/CD-Hit.fasta -c {identity} -n {n}') #Run CD-Hit
    except:
        raise RuntimeError("There was an error running CD-Hit.")
        print("There was an error running CD-Hit.")
    return("CD-Hit.fasta")    

def MAFFT(dirname,finname,sequence): #Run MAFFT
    #Consider using subprocesses instead of os, because we can redirect stdout
    with open(f"{dirname}/{finname}",'r') as fin: #using read-only mode
        lines=fin.readlines()
    did_CDHit_remove=True #Default behavior is to add User_Sequence. CD-Hit can remove it, and that causes problems later.
    for line in lines:
        if line=='>User_Sequence':
            did_CDHit_remove = False
            break
    if did_CDHit_remove:
        with open(f"{dirname}/{finname}",'a') as fout: #We use the append mode to avoid overwriting
            fout.write(f'>User_Sequence\n{sequence}\n')
    try:
        #MAFFT section
        #os.system(str) Executes command of str in the directory from which the python script is executed.
        os.system(f'mafft {dirname}/{finname} > {dirname}/Mafft_Alignment.fasta') #align the sequences passed into the function
    except:
        raise RuntimeError("There was an error creating the sequence alignemnt.")
        print("There was an error creating the sequence alignemnt.")
    return("Mafft_Alignment.fasta")

def Trim_N_C_Termini_Percent (fasta_dict,terminus_cutoff=0.20): #Trim termini based on % of alignment that has a sequence. Currently not used.
    sequences_list = [i for i in fasta_dict.values()]
    gap_percentage_by_pos=[]#Represents the % of each position that is a gap.
    for position in range(len(sequences_list[0])): #For every position in the sequences,
        pos_gap_tally=0
        for i in range(len(fasta_dict)): #For every sequence
            if sequences_list[i][position] =='-': #If that position has a gap
                pos_gap_tally+=1
        gap_percentage_by_pos.append(pos_gap_tally/len(sequences_list))
    #Count N-terminus trim length
    len_to_remove_N=0
    len_to_remove_C=0
    for pos in gap_percentage_by_pos:
        if pos > terminus_cutoff:
            len_to_remove_N+=1
        else:
            break
    #Count C-terminus trim length 
    for i in range(len(gap_percentage_by_pos)-1,-1,-1):
        if gap_percentage_by_pos[i]>terminus_cutoff:
            len_to_remove_C+=1
        else:
            break
    for name,seq in fasta_dict.items():
        fasta_dict.update({name : (seq[len_to_remove_N:-(len_to_remove_C)]) })

def Trim_N_C_Termini (fasta_dict): #Trim termini based on User_Sequence
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

def Remove_Insertions (fasta_dict,User_Sequence,deletion_percentage=0.025,termini_length=0.1): #Remove sequnces that cause insertions in the alighment
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

def Remove_Deletions (fasta_dict,User_Sequence,deletion_percentage=0.025,termini_length=0.1): #Remove sequnces that cause insertions in the alighment
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
            
def Post_MAFFT_processing(dirname,finname,deletion_percentage=0.01,misalignment_cutoff=0.2): #Modifications after the alignment, mostly having to do with gaps.
    fasta_dict = fasta2dict(f"{dirname}/{finname}")
    # These functions **MODIFY** the fasta_dict by doing what their names say they do.
    Trim_N_C_Termini(fasta_dict)
    old_user_length=0
    User_Sequence=fasta_dict.get("User_Sequence")
    while len(User_Sequence)!=old_user_length: #Each of these functions mutate the alignment, so it's important to repeat until we've finsihed finding all misalignments.
        old_user_length=len(User_Sequence)#Store length
        Remove_Insertions(fasta_dict,User_Sequence,deletion_percentage)#Insertions
        Clean_all_gaps(fasta_dict)#Clean
        User_Sequence=fasta_dict.get("User_Sequence")#Update
        Remove_Deletions(fasta_dict,User_Sequence,deletion_percentage)#Deletions
        Clean_all_gaps(fasta_dict)#Clean
        User_Sequence=fasta_dict.get("User_Sequence")#Update
    dict2fasta(fasta_dict,f'{dirname}/Post_MAFFT_Cleaned.fasta')
    '''
        Boolean_Gap_Positions=[]
        User_Sequence=fasta_dict.get("User_Sequence")
        #Tally positions with and without gaps in the provided sequence in the alignment
        for pos in User_Sequence:
            if pos =='-':
                Boolean_Gap_Positions.append(False)
            else:
                Boolean_Gap_Positions.append(True)
        #Check all the other sequences to see how often they diagree
        for Name,Seq in fasta_dict.items():#For all sequences
            score=0 #starting with a clean slate
            for i,residue in enumerate(Seq):
                if Boolean_Gap_Positions[i]: #If the User sequence has an ammino acid and this sequcne does not (deletion)
                    if (residue == '-'):
                        score+=1 #Mark it down
                elif(residue != '-'): #Or if the User sequence has a gap and this sequnce does not (insertion)
                    score+=1 #Mark it down
            if (score>(misalignment_cutoff*len(Seq))):#If more than the specified % of positions disagree,
                Misaligned_Sequences.append(Name) #Add it to the blacklist
        #Remove the misaligned sequences
        for Name in Misaligned_Sequences:
            fasta_dict.pop(Name)
        #Remove all gaps in the alignment, and re-align.
        for key,sequence in fasta_dict.items(): #Remove all the gaps
            fasta_dict.update({key:''.join([ char for char in sequence if char != '-' ])})
        #store as a file for MAFFT alignment, run the alignment, and remove the temporary file.
        dict2fasta(fasta_dict,'temporary1.fasta')
        try:
            #os.system(str) Executes command of str in the directory from which the python script is executed.
            os.system(f'mafft temporary1.fasta > temporary2.fasta') 
        except:
            raise RuntimeError("There was an error aligning the cleaned sequences.")
            print("There was an error ligning the cleaned sequences.")
        os.remove('temporary1.fasta')
        '''
    '''
        Misaligned_Sequences=[]
        #Clean up the sequences that remain
        gap_tally = [0]*(len(Boolean_Gap_Positions))
        number_of_sequences=len(fasta_dict)
        positions_to_pop=[]
        for key in fasta_dict.keys(): #For each position in the alignment 
            i=0
            for residue in fasta_dict[key]:
                if (residue == '-'): #Tally the numer of sequences that have a gap
                    gap_tally[i]+=1
                i+=1
        i=0
        for position in gap_tally: #For each position in the alignment
            if ((position/number_of_sequences)==1): #If all sequences have a gap,
                positions_to_pop.append(i) #remember to remove that position from the whole alignment.
            i+=1
        for key in fasta_dict.keys(): #For each sequence in the alignment
            #Reassemble the sequence using only positions from the alignment that were not majority-gaps
            fasta_dict.update({key:''.join([ fasta_dict.get(key)[i] for i in range(len(fasta_dict.get(key))) if i not in positions_to_pop ])})
        '''    
    return('Post_MAFFT_Cleaned.fasta')

def IQTree_Phylo(dirname,finname): #finname is the alignment
    if not os.path.isdir(f"{dirname}/IQTree_Phylo"): #Make the directory
        os.mkdir(f"{dirname}/IQTree_Phylo")
    try:
        #IQTree section
        os.system(f'iqtree -s {dirname}/{finname} -pre {dirname}/IQTree_Phylo/Phylo -alrt 20000 -bb 20000 -nt AUTO')
    except:
        raise RuntimeError("There was an error building the phylogeny in IQTree.")
        print("There was an error building the phylogeny in IQTree.")

def IQTree_ASR(dirname,finname): #finname is the alignmen
    model=''
    try: #this can fail because for some reason HPCC tries to run it before the IQTree_Phylo has finished running.
        with open(f"{dirname}/IQTree_Phylo/Phylo.iqtree") as fin: #we need to go to the .iqtree file from the phylogony and find the best model for this data set. This will save a ton of time.
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
        os.system(f'iqtree -s {dirname}/{finname} -te {dirname}/IQTree_Phylo/Phylo.treefile -pre {dirname}/IQTree_ASR/ASR -asr -m {model} -redo -nt AUTO')
    except:
        raise RuntimeError("There was an error conducting ASR in IQTree.")
        print("There was an error conducting ASR in IQTree.")

def Select_Nodes_Suppliment(dirname,finname,cutoff=50):
    #This will evaluate the tree, which is in Newick format, and select nodes of high enough confidence to reconstruct an ancestral library.
    #The .treefile has all of the confidence values, and the node names.
    with open(f'{dirname}/{finname}') as treefin:
        treefile = treefin.readlines()[0]
    stinky_nodes=[] #these will be the nodes with poor values. This isn't quite the same value between *.contree and *.trefile, but close
    nodes=treefile.split(')')#Let's split off just the information for each node, which is stored after every close parenthiesis.
    nodes.pop(0)#The above will split off a first section without a node. This does not cause any nodes to be lost.
    for i,node in enumerate(nodes): #rearanging the newick format a little. Each node in nodes will now be stored as "name,int,float"
        try:
            node_info=node.split(',')[0]
            #nodes[i]=f"{node_info.split('/')[0]},{(node_info.split('/')[1]).split(':')[0]},{node_info.split(':')[1]}"
            nodes[i]=node_info.split('/')[0]
            if int((node_info.split('/')[1]).split(':')[0]) < cutoff: #Select the nodes which have a value above the cutoff.
                #tup = ((node_info.split('/')[0]),(node_info.split('/')[1]).split(':')[0])
                stinky_nodes.append(nodes[i])
        except:
            if node.strip() == 'Node1;':
                pass #The root node doesn't have info - this prevents an error in handling that from causing problems
            else:
                print(f'{node} caused problems :(')
                print("There was an error handling the Newick Tree returned by IQTree - This likely indicates a bug with IQTree.")
                raise ValueError("There was an error handling the Newick Tree returned by IQTree - This likely indicates a bug with IQTree.")
    #return [n.split(',')[0] for n in nodes if n.split(',')[0] not in stinky_nodes] 
    return [n for n in nodes if n in stinky_nodes] #return the list of node names whose value is above the cutoff.

def Select_Ancestor_Nodes_OLD(dirname,finname,cutoff=85):
    #This will evaluate the tree, which is in Newick format, and select nodes of high enough confidence to reconstruct an ancestral library.
    #The .treefile has all of the confidence values, and the node names.
    with open(f'{dirname}/{finname}') as treefin:
        treefile = treefin.readlines()[0]
    Good_Nodes=[] #these will be the nodes with poor values. This isn't quite the same value between *.contree and *.trefile, but close
    nodes=treefile.split(')')#Let's split off just the information for each node, which is stored after every close parenthiesis.
    nodes.pop(0)#The above will split off a first section without a node. This does not cause any nodes to be lost.
    for i,node in enumerate(nodes): #rearanging the newick format a little. Each node in nodes will now be stored as "name,int,float"
        try:
            node_info=node.split(',')[0]
            #nodes[i]=f"{node_info.split('/')[0]},{(node_info.split('/')[1]).split(':')[0]},{node_info.split(':')[1]}"
            nodes[i]=node_info.split('/')[0]
            if int((node_info.split('/')[1]).split(':')[0]) >= cutoff: #Select the nodes which have a value above the cutoff.
                #tup = ((node_info.split('/')[0]),(node_info.split('/')[1]).split(':')[0])
                Good_Nodes.append(nodes[i])
        except:
            if node.strip() == 'Node1;':
                pass #The root node doesn't have info - this prevents an error in handling that from causing problems
            else:
                print(f'{node} caused problems :(')
                print("There was an error handling the Newick Tree returned by IQTree - This likely indicates a bug with IQTree.")
                raise ValueError("There was an error handling the Newick Tree returned by IQTree - This likely indicates a bug with IQTree.")
    #return [n.split(',')[0] for n in nodes if n.split(',')[0] not in stinky_nodes]
    return Good_Nodes #return the list of node names whose value is above the cutoff.

def Select_Ancestor_Nodes(dirname): #Select ancestral nodes which have a high enough confidence to resurect ancestors from them.
    #This will evaluate the tree, which is in Newick format, and select nodes of high enough confidence to reconstruct an ancestral library.
    #The Phlyo.treefile has confidence values as (SH-aLRT/UFB), and the ASR.treefile has node names.
    with open(f'{dirname}/IQTree_ASR/ASR.treefile') as ASRtreefin:
        ASRtreefile = ASRtreefin.read()
    with open(f'{dirname}/IQTree_Phylo/Phylo.treefile') as Phylotreefin:
        Phylotreefile = Phylotreefin.read()
    confident_nodes=[] #these will be the nodes with good values. This isn't quite the same value between *.contree and *.trefile, but close
    ASRnodes=[n.split(',')[0].split(':') for n in ASRtreefile.split(')')]  #Let's split off just the information for each node, which is stored after every close parenthiesis.
    Phylonodes=[n.split(',')[0].split(':') for n in Phylotreefile.split(')')]
    ASRnodes.pop(0) #the first split is not actually a node.
    Phylonodes.pop(0)#the first split is not actually a node.
    for i in range(len(ASRnodes)-1): #We exclude the last node, because that's the root, which will mess up the below lines.
        if (float(Phylonodes[i][0].split('/')[0]) > 80) and (float(Phylonodes[i][0].split('/')[1]) > 95): #If the SH-aLRT >80% and the ultrafast bootstraping is >95%
            confident_nodes.append(ASRnodes[i][0]) #Record the name of the high-confidence nodes.
    return confident_nodes #return the list of node names whose value is above the cutoff.
    
def Statefile_to_Dict(dirname): #{NodeX:[[probs],[probs],[probs],...]}
    statefile_dict={} #This is a dicitonary made out of the statefile - its keys are the node names, 
    #and its values are a list of tuples with the amino acid and list of amino acid distributions at each position.
    node=[]
    #This will parse the .state file for the desired nodes to get their AA distribution
    with open(f'{dirname}/IQTree_ASR/ASR.state') as statefin: #Read in each line, skipping the header.
        statelines=[]
        for i,line in enumerate(statefin):
            if i>8:
                statelines.append(line)
    #Now let's pull the data from each line
    working_node=statelines[0].split()[0] #prime the working node
    for line in statelines: # For every line in the state file
        line_list = line.split() #Break up the line into columns
        if working_node == line_list[0]: #If we're still working on the same node
            node.append(line_list[3:]) #record the probability distribution of the position
        else: #If we've come to the end of a node,
            statefile_dict[working_node]=node #Add a key-value pair to the statefile dictionary that is the node name and the node's list
            working_node = line_list[0] #update the working_node value
            node=[] #Clear the working node list
            node.append(line_list[3:]) #add to the node list a touple of the amino acid and the distribution of amino acids at that position.
    statefile_dict[working_node]=node #Be sure to add the last node into the dictionary too!
    return(statefile_dict)

def Statefile_to_Dict_AAs(dirname,finname): #{NodeX:[AA,AA,AA,...]}
    statefile_dict={} #This is a dicitonary made out of the statefile - its keys are the node names, 
    #and its values are a list of tuples with the amino acid and list of amino acid distributions at each position.
    node=[]
    #This will parse the .state file for the desired nodes to get their AA distribution
    with open(f'{dirname}/{finname}') as statefin: #Read in each line, skipping the header.
        statelines=[]
        for i,line in enumerate(statefin):
            if i>8:
                statelines.append(line)
    #Now let's pull the data from each line
    working_node=statelines[0].split()[0] #prime the working node
    for line in statelines: # For every line in the state file
        line_list = line.split() #Break up the stuff
        if working_node == line_list[0]: #If we're still working on the same node
            node.append(line_list[2]) #record the amino acid at that position
        else: #If we've come to the end of a node,
            statefile_dict[working_node]=node #Add a key-value pair to the statefile dictionary that is the node name and the node's list
            working_node = line_list[0] #update the working_node value
            node=[] #Clear the working node list
            node.append((line_list[2],line_list[3:])) #add to the node list a touple of the amino acid and the distribution of amino acids at that position.
    statefile_dict[working_node]=node #Be sure to add the last node into the dictionary too!
    return(statefile_dict)

def Binary_Gap_Analysis(dirname,finname): #Do a binary ASR to determine where gaps in the ancestral sequences reside
    if not os.path.isdir(f"{dirname}/IQTree_Binary"): #Make the directory
        os.mkdir(f"{dirname}/IQTree_Binary")
    fasta_dict = fasta2dict(f"{dirname}/{finname}")
    binary_dict={}
    for key,seq in fasta_dict.items(): #Make a binary alignment of the input fasta
        binary_dict[key]=(''.join(['0' if aa =='-' else '1' for aa in seq])) #This line is taken from Ben - how to give proper credit?
    dict2fasta(binary_dict,f"{dirname}/IQTree_Binary/Binary_Alignment.fasta")
    try:
        os.system(f'iqtree -s {dirname}/IQTree_Binary/Binary_Alignment.fasta -te {dirname}/IQTree_Phylo/Phylo.treefile -pre {dirname}/IQTree_Binary/Binary -blfix -asr -m GTR2+FO -redo -nt AUTO')
    except:
        print("There was an error determining gaps in the ancestral sequence")
        raise RuntimeError("There was an error determining gaps in the ancestral sequence")
    Binary_Statefile_dict = Statefile_to_Dict(dirname,"IQTree_Binary/Binary.state")
    ASR_Statefile_dict = Statefile_to_Dict(dirname,"IQTree_ASR/ASR.state")
    Consensus_Ancestors_with_Gaps={}
    Pos_with_Gaps={} #dictionary of {NodeX:[list of gaps at NodeX]}
    #Find positions that are actually gaps in the ASR
    for node,item in Binary_Statefile_dict.items():
        gap_pos=[]#list of positions with gaps
        for i,pos in enumerate(item): #each position
            if pos[0]=='0':#If the posotion has majority gap
                gap_pos.append(i)
        Pos_with_Gaps[node]=gap_pos #At each node, record the positions in the ancestral sequence that has majority node.
    #Merge the Sequence ASR with the gap ASR
    for node,item in ASR_Statefile_dict.items():
        consensus_seq=''
        for i,pos in enumerate(item):
            if i in Pos_with_Gaps[node]: #If this position at this node is likely a gap, add a gap to the consensus sequence
                consensus_seq+='-'
            else:
                consensus_seq+=pos[0] #Otherwise, add the amino acid from ASR
        Consensus_Ancestors_with_Gaps[node]=consensus_seq
    Clean_all_gaps(Consensus_Ancestors_with_Gaps)
    dict2fasta(Consensus_Ancestors_with_Gaps,f"{dirname}/Consensus_Ancestors_with_Gaps.fasta")

def Degenerate_Nucleotide_Codon(AA_List, source='EColi'):  # This function takes a list of AAs for ONE POSITION and makes a degenerate codon for them.
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
            elif (len(pos)==2):
                primer_txt+=(AA_Pair_lookup_EColi[ pos[0]+pos[1] ])
            else:
                primer_txt+=(Degenerate_Nucleotide_Codon(pos))
    elif source == 'Human':
        for pos in Primer_Request: # pos is the list of AA at the respective position
            if (len(pos)==1):
                if (pos[0]!='-'): # For a gap in the sequence, we will record nothing.
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

def Make_Uncertianty_Libraries (dirname,statefile_dict,nodes,Cutoff=0.125): #Make Uncertianty Libraries for all ancestral sequences with a high enough confidence at thier node
    if not (os.path.isdir(f"{dirname}/Ancestral_DNA_{Cutoff*100}%_Cutoff")):
        os.mkdir(f"{dirname}/Ancestral_DNA_{Cutoff*100}%_Cutoff")
    Consensus_Ancestors_with_Gaps = fasta2dict(f"{dirname}/Consensus_Ancestors_with_Gaps.fasta") #Be sure to import information about where gaps are in the ancestor sequences.
    for node in nodes: #For every node of sufficently high quality, we're going to make a DNA template with uncertianty cutoffs.
        Ancestor_with_Gaps = Consensus_Ancestors_with_Gaps[node]
        node_seq_request=[] # A request is a list of (list of amino acids needed) at each position
        Positions = Satefile_Dict[node]
        for i,pos in enumerate(Positions): # For every position,
            if Ancestor_with_Gaps[i]!='-': #If this position isn't a gap as determined by the Binary ASR
                pos_AAs=[]
                for j,prob in enumerate(pos): #For every probability at that position,
                    if (float(prob) >= Cutoff): # If the amino acid is above the threshold
                        pos_AAs.append(AA_key[j]) # Record that AA at that position
                if not bool(pos_AAs): # Empty lists evaluate as false
                    # If none of the amino acids have a high enough probability to pass the threshold, we'll just record the most likely one.
                    pos_AAs=AA_key[pos.index(max(pos))]
                node_seq_request.append(pos_AAs)
        Degenerate_DNA = Build_DNA_Sequence(node_seq_request)#Make a DNA sequence with degnerate bases
        with open(f"{dirname}/Ancestral_DNA_{Cutoff*100}%_Cutoff/{node}_DNA_Library.txt",'w+') as fout:
            fout.write(Degenerate_DNA)

HaloTag='SGSAEIGTGFPFDPHYVEVLGERMHYVDVGPRDGTPVLFLHGNPTSSYVWRNIIPHVAPTHRCIAPDLIGMGKSDKPDLGYFFDDHVRFMDAFIEALGLEEVVLVIHDWGSALGFHWAKRNPERVKGIAFMEFIRPIPTWDEWPEFARETFQAFRTTDVGRKLIIDQNVFIEGTLPCGVVRPLTEVEMDHYREPFLNPVDREPLWRFPNELPIAGEPANIVALVEEYMDWLHQSPVPKLLFWGTPGVLIPPAEAARLAKSLPNCKAVDIGPGLNLLQEDNPDLIGSEIARWLSTLEISG'
OATP='-MDQNQHLNKTAEAQPSENKKTR-YCNGLKMFLAALSLSFIAKTLGAIIMKSSIIHIERRFEISSSLVGFIDGSFEIGNLLVIVFVSYFGSKLHRPKLIGIGCFIMGIGGVLTALPHFFMGYYRYSKETNINSSENSTSTLSTCLINQILSLNRASPEIVGKGCLKESGSYMWIYVFMGNMLRGIGETPIVPLGLSYIDDFAKEGHSSLYLGILNAIAMIGPIIGFTLGSLFSKMYVDIGYVDLSTIRITPTDSRWVGAWWLNFLVSGLFSIISSIPFFFLPQTPNKPQKERKA-SLSLHVLETNDEKDQTANLTN--QGKNITK-NVTG-FFQSFKSILTNPLYVMFVLLTLLQVSSYIGAFTYVFKYVEQQYGQPSSKANILLGVITIPIFASGMFLGGYIIKKFKLNTVGIAKFSCFTAVMSLSFYLLYFFILCENKSVAGLTMTYDGNNPVTSHRDV-PLSYCNSDCNCDESQWEPVCGNNGITYISPCLAGCKSSSGNKKP---IVFYNCSCLEVTGLQNRNYSAHLGECPRDDACTRKFYFFVAIQVLNLFFSALGGTSHVMLIVKIVQPELKSLALGFHSMVIRALGGILAPIYFGALIDTTCIKWSTNNCGTRGSCRTYNSTSFSRVYLGLSSMLRVSSLVLYIILIYAMKKKYQEKDINASENG-SVMDEANLESLN-KNKHFVPS--AGADSETHC----------'
Estrogen='SNAKRSKKNSLALSLTADQMVSALLDAEPPILYSEYDPTRPFSEASMMGLLTNLADRELVHMINWAKRVPGFVDLTRHDQVHLLECAWLEILMIGLVWRSMEHPGKLLFAPNLLLDRNQGKCVEGMVEIFDMLLATSSRFRMMNLQGEEFVCLKSIILLNSGVYTFLSSTLKSLEEKDHIHRVLDKITDTLIHLMAKAGLTLQQQHQRLAQLLLILSHIRHMSNKGMEHLYSMKCKNVVPSYDLLLEMLDAHRLHAPT'
directory='HPCC_Run_Oct_23'
sequence=HaloTag.replace('-','')
Blastp_out_name = BlastP(directory,sequence,2000,0.3)
#Blastp_out_name = f"BlastP_Results.fasta"
CDHit_out_name = CDHit(directory,Blastp_out_name)
#CDHit_out_name="CD-Hit.fasta"
MAFFT_out_name = MAFFT(directory,CDHit_out_name,sequence)
#MAFFT_out_name="Mafft_Alignment.fasta"
Post_MAFFT_name = Post_MAFFT_processing(directory,MAFFT_out_name)
#Post_MAFFT_name = 'Post_MAFFT_Cleaned.fasta'
IQTree_Phylo(directory, Post_MAFFT_name)
#print("\n\nPhylogeny Finished\n\n")
IQTree_ASR(directory, Post_MAFFT_name)
#print("\n\nASR Finished\n\n")
Binary_Gap_Analysis(directory, Post_MAFFT_name)
#print("\n\nGaps Evaluated\n\n")
Good_Ancestor_Nodes = Select_Ancestor_Nodes(directory) #Returns a list of nodes
Satefile_Dict = Statefile_to_Dict(directory) #Returns a dictionary out of *.state file
Make_Uncertianty_Libraries(directory,Satefile_Dict,Good_Ancestor_Nodes)

#os.system('afplay /System/Library/Sounds/Glass.aiff')
