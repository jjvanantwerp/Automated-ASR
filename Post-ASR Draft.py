# -*- coding: utf-8 -*-
'''
Post-ASR analysis code for an AutoASR Workflow
This program is meant to be used in a larger bash workflow, and has a few dependancies:
   tqdm, Bio, pickle, urllib, dendropy, and requests
A large section of this code was taken from a project by Patrick Finneran (@author: pfinneran)
The project was suplimented and finished by James VanAntwerp in December 2020
Written for the Woldring Lab, Michigan State University in East Lansing, Michigan, USA.
'''

from dendropy import Tree
import os
from tqdm import tqdm
from statistics import mode
from Bio.Align.Applications import MafftCommandline
from Bio import SeqIO
from Bio.Blast import NCBIWWW
import operator
import seaborn as sns
import pandas as pd
from math import ceil
import pickle
import sys
import io
import requests
import urllib
from datetime import date
import xml.etree.ElementTree as ET
from string import ascii_letters

# Patrick's stuff
def retrieve_ancestors(bali_phy_name,burnin=1000,sampling=1):
    ### Open descendants file and place all the sequence names in the labels list
    Desc_Labels = {}
    Outfiles = {}
    Node_Counts = {}
    Node_Counts_flawed = {}
    SampledTrees = []
    print(f"\n\n{os.listdir('Descendants')}\n\n")
    # On Mac OS, the .DS_Store file in every folder could cause the below loop to fail.
    try:
        os.remove("Descendants/.DS_Store")
    except:
        pass

    # This line will loop through every file in the Descendants folder with the naming convention *_Descendants.*
    for Desc in [x for x in os.listdir('Descendants') if x.split('.')[0].split('_')[-1] == 'Descendants']:
        with open('Descendants/{}'.format(Desc)) as infile:
            Desc2 = Desc.split('.')[0]
            lines = infile.readlines()
            Desc_Labels[Desc2] = [x.split()[0].lstrip('>') for x in lines[0::2]]
            if not os.path.exists('Ancestor_Nodes'):
                os.mkdir('Ancestor_Nodes')
            Outfiles[Desc2] = open('Ancestor_Nodes/{}_Sampling_e{}.fasta'.format(Desc2,sampling),'w')
            Node_Counts[Desc2] = 0
            Node_Counts_flawed[Desc2] = 0

    ### Check for directories made by BAli-Phy
    dirs = [d for d in os.listdir() if d[0:len(bali_phy_name)] == bali_phy_name]
    for directory in dirs:
        print('Parsing {0}'.format(directory))
        with open('{0}/C1.trees'.format(directory)) as infile:
            Trees = infile.readlines()
        ### Initialize values
        get_sequence = False
        Searching = False
        iteration = 0
        Fastas = open('{0}/C1.P1.fastas'.format(directory))
        ### The following line James changed "tqdm" tp "tqdm.tqdm" because of namespace resoloution errors.
        for i, line in enumerate(tqdm.tqdm(Fastas.readlines())):
            if line[0:10] == 'iterations':
                iteration = int(line[12:-1])
                if iteration > burnin and (iteration%sampling == 0):
                    tree = Tree()
                    tree = tree.get(data=Trees[iteration],schema='newick',tree_offset=0)
                    SampledTrees.append(Trees[iteration])
                    Node_Labels = {}
                    for Desc in Desc_Labels.keys():
                        node_label = tree.mrca(taxon_labels=Desc_Labels[Desc]).label
                        Node_Labels[Desc] = node_label
                        Node_Labels[node_label] = Desc
                        leaves = [x.taxon.label for x in tree.find_node_with_label(node_label).leaf_nodes()]
                        new_taxa = set(leaves) - set(Desc_Labels[Desc])
                        if len(new_taxa) > 0:
                            Node_Counts_flawed[Desc] += 1
                        Node_Counts[Desc] += 1
                    Searching = True
                else:
                    Searching = False
            elif Searching:
                
                if line.rstrip().lstrip('>') in list(Node_Labels.values()):
                    get_sequence = []
                    for Desc in Node_Labels.keys():
                        if Node_Labels[Desc] == line.rstrip().lstrip('>'):
                            get_sequence.append(Desc)
                elif line[0] == '>':
                    get_sequence = False
                elif get_sequence:
                    for i in range(len(get_sequence)):
                        Outfiles[get_sequence[i]].write('>{0}|Run_{1}|Iteration_{2}\n'.format(get_sequence[i],directory[len(bali_phy_name)+1:],iteration))
                        Outfiles[get_sequence[i]].write(line.rstrip()+'\n')
            Fastas.close()
    if not os.path.exists('SampledTrees'):
        os.mkdir('SampledTrees')
    with open('SampledTrees/e{}.tree'.format(sampling),'w') as outfile:
        outfile.write(''.join(SampledTrees))
    ### Close all files
    with open('Ancestor_Report_Sampling_e{}.csv'.format(sampling),'w') as outfile:
        outfile.write(','.join(['Ancestor','FlawedNodes','TotalNodes','PercentFlawed'])+'\n')
        for Desc in Outfiles.keys():
            Outfiles[Desc].close()
            perc = Node_Counts_flawed[Desc]/Node_Counts[Desc]
            line = '{},{},{},{:.2%}\n'.format(Desc,Node_Counts_flawed[Desc],Node_Counts[Desc],perc)
            outfile.write(line)

def align_ancestors(mafft_exe='"mafft/mafft.bat"'):
    if not os.path.exists('Ancestor_Nodes_Aligned'):
        os.mkdir('Ancestor_Nodes_Aligned')
    if not os.path.exists('Ancestor_Consensus_Seqs'):
        os.mkdir('Ancestor_Consensus_Seqs')
    # Make sure that on Mac OS the .DS_Store file is removed, or the alignments will fail.
    try:
        os.remove("Ancestor_Nodes/.DS_Store")
    except:
        pass
    for fasta in os.listdir('Ancestor_Nodes'):
        print('Aligning and Analyzing {}'.format(fasta.split('.')[0]))
        in_file = 'Ancestor_Nodes/' + fasta
        outfile = 'Ancestor_Nodes_Aligned/' + fasta.split('.')[0] + '_Aligned.fasta'
        mafft_cline = MafftCommandline(mafft_exe, input=in_file) 
        # Patrick's workflow uses MAFFT, I think we should consider something else, possibly splitting the python code 
        # to use somehting not supported in python like MUSCLE. This would let us keep gaps in the alignment, which BaliPhy acutlally handles.
        stdout, stderr = mafft_cline()
        with open(outfile, "w") as handle:
            handle.write(stdout)
        Seqs = retrieve_consensus(fasta.split('.')[0],outfile)
        fname = 'Ancestor_Consensus_Seqs/{0}/{0} Consensus from BAli-Phy.pdf'.format(fasta.split('.')[0])
        seq_heatmap(Seqs,fname)
        pass

def retrieve_consensus(node_name,fasta):
    D = 'Ancestor_Consensus_Seqs/{}'.format(node_name)
    if not os.path.exists(D):
        os.mkdir(D)
    # James changed this line because he had an 'X' in his alignment which raised problems. He added an 'X' term to the list.
    AAs = ['G', 'A', 'V', 'L', 'I', 'W', 'F', 'Y', 'P',  'M', 'C', 'S', 'T', 'N', 'Q', 'D', 'E', 'H', 'K', 'R', 'x', '-']
    Seqs = []
    infile = open(fasta)  
    seq = ''
    #Parsing the Multifasta - Could use or integrate parse_multifasta
    for line in infile:
        if line[0] == '>':
            if seq != '':
                Seqs.append(seq)
            seq = ''
        else:
            seq = seq + line.rstrip()
    if seq != '':
        Seqs.append(seq)
    infile.close()
    consensus = []
    con_Seq = ''
    con_Perc = []
    bPhy_PPs = [] # Here's a place where Patrick is evaluating quality - use this to attach library design.
    BFacts = open(D+'/{0}_from_BAli-Phy_Colors.txt'.format(node_name),'w')
    for i in range(len(Seqs[0])):
        consensus.append({})
        for aa in AAs:
            consensus[i][aa] = 0
        for seq in Seqs:
            consensus[i][seq[i]] += 1
        for aa in AAs:
            perc = consensus[i][aa]/len(Seqs)
            consensus[i][aa] = perc
        con_AA = max(consensus[i].items(), key=operator.itemgetter(1))[0]
        con_Seq = con_Seq + con_AA
        con_Perc.append(consensus[i][con_AA])
        if aa != '-': #I'm not sure if aa has a meaning as a variable here....
            BFacts.write(str(perc)+'\n')
    with open(D+'/{0}_Consensus.fasta'.format(node_name),'w') as outfile:
        outfile.write('>{0}_Sampled_Consensus\n{1}'.format(node_name,con_Seq))
    Seqs = {node_name:{'seq':con_Seq,'score':con_Perc}}
    with open(D+'/{0}_baliphy_PP.pickle'.format(node_name),'wb') as outfile:
        pickle.dump(Seqs[node_name],outfile)
    with open(D+'/{0}_baliphy_All_PPs.pickle'.format(node_name),'wb') as outfile:
        pickle.dump(consensus,outfile)
    return Seqs

def seq_heatmap(Seqs,fname,n_col=50):
    import matplotlib.pyplot as plt
    ### Seqs is a dictionary where labels are the keys 
    ### and the values are two dictionaries: seq and score
    ### seq is used for the annotation
    ### score is used for the color value
    ### Seqs = {node_name:{'seqs':sequence,'score':scores}}
    ### Sequences must be aligned
    ### Sequence must be a string
    ### Score must be a list of same length of Sequence
    
    Sequences = {}
    Scores = {}
    for label in Seqs.keys():
        seqName = label
        seq_length = len(Seqs[label]['seq'])
        Sequences[label] = []
        Scores[label] = []
        for i in range(len(Seqs[label]['seq'])):
            Sequences[label].append(Seqs[label]['seq'][i])
            Scores[label].append(Seqs[label]['score'][i])
    remainder = seq_length%n_col
    for i in range(n_col-remainder):
        for label in Seqs.keys():
            Sequences[label].append(' ')
            Scores[label].append(1)
    ### Setup DataFrames
    df_Scores = pd.DataFrame(Scores)
    df_Scores = df_Scores.T
    df_Seqs = pd.DataFrame(Sequences)
    df_Seqs = df_Seqs.T
    
    
    ### Create cmap
    # What is cmap???
    cmap = [(1, 0, 0), (1, 0, 0), (1, 0, 0),
         (1, 0, 0), (1, 0, 0), (1, 0, 0),
         (1, 0, 0), (1, 0, 0), (1, 0, 0),
         (1, 0, 0), (1, 0.2, 0), (1, 0.2, 0),
         (1, 0.4, 0), (1, 0.4, 0), (1, 0.6, 0),
         (1, 0.6, 0), (1, 0.8, 0), (1, 0.8, 0),
         (1, 0.9, 0), (1, 1, 0)]

    # fig, ax = plt.subplots(nrows=int(ceil(seq_length/n_col)+1),figsize=(10,int(ceil(seq_length/n_col)+1)))
    fig, ax = plt.subplots(nrows=int(ceil(seq_length/n_col)+1),figsize=(10,int(ceil(seq_length/n_col)+1)/3))
    for i in range(int(ceil(seq_length/n_col))):
        h1 = sns.heatmap(df_Scores[df_Scores.columns[i*n_col:(i+1)*n_col]], square=True,xticklabels=False,vmin=0,vmax=1,annot=df_Seqs[df_Seqs.columns[i*n_col:(i+1)*n_col]],fmt = '',ax=ax[i],cbar=False,cmap=cmap)
        ax[i].set_yticklabels(ax[i].get_yticklabels(),rotation=0)
    cmap_ticks = [0]
    for i in range(5):
        cmap_ticks.append(cmap_ticks[-1]+1/5)
    mappable = h1.get_children()[0]

    plt.colorbar(mappable, cax = ax[-1],ticks=cmap_ticks, orientation = 'horizontal')
    title = '{0} Consensus at Each Position'.format(seqName)
    ax[0].set_title(title,fontsize=14)
    plt.savefig(fname,format='pdf', dpi=1200, bbox_inches='tight')
    plt.close()

def simplify_trees():
    ReqClades = {}
    Clade_Failures = {}
    Keep_Taxa = []
    Rename_Taxa = []
    SampledTrees = {}
    for Clade_fasta in [x for x in os.listdir('Descendants') if x.split('.')[0].split('_')[1] == 'Clade']:
        with open('Descendants/' + Clade_fasta) as infile:
            Clade = Clade_fasta.split('_')[0]
            ReqClades[Clade] = [x.split()[0].lstrip('>') for x in infile.readlines() if x[0] == '>']
            Clade_Failures[Clade] = 0
            Keep_Taxa.append(ReqClades[Clade][0])
            Rename_Taxa.append(Clade)
    for tree_fname in [x for x in os.listdir('SampledTrees') if x.split('.')[-1] == 'tree']:
        with open('SampledTrees/'+tree_fname) as infile:
            trees = [Tree().get(data=tree_str,schema='newick') for tree_str in infile.readlines()]
        for tree in trees:
            keep = True
            for Clade in ReqClades:
                node_label = tree.mrca(taxon_labels=ReqClades[Clade]).label
                leaves = [x.taxon.label for x in tree.find_node_with_label(node_label).leaf_nodes()]
                new_taxa = set(leaves) - set(ReqClades[Clade])
                if len(new_taxa) > 0:
                    keep = False
                    Clade_Failures[Clade] += 1
            if keep:
                All_Taxa = set([x.taxon.label for x in tree.leaf_node_iter()])
                Remove_Taxa = list(All_Taxa - set(Keep_Taxa))
                tree.prune_taxa_with_labels(Remove_Taxa)
                for i in range(len(Keep_Taxa)):
                    node = tree.find_node_with_taxon_label(Keep_Taxa[i])
                    node.taxon.label = Rename_Taxa[i]
                nodes = tuple(sorted([tuple(sorted([l.taxon.label.split('|')[2].split('_')[0] if '|' in l.taxon.label else l.taxon.label for l in n.leaf_nodes()])) for n in tree.nodes() if not n.is_leaf()]))
                if nodes not in SampledTrees.keys():
                    SampledTrees[nodes] = []
                    SampledTrees[nodes].append(tree.as_string(schema="newick"))
                else:
                    SampledTrees[nodes].append(tree.as_string(schema="newick"))
        if not os.path.exists('SampledTrees/{}_Simplified'.format(tree_fname.split('.')[0])):
            os.mkdir('SampledTrees/{}_Simplified'.format(tree_fname.split('.')[0]))
        counter = 1
        SampledTrees = {k: v for k, v in sorted(SampledTrees.items(), key=lambda item: len(item[1]),reverse=True)}
        summary = {}
        with open('SampledTrees/{}_Simplified/Summary.txt'.format(tree_fname.split('.')[0]),'w') as outfile:
            for Clade in Clade_Failures.keys():
                outfile.write('{} Clade - {:.2%}/n'.format(Clade,1-(Clade_Failures[Clade]/len(trees))))
                pass
            outfile.write('\n')
            for key in SampledTrees.keys():
                perc = len(SampledTrees[key])/len(trees)
                t = Tree().get(data=SampledTrees[key][0],schema='newick').as_ascii_plot()
                outfile.write('Type_{} - {:.2%}\n{}\n\n'.format(counter,perc,t))
                with open('SampledTrees/e{}_Simplified/Type{}.trees'.format(sampling,counter,perc),'w') as outfile2:
                        outfile2.write(''.join(SampledTrees[key]))
                counter += 1
    pass

# Library design stuff
def parse_multifasta (f_path):
    Sequences_list = []
    with open(f_path) as handle:
        for Seq_Rec_Itr in SeqIO.parse(handle, "fasta"):
            Sequences_list.append(str(Seq_Rec_Itr.seq))
    return Sequences_list

def calculate_position_probability(Sequences_list):
    AA_Probability_Matrix = []
    # For each position,
    for position in range(len(Sequences_list[0])):
        # Make a blank list
        position_aa_distribution = [0]*22
        # And tally each Amino Acid
        for seq in range(len(Sequences_list)):
            AA = Sequences_list[seq][position]
            (position_aa_distribution[AA_key.index(AA)])+=1
        # Then add the collumn to the matrix
        AA_Probability_Matrix.append(position_aa_distribution)
    return AA_Probability_Matrix

def Record_Ancestor_Deviations (matrix, name, Num_Sequences): #All of this below is the handling for writing the human readable file.
    Percent_Share_Cuttoff=0.01 #The human readable file will inculde anything with >1% chance of being there
    space = "  "
    name += "_Non-Concensus_Deviations.txt"
    # Handling output directories
    if not os.path.isdir("Library Design/"):
        os.mkdir("Library Design/")
    if not os.path.isdir("Library Design/Ancestral Deviations/"):
        os.mkdir("Library Design/Ancestral Deviations/")
    while os.path.exists(f"Library Design/Ancestral Deviations/{name}_{Percent_Share_Cuttoff}_Cutoff.txt"):
        name+='_'
    with open(f"Library Design/Ancestral Deviations/{name}_{Percent_Share_Cuttoff}_Cutoff.txt","w+") as fout:
        for position in range(len(matrix)): # For every position in the sequence
            output = f"At position {space}{position+1}, "
            for i in range(len(matrix[position])): # For every AA probabilty
                if ((matrix[position][i] / Num_Sequences) >= Percent_Share_Cuttoff): # If the amino acid is above the threshold
                    percent=round((matrix[position][i]/Num_Sequences*100),1)
                    if (percent > 99.8):
                        output += f"{AA_key[i]} apears 100%  of the time;  "
                    elif (percent >= 10):
                        output += f"{AA_key[i]} apears {percent}% of the time;  "
                    else:
                        output += f"{AA_key[i]} apears {percent}%  of the time;  "
            output+="\n"
            fout.write(output)
            output = ""
            if (position==8):
                space=" "
            if (position==98):
                space=""

def Make_AA_List_for_Cutoff (matrix, Percent_Share_Cuttoff, Num_Sequences):
    Primer_Request = [] # This will be a list of positional requests. - Each item in the list is the list of amino acids needed at that position in the primer
    for position in range(len(matrix)): # For every position in the sequence
        positional_request=[]
        for i in range(len(matrix[position])): # For every AA probabilty
           
            if (i==21): # When evaluating gaps
                if ((matrix[position][i] / Num_Sequences) > 0.5): # If there is a >50% chance of a gap
                    positional_request=['-'] #Just record the gap
                else:
                    pass # Otherwise, don't record the gap.
           
            elif ((matrix[position][i] / Num_Sequences) >= Percent_Share_Cuttoff): # If the amino acid is above the threshold
                positional_request.append(AA_key[i]) # Record that AA at that position
        
        if not positional_request: # Empty lists evaluate as false
            # If none of the amino acids have a high enough probability to pass the threshold, we'll just record the most likely one.
            positional_request.append(AA_key[matrix[position].index(max(matrix[position]))]) 
# Possibly raise an error here? Let the user know
# Add amino acids until the threshold is reached?
# Is this the best position to be "spending" diversity? - AutoASR ancestors already have low confidence/high diversity libraries
# Epistasis could be useful here, but implemntation is beyond the scope and capabilities of this project and coder.
            # This mess will find the index of the maxium value in the probability list for the position, and add that AA to our request.
            # This solves a bug where a position in the primer request was blank.
        
        Primer_Request.append(positional_request) # Record list of AAs for current position
    return Primer_Request 

def code_template():    
    #THIS NEEDS TO BE MADE INTO A DATA FOR PRIMER DESIGN, NOT WRITING TO A FILE!!
    for Ancestor_Alignment in tqdm.tqdm([x for x in os.listdir('Ancestor_Nodes_Aligned') if x.split('.')[0].split('_')[-1] == 'Aligned']):
        # Each itteration of this loop is one file.
        Sequences_list = parse_multifasta(f"Ancestor_Nodes_Aligned/{Ancestor_Alignment}")
        # We now have the sequences parsed. Let's make a matrix.
        Num_Sequences = len(Sequences_list)
        AA_Probability_Matrix = calculate_position_probability(Sequences_list)
        # Now let's write the output.
        write_ancestor_deviations(AA_Probability_Matrix, Ancestor_Alignment)
    print("Concensus deviations have been determined.")

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
    # Tested on Dec 21st - no problems found yet

def Is_Valid_Codon(codon):
    dna_leters = ['a','c','t','g','r','y','m','k','s','w','h','b','v','d','n']
    if isinstance(codon,str):
        return (((len(codon)%3)==0) and all(i in dna_leters for i in codon))
    if isinstance(codon,list):
        return all(((len(c)==3) and all(i in dna_leters for i in c)) for c in codon)
    # Tested on Dec 21st - no problems found yet

def degenerate_codon(AA_List, source='Human'): # AA_List is a list of amino acids to make one degenerate codon for
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
    codon_list=[]
    if (source=='Human'):
        for AA in AA_List:
            codon_list.append(AA_to_Codon_Human[AA])
    elif (source=='EColi'):
        for AA in AA_List:
            codon_list.append(AA_to_Codon_Ecoli[AA])
    else:
        raise ValueError("Please Specify EColi or Human as an expression organism.")
    degenerate_codon = ''
    
    for pos in range(3):
        bases_at_pos=''
        i=0
        #print(f'length of codon list{len(codon_list)}')
        while (i<len(codon_list)) and ((codon_list[i])[pos] not in bases_at_pos) :
            #print(f"i:{i}")
            bases_at_pos+=codon_list[i][pos]
            i+=1
        degenerate_codon+=Mixed_Bases_lookup[bases_at_pos]
    #print(degenerate_codon)
    return degenerate_codon

def Build_Primer (Primer_Request, source='Human'):
    primer_txt = ''
    for pos in Primer_Request: # pos is the list of AA at the respective position
        if (len(pos)==1):
            if (pos[0]!='-'): # For a gap in the sequence, we will record nothing.
                if (source=='Human'):
                    primer_txt+=(AA_to_Codon_Human[pos[0]])
                elif (source=='EColi'):
                    primer_txt+=(AA_to_Codon_Ecoli[pos[0]])
                else:
                    raise ValueError("Please Specify EColi or Human as an expression organism.")
        elif (len(pos)==2):
            primer_txt+=(AA_Pair_lookup[ pos[0]+pos[1] ])
        else:
            primer_txt+=(degenerate_codon(pos))
    if not Is_Valid_Codon(primer_txt):
        raise ValueError ("The primer generated was not a valid DNA sequence. No clue how that happened. If you're seeing this error, it's proabbly caused by a bug.")
    return (primer_txt)

def Build_Libraries(Multifasta,Cutoffs=[0.9],Name=''): # A Multifasta must be provided. The list of cutoff thresholds has a default 90%
    
    if (Name==''):
        print(f"No name provided. Default name: {Multifasta[23:-38]}")
        Name=Multifasta[23:-38] # The default name is the begining of the Multifasta File
    
    Sequences_list = parse_multifasta(Multifasta)
    AA_Probability_Matrix=calculate_position_probability(Sequences_list)
    Num_Sequences = len(Sequences_list)
    Record_Ancestor_Deviations(AA_Probability_Matrix,Name,Num_Sequences)
    
    for Threshold in Cutoffs:
        
        Primer_Request = Make_AA_List_for_Cutoff(AA_Probability_Matrix,Threshold,Num_Sequences)

        # Debugging
        with open(f"Library Design/Debugging/{Name} {Threshold}.txt",'w+') as fout:
            fout_string=''
            for item in Primer_Request:
                fout_string+=(" ".join(item)+"\n")
            fout.write(fout_string)

        if not Is_Valid_AA(Primer_Request):
            raise ValueError ("The primer request is not made of valid amino acids. If you're seeing this error, it's proabbly caused by a bug.")
        Primer_Text = Build_Primer(Primer_Request)
        
        #Now we need to record that primer we've created
        if not os.path.isdir("Library Design/Primers/"):
            os.mkdir("Library Design/Primers/")
        while os.path.exists(f"Library Design/Primers/{Name}_{Threshold*100}%_Cutoff.fasta"):
            Name+='_'
        with open(f"Library Design/Primers/{Name}_{Threshold*100}%_Cutoff.fasta",'w+') as fout:
            fout.write(f">DNA sequence for the library of {Name} with a {Threshold*100}% Cutoff:\n")
            fout.write(Primer_Text)

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

    'act':'h', 'atc':'h', 'cat':'h', 'cta':'h', 'tac':'h',' tca':'h',
    'gct':'b', 'gtc':'b', 'ctg':'b', 'cgt':'b', 'tgc':'b', 'tcg':'b',
    'acg':'v', 'agc':'v', 'cag':'v', 'cga':'v', 'gca':'v', 'gac':'v',
    'agt':'d', 'atg':'d', 'gat':'d', 'gta':'d', 'tga':'d', 'tag':'d',

    'acgt':'n','actg':'n','agct':'n','agtc':'n','atcg':'n','atgc':'n',
    'cagt':'n','catg':'n','gact':'n','gatc':'n','tacg':'n','tagc':'n',
    'cgat':'n','ctag':'n','gcat':'n','gtac':'n','tcag':'n','tgac':'n',
    'cgta':'n','ctga':'n','gcta':'n','gtca':'n','tcga':'n','tgca':'n'
    }
# This dictionary provides the best degenerate codon for every combination of two amino acids
AA_Pair_lookup = { 
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
# This dictionary provides a key for the amino acids to be used by probabilty functions. (Alphabetical)
AA_key = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y','X','-'] # NO TOUCH

if __name__ == '__main__': 
    cutoffs=[0.5,0.33,0.25,0.2,0.1,0.05] # This will only record amino acids with X probablity of being at that site or higher
    # with a 20% cutoff, the amino acid must have a 20% chance of being at that site to be included in the primer request
    for Ancestor_Alignment in tqdm.tqdm([x for x in os.listdir('Ancestor_Nodes_Aligned') if x.split('.')[0].split('_')[-1] == 'Aligned']):
        Build_Libraries((f"Ancestor_Nodes_Aligned/{Ancestor_Alignment}"),cutoffs)