#These are unused functions, in storage.
# Written by James VanAntwerp in September 2021 - vanant25@msu.edu
# Written by Pattrick Finneran, Menten AI, Palo Alto, California, United States of America
# Written for the Woldring Lab, Michigan State University in East Lansing, Michigan, USA.

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
