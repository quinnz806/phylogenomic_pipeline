#script for my pipeline


# Import needed modules
import os
import sys
import glob
import subprocess
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import Phylo

# Sore path to directory containing input files

in_dir = "/shared/forsythe/BB485/Week06/Brass_CDS_seqs/"
print(in_dir)

#store a path to new folder where to put output

out_dir = "/home/quinnz/novus/BB485/Week06/phy_pipeline_out/"
print(out_dir)

# Get a list of all fasts files in that folder

all_in_files=glob.glob(in_dir+"*fasta")
#print(all_in_files)


# Start a for loop for each file name, indie of loop:
    # create a mafft command for the file of intrest 
    #call mafft ftrom the command line using system call

for file in all_in_files[0:10]:
    print(file)
    new_file_path = file.replace(in_dir, out_dir)
    #print(new_file_path)
    
    # Create a command string (this is what get called using the 'system call'.
    aln_cmd = 'mafft --auto --quiet '+file+' > '+new_file_path
    print(aln_cmd)

    #Run the command using a 'system call'
    os.system(aln_cmd) #uncomment once you've check the command

# create a for loop to sun iqtree
#get a list of all fasfta files in the input folder 
all_aln_files=glob.glob(out_dir+"*fasta")

#loop through algn files and run iqtree on each

#empty list, append values to it


for file in all_aln_files: 
    #Read in the tree and store as phylo object
    temp_tree = Phylo.read(tree, "newick")

#Loop through the tips in the tree to find which one contains Es (the outgroup)
    for tip in temp_tree.get_terminals():
        if "Es_" in tip.name:
            es_tip = tip
		#Stope the loop once we found the correct tip
            break
    
#Root the tree by the outgroup taxon
    temp_tree.root_with_outgroup(es_tip)
    
#Get a list of all terminal (aka tips) branches
    all_terminal_branches = temp_tree.get_terminals()
    
#Loop through the branches and store the names of the tips of each
    for t in all_terminal_branches:
        if "Bs_" in t.name:
            Bs_temp=t 
        elif "Cr_" in t.name:
            Cr_temp=t
        elif "At_" in t.name:
            At_temp=t
        else:
            out_temp=t
        
#Make lists of pairs of branches, so that we can ask which is monophyletic
    P1_and_P2=[Bs_temp, Cr_temp]
    P1_and_P3=[Bs_temp, At_temp]
    P2_and_P3=[Cr_temp, At_temp]
    

#Use series of if/else statements to ask which pair in monophyletic
    if bool(temp_tree.is_monophyletic(P1_and_P2)):
        topo_str = "12top"
    elif bool(temp_tree.is_monophyletic(P1_and_P3)):
        topo_str = "13top"
    elif bool(temp_tree.is_monophyletic(P2_and_P3)):
        topo_str = "23top"
    else:
        topo_str = "Unknown"

    print(topo_str)



#create iq tree for each 
#once working, create a job and run it in that 
#run on 24 threads