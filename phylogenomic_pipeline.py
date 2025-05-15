#script for my pipeline


# Import needed modules
import os
import sys
import glob
import subprocess
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import Phylo
import matplotlib.pyplot as plt

# Sore path to directory containing input files

in_dir = "/shared/forsythe/BB485/Week06/Brass_CDS_seqs/"
print(in_dir)

#store a path to new folder where to put output

out_dir = "/home/quinnz/novus/BB485/Week06/final_phylo_out/"
print(out_dir)

# Get a list of all fasts files in that folder

all_in_files=glob.glob(in_dir+"*fasta")
#print(all_in_files)


# Start a for loop for each file name, indie of loop:
    # create a mafft command for the file of intrest 
    #call mafft ftrom the command line using system call

for file in all_in_files:
    #print(file)
    new_file_path = file.replace(in_dir, out_dir)
    #print(new_file_path)
    
    # Create a command string (this is what get called using the 'system call'.
    aln_cmd = 'mafft --auto --quiet '+file+' > '+new_file_path
    #print(aln_cmd)

    #Run the command using a 'system call'
    os.system(aln_cmd) #uncomment once you've check the command

# create a for loop to sun iqtree
#get a list of all fasfta files in the input folder 
all_aln_files=glob.glob(out_dir+"*fasta")

#print(all_aln_files)

#empty list, append values to it
topology_list= []

#loop through algn files and run iqtree on each
for file in all_aln_files:

    new_file_path = file.replace(in_dir, out_dir)
    #Create the command. -nt 2 means two threads. If running this from within a job submission, you could use more threads to make it go faster.
    tree_command = f"iqtree -s {new_file_path} -m TEST -nt 16"

    #Check the command 
    #print(tree_command)

#Run the command using a 'system call'
    os.system(tree_command) #uncomment once you've check the command


for file in all_aln_files:

    treefile = file + ".treefile"

    temp_tree = Phylo.read(treefile, "newick")

#Loop through the tips in the tree to find which one contains Es (the outgroup)
    for tip in temp_tree.get_terminals():
        if "Es_" in tip.name:
            es_tip = tip 

            #stop once we found correct tip

            break
    #root tree by the outgroup 
    temp_tree.root_with_outgroup(es_tip)

    #get a list of terminal (aka tips) branches 
    all_terminal_branches = temp_tree.get_terminals()

    #loop through branches and store the names of the tips of each

    Bs_temp = Cr_temp = At_temp = None

    for t in all_terminal_branches:
        if "Bs_" in t.name:
            Bs_temp=t
        elif "Cr_" in t.name:
            Cr_temp=t
        elif "At_" in t.name:
            At_temp=t
        else:
            out_temp=t

    if not all([Bs_temp, Cr_temp, At_temp]):
        print(f"ERROR: Could not find all expected tips")
        continue
    

    #Make lists of pairs of branches, so that we can ask which is monophyletic
    P1_and_P2=[Bs_temp, Cr_temp]
    P1_and_P3=[Bs_temp, At_temp]
    P2_and_P3=[Cr_temp, At_temp]


    #use a series of if/else statements to ask which pair is monophyletic 
    if bool(temp_tree.is_monophyletic(P1_and_P2)):
        topo_str = "12top"
    elif bool(temp_tree.is_monophyletic(P1_and_P3)):
        topo_str = "13top"
    elif bool(temp_tree.is_monophyletic(P2_and_P3)):
        topo_str = "23top"
    else:
        topo_str= "unknown"
        
    print(topo_str)

    topology_list.append((os.path.basename(treefile), topo_str))

    print(topology_list)

    #using a counter to find topology counts
    # Create an empty dictionary to hold counts
    topo_counts = {}

# Loop through each topology result
    for tree_name, topo in topology_list:
        if topo not in topo_counts:
            topo_counts[topo] = 1
        else:
            topo_counts[topo] += 1

# Print the results
    print("\nTopology Summary:")
    for topo in topo_counts:
        print(f"{topo}: {topo_counts[topo]}")


    
# Data for the pie chart
labels = list(topo_counts.keys())
sizes = list(topo_counts.values())

# Optional: Colors and explode settings
colors = ['#66c2a5', '#fc8d62', '#8da0cb', '#e78ac3']
explode = [0.05] * len(labels)  # Slightly separate each slice

# Create the pie chart
plt.figure(figsize=(6,6))
plt.pie(sizes, labels=labels, autopct='%1.1f%%', startangle=140,
        colors=colors[:len(labels)], explode=explode, shadow=True)
plt.title("Topology Distribution")
plt.axis('equal')  # Makes the pie chart round
plt.tight_layout()

# Show the plot
plt.savefig("topology_pie.png")


