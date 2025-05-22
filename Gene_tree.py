#!/usr/bin/env python

# Code to align sequences from Gene_search.py (and table_to_fasta.py) and generate phylogenetic gene tree
# usage: python Gene_tree.py /PATH_TO/MAFFT FASTA_FOLDER INPUT.fasta
# example: python Gene_tree.py [path]/mafft example_run/bat/all/fasta bat_all_C_final_data_table.fasta

import os
import sys
import pandas as pd
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("mafft_loc",help="MAFFT location required")
parser.add_argument("input_path", help="Input path required")
parser.add_argument("input_file",help="Input file required")

args = parser.parse_args()

mafft_path = str(args.mafft_loc)
in_path = str(args.input_path)
in_file = str(args.input_file)
tree_path = in_path.replace("fasta", "tree")
print(tree_path)
msa_path = tree_path + "/" + in_file.replace(".fasta","_msa.fasta")
print(msa_path)
output = tree_path + "/" + in_file.replace(".fasta","_tree")
print(output)
os.makedirs(tree_path, exist_ok=True)
    
# use bash os to call MAFFT to align sequences
# command line: mafft --auto input > output
os.system(mafft_path +" --auto "+ in_path +"/"+ in_file+" > " + msa_path)

os.chdir(tree_path)
os.getcwd()

# Check if the file with the '.reduced' extension exists
#if os.path.exists(msa_path + ".reduced"):
#    msa_path = msa_path + ".reduced"
#else:
#    msa_path = msa_path

msa_file = msa_path.replace(tree_path + "/","")

print("Updated path" , msa_path, "\n file", msa_file)

out_file = output.replace(tree_path + "/","")
print(out_file)

# use bash os to call RAxML and save file
os.system("raxmlHPC -s " + msa_file + " -p 12345 -m PROTGAMMAAUTO --no-seq-check --print-identical-sequences -n " + out_file)
