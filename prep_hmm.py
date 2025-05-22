#!/usr.bin/env python

# Code to set up hmmer profiles from input fasta
# Converts fasta to msa to stockholm file to hmmer profile
# Input file : Nucleotide fasta file that contains known sequences

# Make sure to install hmmer (see master README file):
# conda install -c biocore hmmer
# Install StringIO
# Install AlignIO

# Check compatability issues with fasta file using:
# sed '/^>/! s/N/-/g; /^>/! s/W/-/g; /^>/! s/K/-/g ; /^>/! s/Y/-/g ; /^>/! s/S/-/g ; /^>/! s/R/-/g; /^>/! s/-//g; s/|/_/g; s/\[//g; s/\]//g; s/\.\./-/g; s/->/-/g;s/\//_/g; /^>/ s/ /_/g' INPUT.fasta | sed '/^[[:space:]]*$/d'
# If there are differences save to a new fasta and run

from io import StringIO
from Bio import AlignIO
import os
import sys
import argparse
import pandas as pd

pd.options.mode.chained_assignment = None  # default='warn'

# Suggested input file name structure: REFERENCE-SPECIES_DOMAIN-TYPE
# prep_hmm.py usage:
# pythin prep_hmm.py REFERENCE-SPECIES_DOMAIN-TYPE.fasta REFERENCE-SPECIES_DOMAIN-TYPE.msa REFERENCE-SPECIES_DOMAIN-TYPE.sto REFERENCE-SPECIES_DOMAIN-TYPE.hmm
# python prep_hmm.py hmmer_profiles_files/6bat_nucl_all.fasta hmmer_profiles_files/6bat_all_msa.fasta hmmer_profiles_files/6bat_all.sto hmmer_profiles_files/primates_single.hmm /titan/brenda/mafft

parser = argparse.ArgumentParser()
parser.add_argument("seq_input_file", help="Input fasta file with known A3s required")
parser.add_argument("msa_file_name", help="MSA file name required")
parser.add_argument("sto_file_name", help="Stolkholm file name required")
parser.add_argument("hmmer_file_name", help="HMMER file name required")
parser.add_argument("mafft_loc",help="MAFFT location required")
args = parser.parse_args()

seq_input = str(args.seq_input_file)
msa_fasta = str(args.msa_file_name)
sto_file = str(args.sto_file_name)
hmmer_file = str(args.hmmer_file_name)
mafft_path = str(args.mafft_loc)

# use bash os to call MAFFT to align sequences
# if using bash directly the command line code is:
# ~/msafiles/mafft-linux64/mafft.bat 6bat_nucl_all.fasta > 6bat_nucl_all_msa.fasta
# ie
# [MAFFT location] seq_input > [aligned sequence output file name]
# change seq_input > [aligned sequence output file name] to generate MSA file
os.system(mafft_path +" "+ seq_input +" > " + msa_fasta)
#os.system(f"{mafft_path} {seq_input} > {msa_fasta}")

# converts aligned file to stockholm file for hmmer set up
# change [aligned sequence output file name],...,[name.sto]
AlignIO.convert(msa_fasta,"fasta",sto_file,"stockholm","DNA")

# use bash os to call hmmer to create profile
# if using bash directly the command line code is:
# hmmbuild output.hmm input.sto
# HMMER docs: http://hmmer.org/documentation.html
# http://eddylab.org/software/hmmer/Userguide.pdf 
# change [hmm_profile name.sto]
os.system(f"hmmbuild {hmmer_file} {sto_file}")
# output is hmm_profile

quit()
