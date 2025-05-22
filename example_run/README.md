# Example Run Guide

Here we briefly describe an example usage of ExTRaCT. For more details go to ExTRaCT/example_Z_domain_bat_all_C.sh.

Activate the environment with all required dependencies.

```#!/bin/sh

conda activate a3pipeline

# update paths
export PATH=/$path_to/emboss.open-bio.org/pub/EMBOSS/EMBOSS-6.6.0/emboss/:$PATH
export PATH=/$path_to/standard-RAxML-master:$PATH
export PATH=/$path_to/scipio-1.4/:$PATH
export PATH=/$path_to/blat/blat/:$PATH
export PATH=/$path_to/bin/:$PATH

```

To run the following example, load a sample (Myotis myotis) genome from GenBank # GCA_014108235.1. 
The following will extract Z-domains based on 
* bat reference sequences
* all (single and double) domain A3
* C motif structure

```#!/bin/sh

cd ExTRaCT

# prepare hmmer profiles
# update mafft path
python prep_hmm.py hmmer_profiles_files/6bat_nucl_all.fasta hmmer_profiles_files/6bat_all_msa.fasta hmmer_profiles_files/6bat_all.sto hmmer_profiles_files/6bat_all.hmm /$path_to/mafft 

# run main algorithm
python Gene_search.py species_names.csv hmmer_profiles_files/6bat_all.hmm example_run/bat/all Zdomains_C.csv bat_all_C

# convert results table to fasta
python table_to_fasta.py example_run/bat/all bat_all_C_final_data_table.csv bat_all_C_final_data_table.fasta Zdomain_AA

# build phylogenetic tree
# update mafft path
python Gene_tree.py /$path_to/mafft example_run/bat/all/fasta bat_all_C_final_data_table.fasta

# update start and end locations and pull updated nucleotide sequences
# copy and paste this code if running batches in parallel and change "run1"
python scipio_run.py example_run/bat/all/bat_all_C_final_data_table.csv example_run/bat/all/scipio/ run1/

# Combine all csv output to generate meta table
awk -F, 'FNR > 1' example_run/bat/all/scipio/run1/*output_sci.csv > example_run/bat/all/scipio/run1_all_scipio.csv
#awk -F, 'FNR > 1' example_run/bat/all/scipio/*/*output_sci.csv > example_run/bat/all/scipio/all_scipio.csv ## for multiple run folders

# Create file with headers
cat example_run/bat/all/scipio/run1/*output_sci.csv | head -n 1 > example_run/bat/all/scipio/headers.csv
cat example_run/bat/all/scipio/headers.csv example_run/bat/all/scipio/run1_all_scipio.csv > example_run/bat/all/scipio/run1_all_scipio_wheader.csv

```

To generate all 27 test cases, copy and paste the code above and make changes to match the species, domain type and motif structure file names. 
Suggested structure: Reference-species_Domain-type_Motif
Where the reference species are bat, laurasiatheria or primate; the domain types are all, single or double; and the motif options are A, B or C.

See example_Z_domain_bat_all_C.sh for more help.


