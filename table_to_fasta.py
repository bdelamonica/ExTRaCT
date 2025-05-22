#!/usr/bin/env python

# Takes Z domain search table and converts sequence lines to fasta file 
# required input: table, output file path/name, sequence column
# usage: python table_to_fasta.py INFILE OUTFILE SEQCOLUMN
# example : python table_to_fasta.py example_run/bat/all bat_all_C_final_data_table.csv bat_all_C_final_data_table.fasta Zdomain_AA

import pandas as pd
import argparse
from Bio import SeqIO
import os
 
def create_fasta(input_path, input_file, output_file, sequence_column):
    # Read the CSV file into a pandas DataFrame
    df = pd.read_csv(input_path + "/" + input_file,header=None)
    header = ['chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand','PBT_NC_sequence','species','fasta','Ztypelist','Ztype','Zcount','namesearch','Zdomain_AA'] ### order here matters, if changes are made to Gene_search.py any changes to the final table columns must be reflected here
    df.columns = header[:len(df.columns)]
    print(df.head())
    print(df.shape)
    df = df.loc[df['Ztype'] != 'x']
    df = df.loc[df['Ztype'] != 'Partial']
    print(df.shape)
    #df = df.drop_duplicates(subset=['species','Zdomain_AA']) # un comment if you do not need duplicates, not recommended since HMMER finds non-overlapping sequences and duplicate sequences may be important for evolutionary analysis
    print(df.shape)
    print(df.head())
    df.to_csv(input_path + "/" + input_file.replace(".csv","_filtered.csv"), index=False)
    # Open the output FASTA file for writing
    with open(output_file, 'w') as fasta_out:
        # Loop through each row in the DataFrame
        for index, row in df.iterrows():
            # Extract the necessary fields for metadata and sequence
            #### If file has different column names make changes here too ####
            species = row['species']
            Ztype = row['Ztype']
            Zcount = row['Zcount']
            chrom = row['chrom']
            chromStart = row['chromStart']
            chromEnd = row['chromEnd']
            sequence = row[sequence_column]  # Use the specified sequence column from input
            
            # Construct the metadata fields to save the fasta sequences
            sequence_id = f"{species}_{Ztype}_{Zcount}"
            sequence_name = f"{species}_{Ztype}_{Zcount}_{chrom}_{chromStart}-{chromEnd}"
            sequence_desc = f"{chrom}_{chromStart}-{chromEnd}"
            
            # Write the sequences to the FASTA file
            fasta_out.write(f">{sequence_name} {sequence_desc}\n")
            fasta_out.write(f"{sequence}\n")

def main():
    # Set up command-line argument parsing
    parser = argparse.ArgumentParser(description="Generate a FASTA file from a CSV input table")
    
    # Define the arguments
    parser.add_argument("input_path", help="Path to the input CSV file")
    parser.add_argument("input_file", help="Input CSV file name")
    parser.add_argument("output_file", help="Output FASTA file name")
    parser.add_argument("sequence_column", help="Name of the sequence column (e.g., 'Zdomain_AA' or 'Zdomain_NC')")
    
    # Parse the arguments
    args = parser.parse_args()
    in_path = str(args.input_path)
    in_file = str(args.input_file)
    
    out_file = str(args.output_file)
    seq_col = str(args.sequence_column)
    
    # Ensure the output directory exists
    out_path = in_path + "/fasta"
    print(out_path)
    os.makedirs(out_path, exist_ok=True)
    
    # Call the function to create the FASTA file
    print(out_path +"/"+ out_file)
    create_fasta(in_path ,in_file, out_path + "/" + out_file, seq_col)


if __name__ == "__main__":
    main()

