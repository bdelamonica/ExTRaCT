#!/usr.bin/env python 

# Extract accurate gene location : find Z domain location from Z domain hmmer profile using scipio

# input : table with headers [chrom,chromStart,chromEnd,name,score,strand,species,fasta,Ztype,Zcount,namesearch,Zdomain_AA] from Gene_search.py
# output : new table with non-empty Z-domains and updated start and stop locations
# usage : python scipio_run.py INPUT_TABLE SCIPIO_PATH SUBFOLDER
# example : python scipio_run.py example_run/bat/all/bat_all_C_final_data_table.csv example_run/bat/all/scipio/ run1/
# SCIPIO may skip non-overlapping repeated sequences, more than 1 YAML entry will be returned
# To get more than one location for the same sequence go to the /scipio/ folder 

import sys
import os
import re
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pybedtools
import argparse

def unzip_fasta(fasta):
    """Unzips the gzipped FASTA file if necessary."""
    if fasta.endswith('.gz'):
        fastafile = re.sub(".gz", "", fasta)  # Remove .gz extension for the unzipped file
        if os.path.exists(fasta):
            print("Unzipping: " + fasta)
            os.system("gunzip " + fasta)
        else:
            print("File " + fasta + " not found. Skipping unzipping.")
    else:
        fastafile = fasta
    return fastafile

def process_fasta(fasta, filtable, savedin, species_name):
    """Process FASTA file and generate output."""

    #print("Species: " ,species_name)
    
    outfile_nc = open(savedin + species_name + "_bat_Z_domain_NC.fasta", "w")
    for _, row in filtable.iterrows():
        #print(filtable.columns)
        
        rname = row[14]
        #print("Temp seq ID: ", row[6] + "_" + row[8] + "_" + str(row[9]))
        #print(rname)

        # Create temp sequence
        temp_seq = SeqRecord(
            Seq(row[13]),
            id=row[7] + "_" + row[10] + "_" + str(row[11]),
            name=rname,
            description=row[12]
        )
        print(temp_seq)
        
        # Save sequence
        tempfasta = savedin + rname + ".fasta"
        SeqIO.write(temp_seq, tempfasta, "fasta")

        # Run Scipio pipeline
        run_scipio_pipeline(tempfasta, savedin, rname, fasta)
        
        # Process the result and update the DataFrame
        update_dataframe_with_scipio_output(filtable, _, savedin, rname, fasta, temp_seq, outfile_nc)
        
    outfile_nc.close()
    
    return filtable

def run_scipio_pipeline(tempfasta, savedin, rname, fastafile):
    """Run Scipio pipeline and convert results to BED format."""
    pslf = savedin + rname + ".psl"
    outf = savedin + rname + ".out"
    scigff = savedin + rname + "_sci.gff"
    scibed = savedin + rname + "_sci.bed"
    #print(scibed)
    os.system("perl /titan/brenda/scipio-1.4/scipio.1.4.1.pl --blat_output="+ pslf + " " + fastafile + " " + tempfasta + " > " + outf )
    os.system("perl /titan/brenda/scipio-1.4/yaml2gff.1.4.pl " + outf + " > " + scigff)
    os.system("/titan/brenda/bin/convert2bed --input=gff < " + scigff + " > " + scibed)


def update_dataframe_with_scipio_output(filtable, idx, savedin, rname, fasta, temp_seq, outfile_nc):
    """Update the DataFrame with results from Scipio."""
    #print("Looking for: " + temp_seq.id)
    scibed = savedin + rname + "_sci.bed"
    try:
        scibedfile = pd.read_csv(scibed, sep="\t", header=None)
        #print(scibedfile)

        # Check if the file is empty
        if not scibedfile.empty:
            chrom = filtable.loc[idx, 'chrom'] = scibedfile.iloc[0, 0]
            start = filtable.loc[idx, 'chromStart'] = scibedfile.iloc[0, 1]
            end = filtable.loc[idx, 'chromEnd'] = scibedfile.iloc[0, 2]
            strand = filtable.loc[idx, 'strand'] = scibedfile.iloc[0, 5]
                #filtable.loc[idx, 12]= str(seq_rd.seq)
                #seq_rd.description = seq_rd.id
                #seq_rd.id = temp_seq.id
                #SeqIO.write(seq_rd, outfile_nc, "fasta")
        else:
            print(f"File {scibed} is empty. Skipping processing for {temp_seq.id}.")
    except Exception as e:
        print(f"Error reading {scibed}: {e}")

def main():
    # Parse arguments
    pd.options.mode.chained_assignment = None  # default='warn'
    parser = argparse.ArgumentParser()
    parser.add_argument("input_table", help="Input csv file with format: [chrom, chromStart, chromEnd, name, score, strand, species, fasta, Ztypelist , Ztype, Zcount, namesearch, Zdomain_AA] required")
    parser.add_argument("path", help="Output Path required")
    parser.add_argument("prefix", help="Output prefix required")
    args = parser.parse_args()

    # Load input data
    itable = pd.read_csv(args.input_table)
    path_Str = str(args.path)
    os.makedirs(path_Str, exist_ok=True)
    prefix_Str = str(args.prefix)
    savedin = path_Str + prefix_Str
    os.makedirs(savedin, exist_ok=True)

    header = ['chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand','PBT_NC_sequence','species','fasta','Ztypelist','Ztype','Zcount','namesearch','Zdomain_AA'] ### ensure headers match
    itable.columns = header[:len(itable.columns)]
    
    # Filter out Zdomain_AA = x rows and save for reference
    df = itable.loc[itable['Zdomain_AA'] != 'x']
    df = df.loc[df['Zdomain_AA'] != 'Partial']
    df['rowname'] = df['species'] + "_" + df['Ztype'] + "_" + df['Zcount'].astype(str) + "_" + df['chrom'].astype(str) + "_" + df['chromStart'].astype(str) + "-" + df['chromEnd'].astype(str)
    df['Sci_Zdomain_NC'] = 'x'
    df.to_csv(savedin + "input_sci.csv", index=False)
    
    #print(df)
    #print(prefix_Str)

    # Identify unique fasta files
    fasta_files = df['fasta'].unique()

    # Loop through each fasta file and process it
    for fasta in fasta_files:
        fastafile = unzip_fasta(fasta)
        #print(fastafile)

        # Get a subset of the table for the current FASTA file
        filtable = df[df['fasta'] == fasta]
         
        species_name = filtable['species'].unique()[0]
        newtable = process_fasta(fastafile, filtable, savedin,species_name)
        print("Species: " ,species_name)
        print("ORIGINAL TABLE: \n", filtable)
        print("Processed: \n", newtable)

        # First, identify rows that have unique combinations of 'chrom' and 'chromStart'
        unique_rows = newtable[~newtable.duplicated(subset=['chrom', 'chromStart', 'chromEnd', 'strand'], keep='first')]
        print("filtered: \n" ,unique_rows)
        newtable.to_csv(savedin + species_name + "_output_sci.csv", index=False)
        unique_rows.to_csv(savedin + species_name + "_output_sci_filtered.csv", index=False)
        
        ### Re-gzip the fastafile if required
        #print("Zip fasta")
        #os.system(f"gzip {fastafile}")
        

if __name__ == "__main__":
    main()



