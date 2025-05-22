#!/usr/bin/env python
# apobec 3: find Z domains from hmmer profile
# input = species table with genome fasta file names
# output = hmm hits and Z domain as fasta
# usage = python Gene_search.py species-input-table.csv reference-species_domain-type.hmm output-folder/reference-species/domain-type motif-file.csv file_start_name
# example run = python Gene_search.py species_names.csv 6bat_all.hmm example_run/bat/all motif_C.csv bat_all_C

import sys
import os
import re
from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd
import pybedtools
import argparse
#import gzip
from datetime import datetime

pd.options.mode.chained_assignment = None  # default='warn'

# Function to check if genome file must be gunzipped
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
        fastafile = fasta # Otherwise, use genome file path and name as is
    return fastafile

# Function to process the HMMER section
def run_hmmer(fastafile, hmmprofile, hmmerfile):
    # Run HMMER to find matches
    hmmercmd = "nhmmer --tblout " + hmmerfile + " " + hmmprofile + " " + fastafile
    os.system(hmmercmd)

    # Edit HMMER output to create BED file
    hmmerbed = hmmerfile + ".bed"
    editcmd1 = "cat " + hmmerfile + " | awk '($1!~\"#\"){print $1,$9,$10,$3,$13,$12}' | sed 's/ /\t/g' > " + hmmerbed
    os.system(editcmd1)

    # Read the bed file into a DataFrame
    bedfile = pd.read_csv(hmmerbed, sep='\t', header=None)

    # Fix bed file if necessary (adjust start and end value order)
    for hit in bedfile.index:
        if bedfile.iloc[hit][5] == "-":
            start = bedfile.iloc[hit][2]
            stop = bedfile.iloc[hit][1]
            bedfile.at[hit, 1] = start
            bedfile.at[hit, 2] = stop

    # Sort bedfile by chrom, start, and end
    temp_bed = bedfile.sort_values(by=[0, 1, 2])
    test_df = temp_bed.reset_index(drop=True)
    print("Run HMMER \n" ) #, test_df)
    test_df.to_csv(hmmerbed,sep=',',header=False, index = False)
    os.system(f"rm {hmmerfile}")
    return hmmerbed

# Function to process the PyBedTools section
# create NC (nucleotide) files required here for getorf
def run_pybedtools(hmmerbed, fastafile, hmmer_nc_fasta,pybt_tbl): # here pybt = pybedtools
    # Read the bed file into a DataFrame
    bedfile = pd.read_csv(hmmerbed, sep=',', header=None) 

    # Convert bedfile to BedTool object
    apobec_bed = pybedtools.BedTool.from_dataframe(bedfile)

    # Get the sequences from the genome file
    a = apobec_bed.sequence(fi=fastafile, s=True)

    # Save sequences here for table 
    sequences = []

    # Run pybedtools command and save NC sequences to FASTA
    with open(hmmer_nc_fasta, "w") as outfile_nc: 
        for seq_record in SeqIO.parse(open(a.seqfn), "fasta"):
            sequences.append(str(seq_record.seq))
            SeqIO.write(seq_record, outfile_nc, "fasta")

    # Add the nucleotide sequences as a new column to the original bedfile
    bedfile['PBT_NC_sequence'] = sequences
    print("Run PyBedTools \n") #, bedfile)
    # Save the updated bedfile with the new sequence column
    
    bedfile.to_csv(pybt_tbl, sep=',', header=False, index=False)

    return pybt_tbl


# Function to process the GetORF section
# pass NC files from pybedtools required here for GetORF
def run_getorf(hmmer_nc,hmmer_nc_2,orffile):
    # Prepare and run GetORF
    os.system(f"cat {hmmer_nc} | sed 's/:/_/' > {hmmer_nc_2}") # if hmmer_nc_2 exists, then this function will be skipped below
    os.system(f"getorf {hmmer_nc_2} -outseq {orffile} -minsize 234 -reverse no")
    # Clean up by removing the temporary files
    os.system(f"rm {hmmer_nc}") # also remove other hmmer_nc files to save space, can comment out for inspection
    print("Run GetORF \n")
    return orffile 

# Function to save only Z domains to a fasta file
def save_z_domains(orffile, test_df, Zdomain_folder, speciesname,Zdomain_aa_file_name, Zdomain_partial_aa_file_name, Zdomain_csv):
    ## Load the Z domain patterns and types from the CSV file 
    Zdomain_data = pd.read_csv(Zdomain_csv)
    Zdomain_patterns = {row['Zdomain_structure']: row['Ztype'] for _, row in Zdomain_data.iterrows()}

    outfile = open(Zdomain_aa_file_name,"w") #output fasta with final Z domains
    outfile_partial= open(Zdomain_partial_aa_file_name,"w") #output fasta with partial Z domains
    Zcounter=0
    #Iterate over sequences in orffile from function run_getorf 
    for orf_record in SeqIO.parse(orffile, "fasta"):
        #print("ORF",Zcounter,"\n",orf_record.description, orf_record.seq)
        ival = test_df[test_df['namesearch'] == orf_record.description.split('(')[0]].index[0]

        # Initialize a list to store the Ztypes for a given sequence (in case there are multiple matches)
        matched_ztypes = []
        matched_domains = [] 

        ## Case when multiple matches are expected

        # Check if the sequence matches any Z domain pattern in Zdomain_patterns
        for pattern, Ztype in Zdomain_patterns.items():
            # Find all matches for the current pattern
            matches = re.findall(pattern, str(orf_record.seq))

            for match in matches:
                print("for loop match = ", match)
                # For each match, increment the counter and process it
                Zcounter += 1
                orf_record.id = f"{speciesname}_{Ztype}_{Zcounter}"
                matched_ztypes.append(Ztype)  # Store the Ztype for this match
                matched_domains.append(match)  # Store the actual matched domain sequence
                test_df.loc[test_df.index[ival], 'Ztypelist'] = 'or'.join(matched_ztypes)
                test_df.loc[test_df.index[ival], 'Ztype'] = Ztype
                test_df.loc[test_df.index[ival], 'Zcount'] = Zcounter
                test_df.loc[test_df.index[ival], 'Zdomain_AA'] = match
                orf_record.seq = Seq(match)  # Update the sequence to the matched Z domain
                SeqIO.write(orf_record, outfile, "fasta")

        # If no match was found, save the sequence as a partial Z domain
        if not matched_ztypes:
            Zcounter += 1
            orf_record.id = f"{speciesname}_PartialZ_{Zcounter}" #speciesname+ "_Partial_" + str(Zcounter)
            SeqIO.write(orf_record, outfile_partial, "fasta")
            test_df.loc[test_df.index[ival], 'Ztype'] = "Partial"
            test_df.loc[test_df.index[ival], 'Zcount']= Zcounter
            test_df.loc[test_df.index[ival], 'Zdomain_AA'] = str(orf_record.seq)

    # Close the output files
    outfile.close() #file with Z domains
    outfile_partial.close() #file with partials saved for additional analysis

    # Save the updated DataFrame to a CSV file 
    test_df = test_df.loc[test_df['Ztype']!='x'].reset_index(drop=True)
    test_df.to_csv(Zdomain_folder + "/" + speciesname + "_data.csv", index=False,header=False) 
    print("Run Save Z: \n" )# ,test_df) 
    # Return the total number of Z domains
    return test_df #f"Total Z domains: {Zcounter}" #return "Total Z domains:" + str(Zcounter)

# Main function
def main():
    # Parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("input_file", help="Input csv file with species name and genome file required")
    parser.add_argument("hmm_profile", help="Input hmm profile")
    parser.add_argument("prefix", help="Output prefix required")
    parser.add_argument("Zdomain_file", help="Gene structure required")
    parser.add_argument("final_table_name", help="Final table name required") # note this is a "prefix" and not the full file name, example: bat_all_C
    args = parser.parse_args()

    # Load input CSV
    excel_file = pd.read_csv(args.input_file)
    hmmprofile = str(args.hmm_profile)
    prefix_str = str(args.prefix) #rename this variable to path
    Zdomain_csv = str(args.Zdomain_file)
    catfile = str(args.final_table_name)
    # Ensure the output directory exists
    output_folder = prefix_str  # Use the prefix argument as the output folder
    os.makedirs(output_folder, exist_ok=True)
    hmmer_folder = prefix_str + "/hmmer"
    os.makedirs(hmmer_folder, exist_ok=True)
    getorf_folder = prefix_str + "/getorf"
    os.makedirs(getorf_folder, exist_ok=True)
    Zdomain_folder = prefix_str+ "/Zdomains"
    os.makedirs(Zdomain_folder, exist_ok=True)
    PyBEDTools_folder = prefix_str+ "/PyBEDTools"
    os.makedirs(PyBEDTools_folder, exist_ok=True)
    # Initialize output folders complete
    
    df_list = []
    Zcounter = 0

    # Process each species
    for row in excel_file.index:

        speciesname = re.sub(" ", "_", excel_file['Species'][row])
        print(speciesname)
        
        # Unzip the fastafile
        fasta = excel_file['Fasta_file'][row]
        fastafile = unzip_fasta(fasta)
        print(fastafile)

        # HMMER, PyBedTools and GetORF processing
        hmmerfile = f"{hmmer_folder}/{speciesname}_tbl"
        hmmerbed = hmmerfile + ".bed" 
        pybt_tbl = f"{PyBEDTools_folder}/{speciesname}_pybt_tbl.csv"
        hmmer_nc_0 = f"{PyBEDTools_folder}/{speciesname}_hmmer_nc_0.fasta"
        hmmer_nc = f"{PyBEDTools_folder}/{speciesname}_hmmer_nc.fasta"
        orffile = f"{getorf_folder}/{speciesname}_orf.fasta"
        # Check if hmmerbed exists and is non-empty
        print(datetime.now())
        if not (os.path.exists(hmmerbed) and os.path.getsize(hmmerbed) > 0) :
            run_hmmer(fastafile, hmmprofile, hmmerfile)
            print(datetime.now())
        # Check if hmmer_nc_0 exists and is non-empty
        if not (os.path.exists(pybt_tbl) and os.path.getsize(pybt_tbl) > 0):
            run_pybedtools(hmmerbed, fastafile, hmmer_nc_0, pybt_tbl)
            print(datetime.now())
        # Check if hmmer_nc and orffile exist and are non-empty
        if not (os.path.exists(orffile) and os.path.getsize(orffile) > 0):
            run_getorf(hmmer_nc_0,hmmer_nc,orffile)
            print(datetime.now())

        # Load HMMER hits into dataframe
        bedfile = pd.read_csv(pybt_tbl, sep=',', header=None)  
        #temp_bed = bedfile.sort_values(by=[0,1,2])
        header = ['chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand','PBT_NC_sequence']
        bedfile.columns = header[:len(bedfile.columns)]
        test_df = bedfile # temp_bed.reset_index(drop=True) change row index

        
        #print(test_df.head())
        test_df['species'] = speciesname
        test_df['fasta'] = fastafile
        test_df['Ztypelist'] = 'x'
        test_df['Ztype'] = 'x'
        test_df['Zcount'] = 0
        test_df['namesearch'] = test_df['chrom'].astype(str) + "_" + test_df['chromStart'].astype(str) + "-" + test_df['chromEnd'].astype(str)
        test_df['Zdomain_AA'] = 'x'

        # Save Z domains to output
        Zdomain_aa_file_name = f"{Zdomain_folder}/{speciesname}_Zdomains.fasta"
        Zdomain_partial_aa_file_name = f"{Zdomain_folder}/{speciesname}_partial_hits.fasta"
        save_z_domains(orffile, test_df, Zdomain_folder, speciesname, Zdomain_aa_file_name, Zdomain_partial_aa_file_name, Zdomain_csv) 
        print(datetime.now())

    #outside for loop :
    # add step where all csv files are combined and save at root folder path
    os.system("cat "+ Zdomain_folder + "/" +"*_data.csv > " + prefix_str + "/" + catfile + "_final_data_table.csv")
    
    

if __name__ == "__main__":
    main()
