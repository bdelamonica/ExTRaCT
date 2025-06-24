# ExTRaCT 

The following includes environment setup, installation requirements and example run instructions.

## Environment Setup

Anaconda is recommended but not required. The current code uses Python version 3.7 or greater and conda version 4.12.0.

```#!/bin/sh

conda create -n a3pipeline
conda activate a3pipeline

git clone https://github.com/bdelamonica/ExTRaCT.git

# alternatively
wget https://github.com/bdelamonica/ExTRaCT/archive/refs/heads/main.zip
unzip main.zip

cd ExTRaCT

```

## Instructions to download dependencies

1. MAFFT

Link: https://mafft.cbrc.jp/alignment/software/source.html

Follow instructions in link to untar, compile and install.

`wget mafft-7.490-without-extensions-src.tgz`

Make sure to save the file path to add to example code: /path_to/mafft

2. HMMER

Link : http://hmmer.org/documentation.html

Follow the documentation to install and run HMMER version hmmer-3.1b2. 

Option for conda install:

`conda install -c biocore hmmer`

3. SCIPIO

Link: https://www.webscipio.org/webscipio/download_scipio

`wget https://www.webscipio.org/webscipio/download_scipio_1_4`

SCIPIO has various dependencies. Below we list crucial documentation and tips to download all

* Convert2bed

Link: https://bedops.readthedocs.io/en/latest/content/installation.html#linux

* BioPerl 

`conda install bioconda::perl-bioperl`

or

`conda install -c bioconda perl-bioperl`

Check version

`perl -MBio::SeqIO -e 'printf "%vd ", $Bio::SeqIO::VERSION'`

Alternative load

`conda install bioconda::perl-bio-featureio`

* YAML

`conda install bioconda::perl-yaml`

Check version

`perl -MYAML -e 'print $YAML::VERSION ."\n";'`

* BLAT

Link: https://genome.ucsc.edu/FAQ/FAQblat.html#blat9 

`rsync -aP hgdownload.soe.ucsc.edu::genome/admin/exe/linux.x86_64/ ./blat`

Save bat = /path_to_blat/blat/blat/

`export PATH=/path_to_blat/blat/blat/:$PATH`

* Once all dependencies are saved test

Check SCIPIO runs

`perl path/to/scipio.1.4.pl`


5. Python requirements

StringIO

AlignIO

Pybedtools

`conda install bioconda::pybedtools`

Openpyxl

`conda install anaconda::openpyxl`

Orffinder

`pip3 install orffinder`

6. EMBOSS

Link : http://emboss.open-bio.org/html/adm/ch01s01.html,http://emboss.open-bio.org/html/adm/ch01s03.html or https://emboss.sourceforge.net/download/#Stable

`wget -m 'ftp://emboss.open-bio.org/pub/EMBOSS/'`

Follow gunzip, tar, compile and make instructions

Test:

`/path_to/emboss/getorf`

`export PATH=/path/to/emboss.open-bio.org/pub/EMBOSS/EMBOSS-6.6.0/emboss/:$PATH`

7. RAxML



## Run Example

Follow the instructions in example_Z_domain_bat_all_C.sh and be sure to make the appropriate file path changes before running.

```#!/bin/sh
Conda activate a3pipeline

export PATH=[location of emboss]/emboss.open-bio.org/pub/EMBOSS/EMBOSS-6.6.0/emboss/:$PATH

export PATH=[location of raxml]/standard-RAxML-master:$PATH

export PATH=[location of sciptio]/scipio-1.4/:$PATH

export PATH=[location of blat]/blat/blat/:$PATH

bash example_Z_domain_bat_all_C.sh

```

The example code will generate results for a single test case on a single target genome. To test all 27 use cases follow Z_domain_bat_all_C.sh directions in the comments.

To run your own experiment make changes to crucial input files and be sure to maintain the correct structure:
1. Input reference species gene FASTA file (example: hmmer_profiles_files/6bat_nucl_all.fasta)
2. Input target genome CSV table (example: species_names.csv) and Save target genomes in folder (example: Bat1kgenomes/..)
3. Input motif files (example: Zdomain_C.csv)

Results are saved under example_run and subfolders.
