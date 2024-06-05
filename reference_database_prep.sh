#!/usr/bin/env bash

# other database ideas: Anvio, NCBI, hmm suites, 

# Commands returning a non-zero exit code will cause an immediate exit
set -e

# generate main Amethyst database
#mkdir /Amethyst/dbs/ && cd /Amethyst/dbs/

# download Sourmash GTDB R08-RS214 genomic representatives database (k31)
wget https://farm.cse.ucdavis.edu/~ctbrown/sourmash-db/gtdb-rs214/gtdb-rs214-reps.k31.zip

# download Sourmash GTDB R08-RS214 all genomes database (k31)
wget https://farm.cse.ucdavis.edu/~ctbrown/sourmash-db/gtdb-rs214/gtdb-rs214-k31.zip

# download Sourmash taxa sheet (note: link broken)
# wget https://osf.io/p6z3

# download CheckM1 database
wget https://zenodo.org/records/7401545/files/checkm_data_2015_01_16.tar.gz

# download CheckM2 database
wget https://zenodo.org/records/4626519/files/uniref100.KO.v1.dmnd.gz

# download GTDB R220 database
wget https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/auxillary_files/gtdbtk_package/full_package/gtdbtk_data.tar.gz

# download Kaiju databases
wget https://kaiju-idx.s3.eu-central-1.amazonaws.com/2023/kaiju_db_nr_euk_2023-05-10.tgz
wget https://kaiju-idx.s3.eu-central-1.amazonaws.com/2023/kaiju_db_nr_2023-05-10.tgz
wget https://kaiju-idx.s3.eu-central-1.amazonaws.com/2023/kaiju_db_refseq_2023-05-23.tgz
wget https://kaiju-idx.s3.eu-central-1.amazonaws.com/2023/kaiju_db_progenomes_2023-05-25.tgz
wget https://kaiju-idx.s3.eu-central-1.amazonaws.com/2023/kaiju_db_fungi_2023-05-26.tgz
wget https://kaiju-idx.s3.eu-central-1.amazonaws.com/2023/kaiju_db_viruses_2023-05-26.tgz
wget https://kaiju-idx.s3.eu-central-1.amazonaws.com/2023/kaiju_db_plasmids_2023-05-26.tgz
wget https://kaiju-idx.s3.eu-central-1.amazonaws.com/2023/kaiju_db_rvdb_2023-05-26.tgz