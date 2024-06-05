#!/usr/bin/env bash

# Commands returning a non-zero exit code will cause an immediate exit
set -e

# generate main Amethyst database
mkdir /Amethyst/dbs/ && cd /Amethyst/dbs/

# download Sourmash GTDB R06-RS202 genomic representatives database
wget https://osf.io/nqmau/download

# download Sourmash GTDB R06-RS202 all genomes database
wget https://osf.io/94mzh/download

# download Sourmash taxa sheet (note: link broken)
wget https://osf.io/p6z3w

# download CheckM1 data
wget https://zenodo.org/records/7401545/files/checkm_data_2015_01_16.tar.gz

# download GTDB R214 data
wget https://data.gtdb.ecogenomic.org/releases/release214/214.0/auxillary_files/gtdbtk_r214_data.tar.gz
