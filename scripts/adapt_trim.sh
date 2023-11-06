#!/bin/bash
cd ./00_data/fastq
mkdir ../../01_qc/trimmed_reads
mkdir ../../01_qc/trimmed_reads/R2
mkdir ../../01_qc/trimmed_reads/R1
cd R1
for i in *R1.fastq*
do
python3 ../../../scripts/multitrim.py -u $i --max -o  ../../../01_qc/trimmed_reads/R1
done

cd ../R2
for i in *R2.fastq*
do
python3 ../../../scripts/multitrim.py -u $i --max -o ../../../01_qc/trimmed_reads/R2
done
