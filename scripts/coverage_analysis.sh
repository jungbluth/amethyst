#!/bin/bash
mkdir ./01_qc/trimmed_reads/R1/nonpareil
cd ./01_qc/trimmed_reads/R1/nonpareil/
nonpareil -s ../../../../00_data/fastq/R1/sub_read1.fq -b out_prefix -d 0.7 -t 20 -R 9500 -X 100 -t 4 -T alignment
