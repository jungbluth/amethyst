#!/bin/bash
cd ./00_data/fastq
mkdir fastqc-R1
mkdir fastqc-R2
cd R1
for i in *R1.fastq*
do
    fastqc $i -o ../fastqc-R1/
done
cd ..
cd R2
for i in *R2.fastq*
do
    fastqc $i -o ../fastqc-R2/
done
cd ..
cd fastqc-R1
multiqc --export .
cd ..
cd fastqc-R2
multiqc --export .
