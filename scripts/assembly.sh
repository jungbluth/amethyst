#!/bin/bash
cd /amethyst/01_qc/trimmed_reads/R1/
gunzip *.gz
cd ..
cd /R2/
gunzip *.gz
cd ..
cd ..
cd ..


mkdir -p 02_assembly
R1_files=$(find 01_qc/trimmed_reads/R1/ -name "*R1.fq")
R2_files=$(find 01_qc/trimmed_reads/R2/ -name "*R2.fq")

if [ -z "$R1_files" ]; then
    echo "No R1 input files found. Exiting script."
    exit 1
fi

if [ -z "$R2_files" ]; then
    echo "No R2 input files found. Exiting script."
    exit 1
fi

for f in $R1_files; do
    base=$(echo "$f" | sed 's/.*\///' | sed 's/_R1.fq/_R2.fq/')
    f2=$(echo "$f" | sed 's/01_qc\/trimmed_reads\/R1/01_qc\/trimmed_reads\/R2/' | sed 's/R1/R2/' $
    output_dir="02_assembly/${base}"

    echo "Processing $f and $f2..."
    echo "Output directory: $output_dir"

    if [ -f "$f" ] && [ -f "$f2" ]; then
       # if [ -d "$output_dir" ]; then
         #   echo "Output directory already exists. Removing existing files."
        #    rm -r "$output_dir"/*
       # else
          #  echo "Creating output directory."
         #   mkdir -p "$output_dir"
        #fi
        echo "Input files found. Running megahit..."
        megahit -1 "$f" -2 "$f2" -m 0.85 -t 12 --min-contig-len 20 --output-prefix "assembly" -o $
        exit_code=$?
        if [ $exit_code -eq 0 ]; then
            echo "megahit executed successfully."
        else
            echo "megahit encountered an error. Exit code: $exit_code"
        fi
