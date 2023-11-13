# Requires conda environments as outlined here: https://github.com/rckarns8/amethyst/wiki/0.-Quick-Start-Guide
# GLOBALS
# RULEORDER DIRECTIVES
# PSEUDO-RULES: 
rule read_qc: 
    input:
        "./00_data/fastq/R1",
        "./00_data/fastq/R2"
    output:
        "./00_data/fastq/fastqc-R1/",
        "./00_data/fastq/fastqc-R2/"
    conda: 
        "mg-qc"
    shell: 
        "./scripts/read_qc.sh"
rule adapt_trim:
    input:
        "./00_data/fastq/R1",
        "./00_data/fastq/R2"
    output:
        "./01_qc/trimmed_reads/R1",
        "./01_qc/trimmed_reads/R2"
    conda:
        "multitrim"
    shell:
        "./scripts/adapt_trim.sh"
rule assembly:
    input:
        "./01_qc/trimmed_reads/R1",
        "./01_qc/trimmed_reads/R2"
    output:
        "./02_assembly/"
    conda:
        "mg-assembly"
    shell: 
        "./scripts/assembly.sh"
