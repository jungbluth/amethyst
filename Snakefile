# Requires X conda environments as outlined here: https://github.com/rckarns8/amethyst

# GLOBALS

#configfile: "config.yaml"

# RULE ORDER DIRECTIVES

# PSEUDO-RULES: LARGER RULES COMPOSED OF MANY SMALLER RULES

rule read_qc: 
    input:
        "00_data/fastq/"
    output:
	"00_data/fastq/fastqc-R1/"
	"00_data/fastq/fastqc-R2/"
    conda: 
        "mg-qc"
    shell: 
        "scripts/read_qc.sh"

rule adapt_trim:
    input:
	"00_data/fastq"
    output:
	"01_qc/trimmed_reads"
    conda:
        "multitrim"
    shell:
        "scripts/adapt_trim.sh"

rule assembly:
    input:
	"01_qc/trimmed_reads/"
    output:
    conda:
	"mg-assembly"
    shell: 
        "scripts/assembly.sh"

        

