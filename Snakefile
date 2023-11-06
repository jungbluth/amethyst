# Requires X conda environments as outlined here: https://github.com/rckarns8/amethyst

# GLOBALS

#configfile: "config.yaml"

# RULE ORDER DIRECTIVES

# PSEUDO-RULES: LARGER RULES COMPOSED OF MANY SMALLER RULES

rule read_qc: 
    input:
        "00_data/fastq/"
    output:
    conda: 
        "mg-qc"
    shell: 
        "scripts/read_qc.sh"

rule adapt_trim:
    input:
	"/00_data/fastq"
    output:
    conda:
        "mg-trim"
    shell:
        "scripts/adapt_trim.sh"


        

