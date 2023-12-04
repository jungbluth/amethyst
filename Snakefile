# Grab a list of all the samples
# Assumes all the file names have the format R[12]_{sample}_R[12].fastq
SAMPLES = glob_wildcards("00_data/fastq/R1/R1_{sample}_R1.fastq").sample

# Master rule that snakemake uses to determine which files need to be 
# generated. 
rule all:
    input:
        expand("00_data/fastq/fastqc-R1/R1_{sample}_R1_fastqc.html", 
            sample=SAMPLES),
        expand("00_data/fastq/fastqc-R2/R2_{sample}_R2_fastqc.html", 
            sample=SAMPLES),
        "00_data/fastq/fastqc-R1/multiqc_report.html",
        "00_data/fastq/fastqc-R2/multiqc_report.html",
        expand("01_qc/trimmed_reads/R1/unpaired.post_trim_R1_{sample}_R1.fq.gz",
            sample=SAMPLES),
        expand("01_qc/trimmed_reads/R2/unpaired.post_trim_R2_{sample}_R2.fq.gz",
            sample=SAMPLES),
        expand("02_assembly/{sample}/{sample}.contigs.fa",
            sample=SAMPLES)

# Run all the samples through FastQC 
rule fastqc: 
    conda: 
        "amethyst-testing-qc"
    input:
        r1 = "00_data/fastq/R1/R1_{sample}_R1.fastq",
        r2 = "00_data/fastq/R2/R2_{sample}_R2.fastq"
    output:
        "00_data/fastq/fastqc-R1/R1_{sample}_R1_fastqc.html",
        "00_data/fastq/fastqc-R2/R2_{sample}_R2_fastqc.html"
    params:
        outfolder1 = "00_data/fastq/fastqc-R1",
        outfolder2 = "00_data/fastq/fastqc-R2"
    threads: 5 
    log:
        "logs/fastqc/{sample}.log"
    benchmark:
        "benchmarks/fastqc/{sample}.txt"
    shell:
        """
        fastqc -t {threads} -o {params.outfolder1} {input.r1}
        fastqc -t {threads} -o {params.outfolder2} {input.r2}
        """

# Run MultiQC on the FastQC reports
rule multiqc:
    conda:
        "amethyst-testing-qc"
    output:
        "00_data/fastq/fastqc-R1/multiqc_report.html",
        "00_data/fastq/fastqc-R2/multiqc_report.html"
    params:
        outfolder1 = "00_data/fastq/fastqc-R1",
        outfolder2 = "00_data/fastq/fastqc-R2"
    log:
        "logs/multiqc/multiqc.log"
    benchmark:
        "benchmarks/multiqc/multiqc.txt"
    shell:
        """
        multiqc --export -o {params.outfolder1} {params.outfolder1}
        multiqc --export -o {params.outfolder2} {params.outfolder2}
        """

# Run multitrim
# multitrim generates files with the same name so we have to create a 
# temp folder for each run. One thing to look at in the future is using
# its --prefix option.
rule multitrim:
    conda:
        "amethyst-testing-mt"
    input:
        r1 = "00_data/fastq/R1/R1_{sample}_R1.fastq",
        r2 = "00_data/fastq/R2/R2_{sample}_R2.fastq"
    output:
        "01_qc/trimmed_reads/R1/unpaired.post_trim_R1_{sample}_R1.fq.gz",
        "01_qc/trimmed_reads/R2/unpaired.post_trim_R2_{sample}_R2.fq.gz"
    params:
        outfolder1 = "01_qc/trimmed_reads/R1",
        outfolder2 = "01_qc/trimmed_reads/R2",
        tempfolder1 = "01_qc/trimmed_reads/R1/{sample}",
        tempfolder2 = "01_qc/trimmed_reads/R2/{sample}",
    threads: 5
    log:
        "logs/multitrim/{sample}.log"
    benchmark:
        "benchmarks/multitrim/{sample}.txt"
    shell:
        """
        mkdir -p {params.tempfolder1}
        multitrim -u {input.r1} -o {params.tempfolder1} -t {threads}  
        mv {params.tempfolder1}/* {params.outfolder1}
        rmdir {params.tempfolder1}
        mkdir -p {params.tempfolder2}
        multitrim -u {input.r2} -o {params.tempfolder2} -t {threads}  
        mv {params.tempfolder2}/* {params.outfolder2}
        rmdir {params.tempfolder2}
        """

# Run megahit
# snakemake will create the output folders since that is the location of the 
# output files we specify. megahit refuses to run if its output folder already
# exists, so because of this, we have to remove the folder snakemake creates
# before we do anything.
# Right now megahit is set to use all the cores and 0.85% of the machine's
# memory. This will probably need to be adjusted when used under other
# situations.
rule megahit:
    conda:
        "amethyst-testing-as"
    input:
        r1 = "01_qc/trimmed_reads/R1/unpaired.post_trim_R1_{sample}_R1.fq.gz",
        r2 = "01_qc/trimmed_reads/R2/unpaired.post_trim_R2_{sample}_R2.fq.gz"
    output:
        "02_assembly/{sample}/{sample}.contigs.fa"
    params:
        r1 = "02_assembly/{sample}_R1.fq",
        r2 = "02_assembly/{sample}_R2.fq",
        sample = "{sample}",
        outfolder = "02_assembly/{sample}"
    threads: 20
    log:
        "logs/megahit/{sample}.log"
    benchmark:
        "benchmarks/megahit/{sample}.txt"
    shell:
        """
        rm -rf {params.outfolder}
        zcat {input.r1} > {params.r1}
        zcat {input.r2} > {params.r2}
        megahit -1 {params.r1} -2 {params.r2} -m 0.85 -t {threads} \
            --min-contig-len 20 --out-prefix {params.sample} \
            -o {params.outfolder}
        rm {params.r1} {params.r2}
        """
