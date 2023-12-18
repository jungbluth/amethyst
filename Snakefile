# Grab a list of all the samples
# Assumes all the file names have the format R[12]_{sample}_R[12].fastq
SAMPLES, = glob_wildcards("00_data/fastq/R1/R1_{sample}_R1.fastq")

# Master rule that snakemake uses to determine which files need to be 
# generated. 
rule all:
    input:
        expand("00_data/fastq/R1/R1_{sample}_R1.fastq", 
            sample=SAMPLES),
        expand("00_data/fastq/R2/R2_{sample}_R2.fastq", 
            sample=SAMPLES),
        expand("00_data/fastq/fastqc-R1/R1_{sample}_R1_fastqc.html", 
            sample=SAMPLES),
        expand("00_data/fastq/fastqc-R2/R2_{sample}_R2_fastqc.html", 
            sample=SAMPLES),
        "00_data/fastq/fastqc-R1/multiqc_report.html",
        "00_data/fastq/fastqc-R2/multiqc_report.html",
        expand("01_qc/trimmed_reads/test/{sample}_1.fq.gz",
            sample=SAMPLES),
        expand("01_qc/trimmed_reads/test/{sample}_2.fq.gz", 
            sample=SAMPLES),
        expand("01_qc/normalized/{sample}_normalized.fq.gz", 
            sample=SAMPLES)
     #  expand("02_assembly/{sample}/{sample}.contigs.fa", 
       #     sample=SAMPLES)

# Run all the samples through FastQC 
rule fastqc: 
    conda: 
        "mg-qc"
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
        "mg-qc"
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

# Run fastp
rule fastp:
    conda: 
        "mg-trim"
    input:
        r1 = "00_data/fastq/R1/R1_{sample}_R1.fastq",
        r2 = "00_data/fastq/R2/R2_{sample}_R2.fastq"
    output:
        o1 = "01_qc/trimmed_reads/test/{sample}_1.fq.gz",
        o2 = "01_qc/trimmed_reads/test/{sample}_2.fq.gz"
    params:
        sample=SAMPLES
    log:
        "logs/fastp/{sample}.log"
    benchmark:
        "benchmarks/fastp/{sample}.txt"
    shell:
        """
        fastp \
        -i {input.r1} \
        -o {output.o1} \
        -I {input.r2} \
        -O {output.o2} \
        --detect_adapter_for_pe \
        -g -l 50 -W 4 -M 20 -w 16 \
        --cut_front \
        -R {params.sample}_test_report
        """


# Run bbnorm

rule bbnorm:
    conda:
        "mg-norm"
    input:
        r1 = "01_qc/trimmed_reads/test/{sample}_1.fq.gz",
        r2 = "01_qc/trimmed_reads/test/{sample}_2.fq.gz" 
    output:
        o1 = "01_qc/normalized/{sample}_normalized.fq.gz"
    params:
        r1 = "01_qc/trimmed_reads/test/{sample}_R1.fq.gz",
        r2 = "01_qc/trimmed_reads/test/{sample}_R2.fq.gz",
        sample=SAMPLES,
        outfolder="01_qc/normalized",
        script="bbmap/bbnorm.sh"
    threads: 20
    log:
        "logs/bbnorm/{sample}.log"
    benchmark:
        "benchmarks/bbnorm/{sample}.txt"
    shell:
        """
        for i in {input.r1};
        do
          bbmap/bbnorm.sh in={input.r1} in2={input.r2} out={params.outfolder}/${{i}}_normalized.fq.gz target=100 min=5 interleaved=FALSE -Xmx40g;
        params.sample=${{i}}
        done
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
        "mg-assembly"
    input:
        r1 = "01_qc/normalized/{params.sample}_normalized.fq",
    output:
        "02_assembly/{params.sample}/{params.sample}.contigs.fa"
    params:
        r1 = "02_assembly/{params.sample}_R1.fq",
        r2 = "02_assembly/{params.sample}_R2.fq",
        sample=SAMPLES,
        outfolder = "02_assembly/{params.sample}"
    threads: 20
    log:
        "logs/megahit/{params.sample}.log"
    benchmark:
        "benchmarks/megahit/{params.sample}.txt"
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