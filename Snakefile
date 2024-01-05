# Grab a list of all the samples
# Assumes all the file names have the format R[12]_{sample}_R[12].fastq
import glob
import re

pattern = re.compile(r'00_data/fastq/R1/R1_(.*)_R1.fastq')
files = glob.glob('00_data/fastq/R1/R1_*.fastq')

SAMPLES = []
for file in files:
    match = pattern.match(file)
    if match:
        SAMPLES.append(match.group(1))

# Master rule that snakemake uses to determine which files need to be 
# generated.
rule all:
    input:
        expand("00_data/fastq/fastqc-R1/R1_{sample}_R1_fastqc.html", sample=SAMPLES),
        expand("00_data/fastq/fastqc-R2/R2_{sample}_R2_fastqc.html", sample=SAMPLES),
        "00_data/fastq/fastqc-R1/multiqc_report.html",
        "00_data/fastq/fastqc-R2/multiqc_report.html",
        expand("01_qc/trimmed_reads/test/{sample}_1.fq.gz", sample=SAMPLES),
        expand("01_qc/trimmed_reads/test/{sample}_2.fq.gz", sample=SAMPLES),
        expand("01_qc/{sample}_normalized.fq.gz", sample=SAMPLES),
        expand("02_assembly/{sample}/{sample}.contigs.fa", sample=SAMPLES),
        expand("02_assembly/{sample}.1.bt2", sample=SAMPLES),
        expand("02_assembly/{sample}/{sample}.sam", sample=SAMPLES),
        expand("02_assembly/{sample}/prodigal/{sample}_contig_cords.gbk", sample=SAMPLES),
        expand("02_assembly/{sample}/prodigal/{sample}_contig_orfs.faa", sample=SAMPLES),
        expand("02_assembly/{sample}/prodigal/{sample}_contig_orfs.fna", sample=SAMPLES),
        expand("02_assembly/{sample}/prodigal/{sample}.gff", sample=SAMPLES),
        expand("02_assembly/{sample}/prodigal/{sample}.faa", sample=SAMPLES)

# Run all the samples through FastQC 
rule fastqc: 
    priority: 1
    conda: 
        "mg-qc"
    input:
        r1 = "00_data/fastq/R1/R1_{sample}_R1.fastq",
        r2 = "00_data/fastq/R2/R2_{sample}_R2.fastq"
    output:
        o1 = "00_data/fastq/fastqc-R1/R1_{sample}_R1_fastqc.html", 
        o2 = "00_data/fastq/fastqc-R2/R2_{sample}_R2_fastqc.html"
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
    priority: 1
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
        cd 00_data/fastq/fastqc-R1
        multiqc --export . 
        cd ..
        cd fastqc-R2
        multiqc --export . 
        """

# Run fastp
rule fastp:
    priority: 1
    conda: 
        "multitrim"
    input:
        r1 = "00_data/fastq/R1/R1_{sample}_R1.fastq",
        r2 = "00_data/fastq/R2/R2_{sample}_R2.fastq"
    output:
        o1 = "01_qc/trimmed_reads/test/{sample}_1.fq.gz",
        o2 = "01_qc/trimmed_reads/test/{sample}_2.fq.gz",
        o3 = "01_qc/trimmed_reads/test/{sample}_test_report.html"
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
        -h {output.o3}
        """
# Run bbnorm
rule bbnorm:
    priority: 1
    conda:
        "mg-norm"
    input:
        r1 = "01_qc/trimmed_reads/test/{sample}_1.fq.gz",
        r2 = "01_qc/trimmed_reads/test/{sample}_2.fq.gz" 
    output:
        o1 = "01_qc/{sample}_normalized.fq.gz"
    params:
        r1 = "01_qc/trimmed_reads/test/{sample}_1.fq.gz",
        r2 = "01_qc/trimmed_reads/test/{sample}_2.fq.gz"
    log:
        "logs/bbnorm/{sample}.log"
    benchmark:
        "benchmarks/bbnorm/{sample}.txt"
    shell:
        """
         /work/hpc/users/rstoro/git/amethyst/bbmap/bbnorm.sh in={input.r1} in2={input.r2} out={output.o1} target=100 min=5 interleaved=FALSE 
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
    priority: 1
    conda:
        "mg-assembly"
    input:
        r1 = "01_qc/trimmed_reads/test/{sample}_1.fq.gz",
        r2 = "01_qc/trimmed_reads/test/{sample}_2.fq.gz"
    output:
        o1 = "02_assembly/{sample}/{sample}.contigs.fa"
    params:
        r1 = "02_assembly/{sample}_R1.fq",
        r2 = "02_assembly/{sample}_R2.fq",
        outfolder = "02_assembly/{sample}",
        prefix = "{sample}"
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
            --min-contig-len 20 --out-prefix {params.prefix} \
            -o {params.outfolder}
        rm {params.r1} {params.r2}
        """
#build db
rule bbdb:
    priority: 2
    conda:
        "mg-binning"
    input:
        seq = "02_assembly/{sample}/{sample}.contigs.fa"
    output:
        o1 = "02_assembly/{sample}.1.bt2",
        o2 = "02_assembly/{sample}.2.bt2",
        o3 = "02_assembly/{sample}.3.bt2",
        o4 = "02_assembly/{sample}.4.bt2",
        o5 = "02_assembly/{sample}.rev.1.bt2",
        o6 = "02_assembly/{sample}.rev.2.bt2"
    params:
        basename="02_assembly/{sample}"
    threads: 32
    log:
        "logs/bbmap/{sample}.log"
    benchmark:
        "benchmarks/bbmap/{sample}.txt"
    shell:
        """
        bowtie2-build --threads 32 {input.seq} {params.basename}
        """
#map and make sam file
rule bbmap:
    conda:
        "mg-binning"
    input:
        r1 = "01_qc/trimmed_reads/test/{sample}_1.fq.gz",
        r2 = "01_qc/trimmed_reads/test/{sample}_2.fq.gz",
        o1 = "02_assembly/{sample}.1.bt2",
        o2 = "02_assembly/{sample}.2.bt2",
        o3 = "02_assembly/{sample}.3.bt2",
        o4 = "02_assembly/{sample}.4.bt2",
        o5 = "02_assembly/{sample}.rev.1.bt2",
        o6 = "02_assembly/{sample}.rev.2.bt2"
    output:
        o1 = "02_assembly/{sample}/{sample}.sam",
        log = "02_assembly/{sample}.bowtie2.log"
    params:
        o2 = "02_assembly/{sample}"
    threads: 32
    log:
        "logs/bbmap/{sample}.log"
    benchmark:
        "benchmarks/bbmap/{sample}.txt"
    shell:
        """
        bowtie2 --threads 32 -x {params.o2} -1 {input.r1} \
        -2 {input.r2} -S {output.o1} > {output.log}
        """
rule prodigal:
    conda:
        "mg-assembly"
    input:
        r1 = "02_assembly/{sample}/{sample}.contigs.fa"
    output:
        o1 = "02_assembly/{sample}/prodigal/{sample}_contig_cords.gbk",
        o2 = "02_assembly/{sample}/prodigal/{sample}_contig_orfs.faa",
        o3 = "02_assembly/{sample}/prodigal/{sample}_contig_orfs.fna"
    threads: 32
    log:
        "logs/prodigal/{sample}.log"
    benchmark:
        "benchmarks/prodigal/{sample}.txt"
    shell:
        """
        prodigal -i {input.r1} -o {output.o1} -a {output.o2} -d {output.o3}
        """
rule prokka:
    conda:
        "mg-assembly2"
    input:
        r1 = "02_assembly/{sample}/{sample}.contigs.fa"
    output:
        o1 = "02_assembly/{sample}/prodigal/{sample}.gff",
        o2 = "02_assembly/{sample}/prodigal/{sample}.faa"
    params:
        outfolder = "02_assembly/{sample}/prodigal",
        prefix = "{sample}" 
    threads: 32
    log:
        "logs/prokka/{sample}.log"
    benchmark:
        "benchmarks/prokka/{sample}.txt"
    shell:
        """
        prokka {input.r1} --outdir {params.outfolder} --prefix {params.prefix}
        """

#####FROM THIS POINT ON IS PSEUDO CODE AND WILL NOT RUN CORRECTLY
#1/5/2024

rule sourmash:
    conda:
        "mg-diversity"
    input:
        r1 = "01_qc/trimmed_reads/test/{sample}_1.fq.gz",
        r2 = "01_qc/trimmed_reads/test/{sample}_2.fq.gz",
        r3 = "01_qc/interleaved/{sample}_interleaved.fq.gz" 
    output:
        o1 = "01_qc/interleaved/{sample}_interleaved.fq.gz",
        o2 = "{params.sample}_reads.sig",
        o3 = "{params.sample}_sourmash_gather_out.csv",
        o4 = "{params.sample}_sourmash_tax"
    params:
        sample = SAMPLES,
        outfolder = "02_assembly/sourmash",
        outfolder2 = "02_assembly/sourmash/tax_out",
        db = "/phodnet/OCED/databases/sourmash/GTDB_k31.sbt"
    threads: 32
    log:
        "logs/sourmash/{params.sample}.log"
    benchmark:
        "benchmarks/sourmash/{params.sample}.txt"
    shell:
        """
        /work/hpc/users/rstoro/git/amethyst/bbmap/bbnorm.sh in={input.r1} in2={input.r2} out={output.o1} 
        sourmash sketch dna {input.r3} -o {params.outfolder}/{output.o2}
        sourmash gather {output.o2} {db} -o {params.outfolder}/{output.o3}
        sourmash tax metagenome -g {output.o3} -t gtdb-rs202.taxonomy.v2.with-strain.csv -o {output.o4} --output-dir {params.outfolder2} --output-format krona --rank species
        sourmash tax metagenome -g {output.o3} -t gtdb-rs202.taxonomy.v2.with-strain.csv -o {output.o4} --output-dir {params.outfolder2} --output-format csv_summary
        """
rule maxbin2:
    conda:
        "mg-binning"
    input:
        r1 = "02_assembly/{params.sample}/{params.sample}.contigs.fa",
        r2 = "01_qc/trimmed_reads/test/{sample}_1.fq.gz",
        r3 = "01_qc/trimmed_reads/test/{sample}_2.fq.gz"
    output:
        o1 = "02_assembly/{params.sample}_MaxBin"
    params:
        sample = SAMPLES
    threads: 32
    log:
        "logs/maxbin/{params.sample}.log"
    benchmark:
        "benchmarks/maxbin/{params.sample}.txt"
    shell:
        """
        run_MaxBin.pl -contig {input.r1} -min_contig_length 1000 \
        -reads {input.r2} -reads2 {input.r3} \
        -out {output.o1} -thread 32
        """
rule checkm:
    conda:
        "mg-checkm"
    input:
        r1 = "02_assembly/{params.sample}_MaxBin",
        r2 = "02_assembly/{params.sample}_temp_folder"
    output:
        o1 = "02_assembly/{params.sample}_checkm"
    params:
        sample = SAMPLES
    threads: 32
    log:
        "logs/checkm/{params.sample}.log"
    benchmark:
        "benchmarks/checkm/{params.sample}.txt"
    shell:
        """
        checkm lineage_wf -t 32 -x fa --tmp {input.r2} {input.r1} {output.o1} >output.log
        """
rule dRep:
    conda:
        "mg-binning"
    input:
        r1 = "02_assembly/{params.sample}_checkm",
        r2 = "dRep_MAGs_map.csv"
    output:
        o1 = "02_assembly/dRep_out"
    params:
        sample = SAMPLES
    threads: 32
    log:
        "logs/dRep/{params.sample}.log"
    benchmark:
        "benchmarks/dRep/{params.sample}.txt"
    shell:
        """
        dRep dereplicate {output.o1} -g {input.r1} --genomeInfo {input.r2} -sa 0.95 -nc 0.2 --ignoreGenomeQuality
        """
rule GTDBtk:
    conda:
        "mg-binning"
    input:
        r1 = "02_assembly/dRep_out"
    output:
        o1 = "03_assignment/GTDBtk"
    params:
        sample = SAMPLES
    threads: 32
    log:
        "logs/GTDBtk/{params.sample}.log"
    benchmark:
        "benchmarks/GTDBtk/{params.sample}.txt"
    shell:
        """
        gtdbtk classify_wf --genome_dir {input.r1} --out_dir {output.o1}
        """
