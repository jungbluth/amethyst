FROM ubuntu:22.04

# docker build --platform=linux/amd64 -t amethyst:latest .
# docker run --platform=linux/amd64 -it amethyst:latest /bin/bash

# source activate <conda-env> to test

RUN apt-get -y -m update && apt-get install -y build-essential wget autoconf git unzip python3 python3-pip

# Set locale environment variables
ENV LANG=en_US.UTF-8 \
    LANGUAGE=en_US:en \
    LC_ALL=en_US.UTF-8

# Install Miniconda
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /miniconda.sh \
    && bash /miniconda.sh -b -p /miniconda \
    && rm /miniconda.sh

# Add Miniconda to PATH
ENV PATH="/miniconda/bin:${PATH}"

# Install Mamba from Conda-Forge
RUN conda install -c conda-forge mamba

# Install Snakemake
RUN mamba install -c conda-forge -c bioconda snakemake

# Install Amethyst
RUN git clone https://github.com/aomlomics/amethyst

# Install Amethyst mg-qc module
RUN mamba create -n mg-qc -c bioconda multiqc fastqc

# Install Amethyst mg-norm module
RUN mamba create -n mg-norm -c bioconda bbmap 

# Install Amethyst Sequence Trimming module
RUN git clone https://github.com/KGerhardt/multitrim && \
  cd multitrim && \
  conda env create -n multitrim -f multitrim.yml

# Install Amethyst mg-assembly module
RUN mamba create -n mg-assembly -c bioconda megahit seqkit prodigal 

# Install Amethyst mg-assembly2 module
RUN mamba create -n mg-assembly2 -c conda-forge -c bioconda -c defaults prokka

# Install Amethyst mg-binning module
RUN mamba create -n mg-binning -c bioconda bowtie2

# Install Amethyst mg-binning2 module
RUN mamba create -n mg-binning2 -c bioconda bowtie2 minimap2 maxbin2 metabat2 seqkit

# Install Amethyst mg-binning3 module
RUN mamba create -n mg-binning3 -c bioconda drep gtdbtk checkm-genome 

# Install Amethyst mg-checkm module
# Skipping because already installed fine in previous module
#RUN mamba create -n mg-checkm -c bioconda checkm-genome

# Install Amethyst mg-diversity module
RUN mamba create -n mg-diversity -c bioconda -c biobakery sourmash krona humann metaphlan seqkit


WORKDIR "/tmp"
