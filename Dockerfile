# Dockerfile for UDN RNA-seq pipeline (adapted from GTEx)
FROM ubuntu:14.04
MAINTAINER Alden Huang

RUN apt-get update && apt-get install -y software-properties-common && add-apt-repository -y ppa:openjdk-r/ppa && \
    apt-get update && apt-get install -y \
        build-essential \
        cmake \
        curl \
        libboost-all-dev \
        libbz2-dev \
        libcurl3-dev \
        liblzma-dev \
        libncurses5-dev \
        libssl-dev \
        openjdk-8-jre-headless \
        unzip \
        wget \
        zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*

#-----------------------------
# Pipeline components
#-----------------------------

# htslib
RUN cd /opt && \
    wget --no-check-certificate https://github.com/samtools/htslib/releases/download/1.6/htslib-1.6.tar.bz2 && \
    tar -xf htslib-1.6.tar.bz2 && rm htslib-1.6.tar.bz2 && cd htslib-1.6 && make && make install && make clean

# samtools
RUN cd /opt && \
    wget --no-check-certificate https://github.com/samtools/samtools/releases/download/1.6/samtools-1.6.tar.bz2 && \
    tar -xf samtools-1.6.tar.bz2 && rm samtools-1.6.tar.bz2 && cd samtools-1.6 && \
    ./configure --with-htslib=/opt/htslib-1.6 && make && make install && make clean

# Picard tools
RUN mkdir /opt/picard-tools && \
    wget --no-check-certificate -P /opt/picard-tools/ https://github.com/broadinstitute/picard/releases/download/2.9.0/picard.jar

# STAR v2.6.0a
RUN cd /opt && \
    wget --no-check-certificate https://github.com/alexdobin/STAR/archive/2.6.0a.tar.gz && \
    tar -xf 2.6.0a.tar.gz && rm 2.6.0a.tar.gz && \
ENV PATH /opt/STAR-2.6.0a/bin/Linux_x86_64:$PATH

# RSEM v1.3.0
RUN cd /opt && \
    wget --no-check-certificate https://github.com/deweylab/RSEM/archive/v1.3.0.tar.gz && \
    tar -xvf v1.3.0.tar.gz && rm v1.3.0.tar.gz && cd RSEM-1.3.0 && make
ENV PATH /opt/RSEM-1.3.0:$PATH

# RNA-SeQC
RUN cd /opt && \
    wget --no-check-certificate https://github.com/francois-a/rnaseqc/releases/download/v1.1.9/RNA-SeQC_1.1.9.zip && \
    unzip RNA-SeQC_1.1.9.zip -d RNA-SeQC_1.1.9 && rm RNA-SeQC_1.1.9.zip

# clean up
RUN apt-get clean && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* && \
    apt-get autoclean && \
    apt-get autoremove -y && \
    rm -rf /var/lib/{apt,dpkg,cache,log}/

####
# Retrieve and build references
# Note: I am unable to test the code below. The cluster I have access to is unable to run docker; my desktop does not have enough memory. It is only provided as a reference.

# hs37d5
RUN mkdir /reference && mkdir /reference/hs37d5 && cd /reference/hs37d5 && \   
    wget --no-check-certificate http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz && \
    zcat hs37d5.fa.gz | head -n 51696830 > hs37d5.nodecoy.fa && \
    rm hs37d5.fa.gz

# gencode19
RUN mkdir /reference/gencode19 && cd /reference/gencode19 && \
    wget --no-check-certificate ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz && \
    zcat gencode.v19.annotation.gtf.gz | sed 's/^chrM/chrMT/1;s/^chr//1' > gencode.v19.annotation.hs37d5.gtf && \
    rm gencode.v19.annotation.gtf.gz

# star_oh74
RUN mkdir /reference/star_oh74 &&
    STAR --runMode genomeGenerate \
    --genomeDir /reference/star_oh74/ \
    --genomeFastaFiles /reference/hs37d5/hs37d5_nodecoy.fa \
    --sjdbGTFfile /reference/gencode19/gencode.v19.annotation.hs37d5.gtf \
    --sjdbOverhang 74 \
    --runThreadN 8


