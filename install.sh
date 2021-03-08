#!/bin/bash
# make sure conda is installed
conda --version

# crate a new env
# conda create -n <env_name>
conda create -n rNGS

# activate existing conda env
# conda activate <env_name>
conda activate rNGS

# install dependencies
conda install -c bioconda \
fastqc trimmomatic \
bowtie2 samtools bedtools \
sambamba picard
