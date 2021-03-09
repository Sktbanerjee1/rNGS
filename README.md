# rNGS
Bash scripts for reproducible NGS data processing.

[DOI](https://zenodo.org/badge/345268065.svg)

## Getting Started

```bash
#!/bin/bash
# make sure git is installed
git --version
# clone this repository
git clone https://github.com/Sktbanerjee1/rNGS.git \
&& cd rNGS/ \
&& ls -l
```

## Required dependencies

| Script        | Reuires       |
| ------------- |:-------------:|
| scripts/fastqc.sh      | [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) |
| scripts/trim_fastq.sh      | [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)|
|scripts/genome_generate.sh | [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)|
|scripts/align.sh | [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)|
|scripts/genomeSize.sh | [Samtools](http://www.htslib.org/)|
|scripts/sam2bed.sh | [Samtools](http://www.htslib.org/), [Picard](https://broadinstitute.github.io/picard/), [Sambamba](https://lomereiter.github.io/sambamba/), [Bedtools](https://bedtools.readthedocs.io/en/latest/) |

By default, the scripts expect these dependencies to be available in  the environment `${PATH}`. Existining installation can be usually checked `<tool_name> --help` or `<tool_name> --version` commands in bash.

For an easier installation, here we use [conda](https://docs.conda.io/en/latest/), an environment manager that provide these softwares as packages.

Create a new file, `install.sh` in the current working directory and copy paste the following code.

```bash
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

# deactivate the environment
conda deactivate rNGS

# finish
echo "installation complete!"
echo "run conda activate rNGS"
```

Run the installer using the following command:

```bash
bash install.sh
```

## Running analysis

### 1. Sequence QC

```bash
bash scripts/fastqc.sh

# Usage: scripts/fastqc.sh -i target_dir -o out_dir -m run_mode
# 	-i path to directory with fastq files
# 	-o Path to the directory where the outputs will be written
# 	-m value should be one of Trimmed/Untrimmed
```

### 2. Read Trimming

```bash
bash scripts/trim_fastq.sh

# Usage: scripts/trim_fastq.sh -i target_dir -o out_dir -s sample_name -a adapter_file -w window -m min_len
# 	-i path to directory with fastq files
# 	-o Path to the directory where the outputs will be written
# 	-s Name of the sample
# 	-a Path to the adapter fasta file
# 	-w Sliding Window
# 	-m Minimum read length
```

### 3. Generate ref genome coordinates

```bash
bash scripts/genome_generate.sh

# Usage: scripts/genome_generate.sh -i ref -o out_dir -p threads
# 	-i path to reference genome
# 	-o Path to the directory where the outputs will be written
# 	-p threads
```

### 4. Align reads

```bash
bash scripts/align.sh

# Usage: scripts/align.sh -i idx -f ref -r rev -o out_dir -s sample
# 	-i generated genome idx
# 	-f paired fwd fastq/fastq.gz
# 	-r paired rev fastq/fastq.gz
# 	-o Path to the directory where the outputs will be written
# 	-s Name of the sample

```

### 5. Convert SAM to BED

```bash
# use scripts/genomeSize.sh
# Generates chromosome sizes
# from the alignment map

bash scripts/sam2bed.sh 

# Usage: scripts/sam2bed.sh -i target -o out_dir -c chr_sizes
# 	-i input sam file
# 	-s sample_name
# 	-o Path to the directory where the outputs will be written
# 	-c sizes of the chromosomes
```
