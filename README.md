# GATA6_bulk_ATAC-seq-analyisis

This repository contains the workflow and scripts used for processing bulk ATAC-seq data for the GATA6 perturbation experiment. The pipeline consists of read trimming and QC, reference indexing, read alignment with QuasR/Rbowtie, peak calling, and differential accessibility analysis.

All steps are run on the DMLS cluster.  
Environments are managed with **micromamba** and R project reproducibility is controlled with **renv**.

---

## Environment Setup

Load micromamba in your shell (add to `~/.bashrc` if not already present):

```bash
eval "$(micromamba shell hook --shell bash)"
```


1) Trimming + QC Environment (atac-trim)


```bash
micromamba create -n atac-trim -c bioconda -c conda-forge \
  trimmomatic fastqc multiqc pigz openjdk=17 -y
```


Activate:

```bash
micromamba activate atac-trim
```

Check:

```bash
trimmomatic -version
fastqc --version
multiqc --version
```

The Nextera adapter file used in trimming is automatically provided:
${CONDA_PREFIX}/share/trimmomatic/adapters/NexteraPE-PE.fa


2) Reference Preparation Environment (refprep)

Used in: bowtie_index.sbatch.sh
Tools required: bowtie-build + samtools

```bash
micromamba create -n refprep -c bioconda -c conda-forge \
  bowtie samtools -y
```

Activate:
```bash
micromamba activate refprep
```


Check:

```bash
which bowtie-build
samtools --version
```

3) R + QuasR Alignment Environment (renv)

This env provides R and system libraries; R packages are tracked inside the project via renv.

```bash
micromamba create -p /home/mchere/data/conda/envs/renv \
  -c conda-forge -c bioconda \
  r-base=4.4 r-devtools r-xml r-data.table r-rcpp r-rlang r-cli \
  bowtie samtools -y

```

Activate:
```bash
micromamba activate /home/mchere/data/conda/envs/renv
```

Inside R (first setup only):
```bash
install.packages("BiocManager")
BiocManager::install(c("QuasR","Rbowtie","GenomicRanges","rtracklayer","BiocParallel"))
install.packages("renv")
renv::init()      # or renv::restore() if lockfile already exists
renv::snapshot()
```

## Directory Structure
```bash
GATA6_bulk-ATAC-seq/
├─ scripts/                     # sbatch scripts and helper scripts
│  ├─ run_trimmomatic.sbatch.sh
│  ├─ quasr_align_mm10.sbatch.sh
│  └─ bowtie_index.sbatch.sh
├─ trimmed_fastqs/              # output of trimming
├─ bams/                        # aligned reads
├─ peaks/                       # MACS2 outputs
├─ qc/                          # FastQC / MultiQC reports
├─ R/                           # R analysis scripts
│  ├─ 01_quasr_align.R
│  ├─ 02_count_and_merge.R
│  ├─ 03_deseq2_analysis.R
│  └─ renv.lock
└─ README.md
```


## Workflow Summary
1) Trim reads
```bash
sbatch scripts/run_trimmomatic.sbatch.sh
```
sbatch scripts/run_trimmomatic.sbatch.sh

2) Prepare bowtie index (if needed)
```bash
sbatch scripts/bowtie_index.sbatch.sh /path/to/genome.fa
```

3) Align reads with QuasR/Rbowtie

Edit quasr_samples_full.tsv to list sample names and paired FASTQs.
Run:

```bash
sbatch scripts/quasr_align_mm10.sbatch.sh
```

This generates:

proj.rds alignment project

alignment QC report (qAlign_QC_mm10.pdf)


4) Peak calling, merging, filtering, DE analysis

Performed in R using Science Apps - scripts will be added later


