# GATA6_bulk_ATAC-seq-analyisis

This repository contains the workflow and scripts used for processing bulk ATAC-seq data for the GATA6 perturbation experiment. The pipeline consists of read trimming and QC, reference indexing, read alignment with QuasR/Rbowtie, peak calling, and differential accessibility analysis.

All steps are run on the DMLS cluster.  
Environments are managed with **micromamba** and R project reproducibility is controlled with **renv**.

---

## Environment Setup

Before starting to set up the environments request the interactive session on cluster

```bash
srun --pty -n 1 -c 2 --time=02:00:00 --mem=32G bash -l
```

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

Files location 

```bash
/shares/domcke.dmls.uzh/external/data/Domcke_bulkATAC-seq
```


1) Trim reads
```bash
sbatch scripts/run_trimmomatic.sbatch.sh
```

*This script also does FastQC before and after trimming:*
```bash
sbatch run_trimmomatic_fastqc.sbatch.sh 
```

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


4) Peak calling,


5) merging, filtering, DE analysis

Performed in R using Science Apps - scripts will be added later

## Specific parameters used

1) For Trimmomatic:

```bash
# ---- Adapters from the environment ----
ADAPTERS="${CONDA_PREFIX}/share/trimmomatic/adapters/NexteraPE-PE.fa"

# ---- Threads ----
THREADS="${SLURM_CPUS_PER_TASK:-8}"

# ---- Run Trimmomatic ----
for f1 in "${INDIR}"/bulk*R1.fastq.gz; do
    [[ -e "$f1" ]] || { echo "No input files matched. Exiting."; exit 1; }

    echo "Processing ${f1}"
    f2="${f1/_R1/_R2}"

    # Safe basenames and stems
    b1="$(basename "$f1")"
    b2="$(basename "$f2")"
    s1="${b1%.fastq.gz}"
    s2="${b2%.fastq.gz}"

    # Output files
    o1="${OUTDIR}/${s1}.trim.fastq.gz"
    u1="${OUTDIR}/${s1}.trim.unpaired.fastq.gz"
    o2="${OUTDIR}/${s2}.trim.fastq.gz"
    u2="${OUTDIR}/${s2}.trim.unpaired.fastq.gz"

    trimmomatic PE -threads "${THREADS}" \
        "$f1" "$f2" \
        "$o1" "$u1" \
        "$o2" "$u2" \
        ILLUMINACLIP:"${ADAPTERS}":2:30:10:1:true \
        TRAILING:3 SLIDINGWINDOW:4:10 MINLEN:20
done
```


```bash

```




2) For alignment:

```bash
Rscript -e '
library(QuasR)
proj <- qAlign(
  sampleFile = "quasr_samples_full.tsv",
  genome     = "/shares/domcke.dmls.uzh/genomes/mm10/mm10.fa",
  aligner    = "Rbowtie",
  paired     = "fr"
)
qQCReport(proj, pdfFile = "/home/mchere/scratch/GATA6_bulk-ATAC-seq_MCLW/qAlign_QC_mm10.pdf")
'
```

*Genome used:*

```bash
/shares/domcke.dmls.uzh/genomes/mm10
```


2) For peak calling:

```bash
# --- MACS3 settings for ATAC-seq (paired-end) ---
# For paired-end ATAC, use -f BAMPE and no shifting (MACS uses fragment spans).
# Genome: mm10 -> use -g mm
QVAL=0.01

# Loop over BAMs and call peaks individually
for BAM in "${BAM_DIR}"/*.bam; do
  # skip index files just in case
  [[ "${BAM}" == *.bam ]] || continue

  SAMPLE="$(basename "${BAM}" .bam)"
  OUTDIR="${OUT_ROOT}/${SAMPLE}"
  mkdir -p "${OUTDIR}"

  echo "[$(date)] Calling peaks for ${SAMPLE}"
  macs3 callpeak \
    -t "${BAM}" \
    -f BAMPE \
    -g mm \
    -n "${SAMPLE}" \
    --outdir "${OUTDIR}" \
    --keep-dup all \
    -q "${QVAL}"

  # Produce a narrowPeak BED sorted and indexed
  if [ -f "${OUTDIR}/${SAMPLE}_peaks.narrowPeak" ]; then
    sort -k1,1 -k2,2n "${OUTDIR}/${SAMPLE}_peaks.narrowPeak" > "${OUTDIR}/${SAMPLE}.narrowPeak.sorted.bed"
    bgzip -f "${OUTDIR}/${SAMPLE}.narrowPeak.sorted.bed"
    tabix -p bed "${OUTDIR}/${SAMPLE}.narrowPeak.sorted.bed.gz"
  fi
```
