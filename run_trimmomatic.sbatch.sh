#!/usr/bin/env bash
#SBATCH --job-name=trimmomatic_only
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --output=logs/trimmomatic_%j.out
#SBATCH --error=logs/trimmomatic_%j.err

set -euo pipefail

# ---- Activate environment ----
eval "$(micromamba shell hook --shell bash)"
micromamba activate atac-trim    # or: trim_qc_env

# ---- IO ----
INDIR="/shares/domcke.dmls.uzh/external/data/Domcke_bulkATAC-seq"
OUTDIR="/home/mchere/scratch/GATA6_bulk-ATAC-seq_MCLW/outs"
mkdir -p logs "${OUTDIR}"

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
