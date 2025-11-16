#!/usr/bin/env bash
#SBATCH --job-name=trim_fastqc
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=logs/trim_fastqc_%j.out
#SBATCH --error=logs/trim_fastqc_%j.err

set -euo pipefail

# --- Settings (toggle if needed) ---
RUN_PRETRIM_QC=false          # set to true if you want FastQC on raw reads, too

# ---- Activate environment ----
eval "$(micromamba shell hook --shell bash)"
micromamba activate atac-trim   # environment must contain trimmomatic + fastqc

# ---- IO ----
INDIR="/shares/domcke.dmls.uzh/external/data/Domcke_bulkATAC-seq"
OUTDIR="/home/mchere/scratch/20251111_GATA6_bulk-ATAC-seq_MC_V2/trimmed_fastqs"
QCDIR="${OUTDIR}/qc"
mkdir -p logs "${OUTDIR}" "${QCDIR}" "${QCDIR}/trimmed"

# Create raw QC dir only if we plan to run pre-trim QC
if [[ "${RUN_PRETRIM_QC}" == "true" ]]; then
  mkdir -p "${QCDIR}/raw"
fi

# ---- Adapters from the environment ----
ADAPTERS="${CONDA_PREFIX}/share/trimmomatic/adapters/NexteraPE-PE.fa"
[[ -s "${ADAPTERS}" ]] || { echo "Adapter file not found: ${ADAPTERS}"; exit 1; }

# ---- Threads ----
THREADS="${SLURM_CPUS_PER_TASK:-8}"

# ---- Optional: FastQC on raw reads ----
if [[ "${RUN_PRETRIM_QC}" == "true" ]]; then
  echo "[INFO] Running FastQC on RAW reads..."
  fastqc -t "${THREADS}" --outdir "${QCDIR}/raw" "${INDIR}"/bulk*R{1,2}.fastq.gz
fi

# ---- Trimming + FastQC on trimmed mates ----
echo "[INFO] Trimming with Trimmomatic and running FastQC on trimmed reads..."
for f1 in "${INDIR}"/bulk*R1.fastq.gz; do
    [[ -e "$f1" ]] || { echo "No input files matched. Exiting."; exit 1; }

    f2="${f1/_R1/_R2}"
    [[ -e "$f2" ]] || { echo "Missing mate for ${f1} -> expected ${f2}"; exit 1; }

    echo "[INFO] Processing: $(basename "$f1") and $(basename "$f2")"

    b1="$(basename "$f1")"; b2="$(basename "$f2")"
    s1="${b1%.fastq.gz}";  s2="${b2%.fastq.gz}"

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

    # FastQC on the paired trimmed reads (skip unpaired)
    fastqc -t "${THREADS}" --outdir "${QCDIR}/trimmed" "$o1" "$o2"
done

echo "[DONE] Trimming + FastQC complete. Results:"
echo "  Trimmed FASTQs : ${OUTDIR}"
echo "  FastQC reports : ${QCDIR}/trimmed"
[[ "${RUN_PRETRIM_QC}" == "true" ]] && echo "  Raw FastQC     : ${QCDIR}/raw"
