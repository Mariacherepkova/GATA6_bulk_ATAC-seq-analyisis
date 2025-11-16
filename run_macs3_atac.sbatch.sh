#!/bin/bash
#SBATCH -J macs3_atac
#SBATCH -p standard                # or: lowprio
#SBATCH -t 12:00:00                # walltime; adjust if needed
#SBATCH -c 4                       # CPU cores
#SBATCH --mem=32G                  # RAM per node
#SBATCH -o logs/macs3_%x_%j.out
#SBATCH -e logs/macs3_%x_%j.err

# --- environment ---

set -euo pipefail

module load mamba
eval "$(micromamba shell hook --shell bash)"
micromamba activate macs3-env

# --- paths ---
PROJ_DIR="$PWD"                    # submit from project root
BAM_DIR="${PROJ_DIR}/bams"
OUT_ROOT="${PROJ_DIR}/peaks"
LOG_DIR="${PROJ_DIR}/logs"

mkdir -p "${OUT_ROOT}" "${LOG_DIR}"

# --- sanity checks ---
command -v macs3 >/dev/null || { echo "MACS3 not found in env"; exit 1; }
[ -d "${BAM_DIR}" ] || { echo "Missing BAM dir: ${BAM_DIR}"; exit 1; }

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

  # Optional: produce a narrowPeak BED sorted and indexed
  if [ -f "${OUTDIR}/${SAMPLE}_peaks.narrowPeak" ]; then
    sort -k1,1 -k2,2n "${OUTDIR}/${SAMPLE}_peaks.narrowPeak" > "${OUTDIR}/${SAMPLE}.narrowPeak.sorted.bed"
    bgzip -f "${OUTDIR}/${SAMPLE}.narrowPeak.sorted.bed"
    tabix -p bed "${OUTDIR}/${SAMPLE}.narrowPeak.sorted.bed.gz"
  fi

done

echo "[$(date)] All MACS3 peak calls finished."
