#!/usr/bin/env bash
#SBATCH --job-name=quasr_align_bulk_mm10
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=/home/mchere/scratch/20251111_GATA6_bulk-ATAC-seq_MC_V2/logs/quasr_align_bulk_mm10_%j.out
#SBATCH --error=/home/mchere/scratch/20251111_GATA6_bulk-ATAC-seq_MC_V2/logs/quasr_align_bulk_mm10_%j.err

set -euo pipefail

module load mamba
eval "$(micromamba shell hook --shell bash)"
#source /apps/opt/spack/linux-ubuntu20.04-x86_64/gcc-9.3.0/mamba-24.9.0-0-vjbtvwqisg5pgjnsdx6hv5beds3tw3jy/etc/profile.d/conda.sh
micromamba activate R_env

Rscript -e '
library(QuasR)
proj <- qAlign(
  sampleFile = "quasr_samples_full.tsv",
  genome     = "/shares/domcke.dmls.uzh/genomes/mm10/mm10.fa",
  aligner    = "Rbowtie",
  paired     = "fr"
)
qQCReport(proj, pdfFile = "/home/mchere/scratch/20251111_GATA6_bulk-ATAC-seq_MC_V2/qAlign_QC_mm10.pdf")
'
