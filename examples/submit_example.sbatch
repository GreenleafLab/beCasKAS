#!/bin/bash
#
# Submission job for the beCasKAS Nextflow pipeline on a SLURM cluster.
#
#SBATCH --job-name=beCasKAS
#SBATCH --partition=sfgf,wjg,biochem,normal   # Customize partitions
#SBATCH --time=48:00:00                       # Customize walltime
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err

set -euo pipefail

# Modify to load java, nextflow, and apptainer on your cluster
module load java
module load biology nextflow

RUN_DIR="${SCRATCH:-$PWD}/beCasKAS/${SLURM_JOB_ID}"
mkdir -p "${RUN_DIR}/logs"
cd "${RUN_DIR}"

# Modify to match your input/output structure
PIPELINE_DIR="/path/to/beCasKAS"
SAMPLESHEET="/path/to/samplesheet.csv"
REFERENCE_DIR="/path/to/refs/hg38"

nextflow run "${PIPELINE_DIR}/main.nf" \
    -profile slurm,apptainer \
    --samplesheet "${SAMPLESHEET}" \
    --refdir "${REFERENCE_DIR}" \
    --outdir results \
    -with-trace \
    -resume
