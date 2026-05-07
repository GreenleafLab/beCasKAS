#!/bin/bash
#
# Driver job for the beCasKAS Nextflow pipeline on a SLURM cluster.
#
# This script does NOT run the pipeline's compute. It runs the Nextflow
# controller, which then submits each pipeline process as its own SLURM job
# via `-profile slurm`. Keep the driver's resource request small: it spends
# most of its time idle, polling squeue and staging files.
#
# Per-process resource requests (CPUs, memory, walltime) are defined in
# nextflow/configs/process.config; do NOT try to size this script to match
# the heaviest process.
#
# EDIT BEFORE SUBMITTING -- the lines marked CUSTOMIZE must match your
# cluster's partitions, modules, and filesystem layout.

#SBATCH --job-name=beCasKAS
#SBATCH --partition=normal              # CUSTOMIZE: partition name on your cluster
#SBATCH --time=48:00:00                 # CUSTOMIZE: walltime ceiling for the whole run
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err

set -euo pipefail

# CUSTOMIZE: load the modules your cluster provides for these tools.
# The driver itself only needs Java + Nextflow; the container runtime
# (singularity/apptainer) is needed because compute jobs use it.
module load java
module load nextflow
module load singularity

# CUSTOMIZE: point this at a fast, large-quota filesystem (e.g. scratch).
# Nextflow writes work/ here; it grows quickly and benefits from parallel I/O.
RUN_DIR="${SCRATCH:-$PWD}/beCasKAS-runs/${SLURM_JOB_ID}"
mkdir -p "${RUN_DIR}"
cd "${RUN_DIR}"

# CUSTOMIZE: paths to the pipeline checkout, your samplesheet, and the
# reference directory built per the README.
PIPELINE_DIR="/path/to/beCasKAS-nf"
SAMPLESHEET="/path/to/samplesheet.csv"
REFERENCE_DIR="/path/to/refs/hg38"

nextflow run "${PIPELINE_DIR}/main.nf" \
    -profile slurm,singularity \
    --samplesheet "${SAMPLESHEET}" \
    --reference_dir "${REFERENCE_DIR}" \
    --outdir results \
    -resume
