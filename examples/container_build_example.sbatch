#!/bin/bash
#SBATCH --job-name=becaskas-sif-build
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --partition=normal
#SBATCH --output=becaskas-sif-build-%j.out
#SBATCH --error=becaskas-sif-build-%j.err
set -euo pipefail

export APPTAINER_TMPDIR=$GROUP_SCRATCH/$USER/apptainer-tmp
export APPTAINER_CACHEDIR=$GROUP_SCRATCH/$USER/apptainer-cache
mkdir -p "$APPTAINER_TMPDIR" "$APPTAINER_CACHEDIR"

# Path to beCasKAS repository
cd $HOME/software/beCasKAS
apptainer build --fakeroot \
    "$GROUP_SCRATCH/$USER/becaskas.sif" \
    env/becaskas.def
