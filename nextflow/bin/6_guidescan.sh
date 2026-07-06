#!/usr/bin/env bash
#
# =============================================================================
# Job Script
# =============================================================================
#
#SBATCH -J guidescan
#SBATCH --export=ALL
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --no-requeue
#SBATCH --partition=sfgf,biochem,wjg
#SBATCH --time=04:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=twang518@stanford.edu
#SBATCH -o guidescan-log.out
#SBATCH -e guidescan-log.err

<<COMMENT
7/2/26

Most-likely gRNA binding site per differential peak, via GuideScan2.
Runs 6A_guidescan.py on a BED of DiffBind peaks (see 5_DiffBind.sh output).
Edit the parameters below, then: sbatch 6_guidescan.sh

COMMENT

set -euo pipefail

# ---------------------------------------------------------------------------
# Customizable parameters
# ---------------------------------------------------------------------------
GRNA="CCCTACCTGTCACCAGGACC"        # protospacer (no PAM)
PAM="NGG"                          # SpCas9 PAM (N = wildcard)
MISMATCHES=10                      # max mismatches in the spacer
RNA_BULGES=2                       # max 1-nt RNA bulges (can stack adjacently up to X nt)
DNA_BULGES=2                       # max 1-nt DNA bulges (can stack adjacently up to X nt)
CHUNK_SIZE=15                      # peaks per parallel unit

SCRIPTDIR="/oak/stanford/groups/wjg/twang518/data/260618_beCasKAS_donor_revs/Scripts"
BED="${SCRIPTDIR}/ABE-diffbind_rep_cd7.bed"                 # DiffBind differential peaks
OUTFILE="${SCRIPTDIR}/cd7_binding_sites.tsv"                # output TSV (first 3 cols = BED)
GENOME="/oak/stanford/groups/wjg/twang518/software/reference_genome/hg38/hg38_Bowtie2Index/hg38.fasta"
GUIDESCAN="/oak/stanford/groups/wjg/twang518/data/260702_guide_align/beCasKAS/scratch/guidescan-cli/build/bin/guidescan"
# ---------------------------------------------------------------------------

# Guard: refuse to run until the gRNA is filled in.
if [[ "${GRNA}" == CHANGEME* || -z "${GRNA}" ]]; then
  echo "ERROR: set GRNA to the CD7 spacer sequence (20 nt, no PAM) before submitting." >&2
  exit 1
fi

# Sherlock environment: python + samtools + curl (libcurl needed by samtools/guidescan).
ml python/3.12.1 biology samtools/1.16.1 system curl/8.17.0
# guidescan was built with gcc/14.2.0; python/samtools pull gcc/12.4.0 (older libstdc++,
# missing GLIBCXX_3.4.32). Load gcc/14.2.0 LAST so its libstdc++ wins on LD_LIBRARY_PATH.
ml gcc/14.2.0

python3 "${SCRIPTDIR}/6A_guidescan.py" \
  --bed "${BED}" \
  --genome "${GENOME}" \
  --grna "${GRNA}" --pam "${PAM}" \
  --mismatches "${MISMATCHES}" \
  --rna-bulges "${RNA_BULGES}" --dna-bulges "${DNA_BULGES}" \
  --chunk-size "${CHUNK_SIZE}" --workers "${SLURM_CPUS_PER_TASK}" \
  --guidescan "${GUIDESCAN}" \
  --samtools "$(command -v samtools)" \
  --outfile "${OUTFILE}"

echo "=== OUTPUT: ${OUTFILE} ==="
column -t "${OUTFILE}"
