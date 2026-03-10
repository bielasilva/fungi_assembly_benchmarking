#!/usr/bin/env bash
set -euo pipefail

source /scratch/gas0042/nanopore_benchmark/real_data/scripts/sourceMe_configs.sh

OLD_DIR="/scratch/gas0042/nanopore_benchmark/real_data/data/filtered/original"
NEW_DIR="/scratch/gas0042/nanopore_benchmark/real_data/data/subsampled"

# How many parallel jobs (default: number of CPUs)
NPROC="20"

# Use pigz if available for faster decompression
DECOMP="gunzip -c"
if command -v pigz >/dev/null 2>&1; then
    DECOMP="pigz -dc"
fi

mkdir -p "$NEW_DIR"

printf '%s\0' "${SAMPLES_FULL[@]}" | xargs -0 -I{} -P "$NPROC" bash -c '
    set -euo pipefail
    OLD_DIR="$1"
    NEW_DIR="$2"
    SAMPLE="$3"

    OUT_DIR="${NEW_DIR}/${SAMPLE}"
    mkdir -p "$OUT_DIR"

    '"$DECOMP"' "${OLD_DIR}/${SAMPLE}/${SAMPLE}_trimmed_R1.fastq.gz" | sed -E "s/\s.*$/\/1/" > "${OUT_DIR}/${SAMPLE}_OG.R1.fq"
    '"$DECOMP"' "${OLD_DIR}/${SAMPLE}/${SAMPLE}_trimmed_R2.fastq.gz" | sed -E "s/\s.*$/\/2/" > "${OUT_DIR}/${SAMPLE}_OG.R2.fq"
' _ "$OLD_DIR" "$NEW_DIR" {}


