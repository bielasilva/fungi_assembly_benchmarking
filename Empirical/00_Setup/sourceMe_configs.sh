#!/bin/env bash

RED='\033[1;31m'
GREEN='\033[1;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

echo_red() {
    echo -e "${RED}$1${NC}"
}
echo_green() {
    echo -e "${GREEN}$1${NC}"
}
echo_yellow() {
    echo -e "${YELLOW}$1${NC}"
}

echo_overwrite() {
    echo -ne "\r\033[2K$1"
}

echo_overwrite_2() {
    echo -ne "\r\033[2K$1\033[B\r\033[2K"
}

load_env () {
    source /home/gas0042/miniforge3/bin/activate $1
}
# -------- Progress bar helpers (ADD THIS BLOCK) --------
fmt_time() {  # HH:MM:SS
    local T=$1
    printf "%02d:%02d:%02d" $((T/3600)) $((T%3600/60)) $((T%60))
}

progress_bar() {  # progress_bar <current> <total> <start_epoch> [sample_number] [extra_message]
    local current=$1 total=$2 start=$3 sample_number=${4:-0} extra_message=${5:-""}
    local max_width=$(tput cols)
    local width=$(( 60 * max_width / 100 ))
    local percent=$(( current * 100 / total ))
    local filled=$(( percent * width / 100 ))
    local empty=$(( width - filled ))
    local elapsed=$(( $(date +%s) - start ))

    local _fill=$(printf "%${filled}s")
    local _empty=$(printf "%${empty}s")

    printf "\r\033[2KSample $sample_number | $extra_message \n"
    printf "\033[2K[${_fill// /=}${_empty// / }] %3d%% (%d/%d) | elapsed %s \033[A" \
        "$percent" "$current" "$total" "$(fmt_time "$elapsed")"
}

start_timer() {
    local __var=$1
    printf -v "$__var" '%s' "$(date +%s)"
}

current_check=0
sample_number=0
new_jobs=0

start_time=$(date +%s)

tput civis 2>/dev/null || true
trap 'tput cnorm 2>/dev/null || true; echo' EXIT
# ------------------------------------------------------

ROOTDIR="/scratch/gas0042/nanopore_benchmark/real_data"

# Get the list of SAMPLEs
SAMPLES_FULL=($(ls ${ROOTDIR}/data/raw/*_R1_001.fastq.gz | sed 's/.*\///; s/_S.*_R1_001.fastq.gz//'))

readarray -t SAMPLES < "${ROOTDIR}/data/samples_list.txt"

# SAMPLES=($(ls ${ROOTDIR}/data/raw/*_R1_001.fastq.gz | sed 's/.*\///; s/_S.*_R1_001.fastq.gz//'))

# Depth in NNx format
depthX=("10x" "15x" "20x" "25x" "30x" "35x" "40x" "50x" "60x")
depthX_NP=("10x" "15x" "20x" "25x" "30x" "35x" "40x" "50x" "60x")
depthX_IL=("10x" "15x" "20x" "25x" "30x" "35x" "40x" "50x" "60x")

# depthX=("OG" "10x" "15x" "20x" "25x" "30x" "35x" "40x" "50x" "60x")
# depthX_NP=("OG" "10x" "15x" "20x" "25x" "30x" "35x" "40x" "50x" "60x")
# depthX_IL=("OG" "10x" "15x" "20x" "25x" "30x" "35x" "40x" "50x" "60x")

# depthX=("OG")
# depthX_NP=("OG")
# depthX_IL=("OG")

ERROR_PROGRAMS=("masurca" "hifiasm" "shasta")

COMPLETED_PROGRAMS=(
    # "nextdenovo"
    # "flye"
    # "spades_short"
    # "spades_hybrid"
    # "abyss_short"
    # "pilon"
    # "raven"
    # "racon"
    # "reference"
)

PROGRAMS=(
    "spades_short"
    "abyss_short"
    "spades_hybrid"
    "polypolish"
    "flye-OVL1000"
    "flye-OVL1500"
    "flye-OVL2000"
    "flye-OVL2500"
    "flye-OVL3000"
)