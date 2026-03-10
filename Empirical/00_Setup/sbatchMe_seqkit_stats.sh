#!/bin/env bash

#SBATCH --job-name=Seqkit_stats
#SBATCH --output=logs/seqkit_stats/seqkit_stats.out
#SBATCH --time=1-00
#SBATCH --ntasks=15
#SBATCH --mem=100GB
#SBATCH --partition=general,jrw0107_std,nova_20,nova_28,nova_ff
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=gabriel.silva@auburn.edu

source /scratch/gas0042/nanopore_benchmark/real_data/scripts/sourceMe_configs.sh

load_env Seqkit

ROOTDIR="/scratch/gas0042/nanopore_benchmark/real_data"

mkdir -p ${ROOTDIR}/results/seqkit_stats

seqkit stats -j 10 -T -a ${ROOTDIR}/data/subsampled/*/*.fq > ${ROOTDIR}/results/seqkit_stats/subsampled_seqkit_stats.tsv