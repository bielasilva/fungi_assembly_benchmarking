#!/bin/env bash

#PBS -N Benchmark_R_simulated
#PBS -j oe
#PBS -o Benchmark_R_simulated.out
#PBS -l walltime=5:00:00
#PBS -l nodes=1
#PBS -l ncpus=10
#PBS -l mem=50gb
#PBS -q large

source /home/aubgxs001/miniforge3/bin/activate Benchmark_R

ROOTDIR="/home/aubgxs001/scratch/github_fungi_assembly_benchmarking"

cd $ROOTDIR/Simulated/02_Processing/R

Rscript 00_orchestrator.R

