#!/bin/env bash

#PBS -N Benchmark_R_empirical
#PBS -j oe
#PBS -o Benchmark_R_empirical.out
#PBS -l walltime=1:00:00
#PBS -l nodes=1
#PBS -l ncpus=10
#PBS -l mem=15gb
#PBS -q medium

source /home/aubgxs001/miniforge3/bin/activate Benchmark_R

ROOTDIR="/home/aubgxs001/scratch/github_fungi_assembly_benchmarking"

cd $ROOTDIR/Empirical/02_Processing/R

Rscript 00_orchestrator.R

