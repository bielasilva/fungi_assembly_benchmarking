#!/bin/env bash

source /home/gas0042/miniforge3/bin/activate Rmarkdown

ROOTDIR="/scratch/gas0042/nanopore_benchmark"

Rscript /scratch/gas0042/nanopore_benchmark/real_data/scripts/tools/quast_busco_merger.R > logs/quast_merge.log 2>&1