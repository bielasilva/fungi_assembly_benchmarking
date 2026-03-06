#!/bin/bash

#SBATCH --job-name=Download_genomes
#SBATCH --output=logs/genomes_download.log
#SBATCH --ntasks=10
#SBATCH --time=2-00:00
#SBATCH --mem=6G
#SBATCH -p general,nova,nova_ff

ROOTDIR="/scratch/gas0042/nanopore_benchmark"

# Create necessary directories
mkdir -p "$ROOTDIR/genomes"

# Download assembly summary file
wget -q ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/fungi/assembly_summary.txt -O "$ROOTDIR/genomes/assembly_summary_$(date +%Y%m%d).txt"

# Download genomes from assembly summary file
awk -F '\t' '$12=="Complete Genome" {print $20}' $ROOTDIR/genomes/assembly_summary*.txt | while read -r url; do
    filename=$(basename "$url")
    fna_file="${filename}_genomic.fna.gz"

    echo "Downloading $fna_file from $url"
    wget -q $url/$fna_file -O $ROOTDIR/genomes/$fna_file
done

# Unzip downloaded genomes
find $ROOTDIR/genomes/ -name '*.fna.gz' -print0 | xargs -P 5 -0 -I{} -n1 sh -c 'gunzip -f "$1"' _ {}