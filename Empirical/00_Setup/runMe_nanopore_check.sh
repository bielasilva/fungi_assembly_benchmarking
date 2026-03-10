#!/bin/env bash

source /scratch/gas0042/nanopore_benchmark/real_data/scripts/sourceMe_configs.sh

for SAMPLE in "${SAMPLES_FULL[@]}"; do

    # echo "Checking sample ${SAMPLE} for nanopore data existence."
    # ls /scratch/acicola_project/fastqs/${SAMPLE}*

    # Check if there are any .fastq or .fastq.gz files for the sample
    if ls /scratch/acicola_project/fastqs/${SAMPLE}*.fastq.gz 1> /dev/null 2>&1; then
        echo "Nanopore data exists for sample ${SAMPLE}."
        # Count number of files
        FILE_COUNT=$(ls /scratch/acicola_project/fastqs/${SAMPLE}*.fastq.gz | wc -l)
        # Get biggest file
        BIGGEST_FILE=$(ls -S /scratch/acicola_project/fastqs/${SAMPLE}*.fastq.gz | head -n 1)
        NEW_DIR="/scratch/gas0042/nanopore_benchmark/real_data/data/filtered/original/${SAMPLE}"
        OUT_DIR="/scratch/gas0042/nanopore_benchmark/real_data/data/subsampled/${SAMPLE}"
        
        # cp "$BIGGEST_FILE" "${NEW_DIR}/${SAMPLE}_trimmed_nanopore.fastq.gz"
        # mv "${NEW_DIR}/${SAMPLE}_trimmed_nanopore.fastq" "${OUT_DIR}/${SAMPLE}_OG.nanopore.fq"
        # sed -i -E "s/\t|\s/\s/g" "${OUT_DIR}/${SAMPLE}_OG.nanopore.fq"
        
        gunzip -c "$BIGGEST_FILE" |
        sed -e "s/\t|\s/|/g" -e "s/400bps_sup@v5.0.0/400bps_sup_v5.0.0/g" > "${OUT_DIR}/${SAMPLE}_OG.nanopore.fq" # Replace tabs with pipes in FASTQ headers (for compatibility)

        # if [[ $FILE_COUNT -gt 1 ]]; then
        #     echo "Sample ${SAMPLE} has ${FILE_COUNT} nanopore files. Biggest file: ${BIGGEST_FILE} (${FILE_SIZE})."
        #     COMBINED_FILE=$(ls -S /scratch/acicola_project/fastqs/combined/${SAMPLE}*.fastq.gz | head -n 1)
        #     COMBINED_SIZE=$(du -h "$COMBINED_FILE" | cut -f1)
        #     echo "Combined nanopore file: ${COMBINED_FILE} (${COMBINED_SIZE})."
        #     echo ""
        # fi
    else
        echo "No nanopore data found for sample ${SAMPLE}."
    fi

done

# Not found
# BS01_360_046-X -> BS01_360_046-1
# BS16_CNT_131-4 -> BS16_CNT_131-x
# CAN7LA 
