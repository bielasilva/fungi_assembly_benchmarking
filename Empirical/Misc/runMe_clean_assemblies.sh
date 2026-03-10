#!/bin/env bash
source /scratch/gas0042/nanopore_benchmark/real_data/scripts/sourceMe_configs.sh

total_checks=$(( ${#SAMPLES[@]} * ${#depthX[@]} * ${#depthX[@]} * ${#PROGRAMS[@]} ))

echo "Total checks to consider: ${total_checks}"

# rm /scratch/gas0042/nanopore_benchmark/scripts/busco_*.log

OVERLAPS=("500" "1000" "1500" "2000" "2500" "3000")

for SAMPLE in "${SAMPLES[@]}"; do
    sample_number=$(( sample_number + 1 ))
    mkdir -p ${ROOTDIR}/results/assemblies_results/${SAMPLE}

    for DEPTH_NP in "${depthX[@]}"; do
        for DEPTH_IL in "${depthX[@]}"; do
            for PROGRAM in "${PROGRAMS[@]}"; do
                current_check=$(( current_check + 1 ))
                progress_bar "$current_check" "$total_checks" "$start_time" "$sample_number" "Cleaned files: $new_jobs"

                case $PROGRAM in
                    "flye-OVL1000" | "flye-OVL1500" | "flye-OVL2000" | "flye-OVL2500" | "flye-OVL3000")
                        ASSEMBLY_DIR="${ROOTDIR}/results/assemblies/${SAMPLE}/${PROGRAM}/NP${DEPTH_NP}"
                        PREFIX="${PROGRAM}_${SAMPLE}_NP${DEPTH_NP}"
                        ASSEMBLY_FASTA="${ASSEMBLY_DIR}/${PREFIX}.fasta"
                        QUAST_DIR="${ROOTDIR}/results/quast/${SAMPLE}/${PREFIX}"
                        BUSCO_DIR="${ROOTDIR}/results/busco/${SAMPLE}/${PREFIX}"
                        BUSCO_JSON="${BUSCO_DIR}/short_summary.specific.fungi_odb10.${PREFIX}.json"
                        BUSCO_TXT="${BUSCO_DIR}/short_summary.specific.fungi_odb10.${PREFIX}.txt"
                        ;;
                    "flye" | "hifiasm" | "canu" | "raven" | "shasta" | "nextdenovo" | "miniasm")
                        ASSEMBLY_DIR="${ROOTDIR}/results/assemblies/${SAMPLE}/${PROGRAM}/NP${DEPTH_NP}"
                        PREFIX="${PROGRAM}_${SAMPLE}_NP${DEPTH_NP}"
                        ASSEMBLY_FASTA="${ASSEMBLY_DIR}/${PREFIX}.fasta"
                        QUAST_DIR="${ROOTDIR}/results/quast/${SAMPLE}/${PREFIX}"
                        BUSCO_DIR="${ROOTDIR}/results/busco/${SAMPLE}/${PREFIX}"
                        BUSCO_JSON="${BUSCO_DIR}/short_summary.specific.fungi_odb10.${PREFIX}.json"
                        BUSCO_TXT="${BUSCO_DIR}/short_summary.specific.fungi_odb10.${PREFIX}.txt"
                        ;;
                    "spades_short" | "abyss_short")
                        ASSEMBLY_DIR="${ROOTDIR}/results/assemblies/${SAMPLE}/${PROGRAM}/IL${DEPTH_IL}"
                        PREFIX="${PROGRAM}_${SAMPLE}_IL${DEPTH_IL}"
                        ASSEMBLY_FASTA="${ASSEMBLY_DIR}/${PREFIX}.fasta"
                        QUAST_DIR="${ROOTDIR}/results/quast/${SAMPLE}/${PREFIX}"
                        BUSCO_DIR="${ROOTDIR}/results/busco/${SAMPLE}/${PREFIX}"
                        BUSCO_JSON="${BUSCO_DIR}/short_summary.specific.fungi_odb10.${PREFIX}.json"
                        BUSCO_TXT="${BUSCO_DIR}/short_summary.specific.fungi_odb10.${PREFIX}.txt"
                        ;;
                    "spades_hybrid" | "masurca" | "polypolish" | "abyss-hybrid" | "pilon" | "racon" | "racon-mp2")
                        ASSEMBLY_DIR="${ROOTDIR}/results/assemblies/${SAMPLE}/${PROGRAM}/NP${DEPTH_NP}_IL${DEPTH_IL}"
                        PREFIX="${PROGRAM}_${SAMPLE}_NP${DEPTH_NP}_IL${DEPTH_IL}"
                        ASSEMBLY_FASTA="${ASSEMBLY_DIR}/${PREFIX}.fasta"
                        QUAST_DIR="${ROOTDIR}/results/quast/${SAMPLE}/${PREFIX}"
                        BUSCO_DIR="${ROOTDIR}/results/busco/${SAMPLE}/${PREFIX}"
                        BUSCO_JSON="${BUSCO_DIR}/short_summary.specific.fungi_odb10.${PREFIX}.json"
                        BUSCO_TXT="${BUSCO_DIR}/short_summary.specific.fungi_odb10.${PREFIX}.txt"
                        ;;
                    "reference")
                        ASSEMBLY_DIR=""
                        ASSEMBLY_FASTA=""
                        PREFIX="reference_${SAMPLE}"
                        QUAST_DIR=""
                        BUSCO_DIR="${ROOTDIR}/results/busco/${SAMPLE}/${PREFIX}"
                        BUSCO_JSON="${BUSCO_DIR}/short_summary.specific.fungi_odb10.${PREFIX}.json"
                        BUSCO_TXT="${BUSCO_DIR}/short_summary.specific.fungi_odb10.${PREFIX}.txt"
                        ;;
                esac

                # Check if the assembly file exists
                if [[ -s "${ASSEMBLY_FASTA}" ]]; then
                    if squeue --me --format "%.100j" | grep -q "^\s*${PREFIX}"; then
                        # echo_overwrite "Assembly for ${SAMPLE} with ${PREFIX} is still in progress. Skipping."
                        continue
                    fi
                    # echo_overwrite "Assembly for ${SAMPLE} with ${PREFIX} already exists. Copying."
                    cp ${ASSEMBLY_FASTA} ${ROOTDIR}/results/assemblies_results/${SAMPLE}/
                    # Make sure the copy was successful and files are the same
                    if cmp -s "${ASSEMBLY_FASTA}" "${ROOTDIR}/results/assemblies_results/${SAMPLE}/${PREFIX}.fasta"; then
                        # echo_overwrite "Assembly for ${SAMPLE} with ${PREFIX} copied successfully."
                        new_jobs=$((new_jobs + 1))
                        rm -rf $ASSEMBLY_DIR
                    else
                        # echo "Assembly for ${SAMPLE} with ${PREFIX} copy failed. Please check."
                        rm -f "${ROOTDIR}/results/assemblies_results/${SAMPLE}/${PREFIX}.fasta"
                    fi
                fi
                
                if [[ -s "${QUAST_DIR}/report.tsv" ]]; then
                    if squeue --me --format "%.100j" | grep -q "quast_${PREFIX}"; then
                        # echo_overwrite "Quast for ${PREFIX} is still in progress. Skipping."
                        continue
                    fi
                    # echo_overwrite "Quast for ${PREFIX} already completed. Copying."
                    mkdir -p ${ROOTDIR}/results/quast_results/${SAMPLE}
                    cp ${QUAST_DIR}/report.tsv ${ROOTDIR}/results/quast_results/${SAMPLE}/${PREFIX}_report.tsv
                    cp ${QUAST_DIR}/transposed_report.tsv ${ROOTDIR}/results/quast_results/${SAMPLE}/${PREFIX}_transposed_report.tsv

                    if cmp -s ${QUAST_DIR}/transposed_report.tsv ${ROOTDIR}/results/quast_results/${SAMPLE}/${PREFIX}_transposed_report.tsv; then
                        # echo_overwrite "Quast report copied successfully for ${PREFIX}."
                        new_jobs=$((new_jobs + 1))
                        rm -r $QUAST_DIR
                    else
                        # echo "Failed to copy Quast report for ${PREFIX}. Exiting."
                        rm -f ${ROOTDIR}/results/quast_results/${SAMPLE}/${PREFIX}_report.tsv
                        rm -f ${ROOTDIR}/results/quast_results/${SAMPLE}/${PREFIX}_transposed_report.tsv
                    fi
                fi

                # 
                if [[ -s "${BUSCO_JSON}" ]]; then
                    if squeue --me --format "%.100j" | grep -q "busco_${PREFIX}"; then
                        # echo "Busco for ${PREFIX} is still in progress. Skipping."
                        continue
                    fi
                    # echo_overwrite "Busco for ${PREFIX} already completed. Copying."
                    mkdir -p ${ROOTDIR}/results/busco_results/${SAMPLE}
                    cp ${BUSCO_JSON} ${ROOTDIR}/results/busco_results/${SAMPLE}/
                    cp ${BUSCO_TXT} ${ROOTDIR}/results/busco_results/${SAMPLE}/
                    if cmp -s ${BUSCO_JSON} ${ROOTDIR}/results/busco_results/${SAMPLE}/$(basename ${BUSCO_JSON}); then
                        # echo_overwrite "Busco results copied successfully for ${PREFIX}."
                        new_jobs=$((new_jobs + 1))
                        rm -rf $BUSCO_DIR
                    else
                        # echo "Failed to copy Busco results for ${PREFIX}. Exiting."
                        rm -f ${ROOTDIR}/results/busco_results/${SAMPLE}/$(basename ${BUSCO_JSON})
                        rm -f ${ROOTDIR}/results/busco_results/${SAMPLE}/$(basename ${BUSCO_TXT})
                    fi
                fi
            done
        done
    done
done


echo_overwrite_2 "Done cleanning the space. Cleaned ${new_jobs} files"