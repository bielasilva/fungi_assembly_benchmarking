#!/bin/env bash
# Count how many assemblies (and their QC outputs) are completed or pending per assembler/config.


source /scratch/gas0042/nanopore_benchmark/real_data/scripts/sourceMe_configs.sh

# Pipelines to check
PROGRAMS=(
    "flye"
    "canu"
    "miniasm"
    "hifiasm"
    "nextdenovo"
    "spades_short"
    "abyss_short"
    "raven"
    "pilon"
    "spades_hybrid"
    "polypolish"
    "racon"
    "racon-mp2"
    "masurca"
    "shasta"
)

# Iterate per program and compute completed vs pending counts
for PROGRAM in "${PROGRAMS[@]}"; do
    echo_red "Checking program: ${PROGRAM}"
    case $PROGRAM in
        "flye" | "hifiasm" | "canu" | "raven" | "shasta" | "nextdenovo" | "miniasm")
            # Long-read-only assemblers; one axis per SAMPLE and ONT depth
            tot_done=0          # assemblies present
            to_do=0             # assemblies missing
            tot_done_quast=0    # QUAST reports present
            to_do_quast=0       # QUAST reports missing
            tot_done_busco=0    # BUSCO summaries present
            to_do_busco=0       # BUSCO summaries missing
            total_jobs=$((${#SAMPLES[@]} * ${#depthX[@]}))  # expected assemblies
            # Iterate over each SAMPLE and depth
            for SAMPLE in "${SAMPLES[@]}"; do
                for DEPTH_NP in "${depthX[@]}"; do
                    PREFIX="${PROGRAM}_${SAMPLE}_NP${DEPTH_NP}"
                    ASSEMBLY_FASTA="${ROOTDIR}/results/assemblies_results/${SAMPLE}/${PREFIX}.fasta"
                    ASSEMBLY_GZ="${ROOTDIR}/results/assemblies_results/${SAMPLE}/${PREFIX}.fasta.gz"
                    # Check if assembly exists and size > 0
                    if [[ -s $ASSEMBLY_FASTA || -s $ASSEMBLY_GZ ]]; then 
                        tot_done=$((tot_done + 1))
                        # Check QUAST output for this assembly
                        if [[ -s "${ROOTDIR}/results/quast_results/${SAMPLE}/${PREFIX}_report.tsv" ]]; then
                            tot_done_quast=$((tot_done_quast + 1))
                        else
                            to_do_quast=$((to_do_quast + 1))
                        fi
                        # Check BUSCO summary for this assembly
                        if [[ -s "${ROOTDIR}/results/busco_results/${SAMPLE}/short_summary.specific.fungi_odb10.${PREFIX}.json" ]]; then
                            tot_done_busco=$((tot_done_busco + 1))
                        else
                            to_do_busco=$((to_do_busco + 1))
                        fi
                    else
                        to_do=$((to_do + 1)) # assembly missing
                    fi
                done
            done
            ;;
        "spades_short" | "abyss_short" )
            # Short-read-only SPAdes; one axis per SAMPLE and Illumina depth
            tot_done=0          # assemblies present
            to_do=0             # assemblies missing
            tot_done_quast=0    # QUAST reports present
            to_do_quast=0       # QUAST reports missing
            tot_done_busco=0    # BUSCO summaries present
            to_do_busco=0       # BUSCO summaries missing
            total_jobs=$((${#SAMPLES[@]} * ${#depthX[@]}))  # expected assemblies
            # Iterate over each SAMPLE and depth
            for SAMPLE in "${SAMPLES[@]}"; do
                for DEPTH_IL in "${depthX[@]}"; do
                    PREFIX="${PROGRAM}_${SAMPLE}_IL${DEPTH_IL}"
                    ASSEMBLY_FASTA="${ROOTDIR}/results/assemblies_results/${SAMPLE}/${PREFIX}.fasta"
                    ASSEMBLY_GZ="${ROOTDIR}/results/assemblies_results/${SAMPLE}/${PREFIX}.fasta.gz"
                    # Check if assembly exists and size > 0
                    if [[ -s $ASSEMBLY_FASTA || -s $ASSEMBLY_GZ ]]; then
                        tot_done=$((tot_done + 1))
                        # Check QUAST output for this assembly
                        if [[ -s "${ROOTDIR}/results/quast_results/${SAMPLE}/${PREFIX}_report.tsv" ]]; then
                            tot_done_quast=$((tot_done_quast + 1))
                        else
                            to_do_quast=$((to_do_quast + 1))
                        fi
                        # Check BUSCO summary for this assembly
                        if [[ -s "${ROOTDIR}/results/busco_results/${SAMPLE}/short_summary.specific.fungi_odb10.${PREFIX}.json" ]]; then
                            tot_done_busco=$((tot_done_busco + 1))
                        else
                            to_do_busco=$((to_do_busco + 1))
                        fi
                    else
                        to_do=$((to_do + 1)) # assembly missing
                    fi
                done
            done
            ;;
        "spades_hybrid" | "masurca" | "polypolish" | "abyss-hybrid" | "racon" | "racon-mp2")
            # Hybrid SPAdes; cross-product of ONT and Illumina depths
            tot_done=0          # assemblies present
            to_do=0             # assemblies missing
            tot_done_quast=0    # QUAST reports present
            to_do_quast=0       # QUAST reports missing
            tot_done_busco=0    # BUSCO summaries present
            to_do_busco=0       # BUSCO summaries missing
            total_jobs=$((${#SAMPLES[@]} * ${#depthX_NP[@]} * ${#depthX_IL[@]})) # expected assemblies
            # Iterate over each SAMPLE, ONT depth, and Illumina depth
            for SAMPLE in "${SAMPLES[@]}"; do
                for DEPTH_NP in "${depthX_NP[@]}"; do
                    for DEPTH_IL in "${depthX_IL[@]}"; do
                        PREFIX="${PROGRAM}_${SAMPLE}_NP${DEPTH_NP}_IL${DEPTH_IL}"
                        ASSEMBLY_FASTA="${ROOTDIR}/results/assemblies_results/${SAMPLE}/${PREFIX}.fasta"
                        ASSEMBLY_GZ="${ROOTDIR}/results/assemblies_results/${SAMPLE}/${PREFIX}.fasta.gz"
                        # Check if assembly exists and size > 0
                        if [[ -s $ASSEMBLY_FASTA || -s $ASSEMBLY_GZ ]]; then
                            tot_done=$((tot_done + 1))
                            # Check QUAST output for this assembly
                            if [[ -s "${ROOTDIR}/results/quast_results/${SAMPLE}/${PREFIX}_report.tsv" ]]; then
                                tot_done_quast=$((tot_done_quast + 1))
                            else
                                to_do_quast=$((to_do_quast + 1))
                            fi
                            # Check BUSCO summary for this assembly
                            if [[ -s "${ROOTDIR}/results/busco_results/${SAMPLE}/short_summary.specific.fungi_odb10.${PREFIX}.json" ]]; then
                                tot_done_busco=$((tot_done_busco + 1))
                            else
                                to_do_busco=$((to_do_busco + 1))
                            fi
                        else
                            to_do=$((to_do + 1)) # assembly missing
                        fi
                    done
                done
            done
            ;;
        "pilon")
            # Pilon polishing of Flye assemblies using Illumina reads; cross-product of ONT and Illumina depths
            tot_done=0          # assemblies present
            to_do=0             # assemblies missing
            tot_done_quast=0    # QUAST reports present
            to_do_quast=0       # QUAST reports missing
            tot_done_busco=0    # BUSCO summaries present
            to_do_busco=0       # BUSCO summaries missing
            total_jobs=$((${#SAMPLES[@]} * ${#depthX_NP[@]} * ${#depthX_IL[@]})) # expected assemblies
            # Iterate over each SAMPLE, ONT depth, and Illumina depth
            for SAMPLE in "${SAMPLES[@]}"; do
                for DEPTH_NP in "${depthX_NP[@]}"; do
                    for DEPTH_IL in "${depthX_IL[@]}"; do
                        PREFIX="pilon_${SAMPLE}_FLYE${DEPTH_NP}_IL${DEPTH_IL}"
                        ASSEMBLY_FASTA="${ROOTDIR}/results/assemblies_results/${SAMPLE}/${PREFIX}.fasta"

                        if [[ -s $ASSEMBLY_FASTA ]]; then
                            tot_done=$((tot_done + 1))

                            if [[ -s "${ROOTDIR}/results/quast_results/${SAMPLE}/${PREFIX}_report.tsv" ]]; then
                                tot_done_quast=$((tot_done_quast + 1))
                            else
                                to_do_quast=$((to_do_quast + 1))
                            fi

                            if [[ -s "${ROOTDIR}/results/busco_results/${SAMPLE}/short_summary.specific.fungi_odb10.${PREFIX}.json" ]]; then
                                tot_done_busco=$((tot_done_busco + 1))
                            else
                                to_do_busco=$((to_do_busco + 1))
                            fi
                        else
                            to_do=$((to_do + 1))
                        fi
                    done
                done
            done
            ;;
        *)
            # Safety net for unexpected program names
            echo "Unknown program: ${PROGRAM}. Skipping."
            ;;
    esac
    
    if [[ $to_do -eq 0 && $to_do_quast -eq 0 && $to_do_busco -eq 0 ]]; then
        echo_green "All jobs for ${PROGRAM} are completed!"
        echo "--------------------------------------------------------------------------------"
    else
        if [[ $to_do -eq 0 ]]; then
            echo_green "All assemblies for ${PROGRAM} are completed!"
        else
            echo_yellow "${tot_done} completed, ${to_do} to do, ${total_jobs} total jobs."
        fi

        if [[ $to_do_quast -eq 0 ]]; then
            echo_green "No QUAST reports for ${PROGRAM} left at the moment! Done ${tot_done_quast}"
        else
            echo_yellow "Quast: ${tot_done_quast} completed, ${to_do_quast} to do."
        fi

        if [[ $to_do_busco -eq 0 ]]; then
            echo_green "No BUSCO reports for ${PROGRAM} left at the moment! Done ${tot_done_busco}"
        else
            echo_yellow "Busco: ${tot_done_busco} completed, ${to_do_busco} to do."
        fi
        echo "--------------------------------------------------------------------------------"
    fi
done
# Notes:
# - Paths assume a specific directory structure under ${ROOTDIR}/results.
# - QUAST and BUSCO checks rely on conventional output filenames.
# - Consider enabling 'set -Eeuo pipefail' for stricter error handling.
# - Consider quoting globs and using find instead of ls | sed to avoid edge cases.