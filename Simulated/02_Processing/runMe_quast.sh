#!/bin/env bash

source /home/gas0042/miniforge3/bin/activate Quast

source /scratch/gas0042/nanopore_benchmark/scripts/sourceMe_configs.sh

# Create output directories
if [[ ! -d logs/quast || ! -d logs/slurm_jobid ]]; then
    mkdir -p logs/quast
    mkdir -p logs/slurm_jobid
    echo -e "JobName,JobID,SubmissionTime" > logs/slurm_jobid/quast_slurm_jobid.csv
else
    cp logs/slurm_jobid/quast_slurm_jobid.csv logs/slurm_jobid/quast_slurm_jobid.csv.bak
fi

total_checks=$(( ${#SAMPLES[@]} * ${#depthX[@]} * ${#depthX[@]} * ${#PROGRAMS[@]} ))

echo "Total checks to consider: ${total_checks}"

for SAMPLE in "${SAMPLES[@]}"; do
    sample_number=$(( sample_number + 1 ))
    ASSEMBLY_DIR="${ROOTDIR}/results/assemblies_results/${SAMPLE}"
    for DEPTH_NP in "${depthX[@]}"; do
        for DEPTH_IL in "${depthX[@]}"; do

            REF_GENOME="${ROOTDIR}/genomes/${SAMPLE}_genomic.fna"
            
            SUBSAMPLE_DIR="${ROOTDIR}/simulated_data/subsampled/${SAMPLE}"
            NANOPORE_FQ="${SUBSAMPLE_DIR}/${SAMPLE}_${DEPTH_NP}.nanopore.fq"
            ILLUMINA1_FQ="${SUBSAMPLE_DIR}/${SAMPLE}_${DEPTH_IL}.R1.fq"
            ILLUMINA2_FQ="${SUBSAMPLE_DIR}/${SAMPLE}_${DEPTH_IL}.R2.fq"

            for PROGRAM in "${PROGRAMS[@]}"; do

                current_check=$(( current_check + 1 ))
                progress_bar "$current_check" "$total_checks" "$start_time" "$sample_number" "New jobs: $new_jobs"
                
                case $PROGRAM in
                    "flye" | "hifiasm" | "canu" | "raven" | "shasta" | "nextdenovo" | "miniasm")
                        ASSEMBLY_FASTA="${ASSEMBLY_DIR}/${PROGRAM}_${SAMPLE}_NP${DEPTH_NP}.fasta"
                        PREFIX="${PROGRAM}_${SAMPLE}_NP${DEPTH_NP}"
                        # READS="--nanopore ${NANOPORE_FQ}"
                        ;;
                    "spades_short" | "abyss-short")
                        ASSEMBLY_FASTA="${ASSEMBLY_DIR}/${PROGRAM}_${SAMPLE}_IL${DEPTH_IL}.fasta"
                        PREFIX="${PROGRAM}_${SAMPLE}_IL${DEPTH_IL}"
                        # READS="--pe1 ${ILLUMINA1_FQ} --pe2 ${ILLUMINA2_FQ}"
                        ;;
                    "spades_hybrid" | "masurca" | "polypolish" | "abyss-hybrid" | "racon" | "racon-mp2")
                        PREFIX="${PROGRAM}_${SAMPLE}_NP${DEPTH_NP}_IL${DEPTH_IL}"
                        ASSEMBLY_FASTA="${ASSEMBLY_DIR}/${PREFIX}.fasta"
                        # READS="--pe1 ${ILLUMINA1_FQ} --pe2 ${ILLUMINA2_FQ} --nanopore ${NANOPORE_FQ}"
                        ;;
                    "pilon")
                        ASSEMBLY_FASTA="${ASSEMBLY_DIR}/pilon_${SAMPLE}_FLYE${DEPTH_NP}_IL${DEPTH_IL}.fasta"
                        PREFIX="pilon_${SAMPLE}_FLYE${DEPTH_NP}_IL${DEPTH_IL}"
                        # READS="--pe1 ${ILLUMINA1_FQ} --pe2 ${ILLUMINA2_FQ} --nanopore ${NANOPORE_FQ}"
                        ;;
                    "reference")
                        ASSEMBLY_FASTA="${ROOTDIR}/genomes/${SAMPLE}_genomic.fna"
                        PREFIX="reference_${SAMPLE}"
                        ;;
                esac

                OUT_DIR="${ROOTDIR}/results/quast/${SAMPLE}/${PREFIX}"
                # Check if the assembly file exists and if the job is already running
                if [[ ! -s "${ASSEMBLY_FASTA}" ]]; then
                    # echo_overwrite "No ${PREFIX} assembly found. Skipping"
                    continue
                elif [[ -s "${OUT_DIR}/report.tsv" ]]; then
                    # echo_overwrite "Quast for ${PREFIX} already completed. Copying."
                    continue
                elif [[ -s "${ROOTDIR}/results/quast_results/${SAMPLE}/${PREFIX}_report.tsv" ]]; then
                    # echo_overwrite "Quast for ${PREFIX} already copied. Skipping."
                    continue
                else
                    # if [[ $(squeue --me | wc -l) -ge 5000 ]]; then
                    #     echo_overwrite_2 "You have reached the maximum number of jobs (5000). Exiting."
                    #     exit 1
                    if squeue --me --format "%.100j" | grep -q "quast_${PREFIX}"; then
                        continue
                    else
                        new_jobs=$((new_jobs + 1))
                        sed -i "/quast_${PREFIX}/d" logs/slurm_jobid/quast_slurm_jobid.csv
                        # rm -rf ${OUT_DIR}
                    fi
                fi

sbatch <<- EOF | sed -e "s/Submitted batch job /quast_${PREFIX},/" -e "s/$/,$(date +'%Y-%m-%d %H:%M:%S')/" >> logs/slurm_jobid/quast_slurm_jobid.csv
#!/bin/env bash

#SBATCH --job-name=quast_${PREFIX}
#SBATCH --output=logs/quast/quast_${PREFIX}_%j.out
#SBATCH --time=1-00
#SBATCH --ntasks=5
#SBATCH --mem=20GB
#SBATCH --partition=nova_20,nova_28,nova_ff
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=gabriel.silva@auburn.edu

# Run QUAST
mkdir -p ${OUT_DIR}

quast.py \
    --threads 10 \
    --fungus \
    --min-contig 0 \
    -o ${OUT_DIR} \
    -r ${REF_GENOME} \
    ${READS} \
    --labels "${PREFIX}" \
    ${ASSEMBLY_FASTA}

if [[ ! -s ${OUT_DIR}/report.tsv ]]; then
    echo "Quast failed for ${PREFIX}. Exiting."
    exit 1
else
    echo "Quast completed successfully for ${PREFIX}."
    ln -sf ${OUT_DIR}/report.tsv ${OUT_DIR}/${PREFIX}_report.tsv
    ln -sf ${OUT_DIR}/transposed_report.tsv ${OUT_DIR}/${PREFIX}_transposed_report.tsv
fi
EOF

            done
        done
    done
done

echo_overwrite_2 "Done Running Quast on all samples. Submitted ${new_jobs} new jobs."