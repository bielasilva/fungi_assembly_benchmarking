#!/bin/env bash

source /scratch/gas0042/nanopore_benchmark/real_data/scripts/sourceMe_configs.sh

load_env Busco

# Create output directories
if [[ ! -d logs/busco || ! -d logs/slurm_jobid ]]; then
    mkdir -p logs/busco
    mkdir -p logs/slurm_jobid
    echo -e "JobName,JobID,SubmissionTime" > logs/slurm_jobid/busco_slurm_jobid.csv
else
    cp logs/slurm_jobid/busco_slurm_jobid.csv logs/slurm_jobid/busco_slurm_jobid.csv.bak
fi

total_checks=$(( ${#SAMPLES[@]} * ${#depthX[@]} * ${#depthX[@]} * ${#PROGRAMS[@]} ))
echo "Total checks to consider: ${total_checks}"

rm /scratch/gas0042/nanopore_benchmark/real_data/scripts/busco_*.log

for SAMPLE in "${SAMPLES[@]}"; do
    sample_number=$(( sample_number + 1 ))
    ASSEMBLY_DIR="${ROOTDIR}/results/assemblies_results/${SAMPLE}"
    for DEPTH_NP in "${depthX[@]}"; do
        for DEPTH_IL in "${depthX[@]}"; do
            for PROGRAM in "${PROGRAMS[@]}"; do

                current_check=$(( current_check + 1 ))
                progress_bar "$current_check" "$total_checks" "$start_time" "$sample_number" "New jobs: $new_jobs"
                
                case $PROGRAM in
                    "flye-OVL1000" | "flye-OVL1500" | "flye-OVL2000" | "flye-OVL2500" | "flye-OVL3000")
                        ASSEMBLY_FASTA="${ASSEMBLY_DIR}/${PROGRAM}_${SAMPLE}_NP${DEPTH_NP}.fasta"
                        PREFIX="${PROGRAM}_${SAMPLE}_NP${DEPTH_NP}"
                        ;;
                    "flye" | "hifiasm" | "canu" | "raven" | "shasta" | "nextdenovo" | "miniasm")
                        ASSEMBLY_FASTA="${ASSEMBLY_DIR}/${PROGRAM}_${SAMPLE}_NP${DEPTH_NP}.fasta"
                        PREFIX="${PROGRAM}_${SAMPLE}_NP${DEPTH_NP}"
                        ;;
                    "spades_short" | "abyss_short")
                        ASSEMBLY_FASTA="${ASSEMBLY_DIR}/${PROGRAM}_${SAMPLE}_IL${DEPTH_IL}.fasta"
                        PREFIX="${PROGRAM}_${SAMPLE}_IL${DEPTH_IL}"
                        ;; 
                    "spades_hybrid" | "masurca" | "polypolish" | "abyss-hybrid" | "racon" | "racon-mp2" | "pilon")
                        PREFIX="${PROGRAM}_${SAMPLE}_NP${DEPTH_NP}_IL${DEPTH_IL}"
                        ASSEMBLY_FASTA="${ASSEMBLY_DIR}/${PREFIX}.fasta"
                        ;;
                    "reference")
                        ASSEMBLY_FASTA="${ROOTDIR}/genomes/${SAMPLE}_genomic.fna"
                        PREFIX="reference_${SAMPLE}"
                        ;;
                esac

            # Check if the assembly file exists and if the job is already running
            if [[ ! -s "${ASSEMBLY_FASTA}" ]]; then
                # echo "No ${PREFIX} assembly found. Skipping"
                continue
            elif [[ -s "${ROOTDIR}/results/busco/${SAMPLE}/${PREFIX}/short_summary.specific.fungi_odb10.${PREFIX}.json" ]]; then
                # echo "busco for ${PREFIX} already completed. Skipping."
                continue
            elif [[ -s "${ROOTDIR}/results/busco_results/${SAMPLE}/short_summary.specific.fungi_odb10.${PREFIX}.json" ]]; then
                # echo "Busco for ${PREFIX} already completed and copied. Skipping."
                continue
            else
                if squeue --me --format "%.100j" | grep -q "busco_${PREFIX}"; then
                    # echo_overwrite "Busco for ${PREFIX} is already running. Skipping."
                    continue
                elif [[ $(squeue --me | wc -l) -ge 5000 ]]; then
                    echo_overwrite_2 "You have reached the maximum number of jobs (5000). Exiting."
                    exit 1
                else
                    new_jobs=$((new_jobs + 1))
                    # echo_overwrite "Submitting Busco for ${PREFIX}."
                    sed -i "/busco_${PREFIX}/d" logs/slurm_jobid/busco_slurm_jobid.csv
                    rm -rf ${ROOTDIR}/results/busco/${SAMPLE}/${PREFIX}
                fi
            fi

sbatch <<- EOF | sed -e "s/Submitted batch job /busco_${PREFIX},/" -e "s/$/,$(date +'%Y-%m-%d %H:%M:%S')/" >> logs/slurm_jobid/busco_slurm_jobid.csv
#!/bin/env bash

#SBATCH --job-name=busco_${PREFIX}
#SBATCH --output=logs/busco/busco_${PREFIX}.out
#SBATCH --time=1-00
#SBATCH --ntasks=10
#SBATCH --mem=50GB
#SBATCH --partition=nova_20,nova_28,nova_ff
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=gabriel.silva@auburn.edu

# Run busco
mkdir -p ${ROOTDIR}/results/busco/${SAMPLE}/${PREFIX}

busco \
    --in $ASSEMBLY_FASTA \
    --out ../results/busco/${SAMPLE}/${PREFIX} \
    --mode genome \
    --lineage_dataset /scratch/gas0042/nanopore_benchmark/busco_downloads/fungi_odb10 \
    --download_path /scratch/gas0042/nanopore_benchmark/busco_downloads/ \
    --offline \
    --cpu 10 --restart

if [[ ! -s ${ROOTDIR}/results/busco/${SAMPLE}/${PREFIX}/short_summary.specific.fungi_odb10.${PREFIX}.json ]]; then
    echo "Busco failed for ${PREFIX}. Exiting."
    exit 1
else
    echo "busco completed successfully for ${PREFIX}."
fi

EOF

            done
        done
    done
done

echo_overwrite_2 "Done running Busco on all samples. Submitted ${new_jobs} new jobs."