#!/bin/env bash

source /home/gas0042/miniforge3/bin/activate Spades

source /scratch/gas0042/nanopore_benchmark/scripts/sourceMe_configs.sh

# Create output directories
if [[ ! -d logs/spades_hybrid || ! -d logs/slurm_jobid ]]; then
    mkdir -p logs/spades_hybrid
    mkdir -p logs/slurm_jobid
    echo -e "JobName,JobID,SubmissionTime" > logs/slurm_jobid/spades_hybrid_slurm_jobid.csv
else
    cp logs/slurm_jobid/spades_hybrid_slurm_jobid.csv logs/slurm_jobid/spades_hybrid_slurm_jobid.csv.bak
fi

total_checks=$(( ${#SAMPLES[@]} * ${#depthX_NP[@]} * ${#depthX_IL[@]} ))
echo "Total checks to consider: ${total_checks}"

# Loop through each sample and depth
for SAMPLE in "${SAMPLES[@]}"; do
    sample_number=$(( sample_number + 1 ))
    for DEPTH_NP in "${depthX_NP[@]}"; do
        for DEPTH_IL in "${depthX_IL[@]}"; do
        
            SUBSAMPLE_DIR="${ROOTDIR}/simulated_data/subsampled/${SAMPLE}"
            OUT_DIR="${ROOTDIR}/results/assemblies/${SAMPLE}/spades_hybrid/NP${DEPTH_NP}_IL${DEPTH_IL}"
            ILLUMINA1_FQ="${SUBSAMPLE_DIR}/${SAMPLE}_${DEPTH_IL}.R1.fq"
            ILLUMINA2_FQ="${SUBSAMPLE_DIR}/${SAMPLE}_${DEPTH_IL}.R2.fq"
            NANOPORE_FQ="${SUBSAMPLE_DIR}/${SAMPLE}_${DEPTH_NP}.nanopore.fq"

            current_check=$(( current_check + 1 ))
            progress_bar "$current_check" "$total_checks" "$start_time" "$sample_number" "New jobs: $new_jobs"

            # Check if it has been run before
            if [[ ! -s ${ILLUMINA1_FQ} || ! -s ${ILLUMINA2_FQ} || ! -s ${NANOPORE_FQ} ]]; then
                echo "SPADES for ${SAMPLE} at NP${DEPTH_NP}_IL${DEPTH_IL} cannot run because the input files do not exist. Skipping."
                continue
            elif [[ -s "${OUT_DIR}/spades_hybrid_${SAMPLE}_NP${DEPTH_NP}_IL${DEPTH_IL}.fasta" ]]; then
                # echo_overwrite "SPADES hybridfor ${SAMPLE} at NP${DEPTH_NP}_IL${DEPTH_IL} already completed. Skipping."
                continue
            elif [[ -s "${ROOTDIR}/results/assemblies_results/${SAMPLE}/spades_hybrid_${SAMPLE}_NP${DEPTH_NP}_IL${DEPTH_IL}.fasta" ]]; then
                # echo_overwrite "SPADES hybrid for ${SAMPLE} at NP${DEPTH_NP}_IL${DEPTH_IL} already completed. Skipping."
                continue
            else
                if squeue --me --format "%.100j" | grep -q spades_hybrid_${SAMPLE}_NP${DEPTH_NP}_IL${DEPTH_IL} ; then
                    # echo_overwrite "SPADES hybrid for ${SAMPLE} at NP${DEPTH_NP}_IL${DEPTH_IL} is already running. Skipping."
                    continue
                elif [[ $(squeue --me | wc -l) -ge 5000 ]]; then
                        echo_overwrite_2 "You have reached the maximum number of jobs (5000). Exiting."
                        exit 1
                else
                    new_jobs=$((new_jobs + 1))
                    sed -i "/spades_hybrid_${SAMPLE}_NP${DEPTH_NP}_IL${DEPTH_IL}/d" logs/slurm_jobid/spades_hybrid_slurm_jobid.csv
                    rm -rf ${OUT_DIR}
                    mkdir -p ${OUT_DIR}
                    # echo_overwrite_2 "Running SPADES Hybrid on ${SAMPLE} at NP${DEPTH_NP}_IL${DEPTH_IL}"
                fi
            fi

sbatch <<- EOF | sed -e "s/Submitted batch job /spades_hybrid_${SAMPLE}_NP${DEPTH_NP}_IL${DEPTH_IL},/g" -e "s/$/,$(date +'%Y-%m-%d %H:%M:%S')/" >> logs/slurm_jobid/spades_hybrid_slurm_jobid.csv
#!/bin/env bash

#SBATCH --job-name=spades_hybrid_${SAMPLE}_NP${DEPTH_NP}_IL${DEPTH_IL}
#SBATCH --output=logs/spades_hybrid/spades_hybrid_${SAMPLE}_NP${DEPTH_NP}_IL${DEPTH_IL}_%j.out
#SBATCH --time=1-00
#SBATCH --ntasks=10
#SBATCH --mem=100GB
#SBATCH --partition=general,jrw0107_std
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=gabriel.silva@auburn.edu
#SBATCH --no-requeue

# Run Spades Hybrid
spades.py \
    --threads 10 \
    -1 ${ILLUMINA1_FQ} \
    -2 ${ILLUMINA2_FQ} \
    --nanopore ${NANOPORE_FQ} \
    -o ${OUT_DIR}



if [[ ! -s ${OUT_DIR}/scaffolds.fasta ]]; then
    echo "SPADES Hybrid failed for ${SAMPLE} at NP${DEPTH_NP}_IL${DEPTH_IL}. Exiting."
    exit 1
else
    ln -sf ${OUT_DIR}/scaffolds.fasta ${OUT_DIR}/spades_hybrid_${SAMPLE}_NP${DEPTH_NP}_IL${DEPTH_IL}.fasta
fi
EOF
        done
    done
done

echo_overwrite_2 "Done submitting SPADES Hybrid for all samples. Submitted ${new_jobs} new jobs."
