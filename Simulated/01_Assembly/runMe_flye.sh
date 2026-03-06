#!/bin/env bash

source /home/gas0042/miniforge3/bin/activate Flye

source /scratch/gas0042/nanopore_benchmark/scripts/sourceMe_configs.sh

# Create output directories
if [[ ! -d logs/flye || ! -d logs/slurm_jobid ]]; then
    mkdir -p logs/flye
    mkdir -p logs/slurm_jobid
    echo -e "JobName,JobID,SubmissionTime" > logs/slurm_jobid/flye_slurm_jobid.csv
else
    cp logs/slurm_jobid/flye_slurm_jobid.csv logs/slurm_jobid/flye_slurm_jobid.csv.bk
fi

total_checks=$(( ${#SAMPLES[@]} * ${#depthX[@]} ))
echo "Total checks to consider: ${total_checks}"

# Loop through each sample and depth
for SAMPLE in "${SAMPLES[@]}"; do
    sample_number=$(( sample_number + 1 ))
    for DEPTH_NP in "${depthX[@]}"; do

        current_check=$(( current_check + 1 ))
        progress_bar "$current_check" "$total_checks" "$start_time" 60 "$sample_number" "New jobs: $new_jobs"

        OUT_DIR="${ROOTDIR}/results/assemblies/${SAMPLE}/flye/${DEPTH_NP}"
        NANOPORE_FQ="${ROOTDIR}/simulated_data/subsampled/${SAMPLE}/${SAMPLE}_${DEPTH_NP}.nanopore.fq"

        # Check if it has been run before
        if [[ ! -s ${NANOPORE_FQ} ]]; then
            echo "FLYE for ${SAMPLE} at ${DEPTH_NP} cannot run because the input file does not exist. Skipping."
            continue
        elif [[ -s "${OUT_DIR}/flye_${SAMPLE}_NP${DEPTH_NP}.fasta" ]]; then
            # echo "FLYE for ${SAMPLE} at ${DEPTH_NP} already completed. Skipping.
            continue
        elif [[ -s "${ROOTDIR}/results/assemblies_results/${SAMPLE}/flye_${SAMPLE}_NP${DEPTH_NP}.fasta" ]]; then
            # echo "FLYE for ${SAMPLE} at ${DEPTH_NP} already completed. Skipping.
            continue
        else
            if squeue --me --format "%.100j" | grep -q "flye_${SAMPLE}_NP${DEPTH_NP}"; then
                # echo "FLYE for ${SAMPLE} at ${DEPTH_NP} is already running. Skipping."
                continue
            elif [[ $(squeue --me | wc -l) -ge 5000 ]]; then
                    echo_overwrite_2 "You have reached the maximum number of jobs (5000). Exiting."
                    exit 1
            else
                new_jobs=$((new_jobs + 1))
                # echo "Submitting FLYE for ${SAMPLE} at ${DEPTH_NP}"
                sed -i "/flye_${SAMPLE}_NP${DEPTH_NP}/d" logs/slurm_jobid/flye_slurm_jobid.csv
                rm -rf ${OUT_DIR}
                mkdir -p ${OUT_DIR}
            fi
        fi

sbatch <<- EOF | sed -e "s/Submitted batch job /flye_${SAMPLE}_NP${DEPTH_NP},/" -e "s/$/,$(date +'%Y-%m-%d %H:%M:%S')/" >> logs/slurm_jobid/flye_slurm_jobid.csv
#!/bin/env bash

#SBATCH --job-name=flye_${SAMPLE}_NP${DEPTH_NP}
#SBATCH --output=logs/flye/flye_${SAMPLE}_NP${DEPTH_NP}_%j.out
#SBATCH --time=1-00
#SBATCH --ntasks=10
#SBATCH --mem=50GB
#SBATCH --partition=$(shuf -n 1 -e general nova nova_ff)
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=gabriel.silva@auburn.edu
#SBATCH --no-requeue

# Run FLYE
flye \
    --threads 10 \
    --nano-hq ${NANOPORE_FQ} \
    --out-dir ${OUT_DIR}

if [[ ! -s ${OUT_DIR}/assembly.fasta ]]; then
    echo "FLYE failed for ${SAMPLE} at ${DEPTH_NP}. Exiting."
    exit 1
else
    ln -sf ${OUT_DIR}/assembly.fasta ${OUT_DIR}/flye_${SAMPLE}_NP${DEPTH_NP}.fasta
fi

EOF
    done
done


echo_overwrite_2 "Done submitting Flye on all samples. Submitted ${new_jobs} new jobs."