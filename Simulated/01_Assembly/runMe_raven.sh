#!/bin/env bash

source /scratch/gas0042/nanopore_benchmark/scripts/sourceMe_configs.sh

load_env Raven

# Create output directories
if [[ ! -s logs/slurm_jobid/raven_slurm_jobid.csv ]]; then
    mkdir -p logs/raven
    mkdir -p logs/slurm_jobid
    echo -e "JobName,JobID,SubmissionTime" > logs/slurm_jobid/raven_slurm_jobid.csv
else
    cp logs/slurm_jobid/raven_slurm_jobid.csv logs/slurm_jobid/raven_slurm_jobid.csv.bk
fi

total_checks=$(( ${#SAMPLES[@]} * ${#depthX[@]} ))
echo "Total checks to consider: ${total_checks}"

# Loop through each sample and depth
for SAMPLE in "${SAMPLES[@]}"; do
    sample_number=$(( sample_number + 1 ))

    for DEPTH_NP in "${depthX[@]}"; do
        
        current_check=$(( current_check + 1 ))
        progress_bar "$current_check" "$total_checks" "$start_time" "$sample_number" "New jobs: $new_jobs"

        OUT_DIR="${ROOTDIR}/results/assemblies/${SAMPLE}/raven/${DEPTH_NP}"
        NANOPORE_FQ="${ROOTDIR}/simulated_data/subsampled/${SAMPLE}/${SAMPLE}_${DEPTH_NP}.nanopore.fq"
        
        # Check if it has been run before
        if [[ ! -s ${NANOPORE_FQ} ]]; then
            echo "raven for ${SAMPLE} at ${DEPTH_NP} cannot run because the input file does not exist. Skipping."
            continue
        elif [[ -s "${OUT_DIR}/raven_${SAMPLE}_NP${DEPTH_NP}.fasta" ]]; then
            # echo "raven for ${SAMPLE} at ${DEPTH_NP} already completed. Skipping."
            continue
        elif [[ -s "${ROOTDIR}/results/assemblies_results/${SAMPLE}/raven_${SAMPLE}_NP${DEPTH_NP}.fasta" ]]; then
            # echo "raven for ${SAMPLE} at ${DEPTH_NP} already completed. Skipping.
            continue
        else
            if squeue --me --format "%.100j" | grep -q "raven_${SAMPLE}_NP${DEPTH_NP}"; then
                # echo "raven for ${SAMPLE} at ${DEPTH_NP} is already running. Skipping."
                continue
            elif [[ $(squeue --me | wc -l) -ge 5000 ]]; then
                    echo_overwrite_2 "You have reached the maximum number of jobs (5000). Exiting."
                    exit 1
            else
                # echo "Submitting raven for ${SAMPLE} at ${DEPTH_NP}"
                new_jobs=$((new_jobs + 1))
                sed -i "/raven_${SAMPLE}_NP${DEPTH_NP}/d" logs/slurm_jobid/raven_slurm_jobid.csv
                rm -rf ${OUT_DIR}
                mkdir -p ${OUT_DIR}
            fi
        fi

sbatch <<- EOF | sed -e "s/Submitted batch job /raven_${SAMPLE}_NP${DEPTH_NP},/" -e "s/$/,$(date +'%Y-%m-%d %H:%M:%S')/" >> logs/slurm_jobid/raven_slurm_jobid.csv
#!/bin/env bash

#SBATCH --job-name=raven_${SAMPLE}_NP${DEPTH_NP}
#SBATCH --output=logs/raven/raven_${SAMPLE}_NP${DEPTH_NP}.out
#SBATCH --time=5-00
#SBATCH --ntasks=10
#SBATCH --mem=100GB
#SBATCH --partition=general
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=gabriel.silva@auburn.edu

cd ${OUT_DIR}

# Run raven
raven \
    --threads 10 \
    ${NANOPORE_FQ} > ${OUT_DIR}/raven_${SAMPLE}_NP${DEPTH_NP}.fasta

if [[ ! -s${OUT_DIR}/raven_${SAMPLE}_NP${DEPTH_NP}.fasta ]]; then
    echo "raven failed for ${SAMPLE} at ${DEPTH_NP}. Exiting."
    exit 1
fi
EOF
    done
done

echo_overwrite_2 "Done submitting Raven on all samples. Submitted ${new_jobs} new jobs."