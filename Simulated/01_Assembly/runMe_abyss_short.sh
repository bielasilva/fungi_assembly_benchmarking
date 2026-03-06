#!/bin/env bash

source /scratch/gas0042/nanopore_benchmark/scripts/sourceMe_configs.sh

load_env Abyss_man

# Create output directories
if [[ ! -d logs/abyss-short || ! -d logs/slurm_jobid ]]; then
    mkdir -p logs/abyss-short
    mkdir -p logs/slurm_jobid
    echo -e "JobName,JobID,SubmissionTime" > logs/slurm_jobid/abyss-short_slurm_jobid.csv
else
    cp logs/slurm_jobid/abyss-short_slurm_jobid.csv logs/slurm_jobid/abyss-short_slurm_jobid.csv.bak
fi

total_checks=$(( ${#SAMPLES[@]} * ${#depthX[@]} ))
echo "Total checks to consider: ${total_checks}"

for SAMPLE in "${SAMPLES[@]}"; do
    sample_number=$(( sample_number + 1 ))
    for DEPTH_IL in "${depthX[@]}"; do

    current_check=$(( current_check + 1 ))
    progress_bar "$current_check" "$total_checks" "$start_time" "$sample_number" "New jobs: $new_jobs"

    SUBSAMPLE_DIR="${ROOTDIR}/simulated_data/subsampled/${SAMPLE}"
    ILLUMINA1_FQ="${SUBSAMPLE_DIR}/${SAMPLE}_${DEPTH_IL}.R1.fq"
    ILLUMINA2_FQ="${SUBSAMPLE_DIR}/${SAMPLE}_${DEPTH_IL}.R2.fq"

    OUT_DIR="${ROOTDIR}/results/assemblies/${SAMPLE}/abyss-short/IL${DEPTH_IL}"

    # Check if it has been run before
    if [[ ! -s ${ILLUMINA1_FQ} || ! -s ${ILLUMINA2_FQ} ]]; then
        echo "Abyss Short for ${SAMPLE} at IL${DEPTH_IL} cannot run because the input files do not exist. Skipping."
        continue
    elif [[ -s "${OUT_DIR}/abyss-short_${SAMPLE}_IL${DEPTH_IL}.fasta" ]]; then
        continue
    elif [[ -s "${ROOTDIR}/results/assemblies_results/${SAMPLE}/abyss-short_${SAMPLE}_IL${DEPTH_IL}.fasta" ]]; then
        continue
    else
        if squeue --me --format "%.100j" | grep -q abyss-short_${SAMPLE}_IL${DEPTH_IL} ; then
            continue
        elif [[ $(squeue --me | wc -l) -ge 5000 ]]; then
                echo_overwrite_2 "You have reached the maximum number of jobs (5000). Exiting."
                exit 1
        else
            new_jobs=$((new_jobs + 1))
            sed -i "/abyss-short_${SAMPLE}_IL${DEPTH_IL}/d" logs/slurm_jobid/abyss-short_slurm_jobid.csv
            mkdir -p ${OUT_DIR}
        fi
    fi

sbatch <<- EOF | sed -e "s/Submitted batch job /abyss-short_${SAMPLE}_IL${DEPTH_IL},/g" -e "s/$/,$(date +'%Y-%m-%d %H:%M:%S')/" >> logs/slurm_jobid/abyss-short_slurm_jobid.csv
#!/bin/env bash

#SBATCH --job-name=abyss-short_${SAMPLE}_IL${DEPTH_IL}
#SBATCH --output=logs/abyss-short/abyss-short_${SAMPLE}_IL${DEPTH_IL}.out
#SBATCH --time=1-00
#SBATCH --ntasks=10
#SBATCH --mem=100GB
#SBATCH --partition=general
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=gabriel.silva@auburn.edu

# Run
abyss-pe -C ${OUT_DIR} \
    j=10 k=82 \
    name=abyss-short_${SAMPLE}_IL${DEPTH_IL} \
    in='${ILLUMINA1_FQ} ${ILLUMINA2_FQ}'

if [[ -s ${OUT_DIR}/abyss-short_${SAMPLE}_IL${DEPTH_IL}-unitigs.fa ]]; then
    mv ${OUT_DIR}/abyss-short_${SAMPLE}_IL${DEPTH_IL}-unitigs.fa ${OUT_DIR}/abyss-short_${SAMPLE}_IL${DEPTH_IL}.fasta
else
    echo "Abyss Short for ${SAMPLE} at IL${DEPTH_IL} failed. No output file generated."
    exit 1
fi

EOF
    done
done


echo_overwrite_2 "Done Running Abyss-short on all samples. Submitted ${new_jobs} new jobs."