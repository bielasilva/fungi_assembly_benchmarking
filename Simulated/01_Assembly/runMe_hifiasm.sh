#!/bin/env bash

source /home/gas0042/miniforge3/bin/activate Hifiasm

source /scratch/gas0042/nanopore_benchmark/scripts/sourceMe_configs.sh

# Create output directories
if [[ ! -d logs/hifiasm || ! -d logs/slurm_jobid ]]; then
    mkdir -p logs/hifiasm
    mkdir -p logs/slurm_jobid
    echo -e "JobName,JobID,SubmissionTime" > logs/slurm_jobid/hifiasm_slurm_jobid.csv
else
    cp logs/slurm_jobid/hifiasm_slurm_jobid.csv logs/slurm_jobid/hifiasm_slurm_jobid.csv.bk
fi

total_checks=$(( ${#SAMPLES[@]} * ${#depthX[@]} ))

echo "Total checks to consider: ${total_checks}"

# Loop through each sample and depth
for SAMPLE in "${SAMPLES[@]}"; do
    sample_number=$(( sample_number + 1 ))
    for DEPTH_NP in "${depthX[@]}"; do

        current_check=$(( current_check + 1 ))
        progress_bar "$current_check" "$total_checks" "$start_time" "$sample_number" "New jobs: $new_jobs"

        OUT_DIR="${ROOTDIR}/results/assemblies/${SAMPLE}/hifiasm/${DEPTH_NP}"

        SUBSAMPLE_DIR="${ROOTDIR}/simulated_data/subsampled/${SAMPLE}"
        NANOPORE_FQ="${SUBSAMPLE_DIR}/${SAMPLE}_${DEPTH_NP}.nanopore.fq"

        # Check if it has been run before
        if [[ ! -s ${NANOPORE_FQ} ]]; then
            echo "Hifiasm for ${SAMPLE} at ${DEPTH_NP} cannot run because the input file does not exist. Skipping."
            continue
        elif [[ -s "${OUT_DIR}/hifiasm_${SAMPLE}_NP${DEPTH_NP}.fasta" ]]; then
            # echo "Hifiasm for ${SAMPLE} at ${DEPTH_NP} already completed. Skipping."
            continue
        elif [[ -s "${ROOTDIR}/results/assemblies_results/${SAMPLE}/hifiasm_${SAMPLE}_NP${DEPTH_NP}.fasta" ]]; then
            # echo "Hifiasm for ${SAMPLE} at ${DEPTH_NP} already completed. Skipping.
            continue
        else
            if squeue --me --format "%.60j" | grep -q "hifiasm_${SAMPLE}_NP${DEPTH_NP}"; then
                # echo "Hifiasm for ${SAMPLE} at ${DEPTH_NP} is already running. Skipping."
                continue
            elif [[ $(squeue --me | wc -l) -ge 5000 ]]; then
                    echo_overwrite_2 "You have reached the maximum number of jobs (5000). Exiting."
                    exit 1
            else
                new_jobs=$((new_jobs + 1))
                # echo "Submitting Hifiasm for ${SAMPLE} at ${DEPTH_NP}"
                sed -i "/hifiasm_${SAMPLE}_NP${DEPTH_NP}/d" logs/slurm_jobid/hifiasm_slurm_jobid.csv
                rm -rf ${OUT_DIR}
                mkdir -p ${OUT_DIR}
            fi
        fi


sbatch <<- EOF | sed -e "s/Submitted batch job /hifiasm_${SAMPLE}_NP${DEPTH_NP},/" -e "s/$/,$(date +'%Y-%m-%d %H:%M:%S')/" >> logs/slurm_jobid/hifiasm_slurm_jobid.csv
#!/bin/env bash

#SBATCH --job-name=hifiasm_${SAMPLE}_NP${DEPTH_NP}
#SBATCH --output=logs/hifiasm/hifiasm_${SAMPLE}_NP${DEPTH_NP}_%j.out
#SBATCH --time=10-00
#SBATCH --ntasks=10
#SBATCH --mem=100GB
#SBATCH --partition=general
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=gabriel.silva@auburn.edu
#SBATCH --no-requeue

# Run Hifiasm
hifiasm \
    -t 10 \
    --ont \
    -o ${OUT_DIR}/hifiasm_${SAMPLE}_NP${DEPTH_NP} \
    ${NANOPORE_FQ}

if [[ ! -s ${OUT_DIR}/hifiasm_${SAMPLE}_NP${DEPTH_NP}.bp.p_utg.gfa ]]; then
    echo "Hifiasm failed for ${SAMPLE} at ${DEPTH_NP}. Exiting."
    exit 1
else
    awk '/^S/{print ">"\$2;print \$3}' ${OUT_DIR}/hifiasm_${SAMPLE}_NP${DEPTH_NP}.bp.p_utg.gfa > ${OUT_DIR}/hifiasm_${SAMPLE}_NP${DEPTH_NP}.fasta
fi
EOF
    done
done


echo_overwrite_2 "Done submitting Hifiasm on all samples. Submitted ${new_jobs} new jobs."