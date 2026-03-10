#!/bin/env bash

source /scratch/gas0042/nanopore_benchmark/real_data/scripts/sourceMe_configs.sh

load_env Flye

# Create output directories
if [[ ! -d logs/flye || ! -d logs/slurm_jobid ]]; then
    mkdir -p logs/flye
    mkdir -p logs/slurm_jobid
    echo -e "JobName,JobID,SubmissionTime" > logs/slurm_jobid/flye_slurm_jobid.csv
else
    cp logs/slurm_jobid/flye_slurm_jobid.csv logs/slurm_jobid/flye_slurm_jobid.csv.bk
fi

# OVERLAPS=("500" "1000" "1500" "2000" "2500" "3000")
OVERLAPS=("1000" "2000" "1500" "2500" "3000")

total_checks=$(( ${#SAMPLES_FULL[@]} * ${#depthX[@]} * ${#OVERLAPS[@]} ))
echo "Total checks to consider: ${total_checks}"

# Loop through each sample and depth
for SAMPLE in "${SAMPLES_FULL[@]}"; do
    sample_number=$(( sample_number + 1 ))
    for DEPTH_NP in "${depthX[@]}"; do
        for OVERLAP in "${OVERLAPS[@]}"; do

        current_check=$(( current_check + 1 ))
        progress_bar "$current_check" "$total_checks" "$start_time" "$sample_number" "New jobs: $new_jobs"

        SUBSAMPLE_DIR="${ROOTDIR}/data/subsampled/${SAMPLE}"
        NANOPORE_FQ="${SUBSAMPLE_DIR}/${SAMPLE}_${DEPTH_NP}.chopper.fq"

        OUT_DIR="${ROOTDIR}/results/assemblies/${SAMPLE}/flye-OVL${OVERLAP}/NP${DEPTH_NP}"

        # Check if it has been run before
        if [[ ! -s ${NANOPORE_FQ} ]]; then
            # echo "FLYE for ${SAMPLE} at ${DEPTH_NP} cannot run because the input file does not exist. Skipping."
            continue
        elif [[ -s "${OUT_DIR}/flye-OVL${OVERLAP}_${SAMPLE}_NP${DEPTH_NP}.fasta" ]]; then
            # echo "FLYE for ${SAMPLE} at ${DEPTH_NP} already completed. Skipping."
            continue
        elif [[ -s "${ROOTDIR}/results/assemblies_results/${SAMPLE}/flye-OVL${OVERLAP}_${SAMPLE}_NP${DEPTH_NP}.fasta" || -s "${ROOTDIR}/results/assemblies_results/${SAMPLE}/flye-OVL${OVERLAP}_${SAMPLE}_NP${DEPTH_NP}.fasta.gz" ]]; then
            # echo "FLYE for ${SAMPLE} at ${DEPTH_NP} already completed. Skipping."
            continue
        else
            if squeue --me --format "%.100j" | grep -q "flye-OVL${OVERLAP}_${SAMPLE}_NP${DEPTH_NP}"; then
                # echo "FLYE for ${SAMPLE} at ${DEPTH_NP} is already running. Skipping."
                continue
            elif [[ $(squeue --me | wc -l) -ge 5000 ]]; then
                    echo_overwrite_2 "You have reached the maximum number of jobs (5000). Exiting."
                    exit 1
            else
                new_jobs=$((new_jobs + 1))
                # echo "Submitting FLYE for ${SAMPLE} at ${DEPTH_NP}"
                sed -i "/flye-OVL${OVERLAP}_${SAMPLE}_NP${DEPTH_NP}/d" logs/slurm_jobid/flye_slurm_jobid.csv
                rm -rf ${OUT_DIR}
                mkdir -p ${OUT_DIR}
            fi
        fi

sbatch <<- EOF | sed -e "s/Submitted batch job /flye-OVL${OVERLAP}_${SAMPLE}_NP${DEPTH_NP},/" -e "s/$/,$(date +'%Y-%m-%d %H:%M:%S')/" >> logs/slurm_jobid/flye_slurm_jobid.csv
#!/bin/env bash

#SBATCH --job-name=flye-OVL${OVERLAP}_${SAMPLE}_NP${DEPTH_NP}
#SBATCH --output=logs/flye/flye-OVL${OVERLAP}_${SAMPLE}_NP${DEPTH_NP}.out
#SBATCH --time=5-00
#SBATCH --ntasks=10
#SBATCH --mem=280GB
#SBATCH --partition=bigmem2,bigmem4

#/SBATCH --mem=180GB
#/SBATCH --partition=general,jrw0107_std

#/]SBATCH --mail-type=FAIL
#/]SBATCH --mail-user=gabriel.silva@auburn.edu

# Run FLYE
flye \
    --threads 10 \
    --nano-hq ${NANOPORE_FQ} \
    --min-overlap ${OVERLAP} \
    --out-dir ${OUT_DIR}

if [[ ! -s ${OUT_DIR}/assembly.fasta ]]; then
    echo "FLYE failed for ${SAMPLE} at ${DEPTH_NP}. Exiting."
    exit 1
else
    ln -sf ${OUT_DIR}/assembly.fasta ${OUT_DIR}/flye-OVL${OVERLAP}_${SAMPLE}_NP${DEPTH_NP}.fasta
fi

EOF
        done
    done
done


echo_overwrite_2 "Done submitting Flye on all samples. Submitted ${new_jobs} new jobs."