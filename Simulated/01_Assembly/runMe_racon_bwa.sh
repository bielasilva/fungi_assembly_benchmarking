#!/bin/env bash

source /home/gas0042/miniforge3/bin/activate Racon

source /scratch/gas0042/nanopore_benchmark/scripts/sourceMe_configs.sh

# Create output directories
if [[ ! -d logs/racon || ! -d logs/slurm_jobid ]]; then
    mkdir -p logs/racon
    mkdir -p logs/slurm_jobid
    echo -e "JobName,JobID,SubmissionTime" > logs/slurm_jobid/racon_slurm_jobid.csv
else
    cp logs/slurm_jobid/racon_slurm_jobid.csv logs/slurm_jobid/racon_slurm_jobid.csv.bak
fi

total_checks=$(( ${#SAMPLES[@]} * ${#depthX_NP[@]} * ${#depthX_IL[@]} ))
echo "Total checks to consider: ${total_checks}"

# Loop through each sample and depth
for SAMPLE in "${SAMPLES[@]}"; do
    sample_number=$(( sample_number + 1 ))
    for DEPTH_NP in "${depthX_NP[@]}"; do
        for DEPTH_IL in "${depthX_IL[@]}"; do
    
            SUBSAMPLE_DIR="${ROOTDIR}/simulated_data/subsampled/${SAMPLE}"
            NANOPORE_FQ="${SUBSAMPLE_DIR}/${SAMPLE}_${DEPTH_NP}.nanopore.fq"
            ILLUMINA1_FQ="${SUBSAMPLE_DIR}/${SAMPLE}_${DEPTH_IL}.R1.fq"
            ILLUMINA2_FQ="${SUBSAMPLE_DIR}/${SAMPLE}_${DEPTH_IL}.R2.fq"

            FLYE_DIR="${ROOTDIR}/results/assemblies_results/${SAMPLE}"
            FLYE_ASM="${FLYE_DIR}/flye_${SAMPLE}_NP${DEPTH_NP}.fasta"

            OUT_DIR="${ROOTDIR}/results/assemblies/${SAMPLE}/racon/NP${DEPTH_NP}_IL${DEPTH_IL}"
            BWA_IDX="${OUT_DIR}/${SAMPLE}_NP${DEPTH_NP}_IL${DEPTH_IL}.bwa_idx"

            BWA_SAM="${OUT_DIR}/${SAMPLE}_NP${DEPTH_NP}_IL${DEPTH_IL}.sorted.sam"
            RACON_OUT="${OUT_DIR}/racon_${SAMPLE}_NP${DEPTH_NP}_IL${DEPTH_IL}"

            current_check=$(( current_check + 1 ))
            progress_bar "$current_check" "$total_checks" "$start_time" "$sample_number" "New jobs: $new_jobs"

            if [[ ! -s "${FLYE_ASM}" ]]; then
                # echo "FLYE assembly for ${SAMPLE} not found. Skipping Racon."
                continue
            elif [[ -s "${RACON_OUT}.fasta" ]]; then
                # echo_overwrite "Racon for ${SAMPLE} at NP${DEPTH_NP}_IL${DEPTH_IL} is already completed. Skipping."
                continue
            elif [[ -s "${ROOTDIR}/results/assemblies_results/${SAMPLE}/racon_${SAMPLE}_NP${DEPTH_NP}_IL${DEPTH_IL}.fasta" ]]; then
                # echo_overwrite "Racon for ${SAMPLE} at NP${DEPTH_NP}_IL${DEPTH_IL} is already completed. Skipping."
                continue
            else
                if [[ $(squeue --me | wc -l) -ge 5000 ]]; then
                    echo_overwrite_2 "You have reached the maximum number of jobs (5000). Exiting."
                    exit 1
                elif squeue --me --format "%.100j" | grep -q "racon_${SAMPLE}_NP${DEPTH_NP}_IL${DEPTH_IL}"; then
                    continue
                else
                    new_jobs=$((new_jobs + 1))
                    sed -i "/racon_${SAMPLE}_NP${DEPTH_NP}_IL${DEPTH_IL}/d" logs/slurm_jobid/racon_slurm_jobid.csv
                    rm -rf ${OUT_DIR}
                fi
            fi

            mkdir -p ${OUT_DIR}

sbatch <<- EOF | sed -e "s/Submitted batch job /racon_${SAMPLE}_NP${DEPTH_NP}_IL${DEPTH_IL},/g" -e "s/$/,$(date +'%Y-%m-%d %H:%M:%S')/" >> logs/slurm_jobid/racon_slurm_jobid.csv
#!/bin/env bash

#SBATCH --job-name=racon_${SAMPLE}_NP${DEPTH_NP}_IL${DEPTH_IL}
#SBATCH --output=logs/racon/racon_${SAMPLE}_NP${DEPTH_NP}_IL${DEPTH_IL}.out
#SBATCH --time=0-10
#SBATCH --ntasks=10
#SBATCH --mem=100GB
#SBATCH --partition=jrw0107_std,general
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=gabriel.silva@auburn.edu
#SBATCH --no-requeue

cd ${OUT_DIR}

cat $ILLUMINA1_FQ $ILLUMINA2_FQ | awk '{ sub(/\/[12]$/, "_" substr(\$0, length, 1)); print }' > ${SAMPLE}_${DEPTH_IL}.combined.fq

# Index the FLYE assembly with BWA
bwa-mem2 index \
    -p ${BWA_IDX} \
    ${FLYE_ASM}

# Align the Illumina reads to the FLYE assembly and sort the BAM file
bwa-mem2 mem \
    -t 10 \
    ${BWA_IDX} \
    ${SAMPLE}_${DEPTH_IL}.combined.fq > ${BWA_SAM}

# Run Racon to polish the FLYE assembly using the Illumina reads
racon \
    -t 10 \
    ${SAMPLE}_${DEPTH_IL}.combined.fq \
    ${BWA_SAM} \
    ${FLYE_ASM} \
    > ${RACON_OUT}.fasta

if [[ ! -s ${RACON_OUT}.fasta ]]; then
    echo "Racon failed for ${SAMPLE} at FLYE${DEPTH_NP}_IL${DEPTH_IL}. Exiting."
    exit 1
fi
EOF
        done
    done
done

echo_overwrite_2 "Done Running Racon on all samples. Submitted ${new_jobs} new jobs."