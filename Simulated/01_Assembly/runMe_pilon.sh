#!/bin/env bash

source /home/gas0042/miniforge3/bin/activate Pilon

source /scratch/gas0042/nanopore_benchmark/scripts/sourceMe_configs.sh

# Create output directories
if [[ ! -d logs/pilon || ! -d logs/slurm_jobid ]]; then
    mkdir -p logs/pilon
    mkdir -p logs/slurm_jobid
    echo -e "JobName,JobID,SubmissionTime" > logs/slurm_jobid/pilon_slurm_jobid.csv
else 
    cp logs/slurm_jobid/pilon_slurm_jobid.csv logs/slurm_jobid/pilon_slurm_jobid.csv.bak
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

            FLYE_DIR="${ROOTDIR}/results/assemblies_results/${SAMPLE}/"
            FLYE_ASM="${FLYE_DIR}/flye_${SAMPLE}_NP${DEPTH_NP}.fasta"

            OUT_DIR="${ROOTDIR}/results/assemblies/${SAMPLE}/pilon/FLYE${DEPTH_NP}_IL${DEPTH_IL}"
            BWA_IDX="${OUT_DIR}/${SAMPLE}_FLYE${DEPTH_NP}_IL${DEPTH_IL}.bwa_idx"

            BWA_BAM="${OUT_DIR}/${SAMPLE}_FLYE${DEPTH_NP}_IL${DEPTH_IL}.sorted.bam"
            PILON_OUT="${OUT_DIR}/pilon_${SAMPLE}_FLYE${DEPTH_NP}_IL${DEPTH_IL}"
            
            current_check=$(( current_check + 1 ))
            progress_bar "$current_check" "$total_checks" "$start_time" "$sample_number" "New jobs: $new_jobs"

            if [[ ! -s "${FLYE_ASM}" ]]; then
                # echo "FLYE assembly for ${SAMPLE} not found. Skipping Pilon."
                continue
            elif [[ -s "${PILON_OUT}.fasta" ]]; then
                # echo_overwrite "Pilon for ${SAMPLE} at NP${DEPTH_NP}_IL${DEPTH_IL} is already completed. Skipping."
                continue
            elif [[ -s "${ROOTDIR}/results/assemblies_results/${SAMPLE}/pilon_${SAMPLE}_FLYE${DEPTH_NP}_IL${DEPTH_IL}.fasta" ]]; then
                # echo_overwrite "Pilon for ${SAMPLE} at NP${DEPTH_NP}_IL${DEPTH_IL} is already completed. Skipping."
                continue
            else
                if [[ $(squeue --me | wc -l) -ge 5000 ]]; then
                    echo_overwrite_2 "You have reached the maximum number of jobs (5000). Exiting."
                    exit 1
                elif squeue --me --format "%.100j" | grep -q "pilon_${SAMPLE}_FLYE${DEPTH_NP}_IL${DEPTH_IL}"; then
                    # echo_overwrite "pilon for ${SAMPLE} at NP${DEPTH_NP}_IL${DEPTH_IL} is already running. Skipping."
                    continue
                else
                    new_jobs=$((new_jobs + 1))
                    # echo_overwrite "Running pilon on ${SAMPLE}_FLYE${DEPTH_NP}_IL${DEPTH_IL}"
                    sed -i "/pilon_${SAMPLE}_FLYE${DEPTH_NP}_IL${DEPTH_IL}/d" logs/slurm_jobid/pilon_slurm_jobid.csv
                    rm -rf ${OUT_DIR}
                    mkdir -p ${OUT_DIR}
                fi
            fi

sbatch <<- EOF | sed -e "s/Submitted batch job /pilon_${SAMPLE}_FLYE${DEPTH_NP}_IL${DEPTH_IL},/g" -e "s/$/,$(date +'%Y-%m-%d %H:%M:%S')/" >> logs/slurm_jobid/pilon_slurm_jobid.csv
#!/bin/env bash

#SBATCH --job-name=pilon_${SAMPLE}_FLYE${DEPTH_NP}_IL${DEPTH_IL}
#SBATCH --output=logs/pilon/pilon_${SAMPLE}_FLYE${DEPTH_NP}_IL${DEPTH_IL}_%j.out
#SBATCH --time=5-00
#SBATCH --ntasks=10
#SBATCH --mem=50GB
#SBATCH --partition=general
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=gabriel.silva@auburn.edu
#SBATCH --no-requeue

# Index the FLYE assembly with BWA
bwa-mem2 index \
    -p ${BWA_IDX} \
    ${FLYE_ASM}

# Align the Illumina reads to the FLYE assembly and sort the BAM file
bwa-mem2 mem \
    -t 10 \
    ${BWA_IDX} \
    ${ILLUMINA1_FQ} \
    ${ILLUMINA2_FQ} \
    | samtools view -b - | samtools sort > ${BWA_BAM}

# Index the BAM file
samtools index ${BWA_BAM}

# Run Pilon to polish the FLYE assembly using the Illumina reads
pilon \
    --changes \
    --genome ${FLYE_ASM} \
    --frags ${BWA_BAM} \
    --output ${PILON_OUT}

if [[ ! -s ${PILON_OUT}.fasta ]]; then
    echo "Pilon failed for ${SAMPLE} at FLYE${DEPTH_NP}_IL${DEPTH_IL}. Exiting."
    exit 1
fi
EOF
        done
    done
done

echo_overwrite_2 "Done Running Pilon on all samples. Submitted ${new_jobs} new jobs."