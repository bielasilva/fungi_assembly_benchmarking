#!/bin/env bash

source /scratch/gas0042/nanopore_benchmark/real_data/scripts/sourceMe_configs.sh

load_env Polypolish

# Create output directories
if [[ ! -d logs/polypolish || ! -d logs/slurm_jobid ]]; then
    mkdir -p logs/polypolish
    mkdir -p logs/slurm_jobid
    echo -e "JobName,JobID,SubmissionTime" > logs/slurm_jobid/polypolish_slurm_jobid.csv
else
    cp logs/slurm_jobid/polypolish_slurm_jobid.csv logs/slurm_jobid/polypolish_slurm_jobid.csv.bak
fi

total_checks=$(( ${#SAMPLES[@]} * ${#depthX_NP[@]} * ${#depthX_IL[@]} ))
echo "Total checks to consider: ${total_checks}"

# Loop through each sample and depth
for SAMPLE in "${SAMPLES[@]}"; do
    sample_number=$(( sample_number + 1 ))
    for DEPTH_NP in "${depthX_NP[@]}"; do
        for DEPTH_IL in "${depthX_IL[@]}"; do
    
            SUBSAMPLE_DIR="${ROOTDIR}/data/subsampled/${SAMPLE}"
            ILLUMINA1_FQ="${SUBSAMPLE_DIR}/${SAMPLE}_${DEPTH_IL}.R1.fq"
            ILLUMINA2_FQ="${SUBSAMPLE_DIR}/${SAMPLE}_${DEPTH_IL}.R2.fq"
            NANOPORE_FQ="${SUBSAMPLE_DIR}/${SAMPLE}_${DEPTH_NP}.chopper.fq"

            FLYE_DIR="${ROOTDIR}/results/assemblies_results/${SAMPLE}/"
            FLYE_ASM="${FLYE_DIR}/flye-OVL1000_${SAMPLE}_NP${DEPTH_NP}.fasta"

            OUT_DIR="${ROOTDIR}/results/assemblies/${SAMPLE}/polypolish/NP${DEPTH_NP}_IL${DEPTH_IL}"
            BWA_IDX="${OUT_DIR}/${SAMPLE}_NP${DEPTH_NP}_IL${DEPTH_IL}.bwa_idx"

            BWA_SAM1="${OUT_DIR}/${SAMPLE}_NP${DEPTH_NP}_IL${DEPTH_IL}.sorted1.sam"
            BWA_SAM2="${OUT_DIR}/${SAMPLE}_NP${DEPTH_NP}_IL${DEPTH_IL}.sorted2.sam"
            POLYPOLISH_OUT="${OUT_DIR}/polypolish_${SAMPLE}_NP${DEPTH_NP}_IL${DEPTH_IL}"

            current_check=$(( current_check + 1 ))
            progress_bar "$current_check" "$total_checks" "$start_time" "$sample_number" "New jobs: $new_jobs"

            if [[ ! -s "${FLYE_ASM}" ]]; then
                # echo "FLYE assembly for ${SAMPLE} not found. Skipping Polypolish."
                continue
            elif [[ -s "${POLYPOLISH_OUT}.fasta" ]]; then
                # echo_overwrite "Polypolish for ${SAMPLE} at NP${DEPTH_NP}_IL${DEPTH_IL} is already completed. Skipping."
                continue
            elif [[ -s "${ROOTDIR}/results/assemblies_results/${SAMPLE}/polypolish_${SAMPLE}_NP${DEPTH_NP}_IL${DEPTH_IL}.fasta" || -s "${ROOTDIR}/results/assemblies_results/${SAMPLE}/polypolish_${SAMPLE}_NP${DEPTH_NP}_IL${DEPTH_IL}.fasta.gz" ]]; then
                # echo_overwrite "Polypolish for ${SAMPLE} at NP${DEPTH_NP}_IL${DEPTH_IL} is already completed. Skipping."
                continue
            else
                if [[ $(squeue --me | wc -l) -ge 5000 ]]; then
                    echo_overwrite_2 "You have reached the maximum number of jobs (5000). Exiting."
                    exit 1
                elif squeue --me --format "%.100j" | grep -q "polypolish_${SAMPLE}_NP${DEPTH_NP}_IL${DEPTH_IL}"; then
                    # echo_overwrite "Polypolish for ${SAMPLE} at NP${DEPTH_NP}_IL${DEPTH_IL} is already running. Skipping."
                    continue
                else
                    new_jobs=$((new_jobs + 1))
                    # echo_overwrite "Running Polypolish on ${SAMPLE}_NP${DEPTH_NP}_IL${DEPTH_IL}"
                    sed -i "/polypolish_${SAMPLE}_NP${DEPTH_NP}_IL${DEPTH_IL}/d" logs/slurm_jobid/polypolish_slurm_jobid.csv
                    rm -rf ${OUT_DIR}
                    # continue
                fi
            fi

            mkdir -p ${OUT_DIR}

sbatch <<- EOF | sed -e "s/Submitted batch job /polypolish_${SAMPLE}_NP${DEPTH_NP}_IL${DEPTH_IL},/g" -e "s/$/,$(date +'%Y-%m-%d %H:%M:%S')/" >> logs/slurm_jobid/polypolish_slurm_jobid.csv
#!/bin/env bash

#SBATCH --job-name=polypolish_${SAMPLE}_NP${DEPTH_NP}_IL${DEPTH_IL}
#SBATCH --output=logs/polypolish/polypolish_${SAMPLE}_NP${DEPTH_NP}_IL${DEPTH_IL}.out
#SBATCH --time=5-00
#SBATCH --ntasks=10
#SBATCH --mem=180GB
#SBATCH --partition=general,jrw0107_std
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=gabriel.silva@auburn.edu

cd ${OUT_DIR}

# Index the FLYE assembly with BWA
bwa-mem2 index \
    -p ${BWA_IDX} \
    ${FLYE_ASM}

# Align the Illumina reads to the FLYE assembly and sort the BAM file
bwa-mem2 mem \
    -a \
    -t 10 \
    ${BWA_IDX} \
    ${ILLUMINA1_FQ} \
    > ${BWA_SAM1}

bwa-mem2 mem \
    -a \
    -t 10 \
    ${BWA_IDX} \
    ${ILLUMINA2_FQ} \
    > ${BWA_SAM2}

# Run Polypolish to polish the FLYE assembly using the Illumina reads
polypolish filter --in1 ${BWA_SAM1} --in2 ${BWA_SAM2} --out1 filtered_1.sam --out2 filtered_2.sam

polypolish polish ${FLYE_ASM} filtered_1.sam filtered_2.sam > ${POLYPOLISH_OUT}.fasta

if [[ ! -s ${POLYPOLISH_OUT}.fasta ]]; then
    echo "Polypolish failed for ${SAMPLE} at NP${DEPTH_NP}_IL${DEPTH_IL}. Exiting."
    exit 1
fi
EOF
        done
    done
done

echo_overwrite_2 "Done Running Polypolish on all samples. Submitted ${new_jobs} new jobs."