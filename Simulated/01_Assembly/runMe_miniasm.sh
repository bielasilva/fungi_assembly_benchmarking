#!/bin/env bash

source /home/gas0042/miniforge3/bin/activate Miniasm

source /scratch/gas0042/nanopore_benchmark/scripts/sourceMe_configs.sh

# Create output directories
if [[ ! -s logs/slurm_jobid/miniasm_slurm_jobid.csv ]]; then
    mkdir -p logs/miniasm
    mkdir -p logs/slurm_jobid
    echo -e "JobName,JobID,SubmissionTime" > logs/slurm_jobid/miniasm_slurm_jobid.csv
else
    cp logs/slurm_jobid/miniasm_slurm_jobid.csv logs/slurm_jobid/miniasm_slurm_jobid.csv.bk
fi

total_checks=$(( ${#SAMPLES[@]} * ${#depthX[@]} ))
echo "Total checks to consider: ${total_checks}"

# Loop through each sample and depth
for SAMPLE in "${SAMPLES[@]}"; do
    sample_number=$(( sample_number + 1 ))
    # Get the genome size
    GENOME_FASTA="${ROOTDIR}/genomes/${SAMPLE}_genomic.fna"
    GENOME_SIZE=$(grep -v ">" "$GENOME_FASTA" | wc -c |  awk '{ printf "%.0f", $1/1e6 }')

    for DEPTH_NP in "${depthX[@]}"; do
        
        current_check=$(( current_check + 1 ))
        progress_bar "$current_check" "$total_checks" "$start_time" "$sample_number" "New jobs: $new_jobs"

        OUT_DIR="${ROOTDIR}/results/assemblies/${SAMPLE}/miniasm/${DEPTH_NP}"
        NANOPORE_FQ="${ROOTDIR}/simulated_data/subsampled/${SAMPLE}/${SAMPLE}_${DEPTH_NP}.nanopore.fq"
        
        # Check if it has been run before
        if [[ ! -s ${NANOPORE_FQ} ]]; then
            echo "miniasm for ${SAMPLE} at ${DEPTH_NP} cannot run because the input file does not exist. Skipping."
            continue
        elif [[ -s "${OUT_DIR}/miniasm_${SAMPLE}_NP${DEPTH_NP}.fasta" ]]; then
            # echo "miniasm for ${SAMPLE} at ${DEPTH_NP} already completed. Skipping."
            continue
        elif [[ -s "${ROOTDIR}/results/assemblies_results/${SAMPLE}/miniasm_${SAMPLE}_NP${DEPTH_NP}.fasta" ]]; then
            # echo "FLYE for ${SAMPLE} at ${DEPTH_NP} already completed. Skipping.
            continue
        else
            if squeue --me --format "%.100j" | grep -q "miniasm_${SAMPLE}_NP${DEPTH_NP}"; then
                # echo "miniasm for ${SAMPLE} at ${DEPTH_NP} is already running. Skipping."
                continue
            elif [[ $(squeue --me | wc -l) -ge 5000 ]]; then
                    echo_overwrite_2 "You have reached the maximum number of jobs (5000). Exiting."
                    exit 1
            else
                # echo "Submitting miniasm for ${SAMPLE} at ${DEPTH_NP}"
                new_jobs=$((new_jobs + 1))
                sed -i "/miniasm_${SAMPLE}_NP${DEPTH_NP}/d" logs/slurm_jobid/miniasm_slurm_jobid.csv
                rm -rf ${OUT_DIR}
                mkdir -p ${OUT_DIR}
            fi
        fi

sbatch <<- EOF | sed -e "s/Submitted batch job /miniasm_${SAMPLE}_NP${DEPTH_NP},/" -e "s/$/,$(date +'%Y-%m-%d %H:%M:%S')/" >> logs/slurm_jobid/miniasm_slurm_jobid.csv
#!/bin/env bash

#SBATCH --job-name=miniasm_${SAMPLE}_NP${DEPTH_NP}
#SBATCH --output=logs/miniasm/miniasm_${SAMPLE}_NP${DEPTH_NP}.out
#SBATCH --time=5-00
#SBATCH --ntasks=10
#SBATCH --mem=120GB
#SBATCH --partition=general
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=gabriel.silva@auburn.edu

# Run minimap2
minimap2 -x ava-ont -t 10 ${NANOPORE_FQ} ${NANOPORE_FQ} | gzip -1 > ${OUT_DIR}/miniasm_${SAMPLE}_NP${DEPTH_NP}.paf.gz

# Run miniasm
miniasm -f \
${NANOPORE_FQ} \
${OUT_DIR}/miniasm_${SAMPLE}_NP${DEPTH_NP}.paf.gz > ${OUT_DIR}/miniasm_${SAMPLE}_NP${DEPTH_NP}.gfa

# Convert GFA to FASTA
awk '/^S/{print ">"\$2"\n"\$3}' ${OUT_DIR}/miniasm_${SAMPLE}_NP${DEPTH_NP}.gfa > ${OUT_DIR}/miniasm_${SAMPLE}_NP${DEPTH_NP}.fasta

if [[ ! -s ${OUT_DIR}/miniasm_${SAMPLE}_NP${DEPTH_NP}.fasta ]]; then
    echo "miniasm failed for ${SAMPLE} at ${DEPTH_NP}. Exiting."
    exit 1
fi
EOF
    done
done

echo_overwrite_2 "Done submitting Miniasm on all samples. Submitted ${new_jobs} new jobs."