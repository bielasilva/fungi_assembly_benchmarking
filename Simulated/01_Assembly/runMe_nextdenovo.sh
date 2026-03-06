#!/bin/env bash

source /home/gas0042/miniforge3/bin/activate Nextdenovo

source /scratch/gas0042/nanopore_benchmark/scripts/sourceMe_configs.sh

# Create output directories
if [[ ! -d logs/nextdenovo || ! -d logs/slurm_jobid ]]; then
    mkdir -p logs/nextdenovo
    mkdir -p logs/slurm_jobid
    echo -e "JobName,JobID,SubmissionTime" > logs/slurm_jobid/nextdenovo_slurm_jobid.csv
else
    cp logs/slurm_jobid/nextdenovo_slurm_jobid.csv logs/slurm_jobid/nextdenovo_slurm_jobid.csv.bk
fi

depthX=("20x" "25x" "30x" "35x" "40x" "50x" "60x" "75x" "100x")

total_checks=$(( ${#SAMPLES[@]} * ${#depthX[@]} ))
echo "Total checks to consider: ${total_checks}"

# Loop through each sample and depth
for SAMPLE in "${SAMPLES[@]}"; do
    sample_number=$(( sample_number + 1 ))
    GENOME_FASTA="${ROOTDIR}/genomes/${SAMPLE}_genomic.fna"
    GENOME_SIZE=$(grep -v ">" "$GENOME_FASTA" | wc -c |  awk '{ printf "%.0f", $1/1e6 }')
    for DEPTH_NP in "${depthX[@]}"; do

        current_check=$(( current_check + 1 ))
        progress_bar "$current_check" "$total_checks" "$start_time" 60 "$sample_number" "New jobs: $new_jobs"

        OUT_DIR="${ROOTDIR}/results/assemblies/${SAMPLE}/nextdenovo/${DEPTH_NP}"
        NANOPORE_FQ="${ROOTDIR}/simulated_data/subsampled/${SAMPLE}/${SAMPLE}_${DEPTH_NP}.nanopore.fq"

        # Check if it has been run before
        if [[ ! -s ${NANOPORE_FQ} ]]; then
            echo "FLYE for ${SAMPLE} at ${DEPTH_NP} cannot run because the input file does not exist. Skipping."
            continue
        elif [[ -s "${OUT_DIR}/nextdenovo_${SAMPLE}_NP${DEPTH_NP}.fasta" ]]; then
            # echo "FLYE for ${SAMPLE} at ${DEPTH_NP} already completed. Skipping.
            continue
        elif [[ -s "${ROOTDIR}/results/assemblies_results/${SAMPLE}/nextdenovo_${SAMPLE}_NP${DEPTH_NP}.fasta" ]]; then
            # echo "FLYE for ${SAMPLE} at ${DEPTH_NP} already completed. Skipping.
            continue
        else
            if squeue --me --format "%.100j" | grep -q "nextdenovo_${SAMPLE}_NP${DEPTH_NP}"; then
                # echo "FLYE for ${SAMPLE} at ${DEPTH_NP} is already running. Skipping."
                continue
            elif [[ $(squeue --me | wc -l) -ge 5000 ]]; then
                    echo_overwrite_2 "You have reached the maximum number of jobs (5000). Exiting."
                    exit 1
            else
                new_jobs=$((new_jobs + 1))
                # echo "Submitting FLYE for ${SAMPLE} at ${DEPTH_NP}"
                sed -i "/nextdenovo_${SAMPLE}_NP${DEPTH_NP}/d" logs/slurm_jobid/nextdenovo_slurm_jobid.csv
                rm -rf ${OUT_DIR}
                mkdir -p ${OUT_DIR}
            fi
        fi

sbatch <<- EOF | sed -e "s/Submitted batch job /nextdenovo_${SAMPLE}_NP${DEPTH_NP},/" -e "s/$/,$(date +'%Y-%m-%d %H:%M:%S')/" >> logs/slurm_jobid/nextdenovo_slurm_jobid.csv
#!/bin/env bash

#SBATCH --job-name=nextdenovo_${SAMPLE}_NP${DEPTH_NP}
#SBATCH --output=logs/nextdenovo/nextdenovo_${SAMPLE}_NP${DEPTH_NP}.out
#SBATCH --time=1-00
#SBATCH --ntasks=10
#SBATCH --mem=300GB
#SBATCH --partition=bigmem2,bigmem4
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=gabriel.silva@auburn.edu
#SBATCH --no-requeue


cd ${OUT_DIR}

ls ${NANOPORE_FQ} > input.fofn

cat > run.cfg << EOL
job_type = local
job_prefix = nextdenovo_${SAMPLE}_NP${DEPTH_NP}
task = all
rewrite = yes
deltmp = yes
parallel_jobs = 10
input_type = raw
read_type = ont
input_fofn = input.fofn
seed_depth = ${DEPTH_NP}
genome_size = ${GENOME_SIZE}m
EOL

# Run nextdenovo
nextDenovo run.cfg

if [[ ! -s ${OUT_DIR}/03.ctg_graph/nd.asm.fasta ]]; then
    echo "nextdenovo failed for ${SAMPLE} at ${DEPTH_NP}. Exiting."
    exit 1
else
    ln -sf ${OUT_DIR}/03.ctg_graph/nd.asm.fasta ${OUT_DIR}/nextdenovo_${SAMPLE}_NP${DEPTH_NP}.fasta
fi

EOF
    done
done


echo_overwrite_2 "Done submitting Nextdenovo on all samples. Submitted ${new_jobs} new jobs."