#!/bin/env bash

source /scratch/gas0042/nanopore_benchmark/scripts/sourceMe_configs.sh

load_env Abyss_man

# Create output directories
if [[ ! -d logs/abyss-hybrid || ! -d logs/slurm_jobid ]]; then
    mkdir -p logs/abyss-hybrid
    mkdir -p logs/slurm_jobid
    echo -e "JobName,JobID,SubmissionTime" > logs/slurm_jobid/abyss-hybrid_slurm_jobid.csv
else
    cp logs/slurm_jobid/abyss-hybrid_slurm_jobid.csv logs/slurm_jobid/abyss-hybrid_slurm_jobid.csv.bak
fi

total_checks=$(( ${#SAMPLES[@]} * ${#depthX_NP[@]} * ${#depthX_IL[@]} ))
echo "Total checks to consider: ${total_checks}"

for SAMPLE in "${SAMPLES[@]}"; do
    sample_number=$(( sample_number + 1 ))
    for DEPTH_NP in "${depthX_NP[@]}"; do
        for DEPTH_IL in "${depthX_IL[@]}"; do

        current_check=$(( current_check + 1 ))
        progress_bar "$current_check" "$total_checks" "$start_time" "$sample_number" "New jobs: $new_jobs"

        SUBSAMPLE_DIR="${ROOTDIR}/simulated_data/subsampled/${SAMPLE}"
        NANOPORE_FQ="${SUBSAMPLE_DIR}/${SAMPLE}_${DEPTH_NP}.nanopore.fq"
        ILLUMINA1_FQ="${SUBSAMPLE_DIR}/${SAMPLE}_${DEPTH_IL}.R1.fq"
        ILLUMINA2_FQ="${SUBSAMPLE_DIR}/${SAMPLE}_${DEPTH_IL}.R2.fq"

        OUT_DIR="${ROOTDIR}/results/assemblies/${SAMPLE}/abyss-hybrid/NP${DEPTH_NP}_IL${DEPTH_IL}"

        # Check if it has been run before
        if [[ ! -s ${ILLUMINA1_FQ} || ! -s ${ILLUMINA2_FQ} || ! -s ${NANOPORE_FQ} ]]; then
            echo "abyss-hybrid for ${SAMPLE} at NP${NANOPORE_FQ} IL${DEPTH_IL} cannot run because the input files do not exist. Skipping."
            continue
        elif [[ -s "${OUT_DIR}/abyss-hybrid_${SAMPLE}_NP${DEPTH_NP}_IL${DEPTH_IL}.fasta" ]]; then
            continue
        elif [[ -s "${ROOTDIR}/results/assemblies_results/${SAMPLE}/abyss-hybrid_${SAMPLE}_NP${DEPTH_NP}_IL${DEPTH_IL}.fasta" ]]; then
            continue
        else
            if squeue --me --format "%.100j" | grep -q abyss-hybrid_${SAMPLE}_NP${DEPTH_NP}_IL${DEPTH_IL} ; then
                continue
            elif [[ $(squeue --me | wc -l) -ge 5000 ]]; then
                    echo_overwrite_2 "You have reached the maximum number of jobs (5000). Exiting."
                    exit 1
            else
                new_jobs=$((new_jobs + 1))
                sed -i "/abyss-hybrid_${SAMPLE}_NP${DEPTH_NP}_IL${DEPTH_IL}/d" logs/slurm_jobid/abyss-hybrid_slurm_jobid.csv
                rm -rf ${OUT_DIR}
                mkdir -p ${OUT_DIR}
            fi
        fi

sbatch <<- EOF | sed -e "s/Submitted batch job /abyss-hybrid_${SAMPLE}_NP${DEPTH_NP}_IL${DEPTH_IL},/g" -e "s/$/,$(date +'%Y-%m-%d %H:%M:%S')/" >> logs/slurm_jobid/abyss-hybrid_slurm_jobid.csv
#!/bin/env bash

#SBATCH --job-name=abyss-hybrid_${SAMPLE}_NP${DEPTH_NP}_IL${DEPTH_IL}
#SBATCH --output=logs/abyss-hybrid/abyss-hybrid_${SAMPLE}_NP${DEPTH_NP}_IL${DEPTH_IL}.out
#SBATCH --time=20-00
#SBATCH --ntasks=10
#SBATCH --mem=100GB
#SBATCH --partition=general
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=gabriel.silva@auburn.edu

# Run
abyss-pe -C ${OUT_DIR} \
    j=10 k=82 \
    name=abyss-hybrid_${SAMPLE}_NP${DEPTH_NP}_IL${DEPTH_IL} \
    lib='pea' long='lra' \
    pea='${ILLUMINA1_FQ} ${ILLUMINA2_FQ}' \
    lra='${NANOPORE_FQ}'

if [[ -s ${OUT_DIR}/abyss-hybrid_${SAMPLE}_NP${DEPTH_NP}_IL${DEPTH_IL}-unitigs.fa ]]; then
    mv ${OUT_DIR}/abyss-hybrid_${SAMPLE}_NP${DEPTH_NP}_IL${DEPTH_IL}-unitigs.fa ${OUT_DIR}/abyss-hybrid_${SAMPLE}_NP${DEPTH_NP}_IL${DEPTH_IL}.fasta
else
    echo "Abyss Hybrid for ${SAMPLE} at IL${DEPTH_IL} failed. No output file generated."
    exit 1
fi

EOF
        done
    done
done

echo_overwrite_2 "Done Running Abyss-hybrid on all samples. Submitted ${new_jobs} new jobs."