#!/bin/env bash

source /scratch/gas0042/nanopore_benchmark/scripts/sourceMe_configs.sh

# load_env Masurca

# export LD_LIBRARY_PATH="/tools/subread-2.0.6/lib:/usr/lib64:/usr/lib:/cm/shared/apps/slurm/current/lib64/slurm:/cm/shared/apps/slurm/current/lib64"
# export PATH="/tools/git-2.31.0/libexec:/tools/git-2.31.0/bin:/tools/subread-2.0.6/bin:/tools/subread-2.0.6:/cm/local/apps/environment-modules/4.2.1/bin:/usr/bin:/cm/shared/apps/slurm/current/sbin:/cm/shared/apps/slurm/current/bin:/usr/local/bin:/usr/local/sbin:/usr/sbin:/opt/dell/srvadmin/bin:/sbin:/usr/sbin:/cm/local/apps/environment-modules/4.2.1/bin"

# Create output directories
if [[ ! -d logs/masurca || ! -d logs/slurm_jobid ]]; then
    mkdir -p logs/masurca
    mkdir -p logs/slurm_jobid
    echo -e "JobName,JobID,SubmissionTime" > logs/slurm_jobid/masurca_slurm_jobid.csv
else
    cp logs/slurm_jobid/masurca_slurm_jobid.csv logs/slurm_jobid/masurca_slurm_jobid.csv.bak
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

            OUT_DIR="${ROOTDIR}/results/assemblies/${SAMPLE}/masurca/NP${DEPTH_NP}_IL${DEPTH_IL}"
            
            if [[ ! -s ${NANOPORE_FQ} || ! -s ${ILLUMINA1_FQ} || ! -s ${ILLUMINA2_FQ} ]]; then
                # echo_overwrite "At least one of the input files for ${SAMPLE} at NP${DEPTH_NP} IL${DEPTH_IL} is missing. Skipping."
                continue
            elif [[ -s "${OUT_DIR}/masurca_${SAMPLE}_NP${DEPTH_NP}_IL${DEPTH_IL}.fasta" ]]; then
                # echo_overwrite "Masurca for ${SAMPLE} at NP${DEPTH_NP}_IL${DEPTH_IL} already completed. Skipping."
                continue
            elif [[ -s "${ROOTDIR}/results/assemblies_results/${SAMPLE}/masurca_${SAMPLE}_NP${DEPTH_NP}_IL${DEPTH_IL}.fasta" ]]; then
                # echo_overwrite "Pilon for ${SAMPLE} already completed. Skipping."
                continue
            else
                if [[ $(squeue --me | wc -l) -ge 5000 ]]; then
                    echo_overwrite_2 "You have reached the maximum number of jobs (5000). Exiting."
                    exit 1
                # elif [[ $(squeue --me | grep -c masurca) -ge 10 ]]; then
                #     echo_overwrite_2 "More than 10 instances of masurca running. Exiting."
                #     exit 1
                elif squeue --me --format "%.100j" | grep -q "masurca_${SAMPLE}_NP${DEPTH_NP}_IL${DEPTH_IL}"; then
                    # echo_overwrite "Masurca for ${SAMPLE} NP${DEPTH_NP}_IL${DEPTH_IL}. Skipping."
                    continue
                else
                    new_jobs=$((new_jobs + 1))
                    # echo_overwrite "Submitting Masurca for ${SAMPLE} NP${DEPTH_NP}_IL${DEPTH_IL}."
                    sed -i "/masurca_${SAMPLE}_NP${DEPTH_NP}_IL${DEPTH_IL}/d" logs/slurm_jobid/masurca_slurm_jobid.csv
                    rm -rf ${OUT_DIR}
                fi
            fi

            mkdir -p ${OUT_DIR}

sbatch <<- EOF | sed -e "s/Submitted batch job /masurca_${SAMPLE}_NP${DEPTH_NP}_IL${DEPTH_IL},/g" -e "s/$/,$(date +'%Y-%m-%d %H:%M:%S')/" >> logs/slurm_jobid/masurca_slurm_jobid.csv
#!/bin/env bash

#SBATCH --job-name=masurca_${SAMPLE}_NP${DEPTH_NP}_IL${DEPTH_IL}
#SBATCH --output=logs/masurca/masurca_${SAMPLE}_NP${DEPTH_NP}_IL${DEPTH_IL}.out
#SBATCH --time=5-00
#SBATCH --cpus-per-task=10
#SBATCH --mem=600GB
#SBATCH --partition=bigmem4

cd ${OUT_DIR}

/home/gas0042/tools/MaSuRCA-4.1.4/bin/masurca -t 10 -i ${ILLUMINA1_FQ},${ILLUMINA2_FQ} -r ${NANOPORE_FQ}

if compgen -G "${OUT_DIR}/CA*/primary.genome.scf.fasta" > /dev/null; then
    echo "Masurca completed successfully for ${SAMPLE} NP${DEPTH_NP} IL${DEPTH_IL}"
    ln -s ${OUT_DIR}/CA*/primary.genome.scf.fasta ${OUT_DIR}/masurca_${SAMPLE}_NP${DEPTH_NP}_IL${DEPTH_IL}.fasta
else
    echo "Masurca failed for ${SAMPLE} NP${DEPTH_NP} IL${DEPTH_IL}"
    # rm *
    touch ${OUT_DIR}/masurca_${SAMPLE}_NP${DEPTH_NP}_IL${DEPTH_IL}.fail
    exit 1
fi

EOF
        done
    done
done

echo_overwrite_2 "Done submitting Masurca on all samples. Submitted ${new_jobs} new jobs."