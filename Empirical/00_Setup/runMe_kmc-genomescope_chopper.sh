#!/bin/bash

source /scratch/gas0042/nanopore_benchmark/real_data/scripts/sourceMe_configs.sh

#* Make sure that the fastq files do not have tabs or whitespaces in the headers 

if [[ ! -d logs/kmc-genomescope_chopperk57p1 || ! -d logs/slurm_jobid ]]; then
    mkdir -p logs/kmc-genomescope_chopperk57p1
    mkdir -p logs/slurm_jobid
    echo -e "JobName,JobID,SubmissionTime" > logs/slurm_jobid/kmc-genomescope_chopperk57p1_slurm_jobid.csv
else
    cp logs/slurm_jobid/kmc-genomescope_chopperk57p1_slurm_jobid.csv logs/slurm_jobid/kmc-genomescope_chopperk57p1_slurm_jobid.csv.bak
fi

for SAMPLE in "${SAMPLES_FULL[@]}"; do

    # Set the output and input
    OUTPUT_DIR="${ROOTDIR}/results/kmc-genomescope_chopperk57p1/${SAMPLE}"
    FASTQ_IN="${ROOTDIR}/data/subsampled/${SAMPLE}/${SAMPLE}_OG.chopper.fq"

    # Check if it has been run before
    if [[ -s "${OUTPUT_DIR}/summary.txt" ]]; then
        # echo "KMC and Genomescope for ${SAMPLE} already completed. Skipping."
        continue
    else
        if squeue --me --format "%.100j" | grep -q "kmc-genomescope_chopperk57p1_${SAMPLE}"; then
            # echo "canu for ${SAMPLE} at ${DEPTH_chopperk13p1} is already running. Skipping."
            continue
        elif [[ $(squeue --me | wc -l) -ge 5000 ]]; then
                echo_overwrite_2 "You have reached the maximum number of jobs (5000). Exiting."
                exit 1
        else
            # echo "Submitting canu for ${SAMPLE} at ${DEPTH_chopperk13p1}"
            sed -i "/kmc-genomescope_chopperk57p1_${SAMPLE}/d" logs/slurm_jobid/kmc-genomescope_chopperk57p1_slurm_jobid.csv
            mkdir -p ${OUTPUT_DIR}
        fi
    fi

sbatch <<- EOF | sed -e "s/Submitted batch job /kmc-genomescope_chopperk57p1_${SAMPLE},/" -e "s/$/,$(date +'%Y-%m-%d %H:%M:%S')/" >> logs/slurm_jobid/kmc-genomescope_chopperk57p1_slurm_jobid.csv
#!/bin/env bash

#SBATCH --job-name=kmc-genomescope_chopperk57p1_${SAMPLE}
#SBATCH --output=logs/kmc-genomescope_chopperk57p1/kmc-genomescope_chopperk57p1_${SAMPLE}.log
#SBATCH --ntasks=10
#SBATCH --time=1-00
#SBATCH --mem=50G
#SBATCH -p nova_20,nova_28,nova_ff
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=gabriel.silva@auburn.edu

ls $FASTQ_IN > ${OUTPUT_DIR}/readfile.fof

source /home/gas0042/miniforge3/bin/activate KMC
# Run
kmc \
    -k57 -m50g -t20 -ci1 \
    @${OUTPUT_DIR}/readfile.fof \
    ${OUTPUT_DIR}/${SAMPLE} \
    /tmp/

kmc_tools transform ${OUTPUT_DIR}/${SAMPLE} histogram ${OUTPUT_DIR}/${SAMPLE}.histo

source /home/gas0042/miniforge3/bin/activate Genomescope

genomescope2 -i ${OUTPUT_DIR}/${SAMPLE}.histo -o ${OUTPUT_DIR} -k 57 -p 1

EOF
done