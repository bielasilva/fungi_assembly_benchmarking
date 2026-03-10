#!/bin/env bash

source /scratch/gas0042/nanopore_benchmark/real_data/scripts/sourceMe_configs.sh

load_env Chopper

ROOTDIR="/scratch/gas0042/nanopore_benchmark/real_data"

# Create output directories
if [[ ! -d logs/chopper || ! -d logs/slurm_jobid ]]; then
    mkdir -p logs/chopper
    mkdir -p logs/slurm_jobid
    echo -e "JobName,JobID,SubmissionTime" > logs/slurm_jobid/chopper_slurm_jobid.csv
else
    cp logs/slurm_jobid/chopper_slurm_jobid.csv logs/slurm_jobid/chopper_slurm_jobid.csv.bak
fi


for SAMPLE in "${SAMPLES_FULL[@]}"; do
    
    SUBSAMPLE_DIR="${ROOTDIR}/data/subsampled/${SAMPLE}"
    NANOPORE_FQ="${SUBSAMPLE_DIR}/${SAMPLE}_OG.nanopore.fq"

    mkdir -p ${ROOTDIR}/results/chopper/${SAMPLE}

    if [[ ! -f ${NANOPORE_FQ} ]]; then
        echo "One or more input files for sample ${SAMPLE} do not exist. Skipping..."
        continue
    elif [[ -f ${SUBSAMPLE_DIR}/${SAMPLE}_OG.chopper.fq ]]; then
        echo "Chopper output file for sample ${SAMPLE} already exist. Skipping..."
        continue
    fi

sbatch <<- EOF | sed -e "s/Submitted batch job /chopper_${SAMPLE}/" -e "s/$/,$(date +'%Y-%m-%d %H:%M:%S')/" >> logs/slurm_jobid/chopper_slurm_jobid.csv
#!/bin/env bash

#SBATCH --job-name=chopper_${SAMPLE}
#SBATCH --output=logs/chopper/chopper_${SAMPLE}.out
#SBATCH --time=1-00
#SBATCH --ntasks=10
#SBATCH --mem=50GB
#SBATCH --partition=general,jrw0107_std,nova_20,nova_28,nova_ff
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=gabriel.silva@auburn.edu

chopper --threads 10 --minlength 1000 -i ${NANOPORE_FQ} > ${SUBSAMPLE_DIR}/${SAMPLE}_OG.chopper.fq

if [ \$? -ne 0 ]; then
    echo "chopper failed for sample ${SAMPLE}"
    exit 1
fi
EOF

done