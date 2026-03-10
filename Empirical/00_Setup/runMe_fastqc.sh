#!/bin/env bash

source /scratch/gas0042/nanopore_benchmark/real_data/scripts/sourceMe_configs.sh

load_env Fastqc

ROOTDIR="/scratch/gas0042/nanopore_benchmark/real_data"

# Create output directories
if [[ ! -d logs/fastqc || ! -d logs/slurm_jobid ]]; then
    mkdir -p logs/fastqc
    mkdir -p logs/slurm_jobid
    echo -e "JobName,JobID,SubmissionTime" > logs/slurm_jobid/fastqc_slurm_jobid.csv
else
    cp logs/slurm_jobid/fastqc_slurm_jobid.csv logs/slurm_jobid/fastqc_slurm_jobid.csv.bak
fi

DEPTH_NP="OG"
DEPTH_IL="OG"

SAMPLES_FULL=("BS00_120_204-2" "BS00_240_163-1")

if [[ ! -d ${ROOTDIR}/results/fastqc ]]; then
    mkdir -p ${ROOTDIR}/results/fastqc
    echo -e "sample\tdepth\ttotal_bases" > ${ROOTDIR}/results/fastqc/total_bases_summary.tsv
fi


for SAMPLE in "${SAMPLES_FULL[@]}"; do
    
    SUBSAMPLE_DIR="${ROOTDIR}/data/subsampled/${SAMPLE}"
    NANOPORE_FQ="${SUBSAMPLE_DIR}/${SAMPLE}_${DEPTH_NP}.nanopore.fq"
    ILLUMINA1_FQ="${SUBSAMPLE_DIR}/${SAMPLE}_${DEPTH_IL}.R1.fq"
    ILLUMINA2_FQ="${SUBSAMPLE_DIR}/${SAMPLE}_${DEPTH_IL}.R2.fq"

    mkdir -p ${ROOTDIR}/results/fastqc/${SAMPLE}

    if [[ ! -f ${NANOPORE_FQ} || ! -f ${ILLUMINA1_FQ} || ! -f ${ILLUMINA2_FQ} ]]; then
        echo "One or more input files for sample ${SAMPLE} do not exist. Skipping..."
        continue
    elif [[ -f ${ROOTDIR}/results/fastqc/${SAMPLE}/${SAMPLE}_${DEPTH_NP}.nanopore_fastqc.zip && \
            -f ${ROOTDIR}/results/fastqc/${SAMPLE}/${SAMPLE}_${DEPTH_IL}.R1_fastqc.zip && \
            -f ${ROOTDIR}/results/fastqc/${SAMPLE}/${SAMPLE}_${DEPTH_IL}.R2_fastqc.zip ]]; then
        echo "FastQC output files for sample ${SAMPLE} already exist. Skipping..."
        continue
    fi

sbatch <<- EOF | sed -e "s/Submitted batch job /fastqc_${SAMPLE}/" -e "s/$/,$(date +'%Y-%m-%d %H:%M:%S')/" >> logs/slurm_jobid/fastqc_slurm_jobid.csv
#!/bin/env bash

#SBATCH --job-name=fastqc_${SAMPLE}
#SBATCH --output=logs/fastqc/fastqc_${SAMPLE}.out
#SBATCH --time=1-00
#SBATCH --ntasks=3
#SBATCH --mem=10GB
#SBATCH --partition=general,jrw0107_std,nova_20,nova_28,nova_ff
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=gabriel.silva@auburn.edu

fastqc --threads 3 --memory 2048 --outdir ${ROOTDIR}/results/fastqc/${SAMPLE} ${NANOPORE_FQ} ${ILLUMINA1_FQ} ${ILLUMINA2_FQ}

if [ \$? -ne 0 ]; then
    echo "FastQC failed for sample ${SAMPLE}"
    exit 1
fi

unzip -p ${ROOTDIR}/results/fastqc/${SAMPLE}/${SAMPLE}_${DEPTH_NP}.nanopore_fastqc.zip ${SAMPLE}_${DEPTH_NP}.nanopore_fastqc/fastqc_data.txt | \
    grep "Total Bases" | cut -f2 | (echo -ne "${SAMPLE}\t${DEPTH_NP}_NP\t" && cat) >> ${ROOTDIR}/results/fastqc/total_bases_summary.tsv

unzip -p ${ROOTDIR}/results/fastqc/${SAMPLE}/${SAMPLE}_${DEPTH_IL}.R1_fastqc.zip ${SAMPLE}_${DEPTH_IL}.R1_fastqc/fastqc_data.txt | \
    grep "Total Bases" | cut -f2 | (echo -ne "${SAMPLE}\t${DEPTH_IL}_R1\t" && cat) >> ${ROOTDIR}/results/fastqc/total_bases_summary.tsv

unzip -p ${ROOTDIR}/results/fastqc/${SAMPLE}/${SAMPLE}_${DEPTH_IL}.R2_fastqc.zip ${SAMPLE}_${DEPTH_IL}.R2_fastqc/fastqc_data.txt | \
    grep "Total Bases" | cut -f2 | (echo -ne "${SAMPLE}\t${DEPTH_IL}_R2\t" && cat) >> ${ROOTDIR}/results/fastqc/total_bases_summary.tsv
EOF

done