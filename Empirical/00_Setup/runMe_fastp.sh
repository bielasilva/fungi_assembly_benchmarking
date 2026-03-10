#!/bin/bash

source /home/gas0042/miniforge3/bin/activate Fastp

ROOTDIR="/scratch/gas0042/nanopore_benchmark/real_data"

# Create output directories
if [[ ! -d logs/fastp || ! -d logs/slurm_jobid ]]; then
    mkdir -p logs/fastp
    mkdir -p logs/slurm_jobid
    echo -e "JobName,JobID,SubmissionTime" > logs/slurm_jobid/fastp_slurm_jobid.csv
else
    cp logs/slurm_jobid/fastp_slurm_jobid.csv logs/slurm_jobid/fastp_slurm_jobid.csv.bak
fi

ADAPTERS="/home/gas0042/miniforge3/envs/Trimmomatic/share/trimmomatic-0.39-2/adapters/all.fa"

mkdir -p ${ROOTDIR}/data/filtered/original

# Get the Illumina reads SAMPLEes
fastqs=($(ls ${ROOTDIR}/data/raw/*_R1_001.fastq.gz))

for fastq in "${fastqs[@]}"; do
    SAMPLE=$(basename ${fastq})
    SAMPLE=${SAMPLE%_S*_R1_001.fastq.gz}

    outdir=${ROOTDIR}/data/filtered/original/${SAMPLE}

    echo "Processing sample: ${SAMPLE}"

    mkdir -p ${outdir}

sbatch <<- EOF # | sed -e "s/Submitted batch job /fastp_${SAMPLE},/" -e "s/$/,$(date +'%Y-%m-%d %H:%M:%S')/" >> logs/slurm_jobid/fastp_slurm_jobid.csv
#!/bin/env bash
#SBATCH --job-name=fastp_${SAMPLE}
#SBATCH --output=logs/fastp2/fastp_${SAMPLE}.log
#SBATCH --ntasks=2
#SBATCH --time=1-00
#SBATCH --mem=4G
#SBATCH -p nova_20,nova_28,nova_ff
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=gabriel.silva@auburn.edu

fastp \
    --thread 10 \
    --correction \
    --detect_adapter_for_pe \
    --trim_front1 20 \
    --trim_front2 20 \
    --trim_tail1 2 \
    --trim_tail2 2 \
    --in1 ${ROOTDIR}/data/raw/${SAMPLE}*_R1_001.fastq.gz \
    --in2 ${ROOTDIR}/data/raw/${SAMPLE}*_R2_001.fastq.gz \
    --out1 ${outdir}/${SAMPLE}_trimmed_R1.fastq.gz \
    --out2 ${outdir}/${SAMPLE}_trimmed_R2.fastq.gz \
    --unpaired1 ${outdir}/${SAMPLE}_unpaired.fastq.gz \
    --unpaired2 ${outdir}/${SAMPLE}_unpaired.fastq.gz \
    --json ${outdir}/${SAMPLE}_fastp.json \
    --html ${outdir}/${SAMPLE}_fastp.html

# mv ${outdir}/${SAMPLE}_trimmed_R1.fastq.gz ${outdir}/${SAMPLE}_trimmed_R1.fastq.gz.tmp
# mv ${outdir}/${SAMPLE}_trimmed_R2.fastq.gz ${outdir}/${SAMPLE}_trimmed_R2.fastq.gz.tmp
# mv ${outdir}/${SAMPLE}_unpaired.fastq.gz ${outdir}/${SAMPLE}_unpaired.fastq.gz.tmp

# gunzip -c ${outdir}/${SAMPLE}_trimmed_R1.fastq.gz.tmp | \
#     sed -E "s/\t|\s/|/g" | gzip > ${outdir}/${SAMPLE}_trimmed_R1.fastq.gz

# gunzip -c ${outdir}/${SAMPLE}_trimmed_R1.fastq.gz.tmp | \
#     sed -E 's/\s.*$/\/1/ | gzip > ${outdir}/${SAMPLE}_trimmed_R1.fastq.gz

# rm ${outdir}/${SAMPLE}_trimmed_R1.fastq.gz.tmp

# gunzip -c ${outdir}/${SAMPLE}_trimmed_R2.fastq.gz.tmp | \
#     sed -E "s/\t|\s/|/g" | gzip > ${outdir}/${SAMPLE}_trimmed_R2.fastq.gz

# rm ${outdir}/${SAMPLE}_trimmed_R2.fastq.gz.tmp

# gunzip -c ${outdir}/${SAMPLE}_unpaired.fastq.gz.tmp | \
#     sed -E "s/\t|\s/|/g" | gzip > ${outdir}/${SAMPLE}_unpaired.fastq.gz

# rm ${outdir}/${SAMPLE}_unpaired.fastq.gz.tmp
EOF
done