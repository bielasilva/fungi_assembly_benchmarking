#!/bin/env bash

source /scratch/gas0042/nanopore_benchmark/real_data/scripts/sourceMe_configs.sh

load_env Coverage

# Create output directories
if [[ ! -d logs/coverage || ! -d logs/slurm_jobid ]]; then
    mkdir -p logs/coverage
    mkdir -p logs/slurm_jobid
    echo -e "JobName,JobID,SubmissionTime" > logs/slurm_jobid/coverage_slurm_jobid.csv
else
    cp logs/slurm_jobid/coverage_slurm_jobid.csv logs/slurm_jobid/coverage_slurm_jobid.csv.bak
fi

for SAMPLE in "${SAMPLES_FULL[@]}"; do
    sample_number=$(( sample_number + 1 ))
    ASSEMBLY_DIR="${ROOTDIR}/results/assemblies_results/${SAMPLE}"
    for PROGRAM in "${PROGRAMS[@]}"; do
        case $PROGRAM in
            "flye-OVL1000" | "flye-OVL1500" | "flye-OVL2000" | "flye-OVL2500" | "flye-OVL3000")
                PREFIX="${PROGRAM}_${SAMPLE}_NPOG"
                ASSEMBLY_FASTA="${ASSEMBLY_DIR}/${PREFIX}.fasta"
                OUTDIR="${ROOTDIR}/results/original_coverage/${SAMPLE}/${PREFIX}"
                BASE_OUT="${OUTDIR}/${PREFIX}"
                mkdir -p ${OUTDIR}

                FASTQ="${ROOTDIR}/data/subsampled/${SAMPLE}/${SAMPLE}_OG.chopper.fq"
                RUNNER="minimap2 -t 10 -ax lr:hq ${ASSEMBLY_FASTA} ${FASTQ}"
            ;;
            "flye" | "hifiasm" | "canu" | "raven" | "shasta" | "nextdenovo" | "miniasm")
                PREFIX="${PROGRAM}_${SAMPLE}_NPOG"
                ASSEMBLY_FASTA="${ASSEMBLY_DIR}/${PREFIX}.fasta"
                OUTDIR="${ROOTDIR}/results/original_coverage/${SAMPLE}/${PREFIX}"
                BASE_OUT="${OUTDIR}/${PREFIX}"
                mkdir -p ${OUTDIR}

                FASTQ="${ROOTDIR}/data/subsampled/${SAMPLE}/${SAMPLE}_OG.chopper.fq"
                RUNNER="minimap2 -t 10 -ax lr:hq ${ASSEMBLY_FASTA} ${FASTQ}"
            ;;
            "spades_short" | "abyss_short")
                PREFIX="${PROGRAM}_${SAMPLE}_ILOG"
                ASSEMBLY_FASTA="${ASSEMBLY_DIR}/${PREFIX}.fasta"
                FASTQ="${ROOTDIR}/data/subsampled/${SAMPLE}/${SAMPLE}_OG.R*.fq"
                
                OUTDIR="${ROOTDIR}/results/original_coverage/${SAMPLE}/${PREFIX}"
                BASE_OUT="${OUTDIR}/${PREFIX}"
                mkdir -p ${OUTDIR}/bwa
                RUNNER="minimap2 -t 10 -ax sr ${ASSEMBLY_FASTA} ${FASTQ}"
                # ln -sf ${ASSEMBLY_FASTA} ${BASE_OUT}.fasta
                # FASTA_IDX="${BASE_OUT}.fasta"
                # RUNNER="bwa-mem2 index ${FASTA_IDX} && bwa-mem2 mem -t 10 ${FASTA_IDX} ${FASTQ}"
            ;;
        esac

        # Check if the assembly file exists and if the job is already running
        if [[ -f "${BASE_OUT}.bam" ]]; then
            # echo "Quast for ${SAMPLE}-${program} already completed. Skipping."
            continue
        elif [[ ! -f "${ASSEMBLY_FASTA}" ]]; then
            # echo "Assembly file ${ASSEMBLY_FASTA} does not exist. Skipping coverage calculation for ${PREFIX}."
            continue
        elif squeue --me --format "%.100j" | grep -q "Coverage_${PREFIX}"; then
            # echo "Quast for ${SAMPLE}-${program} is already running. Skipping."
            continue
        else
            echo "Submitting coverage calculator for ${PREFIX} "
        fi

        

sbatch <<- EOF | sed "s/Submitted batch job /Coverage_${PREFIX}\t/g" >> logs/slurm_jobid/coverage_slurm_jobid.csv
#!/bin/env bash

#SBATCH --job-name=Coverage_${PREFIX}
#SBATCH --output=logs/coverage/coverage_${PREFIX}.out
#SBATCH --time=1-00
#SBATCH --ntasks=15
#SBATCH --mem=100GB
#SBATCH --partition=jrw0107_std,general,nova_20,nova_28,nova_ff
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=gabriel.silva@auburn.edu

echo "Mapping and sorting BAM file"
$RUNNER | \
samtools sort -@ 5 -o ${BASE_OUT}.bam -

echo "Indexing BAM file"
samtools index --threads 10    ${BASE_OUT}.bam

echo "Calculating samtools depth"
samtools depth    --threads 10 ${BASE_OUT}.bam > ${BASE_OUT}.depth.txt

echo "Calculating samtools flagstat"
samtools flagstat --threads 10 ${BASE_OUT}.bam > ${BASE_OUT}.flagstat.txt

echo "Calculating samtools stats"
samtools stats    --threads 10 ${BASE_OUT}.bam > ${BASE_OUT}.stats.txt

echo "Calculating samtools idxstats"
samtools idxstats --threads 10 ${BASE_OUT}.bam > ${BASE_OUT}.idxstats.txt

echo "Calculating samtools coverage"
samtools coverage              ${BASE_OUT}.bam > ${BASE_OUT}.coverage.txt

EOF

    done
done


