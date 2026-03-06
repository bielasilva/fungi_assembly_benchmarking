#!/bin/env bash

source /home/gas0042/miniforge3/bin/activate Seqtk

ROOTDIR="/scratch/gas0042/nanopore_benchmark"

# Get the list of SAMPLEs
SAMPLES=($(ls ${ROOTDIR}/genomes/*.fna | sed 's/.*\///; s/_genomic.fna//'))

# Define the reference genomes directory
REF_DIR="${ROOTDIR}/genomes"

# Depth in NNx format
# depthX=("50x" "30x" "15x" "10x" "05x")
depthX=("10x" "15x" "20x" "25x" "30x" "35x" "40x" "50x" "60x" "75x" "100x")

# Depth fraction to supply to samtools assuming 100x is the original depth
depthP=(0.10 0.15 0.20 0.25 0.30 0.35 0.40 0.50 0.60 0.75 1.00)

# Create output directories
if [[ ! -d logs/subsampler || ! -d logs/slurm_jobid ]]; then
    mkdir -p logs/subsampler
    mkdir -p logs/slurm_jobid
    echo -e "JobName,JobID,SubmissionTime" > logs/slurm_jobid/subsampler_slurm_jobid.csv
else
    cp logs/slurm_jobid/subsampler_slurm_jobid.csv logs/slurm_jobid/subsampler_slurm_jobid.csv.bak
fi

for SAMPLE in "${SAMPLES[@]}"; do
    for i in {0..9}; do
        DEPTH=${depthX[$i]}
        DEPTH_FRACTION=${depthP[$i]}

        # NANOPORE_FQ="${ROOTDIR}/simulated_data/original/${SAMPLE}/nanopore/${SAMPLE}_aligned_reads.fastq"
        # ILLUMINA1_FQ="${ROOTDIR}/simulated_data/original/${SAMPLE}/illumina/${SAMPLE}_R1.fastq"
        # ILLUMINA2_FQ="${ROOTDIR}/simulated_data/original/${SAMPLE}/illumina/${SAMPLE}_R2.fastq"
        
        OUT_DIR="${ROOTDIR}/simulated_data/subsampled/${SAMPLE}"

        NANOPORE_FQ="${OUT_DIR}/${SAMPLE}_100x.nanopore.fq"
        ILLUMINA1_FQ="${OUT_DIR}/${SAMPLE}_100x.R1.fq"
        ILLUMINA2_FQ="${OUT_DIR}/${SAMPLE}_100x.R2.fq"


        OUT_NANOPORE="${OUT_DIR}/${SAMPLE}_${DEPTH}.nanopore.fq"
        OUT_ILLUMINA1="${OUT_DIR}/${SAMPLE}_${DEPTH}.R1.fq"
        OUT_ILLUMINA2="${OUT_DIR}/${SAMPLE}_${DEPTH}.R2.fq"

        # Create output directory if it doesn't exist
        mkdir -p "$OUT_DIR"

        # Check if the assembly file exists and if the job is already running
        if [[ -s $OUT_NANOPORE && -s $OUT_ILLUMINA1 && -s $OUT_ILLUMINA2 ]]; then
            # echo "Subsampling for ${SAMPLE} at ${DEPTH} already completed. Skipping."
            continue
        elif [[ $DEPTH == "100x" ]]; then
            echo "Linking original files for ${SAMPLE} at ${DEPTH}"
            cp ${NANOPORE_FQ} ${OUT_NANOPORE}
            cp ${ILLUMINA1_FQ} ${OUT_ILLUMINA1}
            cp ${ILLUMINA2_FQ} ${OUT_ILLUMINA2}
            if cmp --silent ${NANOPORE_FQ} ${OUT_NANOPORE} && cmp --silent ${ILLUMINA1_FQ} ${OUT_ILLUMINA1} && cmp --silent ${ILLUMINA2_FQ} ${OUT_ILLUMINA2}; then
                # echo "Files for ${SAMPLE} at ${DEPTH} are linked successfully."
                continue
            else
                echo "Error linking files for ${SAMPLE} at ${DEPTH}. Exiting."
            fi
            continue
        elif [[ ! -s $NANOPORE_FQ || ! -s $ILLUMINA1_FQ || ! -s $ILLUMINA2_FQ ]]; then
            echo "Not all files for ${SAMPLE} exist. Skipping."
            continue
        elif squeue --me --format "%.100j" | grep -q "Subsample_${SAMPLE}_${DEPTH}"; then
            # echo "Coverage assessment for ${SAMPLE}-${program} is already running. Skipping."
            continue
        else
            echo "Submitting subsampler for ${SAMPLE} at ${DEPTH}"
            # continue
        fi

sbatch <<- EOF | sed -e "s/Submitted batch job /Subsample_${SAMPLE}_${DEPTH},/" -e "s/$/,$(date +'%Y-%m-%d %H:%M:%S')/" >> logs/slurm_jobid/subsampler_slurm_jobid.csv
#!/bin/env bash

#SBATCH --job-name=Subsample_${SAMPLE}_${DEPTH}
#SBATCH --output=logs/subsampler/subsampler_${SAMPLE}_${DEPTH}_%j.out
#SBATCH --time=1-00
#SBATCH --ntasks=2
#SBATCH --mem=5GB
#SBATCH --partition=jrw0107_std,general,nova,nova_ff
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=gabriel.silva@auburn.edu

seqtk sample -s 10 ${NANOPORE_FQ}  ${DEPTH_FRACTION} > ${OUT_NANOPORE}
seqtk sample -s 10 ${ILLUMINA1_FQ} ${DEPTH_FRACTION} > ${OUT_ILLUMINA1}
seqtk sample -s 10 ${ILLUMINA2_FQ} ${DEPTH_FRACTION} > ${OUT_ILLUMINA2}

READS_COUNT_NP=\$(grep -c "^@NC-" ${OUT_NANOPORE})
READS_COUNT1=\$(grep -c "^@NC_" ${OUT_ILLUMINA1})
READS_COUNT2=\$(grep -c "^@NC_" ${OUT_ILLUMINA2})

if [[ ! \${READS_COUNT1} -eq \${READS_COUNT2} ]]; then
    echo "Subsampling failed for ${SAMPLE} at ${DEPTH} with \${READS_COUNT1} reads in R1 and \${READS_COUNT2} reads in R2. Exiting.Exiting."
    exit 1
elif [[ \${READS_COUNT_NP} -eq 0 || \${READS_COUNT1} -eq 0 || \${READS_COUNT2} -eq 0 ]]; then
    echo "No reads found in one or more files after subsampling for ${SAMPLE} at ${DEPTH}. Exiting."
    exit 1
else
    echo "Subsampling completed for ${SAMPLE} at ${DEPTH} with \${READS_COUNT1} reads in R1 and \${READS_COUNT2} reads in R2."
fi

EOF

    done
done