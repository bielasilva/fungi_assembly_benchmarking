#!/bin/env bash

source /home/gas0042/miniforge3/bin/activate Seqtk

ROOTDIR="/scratch/gas0042/nanopore_benchmark/real_data"

# Create output directories
if [[ ! -d logs/subsampler || ! -d logs/slurm_jobid ]]; then
    mkdir -p logs/subsampler
    mkdir -p logs/slurm_jobid
    echo -e "JobName,JobID,SubmissionTime" > logs/slurm_jobid/subsampler_slurm_jobid.csv
else
    cp logs/slurm_jobid/subsampler_slurm_jobid.csv logs/slurm_jobid/subsampler_slurm_jobid.csv.bak
fi
SAMPLE_LIST="${ROOTDIR}/data/samples_list.txt"

while IFS=, read -r sample depth technology target_coverage coverage_factor
do
        DEPTH=${target_coverage}
        DEPTH_FRACTION=${coverage_factor}
        SAMPLE=${sample}
        
        if [[ $DEPTH == "100x" || $DEPTH == "75x" ]]; then
            # Skipping 100x and 75x 
            continue
        elif [[ $DEPTH_FRACTION > 1 || $DEPTH_FRACTION == "NA" ]]; then
            # Skipping low depth fractions
            continue
        elif grep -q $SAMPLE ${SAMPLE_LIST}; then
            :
        else
            # echo "Sample ${SAMPLE} not found in sample list. Skipping."
            continue
        fi

        OUT_DIR="${ROOTDIR}/data/subsampled/${SAMPLE}"
        
        ILLUMINA1_FQ=`realpath ${ROOTDIR}/data/subsampled/${SAMPLE}/${SAMPLE}_OG.R1.fq`
        ILLUMINA2_FQ=`realpath ${ROOTDIR}/data/subsampled/${SAMPLE}/${SAMPLE}_OG.R2.fq`
        NANOPORE_FQ=`realpath ${ROOTDIR}/data/subsampled/${SAMPLE}/${SAMPLE}_OG.nanopore.fq`
        CHOPPER_FQ=`realpath ${ROOTDIR}/data/subsampled/${SAMPLE}/${SAMPLE}_OG.chopper.fq`

        OUT_ILLUMINA1="${OUT_DIR}/${SAMPLE}_${DEPTH}.R1.fq"
        OUT_ILLUMINA2="${OUT_DIR}/${SAMPLE}_${DEPTH}.R2.fq"
        OUT_NANOPORE="${OUT_DIR}/${SAMPLE}_${DEPTH}.nanopore.fq"
        OUT_CHOPPER="${OUT_DIR}/${SAMPLE}_${DEPTH}.chopper.fq"

        # Create output directory if it doesn't exist
        mkdir -p "$OUT_DIR"

        if [[ $technology == "illumina" ]]; then

                # Check if the assembly file exists and if the job is already running
            if [[ -s $OUT_ILLUMINA1 && -s $OUT_ILLUMINA2 ]]; then
                # echo "Subsampling for ${SAMPLE} at ${DEPTH} already completed. Skipping."
                continue
            elif [[ ! -s $ILLUMINA1_FQ || ! -s $ILLUMINA2_FQ ]]; then
                echo "Not all files for ${SAMPLE} exist. Skipping."
                echo "Missing files: "
                [[ ! -s $ILLUMINA1_FQ ]] && echo "  - ${ILLUMINA1_FQ}"
                [[ ! -s $ILLUMINA2_FQ ]] && echo "  - ${ILLUMINA2_FQ}"
                continue
            elif squeue --me --format "%.100j" | grep -q "Subsample_${SAMPLE}_${DEPTH}"; then
                # echo "Coverage assessment for ${SAMPLE}-${program} is already running. Skipping."
                continue
            else
                echo "Submitting subsampler for ${SAMPLE} at IL${DEPTH}"
            fi
sbatch <<- EOF | sed -e "s/Submitted batch job /Subsample_${SAMPLE}_IL${DEPTH},/" -e "s/$/,$(date +'%Y-%m-%d %H:%M:%S')/" >> logs/slurm_jobid/subsampler_slurm_jobid.csv
#!/bin/env bash

#SBATCH --job-name=Subsample_${SAMPLE}_IL${DEPTH}
#SBATCH --output=logs/subsampler/subsampler_${SAMPLE}_IL${DEPTH}.out
#SBATCH --time=1-00
#SBATCH --ntasks=2
#SBATCH --mem=5GB
#SBATCH --partition=nova_20,nova_28,nova_ff
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=gabriel.silva@auburn.edu

seqtk sample -s 10 ${ILLUMINA1_FQ} ${DEPTH_FRACTION} > ${OUT_ILLUMINA1}
seqtk sample -s 10 ${ILLUMINA2_FQ} ${DEPTH_FRACTION} > ${OUT_ILLUMINA2}

READS_COUNT1=\$(grep -c "^@AV" ${OUT_ILLUMINA1})
READS_COUNT2=\$(grep -c "^@AV" ${OUT_ILLUMINA2})

if [[ ! \${READS_COUNT1} -eq \${READS_COUNT2} ]]; then
    echo "Subsampling failed for ${SAMPLE} at ${DEPTH} with \${READS_COUNT1} reads in R1 and \${READS_COUNT2} reads in R2. Exiting.Exiting."
    exit 1
elif [[ \${READS_COUNT1} -eq 0 || \${READS_COUNT2} -eq 0 ]]; then
    echo "No reads found in one or more files after subsampling for ${SAMPLE} at ${DEPTH}. Exiting."
    exit 1
else
    echo "Subsampling completed for ${SAMPLE} at ${DEPTH} with \${READS_COUNT1} reads in R1 and \${READS_COUNT2} reads in R2."
fi
EOF

        elif [[ $technology == "np" ]]; then
                # Check if the assembly file exists and if the job is already running
            if [[ -s $OUT_NANOPORE ]]; then
                # echo "Subsampling for ${SAMPLE} at ${DEPTH} already completed. Skipping."
                # sed -e "s/400bps_sup@v5.0.0/400bps_sup_v5.0.0/g" -i ${OUT_NANOPORE}
                continue
            elif [[ ! -s $NANOPORE_FQ ]]; then
                echo "Nanopore file for ${SAMPLE} does not exist. Skipping."
                continue
            elif squeue --me --format "%.100j" | grep -q "Subsample_${SAMPLE}_NP${DEPTH}"; then
                # echo "Coverage assessment for ${SAMPLE}-${program} is already running. Skipping."
                continue
            else
                echo "Submitting subsampler for ${SAMPLE} at NP${DEPTH}"
            fi
sbatch <<- EOF | sed -e "s/Submitted batch job /Subsample_${SAMPLE}_NP${DEPTH},/" -e "s/$/,$(date +'%Y-%m-%d %H:%M:%S')/" >> logs/slurm_jobid/subsampler_slurm_jobid.csv
#!/bin/env bash

#SBATCH --job-name=Subsample_${SAMPLE}_NP${DEPTH}
#SBATCH --output=logs/subsampler/subsampler_${SAMPLE}_NP${DEPTH}.out
#SBATCH --time=1-00
#SBATCH --ntasks=2
#SBATCH --mem=5GB
#SBATCH --partition=nova_20,nova_28,nova_ff
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=gabriel.silva@auburn.edu

seqtk sample -s 10 ${NANOPORE_FQ} ${DEPTH_FRACTION} > ${OUT_NANOPORE}

READS_COUNT_NP=\$(grep -c "^@" ${OUT_NANOPORE})

if [[ \${READS_COUNT_NP} -eq 0 ]]; then
    echo "No reads found in nanopore file after subsampling for ${SAMPLE} at ${DEPTH}. Exiting."
    exit 1
else
    echo "Subsampling completed for ${SAMPLE} at ${DEPTH} with \${READS_COUNT_NP} reads in nanopore file."
fi
EOF
        elif [[ $technology == "chopper" ]]; then
                # Check if the assembly file exists and if the job is already running
            if [[ -s $OUT_CHOPPER ]]; then
                # echo "Subsampling for ${SAMPLE} at ${DEPTH} already completed. Skipping."
                # sed -e "s/400bps_sup@v5.0.0/400bps_sup_v5.0.0/g" -i ${OUT_CHOPPER}
                continue
            elif [[ ! -s $CHOPPER_FQ ]]; then
                echo "Nanopore file for ${SAMPLE} does not exist. Skipping."
                continue
            elif squeue --me --format "%.100j" | grep -q "Subsample_${SAMPLE}_CH${DEPTH}"; then
                # echo "Coverage assessment for ${SAMPLE}-${program} is already running. Skipping."
                continue
            else
                echo "Submitting subsampler for ${SAMPLE} at CH${DEPTH}"
            fi
sbatch <<- EOF | sed -e "s/Submitted batch job /Subsample_${SAMPLE}_CH${DEPTH},/" -e "s/$/,$(date +'%Y-%m-%d %H:%M:%S')/" >> logs/slurm_jobid/subsampler_slurm_jobid.csv
#!/bin/env bash

#SBATCH --job-name=Subsample_${SAMPLE}_CH${DEPTH}
#SBATCH --output=logs/subsampler/subsampler_${SAMPLE}_CH${DEPTH}.out
#SBATCH --time=1-00
#SBATCH --ntasks=2
#SBATCH --mem=5GB
#SBATCH --partition=nova_20,nova_28,nova_ff
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=gabriel.silva@auburn.edu

seqtk sample -s 10 ${CHOPPER_FQ} ${DEPTH_FRACTION} > ${OUT_CHOPPER}

READS_COUNT_NP=\$(grep -c "^@" ${OUT_CHOPPER})

if [[ \${READS_COUNT_NP} -eq 0 ]]; then
    echo "No reads found in nanopore file after subsampling for ${SAMPLE} at ${DEPTH}. Exiting."
    exit 1
else
    echo "Subsampling completed for ${SAMPLE} at ${DEPTH} with \${READS_COUNT_NP} reads in nanopore file."
fi
EOF
    fi
done < <(sed 1d /scratch/gas0042/nanopore_benchmark/real_data/results/downsampling_factors.csv)