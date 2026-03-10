#!/bin/env bash

source /scratch/gas0042/nanopore_benchmark/real_data/scripts/sourceMe_configs.sh

load_env Nanoplot

ROOTDIR="/scratch/gas0042/nanopore_benchmark/real_data"

# Create output directories
if [[ ! -d logs/nanoplot || ! -d logs/slurm_jobid ]]; then
    mkdir -p logs/nanoplot
    mkdir -p logs/slurm_jobid
    echo -e "JobName,JobID,SubmissionTime" > logs/slurm_jobid/nanoplot_slurm_jobid.csv
else
    cp logs/slurm_jobid/nanoplot_slurm_jobid.csv logs/slurm_jobid/nanoplot_slurm_jobid.csv.bak
fi

DEPTH_NP="OG"
DEPTH_IL="OG"

# if [[ ! -d ${ROOTDIR}/results/nanoplot ]]; then
#     mkdir -p ${ROOTDIR}/results/nanoplot
#     echo -e "sample\tdepth\ttotal_bases" > ${ROOTDIR}/results/nanoplot/total_bases_summary.tsv
# fi


for SAMPLE in "${SAMPLES_FULL[@]}"; do
    
    SUBSAMPLE_DIR="${ROOTDIR}/data/subsampled/${SAMPLE}"
    NANOPORE_FQ="${SUBSAMPLE_DIR}/${SAMPLE}_${DEPTH_NP}.nanopore.fq"

    mkdir -p ${ROOTDIR}/results/nanoplot/${SAMPLE}

    if [[ ! -s ${NANOPORE_FQ} ]]; then
        echo "One or more input files for sample ${SAMPLE} do not exist. Skipping..."
        continue
    elif [[ -s ${ROOTDIR}/results/nanoplot/${SAMPLE}/${SAMPLE}_${DEPTH_NP}.NanoStats.txt ]]; then
        echo "nanoplot output files for sample ${SAMPLE} already exist. Skipping..."
        continue
    fi

sbatch <<- EOF | sed -e "s/Submitted batch job /nanoplot_${SAMPLE}/" -e "s/$/,$(date +'%Y-%m-%d %H:%M:%S')/" >> logs/slurm_jobid/nanoplot_slurm_jobid.csv
#!/bin/env bash

#SBATCH --job-name=nanoplot_${SAMPLE}
#SBATCH --output=logs/nanoplot/nanoplot_${SAMPLE}.out
#SBATCH --time=1-00
#SBATCH --ntasks=3
#SBATCH --mem=10GB
#SBATCH --partition=general,jrw0107_std,nova_20,nova_28,nova_ff
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=gabriel.silva@auburn.edu

NanoPlot --threads 3 --tsv_stats --outdir ${ROOTDIR}/results/nanoplot/${SAMPLE} --fastq ${NANOPORE_FQ} --prefix ${SAMPLE}_${DEPTH_NP}.

if [ \$? -ne 0 ]; then
    echo "nanoplot failed for sample ${SAMPLE}"
    exit 1
fi

EOF

done