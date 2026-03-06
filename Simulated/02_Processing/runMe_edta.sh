#!/bin/bash

source /home/gas0042/miniforge3/bin/activate Edta

ROOTDIR="/scratch/gas0042/nanopore_benchmark"
GENOMES_DIR="$ROOTDIR/genomes"

# Get the list of SAMPLEs
SAMPLES=($(ls ${ROOTDIR}/genomes/*.fna | sed 's/.*\///; s/_genomic.fna//'))

SAMPLES=(
    "GCF_000002515.2_ASM251v1"
    "GCF_000026945.1_ASM2694v1"
)

for SAMPLE in "${SAMPLES[@]}"; do
    genome="${GENOMES_DIR}/${SAMPLE}_genomic.fna"
    work_dir="$ROOTDIR/results/ref_edta/${SAMPLE}"

    if [[ -f "$work_dir/edta.gtf" ]]; then
        # echo "Skipping ${SAMPLE}, already processed."
        continue
    elif squeue --me --format "%.60j" | grep -q "edta_${SAMPLE}"; then
        # echo "Skipping ${SAMPLE}, job already in queue."
        continue
    else
        echo "Processing genome: ${SAMPLE}"
    fi

sbatch <<- EOF
#!/bin/bash

#SBATCH --job-name=edta_${SAMPLE}
#SBATCH --output=logs/edta/edta_${SAMPLE}.log
#SBATCH --ntasks=10
#SBATCH --time=1-00
#SBATCH --mem=10GB
#SBATCH --partition=nova_20,nova_28,nova_ff
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=gabriel.silva@auburn.edu

# Create output directory for the genome
mkdir -p $work_dir

cd $work_dir

# Run EDTA

EDTA.pl --genome $genome \
    --species others \
    --step all \
    --sensitive 1 \
    --anno 1 \
    --evaluate 1 \
    --threads 10

# ln -s $work_dir/.mod.EDTA.TEanno.gff3 $work_dir/edta_${SAMPLE}.gtf

EOF
done