#!/bin/bash

source /home/gas0042/miniforge3/bin/activate Funannotate

ROOTDIR="/scratch/gas0042/nanopore_benchmark"
GENOMES_DIR="$ROOTDIR/genomes"

gm_path="/home/gas0042/tools/gmes/gmes_linux_64_4/gmes_petap.pl"

# Get the list of SAMPLEs
SAMPLES=($(ls ${ROOTDIR}/genomes/*.fna | sed 's/.*\///; s/_genomic.fna//'))

for SAMPLE in "${SAMPLES[@]}"; do
    genome="${GENOMES_DIR}/${SAMPLE}_genomic.fna"
    work_dir="$ROOTDIR/results/ref_genemark/${SAMPLE}"

    if [[ -f "$work_dir/genemark.gtf" ]]; then
        # echo "Skipping ${SAMPLE}, already processed."
        continue
    elif squeue --me --format "%.60j" | grep -q "genemark_${SAMPLE}"; then
        # echo "Skipping ${SAMPLE}, job already in queue."
        continue
    else
        echo "Processing genome: ${SAMPLE}"
        continue
    fi

sbatch <<- EOF
#!/bin/bash

#SBATCH --job-name=genemark_${SAMPLE}
#SBATCH --output=logs/genemark/genemark_${SAMPLE}.log
#SBATCH --ntasks=10
#SBATCH --time=1-00
#SBATCH --mem=10GB
#SBATCH --partition=nova_20,nova_28,nova_ff
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=gabriel.silva@auburn.edu

# Create output directory for the genome
mkdir -p $work_dir

# Run GeneMark ES Fungal

$gm_path --ES --sequence $genome --fungus --work_dir $work_dir --cores 10

ln -s $work_dir/genemark.gtf $work_dir/genemark_${SAMPLE}.gtf

EOF
done