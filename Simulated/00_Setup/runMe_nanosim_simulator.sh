#!/bin/bash

source /home/gas0042/miniforge3/bin/activate Nanosim

ROOTDIR="/scratch/gas0042/nanopore_benchmark"
GENOMES_DIR="$ROOTDIR/genomes"

for genome in "$GENOMES_DIR"/*.fna; do
    genome_name=$(basename "$genome" _genomic.fna)

    if [[ -f "$ROOTDIR/simulated_data/$genome_name/nanopore/${genome_name}_aligned_reads.fastq" ]]; then
        # echo "Skipping $genome_name, already processed."
        continue
    elif squeue --me --format "%.60j" | grep -q "nanosim_$genome_name"; then
        # echo "Skipping $genome_name, job already in queue."
        continue
    else
        echo "Processing genome: $genome_name"
    fi


sbatch <<- EOF
#!/bin/bash

#SBATCH --job-name=nanosim_$genome_name
#SBATCH --output=logs/nanosim/simulator_$genome_name.log
#SBATCH --ntasks=10
#SBATCH --time=1-00
#SBATCH --mem=10GB
#SBATCH --partition=jrw0107_std,general,nova,nova_ff
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=gabriel.silva@auburn.edu

# Create output directory for the genome
mkdir -p "$ROOTDIR/simulated_data/$genome_name/nanopore"

# Run the nanosim simulator
simulator.py genome \
    --model_prefix "$ROOTDIR/models/human_giab_hg002_sub1M_kitv14_dorado_v3.2.1/training" \
    --ref_g "$genome" \
    --output "$ROOTDIR/simulated_data/$genome_name/nanopore/$genome_name" \
    --coverage 100 \
    --fastq \
    --num_threads 10
EOF
done