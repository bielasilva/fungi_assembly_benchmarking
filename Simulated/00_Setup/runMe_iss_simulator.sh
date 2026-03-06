#!/bin/bash

source /home/gas0042/miniforge3/bin/activate Insilicoseq

ROOTDIR="/scratch/gas0042/nanopore_benchmark"
GENOMES_DIR="$ROOTDIR/genomes"

for genome in "$GENOMES_DIR"/*.fna; do
    genome_name=$(basename "$genome" _genomic.fna)
    genome_size=$(grep -v ">" "$genome" | wc -c)
    reads_number=$(echo "scale=0; $genome_size / 150 * 100" | bc)

    if [[ -f "$ROOTDIR/simulated_data/${genome_name}/illumina/${genome_name}_R1.fastq" ]]; then
        # echo "Skipping ${genome_name}, already processed."
        continue
    elif squeue --me --format "%.60j" | grep -q "Insilicoseq_${genome_name}"; then
        # echo "Skipping ${genome_name}, job already in queue."
        continue
    else
        echo "Processing genome: ${genome_name}"
    fi

sbatch <<- EOF
#!/bin/bash

#SBATCH --job-name=Insilicoseq_${genome_name}
#SBATCH --output=logs/insilicoseq/iss_simulator_${genome_name}.log
#SBATCH --ntasks=10
#SBATCH --time=1-00
#SBATCH --mem=10GB
#SBATCH --partition=jrw0107_std,general,nova,nova_ff
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=gabriel.silva@auburn.edu

# Create output directory for the genome
mkdir -p "$ROOTDIR/simulated_data/${genome_name}/illumina"

# Run the Insilicoseq simulator
iss generate --cpus 10 \
    --draft $genome \
    --model novaseq \
    --n_reads $reads_number \
    --output $ROOTDIR/simulated_data/${genome_name}/illumina/${genome_name}
EOF
done