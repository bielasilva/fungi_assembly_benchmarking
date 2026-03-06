# Nanopore Benchmark project

## Software versions

- **Nanosim**: 3.2.3
- **InSilicoSeq**: 2.0.1
- **BWA-MEM2**: 2.2.1
- **Minimap2**: 2.30-r1287
- **Samtools**: 1.22
- **Flye**: 2.9.6
- **Hifiasm**: 0.25.0
- **Medaka**: 2.1.0
- **Pilon**: 1.24
- **MaSuRCA**: 4.1.4
- **Canu**: 2.3
- **SPAdes**: 4.2.0
- **Abyss**: 2.3.10
- **Raven**: 1.8.3

## Data description

Downloaded 'assembly_summary.txt' from NCBI Refseq on 2025-07-29 then downloaded the assemblies for all Complete Genomes.

## Assembly strategies

Genome Assembly Approaches:

Depths to be explored = ("100x" "75x" "50x" "35x" "10x")
Depths to be explored = ("100x" "90x" "80x" "70x" "60x" "50x" "40x" "30x" "20x" "10x")

- **Illumina-only Assemblers**: SPAdes, ABySS
- **Nanopore-only Assemblers**: Flye, Canu, Hifiasm, Raven
- **Hybrid Assemblers**: ABySS Hybrid, MaSuRCA, hybridSPAdes
- **Nanopore +  Illumina polishing**: Best Nanopore-only Assembler + Pilon, Polypolish, Racon

## Scripts

All scripts are under the folder `scripts/`.

The scripts to prepare the reads are run as follows:

1. `runMe_download_genomes.sh`: Downloads the assemblies from NCBI Refseq.
2. `runMe_iss_simulator.sh`: Simulates Illumina reads using InSilicoSeq.
3. `runMe_nanosim_simulator.sh`: Simulates Nanopore reads using Nanosim.
4. `runMe_coverage_calc.sh`: Calculates the coverage of the simulated reads.
5. `runMe_subsampler.sh`: Subsamples the simulated reads.

The following scripts are used to run the assembly strategies and can be run in any order:

- `runMe_flye.sh`: Assembles the subsampled reads using Flye.
- `runMe_hifiasm.sh`: Assembles the subsampled reads using Hifiasm.
- `runMe_spades_hybrid.sh`: Assembles the subsampled reads using hybridSPAdes.
- `runMe_spades_short.sh`: Assembles the Illumina reads using SPAdes.
- `runMe_canu.sh`: Assembles the subsampled reads using Canu.
- `runMe_pilon.sh`: Polishes the Flye assemblies using Pilon.
- `runMe_racon.sh`: Polishes the Flye assemblies using Racon.
- `runMe_abyss.sh`: Assembles the Illumina reads using ABySS.
- `runMe_velvet.sh`: Assembles the Illumina reads using Velvet.
- `runMe_masurca.sh`: Assembles the subsampled reads using MaSuRCA.
- `runMe_unicycler.sh`: Assembles the subsampled reads using Unicycler.
- `runMe_medaka.sh`: Polishes the Flye assemblies using Medaka.
- `runMe_autocycler.sh`: Assembles the subsampled reads using Autocycler.

Finally, the script `runMe_quast.sh` is used to evaluate the assemblies using QUAST.

The script `runMe_sacct.sh` is used to track the jobs and resources used during the assembly process. It collects job IDs and their associated information from SLURM.

## Simplified code

### Jobs resources tracking

```bash
sacct --format=JobID,JobName%100,Start,End,Elapsed,ReqMem,MaxRSS,ReqCPUS,NCPUS,CPUTime,CPUTimeRAW,ExitCode -P --delimiter="," -j $JOBID_LIST > ${JOBID_FILE%_jobid.csv}_sacct.csv

/usr/bin/time --format="%e,%U,%S,%P,%M,%t" -o logs/time_flye.txt flye --threads 10 --nano-hq [NANOPORE_FASTQ] --out-dir [OUTPUT_DIR]
```

### InSilicoSeq

```bash
iss generate --cpus 10 \
    --draft $genome \
    --model novaseq \
    --n_reads $reads_number \
    --output $ROOTDIR/simulated_data/${genome_name}/illumina/${genome_name}
```

### Nanosim

```bash
simulator.py genome \
    --model_prefix "$ROOTDIR/models/human_giab_hg002_sub1M_kitv14_dorado_v3.2.1/training" \
    --ref_g "$genome" \
    --output "$ROOTDIR/simulated_data/$genome_name/nanopore/$genome_name" \
    --coverage 100 \
    --fastq \
    --num_threads 10
```

### Verify coverage of simulated reads

```bash

# For Illumina reads
RUNNER="bwa-mem2 index [REFERENCE_FASTA] && bwa-mem2 mem -t 10 [REFERENCE_FASTA_IDX] [ORIGINAL_FASTQ]"
# For Nanopore reads
RUNNER="minimap2 -t 10 -ax map-ont [REFERENCE_FASTA] [ORIGINAL_FASTQ]"

# Run the coverage calculation
$RUNNER | \
samtools sort -@ 5 -o [ALIGNED_BAM] -

# Indexing BAM file
samtools index --threads 10    [ALIGNED_BAM]

# Calculating samtools depth
samtools depth    --threads 10 [ALIGNED_BAM] > [OUTPUT_DEPTH]

# Calculating samtools flagstat
samtools flagstat --threads 10 [ALIGNED_BAM] > [OUTPUT_FLAGSTAT]

# Calculating samtools stats
samtools stats    --threads 10 [ALIGNED_BAM] > [OUTPUT_STATS]

# Calculating samtools idxstats
samtools idxstats --threads 10 [ALIGNED_BAM] > [OUTPUT_IDXSTATS]

# Calculating samtools coverage
samtools coverage              [ALIGNED_BAM] > [OUTPUT_COVERAGE]
```

### Subsampler

```bash
seqtk sample -s 10 [ORIGINAL_FASTQ] [DEPTH_FRACTION] > [SUBSAMPLED_FASTQ]
```

### Flye

```bash
flye \
    --threads 10 \
    --nano-hq [NANOPORE_FASTQ] \
    --out-dir [OUTPUT_DIR] 
```

### Hifiasm

```bash
hifiasm \
    -t 10 \
    --ont \
    -o [OUTPUT] \
    [NANOPORE_FASTQ]
```

#### Notes

Hifiasm_GCF_000687475.1_ASM68747v2_NP10xHifiasm_GCF_000687475.1_ASM68747v2_NP10x could not be run. Gives the error:

```text
/cm/local/apps/slurm/var/spool/job3160348/slurm_script: line 13:  4712 Floating point exception hifiasm -t 10 --ont -o /scratch/gas0042/nanopore_benchmark/results/assemblies/GCF_000687475.1_ASM68747v2/hifiasm/10x/hifiasm_GCF_000687475.1_ASM68747v2_NP10x /scratch/gas0042/nanopore_benchmark/simulated_data/subsampled/GCF_000687475.1_ASM68747v2/GCF_000687475.1_ASM68747v2_10x.nanopore.fq
```

### Pilon

```bash
# Index the FLYE assembly with BWA
bwa-mem2 index \
    -p [IDX_PREFIX] \
    [FLYE_ASSEMBLY]

# Align the Illumina reads to the FLYE assembly and sort the BAM file
bwa-mem2 mem \
    -t 10 \
    [IDX_PREFIX] \
    [ILLUMINA_FASTQ_1] \
    [ILLUMINA_FASTQ_2] | \
    samtools view -b - | samtools sort > [ALIGNED_BAM]

# Index the BAM file
samtools index [ALIGNED_BAM]

# Run Pilon to polish the FLYE assembly using the Illumina reads
pilon \
    --changes \
    --genome [FLYE_ASSEMBLY] \
    --frags [ALIGNED_BAM] \
    --output [OUTPUT_PREFIX]
```

### Racon

### Notes

Racon needs the Illumina reads to be in a single file, so the paired-end reads need to be combined first. The following command combines the paired-end reads by changing the suffix `/1` and `/2` to `_1` and `_2`, respectively. Therefore while the reads are still technically paired-end, they are treated as single-end reads by Racon.

```bash
cat [ILLUMINA_FASTQ_1] [ILLUMINA_FASTQ_2] | awk '{ sub(/\/[12]$/, "_" substr(\$0, length, 1)); print }' > [ILLUMINA_FASTQ_COMBINED]

# Index the FLYE assembly with BWA
bwa-mem2 index \
    -p [IDX_PREFIX] \
    [FLYE_ASSEMBLY]

# Align the Illumina reads to the FLYE assembly and sort the BAM file
bwa-mem2 mem \
    -t 10 \
    ${BWA_IDX} \
    ${SAMPLE}_${DEPTH_IL}.combined.fq > ${BWA_SAM}

# Run Racon to polish the FLYE assembly using the Illumina reads
racon \
    -t 10 \
    ${SAMPLE}_${DEPTH_IL}.combined.fq \
    ${BWA_SAM} \
    ${FLYE_ASM} \
    > ${RACON_OUT}.fasta
```

### Canu

Canu needs genome size estimation, so I used the reported genome size from NCBI Refseq. In real scenarios, you might want to use a tool like `jellyfish` or `kmergenie` to estimate the genome size from the reads.

Canu's minimum coverage is set to 10x, it was lowered to 5x for this benchmark since some of the 10x would be reported as "not enough coverage".

```bash
canu \
    -p [OUTPUT_PREFIX] \
    -d [OUTPUT_DIR] \
    genomeSize=[GENOME_SIZE] \
    -nanopore [NANOPORE_FASTQ] \
    useGrid=false \
    minInputCoverage=5 \
    maxThreads=10 \
    maxMemory=120
```

### AbySS

Running Abyss installed from Bioconda was raising errors with the MPI library, then it would get stuck in `Finding adjacent k-mer...`. The error was:

```text
WARNING: There was an error initializing an OpenFabrics device.

  Local host:   node026
  Local device: mlx5_0
```

So I installed most dependencies with Mamba and then compiled Abyss from source.

#### Instalation

```bash
#Mamba env preparation
mamba create -n Abyss_man
mamba activate Abyss_man
mamba install conda-forge::boost conda-forge::openmpi bioconda::google-sparsehash bioconda::btllib bioconda::samtools conda-forge::compilers conda-forge::pigz conda-forge::zsh conda-forge::autoconf conda-forge::automake bioconda::bwa

#Abyss compilation
git clone git@github.com:bcgsc/abyss.git #Version 2.3.10
cd abyss
./autogen.sh
mkdir build
cd build
../configure --prefix=/mmfs1/home/gas0042/miniforge3/envs/Abyss_man --enable-maxk=96 --disable-werror
make
make install
```

Running Abyss-short:

### Racon

### Shasta

#### Notes

Shasta's logs states 'This run used options "--memoryBacking 4K --memoryMode anonymous". This could result in longer run time.
For faster assembly, use "--memoryBacking 2M --memoryMode filesystem" (root privilege via sudo required).
Therefore the results of this run should not be used for the purpose of benchmarking assembly time.
However the memory options don't affect assembly results in any way.'

Also, the cutoff for read size is set to 10kb by default.
