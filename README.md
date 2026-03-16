# Fungi Assembly Benchmarking

This repository stores the scripts used for two benchmarking tracks:

- **Empirical**: assembly benchmarking using real sequencing data.
- **Simulated**: assembly benchmarking using simulated short_reads/Nanopore reads.

## Software versions

- **Nanosim**: 3.2.3
- **InSilicoSeq**: 2.0.1
- **BWA-MEM2**: 2.2.1
- **Minimap2**: 2.30-r1287
- **Samtools**: 1.22
- **Flye**: 2.9.6
- **Hifiasm**: 0.25.0
- **Pilon**: 1.24
- **MaSuRCA**: 4.1.4
- **Canu**: 2.3
- **SPAdes**: 4.2.0
- **ABySS**: 2.3.10
- **Raven**: 1.8.3

## Data source summary

- `assembly_summary.txt` from NCBI RefSeq fungi was downloaded on **2025-07-29**.
- Complete fungal reference genomes were then downloaded and used for simulation/benchmarking.

## Pipeline layout

- `Empirical/00_Setup`: real-read preprocessing, QC, and subsampling.
- `Empirical/01_Assembly`: empirical assemblies and polishing runs.
- `Empirical/02_Processing`: empirical QUAST/BUSCO and aggregation.
- `Empirical/Misc`: helper scripts to consolidate outputs and track completion.
- `Simulated/00_Setup`: download references and simulate/subsample reads.
- `Simulated/01_Assembly`: simulated-read assembly and polishing runs.
- `Simulated/02_Processing`: simulated QUAST/BUSCO + reference annotation helpers.

---

## File index and purpose

### Empirical

#### Empirical/00_Setup

- `runMe_chopper.sh`: filters raw Nanopore reads with Chopper (minimum length) to produce cleaned long-read inputs.
- `runMe_coverage_calc.sh`: maps reads back to references/assemblies and computes coverage/depth/statistics (samtools/minimap2/BWA workflows).
- `runMe_fastp.sh`: trims and filters empirical short_reads paired-end reads using fastp.
- `runMe_fastqc.sh`: runs FastQC on subsampled read sets and summarizes total bases per sample/depth.
- `runMe_get_info.R`: parses fastp/KMC-GenomeScope outputs and generates sample-level depth and downsampling-factor summary tables.
- `runMe_kmc-genomescope_chopper.sh`: runs KMC + GenomeScope on Chopper-filtered Nanopore reads for genome-size/error modeling.
- `runMe_kmc-genomescope_short.sh`: runs KMC + GenomeScope on short_reads reads.
- `runMe_kmc-genomescope_np.sh`: runs KMC + GenomeScope on Nanopore reads.
- `runMe_nanoplot.sh`: runs NanoPlot to summarize Nanopore read quality/length distributions.
- `runMe_nanopore_check.sh`: checks for Nanopore FASTQ availability per sample and extracts/normalizes selected FASTQ inputs.
- `runMe_subsampler.sh`: subsamples short_reads, Nanopore, and Chopper FASTQs to target depths using seqtk.
- `runMe_uncompresser.sh`: decompresses trimmed short_reads FASTQ files to benchmark-ready uncompressed `.fq` files in parallel.
- `sbatchMe_seqkit_stats.sh`: compute seqkit stats for all subsampled FASTQ files.
- `sourceMe_configs.sh`: shared configuration/utilities (sample lists, depth arrays, program lists, SLURM helpers, progress bar, env loader).

#### Empirical/01_Assembly

- `runMe_abyss_short.sh`: submits short_reads-only ABySS assemblies across samples/depths.
- `runMe_flye.sh`: submits Flye long-read assemblies (multiple overlap parameter values) from Chopper-filtered reads.
- `runMe_polypolish.sh`: polishes Flye assemblies using short_reads alignments and Polypolish.
- `runMe_spades_hybrid.sh`: runs hybrid SPAdes with paired short_reads + Nanopore inputs.
- `runMe_spades_short.sh`: runs short_reads-only SPAdes assemblies.

#### Empirical/02_Processing

- `quast_busco_merger.R`: incrementally merges per-assembly QUAST and BUSCO outputs into unified summary tables.
- `runMe_busco.sh`: submits BUSCO completeness analysis jobs for assembled genomes/reference inputs.
- `runMe_quast_busco_merger.sh`: executes the R merger script for QUAST/BUSCO consolidated outputs.
- `runMe_quast.sh`: submits QUAST jobs for all configured assemblies and depth combinations.
- `runMe_sacct.sh`: extracts SLURM accounting metrics for submitted jobs and produces combined sacct tables.
- `stats_figures.qmd`: generates statistical analyses and figures for the empirical assembly benchmarking results.

#### Empirical/Misc

- `runMe_clean_assemblies.sh`: copies finalized assembly/QC outputs into centralized results folders and optionally compresses/removes intermediates.
- `runMe_count_complete.sh`: counts completed vs pending assemblies and associated QUAST/BUSCO outputs by program.

### Simulated

#### Simulated/00_Setup

- `runMe_download_genomes.sh`: downloads RefSeq fungal assembly summary and reference genome FASTA files.
- `runMe_iss_simulator.sh`: simulates short_reads reads with InSilicoSeq.
- `runMe_nanosim_simulator.sh`: simulates Nanopore reads with NanoSim.
- `runMe_subsampler.sh`: subsamples simulated short_reads and Nanopore FASTQs to target coverage fractions.

#### Simulated/01_Assembly

- `runMe_abyss_short.sh`: runs short_reads-only ABySS (`abyss-short`) assemblies.
- `runMe_abyss-hybrid.sh`: runs ABySS hybrid assemblies with short_reads + Nanopore inputs.
- `runMe_canu.sh`: runs Canu long-read assemblies.
- `runMe_flye.sh`: runs Flye long-read assemblies.
- `runMe_hifiasm.sh`: runs Hifiasm in ONT mode.
- `runMe_masurca.sh`: runs MaSuRCA hybrid assemblies.
- `runMe_miniasm.sh`: runs Minimap2 + Miniasm long-read assembly workflow.
- `runMe_nextdenovo.sh`: runs NextDenovo assemblies using generated run configs.
- `runMe_pilon.sh`: polishes Flye assemblies with Pilon using short_reads alignments.
- `runMe_polypolish.sh`: polishes Flye assemblies with Polypolish using paired short_reads alignments.
- `runMe_racon_bwa.sh`: polishes Flye assemblies with Racon using BWA-based read alignments.
- `runMe_racon_minimap.sh`: polishes Flye assemblies with Racon using Minimap2-based alignments.
- `runMe_raven.sh`: runs Raven long-read assemblies.
- `runMe_spades_hybrid.sh`: runs hybrid SPAdes assemblies.
- `runMe_spades_short.sh`: runs short_reads-only SPAdes assemblies.

#### Simulated/02_Processing

- `runMe_busco.sh`: submits BUSCO analyses for simulated assembly outputs.
- `runMe_edta.sh`: runs EDTA transposable-element annotation on selected reference genomes.
- `runMe_genemark.sh`: runs GeneMark-ES/Fungal gene prediction on reference genomes.
- `runMe_quast.sh`: submits QUAST evaluation for simulated assemblies/reference genomes.
- `runMe_sacct.sh`: compiles SLURM resource/accounting summaries for simulated workflow jobs.
- `stats_figures.qmd`: generates statistical figures and analyses for the simulated data.

#### Simulated/02_Processing/Summarizer

- `quast_busco_merger.R`: merges simulated QUAST and BUSCO outputs into unified reports.
- `runMe_quast_busco_merger.sh`: launches the simulated QUAST/BUSCO merger R script.

---

## Notes

- Most scripts are SLURM submission wrappers designed for HPC execution.
- Paths in scripts are currently hard-coded to the original project environment under `/scratch/gas0042/...`.
- Shared helper functions and core arrays are defined in each pipeline’s `sourceMe_configs.sh` file (empirical file is included in this repository).
