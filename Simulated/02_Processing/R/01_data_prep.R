# ============================================================================
# 01_data_prep.R
# ----------------------------------------------------------------------------
# Simulated data extraction, cleaning, and computation of the per-isolate
# genome-size covariate used in the mixed models.
#
# This script is sourced by:
#   - Simulated/02_Processing/stats_figures.qmd
#   - Simulated/R/02_fit_models.R
#
# Inputs (all under Simulated/02_Processing/files/):
#   - all_sacct.csv              SLURM accounting records
#   - quast_busco_results.tsv    Combined QUAST + BUSCO metrics
#   - ref_genome_metrics.csv     Metrics for the original genomes
#
# Outputs:
#   - data_simulated (in-memory)      Cleaned model-ready tibble
#   - reference_metrics (in-memory)   Per-isolate genome-size estimates
# ============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(janitor)
})

# ----------------------------------------------------------------------------
# Helper functions
# ----------------------------------------------------------------------------

`%notin%` <- Negate(`%in%`)

CI95 <- function(x) {
  x <- x[!is.na(x)]
  n <- length(x)
  alpha <- 0.05
  df    <- n - 1
  t_score <- qt(p = alpha / 2, df = df, lower.tail = FALSE)
  t_score * (sd(x) / sqrt(n))
}

# Parse SLURM Elapsed strings to seconds (handles both "DD-HH:MM:SS" and "HH:MM:SS")
parse_elapsed <- function(x) {
  x <- as.character(x) %>% stringr::str_trim()
  ifelse(
    is.na(x) | x == "NA",
    NA_real_,
    ifelse(
      stringr::str_detect(x, "^[0-9]+-"),
      {
        m <- stringr::str_match(x, "^(\\d+)-(\\d{1,2}):(\\d{2}):(\\d{2})")
        d  <- as.numeric(m[, 2])
        h  <- as.numeric(m[, 3])
        mi <- as.numeric(m[, 4])
        s  <- as.numeric(m[, 5])
        d * 86400 + h * 3600 + mi * 60 + s
      },
      {
        m <- stringr::str_match(x, "^(\\d{1,2}):(\\d{2}):(\\d{2})")
        h  <- as.numeric(m[, 2])
        mi <- as.numeric(m[, 3])
        s  <- as.numeric(m[, 4])
        h * 3600 + mi * 60 + s
      }
    )
  )
}

# Combine two assemblies (e.g., Flye + Polypolish) into one polished row.
# Joins on (sample, depth_np); takes most metrics from prog1 (Polypolish)
# and combines computational metrics (RAM via pmax, time via sum).
combine_programs <- function(df, prog1, prog2) {
  combo_name <- paste0(prog1, "_", prog2)
  
  p1 <- df %>% filter(!is.na(n50), program == prog1)
  p2 <- df %>% filter(!is.na(n50), program == prog2)
  
  p2 %>%
    right_join(p1, by = c("sample", "depth_np")) %>%
    transmute(
      strategy               = "polished",
      program                = combo_name,
      sample,
      depth_np,
      depth_il               = depth_il.x,
      n50                    = n50.x,
      busco_complete         = busco_complete.x,
      busco_single_copy      = busco_single_copy.x,
      busco_multi_copy       = busco_multi_copy.x,
      genome_fraction_percent= genome_fraction_percent.x,
      total_length           = total_length.x,
      mism_100kbp            = mism_100kbp.x,
      indel_100kbp           = indel_100kbp.x,
      length_err_pct         = length_err_pct.x,
      gc_err_abs             = gc_err_abs.x,
      busco_complete_diff    = busco_complete_diff.x,
      busco_single_copy_diff = busco_single_copy_diff.x,
      busco_multi_copy_diff  = busco_multi_copy_diff.x,
      maxrss_gb              = pmax(maxrss_gb.x, maxrss_gb.y),
      ram_per_mb             = pmax(ram_per_mb.x, ram_per_mb.y),
      time_per_mb            = time_per_mb.x + time_per_mb.y,
      elapsed_sec            = elapsed_sec.x + elapsed_sec.y
    )
}

colors_programs <- c("flye"        ="#59A14F",
                     "hifiasm"     ="#1F83B4",
                     "canu"        ="#C7519C",
                     "raven"       ="#A52A2A",
                     "nextdenovo"  ="#8C564B",
                     "miniasm"     ="#7B68EE",
                     "spades_short"="#EDC948",
                     "abyss_short" ="#FF7F0E",
                     
                     "abyss_hybrid" ="#FFAA0E",
                     "masurca"      ="#8B4513",
                     "spades_hybrid"="#CD1076",
                     
                     "pilon"          ="#B22222",
                     "polypolish"     ="#00688B",
                     "racon"          ="#008B00",
                     "racon-mp2"      ="#FF7F0E",
                     "flye_pilon"     ="#B22222",
                     "flye_polypolish"="#00688B",
                     "flye_racon"     ="#008B00",
                     "flye_racon-mp2" ="#FF7F0E")

# ----------------------------------------------------------------------------
# Main extraction function
# ----------------------------------------------------------------------------

extract_data_simulated <- function(sacct_file, metrics_file, genome_sizes) {
  
  # Define strategies
  long_reads_strategy  <- c("flye", "hifiasm", "canu", "raven", "shasta", "nextdenovo", "miniasm")
  short_reads_strategy <- c("spades_short", "abyss", "abyss_short")
  hybrid_strategy      <- c("spades_hybrid", "masurca", "abyss_hybrid")
  polished_strategy    <- c("pilon", "polypolish", "racon", "racon-mp2")
  
  rename_pattern <- c("Hifiasm"  = "hifiasm",
                      "abyss-" = "abyss_",
                      "_FLYE"  ="_NP",
                      "_10x"   ="_NP10x",
                      "_15x"   ="_NP15x",
                      "_20x"   ="_NP20x",
                      "_25x"   ="_NP25x",
                      "_30x"   ="_NP30x",
                      "_35x"   ="_NP35x",
                      "_40x"   ="_NP40x",
                      "_50x"   ="_NP50x",
                      "_60x"   ="_NP60x",
                      "_75x"   ="_NP75x",
                      "_100x"  ="_NP100x")

  raw_sacct <- read_csv(sacct_file, col_types = cols()) %>% 
    mutate(JobName = str_replace_all(JobName, pattern = rename_pattern))
  
  raw_metrics <- read_tsv(metrics_file, col_types = cols()) %>% 
    rename(JobName = Assembly) %>%
    mutate(JobName = str_replace_all(JobName, pattern = rename_pattern))
  
  genome_sizes_df <- read_csv(genome_sizes, col_types=cols()) %>% 
    dplyr::select(Sample, RefSizeMbp=size_mbp, RefSize=size_bp, RefContigs=contigs)
  
  raw_metrics_reference <- raw_metrics %>%
    filter(str_detect(JobName, "reference")) %>% 
    select(JobName,
           busco_complete_ref     =busco_complete,
           busco_single_copy_ref  =busco_single_copy,
           busco_multi_copy_ref   =busco_multi_copy,
           busco_fragmented_ref   =busco_fragmented,
           busco_missing_busco_ref=busco_missing_busco,
           busco_avg_identity_ref =busco_avg_identity,
           busco_stop_codon_ref   =busco_stop_codon) %>%
    extract( # Extract relevant columns from JobName
      col  =JobName,
      into =c("Program", "Sample"),
      regex="^(.+?)(?=_GCF)_(GCF.+?)?$"
    ) %>% select(! Program )
  
  # Split the JobID into two parts: the main job ID and the batch job ID
  main <- raw_sacct %>% filter(!str_detect(JobID, "\\.batch$"))
  
  batch <- raw_sacct %>% 
    filter(str_detect(JobID, "\\.batch$")) %>%
    mutate(JobID = str_remove(JobID, "\\.batch$")) %>%
    select(JobID, MaxRSS_batch = MaxRSS)
  
  # Left-join the batch MaxRSS onto the main jobs
  unified <- main %>%
    left_join(batch, by = "JobID") %>%
    mutate(
      MaxRSS    = coalesce(MaxRSS_batch, MaxRSS),
      MaxRSS    = parse_number(MaxRSS),
      MaxRSS_GB = MaxRSS / (1024^2)
    ) %>%
    select(-MaxRSS_batch) %>%
    group_by(JobName) %>%
    slice_max(Start, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    full_join(raw_metrics, by = "JobName") %>% # Join with metrics
    filter(! str_detect(JobName, "reference")) %>%
    extract( # Extract relevant columns from JobName
      col   = JobName,
      into  = c("Program", "Sample", "depthNP", "depthIL"),
      regex = "^(.+?)(?=_GCF)_(GCF.+?)(?:_NP(\\d+x))?(?:_IL(\\d+x))?$"
    ) %>%
    full_join(genome_sizes_df, by="Sample") %>% # Join with genome sizes
    full_join(raw_metrics_reference, by="Sample") %>% # Join with reference metrics
    clean_names() %>% ## filter(start != "None") %>% # Remove jobs that haven't ended
    mutate(
      program     = str_to_lower(program),
      depth_np    = replace_na(parse_number(depth_np), 0),
      depth_il    = replace_na(parse_number(depth_il), 0),
      strategy    = case_when(
        program %in% long_reads_strategy  ~ "long_reads",
        program %in% short_reads_strategy ~ "short_reads",
        program %in% hybrid_strategy      ~ "hybrid",
        program %in% polished_strategy    ~ "polished"
      ),
      elapsed_sec = dplyr::coalesce(suppressWarnings(as.numeric(elapsed)), parse_elapsed(elapsed)),
      elapsed_min    = elapsed_sec / 60,
      max_rss_per_mb = max_rss_gb / ref_size_mbp, # Calculate MaxRSS per Mb of genome size
      elapsed_per_mb = elapsed_sec / ref_size_mbp, # Calculate Elapsed time per Mb of genome size
    ) %>%
    mutate(
      strategy  = factor(strategy, levels = c("short_reads", "long_reads", "hybrid", "polished")),
      program   = factor(program),
      sample    = factor(sample),
      ng50      =  as.numeric(ng50),
      depth_eff = case_when(
        strategy == "short_reads" ~ depth_il,
        strategy == "long_reads"  ~ depth_np,
        strategy == "hybrid"      ~ depth_il + depth_np,
        strategy == "polished"    ~ depth_il + depth_np
      ),
      # Fidelity to truth
      length_err_pct = 100 * (total_length - reference_length) / reference_length,
      gc_err_abs     = abs(gc_percent - reference_gc_percent),
      busco_complete_diff   = busco_complete - busco_complete_ref,
      busco_single_copy_diff= busco_single_copy - busco_single_copy_ref,
      busco_multi_copy_diff = busco_multi_copy - busco_multi_copy_ref
    ) %>%
    rename(
      maxrss_gb    = max_rss_gb,
      time_per_mb  = elapsed_per_mb,
      ram_per_mb   = max_rss_per_mb,
      mism_100kbp  = number_mismatches_per_100_kbp,
      indel_100kbp = number_indels_per_100_kbp,
      misassemblies= number_misassemblies )
  
  
  df_flye_pilon <- combine_programs(unified, "flye", "pilon")
  df_flye_polypolish <- combine_programs(unified, "flye", "polypolish")
  df_flye_racon <- combine_programs(unified, "flye", "racon")
  df_flye_racon_mp2 <- combine_programs(unified, "flye", "racon-mp2")
  
  unified_full <- bind_rows(unified, df_flye_pilon, df_flye_polypolish, df_flye_racon, df_flye_racon_mp2) %>%
    mutate(program = factor(program)) %>%
    filter(!is.na(n50))
  
  return(unified_full)
}

# ----------------------------------------------------------------------------
# Load reference genome metrics for simulated data
# ----------------------------------------------------------------------------

compute_reference_metrics <- function(reference_metrics) {
  
  read_csv(reference_metrics, col_types=cols()) %>% 
    mutate(CDS_mb = CDS / size_mbp,
           gene_density = gene / size_mbp,
           cds_per_gene = CDS / gene,
           intron_mb = intron / size_mbp,
           mRNA_mb = mRNA / size_mbp,
           start_codon_mb = start_codon / size_mbp,
           stop_codon_mb = stop_codon / size_mbp) %>% 
    janitor::clean_names() %>%
    mutate(across(c(size_mbp, gc_content, gene_density, cds_per_gene),
                  ~ as.numeric(scale(.x)), .names = "{.col}_s"))
}

# ----------------------------------------------------------------------------
# Run extraction (when sourced)
# ----------------------------------------------------------------------------

files_dir <- here::here("Simulated", "02_Processing", "R", "files")

data_simulated <- extract_data_simulated(
  sacct_file   = file.path(files_dir, "all_sacct.csv"),
  metrics_file = file.path(files_dir, "quast_busco_results.tsv"),
  genome_sizes = file.path(files_dir, "ref_genome_metrics.csv")
) %>% # Remove two polishing runs with corrupt QUAST/BUSCO outputs caused by failed jobs.
  filter(!(sample == "GCF_009017415.1_ASM901741v1" & depth_np == 100 & 
             depth_il == 25 & program %in% c("flye_polypolish", "polypolish")),
         !(sample == "GCF_000226095.1_ASM22609v1" & depth_np == 50 & 
             depth_il == 75 & program %in% c("flye_pilon", "pilon")),
         !(sample == "GCF_036810415.1_ASM3681041v1" & depth_np == 75 & 
             depth_il == 15 & program %in% c("flye_pilon", "pilon")))

# Compute reference metrics summary (for scaling covariates in models and plotting)
reference_metrics <- compute_reference_metrics(
  file.path(files_dir, "ref_genome_metrics.csv"))
write_csv(reference_metrics,
          file.path(results_dir, "reference_metrics.csv"))

reference_metrics_summary <- reference_metrics %>%
  summarise(
    across(
      c(size_mbp, contigs, gc_content, gene, cds_mb, gene_density,
        cds_per_gene, intron_mb, m_rna_mb, start_codon_mb, stop_codon_mb),
      list(
        min = ~ min(.x, na.rm = TRUE),
        max = ~ max(.x, na.rm = TRUE),
        mean = ~ mean(.x, na.rm = TRUE),
        median = ~ median(.x, na.rm = TRUE),
        sd = ~ sd(.x, na.rm = TRUE)
      )
    )
  ) %>%
  pivot_longer(
    cols = everything(),
    names_to = c("variable", "stat"),
    names_sep = "_(?=[^_]+$)",
    values_to = "value"
  ) %>%
  pivot_wider(names_from = variable, values_from = value) %>%
  mutate(stat = factor(stat, levels = c("min", "max", "mean", "median", "sd"))) %>%
  arrange(stat)
write_csv(reference_metrics_summary,
          file.path(results_dir, "reference_metrics_summary.csv"))

# Compute general stats (mean and CI)
simulated_stats_single <- data_simulated %>%  filter(!is.na(n50)) %>%
  filter(strategy %notin% c("hybrid", "polished")) %>% 
  group_by(strategy, program, depth_eff) %>% 
  summarise(
    n=n(),
    med_n50        =median(n50, na.rm=TRUE),
    med_busco      =median(busco_complete, na.rm=TRUE),
    med_busco_diff =median(busco_complete_diff, na.rm=TRUE),
    med_genfrac    =median(genome_fraction_percent, na.rm=TRUE),
    med_mism       =median(mism_100kbp, na.rm=TRUE),
    med_indel      =median(indel_100kbp, na.rm=TRUE),
    med_len_err    =median(length_err_pct, na.rm=TRUE),
    med_maxrss_gb  =median(maxrss_gb, na.rm=TRUE),
    med_ram_mb     =median(ram_per_mb, na.rm=TRUE),
    med_time_mb    =median(time_per_mb, na.rm=TRUE),
    mean_n50       =mean(n50, na.rm=TRUE),
    mean_ng50      =mean(ng50, na.rm=TRUE),
    mean_busco     =mean(busco_complete, na.rm=TRUE),
    mean_busco_diff=mean(busco_complete_diff, na.rm=TRUE),
    mean_genfrac   =mean(genome_fraction_percent, na.rm=TRUE),
    mean_mism      =mean(mism_100kbp, na.rm=TRUE),
    mean_indel     =mean(indel_100kbp, na.rm=TRUE),
    mean_len_err   =mean(length_err_pct, na.rm=TRUE),
    mean_maxrss_gb =mean(maxrss_gb, na.rm=TRUE),
    mean_ram_mb    =mean(ram_per_mb, na.rm=TRUE),
    mean_time_mb   =mean(time_per_mb, na.rm=TRUE),
    mean_time      =mean(elapsed_sec, na.rm=TRUE),
    CI_n50         =CI95(n50),
    CI_busco       =CI95(busco_complete),
    CI_busco_diff  =CI95(busco_complete_diff),
    CI_genfrac     =CI95(genome_fraction_percent),
    CI_mism        =CI95(mism_100kbp),
    CI_indel       =CI95(indel_100kbp),
    CI_len_err     =CI95(length_err_pct),
    CI_maxrss_gb   =CI95(maxrss_gb),
    CI_ram_mb      =CI95(ram_per_mb),
    CI_time_mb     =CI95(time_per_mb),
    CI_time        =CI95(elapsed_sec),
    .groups="drop"
  ) %>%
  mutate(color=colors_programs[as.character(program)]) %>% 
  arrange(depth_eff, program)

write_csv(simulated_stats_single,
          file.path(results_dir, "simulated_stats_single.csv"))

simulated_stats_multi <- data_simulated %>% filter(!is.na(n50)) %>% 
  group_by(strategy, program, depth_np, depth_il) %>%
  summarise(
    n=n(),
    med_n50        =median(n50, na.rm=TRUE),
    med_busco      =median(busco_complete, na.rm=TRUE),
    med_busco_diff =median(busco_complete_diff, na.rm=TRUE),
    med_genfrac    =median(genome_fraction_percent, na.rm=TRUE),
    med_mism       =median(mism_100kbp, na.rm=TRUE),
    med_indel      =median(indel_100kbp, na.rm=TRUE),
    med_len_err    =median(length_err_pct, na.rm=TRUE),
    med_maxrss_gb  =median(maxrss_gb, na.rm=TRUE),
    med_ram_mb     =median(ram_per_mb, na.rm=TRUE),
    med_time_mb    =median(time_per_mb, na.rm=TRUE),
    mean_n50       =mean(n50, na.rm=TRUE),
    mean_ng50      =mean(ng50, na.rm=TRUE),
    mean_busco     =mean(busco_complete, na.rm=TRUE),
    mean_busco_diff=mean(busco_complete_diff, na.rm=TRUE),
    mean_genfrac   =mean(genome_fraction_percent, na.rm=TRUE),
    mean_mism      =mean(mism_100kbp, na.rm=TRUE),
    mean_indel     =mean(indel_100kbp, na.rm=TRUE),
    mean_len_err   =mean(length_err_pct, na.rm=TRUE),
    mean_maxrss_gb =mean(maxrss_gb, na.rm=TRUE),
    mean_ram_mb    =mean(ram_per_mb, na.rm=TRUE),
    mean_time_mb   =mean(time_per_mb, na.rm=TRUE),
    mean_time      =mean(elapsed_sec, na.rm=TRUE),
    CI_n50         =CI95(n50),
    CI_busco       =CI95(busco_complete),
    CI_busco_diff  =CI95(busco_complete_diff),
    CI_genfrac     =CI95(genome_fraction_percent),
    CI_mism        =CI95(mism_100kbp),
    CI_indel       =CI95(indel_100kbp),
    CI_len_err     =CI95(length_err_pct),
    CI_maxrss_gb   =CI95(maxrss_gb),
    CI_ram_mb      =CI95(ram_per_mb),
    CI_time_mb     =CI95(time_per_mb),
    CI_time        =CI95(elapsed_sec),
    .groups="drop"
  ) %>% 
  mutate(color=colors_programs[as.character(program)]) %>% 
  filter(strategy %in% c("hybrid", "polished"), n > 20) %>% # For plotting only, omit depth/program combinations with very low replication.
  arrange(depth_np, depth_il, program)                      # Models retain all available observations.

write_csv(simulated_stats_multi,
          file.path(results_dir, "simulated_stats_multi.csv"))

message("Data preparation complete. Outputs in: ", results_dir)
save.image(paste(results_dir, "01_data_prep.RData", sep = "/"))