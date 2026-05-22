# ============================================================================
# 01_data_prep.R
# ----------------------------------------------------------------------------
# Empirical data extraction, cleaning, and computation of the per-isolate
# genome-size covariate used in the mixed models.
#
# This script is sourced by:
#   - Empirical/02_Processing/stats_figures.qmd
#   - Empirical/R/02_fit_models.R
#
# Inputs (all under Empirical/02_Processing/files/):
#   - all_sacct_empirical.csv              SLURM accounting records
#   - quast_busco_results_empirical.tsv    Combined QUAST + BUSCO metrics
#   - empirical_samples_info.tsv           Per-sample read and KMC summaries
#
# Outputs:
#   - data_empirical (in-memory)           Cleaned model-ready tibble
#   - genome_size_proxy (in-memory)   Per-isolate genome-size estimates
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
      strategy           = "polished",
      program            = combo_name,
      sample,
      depth_np,
      depth_il           = depth_il.x,
      n50                = n50.x,
      busco_complete     = busco_complete.x,
      busco_single_copy  = busco_single_copy.x,
      busco_multi_copy   = busco_multi_copy.x,
      total_length       = total_length.x,
      maxrss_gb          = pmax(maxrss_gb.x, maxrss_gb.y),
      ram_mb         = pmax(ram_mb.x, ram_mb.y),
      time_mb        = time_mb.x + time_mb.y,
      elapsed_sec        = elapsed_sec.x + elapsed_sec.y
    )
}

# ----------------------------------------------------------------------------
# Main extraction function
# ----------------------------------------------------------------------------

extract_data_empirical <- function(sacct_file, metrics_file, reads_info) {
  
  # Define strategies
  long_reads_strategy  <- c("flye", "flye-ovl1000", "flye-ovl1500", "flye-ovl2000",
                            "flye-ovl2500", "flye-ovl3000", "hifiasm", "canu",
                            "raven", "nextdenovo", "miniasm")
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
  
  reads_info_df <- read_tsv(reads_info, col_types = cols()) %>%
    rename(Sample = sample)
  
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
      regex = "^(.+?)(?=_BS|_LLP)_(BS.+?|LLP.+?)(?:_NP(.[^_]+))?(?:_IL(.+))?$"
    ) %>%
    clean_names() %>%
    filter(!is.na(total_length)) %>%
    mutate(
      program     = str_to_lower(program),
      depth_np    = str_replace_all(depth_np, "OG", "99x"),
      depth_il    = str_replace_all(depth_il, "OG", "99x"),
      depth_np    = replace_na(parse_number(depth_np), 0),
      depth_il    = replace_na(parse_number(depth_il), 0),
      strategy    = case_when(
        program %in% long_reads_strategy  ~ "long_reads",
        program %in% short_reads_strategy ~ "short_reads",
        program %in% hybrid_strategy      ~ "hybrid",
        program %in% polished_strategy    ~ "polished"
      ),
      elapsed_sec = dplyr::coalesce(suppressWarnings(as.numeric(elapsed)),
                                    parse_elapsed(elapsed)),
      elapsed_min    = elapsed_sec / 60,
      max_rss_per_mb = max_rss_gb / (total_length * 1e-6), # Calculate MaxRSS per Mb of assembly size
      elapsed_per_mb = elapsed_sec / (total_length * 1e-6) # Calculate Elapsed time per Mb of assembly size
    ) %>%
    mutate(
      strategy = factor(strategy, levels = c("short_reads", "long_reads", "hybrid", "polished")),
      program  = factor(program),
      sample   = factor(sample),
      depth_eff = case_when(
        strategy == "short_reads" ~ depth_il,
        strategy == "long_reads"  ~ depth_np,
        strategy == "hybrid"      ~ depth_il + depth_np,
        strategy == "polished"    ~ depth_il + depth_np
      )
    ) %>%
    rename(
      maxrss_gb   = max_rss_gb,
      ram_mb  = max_rss_per_mb,
      time_mb = elapsed_per_mb,
      mism_100kbp = number_mismatches_per_100_kbp
    )
  
  # Build polished combination rows (Flye + Polypolish)
  df_flye_polypolish <- combine_programs(unified, "flye-ovl1000", "polypolish")
  unified_full <- bind_rows(unified, df_flye_polypolish) %>%
    mutate(program = factor(program)) %>%
    filter(!is.na(n50))
  
  return(unified_full)
}

# ----------------------------------------------------------------------------
# Compute per-isolate genome-size covariate
# ----------------------------------------------------------------------------
# We use the mean assembly length across the five Flye min-overlap settings
# (1000, 1500, 2000, 2500, 3000 bp) at full sequencing depth (depth_np = 99).
# This estimate is computed once per isolate and applied across all rows for
# that isolate, providing a per-sample genome-size covariate that does not
# vary with program or depth_eff.

compute_genome_size_proxy <- function(unified) {
  flye_variants <- c("flye-ovl1000", "flye-ovl1500", "flye-ovl2000",
                     "flye-ovl2500", "flye-ovl3000")
  
  unified %>%
    filter(program %in% flye_variants, depth_np == 99) %>%
    group_by(sample) %>%
    summarise(
      genome_size_est = mean(total_length, na.rm = TRUE),
      n_flye_variants = sum(!is.na(total_length)),
      .groups = "drop"
    )
}

# ----------------------------------------------------------------------------
# Run extraction (when sourced)
# ----------------------------------------------------------------------------

files_dir <- here("Empirical", "02_Processing", "R", "files")

data_empirical_unfiltered <- extract_data_empirical(
  sacct_file   = file.path(files_dir, "all_sacct_empirical.csv"),
  metrics_file = file.path(files_dir, "quast_busco_results_empirical.tsv"),
  reads_info   = file.path(files_dir, "empirical_samples_info.tsv")
)

# Only keep isolates with sufficient data coverage (matching original logic)
max_assemblies_per_sample <- data_empirical_unfiltered %>%
  filter(strategy %notin% c("hybrid", "polished")) %>%
  group_by(sample) %>%
  summarise(n = n(), .groups = "drop") %>%
  pull(n) %>%
  max()

complete_samples <- data_empirical_unfiltered %>%
  filter(strategy %notin% c("hybrid", "polished")) %>%
  group_by(sample) %>%
  summarise(n = n(), .groups = "drop") %>%
  filter(n >= 78) %>%
  pull(sample)

data_empirical <- data_empirical_unfiltered %>%
  filter(sample %in% complete_samples, program != "flye") %>%
  filter_out(sample=="BS08_360_127-6" & program=="flye-ovl1500" & depth_np=="50")

metrics_empirical <- read_tsv(file.path(files_dir, "empirical_samples_info.tsv"), col_types=cols()) %>% 
  filter(sample %in% complete_samples) %>%
  mutate(program= case_when(
    technology == "illumina" ~ "kmc_il",
    technology == "chopper"  ~ "kmc_chopper",
    TRUE ~ "ukn"
  ))

# Compute genome-size covariate and join into modeling data
genome_size_proxy <- compute_genome_size_proxy(data_empirical)

data_empirical <- data_empirical %>%
  left_join(genome_size_proxy, by = "sample") %>%
  mutate(
      ram_mb_p  = maxrss_gb/(genome_size_est * 1e-6), # Recalculate RAM per Mb using genome size proxy  
      time_mb_p = elapsed_sec/(genome_size_est * 1e-6)  # Recalculate time per Mb using genome size proxy
      ) %>%
  mutate(genome_size_est_s = as.numeric(scale(genome_size_est)))

# Sanity check: every row should have a genome_size_est value
n_missing <- sum(is.na(data_empirical$genome_size_est))
if (n_missing > 0) {
  warning(sprintf(
    "%d rows in data_empirical have no genome_size_est. Check that all isolates have Flye OVL assemblies at OG depth.",
    n_missing
  ))
}

empirical_stats_single <- data_empirical %>%  filter(!is.na(n50)) %>%
  filter(strategy %notin% c("hybrid", "polished"),
         sample %in% complete_samples) %>%
  # filter(depth_np != 99, depth_il != 99) %>%
  mutate(number_contigs_small = number_contigs_0_bp - number_contigs_1000_bp) %>%
  group_by(strategy, program, depth_eff) %>% 
  summarise(
    n=n(),
    med_n50        =median(n50, na.rm=TRUE),
    med_busco      =median(busco_complete, na.rm=TRUE),
    med_busco_copy =median(busco_multi_copy, na.rm=TRUE),
    med_maxrss_gb  =median(maxrss_gb, na.rm=TRUE),
    med_ram_mb     =median(ram_mb, na.rm=TRUE),
    med_time_mb    =median(time_mb, na.rm=TRUE),
    med_sht_contig =median(number_contigs_small, na.rm=TRUE),
    med_length     =median(total_length, na.rm=TRUE),
    med_cov_depth  =median(avg_coverage_depth, na.rm=TRUE),
    mean_n50       =mean(n50, na.rm=TRUE),
    mean_busco     =mean(busco_complete, na.rm=TRUE),
    mean_busco_copy=mean(busco_multi_copy, na.rm=TRUE),
    mean_maxrss_gb =mean(maxrss_gb, na.rm=TRUE),
    mean_ram_mb    =mean(ram_mb, na.rm=TRUE),
    mean_time_mb   =mean(time_mb, na.rm=TRUE),
    mean_ram_mb_p  =mean(ram_mb_p, na.rm=TRUE),
    mean_time_mb_p =mean(time_mb_p, na.rm=TRUE),
    mean_time      =mean(elapsed_sec, na.rm=TRUE),
    mean_sht_contig=mean(number_contigs_small, na.rm=TRUE),
    mean_length    =mean(total_length, na.rm=TRUE),
    mean_cov_depth =mean(avg_coverage_depth, na.rm=TRUE),
    CI_n50         =CI95(n50),
    CI_busco       =CI95(busco_complete),
    CI_busco_copy  =CI95(busco_multi_copy),
    CI_maxrss_gb   =CI95(maxrss_gb),
    CI_ram_mb      =CI95(ram_mb),
    CI_time_mb     =CI95(time_mb),
    CI_ram_mb_p    =CI95(ram_mb_p),
    CI_time_mb_p   =CI95(time_mb_p),
    CI_time        =CI95(elapsed_sec),
    CI_sht_contig  =CI95(number_contigs_small),
    CI_length      =CI95(total_length),
    CI_cov_depth   =CI95(avg_coverage_depth),
    .groups="drop"
  ) %>%
  arrange(depth_eff, program)

write_csv(empirical_stats_single,
          file.path(results_dir, "empirical_stats_single.csv"))

empirical_stats_multi <- data_empirical %>% filter(!is.na(mism_100kbp)) %>% 
  filter(strategy %in% c("hybrid", "polished"),
         sample %in% complete_samples) %>% 
  mutate(number_contigs_small = number_contigs_0_bp - number_contigs_1000_bp) %>%
  group_by(strategy, program, depth_np, depth_il) %>%
  summarise(
    n=n(),
    med_n50        =median(n50, na.rm=TRUE),
    med_busco      =median(busco_complete, na.rm=TRUE),
    med_busco_copy =median(busco_multi_copy, na.rm=TRUE),
    med_maxrss_gb  =median(maxrss_gb, na.rm=TRUE),
    med_ram_mb     =median(ram_mb, na.rm=TRUE),
    med_time_mb    =median(time_mb, na.rm=TRUE),
    med_sht_contig =median(number_contigs_small, na.rm=TRUE),
    med_length     =median(total_length, na.rm=TRUE),
    med_mism       =median(mism_100kbp, na.rm=TRUE),
    mean_n50       =mean(n50, na.rm=TRUE),
    mean_busco     =mean(busco_complete, na.rm=TRUE),
    mean_busco_copy=mean(busco_multi_copy, na.rm=TRUE),
    mean_maxrss_gb =mean(maxrss_gb, na.rm=TRUE),
    mean_ram_mb    =mean(ram_mb, na.rm=TRUE),
    mean_time_mb   =mean(time_mb, na.rm=TRUE),
    mean_ram_mb_p  =mean(ram_mb_p, na.rm=TRUE),
    mean_time_mb_p =mean(time_mb_p, na.rm=TRUE),
    mean_time      =mean(elapsed_sec, na.rm=TRUE),
    mean_sht_contig=mean(number_contigs_small, na.rm=TRUE),
    mean_length    =mean(total_length, na.rm=TRUE),
    mean_mism      =mean(mism_100kbp, na.rm=TRUE),
    CI_n50         =CI95(n50),
    CI_busco       =CI95(busco_complete),
    CI_busco_copy  =CI95(busco_multi_copy),
    CI_maxrss_gb   =CI95(maxrss_gb),
    CI_ram_mb      =CI95(ram_mb),
    CI_time_mb     =CI95(time_mb),
    CI_ram_mb_p    =CI95(ram_mb_p),
    CI_time_mb_p   =CI95(time_mb_p),
    CI_time        =CI95(elapsed_sec),
    CI_sht_contig  =CI95(number_contigs_small),
    CI_length      =CI95(total_length),
    CI_mism        =CI95(mism_100kbp),
    .groups="drop"
  ) %>% 
  filter(strategy %in% c("hybrid", "polished")) %>% 
  arrange(depth_np, depth_il, program)

write_csv(empirical_stats_multi,
          file.path(results_dir, "empirical_stats_multi.csv"))
message("Data preparation complete. Outputs in: ", results_dir)

save.image(file.path(results_dir, "01_data_prep.RData"))
