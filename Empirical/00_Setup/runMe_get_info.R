#!/usr/bin/env Rscript

library(tidyverse)
library(jsonlite)

dir <- "/scratch/gas0042/nanopore_benchmark/real_data/data/filtered/original"
samples <- list.dirs(dir, recursive = FALSE, full.names = TRUE) %>%
                basename()

parse_genomescope <- function(sample, tech) {
    dir <- paste0("/scratch/gas0042/nanopore_benchmark/real_data/results/kmc-genomescope_", tech)
    file <- file.path(dir, sample, "summary.txt")

    if (!file.exists(file)) {
        tibble(
        min_bp          = NA,
        max_bp          = NA,
        mean_bp         = NA,
        min_fit         = NA,
        max_fit         = NA,
        mean_fit        = NA,
        min_error_rate  = NA,
        max_error_rate  = NA,
        mean_error_rate = NA
    )
    } else {
        lines <- readLines(file)

        # extract numbers (strip commas & bp)
        nums <- lines[grepl("Genome Haploid Length", lines)] %>%
            str_extract_all("[0-9,]+") %>% 
            unlist() %>% 
            str_replace_all(",", "") %>% 
            as.numeric()
        
        fit <- lines[grepl("Model Fit", lines)] %>%
            str_extract_all("[0-9.]+") %>% 
            unlist() %>% 
            as.numeric()
        
        error <- lines[grepl("Read Error Rate", lines)] %>%
            str_extract_all("[0-9.]+") %>% 
            unlist() %>% 
            as.numeric()

        tibble(
            min_bp = nums[1],
            max_bp = nums[2],
            mean_bp = mean(nums[1:2]),
            min_fit = fit[1],
            max_fit = fit[2],
            mean_fit = mean(fit[1:2]),
            min_error_rate = error[1],
            max_error_rate = error[2],
            mean_error_rate = mean(error[1:2])
        )
    }
}

parse_fastp <- function(sample) {
    dir <- "/scratch/gas0042/nanopore_benchmark/real_data/data/filtered/original"
    file <- file.path(dir, sample, paste0(sample,"_fastp.json"))
    json_data <- fromJSON(file)

    total_bases <- json_data$summary$after_filtering$total_bases
    mean_length <- json_data$summary$after_filtering$read1_mean_length
    
    return(c(total_bases, mean_length))
}

parse_fastq <- function(sample,tech) {
    dir <- "/scratch/gas0042/nanopore_benchmark/real_data/data/subsampled/"
    if (tech == "np") {
        type <- "nanopore"
    } else if ( str_detect(tech, "chopper") ) {
        type <- "chopper"
    }
    file <- file.path(dir, sample, paste0(sample,"_OG.",type,".fq"))

    if (!file.exists(file)) {
        return(c(NA, NA))
    } else {
        # Calculate total bases and mean read length
        total_bases <- system(paste("awk 'NR%4==2' ", file, " | tr -d '\n' | wc -c"), intern=TRUE) %>% as.numeric()
        # total_reads <- system(paste("grep -c '^@AV' ", file), intern=TRUE) %>% as.numeric()
        total_reads <- system(paste("expr $(cat", file, "| wc -l) / 4"), intern=TRUE) %>% as.numeric()
        mean_length <- total_bases / total_reads

        return(c(total_bases, mean_length))
    }
}

parse_file <- function(sample, tech) {
    genomescope_nums <- parse_genomescope(sample, tech)

    if (tech == "illumina") {
        fastp_nums  <- parse_fastp(sample)
        total_bases <- fastp_nums[1]
        mean_length <- fastp_nums[2]
        depth       <- total_bases / genomescope_nums$mean_bp
    } else if (tech == "np" || str_detect(tech, "chopper")) {
        fastq_nums  <- parse_fastq(sample, tech)
        total_bases <- fastq_nums[1]
        mean_length <- fastq_nums[2]
        depth       <- total_bases / genomescope_nums$mean_bp
    }
    

    tibble(
        sample          = sample,
        technology      = tech,
        min_bp          = genomescope_nums$min_bp,
        max_bp          = genomescope_nums$max_bp,
        mean_bp         = genomescope_nums$mean_bp,
        min_fit         = genomescope_nums$min_fit,
        max_fit         = genomescope_nums$max_fit,
        mean_fit        = genomescope_nums$mean_fit,
        min_error_rate  = genomescope_nums$min_error_rate,
        max_error_rate  = genomescope_nums$max_error_rate,
        mean_error_rate = genomescope_nums$mean_error_rate,
        total_bases     = total_bases,
        mean_length     = mean_length,
        depth           = depth
    )
}

parse_technologies <- function(sample) {
    res <- tibble()
    for (tech in c("illumina", "np", "chopper", "chopperk101", "chopperk13", "chopperk13p1", "chopperk21p1", "chopperk57p1")) {
        parse_file(sample, tech) %>%
            bind_rows(res) -> res
    }
    return(res)
}

results <- map_dfr(samples, parse_technologies)

# Save results to a TSV file
write_tsv(results, "/scratch/gas0042/nanopore_benchmark/real_data/results/real_samples_info.tsv")

# Calculate downsampling factors for various target coverages
results %>% 
  select(sample, depth, technology) %>%
  arrange(depth) %>%
  mutate(`10x`  = 10  / depth,
         `15x`  = 15  / depth,
         `20x`  = 20  / depth,
         `25x`  = 25  / depth,
         `30x`  = 30  / depth,
         `35x`  = 35  / depth,
         `40x`  = 40  / depth,
         `50x`  = 50  / depth,
         `60x`  = 60  / depth,
         `75x`  = 75  / depth,
         `100x` = 100 / depth ) %>%
  pivot_longer(cols=ends_with("x"),
               names_to="target_coverage",
               values_to="coverage_factor") %>%
  write_csv("/scratch/gas0042/nanopore_benchmark/real_data/results/downsampling_factors.csv")
