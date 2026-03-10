#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(tidyverse)
    library(jsonlite)
})

message("Starting incremental aggregator")

# ------------ Paths ------------
root_dir <- "/scratch/gas0042/nanopore_benchmark/real_data"
results_quast <- file.path(root_dir, "results", "quast_results")
results_busco <- file.path(root_dir, "results", "busco_results")

# Unified outputs
quast_out   <- file.path(results_quast, "final_unified_report.tsv")
busco_out   <- file.path(results_busco, "final_unified_summary.tsv")
merged_out  <- file.path(root_dir, "results", "quast_busco_results_real.tsv")

# Manifests (hidden files alongside outputs)
quast_manifest <- file.path(results_quast, ".quast_manifest.tsv")
busco_manifest <- file.path(results_busco, ".busco_manifest.tsv")

# ------------ Helpers ------------
scan_manifest <- function(paths) {
    # Current snapshot of files
    tibble(
        path = paths,
        mtime = as.numeric(file.info(paths)$mtime),
        size  = as.numeric(file.info(paths)$size)
    ) %>% drop_na()
}

load_manifest <- function(manifest_path) {
    if (file.exists(manifest_path)) {
        suppressWarnings(read_tsv(manifest_path, show_col_types = FALSE,
                                    col_types = cols(path="c", mtime="d", size="d")))
    } else {
        tibble(path = character(), mtime = numeric(), size = numeric())
    }
}

save_manifest <- function(df, manifest_path) {
    write_tsv(df, manifest_path)
}

load_unified <- function(path) {
    if (file.exists(path)) {
        suppressWarnings(read_tsv(path, show_col_types = FALSE, col_types = cols(.default = "c")))
    } else {
        tibble(Assembly = character())
    }
}

# Parse QUAST file -> tibble with Assembly column
read_quast_file <- function(p) {
    df <- suppressWarnings(read_tsv(p, show_col_types = FALSE, col_types = cols(.default = "c")))
    if (!"Assembly" %in% names(df)) {
        # Fallback: derive Assembly from filename ".../<PREFIX>_transposed_report.tsv"
        asm <- sub("_transposed_report\\.tsv$", "", basename(p))
        df <- df %>% mutate(Assembly = asm, .before = 1)
    }
    df
}

# Parse BUSCO JSON -> tibble with Assembly column
read_busco_json <- function(jpath) {
    x <- fromJSON(jpath)
    # if (is.null(x$results$avg_identity)) {
    #     message("Skipping empty BUSCO file: ", jpath)
    #     return(tibble(Assembly = character()))
    # } else {
    assembly <- str_split(x$parameters$out, "/")[[1]] %>% dplyr::last()
    as_tibble(x$results) %>% mutate_all(as.character) %>%
        mutate(Assembly = assembly, .before = 1) %>%
        mutate(busco_avg_identity  = ifelse(exists("avg_identity"), avg_identity, "0")) %>%
        rename(
            busco_complete      = `Complete percentage`,
            busco_single_copy   = `Single copy percentage`,
            busco_multi_copy    = `Multi copy percentage`,
            busco_fragmented    = `Fragmented percentage`,
            busco_missing_busco = `Missing percentage`,
            busco_stop_codon    = internal_stop_codon_percent
        ) %>%
        select(Assembly, busco_complete, busco_single_copy, busco_multi_copy,
                busco_fragmented, busco_missing_busco, busco_avg_identity, busco_stop_codon)
    # }
}

# Generic incremental builder:
# - find files matching pattern
# - compare to manifest
# - process only new/changed
# - drop/replace rows for updated assemblies in the unified table
incremental_build <- function(
    files,
    manifest_path,
    unified_path,
    reader_fun,         # function(path) -> tibble (MUST include "Assembly")
    key_col = "Assembly"
) {
    # Current file snapshot
    current <- scan_manifest(files)

    # Prior state
    prev_manifest <- load_manifest(manifest_path)
    unified_prev  <- load_unified(unified_path)

    # Figure out which files are new or modified
    to_join <- current %>% left_join(prev_manifest, by = "path", suffix = c("", ".old"))
    todo <- to_join %>% filter(is.na(mtime.old) | mtime > mtime.old | size != size.old)

    if (nrow(todo) == 0) {
        message("No new/modified files for: ", unified_path)
        # Ensure output exists even on first run
        if (!file.exists(unified_path)) write_tsv(unified_prev, unified_path)
        # Update manifest to current snapshot (also covers deletions/renames if you later want to clean)
        save_manifest(current, manifest_path)
        return(unified_prev)
    }

    message("Processing ", nrow(todo), " new/modified files for: ", unified_path)

    # Read only the TODO set
    new_rows <- map(todo$path, reader_fun) %>% list_rbind()

    # Which assemblies are being updated?
    updated_keys <- new_rows %>% distinct(.data[[key_col]])

    # Drop old rows for those assemblies, then append fresh rows
    unified_next <- unified_prev %>%
        anti_join(updated_keys, by = key_col) %>%
        bind_rows(new_rows)

    # Persist unified + manifest
    write_tsv(unified_next, unified_path)
    save_manifest(current, manifest_path)

    unified_next
}

# ------------ QUAST (incremental) ------------
message("QUAST: scanning")
quast_files <- list.files(
    path = results_quast,
    full.names = TRUE,
    pattern = "_transposed_report\\.tsv$",
    recursive = TRUE
)

quast_unified <- incremental_build(
    files        = quast_files,
    manifest_path = quast_manifest,
    unified_path  = quast_out,
    reader_fun    = read_quast_file,
    key_col       = "Assembly"
)

# ------------ BUSCO (incremental) ------------
message("BUSCO: scanning")
busco_files <- list.files(
    path = results_busco,
    full.names = TRUE,
    pattern = "short_summary.*\\.json$",
    recursive = TRUE
)

busco_unified <- incremental_build(
    files         = busco_files,
    manifest_path = busco_manifest,
    unified_path  = busco_out,
    reader_fun    = read_busco_json,
    key_col       = "Assembly"
)

# ------------ Merge (cheap) ------------
message("Merging QUAST + BUSCO")
final_results <- quast_unified %>%
    full_join(busco_unified, by = "Assembly") %>%
    select(Assembly, everything())

write_tsv(final_results, merged_out)

message("Done!")