# ============================================================================
# 03_figure.R
# ----------------------------------------------------------------------------
# Creates individual plots to be combined into publication ready figures
#
# Provides:
#   - Individual Plots for all the responses
#   - Forest plots for the fixed effects and slopes
#
# Sources:
#   - 01_fit_models.R  (uses statistical summaries from all the points)
#   - 02_fit_models.R  (uses fitted model objects in single_fits / multi_fits)
# ============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggdist)
  library(ggpmisc)
  library(ggExtra)
  library(see)
  library(ggbreak)
  library(paletteer)
  library(ggpubr)
  library(scales)
  library(ggtext)
})
#"janitor", "scales", "rlang", "ggrepel","ggpp", "ggpubr", "ggExtra", "ggpmisc",
color_palette <- paletteer_d("ggthemes::few_Dark")[-1]

response_labels <- c(
  "n50"                     = "N50",
  "busco_complete_diff"     = "BUSCO difference",
  "busco_complete"          = "BUSCO completeness",
  "genome_fraction_percent" = "Genome recovery",
  "length_err_pct"          = "Length error",
  "maxrss_gb"               = "Max RAM",
  "ram_per_mb"              = "Max RAM (/Mbp)",
  "elapsed_sec"             = "Run Time",
  "time_per_mb"             = "Run Time(/Mbp)",
  "mism_100kbp"             = "Mismatches (/100kb)",
  "indel_100kbp"            = "Indels (/100kb)"
)

# Set labels for plotting
program_labels <- c("flye"="Flye",
                    "hifiasm"="Hifiasm",
                    "canu"="Canu",
                    "raven"="Raven",
                    "nextdenovo"="NextDenovo",
                    "miniasm"="Miniasm",
                    "spades_short"="SPAdes Short",
                    "abyss_short"="ABySS Short",
                    "abyss_hybrid"="ABySS Hybrid",
                    "masurca"="MaSuRCA",
                    "spades_hybrid"="SPAdes Hybrid",
                    "pilon"="Pilon",
                    "polypolish"="Polypolish",
                    "racon-mp2"="Racon (Minimap2)",
                    "racon"="Racon (BWA)",
                    "flye_pilon"="Flye + Pilon",
                    "flye_polypolish"="Flye + Polypolish",
                    "flye_racon"="Flye + Racon (BWA)",
                    "flye_racon-mp2"="Flye + Racon (Minimap2)")

# Set color palette
colors_programs <- c("flye"="#59A14F",
                     "hifiasm"="#1F83B4",
                     "canu"="#C7519C",
                     "raven"="#A52A2A",
                     "nextdenovo"="#8C564B",
                     "miniasm"="#7B68EE",
                     "spades_short"="#EDC948",
                     "abyss_short"="#FF7F0E",
                     
                     "abyss_hybrid"="#FFAA0E",
                     "masurca"="#8B4513",
                     "spades_hybrid"="#CD1076",
                     
                     "pilon"="#B22222",
                     "polypolish"="#00688B",
                     "racon"="#008B00",
                     "racon-mp2"="#FF7F0E",
                     "flye_pilon"="#B22222",
                     "flye_polypolish"="#00688B",
                     "flye_racon"="#008B00",
                     "flye_racon-mp2"="#FF7F0E")

figures_dir <- file.path(results_dir, "figures")
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

save_plot <- function(plt, name , size_w, size_h) {
  formats <- c("png", "pdf", "tiff", "jpeg")
  
  for (ext in formats) {
    ggsave(filename = paste0(figures_dir, "/", name, ".", ext),
           plot = plt, device = ext,
           width = size_w, height = size_h,
           units = "in")
  }
}


# ----------------------------------------------------------------------------
# Figures for Single Seq Assemblies
# ----------------------------------------------------------------------------

program_order <- c("abyss_short", "spades_short", "flye", "hifiasm",
                   "canu", "raven", "nextdenovo", "miniasm")

simulated_stats_single_plt <- simulated_stats_single %>%
  mutate(depth_eff=factor(depth_eff),
         program=factor(program, levels=program_order)) 

#
# Simulated Single Sequencing Metrics
#
# Figure 1

{
# N50 plot
plt_n50_single <- simulated_stats_single_plt %>%
  arrange(CI_n50, program) %>% 
  ggplot(aes(x=depth_eff, y=mean_n50, color=program)) +
  geom_line(aes(group=program), position = "identity") +
  geom_point(size=2, position = "identity") +
  geom_linerange(aes(ymin=(mean_n50 - CI_n50), ymax=(mean_n50 + CI_n50)),
                 linewidth=0.5, show.legend = F, position = "identity") +
  scale_y_continuous(labels=scales::label_number(scale=1e-6), breaks = 1e6*c(0,1,2,3,4)) +
  scale_color_manual(values=colors_programs, label=program_labels, name="Program") +
  labs(x="Depth (x)", y="N50 (Mbp)") +
  theme_minimal(base_size=14) +
  theme(panel.grid.minor=element_blank(),
        panel.grid.major.x=element_blank()) 

# Genome Fraction plot
plt_genfrac_single <- simulated_stats_single_plt %>% 
  arrange(CI_genfrac, program) %>%
  ggplot(aes(x=depth_eff, y=mean_genfrac, color=program)) +
  geom_line(aes(group=program), position = "identity") +
  geom_point(size=2, position = "identity") +
  geom_linerange(aes(ymin=(mean_genfrac - CI_genfrac), ymax=(mean_genfrac + CI_genfrac)),
                 linewidth=0.5, show.legend = F, position = "identity") +
  scale_color_manual(values=colors_programs, label=program_labels, name="Program") +
  labs(x="Depth (x)", y="Recovered (%)") +
  theme_minimal(base_size=14)+
  theme(panel.grid.minor=element_blank(),
        panel.grid.major.x=element_blank())

# Length error plot
plt_len_err_single <- simulated_stats_single_plt %>%
  arrange(CI_len_err, program) %>%
  ggplot(aes(x=depth_eff, y=mean_len_err, color=program)) +
  geom_hline(yintercept=0, linetype="dashed", color="grey50") +
  geom_line(aes(group=program), position = "identity") +
  geom_point(size=2, position = "identity") +
  geom_linerange(aes(ymin=(mean_len_err - CI_len_err), ymax=(mean_len_err + CI_len_err)),
                 linewidth=0.5, show.legend = F, position = "identity") +
  scale_color_manual(values=colors_programs, label=program_labels, name="Program") +
  labs(x="Depth (x)", y="Length error (%)") +
  theme_minimal(base_size=14) +
  # coord_cartesian(ylim=c(-5,5)) +
  theme(panel.grid.minor=element_blank(),
        panel.grid.major.x=element_blank())

# BUSCO difference plot
plt_busco_diff_single <- simulated_stats_single_plt %>%
  arrange(CI_busco_diff, program) %>%
  ggplot(aes(x=depth_eff, y=mean_busco_diff, color=program)) +
  geom_hline(yintercept=0, linetype="dashed", color="grey50") +
  geom_line(aes(group=program), position = "identity") +
  geom_point(size=2, position = "identity") +
  geom_linerange(aes(ymin=(mean_busco_diff - CI_busco_diff), ymax=(mean_busco_diff + CI_busco_diff)),
                 linewidth=0.5, show.legend = F, position = "identity") +
  scale_color_manual(values=colors_programs, label=program_labels, name="Program") +
  labs(x="Depth (x)", y="BUSCO difference") +
  theme_minimal(base_size=14) +
  theme(panel.grid.minor=element_blank(),
        panel.grid.major.x=element_blank())

# Mismatches plot
plt_mism_single <- simulated_stats_single_plt %>% 
  arrange(CI_mism, program) %>%
  ggplot(aes(x=depth_eff, y=mean_mism, color=program)) +
  geom_line(aes(group=program), position = "identity") +
  geom_point(size=2, position = "identity") +
  geom_linerange(aes(ymin=(mean_mism - CI_mism), ymax=(mean_mism + CI_mism)),
                 linewidth=0.5, show.legend = F, position = "identity") +
  scale_color_manual(values=colors_programs, label=program_labels, name="Program") +
  scale_y_log10(labels = label_number()) +
  labs(x="Depth (x)", y="Mismatches (/100kb)<p style='font-size:10pt'>log10 scaled</p>") +
  theme_minimal(base_size=14)+
  theme(panel.grid.minor=element_blank(),
        panel.grid.major.x=element_blank(),
        axis.title.y = ggtext::element_markdown())

# Indels plot 
plt_indel_single <- simulated_stats_single_plt %>% 
  arrange(CI_indel, program) %>%
  ggplot(aes(x=depth_eff, y=mean_indel, color=program)) +
  geom_line(aes(group=program), position = "identity") +
  geom_point(size=2, position = "identity") +
  geom_linerange(aes(ymin=(mean_indel - CI_indel), ymax=(mean_indel + CI_indel)),
                 linewidth=0.5, show.legend = F, position = "identity") +
  scale_color_manual(values=colors_programs, label=program_labels, name="Program") +
  scale_y_log10(labels = label_number()) +
  labs(x="Depth (x)", y="Indels (/100kb)<p style='font-size:10pt'>log10 scaled</p>") +
  theme_minimal(base_size=14)+
  theme(panel.grid.minor=element_blank(),
        panel.grid.major.x=element_blank(),
        axis.title.y = ggtext::element_markdown())

ggarrange(plt_n50_single + guides(color=guide_legend(nrow=1)),
          plt_genfrac_single, plt_busco_diff_single,
          plt_len_err_single, plt_mism_single, plt_indel_single,
          ncol=2, nrow=3, labels=c("A", "B", "C", "D", "E", "F"),
          common.legend=TRUE, legend="bottom")
} -> plt_single_metrics

save_plot(
  plt_single_metrics,
  "fig01_simulated_single_metrics",
  11,11
)


#
# Simulated Single Sequencing Resources
#
# Figure 2

{ 
# Max RAM per Mbp plot
plt_ram_mb_single <- simulated_stats_single_plt %>%
  arrange(CI_ram_mb, program) %>%
  ggplot(aes(x=depth_eff, y=mean_ram_mb, color=program)) +
  geom_line(aes(group=program), position = "identity") +
  geom_point(size=2, position = "identity") +
  geom_linerange(aes(ymin=(mean_ram_mb - CI_ram_mb), ymax=(mean_ram_mb + CI_ram_mb)),
                 linewidth=0.5, show.legend = F, position = "identity") +
  scale_color_manual(values=colors_programs, label=program_labels, name="Program") +
  labs(x="Depth (x)", y="Max RAM (GB/Mbp)") +
  theme_minimal(base_size=14)+
  theme(panel.grid.minor=element_blank(),
        panel.grid.major.x=element_blank())

# Time per Mbp plot
plt_time_mb_single <- simulated_stats_single_plt %>%
  arrange(CI_time_mb, program) %>%
  ggplot(aes(x=depth_eff, y=mean_time_mb, color=program)) +
  geom_line(aes(group=program), position = "identity") +
  geom_linerange(aes(ymin=(mean_time_mb - CI_time_mb), ymax=(mean_time_mb + CI_time_mb)),
                 linewidth=0.5, show.legend = F, position = "identity") +
  geom_point(size=2, position = "identity") +
  scale_color_manual(values=colors_programs, label=program_labels, name="Program") +
  scale_y_log10(labels = label_log(digits = 2)) +
  labs(x="Depth (x)", y="Time (s/Mbp)<p style='font-size:10pt'>log10 scaled</p>") +
  theme_minimal(base_size=14)+
  theme(panel.grid.minor=element_blank(),
        panel.grid.major.x=element_blank(),
        axis.title.y = ggtext::element_markdown())

plotA <- plt_ram_mb_single + guides(color=guide_legend(nrow=1))
plotB <- plt_time_mb_single

ggarrange(plotA, plotB,
          ncol=2, nrow=1, labels=c("A", "B"),
          common.legend=TRUE, legend="bottom")
} -> plt_single_resources

save_plot(
  plt_single_resources,
  "fig02_simulated_single_resources",
  11,5
)


# ----------------------------------------------------------------------------
# Figures for Polished Assemblies
# ----------------------------------------------------------------------------

flye_polished <- full_join(simulated_stats_single, simulated_stats_multi) %>% 
  mutate(depth_np=ifelse(is.na(depth_np), depth_eff, depth_np),
         depth_il=ifelse(is.na(depth_il), 100, depth_il)) %>%
  filter(program %in% c("flye", "flye_polypolish", "flye_pilon", "flye_racon", "flye_racon-mp2"), depth_np <= 100, depth_il <= 100) %>%
  mutate(program=factor(program, levels=c("flye_racon", "flye_racon-mp2", "flye_polypolish", "flye_pilon","flye"))) %>% 
  arrange(program)

flye_polished_split <- full_join(simulated_stats_single, simulated_stats_multi) %>% 
  mutate(depth_np=ifelse(is.na(depth_np), depth_eff, depth_np),
         depth_il=ifelse(is.na(depth_il), 0, depth_il)) %>%
  filter(depth_np <= 100, depth_il <= 100) %>%
  mutate(program=factor(program)) %>% 
  arrange(program)

colors_tmp <- c(colors_programs)
colors_tmp[["flye"]] <- "#000000"

#
# Simulated Polished Metrics
#
# Figure 3

{
plt_mism_polished <- flye_polished %>%
  mutate(depth_il=as.numeric(depth_il), depth_np=factor(depth_np)) %>%
  ggplot(aes(x=depth_np, y=mean_mism, color=program, alpha=depth_il)) +
  geom_line(aes(group=interaction(program, depth_il)), position = "identity") +
  geom_point(size=2, position = "identity") +
  geom_linerange(aes(ymin=(mean_mism - CI_mism), ymax=(mean_mism + CI_mism)),
                 linewidth=0.5, show.legend = F, position = "identity") +
  scale_alpha_continuous(name="Depth short-reads (x)", range=c(0.25,1), breaks=c(10,30,50,100)) +
  scale_color_manual(values=colors_tmp, label=program_labels, name="Program") +
  scale_y_log10(labels = label_number(), breaks = c(0.3, 1, 3, 6)) +
  labs(x="Depth long-reads (x)", y="Mismatches (/100kb)<p style='font-size:10pt'>log10 scaled</p>") +
  theme_minimal(base_size=14)+
  theme(panel.grid.minor=element_blank(),
        panel.grid.major.x=element_blank(),
        legend.box="vertical", legend.margin=margin(),
        axis.title.y = ggtext::element_markdown())

plt_indel_polished <- flye_polished %>%
  mutate(depth_il=as.numeric(depth_il), depth_np=factor(depth_np)) %>%
  ggplot(aes(x=depth_np, y=mean_indel, color=program, alpha=depth_il)) +
  geom_line(aes(group=interaction(program, depth_il)), position = "identity") +
  geom_point(size=2, position = "identity") +
  geom_linerange(aes(ymin=(mean_indel - CI_indel), ymax=(mean_indel + CI_indel)),
                 linewidth=0.5, show.legend = F, position = "identity") +
  scale_alpha_continuous(name="Depth short-reads (x)", range=c(0.25,1), breaks=c(10,30,50,100)) +
  scale_color_manual(values=colors_tmp, label=program_labels, name="Program") +
  scale_y_log10(labels = label_number(), breaks = c(0.3, 1, 3, 6)) +
  labs(x="Depth long-reads (x)", y="Indels (/100kb)<p style='font-size:10pt'>log10 scaled</p>") +
  theme_minimal(base_size=14)+
  theme(panel.grid.minor=element_blank(),
        panel.grid.major.x=element_blank(),
        legend.box="vertical", legend.margin=margin(),
        axis.title.y = ggtext::element_markdown())

ggarrange(plt_mism_polished + guides(color=guide_legend(nrow=1, reverse=T)),
          plt_indel_polished,
          ncol=2, labels=c("A","B"),
          common.legend=TRUE, legend="bottom")
} -> plt_polished_metrics

save_plot(
  plt_polished_metrics,
  "fig03_simulated_polished_metrics",
  10,5
)

#
# Simulated Polished Resources RAM
#
# Figure 4
{
plt_polished_split_ram_mb <- function(program_name,title){
  flye_program <- paste0("flye_", program_name)
  
  flye_polished_split %>% filter(program %in% c("flye", !!flye_program, !!program_name)) %>% 
    mutate(program=factor(program, levels = c(!!flye_program, !!program_name, "flye")),
           depth_il=factor(depth_il), depth_np=as.numeric(depth_np)) %>%
    arrange(program) %>%
    ggplot(aes(x=depth_il, y=mean_ram_mb, color=program, alpha=depth_np)) +
    geom_line(aes(group=interaction(program, depth_np))) +
    geom_point(size=2) +
    geom_linerange(aes(ymin=(mean_ram_mb - CI_ram_mb), ymax=(mean_ram_mb + CI_ram_mb)),
                   linewidth=0.5, show.legend = F) +
    scale_alpha_continuous(name="Depth long-reads (x)", range=c(0.25,1), breaks=c(10,30,50,100)) +
    scale_color_manual(values=c("olivedrab4","firebrick3","black"), label=c("Flye + Polisher", "Polisher", "Flye"), name="Program") +
    scale_y_log10(labels = label_number(), breaks = c(0.05, 0.1, 0.2, 0.5, 0.7), limits = c(0.04,0.8)) +
    labs(x="Depth short-reads (x)", title=title,
         y="Max RAM (GB/Mbp)<p style='font-size:10pt'>log10 scaled</p>") +
    theme_minimal(base_size=14)+
    theme(panel.grid.minor=element_blank(),
          panel.grid.major.x=element_blank(),
          axis.text.x = element_text(vjust = 5),
          plot.title = element_text(hjust = 0.5, vjust = -1),
          legend.box="vertical", legend.margin=margin(),
          axis.title.y = ggtext::element_markdown()) +
    guides(color=guide_legend(reverse=T))
}

plt_polished_ram_mb_racon      <- plt_polished_split_ram_mb("racon", "Racon (BWA)")
plt_polished_ram_mb_pilon      <- plt_polished_split_ram_mb("pilon", "Pilon")
plt_polished_ram_mb_polypolish <- plt_polished_split_ram_mb("polypolish", "Polypolish")
plt_polished_ram_mb_racon_mp2  <- plt_polished_split_ram_mb("racon-mp2", "Racon (Minimap2)")

ggarrange(plt_polished_ram_mb_racon + rremove("xlab"),
          plt_polished_ram_mb_racon_mp2 + rremove("xlab") + rremove("ylab"),
          plt_polished_ram_mb_pilon,
          plt_polished_ram_mb_polypolish  + rremove("ylab"),
          ncol=2,nrow=2, labels=c("A", "B", "C", "D"),
          common.legend=TRUE, legend="bottom", hjust = c(-3,-0.5,-3,-0.5),
          widths = c(1.05,1), heights = c(1,1.02))
} -> plt_polished_ram

save_plot(
  plt_polished_ram,
  "fig04_simulated_polished_ram",
  10,7
)

#
# Simulated Polished Resources Time
#
# Figure 5
{
plt_polished_split_time_mb <- function(program_name,title){
  flye_program <- paste0("flye_", program_name)
  
  flye_polished_split %>% filter(program %in% c("flye", !!flye_program, !!program_name)) %>% 
    mutate(program=factor(program, levels = c(!!flye_program, !!program_name, "flye")),
           depth_il=factor(depth_il), depth_np=as.numeric(depth_np)) %>%
    arrange(program) %>%
    ggplot(aes(x=depth_il, y=mean_time_mb, color=program, alpha=depth_np)) +
    geom_line(aes(group=interaction(program, depth_np))) +
    geom_point(size=2) +
    geom_linerange(aes(ymin=(mean_time_mb - CI_time_mb), ymax=(mean_time_mb + CI_time_mb)),
                   linewidth=0.5, show.legend = F) +
    scale_alpha_continuous(name="Depth long-reads (x)", range=c(0.25,1), breaks=c(10,30,50,100)) +
    #scale_y_continuous(labels=scales::label_number(scale=1/60), breaks=c(0,30,60,90,120)) +
    scale_color_manual(values=c("cyan4","firebrick3","black"), label=c("Flye + Polisher", "Polisher", "Flye"), name="Program") +
    scale_y_log10(labels = label_number(), breaks = c(2,5,10,25,50,100,150), limits = c(1,200) ) +
    labs(x="Depth short-reads (x)", title=title,
         y="Time (s/Mbp)<p style='font-size:10pt'>log10 scaled</p>") +
    theme_minimal(base_size=14)+
    theme(panel.grid.minor=element_blank(),
          panel.grid.major.x=element_blank(),
          axis.text.x = element_text(vjust = 5),
          plot.title = element_text(hjust = 0.5, vjust = -1),
          legend.box="vertical", legend.margin=margin(),
          axis.title.y = ggtext::element_markdown()) +
    guides(color=guide_legend(reverse=T))
}

plt_polished_time_mb_racon      <- plt_polished_split_time_mb("racon", "Racon (BWA)")
plt_polished_time_mb_pilon      <- plt_polished_split_time_mb("pilon", "Pilon")
plt_polished_time_mb_polypolish <- plt_polished_split_time_mb("polypolish", "Polypolish")
plt_polished_time_mb_racon_mp2  <- plt_polished_split_time_mb("racon-mp2", "Racon (Minimap2)")

ggarrange(plt_polished_time_mb_racon + rremove("xlab"),
          plt_polished_time_mb_racon_mp2+ rremove("xlab") + rremove("ylab"),
          plt_polished_time_mb_pilon,
          plt_polished_time_mb_polypolish+ rremove("ylab"),
          ncol=2,nrow=2, labels=c("A", "B", "C", "D"),
          common.legend=TRUE, legend="bottom", hjust = c(-3,-0.5,-3,-0.5),
          widths = c(1.05,1), heights = c(1,1.02))
} -> plt_polished_time

save_plot(
  plt_polished_time,
  "fig05_simulated_polished_time",
  10,7
)

# ----------------------------------------------------------------------------
# Figures for Hybrid Assemblies
# ----------------------------------------------------------------------------

stats_hybrid_plot <- simulated_stats_multi %>%
  filter(strategy %in% c("hybrid", "polished"),
         program %notin% c("pilon", "racon", "polypolish","racon-mp2"),
         program %notin% c("flye_pilon", "flye_racon","flye_racon-mp2"),
         depth_np <= 100 ,depth_il <= 100) %>% 
  mutate(program=factor(program, levels=c("abyss_hybrid", "spades_hybrid", "masurca", "flye_polypolish")))

#
# Simulated Hybrid Metric
#
# Figure 6
{
plt_n50_hybrid_xil <- stats_hybrid_plot %>% #mutate(CI_n50 = (CI_n50 * 1e-6), mean_n50 = (mean_n50 * 1e-6)) %>% 
  arrange(CI_n50, program) %>%
  mutate(depth_np=as.numeric(depth_np), depth_il=factor(depth_il)) %>%
  ggplot(aes(x=depth_il, y=mean_n50, color=program, alpha=depth_np)) +
  geom_line(aes(group=interaction(program, depth_np))) +
  geom_point(size=2) +
  geom_linerange(aes(ymin=(mean_n50 - CI_n50), ymax=(mean_n50 + CI_n50)),
                 linewidth=0.5, show.legend = F) +
  scale_y_continuous(labels=scales::label_number(scale=1e-6)) +
  scale_alpha_continuous(name="Depth long-reads (x)", range=c(0.25,1), breaks=c(10,30,50,100)) +
  scale_color_manual(values=colors_programs, label=program_labels, name="Program") +
  labs(x="Depth short-reads (x)", y="N50 (Mbp)") +
  # scale_y_log10(labels = label_number(), breaks = c(0.01, 0.1, 1, 2, 4)) +
  theme_minimal(base_size=14)+
  theme(panel.grid.minor=element_blank(),
        panel.grid.major.x=element_blank(),
        legend.box="vertical", legend.margin=margin())

plt_genfrac_hybrid_xil <- stats_hybrid_plot %>%
  arrange(CI_genfrac, program) %>%
  mutate(depth_np=as.numeric(depth_np), depth_il=factor(depth_il)) %>%
  ggplot(aes(x=depth_il, y=mean_genfrac, color=program, alpha=depth_np)) +
  geom_line(aes(group=interaction(program, depth_np))) +
  geom_point(size=2) +
  geom_linerange(aes(ymin=(mean_genfrac - CI_genfrac), ymax=(mean_genfrac + CI_genfrac)),
                 linewidth=0.5, show.legend = F) +
  scale_alpha_continuous(name="Depth long-reads (x)", range=c(0.25,1), breaks=c(10,30,50,100)) +
  scale_color_manual(values=colors_programs, label=program_labels, name="Program") +
  scale_y_continuous(limits=c(90,100)) +
  labs(x="Depth short-reads (x)", y="Recovered (%)") +
  theme_minimal(base_size=14) +
  theme(panel.grid.minor=element_blank(),
        panel.grid.major.x=element_blank(),
        legend.box="vertical", legend.margin=margin())

plt_busco_diff_hybrid_xil <- stats_hybrid_plot %>%
  arrange(CI_busco_diff, program) %>%
  mutate(depth_np=as.numeric(depth_np), depth_il=factor(depth_il)) %>%
  ggplot(aes(x=depth_il, y=mean_busco_diff, color=program, alpha=depth_np,)) +
  geom_hline(yintercept=0, linetype="dashed", color="grey50") +
  geom_line(aes(group=interaction(program,depth_np))) +
  geom_point(size=2) +
  geom_linerange(aes(ymin=(mean_busco_diff - CI_busco_diff), ymax=(mean_busco_diff + CI_busco_diff)),
                 linewidth=0.5, show.legend = F) +
  scale_alpha_continuous(name="Depth long-reads (x)", range=c(0.25,1), breaks=c(10,30,50,100)) +
  scale_color_manual(values=colors_programs, label=program_labels, name="Program") +
  # scale_y_log10(labels = label_number()) +
  scale_y_continuous(trans = "pseudo_log", breaks = c(-1, -2, -5, -15, -40)) +
  labs(x="Depth short-reads (x)",
       y="BUSCO difference<p style='font-size:10pt'>pseudo_log scaled</p>") +
  theme_minimal(base_size=14) +
  theme(panel.grid.minor=element_blank(),
        panel.grid.major.x=element_blank(),
        legend.box="vertical", legend.margin=margin(),
        axis.title.y = ggtext::element_markdown())

plt_len_err_hybrid_xil <- stats_hybrid_plot %>%
  arrange(CI_len_err, program) %>%
  mutate(depth_np=as.numeric(depth_np), depth_il=factor(depth_il)) %>%
  ggplot(aes(x=depth_il, y=mean_len_err, color=program, alpha=depth_np)) +
  geom_hline(yintercept=0, linetype="dashed", color="grey50") +
  geom_line(aes(group=interaction(program, depth_np))) +
  geom_point(size=2) +
  geom_linerange(aes(ymin=(mean_len_err - CI_len_err), ymax=(mean_len_err + CI_len_err)),
                 linewidth=0.5, show.legend = F) +
  scale_alpha_continuous(name="Depth long-reads (x)", range=c(0.25,1), breaks=c(10,30,50,100)) +
  scale_color_manual(values=colors_programs, label=program_labels, name="Program") +
  labs(x="Depth short-reads (x)", y="Length error (%)") +
  theme_minimal(base_size=14) +
  theme(panel.grid.minor=element_blank(),
        panel.grid.major.x=element_blank(),
        legend.box="vertical", legend.margin=margin())

plt_mism_hybrid_xil <- stats_hybrid_plot %>%
  arrange(CI_mism, program) %>%
  mutate(depth_np=as.numeric(depth_np), depth_il=factor(depth_il)) %>%
  ggplot(aes(x=depth_il, y=mean_mism, color=program, alpha=depth_np)) +
  geom_line(aes(group=interaction(program, depth_np))) +
  geom_point(size=2) +
  geom_linerange(aes(ymin=(mean_mism - CI_mism), ymax=(mean_mism + CI_mism)),
                 linewidth=0.5, show.legend = F) +
  scale_alpha_continuous(name="Depth long-reads (x)", range=c(0.25,1), breaks=c(10,30,50,100)) +
  scale_color_manual(values=colors_programs, label=program_labels, name="Program") +
  scale_y_log10(labels = label_number()) +
  labs(x="Depth short-reads (x)",
       y="Mismatches (/100kb)<p style='font-size:10pt'>log10 scaled</p>") +
  theme_minimal(base_size=14) +
  theme(panel.grid.minor=element_blank(),
        panel.grid.major.x=element_blank(),
        legend.box="vertical", legend.margin=margin(),
        axis.title.y = ggtext::element_markdown())

plt_indel_hybrid_xil <- stats_hybrid_plot %>%
  arrange(CI_indel, program) %>%
  mutate(depth_np=as.numeric(depth_np), depth_il=factor(depth_il)) %>%
  ggplot(aes(x=depth_il, y=mean_indel, color=program, alpha=depth_np)) +
  geom_line(aes(group=interaction(program, depth_np))) +
  geom_point(size=2) +
  geom_linerange(aes(ymin=(mean_indel - CI_indel), ymax=(mean_indel + CI_indel)),
                 linewidth=0.5, show.legend = F) +
  scale_alpha_continuous(name="Depth long-reads (x)", range=c(0.25,1), breaks=c(10,30,50,100)) +
  scale_color_manual(values=colors_programs, label=program_labels, name="Program") +
  scale_y_log10(labels = label_number()) +
  labs(x="Depth short-reads (x)",
       y="Indels (/100kb)<p style='font-size:10pt'>log10 scaled</p>") +
  theme_minimal(base_size=14) +
  theme(panel.grid.minor=element_blank(),
        panel.grid.major.x=element_blank(),
        legend.box="vertical", legend.margin=margin(),
        axis.title.y = ggtext::element_markdown())

ggarrange(plt_n50_hybrid_xil + rremove("xlab"),
          plt_genfrac_hybrid_xil + rremove("xlab"),
          plt_busco_diff_hybrid_xil + rremove("xlab"),
          plt_len_err_hybrid_xil + rremove("xlab"),
          plt_mism_hybrid_xil,
          plt_indel_hybrid_xil,
          ncol=2, nrow=3, labels=c("A", "B", "C", "D", "E", "F"),
          common.legend=TRUE, legend="bottom",
          heights = c(1,1,1.05))
} -> plt_hybrid_metrics

save_plot(
  plt_hybrid_metrics,
  "fig06_simulated_hybrid_metrics",
  10,10
)

#
# Simulated Hybrid Resource
#
# Figure 7
{
plt_time_mp_hybrid_xil <- stats_hybrid_plot %>%
    # filter_out(program == "abyss_hybrid") %>% 
    mutate(CI_time_mb = if_else(program == "abyss_hybrid", 0, CI_time_mb)) %>% 
  arrange(CI_time_mb, program) %>%
  mutate(depth_np=as.numeric(depth_np), depth_il=factor(depth_il)) %>%
  filter(depth_np %in% c(10,30,50,100)) %>%
  ggplot(aes(x=depth_il, y=mean_time_mb, color=program, alpha=depth_np)) +
  geom_line(aes(group=interaction(program, depth_np))) +
  geom_point(size=2) +
  geom_linerange(aes(ymin=(mean_time_mb - CI_time_mb), ymax=(mean_time_mb + CI_time_mb)),
                 linewidth=0.5, show.legend = F) +
  scale_alpha_continuous(name="Depth long-reads (x)", range=c(0.25,1), breaks=c(10,30,50,100)) +
  scale_color_manual(values=colors_programs, label=program_labels, name="Program") +
  scale_y_log10(breaks = c(20,35,60,120,200,300)) +
  labs(x="Depth short-reads (x)",
       y="Time (s/Mbp)<p style='font-size:10pt'>log10 scaled</p>") +
  theme_minimal(base_size=14) +
  theme(panel.grid.minor=element_blank(),
        panel.grid.major.x=element_blank(),
        legend.box="vertical", legend.margin=margin(),
        axis.title.y = ggtext::element_markdown()) +
    annotate("text", x = 8.5, y = 320, label = "*ABySS CI are removed",
             colour = colors_programs["abyss_hybrid"], size = 3)

plt_ram_mb_hybrid_xil <- stats_hybrid_plot %>%
  arrange(CI_ram_mb, program) %>%
  mutate(depth_np=as.numeric(depth_np), depth_il=factor(depth_il)) %>%
  filter(depth_np %in% c(10,30,50,100)) %>%
  ggplot(aes(x=depth_il, y=mean_ram_mb, color=program, alpha=depth_np)) +
  geom_line(aes(group=interaction(program, depth_np))) +
  geom_point(size=2) +
  geom_linerange(aes(ymin=(mean_ram_mb - CI_ram_mb), ymax=(mean_ram_mb + CI_ram_mb)),
                 linewidth=0.5, show.legend = F) +
  scale_alpha_continuous(name="Depth long-reads (x)", range=c(0.25,1), breaks=c(10,30,50,100)) +
  scale_color_manual(values=colors_programs, limits = unique(stats_hybrid_plot$program), label=program_labels, name="Program") +
  labs(x="Depth short-reads (x)", y="Max RAM (GB/Mbp)") +
  theme_minimal(base_size=14)+
  theme(panel.grid.minor=element_blank(),
        panel.grid.major.x=element_blank(),
        legend.box="vertical", legend.margin=margin()) + guides(color=guide_legend(nrow=1, reverse=T))

ggarrange(plt_ram_mb_hybrid_xil,
          plt_time_mp_hybrid_xil,
          ncol=2, labels=c("A", "B") ,
          common.legend=TRUE, legend="bottom")
} -> plt_hybrid_resources

save_plot(
  plt_hybrid_resources,
  "fig07_simulated_hybrid_resources",
  10,5
)

# ----------------------------------------------------------------------------
# Supplementary
# ----------------------------------------------------------------------------

#
# Simulated - Reference Genomes Characteristics
#
# Figure S1
{
  colors_ref <- data.frame (
    light  = c("lightblue", "darkolivegreen2", "coral" , "mediumorchid2", "turquoise"),
    dark   = c("darkblue" , "darkolivegreen4", "coral4", "magenta4"     , "turquoise4")
  )
  
  plot_metric <- function(metric, y_label, colors) {
    box_fill <- colors$light
    jitter_color <- colors$dark
    
    reference_metrics %>% 
      ggplot(aes(x="", y=!!sym(metric))) +
      # geom_violin(width=1, fill=box_fill) +
      geom_violinhalf(color=jitter_color, fill=box_fill, alpha=0.7) +
      geom_boxplot(width=0.3, fill=jitter_color, alpha=0.9, outliers = T) +
      # geom_jitter(width=0.1, alpha=0.5, color=jitter_color) +
      labs(x="", y=y_label) +
      theme_minimal(base_size=10) +
      theme(panel.grid.minor=element_blank(),
            panel.grid.major.x=element_blank())
  }
  
  order_metrics <- c("cds"= "CDS",
                     "intron"="Introns",
                     "gene"= "Genes")
  reference_metrics %>% 
    select(sample, size_mbp, gene, cds, intron, m_rna) %>%
    pivot_longer(cols=c(gene, cds, intron),
                 names_to="feature", values_to="count") %>%
    group_by(sample) %>% mutate(feature=factor(feature, levels=c("cds", "intron", "gene"))) %>%
    ggplot(aes(x=size_mbp, y=count, fill = feature, color = feature)) +
    geom_point(alpha=0.7) +
    stat_poly_line() +
    stat_poly_eq(label.y = c(0.9,0.85,0.8), label.x = 0.05) +
    labs(x="", y="# Features (x1000)") +
    scale_y_continuous(labels=scales::label_number(scale=1/1000),
                       breaks = seq(0, 150e3, by = 25e3)) +
    scale_color_manual(values=colors_ref$dark, name="Feature", labels=order_metrics) +
    scale_fill_manual(values=colors_ref$light, name="Feature", labels=order_metrics) +
    theme_minimal(base_size=10) +
    theme(panel.grid.minor=element_blank(),
          panel.grid.major.x=element_blank(),
          axis.text.x = element_blank(),
          plot.margin=unit(c(0,0.25,-0.5,0.25), "cm"),
          legend.position="top",
          legend.box="vertical", legend.margin=margin()) -> plot_reference_features
  
  reference_metrics %>% 
    select(sample, size_mbp) %>% 
    ggplot(aes(x=size_mbp, y=" ")) +
    geom_boxplot(fill="darkgoldenrod1", color="darkorange4", width=0.75, outlier.color = "darkorange4") +
    labs(x="Size (Mbp)", y=" ") +
    theme_minimal(base_size=10) +
    theme(panel.grid.minor=element_blank(),
          panel.grid.major=element_blank(),
          plot.margin=unit(c(0,0.25,0.25,0.70), "cm")) -> plot_reference_size
  
  ggarrange(
    ggarrange(plot_reference_features, plot_reference_size,
              ncol=1, nrow=2, heights = c(5,1)),
    plot_metric("contigs", "# Contigs", colors_ref[4,]),
    plot_metric("gc_content", "GC Content (%)", colors_ref[5,]),
    ncol=3, nrow=1, labels=c("A", "B", "C"), widths = c(5,1,1))
  
} -> plot_reference_characteristics

save_plot(
  plot_reference_characteristics,
  "figS01_simulated_reference_characteristics",
  8,4
)

#
# Simulated Sup Polished Metrics
#
# Figure S2
{
plt_n50_polished <- flye_polished %>%
  mutate(depth_il=as.numeric(depth_il), depth_np=factor(depth_np)) %>%
  ggplot(aes(x=depth_np, y=mean_n50, color=program, alpha=depth_il)) +
  geom_line(aes(group=interaction(program, depth_il)), position = "identity") +
  geom_point(size=2, position = "identity") +
  geom_linerange(aes(ymin=(mean_n50 - CI_n50), ymax=(mean_n50 + CI_n50)),
                 linewidth=0.5, show.legend = F, position = "identity") +
  scale_y_continuous(labels=scales::label_number(scale=1e-6)) +
  scale_alpha_continuous(name="Depth short-reads (x)", range=c(0.25,1), breaks=c(10,30,50,100)) +
  scale_color_manual(values=colors_tmp, label=program_labels, name="Program") +
  labs(x="Depth long-reads (x)", y="N50 (Mbp)") +
  theme_minimal(base_size=14)+
  theme(panel.grid.minor=element_blank(),
        panel.grid.major.x=element_blank(),
        legend.box="vertical", legend.margin=margin())

plt_len_err_polished <- flye_polished %>%
  mutate(depth_il=as.numeric(depth_il), depth_np=factor(depth_np)) %>%
  ggplot(aes(x=depth_np, y=mean_len_err, color=program, alpha=depth_il)) +
  geom_hline(yintercept=0, linetype="dashed", color="grey50") +
  geom_line(aes(group=interaction(program, depth_il)), position = "identity") +
  geom_point(size=2, position = "identity") +
  geom_linerange(aes(ymin=(mean_len_err - CI_len_err), ymax=(mean_len_err + CI_len_err)),
                 linewidth=0.5, show.legend = F, position = "identity") +
  scale_alpha_continuous(name="Depth short-reads (x)", range=c(0.25,1), breaks=c(10,30,50,100)) +
  scale_color_manual(values=colors_tmp, label=program_labels, name="Program") +
  labs(x="Depth long-reads (x)", y="Length error (%)") +
  theme_minimal(base_size=14) +
  theme(panel.grid.minor=element_blank(),
        panel.grid.major.x=element_blank(),
        legend.box="vertical", legend.margin=margin())

plt_busco_diff_polished <- flye_polished %>%
  mutate(depth_il=as.numeric(depth_il), depth_np=factor(depth_np)) %>%
  ggplot(aes(x=depth_np, y=mean_busco_diff, color=program, alpha=depth_il)) +
  geom_hline(yintercept=0, linetype="dashed", color="grey50") +
  geom_line(aes(group=interaction(program, depth_il)), position = "identity") +
  geom_point(size=2, position = "identity") +
  geom_linerange(aes(ymin=(mean_busco_diff - CI_busco_diff), ymax=(mean_busco_diff + CI_busco_diff)),
                 linewidth=0.5, show.legend = F, position = "identity") +
  scale_alpha_continuous(name="Depth short-reads (x)", range=c(0.25,1), breaks=c(10,30,50,100)) +
  scale_color_manual(values=colors_tmp, label=program_labels, name="Program") +
  labs(x="Depth long-reads (x)", y="BUSCO difference") +
  theme_minimal(base_size=14) +
  theme(panel.grid.minor=element_blank(),
        panel.grid.major.x=element_blank(),
        legend.box="vertical", legend.margin=margin())

plt_genfrac_polished <- flye_polished %>%
  mutate(depth_il=as.numeric(depth_il), depth_np=factor(depth_np)) %>%
  ggplot(aes(x=depth_np, y=mean_genfrac, color=program, alpha=depth_il)) +
  geom_line(aes(group=interaction(program, depth_il)), position = "identity") +
  geom_point(size=2, position = "identity") +
  geom_linerange(aes(ymin=(mean_genfrac - CI_genfrac), ymax=(mean_genfrac + CI_genfrac)),
                 linewidth=0.5, show.legend = F, position = "identity") +
  scale_alpha_continuous(name="Depth short-reads (x)", range=c(0.25,1), breaks=c(10,30,50,100)) +
  scale_color_manual(values=colors_tmp, label=program_labels, name="Program") +
  labs(x="Depth long-reads (x)", y="Recovered (%)") +
  theme_minimal(base_size=14) +
  theme(panel.grid.minor=element_blank(),
        panel.grid.major.x=element_blank(),
        legend.box="vertical", legend.margin=margin())

ggarrange(plt_n50_polished + guides(color=guide_legend(nrow=1, reverse=T)),
          plt_len_err_polished,
          plt_busco_diff_polished, plt_genfrac_polished,
          ncol=2, nrow=2, labels=c("A", "B", "C", "D"),
          common.legend=TRUE, legend="bottom") 
} -> plt_polished_metrics_sup

save_plot(
  plt_polished_metrics_sup,
  "figS02_simulated_polished_metrics_sup",
  10,10
)

#
# Simulated Sup Hybrid Resources RAM
#
# Figure S3
{
plt_time_mp_hybrid_xil_full <- stats_hybrid_plot %>%
  arrange(CI_time_mb, program) %>%
  mutate(depth_np=as.numeric(depth_np), depth_il=factor(depth_il)) %>%
  ggplot(aes(x=depth_il, y=mean_time_mb, color=program, alpha=depth_np)) +
  geom_line(aes(group=interaction(program, depth_np))) +
  geom_point(size=2) +
  geom_linerange(aes(ymin=(mean_time_mb - CI_time_mb), ymax=(mean_time_mb + CI_time_mb)),
                 linewidth=0.5, show.legend = F) +
  scale_alpha_continuous(name="Depth long-reads (x)", range=c(0.25,1), breaks=c(10,30,50,100)) +
  scale_color_manual(values=colors_programs, label=program_labels, name="Program") +
  labs(x="Depth short-reads (x)", y="Time (s/Mbp)") +
  theme_minimal(base_size=14)+
  facet_wrap(~depth_np, scales = "free",
               labeller = labeller(depth_np = ~ paste0("Long-reads ", .x, "x")),
               nrow = 4,  ncol = 3) +
  theme(panel.grid.minor=element_blank(),
        panel.grid.major.x=element_blank(),
        legend.box="horizontal", legend.margin=margin(),
        legend.position = c(0.9, 0.08),
        legend.justification = c(1, 0)) +
    guides(alpha="none",
           color=guide_legend(ncol=1,byrow=F))

egg::tag_facet(plt_time_mp_hybrid_xil_full,
                                      x = -Inf,
                                      open = "", close = "",
                                      tag_pool = LETTERS) + theme(strip.text = element_text())
} -> plt_hybrid_time_sup
save_plot(
  plt_hybrid_time_sup,
  "figS03_simulated_hybrid_time_sup",
  15,15
)

#
# Simulated Sup Hybrid Resources Time
#
# Figure S4
{
plt_ram_mb_hybrid_xil_full <- stats_hybrid_plot %>%
  arrange(CI_ram_mb, program) %>%
  mutate(depth_np=as.numeric(depth_np), depth_il=factor(depth_il)) %>%
  ggplot(aes(x=depth_il, y=mean_ram_mb, color=program, alpha=depth_np)) +
  geom_line(aes(group=interaction(program, depth_np))) +
  geom_point(size=2) +
  geom_linerange(aes(ymin=(mean_ram_mb - CI_ram_mb), ymax=(mean_ram_mb + CI_ram_mb)),
                 linewidth=0.5, show.legend = F) +
  scale_alpha_continuous(name="Depth long-reads (x)", range=c(0.25,1), breaks=c(10,30,50,100)) +
  scale_color_manual(values=colors_programs, label=program_labels, name="Program") +
  labs(x="Depth short-reads (x)", y="Max RAM (GB/Mbp)") +
  theme_minimal(base_size=14) +
  facet_wrap(~depth_np, 
             labeller = labeller(depth_np = ~ paste0("Long-reads ", .x, "x")),
             nrow = 4,  ncol = 3) +
    theme(panel.grid.minor=element_blank(),
          panel.grid.major.x=element_blank(),
          legend.box="horizontal", legend.margin=margin(),
          legend.position = c(0.9, 0.08),
          legend.justification = c(1, 0)) +
    guides(alpha="none",
           color=guide_legend(ncol=1,byrow=F))

egg::tag_facet(plt_ram_mb_hybrid_xil_full,
                                     x = -Inf, y = 0.90,
                                     open = "", close = "",
                                     tag_pool = LETTERS) + theme(strip.text = element_text())
} -> plt_hybrid_ram_sup
save_plot(
  plt_hybrid_ram_sup,
  "figS04_simulated_hybrid_ram_sup",
  15,15
)

# ----------------------------------------------------------------------------
# Figures for Model 
# ----------------------------------------------------------------------------

filt_responses <- c( "n50",
                     "busco_complete_diff",
                     "genome_fraction_percent",
                     "length_err_pct",
                     "mism_100kbp",
                     "indel_100kbp",
                     "ram_per_mb",
                     "time_per_mb"
)

response_labels <- c(
  "n50"                     = "N50",
  "busco_complete_diff"     = "BUSCO difference",
  "genome_fraction_percent" = "Genome recovery",
  "length_err_pct"          = "Length error",
  "ram_per_mb"              = "Max RAM (/Mbp)",
  "time_per_mb"             = "Run Time(/Mbp)",
  "mism_100kbp"             = "Mismatches (/100kb)",
  "indel_100kbp"            = "Indels (/100kb)"
)


forest_full_slopes <- function(fit_results) {

  depth_df <- map_dfr(fit_results, \(f) f$depth_slopes) %>%
                rename(covariate=slope_type)
  
  df <- bind_rows(depth_df,
    map_dfr(fit_results, \(f) f$genome_slopes)) %>% 
    mutate(signif = case_when(
      is.na(p_value)  ~ "",
      p_value < 0.001 ~ "< 0.001 ",
      p_value < 0.05  ~ "< 0.05",
      p_value >= 0.05  ~ "≥ 0.05",
      TRUE      ~ ""),
      covariate = str_replace_all(covariate,
                             c("program" = "Program",
                               "eff" = "Depth",
                               "np" = "Depth LR",
                               "il" = "Depth SR",
                               "size_mbp_s" = "Genome Size",
                               "gene_density_s" = "Gene Density",
                               "gc_content_s" = "GC Content",
                               "cds_per_gene_s" = "CDS/gene",
                               ":" = " x ")),
      program = fct_inorder(program),
      response = str_replace_all(response, response_labels),
      response = fct_relevel(response, rev(response_labels))
    )
  
  tmp_ploter <- function(df, title) {
    ggplot(df, aes(x = slope, y = response , color = program, shape = signif)) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "grey20") +
      geom_pointrange(
        aes(xmin = lower.CL, xmax = upper.CL),
        size  = 0.4,
        position = position_dodge(width = 0.5)
      ) +
      scale_y_discrete(labels = response_labels) +
      scale_color_manual(values=colors_programs, label=program_labels, name="Program") +
      labs(
        x     = paste0("Slope (", title, ")"),
        y     = NULL,
        color = "Response",
        shape = expression(italic(P)~"val")
      ) +
      theme_minimal(base_size = 12) +
      theme(
        panel.grid.minor   = element_blank(),
        panel.grid.major.y = element_line(color = "gray92", linewidth = 12),
        panel.grid.major.x = element_line(color = "gray", linewidth = 0.75),
        axis.text.y.right = element_blank(),
        axis.line.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        legend.position = "bottom"
      )
  }
  
  plt_size <- tmp_ploter(df %>% filter(covariate == "Genome Size"), "Genome Size")
  plt_gene_density <- tmp_ploter(df %>% filter(covariate == "Gene Density"), "Gene Density")
  plt_gc <- tmp_ploter(df %>% filter(covariate == "GC Content"), "GC Content")
  plt_cds <- tmp_ploter(df %>% filter(covariate == "CDS/gene"), "CDS/gene")

  if ("Depth LR" %in% unique(df$covariate) && "Depth SR" %in% unique(df$covariate)) {
    plt_lr <- tmp_ploter(df %>% filter(covariate == "Depth LR"), "Depth Long-reads") +
                  rremove("legend")
    
    plt_sr <- tmp_ploter(df %>% filter(covariate == "Depth SR"), "Depth Short-reads") +
                  ggbreak::scale_x_break(c(0.03, 0.24), scales = "fixed") +
                  rremove("y.text") + rremove("legend")

    plt_top <- ggarrange(plt_lr, print(plt_sr),
                          ncol=2, nrow=1, labels=c("A", "B"),
                          widths = c(1.2,1)) + rremove("legend")

    plt_bottom <- ggarrange(plt_size, plt_gene_density + rremove("y.text"),
                            plt_gc + rremove("y.text"), plt_cds + rremove("y.text"),
                            ncol=4, nrow=1, labels=c("C", "D", "E", "F"),
                            common.legend = TRUE, legend="bottom",
                            widths = c(1.5,1,1,1))

    ggarrange(plt_top, plt_bottom,
              ncol=1, nrow=2, labels=c("", ""),
              common.legend = TRUE, legend="bottom",
              heights = c(1.2,1))
  } else {
    plt_top <- tmp_ploter(df %>% filter(covariate == "Depth"), "Depth") +
    ggbreak::scale_x_break(c(0.1, 0.2), scales = "fixed") +
    guides(shape=guide_legend(nrow=3,byrow=F),
            color=guide_legend(nrow=2,byrow=F)) + rremove("legend")

    plt_bottom <- ggarrange(plt_size, plt_gene_density + rremove("y.text"),
                            plt_gc + rremove("y.text"), plt_cds + rremove("y.text"),
                            ncol=4, nrow=1, labels=c("B", "C", "D", "E"),
                            common.legend = TRUE, legend="bottom",
                            widths = c(1.5,1,1,1))

    ggarrange(print(plt_top), plt_bottom,
              ncol=1, nrow=2, labels=c("A", ""),
              common.legend = TRUE, legend="bottom",
              heights = c(1.2,1))
  }
}

plt_forest_full_slopes_single <- forest_full_slopes(single_fits[filt_responses])
save_plot(
  plt_forest_full_slopes_single,
  "figS05_simulated_forest_full_slopes_single",
  15,10
)

plt_forest_full_slopes_multi <- forest_full_slopes(multi_fits[filt_responses])
save_plot(
  plt_forest_full_slopes_multi,
  "figS06_simulated_forest_full_slopes_multi",
  15,10
)

# forest_all_slopes <- function(fit_results) {

#   df <- map_dfr(fit_results, function(f) {
#     f$depth_slopes
#   }) %>% 
#     mutate(signif = case_when(
#       is.na(p_value)  ~ "",
#       p_value < 0.001 ~ "< 0.001 ",
#       p_value < 0.05  ~ "< 0.05",
#       p_value >= 0.05  ~ "≥ 0.05",
#       TRUE      ~ ""),
#       # slope = case_when(
#       #   response == "mism_100kbp" ~ slope * -1,
#       #   response == "indel_100kbp" ~ slope * -1,
#       #   TRUE ~ slope),
#       # lower.CL = case_when(
#       #   response == "mism_100kbp" ~ lower.CL * -1,
#       #   response == "indel_100kbp" ~ lower.CL * -1,
#       #   TRUE ~ lower.CL),
#       # upper.CL = case_when(
#       #   response == "mism_100kbp" ~ upper.CL * -1,
#       #   response == "indel_100kbp" ~ upper.CL * -1,
#       #   TRUE ~ upper.CL),
#       program = fct_inorder(program),
#       response = str_replace_all(response, response_labels),
#       response = fct_relevel(response, rev(response_labels))
#     )
  
#   tmp_ploter <- function(df, title) {
#   ggplot(df, aes(x = slope, y = response , color = program, shape = signif)) +
#     geom_vline(xintercept = 0, linetype = "dashed", color = "grey20") +
#     geom_pointrange(
#       aes(xmin = lower.CL, xmax = upper.CL),
#       size  = 0.4,
#       position = position_dodge(width = 0.5)
#     ) +
#     scale_y_discrete(labels = response_labels) +
#     scale_color_manual(values=colors_programs, label=program_labels, name="Program") +
#     labs(
#       x     = "Slope (change per unit  of depth)",
#       y     = NULL,
#       color = "Response",
#       shape = expression(italic(P)~"val"),
#       title = title
#     ) +
#     theme_minimal(base_size = 12) +
#     theme(
#       panel.grid.minor   = element_blank(),
#       panel.grid.major.y = element_line(color = "gray92", linewidth = 18),
#       panel.grid.major.x = element_line(color = "gray", linewidth = 0.75),
#       axis.text.y.right = element_blank(),
#       axis.line.y.right = element_blank(),
#       axis.ticks.y.right = element_blank(),
#       plot.title = element_text(hjust = 0.5, vjust = -1),
#       legend.position = "top"
#     ) +
#     guides(shape=guide_legend(nrow=3,byrow=F),
#            color=guide_legend(nrow=3,byrow=F)
#     )}

#   if (length(unique(df$slope_type)) == 2){
#     plt_lr <- tmp_ploter(df %>% filter(slope_type == "np"), "Depth Long-reads") +
#     guides(shape=guide_legend(nrow=3,byrow=F),
#             color=guide_legend(nrow=2,byrow=F))
    
#     plt_sr <- tmp_ploter(df %>% filter(slope_type == "il"), "Depth Short-reads") +
#                   ggbreak::scale_x_break(c(0.03, 0.24), scales = "fixed") +
#                   rremove("y.text") + rremove("legend")

#     ggarrange(plt_lr, print(plt_sr),
#               ncol=2, nrow=1, labels=c("A", "B"),
#               common.legend = TRUE, legend="top",
#               widths = c(1.2,1))
#   } else {
#     tmp_ploter(df, "Depth") +
#     ggbreak::scale_x_break(c(0.1, 0.2), scales = "fixed") +
#     guides(shape=guide_legend(nrow=3,byrow=F),
#             color=guide_legend(nrow=2,byrow=F)) +
#     theme(legend.position = "top")
#   }
# }

# plt_forest_all_slopes_single <- forest_all_slopes(single_fits[filt_responses])
# plt_forest_all_slopes_multi <- forest_all_slopes(multi_fits[filt_responses])

# save_plot(
#   plt_forest_all_slopes_single,
#   "simulated_forest_all_slopes_single",
#   15,7
# )

# save_plot(
#   plt_forest_all_slopes_multi,
#   "simulated_forest_all_slopes_multi",
#   15,7
# )

# forest_all_fixed <- function(fit_results) {
#   df <- map_dfr(fit_results, function(f) {
#     f$fixed_effects %>%
#       filter(term != "(Intercept)") %>%
#       mutate(response = f$fit$response)
#   }) %>% 
#     mutate(term = str_replace_all(term, c("program" = "",
#                                           "depth_eff" = "Depth",
#                                           "depth_np" = "Depth LR",
#                                           "depth_il" = "Depth SR",
#                                           "size_mbp_s" = "Genome Size",
#                                           "gene_density_s" = "Gene Density",
#                                           "gc_content_s" = "GC Content",
#                                           "cds_per_gene_s" = "CDS/gene",
#                                           ":" = " x ")),
#            term = str_replace_all(term, coll(program_labels)),
#            response = str_replace_all(response, response_labels),
#            response = fct_relevel(response, response_labels)
#     ) %>% group_by(response) %>% 
#     arrange(estimate, .by_group = T) %>% 
#     mutate()
  
#   ggplot(df, aes(x = estimate, y = term, color = response)) +
#     geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
#     geom_pointrange(
#       aes(xmin = conf.low, xmax = conf.high),
#       size = 0.3,
#       position = position_dodge(width = 0.5)
#     ) +
#     scale_x_continuous(trans = "pseudo_log") +
#     scale_color_manual(values=color_palette,,  name="Response") +
#     labs(
#       x     = "Estimate<p style='font-size:10pt'>log10 scaled</p>",
#       y     = NULL
#     ) +
#     theme_minimal(base_size = 12) +
#     theme(
#       panel.grid.minor = element_blank(),
#       axis.title.x = ggtext::element_markdown(),
#       legend.position  = "top"
#     )
# }

# plt_forest_all_fixed_single <- forest_all_fixed(single_fits[filt_responses])
# save_plot(
#   plt_forest_all_fixed_single,
#   "simulated_forest_all_fixed_single",
#   15,10
# )

# plt_forest_all_fixed_multi <- forest_all_fixed(multi_fits[filt_responses])
# save_plot(
#   plt_forest_all_fixed_multi,
#   "simulated_forest_all_fixed_multi",
#   15,10
# ) 

# forest_all_global <- function(fit_results) {
#   df <- map_dfr(fit_results, function(f) {
#     f$global_effects %>%
#       filter(term != "(Intercept)") %>%
#       mutate(response = f$fit$response)
#   }) %>% 
#     mutate(signif = case_when(
#       is.na(p_value)  ~ "",
#       p_value < 0.001 ~ "< 0.001 ",
#       p_value < 0.05  ~ "< 0.05",
#       p_value >= 0.05  ~ "≥ 0.05",
#       TRUE      ~ ""),
#       term = str_replace_all(term, c("program" = "Program",
#                                      "depth_eff" = "Depth",
#                                      "depth_np" = "Depth LR",
#                                      "depth_il" = "Depth SR",
#                                      "size_mbp_s" = "Genome Size",
#                                      "gene_density_s" = "Gene Density",
#                                      "gc_content_s" = "GC Content",
#                                      "cds_per_gene_s" = "CDS/gene",
#                                      ":" = " x ")),
#       response = str_replace_all(response, response_labels),
#       response = fct_relevel(response, response_labels)
#     )
  
#   ggplot(df, aes(x = `chi_sq`, y = term, color = response, shape = signif)) +
#     geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
#     geom_point(size = 2,
#                position = position_dodge(width = 0.5)) +
#     scale_x_continuous(trans = "pseudo_log", breaks = c(0,1,10, 100, 1000, 10000, 50000)) +
#     labs(
#       title = NULL,
#       x     = expression(chi^2),
#       y     = NULL,
#       shape = expression(italic(P)~"val")
#     ) +
#     theme_minimal(base_size = 12) +
#     scale_color_manual(values=color_palette,  name="Response") +
#     theme(
#       panel.grid.minor = element_blank(),
#       legend.position  = "top",
#       panel.grid.major.x = element_line(color = "gray", size = 0.75),
#       plot.title = element_text(hjust = 0.5, vjust = -1)
#     ) +
#     guides(shape=guide_legend(nrow=3,byrow=F),
#            color=guide_legend(nrow=3,byrow=F)
#            ) -> tmp_plt
  
#   if ("Depth LR" %in% df$term){
#     tmp_plt + 
#     theme(
#       panel.grid.major.y = element_line(color = "gray92", linewidth = 8),
#     ) 
#   } else {
#     tmp_plt + theme(
#       panel.grid.major.y = element_line(color = "gray92", linewidth = 12),
#     )
#   }
# }

# plt_forest_all_global_single <- forest_all_global(single_fits[filt_responses])
# save_plot(
#   plt_forest_all_global_single,
#   "simulated_forest_all_global_single",
#   15,10
# )

# plt_forest_all_global_multi <- forest_all_global(multi_fits[filt_responses])
# save_plot(
#   plt_forest_all_global_multi,
#   "simulated_forest_all_global_multi",
#   15,10
# )




# forest_all_genome_slopes <- function(fit_results) {
#   df <- map_dfr(fit_results, function(f) {
#     f$genome_slopes
#   }) %>% 
#     mutate(signif = case_when(
#       is.na(p_value)  ~ "",
#       p_value < 0.001 ~ "< 0.001 ",
#       p_value < 0.05  ~ "< 0.05",
#       p_value >= 0.05  ~ "≥ 0.05",
#       TRUE      ~ ""),
#       covariate = str_replace_all(covariate,
#                              c("program" = "Program",
#                                "depth_eff" = "Depth",
#                                "depth_np" = "Depth LR",
#                                "depth_il" = "Depth SR",
#                                "size_mbp_s" = "Genome Size",
#                                "gene_density_s" = "Gene Density",
#                                "gc_content_s" = "GC Content",
#                                "cds_per_gene_s" = "CDS/gene",
#                                ":" = " x ")),
#       program = fct_inorder(program),
#       response = str_replace_all(response, response_labels),
#       response = fct_relevel(response, rev(response_labels))
#     )
  
#   ggplot(df, aes(x = slope, y = response , color = program, shape = signif)) +
#     geom_vline(xintercept = 0, linetype = "dashed", color = "grey20") +
#     geom_pointrange(
#       aes(xmin = lower.CL, xmax = upper.CL),
#       size  = 0.4,
#       position = position_dodge(width = 0.5)
#     ) +
#     scale_y_discrete(labels = response_labels) +
#     scale_color_manual(values=colors_programs, label=program_labels, name="Program") +
#     labs(
#       x     = "Slope",
#       y     = NULL,
#       color = "Response",
#       shape = expression(italic(P)~"val")
#     ) +
#     theme_minimal(base_size = 12) +
#     theme(
#       panel.grid.minor   = element_blank(),
#       panel.grid.major.y = element_line(color = "gray92", linewidth = 18),
#       panel.grid.major.x = element_line(color = "gray", linewidth = 0.75),
#       axis.text.y.right = element_blank(),
#       axis.line.y.right = element_blank(),
#       axis.ticks.y.right = element_blank(),
#       legend.position = "top"
#     ) +
#     facet_wrap(~ covariate , scales = "free_x", ncol = 4) +
#     guides(shape=guide_legend(nrow=3,byrow=F),
#            color=guide_legend(nrow=3,byrow=F)
#     )

    
# }

# plt_forest_all_genome_slopes_single <- forest_all_genome_slopes(single_fits[filt_responses])
# save_plot(
#   plt_forest_all_genome_slopes_single,
#   "simulated_forest_all_genome_slopes_single",
#   15,10
# )

# plt_forest_all_genome_slopes_multi <- forest_all_genome_slopes(multi_fits[filt_responses])
# save_plot(
#   plt_forest_all_genome_slopes_multi,
#   "simulated_forest_all_genome_slopes_multi",
#   15,10
# )

save.image(file.path(results_dir, "03_figures.RData"))