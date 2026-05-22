# ============================================================================
# 03_figures.R
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
  library(paletteer)
  library(ggpubr)
  library(scales)
  library(ggtext)
})

color_palette <- paletteer_d("ggsci::default_locuszoom")

response_labels <- c(
  "n50"             = "N50",
  "busco_complete"  = "BUSCO completeness",
  "maxrss_gb"       = "Max RAM",
  "ram_mb"          = "Max RAM (/Mbp)",
  "time_mb"         = "Run Time(/Mbp)",
  "ram_mb_p"        = "Max RAM (/Mbp)",
  "time_mb_p"       = "Run Time(/Mbp)",
  "elapsed_sec"     = "Run Time",
  "mism_100kbp"     = "Mismatches (/100kb)",
  "indel_100kbp"    = "Indels (/100kb)"
)

# Set color palette
colors_programs <- c("flye"           = "#2CA02C",
                     "flye-ovl1000"   = "#59A14F",
                     "flye-ovl1500"   = "#1F83B4",
                     "flye-ovl2000"   = "#A52A2A",
                     "flye-ovl2500"   = "#C7519C",
                     "flye-ovl3000"   = "#8C564B",
                     "kmc_np"         = "#CD1076",
                     "kmc_il"         = "#68228B",
                     "kmc_chopper"    = "#CD1076",
                     "spades_short"   = "#EDC948",
                     "abyss_short"    = "#FF7F0E",
                     "spades_hybrid"  = "#EDC948",
                     "abyss_hybrid"   = "#FF7F0E",
                     "polypolish"     = "#00688B",
                     "flye-ovl1000_polypolish"= "#00688B"
)

program_labels <- c("flye"           ="Flye",
                    "flye-ovl1000"   ="Flye Min Overlap 1000",
                    "flye-ovl1500"   ="Flye Min Overlap 1500",
                    "flye-ovl2000"   ="Flye Min Overlap 2000",
                    "flye-ovl2500"   ="Flye Min Overlap 2500",
                    "flye-ovl3000"   ="Flye Min Overlap 3000",
                    "kmc_np"         ="Estimated Nanopore",
                    "kmc_il"         ="Estimated Short Reads",
                    "kmc_chopper"    ="Estimated Long Reads",
                    "spades_short"   ="SPAdes Short",
                    "abyss_short"    ="ABySS Short",
                    "abyss_hybrid"   ="ABySS Hybrid",
                    "spades_hybrid"  ="SPAdes Hybrid",
                    "polypolish"     ="Polypolish",
                    "flye-ovl1000_polypolish"="Flye Min Ovl 1000 + Polypolish")

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
#
# Empirical Single Sequencing
#
# Figure 8

# ----------------------------------------------------------------------------
# Figures for Single Seq Assemblies
# ----------------------------------------------------------------------------

program_order <- c("spades_short", "abyss_short", "flye-ovl1000",
                   "flye-ovl1500", "flye-ovl2000", "flye-ovl2500",
                   "flye-ovl3000")

x_labels <- c("10x", "15x", "20x", "25x", "30x", "35x", "40x", "50x", "60x", "OG")

empirical_stats_single_plt <- empirical_stats_single %>%
  filter(program %in% program_order) %>% 
  mutate(depth_eff=factor(depth_eff),
         program=factor(program, levels=program_order))

{
# N50 plot
plt_n50_single <- empirical_stats_single_plt %>%
  ggplot(aes(x=depth_eff, y=mean_n50, color=program)) +
  geom_line(aes(group=program)) +
  geom_point(size=2) +
  geom_linerange(aes(ymin=(mean_n50 - CI_n50), ymax=(mean_n50 + CI_n50)),
                 linewidth=0.5, show.legend = F) +
  scale_y_continuous(labels=scales::label_number(scale=1e-6)) +
  scale_x_discrete(labels=x_labels) +
  scale_color_manual(values=colors_programs, label=program_labels, name="Program") +
  labs(x="Depth (x)", y="N50 (Mbp)", title=element_blank()) +
  theme_minimal(base_size=10) +
  theme(panel.grid.minor=element_blank(),
        panel.grid.major.x=element_blank()) 

# Length plot
plt_length_single <- empirical_stats_single_plt %>%
  ggplot(aes(x=depth_eff, y=mean_length, color=program)) +
  geom_line(aes(group=program)) +
  geom_point(size=2) +
  geom_linerange(aes(ymin=(mean_length - CI_length), ymax=(mean_length + CI_length)),
                 linewidth=0.5, show.legend = F) +
  scale_y_continuous(labels=scales::label_number(scale=1e-6)) +
  scale_x_discrete(labels=x_labels) +
  scale_color_manual(values=colors_programs, label=program_labels, name="Program") +
  labs(x="Depth (x)", y="Assembly Size (Mbp)", title=element_blank()) +
  theme_minimal(base_size=10) +
  theme(panel.grid.minor=element_blank(),
        panel.grid.major.x=element_blank())

# BUSCO complete plot
plt_busco_single <- empirical_stats_single_plt %>%
  ggplot(aes(x=depth_eff, y=mean_busco, color=program)) +
  geom_line(aes(group=program)) +
  geom_point(size=2) +
  geom_linerange(aes(ymin=(mean_busco - CI_busco), ymax=(mean_busco + CI_busco)),
                 linewidth=0.5, show.legend = F) +
  scale_x_discrete(labels=x_labels) +
  scale_color_manual(values=colors_programs, label=program_labels, name="Program") +
  labs(x="Depth (x)", y="BUSCO %", title=element_blank()) +
  theme_minimal(base_size=10) +
  # coord_cartesian(ylim=c(-5,5)) +
  theme(panel.grid.minor=element_blank(),
        panel.grid.major.x=element_blank())

# BUSCO Copy plot
plt_busco_copy_single <- empirical_stats_single_plt %>%
  ggplot(aes(x=depth_eff, y=mean_busco_copy, color=program)) +
  geom_line(aes(group=program)) +
  geom_point(size=2) +
  geom_linerange(aes(ymin=(mean_busco_copy - CI_busco_copy), ymax=(mean_busco_copy + CI_busco_copy)),
                 linewidth=0.5, show.legend = F) +
  scale_x_discrete(labels=x_labels) +
  scale_color_manual(values=colors_programs, label=program_labels, name="Program") +
  labs(x="Depth (x)", y="BUSCO % Multiple Copy", title=element_blank()) +
  theme_minimal(base_size=10) +
  # coord_cartesian(ylim=c(-5,5)) +
  theme(panel.grid.minor=element_blank(),
        panel.grid.major.x=element_blank())

# Max RAM per Mbp plot
plt_ram_mb_single <- empirical_stats_single_plt %>%
  ggplot(aes(x=depth_eff, y=mean_ram_mb_p, color=program)) +
  geom_line(aes(group=program)) +
  geom_point(size=2) +
  geom_linerange(aes(ymin=(mean_ram_mb_p - CI_ram_mb_p), ymax=(mean_ram_mb_p + CI_ram_mb_p)),
                 linewidth=0.5, show.legend = F) +
  scale_x_discrete(labels=x_labels) +
  scale_color_manual(values=colors_programs, label=program_labels, name="Program") +
  labs(x="Depth (x)", y="Max RAM (GB/Mbp)", title=element_blank()) +
  theme_minimal(base_size=10)+
  theme(panel.grid.minor=element_blank(),
        panel.grid.major.x=element_blank()) + coord_cartesian(ylim=c(0,2))

# Time per Mbp plot
plt_time_mb_single <- empirical_stats_single_plt %>%
  ggplot(aes(x=depth_eff, y=mean_time_mb_p, color=program)) +
  geom_line(aes(group=program)) +
  geom_linerange(aes(ymin=(mean_time_mb_p - CI_time_mb_p), ymax=(mean_time_mb_p + CI_time_mb_p)),
                 linewidth=0.5, show.legend = F) +
  geom_point(size=2) +
  scale_x_discrete(labels=x_labels) +
  #scale_y_continuous(labels=scales::label_number(scale=1/60)) +
  scale_color_manual(values=colors_programs, label=program_labels, name="Program") +
  labs(x="Depth (x)", y="Time (s/Mbp)", title=element_blank()) +
  theme_minimal(base_size=10)+
  theme(panel.grid.minor=element_blank(),
        panel.grid.major.x=element_blank()) + coord_cartesian(ylim=c(0,350))

ggarrange(plt_n50_single, plt_length_single, plt_ram_mb_single, 
          plt_busco_single, plt_busco_copy_single, plt_time_mb_single,
          ncol=3, nrow=2, labels=c("A", "B", "C", "D", "E", "F"),
          common.legend=TRUE, legend="bottom") 
} -> plt_empirical_single

save_plot(
  plt_empirical_single,
  "fig08_empirical_single",
  10,6
)

#
# Empirical Mismatched Polished 
#
# Figure S5
empirical_stats_multi_plt <- empirical_stats_multi %>%
  filter(program %in% program_order) %>%
  mutate(depth_np=factor(depth_np),
          depth_il=factor(depth_il),
          program=factor(program, levels=program_order))

plt_polished_mismatches <- empirical_stats_multi_plt %>%
  ggplot(aes(x=depth_il, y=mean_mism, color=program, alpha=depth_np)) +
  geom_line(aes(group=interaction(program, depth_np))) +
  geom_point(size=2) +
  geom_linerange(aes(ymin=(mean_mism - CI_mism), ymax=(mean_mism + CI_mism)),
                  linewidth=0.5, show.legend = F) +
  #scale_x_continuous(breaks=c(10,15,20,25,30,35,40,50,60)) +
  scale_alpha_ordinal(name="Depth Long-reads (x)", range=c(0.25,1), breaks=c(10,30,60))+
  scale_color_manual(values=colors_programs, label=program_labels, name="Program") +
  labs(x="Depth Short-reads (x)", y="Mismatches(/100kb)", title=element_blank()) +
  theme_minimal(base_size=14)+
  theme(panel.grid.minor=element_blank(),
        panel.grid.major.x=element_blank(),
        legend.position = "bottom", legend.margin=margin())

save_plot(
  plt_polished_mismatches,
  "figS05_mismatches_polished",
  10,6
)

#
# Empirical Size Estimation And Comparisson 
#
# Figure S6


plot_only <- c("spades_short", "spades_hybrid", "abyss_short", "kmc_il", "kmc_chopper",
               "flye-ovl1000","flye-ovl1500", "flye-ovl2000", "flye-ovl2500",
               "flye-ovl3000","flye-ovl1000_polypolish")

length_df <- data_empirical %>% unite("depth", depth_np, depth_il, sep="") %>%
  filter(str_detect(depth, "99")) %>% 
  select(sample, program, total_length) 

length_df <- metrics_empirical %>%
  mutate(program= case_when(
    technology == "illumina" ~ "kmc_il",
    technology == "chopper"  ~ "kmc_chopper",
    TRUE ~ "ukn"
  )) %>% 
  select(sample, min_bp, max_bp, mean_bp, program) %>% 
  rename(total_length=mean_bp) %>%
  bind_rows(length_df) %>% 
  filter(program %in% plot_only) %>% 
  mutate(program = factor(program, levels = plot_only))

length_df %>% filter(program %in% plot_only) %>%
  group_by(program) %>%
  mutate(total_length = total_length * 1e-6) %>%
  summarise(med_length=median(total_length), 
            mean_length=mean(total_length),
            CI_length=CI95(total_length),
            .groups="drop") %>%
  arrange(med_length) -> length_summary
  
size_estimation_plot <- length_df %>%
  filter(program %in% plot_only) %>%
  ggplot(aes(x=program, y=total_length, color=program)) +
  geom_point(size=3, alpha=0.15) +
  geom_line(aes(group=sample), alpha=0.15) +
  geom_boxplot() +
  labs(x="", y="Assembly size (Mbp)") +
  scale_y_continuous(labels=scales::label_number(scale=1e-06)) +
  scale_x_discrete(labels=program_labels) +
  theme_minimal(base_size=14) +
  theme(panel.grid.minor=element_blank(),
        panel.grid.major.x=element_blank(),
        legend.position="none",
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))


save_plot(
  size_estimation_plot,
  "figS06_size_estimation_plot",
  10,6
)

# ----------------------------------------------------------------------------
# Figures for Model 
# ----------------------------------------------------------------------------

filt_responses <- c(
                    "n50",
                    "busco_complete",
                    "# mism_100kbp",
                    "# indel_100kbp,
                    maxrss_gb,
                    ram_mb,
                    time_mb,
                    ram_mb_p,
                    time_mb_p,
                    elapsed_sec"
                  )

response_labels <- c(
  "n50"             = "N50",
  "busco_complete"  = "BUSCO completeness",
  "maxrss_gb"       = "Max RAM",
  "ram_mb_p"        = "Max RAM (/Mbp)\nEstimated Size",
  "time_mb_p"       = "Run Time (/Mbp)\nEstimated Size",
  "ram_mb"          = "Max RAM (/Mbp)\nAssembly Size",
  "time_mb"         = "Run Time (/Mbp)\nAssembly Size",
  "elapsed_sec"     = "Run Time",
  "mism_100kbp"     = "Mismatches (/100kb)",
  "indel_100kbp"    = "Indels (/100kb)"
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
                               "genome_size_est" = "Genome Size Estimation",
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
        position = position_dodge(width = 0.7)
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
        panel.grid.major.y = element_line(color = "gray92", linewidth = 20),
        panel.grid.major.x = element_line(color = "gray", linewidth = 0.75),
        axis.text.y.right = element_blank(),
        axis.line.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        legend.position = "bottom"
      )
  }
  
  plt_size <- tmp_ploter(df %>% filter(covariate == "Genome Size Estimation"), "Genome Size Estimation") +
                rremove("y.text")
    

  if ("Depth LR" %in% unique(df$covariate) && "Depth SR" %in% unique(df$covariate)) {
    plt_lr <- tmp_ploter(df %>% filter(covariate == "Depth LR"), "Depth Long-reads") +
                  rremove("legend")
    
    plt_sr <- tmp_ploter(df %>% filter(covariate == "Depth SR"), "Depth Short-reads") +
                  # ggbreak::scale_x_break(c(0.03, 0.24), scales = "fixed") +
                  rremove("y.text") + rremove("legend")

    ggarrange(plt_lr, plt_sr, plt_size,
              ncol=3, nrow=1, labels=c("A", "B", "C"),
              common.legend = TRUE, legend="bottom",
              widths = c(1.2,1,1))
  } else {    
    plt_depth <- tmp_ploter(df %>% filter(covariate == "Depth"), "Depth") +
            guides(shape=guide_legend(nrow=1,byrow=F),
                    color=guide_legend(nrow=2,byrow=F)) 
    
    ggarrange(plt_depth, plt_size,
              ncol=2, nrow=1, labels=c("A", "B"),
              common.legend = TRUE, legend="bottom",
              heights = c(1.2,1))
  }
}

plt_forest_full_slopes_single <- forest_full_slopes(single_fits)
save_plot(
  plt_forest_full_slopes_single,
  "figS07_simulated_forest_full_slopes_single",
  15,7
)

plt_forest_full_slopes_multi <- forest_full_slopes(multi_fits)
save_plot(
  plt_forest_full_slopes_multi,
  "figS08_simulated_forest_full_slopes_multi",
  15,7
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
#       panel.grid.major.y = element_line(color = "gray92", linewidth = 15),
#       panel.grid.major.x = element_line(color = "gray", linewidth = 0.75),
#       axis.text.y.right = element_blank(),
#       axis.line.y.right = element_blank(),
#       axis.ticks.y.right = element_blank(),
#       plot.title = element_text(hjust = 0.5, vjust = -1),
#       legend.position = "bottom"
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
#               common.legend = TRUE, legend="bottom",
#               widths = c(1.2,1))
#   } else {
#     tmp_ploter(df, null)+
#     ggbreak::scale_x_break(c(0.1, 0.2), scales = "fixed") +
#       guides(shape=guide_legend(nrow=3,byrow=F),
#              color=guide_legend(nrow=2,byrow=F))
#   }
# }

# plt_forest_all_slopes_single <- forest_all_slopes(single_fits)
# plt_forest_all_slopes_multi <- forest_all_slopes(multi_fits)
# save_plot(
#   plt_forest_all_slopes_single,
#   "simulated_forest_all_slopes_single",
#   15,10
# )
# save_plot(
#   plt_forest_all_slopes_multi,
#   "simulated_forest_all_slopes_multi",
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
#       title = "Global effects (Type II ANOVA)",
#       x     = expression(chi^2),
#       y     = NULL,
#       shape = expression(italic(P)~"val")
#     ) +
#     theme_minimal(base_size = 12) +
#     scale_color_manual(values=color_palette,  name="Response") +
#     theme(
#       panel.grid.minor = element_blank(),
#       legend.position  = "bottom",
#       plot.title = element_text(hjust = 0.5, vjust = -1),
#       panel.grid.major.x = element_line(color = "gray", size = 0.75),
#     ) +
#     guides(shape=guide_legend(nrow=3,byrow=F),
#            color=guide_legend(nrow=3,byrow=F)
#            ) -> tmp_plt
  
#   if ("Depth LR" %in% df$term){
#     tmp_plt + theme(
#       panel.grid.major.y = element_line(color = "gray92", size = 8),
#     ) 
#   } else {
#     tmp_plt + theme(
#       panel.grid.major.y = element_line(color = "gray92", size = 12),
#     )
#   }
# }

# plt_forest_all_global_single <- forest_all_global(single_fits)
# save_plot(
#   plt_forest_all_global_single,
#   "simulated_forest_all_global_single",
#   15,10
# )

# plt_forest_all_global_multi <- forest_all_global(multi_fits)
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
#       response = fct_relevel(response, response_labels)
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
#       panel.grid.major.y = element_line(color = "gray92", linewidth = 15),
#       panel.grid.major.x = element_line(color = "gray", linewidth = 0.75),
#       axis.text.y.right = element_blank(),
#       axis.line.y.right = element_blank(),
#       axis.ticks.y.right = element_blank(),
#       legend.position = "bottom"
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
#       legend.position  = "bottom"
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