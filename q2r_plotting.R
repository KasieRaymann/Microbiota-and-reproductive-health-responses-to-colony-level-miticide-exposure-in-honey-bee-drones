# QIIME2R PLOTTING PIPELINE

# Author:        Layla Jeries, with some ChatGPT assistance
# Date created:  2026-09-17
# Last modified: 2026-11-06

########################################## AI USAGE DISCLAIMER ###########################################
## This disclaimer has been included because some or all parts of this script were written with ChatGPT ##
## assistance. The author accepts full responsibility for the entirety of this script. This disclaimer, ##
## in-line comments, and all other documentation are written by the author and without the aid of       ##
## ChatGPT or a similar AI/chatbot.                                                                     ##
##########################################################################################################

#### BUILD PHYLOSEQ OBJECT ####
# Set correct working directory
setwd("path/to/analysis-folder")

##### Import libraries ####

library(tidyverse)    # Contains packages for core data analyses
library(phyloseq)     # For construction and analyses of phyloseq objects
library(qiime2R)      # For read-in of qiime2 .qza files
library(ggsignif)     # For plotting significance bars on alpha diversity boxplots

##### Generate phyloseq object ####

###### qza to phyloseq #####
# Confirm your filenames!
physeq <- qza_to_phyloseq(
  features = "feature_table.qza",
  tree     = "rooted-tree.qza",
  taxonomy = "taxonomy.qza",
  metadata = "mapping.tsv"     # Must be a qiime2 metadata file in .tsv format with "#q2:types" line!!!
)

###### double check phyloseq object #####
physeq
sample_names(physeq)
sample_variables(physeq)
rank_names(physeq)
ntaxa(physeq)
nsamples(physeq)

##### Organize metadata ####

###### load metadata as a dataframe #####
metadata_df <- data.frame(sample_data(physeq)) |>
  tibble::rownames_to_column("SampleID")

###### format metadata for consistently ordered legends #####
metadata_df$C_T <- factor(
  metadata_df$C_T,
  levels = c("Control", "Treatment")
)
metadata_df$Best_Worst <- factor(
  metadata_df$Best_Worst,
  levels = c("Best", "Average", "Worst")
)
metadata_df$CT_BestWorst <- factor(
  metadata_df$CT_BestWorst,
  levels = c("ControlBest", "ControlAverage", "ControlWorst", "TreatmentBest", "TreatmentAverage", "TreatmentWorst")
)

#### PLOT ALPHA DIVERSITY ####
# Set correct working directory
setwd("~/path/to/core_diversity")

##### Prepare data ####

###### load in qiime2 vectors #####
shannon_df <- read_qza("shannon_vector.qza")$data |>
  tibble::rownames_to_column("SampleID")
faithpd_df <- read_qza("faith_pd_vector.qza")$data |>
  tibble::rownames_to_column("SampleID")

###### create alpha diversity dataframe #####
# Other alpha diversity metrics may be added.
alpha_df <- shannon_df |>
  left_join(faithpd_df, by = "SampleID") |>
  left_join(metadata_df, by = "SampleID")

##### Function configurations ####

###### define metric labels #####
# Do NOT change the name of the alpha_labels variable.
# Other alpha diversity metrics may be added.
alpha_labels <- c(
  shannon_entropy   = "Shannon Diversity",
  faith_pd          = "Faith's Phylogenetic Diversity"
)

###### define groups as x-axis labels #####
# Do NOT change the name of the group_labels variable.
group_labels <- c(
  C_T          = "Treatment Group",
  Best_Worst   = "Score",
  CT_BestWorst = "Treatment Group x Score"
)

###### define color palettes #####
# Do NOT change the name of the palettes variable.
palettes <- list(
  
  C_T         = c(
    Control   = "#67A9CF",
    Treatment = "#F1A340"
  ),
  
  Best_Worst = c(
    Best    = "#0067A5",
    Average = "#FFBF00",
    Worst   = "#F04923"
  ),
  
  CT_BestWorst       = c(
    ControlBest      = "#D1E5F0",
    ControlAverage   = "#67A9CF",
    ControlWorst     = "#2166AC",
    TreatmentBest    = "#FFCE8A",
    TreatmentAverage = "#F1A340",
    TreatmentWorst   = "#B35806"
  )
)

##### Function ####

plot_alpha <- function(data, group, metric, title_add, sig_file) {
  p <- ggplot(data,aes(x = .data[[group]],
                       y = .data[[metric]],
                       fill = .data[[group]]
  )
  ) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.15) +
    scale_fill_manual(
      values = palettes[[group]]
    ) +
    labs(
      title = paste(alpha_labels[[metric]], title_add, sep = " | "),
      x = group_labels[[group]],
      y = alpha_labels[[metric]]
    ) +
    theme_bw() +
    theme(
      plot.title  = element_text(face = "bold", hjust = 0.5),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  
  # Read significance file
  sig <- read.csv(
    sig_file,
    check.names = FALSE
  )
  
  annotation_df <- sig %>%
    filter(`q-value` < 0.05) %>%
    mutate(
      group1 = sub(" \\(n=.*\\)$", "", `Group 1`),
      group2 = sub(" \\(n=.*\\)$", "", `Group 2`),
      
      groups = purrr::map2(
        group1,
        group2,
        c
      ),
      
      label = case_when(
        `q-value` < 0.001 ~ "***",
        `q-value` < 0.01  ~ "**",
        `q-value` < 0.05  ~ "*"
      )
    )
  
  # Make significance comparisons
  if (nrow(annotation_df) > 0) {
    
    y_max <- max(
      data[[metric]],
      na.rm = TRUE
    )
    
    y_range <- diff(
      range(
        data[[metric]],
        na.rm = TRUE
      )
    )
    
    annotation_df$y_position <- y_max +
      seq(
        0.10 * y_range,
        0.10 * y_range * nrow(annotation_df),
        by = 0.10 * y_range
      )
    
    p <- p +
      ggsignif::geom_signif(
        comparisons = annotation_df$groups,
        annotations = annotation_df$label,
        y_position = annotation_df$y_position,
        textsize = 6,
        vjust = 0.75,
        tip_length = 0.01
      )
  }
  
  return(p)
}

##### Plotting example for treatment group ####
###### for C_T #####
shannon_plot_CT <- plot_alpha(
  alpha_df,
  group  = "C_T",
  metric = "shannon_entropy",
  title_add = "Guts",
  sig_file = "sig_csvs_shannon/kruskal-wallis-pairwise-C_T.csv"
)
shannon_plot_CT

faithpd_plot_CT <- plot_alpha(
  alpha_df,
  group  = "C_T",
  metric = "faith_pd",
  title_add = "Guts",
  sig_file = "sig_csvs_faith_pd/kruskal-wallis-pairwise-C_T.csv"
)
faithpd_plot_CT

#### PLOT BETA DIVERSITY ####
##### Prepare data ####

###### load in qiime2 vectors #####
bray_pcoa <- read_qza("bray_curtis_pcoa_results.qza")
weighted_pcoa <- read_qza("weighted_unifrac_pcoa_results.qza")
unweighted_pcoa <- read_qza("unweighted_unifrac_pcoa_results.qza")

###### create separate beta diversity dataframes #####
bray_df <- bray_pcoa$data$Vectors |>
  left_join(metadata_df, by = "SampleID")

weighted_df <- weighted_pcoa$data$Vectors |>
  left_join(metadata_df, by = "SampleID")

unweighted_df <- unweighted_pcoa$data$Vectors |>
  left_join(metadata_df, by = "SampleID")

##### Function configurations ####

###### define metric labels #####
# Do NOT change the name of the beta_labels variable.
# Other beta diversity metrics may be added.
beta_labels <- c(
  bray = "Bray-Curtis",
  weighted = "Weighted UniFrac",
  unweighted = "Unweighted UniFrac"
)

###### define shape labels (optional) #####
# Other shapes may be added.
# You must match names of the variables storing each shape vector c() to how they appear in the metadata file.
shape_palettes <- list(
  
  Sample_Type = c(
    Gut   = 16,
    Repro = 17
  ),
  
  C_T = c(
    Control = 1,
    Treatment = 16
  )
  
)

###### define color palettes #####
# Color palettes should have already been defined in script "2_plotalpha."
# Can be rerun now if need be.
# OR can define a new set of color palettes for beta diversity plots (store in variable ""palettes")

##### Function ####

plot_pcoa <- function(data, ordination, group, title, shape_group = NULL, ellipse_group = NULL) {
  
  # Create PC labels
  pc1_label <- paste0(
    "PC1 (",
    round(ordination$data$ProportionExplained[1, "PC1"] * 100, 1),
    "%)"
  )
  pc2_label <- paste0(
    "PC2 (",
    round(ordination$data$ProportionExplained[1, "PC2"] * 100, 1),
    "%)"
  )
  
  # Plot
  p <- ggplot(
    data,
    aes(
      x = PC1,
      y = PC2,
      color = .data[[group]]
    )
  )
  
  # Add ellipses if requested
  if (!is.null(ellipse_group)) {
    
    p <- p +
      stat_ellipse(
        aes(
          group = .data[[ellipse_group]]
        ),
        type = "norm",
        linewidth = 0.8
      )
    
  }
  
  # Add shapes if requested
  if (!is.null(shape_group)) {
    
    p <- p +
      geom_point(
        aes(shape = .data[[shape_group]]),
        size = 3
      ) +
      scale_shape_manual(
        values = shape_palettes[[shape_group]]
      )
    
  } else {
    
    p <- p +
      geom_point(size = 3)
    
  }
  
  p <- p +
    scale_color_manual(
      values = palettes[[group]]
    ) +
    labs(
      title = title,
      x = pc1_label,
      y = pc2_label,
      color = group_labels[[group]]
    ) +
    coord_fixed(ratio = 1) +
    theme_bw() +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      aspect.ratio = 1
    )
  
  return(p)
  
}

##### Plotting example for treatment group ####
###### for C_T ######
bray_plot_CT <- plot_pcoa(
  bray_df,
  bray_pcoa,
  group = "C_T",
  shape_group = "Sample_Type",
  ellipse_group = "C_T",
  title = "Bray-Curtis for Treatment Group | Guts"
)
bray_plot_CT

wt_plot_CT <- plot_pcoa(
  weighted_df,
  weighted_pcoa,
  group = "C_T",
  shape_group = "Sample_Type",
  ellipse_group = "C_T",
  title = "Weighted UniFrac for Treatment Group | Guts"
)
wt_plot_CT

uwt_plot_CT <- plot_pcoa(
  unweighted_df,
  unweighted_pcoa,
  group = "C_T",
  shape_group = "Sample_Type",
  ellipse_group = "C_T",
  title = "Unweighted UniFrac for Treatment Group | Guts"
)
uwt_plot_CT