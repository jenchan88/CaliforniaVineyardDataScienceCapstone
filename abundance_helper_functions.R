library(here)
library(tidyverse)
library(paletteer)
library(broom)
library(vegan)
library(ggplot2)
library(patchwork)
library(stringr)

get_correlation_matrix <- function(data, taxa, treatment){
  
  cb_trt <- data |>
    filter(Treatment == treatment) |>
    select(all_of(taxa), Sample_ID, Count) |>
    pivot_wider(names_from = Sample_ID, values_from = Count, values_fill = 0)
  
  taxa_labels <- cb_trt[[taxa]]
  
  
  use <- cb_trt |>
    select(-Class)
  
  filtered <- use[, apply(use, 2, sd) != 0]
  
  # Compute the correlation matrix if no constant species remain
  correlation_matrix <- cor(t(filtered))
  
  colnames(correlation_matrix) <- taxa_labels
  rownames(correlation_matrix) <- taxa_labels
  
  return(correlation_matrix)
}

get_correlation_heatmap_2 <- function(matrix, taxa){
  heatmap <- matrix[taxa, taxa]
  heatmap[is.na(heatmap)] <- 0
  library(reshape2)
  
  heatmap_long <- melt(heatmap)
  
  ggplot(heatmap_long, aes(x = Var1, y = Var2, fill = value)) +
    geom_tile() +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, limits = c(-1, 1)) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 20, hjust = 0)) +
    xlab("") + 
    ylab("") +
    scale_y_discrete(limits=rev) +
    scale_x_discrete(position = "top")
}

get_correlation_heatmap_differential <- function(matrix, taxa){
  heatmap <- matrix[taxa, taxa]
  heatmap[is.na(heatmap)] <- 0
  library(reshape2)
  
  heatmap_long <- melt(heatmap)
  
  max_val <- max(abs(heatmap_long$value), na.rm = TRUE)
  
  ggplot(heatmap_long, aes(x = Var1, y = Var2, fill = value)) +
    geom_tile() +
    scale_fill_gradientn(colors = c("blue", "white", "red"), 
                         limits = c(-2, 2), 
                         values = scales::rescale(c(-2, -0.5, 0, 0.5, 2)),
                         oob=scales::squish) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 20, hjust = 0)) +
    xlab("") + 
    ylab("") +
    scale_y_discrete(limits=rev) +
    scale_x_discrete(position = "top")
}

get_correlation_difference_matrix <- function(t1_correlations, t2_correlations) {
  # Standardize correlations
  std_function <- function(r) {
    if (is.na(r)) {
      return(NA)
    } else if (abs(r) >= 1) {
      return(NA)
    } else {
      return((abs(r) / sqrt((1 - abs(r)) * abs(r))) * sign(r))
    }
  }
  
  t1_std <- apply(t1_correlations, c(1, 2), std_function)
  t2_std <- apply(t2_correlations, c(1, 2), std_function)
  
  # Compute the difference matrix
  diff_std <- t1_std - t2_std
  return(diff_std)
}

get_correlation_distance_matrix <- function(t1_correlations, t2_correlations) {
  # Standardize correlations
  std_function <- function(r) {
    if (is.na(r)) {
      return(NA)
    } else if (abs(r) >= 1) {
      return(NA)
    } else {
      return((abs(r) / sqrt((1 - abs(r)) * abs(r))) * sign(r))
    }
  }
  
  t1_std <- apply(t1_correlations, c(1, 2), std_function)
  t2_std <- apply(t2_correlations, c(1, 2), std_function)
  
  # Compute the difference matrix
  diff_std <- t1_std - t2_std
  
  # Compute distances
  inverse_function <- function(x) {
    if (is.na(x) || x == 0) {
      return(Inf)  # Prevent division by zero
    } else {
      return(1 / abs(x))
    }
  }
  
  dist_matrix <- apply(diff_std, c(1, 2), inverse_function)
  
  return(dist_matrix)
}

plot_relative_abundances <- function(data, taxa, taxa_list, experiment, t1, t2) {
  df <- data |>
    filter({{taxa}} %in% taxa_list) |>
    filter(Experiment == experiment) |>
    filter(Treatment %in% c(t1, t2)) |>
    group_by({{taxa}}, Treatment) |>
    summarize(Total_Count = sum(Count), .groups = 'drop')
  
  p <- ggplot(df, aes(x = {{taxa}}, y = Total_Count, fill = Treatment)) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(x = as_label(enquo(taxa)), y = "Total Abundance", title = "Relative Abundance by Treatment") +
    theme_minimal()
  
  return(p)
}

plot_std_counts_over_time <- function(data, taxa, taxa_list, experiment, t1, t2) {
  
  data$Sample <- paste(data$Season, data$Replicate, data$Year, sep = "_")
  
  sample_levels <- c(
    "Bud_A_2016", "Bud_B_2016", "Bud_C_2016",
    "Bloom_A_2016", "Bloom_B_2016", "Bloom_C_2016",
    "Veraison_A_2016", "Veraison_B_2016", "Verasion_C_2016",
    "Harvest_A_2016", "Harvest_B_2016", "Harvest_C_2016",
    "Bud_A_2017", "Bud_B_2017", "Bud_C_2017",
    "Bloom_A_2017", "Bloom_B_2017", "Bloom_C_2017",
    "Veraison_A_2017", "Veraison_B_2017", "Verasion_C_2017",
    "Harvest_A_2017", "Harvest_B_2017", "Harvest_C_2017",
    "Bud_A_2018", "Bud_B_2018", "Bud_C_2018",
    "Bloom_A_2018", "Bloom_B_2018", "Bloom_C_2018",
    "Veraison_A_2018", "Veraison_B_2018", "Verasion_C_2018",
    "Harvest_A_2018", "Harvest_B_2018", "Harvest_C_2018"
  )
  
  df <- data |>
    filter({{taxa}} %in% taxa_list, Experiment == experiment) |>
    select({{taxa}}, Sample, Treatment, Count) |>
    group_by({{taxa}}, Sample, Treatment) |>
    summarize(Count = sum(Count), .groups = "drop") |>
    group_by({{taxa}}) |>
    mutate(Scaled_Count = (Count - mean(Count)) / sd(Count)) |>
    ungroup() |>
    mutate(Sample = factor(Sample, levels = sample_levels, ordered = TRUE)) |>
    filter(!is.na(Sample))
  
  T1_counts <- df |> 
    filter(Treatment == t1)
  
  T2_counts <- df |> 
    filter(Treatment == t2)
  
  plot_T1 <- ggplot(T1_counts, aes(x = Sample, y = Scaled_Count, color = {{taxa}}, group = {{taxa}})) +
    geom_line() +
    ggtitle(t1) +
    theme_minimal() +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  
  plot_T2 <- ggplot(T2_counts, aes(x = Sample, y = Scaled_Count, color = {{taxa}}, group = {{taxa}})) +
    geom_line() +
    ggtitle(t2) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  combined_plot <- (plot_T1 / plot_T2) + plot_layout(guides = "collect") & theme(legend.position = "right")
  
  return(combined_plot)
}

plot_raw_counts_by_treatment <- function(data, taxa, taxa_list, experiment, t1, t2) {
  
  # Create Sample name
  data$Sample <- paste(data$Season, data$Replicate, data$Year, sep = "_")
  
  # Optional: preserve consistent sample order
  sample_levels <- unique(data$Sample)
  
  df <- data |>
    filter({{taxa}} %in% taxa_list, Experiment == experiment, Treatment %in% c(t1, t2)) |>
    select({{taxa}}, Sample, Treatment, Count) |>
    group_by({{taxa}}, Sample, Treatment) |>
    summarize(Count = sum(Count), .groups = "drop") |>
    filter(!is.na(Sample)) |>
    mutate(Sample = factor(Sample, levels = sample_levels, ordered = TRUE))
  
  # Split by treatment
  T1_df <- df |> filter(Treatment == t1)
  T2_df <- df |> filter(Treatment == t2)
  
  plot_T1 <- ggplot(T1_df, aes(x = Sample, y = Count, color = {{taxa}})) +
    geom_point(size = 2, alpha = 0.8, position = position_jitter(width = 0.2)) +
    ggtitle(t1) +
    theme_minimal() +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
    ylab("Raw Count")
  
  plot_T2 <- ggplot(T2_df, aes(x = Sample, y = Count, color = {{taxa}})) +
    geom_point(size = 2, alpha = 0.8, position = position_jitter(width = 0.2)) +
    ggtitle(t2) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ylab("Raw Count") +
    xlab("Sample")
  
  # Combine plots vertically with shared legend
  combined_plot <- (plot_T1 / plot_T2) + patchwork::plot_layout(guides = "collect") &
    theme(legend.position = "right")
  
  return(combined_plot)
}


