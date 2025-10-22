corr_uv_df <- function(data, taxa, t1, t2, experiment, treatment1, treatment2){
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  
  data$Sample <- paste(data$Season, data$Replicate, data$Year, sep = "_")
  
  df <- data |>
    filter({{taxa}} %in% c(t1, t2)) |>
    filter(Experiment == experiment) |>
    filter(Treatment %in% c(treatment1, treatment2)) |>
    select({{taxa}}, Sample, Treatment, Count, Year, Season) |>
    group_by({{taxa}}) |>
    mutate(Scaled_Count = (Count - mean(Count)) / sd(Count)) |>
    ungroup() |>
    filter(!is.na(Sample))
  
  df_t1 <- df |> 
    filter({{taxa}} == t1) |> 
    rename(Scaled_t1 = Scaled_Count)
  df_t2 <- df |> 
    filter({{taxa}} == t2) |> 
    rename(Scaled_t2 = Scaled_Count)
  
  df_uivi <- df_t1 |> 
    inner_join(df_t2, by = c("Sample", "Treatment", "Year", "Season")) |> 
    mutate(uivi = Scaled_t1 * Scaled_t2) |>
    mutate(Year = as.character(Year)) |>
    mutate(Season = factor(Season, levels= c("Bud", "Bloom", "Veraison", "Harvest"), ordered=TRUE)) |>
    group_by(Treatment, Season, Year) |>
    summarize(uivi = sum(uivi)) |>
    mutate(Month = case_when(Season == "Bud" ~ "03",
                             Season == "Bloom" ~ "05",
                             Season == "Veraison" ~ "07",
                             Season == "Harvest" ~ "09"
    ),
    date = make_date(Year,Month,1))
  
  return(data = df_uivi)
}

plot_rainfall_overlay <- function(rainfall_df, uivi_df, tax_names) {
  
  taxa_string <- if (length(tax_names) == 2) {
    paste(tax_names[1], tax_names[2], sep = " and ")
  } else {
    paste(paste(tax_names[-length(tax_names)], sep = ", "), "and", tax_names[length(tax_names)])
  }
  plot_title <- paste("Correlation Contribution plot for", taxa_string)
  
  ggplot() +
    geom_bar(data = rainfall_df, 
             aes(x = date, y = cumulative_rainfall), 
             stat = "identity", fill = "skyblue") +
    
    geom_point(data = uivi_df, 
               aes(x = date, y = uivi, color = Treatment), 
               size = 2) +
    
    geom_line(data = uivi_df, 
              aes(x = date, y = uivi, group = date), 
              color = "gray50", linewidth = 0.5) +
    
    scale_x_date(date_labels = "%b %Y", date_breaks = "1 month") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = plot_title, 
         x = "Date", y = "")
}

plot_CN_overlay <- function(CN_df, uivi_df, tax_names, experiment, t1, t2) {
  
  # Format taxa for title
  taxa_string <- if (length(tax_names) == 2) {
    paste(tax_names, collapse = " and ")
  } else {
    paste(paste(tax_names[-length(tax_names)], collapse = ", "), "and", tax_names[length(tax_names)])
  }
  
  plot_title_prefix <- paste("Overlay for", taxa_string, "in", experiment)
  
  CN_exp <- CN_df |>
    filter(Treatment %in% c(t1, t2)) |>
    select(Treatment, Season, Year, `%C (average)`, `%N (average)`, `C/N (average)`)
  
  CN_overlay <- merge(uivi_df, CN_exp) |>
    select(-c(Month, date)) |>
    mutate(Sample = factor(
      paste0(Season, "_", Year),
      levels = c("Bud_2016", "Bloom_2016", "Veraison_2016", "Harvest_2016", 
                 "Bud_2017", "Bloom_2017", "Veraison_2017", "Harvest_2017", 
                 "Bud_2018", "Bloom_2018", "Veraison_2018", "Harvest_2018"),
      ordered = TRUE))
  
  make_plot <- function(y_var, fill_color, label) {
    ggplot() +
      geom_col(data = CN_overlay, aes(x = Sample, y = .data[[y_var]]), 
               position = position_dodge(width = 0.9), fill = fill_color) +
      geom_point(data = CN_overlay, aes(x = Sample, y = uivi, color = Treatment), size = 2) +
      geom_line(data = CN_overlay, aes(x = Sample, y = uivi, group = Sample), 
                color = "gray50", linewidth = 0.5) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            axis.ticks.y = element_blank(), 
            axis.text.y = element_blank()) +
      labs(title = paste(plot_title_prefix, "-", label),
           x = "Sample", y = "")
  }
  
  list(
    plot_percent_C = make_plot("%C (average)", "orange3", "%C"),
    plot_percent_N = make_plot("%N (average)", "darkolivegreen3", "%N"),
    plot_CN_ratio  = make_plot("C/N (average)", "orange3", "C/N Ratio")
  )
}

plot_CN_overlay_scaled <- function(CN_df, uivi_df, tax_names, experiment, t1, t2) {
  
  # Format taxa for title
  taxa_string <- if (length(tax_names) == 2) {
    paste(tax_names, collapse = " and ")
  } else {
    paste(paste(tax_names[-length(tax_names)], collapse = ", "), "and", tax_names[length(tax_names)])
  }
  
  plot_title_prefix <- paste("Overlay for", taxa_string, "in", experiment)
  
  CN_exp <- CN_df |>
    filter(Treatment %in% c(t1, t2)) |>
    select(Treatment, Season, Year, `%C (average)`, `%N (average)`, `C/N (average)`)
  
  CN_overlay <- merge(uivi_df, CN_exp) |>
    select(-c(Month, date)) |>
    mutate(Sample = factor(
      paste0(Season, "_", Year),
      levels = c("Bud_2016", "Bloom_2016", "Veraison_2016", "Harvest_2016", 
                 "Bud_2017", "Bloom_2017", "Veraison_2017", "Harvest_2017", 
                 "Bud_2018", "Bloom_2018", "Veraison_2018", "Harvest_2018"),
      ordered = TRUE)) |>
    group_by(Treatment) |>
    mutate(uivi_scaled = scales::rescale(uivi, to = c(0, 10))) |>  # adjust this range
    ungroup()
  
  make_plot <- function(y_var, fill_color, label, ylim = NULL) {
    p <- ggplot() +
      geom_col(data = CN_overlay, aes(x = Sample, y = .data[[y_var]]), 
               position = position_dodge(width = 0.9), fill = fill_color) +
      geom_point(data = CN_overlay, aes(x = Sample, y = uivi_scaled, color = Treatment), size = 2) +
      geom_line(data = CN_overlay, aes(x = Sample, y = uivi_scaled, group = Treatment), 
                color = "gray50", linewidth = 0.5) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), 
            axis.ticks.y = element_blank(), 
            axis.text.y = element_blank()) +
      labs(title = paste(plot_title_prefix, "-", label),
           x = "Sample", y = "")
    
    if (!is.null(ylim)) {
      p <- p + coord_cartesian(ylim = ylim)
    }
    
    return(p)
  }
  
  list(
    plot_percent_C = make_plot("%C (average)", "gray50", "%C"),
    plot_percent_N = make_plot("%N (average)", "darkolivegreen3", "%N"),
    plot_CN_ratio  = make_plot("C/N (average)", "gray50", "C/N Ratio")
  )
}