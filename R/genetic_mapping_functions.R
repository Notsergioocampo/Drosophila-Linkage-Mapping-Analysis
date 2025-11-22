#' Build Triple Point Map
#'
#' This function takes a dataframe of Drosophila class counts and phenotypes,
#' classifies them into Parental (P), Single Crossover (SCO), and Double Crossover (DCO)
#' classes, and calculates recombination frequencies and map distances (in cM).
#'
#' @param df A dataframe containing at least 'phenotype' and 'count' columns.
#'           It is expected to have specific phenotype strings:
#'           "all wt", "eye-", "wing-", "body-", "eye- wing-", "eye- body-",
#'           "body- wing-", "eye- wing- body-"
#'
#' @return A list containing:
#'         - data: The original dataframe with an added 'class' column.
#'         - N: Total number of offspring.
#'         - distances: A named vector of map distances (EW, WB, EB).
#'         - map_summary: A text summary of the map.
#' @export
build_triple_map <- function(df) {
  # 1. Total offspring
  N <- sum(df$count)

  # 2. Classify offspring into P, SCO_EW, SCO_WB, DCO
  # Gene order: eye - wing - body

  df$class <- NA_character_

  # Parental: all wt (000), eye- wing- body- (111)
  df$class[df$phenotype == "all wt"] <- "P"
  df$class[df$phenotype == "eye- wing- body-"] <- "P"

  # DCO: wing- (010), eye- body- (101)
  # (Middle gene 'wing' switches)
  df$class[df$phenotype == "wing-"] <- "DCO"
  df$class[df$phenotype == "eye- body-"] <- "DCO"

  # SCO_EW: eye- (100), body- wing- (011)
  # (Crossover between eye and wing)
  df$class[df$phenotype == "eye-"] <- "SCO_EW"
  df$class[df$phenotype == "body- wing-"] <- "SCO_EW"

  # SCO_WB: body- (001), eye- wing- (110)
  # (Crossover between wing and body)
  df$class[df$phenotype == "body-"] <- "SCO_WB"
  df$class[df$phenotype == "eye- wing-"] <- "SCO_WB"

  # 3. Summarize counts by recombination class
  counts_by_class <- tapply(df$count, df$class, sum)

  # Handle potential missing classes (set to 0 if NA)
  if (is.na(counts_by_class["SCO_EW"])) counts_by_class["SCO_EW"] <- 0
  if (is.na(counts_by_class["SCO_WB"])) counts_by_class["SCO_WB"] <- 0
  if (is.na(counts_by_class["DCO"])) counts_by_class["DCO"] <- 0

  SCO_EW <- as.numeric(counts_by_class["SCO_EW"])
  SCO_WB <- as.numeric(counts_by_class["SCO_WB"])
  DCO <- as.numeric(counts_by_class["DCO"])

  # 4. Compute recombination frequencies & map distances

  # 4a. eye-wing interval
  recomb_EW <- (SCO_EW + DCO) / N
  dist_EW_cM <- recomb_EW * 100

  # 4b. wing-body interval
  recomb_WB <- (SCO_WB + DCO) / N
  dist_WB_cM <- recomb_WB * 100

  # 4c. outer interval eye-body (2-point)
  # For outer genes, recombinants are where eye and body differ (10 or 01)
  outer_recomb_classes <- c("eye-", "body-", "eye- wing-", "body- wing-")
  recomb_EB <- sum(df$count[df$phenotype %in% outer_recomb_classes]) / N
  dist_EB_cM <- recomb_EB * 100

  # Total Map Length (sum of intervals)
  total_map_dist <- dist_EW_cM + dist_WB_cM

  # 5. Create summary text
  summary_text <- paste0(
    "Total N = ", N, "\n\n",
    "Gene order (from 3-point cross): eye - wing - body\n\n",
    "Map distances:\n",
    sprintf("  eye  - wing  : %5.2f cM\n", dist_EW_cM),
    sprintf("  wing - body  : %5.2f cM\n", dist_WB_cM),
    "  -------------------------\n",
    sprintf("  Total Map Length: %5.2f cM (eye - body)\n\n", total_map_dist),
    sprintf("  (Uncorrected 2-point fraction: %5.2f cM)\n", dist_EB_cM),
    "   *The 2-point fraction is smaller because it misses double crossovers.\n"
  )

  list(
    data = df,
    N = N,
    distances = c(EW = as.numeric(dist_EW_cM), WB = as.numeric(dist_WB_cM), EB = as.numeric(dist_EB_cM), Total = as.numeric(total_map_dist)),
    map_summary = summary_text
  )
}

#' Draw Genetic Map
#'
#' Creates a publication-quality visualization of the genetic map using ggplot2.
#'
#' @param map_results The list output from build_triple_map().
#' @return A ggplot object.
#' @export
draw_genetic_map <- function(map_results) {
  require(ggplot2)

  # Extract distances
  d_EW <- map_results$distances["EW"]
  d_WB <- map_results$distances["WB"]

  # Define positions (cumulative distance)
  gene_pos <- data.frame(
    gene = c("Crimson (eye)", "No Take Off (wing)", "Spray Tan (body)"),
    position = c(0, d_EW, d_EW + d_WB),
    type = "gene"
  )

  total_len <- max(gene_pos$position)

  # Create plot
  p <- ggplot(gene_pos, aes(x = position, y = 0)) +
    # Chromosome rod (thicker, light grey base)
    annotate("segment",
      x = -5, xend = total_len + 5, y = 0, yend = 0,
      linewidth = 4, lineend = "round", color = "grey90"
    ) +

    # Active chromosome segment (darker)
    annotate("segment",
      x = 0, xend = total_len, y = 0, yend = 0,
      linewidth = 4, lineend = "round", color = "grey70"
    ) +

    # Gene markers (points)
    geom_point(aes(fill = gene), size = 6, shape = 21, color = "white", stroke = 1.5) +

    # Gene labels (above)
    geom_text(aes(label = gene, color = gene), y = 0.6, fontface = "bold", size = 5) +

    # Position labels (below markers)
    geom_text(aes(label = sprintf("%.1f cM", position)), y = -0.3, size = 4, color = "grey40") +

    # Interval 1: Eye - Wing
    annotate("segment",
      x = 0.5, xend = d_EW - 0.5, y = -0.7, yend = -0.7,
      arrow = arrow(length = unit(0.03, "inches"), ends = "both"), color = "grey60", linewidth = 0.5
    ) +
    annotate("text",
      x = d_EW / 2, y = -0.9,
      label = sprintf("%.1f cM", d_EW), fontface = "italic", color = "grey50", size = 4
    ) +

    # Interval 2: Wing - Body
    annotate("segment",
      x = d_EW + 0.5, xend = (d_EW + d_WB) - 0.5, y = -0.7, yend = -0.7,
      arrow = arrow(length = unit(0.03, "inches"), ends = "both"), color = "grey60", linewidth = 0.5
    ) +
    annotate("text",
      x = d_EW + (d_WB / 2), y = -0.9,
      label = sprintf("%.1f cM", d_WB), fontface = "italic", color = "grey50", size = 4
    ) +

    # Scales and Theme
    scale_x_continuous(limits = c(-10, total_len + 10)) +
    scale_y_continuous(limits = c(-2, 2)) +
    scale_fill_manual(values = c("#e74c3c", "#3498db", "#f1c40f")) + # Professional flat colors
    scale_color_manual(values = c("#c0392b", "#2980b9", "#f39c12")) +
    labs(
      title = "Genetic Map of Chromosome 2",
      subtitle = "Relative positions of cn, nto, and sp loci"
    ) +
    theme_publication() +
    theme(
      axis.line = element_blank(),
      axis.text = element_blank(),
      axis.title = element_blank(),
      axis.ticks = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "none"
    )

  return(p)
}

#' Plot Phenotype Counts
#'
#' Creates a bar chart of counts for each phenotype, colored by class.
#'
#' @param df The dataframe from build_triple_map() containing 'phenotype', 'count', and 'class'.
#' @return A ggplot object.
#' @export
plot_phenotype_counts <- function(df) {
  require(ggplot2)
  require(dplyr)

  # Order phenotypes by count
  df_plot <- df %>%
    mutate(phenotype = reorder(phenotype, -count))

  ggplot(df_plot, aes(x = phenotype, y = count, fill = class)) +
    geom_bar(stat = "identity", color = "black", width = 0.7) +
    scale_fill_brewer(palette = "Pastel1", name = "Recombination Class") +
    labs(
      title = "Distribution of F2 Phenotypes",
      x = "Phenotype",
      y = "Count (Number of Flies)"
    ) +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      legend.position = "top"
    )
}

#' Plot Class Distribution
#'
#' Creates a summary bar chart of total counts per recombination class.
#'
#' @param df The dataframe from build_triple_map().
#' @return A ggplot object.
#' @export
plot_class_distribution <- function(df) {
  require(ggplot2)
  require(dplyr)

  # Summarize by class
  df_summary <- df %>%
    group_by(class) %>%
    summarise(total_count = sum(count)) %>%
    mutate(class = factor(class, levels = c("P", "SCO_EW", "SCO_WB", "DCO")))

  ggplot(df_summary, aes(x = class, y = total_count, fill = class)) +
    geom_bar(stat = "identity", color = "black", width = 0.6) +
    geom_text(aes(label = total_count), vjust = -0.5, size = 5) +
    scale_fill_brewer(palette = "Pastel1") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) + # Add 15% space at top
    labs(
      title = "Total Counts by Recombination Class",
      subtitle = "Demonstrating Parental Excess and Double Crossover Rarity",
      x = "Recombination Class",
      y = "Total Count"
    ) +
    guides(fill = "none") + # Hide legend as x-axis is sufficient
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = element_text(hjust = 0.5, face = "italic"),
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 13)
    )
}

#' Calculate Chi-Square Statistics
#'
#' Calculates the Chi-Square statistic for the null hypothesis of independent assortment (1:1:1:1:1:1:1:1 ratio).
#'
#' @param df The dataframe from build_triple_map().
#' @return A list containing the Chi-Square table and the test result.
#' @export
calculate_chi_square <- function(df) {
  require(dplyr)

  N <- sum(df$count)
  expected_count <- N / 8

  chi_sq_table <- df %>%
    select(phenotype, count, class) %>%
    mutate(
      Expected = expected_count,
      Difference = count - Expected,
      Diff_Sq = Difference^2,
      Chi_Sq_Component = Diff_Sq / Expected
    ) %>%
    arrange(desc(count))

  total_chi_sq <- sum(chi_sq_table$Chi_Sq_Component)
  df_degrees <- 8 - 1
  p_value <- pchisq(total_chi_sq, df = df_degrees, lower.tail = FALSE)

  # Format for display
  display_table <- chi_sq_table %>%
    mutate(
      Expected = round(Expected, 2),
      Difference = round(Difference, 2),
      Diff_Sq = round(Diff_Sq, 2),
      Chi_Sq_Component = round(Chi_Sq_Component, 2)
    ) %>%
    rename(
      "Phenotype" = phenotype,
      "Observed" = count,
      "Class" = class,
      "Expected" = Expected,
      "O - E" = Difference,
      "(O - E)^2" = Diff_Sq,
      "Chi-Square" = Chi_Sq_Component
    )

  return(list(
    table = display_table,
    statistic = total_chi_sq,
    df = df_degrees,
    p_value = p_value,
    chi_sq_data = chi_sq_table # Return raw data for plotting
  ))
}

plot_chi_square_contribution <- function(df) {
  # Recalculate chi-square data locally to ensure we have the components
  N <- sum(df$count)
  expected_per_class <- N / 8

  chi_sq_data <- df %>%
    mutate(
      Expected = expected_per_class,
      Chi_Sq_Component = ((count - Expected)^2) / Expected
    ) %>%
    arrange(desc(Chi_Sq_Component))

  # Create the plot
  p <- ggplot(chi_sq_data, aes(x = reorder(phenotype, Chi_Sq_Component), y = Chi_Sq_Component, fill = class)) +
    geom_bar(stat = "identity", color = "black", width = 0.7) +
    coord_flip() +
    scale_fill_brewer(palette = "Pastel1", name = "Recombination Class") +
    labs(
      title = "Contribution to Chi-Square Statistic",
      subtitle = "Parental classes drive the deviation from independent assortment",
      x = "Phenotype",
      y = "Chi-Square Contribution ((O-E)^2 / E)"
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      axis.text = element_text(size = 10)
    )

  return(p)
}


#' Publication Quality Theme
#'
#' A clean, professional theme for ggplot2 figures.
#' @export
theme_publication <- function(base_size = 12, base_family = "sans") {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      plot.title = element_text(face = "bold", size = rel(1.2), hjust = 0.5),
      text = element_text(),
      panel.background = element_rect(colour = NA),
      plot.background = element_rect(colour = NA),
      panel.border = element_rect(colour = NA),
      axis.title = element_text(face = "bold", size = rel(1)),
      axis.title.y = element_text(angle = 90, vjust = 2),
      axis.title.x = element_text(vjust = -0.2),
      axis.text = element_text(),
      axis.line = element_line(colour = "black"),
      axis.ticks = element_line(),
      panel.grid.major = element_line(colour = "#f0f0f0"),
      panel.grid.minor = element_blank(),
      legend.key = element_rect(colour = NA),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.key.size = unit(0.5, "cm"),
      legend.title = element_text(face = "italic"),
      plot.margin = unit(c(10, 5, 5, 5), "mm"),
      strip.background = element_rect(colour = "#f0f0f0", fill = "#f0f0f0"),
      strip.text = element_text(face = "bold")
    )
}

#' Plot Phenotype Counts
#'
#' @param df Dataframe with 'phenotype' and 'count' columns.
#' @return A ggplot object.
#' @export
plot_phenotype_counts <- function(df) {
  ggplot(df, aes(x = reorder(phenotype, -count), y = count)) +
    geom_bar(stat = "identity", fill = "#2c3e50", width = 0.7) + # Professional dark blue
    labs(
      title = "Distribution of F2 Phenotypes",
      x = "Phenotype",
      y = "Count"
    ) +
    theme_publication() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

#' Plot Class Distribution
#'
#' @param df Dataframe with 'class' and 'count' columns.
#' @return A ggplot object.
#' @export
plot_class_distribution <- function(df) {
  # Define professional colors for classes
  class_colors <- c(
    "Parental" = "#2c3e50", # Dark Blue
    "SCO_EW" = "#e67e22", # Muted Orange
    "SCO_WB" = "#27ae60", # Muted Green
    "DCO" = "#c0392b" # Muted Red
  )

  ggplot(df, aes(x = class, y = count, fill = class)) +
    geom_bar(stat = "identity", width = 0.6) +
    scale_fill_manual(values = class_colors) +
    labs(
      title = "Crossover Class Frequency",
      x = "Crossover Class",
      y = "Number of Offspring"
    ) +
    theme_publication() +
    theme(legend.position = "none")
}

#' Plot Chi-Square Contribution
#'
#' @param df Dataframe with 'phenotype', 'observed', 'expected' columns.
#' @return A ggplot object.
#' @export
plot_chi_square_contribution <- function(df) {
  # Calculate contribution
  N <- sum(df$count)
  expected_val <- N / 8

  chi_sq_data <- df %>%
    mutate(
      observed = count,
      expected = expected_val,
      Chi_Sq_Component = ((observed - expected)^2) / expected,
      class = case_when(
        phenotype %in% c("Wild Type", "Crimson, No Take Off, Spray Tan") ~ "Parental",
        phenotype %in% c("Crimson", "No Take Off, Spray Tan") ~ "SCO (C-N)",
        phenotype %in% c("Spray Tan", "Crimson, No Take Off") ~ "SCO (N-S)",
        TRUE ~ "DCO"
      )
    )

  # Professional palette
  class_colors <- c(
    "Parental" = "#2c3e50",
    "SCO (C-N)" = "#e67e22",
    "SCO (N-S)" = "#27ae60",
    "DCO" = "#c0392b"
  )

  ggplot(chi_sq_data, aes(x = reorder(phenotype, Chi_Sq_Component), y = Chi_Sq_Component, fill = class)) +
    geom_bar(stat = "identity", width = 0.7) +
    scale_fill_manual(values = class_colors, name = "Class") +
    coord_flip() +
    labs(
      title = "Contribution to Chi-Square Statistic",
      subtitle = "Higher values indicate greater deviation from independent assortment",
      x = "Phenotype",
      y = "Chi-Square Contribution"
    ) +
    theme_publication()
}

#' Plot Genetic Map Comparison
#'
#' @param map_results List output from build_triple_map
#' @return A ggplot object
#' @export
plot_map_comparison <- function(map_results) {
  # Extract calculated distances
  dist_EW <- map_results$distances["EW"]
  dist_WB <- map_results$distances["WB"]

  # Create data for plotting
  # Calculated Map (Cumulative positions)
  calc_map <- tibble(
    Map = "Calculated (F2)",
    Gene = c("Crimson", "No Take Off", "Spray Tan"),
    Position = c(0, dist_EW, dist_EW + dist_WB),
    Label = c(paste0("cn (0.0)"), paste0("nto (", round(dist_EW, 1), ")"), paste0("sp (", round(dist_EW + dist_WB, 1), ")"))
  )

  # Standard Map (Approximate locations on Chr 2R)
  # cn = 57.5, vg = 67.0, bw = 104.5 (using standard markers as proxies for this example if exact loci unknown,
  # but let's stick to the relative distances for comparison or use the provided standard values if any.
  # Assuming standard values from literature for these specific mutants:
  # Let's assume for the visual comparison we align them at 0.
  # If we don't have exact standard loci, we can omit the standard map or use a placeholder.
  # For this report, let's visualize the CALCULATED map clearly.

  # Let's just plot the calculated map to be safe and clean, or compare to a "Theoretical" if we had one.
  # The previous code had a "Standard" map hardcoded. I will keep a simplified version.

  # Re-using the structure but making it cleaner.

  p <- ggplot(calc_map, aes(x = Position, y = Map)) +
    geom_line(linewidth = 2, color = "gray80") +
    geom_point(size = 5, color = "#2c3e50") +
    geom_text(aes(label = Label), vjust = -1.5, size = 4, fontface = "bold") +
    scale_x_continuous(limits = c(-5, max(calc_map$Position) + 10)) +
    labs(
      title = "Genetic Map of Chromosome 2",
      subtitle = "Relative positions based on recombination frequencies (cM)",
      x = "Map Distance (cM)",
      y = ""
    ) +
    theme_publication() +
    theme(
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.y = element_blank(), # Hide y axis text as we only have one map
      panel.grid.major.y = element_blank()
    )

  return(p)
}

#' Calculate Interference
#'
#' Calculates the Coefficient of Coincidence (CoC) and Interference (I).
#'
#' @param map_results The list output from build_triple_map().
#' @return A dataframe summary of the interference analysis.
#' @export
calculate_interference <- function(map_results) {
  N <- map_results$N

  # Get recombination fractions (decimal)
  rf_EW <- map_results$distances["EW"] / 100
  rf_WB <- map_results$distances["WB"] / 100

  # Observed DCOs
  obs_DCO_count <- sum(map_results$data$count[map_results$data$class == "DCO"])
  obs_DCO_freq <- obs_DCO_count / N

  # Expected DCOs (product rule)
  exp_DCO_freq <- rf_EW * rf_WB
  exp_DCO_count <- exp_DCO_freq * N

  # Metrics
  CoC <- obs_DCO_freq / exp_DCO_freq
  Interference <- 1 - CoC

  # Create summary table
  tibble(
    Metric = c("Observed DCOs", "Expected DCOs", "Coeff. of Coincidence (CoC)", "Interference (I)"),
    Value = c(obs_DCO_count, round(exp_DCO_count, 1), round(CoC, 3), round(Interference, 3)),
    Description = c(
      "Actual double crossovers found",
      "Predicted if intervals were independent",
      "Observed / Expected",
      "1 - CoC (Positive means suppression)"
    )
  )
}

#' Calculate Map Confidence Intervals

#'
#' Calculates 95% Confidence Intervals for map distances using the standard error of a proportion.
#'
#' @param map_results The list output from build_triple_map().
#' @return A dataframe with distances and CIs.
#' @export
calculate_map_confidence_intervals <- function(map_results) {
  require(dplyr)

  N <- map_results$N

  # Extract distances (in cM)
  dist_EW <- map_results$distances["EW"]
  dist_WB <- map_results$distances["WB"]
  dist_Total <- map_results$distances["Total"]

  # Convert to proportions (p)
  p_EW <- dist_EW / 100
  p_WB <- dist_WB / 100
  p_Total <- dist_Total / 100

  # Calculate Standard Error (SE) = sqrt(p(1-p)/N)
  se_EW <- sqrt((p_EW * (1 - p_EW)) / N)
  se_WB <- sqrt((p_WB * (1 - p_WB)) / N)
  se_Total <- sqrt((p_Total * (1 - p_Total)) / N)

  # Calculate 95% CI = p +/- 1.96 * SE
  # Convert back to cM (* 100)
  create_ci_row <- function(name, dist, se) {
    lower <- (dist / 100 - 1.96 * se) * 100
    upper <- (dist / 100 + 1.96 * se) * 100

    tibble(
      Interval = name,
      `Distance (cM)` = round(dist, 1),
      `Standard Error` = round(se * 100, 2), # SE in cM
      `95% CI Lower` = round(lower, 1),
      `95% CI Upper` = round(upper, 1),
      `Formatted CI` = sprintf("%.1f (%.1f–%.1f)", dist, lower, upper)
    )
  }

  bind_rows(
    create_ci_row("Crimson – No Take Off", dist_EW, se_EW),
    create_ci_row("No Take Off – Spray Tan", dist_WB, se_WB),
    create_ci_row("Crimson – Spray Tan (Total)", dist_Total, se_Total)
  )
}

#' Plot Map Comparison Dumbbell
#'
#' Creates a dumbbell plot comparing Calculated vs Standard map distances.
#'
#' @param map_results The list output from build_triple_map().
#' @return A ggplot object.
#' @export
plot_map_comparison_dumbbell <- function(map_results) {
  require(ggplot2)
  require(dplyr)
  require(tidyr)

  # 1. Get Calculated Distances
  calc_EW <- map_results$distances["EW"]
  calc_WB <- map_results$distances["WB"]

  # 2. Define Standard Distances (from literature/text)
  # Brown (104.5) - Vestigial (67.0) = 37.5 cM (Corresponds to Crimson-NTO)
  # Vestigial (67.0) - Black (48.5) = 18.5 cM (Corresponds to NTO-Spray Tan)
  std_EW <- 37.5
  std_WB <- 18.5

  # 3. Create Dataframe
  plot_data <- tibble(
    Interval = c("Crimson – No Take Off", "No Take Off – Spray Tan"),
    Calculated = c(calc_EW, calc_WB),
    Standard = c(std_EW, std_WB)
  ) %>%
    pivot_longer(
      cols = c("Calculated", "Standard"),
      names_to = "Source",
      values_to = "Distance"
    )

  # 4. Create Dumbbell Plot
  ggplot(plot_data, aes(x = Distance, y = Interval)) +
    # Line connecting the dots
    geom_line(aes(group = Interval), color = "grey60", linewidth = 1.5) +
    # Points
    geom_point(aes(color = Source), size = 5) +
    # Colors
    scale_color_manual(values = c("Calculated" = "#e74c3c", "Standard" = "#2c3e50")) +
    # Labels
    geom_text(aes(label = round(Distance, 1)), vjust = -1.5, size = 3.5) +
    # Scales
    scale_x_continuous(limits = c(0, 50)) +
    # Titles
    labs(
      title = "Map Distance Comparison",
      subtitle = "Experimental Results vs. Standard Literature Values",
      x = "Map Distance (cM)",
      y = NULL,
      color = "Data Source"
    ) +
    theme_publication() +
    theme(
      panel.grid.major.y = element_blank(),
      legend.position = "top"
    )
}
