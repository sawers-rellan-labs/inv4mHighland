#' Seedling Root and Shoot Analysis for Chamber Experiments
#'
#' This script analyzes 3-day seedling root and shoot measurements from
#' controlled chamber experiments, testing inv4m effects across different
#' genetic backgrounds (donors). Part of the developmental axis analysis
#' from seedling to SAM to field stages.
#'
#' @author Francisco Rodriguez
#' @date 2025-08-04

# Load required libraries -----------------------------------------------
library(dplyr)
library(ggplot2)
library(ggbeeswarm)
library(ggpubr)
library(rstatix)
library(tidyr)
library(ggtext)

# Configuration ----------------------------------------------------------
DATA_FILE <- "~/Desktop/inv4m_root_phenotype.csv"
OUTPUT_DIR <- "~/Desktop"

# Color schemes
GENOTYPE_COLORS <- c("CTRL" = "gold", "Inv4m" = "purple4")
DONOR_COLORS <- c("TMEX" = "royalblue", "MI21" = "tomato")

# Helper functions -------------------------------------------------------

#' Validate and load seedling data
#'
#' @param file_path Path to the seedling data CSV file
#' @return Data frame with validated seedling measurements
#' @export
load_seedling_data <- function(file_path) {
  if (!file.exists(file_path)) {
    stop("Data file not found: ", file_path)
  }
  
  # Load data
  data <- read.csv(file_path, stringsAsFactors = FALSE)
  
  # Validate required columns
  required_cols <- c("Root_Length_cm", "Shoot_Length_cm", "Date", "Donor", "Genotype")
  missing_cols <- setdiff(required_cols, colnames(data))
  
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  
  # Standardize factor levels
  data$Donor <- factor(data$Donor, levels = c("TMEX", "MI21"))
  data$Genotype <- factor(data$Genotype, levels = c("CTRL", "Inv4m"))
  data$Date <- as.factor(data$Date)
  
  # Add combined treatment variable for easier analysis
  data$Treatment <- interaction(data$Donor, data$Genotype, sep = "_")
  
  # Validate data ranges
  if (any(data$Root_Length_cm < 0, na.rm = TRUE) || 
      any(data$Shoot_Length_cm < 0, na.rm = TRUE)) {
    warning("Negative length measurements detected")
  }
  
  return(data)
}

#' Perform statistical analysis for a single trait
#'
#' @param data Seedling data frame
#' @param trait_name Name of the trait column to analyze
#' @param formula_rhs Right-hand side of the linear model formula
#' @return List containing model summary and effects
#' @export
analyze_trait <- function(data, trait_name, 
                         formula_rhs = "Date + Donor * Genotype") {
  # Validate inputs
  if (!trait_name %in% colnames(data)) {
    stop("Trait '", trait_name, "' not found in data")
  }
  
  # Check for sufficient data
  valid_obs <- sum(!is.na(data[[trait_name]]))
  if (valid_obs < 10) {
    warning("Insufficient observations for trait: ", trait_name)
    return(NULL)
  }
  
  # Fit linear model
  formula_str <- paste(trait_name, "~", formula_rhs)
  model <- lm(as.formula(formula_str), data = data)
  
  # Extract results
  model_summary <- summary(model)
  
  # Create effects data frame
  effects <- data.frame(
    term = rownames(model_summary$coefficients),
    estimate = model_summary$coefficients[, "Estimate"],
    std_error = model_summary$coefficients[, "Std. Error"],
    t_value = model_summary$coefficients[, "t value"],
    p_value = model_summary$coefficients[, "Pr(>|t|)"],
    trait = trait_name,
    stringsAsFactors = FALSE
  )
  
  return(list(
    model = model,
    summary = model_summary,
    effects = effects,
    r_squared = model_summary$r.squared,
    adj_r_squared = model_summary$adj.r.squared
  ))
}

#' Perform pairwise t-tests for treatment comparisons
#'
#' @param data Seedling data frame
#' @param traits Vector of trait names to analyze
#' @return List of t-test results for different comparison types
#' @export
perform_pairwise_tests <- function(data, traits) {
  # Input validation
  if (!is.character(traits) || length(traits) == 0) {
    stop("traits must be a non-empty character vector")
  }
  
  missing_traits <- setdiff(traits, colnames(data))
  if (length(missing_traits) > 0) {
    warning("Missing traits: ", paste(missing_traits, collapse = ", "))
    traits <- intersect(traits, colnames(data))
  }
  
  results <- list()
  
  for (trait in traits) {
    cat("Performing t-tests for", trait, "...\n")
    
    # Inversion effect within each donor
    inversion_tests <- data %>%
      group_by(Donor) %>%
      t_test(as.formula(paste(trait, "~ Genotype"))) %>%
      add_significance() %>%
      add_xy_position(x = "Genotype") %>%
      mutate(trait = trait, comparison_type = "inversion_by_donor")
    
    # Donor effect within each genotype
    donor_tests <- data %>%
      group_by(Genotype) %>%
      t_test(as.formula(paste(trait, "~ Donor"))) %>%
      add_significance() %>%
      add_xy_position(x = "Donor") %>%
      mutate(trait = trait, comparison_type = "donor_by_genotype")
    
    results[[trait]] <- list(
      inversion = inversion_tests,
      donor = donor_tests
    )
  }
  
  return(results)
}

#' Calculate effect sizes (Cohen's d) for treatment comparisons
#'
#' @param data Seedling data frame
#' @param traits Vector of trait names to analyze
#' @return Data frame with effect size estimates
#' @export
calculate_effect_sizes <- function(data, traits) {
  effect_sizes <- list()
  
  for (trait in traits) {
    # Inversion effect within each donor
    inversion_d <- data %>%
      group_by(Donor) %>%
      cohens_d(as.formula(paste(trait, "~ Genotype"))) %>%
      mutate(trait = trait, comparison_type = "inversion_by_donor")
    
    # Donor effect within each genotype
    donor_d <- data %>%
      group_by(Genotype) %>%
      cohens_d(as.formula(paste(trait, "~ Donor"))) %>%
      mutate(trait = trait, comparison_type = "donor_by_genotype")
    
    effect_sizes[[trait]] <- bind_rows(inversion_d, donor_d)
  }
  
  return(bind_rows(effect_sizes))
}

#' Create standardized plot for a single trait
#'
#' @param data Seedling data frame
#' @param trait_name Name of the trait to plot
#' @param plot_type Type of plot ("inversion" or "donor")
#' @param test_results T-test results for statistical annotations
#' @param y_limits Optional y-axis limits
#' @return ggplot object
#' @export
create_trait_plot <- function(data, trait_name, plot_type, test_results, 
                             y_limits = NULL) {
  # Input validation
  plot_type <- match.arg(plot_type, c("inversion", "donor"))
  
  if (!trait_name %in% colnames(data)) {
    stop("Trait '", trait_name, "' not found in data")
  }
  
  # Create plot title
  plot_title <- switch(trait_name,
    "Root_Length_cm" = ifelse(plot_type == "inversion", 
                             "Root Length: Inversion Effect by Donor",
                             "Root Length: Donor Effect by Genotype"),
    "Shoot_Length_cm" = ifelse(plot_type == "inversion",
                              "Shoot Length: Inversion Effect by Donor", 
                              "Shoot Length: Donor Effect by Genotype"),
    paste(trait_name, "Analysis")
  )
  
  # Create y-axis label
  y_label <- switch(trait_name,
    "Root_Length_cm" = "Root Length (cm)",
    "Shoot_Length_cm" = "Shoot Length (cm)",
    trait_name
  )
  
  # Set up plot aesthetics based on plot type
  if (plot_type == "inversion") {
    # Inversion effect plot
    p <- ggplot(data, aes(x = Genotype, y = .data[[trait_name]], color = Genotype)) +
      facet_wrap(~ Donor, scales = "free_x") +
      scale_color_manual(values = GENOTYPE_COLORS) +
      stat_pvalue_manual(
        test_results,
        label = "p.signif",
        tip.length = 0.02,
        size = 4
      )
  } else {
    # Donor effect plot
    p <- ggplot(data, aes(x = Donor, y = .data[[trait_name]], color = Donor)) +
      facet_wrap(~ Genotype, scales = "free_x") +
      scale_color_manual(values = DONOR_COLORS) +
      stat_pvalue_manual(
        test_results,
        label = "p.signif",
        tip.length = 0.02,
        size = 4
      )
  }
  
  # Add common plot elements
  p <- p +
    geom_boxplot(
      fill = "white",
      outlier.shape = NA,
      width = 0.5
    ) +
    geom_quasirandom(
      size = 2.5,
      alpha = 0.8,
      width = 0.3
    ) +
    labs(
      title = plot_title,
      y = y_label,
      color = ifelse(plot_type == "inversion", "Genotype", "Donor")
    ) +
    theme_classic(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      legend.position = "right",
      strip.background = element_blank(),
      strip.text = element_text(face = "bold")
    )
  
  # Apply y-axis limits if provided
  if (!is.null(y_limits)) {
    p <- p + coord_cartesian(ylim = y_limits)
  }
  
  return(p)
}

#' Create multi-panel figure with all trait comparisons
#'
#' @param data Seedling data frame
#' @param traits Vector of trait names
#' @param test_results T-test results
#' @return ggplot object
#' @export
create_multipanel_figure <- function(data, traits, test_results) {
  plots <- list()
  
  for (trait in traits) {
    # Create inversion effect plots
    inversion_plot <- create_trait_plot(
      data, trait, "inversion", 
      test_results[[trait]]$inversion
    ) + theme(legend.position = "none", plot.title = element_blank())
    
    # Create donor effect plots  
    donor_plot <- create_trait_plot(
      data, trait, "donor",
      test_results[[trait]]$donor
    ) + theme(legend.position = "none", plot.title = element_blank())
    
    plots[[paste0(trait, "_inversion")]] <- inversion_plot
    plots[[paste0(trait, "_donor")]] <- donor_plot
  }
  
  # Extract legends
  legend_genotype <- get_legend(
    create_trait_plot(data, traits[1], "inversion", test_results[[traits[1]]]$inversion) + 
      theme(legend.position = "bottom")
  )
  
  legend_donor <- get_legend(
    create_trait_plot(data, traits[1], "donor", test_results[[traits[1]]]$donor) + 
      theme(legend.position = "bottom")
  )
  
  # Arrange plots
  if (length(traits) == 2) {
    # Standard layout for root and shoot
    left_panels <- ggarrange(
      plots[[paste0(traits[2], "_inversion")]],  # Shoot inversion
      plots[[paste0(traits[1], "_inversion")]],  # Root inversion
      ncol = 1,
      labels = c("A", "C"),
      label.x = 0.02,
      label.y = 0.98,
      font.label = list(size = 16, face = "bold")
    )
    
    right_panels <- ggarrange(
      plots[[paste0(traits[2], "_donor")]],      # Shoot donor
      plots[[paste0(traits[1], "_donor")]],      # Root donor
      ncol = 1,
      labels = c("B", "D"),
      label.x = 0.02,
      label.y = 0.98,
      font.label = list(size = 16, face = "bold")
    )
    
    main_figure <- ggarrange(
      left_panels,
      right_panels,
      ncol = 2,
      widths = c(1, 1)
    )
  } else {
    # Flexible layout for other numbers of traits
    main_figure <- ggarrange(
      plotlist = plots,
      ncol = 2,
      labels = LETTERS[1:length(plots)],
      font.label = list(size = 16, face = "bold")
    )
  }
  
  # Add legends
  final_figure <- ggarrange(
    main_figure,
    ggarrange(legend_genotype, legend_donor, ncol = 2),
    nrow = 2,
    heights = c(1, 0.1)
  )
  
  return(final_figure)
}

#' Generate descriptive statistics summary
#'
#' @param data Seedling data frame
#' @param traits Vector of trait names
#' @return Data frame with summary statistics
#' @export
generate_summary_stats <- function(data, traits) {
  summary_stats <- data %>%
    select(Genotype, Donor, all_of(traits)) %>%
    pivot_longer(
      cols = all_of(traits),
      names_to = "trait",
      values_to = "value"
    ) %>%
    group_by(trait, Genotype, Donor) %>%
    summarise(
      n = sum(!is.na(value)),
      mean = mean(value, na.rm = TRUE),
      sd = sd(value, na.rm = TRUE),
      se = sd / sqrt(n),
      median = median(value, na.rm = TRUE),
      min = min(value, na.rm = TRUE),
      max = max(value, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(trait, Genotype, Donor)
  
  return(summary_stats)
}

# Main analysis pipeline -------------------------------------------------

#' Main function to run complete seedling analysis
#'
#' @export
run_seedling_analysis <- function() {
  cat("Starting seedling root and shoot analysis...\n")
  
  # 1. Load and validate data
  cat("Loading seedling data...\n")
  data_with_factors <- load_seedling_data(DATA_FILE)
  
  cat("Loaded", nrow(data_with_factors), "observations\n")
  
  # 2. Define traits for analysis
  traits <- c("Root_Length_cm", "Shoot_Length_cm")
  
  # 3. Perform linear model analysis
  cat("Performing linear model analysis...\n")
  trait_models <- list()
  
  for (trait in traits) {
    cat("Analyzing", trait, "...\n")
    trait_models[[trait]] <- analyze_trait(data_with_factors, trait)
    
    if (!is.null(trait_models[[trait]])) {
      cat("  R-squared:", round(trait_models[[trait]]$r_squared, 3), "\n")
    }
  }
  
  # 4. Perform pairwise t-tests
  cat("Performing pairwise t-tests...\n")
  test_results <- perform_pairwise_tests(data_with_factors, traits)
  
  # 5. Calculate effect sizes
  cat("Calculating effect sizes...\n")
  effect_sizes <- calculate_effect_sizes(data_with_factors, traits)
  
  # 6. Generate summary statistics
  cat("Generating summary statistics...\n")
  summary_stats <- generate_summary_stats(data_with_factors, traits)
  
  # 7. Create visualizations
  cat("Creating visualizations...\n")
  
  # Individual plots
  individual_plots <- list()
  for (trait in traits) {
    individual_plots[[paste0(trait, "_inversion")]] <- create_trait_plot(
      data_with_factors, trait, "inversion", 
      test_results[[trait]]$inversion
    )
    
    individual_plots[[paste0(trait, "_donor")]] <- create_trait_plot(
      data_with_factors, trait, "donor",
      test_results[[trait]]$donor
    )
  }
  
  # Multi-panel figure
  multipanel_figure <- create_multipanel_figure(
    data_with_factors, traits, test_results
  )
  
  # 8. Save results
  cat("Saving results...\n")
  
  # Save summary statistics
  write.csv(
    summary_stats,
    file.path(OUTPUT_DIR, "seedling_summary_statistics.csv"),
    row.names = FALSE
  )
  
  # Save effect sizes
  write.csv(
    effect_sizes,
    file.path(OUTPUT_DIR, "seedling_effect_sizes.csv"),
    row.names = FALSE
  )
  
  # Save t-test results
  all_tests <- bind_rows(
    lapply(names(test_results), function(trait) {
      bind_rows(
        test_results[[trait]]$inversion,
        test_results[[trait]]$donor
      )
    })
  )
  
  write.csv(
    all_tests,
    file.path(OUTPUT_DIR, "seedling_statistical_tests.csv"),
    row.names = FALSE
  )
  
  # Save model summaries
  model_effects <- bind_rows(lapply(trait_models, function(x) {
    if (!is.null(x)) x$effects else NULL
  }))
  
  write.csv(
    model_effects,
    file.path(OUTPUT_DIR, "seedling_model_effects.csv"),
    row.names = FALSE
  )
  
  # Save plots
  ggsave(
    file.path(OUTPUT_DIR, "seedling_analysis_multipanel.pdf"),
    multipanel_figure,
    width = 12, height = 10, dpi = 300
  )
  
  ggsave(
    file.path(OUTPUT_DIR, "seedling_analysis_multipanel.png"),
    multipanel_figure,
    width = 12, height = 10, dpi = 300
  )
  
  cat("Analysis complete!\n")
  cat("Results saved to:\n")
  cat("- seedling_summary_statistics.csv: Descriptive statistics by group\n")
  cat("- seedling_effect_sizes.csv: Cohen's d effect sizes\n")
  cat("- seedling_statistical_tests.csv: T-test results\n")
  cat("- seedling_model_effects.csv: Linear model coefficients\n")
  cat("- seedling_analysis_multipanel.pdf/png: Multi-panel figure\n")
  
  return(list(
    data = data_with_factors,
    models = trait_models,
    tests = test_results,
    effect_sizes = effect_sizes,
    summary_stats = summary_stats,
    plots = individual_plots,
    multipanel = multipanel_figure
  ))
}

# Execute analysis if script is run directly
if (!interactive()) {
  results <- run_seedling_analysis()
}
