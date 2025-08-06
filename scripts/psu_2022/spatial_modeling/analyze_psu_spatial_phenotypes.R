#' Spatial Phenotype Analysis for PSU Field Experiment
#'
#' This script implements spherical spatial correlation modeling for field
#' experiment analysis, representing a core methodological innovation for
#' detecting consistent inv4m effects across environments while accounting
#' for spatial field variation.
#'
#' @author Francisco Rodriguez
#' @date 2025-08-04

# Load required libraries -----------------------------------------------
library(dplyr)
library(ggplot2)
library(nlme)          # Mixed-effects models with spatial correlation
library(tidyr)         # Data reshaping
library(rstatix)       # Statistical tests
library(FactoMineR)    # PCA analysis
library(factoextra)    # PCA visualization
library(emmeans)       # Estimated marginal means
library(ggtext)        # Enhanced text formatting
library(ggfx)          # Enhanced graphics effects
library(ggbeeswarm)    # Beeswarm plots
library(ggpubr)        # Publication-ready plots
library(growthcurver)  # Growth curve modeling

# Configuration ----------------------------------------------------------
DATA_DIR <- "../../../data"
OUTPUT_DIR <- "results_psu_2022"

# Color palette
PALETTE <- c("gold", "#4a0f82")

# Mineral analysis configuration
MINERAL_NAMES <- c(
  "Al", "B", "Ca", "Fe", "K", "Mg", "Mn", "Mo",
  "Ni", "P", "S", "Zn", "Na", "Cu"
)

RELIABLE_MINERALS <- c("Ca", "Fe", "K", "Mg", "Mn", "P", "S", "Zn")
UNRELIABLE_MINERALS <- c("Al", "B", "Mo", "Ni", "Na", "Cu")

# Helper functions -------------------------------------------------------

#' Validate input file paths
#'
#' @param file_paths Named vector of file paths to validate
validate_file_paths <- function(file_paths) {
  missing_files <- file_paths[!file.exists(file_paths)]
  if (length(missing_files) > 0) {
    stop("Missing files: ", paste(names(missing_files), collapse = ", "))
  }
}

#' Load and prepare ionome data
#'
#' @param mineral_csv Path to mineral data CSV
#' @return Processed mineral data frame
load_ionome_data <- function(mineral_csv) {
  validate_file_paths(c("mineral" = mineral_csv))
  
  raw_data <- read.csv(mineral_csv, na.strings = c("", "#N/A", "NA"))
  
  # Standardize tissue naming
  raw_data$tissue[raw_data$tissue == "stalk"] <- "stover"
  
  return(raw_data)
}

#' Load and prepare field phenotype data
#'
#' @param plant_csv Path to plant phenotype CSV
#' @param ear_csv Path to ear phenotype CSV
#' @return Processed phenotype data frame
load_field_phenotypes <- function(plant_csv, ear_csv) {
  file_paths <- c("plant" = plant_csv, "ear" = ear_csv)
  validate_file_paths(file_paths)
  
  # Load ear data
  psu_ear <- read.csv(ear_csv, na.strings = c("", "n/a", "NA"), skip = 1) %>%
    select(-description, -RK, -CC, -NIR, ear_rep = rep) %>%
    rename(rowid = row) %>%
    arrange(rowid) %>%
    group_by(rowid) %>%
    select(-ear_rep) %>%
    summarise_all(mean, na.rm = TRUE) %>%
    droplevels()
  
  ear_cols <- colnames(psu_ear)[-1]
  
  # Load plant data
  psu_plant <- read.csv(plant_csv) %>%
    rename(Genotype = Who.What, rowid = P22.) %>%
    filter(rowid >= 3004, rowid <= 4192) %>%  # Inv4m experiment rows
    mutate(Treatment = factor(Treatment, levels = c("HighP", "LowP"))) %>%
    filter(Genotype %in% c("CTRL", "INV4M")) %>%
    mutate(Genotype = factor(Genotype, levels = c("CTRL", "INV4M"))) %>%
    droplevels()
  
  # Standardize factor levels
  psu_plant$Rep <- as.factor(psu_plant$Rep)
  levels(psu_plant$Genotype) <- c("CTRL", "Inv4m")
  levels(psu_plant$Treatment) <- c("+P", "-P")
  
  # Standardize column names
  short_names <- c("PH", "STW40", "STW50", "STW60", "STWHV")
  long_names <- c(
    "Height_Anthesis", "X40_DAP_dw", "X50_DAP_dw", 
    "X60_DAP_dw", "harvest_dw"
  )
  colnames(psu_plant)[colnames(psu_plant) %in% long_names] <- short_names
  
  return(list(plant = psu_plant, ear = psu_ear, ear_cols = ear_cols))
}

#' Combine phenotype and ionome data
#'
#' @param plant_data Plant phenotype data
#' @param ear_data Ear phenotype data  
#' @param ionome_data Ionome data
#' @return Combined phenotype data frame
combine_phenotype_data <- function(plant_data, ear_data, ionome_data) {
  # Combine phenotype data
  pheno <- plant_data %>%
    select(
      rowid, Plot, Rep, Plot_Row, Plot_Column, 
      Treatment, Genotype, PH:DTS, STW40:STWHV
    ) %>%
    left_join(ear_data, by = "rowid") %>%
    mutate(HI = TKW / STWHV) %>%
    mutate(
      Plot_Column = case_when(
        Plot == "PlotVIII" ~ Plot_Column + 10,
        Plot == "PlotVI" ~ Plot_Column,
        TRUE ~ NA_real_
      )
    ) %>%
    filter(!is.na(PH)) %>%
    droplevels()
  
  # Add ionome data
  if (!is.null(ionome_data)) {
    ionome_wide <- ionome_data %>%
      select(rowid, tissue, all_of(RELIABLE_MINERALS)) %>%
      pivot_wider(
        names_from = tissue,
        values_from = all_of(RELIABLE_MINERALS)
      )
    
    pheno <- pheno %>%
      left_join(ionome_wide, by = "rowid")
  }
  
  return(pheno)
}

#' Fit spatial mixed-effects model for single trait
#'
#' @param data Data frame containing the phenotype data
#' @param trait_name Name of the trait column to analyze
#' @param formula_rhs Right-hand side of model formula (default: "Genotype*Treatment")
#' @return Fitted lme model object
fit_spatial_model <- function(data, trait_name, 
                             formula_rhs = "Genotype*Treatment") {
  if (!trait_name %in% colnames(data)) {
    stop("Trait '", trait_name, "' not found in data")
  }
  
  if (sum(!is.na(data[[trait_name]])) < 10) {
    warning("Insufficient data for trait: ", trait_name)
    return(NULL)
  }
  
  # Construct formula
  formula_str <- paste(trait_name, "~", formula_rhs)
  
  tryCatch({
    model <- lme(
      fixed = as.formula(formula_str),
      random = ~ 1 | Rep,
      correlation = corSpher(form = ~ Plot_Row + Plot_Column),
      na.action = na.exclude,
      data = data
    )
    
    return(model)
  }, error = function(e) {
    warning("Model fitting failed for ", trait_name, ": ", e$message)
    return(NULL)
  })
}

#' Extract effects from spatial model
#'
#' @param model Fitted lme model
#' @param trait_name Name of the trait
#' @return Data frame with model coefficients and statistics
extract_spatial_effects <- function(model, trait_name) {
  if (is.null(model)) {
    return(NULL)
  }
  
  effects <- coef(summary(model)) %>%
    as.data.frame() %>%
    tibble::rownames_to_column("predictor")
  
  colnames(effects)[3] <- "Std.Error"
  colnames(effects)[6] <- "p.value"
  
  effects$response <- trait_name
  effects <- effects %>%
    select(response, everything()) %>%
    filter(!grepl("Intercept", predictor))
  
  return(effects)
}

#' Analyze multiple traits with spatial models
#'
#' @param data Combined phenotype data
#' @param trait_list Vector of trait names to analyze
#' @param formula_rhs Right-hand side of model formula
#' @return List containing effects and models
analyze_traits_spatial <- function(data, trait_list, 
                                  formula_rhs = "Genotype*Treatment") {
  effects_list <- list()
  models_list <- list()
  
  for (trait in trait_list) {
    cat("Analyzing trait:", trait, "\n")
    
    model <- fit_spatial_model(data, trait, formula_rhs)
    
    if (!is.null(model)) {
      effects <- extract_spatial_effects(model, trait)
      
      if (!is.null(effects)) {
        effects_list[[trait]] <- effects
        models_list[[trait]] <- model
      }
    }
  }
  
  # Combine effects
  combined_effects <- bind_rows(effects_list) %>%
    mutate(p.adjust = p.adjust(p.value, method = "fdr")) %>%
    arrange(p.adjust)
  
  return(list(
    effects = combined_effects,
    models = models_list
  ))
}

#' Create standardized effect sizes
#'
#' @param data Original data
#' @param trait_list Vector of trait names
#' @param significant_effects Significant effects from unscaled analysis
#' @return Data frame with standardized effects
calculate_standardized_effects <- function(data, trait_list, significant_effects) {
  # Scale the data
  scaled_data <- data %>%
    select(rowid:Genotype) %>%
    bind_cols(
      data %>%
        select(all_of(trait_list)) %>%
        scale() %>%
        as.data.frame()
    )
  
  # Re-analyze with scaled data
  scaled_results <- analyze_traits_spatial(scaled_data, trait_list)
  
  # Filter to significant terms from unscaled analysis
  significant_scaled <- scaled_results$effects %>%
    inner_join(
      significant_effects %>% select(response, predictor),
      by = c("response", "predictor")
    ) %>%
    mutate(
      CI.upper = Value + qnorm(0.975) * Std.Error,
      CI.lower = Value - qnorm(0.975) * Std.Error
    )
  
  return(significant_scaled)
}

#' Generate publication-ready effect plot
#'
#' @param effects_data Standardized effects data
#' @param trait_labels Named vector of trait labels for plotting
#' @return ggplot object
create_effects_plot <- function(effects_data, trait_labels) {
  # Prepare data for plotting
  plot_data <- effects_data %>%
    mutate(
      predictor = factor(
        predictor, 
        levels = c("Treatment-P", "GenotypeInv4m", "GenotypeInv4m:Treatment-P")
      )
    ) %>%
    group_by(response) %>%
    mutate(max_effect = max(abs(Value))) %>%
    ungroup() %>%
    arrange(max_effect) %>%
    mutate(response = forcats::fct_inorder(response))
  
  # Update predictor labels
  levels(plot_data$predictor) <- c("-P", "Inv4m", "Inv4m:-P")
  
  # Add trait group classification
  plot_data$trait_group <- if_else(
    grepl("_seed|_stover", plot_data$response, perl = TRUE),
    "Ionome", "Direct"
  )
  
  # Create trait labels
  plot_data$trait_label <- factor(trait_labels[as.character(plot_data$response)])
  plot_data$trait_label <- forcats::fct_reorder(plot_data$trait_label, plot_data$max_effect)
  
  # Create the plot
  effects_plot <- plot_data %>%
    mutate(predictor = gsub("Inv4m", "<i>Inv4m</i>", predictor)) %>%
    mutate(predictor = gsub("<i>Inv4m</i>:", "<i>Inv4m</i> × ", predictor)) %>%
    ggplot(aes(x = Value, y = trait_label, group = trait_group)) +
    xlab("Standardized Effect") +
    ylab("Trait") +
    geom_vline(xintercept = 0, lty = 2) +
    geom_point(position = position_dodge(0.4), size = 4) +
    geom_errorbar(
      aes(xmin = CI.lower, xmax = CI.upper),
      position = position_dodge(0.4),
      width = 0.2,
      size = 0.7
    ) +
    facet_wrap(~ predictor, ncol = 3) +
    guides(shape = guide_legend(title = NULL, reverse = TRUE)) +
    theme_classic2(base_size = 25) +
    theme(
      strip.background = element_blank(),
      strip.text = element_markdown(face = "bold", size = 25),
      axis.title.y = element_blank(),
      axis.text.y = element_markdown(hjust = 1, color = "black"),
      panel.grid.major.y = element_line(color = "grey90"),
      plot.caption = element_text(hjust = 0)
    )
  
  return(effects_plot)
}

#' Perform t-tests and add statistical annotations
#'
#' @param data Phenotype data
#' @param significant_traits Vector of significant trait names
#' @return List of t-test results for different comparisons
perform_trait_tests <- function(data, significant_traits) {
  # Reshape data for testing
  test_data <- data %>%
    select(rowid, Genotype, Treatment, all_of(significant_traits)) %>%
    pivot_longer(
      cols = all_of(significant_traits),
      names_to = "trait",
      values_to = "value"
    ) %>%
    group_by(Treatment) %>%
    filter(trait %in% significant_traits)
  
  # Genotype effect tests (within each treatment)
  genotype_tests <- test_data %>%
    group_by(Treatment, trait) %>%
    t_test(value ~ Genotype, ref.group = "CTRL") %>%
    add_significance() %>%
    add_x_position(x = "Genotype") %>%
    add_y_position(scales = "free") %>%
    mutate(group2 = gsub("Inv4m", "<i>Inv4m</i>", group2))
  
  # Treatment effect tests (within each genotype)
  treatment_tests <- test_data %>%
    group_by(Genotype, trait) %>%
    t_test(value ~ Treatment, ref.group = "+P") %>%
    add_significance() %>%
    add_x_position(x = "Treatment") %>%
    add_y_position(scales = "free") %>%
    mutate(y.position = 1.05 * y.position) %>%
    mutate(Genotype = gsub("Inv4m", "<i>Inv4m</i>", Genotype))
  
  return(list(
    genotype = genotype_tests,
    treatment = treatment_tests,
    test_data = test_data
  ))
}

#' Create individual trait plots with statistical annotations
#'
#' @param data Phenotype data
#' @param trait_name Name of the trait to plot
#' @param test_results Results from perform_trait_tests
#' @param plot_title Title for the plot
#' @param y_label Y-axis label
#' @param y_limits Y-axis limits
#' @return ggplot object
create_trait_plot <- function(data, trait_name, test_results, 
                             plot_title, y_label, y_limits = NULL) {
  # Filter relevant test results
  genotype_test <- test_results$genotype %>%
    filter(trait == trait_name)
  
  treatment_test <- test_results$treatment %>%
    filter(trait == trait_name)
  
  # Determine which test to use based on significant effects
  use_genotype_test <- nrow(genotype_test) > 0
  use_treatment_test <- nrow(treatment_test) > 0
  
  # Base plot
  plot_data <- data %>%
    mutate(Genotype = gsub("Inv4m", "<i>Inv4m</i>", Genotype)) %>%
    mutate(Genotype = factor(Genotype, levels = c("CTRL", "<i>Inv4m</i>")))
  
  if (use_genotype_test) {
    # Genotype effect plot
    p <- plot_data %>%
      ggplot(aes(x = Genotype, y = .data[[trait_name]], col = Genotype)) +
      ggtitle(plot_title) +
      ylab(y_label) +
      geom_boxplot(width = 0.25, linewidth = 1.5, alpha = 0) %>%
      with_shadow(colour = "black", x_offset = 0, y_offset = 0, sigma = 1) +
      geom_quasirandom(size = 2) %>%
      with_shadow(colour = "black", x_offset = 0, y_offset = 0, sigma = 1) +
      stat_pvalue_manual(
        genotype_test %>% mutate(y.position = y.position + max(plot_data[[trait_name]], na.rm = TRUE) * 0.02),
        size = 10,
        bracket.size = 0.8,
        hide.ns = TRUE
      ) +
      scale_color_manual(values = PALETTE) +
      facet_wrap(~ Treatment) +
      theme_classic2(base_size = 20) +
      theme(
        plot.title = element_markdown(hjust = 0.5, face = "bold"),
        axis.title.x = element_blank(),
        axis.text.x = element_markdown(face = "bold", color = "black", size = 20, angle = 45, hjust = 1),
        strip.background = element_blank(),
        strip.text = element_text(size = 25),
        legend.position = "none"
      )
  } else if (use_treatment_test) {
    # Treatment effect plot
    p <- plot_data %>%
      ggplot(aes(x = Treatment, y = .data[[trait_name]], col = Treatment)) +
      ggtitle(plot_title) +
      ylab(y_label) +
      geom_boxplot(width = 0.25, linewidth = 1, alpha = 0) %>%
      with_shadow(colour = "black", x_offset = 0, y_offset = 0, sigma = 1) +
      geom_quasirandom(size = 2) %>%
      with_shadow(colour = "black", x_offset = 0, y_offset = 0, sigma = 1) +
      stat_pvalue_manual(
        treatment_test,
        size = 10,
        bracket.size = 0.8,
        hide.ns = TRUE
      ) +
      scale_color_manual(values = PALETTE) +
      facet_wrap(~ factor(Genotype, levels = c("CTRL", "<i>Inv4m</i>"))) +
      theme_classic2(base_size = 20) +
      theme(
        plot.title = element_markdown(hjust = 0.5, face = "bold"),
        axis.title.y = element_markdown(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 25, face = "bold", color = "black"),
        strip.background = element_blank(),
        strip.text = element_markdown(size = 20),
        legend.position = "none"
      )
  } else {
    # Basic plot without statistical tests
    p <- plot_data %>%
      ggplot(aes(x = Genotype, y = .data[[trait_name]], col = Genotype)) +
      ggtitle(plot_title) +
      ylab(y_label) +
      geom_boxplot(width = 0.25, linewidth = 1.5, alpha = 0) +
      geom_quasirandom(size = 2) +
      scale_color_manual(values = PALETTE) +
      facet_wrap(~ Treatment) +
      theme_classic2(base_size = 20) +
      theme(
        plot.title = element_markdown(hjust = 0.5, face = "bold"),
        axis.title.x = element_blank(),
        axis.text.x = element_markdown(face = "bold", color = "black", size = 20),
        strip.background = element_blank(),
        strip.text = element_text(size = 25),
        legend.position = "none"
      )
  }
  
  if (!is.null(y_limits)) {
    p <- p + coord_cartesian(ylim = y_limits)
  }
  
  return(p)
}

# Main analysis pipeline -------------------------------------------------

#' Main function to run spatial phenotype analysis
#'
#' @export
run_spatial_analysis <- function() {
  cat("Starting spatial phenotype analysis...\n")
  
  # File paths
  mineral_csv <- file.path(DATA_DIR, "PSU_inv4m_ionome_all.csv")
  plant_csv <- file.path(DATA_DIR, "22_NCS_PSU_LANGEBIO_FIELDS_PSU_P_field.csv")
  ear_csv <- file.path(DATA_DIR, "22_NCS_PSU_LANGEBIO_FIELDS_PSU_P_field_ear_pheno.csv")
  
  # 1. Load data
  cat("Loading phenotype and ionome data...\n")
  ionome_data <- load_ionome_data(mineral_csv)
  field_data <- load_field_phenotypes(plant_csv, ear_csv)
  
  # 2. Combine data
  cat("Combining phenotype data...\n")
  pheno <- combine_phenotype_data(
    field_data$plant,
    field_data$ear,
    ionome_data
  )
  
  cat("Final dataset contains", nrow(pheno), "observations\n")
  
  # 3. Define traits for analysis
  trait_vars <- colnames(pheno %>% select(PH:last_col(), -contains("_")))
  ionome_vars <- colnames(pheno %>% select(contains("_seed"), contains("_stover")))
  all_vars <- c(trait_vars, ionome_vars)
  
  cat("Analyzing", length(all_vars), "traits\n")
  
  # 4. Spatial analysis
  cat("Performing spatial mixed-effects modeling...\n")
  spatial_results <- analyze_traits_spatial(pheno, all_vars)
  
  # 5. Identify significant effects
  significant_effects <- spatial_results$effects %>%
    filter(p.adjust < 0.05, !grepl("Intercept", predictor))
  
  cat("Found", nrow(significant_effects), "significant effects\n")
  
  # 6. Calculate standardized effects
  cat("Calculating standardized effect sizes...\n")
  standardized_effects <- calculate_standardized_effects(
    pheno, all_vars, significant_effects
  )
  
  # 7. Create trait labels for plotting
  trait_labels <- c(
    Ca_seed = "Seed Calcium",
    Mg_seed = "Seed Magnesium", 
    S_stover = "Stover Sulphur",
    Ca_stover = "Stover Calcium",
    Zn_stover = "Stover Zinc",
    P_seed = "Seed Phosphorus",
    P_stover = "Stover Phosphorus",
    KW50 = "50 Kernel Weight",
    CD = "<b>Cob Diameter</b>",
    PH = "<b>Plant Height at Anthesis</b>",
    DTA = "<b>Days to Anthesis</b>",
    STWHV = "Stover Weight at Harvest",
    DTS = "<b>Days to Silking</b>",
    STW40 = "Stover Weight at 40 DAP",
    STW60 = "Stover Weight at 60 DAP",
    STW50 = "Stover Weight at 50 DAP"
  )
  
  # 8. Generate main effects plot
  cat("Creating effects visualization...\n")
  effects_plot <- create_effects_plot(standardized_effects, trait_labels)
  
  # 9. Statistical testing for individual traits
  significant_traits <- unique(significant_effects$response)
  test_results <- perform_trait_tests(pheno, significant_traits)
  
  # 10. Save results
  cat("Saving results...\n")
  write.csv(
    spatial_results$effects,
    file.path(OUTPUT_DIR, "phenotypic_effects.csv"),
    row.names = FALSE
  )
  
  write.csv(
    standardized_effects,
    file.path(OUTPUT_DIR, "standardized_phenotypic_effects.csv"),
    row.names = FALSE
  )
  
  # Save main plot
  ggsave(
    file.path(OUTPUT_DIR, "spatial_effects_plot.pdf"),
    effects_plot,
    width = 12, height = 7
  )
  
  cat("Analysis complete!\n")
  cat("Results saved to:\n")
  cat("- phenotypic_effects.csv: Raw spatial model effects\n")
  cat("- standardized_phenotypic_effects.csv: Standardized effect sizes\n")
  cat("- spatial_effects_plot.pdf: Main effects visualization\n")
  
  return(list(
    effects = spatial_results$effects,
    standardized_effects = standardized_effects,
    models = spatial_results$models,
    test_results = test_results,
    phenotype_data = pheno
  ))
}

# Execute analysis if script is run directly
if (!interactive()) {
  results <- run_spatial_analysis()
}
