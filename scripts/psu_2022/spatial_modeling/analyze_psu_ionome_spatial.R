#' Spatial Ionome Analysis for PSU Field Experiment
#'
#' This script analyzes ionome data from the PSU field experiment using
#' spatial mixed-effects models. It applies spherical spatial correlation
#' modeling to detect consistent inv4m effects on mineral concentrations
#' in seed and stover tissues.
#'
#' @author Francisco Rodriguez
#' @date 2025-08-04

# Load required libraries -----------------------------------------------
library(dplyr)
library(janitor)
library(ggplot2)
library(FactoMineR)
library(factoextra)
library(nlme)         # Mixed-effects models with spatial correlation
library(tidyr)        # Data reshaping
library(rstatix)      # Statistical tests
library(ggbeeswarm)   # Beeswarm plots
library(ggpubr)       # Publication-ready plots
library(ggfx)         # Enhanced graphics effects
library(ggtext)       # Enhanced text formatting
library(tibble)       # Modern data frames
library(forcats)      # Factor manipulation

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

#' Load and process ionome data
#'
#' @param mineral_csv Path to ionome data CSV file
#' @return Processed ionome data frame
#' @export
load_ionome_data <- function(mineral_csv) {
  if (!file.exists(mineral_csv)) {
    stop("Ionome data file not found: ", mineral_csv)
  }
  
  raw_data <- read.csv(mineral_csv, na.strings = c("", "#N/A", "NA"))
  
  # Standardize tissue naming
  raw_data$tissue[raw_data$tissue == "stalk"] <- "stover"
  
  # Validate required columns
  required_cols <- c("rowid", "tissue", "rack")
  missing_cols <- setdiff(required_cols, colnames(raw_data))
  
  if (length(missing_cols) > 0) {
    warning("Missing columns in ionome data: ", 
            paste(missing_cols, collapse = ", "))
  }
  
  return(raw_data)
}

#' Load and process field phenotype data
#'
#' @param plant_csv Path to plant phenotype CSV
#' @param ear_csv Path to ear phenotype CSV
#' @return List containing processed phenotype data
#' @export
load_field_phenotypes <- function(plant_csv, ear_csv) {
  # Validate file paths
  files <- c("plant" = plant_csv, "ear" = ear_csv)
  missing_files <- files[!file.exists(files)]
  
  if (length(missing_files) > 0) {
    stop("Missing files: ", paste(names(missing_files), collapse = ", "))
  }
  
  # Load ear phenotype data
  psu_ear <- read.csv(ear_csv, na.strings = c("", "n/a", "NA"), skip = 1) %>%
    select(-description, -RK, -CC, -NIR, ear_rep = rep) %>%
    rename(rowid = row) %>%
    arrange(rowid) %>%
    group_by(rowid) %>%
    select(-ear_rep) %>%
    summarise_all(mean, na.rm = TRUE) %>%
    droplevels()
  
  ear_cols <- colnames(psu_ear)[-1]
  
  # Load plant phenotype data
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
  
  return(list(
    plant = psu_plant,
    ear = psu_ear,
    ear_cols = ear_cols
  ))
}

#' Combine phenotype and ionome data
#'
#' @param field_data List containing plant and ear data
#' @param ionome_data Raw ionome data
#' @return Combined data frame with phenotype and ionome measurements
#' @export
combine_phenotype_ionome <- function(field_data, ionome_data) {
  # Combine phenotype data
  pheno <- field_data$plant %>%
    select(
      rowid, Plot, Rep, Plot_Row, Plot_Column,
      Treatment, Genotype, PH:DTS, STW40:STWHV
    ) %>%
    left_join(field_data$ear, by = "rowid") %>%
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
  
  # Add ionome data in wide format
  ionome_wide <- ionome_data %>%
    select(rowid, tissue, all_of(RELIABLE_MINERALS)) %>%
    pivot_wider(
      names_from = tissue,
      values_from = all_of(RELIABLE_MINERALS)
    )
  
  combined_data <- pheno %>%
    left_join(ionome_wide, by = "rowid")
  
  return(combined_data)
}

#' Prepare ionome data for spatial analysis
#'
#' @param combined_data Data frame with phenotype and ionome data
#' @return Data frame formatted for spatial modeling
#' @export
prepare_ionome_for_spatial <- function(combined_data) {
  # Create tissue-specific variable names
  reliable_by_tissue <- expand.grid(
    mineral = RELIABLE_MINERALS,
    tissue = c("seed", "stover")
  ) %>%
    mutate(var_name = paste(mineral, tissue, sep = "_")) %>%
    pull(var_name)
  
  # Calculate mineral ratios (seed:stover)
  ratio_vars <- paste0(RELIABLE_MINERALS, "_ratio")
  
  ionome_formatted <- combined_data
  
  # Add ratio calculations
  for (mineral in RELIABLE_MINERALS) {
    seed_col <- paste0(mineral, "_seed")
    stover_col <- paste0(mineral, "_stover")
    ratio_col <- paste0(mineral, "_ratio")
    
    if (all(c(seed_col, stover_col) %in% colnames(ionome_formatted))) {
      ionome_formatted[[ratio_col]] <- 
        ionome_formatted[[seed_col]] / ionome_formatted[[stover_col]]
    }
  }
  
  # Ensure factors are properly set
  ionome_formatted$Rep <- as.factor(ionome_formatted$Rep)
  ionome_formatted$Treatment <- as.factor(ionome_formatted$Treatment)
  ionome_formatted$Genotype <- as.factor(ionome_formatted$Genotype)
  
  return(ionome_formatted)
}

#' Fit spatial mixed-effects model for ionome traits
#'
#' @param data Prepared ionome data
#' @param traits Vector of trait names to analyze
#' @param formula_rhs Right-hand side of model formula
#' @return Data frame with model effects
#' @export
fit_ionome_spatial_models <- function(data, traits, 
                                     formula_rhs = "Treatment*Genotype") {
  # Validate spatial coordinates
  required_coords <- c("Plot_Row", "Plot_Column")
  missing_coords <- setdiff(required_coords, colnames(data))
  
  if (length(missing_coords) > 0) {
    stop("Missing spatial coordinates: ", 
         paste(missing_coords, collapse = ", "))
  }
  
  effects_list <- list()
  
  for (trait in traits) {
    if (!trait %in% colnames(data)) {
      warning("Trait not found: ", trait)
      next
    }
    
    # Check for sufficient data
    valid_obs <- sum(!is.na(data[[trait]]))
    if (valid_obs < 10) {
      warning("Insufficient data for trait: ", trait)
      next
    }
    
    cat("Fitting spatial model for:", trait, "\n")
    
    # Construct formula
    formula_str <- paste(trait, "~", formula_rhs)
    
    # Fit spatial mixed-effects model
    model <- tryCatch({
      lme(
        fixed = as.formula(formula_str),
        random = ~ 1 | Rep,
        correlation = corSpher(form = ~ Plot_Row + Plot_Column),
        na.action = na.exclude,
        data = data
      )
    }, error = function(e) {
      warning("Model fitting failed for ", trait, ": ", e$message)
      NULL
    })
    
    if (!is.null(model)) {
      # Extract effects
      effects <- coef(summary(model)) %>%
        as.data.frame() %>%
        rownames_to_column("predictor")
      
      colnames(effects)[3] <- "Std.Error"
      colnames(effects)[6] <- "p.value"
      effects$response <- trait
      
      effects <- effects %>%
        select(response, everything()) %>%
        filter(!grepl("Intercept", predictor))
      
      effects_list[[trait]] <- effects
    }
  }
  
  # Combine all effects
  if (length(effects_list) > 0) {
    combined_effects <- bind_rows(effects_list) %>%
      mutate(p.adjust = p.adjust(p.value, method = "fdr")) %>%
      arrange(p.adjust)
    
    return(combined_effects)
  } else {
    stop("No models successfully fitted")
  }
}

#' Calculate standardized effect sizes for ionome traits
#'
#' @param data Original ionome data
#' @param traits Vector of trait names
#' @param significant_effects Data frame of significant effects
#' @return Data frame with standardized effects
#' @export
calculate_ionome_standardized_effects <- function(data, traits, significant_effects) {
  # Scale the data
  scaled_data <- data %>%
    select(rowid:Genotype) %>%
    bind_cols(
      data %>%
        select(all_of(traits)) %>%
        scale() %>%
        as.data.frame()
    )
  
  # Re-fit models with scaled data
  scaled_effects <- fit_ionome_spatial_models(scaled_data, traits)
  
  # Filter to significant terms from unscaled analysis
  standardized_effects <- scaled_effects %>%
    inner_join(
      significant_effects %>% select(response, predictor),
      by = c("response", "predictor")
    ) %>%
    mutate(
      CI.upper = Value + qnorm(0.975) * Std.Error,
      CI.lower = Value - qnorm(0.975) * Std.Error
    )
  
  return(standardized_effects)
}

#' Create ionome effects visualization
#'
#' @param effects_data Standardized effects data
#' @return ggplot object
#' @export
create_ionome_effects_plot <- function(effects_data) {
  # Prepare data for plotting
  plot_data <- effects_data %>%
    mutate(
      response = gsub("_", " ", response),
      predictor = factor(
        predictor,
        levels = c("Treatment-P", "GenotypeInv4m", "GenotypeInv4m:Treatment-P")
      )
    ) %>%
    group_by(response) %>%
    mutate(max_effect = max(abs(Value))) %>%
    ungroup() %>%
    group_by(predictor) %>%
    arrange(predictor) %>%
    mutate(.r = row_number()) %>%
    ungroup() %>%
    mutate(response = fct_reorder(response, max_effect))
  
  # Update predictor labels
  levels(plot_data$predictor) <- c("-P", "Inv4m", "Inv4m:-P")
  
  # Create the plot
  effects_plot <- plot_data %>%
    ggplot(aes(x = Value, y = response, shape = predictor)) +
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
    guides(shape = guide_legend(title = NULL)) +
    scale_shape_manual(values = c(19, 18, 13)) +
    theme_classic2(base_size = 20) +
    grids() +
    theme(
      legend.position = "top",
      axis.title.y = element_blank(),
      axis.text.y = element_text(size = 25, hjust = 0, face = "bold"),
      plot.caption = element_text(size = 25, hjust = 0)
    )
  
  return(effects_plot)
}

# Main analysis pipeline -------------------------------------------------

#' Main function to run spatial ionome analysis
#'
#' @export
run_ionome_spatial_analysis <- function() {
  cat("Starting spatial ionome analysis...\n")
  
  # File paths
  mineral_csv <- file.path(DATA_DIR, "PSU_inv4m_ionome_all.csv")
  plant_csv <- file.path(DATA_DIR, "22_NCS_PSU_LANGEBIO_FIELDS_PSU_P_field.csv")
  ear_csv <- file.path(DATA_DIR, "22_NCS_PSU_LANGEBIO_FIELDS_PSU_P_field_ear_pheno.csv")
  
  # 1. Load data
  cat("Loading ionome and phenotype data...\n")
  ionome_data <- load_ionome_data(mineral_csv)
  field_data <- load_field_phenotypes(plant_csv, ear_csv)
  
  # 2. Combine and prepare data
  cat("Combining and preparing data for spatial analysis...\n")
  combined_data <- combine_phenotype_ionome(field_data, ionome_data)
  ionome_spatial <- prepare_ionome_for_spatial(combined_data)
  
  # 3. Define traits for analysis
  reliable_by_tissue <- expand.grid(
    mineral = RELIABLE_MINERALS,
    tissue = c("seed", "stover")
  ) %>%
    mutate(var_name = paste(mineral, tissue, sep = "_")) %>%
    pull(var_name)
  
  ratio_vars <- paste0(RELIABLE_MINERALS, "_ratio")
  all_traits <- c(reliable_by_tissue, ratio_vars)
  
  # Filter to available traits
  available_traits <- intersect(all_traits, colnames(ionome_spatial))
  cat("Analyzing", length(available_traits), "ionome traits\n")
  
  # 4. Fit spatial models
  cat("Fitting spatial mixed-effects models...\n")
  spatial_effects <- fit_ionome_spatial_models(ionome_spatial, available_traits)
  
  # 5. Identify significant effects
  significant_effects <- spatial_effects %>%
    filter(p.adjust < 0.05, !grepl("Intercept", predictor))
  
  cat("Found", nrow(significant_effects), "significant effects\n")
  
  # 6. Calculate standardized effects
  if (nrow(significant_effects) > 0) {
    cat("Calculating standardized effect sizes...\n")
    standardized_effects <- calculate_ionome_standardized_effects(
      ionome_spatial, available_traits, significant_effects
    )
    
    # 7. Create visualization
    cat("Creating effects visualization...\n")
    effects_plot <- create_ionome_effects_plot(standardized_effects)
  } else {
    standardized_effects <- NULL
    effects_plot <- NULL
    warning("No significant effects found")
  }
  
  # 8. Save results
  cat("Saving results...\n")
  
  # Save raw effects
  write.csv(
    spatial_effects,
    file.path(OUTPUT_DIR, "ionome_spatial_effects.csv"),
    row.names = FALSE
  )
  
  # Save significant effects
  if (!is.null(standardized_effects)) {
    write.csv(
      standardized_effects,
      file.path(OUTPUT_DIR, "ionome_standardized_effects.csv"),
      row.names = FALSE
    )
    
    # Save plot
    ggsave(
      file.path(OUTPUT_DIR, "ionome_effects_plot.pdf"),
      effects_plot,
      width = 10, height = 8
    )
  }
  
  cat("Analysis complete!\n")
  cat("Results saved to:\n")
  cat("- ionome_spatial_effects.csv: Raw spatial model effects\n")
  if (!is.null(standardized_effects)) {
    cat("- ionome_standardized_effects.csv: Standardized effect sizes\n")
    cat("- ionome_effects_plot.pdf: Effects visualization\n")
  }
  
  return(list(
    effects = spatial_effects,
    standardized_effects = standardized_effects,
    plot = effects_plot,
    data = ionome_spatial
  ))
}

# Execute analysis if script is run directly
if (!interactive()) {
  results <- run_ionome_spatial_analysis()
}
