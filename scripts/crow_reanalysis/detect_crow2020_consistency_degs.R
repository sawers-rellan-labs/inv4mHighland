#' Detect Consistently Differentially Expressed Genes in Crow 2020 Data
#'
#' This script performs differential expression analysis on the Crow 2020 dataset
#' using a consistency-focused approach rather than mashr. The goal is to identify
#' genes with robust inv4m effects across temperature treatments, specifically
#' targeting PCNA2 and MCM5 temperature responses in proliferating tissues.
#'
#' Based on the analysis framework from detect_rnaseq_fdr_degs.R but adapted
#' for the Crow 2020 experimental design with temperature × genotype interactions.
#'
#' @author Francisco Rodriguez
#' @date 2025-08-06

# Load required libraries and configuration ------------------------------
source("common_config.R")

# Additional libraries specific to this analysis
optional_packages <- c("ggtext", "ggpubr", "corrplot", "pheatmap")
for (pkg in optional_packages) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    suppressMessages(library(pkg, character.only = TRUE))
  } else {
    warning("Optional package ", pkg, " not available. Some visualizations may not work.")
  }
}

# Analysis-specific configuration ----------------------------------------
MIN_CPM <- ANALYSIS_SETTINGS$min_cpm_threshold
MIN_SAMPLES <- ANALYSIS_SETTINGS$min_count_samples
FDR_THRESHOLD <- ANALYSIS_SETTINGS$fdr_threshold
LOGFC_THRESHOLD <- ANALYSIS_SETTINGS$logfc_threshold

# Temperature effect thresholds (different from general effects)
TEMP_LOGFC_THRESHOLD <- 0.5  # More sensitive for temperature effects

# Helper functions -------------------------------------------------------

#' Load Kallisto quantification data
#'
#' @param kallisto_dir Directory containing sample subdirectories with abundance.tsv
#' @param sample_metadata Data frame with sample information  
#' @return List containing count matrix and sample information
load_kallisto_data <- function(kallisto_dir, sample_metadata) {
  cat("Loading Kallisto quantification data...\n")
  
  if (!dir.exists(kallisto_dir)) {
    stop("Kallisto directory not found: ", kallisto_dir)
  }
  
  # Get sample directories
  sample_dirs <- list.dirs(kallisto_dir, recursive = FALSE, full.names = TRUE)
  sample_names <- basename(sample_dirs)
  
  if (length(sample_dirs) == 0) {
    stop("No sample directories found in: ", kallisto_dir)
  }
  
  cat("Found", length(sample_dirs), "sample directories\n")
  
  # Load transcript to gene mapping if available
  tx2gene <- NULL
  if (file.exists(file.path(DATA_DIR, "tx2gene.csv"))) {
    tx2gene <- read.csv(file.path(DATA_DIR, "tx2gene.csv"))
  } else {
    cat("No tx2gene mapping found, using transcript-level analysis\n")
  }
  
  # Load quantifications
  txi <- load_kallisto_results(sample_dirs, tx2gene)
  
  # Match with sample metadata
  run_ids <- intersect(sample_names, sample_metadata$Run)
  
  if (length(run_ids) == 0) {
    stop("No matching samples found between Kallisto results and metadata")
  }
  
  cat("Matched", length(run_ids), "samples with metadata\n")
  
  # Filter data to matched samples
  sample_idx <- match(run_ids, sample_names)
  metadata_idx <- match(run_ids, sample_metadata$Run)
  
  filtered_counts <- txi$counts[, sample_idx, drop = FALSE]
  filtered_metadata <- sample_metadata[metadata_idx, ]
  
  # Ensure sample order matches
  colnames(filtered_counts) <- run_ids
  rownames(filtered_metadata) <- run_ids
  
  return(list(
    counts = filtered_counts,
    abundance = txi$abundance[, sample_idx, drop = FALSE],
    length = txi$length[, sample_idx, drop = FALSE],
    sample_info = filtered_metadata,
    tx2gene = tx2gene
  ))
}

#' Prepare sample metadata for analysis
#'
#' @param raw_metadata Raw sample metadata from SRA
#' @return Cleaned and parsed sample metadata
prepare_sample_metadata <- function(raw_metadata) {
  cat("Preparing sample metadata...\n")
  
  if (nrow(raw_metadata) == 0) {
    stop("No sample metadata provided")
  }
  
  # Parse sample information (based on check_sra_metadata.R patterns)
  metadata <- raw_metadata %>%
    # Extract genotype information
    mutate(
      genotype = case_when(
        grepl("PT_NIL", old_line, ignore.case = TRUE) ~ "INV4",
        grepl("Mi21_NIL", old_line, ignore.case = TRUE) ~ "INV4", 
        grepl("B73|Mo17", old_line, ignore.case = TRUE) ~ "CTRL",
        TRUE ~ "UNKNOWN"
      ),
      # Extract temperature from inv_temp column
      temperature = case_when(
        grepl("16|18|20", inv_temp) ~ "cold",
        grepl("22|24|25", inv_temp) ~ "optimal",
        grepl("28|30|32", inv_temp) ~ "warm",
        TRUE ~ "unknown"
      ),
      # Classify tissue types
      tissue_type = case_when(
        grepl("SAM", tissue, ignore.case = TRUE) ~ "proliferating",
        grepl("root", tissue, ignore.case = TRUE) ~ "proliferating",
        grepl("leaf", tissue, ignore.case = TRUE) ~ "mixed",
        TRUE ~ "other"
      ),
      # Create interaction term
      geno_temp = interaction(genotype, temperature),
      # Convert to factors with appropriate levels
      genotype = factor(genotype, levels = c("CTRL", "INV4")),
      temperature = factor(temperature, levels = c("cold", "optimal", "warm")),
      tissue_type = factor(tissue_type),
      old_line = factor(old_line)
    ) %>%
    # Filter for proliferating tissues (primary focus)
    filter(
      tissue_type == "proliferating",
      genotype != "UNKNOWN",
      temperature != "unknown"
    )
  
  cat("Filtered to", nrow(metadata), "proliferating tissue samples\n")
  cat("Genotypes:", table(metadata$genotype), "\n")
  cat("Temperatures:", table(metadata$temperature), "\n")
  
  return(metadata)
}

#' Create filtered DGEList object
#'
#' @param counts_data Kallisto data object
#' @param min_cpm Minimum CPM threshold
#' @param min_samples Minimum number of samples
#' @return Filtered and normalized DGEList object
create_filtered_dgelist <- function(counts_data, min_cpm = MIN_CPM, min_samples = MIN_SAMPLES) {
  cat("Creating filtered DGEList object...\n")
  
  # Create DGEList
  y <- DGEList(
    counts = round(counts_data$counts),
    samples = counts_data$sample_info
  )
  
  # Add group information
  y$samples$group <- interaction(
    y$samples$genotype,
    y$samples$temperature
  )
  
  cat("Initial genes:", nrow(y), "\n")
  cat("Initial samples:", ncol(y), "\n")
  
  # Filter lowly expressed genes
  # Require at least min_cpm CPM in at least min_samples samples
  keep <- rowSums(cpm(y) > min_cpm) >= min_samples
  y_filtered <- y[keep, ]
  
  cat("Genes after filtering:", nrow(y_filtered), "\n")
  
  # Check library sizes and flag low quality samples
  y_filtered$samples$lib_size_million <- y_filtered$samples$lib.size / 1e6
  low_lib_threshold <- quantile(y_filtered$samples$lib.size, 0.1)
  y_filtered$samples$low_lib_size <- y_filtered$samples$lib.size < low_lib_threshold
  
  # Normalize for composition bias
  y_normalized <- calcNormFactors(y_filtered)
  
  cat("Library sizes (millions):\n")
  print(summary(y_normalized$samples$lib_size_million))
  
  return(y_normalized)
}

#' Perform exploratory data analysis
#'
#' @param dge_list DGEList object
#' @param output_dir Output directory for plots
#' @return List with MDS data and plots
perform_eda <- function(dge_list, output_dir) {
  cat("Performing exploratory data analysis...\n")
  
  # Calculate MDS
  mds <- plotMDS(dge_list, ndim = 4, plot = FALSE)
  
  # Prepare data for plotting
  mds_data <- dge_list$samples %>%
    mutate(
      dim1 = mds$x,
      dim2 = mds$y,
      dim3 = mds$cmdscale.out[, 3],
      dim4 = mds$cmdscale.out[, 4]
    )
  
  # Create genotype separation plot (dimensions 1-2)
  genotype_plot <- mds_data %>%
    ggplot(aes(x = dim1, y = dim2, fill = genotype, shape = temperature)) +
    geom_point(size = 4, alpha = 0.8) +
    scale_fill_manual(
      values = COLOR_SCHEMES$genotype,
      labels = c("CTRL" = "Control", "INV4" = "*Inv4m*")
    ) +
    scale_shape_manual(values = c(21, 22, 24)) +
    labs(
      x = paste0("Dimension 1 (", round(100 * mds$var.explained[1], 1), "%)"),
      y = paste0("Dimension 2 (", round(100 * mds$var.explained[2], 1), "%)"),
      title = "Sample Clustering by Genotype"
    ) +
    guides(
      fill = guide_legend(
        title = "Genotype",
        override.aes = list(shape = 21, size = 5)
      ),
      shape = guide_legend(
        title = "Temperature",
        override.aes = list(size = 5)
      )
    ) +
    theme_classic(base_size = PLOT_SETTINGS$base_size) +
    theme(
      legend.text = element_markdown(),
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
  
  # Create temperature effect plot (dimensions 3-4)  
  temperature_plot <- mds_data %>%
    ggplot(aes(x = dim3, y = dim4, fill = temperature, shape = genotype)) +
    geom_point(size = 4, alpha = 0.8) +
    scale_fill_manual(values = COLOR_SCHEMES$temperature) +
    scale_shape_manual(values = c(21, 24)) +
    labs(
      x = paste0("Dimension 3 (", round(100 * mds$var.explained[3], 1), "%)"),
      y = paste0("Dimension 4 (", round(100 * mds$var.explained[4], 1), "%)"),
      title = "Sample Clustering by Temperature"
    ) +
    guides(
      fill = guide_legend(
        title = "Temperature",
        override.aes = list(shape = 21, size = 5)
      ),
      shape = guide_legend(
        title = "Genotype", 
        override.aes = list(size = 5)
      )
    ) +
    theme_classic(base_size = PLOT_SETTINGS$base_size) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  
  # Save plots
  plot_files <- list(
    genotype = file.path(output_dir, "plots", "crow2020_mds_genotype.pdf"),
    temperature = file.path(output_dir, "plots", "crow2020_mds_temperature.pdf")
  )
  
  ggsave(plot_files$genotype, genotype_plot, 
         width = PLOT_SETTINGS$width, height = PLOT_SETTINGS$height, 
         dpi = PLOT_SETTINGS$dpi)
  
  ggsave(plot_files$temperature, temperature_plot,
         width = PLOT_SETTINGS$width, height = PLOT_SETTINGS$height,
         dpi = PLOT_SETTINGS$dpi)
  
  cat("EDA plots saved to:", dirname(plot_files$genotype), "\n")
  
  return(list(
    mds_data = mds_data,
    plots = list(genotype = genotype_plot, temperature = temperature_plot),
    files = plot_files
  ))
}

#' Perform differential expression analysis
#'
#' @param dge_list DGEList object
#' @return List containing voom object and fitted model
perform_differential_analysis <- function(dge_list) {
  cat("Performing differential expression analysis...\n")
  
  # Create design matrix
  # Focus on genotype × temperature interaction effects
  design <- model.matrix(
    ~ old_line + genotype * temperature,
    data = dge_list$samples
  )
  
  cat("Design matrix dimensions:", dim(design), "\n")
  cat("Design matrix columns:\n")
  print(colnames(design))
  
  # Voom transformation with quality weights
  voom_result <- voom(dge_list, design = design, plot = FALSE)
  
  # Fit linear model
  fit <- lmFit(voom_result, design)
  ebayes_fit <- eBayes(fit, robust = TRUE)
  
  # Save intermediate results
  saveRDS(voom_result, file.path(OUTPUT_DIR, "tables", "crow2020_voom_object.rds"))
  saveRDS(ebayes_fit, file.path(OUTPUT_DIR, "tables", "crow2020_ebayes_fit.rds"))
  
  return(list(
    voom = voom_result,
    fit = ebayes_fit,
    design = design
  ))
}

#' Extract differential expression results
#'
#' @param ebayes_fit eBayes fitted model
#' @param gene_annotation Gene annotation data
#' @return Data frame with DE results
extract_de_results <- function(ebayes_fit, gene_annotation = NULL) {
  cat("Extracting differential expression results...\n")
  
  # Define predictors of interest
  predictors <- c(
    "genotypeINV4",                          # Main genotype effect
    "temperatureoptimal",                    # Temperature: optimal vs cold
    "temperaturewarm",                       # Temperature: warm vs cold  
    "genotypeINV4:temperatureoptimal",       # Interaction: INV4 × optimal
    "genotypeINV4:temperaturewarm"           # Interaction: INV4 × warm
  )
  
  # Check which predictors exist in the model
  available_predictors <- intersect(predictors, colnames(ebayes_fit$coefficients))
  
  if (length(available_predictors) == 0) {
    stop("No target predictors found in model. Available: ",
         paste(colnames(ebayes_fit$coefficients), collapse = ", "))
  }
  
  cat("Extracting results for predictors:", paste(available_predictors, collapse = ", "), "\n")
  
  # Extract results for each predictor
  results_list <- list()
  
  for (predictor in available_predictors) {
    result <- topTable(
      ebayes_fit,
      coef = predictor,
      sort.by = "none",
      n = Inf
    )
    
    result$predictor <- predictor
    result$gene <- rownames(result)
    
    # Calculate confidence intervals
    conf_interval <- qt(0.975, ebayes_fit$df.residual + ebayes_fit$df.prior) *
      ebayes_fit$stdev.unscaled[, predictor] * sqrt(ebayes_fit$s2.post)
    
    result$logFC_lower <- result$logFC - conf_interval
    result$logFC_upper <- result$logFC + conf_interval
    
    results_list[[predictor]] <- result
  }
  
  # Combine all results
  combined_results <- bind_rows(results_list)
  
  # Add gene annotation if provided
  if (!is.null(gene_annotation) && nrow(gene_annotation) > 0) {
    combined_results <- combined_results %>%
      left_join(gene_annotation, by = c("gene" = "gene_id"))
  }
  
  return(combined_results)
}

#' Classify differential expression results
#'
#' @param de_results DE results data frame
#' @param target_genes List of target genes to highlight
#' @return Classified DE results
classify_de_results <- function(de_results, target_genes = TARGET_GENES) {
  cat("Classifying differential expression results...\n")
  
  results_classified <- de_results %>%
    mutate(
      # Basic significance
      is_significant = adj.P.Val < FDR_THRESHOLD,
      neglogP = -log10(adj.P.Val),
      
      # Different thresholds for different effect types
      logfc_threshold = case_when(
        grepl("temperature", predictor) ~ TEMP_LOGFC_THRESHOLD,
        TRUE ~ LOGFC_THRESHOLD
      ),
      
      # Direction classification
      upregulated = (logFC > logfc_threshold) & is_significant,
      downregulated = (logFC < -logfc_threshold) & is_significant,
      
      regulation = case_when(
        upregulated ~ "Upregulated",
        downregulated ~ "Downregulated", 
        is_significant ~ "Significant (small effect)",
        TRUE ~ "Not significant"
      ),
      
      # Target gene classification
      is_target_gene = gene %in% unlist(target_genes),
      target_gene_name = case_when(
        gene %in% target_genes$PCNA2 ~ "PCNA2",
        gene %in% target_genes$MCM5 ~ "MCM5",
        gene %in% target_genes$ZmCCT ~ "ZmCCT",
        gene %in% target_genes$ZmCCT9 ~ "ZmCCT9",
        TRUE ~ NA_character_
      ),
      
      # Effect type classification
      effect_type = case_when(
        predictor == "genotypeINV4" ~ "Main genotype effect",
        grepl("temperature", predictor) & !grepl(":", predictor) ~ "Main temperature effect",
        grepl(":", predictor) ~ "Genotype × Temperature interaction",
        TRUE ~ "Other"
      )
    )
  
  # Summary statistics
  cat("Classification summary:\n")
  print(table(results_classified$regulation, results_classified$effect_type))
  
  if (any(results_classified$is_target_gene)) {
    cat("\nTarget gene results:\n")
    target_summary <- results_classified %>%
      filter(is_target_gene) %>%
      select(target_gene_name, predictor, logFC, adj.P.Val, regulation)
    print(target_summary)
  }
  
  return(results_classified)
}

#' Create target gene heatmap
#'
#' @param voom_data Voom-normalized expression data
#' @param sample_info Sample metadata
#' @param target_genes List of target genes
#' @param output_file Output file path
create_target_gene_heatmap <- function(voom_data, sample_info, target_genes, output_file) {
  cat("Creating target gene heatmap...\n")
  
  # Get target gene IDs
  target_gene_ids <- unlist(target_genes)
  available_targets <- intersect(target_gene_ids, rownames(voom_data))
  
  if (length(available_targets) == 0) {
    warning("No target genes found in expression data")
    return(NULL)
  }
  
  cat("Found", length(available_targets), "target genes in data\n")
  
  # Extract expression data for target genes
  target_expr <- voom_data[available_targets, , drop = FALSE]
  
  # Create annotation for samples
  sample_annotation <- sample_info %>%
    select(genotype, temperature, tissue) %>%
    as.data.frame()
  rownames(sample_annotation) <- rownames(sample_info)
  
  # Create annotation colors
  annotation_colors <- list(
    genotype = COLOR_SCHEMES$genotype,
    temperature = COLOR_SCHEMES$temperature
  )
  
  # Create gene labels
  gene_labels <- sapply(rownames(target_expr), function(gene_id) {
    for (gene_name in names(target_genes)) {
      if (gene_id %in% target_genes[[gene_name]]) {
        return(gene_name)
      }
    }
    return(gene_id)
  })
  
  # Create heatmap
  pheatmap(
    target_expr,
    annotation_col = sample_annotation,
    annotation_colors = annotation_colors,
    labels_row = gene_labels,
    scale = "row",
    clustering_distance_rows = "correlation",
    clustering_distance_cols = "correlation", 
    filename = output_file,
    width = 12,
    height = 8
  )
  
  cat("Target gene heatmap saved to:", output_file, "\n")
  return(target_expr)
}

# Main analysis pipeline -------------------------------------------------

#' Main function to run the complete analysis
#'
#' @export
run_crow2020_analysis <- function() {
  cat("Starting Crow 2020 consistency-focused DEG analysis...\n")
  cat("Analysis based on:", length(TARGET_GENES), "target genes\n")
  cat("Output directory:", OUTPUT_DIR, "\n")
  
  # 1. Load sample metadata
  cat("\n=== Step 1: Loading sample metadata ===\n")
  raw_metadata <- load_sra_metadata()
  
  if (nrow(raw_metadata) == 0) {
    stop("No sample metadata available. Please run check_sra_metadata.R first.")
  }
  
  # For now, create sample metadata from the check_sra_metadata.R output
  # This would typically be loaded from the processed metadata file
  if (file.exists(file.path(DATA_DIR, "crow2020_apical_tissue_samples.tab"))) {
    sample_metadata <- read.table(
      file.path(DATA_DIR, "crow2020_apical_tissue_samples.tab"),
      header = TRUE,
      sep = "\t"
    )
  } else {
    warning("Crow 2020 sample metadata not found, using simulated metadata for testing")
    # Create minimal test metadata
    sample_metadata <- data.frame(
      old_line = c("PT_NIL", "PT_NIL", "Mi21_NIL", "Mi21_NIL"),
      inv_temp = c("cold_16", "warm_28", "cold_16", "warm_28"),
      tissue = c("SAM", "SAM", "root", "root"),
      samp = c("test1", "test2", "test3", "test4"),
      Run = c("SRR001", "SRR002", "SRR003", "SRR004")
    )
  }
  
  prepared_metadata <- prepare_sample_metadata(sample_metadata)
  
  # 2. Load Kallisto quantification data
  cat("\n=== Step 2: Loading expression data ===\n")
  kallisto_dir <- file.path(OUTPUT_DIR, "kallisto")
  
  if (!dir.exists(kallisto_dir)) {
    warning("Kallisto results directory not found: ", kallisto_dir)
    warning("Please run process_crow2020_rnaseq.sh first to quantify samples")
    
    # For testing, create simulated expression data
    cat("Creating simulated expression data for testing...\n")
    n_genes <- 1000
    n_samples <- nrow(prepared_metadata)
    
    counts_matrix <- matrix(
      rpois(n_genes * n_samples, lambda = 100),
      nrow = n_genes,
      ncol = n_samples
    )
    rownames(counts_matrix) <- paste0("gene_", 1:n_genes)
    colnames(counts_matrix) <- prepared_metadata$Run
    
    # Add target genes with some signal
    target_gene_ids <- unlist(TARGET_GENES)[1:min(length(unlist(TARGET_GENES)), 10)]
    if (length(target_gene_ids) > 0) {
      # Replace some genes with target gene IDs
      rownames(counts_matrix)[1:length(target_gene_ids)] <- target_gene_ids
      
      # Add some differential expression signal for testing
      inv4_samples <- prepared_metadata$genotype == "INV4"
      if (any(inv4_samples)) {
        counts_matrix[1:length(target_gene_ids), inv4_samples] <- 
          counts_matrix[1:length(target_gene_ids), inv4_samples] * 2
      }
    }
    
    expression_data <- list(
      counts = counts_matrix,
      abundance = counts_matrix,
      length = matrix(1000, nrow = nrow(counts_matrix), ncol = ncol(counts_matrix)),
      sample_info = prepared_metadata,
      tx2gene = NULL
    )
  } else {
    expression_data <- load_kallisto_data(kallisto_dir, prepared_metadata)
  }
  
  # 3. Create filtered DGEList
  cat("\n=== Step 3: Filtering and normalization ===\n")
  dge_list <- create_filtered_dgelist(expression_data)
  
  # 4. Exploratory data analysis
  cat("\n=== Step 4: Exploratory data analysis ===\n")
  eda_results <- perform_eda(dge_list, OUTPUT_DIR)
  
  # 5. Differential expression analysis
  cat("\n=== Step 5: Differential expression analysis ===\n")
  de_analysis <- perform_differential_analysis(dge_list)
  
  # 6. Extract and classify results
  cat("\n=== Step 6: Extracting and classifying results ===\n")
  de_results <- extract_de_results(de_analysis$fit, genes)
  classified_results <- classify_de_results(de_results, TARGET_GENES)
  
  # 7. Create target gene heatmap
  cat("\n=== Step 7: Creating visualizations ===\n")
  heatmap_file <- file.path(OUTPUT_DIR, "plots", "crow2020_target_genes_heatmap.pdf")
  target_expr <- create_target_gene_heatmap(
    de_analysis$voom$E,
    dge_list$samples,
    TARGET_GENES,
    heatmap_file
  )
  
  # 8. Save results
  cat("\n=== Step 8: Saving results ===\n")
  
  # All results
  write.csv(
    classified_results,
    file.path(OUTPUT_DIR, "tables", "crow2020_all_de_results.csv"),
    row.names = FALSE
  )
  
  # Significant DEGs only
  sig_degs <- classified_results %>%
    filter(is_significant) %>%
    arrange(predictor, adj.P.Val)
  
  write.csv(
    sig_degs,
    file.path(OUTPUT_DIR, "tables", "crow2020_significant_degs.csv"),
    row.names = FALSE
  )
  
  # Target gene results
  target_results <- classified_results %>%
    filter(is_target_gene) %>%
    arrange(target_gene_name, predictor, adj.P.Val)
  
  write.csv(
    target_results,
    file.path(OUTPUT_DIR, "tables", "crow2020_target_gene_results.csv"),
    row.names = FALSE
  )
  
  # Generate summary report
  cat("\n=== Analysis Summary ===\n")
  cat("Total genes analyzed:", nrow(de_results), "\n")
  cat("Significant DEGs:", nrow(sig_degs), "\n")
  cat("Target genes found:", sum(classified_results$is_target_gene), "\n")
  
  if (nrow(target_results) > 0) {
    cat("\nTarget gene summary:\n")
    target_summary <- target_results %>%
      group_by(target_gene_name, effect_type) %>%
      summarise(
        n_significant = sum(is_significant),
        mean_logFC = mean(logFC),
        .groups = "drop"
      )
    print(target_summary)
  }
  
  # Effect type summary
  cat("\nEffect type summary:\n")
  effect_summary <- sig_degs %>%
    group_by(effect_type, regulation) %>%
    summarise(count = n(), .groups = "drop")
  print(effect_summary)
  
  cat("\nResults saved to:\n")
  cat("- All results:", file.path(OUTPUT_DIR, "tables", "crow2020_all_de_results.csv"), "\n")
  cat("- Significant DEGs:", file.path(OUTPUT_DIR, "tables", "crow2020_significant_degs.csv"), "\n") 
  cat("- Target genes:", file.path(OUTPUT_DIR, "tables", "crow2020_target_gene_results.csv"), "\n")
  cat("- Plots:", file.path(OUTPUT_DIR, "plots/"), "\n")
  
  return(list(
    expression_data = expression_data,
    dge_list = dge_list,
    de_analysis = de_analysis,
    results = classified_results,
    target_results = target_results,
    eda = eda_results
  ))
}

# Execute analysis if script is run directly
if (!interactive()) {
  results <- run_crow2020_analysis()
}