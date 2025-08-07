#' Detect Consistently Differentially Expressed Genes using FDR
#'
#' This script performs differential expression analysis for the PSU field
#' experiment using a standard FDR-based approach rather than mashr.
#' The goal is to identify genes with consistent inv4m effects across
#' phosphorus treatments.
#'
#' @author Fausto Rodriguez
#' @date 2025-08-04

# Load required libraries -----------------------------------------------
library(edgeR)         # Bioconductor: differential expression analysis
library(limma)         # Bioconductor: linear models for microarray/RNA-seq
library(rtracklayer)   # Bioconductor: genomic annotation handling
library(GenomicRanges) # Bioconductor: genomic ranges operations
library(dplyr)         # CRAN: data manipulation
library(ggplot2)       # CRAN: plotting
library(ggpubr)        # CRAN: publication ready plots
library(ggtext)        # CRAN: formatted text in plots

# Configuration ----------------------------------------------------------
DATA_DIR <- "data"
OUTPUT_DIR <- "results"
MIN_LIB_SIZE <- 2e7  # Minimum library size threshold
LOGFC_THRESHOLD_GENERAL <- 2     # LogFC threshold for most predictors
LOGFC_THRESHOLD_TISSUE <- 0.7    # LogFC threshold for tissue effects
FDR_THRESHOLD <- 0.05

# Helper functions -------------------------------------------------------

#' Load and prepare expression data
#'
#' @param counts_file Path to counts file
#' @param metadata_file Path to metadata file
#' @return List containing processed DGEList and sample information
load_expression_data <- function(counts_file, metadata_file) {
  # Validate inputs
  if (!file.exists(counts_file)) {
    stop("Counts file not found: ", counts_file)
  }
  if (!file.exists(metadata_file)) {
    stop("Metadata file not found: ", metadata_file)
  }
  
  # Load data
  counts <- read.csv(counts_file)
  sample_info <- read.csv(metadata_file)
  
  # Prepare count matrix
  genes <- data.frame(gene = counts[, 2])
  count_matrix <- as.matrix(counts[, -c(1:2)])
  rownames(count_matrix) <- genes$gene
  
  # Match samples with metadata
  tag <- sample_info$side_tag
  names(tag) <- sample_info$library
  
  sample_names <- tag[colnames(count_matrix)]
  colnames(count_matrix) <- sample_names
  
  # Ensure sample order matches
  if (!all(sample_names %in% sample_info$side_tag)) {
    stop("Sample name mismatch between counts and metadata")
  }
  
  sample_info <- sample_info[match(sample_names, sample_info$side_tag), ]
  
  return(list(
    counts = count_matrix,
    sample_info = sample_info,
    genes = genes
  ))
}

#' Create DGEList object with quality filtering
#'
#' @param counts Count matrix
#' @param sample_info Sample metadata
#' @param min_lib_size Minimum library size threshold
#' @return Filtered DGEList object
create_filtered_dgelist <- function(counts, sample_info, min_lib_size) {
  # Create DGEList
  y <- DGEList(counts = counts, samples = sample_info)
  y$group <- interaction(y$samples$Treatment, y$samples$Genotype)
  
  # Filter lowly expressed genes
  keep <- filterByExpr(y, group = y$group)
  y_filtered <- y[keep, ]
  
  # Identify low quality libraries
  y_filtered$samples$lowCount <- y_filtered$samples$lib.size < min_lib_size
  
  # Filter out low quality libraries
  y_final <- y_filtered[, !y_filtered$samples$lowCount]
  
  # Normalize
  y_final <- calcNormFactors(y_final)
  
  y_final
}

#' Perform MDS analysis and generate plots
#'
#' @param dge_list DGEList object
#' @param output_dir Directory for saving plots
#' @return MDS coordinates data frame
perform_mds_analysis <- function(dge_list, output_dir) {
  # Calculate MDS
  mds <- plotMDS(dge_list, pch = 21, plot = FALSE)
  mds2 <- plotMDS(dge_list, pch = 21, dim.plot = c(3, 4), plot = FALSE)
  
  # Prepare data for plotting
  mds_data <- dge_list$samples
  mds_data$dim1 <- mds$x
  mds_data$dim2 <- mds$y
  mds_data$dim3 <- mds2$x
  mds_data$dim4 <- mds2$y
  
  # Format factor levels
  mds_data$Treatment <- factor(mds_data$Treatment)
  levels(mds_data$Treatment) <- c("+P", "-P")
  mds_data$Genotype <- factor(mds_data$Genotype)
  mds_data$RNA_Plant <- factor(mds_data$RNA_Plant)
  
  # Create MDS plots
  generate_mds_plots(mds_data, mds, mds2, output_dir)
  
  mds_data
}

#' Generate MDS plots
#'
#' @param mds_data MDS data with coordinates
#' @param mds MDS object for dimensions 1-2
#' @param mds2 MDS object for dimensions 3-4
#' @param output_dir Output directory
generate_mds_plots <- function(mds_data, mds, mds2, output_dir) {
  # Genotype plot (dimensions 3-4)
  labels <- c("CTRL", "*Inv4m*")
  names(labels) <- c("CTRL", "INV4")
  
  genotype_plot <- mds_data %>%
    ggplot(aes(x = dim3, y = dim4, fill = Genotype, shape = Treatment)) +
    xlab(paste0("dim3 (", round(100 * mds2$var.explained[3]), "%)")) +
    ylab(paste0("dim4 (", round(100 * mds2$var.explained[4]), "%)")) +
    geom_point(size = 4) +
    scale_fill_viridis_d(direction = -1, labels = labels) +
    scale_shape_manual(values = c(24, 21)) +
    guides(
      shape = "none",
      fill = guide_legend(
        title = "Genotype",
        order = 2,
        override.aes = list(
          geom = "point",
          shape = 22,
          size = 7,
          reverse = TRUE
        )
      )
    ) +
    theme_classic2(base_size = 25) +
    theme(
      legend.position = c(0.89, 0.9),
      legend.text = element_markdown(),
      legend.spacing = unit(0, "line"),
      legend.box.spacing = unit(0, "line")
    )
  
  # Leaf tissue plot (dimensions 1-2)
  tissue_plot <- mds_data %>%
    mutate(leaf = factor(leaf_tissue)) %>%
    ggplot(aes(x = dim1, y = dim2)) +
    xlab(paste0("dim1 (", round(100 * mds$var.explained[1]), "%)")) +
    ylab(paste0("dim2 (", round(100 * mds$var.explained[2]), "%)")) +
    geom_point(aes(fill = leaf, shape = Treatment), size = 4) +
    scale_fill_viridis_d() +
    scale_shape_manual(values = c(24, 21)) +
    guides(
      shape = guide_legend(
        title = "Treatment",
        order = 1,
        override.aes = list(size = 7)
      ),
      fill = guide_legend(
        title = "Leaf",
        order = 2,
        override.aes = list(geom = "point", shape = 22, size = 7)
      )
    ) +
    theme_classic2(base_size = 25) +
    theme(
      legend.box = "horizontal",
      legend.spacing = unit(0, "line"),
      legend.box.spacing = unit(0, "in"),
      legend.position = c(0.75, 0.17)
    )
  
  # Save plots
  png(file.path(output_dir, "inv4m_expression_MDS.png"),
      width = 7, height = 7, units = "in", res = 300)
  print(tissue_plot)
  dev.off()
  
  return(list(genotype = genotype_plot, tissue = tissue_plot))
}

#' Perform differential expression analysis
#'
#' @param dge_list DGEList object
#' @param design_formula Design formula for the model
#' @return List containing fitted model and results
perform_de_analysis <- function(dge_list, design_formula) {
  # Create design matrix
  design <- model.matrix(design_formula, data = dge_list$samples)
  
  # Perform voom transformation
  voom_result <- voom(dge_list, design = design, plot = FALSE)
  
  # Fit linear model
  fit <- lmFit(voom_result)
  ebayes_fit <- eBayes(fit, robust = TRUE)
  
  list(
    voom = voom_result,
    fit = ebayes_fit,
    design = design
  )
}

#' Extract differential expression results
#'
#' @param ebayes_fit eBayes fitted object
#' @param predictors Vector of predictor names to extract
#' @param gene_annotation Gene annotation data frame
#' @return Data frame with DE results
get_DE_effects <- function(ebayes_fit, predictors, gene_annotation) {
  results_list <- list()
  
  for (predictor in predictors) {
    # Extract results for this predictor
    result <- topTable(ebayes_fit, coef = predictor, sort.by = "none", n = Inf)
    result$predictor <- predictor
    result <- tibble::rownames_to_column(result, "Response")
    
    # Calculate confidence intervals
    conf_interval <- qt(0.975, ebayes_fit$df.residual + ebayes_fit$df.prior) *
      ebayes_fit$stdev.unscaled[, predictor] * sqrt(ebayes_fit$s2.post)
    
    result$upper <- result$logFC + conf_interval
    result$lower <- result$logFC - conf_interval
    
    results_list[[predictor]] <- result
  }
  
  # Combine results
  combined_results <- dplyr::bind_rows(results_list)
  
  # Add gene annotation
  if (!is.null(gene_annotation)) {
    combined_results <- combined_results %>%
      left_join(gene_annotation, by = c("Response" = "gene_model"),
                relationship = "many-to-many")
  }
  
  return(combined_results)
}

#' Add statistical classifications to DE results
#'
#' @param effects DE results data frame
#' @param logfc_thresholds Named vector of logFC thresholds per predictor type
#' @param fdr_threshold FDR significance threshold
#' @return DE results with classification columns
classify_de_results <- function(effects, logfc_thresholds, fdr_threshold) {
  effects <- effects %>%
    mutate(
      predictor = factor(predictor, levels = names(logfc_thresholds)),
      P = adj.P.Val,
      neglogP = -log10(adj.P.Val),
      is_significant = adj.P.Val < fdr_threshold
    )
  
  # Apply different thresholds based on predictor type
  effects <- effects %>%
    mutate(
      logfc_threshold = case_when(
        predictor == "leaf_tissue" ~ logfc_thresholds["leaf_tissue"],
        TRUE ~ logfc_thresholds["general"]
      ),
      upregulated = (logFC > logfc_threshold) & is_significant,
      downregulated = (logFC < -logfc_threshold) & is_significant,
      regulation = case_when(
        is_significant & upregulated ~ "Upregulated",
        is_significant & downregulated ~ "Downregulated",
        TRUE ~ "Unregulated"
      ),
      is_DEG = is_significant & regulation != "Unregulated"
    )
  
  return(effects)
}

#' Load genomic coordinates and define inversion regions
#'
#' @param gff_path Path to GFF file
#' @return List with gene coordinates and inversion boundaries
load_genomic_coordinates <- function(gff_path) {
  # Load GFF
  v5_gff <- rtracklayer::import(gff_path) %>%
    subset(type == "gene" & seqnames %in% 1:10)
  
  genes <- as.data.frame(v5_gff)
  genes$ID <- gsub("gene:", "", genes$ID)
  
  # Define inversion boundaries (from MCScan analysis)
  inv4m_start <- genes[genes$ID == "Zm00001eb190470", "start"]
  inv4m_end <- genes[genes$ID == "Zm00001eb194800", "end"]
  
  # Define introgression boundaries (from genotyping data)
  introgression_start <- 157012149
  introgression_end <- 195900523
  
  # Identify gene sets
  inv4m_genes <- genes %>%
    filter(seqnames == 4, start >= inv4m_start, end <= inv4m_end) %>%
    pull(ID)
  
  shared_introgression_genes <- genes %>%
    filter(
      seqnames == 4,
      start >= introgression_start,
      end <= introgression_end
    ) %>%
    pull(ID)
  
  flanking_genes <- setdiff(shared_introgression_genes, inv4m_genes)
  
  return(list(
    genes = genes,
    inv4m_start = inv4m_start,
    inv4m_end = inv4m_end,
    introgression_start = introgression_start,
    introgression_end = introgression_end,
    inv4m_genes = inv4m_genes,
    shared_introgression_genes = shared_introgression_genes,
    flanking_genes = flanking_genes
  ))
}

#' Add Mahalanobis outlier detection
#'
#' @param data DE results data frame
#' @param distance_quantile Quantile threshold for outlier detection
#' @param fdr_threshold FDR threshold
#' @return Data frame with outlier classifications
add_mahalanobis_outliers <- function(data, distance_quantile = 0.05,
                                     fdr_threshold = 0.05) {
  data_with_outliers <- lapply(
    split(data, factor(data$predictor)),
    function(per_predictor) {
      bivariate <- per_predictor %>% select(logFC, neglogP)
      per_predictor$mahalanobis <- mahalanobis(
        x = bivariate,
        center = colMeans(bivariate),
        cov = cov(bivariate)
      )
      per_predictor
    }
  ) %>%
    bind_rows()
  
  cutoff <- qchisq(p = 1 - distance_quantile, df = 2)
  
  data_with_outliers$is_mh_outlier <- (data_with_outliers$adj.P.Val < fdr_threshold) &
    (data_with_outliers$mahalanobis > cutoff)
  
  return(
    data_with_outliers %>%
      ungroup() %>%
      group_by(predictor, regulation) %>%
      arrange(-mahalanobis, .by_group = TRUE) %>%
      ungroup()
  )
}

# Main analysis pipeline -------------------------------------------------


cat("Starting differential expression analysis...\n")
  
# 1. Load and prepare data
cat("Loading expression data...\n")
expr_data <- load_expression_data(
  file.path(DATA_DIR, "inv4mRNAseq_gene_sample_exp.csv"),
  file.path(DATA_DIR, "PSU-PHO22_Metadata.csv")
)

# 2. Create filtered DGEList
cat("Creating filtered DGEList...\n")
dge_list <- create_filtered_dgelist(
  expr_data$counts,
  expr_data$sample_info,
  MIN_LIB_SIZE
)

cat("Retained", ncol(dge_list), "samples after quality filtering\n")

# 3. Perform MDS analysis
cat("Performing MDS analysis...\n")
mds_data <- perform_mds_analysis(dge_list, OUTPUT_DIR)

# 4. Perform differential expression analysis
cat("Performing differential expression analysis...\n")
design_formula <- ~ Plot_Column + Plot_Row + leaf_tissue + Treatment * Genotype
de_results <- perform_de_analysis(dge_list, design_formula)

# Save normalized expression data
saveRDS(de_results$voom$E, file.path(OUTPUT_DIR, "normalized_expression_logCPM.rda"))
saveRDS(de_results$voom, file.path(OUTPUT_DIR, "normalized_expression_voom_object.rda"))

# 5. Load gene annotation
cat("Loading gene annotation...\n")
gene_symbol <- read.table(
  file.path(DATA_DIR, "gene_symbol.tab"),
  quote = "", header = TRUE, sep = "\t", na.strings = ""
)

pannzer_path <- "~/Desktop/PANNZER/PANNZER_DESC.tab"
if (file.exists(pannzer_path)) {
  pannzer <- read.table(
    pannzer_path,
    quote = "", header = TRUE, sep = "\t", na.strings = ""
  ) %>%
    group_by(gene_model) %>%
    slice(1) %>%
    select(gene_model, desc)
} else {
  warning("PANNZER annotation file not found, proceeding without descriptions")
  pannzer <- data.frame(gene_model = character(), desc = character())
}

gene_annotation <- gene_symbol %>%
  left_join(pannzer, by = "gene_model")

# 6. Extract and classify results
cat("Extracting differential expression results...\n")
predictors_of_interest <- c(
  "leaf_tissue",
  "TreatmentLow_P",
  "GenotypeINV4",
  "TreatmentLow_P:GenotypeINV4"
)

effects <- get_DE_effects(
  de_results$fit,
  predictors_of_interest,
  gene_annotation
)

# 7. Load genomic coordinates
cat("Loading genomic coordinates...\n")
genomic_data <- load_genomic_coordinates(
  "~/ref/zea/Zea_mays.Zm-B73-REFERENCE-NAM-5.0.59.chr.gff3"
)

# 7. Apply statistical classifications
cat("Applying statistical classifications...\n")
logfc_thresholds <- c(
  "general" = LOGFC_THRESHOLD_GENERAL,
  "leaf_tissue" = LOGFC_THRESHOLD_TISSUE
)

effects <- classify_de_results(effects, logfc_thresholds, FDR_THRESHOLD)

# 8. Merge effects with gene annotation
effects <- effects %>%
  mutate(desc_merged = coalesce(locus_name, desc)) %>%
  select(predictor, regulation, Response, locus_symbol, desc_merged, everything()) %>%
  rename(gene = Response) %>%
  inner_join(
    genomic_data$genes %>%
      as.data.frame() %>%
      select(
        gene = gene_id,
        CHR = seqnames,
        BP = start
      ) %>%
      mutate(CHR = as.integer(as.character(CHR)))
  )


# 9. Add genomic location classifications
effects <- effects %>%
  mutate(
    in_shared = CHR == 4 & 
      (BP >= genomic_data$introgression_start) & 
      (BP <= genomic_data$introgression_end),
    in_inv4m = CHR == 4 & 
      (BP >= genomic_data$inv4m_start) & 
      (BP <= genomic_data$inv4m_end),
    in_flanking = in_shared & !in_inv4m,
    region = case_when(
      in_inv4m ~ "inv4m",
      in_flanking ~ "flanking",
      TRUE ~ "outside"
    )
  )

# 10. Add Mahalanobis outliers
cat("Adding Mahalanobis outlier detection...\n")
effects_with_outliers <- add_mahalanobis_outliers(effects)

# 11. Generate summary outputs
cat("Generating summary outputs...\n")

# DEG summary
deg_summary <- effects_with_outliers %>%
  filter(is_significant & regulation != "Unregulated") %>%
  mutate(
    in_cis = gene %in% genomic_data$shared_introgression_genes,
    in_trans = !in_cis,
    in_Inv4m = gene %in% genomic_data$inv4m_genes
  ) %>%
  select(
    predictor, regulation, gene, locus_symbol,
    description = desc_merged, logFC, neglogP,
    in_cis, in_trans, in_Inv4m
  ) %>%
  arrange(regulation, -neglogP, .by_group = TRUE) %>%
  group_by(predictor, regulation) %>%
  arrange(regulation, -neglogP, .by_group = TRUE)

# Save results
write.csv(
  effects_with_outliers,
  file.path(OUTPUT_DIR, "predictor_effects.csv"),
  row.names = FALSE
)

write.csv(
  deg_summary %>% filter(predictor == "GenotypeINV4"),
  file.path(OUTPUT_DIR, "PSU_Inv4m_DEGs.csv"),
  row.names = FALSE
)

write.csv(
  deg_summary %>% filter(predictor == "Treatment-P"),
  file.path(OUTPUT_DIR, "PSU_phosphorus_DEGs.csv"),
  row.names = FALSE
)

cat("Analysis complete!\n")
cat("Results saved to:\n")
cat("- predictor_effects.csv: All effects with statistical classifications\n")
cat("- PSU_Inv4m_DEGs.csv: Inv4m differential expression results\n")
cat("- PSU_phosphorus_DEGs.csv: Phosphorus differential expression results\n")
cat("- normalized_expression_logCPM.rda: Normalized log-CPM expression data\n")


