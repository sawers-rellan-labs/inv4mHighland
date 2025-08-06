#!/usr/bin/env Rscript

#' Differential Expression Analysis Pipeline
#' 
#' This script performs differential expression analysis on RNA-seq data using edgeR and limma.
#' It includes quality control, normalization, and statistical testing steps.

# Load required libraries
library(edgeR)      # Bioconductor: differential expression analysis
library(limma)      # Bioconductor: linear models for microarray/RNA-seq data
library(rtracklayer) # Bioconductor: genomic annotation handling
library(GenomicRanges) # Bioconductor: genomic ranges operations
library(dplyr)      # CRAN: data manipulation
library(ggplot2)    # CRAN: plotting
library(ggpubr)     # CRAN: publication ready plots
library(ggtext)     # CRAN: formatted text in plots


#' Read and Process Expression Data
#' 
#' @param counts_file Path to counts data CSV
#' @param metadata_file Path to sample metadata CSV
#' @return DGEList object with filtered data
read_expression_data <- function(counts_file, metadata_file) {
  # Read input files
  counts <- read.csv(counts_file)
  sampleInfo <- read.csv(metadata_file)
  
  # Process sample tags
  tag <- sampleInfo$side_tag
  names(tag) <- sampleInfo$library
  
  # Format counts matrix
  genes <- data.frame(gene = counts[,2])
  counts <- as.matrix(counts[,-c(1:2)])
  rownames(counts) <- genes$gene
  
  # Match sample names
  sampleNames <- tag[colnames(counts)]
  colnames(counts) <- sampleNames
  sampleInfo <- sampleInfo[match(sampleNames, sampleInfo$side_tag),]
  
  # Create DGEList object
  y <- DGEList(counts = counts, samples = sampleInfo)
  y$group <- interaction(y$samples$Treatment, y$samples$Genotype)
  
  # Filter low expression genes
  keep <- filterByExpr(y, group = y$group)
  y_filtered <- y[keep,]
  
  # Filter low quality samples
  y_filtered$samples$lowCount <- y_filtered$samples$lib.size < 2e7
  y_filtered <- y_filtered[,!y_filtered$samples$lowCount]
  
  return(y_filtered)
}

#' Perform MDS Analysis
#' 
#' @param y_filtered Filtered DGEList object
#' @return List containing MDS coordinates and plot data
perform_mds_analysis <- function(y_filtered) {
  mds <- plotMDS(y_filtered, pch=21, plot=FALSE)
  mds2 <- plotMDS(y_filtered, pch=21, dim.plot = c(3,4), plot=FALSE)
  
  plot_data <- y_filtered$samples
  plot_data$x <- mds2$x
  plot_data$y <- mds2$y
  plot_data$Treatment <- factor(plot_data$Treatment)
  plot_data$Genotype <- factor(plot_data$Genotype)
  
  return(list(mds=mds, mds2=mds2, plot_data=plot_data))
}

#' Run Differential Expression Analysis
#' 
#' @param y_filtered Normalized DGEList object
#' @param design Design matrix for the model
#' @return List containing voom and eBayes results
run_de_analysis <- function(y_filtered, metadata) {
  # Normalize
  y_filtered <- calcNormFactors(y_filtered)
  
  # Create design matrix
  design <- model.matrix(
    ~ Plot_Column + Plot_Row + leaf_tissue + Treatment*Genotype,
    metadata
  )
  
  # Voom transformation
  voom_result <- voom(y_filtered, design=design, plot=FALSE)
  
  # Fit models
  fit <- lmFit(voom_result)
  eb_fit <- eBayes(fit, robust=TRUE)
  
  return(list(voom=voom_result, fit=eb_fit))
}

#' Process DE Results
#' 
#' @param eb_fit eBayes fit object
#' @param coefficients Coefficients of interest
#' @param gene_info Gene annotation information
#' @return Data frame of processed results
process_de_results <- function(eb_fit, coefficients, gene_info) {
  results <- list()
  
  for(coef in coefficients) {
    # Get results table
    r <- cbind(
      topTable(eb_fit, coef=coef, sort.by='none', n=Inf),
      data.frame(predictor=coef)
    ) %>% tibble::rownames_to_column("Response")
    
    # Calculate confidence intervals
    cr <- qt(0.975, eb_fit$df.residual + eb_fit$df.prior) * 
      eb_fit$stdev.unscaled[,coef] * sqrt(eb_fit$s2.post)
    r$upper <- r$logFC + cr
    r$lower <- r$logFC - cr
    
    results[[coef]] <- r
  }
  
  # Combine and annotate results
  effects <- results %>% 
    dplyr::bind_rows() %>%
    mutate(
      P = adj.P.Val,
      neglogP = -log10(adj.P.Val),
      is_significant = adj.P.Val < 0.05,
      regulation = case_when(
        is_significant & logFC > 2 ~ "Upregulated",
        is_significant & logFC < -2 ~ "Downregulated",
        TRUE ~ "Unregulated"
      )
    ) %>%
    left_join(gene_info, by=c("Response"="gene_id"))
  
  return(effects)
}

#' Create MDS Plot by Treatment and Leaf
#' 
#' @param mds_data MDS coordinates and metadata
#' @param var_explained Variance explained by each dimension
#' @return ggplot object
plot_mds_treatment_leaf <- function(mds_data, var_explained) {
  mds_data %>%
    mutate(leaf = factor(leaf_tissue)) %>%
    ggplot(aes(x=x, y=y)) + 
    xlab(paste0("dim1 (", round(100*var_explained[1],), "%)")) +
    ylab(paste0("dim2 (", round(100*var_explained[2],), "%)")) +
    geom_point(aes(fill = leaf, shape = Treatment), size=4) +
    scale_fill_viridis_d() +
    scale_shape_manual(values=c(24,21)) +
    guides(
      shape = guide_legend(
        title = "Treatment",
        order = 1,
        override.aes = list(size=7)
      ),
      fill = guide_legend(
        title = "Leaf", 
        order = 2,
        override.aes = list(geom = "point", shape = 22, size=7)
      )
    ) +
    ggpubr::theme_classic2(base_size = 25) +
    theme(
      legend.box = "horizontal",
      legend.spacing = unit(0,"line"),
      legend.box.spacing = unit(0, "in"),
      legend.position = c(0.75,0.17)
    )
}

#' Create MDS Plot by Genotype and Treatment
#' 
#' @param mds_data MDS coordinates and metadata
#' @param var_explained Variance explained by each dimension
#' @param dim_plot Vector of dimensions to plot (e.g., c(3,4))
#' @return ggplot object
plot_mds_genotype <- function(mds_data, var_explained, dim_plot = c(3,4)) {
  labels <- c("CTRL", "*Inv4m*")
  names(labels) <- c("CTRL", "INV4")
  
  mds_data %>%
    ggplot(aes(x=x, y=y, fill=Genotype, shape = Treatment)) + 
    xlab(paste0("dim", dim_plot[1], " (", round(100*var_explained[dim_plot[1]],), "%)")) +
    ylab(paste0("dim", dim_plot[2], " (", round(100*var_explained[dim_plot[2]],), "%)")) +
    geom_point(size=4) +
    scale_fill_viridis_d(direction = -1, labels = labels) +
    scale_shape_manual(values=c(24,21)) +
    guides(
      shape = "none",
      fill = guide_legend(
        title = "Genotype", 
        order = 2,
        override.aes = list(geom = "point", shape = 22, size=7, reverse = TRUE)
      )
    ) +
    ggpubr::theme_classic2(base_size = 25) +
    theme(
      legend.position = c(0.89,0.9),
      legend.text = ggtext::element_markdown(),
      legend.spacing = unit(0,"line"),
      legend.box.spacing = unit(0, "line")
    )
}

#' Create Multiple MDS Plots by Different Variables
#' 
#' @param mds_data MDS coordinates and metadata
#' @param output_file Optional file path to save plots
#' @return List of ggplot objects
plot_mds_variables <- function(mds_data, output_file = NULL) {
  plots <- list()
  
  # Treatment plot
  plots$treatment <- ggplot(mds_data, aes(x=x, y=y)) + 
    geom_point(aes(color = Treatment)) +
    ggpubr::theme_classic2(base_size = 25) +
    theme(
      legend.box = "horizontal",
      legend.spacing = unit(0,"line"),
      legend.box.spacing = unit(0, "in"),
      legend.position = c(0.75,0.2)
    )
  
  # Row and Treatment plot
  plots$row_treatment <- ggplot(mds_data, aes(x=x, y=y)) + 
    geom_point(aes(color = row, shape = Treatment)) +
    ggpubr::theme_classic2(base_size = 25) +
    theme(
      legend.box = "horizontal",
      legend.spacing = unit(0,"line"),
      legend.box.spacing = unit(0, "in"),
      legend.position = c(0.75,0.2)
    )
  
  # Time plot
  plots$time <- ggplot(mds_data, aes(x=x, y=y)) + 
    geom_point(aes(color = decimal_time)) +
    ggpubr::theme_classic2(base_size = 25) +
    theme(
      legend.box = "horizontal",
      legend.spacing = unit(0,"line"),
      legend.box.spacing = unit(0, "in"),
      legend.position = c(0.75,0.2)
    )
  
  # Collector plot
  plots$collector <- ggplot(mds_data, aes(x=x, y=y)) + 
    geom_point(aes(color = COLLECTOR)) +
    ggpubr::theme_classic2(base_size = 25) +
    theme(
      legend.box = "horizontal",
      legend.spacing = unit(0,"line"),
      legend.box.spacing = unit(0, "in"),
      legend.position = c(0.75,0.2)
    )
  
  # Genotype plot
  plots$genotype <- ggplot(mds_data, aes(x=x, y=y)) + 
    geom_point(aes(color = Genotype)) +
    ggpubr::theme_classic2(base_size = 25) +
    theme(
      legend.box = "horizontal",
      legend.spacing = unit(0,"line"),
      legend.box.spacing = unit(0, "in"),
      legend.position = c(0.75,0.2)
    )
  
  if (!is.null(output_file)) {
    pdf(file = output_file)
    for (plot in plots) {
      print(plot)
    }
    dev.off()
  }
  
  return(plots)
}

#' Enhanced MDS Analysis Function
#' 
#' @param y_filtered Filtered DGEList object
#' @param output_dir Directory to save plots
#' @return List containing MDS coordinates and plots
perform_mds_analysis <- function(y_filtered, output_dir = NULL) {
  # Calculate MDS coordinates
  mds <- plotMDS(y_filtered, pch=21, plot=FALSE)
  mds2 <- plotMDS(y_filtered, pch=21, dim.plot = c(3,4), plot=FALSE)
  
  # Prepare plot data
  plot_data <- y_filtered$samples
  plot_data$x <- mds$x
  plot_data$y <- mds$y
  plot_data$Treatment <- factor(plot_data$Treatment)
  plot_data$Genotype <- factor(plot_data$Genotype)
  
  # Create plots
  plots <- list()
  plots$treatment_leaf <- plot_mds_treatment_leaf(plot_data, mds$var.explained)
  plots$genotype <- plot_mds_genotype(plot_data, mds2$var.explained)
  plots$variables <- plot_mds_variables(plot_data)
  
  # Save plots if output directory is provided
  if (!is.null(output_dir)) {
    # Save treatment and leaf plot
    ggsave(
      filename = file.path(output_dir, "mds_treatment_leaf.png"),
      plot = plots$treatment_leaf,
      width = 7,
      height = 7,
      units = "in",
      dpi = 300
    )
    
    # Save genotype plot
    ggsave(
      filename = file.path(output_dir, "mds_genotype.png"),
      plot = plots$genotype,
      width = 7,
      height = 7,
      units = "in",
      dpi = 300
    )
    
    # Save variable plots
    pdf(file = file.path(output_dir, "mds_variables.pdf"))
    for (plot in plots$variables) {
      print(plot)
    }
    dev.off()
  }
  
  return(list(
    mds = mds,
    mds2 = mds2,
    plot_data = plot_data,
    plots = plots
  ))
}

# Update main function to include MDS plotting
main <- function(counts_file, metadata_file, gff_file, output_dir) {
  # Create output directory if it doesn't exist
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Read and process data
  y_filtered <- read_expression_data(counts_file, metadata_file)
  
  # Perform MDS analysis with plots
  mds_results <- perform_mds_analysis(y_filtered, output_dir)
  summary(mds_results)
  # Run DE analysis
  de_results <- run_de_analysis(y_filtered, y_filtered$samples)
  
  # Define coefficients of interest
  coefficients <- c("leaf_tissue", "TreatmentLow_P", "GenotypeINV4", "TreatmentLow_P:GenotypeINV4")
  
  # Get gene annotations
  gene_info <- import(gff_file) %>%
    subset(type=="gene") %>%
    as.data.frame()
  
  # Process results
  effects <- process_de_results(de_results$fit, coefficients, gene_info)
  
  # Save results
  write.csv(effects, file.path(output_dir, "DEG_results.csv"), row.names=FALSE)
  
  return(list(
    effects = effects,
    mds_results = mds_results
  ))
}

quartz()
mds_results$plots$genotype

# Example usage:
if (FALSE) {
  results <- main(
    counts_file = "data/inv4mRNAseq_gene_sample_exp.csv",
    metadata_file = "data/PSU-PHO22_Metadata.csv",
    gff_file = "~/ref/zea/Zea_mays.Zm-B73-REFERENCE-NAM-5.0.59.chr.gff3",
    output_dir = "results"
  )
}