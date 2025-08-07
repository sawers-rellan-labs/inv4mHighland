#' Spatial Visualization Functions
#'
#' Functions for creating spatial plots and visualizations of field experiment data.

#' Create spatial distribution plot
#' 
#' Creates a spatial plot showing the distribution of a phenotype across the field
#' 
#' @param data Data frame containing spatial coordinates and phenotype data
#' @param col_var Column name for the phenotype to plot (unquoted)
#' @param plot_title Title for the plot
#' @param legend_name Name for the legend
#' @param x_lab X-axis label (optional)
#' 
#' @return ggplot object
#' @export
#' @import ggplot2 
#' @import dplyr
#' @importFrom rlang ensym
create_spatial_plot <- function(data, col_var, plot_title, legend_name, x_lab = "") {
  col_sym <- ensym(col_var)
  
  data %>%
    filter(!is.na(!!col_sym)) %>%
    ggplot(aes(x = x, y = y)) +
    geom_point(aes(color = !!col_sym), size = 1.2, alpha = 0.8) +
    scale_color_distiller(palette = "RdYlGn", direction = 1, name = legend_name) +
    labs(title = plot_title, x = x_lab, y = "Field Y position") +
    theme_classic2(base_size = 10) +
    theme(legend.position = "right") +
    coord_equal()
}

#' Create legacy spatial heat map
#' 
#' Legacy function for creating spatial heat maps with blocks
#' 
#' @param data Data frame containing spatial data
#' @param trait Column name for the trait to plot
#' @param title Plot title
#' @param midpoint Midpoint for color scale (optional)
#' 
#' @return ggplot object
#' @export
#' @import ggplot2
plot_trait_spatial <- function(data, trait, title, midpoint = NULL) {
  if (is.null(midpoint)) {
    midpoint <- median(data[[trait]], na.rm = TRUE)
  }
  
  ggplot(data, aes(x = x, y = y)) +
    geom_point(aes(color = get(trait)), size = 2) +
    scale_color_gradient2(low = "blue", mid = "white", high = "red",
                          midpoint = midpoint) +
    facet_wrap(~block, labeller = labeller(block = function(x) paste("Block", x))) +
    labs(title = title, color = trait, x = "X Position", y = "Y Position") +
    coord_equal() +
    theme_minimal()
}

#' Show spatial distribution for multiple traits
#' 
#' Creates spatial distribution plots for multiple traits and displays them in groups of three.
#' By default shows only the first 9 traits. For more than 9 traits, user must explicitly 
#' request plotting up to 30 traits maximum.
#' 
#' @param data Dataset containing the phenotypes and spatial coordinates (x, y)
#' @param traits Vector of trait names to plot
#' @param force_all Logical, whether to force plotting more than 9 traits (default: FALSE)
#' @param max_traits Maximum number of traits to plot when force_all=TRUE (default: 30)
#' 
#' @return Invisibly returns list of generated plots
#' @export
#' @import ggplot2
#' @import dplyr
#' @import ggpubr
#' @importFrom rlang sym
show_spatial_distribution <- function(data, traits, force_all = FALSE, max_traits = 30) {
  
  # Check if spatial coordinates are available
  if (!all(c("x", "y") %in% names(data))) {
    cat("No spatial coordinates available for spatial plots\n")
    return(invisible(NULL))
  }
  
  # Validate traits exist in data
  available_traits <- intersect(traits, names(data))
  missing_traits <- setdiff(traits, names(data))
  
  if (length(missing_traits) > 0) {
    cat("Warning: These traits not found in data:", paste(missing_traits, collapse = ", "), "\n")
  }
  
  if (length(available_traits) == 0) {
    cat("No valid traits found for plotting\n")
    return(invisible(NULL))
  }
  
  # Apply trait limits
  if (length(available_traits) > 9 && !force_all) {
    cat("More than 9 traits provided. Showing only first 9 traits.\n")
    cat("Use force_all=TRUE to plot up to", max_traits, "traits.\n")
    available_traits <- available_traits[1:9]
  } else if (length(available_traits) > max_traits) {
    cat("Too many traits. Limiting to first", max_traits, "traits.\n")
    available_traits <- available_traits[1:max_traits]
  }
  
  cat("Creating spatial plots for", length(available_traits), "traits\n")
  
  # Define trait labels for plotting
  trait_labels <- list(
    DTA = "Days to Anthesis", DTS = "Days to Silking", LAE = "Leaves Above Ear",
    PH = "Plant Height", EN = "Ear Number", SL = "Sheath Length", 
    BL = "Blade Length", BW = "Blade Width", EBA = "Estimated Blade Area",
    DTA_GDD = "Anthesis GDD", DTS_GDD = "Silking GDD"
  )
  
  # Create all spatial plots
  spatial_plots <- list()
  for (trait in available_traits) {
    plot_title <- if (trait %in% names(trait_labels)) {
      trait_labels[[trait]]
    } else {
      trait  # Use trait name if no pretty label available
    }
    
    tryCatch({
      spatial_plots[[trait]] <- create_spatial_plot(data, !!sym(trait), plot_title, trait)
      cat("  Created plot for", trait, "\n")
    }, error = function(e) {
      cat("  ERROR creating plot for", trait, ":", e$message, "\n")
    })
  }
  
  if (length(spatial_plots) == 0) {
    cat("No plots could be created\n")
    return(invisible(NULL))
  }
  
  # Display plots in groups of three
  plot_names <- names(spatial_plots)
  group_plots <- list()
  
  for (i in seq(1, length(plot_names), by = 3)) {
    group_plots_current <- list()
    group_names <- plot_names[i:min(i+2, length(plot_names))]
    
    for (j in seq_along(group_names)) {
      plot_name <- group_names[j]
      if (j == 1) {
        group_plots_current[[j]] <- spatial_plots[[plot_name]]
      } else {
        # Remove y-axis labels for plots 2 and 3 in each group
        group_plots_current[[j]] <- spatial_plots[[plot_name]] + 
          theme(axis.title.y = element_blank(), 
                axis.text.y = element_blank(), 
                axis.ticks.y = element_blank())
      }
    }
    
    # Arrange and print the group
    if (length(group_plots_current) > 1) {
      group_arranged <- do.call(ggarrange, c(group_plots_current, list(nrow = 1)))
      print(group_arranged)
      group_plots[[length(group_plots) + 1]] <- group_arranged
    } else {
      print(group_plots_current[[1]])
      group_plots[[length(group_plots) + 1]] <- group_plots_current[[1]]
    }
  }
  
  cat("Displayed", length(spatial_plots), "spatial plots in", length(group_plots), "groups\n")
  
  # Return plots invisibly for further use if needed
  invisible(spatial_plots)
}