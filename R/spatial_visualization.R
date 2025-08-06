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