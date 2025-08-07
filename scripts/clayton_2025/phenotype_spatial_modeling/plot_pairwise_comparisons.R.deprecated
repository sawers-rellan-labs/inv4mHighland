# Load libraries
library(dplyr)
library(ggplot2)
library(ggpubr)
library(rstatix)

# trait plot ---
# Define custom mean ± 95% CI function
mean_ci_95 <- function(x) {
  m <- mean(x, na.rm = TRUE)
  se <- sd(x, na.rm = TRUE) / sqrt(length(na.omit(x)))
  ci <- qt(0.975, df = length(na.omit(x)) - 1) * se
  return(c(y = m, ymin = m - ci, ymax = m + ci))
}

# Define color palette (adjust if needed)
pal <- c("INV4M" = "purple4", "CTRL" = "gold")

# Define generic plot function
plot_trait <- function(data, trait_name, y_label, pal =pal) {
  # Compute t-tests and CI positioning
  stat <- data %>%
    group_by(donor) %>%
    t_test(reformulate("inv4m", response = trait_name)) %>%
    adjust_pvalue(method = "bonferroni") %>%
    add_significance() %>%
    add_xy_position(x = "inv4m") 
  
  # Create plot
  p <- ggplot(data, aes_string(x = "inv4m", y = trait_name, color = "inv4m", group = "inv4m")) +
    geom_point(position = position_jitter(width = 0.25, height = 0.2),
               shape = 21, fill = "white", size = 2, stroke = 0.3, alpha = 0.8) +
    stat_summary(fun.data = mean_ci_95, geom = "pointrange",
                 position = position_dodge(width = 0.3),
                 size = 0.5) +
    stat_pvalue_manual(
      stat,
      label = "p.adj",
      tip.length = 0.02,
      size = 4
    ) +
    facet_wrap(~donor, scales = "free_x") +
    scale_color_manual(values = pal) +
    labs(
      x = "Genotype",
      y = y_label,
      fill = "Genotype"
    ) +
    ggpubr::theme_classic2(base_size = 20) +
    theme(
      strip.background = element_blank(),
      strip.text = element_text(face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  return(p)
}

# Load your data
field_data <- read.csv("~/Desktop/CLY25_Inv4m.csv")
field_data$inv4m <- factor(field_data$inv4m)

# Generate plots
p_dts <- plot_trait(field_data, "DTS", "Days to Silking")
p_dta <- plot_trait(field_data, "DTA", "Days to Anthesis ")
p_ph  <- plot_trait(field_data, "PH",  "Plant Height [cm]")

ggarrange(p_ph, 
          p_dta + coord_cartesian(ylim = c(74, 85)) + scale_y_continuous(breaks = 74:85), 
          p_dts + coord_cartesian(ylim = c(74, 85))  + scale_y_continuous(breaks = 74:85), 
          ncol = 3,align = "hv"
            )


# Generate plots
p_SL <- plot_trait(field_data, "SL", "Sheath Length")
p_BL <- plot_trait(field_data, "BL", "Blade Length")
p_BW  <- plot_trait(field_data, "BW",  "Blade Width")

# Display them (e.g., in quartz or RStudio)
ggarrange(p_ph, 
          p_dta + coord_cartesian(ylim = c(74, 85)) + scale_y_continuous(breaks = 74:85), 
          p_dts + coord_cartesian(ylim = c(74, 85))  + scale_y_continuous(breaks = 74:85), 
          ncol = 3,align = "hv"
)