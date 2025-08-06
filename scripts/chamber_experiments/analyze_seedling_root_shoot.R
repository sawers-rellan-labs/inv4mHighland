# Load necessary libraries
library(dplyr)
library(ggplot2)
library(ggbeeswarm)
library(ggpubr)
library(rstatix)
library(tidyr)

# Create the dataframe
seedling_data <- read.csv("~/Desktop/inv4m_root_phenotype.csv")

rl <- lm(data=seedling_data, Root_Length_cm ~  Date + Donor*Genotype)

summary(rl)

sl <- lm(data=seedling_data, Shoot_Length_cm ~  Date + Donor*Genotype)

summary(sl)

# =============================================================================
# 1. ROOT LENGTH ANALYSIS
# =============================================================================

cat("=== ROOT LENGTH ANALYSIS ===\n\n")

# Statistical tests for inversion effect within each donor
root_inversion_stats <- data_with_factors %>%
  group_by(Donor) %>%
  t_test(Root_Length_cm ~ Genotype) %>%
  add_significance() %>%
  add_xy_position(x = "Genotype")

# Statistical tests for donor effect within each genotype
root_donor_stats <- data_with_factors %>%
  group_by(Genotype) %>%
  t_test(Root_Length_cm ~ Donor) %>%
  add_significance() %>%
  add_xy_position(x = "Donor")

# Root length plot - Inversion effect by Donor
root_inversion_plot <- ggplot(data_with_factors, aes(x = Genotype, y = Root_Length_cm)) +
  geom_boxplot(
    aes(color = Genotype),
    fill = "white",
    outlier.shape = NA,
    width = 0.5
  ) +
  geom_quasirandom(
    aes(color = Genotype),
    size = 2.5,
    alpha = 0.8,
    width = 0.3
  ) +
  stat_pvalue_manual(
    root_inversion_stats,
    label = "p.signif",
    tip.length = 0.02,
    size = 4
  ) +
  facet_wrap(~ Donor, scales = "free_x") +
  scale_color_manual(values = c("CTRL" = "gold", "Inv4m" = "purple4")) +
  scale_y_continuous(limits = c(0, 12), breaks = seq(0, 12, 2)) +
  labs(
    title = "Root Length: Inversion Effect by Donor",
    x = "Inversion Genotype",
    y = "Root Length (cm)",
    color = "Genotype"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "right",
    strip.background = element_blank(),
    strip.text = element_text(face = "bold")
  )

# Root length plot - Donor effect by Genotype
root_donor_plot <- ggplot(data_with_factors, aes(x = Donor, y = Root_Length_cm)) +
  geom_boxplot(
    aes(color = Donor),
    fill = "white",
    outlier.shape = NA,
    width = 0.5
  ) +
  geom_quasirandom(
    aes(color = Donor),
    size = 2.5,
    alpha = 0.8,
    width = 0.3
  ) +
  stat_pvalue_manual(
    root_donor_stats,
    label = "p.signif",
    tip.length = 0.02,
    size = 4
  ) +
  facet_wrap(~ Genotype, scales = "free_x") +
  scale_color_manual(values = c("TMEX" = "royalblue", "MI21" = "tomato")) +
  scale_y_continuous(limits = c(0, 12), breaks = seq(0, 12, 2)) +
  labs(
    title = "Root Length: Donor Effect by Genotype",
    x = "Donor",
    y = "Root Length (cm)",
    color = "Donor"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "right",
    strip.background = element_blank(),
    strip.text = element_text(face = "bold")
  )

# =============================================================================
# 2. SHOOT LENGTH ANALYSIS
# =============================================================================

cat("=== SHOOT LENGTH ANALYSIS ===\n\n")

# Statistical tests for inversion effect within each donor
shoot_inversion_stats <- data_with_factors %>%
  group_by(Donor) %>%
  t_test(Shoot_Length_cm ~ Genotype) %>%
  add_significance() %>%
  add_xy_position(x = "Genotype")

# Statistical tests for donor effect within each genotype
shoot_donor_stats <- data_with_factors %>%
  group_by(Genotype) %>%
  t_test(Shoot_Length_cm ~ Donor) %>%
  add_significance() %>%
  add_xy_position(x = "Donor")

# Shoot length plot - Inversion effect by Donor
shoot_inversion_plot <- ggplot(data_with_factors, aes(x = Genotype, y = Shoot_Length_cm)) +
  geom_boxplot(
    aes(color = Genotype),
    fill = "white",
    outlier.shape = NA,
    width = 0.5
  ) +
  geom_quasirandom(
    aes(color = Genotype),
    size = 2.5,
    alpha = 0.8,
    width = 0.3
  ) +
  stat_pvalue_manual(
    shoot_inversion_stats,
    label = "p.signif",
    tip.length = 0.02,
    size = 4
  ) +
  facet_wrap(~ Donor, scales = "free_x") +
  scale_color_manual(values = c("CTRL" = "gold", "Inv4m" = "purple4")) +
  scale_y_continuous(limits = c(0, 6), breaks = seq(0, 6, 1)) +
  labs(
    title = "Shoot Length: Inversion Effect by Donor",
    x = "Inversion Genotype",
    y = "Shoot Length (cm)",
    color = "Genotype"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "right",
    strip.background = element_blank(),
    strip.text = element_text(face = "bold")
  )

# Shoot length plot - Donor effect by Genotype
shoot_donor_plot <- ggplot(data_with_factors, aes(x = Donor, y = Shoot_Length_cm)) +
  geom_boxplot(
    aes(color = Donor),
    fill = "white",
    outlier.shape = NA,
    width = 0.5
  ) +
  geom_quasirandom(
    aes(color = Donor),
    size = 2.5,
    alpha = 0.8,
    width = 0.3
  ) +
  stat_pvalue_manual(
    shoot_donor_stats,
    label = "p.signif",
    tip.length = 0.02,
    size = 4
  ) +
  facet_wrap(~ Genotype, scales = "free_x") +
  scale_color_manual(values = c("TMEX" = "royalblue", "MI21" = "tomato")) +
  scale_y_continuous(limits = c(0, 6), breaks = seq(0, 6, 1)) +
  labs(
    title = "Shoot Length: Donor Effect by Genotype",
    x = "Donor",
    y = "Shoot Length (cm)",
    color = "Donor"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "right",
    strip.background = element_blank(),
    strip.text = element_text(face = "bold")
  )

# =============================================================================
# 3. DISPLAY PLOTS
# =============================================================================

print("Root Length - Inversion Effect:")
quartz()
print(root_inversion_plot)

print("Root Length - Donor Effect:")
quartz()
print(root_donor_plot)

print("Shoot Length - Inversion Effect:")
print(shoot_inversion_plot)

print("Shoot Length - Donor Effect:")
print(shoot_donor_plot)

# =============================================================================
# 4. STATISTICAL SUMMARIES
# =============================================================================

cat("\n=== STATISTICAL SUMMARIES ===\n\n")

# Linear models with interactions
root_lm <- lm(Root_Length_cm ~ Date + Donor*Genotype, data = data_with_factors)
shoot_lm <- lm(Shoot_Length_cm ~  Date + Donor*Genotype, data = data_with_factors)

cat("ROOT LENGTH MODEL:\n")
print(summary(root_lm))

cat("\nSHOOT LENGTH MODEL:\n")
print(summary(shoot_lm))

# Summary statistics by groups
cat("\n=== DESCRIPTIVE STATISTICS ===\n\n")

summary_stats <- data_with_factors %>%
  group_by(Genotype, Donor) %>%
  summarise(
    Root_Mean = mean(Root_Length_cm),
    Root_SD = sd(Root_Length_cm),
    Root_SE = sd(Root_Length_cm) / sqrt(n()),
    Shoot_Mean = mean(Shoot_Length_cm),
    Shoot_SD = sd(Shoot_Length_cm),
    Shoot_SE = sd(Shoot_Length_cm) / sqrt(n()),
    n = n(),
    .groups = "drop"
  )

print(summary_stats)

# Effect sizes for inversion effect within each donor
cat("\n=== EFFECT SIZES (Cohen's d) ===\n\n")

root_effect_sizes <- data_with_factors %>%
  group_by(Donor) %>%
  cohens_d(Root_Length_cm ~ Genotype)

shoot_effect_sizes <- data_with_factors %>%
  group_by(Donor) %>%
  cohens_d(Shoot_Length_cm ~ Genotype)

cat("Root Length Effect Sizes:\n")
print(root_effect_sizes)

cat("\nShoot Length Effect Sizes:\n")
print(shoot_effect_sizes)

# Add this code after generating all four plots in your existing script
# =============================================================================
# 5. CREATE MULTIPANEL FIGURE
# =============================================================================

# First, we need to remove individual legends from plots and create a common legend
# We'll extract the legend from one plot and use it for all

# Modify plots to remove legends and titles (we'll add a common legend later)
shoot_inversion_plot_no_leg <- shoot_inversion_plot + 
  theme(legend.position = "none", plot.title = element_blank())

shoot_donor_plot_no_leg <- shoot_donor_plot + 
  theme(legend.position = "none", plot.title = element_blank())

root_inversion_plot_no_leg <- root_inversion_plot + 
  theme(legend.position = "none", plot.title = element_blank())

root_donor_plot_no_leg <- root_donor_plot + 
  theme(legend.position = "none", plot.title = element_blank())

# Extract legends from the original plots
# We need both color schemes, so we'll create a combined legend
legend_genotype <- get_legend(
  shoot_inversion_plot + 
    theme(legend.position = "bottom")
)

legend_donor <- get_legend(
  shoot_donor_plot + 
    theme(legend.position = "bottom")
)

# Create the multipanel figure
# Top row: Shoot analysis (inversion effect left, donor effect right)
# Bottom row: Root analysis (inversion effect left, donor effect right)

multipanel_figure <- ggarrange(
  # Top row - Shoot analysis
  shoot_inversion_plot_no_leg, shoot_donor_plot_no_leg,
  # Bottom row - Root analysis  
  root_inversion_plot_no_leg, root_donor_plot_no_leg,
  ncol = 2, 
  nrow = 2,
  labels = c("A", "B", "C", "D"),
  label.x = 0.02,
  label.y = 0.98,
  font.label = list(size = 16, face = "bold"),
  common.legend = FALSE,  # We'll add legends separately
  align = "hv"
)

# Since we have two different color schemes (Genotype and Donor), 
# we need to add both legends. Let's create a version with legends at the bottom

# For a cleaner approach, let's create plots with appropriate legends
# and use a single common legend approach

# Alternative approach: Create figure with genotype colors for left panels
# and donor colors for right panels

left_panels <- ggarrange(
  shoot_inversion_plot + theme(legend.position = "none", plot.title = element_blank()),
  root_inversion_plot + theme(legend.position = "none", plot.title = element_blank()),
  ncol = 1,
  labels = c("A", "C"),
  label.x = 0.02,
  label.y = 0.98,
  font.label = list(size = 16, face = "bold")
)

right_panels <- ggarrange(
  shoot_donor_plot + theme(legend.position = "none", plot.title = element_blank()),
  root_donor_plot + theme(legend.position = "none", plot.title = element_blank()),
  ncol = 1,
  labels = c("B", "D"),
  label.x = 0.02,
  label.y = 0.98,
  font.label = list(size = 16, face = "bold")
)

# Combine left and right panels
main_figure <- ggarrange(
  left_panels,
  right_panels,
  ncol = 2,
  widths = c(1, 1)
)

# Add legends at the bottom
final_figure <- ggarrange(
  main_figure,
  ggarrange(legend_genotype, legend_donor, ncol = 2),
  nrow = 2,
  heights = c(1, 0.1)
)

# Display the figure
print(final_figure)

# Save the figure
ggsave(
  filename = "seedling_analysis_multipanel.pdf",
  plot = final_figure,
  width = 12,
  height = 10,
  dpi = 300
)

ggsave(
  filename = "seedling_analysis_multipanel.png",
  plot = final_figure,
  width = 12,
  height = 10,
  dpi = 300
)

# Alternative: If you want a single unified legend approach
# Create a combined plot with all data for legend extraction

# Create a dummy plot with both color schemes for a unified legend
dummy_data <- data.frame(
  x = c("CTRL", "Inv4m", "TMEX", "MI21"),
  y = c(1, 1, 1, 1),
  Type = c("Genotype", "Genotype", "Donor", "Donor"),
  Category = c("CTRL", "Inv4m", "TMEX", "MI21")
)

dummy_plot <- ggplot(dummy_data, aes(x = x, y = y, color = Category)) +
  geom_point(size = 3) +
  scale_color_manual(
    values = c(
      "CTRL" = "gold", 
      "Inv4m" = "purple4",
      "TMEX" = "royalblue", 
      "MI21" = "tomato"
    ),
    name = "Factor Levels",
    labels = c(
      "CTRL" = "Control (CTRL)",
      "Inv4m" = "Inversion (Inv4m)",
      "TMEX" = "TMEX Donor",
      "MI21" = "MI21 Donor"
    )
  ) +
  theme_classic() +
  theme(legend.position = "bottom") +
  guides(color = guide_legend(nrow = 2, byrow = TRUE))

combined_legend <- get_legend(dummy_plot)

# Create final figure with unified legend
final_figure_unified <- ggarrange(
  multipanel_figure,
  combined_legend,
  nrow = 2,
  heights = c(1, 0.1)
)

# Display and save the unified legend version
quartz()
print(final_figure_unified)

ggsave(
  filename = "seedling_analysis_multipanel_unified.pdf",
  plot = final_figure_unified,
  width = 12,
  height = 10,
  dpi = 300
)

# Print confirmation
cat("\nMultipanel figures saved as:\n")
cat("- seedling_analysis_multipanel.pdf/png (separate legends)\n")
cat("- seedling_analysis_multipanel_unified.pdf (unified legend)\n")
