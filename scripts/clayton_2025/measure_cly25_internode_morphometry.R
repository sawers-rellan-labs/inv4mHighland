# Complete Maize Height Analysis with Residual Bias Detection
# Author: [Your name]
# Date: [Current date]

# Load required libraries ----
library(dplyr)
library(ggplot2)
library(ggpubr)

# Data import and preprocessing ----
field_data <- read.csv("~/Desktop/CLY25_Inv4m.csv", na.strings=c("","NA",NA)) 

# Include spatial coordinates in metadata
metadata <- field_data %>%
  select(plant:group, X_pos, Y_pos)

phytomere <- read.csv("~/Desktop/CLY25_Inv4m_Internode.csv", na.strings=c("","na","NA",NA)) 
table(phytomere$length)

dissection  <- read.csv("~/Desktop/CLY25_Inv4m_dissection.csv")


# Data cleaning and preparation ----
phytomere_processed <- metadata %>%
  inner_join(phytomere, by = "plant") %>%
  rename(genotype = "inv4m_gt") %>%

  filter(!is.na(length)) %>%
  filter(!(grepl("no tassel",notes))) %>%
  group_by(plant) %>%
  arrange(-internode) %>%
  mutate(from_top = row_number()) 

phytomere_processed <- phytomere_processed  %>%
  filter( donor == "MI21" | donor == "TMEX" & from_top <=10)


for_pca <- phytomere_processed %>%
  select(plant,internode,length) %>% 
  tidyr::pivot_wider(names_from = "internode", values_from = "length", names_prefix = "top") %>%
  arrange(plant)

# Outlier detection ----
is_ear_node_outlier <- with(phytomere_processed,
                            from_top < 7 & from_top >= 3 & length > 20
)

is_tassel_outlier <- with(phytomere_processed,
                          length > 35
)

# Identify outlier plants
outliers <- phytomere_processed$plant[is_ear_node_outlier | is_tassel_outlier]
review <- c(10, 15, 55, 69,96,110,134,166,197,228,231,232,242,259,260,271,297,338,356,420)

# Filter data for analysis ----
 filtered <- phytomere_processed 
#  filter(!plant %in% outliers)

# Height comparison analysis ----
compare_height <- filtered %>% 
  summarise(
    donor = first(donor), 
    genotype = first(genotype),
    sum_nodes = sum(length)
  ) %>%
  distinct() %>%
  inner_join(data %>% select(plant:group, PH, X_pos, Y_pos), by = "plant")


# Calculate regression statistics and differences ----
height_model <- lm(sum_nodes ~ PH, data = compare_height)
compare_height$residuals <- residuals(height_model)
compare_height$difference <- compare_height$sum_nodes - compare_height$PH
compare_height$fitted_values <- fitted(height_model)

r_squared <- summary(height_model)$r.squared
slope <- coef(height_model)[2]
intercept <- coef(height_model)[1]

# Create equation text
equation_text <- sprintf("y = %.2fx + %.2f\nR² = %.3f", slope, intercept, r_squared)

# Generate internode length plots ----

# Plot 1: Line plot by internode position above ground
plot_internode_lines <- filtered %>%
  ggplot(aes(x = internode, y = length, color = genotype, group = plant)) +
  geom_line() +
  xlab("Internode position above ground") +
  scale_color_manual(values = c("gold", "purple4")) +
  scale_fill_manual(values = c("gold", "purple4")) +
  scale_y_continuous(breaks = 5 * (1:6)) +
  scale_x_continuous(breaks = 1:16) +
  theme_classic2(base_size = 20) +
  theme(legend.position = "top")

# Plot 2: Boxplot by internode position above ground
plot_internode_boxplot <- filtered %>%
  ggplot(aes(x = factor(internode), y = length, color = genotype)) +
  geom_boxplot() +
  xlab("Internode position above ground") +
  scale_color_manual(values = c("gold", "purple4")) +
  scale_y_continuous(breaks = 5 * (1:6)) +
  theme_classic2(base_size = 20) +
  theme(legend.position = "top")

# Plot 3: Line plot by position from top with ear nodes highlighted
plot_from_top_lines <- filtered %>%
  arrange(genotype, length) %>%
  ggplot(aes(x = factor(from_top), y = length, color = genotype, group = plant)) +
  geom_line() +
  geom_point(data = filtered %>% filter(has_ear == 1),
             aes(x = factor(from_top), y = length, fill = genotype),
             size = 3, 
             shape = 21,
             color = "white",
             stroke = 1) +
  xlab("Internode position from top") +
  scale_x_discrete(limits = rev) +
  scale_color_manual(values = c("gold", "purple4")) +
  scale_fill_manual(values = c("gold", "purple4"))+
  facet_wrap(. ~ donor, nrow=2) +
  theme_classic2(base_size = 20) +
  theme(legend.position = "top",
        strip.background = element_blank())

quartz(height=12)
plot_from_top_lines

# Plot 4: Boxplot by position from top
plot_from_top_boxplot <- filtered %>%
  ggplot(aes(x = factor(from_top), y = length, color = genotype)) +
  geom_boxplot() +
  xlab("Internode position from top") +
  scale_color_manual(values = c("gold", "purple4")) +
  scale_x_discrete(limits = rev) +
  facet_wrap(. ~ donor, nrow=2) +
  theme_classic2(base_size = 20) +
  theme(legend.position = "top",
        strip.background = element_blank())

quartz(height=12)
plot_from_top_boxplot


# Plot 5: Height comparison with regression line and statistics
plot_height_comparison <- compare_height %>%
  ggplot(aes(x = PH, y = sum_nodes, color = genotype)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed") +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "solid", size = 1) +
  annotate("text", 
           x = min(compare_height$PH, na.rm = TRUE) + 
             0.05 * diff(range(compare_height$PH, na.rm = TRUE)),
           y = max(compare_height$sum_nodes, na.rm = TRUE) - 
             0.1 * diff(range(compare_height$sum_nodes, na.rm = TRUE)),
           label = equation_text, 
           hjust = 0, vjust = 1, size = 5) +
  ggtitle("plant height comparison") +
  xlab("Direct height measurement (cm)") +
  ylab("Sum of internode lengths (cm)") +
  scale_color_manual(values = c("gold", "purple4")) +
  theme_classic2(base_size = 20) +
  theme(legend.position = "top")

# RESIDUAL AND DIFFERENCE BIAS ANALYSIS ----

# Identify plants with negative and positive residuals ----
compare_height$under_residuals <- compare_height$residuals < -10
compare_height$over_residuals <- compare_height$residuals > 10

# Identify plants with negative and positive differences ----
compare_height$under_difference <- compare_height$difference < -10
compare_height$over_difference <- compare_height$difference > 10

# Create combined residual categories
compare_height$residual_category <- case_when(
  compare_height$under_residuals ~ "under",
  compare_height$over_residuals ~ "over", 
  TRUE ~ "within"
)

# Create combined difference categories
compare_height$difference_category <- case_when(
  compare_height$under_difference ~ "under",
  compare_height$over_difference ~ "over", 
  TRUE ~ "within"
)

# Create whole field dataset for spatial visualization ----
whole_field <- data %>%
  select(plant, X_pos, Y_pos) %>%
  right_join(compare_height %>% select(plant, residual_category, difference_category, difference), by = "plant") %>%
  mutate(
    residual_category = ifelse(is.na(residual_category), "not_measured", residual_category),
    difference_category = ifelse(is.na(difference_category), "not_measured", difference_category)
  )

# Create residual scatter plot with both under and over ----
plot_residuals_scatter <- compare_height %>%
  ggplot(aes(x = fitted_values, y = residuals)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_hline(yintercept = -10, linetype = "dashed", color = "blue", alpha = 0.7) +
  geom_hline(yintercept = 10, linetype = "dashed", color = "red", alpha = 0.7) +
  geom_point(aes(color = residual_category, 
                 shape = residual_category,
                 size = residual_category), alpha = 0.7) +
  scale_color_manual(values = c("within" = "gray60", 
                                "under" = "blue", 
                                "over" = "red"),
                     name = "Residuals") +
  scale_shape_manual(values = c("within" = 16, 
                                "under" = 25,  # downward triangle
                                "over" = 24),  # upward triangle
                     name = "Residuals") +
  scale_size_manual(values = c("within" = 2, 
                               "under" = 4,
                               "over" = 4),
                    name = "Residuals") +
  labs(title = "Residuals vs Fitted Values",
       x = "Fitted values (sum of internodes, cm)",
       y = "Residuals (cm)") +
  theme_classic2(base_size = 14) +
  theme(legend.position = "top")

# Create spatial plot showing whole field  ----
plot_spatial_residuals <- whole_field %>%
  ggplot(aes(x = X_pos, y = Y_pos)) +
  geom_point(aes(color = residual_category, 
                 shape = residual_category, 
                 size = residual_category), alpha = 0.7) +
  scale_color_manual(values = c("within" = "gray60", 
                                "under" = "blue", 
                                "over" = "red",
                                "not_measured" = "lightgray"),
                     name = "Residuals") +
  scale_shape_manual(values = c("within" = 16, 
                                "under" = 25,  # downward triangle
                                "over" = 24,   # upward triangle
                                "not_measured" = 1),  # open circle
                     name = "Residuals") +
  scale_size_manual(values = c("within" = 2, 
                               "under" = 4,
                               "over" = 4,
                               "not_measured" = 1),
                    name = "Residuals") +
  labs(title = "Spatial Distribution of Residuals",
       x = "Field X position",
       y = "Field Y position") +
  theme_classic2(base_size = 14) +
  theme(legend.position = "right") +
  coord_equal()

# Create difference scatter plot ----
plot_difference_scatter <- compare_height %>%
  ggplot(aes(x = PH, y = difference)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_hline(yintercept = -10, linetype = "dashed", color = "blue", alpha = 0.7) +
  geom_hline(yintercept = 10, linetype = "dashed", color = "red", alpha = 0.7) +
  geom_point(aes(color = difference_category, 
                 shape = difference_category,
                 size = difference_category), alpha = 0.7) +
  scale_color_manual(values = c("within" = "gray60", 
                                "under" = "blue", 
                                "over" = "red"),
                     name = "Difference") +
  scale_shape_manual(values = c("within" = 16, 
                                "under" = 25,  # downward triangle
                                "over" = 24),  # upward triangle
                     name = "Difference") +
  scale_size_manual(values = c("within" = 2, 
                               "under" = 4,
                               "over" = 4),
                    name = "Difference") +
  labs(title = "Difference vs Direct Height",
       x = "Direct height measurement (cm)",
       y = "Difference (sum_nodes - PH, cm)") +
  theme_classic2(base_size = 14) +
  theme(legend.position = "top")

# Create spatial plot showing whole field with differences ----
plot_spatial_difference <- whole_field %>%
  ggplot(aes(x = X_pos, y = Y_pos)) +
  geom_point(aes(color = difference_category, 
                 shape = difference_category, 
                 size = difference_category), alpha = 0.7) +
  scale_color_manual(values = c("within" = "gray60", 
                                "under" = "blue", 
                                "over" = "red",
                                "not_measured" = "lightgray"),
                     name = "Difference") +
  scale_shape_manual(values = c("within" = 16, 
                                "under" = 25,  # downward triangle
                                "over" = 24,   # upward triangle
                                "not_measured" = 1),  # open circle
                     name = "Difference") +
  scale_size_manual(values = c("within" = 2, 
                               "under" = 4,
                               "over" = 4,
                               "not_measured" = 1),
                    name = "Difference") +
  labs(title = "Spatial Distribution of Differences",
       x = "Field X position",
       y = "Field Y position") +
  theme_classic2(base_size = 14) +
  theme(legend.position = "right") +
  coord_equal()
data %>%arrange(-DTA)
# Field layout plot with plant numbers ----
plot_field_layout <- data %>%
  ggplot(aes(x = X_pos, y = Y_pos)) +
  geom_text(aes(label = row_id), size = 2, vjust = -0.8) +
  labs(title = "Field Layout",
       x = "Field X position",
       y = "Field Y position") +
  theme_classic2(base_size = 14) +
  coord_equal()
quartz()
plot_field_layout
compare_height %>%
  arrange(-difference)

# Display all plots ----
print(plot_internode_lines)
print(plot_internode_boxplot)
print(plot_from_top_lines)
print(plot_from_top_boxplot)
print(plot_height_comparison)
print(plot_field_layout)
print(plot_residuals_scatter)
print(plot_spatial_residuals)
print(plot_difference_scatter)
print(plot_spatial_difference)
