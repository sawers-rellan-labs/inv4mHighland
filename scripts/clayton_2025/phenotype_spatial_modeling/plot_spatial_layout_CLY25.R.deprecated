# Spatial Distribution of Flowering Time (DTS and DTA)
# Author: [Your name]
# Date: [Current date]

# Load required libraries ----
library(dplyr)
library(ggplot2)
library(ggpubr)

# Data import ----
data <- read.csv("~/Desktop/CLY25_Inv4m.csv")

# Data preprocessing ----
flowering_data <- data %>%
  select(Plant, X_pos, Y_pos, rep, donor, genotype= inv4m_gt,DTS, DTA,PH) %>%
  filter(!is.na(X_pos), !is.na(Y_pos))  # Remove plants without spatial coordinates

noise <- runif(nrow(flowering_data), min = 0.0, max = 0.01)
flowering_data$Y_pos<- flowering_data$Y_pos + noise

# Create DTS spatial plot ----
plot_dts_spatial <- flowering_data %>%
  filter(!is.na(DTS)) %>%
  ggplot(aes(x = X_pos, y = Y_pos)) +
  geom_point(aes(color = DTS), size = 3, alpha = 0.8) +
  scale_color_distiller(palette = "RdYlGn", direction = 1,
                        name = "DTS\n(days)") +
  labs(title = "Days to Silking (DTS)",
       x = "Field X position",
       y = "Field Y position") +
  theme_classic2(base_size = 14) +
  theme(legend.position = "right") +
  coord_equal()

# Create DTA spatial plot ----
plot_dta_spatial <- flowering_data %>%
  filter(!is.na(DTA)) %>%
  ggplot(aes(x = X_pos, y = Y_pos)) +
  geom_point(aes(color = DTA), size = 3, alpha = 0.8) +
  scale_color_distiller(palette = "RdYlGn", direction = 1,
                        name = "DTA\n(days)") +
  labs(title = "Days to Anthesis (DTA)",
       x = "Field X position",
       y = "Field Y position") +
  theme_classic2(base_size = 14) +
  theme(legend.position = "right") +
  coord_equal()

# Create DTA spatial plot ----
plot_ph_spatial <- flowering_data %>%
  filter(!is.na(PH)) %>%
  ggplot(aes(x = X_pos, y = Y_pos)) +
  geom_point(aes(color = PH), size = 3, alpha = 0.8) +
  scale_color_distiller(palette = "RdYlGn", direction = 1,
                        name = "Height [cm]") +
  labs(title = "Plant Height [cm] ",
       x = "Field X position",
       y = "Field Y position") +
  theme_classic2(base_size = 14) +
  theme(legend.position = "right") +
  coord_equal()

# Display plots ----
quartz()
print(plot_dts_spatial)
quartz()
print(plot_dta_spatial)
quartz()
print(plot_ph_spatial)

# Print basic statistics ----
cat("\n=== DTS SUMMARY ===\n")
cat("N plants with DTS data:", sum(!is.na(flowering_data$DTS)), "\n")
cat("Mean DTS:", round(mean(flowering_data$DTS, na.rm = TRUE), 1), "days\n")
cat("Median DTS:", round(median(flowering_data$DTS, na.rm = TRUE), 1), "days\n")
cat("Range DTS:", round(min(flowering_data$DTS, na.rm = TRUE), 1), "-", 
    round(max(flowering_data$DTS, na.rm = TRUE), 1), "days\n")

cat("\n=== DTA SUMMARY ===\n")
cat("N plants with DTA data:", sum(!is.na(flowering_data$DTA)), "\n")
cat("Mean DTA:", round(mean(flowering_data$DTA, na.rm = TRUE), 1), "days\n")
cat("Median DTA:", round(median(flowering_data$DTA, na.rm = TRUE), 1), "days\n")
cat("Range DTA:", round(min(flowering_data$DTA, na.rm = TRUE), 1), "-", 
    round(max(flowering_data$DTA, na.rm = TRUE), 1), "days\n")

