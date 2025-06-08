rm(list=ls())


# Load necessary libraries
library(ggplot2)
library(dplyr)
library(gridExtra)
library(tidyverse)
library(patchwork)
library(viridis)
library(cowplot)
# library(betareg)
library(lme4)

# Set working directory
# setwd('./fullSummaries')


# Load datasets
full_results.05 = read.csv("./fullSummaries/Full_Results_Summary0.05.csv")
full_results.1 = read.csv("./fullSummaries/Full_Results_Summary0.1.csv")
full_results.2 = read.csv("./fullSummaries/Full_Results_Summary0.2.csv")


# Combine datasets
# bias_results <- bind_rows(bias_results.05.3, bias_results.1.3, bias_results.2.3)
full_results <- bind_rows(
  full_results.05 %>% mutate(Density = "0.05"),
  full_results.1 %>% mutate(Density = "0.1"),
  full_results.2 %>% mutate(Density = "0.2")
)


# Adjust data type
full_results$surveyLength = as.factor(full_results$surveyLength)
full_results$radius = as.factor(full_results$radius)
full_results$nSurveys = as.factor(full_results$nSurveys)
full_results$density = as.factor(full_results$density)

# Rename columns
colnames(full_results)[colnames(full_results) == "ExpDec_DailyBias"] <- "ExponentialDecay_DailyBias"
colnames(full_results)[colnames(full_results) == "ExpDec_SeasonBias"] <- "ExponentialDecay_SeasonBias"
colnames(full_results)[colnames(full_results) == "ExpDec_InstBias"] <- "ExponentialDecay_InstBias"
colnames(full_results)[colnames(full_results) == "Initial_InstantBias"] <- "Initial_InstBias"


# Trim results
full_results_trimmed <- full_results %>%
  select(simulation, surveyLength, intervalLength, nSurveys, radius, density, 
         Initial_InstBias, Initial_DailyBias, Initial_SeasonBias, 
         ExponentialDecay_success, ExponentialDecay_DailyBias, ExponentialDecay_SeasonBias, ExponentialDecay_InstBias,
         Gompertz_success, Gompertz_DailyBias, Gompertz_SeasonBias, Gompertz_InstBias,
         MichaelisMenten_success, MichaelisMenten_DailyBias, MichaelisMenten_SeasonBias, MichaelisMenten_InstBias,
         LogisticGrowth_success, LogisticGrowth_DailyBias, LogisticGrowth_SeasonBias, LogisticGrowth_InstBias,
         AsymptoticRegression_success, AsymptoticRegression_DailyBias, AsymptoticRegression_SeasonBias, AsymptoticRegression_InstBias,
         UpperHingeThreshold_success, UpperHingeThreshold_threshold, UpperHingeThreshold_DailyBias, UpperHingeThreshold_SeasonBias, UpperHingeThreshold_InstBias)

# Pivot longer to make each scenario a row
full_results_long <- full_results_trimmed %>%
  pivot_longer(
    cols = c(
      Initial_InstBias, Initial_DailyBias, Initial_SeasonBias,
      ExponentialDecay_success, ExponentialDecay_DailyBias, ExponentialDecay_SeasonBias, ExponentialDecay_InstBias,
      Gompertz_success, Gompertz_DailyBias, Gompertz_SeasonBias, Gompertz_InstBias,
      MichaelisMenten_success, MichaelisMenten_DailyBias, MichaelisMenten_SeasonBias, MichaelisMenten_InstBias,
      LogisticGrowth_success, LogisticGrowth_DailyBias, LogisticGrowth_SeasonBias, LogisticGrowth_InstBias,
      AsymptoticRegression_success, AsymptoticRegression_DailyBias, AsymptoticRegression_SeasonBias, AsymptoticRegression_InstBias,
      UpperHingeThreshold_success, UpperHingeThreshold_threshold, UpperHingeThreshold_DailyBias, UpperHingeThreshold_SeasonBias, UpperHingeThreshold_InstBias
    ),
    names_to = c("Model", ".value"),  # Split names into "Model" and corresponding value types
    names_sep = "_"  # Separator between model name and variable type
  )

full_results_long$absDailyBias = abs(full_results_long$DailyBias)
full_results_long = full_results_long[full_results_long$Model != "UpperHingeThreshold", ]



################################################################################
# Figure 4
################################################################################
full_results_psiVariation <- full_results |>
  dplyr::filter(surveyLength %in% c(5, 10, 15, 20, 30)) |>
  dplyr::mutate(standardEstimate = plogis(psi))

density_levels <- sort(unique(full_results_psiVariation$density))
density_colors <- setNames(
  c("#0077FF",  
    "#00C853",  
    "#FFFF00"   
  ),
  as.character(density_levels)
)


# Update the radius labels for facets
radius_labels <- c(
  "50" = "50 m radius",
  "100" = "100 m radius"
)

psi_violin_plot <- full_results_psiVariation |>
  ggplot(aes(x = factor(density), y = standardEstimate, fill = factor(density))) +
  geom_violin(
    alpha = 0.5,
    scale = "width",
    trim = TRUE,
    position = position_identity()
  ) +
  # Add IQR as error bars
  stat_summary(
    fun.min = ~ quantile(.x, 0.025),
    fun.max = ~ quantile(.x, 0.975),
    fun = median,
    geom = "errorbar",
    color = "black",
    width = 0.35,
    size = 1.5
  ) +
  # Add median point
  stat_summary(
    fun = median,
    geom = "point",
    shape = 21,
    color = "black",
    fill = "black",
    size = 3,
    stroke = 1
  ) +
  facet_wrap(~ radius, labeller = as_labeller(radius_labels)) +
  scale_fill_manual(values = density_colors) +
  labs(
    x = "Density (thrush / ha)",
    y = "Occupancy estimate",
    title = "Variation in traditional occupancy estimates"
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(size = 12)
    # Grid lines are now back by default with theme_minimal()
  )

print(psi_violin_plot)


# ------------------------------------------------------------------------------
# Compute variance, 95% confidence limits, and range of standard occupancy estimates
# by density and radius
# ------------------------------------------------------------------------------
psi_summary <- full_results_psiVariation %>%
  group_by(density, radius) %>%
  summarise(
    variance     = var(standardEstimate, na.rm = TRUE),
    lower_95cl   = quantile(standardEstimate, 0.025, na.rm = TRUE),
    upper_95cl   = quantile(standardEstimate, 0.975, na.rm = TRUE),
    min_value    = min(standardEstimate, na.rm = TRUE),
    max_value    = max(standardEstimate, na.rm = TRUE),
    .groups = "drop"
  )

# ------------------------------------------------------------------------------
# Print the summary
# ------------------------------------------------------------------------------
print(psi_summary)



#-------------------------------------------------------------------------------
# Figure 5
#-------------------------------------------------------------------------------

# Filter data for Initial model
standard_data <- full_results_long %>%
  filter(Model == "Initial")

# Reorder intervalLength factor to ensure correct order
standard_data$intervalLength <- factor(
  standard_data$intervalLength,
  levels = c("none", "24Hours", "10Days")
)

#------------------------------------------------------------------------------
# Convert data to long format for plotting all bias types together
#------------------------------------------------------------------------------
bias_long <- standard_data %>%
  select(density, radius, intervalLength, nSurveys, surveyLength, DailyBias, InstBias, SeasonBias) %>%
  pivot_longer(cols = c(InstBias, DailyBias, SeasonBias), names_to = "BiasType", values_to = "BiasValue") %>%
  mutate(BiasType = factor(BiasType, levels = c("InstBias", "DailyBias", "SeasonBias"),  # Custom order
                           labels = c("Instantaneous", "Daily", "Seasonal")))  # Rename bias types

#------------------------------------------------------------------------------
# Define blending-friendly colors (RGB-mixed)
#------------------------------------------------------------------------------

# Generate 3 distinct colors from the "plasma" palette
bias_colors <- setNames(viridis(3, alpha = 1, option = "viridis"), c("Instantaneous", "Daily", "Seasonal"))



#------------------------------------------------------------------------------
# Create lookup table for x-axis labels
#------------------------------------------------------------------------------
x_axis_labels <- list(
  "density" = "Density (thrush / ha)",
  "radius" = "Radius (m)",
  "intervalLength" = "Time between surveys",
  "nSurveys" = "# of surveys",
  "surveyLength" = "Survey duration (min)"
)

#------------------------------------------------------------------------------
# Function to create violin plots with fully overlapping bias distributions
#------------------------------------------------------------------------------
create_violin_plot <- function(x_axis, title_prefix) {
  ggplot(bias_long, aes_string(
    x = paste0("factor(", x_axis, ")"),  
    y = "BiasValue",
    fill = "BiasType"
  )) +
    geom_violin(alpha = 0.5, scale = "width", trim = TRUE, position = position_identity()) +  
    geom_hline(yintercept = 0, linetype = "dashed", color = "red", size = 1) +  
    scale_fill_manual(values = bias_colors, name = "Bias type") +  
    ggtitle(x_axis_labels[[x_axis]]) +  # Move x-axis label to the top
    labs(y = "Bias") +  # Keep only y-axis label
    theme_bw() +
    theme(
      plot.title = element_text(size = 12, hjust = 0.5),  # Title centered at the top
      axis.title.x = element_blank(),  # Remove default x-axis label at the bottom
      axis.text.x = element_text(angle = 45, hjust = 1),  
      strip.text = element_text(size = 12)  
    )
}

#------------------------------------------------------------------------------
# Generate violin plots for each methodological parameter
#------------------------------------------------------------------------------
plot_surveyLength <- create_violin_plot("surveyLength", "Survey Duration 5")
plot_nSurveys <- create_violin_plot("nSurveys", "# of Surveys 5")
plot_radius <- create_violin_plot("radius", "Radius 5")
plot_density <- create_violin_plot("density", "Density 5")
plot_intervalLength <- create_violin_plot("intervalLength", "Time Between Surveys 5")

#------------------------------------------------------------------------------
# Combine plots into one layout with a larger global title and **single** legend
#------------------------------------------------------------------------------
combined_plot <- (
  (plot_density | plot_radius) /  # Row 1
    (plot_intervalLength | plot_nSurveys) /  # Row 2
    plot_surveyLength  # Row 3
) +
  plot_layout(guides = "collect") +  # Collect and unify legend
  plot_annotation(
    title = "Effects of density & survey protocol on traditional occupancy estimates",
    theme = theme(
      plot.title = element_text(size = 15, face = "bold", hjust = 0.5),  # Centered global title
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.title = element_text(size = 12, hjust = 1)  # Align legend title to the left
    )
  ) & 
  guides(fill = guide_legend(nrow = 1, title.position = "left"))  # Place legend title to the left

# Print the final combined plot
print(combined_plot)



underestimate_proportions <- bias_long %>%
  group_by(BiasType) %>%
  summarise(
    total = n(),
    underestimates = sum(BiasValue < 0, na.rm = TRUE),
    proportion_under = underestimates / total
  )

underestimate_proportions





# Group by BiasType and compute mean and 95% confidence limits (2.5th and 97.5th percentiles)
bias_summary <- bias_long %>%
  group_by(BiasType) %>%
  summarise(
    Mean = mean(BiasValue, na.rm = TRUE),
    Lower_95CL = quantile(BiasValue, 0.025, na.rm = TRUE),
    Upper_95CL = quantile(BiasValue, 0.975, na.rm = TRUE)
  )

# Print summary
print(bias_summary)


# Filter for DailyBias
daily_bias <- bias_long %>% filter(BiasType == "Daily")

# Function to compute mean and 95% confidence limits
compute_stats <- function(data, group_var) {
  data %>%
    group_by({{group_var}}) %>%
    summarise(
      Mean = mean(BiasValue, na.rm = TRUE),
      Lower_95CL = quantile(BiasValue, 0.025, na.rm = TRUE),
      Upper_95CL = quantile(BiasValue, 0.975, na.rm = TRUE),
      .groups = "drop"
    )
}

# Compute statistics for each parameter separately
density_summary <- compute_stats(daily_bias, density)
radius_summary <- compute_stats(daily_bias, radius)
interval_summary <- compute_stats(daily_bias, intervalLength)
nSurveys_summary <- compute_stats(daily_bias, nSurveys)
surveyLength_summary <- compute_stats(daily_bias, surveyLength)

# Print summaries
print(density_summary)
print(radius_summary)
print(interval_summary)
print(nSurveys_summary)
print(surveyLength_summary)



# Figure 6
# ------------------------------------------------------------------------------
# Filter out "Initial" model
# ------------------------------------------------------------------------------
filtered_results <- full_results_long %>%
  filter(Model != "Initial")

# ------------------------------------------------------------------------------
# Compute Success Rate per Model and Density
# ------------------------------------------------------------------------------
success_rates <- filtered_results %>%
  group_by(Model, density) %>%
  summarise(
    success_rate = mean(success, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    success_label = ifelse(is.na(success_rate), "", paste0(round(success_rate * 100, 1), "%"))
  )

# ------------------------------------------------------------------------------
# Compute MSE per Model and Density
# ------------------------------------------------------------------------------
mse_values <- filtered_results %>%
  group_by(Model, density) %>%
  summarise(
    mse = mean(DailyBias^2, na.rm = TRUE) + var(DailyBias, na.rm = TRUE),
    .groups = "drop"
  )

# Merge Success Rate & MSE for labeling
plot_labels <- success_rates %>%
  left_join(mse_values, by = c("Model", "density"))

# ------------------------------------------------------------------------------
# Calculate the statistical summaries for each group
# ------------------------------------------------------------------------------
box_stats <- filtered_results %>%
  group_by(Model, density) %>%
  summarise(
    q3 = quantile(DailyBias, 0.75, na.rm = TRUE),
    iqr = IQR(DailyBias, na.rm = TRUE),
    upper_whisker = min(max(DailyBias, na.rm = TRUE), q3 + 1.5 * iqr),
    .groups = "drop"
  )

# Merge the box stats with the plot labels
plot_labels <- plot_labels %>%
  left_join(box_stats, by = c("Model", "density"))

# ------------------------------------------------------------------------------
# Define custom labels for density levels
# ------------------------------------------------------------------------------
density_labels <- c(
  "0.2" = "0.2 thrush / ha",
  "0.1" = "0.1 thrush / ha",
  "0.05" = "0.05 thrush / ha"
)

# ------------------------------------------------------------------------------
# Final plot
# ------------------------------------------------------------------------------
ggplot(filtered_results, aes(x = Model, y = DailyBias, fill = Model)) +
  geom_boxplot(alpha = 0.7, show.legend = FALSE, outlier.alpha = 0, outlier.size = 0) +
  scale_fill_viridis_d(option = "viridis", name = NULL) +  # Remove legend title "Model"
  scale_x_discrete(labels = c(
    "AsymptoticRegression" = "Asymptotic regression",
    "ExponentialDecay" = "Exponential decay",
    "Gompertz" = "Gompertz",
    "LogisticGrowth" = "Logistic growth",
    "MichaelisMenten" = "Michaelis-Menten"
  )) +
  facet_grid(
    rows = vars(density),
    labeller = labeller(density = density_labels),
    switch = "y"
  ) +
  labs(
    title = "Comparing daily bias between asymptotic models",
    y = "Daily bias",
    shape = ""
  ) +
  ylim(-0.8, 1.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", size = 1.5) +
  geom_text(
    data = plot_labels,
    aes(
      x = Model,
      y = upper_whisker + 0.3,
      label = success_label
    ),
    vjust = 0,
    fontface = "bold",
    color = "black",
    size = 3
  ) +
  geom_point(
    data = plot_labels,
    aes(
      x = Model,
      y = mse,
      shape = "MSE"
    ),
    color = "black",
    fill = "black",
    size = 3
  ) +
  scale_shape_manual(values = c("MSE" = 24)) +
  guides(fill = "none", shape = guide_legend(override.aes = list(size = 4))) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.x = element_blank(),  # ✅ Removed x-axis label
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none",
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 90, face = "plain"),  # Density strips unbolded
  )


# ------------------------------------------------------------------------------
# Compute Mean Daily Bias & 95% Confidence Limits
# ------------------------------------------------------------------------------
bias_summary <- full_results_long %>%
  group_by(Model, density) %>%
  summarise(
    Mean_Bias = mean(DailyBias, na.rm = TRUE),
    Lower_95CL = quantile(DailyBias, 0.025, na.rm = TRUE),
    Upper_95CL = quantile(DailyBias, 0.975, na.rm = TRUE),
    .groups = "drop"
  )

# ------------------------------------------------------------------------------
# Compute MSE per Model and Density
# ------------------------------------------------------------------------------
mse_summary <- full_results_long %>%
  group_by(Model, density) %>%
  summarise(
    MSE = mean(DailyBias^2, na.rm = TRUE) + var(DailyBias, na.rm = TRUE),  # MSE = Mean(Squared Error) + Variance
    .groups = "drop"
  )

# Merge Bias Summary & MSE
full_summary <- bias_summary %>%
  left_join(mse_summary, by = c("Model", "density"))

# Print Summary
print(full_summary)



# ------------------------------------------------------------------------------
# Compute Mean DailyBias & 95% Confidence Limits by Model
# ------------------------------------------------------------------------------
bias_summary <- full_results_long %>%
  group_by(Model) %>%
  summarise(
    Mean_Bias = mean(DailyBias, na.rm = TRUE),
    Lower_95CL = quantile(DailyBias, 0.025, na.rm = TRUE),
    Upper_95CL = quantile(DailyBias, 0.975, na.rm = TRUE),
    .groups = "drop"
  )

# ------------------------------------------------------------------------------
# Compute MSE per Model (ignoring density)
# ------------------------------------------------------------------------------
mse_summary <- full_results_long %>%
  group_by(Model) %>%
  summarise(
    MSE = mean(DailyBias^2, na.rm = TRUE) + var(DailyBias, na.rm = TRUE),
    .groups = "drop"
  )

# ------------------------------------------------------------------------------
# Merge Bias Summary & MSE
# ------------------------------------------------------------------------------
full_summary <- bias_summary %>%
  left_join(mse_summary, by = "Model")

# ------------------------------------------------------------------------------
# Print Summary
# ------------------------------------------------------------------------------
print(full_summary)

mean(filtered_results$DailyBias, na.rm = TRUE)
quantile(filtered_results$DailyBias, probs = c(0.025, 0.975), na.rm = TRUE)




# Figure S1
# ------------------------------------------------------------------------------
# Filter out "Initial" model
# ------------------------------------------------------------------------------
filtered_results <- full_results_long %>%
  filter(Model != "Initial")

# ------------------------------------------------------------------------------
# Compute Success Rate per Model and Density
# ------------------------------------------------------------------------------
success_rates <- filtered_results %>%
  group_by(Model, density) %>%
  summarise(
    success_rate = mean(success, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    success_label = ifelse(is.na(success_rate), "", paste0(round(success_rate * 100, 1), "%"))
  )

# ------------------------------------------------------------------------------
# Compute MSE per Model and Density
# ------------------------------------------------------------------------------
mse_values <- filtered_results %>%
  group_by(Model, density) %>%
  summarise(
    mse = mean(InstBias^2, na.rm = TRUE) + var(InstBias, na.rm = TRUE),
    .groups = "drop"
  )

# Merge Success Rate & MSE for labeling
plot_labels <- success_rates %>%
  left_join(mse_values, by = c("Model", "density"))

# ------------------------------------------------------------------------------
# Calculate the statistical summaries for each group
# ------------------------------------------------------------------------------
box_stats <- filtered_results %>%
  group_by(Model, density) %>%
  summarise(
    q3 = quantile(InstBias, 0.75, na.rm = TRUE),
    iqr = IQR(InstBias, na.rm = TRUE),
    upper_whisker = min(max(InstBias, na.rm = TRUE), q3 + 1.5 * iqr),
    .groups = "drop"
  )

# Merge the box stats with the plot labels
plot_labels <- plot_labels %>%
  left_join(box_stats, by = c("Model", "density"))

# ------------------------------------------------------------------------------
# Define custom labels for density levels
# ------------------------------------------------------------------------------
density_labels <- c(
  "0.2" = "0.2 thrush / ha",
  "0.1" = "0.1 thrush / ha",
  "0.05" = "0.05 thrush / ha"
)

# ------------------------------------------------------------------------------
# Final plot
# ------------------------------------------------------------------------------
ggplot(filtered_results, aes(x = Model, y = InstBias, fill = Model)) +
  geom_boxplot(alpha = 0.7, show.legend = FALSE, outlier.alpha = 0, outlier.size = 0) +
  scale_fill_viridis_d(option = "viridis", name = NULL) +  # Remove legend title "Model"
  scale_x_discrete(labels = c(
    "AsymptoticRegression" = "Asymptotic regression",
    "ExponentialDecay" = "Exponential decay",
    "Gompertz" = "Gompertz",
    "LogisticGrowth" = "Logistic growth",
    "MichaelisMenten" = "Michaelis-Menten"
  )) +
  facet_grid(
    rows = vars(density),
    labeller = labeller(density = density_labels),
    switch = "y"
  ) +
  labs(
    title = "Comparing instantaneous bias between asymptotic models",
    y = "Instantaneous bias",
    shape = ""
  ) +
  ylim(-0.8, 1.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", size = 1.5) +
  geom_text(
    data = plot_labels,
    aes(
      x = Model,
      y = upper_whisker + 0.1,  # ✅ This adapts per facet
      label = success_label
    ),
    vjust = 0,
    fontface = "bold",
    color = "black",
    size = 3
  ) +
  geom_point(
    data = plot_labels,
    aes(
      x = Model,
      y = mse,
      shape = "MSE"
    ),
    color = "black",
    fill = "black",
    size = 3
  ) +
  scale_shape_manual(values = c("MSE" = 24)) +
  guides(fill = "none", shape = guide_legend(override.aes = list(size = 4))) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.x = element_blank(),  # ✅ Removed x-axis label
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none",
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 90, face = "plain"),  # Density strips unbolded
  )



# ------------------------------------------------------------------------------
# Compute Mean InstBias & 95% Confidence Limits
# ------------------------------------------------------------------------------
bias_summary <- full_results_long %>%
  group_by(Model, density) %>%
  summarise(
    Mean_Bias = mean(InstBias, na.rm = TRUE),
    Lower_95CL = quantile(InstBias, 0.025, na.rm = TRUE),
    Upper_95CL = quantile(InstBias, 0.975, na.rm = TRUE),
    .groups = "drop"
  )

# ------------------------------------------------------------------------------
# Compute MSE per Model and Density for InstBias
# ------------------------------------------------------------------------------
mse_summary <- full_results_long %>%
  group_by(Model, density) %>%
  summarise(
    MSE = mean(InstBias^2, na.rm = TRUE) + var(InstBias, na.rm = TRUE),
    .groups = "drop"
  )

# Merge Bias Summary & MSE
full_summary <- bias_summary %>%
  left_join(mse_summary, by = c("Model", "density"))

# Print Summary
print(full_summary)


mean(filtered_results$InstBias, na.rm = TRUE)
quantile(filtered_results$InstBias, probs = c(0.025, 0.975), na.rm = TRUE)



instbias_summary_by_model <- filtered_results %>%
  group_by(Model) %>%
  summarise(
    mean_bias = mean(InstBias, na.rm = TRUE),
    lower_95cl = quantile(InstBias, 0.025, na.rm = TRUE),
    upper_95cl = quantile(InstBias, 0.975, na.rm = TRUE),
    .groups = "drop"
  )

print(instbias_summary_by_model)


# NOW FOR SEASONAL

# Figure S2
# ------------------------------------------------------------------------------
# Filter out "Initial" model
# ------------------------------------------------------------------------------
filtered_results <- full_results_long %>%
  filter(Model != "Initial")

# ------------------------------------------------------------------------------
# Compute Success Rate per Model and Density
# ------------------------------------------------------------------------------
success_rates <- filtered_results %>%
  group_by(Model, density) %>%
  summarise(
    success_rate = mean(success, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    success_label = ifelse(is.na(success_rate), "", paste0(round(success_rate * 100, 1), "%"))
  )

# ------------------------------------------------------------------------------
# Compute MSE per Model and Density
# ------------------------------------------------------------------------------
mse_values <- filtered_results %>%
  group_by(Model, density) %>%
  summarise(
    mse = mean(SeasonBias^2, na.rm = TRUE) + var(SeasonBias, na.rm = TRUE),
    .groups = "drop"
  )

# Merge Success Rate & MSE for labeling
plot_labels <- success_rates %>%
  left_join(mse_values, by = c("Model", "density"))

# ------------------------------------------------------------------------------
# Calculate the statistical summaries for each group
# ------------------------------------------------------------------------------
box_stats <- filtered_results %>%
  group_by(Model, density) %>%
  summarise(
    q3 = quantile(SeasonBias, 0.75, na.rm = TRUE),
    iqr = IQR(SeasonBias, na.rm = TRUE),
    upper_whisker = min(max(SeasonBias, na.rm = TRUE), q3 + 1.5 * iqr),
    .groups = "drop"
  )

# Merge the box stats with the plot labels
plot_labels <- plot_labels %>%
  left_join(box_stats, by = c("Model", "density"))

# ------------------------------------------------------------------------------
# Define custom labels for density levels
# ------------------------------------------------------------------------------
density_labels <- c(
  "0.2" = "0.2 thrush / ha",
  "0.1" = "0.1 thrush / ha",
  "0.05" = "0.05 thrush / ha"
)

# ------------------------------------------------------------------------------
# Final plot
# ------------------------------------------------------------------------------
ggplot(filtered_results, aes(x = Model, y = SeasonBias, fill = Model)) +
  geom_boxplot(alpha = 0.7, show.legend = FALSE, outlier.alpha = 0, outlier.size = 0) +
  scale_fill_viridis_d(option = "viridis", name = NULL) +  # Remove legend title "Model"
  scale_x_discrete(labels = c(
    "AsymptoticRegression" = "Asymptotic regression",
    "ExponentialDecay" = "Exponential decay",
    "Gompertz" = "Gompertz",
    "LogisticGrowth" = "Logistic growth",
    "MichaelisMenten" = "Michaelis-Menten"
  )) +
  facet_grid(
    rows = vars(density),
    labeller = labeller(density = density_labels),
    switch = "y"
  ) +
  labs(
    title = "Comparing seasonal bias between asymptotic models",
    y = "Seasonal bias",
    shape = ""
  ) +
  ylim(-0.8, 1.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", size = 1.5) +
  geom_text(
    data = plot_labels,
    aes(
      x = Model,
      y = 0.8,  # ✅ This adapts per facet
      label = success_label
    ),
    vjust = 0,
    fontface = "bold",
    color = "black",
    size = 3
  ) +
  geom_point(
    data = plot_labels,
    aes(
      x = Model,
      y = mse,
      shape = "MSE"
    ),
    color = "black",
    fill = "black",
    size = 3
  ) +
  scale_shape_manual(values = c("MSE" = 24)) +
  guides(fill = "none", shape = guide_legend(override.aes = list(size = 4))) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.x = element_blank(),  # ✅ Removed x-axis label
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none",
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 90, face = "plain"),  # Density strips unbolded
  )

# ------------------------------------------------------------------------------
# Compute Mean SeasonBias & 95% Confidence Limits
# ------------------------------------------------------------------------------
bias_summary <- full_results_long %>%
  group_by(Model, density) %>%
  summarise(
    Mean_Bias = mean(SeasonBias, na.rm = TRUE),
    Lower_95CL = quantile(SeasonBias, 0.025, na.rm = TRUE),
    Upper_95CL = quantile(SeasonBias, 0.975, na.rm = TRUE),
    .groups = "drop"
  )

# ------------------------------------------------------------------------------
# Compute MSE per Model and Density for SeasonBias
# ------------------------------------------------------------------------------
mse_summary <- full_results_long %>%
  group_by(Model, density) %>%
  summarise(
    MSE = mean(SeasonBias^2, na.rm = TRUE) + var(SeasonBias, na.rm = TRUE),
    .groups = "drop"
  )

# Merge Bias Summary & MSE
full_summary <- bias_summary %>%
  left_join(mse_summary, by = c("Model", "density"))

# Print Summary
print(full_summary)

# ------------------------------------------------------------------------------
# Overall Summary: Mean and 95% CL for all SeasonBias values
# ------------------------------------------------------------------------------
mean(filtered_results$SeasonBias, na.rm = TRUE)
quantile(filtered_results$SeasonBias, probs = c(0.025, 0.975), na.rm = TRUE)


# ------------------------------------------------------------------------------
# Compute Mean SeasonBias & 95% Confidence Limits by Model
# ------------------------------------------------------------------------------
bias_summary <- full_results_long %>%
  group_by(Model) %>%
  summarise(
    Mean_Bias = mean(SeasonBias, na.rm = TRUE),
    Lower_95CL = quantile(SeasonBias, 0.025, na.rm = TRUE),
    Upper_95CL = quantile(SeasonBias, 0.975, na.rm = TRUE),
    .groups = "drop"
  )

# ------------------------------------------------------------------------------
# Compute MSE per Model (ignoring density)
# ------------------------------------------------------------------------------
mse_summary <- full_results_long %>%
  group_by(Model) %>%
  summarise(
    MSE = mean(SeasonBias^2, na.rm = TRUE) + var(SeasonBias, na.rm = TRUE),
    .groups = "drop"
  )

# ------------------------------------------------------------------------------
# Merge Bias Summary & MSE
# ------------------------------------------------------------------------------
full_summary <- bias_summary %>%
  left_join(mse_summary, by = "Model")

# ------------------------------------------------------------------------------
# Print Summary
# ------------------------------------------------------------------------------
print(full_summary)



#-------------------------------------------------------------------------------
# Figure 7
#-------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# Filter data for Initial and AR models
# ------------------------------------------------------------------------------
combined_data <- full_results_long %>%
  filter(Model %in% c("Initial", "AsymptoticRegression"))

#------------------------------------------------------------------------------
# Define consistent viridis colors with updated model labels
#------------------------------------------------------------------------------
model_colors <- setNames(
  viridis::viridis(2, alpha = 1, option = "viridis"),
  c("Traditional methods", "Asymptotic regression")
)

# Update combined_data to reflect new Model labels for plot legend
combined_data <- combined_data %>%
  mutate(Model = recode(Model, 
                        "Initial" = "Traditional methods",
                        "AsymptoticRegression" = "Asymptotic regression"
  ))

#------------------------------------------------------------------------------
# Define x-axis labels and ensure correct factor ordering BEFORE MSE computation
#------------------------------------------------------------------------------
x_axis_labels <- c(
  "density" = "Density (thrush / ha)",
  "radius" = "Radius (m)",
  "intervalLength" = "Time between surveys",
  "nSurveys" = "# of surveys",
  "surveyLength" = "Survey duration (min)"
)

# Define correct factor levels for categorical variables (to ensure correct x-axis order)
factor_levels <- list(
  "density" = c("0.05", "0.1", "0.2"),
  "radius" = c("50", "100"),  # Correct order: 50m first
  "intervalLength" = c("none", "24Hours", "10Days"),  # Correct order: None → 24 Hours → 10 Days
  "nSurveys" = NULL,  # Keep default order
  "surveyLength" = c("5", "10", "15", "20", "30")  # Correct order: 5 → 10 → 15 → 20 → 30
)

# Apply correct factor levels where needed
for (var in names(factor_levels)) {
  if (!is.null(factor_levels[[var]])) {
    combined_data[[var]] <- factor(combined_data[[var]], levels = factor_levels[[var]])
  }
}

#------------------------------------------------------------------------------
# Compute MSE independently for each grouping in each subplot
#------------------------------------------------------------------------------
mse_values <- combined_data %>%
  pivot_longer(cols = c(density, radius, intervalLength, nSurveys, surveyLength), 
               names_to = "x_axis", values_to = "x_value") %>%
  drop_na(x_value) %>%
  mutate(x_axis = recode(x_axis, !!!x_axis_labels)) %>%
  group_by(Model, x_axis, x_value) %>%
  summarise(
    mse = mean(DailyBias^2, na.rm = TRUE) + var(DailyBias, na.rm = TRUE),  # MSE calculation
    .groups = "drop"
  )

#------------------------------------------------------------------------------
# Updated function to create box plots (with legend for models & MSE points)
#------------------------------------------------------------------------------
create_box_plot <- function(x_axis, title_prefix) {
  ggplot(combined_data, aes_string(
    x = paste0("factor(", x_axis, ")"), 
    y = "DailyBias",
    fill = "Model"
  )) +
    geom_boxplot(alpha = 0.7, position = position_dodge(width = 0.75), 
                 outlier.alpha = 0, outlier.size = 0) +
    scale_fill_manual(values = model_colors, name = "Model") +
    ggtitle(title_prefix) +  # X-axis label at top
    labs(y = "Daily bias", x = NULL) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red", size = 1) +
    
    # Ensure MSE is placed correctly above its respective box plot
    geom_point(
      data = mse_values %>% filter(x_axis == title_prefix),
      aes(x = factor(x_value), y = mse, shape = "MSE", group = Model), 
      color = "black",
      fill = "black",
      size = 3,
      position = position_dodge(width = 0.75)  # Align to box plots
    ) +
    
    scale_shape_manual(values = c("MSE" = 24), name = "") +  # Black triangle legend
    theme_bw() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.text = element_text(size = 12),
      plot.title = element_text(size = 12, hjust = 0.5),
      legend.position = "bottom",
      legend.direction = "horizontal"
    )
}

# ------------------------------------------------------------------------------
# Generate box plots for each methodological parameter with updated legend labels & MSE
# ------------------------------------------------------------------------------
plot_density <- create_box_plot("density", "Density (thrush / ha)")
plot_radius <- create_box_plot("radius", "Radius (m)")
plot_intervalLength <- create_box_plot("intervalLength", "Time between surveys")
plot_nSurveys <- create_box_plot("nSurveys", "# of surveys")
plot_surveyLength <- create_box_plot("surveyLength", "Survey duration (min)")

# ------------------------------------------------------------------------------
# Combine plots into one layout with global title & single horizontal legend
# ------------------------------------------------------------------------------
combined_plot <- (
  (plot_density | plot_radius) /
    (plot_intervalLength | plot_nSurveys) /
    plot_surveyLength
) +
  patchwork::plot_layout(guides = "collect") +
  patchwork::plot_annotation(
    title = "Daily occupancy bias: traditional methods vs. asymptotic regression",
    theme = theme(
      plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.title = element_text(size = 12)
    )
  )

# Display the plot
combined_plot





#------------------------------------------------------------------------------
# Binomial Regression Model
#------------------------------------------------------------------------------
combined_data <- full_results_long %>%
  filter(Model %in% c("Initial", "AsymptoticRegression"))

# Convert to binary outcome (1 if AR is less biased, 0 otherwise)
binomial_data <- combined_data %>%
  select(density, radius, intervalLength, nSurveys, surveyLength, Model, DailyBias, simulation) %>%  
  pivot_wider(names_from = Model, values_from = DailyBias) %>%  # Spread into "Initial" & "AsymptoticRegression"
  mutate(closer_to_zero = ifelse(abs(AsymptoticRegression) < abs(Initial), 1, 0))  # Binary outcome

# Convert categorical variables to factors and set new reference levels
binomial_data <- binomial_data %>%
  mutate(
    intervalLength = relevel(factor(intervalLength), ref = "none"),  # Set intervalLength reference to 24 Hours
    nSurveys = relevel(factor(nSurveys), ref = "2"),  # Set nSurveys reference to 3
    surveyLength = relevel(factor(surveyLength), ref = "5")  # Set surveyLength reference to 15
  )

binomial_model <- glmer(
  closer_to_zero ~ density + radius + intervalLength + nSurveys + surveyLength + (1 | simulation), 
  data = binomial_data, 
  family = binomial, 
  control = glmerControl(optimizer = "bobyqa")
)

# View the model summary
summary(binomial_model)

# Higher survey effort weakens ability. Lower survey effort makes it more valuable

filtered_data = full_results_long[full_results_long$Model ==c("AsymptoticRegression", "Initial"),]
filtered_data = filtered_data[filtered_data$nSurveys == 2,]
filtered_data = filtered_data[filtered_data$surveyLength == 5,]
filtered_data = filtered_data[filtered_data$intervalLength == "none",]


