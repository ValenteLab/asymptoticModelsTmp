rm(list=ls())



#Density to analyze
densToAnalyze = 0.2
# densToAnalyze = 0.1
# densToAnalyze = 0.05

# Record the start time for the script execution
startTime = Sys.time()

# Load required libraries
require(gdata)      # General-purpose R data functions
require(tidyverse)  # Collection of useful packages (ggplot2, dplyr, etc.)
require(minpack.lm) # Nonlinear least squares optimization
require(forcats)    # Handling factors in tidy data
require(chngpt)     # Change-point regression models

# Initialize lists to store data for different radii
data_50m_list <- list()  # For data with 50m radius
data_100m_list <- list() # For data with 100m radius

# Loop through files to read data
for (i in 1:2) { 
  # Construct the file name
  file_name <- paste0("occResults", i, ".csv")
  
  # Read the file
  tmpRes <- read.csv(file_name)
  
  # Logit-transform occupancy values for analysis
  tmpRes$logitInstOcc <- qlogis(tmpRes$instOcc)
  tmpRes$logitDailyOcc <- qlogis(tmpRes$dailyOcc)
  tmpRes$logitSeasonOcc <- qlogis(tmpRes$seasonOcc)
  
  # Calculate occupancy biases
  tmpRes$Initial_InstantBias <- plogis(tmpRes$psi) - tmpRes$instOcc
  tmpRes$Initial_DailyBias <- plogis(tmpRes$psi) - tmpRes$dailyOcc
  tmpRes$Initial_SeasonBias <- plogis(tmpRes$psi) - tmpRes$seasonOcc
  
  # Store data in respective lists based on radius
  data_50m_list[[i]] <- tmpRes %>% filter(radius == 50)
  data_100m_list[[i]] <- tmpRes %>% filter(radius == 100)
}

# Combine data for 50m radius across all files
combined_50m <- bind_rows(
  lapply(seq_along(data_50m_list), function(i) {
    data_50m_list[[i]] %>%
      rename(simulation = iteration) %>%  # Rename 'iteration' to 'simulation'
      select(simulation, everything())   # Move 'simulation' to the first column
  })
)

# Combine data for 100m radius across all files
combined_100m <- bind_rows(
  lapply(seq_along(data_100m_list), function(i) {
    data_100m_list[[i]] %>%
      rename(simulation = iteration) %>%  # Rename 'iteration' to 'simulation'
      select(simulation, everything())   # Move 'simulation' to the first column
  })
)

# Combine both 50m and 100m datasets
combined_df <- bind_rows(combined_50m, combined_100m)

# Remove the 'X' column (unnecessary in the data)
combined_df <- combined_df %>% select(-X)

# Group data by simulation settings and ensure exactly 30 rows per simulation
grouped_data <- combined_df %>%
  group_by(intervalLength, nSurveys, radius, simulation) %>%
  filter(n() == 90) %>%  # Ensure each simulation has exactly 30 rows
  ungroup()

# Split the grouped data into separate data frames for each combination
grouped_data_frames <- grouped_data %>%
  group_split(intervalLength, nSurveys, radius) %>%  # Split by combinations
  set_names(  # Assign meaningful names to the resulting list elements
    grouped_data %>%
      distinct(intervalLength, nSurveys, radius) %>%
      mutate(
        group_name = paste0(
          "interval_", intervalLength, 
          "_nSurveys_", nSurveys,
          "_radius_", radius
        )
      ) %>%
      pull(group_name)  # Extract names for list
  )

# Model fitting function
fit_model <- function(data, sim_index, model_formula, start_list, model_name, dens) {
  
  data = data %>% 
    filter(density==dens)
  
  tryCatch({
    if (model_name == "UpperHingeThreshold") {
      # Fit the Upper Hinge Threshold Regression Model
      fitted_model <- chngptm(
        formula.1 = psi ~ 1,           # Mean occupancy (constant)
        formula.2 = ~ surveyLength,   # Change point depends on survey length
        data = data,
        family = "gaussian",           # Use Gaussian regression
        type = "M10",                  # Model type
        lb.quantile = 0.01,           # Lower quantile bound
        ub.quantile = 0.99,           # Upper quantile bound
        tol = 1e-04,                  # Convergence tolerance
        maxit = 1000                  # Maximum iterations
      )
      
      # Extract asymptote and threshold
      asymptote <- coef(fitted_model)["(Intercept)"]
      threshold <- fitted_model$chngpt
      
      # Predict psi at surveyLength = 0
      yIntercept <- predict(fitted_model, newdata = data.frame(surveyLength = 0))
      
      # Calculate SSE
      residuals <- data$psi - predict(fitted_model, newdata = data)
      sse <- sum(residuals^2, na.rm = TRUE)
      
      # Add results to the data
      data <- data %>%
        mutate(
          UpperHingeThreshold_success = TRUE,
          UpperHingeThreshold_predicted_asymptote = asymptote,
          UpperHingeThreshold_threshold = threshold,
          UpperHingeThreshold_sse = sse,
          UpperHingeThreshold_predicted_yIntercept = yIntercept
        )
    } else {
      # Fit other nonlinear models
      fitted_model <- nlsLM(
        model_formula,               # Nonlinear formula
        data = data,
        start = start_list,          # Starting parameters
        control = nls.lm.control(maxiter = 500)
      )
      
      
      coefs <- coef(fitted_model)
      names(coefs) <- gsub("\\.\\d+%", "", names(coefs)) # This line removes a bug that causes issues during extraction
      
      # Predict psi at surveyLength = 0
      yIntercept <- predict(fitted_model, newdata = data.frame(surveyLength = 0))
      
      # Extract model-specific asymptote
      asymptote <- switch(
        model_name,
        "ExponentialDecay" = coefs["asymptote"],
        "MichaelisMenten" = coefs["Vmax"] + coefs["c"],  # Modified for nonzero intercept
        "AsymptoticRegression" = coefs["a"] + coefs["b"],
        "LogisticGrowth" = coefs["K"],
        "Gompertz" = coefs["a"],
        NA
      )
      
      residuals <- data$psi - predict(fitted_model, newdata = data)
      sse <- sum(residuals^2, na.rm = TRUE)
      
      # Add results to the data
      data <- data %>%
        mutate(
          !!paste0(model_name, "_success") := TRUE,
          !!paste0(model_name, "_predicted_asymptote") := asymptote,
          !!paste0(model_name, "_sse") := sse,
          !!paste0(model_name, "_predicted_yIntercept") := yIntercept
        )
    }
    return(data)
  }, error = function(e) {
    # Handle errors for all models
    if (model_name == "UpperHingeThreshold") {
      data <- data %>%
        mutate(
          UpperHingeThreshold_success = FALSE,
          UpperHingeThreshold_predicted_asymptote = NA,
          UpperHingeThreshold_threshold = NA,
          UpperHingeThreshold_sse = NA,
          UpperHingeThreshold_predicted_yIntercept = NA
        )
    } else {
      data <- data %>%
        mutate(
          !!paste0(model_name, "_success") := FALSE,
          !!paste0(model_name, "_predicted_asymptote") := NA,
          !!paste0(model_name, "_sse") := NA,
          !!paste0(model_name, "_predicted_yIntercept") := NA
        )
    }
    return(data)
  })
}


# Define truncation sizes for subsets of survey data representing various survey lenths (durations)
truncation_sizes <- c(5, 10, 15, 20, 30)

# Initialize an empty list to store summary results
results_summary_list <- list()

# Loop through each grouped data frame (created earlier from simulations)
for (group_name in names(grouped_data_frames)) {  # Opening for group_name
  group_data <- grouped_data_frames[[group_name]]
  
  # Create truncated subsets
  truncated_data_frames <- lapply(truncation_sizes, function(size) {
    group_data %>%
      group_by(simulation) %>%   # Group by simulation
      slice_head(n = size) %>%   # Select the first `size` rows for each simulation
      ungroup()                  # Remove grouping after slicing
  }) %>%
    set_names(paste0("first_", truncation_sizes, "_entries"))  # Name subsets by truncation size
  
  # Fit models for each truncated subset
  for (subset_name in names(truncated_data_frames)) {  
    data_subset <- truncated_data_frames[[subset_name]]
    
    # Initialize model storage
    model_storage_list <- list(
      ExponentialDecay = vector("list", length = length(unique(data_subset$simulation))),
      MichaelisMenten = vector("list", length = length(unique(data_subset$simulation))),
      AsymptoticRegression = vector("list", length = length(unique(data_subset$simulation))),
      LogisticGrowth = vector("list", length = length(unique(data_subset$simulation))),
      Gompertz = vector("list", length = length(unique(data_subset$simulation))),
      UpperHingeThreshold = vector("list", length = length(unique(data_subset$simulation)))  # Add threshold model
    )
    
    # Used to track progress
    print(paste(group_name, "-", subset_name, "initialized")) 
    
    # Fit models for each simulation
    results <- data_subset %>%
      group_by(simulation) %>%  # Group data by simulation
      group_modify(~ {  # Process each group (simulation-specific)
        simulation_id <- .y$simulation[1]  # Extract the simulation ID
        
        # Fit models
        for (model_name in names(model_storage_list)) {  # Iterate over models
          # Fit the current model using the `fit_model` function
          .x <- fit_model(
            .x,
            simulation_id,
            model_formula = switch(  # Specify model formulas dynamically
              model_name,
              "ExponentialDecay" = psi ~ initial_offset * exp(-decay_rate * surveyLength) + asymptote,
              "MichaelisMenten" = psi ~ Vmax * surveyLength / (Km + surveyLength) + c,
              "AsymptoticRegression" = psi ~ a + b * (1 - exp(-c * surveyLength)),
              "LogisticGrowth" = psi ~ K / (1 + exp(-r * (surveyLength - x0))),
              "Gompertz" = psi ~ a * exp(-exp(b - c * surveyLength)),
              "UpperHingeThreshold" = NULL  # Formula handled previously in fit_model
            ),
            start_list = switch( # Provide initial parameter values for models
              model_name,
              "ExponentialDecay" = list(initial_offset = quantile(.x$psi, 0.5, na.rm = TRUE), decay_rate = 0.1, asymptote = quantile(.x$psi, 0.9, na.rm = TRUE)),
              "MichaelisMenten" = list(Vmax = max(.x$psi, na.rm = TRUE), Km = mean(.x$surveyLength, na.rm = TRUE), c = 0.1),
              "AsymptoticRegression" = list(a = min(.x$psi, na.rm = TRUE), b = 0.5, c = 0.1),
              "LogisticGrowth" = list(K = max(.x$psi, na.rm = TRUE), r = 0.1, x0 = mean(.x$surveyLength, na.rm = TRUE)),
              "Gompertz" = list(a = max(.x$psi, na.rm = TRUE), b = 0.1, c = 0.1),
              "UpperHingeThreshold" = NULL  # Start list not applicable
            ),
            model_name = model_name, dens=densToAnalyze
          )
        }  
        .x
      }) %>%  
      ungroup()
    
    # Summarize results for this subset
    results_summary <- results %>%
      group_by(simulation) %>%
      slice_tail(n = 1) %>%  # Keep last row for each simulation to select last minute (5, 10, 15, 20, 30th minutes)
      ungroup() %>%
      mutate(
        # Calculate biases for daily, seasonal, and instantaneous occupancy
        # Daily Biases
        ExpDec_DailyBias = plogis(ExponentialDecay_predicted_asymptote) - dailyOcc,
        MichaelisMenten_DailyBias = plogis(MichaelisMenten_predicted_asymptote) - dailyOcc,
        AsymptoticRegression_DailyBias = plogis(AsymptoticRegression_predicted_asymptote) - dailyOcc,
        LogisticGrowth_DailyBias = plogis(LogisticGrowth_predicted_asymptote) - dailyOcc,
        Gompertz_DailyBias = plogis(Gompertz_predicted_asymptote) - dailyOcc,
        UpperHingeThreshold_DailyBias = plogis(UpperHingeThreshold_predicted_asymptote) - dailyOcc,
        # Season Biases
        ExpDec_SeasonBias = plogis(ExponentialDecay_predicted_asymptote) - seasonOcc,
        MichaelisMenten_SeasonBias = plogis(MichaelisMenten_predicted_asymptote) - seasonOcc,
        AsymptoticRegression_SeasonBias = plogis(AsymptoticRegression_predicted_asymptote) - seasonOcc,
        LogisticGrowth_SeasonBias = plogis(LogisticGrowth_predicted_asymptote) - seasonOcc,
        Gompertz_SeasonBias = plogis(Gompertz_predicted_asymptote) - seasonOcc,
        UpperHingeThreshold_SeasonBias = plogis(UpperHingeThreshold_predicted_asymptote) - seasonOcc,
        # Instantaneous Biases
        ExpDec_InstBias = plogis(ExponentialDecay_predicted_yIntercept) - instOcc,
        MichaelisMenten_InstBias = plogis(MichaelisMenten_predicted_yIntercept) - instOcc,
        AsymptoticRegression_InstBias = plogis(AsymptoticRegression_predicted_yIntercept) - instOcc,
        LogisticGrowth_InstBias = plogis(LogisticGrowth_predicted_yIntercept) - instOcc,
        Gompertz_InstBias = plogis(Gompertz_predicted_yIntercept) - instOcc,
        UpperHingeThreshold_InstBias = plogis(UpperHingeThreshold_predicted_yIntercept) - instOcc
      )
    
    # Store in results list
    results_summary_list[[paste(group_name, subset_name, sep = "_")]] <- results_summary
    print("Processed") # Used to track progress
  }  
}  

# Record the end time and calculate elapsed time
endTime = Sys.time()
elapsedTime = endTime - startTime
print(elapsedTime)


# Combine all data frames in results_summary_list into a single data frame
combined_results_summary <- bind_rows(results_summary_list)  


# Save the combined data frame as a CSV
write.csv(combined_results_summary, file = paste0('Full_Results_Summary', densToAnalyze, '.csv'), row.names = FALSE)






################################################################################
#Lucas, is anything below here necessary?  I chunked it all out for now.
################################################################################



# full_results.05.3 = read.csv("Full_Results_Summary.05.3.csv")
# 
# # Rename columns directly
# colnames(full_results.05.3)[colnames(full_results.05.3) == "ExpDec_DailyBias"] <- "ExponentialDecay_DailyBias"
# colnames(full_results.05.3)[colnames(full_results.05.3) == "ExpDec_SeasonBias"] <- "ExponentialDecay_SeasonBias"
# colnames(full_results.05.3)[colnames(full_results.05.3) == "ExpDec_InstBias"] <- "ExponentialDecay_InstBias"
# colnames(full_results.05.3)[colnames(full_results.05.3) == "Initial_InstantBias"] <- "Initial_InstBias"
# 
# # Trim
# full_results.05.3_trimmed <- full_results.05.3 %>%
#   select(simulation, surveyLength, intervalLength, nSurveys, radius, density, 
#          Initial_InstBias, Initial_DailyBias, Initial_SeasonBias, 
#          ExponentialDecay_success, ExponentialDecay_DailyBias, ExponentialDecay_SeasonBias, ExponentialDecay_InstBias,
#          Gompertz_success, Gompertz_DailyBias, Gompertz_SeasonBias, Gompertz_InstBias,
#          MichaelisMenten_success, MichaelisMenten_DailyBias, MichaelisMenten_SeasonBias, MichaelisMenten_InstBias,
#          LogisticGrowth_success, LogisticGrowth_DailyBias, LogisticGrowth_SeasonBias, LogisticGrowth_InstBias,
#          AsymptoticRegression_success, AsymptoticRegression_DailyBias, AsymptoticRegression_SeasonBias, AsymptoticRegression_InstBias,
#          UpperHingeThreshold_success, UpperHingeThreshold_threshold, UpperHingeThreshold_DailyBias, UpperHingeThreshold_SeasonBias, UpperHingeThreshold_InstBias)
# 
# # Pivot longer to make scenario a row
# full_results.05.3_long <- full_results.05.3_trimmed %>%
#   pivot_longer(
#     cols = c(
#       Initial_InstBias, Initial_DailyBias, Initial_SeasonBias,
#       ExponentialDecay_success, ExponentialDecay_DailyBias, ExponentialDecay_SeasonBias, ExponentialDecay_InstBias,
#       Gompertz_success, Gompertz_DailyBias, Gompertz_SeasonBias, Gompertz_InstBias,
#       MichaelisMenten_success, MichaelisMenten_DailyBias, MichaelisMenten_SeasonBias, MichaelisMenten_InstBias,
#       LogisticGrowth_success, LogisticGrowth_DailyBias, LogisticGrowth_SeasonBias, LogisticGrowth_InstBias,
#       AsymptoticRegression_success, AsymptoticRegression_DailyBias, AsymptoticRegression_SeasonBias, AsymptoticRegression_InstBias,
#       UpperHingeThreshold_success, UpperHingeThreshold_threshold, UpperHingeThreshold_DailyBias, UpperHingeThreshold_SeasonBias, UpperHingeThreshold_InstBias
#     ),
#     names_to = c("Model", ".value"),  # Split names into "Model" and corresponding value types
#     names_sep = "_"  # Separator between model name and variable type
#   )
# 
# 
# # Helper function to calculate mean squared error (MSE)
# calc_mse <- function(values) {
#   values <- na.omit(values)  # Remove NA values to avoid computation errors
#   mean_val <- mean(values)  # Compute the mean of the values
#   variance_val <- var(values)  # Compute the variance of the values
#   mse <- (mean_val^2) + variance_val  # Calculate MSE as squared mean plus variance
#   mse  # Return MSE
# }
# 
# # Function to calculate variance confidence limits
# calc_var_cl <- function(values) {
#   values <- na.omit(values)  # Remove NAs
#   n <- length(values)  # Sample size
#   if (n < 2) return(c(NA, NA))  # Can't calculate variance with <2 samples
#   
#   var_val <- var(values)  # Sample variance
#   lower_cl <- (n - 1) * var_val / qchisq(0.975, df = n - 1)
#   upper_cl <- (n - 1) * var_val / qchisq(0.025, df = n - 1)
#   return(c(lower_cl, upper_cl))
# }
# 
# # Group by scenario and calculate summary statistics
# summary_results <- full_results.05.3_long %>%
#   group_by(nSurveys, intervalLength, surveyLength, radius, Model) %>%
#   summarise(
#     # Bias Means
#     mean_InstBias = mean(InstBias, na.rm = TRUE),
#     mean_DailyBias = mean(DailyBias, na.rm = TRUE),
#     mean_SeasonBias = mean(SeasonBias, na.rm = TRUE),
#     
#     # Bias Variances
#     InstBias_Variance = var(InstBias, na.rm = TRUE),
#     DailyBias_Variance = var(DailyBias, na.rm = TRUE),
#     SeasonBias_Variance = var(SeasonBias, na.rm = TRUE),
#     
#     # Variance Confidence Limits
#     InstBias_Var_LowerCL = calc_var_cl(InstBias)[1],
#     InstBias_Var_UpperCL = calc_var_cl(InstBias)[2],
#     DailyBias_Var_LowerCL = calc_var_cl(DailyBias)[1],
#     DailyBias_Var_UpperCL = calc_var_cl(DailyBias)[2],
#     SeasonBias_Var_LowerCL = calc_var_cl(SeasonBias)[1],
#     SeasonBias_Var_UpperCL = calc_var_cl(SeasonBias)[2],
#     
#     # Bias Confidence Limits
#     InstBias_LowerCL = quantile(InstBias, 0.025, na.rm = TRUE),
#     InstBias_UpperCL = quantile(InstBias, 0.975, na.rm = TRUE),
#     DailyBias_LowerCL = quantile(DailyBias, 0.025, na.rm = TRUE),
#     DailyBias_UpperCL = quantile(DailyBias, 0.975, na.rm = TRUE),
#     SeasonBias_LowerCL = quantile(SeasonBias, 0.025, na.rm = TRUE),
#     SeasonBias_UpperCL = quantile(SeasonBias, 0.975, na.rm = TRUE),
#     
#     # Bias MSEs
#     MSE_InstBias = calc_mse(InstBias),
#     MSE_DailyBias = calc_mse(DailyBias),
#     MSE_SeasonBias = calc_mse(SeasonBias),
#     
#     # MSE Confidence Limits
#     MSE_InstBias_LowerCL = mean(InstBias, na.rm = TRUE)^2 + calc_var_cl(InstBias)[1],
#     MSE_InstBias_UpperCL = mean(InstBias, na.rm = TRUE)^2 + calc_var_cl(InstBias)[2],
#     MSE_DailyBias_LowerCL = mean(DailyBias, na.rm = TRUE)^2 + calc_var_cl(DailyBias)[1],
#     MSE_DailyBias_UpperCL = mean(DailyBias, na.rm = TRUE)^2 + calc_var_cl(DailyBias)[2],
#     MSE_SeasonBias_LowerCL = mean(SeasonBias, na.rm = TRUE)^2 + calc_var_cl(SeasonBias)[1],
#     MSE_SeasonBias_UpperCL = mean(SeasonBias, na.rm = TRUE)^2 + calc_var_cl(SeasonBias)[2],
#     
#     # Threshold statistics (only for UpperHingeThreshold)
#     mean_threshold = ifelse(Model == "UpperHingeThreshold", mean(threshold, na.rm = TRUE), NA),
#     threshold_LowerCL = ifelse(Model == "UpperHingeThreshold", quantile(threshold, 0.025, na.rm = TRUE), NA),
#     threshold_UpperCL = ifelse(Model == "UpperHingeThreshold", quantile(threshold, 0.975, na.rm = TRUE), NA),
#     
#     # Success Rate
#     success_rate = mean(success, na.rm = TRUE)
#   ) %>%
#   ungroup()
# 
# # Remove duplicate rows (where every value is identical)
# summary_results <- summary_results %>%
#   distinct(.keep_all = TRUE)
# 
# # Add density back
# summary_results$density = 0.05
# 
# # change relevant values to factors
# summary_results$intervalLength = as.factor(summary_results$intervalLength)
# summary_results$nSurveys = as.factor(summary_results$nSurveys)
# summary_results$radius = as.factor(summary_results$radius)
# 
# # save
# write.csv(summary_results, file = "bias_results.05.3.csv", row.names = FALSE) 
