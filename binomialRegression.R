rm(list=ls())


# Load necessary libraries
library(dplyr)
library(tidyverse)
library(lme4)

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


#------------------------------------------------------------------------------
# Binomial Regression Model
#------------------------------------------------------------------------------
combined_data <- full_results_long %>%
  filter(Model %in% c("Initial", "AsymptoticRegression")) # "Initial" = traditional SOMs

# Convert to binary outcome (1 if Asymptotic Regression is less biased, 0 otherwise)
binomial_data <- combined_data %>%
  select(density, radius, intervalLength, nSurveys, surveyLength, Model, DailyBias, simulation) %>%  
  pivot_wider(names_from = Model, values_from = DailyBias) %>%  
  mutate(closer_to_zero = ifelse(abs(AsymptoticRegression) < abs(Initial), 1, 0))  # Binary outcome

# Convert categorical variables to factors and set reference levels
binomial_data <- binomial_data %>%
  mutate(
    intervalLength = relevel(factor(intervalLength), ref = "none"), 
    nSurveys = relevel(factor(nSurveys), ref = "2"),
    surveyLength = relevel(factor(surveyLength), ref = "5")
  )

binomial_model <- glmer(
  closer_to_zero ~ density + radius + intervalLength + nSurveys + surveyLength + (1 | simulation), 
  data = binomial_data, 
  family = binomial, 
  control = glmerControl(optimizer = "bobyqa")
)

# View the model summary
summary(binomial_model)

