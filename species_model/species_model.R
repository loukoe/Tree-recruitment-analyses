# Libraries ---------------------------------------------------------------

library(dplyr)
library(nnet)
library(caret)
library(boot)
library(rsample)
library(pscl)
library(purrr)
library(margins)

# Data preparations ---------------------------------------------------------------

species_data <- read.csv("species_model_data_subset.csv")

# Transformations and standardization
# Function to standardize continuous variables
VarStd <- function(x) { (x - mean(x)) / sd(x) }

# Standardizing and transforming variables for NFI method and forest structure
species_data <- species_data %>%
  mutate(
    s_log_n_tree_alive = VarStd(log(n_tree_alive + 1)),
    s_root_ba_sum_alive = VarStd(sqrt(ba_sum_alive)),
    s_ba_skew = VarStd(ba_skew),
    s_log_ba_sum_dead = VarStd(log(ba_sum_dead + 1))
  )

# List of variables to standardize for the five different models
vars_model_1 <- c("MaT", "TaP", "ThHUi", "var33")
vars_model_2 <- c("MweqT", "MwaqT", "MdrqP", "var32")
vars_model_3 <- c("MaDR", "var23", "var37", "SDmT", "SDmP", "SDmR")
vars_model_4 <- c("MAXmT", "MAXmP", "MINmP", "var30")
vars_soil_deposition <- c("RedN", "CEC", "SLTPPT")

# Standardizing variables in model_data efficiently
species_data <- species_data %>%
  mutate(
    across(all_of(vars_model_1), VarStd, .names = "s_{.col}"),
    across(all_of(vars_model_2), VarStd, .names = "s_{.col}"),
    across(all_of(vars_model_3), VarStd, .names = "s_{.col}"),
    across(all_of(vars_model_4), VarStd, .names = "s_{.col}"),
    across(all_of(vars_soil_deposition), VarStd, .names = "s_{.col}")
  )

colnames(species_data)

# Model fitting -----------------------------------------------------------

MN5b <- multinom(ing_sp ~ s_root_ba_sum_alive + s_log_n_tree_alive + s_log_ba_sum_dead +  
                   s_ba_skew  + forest + lead_sp +
                   s_RedN + s_CEC + s_SLTPPT +
                   s_ThHUi + s_var33  + s_MweqT + s_MaDR + 
                   s_SDmT + s_var23 + s_var37 + s_SDmR +
                   s_MAXmT + s_MAXmP + s_MINmP,
                 data=species_data, maxit=10000, trace=T,MaxNWts=10000)
coef(MN5b)

# Model diagnostics -------------------------------------------------------

# Evaluation
# Model AIC and accuracy
AIC(MN5b)
preds5 <- predict(MN5b, type="class", newdata=species_data) # Make model predictions
postResample(species_data$ing_sp,preds5)[1] # Print model accuracy


# Validation
# Set seed for reproducibility
set.seed(123)

# Create 10 folds for cross-validation
cv_folds <- vfold_cv(species_data, v = 10)

# Define a function to train the model and compute classification accuracy
cv_fun_multinom <- function(split) {
  train_data <- analysis(split)  # Training data
  test_data  <- assessment(split)  # Test data
  
  # Fit the multinomial model on training data
  model <- multinom(
    ing_sp ~ s_root_ba_sum_alive + s_log_n_tree_alive + s_log_ba_sum_dead +  
      s_ba_skew  + forest + lead_sp +
      s_RedN + s_CEC + s_SLTPPT + # Soil + deposition +
      s_ThHUi + s_var33  + s_MweqT + s_MaDR + 
      s_SDmT + s_var23 + s_var37 + s_SDmR +
      s_MAXmT + s_MAXmP + s_MINmP,
    data = train_data,
    maxit = 10000, trace = FALSE, MaxNWts = 10000
  )
  
  # Predict on test data
  predictions <- predict(model, newdata = test_data, type = "class")
  
  # Compute classification accuracy
  accuracy <- mean(predictions == test_data$ing_sp)
  
  return(accuracy)
}

# Run Cross-Validation
cv_results <- cv_folds %>%
  mutate(accuracy = map_dbl(splits, cv_fun_multinom))

# Compute Average Accuracy Across Folds
mean(cv_results$accuracy)


# Sensitivity analysis
# Calculate the Average Marginal Effects for the numeric covariates
response_variable <- "ing_sp"
covariates <- c("s_root_ba_sum_alive", "s_log_n_tree_alive","s_log_ba_sum_dead",
                "s_ba_skew",
                "s_RedN","s_CEC" , "s_SLTPPT", 
                "s_ThHUi","s_var33","s_MweqT","s_MaDR", 
                "s_SDmT","s_var23", "s_var37", "s_SDmR", 
                "s_MAXmT", "s_MAXmP", "s_MINmP")

calculate_AME <- function(model, data, covariate, delta = 1e-10) {
  # Predict probabilities for the original data
  original_probs <- predict(model, newdata = data, type = "probs")
  
  # Check if prediction was successful
  if (is.null(original_probs)) {
    warning(paste("Original prediction failed for covariate:", covariate))
    return(NULL)
  }
  
  # Perturb the specified covariate by a small amount (delta)
  data_perturbed <- data
  data_perturbed[[covariate]] <- data_perturbed[[covariate]] + delta
  
  # Predict probabilities for the perturbed data
  perturbed_probs <- predict(model, newdata = data_perturbed, type = "probs")
  
  # Check if perturbed prediction was successful
  if (is.null(perturbed_probs)) {
    warning(paste("Perturbed prediction failed for covariate:", covariate))
    return(NULL)
  }
  
  # Calculate the marginal effect as the change in probability divided by delta
  marginal_effects <- (perturbed_probs - original_probs) / delta
  
  # Average the marginal effects across all observations
  AME <- colMeans(marginal_effects)
  
  return(AME)
}

# Apply the function to the numeric covariates
AMEs_list <- lapply(covariates, function(cov) calculate_AME(MN5b, species_data, cov))

# Combine the results into a data frame for easy viewing
AMEs_df <- do.call(rbind, AMEs_list)
rownames(AMEs_df) <- covariates[!sapply(AMEs_list, is.null)]
colnames(AMEs_df) <- levels(species_data[[response_variable]]) # Replace with actual response variable name

print(AMEs_df) # Average marginal effects per covariate and species_class
rowSums(AMEs_df) # Average marginal effects per covariate
