# Libraries ---------------------------------------------------------------
library(dplyr)
library(glmmTMB)
library(DHARMa)
library(boot)

# Data preparations ---------------------------------------------------------------

count_data <- read.csv("count_model_data_subset.csv")

# Transformations and standardization

# Function to standardize continuous variables
VarStd <- function(x) { (x - mean(x)) / sd(x) }

# Standardizing variables for NFI method and forest structure
count_data <- count_data %>%
  mutate(
    s_l_dbh = VarStd(log(dbh)),
    s_log_n_tree_alive = VarStd(log(n_tree_alive + 1)),
    s_root_ba_sum_alive = VarStd(sqrt(ba_sum_alive)),
    s_ba_skew = VarStd(ba_skew),
    s_log_ba_sum_dead = VarStd(log(ba_sum_dead + 1))
  )

# List of variables to standardize for models
vars_model_1 <- c("MaT", "TaP", "ThHUi", "var33")
vars_model_2 <- c("MweqT", "MwaqT", "MdrqP", "var32")
vars_model_3 <- c("MaDR", "var23", "var37", "SDmT", "SDmP", "SDmR")
vars_model_4 <- c("MAXmT", "MAXmP", "MINmP", "var30")
vars_soil_deposition <- c("RedN", "CEC", "SLTPPT")

# Standardizing variables in model_data efficiently
count_data <- count_data %>%
  mutate(
    across(all_of(vars_model_1), VarStd, .names = "s_{.col}"),
    across(all_of(vars_model_2), VarStd, .names = "s_{.col}"),
    across(all_of(vars_model_3), VarStd, .names = "s_{.col}"),
    across(all_of(vars_model_4), VarStd, .names = "s_{.col}"),
    across(all_of(vars_soil_deposition), VarStd, .names = "s_{.col}")
  )

# Model fitting of model 5 -----------------------------------------------------------

# Negative binomial model
mNB1 <- glmmTMB(number.in ~ s_l_dbh  +  
                  s_root_ba_sum_alive * s_log_n_tree_alive + s_log_ba_sum_dead + 
                  forest + #s_ba_skew + 
                  s_RedN + s_CEC + s_SLTPPT + s_ThHUi  + s_MaDR +  s_MweqT +
                  s_MAXmT + s_MAXmP + s_var23 + s_var37  + # + s_var30
                  offset(log(ing.plot.area))+ offset(log(interval)),
                family="nbinom2",data=count_data,verbose=TRUE,
                control = glmmTMBControl(profile=quote(length(parameters$beta)>=5),parallel = 12,
                                         optimizer = optim, optArgs = list(method="BFGS")))
gc()
summary(mNB1)

# Model diagnostics -------------------------------------------------------

# Evaluation
# Model AIC
AIC(mNB1)

# Compare number of observed plots without recruitment to predictions
sum(count_data$number.in==0) # Observed

P <- simulate(mNB1, 1000) # Simulate 1000 datasets from the model
mean(colSums(P==0)) # Count number of 0s in the simulated datasets


# Validation
# Function to compute deviance residuals via DHARMa
cv_fun_dharma <- function(data, indices) {
  # Bootstrap sample
  train_data <- data[indices, ]
  
  # Fit the glmmTMB model
  model <- glmmTMB(
    number.in ~ s_l_dbh + 
      s_root_ba_sum_alive * s_log_n_tree_alive + 
      s_log_ba_sum_dead + forest + 
      s_RedN + s_CEC + s_SLTPPT + s_ThHUi + 
      s_MaDR + s_MweqT + s_MAXmT + 
      s_MAXmP + s_var23 + s_var37 + 
      offset(log(ing.plot.area)) + 
      offset(log(interval)),
    family = nbinom2, 
    data = train_data,
    control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS"))
  )
  
  # Simulate residuals using DHARMa
  sim_res <- simulateResiduals(fittedModel = model, plot = FALSE)
  
  # Calculate mean deviance residuals
  deviance_res <- mean(residuals(sim_res, type = "deviance"))
  
  return(deviance_res)
}

# Perform Bootstrapped Cross-Validation with 100 resamples
set.seed(123)
cv_results_boot2 <- boot(
  data = count_data, 
  statistic = cv_fun_dharma, 
  R = 100  # Number of bootstrap samples
)

# Print bootstrap results
print(cv_results_boot2)
