# Get bootstrap standard error of estimates in real data

# Root folder of project

source("code/preliminaries.R")
source("code/cmmp_functions.R")
source("code/sim_functions.R")


data <- readRDS("data/clustered_6.RDS")

meth_var <- grep("^meth_", colnames(data), value = TRUE) # meth var names
cna_var <- grep("^cna_", colnames(data), value = TRUE)

Y_names <- meth_var
X_names <- c(cna_var, "RACE", "SEX", "AGE", "STAGE", "cancer_id")
X_string <- paste0(X_names, collapse = " + ")

# CMMP results
cmmp_res <- readRDS("result_matrices/real_data_cmmp_res_pred.RDS")
cluster_predictions <- cmmp_res[[1]]
re_var_list <- cmmp_res[[4]]
error_var_list <- cmmp_res[[5]]
predictions <- cmmp_res[[10]]

# Get the job number and seeds
command_args <- commandArgs(trailingOnly = T) 
curr_task_id <- as.integer(command_args[1])
print(curr_task_id)

# Just one bootstrap iteration in this file

# Create bootstrapped data------------------------------------------------------

set.seed(curr_task_id)

# Assign id to make it easy to track patients in BS samples
data <- data %>%
  mutate(id = 1:nrow(data))

## Generate boostrap sample
# Create dataframe to sample from
# Get fixed effect predictions for both the training and test data
data_joined <- cbind(data[, c(X_names, "cluster", "train", "id")], predictions) # Only data we need

# Bootstrap sample (stratified on train/test)
bootstrap_data <- data_joined %>%
  group_by(train) %>%
  slice_sample(prop = 1, replace = T) %>%
  ungroup()

theta <- data[, Y_names] # Just to get the structure of the data

for (outcome_var in Y_names) {
  # Attach clusters for this outcome variable for the test set
  test_cluster_df <- data.frame(id = data[!data$train, "id"], 
                                test_cluster = cluster_predictions[, outcome_var]) # Some are NA
  bootstrap_data <- bootstrap_data %>%
    left_join(test_cluster_df, by = "id") %>%
    mutate(cluster = ifelse(is.na(cluster), test_cluster, cluster))
  
  # Generate random effect
  bs_rand_ef <- rnorm(6, mean = 0, sd = sqrt(re_var_list[outcome_var]))
  
  # Attach random effects
  # This will give a warning because Population Mean is a level in the clusters but not in bs_rand_ef
  bootstrap_data <- bootstrap_data %>% 
    left_join(data.frame(cluster = as.factor(1:6), rand_ef = bs_rand_ef), by = "cluster") %>%
    mutate(rand_ef = replace_na(rand_ef, 0))
  
  # Generate new data
  error_var <- error_var_list[outcome_var]
  error <- rnorm(nrow(bootstrap_data), mean = 0, sd = sqrt(error_var))
  theta[, outcome_var] <- bootstrap_data[, outcome_var] + bootstrap_data$rand_ef 
  bootstrap_data[, outcome_var] <- theta[, outcome_var] + error
  
  # Remove cluster for test data (and other temporary cols)
  bootstrap_data[!bootstrap_data$train, "cluster"] <- NA
  bootstrap_data$test_cluster <- NULL
  bootstrap_data$rand_ef <- NULL
}

bs_data <- as.data.frame(bootstrap_data)

# Run methods ------------------------------------------------------------------

# Masks that capture which samples are minorities
test_minority_mask <- (!bs_data$train) & (bs_data$race != "White")
mask_2 <- bs_data[!bs_data$train, ]$race != "White" # Just on test
num_test_min <- sum(test_minority_mask)


## 1: Regression prediction -----------------------------------------------------

reg_pred_pred <- reg_pred(bs_data)
reg_pred_mse <- mse_all(theta[!bs_data$train,], reg_pred_pred)
reg_pred_mse_minority <- mse_all(theta[test_minority_mask,], 
                                 reg_pred_pred[mask_2, , drop = F])
reg_pred_res <- rbind(data.frame(metric = "mse", outcome = Y_names, value = reg_pred_mse),
                      data.frame(metric = "mse_minority", outcome = Y_names, value = reg_pred_mse_minority))
reg_pred_res$method <- "reg_pred"

## 2: k-cmmp --------------------------------------------------------------------

# Use CMMP to get predictions on test data
cmmp_res <- cmmp_core(bs_data, predict = T)

k_cmmp_pred <- cmmp_res[[2]]
k_cmmp_mse <- mse_all(theta[!bs_data$train, ], k_cmmp_pred)
k_cmmp_mse_minority <- mse_all(data[test_minority_mask, Y_names], 
                               k_cmmp_pred[mask_2, , drop = F])
k_cmmp_res <- rbind(data.frame(metric = "mse", outcome = Y_names, value = k_cmmp_mse),
                    data.frame(metric = "mse_minority", outcome = Y_names, value = k_cmmp_mse_minority),
                    data.frame(metric = "re_var", outcome = Y_names, value = cmmp_res[[4]]), 
                    data.frame(metric = "sigma_sq", outcome = Y_names, value = cmmp_res[[5]]),
                    data.frame(metric = "snr", outcome = Y_names, value = cmmp_res[[6]]))
k_cmmp_res$method <- "k_cmmp"

## 3: Naive --------------------------------------------------------------------

bar_y_pred <- bar_y(bs_data)
bar_y_mse <- mse_all(theta[!bs_data$train,], bar_y_pred)
bar_y_mse_minority <- mse_all(theta[test_minority_mask,], 
                              bar_y_pred[mask_2, , drop = F])
bar_y_res <- rbind(data.frame(metric = "mse", outcome = Y_names, value = bar_y_mse),
                      data.frame(metric = "mse_minority", outcome = Y_names, value = bar_y_mse_minority))
bar_y_res$method <- "bar_y"

# Compile results -------------------------------------------------------------

result_df <- bind_rows(reg_pred_res, k_cmmp_res, bar_y_res)
saveRDS(result_df, paste0("bootstrap_res/bs_seed_", curr_task_id, ".RDS"))
