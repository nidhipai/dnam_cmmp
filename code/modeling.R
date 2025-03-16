
# Run rp and k_cmmp methods on the data

# Root folder of project
# set wd here!

source("code/preliminaries.R")
source("code/cmmp_functions.R")
source("code/sim_functions.R")


data <- readRDS("data/clustered_6.RDS")


meth_var <- grep("^meth_", colnames(data), value = TRUE) # meth var names
cna_var <- grep("^cna_", colnames(data), value = TRUE)

Y_names <- meth_var
X_names <- c(cna_var, "RACE", "SEX", "AGE", "STAGE", "cancer_id")
X_string <- paste0(X_names, collapse = " + ")

# Masks that capture which samples are minorities
test_minority_mask <- (!data$train) & (data$race != "White")
mask_2 <- data[!data$train, ]$race != "White" # Just on test
num_test_min <- sum(test_minority_mask)

# 1: Regression prediction -----------------------------------------------------

reg_pred_pred <- reg_pred(data)
reg_pred_mse <- mse_all(data[!data$train, Y_names], reg_pred_pred)
reg_pred_mse_minority <- mse_all(data[test_minority_mask, Y_names], 
                                 reg_pred_pred[mask_2, , drop = F])
reg_pred_res <- rbind(data.frame(metric = "mse", outcome = Y_names, value = reg_pred_mse),
                      data.frame(metric = "mse_minority", outcome = Y_names, value = reg_pred_mse_minority))
reg_pred_res$method <- "reg_pred"

# 2: k-cmmp --------------------------------------------------------------------

num_clusters <- 6

# Use CMMP to get predictions on test data
cmmp_res <- cmmp_core(data)

# Return predictions and list of clusters for train, test, plus metrics
k_cmmp_res <- list(cmmp_res[[2]], data[data$train, "cluster"], cmmp_res[[1]],
     cmmp_res[[4]], cmmp_res[[5]], cmmp_res[[6]])

k_cmmp_pred <- k_cmmp_res[[1]]
k_cmmp_mse <- mse_all(data[!data$train, Y_names], k_cmmp_pred)
k_cmmp_mse_minority <- mse_all(data[test_minority_mask, Y_names], 
                               k_cmmp_pred[mask_2, , drop = F])
k_cmmp_res <- rbind(data.frame(metric = "mse", outcome = Y_names, value = k_cmmp_mse),
                    data.frame(metric = "mse_minority", outcome = Y_names, value = k_cmmp_mse_minority),
                    data.frame(metric = "re_var", outcome = Y_names, value = k_cmmp_res[[4]]), 
                    data.frame(metric = "sigma_sq", outcome = Y_names, value = k_cmmp_res[[5]]),
                    data.frame(metric = "snr", outcome = Y_names, value = k_cmmp_res[[6]]))
k_cmmp_res$method <- "k_cmmp"

# Save details for this model as well 
saveRDS(cmmp_res, "result_matrices/real_data_kcmmp_covar_cmmp_res_2.RDS")

result_df <- bind_rows(reg_pred_res, k_cmmp_res)
saveRDS(result_df, "result_matrices/real_data_kcmmp_covar_2.RDS")




















