# Prep ---------------------------------------------------------------

# Root folder of project
# set wd!

# Source files and saved data
source("code/preliminaries.R")
source("code/sim_functions.R")
source("code/cmmp_functions.R")

data <- readRDS("data/bound_data.RDS")
cna_var <- grep("^cna_", colnames(data), value = TRUE)

sim_ids <- readRDS("simulation_results/sim_ids.RDS")

# Read in sim_config and task_id
command_args <- commandArgs(trailingOnly = T)
sim_config <- readRDS(paste0("sim_configs/", command_args[1])) 
curr_task_id <- command_args[2]

# Filter sim_ids for this task_id
task_sim_config <- sim_config %>% filter(task_id == curr_task_id)
sim_id_list <- task_sim_config$sim_id

# Likely fixed for all simulations
train_prop <- .7

# Iterate over sims matching task_id
for (curr_sim_id in sim_id_list) {
  print(paste0("curr_sim_id/total in task: ", match(curr_sim_id, sim_id_list), 
               " / ", length(sim_id_list)))
  
  # Set up sim -----------------------------------------------------------------
  sim_id_row <- sim_ids %>% filter(sim_id == curr_sim_id)
  
  tryCatch({
    # Helper lists
    Y_names <- paste0("Y", 1:sim_id_row$k)
    cna_names <- paste0("cna_", 1:sim_id_row$num_cna)
    theta_names <- paste0("theta_", 1:sim_id_row$k)
    X_names <- c("race", "cancer", "sex", "age", "stage", cna_names)
    X_string <- paste0(X_names, collapse = " + ")
    
    # Set seed to generate sample
    set.seed(sim_id_row$sim_id)
    # Generate sample
    sample_data <- generate_sample(sim_id_row$k, sim_id_row$num_cna, sim_id_row$n, 
                                   sim_id_row$j, sim_id_row$sigma_a_diag, 
                                   sim_id_row$sigma_e_diag, sim_id_row$z_factor,
                                   sim_id_row$bias)
    
    # Run methods ----------------------------------------------------------------
    
    # Split into train and test
    # Stratified so training set includes all races
    train_index <- suppressWarnings(caret::createDataPartition(sample_data$race, p = train_prop, list = F))
    sample_data$train <- 1:nrow(sample_data) %in% train_index

    
    # Fit each model and get mspe for each outcome
    reg_pred_pred <- reg_pred(sample_data)
    reg_pred_mse <- mse_all(sample_data[!sample_data$train, theta_names], reg_pred_pred)
    reg_pred_res <- data.frame(metric = "mse", outcome = Y_names, value = reg_pred_mse, method = "reg_pred")
    
    k_cmmp_res <- k_cmmp(sample_data, sim_id_row$est_clusters)
    k_cmmp_pred <- k_cmmp_res[[1]]
    ari_train <- mclust::adjustedRandIndex(k_cmmp_res[[2]],
                                           sample_data[sample_data$train, "group"])
    ari_test <- mclust::adjustedRandIndex(cmmp_clust_to_cluster(k_cmmp_res[[3]]),
                                          sample_data[!sample_data$train, "group"])
    k_cmmp_mse <- mse_all(sample_data[!sample_data$train, theta_names], k_cmmp_pred)
    k_cmmp_res <- rbind(data.frame(metric = "mse", outcome = Y_names, value = k_cmmp_mse),
                        data.frame(metric = "ari_train", outcome = NA, value = ari_train),
                        data.frame(metric = "ari_test", outcome = NA, value = ari_test)
                        #data.frame(metric = "re_var", outcome = Y_names, value = k_cmmp_res[[4]]), 
                        #data.frame(metric = "sigma_sq", outcome = Y_names, value = k_cmmp_res[[5]]),
                        #data.frame(metric = "snr", outcome = Y_names, value = k_cmmp_res[[6]])
                        )
    k_cmmp_res$method <- "k_cmmp"
    
    cmmp_res <- group_cmmp(sample_data) # CMMP without cluster, just use group
    cmmp_pred <- cmmp_res[[1]]
    ari_train <- mclust::adjustedRandIndex(cmmp_res[[2]],
                                           sample_data[sample_data$train, "group"])
    ari_test <- mclust::adjustedRandIndex(cmmp_clust_to_cluster(cmmp_res[[3]]),
                                          sample_data[!sample_data$train, "group"])
    cmmp_mse <- mse_all(sample_data[!sample_data$train, theta_names], cmmp_pred)
    cmmp_res <- rbind(data.frame(metric = "mse", outcome = Y_names, value = cmmp_mse),
                      data.frame(metric = "ari_train", outcome = NA, value = ari_train),
                      data.frame(metric = "ari_test", outcome = NA, value = ari_test)
                     # data.frame(metric = "re_var", outcome = Y_names, value = cmmp_res[[4]]), 
                      #data.frame(metric = "sigma_sq", outcome = Y_names, value = cmmp_res[[5]]),
                      #data.frame(metric = "snr", outcome = Y_names, value = cmmp_res[[6]])
                     )
    cmmp_res$method <- "cmmp"
  
    
    bar_y_pred <- bar_y(sample_data)
    bar_y_mse <- mse_all(sample_data[!sample_data$train, theta_names], bar_y_pred)
    bar_y_res <- data.frame(metric = "mse", outcome = Y_names, value = bar_y_mse, method = "bar_y")
    bar_y_res$method <- "bar_y"
    
    result_df <- bind_rows(reg_pred_res, k_cmmp_res, cmmp_res, bar_y_res)
    
    # If this sim_id requires mse_minority 
    if (sim_id_row$set == "bias" | sim_id_row$set == "base") {
      # Helper masks
      test_minority_mask <- (!sample_data$train) & (sample_data$race != "White")
      mask_2 <- sample_data[!sample_data$train, ]$race != "White" # Just on test
      num_test_min <- sum(test_minority_mask)
      reg_pred_mse_minority <- mse_all(sample_data[test_minority_mask, theta_names], 
                                       reg_pred_pred[mask_2, , drop = F])
      k_cmmp_mse_minority <- mse_all(sample_data[test_minority_mask, theta_names], 
                                     k_cmmp_pred[mask_2, , drop = F])
      cmmp_mse_minority <- mse_all(sample_data[test_minority_mask, theta_names], 
                                   cmmp_pred[mask_2, ])
      bar_y_mse_minority <- mse_all(sample_data[test_minority_mask, theta_names], 
                                    bar_y_pred[mask_2, ])
      # Bind together results
      minority_result_df <- rbind(
        data.frame(metric = "mse_minority", outcome = Y_names, value = reg_pred_mse_minority, method = "reg_pred"),
        data.frame(metric = "mse_minority", outcome = Y_names, value = k_cmmp_mse_minority, method = "k_cmmp"),
        data.frame(metric = "mse_minority", outcome = Y_names, value = cmmp_mse_minority, method = "cmmp"),
        data.frame(metric = "mse_minority", outcome = Y_names, value = bar_y_mse_minority, method = "bar_y"))
      
      result_df <- bind_rows(result_df, minority_result_df)
    }
  
    result_df$sim_id <- curr_sim_id

    # Save results
    saveRDS(result_df, paste0("sim_results_id/sim_id_", curr_sim_id, ".RDS"))
    
  }, error = function(e) {
    print(e)
  })
}
