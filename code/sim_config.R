# Root folder of project
setwd("/panfs/jay/groups/23/murra484/pai00032/dnam_cluster/dnam_prediction")

# Template ---------------------------------------------------------------------
sim_ids <- data.frame(matrix(nrow = 0, ncol = 10))
colnames(sim_ids) <- c("sim_id", "k","num_cna", "n", "j", "bias", 
                       "sigma_a_diag", "sigma_e_diag", "est_clusters", 
                       "z_factor")
#saveRDS(sim_id_template, "simulation_results/sim_ids.RDS")

base_sim_id <- data.frame(list(NA, 2275, 65, 655, 6, .78, .2, .9, 6, 1))
#base_sim_id <- data.frame(list(NA, 5, 65, 655, 6, .78, .2, .9, 6, 1))
colnames(base_sim_id) <- c("sim_id", "k","num_cna", "n", "j", "bias", 
                           "sigma_a_diag", "sigma_e_diag", "est_clusters", 
                           "z_factor")

# Add variations of interest ("var"= variation)
sigma_a_diag_var <- c(.1, .5, 1, 2, 3)
j_var <- c(1, 2, 4, 8, 10, 16)
est_clusters_var <- c(1, 2, 4, 8, 10, 16, 20)
bias_var <- c(.5, .6, .681, .9)

# first 1 is for the base
total_vars <- 1 + length(sigma_a_diag_var) + length(j_var) + length(est_clusters_var) + length(bias_var)

sim_id_var <- base_sim_id[rep(1, total_vars),]
sim_id_var$set <- "base"
curr <- 2
sim_id_var[curr:(curr + length(sigma_a_diag_var) - 1), c("sigma_a_diag", "set")] <- data.frame(sigma_a_diag_var, "sigma_a_diag")
curr <- curr + length(sigma_a_diag_var)
sim_id_var[curr:(curr + length(j_var) - 1), c("j", "set")] <- data.frame(j_var, "j")
curr <- curr + length(j_var)
sim_id_var[curr:(curr + length(est_clusters_var) - 1), c("est_clusters", "set")] <- data.frame(est_clusters_var, "est_clusters")
curr <- curr + length(est_clusters_var)
sim_id_var[curr:(curr + length(bias_var) - 1), c("bias", "set")] <- data.frame(bias_var, "bias")

num_replicates <- 300
sim_id_df <- sim_id_var[rep(1:total_vars, each = num_replicates), ]

## Add to sim_ids.RDS
# Assign the actual id column
sim_id_df$sim_id <- 1:nrow(sim_id_df)
# Bind to sim_ids and save
sim_ids <- bind_rows(sim_ids, sim_id_df)
saveRDS(sim_ids, "simulation_results/sim_ids.RDS")

# Create sim_config.RDS
sim_config <- data.frame(sim_id = sim_id_df$sim_id)
# Simple task so give 10(0) simulations per job id # 100 for real, 10 for testing
sim_config$task_id <- rep(1:as.integer(nrow(sim_config) / 2), length = nrow(sim_config))
saveRDS(sim_config, "sim_configs/full_sim_config.RDS")


# Simulation 2: Run everything that hasn't been run yet ------------------------

# Get all sim_ids in compiled results

compiled_ids <- readRDS("simulation_results/full_simulation_results.RDS") %>%
  select(sim_id) %>%
  unique() %>% unlist()

# Get all sim_ids that haven't been compiled
rds_files <- list.files(path = "sim_results_id", pattern = "\\.RDS$", full.names = F)
indiv_ids <- gsub("sim_id_(\\d+)\\.RDS", "\\1", rds_files) %>% as.numeric()

# Filter out sim_ids that haven't been run yet and reassign task_ids
remaining_sim_config <- readRDS("simulation_results/sim_ids.RDS") %>%
  filter(!(sim_id %in% indiv_ids)) %>%
  #filter(!(sim_id %in% c(compiled_ids, indiv_ids))) %>%
  select(sim_id) %>% # 5 simulations per task_id
  mutate(task_id = rep(1:as.integer(nrow(.) / 1), length = nrow(.)))

saveRDS(remaining_sim_config, "sim_configs/remaining_sim_config.RDS")

# Test using subset of sim_ids -------------------------------------------------

sim_ids <- readRDS("simulation_results/sim_ids.RDS")

sim_ids_try <- c(1:100, 301:400, 601:700, 901:1000)
sim_config_test <- data.frame(sim_id = sim_ids_try, task_id = 1:length(sim_ids_try))
saveRDS(sim_config_test, "sim_configs/sim_config_test.RDS")




