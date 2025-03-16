library(dplyr)
library(tidyr)
library(ggplot2)

# Root folder of project
# set wd here!

# Compile individual simulation results & move indiv ones to archive -----------

rds_files <- list.files(path = "sim_results_id", pattern = "\\.RDS$", full.names = TRUE)

indiv_dfs <- list() # df from each one
# Loop over each file, read it, and append it to the list
compiled_list <- list()
for (f in 1:length(rds_files)) {
  if (f %% 500 == 0) {
    print(f)
    compiled_temp <- do.call(rbind, indiv_dfs)
    compiled_list[[length(compiled_list) + 1]] <- compiled_temp
    indiv_dfs <- list()
  }
  df <- readRDS(rds_files[f])
  indiv_dfs[[length(indiv_dfs) + 1]] <- df
  
  # Move individual files to archive
  file.copy(rds_files[f], "sim_results_id/archive_dec")
  unlink(rds_files[f])
}
compiled_list[[length(compiled_list) + 1]] <- do.call(rbind, indiv_dfs)

saveRDS(compiled_list, "simulation_results/compiled_list_dec31.RDS")
compiled <- do.call(rbind, compiled_list)

saveRDS(compiled, "simulation_results/full_simulation_results.RDS")


# Load in existing compiled results and append if sim_id not included yet
full_sim_results <- readRDS("simulation_results/full_simulation_results.RDS")
new_compiled <- compiled %>%
  filter(!(sim_id %in% full_sim_results$sim_id))
full_sim_results <- bind_rows(full_sim_results, new_compiled)
#saveRDS(full_sim_results, "simulation_results/full_simulation_results.RDS")


# Analyze simulation set up ----------------------------------------------------

sim_ids <- readRDS("simulation_results/sim_ids.RDS")

# Fix sim_config error
full_sim_results <- full_sim_results %>%
  filter((sim_id < 1501) | (sim_id > 1800))


# Plot 1: sigma_a_diag
sigma_a_diag_ids <- sim_ids %>%
  filter(set == "sigma_a_diag" | set == "base")

sim_res <- full_sim_results %>%
  right_join(sigma_a_diag_ids, by = "sim_id")

sigma_plot <- sim_res %>%
  filter(metric == "mse") %>% 
  mutate(sigma_a_diag = factor(sigma_a_diag),
         method = factor(method)) %>%
  mutate(method = fct_recode(method,"Naive y" = "bar_y", "Oracle CMMP" = "cmmp", 
                                       "k-means + CMMP" = "k_cmmp", "Regression pred." = "reg_pred")) %>%
  ggplot(aes(y = value, x = sigma_a_diag, fill = method)) +
  geom_boxplot(outliers = F) +
  theme_bw() +
  guides(fill = "none") +
  #facet_wrap(~sigma_a_diag, scales = "free") +
  labs(x = "Std. dev. of random effect", y = "MSE", fill = "Method", title = "(a)")

# Plot 2: True clusters

j_ids <- sim_ids %>%
  filter(set == "j" | set == "base")

sim_res <- full_sim_results %>%
  right_join(j_ids, by = "sim_id")

j_plot <- sim_res %>%
  filter(metric == "mse") %>%
  mutate(j = factor(j),
         method = factor(method)) %>%
  mutate(method = fct_recode(method,"Naive y" = "bar_y", "Oracle CMMP" = "cmmp", 
                             "k-means + CMMP" = "k_cmmp", "Regression pred." = "reg_pred")) %>%
  filter(method != "Naive y") %>%
  ggplot(aes(x = j, y = value, fill = method)) +
  geom_boxplot(outliers = F) +
  theme_bw() +
  guides(fill = "none") +
  scale_fill_manual(values = c("#7CAE00", "#00BFC4", "#C77CFF")) +
  labs(x = "Number of clusters", y = "MSE", fill = "Method", title = "(b)")

# Plot 3: Est clusters
est_clusters_ids <- sim_ids %>%
  filter(set == "est_clusters" | set == "base")

sim_res <- full_sim_results %>%
  right_join(est_clusters_ids, by = "sim_id") %>%
  filter((set == "est_clusters" & method == "k_cmmp") | set == "base") %>%
  mutate(est_clusters = ifelse(method == "k_cmmp", est_clusters, "Base"))

est_clusters_plot <- sim_res %>%
  filter(metric == "mse") %>%
  filter(method != "y_bar") %>% # TODO REMOVE
  mutate(est_clusters = factor(est_clusters),
         method = factor(method)) %>%
  mutate(est_clusters = fct_relevel(est_clusters, c("Base", 1, 2, 4, 6, 8, 10, 16, 20))) %>%
  mutate(method = fct_recode(method,"Naive y" = "bar_y", "Oracle CMMP" = "cmmp", 
                             "k-means + CMMP" = "k_cmmp", "Regression pred." = "reg_pred")) %>%
  ggplot(aes(x = est_clusters, y = value, fill = method)) +
  geom_boxplot(outliers = F) +
  theme_bw() +
  guides(fill = "none") +
  labs(x = "Number of estimated clusters", y = "MSE", fill = "Method", title = "(c)")

# Plot 4: Bias
bias_ids <- sim_ids %>%
  filter(set == "bias" | set == "base")

sim_res <- full_sim_results %>%
  right_join(bias_ids, by = "sim_id")

bias_plot <- sim_res %>%
  filter(metric == "mse" | metric == "mse_minority") %>%
  mutate(bias = factor(bias),
         method = factor(method)) %>%
  mutate(method = fct_recode(method,"Naive y" = "bar_y", "Oracle CMMP" = "cmmp", 
                             "k-means + CMMP" = "k_cmmp", "Regression pred." = "reg_pred")) %>%
  ggplot(aes(x = bias, y = value, fill = method)) +
  geom_boxplot(outliers = F) +
  theme_bw() +
  #facet_wrap(~ metric) +
  #guides(fill = "none") +
  labs(x = "Proportion of sample that is White", y = "Minority MSE", 
       fill = "Method", title = "(d)")

# Combine plots

combined_plot <- (sigma_plot + j_plot) / (est_clusters_plot + bias_plot)

ggsave("images/sim_boxplots.png", combined_plot, height = 8, width = 10,
       dpi = 800)


# Level of clustering accuracy

sim_res <- full_sim_results %>%
  right_join(j_ids, by = "sim_id") 

sim_res %>% filter(metric == "ari_train" | metric == "ari_test") %>%
  mutate(j = factor(j),
         method = factor(method)) %>%
  mutate(method = fct_recode(method, "Oracle CMMP" = "cmmp", 
                             "k-means + CMMP" = "k_cmmp",)) %>%
  ggplot(aes(x = j, y = value, fill = method)) +
  geom_boxplot(outliers = F) +
  theme_bw() +
  facet_wrap(~metric) +
  #guides(fill = "none") +
  labs(x = "Number of estimated clusters", y = "MSE", fill = "Method")


# Across sigma_a_diag instead

sim_res <- full_sim_results %>%
  right_join(sigma_a_diag_ids, by = "sim_id") 

ari_plot <- sim_res %>% filter(metric == "ari_train" | metric == "ari_test") %>%
  mutate(sigma_a_diag = factor(sigma_a_diag),
         method = factor(method),
         metric = factor(metric)) %>%
  mutate(method = fct_recode(method, "Oracle CMMP" = "cmmp", 
                             "k-means + CMMP" = "k_cmmp",)) %>%
  mutate(metric = fct_recode(metric, "Test ARI" = "ari_test", 'Train ARI' = "ari_train"),
         metric = fct_relevel(metric, "Train ARI", "Test ARI")) %>%
  ggplot(aes(x = sigma_a_diag, y = value, fill = method)) +
  geom_boxplot(outliers = F) +
  theme_bw() +
  facet_wrap(~metric) +
  #guides(fill = "none") +
  scale_fill_manual(values = c("#7CAE00", "#00BFC4")) +
  theme(strip.background = element_rect(fill = "transparent", color = "transparent")) +
  labs(x = "Std. dev. of random effect", y = "MSE", fill = "Method")


ggsave("images/ari_plot.png", ari_plot, height = 4, width = 6, dpi = 600)


  








