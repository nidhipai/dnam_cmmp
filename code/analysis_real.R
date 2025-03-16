library(dplyr)
library(tidyr)
library(ggplot2)

# Overall results --------------------------------------------------------------

real_df <- readRDS("result_matrices/real_data_kcmmp_covar_2.RDS")
real_data_boxplot <- real_df %>%
  filter(metric == "mse") %>%
  mutate(method = fct_recode(method, "k-means + CMMP" = "k_cmmp", "Regression pred." = "reg_pred")) %>%
  ggplot(aes(y = value, x = method, fill = method)) +
    geom_boxplot() +
    theme_bw() +
    guides(fill = "none") +
    scale_fill_manual(values = c("#00BFC4", "#C77CFF")) +
    labs(x = "Method", y = "MSE")
ggsave("images/real_data_boxplot.png", real_data_boxplot, 
       height = 3, width = 4, dpi = 800)

rp <- real_df %>%
  filter(metric == "mse", method == "reg_pred") %>%
  select(value)
kcmmp <- real_df %>%
  filter(metric == "mse", method == "k_cmmp") %>%
  select(value)

t.test(rp, kcmmp)

# Betas ---------------------------------------------------------------

real_data_res <- readRDS("result_matrices/real_data_kcmmp_covar_cmmp_res_2.RDS")

betas <- data.frame(t(real_data_res[[7]])) %>%
  mutate(outcome = rownames(.))

betas %>%
  select(starts_with("RACE"), outcome) %>%
  select(- RACENHPI) %>%
  pivot_longer(cols = -outcome, names_to = "race", values_to = "beta") %>%
  ggplot(aes(x = outcome, y = beta, color = race)) +
    geom_point()

# Check betas ------------------------------------------------------------------
# Not in paper

cmmp_core_res_full <- readRDS("result_matrices/real_data_kcmmp_covar_cmmp_res_2.RDS")

data.frame(RE = cmmp_core_res_full[[4]], Noise = cmmp_core_res_full[[5]]) %>%
  mutate(oc = Y_names) %>%
  pivot_longer(cols = -oc, names_to = "Component", values_to = "var") %>%
  ggplot(aes(x = var, fill = Component)) +
    geom_histogram(color = "#636363", position = "dodge") +
    labs(x = "Estimated variance", y = "Number of outcomes") +
    theme_bw()

cmmp_core_res_full[[6]]
# Plot histogram of SNR over outcomes
data.frame(snr = cmmp_core_res_full[[6]]) %>%
 ggplot(aes(x = snr)) +
  geom_histogram(color = "#636363", fill = "skyblue", position = "dodge") +
  labs(x = "SNR", y = "Number of outcomes") +
  theme_bw()

# Check variance of X beta for various X and beta standard deviations
design_mat <- model.matrix(as.formula(paste0("~", X_string)), data[data$train, ])

p <- length(cna_var) + 11 # Number of covariates (num_cna + race + cancer + intercept)
beta_test <- matrix(rnorm(p, mean = 0, sd = .1), nrow = p, ncol = 1)
design_mat %*% beta_test %>% var()

sapply(1:2000, function(x) {
  beta_test <- matrix(rnorm(p, mean = 0, sd = .1), nrow = p, ncol = 1)
  design_mat %*% beta_test %>% var()
}) %>%
  data.frame(x = .) %>%
  ggplot(aes(x = x)) +
    geom_histogram(fill = "skyblue", color = "#636363") +
  theme_bw()

# Assignment of random effects in the test set ---------------------------------

data <- readRDS("data/clustered_6.RDS")
meth_var <- grep("^meth_", colnames(data), value = TRUE) # meth var names
cna_var <- grep("^cna_", colnames(data), value = TRUE)

cont_vars <- c("AGE", "FRACTION_GENOME_ALTERED")
factor_vars <- c("RACE", "SEX", "STAGE")

Y_names <- meth_var
X_names <- c(cna_var, "RACE", "SEX", "AGE", "STAGE", "cancer_id")
X_string <- paste0(X_names, collapse = " + ")

# Function to get the most common value, resolving ties randomly
most_common_value <- function(row) {
  freq <- table(as.matrix(row)) # Get frequency of each value
  max_freq <- max(freq) # Find the maximum frequency
  common_values <- names(freq[freq == max_freq]) # Get values with max frequency
  sample(common_values, 1) # Randomly select one in case of tie
}

test_clusters <- data %>%
  select(train, all_of(c(X_names, cont_vars, factor_vars))) %>%
  filter(train == F) %>%
  mutate(RACE = fct_drop(RACE))
test_clusters$mode_assignment <- apply(real_data_res[[1]], 1, most_common_value)

# # Same analysis as clustering

# Cluster vs. cancer type
test_clusters %>%
  slice_sample(n = nrow(test_clusters)) %>% # Shuffle
  mutate(cancer_id = fct_relabel(cancer_id, toupper)) %>%
  mutate(RACE = race_abbrev[RACE]) %>%
  ggplot(aes(x = mode_assignment, y = RACE)) +
  geom_jitter(width = .25, height = .20,  alpha = .4) +
  labs(x = "Mode cluster", y = "Race") +
  theme_bw() +
  scale_y_discrete(limits = rev)

#ggsave("images/clusters_jitter_race_cancer.png", height = 4, width = 6)

# Plot clusters by race
ggplot(test_clusters, aes(x = mode_assignment, fill = RACE)) +
  geom_bar(position = "dodge") +
  theme_bw() +
  ggtitle("Test clusters by race") + 
  labs(x = "Cluster index", y = "Count")

# Test for an association between clinical variables and cluster----

# Cont variables
p_values_cont <- sapply(cont_vars, function(var_name) {
  lm_res <- lm(as.formula(paste(var_name, " ~ mode_assignment")), data = test_clusters)
  return(anova(lm_res)[1, "Pr(>F)"])
})

# Categorical variables
p_values_categorical <- sapply(factor_vars, function(var_name) {
  res <- chisq.test(table(test_clusters[, c(var_name, "mode_assignment")]))
  return(res$p.value)
})

all_p_values <- c(p_values_cont, p_values_categorical)
all_p_values < (.05 / length(all_p_values))

# Plot clinical variable vs. cluster for any interesting ones
# Fraction genome altered, ordered by mean
train_data %>%
  #mutate(cross = paste(cluster, cancer_id)) %>%
  ggplot(aes(x = cluster, y = FRACTION_GENOME_ALTERED)) +
  geom_boxplot(position = "dodge") +
  theme_bw() +
  labs(x = "Cluster", y = "Fraction genome altered")

# Stage
ggplot(train_data, aes(x = cluster, 
                       fill = STAGE)) +
  geom_bar(position = "fill") +
  theme_bw() +
  labs(x = "Cluster", y = "Stage")

# Age - not very interesting
ggplot(train_data, aes(x = cluster, 
                       y = AGE)) +
  geom_boxplot(position = "dodge") +
  theme_bw() + 
  labs(x = "Cluster", y = "Age")

# Magnitude of RE --------------------------------------------------------------

not_cna_var <- c("RACE", "SEX", "AGE", "cancer_id", "STAGE",
                 "CANCER_TYPE_DETAILED", "FRACTION_GENOME_ALTERED", "HISTOLOGICAL_DIAGNOSIS",
                 "TISSUE_SOURCE_SITE")

# colnames meth, cluster, re
random_pred <- real_data_res[[8]] %>%
  mutate(cluster = rownames(.)) %>%
  pivot_longer(cols = -cluster, names_to = "meth", values_to = "re")

re_df <- data %>%
  filter(!train) %>%
  select(all_of(not_cna_var)) %>%
  cbind(real_data_res[[1]]) %>%
  pivot_longer(cols = -not_cna_var, names_to = "meth", values_to = "cluster") %>%
  left_join(random_pred, by = c("meth", "cluster")) %>%
  mutate(RACE = fct_relevel(RACE, c("White", "BAA", "Asian", "NHPI", "AIAN")))

re_df %>%
  filter(cluster != "Population Mean") %>%
  #filter(meth == Y_names[1]) %>%
  ggplot(aes(x = re, y = FRACTION_GENOME_ALTERED)) +
    geom_point(alpha = .08) +
    geom_smooth() +
    #facet_wrap(~ TISSUE_SOURCE_SITE) +
    theme_bw() +
    labs(y = "Dist. of absolute RE in test set")

re_df %>%
  filter(cluster != "Population Mean") %>%
  #filter(meth == Y_names[1]) %>%
  ggplot(aes(x = abs(re))) +
  geom_density() +
  facet_wrap(~ HISTOLOGICAL_DIAGNOSIS) +
  theme_bw() +
  labs(y = "Dist. of absolute RE in test set")

lm(abs(re) ~ FRACTION_GENOME_ALTERED, data = re_df) %>% summary()

# Histogram/density of RE
# Figure in paper
test_RE_plot <- re_df %>%
  filter(cluster != "Population Mean") %>%
  mutate(cluster = paste0("Cluster ", cluster)) %>%
  ggplot(aes(x = re, fill = cluster)) +
    geom_density(alpha = 1) +
    facet_wrap(~cluster) +
    theme_bw() +
    scale_fill_brewer(palette = "Set3") +
    #scale_color_brewer(palette = "Set3") +
    guides(fill = "none") +
    theme(strip.background = element_rect(fill = "transparent", color = "transparent")) +
    labs(x = "Assigned random effect", y = "Density")

ggsave("images/test_RE_plot.png", test_RE_plot, height = 4, width = 6,
       dpi = 800)

# Part R2 ----------------------------------------------------------------------

library(lme4)
library(partR2)

# We only need this on the training data

r2_lme4 <- sapply(meth_var, function(outcome_var) {
  lme4_formula <- as.formula(paste0(outcome_var, "~ ", X_string, "+ (1 | cluster)"))
  model <- lmer(lme4_formula, data = data[data$train, ])
  pr2_res <- partR2(model, data = data, partbatch = list(cna_var, c("RACE", "SEX", "AGE", "STAGE", "cancer_id")), nboot=NULL)
  as.matrix(pr2_res$R2[1:3, "estimate"])
})

saveRDS(r2_lme4, "result_matrices/r2_lme4.RDS")

r2_plot <- t(r2_lme4) %>%
  data.frame() %>%
  setNames(c("Model", "CNA variables", "Clinical variables")) %>%
  mutate(meth = rownames(.)) %>%
  pivot_longer(cols = -meth, names_to = "scope", values_to = "r2") %>%
  ggplot(aes(y = r2, fill = scope)) +
    geom_boxplot() +
    facet_wrap(~scope) +
    guides(fill = "none") +
    theme_bw() +
    theme(strip.background = element_rect(fill = "transparent", color = "transparent"),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank()) +
    labs(y = expression(R^2))

ggsave("images/r2_plot.png", r2_plot, height = 3, width = 4, dpi = 800)  










