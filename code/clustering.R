data <- readRDS("data/bound_data_2.RDS")

meth_var <- grep("^meth_", colnames(data), value = TRUE) # meth var names
cna_var <- data %>%
  select(starts_with("cna")) %>%
  colnames()

# Functions  -------------------------------------------------------------------
# Gini isn't used in the paper

# Calculate gini coefficient according to Gini_max equation in mDNA paper
gini_index <- function(cluster, group) {
  # Remove those with NA race
  cluster <- cluster[!is.na(group)]
  group <- group[!is.na(group)]
  groups <- unique(group)
  clusters <- unique(cluster)
  sums <- lapply(clusters, function(c) { #TODO make sapply
    n_i <- sum(cluster == c) # Number of rows in cluster
    term <- lapply(groups, function(g) {
      n_ij <- sum((group == g) & (cluster == c))
      return(n_ij / n_i * (1 - n_ij / n_i))
    })
    return(sum(unlist(term)))
  })
  gini_stat <- max(unlist(sums))
  return(gini_stat)
}

# Clustering -------------------------------------------------------------------

# Get only methylation variables from training data
train_data <- data[data$train, ]
meth_only <- train_data[, meth_var]

# Choose number of clusters ----------------------------------------------------

k_max <- 15

set.seed(123)
# Get gap statistics to choose optimal clusters
gap_res <- clusGap(meth_only, function(x, k) kmeans(x, k), K.max = k_max, B = 100)
cluster_stat <- as.data.frame(gap_res$Tab[, c(3,4)]) # Extract gap stat, SE

saveRDS(cluster_stat, "result_matrices/cluster_stat_B100.RDS")

# Get gini indicies for race, cancer, and cells
for (k in 1:k_max) {
  kmeans_cluster <- kmeans(meth_only, k)$cluster
  cluster_stat[k, "race_gini"] <- gini_index(kmeans_cluster, train_data$RACE)
  cluster_stat[k, "cancer_gini"] <- gini_index(kmeans_cluster, train_data$cancer_id)
  cluster_stat[k, "race_cancer_gini"] <- gini_index(kmeans_cluster, train_data$race_cancer)
}

# # Plot Gap stat and gini stats for one cancer combo, for k = 1:10
# Data cleaning for plot
se_temp <- cluster_stat$SE.sim
k_plot_data <- cluster_stat[, -which(names(cluster_stat) == "SE.sim")]
k_plot_data$k <- rownames(k_plot_data)
k_plot_data$k <- factor(k_plot_data$k,
                        levels = as.character(seq(1:k_max)))
k_plot_data <- k_plot_data %>% pivot_longer(cols = !k, names_to = "metric", values_to = "value")
gap_rows <- k_plot_data$metric == "gap"
k_plot_data[gap_rows, "SE.sim"] <- se_temp # Save SE for gap stats
# Create lower and upper limits for gap stat
k_plot_data %>%
  rowwise() %>% 
  #filter(metric == "gap") %>%
  mutate(upper = sum(value, SE.sim),
         lower = sum(value, -1 * SE.sim)) %>%
  ggplot(aes(x = k, y = value, group = metric)) +
    geom_line(aes(color = metric)) +
    geom_errorbar(aes(ymin = lower, ymax = upper), width = .2) +
    labs(x = "k (number of clusters)", y = "Value") +
    ggtitle("Gap and gini statistics by number of clusters (kmeans)") +
    theme_bw()

# Find best k based on gap statistic
k_optimal <- maxSE(cluster_stat$gap, cluster_stat$SE.sim, 
                   method = "Tibs2001SEmax") # 6
# Calculate value of each k based on gini (not used for choosing k)
cluster_stat$gap_plus_gini <- cluster_stat$gap + cluster_stat$cancer_gini
cluster_gof <- max(cluster_stat$gap_plus_gini)

set.seed(123)
# Cluster the data according to best number of clusters using gap
data$cluster <- NULL
#kmeans_res <- kmeans(meth_only, k_optimal)
kmeans_res <- kmeans(meth_only, 6)
data[data$train, "cluster"] <- factor(kmeans_res$cluster)
train_data$cluster <- factor(kmeans_res$cluster)

saveRDS(data, "data/clustered_6.RDS")

# Exploring the clusters -------------------------------------------------------

set.seed(123) # for the jitter
# Cluster vs. cancer types
train_data %>%
  slice_sample(n = nrow(train_data)) %>% # Shuffle
  mutate(cancer_id = fct_relabel(cancer_id, toupper)) %>%
  mutate(RACE = race_abbrev[RACE]) %>%
ggplot(aes(x = cluster, y = RACE, color = cancer_id, shape = cancer_id)) +
  geom_jitter(width = .25, height = .20,  alpha = .4) +
  labs(x = "Cluster", y = "Race", color = "Cancer type", shape = "Cancer type") +
  theme_bw() +
  scale_y_discrete(limits = rev)

ggsave("images/clusters_jitter_race_cancer.png", 
       height = 3, width = 6, dpi = 600)

# Plot clusters by race
ggplot(train_data, aes(x = cluster, fill = RACE)) +
  geom_bar(position = "dodge") +
  theme_bw() +
  ggtitle("Clusters by race") + 
  labs(x = "Cluster index", y = "Count")

# Test for an association between clinical variables and cluster----
cont_vars <- c("AGE", "FRACTION_GENOME_ALTERED")
factor_vars <- c("RACE", "SEX", "STAGE", "CANCER_TYPE")
# Cont variables
p_values_cont <- sapply(cont_vars, function(var_name) {
  lm_res <- lm(as.formula(paste(var_name, " ~ cluster")), data = train_data)
  return(anova(lm_res)[1, "Pr(>F)"])
})

# Categorical variables
p_values_categorical <- sapply(factor_vars, function(var_name) {
  res <- chisq.test(table(train_data[, c(var_name, "cluster")]))
  return(res$p.value)
})

all_p_values <- c(p_values_cont, p_values_categorical)
all_p_values < (.05 / length(all_p_values))

all_p_values * length(all_p_values)


# Plot clinical variable vs. cluster for any interesting ones
# Fraction genome altered, ordered by mean
frac_genome_cluster <- train_data %>%
  #mutate(cross = paste(cluster, cancer_id)) %>%
ggplot(aes(x = cluster, y = FRACTION_GENOME_ALTERED, fill = cluster)) +
  geom_boxplot(position = "dodge") +
  theme_bw() +
  guides(fill = "none") +
  labs(x = "Cluster", y = "Fraction genome altered")

ggsave("images/frac_genome_cluster.png", frac_genome_cluster, 
       height = 3.5, width = 6, dpi = 800)

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

# PCs
pc <- prcomp(train_data[, meth_var], center = T, scale. = T)
train_data_pc <- train_data %>%
  cbind(pc$x)

ggplot(train_data_pc, aes(x = PC1, y = PC2, color = cluster)) +
  geom_point() +
  theme_bw() +
  labs(color = "Cluster")
ggsave("images/pc_clusters.png", height = 5, width = 6)

train_data_pc %>%
  select(cluster, starts_with("PC")) %>%
  pivot_longer(cols = -cluster, names_to = "PC", values_to = "value") %>%
  filter(PC %in% paste0("PC", 1:6)) %>%
  ggplot(aes(x = cluster, y = value, fill = cluster)) +
    facet_wrap(~ PC) +
    geom_boxplot() +
    theme_bw() +
    labs(x = "Cluster", y = "PC value") + 
    guides(fill = "none")
ggsave("images/pc_6_clusters.png", height = 6, width = 8)


















