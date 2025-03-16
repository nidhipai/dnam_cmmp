# Preprocessing ---------------------------------------------------------------
set.seed(12345)

all_merged <- readRDS("data/all_merged.RDS")

cancer_ids <- c("cesc", "luad") # cesc has to be first
# Each cancer has a different naming convention for the stage variable
stage_var_names <- list(cesc = "CLINICAL_STAGE", luad = "AJCC_PATHOLOGIC_TUMOR_STAGE")

# Define clinical variables of interest and secondary to interest
clinical_var <- c("AGE", "SEX", "RACE", "STAGE")
clinical_var_2 <- c("AGE", "SEX", "RACE", "STAGE", "CANCER_TYPE_DETAILED", "ETHNICITY",
                    "FRACTION_GENOME_ALTERED", "MUTATION_COUNT", "CANCER_TYPE")

# Variables that should be factors, of interest only (didn't look through all of them)
factor_vars <- c("RACE", "SEX", "STAGE", "CANCER_TYPE_DETAILED", "CANCER_TYPE")
cont_vars <- c("AGE", "FRACTION_GENOME_ALTERED", "MUTATION_COUNT")

# Remove -, *, and / from variable names in cna and meth
for (id in cancer_ids) {
  colnames(all_merged[[id]]) <- gsub("-", "_", colnames(all_merged[[id]]))
  colnames(all_merged[[id]]) <- gsub("/", "_", colnames(all_merged[[id]]))
  colnames(all_merged[[id]]) <- gsub("\\*", "_", colnames(all_merged[[id]]))
}

# Find variables that are in all datasets
# Can only intersect two variables at a time
for (id in cancer_ids) {
  merged <- all_merged[[id]]
  # Need to reconcile any columns with different names by dataset
  # Rename stage column as such
  stage_index <- which(colnames(merged) == stage_var_names[[id]])
  colnames(merged)[stage_index] <- "STAGE"
  
  # Remove patients missing key covariates
  # Drop patients that missing value for stage, age, race, sex
  merged <- merged %>%
    drop_na(RACE, AGE, STAGE, SEX)
  # We use an M-value transform, so need to remove 0 and 1, values which will go NA
  meth_var <- grep("^meth_", colnames(merged), value = TRUE)
  merged[, meth_var] <- sapply(1:length(meth_var), function(x) {
    meth_col <- merged[, meth_var[x]]
    ifelse((meth_col == 0) | (meth_col == 1),
           NA,
           meth_col)
  }) %>%
    setNames(meth_var)
  
  # Remove cols with missing data
  merged <- merged[, colSums(is.na(merged)) == 0] # TODO collect number here
  
  all_merged[[id]] <- merged
  
  if (id == "cesc") { # first one
    common_colnames <- colnames(all_merged[["cesc"]])
  } else {
    id_colnames <- colnames(merged)
    common_colnames <- intersect(id_colnames, common_colnames)
  }
}
# Remove variables that aren't in all datasets
for (id in cancer_ids) {
  all_merged[[id]] <- all_merged[[id]] %>% select(common_colnames)
}

# Processing to be done separately for each cancer
all_processed <- list()
for (id in cancer_ids) {
  merged <- all_merged[[id]]

  # # For filtering variables, we want to filter according to CESC traits
  if (id == "cesc") {
    ## Filter meth variables
    all_meth <- grep("^meth_", colnames(merged), value = TRUE)
    meth_in_merged <- merged[, all_meth] # Subset only methylation columns
    # Drop low variance methylation sites and those with missing data
    meth_stds <- sapply(meth_in_merged, function(x) sd(x, na.rm = TRUE))
    meth_missings <- sapply(meth_in_merged, function(x) sum(is.na(x)))
    paste0("Number of meth cols with missing values: ", sum(meth_missings > 0))
    # .2 is from https://www.nature.com/articles/s41392-019-0081-6
    meth_to_remove <- all_meth[meth_stds < .20 | meth_missings != 0] # 2379 remain
    
    ## Filter through sSCNA
    all_cna <- grep("^cna_", colnames(merged), value = TRUE)
    cna_in_merged <- merged[, all_cna] # Subset CNA columns
    
    cna_var <- all_cna #[!all_cna %in% cna_to_remove_dup] # Remaining CNA variables
    cna_in_merged <- merged[, cna_var]
    # Remove cols with variance lower than percentile threshold
    cna_stds <- sapply(cna_in_merged, function(x) sd(x))
    top_percentile <- quantile(cna_stds, probs = 0.85)
    cna_to_remove_var <- cna_var[cna_stds < top_percentile]
    
    cna_var <- all_cna[!all_cna %in% c(cna_to_remove_var)]
    cna_in_merged <- merged[, cna_var]
    # Remove one gene from pairs with high correlation
    cor_mat <- cor(cna_in_merged)
    # Returns cols by index to remove, min number of cols s.t. all corrs < cutoff
    # was .7 originally - might be able to move up without changing the number of variables
    cols_to_remove_cor <- caret::findCorrelation(cor_mat, cutoff = .70, verbose = F)
    # Get cols in cna_matrix
    cna_to_remove_corr <- colnames(cna_in_merged)[cols_to_remove_cor]
    cna_to_remove <- c(cna_to_remove_corr, cna_to_remove_var)#, cna_to_remove_dup)
    
    cols_to_remove <- c(meth_to_remove, cna_to_remove)
    cols_to_keep <- common_colnames[!common_colnames %in% cols_to_remove]
  }
  # Remove all cols not needed
  merged <- merged[, cols_to_keep]
  
  if (id == "cesc") {
    # Train test split
    total_samples <- nrow(merged)
    train_indicies <- sample(seq_len(total_samples), as.integer(.7 * total_samples),
                             replace = FALSE)
    merged$train <- seq_len(total_samples) %in% train_indicies
  } else {
    merged$train <- TRUE
  }

  # Add variable for cancer type
  merged$cancer_id <- id
  
  all_processed[[id]] <- merged
}

# Stack data frames for different cancers --------------------------------------

data <- bind_rows(all_processed[["cesc"]], all_processed[["luad"]])

# Replace methylation values with m-value transform
meth_var <- grep("^meth_", colnames(data), value = TRUE) # meth var names
data[, meth_var] <- lapply(data[, meth_var], function(x) log2(x / (1 - x)))

# Standardize methylation variables within cancer
meth_var <- grep("^meth_", colnames(data), value = TRUE)
data[, meth_var] <- data %>%
  select(patientId, cancer_id, starts_with("meth")) %>%
  pivot_longer(cols = -c(patientId, cancer_id), names_to = "meth_var_name", values_to = "value") %>%
  group_by(cancer_id, meth_var_name) %>%
  mutate(value = (value - mean(value)) / sd(value)) %>%
  ungroup() %>%
  pivot_wider(names_from = meth_var_name, values_from = value) %>%
  select(-c(patientId, cancer_id))

# Turn clinical variables into factors
data[, factor_vars] <- lapply(data[, factor_vars], factor)
data$cancer_id <- factor(data$cancer_id)

# Consolidate stage variables
data$STAGE <- forcats::fct_collapse(data$STAGE,
                                      StageI = c("Stage I", "Stage IA", "Stage IA1", "Stage IA2", 
                                                 "Stage IB", "Stage IB1", "Stage IB2"),
                                      StageII = c("Stage II", "Stage IIA", "Stage IIA1", 
                                                  "Stage IIA2", "Stage IIB"),
                                      StageIII = c("Stage III", "Stage IIIA", "Stage IIIB"),
                                      StageIV = c("Stage IV", "Stage IVA", "Stage IVB"))

# Abbreviate races and ethnicity
data$RACE <- recode_factor(data$RACE, 
                           `AMERICAN INDIAN OR ALASKA NATIVE` = "AIAN",
                           `ASIAN` = "Asian",
                           `BLACK OR AFRICAN AMERICAN` = "BAA",
                           `WHITE` = "White",
                           `NATIVE HAWAIIAN OR OTHER PACIFIC ISLANDER` = "NHPI")

# Create a variable which is cancer x race
data$race_cancer <- interaction(data$RACE, data$CANCER_TYPE)

saveRDS(data, "data/bound_data_2.RDS")

# Plots ------------------------------------------------------------------------

## Figure 1
# Should match up with race weights in the sim_functions file

race_prop_plot <- data %>%
  filter(cancer_id == "cesc") %>%
  group_by(RACE) %>%
  summarize(n = n()) %>%
  mutate(TCGA = n / sum(n)) %>%
  cbind(U.S. = cesc_rates) %>%
  select(-n) %>%
  mutate(RACE = race_abbrev[RACE]) %>%
  pivot_longer(cols = -RACE, names_to = "set", values_to = "rate") %>%
  ggplot(aes(x = set, y = rate, fill = RACE)) +
    geom_bar(stat = "identity") +
    theme_bw() +
    labs(fill = "Race", x = NULL) +
    theme(axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank()) +
    scale_fill_viridis(discrete = T, option = "magma")

ggsave("images/race_prop_plot.png", race_prop_plot, 
       height = 2.5, width = 3, dpi = 600)


# Extra plots not included in the paper

# Count of stage x cancer
ggplot(data, aes(x = STAGE, fill = cancer_id)) +
  geom_bar() +
  labs(x = "Clinical stage", y = "Count") +
  theme_bw() +
  ggtitle("Count of patients by stage and cancer") +
  #scale_fill_manual(values = colors) +
  guides(fill = guide_legend(title = "Cancer type"))

# Count of race x stage
ggplot(data, aes(x = STAGE, fill = RACE)) +
  geom_bar() +
  labs(x = "Clinical stage", y = "Count") +
  theme_bw() +
  ggtitle("Count of patients by stage and race") +
  #scale_fill_manual(values = colors) +
  guides(fill = guide_legend(title = "Race"))

# Count of cancer x race
# Shows lack of diversity in race
ggplot(data, aes(x = RACE, fill = cancer_id)) +
  geom_bar() +
  labs(x = "Race", y = "Count") +
  theme_bw() +
  ggtitle("Count of patients by race and cancer type") +
  #scale_fill_manual(values = colors) +
  guides(fill = guide_legend(title = "Cancer type"))

# Count of race x ethnicity
ggplot(data, aes(x = RACE, fill = ETHNICITY)) +
  geom_bar() +
  labs(x = "Race", y = "Count") +
  theme_bw() +
  ggtitle("Count of patients by race and ethnicity") + 
  #scale_fill_manual(values = colors) 
  guides(fill = guide_legend(title = "Ethnicity"))

# Get distribution of covariates in real data ----------------------------------

cna_var <- grep("^cna_", colnames(merged), value = TRUE)
cna_data <- data[, cna_var]

# Plot distribution of some cna variables
cna_data %>%
  pivot_longer(cols = everything(), names_to = "cna", values_to = "value") %>%
  filter(cna %in% sample(cna_var, size = 6)) %>%
  ggplot(aes(x = value)) +
  geom_histogram() +
  facet_wrap(~ cna) +
  theme_bw() +
  labs(x = "Value", y = "Count", title = "Distribution of CNA variables")

# Check mean and sd
cna_mean_sd <- sapply(cna_var, function(cna) {
  c(mean(data[, cna], na.rm = T), sd(data[, cna], na.rm = T))
})
cna_mean_sd <- cna_mean_sd %>% t()
colnames(cna_mean_sd) <- c("mean", "sd")


# Covariance between methylation variables in the real data
meth_var <- grep("^meth_", colnames(data), value = TRUE) # meth var names
meth_data <- data[, meth_var]
meth_cor <- cor(meth_data)
meth_cov <- var(meth_data)

corrplot(meth_cor[100:150, 100:150], order = "hclust")

corrplot(meth_cor[1:20, 1:20])

heatmap(meth_cov[1:20, 1:20])
heatmap(meth_cor[1:20, 1:20])











