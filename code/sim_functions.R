# Constants --------------------------------------------------------------------

# Weights of cancer incidence in US population
race_weights <- cesc_rates

# Generate sample ----------------------------------------------------------

generate_sample <- function(k, num_cna, n, j, sigma_a_diag, sigma_e_diag, 
                                z_factor, bias) {
  p <- num_cna + 11 # Number of covariates (num_cna + intercept + covariates)
  
  ## Generate covariates according to empirical distribution of CNA variables
  # Choose num_cna variables to pull from
  model_cna <- sample(cna_var, num_cna, replace = T)
  X_cna <- sapply(model_cna, function(cna) {
    sample(data[, cna], size = n, replace = T)
  })
  colnames(X_cna) <- paste0("X", 1:num_cna)
  
  ## Generate race variables
  # If bias (proportion white) is not specified, keep it default (unbiased)
  bias <- ifelse(is.null(bias), race_weights["White"], bias)
  biased_weights <- race_weights # Create copy
  # Update weights by bias
  for (w in names(race_weights)) {
    if (w == "White") biased_weights[w] <- bias
    else biased_weights[w] <- race_weights[w] * (1 - bias) / (1 - race_weights["White"])
  }
  # Include at least one from each race, then draw according to weights
  race_unshuff <- c(names(race_weights), sample(names(race_weights),
                                                prob = biased_weights,
                                                size = n - length(names(race_weights)),
                                                replace = T))
  
  # Shuffle order and turn into a factor
  race <- factor(sample(race_unshuff))
  
  # Generate cancer variables
  cancer <- factor(sample(data$cancer_id[!is.na(data$cancer_id)],
                          size = n, replace = T))
  
  # Generate sex based on cancer
  sex <- rep("Female", n) # Initiate with all Female
  prop_luad_male <- nrow(data[data$SEX == "Male", ]) / nrow(data[data$cancer_id == "luad", ])
  sex[cancer == "luad"] <- ifelse(rbinom(sum(cancer == "luad"), 1, prop_luad_male), "Male", "Female")
  sex <- factor(sex)
  
  # Generate age
  age <- sample(data$AGE[!is.na(data$AGE)], size = n, replace = T)
  
  # Generate cancer stage, ensuring at least one for each stage
  stage <- c(unique(data$STAGE), sample(data$STAGE[!is.na(data$STAGE)], 
                                        size = n - length(unique(data$STAGE)),
                                        replace = T)) %>% 
    sample() # Shuffle
  
  # Combine X variables
  X_covar <- model.matrix(~ race + cancer + sex + age + stage, 
                                data = data.frame(race, cancer, sex, age, stage))
  X <- cbind(X_covar, X_cna)
  
  # Generate fixed effects (beta)
  # beta_vec <- runif(p * k, min = -.3, max = .3)
  beta_vec <- rnorm(p * k, mean = 0, sd = .03)
  beta <- matrix(beta_vec, nrow = p, ncol = k)
  
  # Generate clusters
  group <- sample(1:j, size = n, replace = T)
  Z <- matrix(0, nrow = n, ncol = j)
  # Add 1 to designated col/group in each row/sample
  for (i in 1:n) {
    Z[i, group[i]] <- 1 * z_factor
  }
  
  # Generate random effects (alpha) - only random intercept for now
  # Original version: .07 all around
  #sigma_a <- matrix(.07, nrow = k, ncol = k)
  
  if (k <= 3) {
    sigma_a <- matrix(.07, nrow = k, ncol = k)
    diag(sigma_a) <- sigma_a_diag
    alpha <- mvrnorm(j, mu = rep(0, k), Sigma = sigma_a) # j x k
  } else { # Clumpy dependence
    # Randomly choose two blocks of k/4 to be correlated
    dep_meth <- sample(1:k, as.integer(k/4) * 2, replace = F) # Evenly size
    block_size <- ceiling(length(dep_meth) / 2)
    block1 <- dep_meth[1:block_size]
    block2 <- dep_meth[(block_size + 1):length(dep_meth)]
    indep_idx <- setdiff(1:k, c(block1, block2))
    
    block_size <- length(block1)
    sigma_dep <- matrix(.07, nrow = block_size, ncol = block_size)
    diag(sigma_dep) <- sigma_a_diag
    
    # Generate RE for each block seperately to avoid numerical issues
    block1_re <- mvrnorm(j, mu = rep(0, block_size), Sigma = sigma_dep) # j x block_size
    block2_re <- mvrnorm(j, mu = rep(0, block_size), Sigma = sigma_dep)
    num_ind_re <- k - 2 * block_size
    indep_re <- matrix(rnorm(j * num_ind_re, mean = 0, sd = sigma_a_diag), nrow = j) 
    
    alpha <- matrix(0, nrow = j, ncol = k)
    alpha[, block1] <- block1_re
    alpha[, block2] <- block2_re
    alpha[, indep_idx] <- indep_re
  }
  
  # Generate errors
  #sigma_e <- matrix(.07, nrow = k, ncol = k)
  sigma_e <- matrix(0, nrow = k, ncol = k)
  diag(sigma_e) <- sigma_e_diag
  epsilon <- mvrnorm(n, mu = rep(0, k), Sigma = sigma_e)
  
  # Calculate y's
  Y <- X %*% beta + Z %*% alpha + epsilon
  theta <- X %*% beta + Z %*% alpha # Mixed effect
  
  # Combine data
  sim_dat <- data.frame(Y, X_cna, race, cancer, sex, age, stage, group, theta)
  colnames(sim_dat) <- c(Y_names, cna_names, "race", 
                         "cancer", "sex", "age", "stage", "group", theta_names)
  
  sim_dat
}

# Methods for estimation of theta-----------------------------------------------

# Each method spits out sum(!data_train) x k predictions

# Method 1: Regression prediction
reg_pred <- function(data) {
  # Returns matrix of preds
  sapply(Y_names, function(Y_name) {
    # Fit model on train
    lm_res <- lm(as.formula(paste0(Y_name, "~", X_string)), data = data[data$train, ])
    # Get predictions on test set
    predict(lm_res, data[!data$train, ])
  })
}

# Method 2: k-means + CMMP
k_cmmp <- function(data, num_clusters) {
  # Get kmeans clusters on training data
  data[data$train, "cluster"] <- kmeans(data[data$train, Y_names], num_clusters)$cluster
  # dist_mat <- dist(data[data$train, Y_names])
  # hclust_res <- hclust(dist_mat, method = "ward.D")
  # data[data$train, "cluster"] <- cutree(hclust_res, k = num_clusters)
  
  # Use CMMP to get predictions on test data
  cmmp_res <- cmmp_core(data)
  
  # Return predictions and list of clusters for train, test, plus metrics
  list(cmmp_res[[2]], data[data$train, "cluster"], cmmp_res[[1]],
       cmmp_res[[4]], cmmp_res[[5]], cmmp_res[[6]])
}

# Method 2.5: CMMP with true clusters
group_cmmp <- function(data) {
  # The cluster is the group
  data[data$train, "cluster"] <- data[data$train, "group"]
  
  # Use CMMP to get predictions on test data
  cmmp_res <- cmmp_core(data)
  
  # Return predictions and list of clusters for train, test, plus metrics
  list(cmmp_res[[2]], data[data$train, "cluster"], cmmp_res[[1]], 
       cmmp_res[[4]], cmmp_res[[5]], cmmp_res[[6]])
}

# Method 3: \bar Y
bar_y <- function(data) {
  data[!data$train, Y_names] # Prediction is just y
}

# Metrics ----------------------------------------------------------------------

# MSE directly on theta (only doable in sim)
mse <- function(truth, pred) {
  mean((truth - pred)^2)
}
mse_all <- function(truth, pred) {
  if (nrow(truth) == 0) {
    return(NA)
  }
  sapply(1:ncol(truth), function(Y_index) {
    mse(truth[, Y_index], pred[, Y_index])
  })
}