# Modified the actual cmmp function to return the lme model

r.in <- as.formula("~ 1 | cluster")

cmmp_core <- function(data) {
  # Iterate over outcomes and get CMMP clusters/pred/interval
  cmmp_res_all <- lapply(Y_names,  function(outcome_var) {
    f.in <<- as.formula(paste0(outcome_var, "~ ", X_string))
    f.in <- as.formula(paste0(outcome_var, "~ ", X_string))
    # TODO THIS (<<-) IS A HACK - see notes
    cmmp_res <- cmmp(f.in, r.in,
         train = data[data$train, ],
         x.new = data[!data$train, X_names],
         y.new = data[!data$train, outcome_var],
         interval = TRUE)
    
    model <- cmmp_res[[4]]
    # Returns model; replace with estimated variance component
    cmmp_res[[4]] <- VarCorr(model)[1,1]
    cmmp_res[[5]] <- model$sigma^2 # Also add variance of the error
    # SNR
    design_mat <- model.matrix(as.formula(paste0("~", X_string)), data[data$train, ])
    xbeta_temp <- design_mat %*% t(as.matrix(coef(model)[1, ])) %>% var()
    cmmp_res[[6]] <- xbeta_temp / cmmp_res[[5]]
    cmmp_res[[7]] <- t(coef(model)[1, ]) # Beta values
    data_with_cmmp_cluster <- data
    # Random effects
    cmmp_res[[8]] <- ranef(model)
    cmmp_res
  })
  
  # Create a list of dataframes and list of estimated var comp.
  group.est_list <- lapply(1:length(cmmp_res_all), function(x) cmmp_res_all[[x]]$group.est)
  cmmp_cluster <- data.frame(do.call(cbind, group.est_list))
  mixed.pred_list <- lapply(1:length(cmmp_res_all), function(x) cmmp_res_all[[x]]$mixed.pred)
  cmmp_est <- data.frame(do.call(cbind, mixed.pred_list))
  random.pred_list <- lapply(1:length(cmmp_res_all), function(x) cmmp_res_all[[x]][[8]])
  random_pred <- data.frame(do.call(cbind, random.pred_list))
  colnames(cmmp_cluster) <- colnames(cmmp_est) <- colnames(random_pred) <- Y_names
  # Don't need the intervals for now and they take space
  #new.interval_list <- lapply(1:length(cmmp_res_all), function(x) cmmp_res_all[[x]]$new.interval)
  #colnames(cmmp_interval) <- paste0(rep(Y_names, each = 2), "_", c("Lower bound", "Upper bound"))
  #cmmp_interval <- data.frame(do.call(cbind, new.interval_list))
  
  var_comp_list <- sapply(1:length(cmmp_res_all), function(x) cmmp_res_all[[x]][[4]])
  var_comp_list <- as.numeric(var_comp_list)
  names(var_comp_list) <- Y_names
  sigma_list <- sapply(1:length(cmmp_res_all), function(x) cmmp_res_all[[x]][[5]])
  sigma_list <- as.numeric(sigma_list)
  names(sigma_list) <- Y_names
  snr_list <- sapply(1:length(cmmp_res_all), function(x) cmmp_res_all[[x]][[6]])
  snr_list <- as.numeric(snr_list)
  names(snr_list) <- Y_names
  # Beta values
  beta_matrix <- sapply(1:length(cmmp_res_all), function(x) cmmp_res_all[[x]][[7]])
  colnames(beta_matrix) <- Y_names
  rownames(beta_matrix) <- rownames(cmmp_res_all[[1]][[7]])
  #R2 values
  r2_all <- sapply(1:length(cmmp_res_all), function(x) cmmp_res_all[[x]][[9]])
  colnames(r2_all) <- Y_names
  rownames(r2_all) <- c("cna", "clinical")
    
  return(list(cmmp_cluster, cmmp_est, NULL, var_comp_list, sigma_list, snr_list,
              beta_matrix, random_pred, r2_all))
}

# CMMP outputs a cluster prediction for each outcome, including "pop mean",
# and we want one prediction per obsevation
cmmp_clust_to_cluster <- function(cmmp_clusters) {
  for (i in 1:nrow(cmmp_clusters)) {
    tab <- table(unlist(cmmp_clusters[i, ]))
    if ((length(tab) == 1) && (names(tab) == "Population Mean")) {
      cmmp_clusters[i, "cluster"] <- -1
    } else {
      # List of cluster names that are a mode, other than "Population Mean"
      tab <- tab[names(tab) != "Population Mean"]
      max_freq <- max(tab)
      max_clusters <- names(which(tab == max_freq))
      cmmp_clusters[i, "cluster"] <- sample(max_clusters, 1)
    }
  }
  return(cmmp_clusters$cluster)
}

# Actual cmmp function ---------------------------------------------------------

cmmp <- function(f.in, r.in, train, x.new, y.new, x.fut = NULL, match.train = NULL, a1 = NULL, n.new = NULL, interval = NULL, ...){
  if (class(f.in) != "formula" || class(r.in) != "formula")
    stop("'f.in and r.in must be formula")
  cmmp.result <- NULL
  if(is.null(match.train)){match.train = FALSE}
  if(is.null(n.new)){ n.new = 1}
  if(is.null(a1)){ a1 = 0.05}
  r.ch <- as.character(r.in)
  f.ch <- as.character(f.in)
  if (length(r.ch) > 1) {r.ch <- r.ch[length(r.ch)]}
  f.in <- as.formula(f.in)
  r.in <- as.formula(r.in)
  group.var0 <- sub("^[^|]*", "", r.ch)
  group.var1 <- gsub("|", "", group.var0, fixed = TRUE)
  group.var <- gsub(" ", "", group.var1, fixed = TRUE)
  
  ## mixed effect model
  obj <- lme(f.in, train, r.in)
  var.e.est <- as.numeric(VarCorr(obj)[1,1])	  	## est of var.e
  var.a.est <- as.numeric(VarCorr(obj)[2,1])	  	## est of var.a
  ran.names <- as.factor(rownames(ranef(obj)))
  I.est <- theta.cmmp <- c()
  if(is.null(dim(x.new))) {
    x.new <- matrix(x.new, length(x.new),1)
    x.name0 <- f.ch[length(f.ch)]
    x.name1 <- sub("^[^+]*", "", x.name0)
    x.name2 <- gsub("+", "", x.name1, fixed = TRUE)
    colnames(x.new) <- gsub(" ", "", x.name2, fixed = TRUE)
  }
  ## apply cmmp method
  for (i in 1:dim(x.new)[1]){
    y.temp <- y.new[i]
    if(is.data.frame(x.new)){
      x.new.df <- data.frame(x.new[i,], g = ran.names, row.names= NULL)
    }else{
      x.new.df <- data.frame(matrix(x.new[i,], byrow= T, nrow = length(ran.names), ncol = dim(x.new)[2]), g = ran.names, row.names= NULL)
    }
    colnames(x.new.df) <- c(colnames(x.new), group.var)
    theta.pred <- predict(obj, x.new.df, level = 0:1)
    if(match.train == TRUE){
      theta.temp <- theta.pred[,dim(theta.pred)[2]]   ## matched version
    }else{
      theta.temp <- c(theta.pred[,dim(theta.pred)[2]],theta.pred[1,2])  ## unmatched version
    }
    sp.temp <- theta.temp^2-2*theta.temp*y.temp
    I.temp <- match(min(sp.temp),sp.temp)
    theta.cmmp[i] <- theta.temp[I.temp]
    I.est[i] <- I.temp
  }
  ## matched or unmatched? group name changes
  if(match.train == TRUE){
    group.all <- as.character(ran.names)
  }else{
    group.all <- c(as.character(ran.names),"Population Mean")
  }
  group.est <- group.all[I.est]
  cmmp.result$group.est <- group.est    ## estimated group name
  cmmp.result$mixed.pred <- theta.cmmp  ## prediction of mixed effect
  
  ## prediction interval of new observations (only for NER model)
  if (!is.null(interval)){
    new.length <- qnorm((1-a1/2))*sqrt(var.e.est/n.new)
    new.interval <- cbind(theta.cmmp - new.length, theta.cmmp + new.length)
    colnames(new.interval) <- c("Lower bound", "Upper bound")
    cmmp.result$new.interval <- new.interval
  }
  
  if(!is.null(x.fut)){
    ## for the prediction of future observation (x.new, y.new) is assumed to be the intermediate data that come from the same group with the future observation
    mu.fut.pred <- predict(obj, x.fut, level = 0)[1]
    alpha.all <- rbind(as.matrix(random.effects(obj)),0)
    if(dim(alpha.all)[2]>1){
      alpha.est <- apply(alpha.all[I.est,],2,mean)
      a.est.temp <- as.matrix(random.effects(obj))
      a.est.temp1 <- rbind(a.est.temp,alpha.est)
      obj2 <- obj
      obj2$coefficients$random$Subject <- a.est.temp1
      x.fut$group <- "alpha.est"
      colnames(x.fut)[length(x.fut)] <- group.var
      y.fut.pred <- as.numeric(predict(obj2, x.fut, Subject = "alpha.est"), level = 1)
    }else{
      alpha.est <- mean(alpha.all[I.est])
      y.fut.pred <- mu.fut.pred + alpha.est
    }
    cmmp.result$fut.pred <- y.fut.pred
    ## prediction interval of future observation (only for NER model)
    if (!is.null(interval)){
      fut.length <- qnorm((1-a1/2))*sqrt(var.e.est*(1+1/n.new))
      cmmp.result$fut.interval <- c(y.fut.pred - fut.length, y.fut.pred + fut.length)
    }
  }
  cmmp.result$obj <- obj
  
  return(cmmp.result)
  
}

