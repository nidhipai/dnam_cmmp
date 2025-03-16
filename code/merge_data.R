# Merge data -------------------------------------------------------------------
raw_data <- readRDS("data/raw_data.RDS")

cancer_ids <- c("cesc", "luad")
all_merged <- list()

# start for loop here
for (id in cancer_ids) {
  cna <- raw_data[[id]][["cna"]]
  meth <- raw_data[[id]][["meth"]]
  clinical <- raw_data[[id]][["clinical"]]
  
  rownames(cna) <- cna[, 2] # Make hugoGeneName the row id
  cna <- cna[, c(4:ncol(cna))] # Lose the meta information
  # B/c there are duplicates in hugoGeneName between the two omics, need to add
  # what set they are from
  rownames(cna) <- paste0("cna_", rownames(cna))
  cna <- t(cna) # Put patients/samples on rows, genes on cols
  
  # Repeat for meth
  rownames(meth) <- meth[, 2] # Make hugoGeneName the row id
  meth <- meth[, c(4:ncol(meth))] # Lose the meta information
  rownames(meth) <- paste0("meth_", rownames(meth))
  meth <- t(meth)
  
  rownames(clinical) <- clinical$sampleId # Make row names sample id
  clinical$sampleId <- NULL # Remove sample id from col
  
  # Merge data frames
  cna_meth <- merge(meth, cna, by = "row.names") # Inner join
  rownames(cna_meth) <- cna_meth$Row.names
  cna_meth$Row.names <- NULL
  merged <- merge(cna_meth, clinical, by = "row.names") # Also inner join
  rownames(merged) <- merged$Row.names
  merged$Row.names <- NULL
  
  # Save to list
  all_merged[[id]] <- merged
}

saveRDS(all_merged, file = "data/all_merged.RDS")
