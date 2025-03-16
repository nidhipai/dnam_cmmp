library(TCGAretriever)

# Download data ----------------------------------------------------------------

cancer_ids <- c("cesc", "luad")
raw_data <- list() # List of lists (cancers) of dataframes (clinical, meth, cna)

for (id in cancer_ids) {
  csid <- paste0(id, "_tcga")
  meth_feat <- paste0(id, "_tcga_methylation_hm450")
  mut_feat <- paste0(id, "_tcga_linear_CNA")
  #mut_feat <- "cesc_tcga_gistic" # discrete (binary 1/-1) version
  cases <- paste0(id, "_tcga_all")
  
  # Get clinical/covariate data
  clinical <- get_clinical_data(csid = csid, case_list_id = cases)
  
  # Get methylation data
  meth <- fetch_all_tcgadata(case_list_id = cases, gprofile_id = meth_feat, mutations = FALSE)
  
  # Get sCNA features
  cna <- fetch_all_tcgadata(case_list_id = cases, gprofile_id = mut_feat, mutations = FALSE)
  
  cancer_data <- list(clinical = clinical, meth = meth, cna = cna)
  raw_data[[id]] <- cancer_data
}

saveRDS(raw_data, "data/raw_data.RDS")