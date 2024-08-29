
build_phylo = function(covariate, fit_TRACE, tax_matrix, Xmodel){
  "----------------------------------------------------------------------------
  Return a phyloseq object to make Krona wheels plot
  - covariate: covariate from w
  - fit_TRACE: output from a TRACE model
  - tax_matrix: data frame with taxonoical tree classification of data 
  - Xmodel: model matrix covariate used for fitting fit_TRACE
  ----------------------------------------------------------------------------"
  
  p = ncol(fit_TRACE$Y)
  n = nrow(fit_TRACE$Y)
  # posterior mean of marginal species occurence
  pi = pnorm(Xmodel %*% fit_TRACE$coefficients)[,1:p]
  # otu_matrix
  otu_mat = as.data.frame(t(pi))
  rownames(otu_mat) = paste0("otu", 1:p)
  colnames(otu_mat) = paste0("s", 1:n)
  otu_mat = as.matrix(otu_mat)
  # construct samples_df
  covariate_fungi = cbind.data.frame(c(1:n), covariate_fungi)
  colnames(covariate_fungi)[1] = "sample"
  covariate_fungi = cbind.data.frame(covariate_fungi, rep("a",NROW(covariate_fungi)))
  colnames(covariate_fungi)[4] = "no_cov"
  samples_df = covariate_fungi %>% 
    tibble::column_to_rownames("sample") 
  rownames(samples_df) = paste0("s", 1:n)
  colnames(samples_df)[1]="site_id"
  # obtain phyloseq object
  OTU = otu_table(otu_mat, taxa_are_rows = TRUE) 
  samples = sample_data(samples_df)
  TAX = tax_table(tax_matrix)
  fungi_phylo = phyloseq(OTU, TAX, samples)
  return(fungi_phylo)
}

sprichnness_site = function(fit_TRACE, Xmodel){
  # compute sample-specific species richness
  p = ncol(fit_TRACE$Y)
  n = nrow(fit_TRACE$Y)
  pi = pnorm(Xmodel %*% fit_TRACE$coefficients)[,1:p]
  sprich = rowSums(pi)
  names(sprich) = NULL
  sprich
}
