#' Extract an ordered list of feature importance for each of the first 50 principal components
#'
#' @param seurat_object The Seurat object you wish to extract PCA feature loadings from
#' @export
TS_extract_ranked_genes_per_PC <- function(seurat_object) {

  require(dplyr)
  require(Seurat)

  temp1 <- Seurat::Loadings(object = seurat_object, reduction = "pca")  # Extract loadings for PCA
  temp2 <- as.data.frame(abs(temp1))  # Make values absolute

  df <- data.frame(NA_col = rep(NA, length(Seurat::VariableFeatures(seurat_object))))     # Create empty dataframe
  for (i in 1:50) {  # Loop to extract VariableFeature loadings for each PC - number of PCs=50
    d <- rownames(dplyr::slice_max(temp2, get(paste0("PC_",i)), n = length(Seurat::VariableFeatures(seurat_object)))) # Extract loadings
    df[i] <- d # Print to output dataframe
    colnames(df)[i] <- paste0("PC_",i) # Correct column name for each component
  }
  return(df)
}
