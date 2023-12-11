#' Extract the 'rank' for each of the Ig-genes
#'
#' @description
#' Takes as input the output from [TS::TS_extract_ranked_genes_per_PC]
#' @param PCA_loadings The output from [TS::TS_extract_ranked_genes_per_PC]
#' @export
TS_extract_ig_ranks_per_PC <- function(PCA_loadings) {
  temp_FUN <- function(x) {grep(pattern = "^IG[HKL][ADEGMVJC]", x = x)}
  temp1 <- sapply(PCA_loadings, temp_FUN)
  temp2 <- as.data.frame(temp1)
  temp2
}
