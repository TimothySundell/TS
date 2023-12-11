#' Extract the mean and median values of Ig-genes from a list of PCA feature loadings
#'
#' @description
#' Takes as input the output from [TS::TS_extract_ranked_genes_per_PC]
#' @param PCA_loadings List, the output from [TS::TS_extract_ranked_genes_per_PC]
#' @export
TS_extract_ig_mean_median_rank_per_PC <- function(PCA_loadings) {
  temp_FUN <- function(x) {list(grep(pattern = "^IG[HKL][ADEGMVJC]", x = x) %>% mean() %>% round(digits = 0), grep(pattern = "^IG[HKL][ADEGMVJC]", x = x) %>% median())}
  temp1 <- sapply(PCA_loadings, temp_FUN)
  temp2 <- as.data.frame(as.list(temp1))
  row.names(temp2) <- c("mean", "median")
  temp2
}
