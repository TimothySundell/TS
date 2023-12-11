#' Compare number of clusters for combinations of number of PCs and resolution
#'
#' @description
#' Old function. Not sure of how well it functions, or what it produces.
#' I think this function counts the number of clusters for each combination of
#' - 'dims = 1:i'
#' - 'resolution = seq(from = 0.1, to = 2, by = 0.1)'
#'
#' @param seurat_object Your input seurat object
#' @param PCs Number of the maximum PC you want to look at
#' @importFrom magrittr  %>%
#' @export
TS_compare_PC_and_resolution <- function(seurat_object, PCs) {

  require(Seurat)
  require(magrittr)

  df <- data.frame(NA_col = rep(NA, 20))
  rownames(df) <- seq(from = 0.1, to = 2, by = 0.1) # Create empty dataframe with 20 rows
  DefaultAssay(seurat_object) <- "SCT"
  counter = 1
  for (i in PCs) {
    temp1 <- Seurat::FindNeighbors(seurat_object, dims = 1:i, verbose = FALSE, nn.method = "rann") %>%
      Seurat::FindClusters(verbose = FALSE, resolution = seq(from = 0.1, to = 2, by = 0.1))
    temp2 <- c(length(table(temp1$SCT_snn_res.0.1)),
               length(table(temp1$SCT_snn_res.0.2)),
               length(table(temp1$SCT_snn_res.0.3)),
               length(table(temp1$SCT_snn_res.0.4)),
               length(table(temp1$SCT_snn_res.0.5)),
               length(table(temp1$SCT_snn_res.0.6)),
               length(table(temp1$SCT_snn_res.0.7)),
               length(table(temp1$SCT_snn_res.0.8)),
               length(table(temp1$SCT_snn_res.0.9)),
               length(table(temp1$SCT_snn_res.1)),
               length(table(temp1$SCT_snn_res.1.1)),
               length(table(temp1$SCT_snn_res.1.2)),
               length(table(temp1$SCT_snn_res.1.3)),
               length(table(temp1$SCT_snn_res.1.4)),
               length(table(temp1$SCT_snn_res.1.5)),
               length(table(temp1$SCT_snn_res.1.6)),
               length(table(temp1$SCT_snn_res.1.7)),
               length(table(temp1$SCT_snn_res.1.8)),
               length(table(temp1$SCT_snn_res.1.9)),
               length(table(temp1$SCT_snn_res.2)))
    df[counter] <- temp2 # Print to output dataframe
    colnames(df)[counter] <- paste0("PC_1:",i) # Correct column name for each component
    counter = counter + 1
  }
  return(df)
}
