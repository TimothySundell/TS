#' Calculate cellnumbers per cluster for each combination of 'number or PCs' and 'resolution'
#'
#' @param seurat_object Input seurat object
#' @param PCs Number of maximum PC you want to look at
#' @importFrom magrittr  %>%
#' @export
TS_cellnumbers_per_PC_and_resolution <- function(seurat_object, PCs) {

  require(Seurat)
  require(magrittr)

  df <- data.frame(NA_col = rep(NA, 20))
  rownames(df) <- seq(from = 0.1, to = 2, by = 0.1) # Create empty dataframe with 20 rows
  DefaultAssay(seurat_object) <- "SCT"
  counter = 1
  for (i in PCs) {
    temp1 <- Seurat::FindNeighbors(seurat_object, dims = 1:i, verbose = FALSE, nn.method = "rann") %>%
      Seurat::FindClusters(verbose = FALSE, resolution = seq(from = 0.1, to = 2, by = 0.1))
    temp2 <- list((table(temp1$SCT_snn_res.0.1)),
                  (table(temp1$SCT_snn_res.0.2)),
                  (table(temp1$SCT_snn_res.0.3)),
                  (table(temp1$SCT_snn_res.0.4)),
                  (table(temp1$SCT_snn_res.0.5)),
                  (table(temp1$SCT_snn_res.0.6)),
                  (table(temp1$SCT_snn_res.0.7)),
                  (table(temp1$SCT_snn_res.0.8)),
                  (table(temp1$SCT_snn_res.0.9)),
                  (table(temp1$SCT_snn_res.1)),
                  (table(temp1$SCT_snn_res.1.1)),
                  (table(temp1$SCT_snn_res.1.2)),
                  (table(temp1$SCT_snn_res.1.3)),
                  (table(temp1$SCT_snn_res.1.4)),
                  (table(temp1$SCT_snn_res.1.5)),
                  (table(temp1$SCT_snn_res.1.6)),
                  (table(temp1$SCT_snn_res.1.7)),
                  (table(temp1$SCT_snn_res.1.8)),
                  (table(temp1$SCT_snn_res.1.9)),
                  (table(temp1$SCT_snn_res.2)))
    df[counter] <- temp2 # Print to output dataframe
    colnames(df)[counter] <- paste0("PC_1:",i) # Correct column name for each component
    counter = counter + 1
  }
  return(df)
}
