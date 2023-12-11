#' Export cell barcodes for each cluster
#'
#' @param seurat_object The Seurat object you wish to retrieve the identities from. Writes to your working directory.
#' @export
TS_export_cluster_ids <- function(seurat_object) {

  require(Seurat)
  require(utils)

  temp1 <- levels(seurat_object)
  temp2 <- (length(temp1)-1)
  for (i in 0:temp2) {
    temp <- Seurat::WhichCells(seurat_object, idents = i)
    utils::write.csv(temp, file = paste("cluster", i, "_ids.csv", sep = ""), row.names = FALSE)
    print(paste("Exported IDs for", i))
  }
}
