#' Export identities for each cluster, which are not numerical
#'
#' @param seurat_object The Seurat object you wish to retrieve the identities from. Writes to working directory.
#' @export
TS_export_named_ids <- function(seurat_object) {

  require(Seurat)
  require(utils)

  temp1 <- levels(seurat_object)
  temp2 <- (as.vector(temp1))
  for (i in temp2) {
    temp <- Seurat::WhichCells(seurat_object, idents = i)
    utils::write.csv(temp, file = paste(i, "_ids.csv", sep = ""), row.names = FALSE)
    print(paste("Exported IDs for", i))
  }
}
