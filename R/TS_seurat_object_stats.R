#' Prints stats about your Seurat object
#'
#' @param seurat_object Your Seurat object
#' @export
TS_seurat_object_stats <- function(seurat_object) {

  require(Seurat)
  require(glue)

  print(glue::glue("Dataset: ", {substitute(seurat_object)}, "\n\n", "Cells per cluster: "))
  print(table(Seurat::Idents(seurat_object)))
  print(glue::glue("\n\nTotal number of cells: ", {sum(table(Seurat::Idents(seurat_object)))}))
}
