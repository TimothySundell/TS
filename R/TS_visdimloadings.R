#' Plot gene contribution to each principal component
#'
#' @param seurat_object The Seurat object you want to plot
#' @param dims What principal components to plot, default = 10
#' @export
TS_visdimloadings <- function(seurat_object = get(default_seurat_object), dims = 1:10) {

  require(Seurat)

  for (i in dims) {
    p1 <- Seurat::VizDimLoadings(object = seurat_object, dims = i)
    print(p1)
  }
}
