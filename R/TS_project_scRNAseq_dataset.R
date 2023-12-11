#' Projects one Seurat object onto another
#'
#' @description
#' A work in progress...
#' Will output resulting DimPlots with data from [Seurat::MapQuery]
#'
#' @param reference_object A Seurat object, will be used as reference.
#' @param query_object A Seurat object, will be used projected onto your reference object.
#' @param ref_name Name of reference object. Used for resulting plots.
#' @param query_name Name of query object. Used for resulting plots.
#' @export
TS_project_scRNAseq_dataset <- function(reference_object, query_object, ref_name, query_name) {

  # !! Change according to what defaultassay is required !!
  # "SCT" for single datasets
  # "integrated" for integrated datasets
  # Add if/else support for writing to file
  # Add auto-detection if reference object should have Defaultassay <- integrated || SCT

  require(Seurat)
  require(ggplot2)

  DefaultAssay(reference_object) <- "integrated" #this has already corrected cc genes
  DefaultAssay(query_object) <- "integrated" #this has already corrected cc genes

  ref <- Seurat::RunUMAP(reference_object, dims = 1:20, reduction = "pca", return.model = T)
  anchors <- Seurat::FindTransferAnchors(reference = ref, query = query_object, dims = 1:20, reference.reduction = "pca")
  query <- Seurat::MapQuery(anchorset = anchors, query = query_object, reference = ref, refdata = list(seurat_clusters = "seurat_clusters"), reference.reduction = "pca", reduction.model = "umap")

  p1 <- Seurat::DimPlot(ref, group.by = "seurat_clusters", label = T) + ggplot2::ggtitle(paste0(ref_name, " - reference"))
  p2 <- Seurat::DimPlot(query, group.by = "predicted.seurat_clusters", reduction = "ref.umap", label = T, repel = T) + ggplot2::ggtitle(paste0(query_name, " - query on ref. umap"))
  p3 <- Seurat::DimPlot(query, group.by = "predicted.seurat_clusters") + ggplot2::ggtitle(paste0(query_name, " - query on its own umap"))

  print(p1)
  print(p2)
  print(p1+p2)
  print(p3)
}
