#' Projects one Seurat object onto another and exports resulting dataframe
#'
#' @description
#' A work in progress...
#' Will output resulting dataframe from [Seurat::MapQuery]
#'
#' @param reference_object A Seurat object, will be used as reference.
#' @param query_object A Seurat object, will be used projected onto your reference object.
#' @param dims Number of Principal Components to use for the UMAP. Defaults to '1:20'
#' @param reference_object_assay Assay to use for reference Seurat object. Defaults to 'integrated'
#' @param query_object_assay Assay to use for query Seurat object. Defaults to 'integrated'
#' @param normalization.method Normalization method used by FindTransferAnchors. Alternatives are 'LogNormalize' or 'SCT'. Defaults to 'SCT'.
#' @param plot_3D Whether to calculate 3 dimensions for the UMAP model. Defaults to 'FALSE'
#' @export
TS_project_scRNAseq_dataset_export_query <- function(reference_object, query_object, dims = 1:20, reference_object_assay = "integrated", query_object_assay = "integrated", normalization.method = "SCT", plot_3D = F) {

  # !! Change according to what defaultassay is required !!
  # "SCT" for single datasets
  # "integrated" for integrated datasets
  # Add if/else support for writing to file
  # Add auto-detection if reference object should have Defaultassay <- integrated || SCT

  library(Seurat)
  library(ggplot2)


  if(plot_3D) {
    cat(paste0("\nProjection based on dimensions ", deparse(dims), " of the reference PCA\n"))

    DefaultAssay(reference_object) <- reference_object_assay
    DefaultAssay(query_object) <- query_object_assay

    ref <- Seurat::RunUMAP(reference_object, dims = dims, reduction = "pca", return.model = T, n.components = 3)
    anchors <- Seurat::FindTransferAnchors(reference = ref, query = query_object, dims = dims, reference.reduction = "pca")
    query <- Seurat::MapQuery(anchorset = anchors, query = query_object, reference = ref, refdata = list(seurat_clusters = "seurat_clusters"), reference.reduction = "pca", reduction.model = "umap")

    return(query)
  }

  else{
    cat(paste0("\nProjection based on dimensions ", deparse(dims), " of the reference PCA\n"))

    DefaultAssay(reference_object) <- reference_object_assay
    DefaultAssay(query_object) <- query_object_assay

    ref <- Seurat::RunUMAP(reference_object, dims = dims, reduction = "pca", return.model = T)
    anchors <- Seurat::FindTransferAnchors(reference = ref, query = query_object, dims = dims, reference.reduction = "pca", normalization.method = normalization.method)
    query <- Seurat::MapQuery(anchorset = anchors, query = query_object, reference = ref, refdata = list(seurat_clusters = "seurat_clusters"), reference.reduction = "pca", reduction.model = "umap")

    return(query)
  }


}
