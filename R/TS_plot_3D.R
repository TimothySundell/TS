#' Plot your Seurat object in 3D.
#'
#' @description
#' Will be coloured according to cluster idents with a [scales::hue_pal()] colouring scheme.
#' You need to have an already analysed Seurat object to be able to run the function.
#'
#' @param seurat_object The Seurat object you want to analyse. Defaults to 'default_seurat_object'
#' @param dims What Principal Components to use for calculating the UMAP.
#' @importFrom magrittr  %>%
#' @export
TS_plot_3D <- function(seurat_object = get(default_seurat_object), dims = 1:20) {

  require(plotly)
  require(scales)
  require(magrittr)
  require(Seurat)

  temp1 <- Seurat::RunUMAP(object = seurat_object, dims = dims, n.components = 3, verbose = F)
  df <- data.frame(umap1 = temp1@reductions$umap@cell.embeddings[,1],
                   umap2 = temp1@reductions$umap@cell.embeddings[,2],
                   umap3 = temp1@reductions$umap@cell.embeddings[,3],
                   cell = Seurat::Idents(temp1))
  plotly::plot_ly(df, x = ~umap1, y = ~umap2, z = ~umap3, color = ~cell, colors = scales::hue_pal()(length(levels(df$cell)))) %>%
    plotly::add_markers() %>% # Don't know if needed
    plotly::layout(scene = list(xaxis = list(title = 'umap 1'),
                                yaxis = list(title = 'umap 2'),
                                zaxis = list(title = 'umap 3')))
}
