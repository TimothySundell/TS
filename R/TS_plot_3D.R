#' Plot your Seurat object in 3D.
#'
#' @description
#' Will be coloured according to cluster idents with a scales::hue_pal() colouring scheme.
#' You need to have an already analysed Seurat object to be able to run the function.
#'
#' @param seurat_object The Seurat object you want to analyse. Defaults to 'default_seurat_object'
#' @param dims What Principal Components to use for calculating the UMAP.
#' @param group.by Variable to group cells by. Default is Seurat::Idents(), alternatives are any available metadata in the object.
#' @param colors A vector of colors to use for the group.by parameter. Default is hue_pal()(length of unique values in your group.by parameter)
#' @importFrom magrittr  %>%
#' @export
TS_plot_3D <- function(seurat_object = get(default_seurat_object),
                       dims = 1:20,
                       group.by = NULL,
                       colors = NULL) {

  require(plotly)
  require(scales)
  require(magrittr)
  require(Seurat)

  message("Calculating UMAP coordinates with 3 components")
  temp1 <- Seurat::RunUMAP(object = seurat_object, dims = dims, n.components = 3, verbose = F)

  if(is.null(group.by)){
    df <- data.frame(umap1 = temp1@reductions$umap@cell.embeddings[,1],
                     umap2 = temp1@reductions$umap@cell.embeddings[,2],
                     umap3 = temp1@reductions$umap@cell.embeddings[,3],
                     group.by = Seurat::Idents(temp1))
  } else{
    df <- data.frame(umap1 = temp1@reductions$umap@cell.embeddings[,1],
                     umap2 = temp1@reductions$umap@cell.embeddings[,2],
                     umap3 = temp1@reductions$umap@cell.embeddings[,3],
                     group.by = as.factor(temp1[[group.by]][,1]))
    }

  if(is.null(colors)){
  message("Drawing 3D plot")
  plotly::plot_ly(df, x = ~umap1, y = ~umap2, z = ~umap3, color = ~group.by, colors = scales::hue_pal()(length(levels(df$group.by)))) %>%
    plotly::add_markers() %>% # Don't know if needed
    plotly::layout(scene = list(xaxis = list(title = 'umap 1'),
                                yaxis = list(title = 'umap 2'),
                                zaxis = list(title = 'umap 3')))
  } else{
    message("Drawing 3D plot with custom colors.
Make sure that the number of supplied colors match the length of your group.by vector")
    plotly::plot_ly(df, x = ~umap1, y = ~umap2, z = ~umap3, color = ~group.by, colors = colors) %>%
      plotly::add_markers() %>% # Don't know if needed
      plotly::layout(scene = list(xaxis = list(title = 'umap 1'),
                                  yaxis = list(title = 'umap 2'),
                                  zaxis = list(title = 'umap 3')))
  }

}
