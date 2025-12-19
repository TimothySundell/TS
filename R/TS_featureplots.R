#' Output one or multiple featureplots
#'
#' @param features Vector of features to plot.
#' @param seurat_object The Seurat object you wish to plot. Defaults to "default_seurat_object".
#' @param wrap_plots Whether to wrap plots into one by patchwork::wrap_plots(). Defaults to TRUE.
#' @param wrap_n_columns If [wrap_plots = T], how many columns do you want for the final plot.
#' @param dimplot_group.by What metadata to group cells by. Defaults to Seurat::Idents(seurat_object)
#' @param add_dimplot If [wrap_plots = T], whether to include Seurat::DimPlot() as the first panel. Defaults to TRUE
#' @param dimplot_cols Colours to use for the Seurat::DimPlot(). Defaults to scales::hue_pal()(length(levels(Seurat::Idents(seurat_object))))
#' @param dimplot_label Whether to label the DimPlot by [dimplot_group.by]. Defaults to TRUE.
#' @param expression_cols Vector of colors to use (low->high), defaults to a 15-step gradient from "lightgray" to "#950028"
#' @param reduction Dimensionality reduction technique to use for plotting, default = "umap".
#' @param pt.size Size of each dot, default = 1.
#' @param order Whether to plot positive cells last, e.g. on top, default = TRUE.
#' @param assay What assay to pull gene expression data from. Defaults to "RNA", as this is the convention.
#' @param reverse_y_scale Logical. Should the Y axis be reversed?
#' @param reverse_x_scale Logical. Should the X axis be reversed?
#' @export
TS_featureplots <- function(
    seurat_object = get(default_seurat_object),
    features,
    wrap_plots = TRUE,
    wrap_n_columns = NULL,
    add_dimplot = TRUE,
    dimplot_group.by = NULL,
    dimplot_cols = NULL,
    dimplot_label = TRUE,
    expression_cols = c(
      "lightgray",
      "#FEFDE2",
      "#FEEAC7",
      "#FDD7AC",
      "#FDC591",
      "#FDB276",
      "#FC9F5B",
      "#FC8C40",
      "#FC8C40",
      "#EB753C",
      "#DA5D38",
      "#C94634",
      "#B72F30",
      "#A6172C",
      "#950028"
    ),
    reduction = "umap",
    pt.size = 1,
    order = TRUE,
    assay = NULL,
    reverse_y_scale = FALSE,
    reverse_x_scale = FALSE) {

  if(!requireNamespace("Seurat", quietly = TRUE)){
    library(Seurat)
  }

  if(is.null(assay)) {
    Seurat::DefaultAssay(seurat_object) <- "RNA"
  } else{
    Seurat::DefaultAssay(seurat_object) <- assay
  }

  plot_list <- list()

  if(add_dimplot){
    p <- Seurat::DimPlot(
      object = seurat_object,
      pt.size = pt.size,
      reduction = reduction,
      group.by = dimplot_group.by,
      cols = dimplot_cols,
      label = dimplot_label,
      combine = TRUE
    ) +
      Seurat::NoLegend() +
      Seurat::NoAxes()
    if(reverse_y_scale) p <- p + ggplot2::scale_y_reverse()
    if(reverse_x_scale) p <- p + ggplot2::scale_x_reverse()
    plot_list[["Dimplot"]] <- p
  }

  for(gene in features){
    p <- Seurat::FeaturePlot(
      object = seurat_object,
      features = gene,
      cols = expression_cols,
      pt.size = pt.size,
      order = order,
      reduction = reduction,
      combine = TRUE
    ) +
      Seurat::NoAxes()
    if(reverse_y_scale) p <- p + ggplot2::scale_y_reverse()
    if(reverse_x_scale) p <- p + ggplot2::scale_x_reverse()
    plot_list[[gene]] <- p
  }

  if(wrap_plots){
    if(!requireNamespace("patchwork", quietly = TRUE)){
      library(patchwork)
    }

    return(patchwork::wrap_plots(plot_list, ncol = wrap_n_columns))
  }

  # If only one plot, return ggplot directly to allow chaining
  if(length(plot_list) == 1){
    return(plot_list[[1]])
  }
  # If plots are not wrapped
  return(plot_list)

}
