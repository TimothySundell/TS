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
#' @export
TS_featureplots <- function(features,
                            seurat_object = get(default_seurat_object),
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
                            assay = NULL) {
  library(Seurat)

  if(is.null(assay)) {
    Seurat::DefaultAssay(seurat_object) <- "RNA"
  } else{
    Seurat::DefaultAssay(seurat_object) <- assay
  }

  if(wrap_plots) {
    suppressPackageStartupMessages(library(patchwork))
    plot_list <- list()
    if(add_dimplot) {
      plot_list[["Dimplot"]] <- Seurat::DimPlot(
          object = seurat_object,
          pt.size = pt.size,
          reduction = reduction,
          group.by = dimplot_group.by,
          cols = dimplot_cols,
          label = dimplot_label,
          combine = TRUE
        ) +
        Seurat::NoLegend() + Seurat::NoAxes()
    } else{}
    for (i in 1:length(features)) {
      plot_list[[features[i]]] <-
        Seurat::FeaturePlot(
          object = seurat_object,
          features = features[i],
          cols = expression_cols,
          pt.size = pt.size,
          order = order,
          reduction = reduction,
        ) + Seurat::NoAxes()
    }
    print(patchwork::wrap_plots(plot_list, ncol = wrap_n_columns, guides = "collect"))
  } else{
    for (i in 1:length(features)) {
      temp_plot <-
        Seurat::FeaturePlot(
          object = seurat_object,
          features = features[i],
          cols = expression_cols,
          pt.size = pt.size,
          order = order,
          reduction = reduction
        )
      print(temp_plot)
    }
  }
}

