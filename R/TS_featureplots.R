#' Output one or multiple featureplots
#'
#' @param features Vector of features to plot.
#' @param seurat_object The Seurat object you wish to plot.
#' @param wrap_plots Whether to wrap plots into one by patchwork::wrap_plots().
#' @param wrap_columns If [wrap_plots = T], how many columns do you want for the final plot.
#' @param cols Vector of colors to use (low->high), defaults to a 15-step gradient from "lightgray" to "#950028"
#' @param reduction Dimensionality reduction technique to use for plotting, default = "umap".
#' @param pt.size Size of each dot, default = 1.
#' @param order Whether to plot positive cells last, e.g. on top, default = TRUE.
#' @export
TS_featureplots <-
  function(features,
           seurat_object = get(default_seurat_object),
           wrap_plots = FALSE,
           wrap_columns = NULL,
           cols = c(
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
           order = T) {
    library(Seurat)

    if (!wrap_plots) {
      for (i in 1:length(features)) {
        temp_plot <-
          Seurat::FeaturePlot(
            object = seurat_object,
            features = features[i],
            cols = cols,
            pt.size = pt.size,
            order = order,
            reduction = reduction
          )
        print(temp_plot)
      }
    } else{
      library(patchwork)
      plot_list <- list()
      for (i in 1:length(features)) {
        plot_list[[i]] <-
          Seurat::FeaturePlot(
            object = seurat_object,
            features = features[i],
            cols = cols,
            pt.size = pt.size,
            order = order,
            reduction = reduction
          )
      }
        print(patchwork::wrap_plots(plot_list, ncol = wrap_columns))

    }

  }
