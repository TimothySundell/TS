#' Output one or multiple featureplots
#'
#' @param features Vector of features to plot
#' @param seurat_object The Seurat object you wish to plot
#' @param pt.size Dot size, default = 1
#' @param order Whether to plot positive cells last, e.g. on top, default = TRUE
#' @param cols Vector of 2 colors to use (low->high), default = c("lightgray", "red")
#' @param reduction Dimensionality reduction technique to use for plotting, default = "umap"
#' @export
TS_featureplots <-
  function(features,
           seurat_object = get(default_seurat_object),
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

    require(Seurat)

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
  }
