#' Export featureplot to pdf in A4 format
#'
#' @param seurat_object The Seurat object you wish to plot
#' @param file_prefix Prefix for the saved file, e.g. sample name
#' @param features Vector of features to plot
#' @param pt.size Dot size, default = 1.5
#' @param reduction Dimensionality reduction technique to use for plotting, default = "umap"
#' @param order Whether to plot positive cells last, e.g. on top, default = 1.5
#' @param cols Vector of 2 colors to use (low->high), default = c("lightgray", "red")
#' @export
TS_featureplot_to_pdf_A4 <- function(seurat_object, file_prefix, features, pt.size = 1.5, reduction = "umap", order = TRUE, cols = c("lightgray", "#FEFDE2", "#FEEAC7", "#FDD7AC", "#FDC591", "#FDB276", "#FC9F5B", "#FC8C40", "#FC8C40", "#EB753C", "#DA5D38", "#C94634", "#B72F30", "#A6172C", "#950028")) {

  require(Seurat)
  require(ggplot2)

  for(i in 1:length(feature_list)) {
    temp_plot <- Seurat::FeaturePlot(object = seurat_object, reduction = reduction, order = order, cols = cols, pt.size = pt.size, features = features[i])
    ggplot2::ggsave(filename = paste(file_prefix,"_", feature_list[i], ".pdf", sep = ""), plot = temp_plot, units = "in", width = 11.69, height = 8.27)
    print(paste("Exported figure for", feature_list[i]))
  }
  print(paste("Finished exporting", i, "figure(s)"))
}
