#' Export featureplot to pdf without legend and axes
#'
#' @param seurat_object The Seurat object you wish to plot
#' @param file_prefix Prefix for the saved file, e.g. sample name
#' @param features Vector of features to plot
#' @param pt.size Dot size, default = 0.5
#' @param order Whether to plot positive cells last, e.g. on top, default = TRUE
#' @param cols Vector of 2 colors to use (low->high), default = c("lightgray", "red")
#' @param reduction Dimensionality reduction technique to use for plotting, default = "umap"
#' @param width Width (inches) of the resulting plot.
#' @param height Height (inches) of the resulting plot.
#' @export
TS_featureplot_to_pdf_no_leg_no_axes <-
  function (seurat_object,
            file_prefix,
            features,
            pt.size = 0.5,
            order = TRUE,
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
            width = 4.4,
            height = 3.78) {

    require(Seurat)
    require(ggplot2)

    for (i in 1:length(features)) {
      temp_plot <- Seurat::FeaturePlot(
        object = seurat_object,
        order = order,
        cols = cols,
        features = features[i],
        reduction = reduction,
        pt.size = pt.size
      ) +
        ggplot2::theme_classic(base_family = "Arial",
                      base_size = 15) +
        ggplot2::theme(
          plot.title = element_text(hjust = 0.5),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.line.x = element_blank(),
          axis.line.y = element_blank(),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.position = "none",
          validate = T
        )
      ggplot2::ggsave(
        filename = paste(file_prefix, "_", features[i], "_", "edited", ".pdf", sep = ""),
        plot = temp_plot,
        units = "in",
        width = width,
        height = height
      )
      print(paste("Exported figure for", features[i]))
    }
    print(paste("Finished exporting", i, "figure(s)"))
  }
