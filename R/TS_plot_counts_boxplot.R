#' Plot counts for a gene in a Seurat object as boxplot
#'
#' @param seurat_object The Seurat object you want to analyse
#' @param slot Slot to plot. Defaults to "counts"
#' @param gene Gene to plot
#' @export
TS_plot_counts_boxplot <- function(seurat_object = get(default_seurat_object), slot = "counts", gene) {

  require(Seurat)
  require(ggplot2)

  temp_plot <- Seurat::FetchData(seurat_object, slot = slot, vars = c(gene, "ident"))
  ggplot2::ggplot(data = temp_plot, aes(y = get(gene), x = ident)) +
    ggplot2::geom_boxplot() +
    ggplot2::theme_bw(base_family = "Arial") +
    ggplot2::labs(title = paste("Counts distribution for", gene, "in", deparse(substitute(seurat_object))), x = "Cluster", y = "Number of counts")
}
