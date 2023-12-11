#' Export any plot as svg
#'
#' @param plot Needs a ggplot as input
#' @param filename Filename for the resulting plot, needs to be within "".
#' @export
TS_svg_from_plot <- function(plot, filename) {

  require(ggplot2)

  ggplot2::ggsave(filename = paste(filename, ".svg", sep = ""), plot = plot, units = "in", width = 8.85, height = 5.99)
  print(paste("Exported ",filename,".svg", sep = ""))
  print("Done")
}
