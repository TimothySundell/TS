#' Create column chart with pattern with ggplot2
#'
#' @param plot_data A dataframe containing the long format data which you want to plot with precomputed proportions.
#' @param x Variable to plot on x-axis.
#' @param y Variable to plot on y-axis.
#' @param fill Variable to fill on (group variable).
#' @param fill_colours Vector. Fill colours to use.
#' @param name String. Legend title. Defaults to "Sample"
#' @param ggplot_theme Optional. A specified ggplot_theme to use.
#' @param ylim  Optional, numeric. Use to specify max y_lim value.
#'
#' @return A tibble with columns: gene, group1, group2, p_value, p_adj
#' @importFrom magrittr %>%
#' @importFrom ggplot2 ggplot geom_col labs aes theme ylab ylim
#' @import ggpattern
#' @export

TS_plot_column_pattern_chart <- function(data,
                                 x,
                                 y,
                                 fill,
                                 fill_colours,
                                 name = "",
                                 ggplot_theme,
                                 ylim_max,
                                 border_colour,
                                 plot_title,
                                 plot_subtitle,
                                 pattern,
                                 pattern_colour = "black",
                                 pattern_fill = "black",
                                 pattern_density = 0.1,
                                 plot_n_sequences){

  count_information <- data %>%
    dplyr::mutate(count = as.integer(count)) %>%
    dplyr::group_by(grouping_variable) %>%
    dplyr::summarise(count = sum(count)) %>%
    dplyr::mutate(label = paste0(grouping_variable, " (n = ", count, ")")) %>%
    dplyr::pull(label)
  message((paste0(count_information, collapse = "\n")))

  if(plot_n_sequences){
    labels <- count_information
  } else{
    labels <- unique(data$grouping_variable)
  }

  plot <- ggplot2::ggplot(data = data, mapping = ggplot2::aes(x = gene_list, y = proportion, fill = grouping_variable)) +
    ggpattern::geom_col_pattern(
      mapping = ggplot2::aes(pattern = as.character(grouping_variable)),
      width = 0.8,
      position = position_dodge(),
      color = border_colour,
      pattern_colour = pattern_colour,
      pattern_fill = pattern_fill,
      pattern_density = pattern_density
      ) +
    {if(is.null(fill_colours)) ggplot2::scale_fill_grey(start = 0.9, end = 0.3, name = name, labels = labels)} +
    {if(!is.null(fill_colours)) ggplot2::scale_fill_manual(values = fill_colours, name = name, labels = labels)} +
    ggplot_theme +
    ggplot2::theme(
      axis.text.x = element_text(angle = 60,
                                 hjust = 1,
                                 vjust = 1,
                                 family = "Helvetica",
                                 color = "black"),
      axis.text.y = element_text(family = "Helvetica",
                                 color = "black"),
      text = element_text(family = "Helvetica",
                          color = "black"),
      axis.title.x = element_blank(),
      axis.ticks = element_line(color = "black"),
      legend.key.size = ggplot2::unit(0.5, "in")) +
    ggplot2::labs(title = plot_title,
                  subtitle = plot_subtitle) +
    ggplot2::ylab("%") +

  # Only set y limits if actually entered
  {if(!is.null(ylim_max)) ggplot2::ylim(0, ylim_max)}

  return(plot)
}
