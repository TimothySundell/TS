#' Create column chart with ggplot2
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
#' @export

TS_plot_column_chart <- function(data,
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
  plot <- ggplot2::ggplot() +
    ggplot2::geom_col(
      data = data,
      mapping = ggplot2::aes(x = gene_list,
                    y = proportion,
                    fill = grouping_variable),
      width = 0.8,
      position = position_dodge(),
      color = border_colour) +
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
      axis.ticks = element_line(color = "black")) +
    ggplot2::labs(title = plot_title,
                  subtitle = plot_subtitle) +
    ggplot2::ylab("%") +

  # Only set y limits if actually entered
  {if(!is.null(ylim_max)) ggplot2::ylim(0, ylim_max)}

  return(plot)
}


# Implement something similar to this?
# To add number of sequences information
# ggplot2::scale_fill_manual(values = fill_colours,
#                            name = "Sample",
#                            labels = c(paste0(sample1, "\n(n = ", intermediate_export[1, 6], ")"), paste0(sample2, "\n(n = ", intermediate_export[1, 7], ")")))
