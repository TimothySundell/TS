#' Adds test results to column charts
#'
#' @description
#' Helper function for TS_plot_IG_usage
#' Adds test results from TS_compute_pairwise_fisher() to TS_plot_column_chart() charts
#'
#' @param plot GGplot2 object. The output from TS_plot_column_chart()
#' @param compare_groups Characer vector with groups to be compared
#' @param gene_list_level_plotting Character vector with factor levels for x_axis categories (genes)
#' @param df_stats Dataframe. The output from TS_compute_pairwise_fisher()
#' @return Returns a ggplot object with geom_segment() and geom_text() added
#'
#' @export
#'
#' @examples
#' plot <- TS_add_sig_segments_and_labels(plot = plot,
#' compare_groups = compare_groups,
#' gene_list_levels_plotting = gene_list_levels_plotting,
#' df_stats = df_stats)

TS_add_sig_segments_and_labels <- function(plot,
                                           compare_groups,
                                           gene_list_levels_plotting,
                                           df_stats,
                                           sig_distance_above_bar,
                                           sig_distance_between_significancies,
                                           sig_label_distance_from_line,
                                           sig_line_thickness,
                                           sig_label_text_size){
  plot_build <- ggplot2::ggplot_build(plot)$data[[1]] %>%
    dplyr::mutate(sample = compare_groups[base::match(fill, base::unique(fill))],
                  gene = base::rep(gene_list_levels_plotting, each = base::length(compare_groups)))

  plot_data_summarised <- plot_build %>%
    dplyr::select(sample, gene, x, y) %>%
    dplyr::group_by(sample, gene) %>%
    summarise(
      x = x[match(TRUE, !is.na(x))],
      y_max = max(y, na.rm = TRUE),
      .groups = "drop")

  segments <- df_stats %>%
    dplyr::left_join(plot_data_summarised, by = c("gene" = "gene", "group1" = "sample")) %>%
    dplyr::rename(x_start = x, y1 = y_max) %>%
    dplyr::left_join(plot_data_summarised, by = c("gene" = "gene", "group2" = "sample")) %>%
    dplyr::rename(x_end = x, y2 = y_max) %>%
    # compute a raw y position based on bar heights
    dplyr::mutate(
      y_seg_raw = pmax(y1, y2, na.rm = TRUE) + sig_distance_above_bar
    ) %>%
    # for each gene, stack brackets regardless of which groups are compared
    dplyr::group_by(gene) %>%
    dplyr::arrange(pmin(x_start, x_end), .by_group = TRUE) %>%
    dplyr::mutate(
      .dup_idx = dplyr::row_number(),
      y_seg = y_seg_raw + (.dup_idx - 1) * sig_distance_between_significancies
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(-.dup_idx, -y_seg_raw)

  plot_return <- plot +
    ggplot2::geom_segment(data = segments, mapping = aes(x = x_start, y = y_seg, xend = x_end, yend = y_seg), linewidth = sig_line_thickness, inherit.aes = F) +
    ggplot2::geom_text(data = segments, mapping = aes(x = ((x_start + x_end) / 2), y = y_seg + sig_label_distance_from_line, label = p_sig), size = sig_label_text_size, inherit.aes = F)

  return(plot_return)
}
