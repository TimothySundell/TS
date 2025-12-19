# Adds test results to column charts

Helper function for TS_plot_IG_usage Adds test results from
TS_compute_pairwise_fisher() to TS_plot_column_chart() charts

## Usage

``` r
TS_add_sig_segments_and_labels(
  plot,
  compare_groups,
  gene_list_levels_plotting,
  df_stats,
  sig_distance_above_bar,
  sig_distance_between_significancies,
  sig_label_distance_from_line,
  sig_line_thickness,
  sig_label_text_size
)
```

## Arguments

- plot:

  GGplot2 object. The output from TS_plot_column_chart()

- compare_groups:

  Characer vector with groups to be compared

- df_stats:

  Dataframe. The output from TS_compute_pairwise_fisher()

- gene_list_level_plotting:

  Character vector with factor levels for x_axis categories (genes)

## Value

Returns a ggplot object with geom_segment() and geom_text() added

## Examples

``` r
plot <- TS_add_sig_segments_and_labels(plot = plot,
compare_groups = compare_groups,
gene_list_levels_plotting = gene_list_levels_plotting,
df_stats = df_stats)
#> Error in UseMethod("ggplot_build"): no applicable method for 'ggplot_build' applied to an object of class "function"
```
