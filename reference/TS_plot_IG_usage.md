# Plots Ig usage between two groups using Fisher's Exact Test

Plot the Ig usage from two groups and perform Fisher's Exact Test to
find significant differences.

## Usage

``` r
TS_plot_IG_usage(
  input_file,
  compare_groups,
  compare_what = "IGHV",
  grouping_variable = "named_clusters",
  export_pdf = F,
  return_only_data = F,
  return_only_sig_data = F,
  file_name_suffix = NULL,
  custom_gene_levels = NULL,
  file_type = "IMGT",
  ylim_max = NULL,
  p.adj.method = "hochberg",
  fill_colours = NULL,
  border_colour = "black",
  ggplot_theme = ggplot2::theme_classic(),
  plot_width = 15,
  plot_height = 5,
  remove_genes_present_in_one_dataset = F,
  plot_zeroes = T,
  plot_subtitle = "% of max",
  custom_plot_title = NULL,
  legend_title = "",
  plot_n_sequences = F,
  add_pattern = F,
  pattern = grouping_variable,
  pattern_colour = "black",
  pattern_fill = "black",
  pattern_density = 0.1,
  add_significancies = T,
  plot_only_significant_comparisons = T,
  sig_distance_above_bar = 1,
  sig_distance_between_significancies = 0.8,
  sig_line_thickness = 1,
  sig_label_distance_from_line = 0.2,
  sig_label_text_size = 6,
  return_usage_proportions = F
)
```

## Arguments

- input_file:

  A dataframe (e.g. from 10X Genomics Cellranger VDJ output) that
  includes a grouping column.

- compare_groups:

  Character vector. Names of group values to compare.

- compare_what:

  Genes to compare (e.g. "IGHV", "IGKV"). Defaults to "IGHV".

- grouping_variable:

  Column used to group data. Defaults to "named_clusters"

- export_pdf:

  Logical. If TRUE, exports the plot to a PDF. Defaults to FALSE

- return_only_data:

  Logical. If TRUE, returns only test results. Defaults to FALSE.

- return_only_sig_data:

  Logical. If TRUE, returns only significant test results. Defaults to
  FALSE.

- file_name_suffix:

  Character string. Optional suffix to be added to the exported pdf.

- custom_gene_levels:

  Character vector. Optional, vector used to set gene list levels.

- file_type:

  Set type of file/annotation input. Defaults to "cellranger". Options
  are "IMGT" and "AIRR".

- ylim_max:

  Numeric. Optional upper limit for y-axis. ggplot2 will adjust the
  y-axis scale to fit the data.

- p.adj.method:

  Method for adjust for multiple tests. Defaults to 'Hochberg'.
  Alternative is 'FDR'.

- fill_colours:

  A character vector of fill colours. Defaults to greyscale.

- border_colour:

  Colour of bar borders. Defaults to "black".

- ggplot_theme:

  Optional ggplot2 theme to apply. Defaults to "theme_classic".

- plot_width:

  Numeric. Width of exported plot in inches. Used if 'export_pdf =
  TRUE'. Defaults to 15.

- plot_height:

  Numeric. Height of exported plot in inches. Used if 'export_pdf =
  TRUE'. Defaults to 5.

- remove_genes_present_in_one_dataset:

  Logical. If TRUE, removes genes not shared by both groups. Defaults to
  FALSE

- plot_zeroes:

  Logical. If FALSE, genes with 0 counts in one group are not plotted,
  resulting in a wider bar for the other group.

- plot_subtitle:

  Character string. Plot subtitle. Defaults to "% of max".

- legend_title:

  Character string. Optional. Title of the legend.

- plot_n_sequences:

  Logical. Should the number of sequences used for each group be plotted
  in the legend? Defaults to FALSE.

- add_pattern:

  Logical. Should the bars be plotted with pattern? Defaults to FALSE.
  Requires `ggpattern`

- pattern:

  What to base the pattern on. Defaults to the grouping_variable.

- pattern_colour:

  String. Colour of the pattern. Defaults to black.

- pattern_fill:

  String. Fill colour of the pattern. Defaults to black.

- pattern_density:

  Numeric. Density of the pattern. Defaults to 0.1

- add_significancies:

  Whether to calculate and plot test results. Defaults to TRUE.

- plot_only_significant_comparisons:

  Logical. Should only significant test results be plotted? Defaults to
  TRUE.

- sig_distance_above_bar:

  Numeric. What distance above the tallest bar the stats results should
  be plotted. Defaults to 1.

- sig_distance_between_significancies:

  Numeric. What distance should be between each test result? Defaults to
  0.8.

- sig_line_thickness:

  Numeric. Thickness of the comparison line. Defaults to 1.

- sig_label_distance_from_line:

  Numeric. Distance of asterisk from the comparison line. Defaults to
  0.2.

- sig_label_text_size:

  Numeric. Size of asterisks. Defaults to 6.

- return_usage_proportions:

  Logical. Whether to return the proportions of usage per group?
  Defaults to FALSE.

## Value

Returns a ggplot object, or a data.frame containing test results.

## Examples

``` r
# TS_plot_IG_usage(input_file = VDJ_data, compare_groups = c("Group1", "Group2", "Group3"), compare_what = "IGHV", export_pdf = TRUE, file_type = "IMGT")
```
