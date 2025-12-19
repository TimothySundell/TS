# Plots Ig usage between two groups using Fisher's Exact Test

Plot the Ig usage from two groups and perform Fisher's Exact Test to
find significant differences. Maximum number of groups = 2

## Usage

``` r
TS_plot_IG_usage(
  input_file,
  sample1,
  sample2,
  compare_what = "IGHV",
  grouping_variable = "named_clusters",
  export_pdf = F,
  return_only_data = F,
  ylim_max = NULL,
  p.adj.method = "hochberg",
  fill_colours = c("gray80", "gray50"),
  border_colour = "black",
  ggplot_theme = ggplot2::theme_classic(),
  plot_width = 15,
  plot_height = 5,
  remove_genes_present_in_one_dataset = F,
  plot_zeroes = T,
  file_type = "cellranger"
)
```

## Arguments

- input_file:

  A dataframe (e.g. from 10X Genomics Cellranger VDJ output) that
  includes a grouping column.

- sample1:

  Grouping variable value for sample 1

- sample2:

  Grouping variable value for sample 2

- compare_what:

  Genes to compare (e.g. "IGHV", "IGKV"). Defaults to "IGHV".

- grouping_variable:

  Column used to group data. Defaults to "named_clusters"

- export_pdf:

  Logical. If TRUE, exports the plot to a PDF. Defaults to FALSE

- return_only_data:

  Logical. If TRUE, returns only the calculated frequencies and the test
  matrices. Defaults to FALSE.

- ylim_max:

  Numeric. Optional upper limit for y-axis. ggplot2 will adjust the
  y-axis scale to fit the data.

- p.adj.method:

  Method for adjust for multiple tests. Defaults to 'Hochberg'.
  Alternative is 'FDR'.

- fill_colours:

  A character vector of fill colours for the two groups. Defaults to
  c("grey80", "grey50").

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

- file_type:

  Set type of file/annotation input. Defaults to "cellranger". Options
  are "IMGT" and "AIRR".

## Value

Returns a ggplot object, or a data.frame containing counts and test
results if ' return_only_data = TRUE'.

## Examples

``` r
# TS_plot_IG_usage(input_file = VDJ_data, sample1 = "group_A", sample2 = "group_B", compare_what = "IGHV")
```
