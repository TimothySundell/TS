# Create column chart with ggplot2

Create column chart with ggplot2

## Usage

``` r
TS_plot_column_chart(
  data,
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
  plot_n_sequences
)
```

## Arguments

- x:

  Variable to plot on x-axis.

- y:

  Variable to plot on y-axis.

- fill:

  Variable to fill on (group variable).

- fill_colours:

  Vector. Fill colours to use.

- name:

  String. Legend title. Defaults to "Sample"

- ggplot_theme:

  Optional. A specified ggplot_theme to use.

- plot_data:

  A dataframe containing the long format data which you want to plot
  with precomputed proportions.

- ylim:

  Optional, numeric. Use to specify max y_lim value.

## Value

A tibble with columns: gene, group1, group2, p_value, p_adj
