# Convert IMGT VDJ data to Cell Ranger–like format

Standardizes IMGT VDJ data to match Cell Ranger expectations. Also
checks for the presence of a grouping column and required group values.

## Usage

``` r
TS_format_IMGT_to_cellranger(imgt_data, grouping_variable, sample1, sample2)
```

## Arguments

- imgt_data:

  Data frame from IMGT.

- grouping_variable:

  Name of the grouping column (character).

- sample1:

  Grouping variable value for sample 1.

- sample2:

  Grouping variable value for sample 2.

## Value

A standardized data frame ready for use in
[`TS_plot_IG_usage()`](https://timothysundell.github.io/TS/reference/TS_plot_IG_usage.md).
