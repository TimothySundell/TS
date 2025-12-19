# Convert IMGT VDJ data to Cell Ranger–like format

Standardises IMGT VDJ data to match Cell Ranger expectations. Also
checks for the presence of a grouping column and required group values.

## Usage

``` r
TS_format_IMGT_to_cellranger(imgt_data)
```

## Arguments

- imgt_data:

  Data frame from IMGT.

## Value

A standardized data frame ready for use in
[`TS_plot_IG_usage()`](https://timothysundell.github.io/TS/reference/TS_plot_IG_usage.md).
