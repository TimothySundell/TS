# Plots IGHV usage between two groups using Fisher's Exact Test

Plot the IGHV usage from two groups and perform Fisher's Exact Test to
find significant differences. Limited to comparing two groups.

## Usage

``` r
TS_IGHV_usage_prop_with_stats(
  input_file,
  sample1,
  sample2,
  export_pdf = F,
  return_only_data = F,
  ylim_max = NULL,
  p.adj.method = "hochberg",
  fill_colours = c("gray80", "gray50"),
  border_colour = "black"
)
```

## Arguments

- input_file:

  Takes as input a dataframe from 10X Genomixs Cellranger VDJ pipeline
  that has been through QC with added sample/grouping information in the
  column "sample"

- sample1:

  Grouping variable value for sample 1

- sample2:

  Grouping variable value for sample 2

- export_pdf:

  Whether to automatically export the plot as PDF in your working
  directory. Defaults to 'FALSE'

- return_only_data:

  Whether to only return the the calculated frequencies and the test
  matrices. 'FALSE' returns the plot in standard out. Defaults to
  'FALSE'

- ylim_max:

  Defaults to NULL. ggplot2 will adjust the y-axis scale to fit the
  data. Set this if you want to keep it consistent between plots.

- p.adj.method:

  Method for adjust for multiple tests. Defaults to
  'Benjamini-Hochberg'. Alternative is 'FDR'.
