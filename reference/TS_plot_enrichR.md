# Plot the output from [TS_run_enrichR](https://timothysundell.github.io/TS/reference/TS_run_enrichR.md)

Plot the output from
[TS_run_enrichR](https://timothysundell.github.io/TS/reference/TS_run_enrichR.md)

## Usage

``` r
TS_plot_enrichR(input, sort_on = "neg_log10", n = 10)
```

## Arguments

- input:

  The output list from
  [TS_run_enrichR](https://timothysundell.github.io/TS/reference/TS_run_enrichR.md)

- sort_on:

  What to sort top n terms on. Defaults to 'neg_log10', can also be
  'combined_score'

- n:

  Number of terms to plot for each database. Defaults to '10'.
