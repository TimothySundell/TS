# Wrapper around the enrichR package.

Performs enrichment analyses for a vector of genes Depends on -enrichR
(and that their servers are up, heh)

- A live Internet connection

- magrittr

- dplyr Returns a list of the results of each database tested for.

## Usage

``` r
TS_run_enrichR(
  gene_list,
  plot_enrichr = F,
  database_to_plot = "GO_Biological_Process_2023",
  databases = c("GO_Biological_Process_2023", "GO_Molecular_Function_2023",
    "GO_Cellular_Component_2023")
)
```

## Arguments

- gene_list:

  Character vector containing genes to test

- plot_enrichr:

  Whether to plot the output using the enrichR function
  [plotEnrich](https://rdrr.io/pkg/enrichR/man/plotEnrich.html).
  Defaults to 'FALSE'

- database_to_plot:

  What database to plot the output from. Defaults to
  '"GO_Biological_Process_2023"'

- databases:

  What databases to test against. Defaults to
  'c("GO_Biological_Process_2023", "GO_Molecular_Function_2023",
  "GO_Cellular_Component_2023")'
