# Find conserved markers within two groups

Not sure this function does what I want it to...

## Usage

``` r
TS_FindConservedMarkers(
  seurat_object,
  ident.1,
  ident.2,
  grouping_variable_name = "tested_population",
  only_pos_markers = T
)
```

## Arguments

- seurat_object:

  Your Seurat object

- ident.1:

  Name for cluster 1

- ident.2:

  Name for cluster 2

- grouping_variable_name:

  Name of the combined clusters, can be anything. Defaults to
  'tested_population'.

- only_pos_markers:

  Whether to only highlight positively upregulated markers between the
  grouped population and the rest. Defaults to 'TRUE'.
