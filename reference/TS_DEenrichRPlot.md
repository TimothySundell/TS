# A wrapper around the DEnrichRPlot from the 'Seurat' package

A wrapper around the DEnrichRPlot from the 'Seurat' package

## Usage

``` r
TS_DEenrichRPlot(
  seurat_object = get(default_seurat_object),
  ident.1 = 0,
  ident.2 = 1,
  test.use = "wilcox",
  max.genes = 500
)
```

## Arguments

- seurat_object:

  Your 'Seurat' object. Defaults to 'default_seurat_object'

- ident.1:

  Name of cluster 1. Defaults to '0'

- ident.2:

  Name of cluster 2. Defaults to '1'

- test.use:

  The test you want to use for DE testing. Defaults to 'Wilcoxons ranked
  sum' test.

- max.genes:

  The maximum number of genes to use. Defaults to 500
