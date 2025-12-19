# Calculate light chain ratios of a Seurat object, requires light chain data in metadata "light.chain"

Calculate light chain ratios of a Seurat object, requires light chain
data in metadata "light.chain"

## Usage

``` r
TS_calculate_ratios(
  seurat_object = get(default_seurat_object),
  group_var = "ident"
)
```

## Arguments

- seurat_object:

  The Seurat object you want to analyse

- group_var:

  What parameter to group by, has to exist as a metadata slot in the
  Seurat object. Defaults to "ident"
