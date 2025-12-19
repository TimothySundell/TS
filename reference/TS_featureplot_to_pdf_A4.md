# Export featureplot to pdf in A4 format

Export featureplot to pdf in A4 format

## Usage

``` r
TS_featureplot_to_pdf_A4(
  seurat_object,
  file_prefix,
  features,
  pt.size = 1.5,
  reduction = "umap",
  order = TRUE,
  cols = c("lightgray", "#FEFDE2", "#FEEAC7", "#FDD7AC", "#FDC591", "#FDB276", "#FC9F5B",
    "#FC8C40", "#FC8C40", "#EB753C", "#DA5D38", "#C94634", "#B72F30", "#A6172C",
    "#950028")
)
```

## Arguments

- seurat_object:

  The Seurat object you wish to plot

- file_prefix:

  Prefix for the saved file, e.g. sample name

- features:

  Vector of features to plot

- pt.size:

  Dot size, default = 1.5

- reduction:

  Dimensionality reduction technique to use for plotting, default =
  "umap"

- order:

  Whether to plot positive cells last, e.g. on top, default = 1.5

- cols:

  Vector of 2 colors to use (low-\>high), default = c("lightgray",
  "red")
