# Output one or multiple featureplots

Output one or multiple featureplots

## Usage

``` r
TS_featureplots(
  features,
  seurat_object = get(default_seurat_object),
  wrap_plots = TRUE,
  wrap_n_columns = NULL,
  add_dimplot = TRUE,
  dimplot_group.by = NULL,
  dimplot_cols = NULL,
  dimplot_label = TRUE,
  expression_cols = c("lightgray", "#FEFDE2", "#FEEAC7", "#FDD7AC", "#FDC591", "#FDB276",
    "#FC9F5B", "#FC8C40", "#FC8C40", "#EB753C", "#DA5D38", "#C94634", "#B72F30",
    "#A6172C", "#950028"),
  reduction = "umap",
  pt.size = 1,
  order = TRUE,
  assay = NULL
)
```

## Arguments

- features:

  Vector of features to plot.

- seurat_object:

  The Seurat object you wish to plot. Defaults to
  "default_seurat_object".

- wrap_plots:

  Whether to wrap plots into one by patchwork::wrap_plots(). Defaults to
  TRUE.

- wrap_n_columns:

  If wrap_plots = T, how many columns do you want for the final plot.

- add_dimplot:

  If wrap_plots = T, whether to include Seurat::DimPlot() as the first
  panel. Defaults to TRUE

- dimplot_group.by:

  What metadata to group cells by. Defaults to
  Seurat::Idents(seurat_object)

- dimplot_cols:

  Colours to use for the Seurat::DimPlot(). Defaults to
  scales::hue_pal()(length(levels(Seurat::Idents(seurat_object))))

- dimplot_label:

  Whether to label the DimPlot by dimplot_group.by. Defaults to TRUE.

- expression_cols:

  Vector of colors to use (low-\>high), defaults to a 15-step gradient
  from "lightgray" to "#950028"

- reduction:

  Dimensionality reduction technique to use for plotting, default =
  "umap".

- pt.size:

  Size of each dot, default = 1.

- order:

  Whether to plot positive cells last, e.g. on top, default = TRUE.

- assay:

  What assay to pull gene expression data from. Defaults to "RNA", as
  this is the convention.
