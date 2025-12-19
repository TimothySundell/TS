# Plot your Seurat object in 3D.

Will be coloured according to cluster idents with a scales::hue_pal()
colouring scheme. You need to have an already analysed Seurat object to
be able to run the function.

## Usage

``` r
TS_plot_3D(
  seurat_object = get(default_seurat_object),
  dims = 1:20,
  group.by = NULL,
  colors = NULL,
  reduction = "umap"
)
```

## Arguments

- seurat_object:

  The Seurat object you want to analyse. Defaults to
  'default_seurat_object'

- dims:

  What Principal Components to use for calculating the UMAP.

- group.by:

  Variable to group cells by. Default is Seurat::Idents(), alternatives
  are any available metadata in the object.

- colors:

  A vector of colors to use for the group.by parameter. Default is
  hue_pal()(length of unique values in your group.by parameter)

- reduction:

  What reduction to use, defaults to "umap". Alternatively, you can use
  "harmony"
