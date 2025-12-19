# Plot counts for a gene in a Seurat object as violinplot

Plot counts for a gene in a Seurat object as violinplot

## Usage

``` r
TS_plot_counts_violin(
  seurat_object = get(default_seurat_object),
  slot = "counts",
  gene
)
```

## Arguments

- seurat_object:

  The Seurat object you want to analyse

- slot:

  Slot to plot. Defaults to "counts"

- gene:

  Gene to plot
