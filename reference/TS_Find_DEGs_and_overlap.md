# Find the upregulated DEGs for each cluster, comparing two clusters. As well as the conserved markers between them.

Returns a list containing:

- Upregulated markers in cluster 1 compared to cluster 2

- Upregulated markers in cluster 2 compared to cluster 1

- Shared markers between them, compared to all other clusters in the
  Seurat object

Depends on:

- Seurat

- magrittr

- tidyverse

- VennDetail

Automatically filters resulting DEG list so that '0 \< p_val_adj \>
0.05.

## Usage

``` r
TS_Find_DEGs_and_overlap(
  seurat_object,
  ident.1,
  ident.2,
  grouping_variable_name = "tested_population",
  only_pos_markers = T,
  test_use = "wilcox",
  plot_venn = F,
  export_venn = F,
  filename_pdf = "Venn_output",
  width = 8.27,
  heigth = 5.83
)
```

## Arguments

- seurat_object:

  Your Seurat object

- ident.1:

  Name of cluster 1

- ident.2:

  Name of cluster 2

- grouping_variable_name:

  Name of the combined group. Could be anything. Defaults to
  'tested_population'.

- only_pos_markers:

  Whether to only return positively upregulated markers for the specific
  cluster. Defaults to 'TRUE'

- test_use:

  Statistical test used for DE testing. Defaults to 'Wilcoxon ranked sum
  test'.

- plot_venn:

  Whether to plot a Venn diagram for the gene lists. Defaults to 'FALSE'

- export_venn:

  Whether to export the Venn diagram to working directory. Defaults to
  'FALSE'

- filename_pdf:

  File name for exported Venn diagram.

- width:

  Width of saved Venn diagram. Defaults to '8.27 inches', which is A5
  landscape format.

- height:

  Height of saved Venn diagram. Defaults to '5.83 inches', which is A5
  landscape format.
