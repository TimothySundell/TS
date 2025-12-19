# Projects one Seurat object onto another and exports resulting dataframe

A work in progress... Will output resulting dataframe from
[Seurat::MapQuery](https://satijalab.org/seurat/reference/MapQuery.html)

## Usage

``` r
TS_project_scRNAseq_dataset_export_query(
  reference_object,
  query_object,
  dims = 1:20,
  reference_object_assay = "integrated",
  query_object_assay = "integrated",
  normalization.method = "SCT",
  plot_3D = F
)
```

## Arguments

- reference_object:

  A Seurat object, will be used as reference.

- query_object:

  A Seurat object, will be used projected onto your reference object.

- dims:

  Number of Principal Components to use for the UMAP. Defaults to '1:20'

- reference_object_assay:

  Assay to use for reference Seurat object. Defaults to 'integrated'

- query_object_assay:

  Assay to use for query Seurat object. Defaults to 'integrated'

- normalization.method:

  Normalization method used by FindTransferAnchors. Alternatives are
  'LogNormalize' or 'SCT'. Defaults to 'SCT'.

- plot_3D:

  Whether to calculate 3 dimensions for the UMAP model. Defaults to
  'FALSE'
