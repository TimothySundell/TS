# Projects one Seurat object onto another

A work in progress... Will output resulting DimPlots with data from
[Seurat::MapQuery](https://satijalab.org/seurat/reference/MapQuery.html)

## Usage

``` r
TS_project_scRNAseq_dataset(
  reference_object,
  query_object,
  ref_name,
  query_name
)
```

## Arguments

- reference_object:

  A Seurat object, will be used as reference.

- query_object:

  A Seurat object, will be used projected onto your reference object.

- ref_name:

  Name of reference object. Used for resulting plots.

- query_name:

  Name of query object. Used for resulting plots.
