# Process CellRanger output files with Seurat

Input the three matrices from CellRanger and analyse them automatically.
Will first run an analysis pipeline with log-normalized counts, allow
you to set parameters on-the-fly, remove unwanted clusters, and then run
SCTransform_v1 on the filtered data.

Based on Seurat, using the workflow published in Sundell et. al., 2022
(https://doi.org/10.1093/bfgp/elac044)

## Usage

``` r
TS_Seurat_from_file(
  input_data,
  sample_ID = NULL,
  remove_ig_genes = T,
  regular_clustering_dims = NULL,
  SCT_clustering_dims = NULL,
  additional_markers = NULL,
  project_name = NULL,
  min.cells = 3,
  min.features = 200,
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  selection.method = "vst",
  nfeatures = 2000,
  nn.method = "rann"
)
```

## Arguments

- input_data:

  File path to the folder containing the three matrices used for Seurat

- sample_ID:

  Required ID-tag for your sample. Will be added as metadata in column
  "sample_ID". Handy in downstream analyses, e.g. integration,
  batch-correction etc.

- remove_ig_genes:

  Filters out Ig-genes prior to clustering. Which can help in finding
  biologically relevant clusters, as shown in Sundell et. al., 2022.
  Defaults to TRUE.

- regular_clustering_dims:

  What PCA components to use for building NN graph and dimensionality
  reduction with regular normalization/scaling. Defaults to NULL.

- SCT_clustering_dims:

  What PCA components to use for building NN graph and dimensionality
  reduction with SCTransformed data. Defaults to NULL.

- additional_markers:

  Additional gene names to plot to aid in exclusion of unwanted cell
  types. Default markers are "percent.mt", "percent.ribo",
  "percent.malat1", "CD19", "CD3E", "CD14", "CD56".

- project_name:

  Required by Seurat. Defaults to the same value as 'sample_ID'. Added
  automatically as metadata in column 'orig.ident'

- min.cells:

  Include features/genes detected in at least this many cells. Defaults
  to 3.

- min.features:

  Keep cells with at least this many features detected. Defaults to 200.

- normalization.method:

  Normalization method to use. Alternatives are "LogNormalize"
  (default), "CLR", "RC".

- scale.factor:

  Scale factor for cell-level normalization. Defaults to 10000

- selection.method:

  Method for finding variable features. Alternatives are "vst"
  (default), "mean.var.plot", "dispersion".

- nn.method:

  Nearest-neighbour method to use for clustering with regular
  normalization/scaling. Alternatives are "rann" or "annoy". Defaults to
  "rann".

- nFeatures:

  Number of variable features to calculate. Defaults to 2000
