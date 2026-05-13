# Plan collection of selected Cell Ranger files

Builds a copy plan for selected Cell Ranger files without changing the
file system. Matrix files are planned from sample-level filtered matrix
output into per-sample `filtered_feature_bc_matrix` directories with
their original filenames, so they remain directly readable by Seurat.
V(D)J files are planned from run-level unfiltered `outs/vdj_b` and
`outs/vdj_t` directories into per-sample `vdj_b` and `vdj_t` directories
with the sample name prepended.

## Usage

``` r
TS_plan_cellranger_file_collection(
  input_dir,
  dest_dir,
  sample_regex = NULL,
  overwrite = FALSE,
  strict = TRUE,
  ...
)
```

## Arguments

- input_dir:

  Character scalar. Project or Cell Ranger output directory to search.

- dest_dir:

  Character scalar. Destination directory for the collected sample
  folders.

- sample_regex:

  Optional regular expression used to extract sample names. See
  [`TS_find_cellranger_samples()`](https://timothysundell.github.io/TS/reference/TS_find_cellranger_samples.md).

- overwrite:

  Logical. If `FALSE`, the function errors when planned target files
  already exist.

- strict:

  Logical. If `TRUE`, existing raw V(D)J directories must contain all
  requested `all_contig*` files. If `FALSE`, missing V(D)J files are
  recorded in the plan and available files can still be copied.

- ...:

  Additional arguments passed to
  [`TS_find_cellranger_samples()`](https://timothysundell.github.io/TS/reference/TS_find_cellranger_samples.md),
  such as `matrix_dir_names` or `vdj_dir_names`.

## Value

A tibble with one row per planned or missing file. Important columns
include `sample_id`, `file_group`, `source_path`, `dest_path`,
`source_exists`, `copy`, and `status`.

## Examples

``` r
project <- file.path(tempdir(), "cellranger_project_plan")
matrix_dir <- file.path(project, "sampleA", "cellranger_outs", "filtered_feature_bc_matrix")
dir.create(matrix_dir, recursive = TRUE, showWarnings = FALSE)
file.create(file.path(matrix_dir, c("barcodes.tsv.gz", "features.tsv.gz", "matrix.mtx.gz")))
#> [1] TRUE TRUE TRUE

plan <- TS_plan_cellranger_file_collection(
  input_dir = project,
  dest_dir = file.path(tempdir(), "cellranger_transfer_plan")
)
plan
#> # A tibble: 3 × 12
#>   sample_id file_group           source_dir source_filename source_path dest_dir
#>   <chr>     <chr>                <chr>      <chr>           <chr>       <chr>   
#> 1 sampleA   filtered_feature_bc… /tmp/Rtmp… barcodes.tsv.gz /tmp/Rtmp9… /tmp/Rt…
#> 2 sampleA   filtered_feature_bc… /tmp/Rtmp… features.tsv.gz /tmp/Rtmp9… /tmp/Rt…
#> 3 sampleA   filtered_feature_bc… /tmp/Rtmp… matrix.mtx.gz   /tmp/Rtmp9… /tmp/Rt…
#> # ℹ 6 more variables: dest_filename <chr>, dest_path <chr>, required <lgl>,
#> #   source_exists <lgl>, copy <lgl>, status <chr>
```
