# Collect selected Cell Ranger files into per-sample folders

Copies selected Cell Ranger files into a transfer-friendly destination
directory. Each sample receives its own folder. Filtered matrix
filenames are kept unchanged for Seurat compatibility, and unfiltered
V(D)J `all_contig*` files are renamed with the sample ID as a prefix.
Per-sample V(D)J folders are intentionally not used because they contain
filtered contig outputs.

## Usage

``` r
TS_collect_cellranger_files(
  input_dir,
  dest_dir,
  sample_regex = NULL,
  overwrite = FALSE,
  confirm = interactive(),
  execute = TRUE,
  strict = TRUE,
  ...
)
```

## Arguments

- input_dir:

  Character scalar. Project or Cell Ranger output directory to search.

- dest_dir:

  Character scalar. Destination directory for collected files.

- sample_regex:

  Optional regular expression used to extract sample names. See
  [`TS_find_cellranger_samples()`](https://timothysundell.github.io/TS/reference/TS_find_cellranger_samples.md).

- overwrite:

  Logical. If `FALSE`, existing destination files cause an error before
  copying starts.

- confirm:

  Logical. If `TRUE`, ask for confirmation before copying.

- execute:

  Logical. If `FALSE`, return the planned copy operations without
  creating directories or copying files.

- strict:

  Logical. If `TRUE`, existing raw V(D)J directories must contain all
  requested `all_contig*` files. If `FALSE`, missing V(D)J files are
  skipped and recorded in the returned plan.

- ...:

  Additional arguments passed to
  [`TS_find_cellranger_samples()`](https://timothysundell.github.io/TS/reference/TS_find_cellranger_samples.md),
  such as `matrix_dir_names` or `vdj_dir_names`.

## Value

Invisibly returns the copy plan tibble.

## Examples

``` r
project <- file.path(tempdir(), "cellranger_project_collect")
matrix_dir <- file.path(project, "sampleA", "cellranger_outs", "filtered_feature_bc_matrix")
dir.create(matrix_dir, recursive = TRUE, showWarnings = FALSE)
file.create(file.path(matrix_dir, c("barcodes.tsv.gz", "features.tsv.gz", "matrix.mtx.gz")))
#> [1] TRUE TRUE TRUE

TS_collect_cellranger_files(
  input_dir = project,
  dest_dir = file.path(tempdir(), "cellranger_transfer_collect"),
  execute = FALSE
)
#> Proposed Cell Ranger file collection:
#> # A tibble: 3 × 5
#>   sample_id file_group                 source_path              dest_path status
#>   <chr>     <chr>                      <chr>                    <chr>     <chr> 
#> 1 sampleA   filtered_feature_bc_matrix /tmp/Rtmp9VnWrG/cellran… /tmp/Rtm… ready 
#> 2 sampleA   filtered_feature_bc_matrix /tmp/Rtmp9VnWrG/cellran… /tmp/Rtm… ready 
#> 3 sampleA   filtered_feature_bc_matrix /tmp/Rtmp9VnWrG/cellran… /tmp/Rtm… ready 
#> Dry run only. No directories were created and no files were copied.
```
