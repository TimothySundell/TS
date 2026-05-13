# Find sample-level Cell Ranger output folders

Recursively searches a project folder for Cell Ranger matrix directories
and infers the corresponding sample names. The function is designed for
common `cellranger multi` layouts such as
`outs/per_sample_outs/<sample>` and for project folders organised as
`project/<sample>/cellranger_outs`.

The filtered gene expression matrix is detected at sample level. V(D)J
directories are detected at run level, for example `outs/vdj_b` and
`outs/vdj_t`, because `cellranger multi` stores the unfiltered
`all_contig*` files there rather than in `per_sample_outs/<sample>`.

## Usage

``` r
TS_find_cellranger_samples(
  input_dir,
  sample_regex = NULL,
  matrix_dir_names = c("filtered_feature_bc_matrix", "sample_filtered_feature_bc_matrix"),
  vdj_dir_names = c("vdj_b", "vdj_t")
)
```

## Arguments

- input_dir:

  Character scalar. Project or Cell Ranger output directory to search.

- sample_regex:

  Optional regular expression used to extract sample names from matrix
  directory paths. If the regular expression contains a capture group,
  the first captured group is used; otherwise the full match is used.

- matrix_dir_names:

  Character vector of matrix directory names to search for. Defaults to
  Cell Ranger's `filtered_feature_bc_matrix` and
  `sample_filtered_feature_bc_matrix`.

- vdj_dir_names:

  Character vector of V(D)J directory names to detect.

## Value

A tibble with one row per detected sample matrix and columns
`sample_id`, `sample_root`, `outs_dir`, `matrix_dir`, `vdj_b_dir`,
`vdj_t_dir`, and `detection_note`.

## Examples

``` r
project <- file.path(tempdir(), "cellranger_project")
matrix_dir <- file.path(project, "sampleA", "cellranger_outs", "filtered_feature_bc_matrix")
dir.create(matrix_dir, recursive = TRUE, showWarnings = FALSE)
file.create(file.path(matrix_dir, c("barcodes.tsv.gz", "features.tsv.gz", "matrix.mtx.gz")))
#> [1] TRUE TRUE TRUE

TS_find_cellranger_samples(project)
#> # A tibble: 1 × 7
#>   sample_id sample_root   outs_dir matrix_dir vdj_b_dir vdj_t_dir detection_note
#>   <chr>     <chr>         <chr>    <chr>      <chr>     <chr>     <chr>         
#> 1 sampleA   /tmp/RtmpBku… /tmp/Rt… /tmp/Rtmp… NA        NA        cellranger_ou…
```
