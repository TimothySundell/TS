# Compare number of clusters for combinations of number of PCs and resolution

Old function. Not sure of how well it functions, or what it produces. I
think this function counts the number of clusters for each combination
of

- 'dims = 1:i'

- 'resolution = seq(from = 0.1, to = 2, by = 0.1)'

## Usage

``` r
TS_compare_PC_and_resolution(seurat_object, PCs)
```

## Arguments

- seurat_object:

  Your input seurat object

- PCs:

  Number of the maximum PC you want to look at
