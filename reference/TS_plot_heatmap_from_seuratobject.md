# Plot a heatmap from a Seurat object

Creates a `ComplexHeatmap` heatmap from a Seurat object using either
cell-level expression values or averaged expression per group. Columns
can be grouped by
[`Seurat::Idents()`](https://satijalab.github.io/seurat-object/reference/Idents.html)
or by any metadata column, and cell-level heatmaps can optionally be
ordered by pseudotime. Only requested features are scaled; a stored
`scale.data` layer is not required.

## Usage

``` r
TS_plot_heatmap_from_seuratobject(
  seurat_object,
  features,
  plot_average_expression = FALSE,
  assay = "RNA",
  layer = "data",
  scale_rows = TRUE,
  group_by = "ident",
  group_levels = NULL,
  pseudotime_var = "pseudotime",
  show_genes = TRUE,
  show_row_dend = FALSE,
  show_row_names = TRUE,
  row_names_fontsize = NULL,
  colors = NULL,
  breaks = NULL,
  color_steps = 11,
  title = NULL,
  border_col = "black",
  split_rows_by = NULL,
  group_bar_fill_colors = NULL,
  rasterise = TRUE,
  raster_by_magick = TRUE,
  legend_name = NULL,
  cluster_rows = FALSE,
  row_clustering_distance = "euclidean",
  row_clustering_method = "complete",
  raster_quality = 2,
  row_names_side = "left",
  row_dend_side = "left",
  row_dend_width = 15,
  highlight_genes = NULL,
  order_pseudotime = FALSE
)
```

## Arguments

- seurat_object:

  A Seurat object containing the expression data to plot.

- features:

  Character vector of genes/features to include in the heatmap.
  Duplicated features are removed, and features not present in the
  object are ignored with a warning.

- plot_average_expression:

  Logical. If `TRUE`, plots average expression per group. If `FALSE`,
  plots cell-level expression.

- assay:

  Character string. Assay to use. Defaults to `"RNA"`.

- layer:

  Character string. Must be `"data"` (the default), containing
  log-normalized expression values.

- scale_rows:

  Logical. Whether to z-score scale each feature across columns before
  plotting. Defaults to `TRUE`. Cell-level values are scaled across
  cells. For average-expression plots, values are first averaged in
  linear space within each group, transformed with
  [`log1p()`](https://rdrr.io/r/base/Log.html), and then scaled across
  groups, giving each group equal weight. Scaling uses
  [`Seurat::ScaleData()`](https://satijalab.org/seurat/reference/ScaleData.html)
  with its default upper cap of 10. Constant features (including a
  single-group average) become zero. Set to `FALSE` to plot
  log-normalized values without scaling.

- group_by:

  Character string. Grouping variable for columns. Use `"ident"` to
  group by
  [`Seurat::Idents()`](https://satijalab.github.io/seurat-object/reference/Idents.html),
  or provide a metadata column name.

- group_levels:

  Optional character vector specifying the order of groups in the
  heatmap. If `NULL`, the existing factor levels are used when
  available.

- pseudotime_var:

  Character string. Metadata column used when `order_pseudotime = TRUE`.
  Defaults to `"pseudotime"`.

- show_genes:

  Logical. Whether gene names should be shown when row labels are
  enabled.

- show_row_dend:

  Logical. Whether to show the row dendrogram.

- show_row_names:

  Logical. Whether to show row names. If `split_rows_by` is supplied,
  row names are instead drawn in the row annotation.

- row_names_fontsize:

  Numeric. Font size for row labels.

- colors:

  Optional color specification for the heatmap scale. Can be either a
  character vector of colors or a palette function such as
  `scales::brewer_pal(palette = "RdBu")`. If `NULL`, defaults are chosen
  based on whether the plotted matrix is scaled or unscaled.

- breaks:

  Optional numeric vector of breakpoints for the heatmap color scale. If
  `NULL`, sensible defaults are chosen from the plotted matrix. If the
  number of supplied colors differs from the number of breaks, the
  colors are distributed evenly across the range of `breaks`.

- color_steps:

  Integer. Number of colors to generate when `colors` is a palette
  function. Defaults to `11`.

- title:

  Optional overall heatmap title. Defaults to `NULL`, meaning no title
  is drawn.

- border_col:

  Border color for the heatmap body. Defaults to `"black"`.

- split_rows_by:

  Optional vector used to split rows into groups. Must have the same
  length as `features` before missing features are removed.

- group_bar_fill_colors:

  Optional vector of colors for the column group annotation. Can be
  unnamed, or a named vector matching the group labels.

- rasterise:

  Logical. Whether to rasterize the heatmap body.

- raster_by_magick:

  Logical. Passed to
  [`ComplexHeatmap::Heatmap()`](https://rdrr.io/pkg/ComplexHeatmap/man/Heatmap.html).

- legend_name:

  Optional legend title. If `NULL`, a default is chosen based on plot
  mode.

- cluster_rows:

  Logical. Whether to cluster rows hierarchically.

- row_clustering_distance:

  Character string specifying the row distance metric. Defaults to
  `"euclidean"`. Also supports `"pearson"` and `"spearman"`.

- row_clustering_method:

  Character string specifying the agglomeration method passed to
  [`fastcluster::hclust()`](https://rdrr.io/pkg/fastcluster/man/hclust.html).
  Defaults to `"complete"`.

- raster_quality:

  Numeric. Raster quality passed to
  [`ComplexHeatmap::Heatmap()`](https://rdrr.io/pkg/ComplexHeatmap/man/Heatmap.html).

- row_names_side:

  Character string. Side for row names. Defaults to `"left"`.

- row_dend_side:

  Character string. Side for the row dendrogram. Defaults to `"left"`.

- row_dend_width:

  Numeric. Width of the row dendrogram in mm.

- highlight_genes:

  Optional character vector of genes to highlight with labels.

- order_pseudotime:

  Logical. If `TRUE`, plots cell-level expression with columns ordered
  by `pseudotime_var`.

## Value

A `ComplexHeatmap` heatmap object. If `title` is supplied, the returned
object is the result of
[`ComplexHeatmap::draw()`](https://rdrr.io/pkg/ComplexHeatmap/man/draw-dispatch.html).

## Examples

``` r
# Average expression grouped by Seurat identities
# TS_plot_heatmap_from_seuratobject(
#   seurat_object = seurat_obj,
#   features = c("CD79A", "MS4A1", "CD74"),
#   plot_average_expression = TRUE
# )

# Cell-level heatmap grouped by a metadata column with custom order
# TS_plot_heatmap_from_seuratobject(
#   seurat_object = seurat_obj,
#   features = c("MKI67", "TOP2A", "PCNA"),
#   group_by = "Phase",
#   group_levels = c("G1", "S", "G2M")
# )

# Order cells by pseudotime
# TS_plot_heatmap_from_seuratobject(
#   seurat_object = seurat_obj,
#   features = c("LEF1", "RAG1", "RAG2"),
#   order_pseudotime = TRUE
# )

# Use a richer Brewer-style diverging palette
# TS_plot_heatmap_from_seuratobject(
#   seurat_object = seurat_obj,
#   features = c("LEF1", "RAG1", "RAG2"),
#   colors = scales::brewer_pal(palette = "RdBu")(11),
#   breaks = c(-1, 0, 1)
# )
```
