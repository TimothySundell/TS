#' Plot a heatmap from a Seurat object
#'
#' @description
#' Creates a `ComplexHeatmap` heatmap from a Seurat object using either
#' cell-level expression values or averaged expression per group.
#' Columns can be grouped by `Seurat::Idents()` or by any metadata column,
#' and cell-level heatmaps can optionally be ordered by pseudotime.
#'
#' @param seurat_object A Seurat object containing the expression data to plot.
#' @param features Character vector of genes/features to include in the heatmap.
#' Duplicated features are removed, and features not present in the object are ignored with a warning.
#' @param plot_average_expression Logical. If `TRUE`, plots average expression per group.
#' If `FALSE`, plots cell-level expression.
#' @param assay Character string. Assay to use. Defaults to `"RNA"`.
#' @param layer Character string. Layer to pull expression values from.
#' This is passed on to Seurat functions directly.
#' @param scale_rows Logical. Whether to z-score scale each feature across columns
#' before plotting. Defaults to `FALSE`. If you want pre-scaled values, you can
#' instead use `layer = "scale.data"`.
#' @param group_by Character string. Grouping variable for columns. Use `"ident"`
#' to group by `Seurat::Idents()`, or provide a metadata column name.
#' @param group_levels Optional character vector specifying the order of groups in the heatmap.
#' If `NULL`, the existing factor levels are used when available.
#' @param pseudotime_var Character string. Metadata column used when `order_pseudotime = TRUE`.
#' Defaults to `"pseudotime"`.
#' @param show_genes Logical. Whether gene names should be shown when row labels are enabled.
#' @param show_row_dend Logical. Whether to show the row dendrogram.
#' @param show_row_names Logical. Whether to show row names.
#' If `split_rows_by` is supplied, row names are instead drawn in the row annotation.
#' @param row_names_fontsize Numeric. Font size for row labels.
#' @param colors Optional color specification for the heatmap scale. Can be either
#' a character vector of colors or a palette function such as
#' `scales::brewer_pal(palette = "RdBu")`. If `NULL`, defaults are chosen based on
#' whether the plotted matrix is scaled or unscaled.
#' @param breaks Optional numeric vector of breakpoints for the heatmap color scale.
#' If `NULL`, sensible defaults are chosen from the plotted matrix. If the number of
#' supplied colors differs from the number of breaks, the colors are distributed
#' evenly across the range of `breaks`.
#' @param color_steps Integer. Number of colors to generate when `colors` is a palette
#' function. Defaults to `11`.
#' @param title Optional overall heatmap title. Defaults to `NULL`, meaning no title is drawn.
#' @param border_col Border color for the heatmap body. Defaults to `"black"`.
#' @param split_rows_by Optional vector used to split rows into groups.
#' Must have the same length as `features` before missing features are removed.
#' @param group_bar_fill_colors Optional vector of colors for the column group annotation.
#' Can be unnamed, or a named vector matching the group labels.
#' @param rasterise Logical. Whether to rasterize the heatmap body.
#' @param raster_by_magick Logical. Passed to `ComplexHeatmap::Heatmap()`.
#' @param legend_name Optional legend title. If `NULL`, a default is chosen based on plot mode.
#' @param cluster_rows Logical. Whether to cluster rows hierarchically.
#' @param row_clustering_distance Character string specifying the row distance metric.
#' Defaults to `"euclidean"`. Also supports `"pearson"` and `"spearman"`.
#' @param row_clustering_method Character string specifying the agglomeration method
#' passed to `fastcluster::hclust()`. Defaults to `"complete"`.
#' @param raster_quality Numeric. Raster quality passed to `ComplexHeatmap::Heatmap()`.
#' @param row_names_side Character string. Side for row names. Defaults to `"left"`.
#' @param row_dend_side Character string. Side for the row dendrogram. Defaults to `"left"`.
#' @param row_dend_width Numeric. Width of the row dendrogram in mm.
#' @param highlight_genes Optional character vector of genes to highlight with labels.
#' @param order_pseudotime Logical. If `TRUE`, plots cell-level expression with columns ordered by `pseudotime_var`.
#'
#' @return A `ComplexHeatmap` heatmap object. If `title` is supplied, the returned
#' object is the result of `ComplexHeatmap::draw()`.
#'
#' @export
#'
#' @examples
#' # Average expression grouped by Seurat identities
#' # TS_plot_heatmap_from_seuratobject(
#' #   seurat_object = seurat_obj,
#' #   features = c("CD79A", "MS4A1", "CD74"),
#' #   plot_average_expression = TRUE
#' # )
#'
#' # Cell-level heatmap grouped by a metadata column with custom order
#' # TS_plot_heatmap_from_seuratobject(
#' #   seurat_object = seurat_obj,
#' #   features = c("MKI67", "TOP2A", "PCNA"),
#' #   group_by = "Phase",
#' #   group_levels = c("G1", "S", "G2M")
#' # )
#'
#' # Order cells by pseudotime
#' # TS_plot_heatmap_from_seuratobject(
#' #   seurat_object = seurat_obj,
#' #   features = c("LEF1", "RAG1", "RAG2"),
#' #   order_pseudotime = TRUE
#' # )
#'
#' # Use a richer Brewer-style diverging palette
#' # TS_plot_heatmap_from_seuratobject(
#' #   seurat_object = seurat_obj,
#' #   features = c("LEF1", "RAG1", "RAG2"),
#' #   colors = scales::brewer_pal(palette = "RdBu")(11),
#' #   breaks = c(-1, 0, 1)
#' # )
TS_plot_heatmap_from_seuratobject = function(
  seurat_object,
  features,
  plot_average_expression = FALSE,
  assay = "RNA",
  layer = "scale.data",
  scale_rows = FALSE,
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
) {
  library(ComplexHeatmap)
  library(scales)
  library(circlize)
  library(dplyr)
  library(Seurat)
  library(grid)

  original_features <- features

  if (!assay %in% Seurat::Assays(seurat_object)) {
    stop(
      "`assay = ", assay,
      "` was not found in the Seurat object. Available assays: ",
      paste(Seurat::Assays(seurat_object), collapse = ", ")
    )
  }

  if (!is.null(split_rows_by) && length(original_features) != length(split_rows_by)) {
    stop(
      "Length of features and split_rows_by are not the same. Most likely due to duplicated genes in your list. Filter input data and try again"
    )
  }

  if (anyDuplicated(original_features)) {
    keep_idx <- !duplicated(original_features)
    features <- original_features[keep_idx]
    if (!is.null(split_rows_by)) {
      split_rows_by <- split_rows_by[keep_idx]
    }
  } else {
    features <- original_features
  }

  assay_data <- Seurat::GetAssayData(
    object = seurat_object,
    assay = assay,
    layer = layer
  )

  features_in_object <- rownames(assay_data)
  features_present <- features[features %in% features_in_object]
  features_not_present <- base::setdiff(features, features_present)

  if (length(features_not_present) > 0) {
    cli::cli_alert_warning(
      "These features were not present: {features_not_present}"
    )
  }

  if (length(features_present) < 1) {
    stop("None of your supplied features are present in the seurat object")
  }

  if (!is.null(split_rows_by)) {
    split_rows_by <- split_rows_by[features %in% features_present]
    names(split_rows_by) <- features_present
  }

  if (isTRUE(order_pseudotime)) {
    plot_mode <- "pseudotime"
  } else if (isTRUE(plot_average_expression)) {
    plot_mode <- "average"
  } else if (identical(plot_average_expression, FALSE)) {
    plot_mode <- "cells"
  } else {
    stop(
      "You tried to run `plot_average_expression = ",
      plot_average_expression,
      "`. This argument only accepts TRUE or FALSE."
    )
  }

  resolve_heatmap_colors <- function(colors, breaks, color_steps) {
    if (is.function(colors)) {
      colors <- colors(color_steps)
    }

    if (!is.character(colors) || length(colors) < 2) {
      stop(
        "`colors` must be either a character vector with at least two colors or a palette function."
      )
    }

    if (!is.numeric(breaks) || length(breaks) < 2) {
      stop("`breaks` must be a numeric vector with at least two values.")
    }

    if (length(colors) != length(breaks)) {
      breaks <- seq(min(breaks), max(breaks), length.out = length(colors))
    }

    list(colors = colors, breaks = breaks)
  }

  get_default_group_levels <- function(group_by) {
    if (identical(group_by, "ident")) {
      return(levels(Seurat::Idents(seurat_object)))
    }

    meta_data <- seurat_object[[]]

    if (!group_by %in% colnames(meta_data)) {
      stop(
        "`group_by = ", group_by,
        "` was not found in the Seurat metadata. Use `ident` or a valid metadata column name."
      )
    }

    group_values <- meta_data[[group_by]]

    if (is.factor(group_values)) {
      return(levels(group_values))
    }

    unique(as.character(group_values))
  }

  grouping_levels <- if (is.null(group_levels)) {
    get_default_group_levels(group_by)
  } else {
    unique(as.character(group_levels))
  }

  make_grouping_factor <- function(x, levels = grouping_levels) {
    x <- as.character(x)
    effective_levels <- c(levels, setdiff(unique(x), levels))
    factor(x, levels = effective_levels)
  }

  make_group_colors <- function(group_levels, provided_colors = NULL) {
    if (is.null(provided_colors)) {
      group_colors <- scales::hue_pal()(length(group_levels))
      names(group_colors) <- group_levels
      return(group_colors)
    }

    if (!is.null(names(provided_colors))) {
      group_colors <- provided_colors[group_levels]
      missing_idx <- is.na(group_colors)

      if (any(missing_idx)) {
        fallback_colors <- scales::hue_pal()(sum(missing_idx))
        group_colors[missing_idx] <- fallback_colors
      }

      names(group_colors) <- group_levels
      return(group_colors)
    }

    if (length(provided_colors) < length(group_levels)) {
      stop(
        "`group_bar_fill_colors` must be NULL, a named vector, or have at least as many colors as grouping levels."
      )
    }

    group_colors <- provided_colors[seq_along(group_levels)]
    names(group_colors) <- group_levels
    group_colors
  }

  make_row_split_annotation <- function(expression_mat, split_rows_by, row_names_fontsize) {
    if (is.null(split_rows_by)) {
      return(NULL)
    }

    ComplexHeatmap::rowAnnotation(
      a = ComplexHeatmap::anno_block(
        gp = grid::gpar(fill = "black"),
        width = grid::unit(0.5, "mm")
      ),
      b = ComplexHeatmap::anno_text(
        rownames(expression_mat),
        location = 1,
        just = "right",
        gp = grid::gpar(fontsize = row_names_fontsize)
      )
    )
  }

  make_gene_highlight_annotation <- function(expression_mat, highlight_genes) {
    if (is.null(highlight_genes)) {
      return(NULL)
    }

    highlight_genes <- highlight_genes[highlight_genes %in% rownames(expression_mat)]

    if (length(highlight_genes) < 1) {
      return(NULL)
    }

    ComplexHeatmap::rowAnnotation(
      foo = ComplexHeatmap::anno_mark(
        at = match(highlight_genes, rownames(expression_mat)),
        labels = highlight_genes
      )
    )
  }

  compute_row_distance <- function(expression_mat, distance_method) {
    if (distance_method %in% c("pearson", "spearman")) {
      correlation_mat <- stats::cor(
        t(expression_mat),
        method = distance_method,
        use = "pairwise.complete.obs"
      )
      correlation_mat[is.na(correlation_mat)] <- 0
      diag(correlation_mat) <- 1
      return(stats::as.dist(1 - correlation_mat))
    }

    stats::dist(expression_mat, method = distance_method)
  }

  row_clustering_fun <- if (isTRUE(cluster_rows)) {
    if (!requireNamespace("fastcluster", quietly = TRUE)) {
      stop(
        "Package `fastcluster` is required when `cluster_rows = TRUE`. Please install it first."
      )
    }

    if (
      identical(row_clustering_method, "ward.D2") &&
      row_clustering_distance %in% c("pearson", "spearman")
    ) {
      warning(
        "`ward.D2` is usually intended for Euclidean distances. ",
        "For correlation-based clustering, `complete` or `average` are usually safer choices."
      )
    }

    function(m) {
      row_distance <- compute_row_distance(
        expression_mat = m,
        distance_method = row_clustering_distance
      )

      fastcluster::hclust(
        d = row_distance,
        method = row_clustering_method
      )
    }
  } else {
    FALSE
  }

  build_average_expression <- function() {
    expression_mat <- Seurat::AverageExpression(
      seurat_object,
      features = features_present,
      assays = assay,
      layer = layer,
      group.by = group_by
    )[[assay]]

    expression_mat <- base::as.matrix(expression_mat)

    present_levels <- grouping_levels[grouping_levels %in% colnames(expression_mat)]
    remaining_levels <- setdiff(colnames(expression_mat), present_levels)
    expression_mat <- expression_mat[, c(present_levels, remaining_levels), drop = FALSE]

    if (isTRUE(scale_rows)) {
      expression_mat <- expression_mat |>
        base::t() |>
        base::scale() |>
        base::t()
    }

    group_levels <- colnames(expression_mat)
    column_groups <- factor(group_levels, levels = group_levels)

    list(
      expression_mat = expression_mat,
      column_groups = column_groups,
      column_split = column_groups,
      top_annotation_mode = "block"
    )
  }

  build_cell_expression <- function(order_by_pseudotime = FALSE) {
    cell_metadata <- seurat_object[[]]

    grouping_values <- if (identical(group_by, "ident")) {
      Seurat::Idents(seurat_object)
    } else {
      if (!group_by %in% colnames(cell_metadata)) {
        stop(
          "`group_by = ", group_by,
          "` was not found in the Seurat metadata. Use `ident` or a valid metadata column name."
        )
      }

      cell_metadata[[group_by]]
    }

    if (isTRUE(order_by_pseudotime) && !pseudotime_var %in% colnames(cell_metadata)) {
      stop(
        "`pseudotime_var = ", pseudotime_var,
        "` was not found in the Seurat metadata."
      )
    }

    expression_mat <- assay_data[features_present, , drop = FALSE] |>
      base::as.matrix() |>
      base::t()

    data <- data.frame(
      group = grouping_values,
      stringsAsFactors = FALSE,
      row.names = rownames(expression_mat)
    )

    if (isTRUE(order_by_pseudotime)) {
      data$pseudotime <- cell_metadata[[pseudotime_var]]
      data <- data[order(data$pseudotime, decreasing = FALSE), , drop = FALSE]
      expression_mat <- expression_mat[rownames(data), , drop = FALSE]
    }

    if (isTRUE(scale_rows)) {
      expression_mat <- base::scale(expression_mat)
    }

    expression_mat <- base::t(expression_mat)
    column_groups <- make_grouping_factor(data$group)

    list(
      expression_mat = base::as.matrix(expression_mat),
      column_groups = column_groups,
      column_split = if (isTRUE(order_by_pseudotime)) NULL else column_groups,
      top_annotation_mode = if (isTRUE(order_by_pseudotime)) "per_cell" else "block"
    )
  }

  plot_data <- switch(
    plot_mode,
    average = build_average_expression(),
    cells = build_cell_expression(order_by_pseudotime = FALSE),
    pseudotime = build_cell_expression(order_by_pseudotime = TRUE)
  )

  expression_mat <- plot_data$expression_mat
  column_groups <- plot_data$column_groups
  column_split <- plot_data$column_split
  group_levels <- levels(column_groups)
  group_colors <- make_group_colors(group_levels, group_bar_fill_colors)

  is_scaled_matrix <- isTRUE(scale_rows) || identical(layer, "scale.data")

  if (is.null(colors)) {
    colors <- if (is_scaled_matrix) {
      c("#FF00FF", "black", "#FFFF00")
    } else {
      c("lightgrey", "#08519C")
    }
  }

  if (is.null(breaks)) {
    expr_range <- range(expression_mat, na.rm = TRUE, finite = TRUE)

    if (is_scaled_matrix) {
      breaks <- c(-2, 0, 2)
    } else if (diff(expr_range) == 0) {
      breaks <- expr_range + c(-0.5, 0.5)
    } else {
      breaks <- expr_range
    }
  }

  heatmap_colors <- resolve_heatmap_colors(
    colors = colors,
    breaks = breaks,
    color_steps = color_steps
  )
  col_fun <- circlize::colorRamp2(
    breaks = heatmap_colors$breaks,
    colors = heatmap_colors$colors
  )

  if (!is.null(split_rows_by)) {
    split_rows_by <- split_rows_by[rownames(expression_mat)]
  }

  presence_of_row_split <- make_row_split_annotation(
    expression_mat = expression_mat,
    split_rows_by = split_rows_by,
    row_names_fontsize = row_names_fontsize
  )

  genes_highlight <- make_gene_highlight_annotation(
    expression_mat = expression_mat,
    highlight_genes = highlight_genes
  )

  if (is.null(legend_name)) {
    legend_name <- if (!is_scaled_matrix) {
      "Normalized expression"
    } else if (plot_mode == "average") {
      "Average expression"
    } else {
      "Expression"
    }
  }

  top_annotation <- if (identical(plot_data$top_annotation_mode, "per_cell")) {
    ComplexHeatmap::HeatmapAnnotation(
      group = column_groups,
      col = list(group = group_colors),
      show_annotation_name = FALSE
    )
  } else {
    ComplexHeatmap::HeatmapAnnotation(
      foo = ComplexHeatmap::anno_block(
        gp = grid::gpar(fill = unname(group_colors[group_levels]))
      )
    )
  }

  p <- ComplexHeatmap::Heatmap(
    expression_mat,
    name = legend_name,

    # Columns
    column_split = column_split,
    cluster_columns = FALSE,
    cluster_column_slices = FALSE,
    show_column_dend = FALSE,
    show_column_names = FALSE,
    column_gap = grid::unit(1, "mm"),
    column_title_rot = 45,

    # Rows
    row_names_side = row_names_side,
    show_row_names = if (is.null(split_rows_by)) {
      isTRUE(show_row_names) && isTRUE(show_genes)
    } else {
      FALSE
    },
    row_names_gp = grid::gpar(fontsize = row_names_fontsize),

    # Hierarchical clustering
    cluster_rows = row_clustering_fun,
    show_row_dend = show_row_dend,
    row_dend_side = row_dend_side,
    row_dend_width = grid::unit(row_dend_width, "mm"),

    # Row annotation
    left_annotation = presence_of_row_split,
    right_annotation = genes_highlight,
    row_split = split_rows_by,
    row_gap = grid::unit(2, "mm"),

    # Colors
    col = col_fun,

    # Top annotation
    top_annotation = top_annotation,

    # Other
    border = border_col,
    use_raster = rasterise,
    raster_quality = raster_quality,
    raster_by_magick = raster_by_magick
  )

  if (is.null(title)) {
    p
  } else {
    ComplexHeatmap::draw(
      p,
      column_title = title,
      column_title_gp = grid::gpar(fontsize = 15)
    )
  }
}
