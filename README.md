
<!-- README.md is generated from README.Rmd. Please edit that file -->

# TS - A collection of QoL functions <a href="https://github.com/TimothySundell/TS"><img src="man/Figures/TS_logo.png" align="right" height="138" /></a>

<!-- badges: start -->

[![](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
[![Project Status: Active - The project has reached a stable, usable
state and is being actively
developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![](https://img.shields.io/github/languages/code-size/TimothySundell/TS.svg)](https://github.com/TimothySundell/TS)
[![](https://img.shields.io/badge/doi-10.1093/bfgp/elac044-green.svg)](https://doi.org/10.1093/bfgp/elac044)
<!-- badges: end -->

TS was written to improve the quality of life for its creator, maybe it
can help you as well?

It contains a collection of functions that may be relevant people
analysing single-cell RNA-sequencing data, and analysing
Ig(BCR)/TCR-sequences.

**Please** give our method article a look at the **DOI** button above.
It shows the importance of removing BCR/TCR genes prior to performing
unsupervised clustering and downstream analyses of scRNA-seq data
containing lymphocytes.

You can find a reference webpage
[*here*](https://timothysundell.github.io/TS/), or under the *About*
section in the top right corner of this webpage.

------------------------------------------------------------------------

## Installation

To install the latest version of the TS package:

``` r
devtools::install_github(TimothySundell/TS)
```

------------------------------------------------------------------------

## Comments from the author

My plan is to keep this repository updated as much as possible. If
something is not working, please post in the *Issues* section.

As of now, the package does not pass *R CMD CHECK* as the tidyverse way
of coding leaves me a lot of unlinked global variables to reference.
Will maybe be fixed in the future.

------------------------------------------------------------------------

## The TS collection of functions can be divided into 4 main categories:

### Functions for plotting figures

- TS_DEenrichRPlot
- TS_featureplot_to_pdf
- TS_featureplot_to_pdf_A4
- TS_featureplot_to_pdf_no_leg_no_axes
- TS_featureplots
- TS_plot_3D
- TS_plot_counts_boxplot
- TS_plot_counts_violin
- TS_plot_enrichR
- TS_project_scRNAseq_dataset
- TS_svg_from_featureplot
- TS_svg_from_plot
- TS_visdimloadings

### Functions for analysing sc-V(D)J-sequencing data

- TS_IGHV_usage_prop_with_stats
- TS_calculate_ratios
- TS_count_VDJ
- TS_count_VDJ_7.0

### Functions for analysing scRNA-seq data

- TS_FindConservedMarkers
- TS_FindConservedMarkers
- TS_Find_DEGS_and_overlap
- TS_cellnumbers_per_PC_and_resolution
- TS_compare_PC_and_resolution
- TS_export_cluster_ids
- TS_export_named_ids
- TS_extract_ig_mean_median_rank_per_PC
- TS_extract_ig_ranks_per_PC
- TS_extract_ranked_genes_per_PC
- TS_filter_marker_list
- TS_project_scRNAseq_dataset_export_query
- TS_run_enrichR
- TS_seurat_object_stats

### General data-wrangling

- TS_dataframe_to_fasta
