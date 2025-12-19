# Package index

## All functions

- [`as.factor()`](https://timothysundell.github.io/TS/reference/IGHD_family_levels.md)
  : Human IGHD families in ascending numeric order
- [`as.factor()`](https://timothysundell.github.io/TS/reference/IGHD_levels.md)
  : Human IGHD genes according to their relative locus location
- [`as.factor()`](https://timothysundell.github.io/TS/reference/IGHJ_levels.md)
  : Human IGHJ genes according to their relative locus location
- [`as.factor()`](https://timothysundell.github.io/TS/reference/IGHV_family_levels.md)
  : Human IGHV families in ascending numeric order
- [`as.factor()`](https://timothysundell.github.io/TS/reference/IGHV_family_levels_grouped.md)
  : Human IGHV families in ascending numeric order with IGHV6 and 7
  grouped
- [`as.factor()`](https://timothysundell.github.io/TS/reference/IGHV_levels.md)
  : Human IGHV genes according to their relative locus location
- [`as.factor()`](https://timothysundell.github.io/TS/reference/IGKJ_levels.md)
  : Human IGKJ genes according to their relative locus location
- [`as.factor()`](https://timothysundell.github.io/TS/reference/IGKV_family_levels.md)
  : Human IGKV families in ascending numeric order
- [`as.factor()`](https://timothysundell.github.io/TS/reference/IGKV_levels.md)
  : Human IGKV genes according to their relative locus location
- [`as.factor()`](https://timothysundell.github.io/TS/reference/IGLC_levels.md)
  : Human IGLC genes according to their relative locus location
- [`as.factor()`](https://timothysundell.github.io/TS/reference/IGLJ_levels.md)
  : Human IGLJ genes according to their relative locus location
- [`as.factor()`](https://timothysundell.github.io/TS/reference/IGLV_family_levels.md)
  : Human IGLV families in ascending numeric order
- [`as.factor()`](https://timothysundell.github.io/TS/reference/IGLV_levels.md)
  : Human IGLV genes according to their relative locus location
- [`TS_DEenrichRPlot()`](https://timothysundell.github.io/TS/reference/TS_DEenrichRPlot.md)
  : A wrapper around the DEnrichRPlot from the 'Seurat' package
- [`TS_FindConservedMarkers()`](https://timothysundell.github.io/TS/reference/TS_FindConservedMarkers.md)
  : Find conserved markers within two groups
- [`TS_Find_DEGs_and_overlap()`](https://timothysundell.github.io/TS/reference/TS_Find_DEGs_and_overlap.md)
  : Find the upregulated DEGs for each cluster, comparing two clusters.
  As well as the conserved markers between them.
- [`TS_IGHV_usage_prop_with_stats()`](https://timothysundell.github.io/TS/reference/TS_IGHV_usage_prop_with_stats.md)
  : Plots IGHV usage between two groups using Fisher's Exact Test
- [`TS_Seurat_from_file()`](https://timothysundell.github.io/TS/reference/TS_Seurat_from_file.md)
  : Process CellRanger output files with Seurat
- [`TS_add_sig_segments_and_labels()`](https://timothysundell.github.io/TS/reference/TS_add_sig_segments_and_labels.md)
  : Adds test results to column charts
- [`TS_calculate_ratios()`](https://timothysundell.github.io/TS/reference/TS_calculate_ratios.md)
  : Calculate light chain ratios of a Seurat object, requires light
  chain data in metadata "light.chain"
- [`TS_cellnumbers_per_PC_and_resolution()`](https://timothysundell.github.io/TS/reference/TS_cellnumbers_per_PC_and_resolution.md)
  : Calculate cellnumbers per cluster for each combination of 'number or
  PCs' and 'resolution'
- [`TS_compare_PC_and_resolution()`](https://timothysundell.github.io/TS/reference/TS_compare_PC_and_resolution.md)
  : Compare number of clusters for combinations of number of PCs and
  resolution
- [`TS_compute_pairwise_fisher()`](https://timothysundell.github.io/TS/reference/TS_compute_pairwise_fisher.md)
  : Compute pairwise Fisher tests for gene usage
- [`TS_count_VDJ()`](https://timothysundell.github.io/TS/reference/TS_count_VDJ.md)
  : Not sure if this function is desirable anymore. So use it at your
  own risk
- [`TS_count_VDJ_7.0()`](https://timothysundell.github.io/TS/reference/TS_count_VDJ_7.0.md)
  : Count VDJ from Cellranger 7.0
- [`TS_dataframe_to_fasta()`](https://timothysundell.github.io/TS/reference/TS_dataframe_to_fasta.md)
  : Converts a dataframe and saves it in Fasta format
- [`TS_export_cluster_ids()`](https://timothysundell.github.io/TS/reference/TS_export_cluster_ids.md)
  : Export cell barcodes for each cluster
- [`TS_export_named_ids()`](https://timothysundell.github.io/TS/reference/TS_export_named_ids.md)
  : Export identities for each cluster, which are not numerical
- [`TS_extract_ig_mean_median_rank_per_PC()`](https://timothysundell.github.io/TS/reference/TS_extract_ig_mean_median_rank_per_PC.md)
  : Extract the mean and median values of Ig-genes from a list of PCA
  feature loadings
- [`TS_extract_ig_ranks_per_PC()`](https://timothysundell.github.io/TS/reference/TS_extract_ig_ranks_per_PC.md)
  : Extract the 'rank' for each of the Ig-genes
- [`TS_extract_ranked_genes_per_PC()`](https://timothysundell.github.io/TS/reference/TS_extract_ranked_genes_per_PC.md)
  : Extract an ordered list of feature importance for each of the first
  50 principal components
- [`TS_featureplot_to_pdf()`](https://timothysundell.github.io/TS/reference/TS_featureplot_to_pdf.md)
  : Export featureplot as pdf to your working directory
- [`TS_featureplot_to_pdf_A4()`](https://timothysundell.github.io/TS/reference/TS_featureplot_to_pdf_A4.md)
  : Export featureplot to pdf in A4 format
- [`TS_featureplot_to_pdf_no_leg_no_axes()`](https://timothysundell.github.io/TS/reference/TS_featureplot_to_pdf_no_leg_no_axes.md)
  : Export featureplot to pdf without legend and axes
- [`TS_featureplots()`](https://timothysundell.github.io/TS/reference/TS_featureplots.md)
  : Output one or multiple featureplots
- [`TS_filter_marker_list()`](https://timothysundell.github.io/TS/reference/TS_filter_marker_list.md)
  : Filter list of markers from FindMarkers function, needs to be output
  into new df
- [`TS_format_IMGT_to_cellranger()`](https://timothysundell.github.io/TS/reference/TS_format_IMGT_to_cellranger.md)
  : Convert IMGT VDJ data to Cell Ranger–like format
- [`TS_plot_3D()`](https://timothysundell.github.io/TS/reference/TS_plot_3D.md)
  : Plot your Seurat object in 3D.
- [`TS_plot_IG_usage()`](https://timothysundell.github.io/TS/reference/TS_plot_IG_usage.md)
  : Plots Ig usage between two groups using Fisher's Exact Test
- [`TS_plot_IG_usage_legacy()`](https://timothysundell.github.io/TS/reference/TS_plot_IG_usage_legacy.md)
  : Plots Ig usage between two groups using Fisher's Exact Test
- [`TS_plot_column_chart()`](https://timothysundell.github.io/TS/reference/TS_plot_column_chart.md)
  : Create column chart with ggplot2
- [`TS_plot_column_pattern_chart()`](https://timothysundell.github.io/TS/reference/TS_plot_column_pattern_chart.md)
  : Create column chart with pattern with ggplot2
- [`TS_plot_counts_boxplot()`](https://timothysundell.github.io/TS/reference/TS_plot_counts_boxplot.md)
  : Plot counts for a gene in a Seurat object as boxplot
- [`TS_plot_counts_violin()`](https://timothysundell.github.io/TS/reference/TS_plot_counts_violin.md)
  : Plot counts for a gene in a Seurat object as violinplot
- [`TS_plot_enrichR()`](https://timothysundell.github.io/TS/reference/TS_plot_enrichR.md)
  : Plot the output from TS_run_enrichR
- [`TS_project_scRNAseq_dataset()`](https://timothysundell.github.io/TS/reference/TS_project_scRNAseq_dataset.md)
  : Projects one Seurat object onto another
- [`TS_project_scRNAseq_dataset_export_query()`](https://timothysundell.github.io/TS/reference/TS_project_scRNAseq_dataset_export_query.md)
  : Projects one Seurat object onto another and exports resulting
  dataframe
- [`TS_run_enrichR()`](https://timothysundell.github.io/TS/reference/TS_run_enrichR.md)
  : Wrapper around the enrichR package.
- [`TS_seurat_object_stats()`](https://timothysundell.github.io/TS/reference/TS_seurat_object_stats.md)
  : Prints stats about your Seurat object
- [`TS_svg_from_featureplot()`](https://timothysundell.github.io/TS/reference/TS_svg_from_featureplot.md)
  : Export featureplot as svg to your working directory
- [`TS_svg_from_plot()`](https://timothysundell.github.io/TS/reference/TS_svg_from_plot.md)
  : Export any plot as svg
- [`TS_visdimloadings()`](https://timothysundell.github.io/TS/reference/TS_visdimloadings.md)
  : Plot gene contribution to each principal component
- [`ram_usage()`](https://timothysundell.github.io/TS/reference/ram_usage.md)
  : RAM usage of objects.
