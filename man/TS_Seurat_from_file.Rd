% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TS_Seurat_from_file.R
\name{TS_Seurat_from_file}
\alias{TS_Seurat_from_file}
\title{Process CellRanger output files with Seurat}
\usage{
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
}
\arguments{
\item{input_data}{File path to the folder containing the three matrices used for Seurat}

\item{sample_ID}{Required ID-tag for your sample. Will be added as metadata in column "sample_ID". Handy in downstream analyses, e.g. integration, batch-correction etc.}

\item{remove_ig_genes}{Filters out Ig-genes prior to clustering. Which can help in finding biologically relevant clusters, as shown in Sundell et. al., 2022. Defaults to TRUE.}

\item{regular_clustering_dims}{What PCA components to use for building NN graph and dimensionality reduction with regular normalization/scaling. Defaults to NULL.}

\item{SCT_clustering_dims}{What PCA components to use for building NN graph and dimensionality reduction with SCTransformed data. Defaults to NULL.}

\item{additional_markers}{Additional gene names to plot to aid in exclusion of unwanted cell types. Default markers are "percent.mt", "percent.ribo", "percent.malat1", "CD19", "CD3E", "CD14", "CD56".}

\item{project_name}{Required by Seurat. Defaults to the same value as \link{sample_ID}. Added automatically as metadata in column "orig.ident"}

\item{min.cells}{Include features/genes detected in at least this many cells. Defaults to 3.}

\item{min.features}{Keep cells with at least this many features detected. Defaults to 200.}

\item{normalization.method}{Normalization method to use. Alternatives are "LogNormalize" (default), "CLR", "RC".}

\item{scale.factor}{Scale factor for cell-level normalization. Defaults to 10000}

\item{selection.method}{Method for finding variable features. Alternatives are "vst" (default), "mean.var.plot", "dispersion".}

\item{nn.method}{Nearest-neighbour method to use for clustering with regular normalization/scaling. Alternatives are "rann" or "annoy". Defaults to "rann".}

\item{nFeatures}{Number of variable features to calculate. Defaults to 2000}
}
\description{
Input the three matrices from CellRanger and analyse them automatically.
Will first run an analysis pipeline with log-normalized counts,
allow you to set parameters on-the-fly, remove unwanted clusters,
and then run SCTransform_v1 on the filtered data.

Based on Seurat, using the workflow published in Sundell et. al., 2022 (https://doi.org/10.1093/bfgp/elac044)
}
