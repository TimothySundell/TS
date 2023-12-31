% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TS_calculate_ratios.R
\name{TS_calculate_ratios}
\alias{TS_calculate_ratios}
\title{Calculate light chain ratios of a Seurat object, requires light chain data in metadata "light.chain"}
\usage{
TS_calculate_ratios(
  seurat_object = get(default_seurat_object),
  group_var = "ident"
)
}
\arguments{
\item{seurat_object}{The Seurat object you want to analyse}

\item{group_var}{What parameter to group by, has to exist as a metadata slot in the Seurat object. Defaults to "ident"}
}
\description{
Calculate light chain ratios of a Seurat object, requires light chain data in metadata "light.chain"
}
