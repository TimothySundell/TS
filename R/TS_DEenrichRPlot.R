#' A wrapper around the DEnrichRPlot from the 'Seurat' package
#'
#' @param seurat_object Your 'Seurat' object. Defaults to 'default_seurat_object'
#' @param ident.1 Name of cluster 1. Defaults to '0'
#' @param ident.2 Name of cluster 2. Defaults to '1'
#' @param test.use The test you want to use for DE testing. Defaults to 'Wilcoxons ranked sum' test.
#' @param max.genes The maximum number of genes to use. Defaults to 500
#' @export
TS_DEenrichRPlot <-
  function(seurat_object = get(default_seurat_object),
           ident.1 = 0,
           ident.2 = 1,
           test.use = "wilcox",
           max.genes = 500) {

    require(enrichR)
    require(Seurat)

    return(
      Seurat::DEenrichRPlot(
        object = seurat_object,
        ident.1 = ident.1,
        ident.2 = ident.2,
        test.use = "negbinom",
        enrich.database = "GO_Biological_Process_2021",
        max.genes = max.genes
      )
    )
  }
