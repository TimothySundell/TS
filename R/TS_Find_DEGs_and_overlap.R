#' Find the upregulated DEGs for each cluster, comparing two clusters. As well as the conserved markers between them.
#'
#' @description
#' Returns a list containing:
#' - Upregulated markers in cluster 1 compared to cluster 2
#' - Upregulated markers in cluster 2 compared to cluster 1
#'  - Shared markers between them, compared to all other clusters in the Seurat object
#'
#'  Depends on:
#'  - Seurat
#'  - magrittr
#'  - tidyverse
#'  - VennDetail
#'
#'  Automatically filters resulting DEG list so that '0 < p_val_adj > 0.05.
#'
#' @param seurat_object Your Seurat object
#' @param ident.1 Name of cluster 1
#' @param ident.2 Name of cluster 2
#' @param grouping_variable_name Name of the combined group. Could be anything. Defaults to 'tested_population'.
#' @param only_pos_markers Whether to only return positively upregulated markers for the specific cluster. Defaults to 'TRUE'
#' @param test_use Statistical test used for DE testing. Defaults to 'Wilcoxon ranked sum test'.
#' @param plot_venn Whether to plot a Venn diagram for the gene lists. Defaults to 'FALSE'
#' @param export_venn Whether to export the Venn diagram to working directory. Defaults to 'FALSE'
#' @param filename_pdf File name for exported Venn diagram.
#' @param width  Width of saved Venn diagram. Defaults to '8.27 inches', which is A5  landscape format.
#' @param height Height of saved Venn diagram. Defaults to '5.83 inches', which is A5 landscape format.
#' @importFrom magrittr  %>%
#' @export
TS_Find_DEGs_and_overlap <- function(seurat_object, ident.1, ident.2, grouping_variable_name = "tested_population", only_pos_markers = T, test_use = "wilcox", plot_venn = F, export_venn = F, filename_pdf = "Venn_output", width = 8.27, heigth = 5.83){

  library(Seurat)
  library(magrittr)
  library(dplyr)
  library(VennDetail)

  name_1 <- paste0("Upregulated_", ident.1)
  name_2 <- paste0("Upregulated_", ident.2)
  name_3 <- paste0("Shared_markers")

  res <- list()

  # Find upregulated genes in ident.1
  cat("\nFinding genes upregulated in cluster", ident.1, "\n")
  res[[name_1]] <-
    Seurat::FindMarkers(
      object = seurat_object,
      ident.1 = ident.1,
      ident.2 = ident.2,
      only.pos = only_pos_markers,
      test.use = test_use
    ) %>% dplyr::filter(dplyr::between(p_val_adj, 0, 0.05)) %>%
    dplyr::mutate(gene = rownames(.))

  cat("\nFinding genes upregulated in cluster", ident.2, "\n")
  # Find upregulated genes in ident.2
  res[[name_2]] <-
    Seurat::FindMarkers(
      object = seurat_object,
      ident.1 = ident.2,
      ident.2 = ident.1,
      only.pos = only_pos_markers,
      test.use = test_use
    ) %>% dplyr::filter(dplyr::between(p_val_adj, 0, 0.05)) %>%
    dplyr::mutate(gene = rownames(.))

  # Find shared markers
  ## Add metadata to allow analyses
  cat("\nFinding genes shared between cluster", ident.1, "and cluster", ident.2, "\n")
  seurat_object <-
    Seurat::AddMetaData(
      object = seurat_object,
      metadata = Seurat::FetchData(object = seurat_object, vars = "ident") %>%
        dplyr::rename(group = ident) %>%
        dplyr::mutate(
          group = dplyr::case_when(
            group == ident.1 ~ grouping_variable_name,
            group == ident.2 ~ grouping_variable_name,
            .default = group
          )
        ),
      col.name = "group"
    )

  res[[name_3]] <-
    Seurat::FindConservedMarkers(
      object = seurat_object,
      ident.1 = ident.1,
      ident.2 = ident.2,
      grouping.var = "group"
    )

  temp_colname1 <- as.name(colnames(res[[3]][5]))

  res[[3]] <- res[[3]] %>%
    dplyr::filter(dplyr::between(!!temp_colname1, 0, 0.05))

  if(only_pos_markers) {
    temp_colname2 <- as.name(colnames(res[[3]][2]))
    res[[3]] <- res[[3]] %>% dplyr::filter(!!temp_colname2 > 0)
  } else{}

  # Calculate Venn diagram information
  cat("Drawing Venn diagram")
  res[["Venn.diagram"]] <-
    VennDetail::venndetail(x = lapply(
      X = res[1:3],
      FUN = function(x) {
        rownames(x)
      }
    ))
  if(plot_venn){
    grid::grid.newpage()
    graphics::plot(res[[4]])
  } else{}
  if(export_venn){
    grDevices::pdf(file = paste0(getwd(), "/", filename_pdf, ".pdf"), width = width, height = height)
    graphics::plot(res[[4]])
    grDevices::dev.off()
  } else{}

  return(res)

}
