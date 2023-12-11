#' Find conserved markers within two groups
#'
#' @description
#' Not sure this function does what I want it to...
#'
#' @param seurat_object Your Seurat object
#' @param ident.1 Name for cluster 1
#' @param ident.2 Name for cluster 2
#' @param grouping_variable_name Name of the combined clusters, can be anything. Defaults to 'tested_population'.
#' @param only_pos_markers Whether to only highlight positively upregulated markers between the grouped population and the rest. Defaults to 'TRUE'.
#' @importFrom magrittr  %>%
#' @export
TS_FindConservedMarkers <- function(seurat_object, ident.1, ident.2, grouping_variable_name = "tested_population", only_pos_markers = T){

  require(Seurat)
  require(magrittr)
  require(dplyr)


  if(only_pos_markers){
    seurat_object <-
      Seurat::AddMetaData(
        object = seurat_object,
        metadata = Seurat::FetchData(object = seurat_object, vars = "ident") %>%
          dplyr::rename(group = ident) %>%
          dplyr::mutate(group = dplyr::case_when(
            group == ident.1 ~ grouping_variable_name,
            group == ident.2 ~ grouping_variable_name,
            .default = group
          )),
        col.name = "group"
      )

    temp <- Seurat::FindConservedMarkers(object = seurat_object, ident.1 = ident.1, ident.2 = ident.2, grouping.var = "group")
    temp_colname <- as.name(colnames(temp[2]))
    temp <- temp %>% dplyr::filter(!!temp_colname > 0)

    return(temp)

  } else{

    seurat_object <-
      Seurat::AddMetaData(
        object = seurat_object,
        metadata = Seurat::FetchData(object = seurat_object, vars = "ident") %>%
          dplyr::rename(group = ident) %>%
          dplyr:: mutate(group = case_when(
            group == ident.1 ~ grouping_variable_name,
            group == ident.2 ~ grouping_variable_name,
            .default = group
          )),
        col.name = "group"
      )

    temp <- Seurat::FindConservedMarkers(object = seurat_object, ident.1 = ident.1, ident.2 = ident.2, grouping.var = "group")
    return(temp)
  }
}
