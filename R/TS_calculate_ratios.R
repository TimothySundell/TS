#' Calculate light chain ratios of a Seurat object, requires light chain data in metadata "light.chain"
#'
#' @param seurat_object The Seurat object you want to analyse
#' @param group_var What parameter to group by, has to exist as a metadata slot in the Seurat object. Defaults to "ident"
#' @importFrom magrittr  %>%
#' @export
TS_calculate_ratios <- function(seurat_object = get(default_seurat_object), group_var = "ident") {
  a <- Seurat::FetchData(object = seurat_object, vars = c("light.chain", group_var))
  b <- a %>% dplyr::group_by(.data[[group_var]]) %>% dplyr::count(light.chain)
  c <- b %>% dplyr::filter(light.chain != "NA") %>% data.frame()
  d <- c %>% dplyr::filter(light.chain == "IGK") %>% dplyr::rename(IGK = n) %>% select(-"light.chain")
  e <- c %>% dplyr::filter(light.chain == "IGL") %>% dplyr::rename(IGL = n) %>% select(-"light.chain")
  f <- merge(d, e, by = group_var)
  f$ratio <- f$IGK/f$IGL
  return(f)
}
