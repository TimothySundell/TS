#' Not sure if this function is desirable anymore. So use it at your own risk
#' @description
#' Not sure what this function does
#'
#' @param vdj_file I guess it takes the input from cellranger vdj
#' @importFrom magrittr  %>%
#' @export
TS_count_VDJ <- function(vdj_file) {

  require(dplyr)
  require(magrittr)

  vdj_file %>%
    dplyr::group_by(barcode) %>%
    dplyr::mutate(
      # Label which IGH contigs are productive, unproductive and sterile
      HC_productive = dplyr::case_when(
        chain == "IGH" & productive == "True" ~ 1),
      HC_unproductive = dplyr::case_when(
        chain == "IGH" & productive == "False" & v_gene != "None" & j_gene != "None" & c_gene != "None" ~ 1),
      HC_sterile = dplyr::case_when(
        chain == "IGH" & productive == "False" & v_gene == "None" & j_gene != "None" & c_gene != "None" ~ 1),

      # Label which IGK contigs are productive, unproductive and sterile
      IGK_productive = dplyr::case_when(
        chain == "IGK" & productive == "True" ~ 1),
      IGK_unproductive = dplyr::case_when(
        chain == "IGK" & productive == "False" & v_gene != "None" & j_gene != "None" & c_gene != "None" ~ 1),
      IGK_sterile = dplyr::case_when(
        chain == "IGK" & productive == "False" & v_gene == "None" & j_gene != "None" & c_gene != "None" ~ 1),

      # Label which IGK contigs are productive, unproductive and sterile
      IGL_productive = dplyr::case_when(
        chain == "IGL" & productive == "True" ~ 1),
      IGL_unproductive = dplyr::case_when(
        chain == "IGL" & productive == "False" & v_gene != "None" & j_gene != "None" & c_gene != "None" ~ 1),
      IGL_sterile = dplyr::case_when(
        chain == "IGL" & productive == "False" & v_gene == "None" & j_gene != "None" & c_gene != "None" ~ 1),

      # Sum how many IGH that are productive, unproductive and sterile for each cell
      HC_productive_sum = sum(HC_productive, na.rm = T),
      HC_unproductive_sum = sum(HC_unproductive, na.rm = T),
      HC_sterile_sum = sum(HC_sterile, na.rm = T),

      # Sum how many IGK that are productive, unproductive and sterile for each cell
      IGK_productive_sum = sum(IGK_productive, na.rm = T),
      IGK_unproductive_sum = sum(IGK_unproductive, na.rm = T),
      IGK_sterile_sum = sum(IGK_sterile, na.rm = T),

      # Sum how many IGK that are productive, unproductive and sterile for each cell
      IGL_productive_sum = sum(IGL_productive, na.rm = T),
      IGL_unproductive_sum = sum(IGL_unproductive, na.rm = T),
      IGL_sterile_sum = sum(IGL_sterile, na.rm = T),

      # Sum how many light chain that are productive, unproductive and sterile for each cell
      LC_productive_sum = sum(c(IGK_productive, IGL_productive), na.rm = T),
      LC_unproductive_sum = sum(c(IGK_unproductive, IGL_unproductive), na.rm = T),
      LC_sterile_sum = sum(c(IGK_sterile, IGL_sterile), na.rm = T))
}
