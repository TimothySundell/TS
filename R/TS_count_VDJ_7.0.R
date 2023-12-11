#' Count VDJ from Cellranger 7.0
#'
#' @description
#' Not sure if function is usable.
#' @param vdj_file VDJ table which is output from the Cellranger suite.
#' @export
TS_count_VDJ_7.0 <- function(vdj_file) {

  require(dplyr)

  vdj_file %>%
    dplyr::group_by(barcode) %>%
    dplyr::mutate(
      # Label which IGH contigs are productive, unproductive and sterile
      HC_productive = dplyr::case_when(
        chain == "IGH" & productive == "TRUE" ~ 1),
      HC_unproductive = dplyr::case_when(
        chain == "IGH" & productive == "FALSE" & is.na(v_gene) == FALSE & is.na(j_gene) == FALSE & is.na(c_gene) == FALSE ~ 1),
      HC_sterile = dplyr::case_when(
        chain == "IGH" & productive == "FALSE" & is.na(v_gene) == TRUE & is.na(j_gene) == FALSE & is.na(c_gene) == FALSE ~ 1),

      # Label which IGK contigs are productive, unproductive and sterile
      IGK_productive = dplyr::case_when(
        chain == "IGK" & productive == "TRUE" ~ 1),
      IGK_unproductive = dplyr::case_when(
        chain == "IGK" & productive == "FALSE" & is.na(v_gene) == FALSE & is.na(j_gene) == FALSE & is.na(c_gene) == FALSE ~ 1),
      IGK_sterile = dplyr::case_when(
        chain == "IGK" & productive == "FALSE" & is.na(v_gene) == TRUE & is.na(j_gene) == FALSE & is.na(c_gene) == FALSE ~ 1),

      # Label which IGK contigs are productive, unproductive and sterile
      IGL_productive = dplyr::case_when(
        chain == "IGL" & productive == "TRUE" ~ 1),
      IGL_unproductive = dplyr::case_when(
        chain == "IGL" & productive == "FALSE" & is.na(v_gene) == FALSE & is.na(j_gene) == FALSE & is.na(c_gene) == FALSE ~ 1),
      IGL_sterile = dplyr::case_when(
        chain == "IGL" & productive == "FALSE" & is.na(v_gene) == TRUE & is.na(j_gene) == FALSE & is.na(c_gene) == FALSE ~ 1),

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
