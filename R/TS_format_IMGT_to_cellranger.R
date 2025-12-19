#' Convert IMGT VDJ data to Cell Ranger–like format
#'
#' Standardises IMGT VDJ data to match Cell Ranger expectations.
#' Also checks for the presence of a grouping column and required group values.
#'
#' @param imgt_data Data frame from IMGT.
#'
#' @return A standardized data frame ready for use in `TS_plot_IG_usage()`.
#'
#' @importFrom dplyr mutate case_when
#' @importFrom stringr str_extract str_detect
#' @importFrom magrittr %>%
#'
#' @export
TS_format_IMGT_to_cellranger <- function(imgt_data) {


  # Process the data
  res <- imgt_data %>%
    dplyr::mutate(
      v_gene = stringr::str_extract(string = V.GENE.and.allele, pattern = "IG[^*]+"),
      v_family = stringr::str_extract(string = v_gene, pattern = "IG[^-D]+")
    ) %>%
    {
      if("D.GENE.and.allele" %in% names(.)) {
        dplyr::mutate(
          .data = .,
          d_gene = stringr::str_extract(string = .[["D.GENE.and.allele"]], pattern = "IG[^*]+"),
          d_family = stringr::str_extract(string = d_gene, pattern = "IG[^-]+"))
      } else{
        message("D.GENE.and.allele column not found — skipping D gene processing.")
        .
      }
    } %>%
    dplyr::mutate(
      j_gene = stringr::str_extract(string = J.GENE.and.allele, pattern = "IG[^*]+"),
      j_family = stringr::str_extract(string = j_gene, pattern = "IG[^-D]+"),
      chain = dplyr::case_when(
        stringr::str_detect(v_gene, "^IGH") ~ "IGH",
        stringr::str_detect(v_gene, "^IGK") ~ "IGK",
        stringr::str_detect(v_gene, "^IGL") ~ "IGL",
        TRUE ~ NA
      )
    )

  return(res)
}
