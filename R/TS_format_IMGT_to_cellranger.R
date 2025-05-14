#' Convert IMGT VDJ data to Cell Ranger–like format
#'
#' Standardizes IMGT VDJ data to match Cell Ranger expectations.
#' Also checks for the presence of a grouping column and required group values.
#'
#' @param imgt_data Data frame from IMGT.
#' @param grouping_variable Name of the grouping column (character).
#' @param sample1 Grouping variable value for sample 1.
#' @param sample2 Grouping variable value for sample 2.
#'
#' @return A standardized data frame ready for use in `TS_plot_IG_usage()`.
#'
#' @export
TS_format_IMGT_to_cellranger <- function(imgt_data, grouping_variable, sample1, sample2) {
  # Check that grouping_variable column exists
  if (!grouping_variable %in% colnames(imgt_data)) {
    stop(paste0("The grouping column '", grouping_variable, "' is not present in the data."))
  }

  # Check that both sample1 and sample2 are present in the grouping column
  group_values <- unique(imgt_data[[grouping_variable]])
  if (!all(c(sample1, sample2) %in% group_values)) {
    stop(paste0("One or both of the specified samples ('", sample1, "', '", sample2,
                "') are not found in the grouping column '", grouping_variable, "'.\n",
                "Found values: ", paste(group_values, collapse = ", ")))
  }

  # Process the data
  res <- imgt_data %>%
    tidyr::separate_wider_delim(
      cols = V.GENE.and.allele,
      delim = ",",
      too_many = "drop",
      names = "v_gene_and_allele",
      cols_remove = FALSE
    ) %>%
    dplyr::mutate(
      v_gene = stringr::str_sub(v_gene_and_allele, start = 8, end = -6),
      v_family = stringr::str_sub(v_gene, start = 1, end = 5)
    ) %>%
    {
      if ("D.GENE.and.allele" %in% names(.)) {
        tidyr::separate_wider_delim(
          data = .,
          cols = D.GENE.and.allele,
          delim = ",",
          too_many = "drop",
          names = "d_gene_and_allele",
          cols_remove = FALSE
        ) %>%
          dplyr::mutate(
            d_gene = stringr::str_sub(d_gene_and_allele, start = 8, end = -6),
            d_family = stringr::str_sub(d_gene, start = 1, end = 5)
          )
      } else {
        message("D.GENE.and.allele column not found — skipping D gene processing.")
        .
      }
    } %>%
    tidyr::separate_wider_delim(
      cols = J.GENE.and.allele,
      delim = ",",
      too_many = "drop",
      names = "j_gene_and_allele",
      cols_remove = FALSE
    ) %>%
    dplyr::mutate(
      j_gene = stringr::str_sub(j_gene_and_allele, start = 8, end = -6),
      chain = dplyr::case_when(
        stringr::str_detect(v_gene, "^IGH") ~ "IGH",
        stringr::str_detect(v_gene, "^IGK") ~ "IGK",
        stringr::str_detect(v_gene, "^IGL") ~ "IGL",
        TRUE ~ NA
      )
    )

  return(res)
}
