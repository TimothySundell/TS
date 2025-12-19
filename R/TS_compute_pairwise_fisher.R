#' Compute pairwise Fisher tests for gene usage
#'
#' @description
#' Given a data.frame of gene counts per group, performs Fisher's Exact Test for all pairwise group combinations for each gene.
#'
#' @param df A data.frame or tibble containing at least two columns: one for genes and one for groups.
#' @param gene_var String; name of the column in `df` that identifies the gene (e.g. "gene_list").
#' @param group_var String; name of the column in `df` that identifies the group (e.g. "named_clusters").
#' @param p_adj_method String; method for multiple-testing correction passed to `p.adjust` (default: "hochberg").
#'
#' @return A tibble with columns: gene, group1, group2, p_value, p_adj, p_sig
#' @importFrom dplyr %>% sym enquo mutate group_by summarise ungroup pull count filter select
#' @importFrom tidyr pivot_wider complete unnest nest
#' @importFrom purrr map map_dfr
#' @importFrom stats fisher.test p.adjust
#' @importFrom tibble deframe
#' @export
TS_compute_pairwise_fisher <- function(df, gene_var, group_var, p_adj_method = "hochberg") {
  # Convert to symbols
  gene_sym <- rlang::sym(gene_var)
  group_sym <- rlang::sym(group_var)

  # Summarise counts per gene × group, filling missing combinations with zero
  # Drop genes that are zero across all groups
  count_df <- df %>%
    dplyr::count(!!gene_sym, !!group_sym) %>%
    tidyr::complete(!!gene_sym, !!group_sym, fill = list(n = 0)) %>%
    dplyr::group_by(!!gene_sym) %>%
    dplyr::filter(sum(n) > 0) %>%
    dplyr::ungroup()

  # Compute total counts per group
  group_totals <- count_df %>%
    dplyr::group_by(!!group_sym) %>%
    dplyr::summarise(total = sum(n), .groups = "drop") %>%
    tibble::deframe()

  # Identify unique group names
  groups <- names(group_totals)

  # Nest by gene
  nested <- count_df %>%
    dplyr::group_by(!!gene_sym) %>%
    tidyr::nest(data = c(!!group_sym, n))

  # For each gene, run all pairwise Fisher tests
  results <- nested %>%
    dplyr::mutate(test_results = purrr::map(data, function(subdf) {
      # Create a named vector of counts
      counts <- setNames(subdf$n, subdf[[group_var]])
      # Use precomputed group totals
      totals <- group_totals
      # Total others per group = total cells minus count for that gene
      others <- totals - counts

      # All pairwise combinations
      combos <- combn(groups, 2, simplify = FALSE)
      purrr::map_dfr(combos, function(pair) {
        g1 <- pair[1]; g2 <- pair[2]

        c1 <- counts[g1]
        c2 <- counts[g2]

        # Skip comparisons where the gene is absent in both groups
        if (c1 == 0 && c2 == 0) {
          return(NULL)
        }

        mat <- matrix(c(
          c1, c2,
          others[g1], others[g2]
        ), nrow = 2)

        pval <- stats::fisher.test(mat)$p.value
        tibble::tibble(
          group1 = g1,
          group2 = g2,
          p_value = pval
        )
      })
    })) %>%
    tidyr::unnest(test_results) %>%
    dplyr::group_by(group1, group2) %>%
    dplyr::mutate(
      p_adj = stats::p.adjust(p_value, method = p_adj_method),
      p_sig = dplyr::case_when(
        p_adj <= .0001 ~ "****",
        p_adj <= .001 ~ "***",
        p_adj <= .01 ~ "**",
        p_adj <= .05 ~ "*",
        p_adj > .05 ~ "ns"
      )
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(
      gene = !!gene_sym,
      group1, group2, p_value, p_adj, p_sig
    )

  return(results)
}
