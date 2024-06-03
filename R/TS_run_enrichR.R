#' Wrapper around the enrichR package.
#'
#' @description
#' Performs enrichment analyses for a vector of genes
#' Depends on
#' -enrichR (and that their servers are up, heh)
#' - A live Internet connection
#' - magrittr
#' - dplyr
#' Returns a list of the results of each database tested for.
#'
#' @param gene_list Character vector containing genes to test
#' @param plot_enrichr Whether to plot the output using the enrichR function [plotEnrich]. Defaults to 'FALSE'
#' @param database_to_plot What database to plot the output from. Defaults to '"GO_Biological_Process_2023"'
#' @param databases What databases to test against. Defaults to 'c("GO_Biological_Process_2023", "GO_Molecular_Function_2023", "GO_Cellular_Component_2023")'
#' @importFrom magrittr  %>%
#' @export
TS_run_enrichR <- function(gene_list, plot_enrichr = F, database_to_plot = "GO_Biological_Process_2023", databases = c("GO_Biological_Process_2023", "GO_Molecular_Function_2023", "GO_Cellular_Component_2023")) {

  library(enrichR)
  library(magrittr)
  library(dplyr)

  temp.enrichr <- enrichR::enrichr(genes = gene_list, databases = databases) %>%
    lapply(X = ., FUN = function(x){
      x <- x %>% dplyr::filter(Adjusted.P.value < 0.05)
    }) %>%
    lapply(X =., FUN = function(x){
      x <- x %>% tidyr::separate_wider_regex(
        cols = Term,
        patterns = c(Term = ".*", "\\ ", 'GO_ID' = "\\(.*"))
    }) %>%
    lapply(X = ., FUN = function(x){
      x <- x %>% dplyr::mutate(GO_ID = base::substring(text = GO_ID, first = 2, last = 11))
    })

  if(plot_enrichr){
    print(enrichr::plotEnrich(df = temp.enrichr[[database_to_plot]]))
  } else{}

  return(temp.enrichr)
}
