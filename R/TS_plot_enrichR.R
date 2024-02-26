#' Plot the output from [TS_run_enrichR]
#'
#' @param input The output list from [TS::TS_run_enrichR]
#' @param sort_on What to sort top n terms on. Defaults to 'neg_log10', can also be 'combined_score'
#' @param n Number of terms to plot for each database. Defaults to '10'.
#' @importFrom magrittr  %>%
#' @export
TS_plot_enrichR <- function(input, sort_on = "neg_log10", n = 10){

  library(dplyr)
  library(ggplot2)
  library(stringr)
  library(magrittr)
  library(forcats)
  library(scales)

  plot_title <- c(substitute(input))
  plot_title <- gsub(pattern = "_", replacement = "\\ ", x = plot_title)

  named_colours <- scales::hue_pal()(3)
  names(named_colours) <- c("BP", "CC", "MF") # Allows colouring to be consistent even when there are missing "Classes"

  if(sort_on == "neg_log10") {

    data <- rbind(

      input$GO_Biological_Process_2023 %>%
        dplyr::mutate(
          Class = "BP",
          neg_log10 = -log10(Adjusted.P.value)) %>%
        dplyr::arrange(dplyr::desc(neg_log10)),

      input$GO_Molecular_Function_2023 %>%
        dplyr::mutate(
          Class = "MF",
          neg_log10 = -log10(Adjusted.P.value)) %>%
        dplyr::arrange(dplyr::desc(neg_log10)),

      input$GO_Cellular_Component_2023 %>%
        dplyr::mutate(
          Class = "CC",
          neg_log10 = -log10(Adjusted.P.value)) %>%
        dplyr::arrange(dplyr::desc(neg_log10))) %>%

      dplyr::mutate(Class = factor(Class, levels = c("BP", "MF", "CC"))) %>%
      tidyr::separate_wider_delim(cols = Genes, delim = ";", names_sep = "", too_few = "align_start") %>%

      {if(!"Genes1" %in% colnames(.)) dplyr::mutate(., "Genes1" = NA) else .} %>% # Added this block to allow matrices where no term has enough genes to run "unite"
      {if(!"Genes2" %in% colnames(.)) dplyr::mutate(., "Genes2" = NA) else .} %>%
      {if(!"Genes3" %in% colnames(.)) dplyr::mutate(., "Genes3" = NA) else .} %>%
      {if(!"Genes4" %in% colnames(.)) dplyr::mutate(., "Genes4" = NA) else .} %>%
      {if(!"Genes5" %in% colnames(.)) dplyr::mutate(., "Genes5" = NA) else .} %>%

      tidyr::unite(col = First_five_genes, Genes1, Genes2, Genes3, Genes4, Genes5, sep = ", ", remove = F, na.rm = T) %>%
      dplyr::mutate(Term = factor(Term, levels = Term)) %>%
      dplyr::group_by(Class) %>%
      dplyr::slice_max(order_by = neg_log10, n = n)

    plot <- ggplot2::ggplot(data = data, aes(x = neg_log10, y = forcats::fct_rev(Term), fill = Class)) +
      ggplot2::geom_vline(xintercept = -log10(0.05)) +
      ggplot2::geom_col() +
      ggplot2::scale_fill_manual(values = named_colours, aesthetics = "fill") + # added this to work with the named colours
      ggplot2::facet_grid(Class ~ ., scales = "free", space = "free") +
      ggplot2::geom_text(aes(x = 0, label = First_five_genes), hjust = 0, check_overlap = F)+
      ggplot2::ggtitle(
        label = stringr::str_sub(string = plot_title, start = 1L, end = -9L),
        subtitle = paste0("Ordered by ", substitute(sort_on))) +
      ggplot2::ylab(label = "Term")

    print(plot)
  }

  else if(sort_on == "combined_score") {

    data <- rbind(

      input$GO_Biological_Process_2023 %>%
        dplyr::mutate(
          Class = "BP",
          neg_log10 = -log10(Adjusted.P.value)) %>%
        dplyr::arrange(desc(Combined.Score)),

      input$GO_Molecular_Function_2023 %>%
        dplyr::mutate(
          Class = "MF",
          neg_log10 = -log10(Adjusted.P.value)) %>%
        dplyr::arrange(desc(Combined.Score)),

      input$GO_Cellular_Component_2023 %>%
        dplyr::mutate(
          Class = "CC",
          neg_log10 = -log10(Adjusted.P.value)) %>%
        dplyr::arrange(desc(Combined.Score))) %>%

      dplyr::mutate(Class = factor(Class, levels = c("BP", "MF", "CC"))) %>%
      tidyr::separate_wider_delim(cols = Genes, delim = ";", names_sep = "", too_few = "align_start") %>%

      {if(!"Genes1" %in% colnames(.)) dplyr::mutate(., "Genes1" = NA) else .} %>% # Added this block to allow matrices where no term has enough genes to run "unite"
      {if(!"Genes2" %in% colnames(.)) dplyr::mutate(., "Genes2" = NA) else .} %>%
      {if(!"Genes3" %in% colnames(.)) dplyr::mutate(., "Genes3" = NA) else .} %>%
      {if(!"Genes4" %in% colnames(.)) dplyr::mutate(., "Genes4" = NA) else .} %>%
      {if(!"Genes5" %in% colnames(.)) dplyr::mutate(., "Genes5" = NA) else .} %>%

      tidyr::unite(col = First_five_genes, Genes1, Genes2, Genes3, Genes4, Genes5, sep = ", ", remove = F, na.rm = T) %>%
      dplyr::mutate(Term = factor(Term, levels = Term)) %>%
      dplyr::group_by(Class) %>%
      dplyr::slice_max(order_by = Combined.Score, n = n)

    data$First_five_genes <- gsub(pattern = "\\,NA", replacement = "", x = data$First_five_genes)

    plot <- ggplot2::ggplot(data = data, aes(x = neg_log10, y = fct_rev(Term), fill = Class)) +
      ggplot2::geom_vline(xintercept = -log10(0.05)) +
      ggplot2::geom_col() +
      ggplot2::scale_fill_manual(values = named_colours, aesthetics = "fill") + # added this to work with the named colors
      ggplot2::facet_grid(Class ~ ., scales = "free", space = "free") +
      ggplot2::geom_text(aes(x = 0, label = First_five_genes), hjust = 0, check_overlap = F) +
      ggplot2::ggtitle(
        label = stringr::str_sub(string = plot_title, start = 1L, end = -9L),
        subtitle = paste0("Ordered by ", substitute(sort_on))) +
      ggplot2::ylab(label = "Term")

    print(plot)
  }
}
