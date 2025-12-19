#' Plots Ig usage between two groups using Fisher's Exact Test
#'
#' @description
#' DEPRACATED! Instead, use TS_plot_IG_usage()
#' Plot the Ig usage from two groups and perform Fisher's Exact Test to find significant differences.
#' Maximum number of groups = 2
#'
#' @param input_file A dataframe (e.g. from 10X Genomics Cellranger VDJ output) that includes a grouping column.
#' @param sample1 Grouping variable value for sample 1
#' @param sample2 Grouping variable value for sample 2
#' @param compare_what Genes to compare (e.g. "IGHV", "IGKV"). Defaults to "IGHV".
#' @param grouping_variable Column used to group data. Defaults to "named_clusters"
#' @param export_pdf Logical. If TRUE, exports the plot to a PDF. Defaults to FALSE
#' @param return_only_data Logical. If TRUE, returns only the calculated frequencies and the test matrices. Defaults to FALSE.
#' @param file_name_suffix Character string. Optional suffix to be added to the exported pdf.
#' @param custom_gene_levels Character vector. Optional, vector used to set gene list levels.
#' @param ylim_max Numeric. Optional upper limit for y-axis. ggplot2 will adjust the y-axis scale to fit the data.
#' @param p.adj.method Method for adjust for multiple tests. Defaults to 'Hochberg'. Alternative is 'FDR'.
#' @param fill_colours A character vector of fill colours for the two groups. Defaults to c("grey80", "grey50").
#' @param border_colour Colour of bar borders. Defaults to "black".
#' @param ggplot_theme Optional ggplot2 theme to apply. Defaults to "theme_classic".
#' @param plot_width Numeric. Width of exported plot in inches. Used if 'export_pdf = TRUE'. Defaults to 15.
#' @param plot_height Numeric. Height of exported plot in inches. Used if 'export_pdf = TRUE'. Defaults to 5.
#' @param remove_genes_present_in_one_dataset Logical. If TRUE, removes genes not shared by both groups. Defaults to FALSE
#' @param plot_zeroes Logical. If FALSE, genes with 0 counts in one group are not plotted, resulting in a wider bar for the other group.
#' @param file_type Set type of file/annotation input. Defaults to "IMGT". Options are "IMGT" and "AIRR".

#'
#' @return Returns a ggplot object, or a data.frame containing counts and test results if ' return_only_data = TRUE'.
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr group_by summarise mutate filter arrange pull rename count n_distinct ungroup case_when select all_of
#' @importFrom tidyr complete pivot_wider
#' @importFrom ggplot2 ggplot aes geom_col scale_fill_manual theme labs ylab ylim
#' @importFrom ggrepel geom_label_repel
#' @importFrom stats fisher.test p.adjust
#' @importFrom purrr pluck map2_dbl
#' @importFrom stringr str_sub
#' @importFrom rlang .data
#' @export
#'
#' @examples
#' # TS_plot_IG_usage(input_file = VDJ_data, sample1 = "group_A", sample2 = "group_B", compare_what = "IGHV")

## FIX:
# Remove legend title

TS_plot_IG_usage_legacy <-
  function(input_file,
           sample1,
           sample2,
           compare_what = "IGHV",
           grouping_variable = "named_clusters",
           export_pdf = F,
           return_only_data = F,
           file_name_suffix = NULL,
           custom_gene_levels = NULL,
           ylim_max = NULL,
           p.adj.method = "hochberg",
           fill_colours = c("gray80", "gray50"),
           border_colour = "black",
           ggplot_theme = ggplot2::theme_classic(),
           plot_width = 15,
           plot_height = 5,
           remove_genes_present_in_one_dataset = F,
           plot_zeroes = T,
           file_type = "IMGT") {

    options(dplyr.summarise.inform = F) # TMI!

    #######################
    ## Conform file_type ##
    #######################

    # IMGT
    if(file_type == "IMGT"){
      input_file <- TS_format_IMGT_to_cellranger(imgt_data = input_file)
    }

    # AIRR
    if(file_type == "AIRR"){
      stop("AIRR file format is not yet supported")
      # input_file <- TS_format_AIRR_to_cellranger()
    }

    ########################
    ## Set variable names ##
    ########################
    total_sample1 <- paste0("total_", sample1)
    total_sample2 <- paste0("total_", sample2)

    count_sample1 <- paste0("count_", sample1)
    count_sample2 <- paste0("count_", sample2)

    frequency_sample1 <- paste0("frequency_", sample1)
    frequency_sample2 <- paste0("frequency_", sample2)

    ######################
    ## Set up variables ##
    ######################

    # <--- Check compare_what is supported ---> #
    # Supported comparisons
    supported_comparisons <- c("IGHV", "IGHV_FAMILY", "IGHV_FAMILY_GROUPED", "IGHD", "IGHD_FAMILY", "IGHJ", "IGKV", "IGKV_FAMILY", "IGKJ", "IGLV", "IGLV_FAMILY", "IGLJ", "IGLC")
    compare_what <- toupper(compare_what) # Shift to upper case

    if(!compare_what %in% supported_comparisons){
      stop(paste0("\nYour choice of compare_what = ", compare_what, " is not supported. \n Options include: \n'IGHV', ' IGHV_family', 'IGHV_family_grouped', 'IGHJ', \n'IGKV', 'IGKV_family', 'IGKJ', \n'IGLV', 'IGLV_family', 'IGLJ', 'IGLC'\nCapital doesn't matter"))
    }

    # <--- IGH ---> #
    # IGHV usage
    if(compare_what == "IGHV"){
      input_file <- input_file %>%
        dplyr::rename(gene_list = v_gene)
      gene_list_levels <- IGHV_levels
      plot_title <- paste0("IGHV usage - ", sample1, " vs ", sample2)
    }

    # IGHV family usage
    if(compare_what == "IGHV_FAMILY"){
      input_file <- input_file %>%
        dplyr::mutate(gene_list = stringr::str_sub(string = v_gene, start = 1, end = 5))
      gene_list_levels <- IGHV_family_levels
      plot_title <- paste0("IGHV family usage - ", sample1, " vs ", sample2)
    }

    # IGHV family usage - with IGHV6 and 7 grouped
    if(compare_what == "IGHV_FAMILY_GROUPED"){
      input_file <- input_file %>%
        dplyr::mutate(gene_list = stringr::str_sub(string = v_gene, start = 1, end = 5),
                      gene_list = case_when(gene_list == "IGHV6" ~ "IGHV6/7",
                                            gene_list == "IGHV7" ~ "IGHV6/7",
                                            .default = gene_list))
      gene_list_levels <- IGHV_family_levels_grouped
      plot_title <- paste0("IGHV family usage - ", sample1, " vs ", sample2)
    }
    # IGHD usage
    if(compare_what == "IGHD"){
      input_file <- input_file %>%
        dplyr::rename(gene_list = d_gene)
      gene_list_levels <- IGHD_levels
      plot_title <- paste0("IGHD usage - ", sample1, " vs ", sample2)
      if(file_type == "cellranger"){
        message("You are comparing IGHD usage with data from Cell Ranger. This is not advisable to do as they usually map <50% of the genes\nAnyway, continuing...")
      }
    }

    # IGHD family usage
    if(compare_what == "IGHD_FAMILY"){
      input_file <- input_file %>%
        dplyr::mutate(gene_list = stringr::str_sub(string = d_gene, start = 1, end = 5))
      gene_list_levels <- IGHD_family_levels
      plot_title <- paste0("IGHD family usage - ", sample1, " vs ", sample2)
      if(file_type == "cellranger"){
        message("You are comparing IGHD usage with data from Cell Ranger. This is not advisable to do as they usually map <50% of the genes\nAnyway, continuing...")
      }
    }
    # IGHJ usage
    if(compare_what == "IGHJ"){
      input_file <- input_file %>%
        dplyr::rename(gene_list = j_gene)
      gene_list_levels <- IGHJ_levels
      plot_title <- paste0("IGHJ usage - ", sample1, " vs ", sample2)
    }

    # <--- IGK ---> #
    # IGKV usage
    if(compare_what == "IGKV"){
      input_file <- input_file %>%
        dplyr::rename(gene_list = v_gene)
      gene_list_levels <- IGKV_levels
      plot_title <- paste0("IGKV usage - ", sample1, " vs ", sample2)
    }

    # IGKV family usage
    if(compare_what == "IGKV_FAMILY"){
      input_file <- input_file %>%
        dplyr::mutate(gene_list = stringr::str_sub(string = v_gene, start = 1, end = 5))
      gene_list_levels <- IGKV_family_levels
      plot_title <- paste0("IGKV family usage - ", sample1, " vs ", sample2)
    }

    # IGKJ usage
    if(compare_what == "IGKJ"){
      input_file <- input_file %>%
        dplyr::rename(gene_list = j_gene)
      gene_list_levels <- IGKJ_levels
      plot_title <- paste0("IGKJ usage - ", sample1, " vs ", sample2)
    }

    # <--- IGL ---> #
    # IGLV usage
    if(compare_what == "IGLV"){
      input_file <- input_file %>%
        dplyr::rename(gene_list = v_gene)
      gene_list_levels <- IGLV_levels
      plot_title <- paste0("IGLV usage - ", sample1, " vs ", sample2)
    }

    # IGLV family usage
    if(compare_what == "IGLV_FAMILY"){
      input_file <- input_file %>%
        dplyr::mutate(gene_list = stringr::str_sub(string = v_gene, start = 1, end = 5))
      gene_list_levels <- IGLV_family_levels
      plot_title <- paste0("IGLV family usage - ", sample1, " vs ", sample2)
    }

    # IGLJ usage
    if(compare_what == "IGLJ"){
      input_file <- input_file %>%
        dplyr::rename(gene_list = j_gene)
      gene_list_levels <- IGLJ_levels
      plot_title <- paste0("IGLJ usage - ", sample1, " vs ", sample2)
    }

    # IGLC usage
    if(compare_what == "IGLC"){

      if(file_type == "IMGT"){
        stop("IMGT does not call constant region genes, therefore it is missing in the data")
      }

      input_file <- input_file %>%
        dplyr::rename(gene_list = c_gene)
      gene_list_levels <- IGLC_levels
      plot_title <- paste0("IGLC usage - ", sample1, " vs ", sample2)
    }

    # If custom gene_list_levels is supplied
    if(!is.null(custom_gene_levels)){
      gene_list_levels <- custom_gene_levels
    }

    # <--- General settings ---> #
    plot_subtitle <- paste0("% of max")


    #############################
    ## Filter input data frame ##
    #############################

    # Retrieve genes present in both groups
    # Only validated for IGHV
    if(remove_genes_present_in_one_dataset){

      present_genes <- input_file %>%
        dplyr::select(gene_list, dplyr::all_of(grouping_variable)) %>%
        dplyr::filter(.data[[grouping_variable]] %in% c(sample1, sample2)) %>%
        dplyr::filter(!is.na(gene_list)) %>%
        dplyr::group_by(.data[[grouping_variable]]) %>%
        dplyr::count(gene_list) %>%
        dplyr::group_by(gene_list) %>%
        dplyr::summarise(n_groups = dplyr::n_distinct(.data[[grouping_variable]]), .groups = "drop") %>%
        dplyr::filter(n_groups == 2) %>%
        dplyr::pull(gene_list)

      # Filter data
      input_file <- input_file %>%
        dplyr::filter(.data[[grouping_variable]] %in% c(sample1, sample2)) %>%
        dplyr::select(gene_list, dplyr::all_of(grouping_variable)) %>%
        dplyr::mutate(
          gene_list = factor(gene_list, levels = gene_list_levels[gene_list_levels %in% present_genes]),
          grouping_variable = factor(.data[[grouping_variable]], levels = c(sample1, sample2))) %>%
        dplyr::filter(!is.na(gene_list) & !is.na(.data[[grouping_variable]])) %>%
        dplyr::group_by(grouping_variable) %>%
        dplyr::count(gene_list)
    } else{

      # Filter data
      input_file <- input_file %>%
        dplyr::filter(.data[[grouping_variable]] %in% c(sample1, sample2)) %>% # Keep only cells from chosen groups
        dplyr::select(gene_list, dplyr::all_of(grouping_variable)) %>% # Remove excess columns
        dplyr::mutate(
          gene_list = factor(gene_list, levels = gene_list_levels),
          grouping_variable = factor(.data[[grouping_variable]], levels = c(sample1, sample2))) %>% # Convert to factors and set levels
        dplyr::filter(!is.na(gene_list) & !is.na(.data[[grouping_variable]])) %>% # Remove NA values in genes and groups
        dplyr::group_by(grouping_variable) %>%
        dplyr::count(gene_list)
    }

    ##########################################
    ### Calculate usage for stats and plot ###
    ##########################################
    # Filter data and calculate stats #
    res <- input_file %>%
      dplyr::group_by(gene_list, grouping_variable) %>%
      dplyr::summarise(count = n) %>%
      dplyr::ungroup() %>%
      dplyr::arrange(grouping_variable) %>%
      dplyr::group_by(grouping_variable) %>%
      dplyr::mutate(total = sum(count),
                    frequency = count / sum(count)) %>%
      tidyr::pivot_wider(
        names_from = grouping_variable,
        values_from = c(count, frequency, total),
        values_fill = 0
      ) %>%
      dplyr::mutate(
        !!total_sample1 := get(total_sample1) %>% max(),
        !!total_sample2 := get(total_sample2) %>% max()
      ) %>%
      {. ->> intermediate_export} %>%
      dplyr::ungroup() %>%
      dplyr::group_by(gene_list) %>%
      dplyr::mutate(
        cm = list(# Not required for the analyses per se, but nice to have if you want to go back and have a look at a specific contingency matrix.
          c(
            get(count_sample1),
            get(count_sample2),
            get(total_sample1) - get(count_sample1),
            get(total_sample2) - get(count_sample2)
          ) %>%
            matrix(nrow = 2)),
        p = c(
          get(count_sample1),
          get(count_sample2),
          get(total_sample1) - get(count_sample1),
          get(total_sample2) - get(count_sample2)
        ) %>%
          matrix(nrow = 2) %>%
          stats::fisher.test() %>%
          purrr::pluck("p.value"),
        .keep = "unused"
      ) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(
        p.adj = stats::p.adjust(p, method = p.adj.method),
        p.sig = dplyr::case_when(
          p.adj <= .0001 ~ "****",
          p.adj <= .001 ~ "***",
          p.adj <= .01 ~ "**",
          p.adj <= .05 ~ "*",
          p.adj >= .05 ~ "ns"
        )
      )
    # Return dataframe if specified
    if (return_only_data)
    {
      return(res)
    }

    # Otherwise plot the resulting IG-usage chart with stats-results denoted
    else{

      ##################
      ### Plot usage ###
      ##################

      # Code to fill NAs with 0.1, so that lines are printed in the final plot (that are still present)
      if(plot_zeroes){
        present_genes <- input_file %>%
          dplyr::pull(gene_list) %>%
          unique()
        gene_list_levels_plotting <- gene_list_levels[gene_list_levels %in% present_genes]
      } else{gene_list_levels_plotting <- gene_list_levels}

      # If custom gene_list_levels is supplied
      if(!is.null(custom_gene_levels)){
        gene_list_levels_plotting <- custom_gene_levels
      }

      if(length(gene_list_levels_plotting) == 0){
        stop("Error: The gene list for plotting has 0 entries")
      }

      plot_data <- input_file %>%
        dplyr::mutate(gene_list = factor(gene_list, levels = gene_list_levels_plotting)) %>%
        dplyr::group_by(grouping_variable, gene_list) %>%
        dplyr::summarise(count = n, .groups = "drop") %>%
        dplyr::mutate(count = as.double(count)) %>%
        tidyr::complete(gene_list, grouping_variable, fill = list(count = 0.01)) %>%
        dplyr::group_by(grouping_variable) %>%
        dplyr::mutate(proportion = count / sum(count) * 100)

      if(!is.null(ylim_max)){
        if(ylim_max < max(plot_data$proportion)){
          message("Warning: You have set your ylim_max to ", ylim_max, ". Increase the value of ylim_max to plot them.
  These values will NOT be plotted: ")
          print(plot_data[plot_data$proportion > ylim_max,])
        }
      }

      plot <- ggplot2::ggplot() +
        ggplot2::geom_col(
          data = plot_data,
          mapping = aes(x = gene_list,
                        y = proportion,
                        fill = grouping_variable),
          width = 0.8,
          position = position_dodge(),
          color = border_colour) +
        ggplot2::scale_fill_manual(
          values = fill_colours,
          name = "Sample",
          labels = c(
            paste0(sample1, "\n(n = ", intermediate_export[1,6], ")"),
            paste0(sample2, "\n(n = ", intermediate_export[1,7], ")"))) +
        ggplot_theme +
        ggplot2::theme(
          axis.text.x = element_text(angle = 60,
                                     hjust = 1,
                                     vjust = 1,
                                     family = "Arial"),
          text = element_text(color = "black")) +
        ggplot2::labs(title = plot_title,
                      subtitle = plot_subtitle) +
        ggplot2::ylab("%") + # NEW

        # Only set y limits if actually entered
        {if(!is.null(ylim_max)) ggplot2::ylim(0, ylim_max)} +

        # Add asterisks for significant differences
        ggrepel::geom_label_repel(
          data = dplyr::filter(res, p.adj < 0.05),
          label.padding = 0.01,
          parse = F,
          size = 5,
          label.size = NA,
          direction = "y",
          nudge_y = 0.01,
          mapping = aes(
            x = gene_list,
            y = purrr::map2_dbl((get(frequency_sample1) * 100), (get(frequency_sample2) * 100), max),
            label = p.sig
          )
        )
      if (export_pdf) {
        if(!is.null(file_name_suffix)){
        export_pdf_filename_suffix <- paste0("_", file_name_suffix)
        } else{export_pdf_filename_suffix <- NULL}
        export_pdf_filename <- paste0(compare_what, "_usage_", sample1, "_vs_", sample2, export_pdf_filename_suffix, "_", file_type, ".pdf")

        ggplot2::ggsave(
          filename = export_pdf_filename,
          plot = plot,
          units = "in",
          width = plot_width,
          height = plot_height
        )
        message(paste0("\n  Plot exported to:\n  ", getwd(), "\n\n  With the filename:\n  ", export_pdf_filename))
      } else{}
      return(plot)
    }
  }
