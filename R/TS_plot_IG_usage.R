#' Plots Ig usage between two groups using Fisher's Exact Test
#'
#' @description
#' Plot the Ig usage from two groups and perform Fisher's Exact Test to find significant differences.
#'
#' @param input_file A dataframe (e.g. from 10X Genomics Cellranger VDJ output) that includes a grouping column.
#' @param compare_groups Character vector. Names of group values to compare.
#' @param compare_what Genes to compare (e.g. "IGHV", "IGKV"). Defaults to "IGHV".
#' @param grouping_variable Column used to group data. Defaults to "named_clusters"
#' @param export_pdf Logical. If TRUE, exports the plot to a PDF. Defaults to FALSE
#' @param return_only_data Logical. If TRUE, returns only test results. Defaults to FALSE.
#' @param return_only_sig_data Logical. If TRUE, returns only significant test results. Defaults to FALSE.
#' @param file_name_suffix Character string. Optional suffix to be added to the exported pdf.
#' @param custom_gene_levels Character vector. Optional, vector used to set gene list levels.
#' @param file_type Set type of file/annotation input. Defaults to "cellranger". Options are "IMGT" and "AIRR".
#' @param ylim_max Numeric. Optional upper limit for y-axis. ggplot2 will adjust the y-axis scale to fit the data.
#' @param p.adj.method Method for adjust for multiple tests. Defaults to 'Hochberg'. Alternative is 'FDR'.
#' @param fill_colours A character vector of fill colours. Defaults to greyscale.
#' @param border_colour Colour of bar borders. Defaults to "black".
#' @param ggplot_theme Optional ggplot2 theme to apply. Defaults to "theme_classic".
#' @param plot_width Numeric. Width of exported plot in inches. Used if 'export_pdf = TRUE'. Defaults to 15.
#' @param plot_height Numeric. Height of exported plot in inches. Used if 'export_pdf = TRUE'. Defaults to 5.
#' @param remove_genes_present_in_one_dataset Logical. If TRUE, removes genes not shared by both groups. Defaults to FALSE
#' @param plot_zeroes Logical. If FALSE, genes with 0 counts in one group are not plotted, resulting in a wider bar for the other group.
#' @param plot_subtitle Character string. Plot subtitle. Defaults to "% of max".
#' @param legend_title Character string. Optional. Title of the legend.
#' @param plot_n_sequences Logical. Should the number of sequences used for each group be plotted in the legend? Defaults to FALSE.
#' @param add_pattern Logical. Should the bars be plotted with pattern? Defaults to FALSE. Requires `ggpattern`
#' @param pattern What to base the pattern on. Defaults to the grouping_variable.
#' @param pattern_colour String. Colour of the pattern. Defaults to black.
#' @param pattern_fill String. Fill colour of the pattern. Defaults to black.
#' @param pattern_density Numeric. Density of the pattern. Defaults to 0.1
#' @param add_significancies Whether to calculate and plot test results. Defaults to TRUE.
#' @param plot_only_significant_comparisons Logical. Should only significant test results be plotted? Defaults to TRUE.
#' @param sig_distance_above_bar Numeric. What distance above the tallest bar the stats results should be plotted. Defaults to 1.
#' @param sig_distance_between_significancies Numeric. What distance should be between each test result? Defaults to 0.8.
#' @param sig_line_thickness Numeric. Thickness of the comparison line. Defaults to 1.
#' @param sig_label_distance_from_line Numeric. Distance of asterisk from the comparison line. Defaults to 0.2.
#' @param sig_label_text_size Numeric. Size of asterisks. Defaults to 6.
#' @param return_usage_proportions Logical. Whether to return the proportions of usage per group? Defaults to FALSE.
#'
#' @return Returns a ggplot object, or a data.frame containing test results.
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
#' # TS_plot_IG_usage(input_file = VDJ_data, compare_groups = c("Group1", "Group2", "Group3"), compare_what = "IGHV", export_pdf = TRUE, file_type = "IMGT")
TS_plot_IG_usage <-
  function(
      input_file,
      compare_groups,
      compare_what = "IGHV",
      grouping_variable = "named_clusters",
      export_pdf = F,
      return_only_data = F,
      return_only_sig_data = F,
      file_name_suffix = NULL,
      custom_gene_levels = NULL,
      file_type = "IMGT",
      ylim_max = NULL,
      p.adj.method = "hochberg",
      fill_colours = NULL,
      border_colour = "black",
      ggplot_theme = ggplot2::theme_classic(),
      plot_width = 15,
      plot_height = 5,
      remove_genes_present_in_one_dataset = F,
      plot_zeroes = T,
      plot_subtitle = "% of max",
      custom_plot_title = NULL,
      legend_title = "",
      plot_n_sequences = F,
      add_pattern = F,
      pattern = grouping_variable,
      pattern_colour = "black",
      pattern_fill = "black",
      pattern_density = 0.1,
      add_significancies = T,
      plot_only_significant_comparisons = T,
      sig_distance_above_bar = 1,
      sig_distance_between_significancies = 0.8,
      sig_line_thickness = 1,
      sig_label_distance_from_line = 0.2,
      sig_label_text_size = 6,
      return_usage_proportions = F) {
    options(dplyr.summarise.inform = F) # TMI!

    #######################
    ## Conform file_type ##
    #######################

    # IMGT
    if (file_type == "IMGT") {
      input_file <- TS_format_IMGT_to_cellranger(imgt_data = input_file)
    }

    # AIRR
    if (file_type == "AIRR") {
      stop("AIRR file format is not yet supported")
      # input_file <- TS_format_AIRR_to_cellranger()
    }

    ######################
    ## Set up variables ##
    ######################

    # <--- Check compare_what is supported ---> #
    # Supported comparisons
    supported_comparisons <- c(
      "IGHV",
      "IGHV_FAMILY",
      "IGHV_FAMILY_GROUPED",
      "IGHD",
      "IGHD_FAMILY",
      "IGHJ",
      "IGKV",
      "IGKV_FAMILY",
      "IGKJ",
      "IGLV",
      "IGLV_FAMILY",
      "IGLJ",
      "IGLC"
    )
    compare_what <- toupper(compare_what) # Shift to upper case

    if (!compare_what %in% supported_comparisons) {
      stop(paste0(
        "\nYour choice of compare_what = ",
        compare_what,
        " is not supported. \n Options include: \n'IGHV', ' IGHV_family', 'IGHV_family_grouped', 'IGHJ', \n'IGKV', 'IGKV_family', 'IGKJ', \n'IGLV', 'IGLV_family', 'IGLJ', 'IGLC'\nCapital doesn't matter"
      ))
    }

    # Automate this next section
    # What is required
    #
    # Rename gene_list = [whatever variable should be used]
    # Recognize the pattern,
    #  --> based on this set the correct gene_levels
    #  --> And the correct plot_title
    #
    # <--- IGH ---> #
    # IGHV usage
    if (compare_what == "IGHV") {
      input_file <- input_file %>%
        dplyr::rename(gene_list = v_gene)
      gene_list_levels <- IGHV_levels
      plot_title <- paste0(
        "IGHV usage",
        paste0(compare_groups, collapse = " vs ")
      )
    }

    # IGHV family usage
    if (compare_what == "IGHV_FAMILY") {
      input_file <- input_file %>%
        dplyr::mutate(
          gene_list = stringr::str_sub(string = v_gene, start = 1, end = 5)
        )
      gene_list_levels <- IGHV_family_levels
      plot_title <- paste0(
        "IGHV family usage - ",
        paste0(compare_groups, collapse = " vs ")
      )
    }

    # IGHV family usage - with IGHV6 and 7 grouped
    if (compare_what == "IGHV_FAMILY_GROUPED") {
      input_file <- input_file %>%
        dplyr::mutate(
          gene_list = stringr::str_sub(string = v_gene, start = 1, end = 5),
          gene_list = case_when(
            gene_list == "IGHV6" ~ "IGHV6/7",
            gene_list == "IGHV7" ~ "IGHV6/7",
            .default = gene_list
          )
        )
      gene_list_levels <- IGHV_family_levels_grouped
      plot_title <- paste0(
        "IGHV family usage - ",
        paste0(compare_groups, collapse = " vs ")
      )
    }
    # IGHD usage
    if (compare_what == "IGHD") {
      input_file <- input_file %>%
        dplyr::rename(gene_list = d_gene)
      gene_list_levels <- IGHD_levels
      plot_title <- paste0("IGHD usage")
      if (file_type == "cellranger") {
        message(
          "You are comparing IGHD usage with data from Cell Ranger. This is not advisable to do as they usually map <50% of the genes\nAnyway, continuing..."
        )
      }
    }

    # IGHD family usage
    if (compare_what == "IGHD_FAMILY") {
      input_file <- input_file %>%
        dplyr::mutate(
          gene_list = stringr::str_sub(string = d_gene, start = 1, end = 5)
        )
      gene_list_levels <- IGHD_family_levels
      plot_title <- paste0(
        "IGHD family usage - ",
        paste0(compare_groups, collapse = " vs ")
      )
      if (file_type == "cellranger") {
        message(
          "You are comparing IGHD usage with data from Cell Ranger. This is not advisable to do as they usually map <50% of the genes\nAnyway, continuing..."
        )
      }
    }
    # IGHJ usage
    if (compare_what == "IGHJ") {
      input_file <- input_file %>%
        dplyr::rename(gene_list = j_gene)
      gene_list_levels <- IGHJ_levels
      plot_title <- paste0(
        "IGHJ usage - ",
        paste0(compare_groups, collapse = " vs ")
      )
    }

    # <--- IGK ---> #
    # IGKV usage
    if (compare_what == "IGKV") {
      input_file <- input_file %>%
        dplyr::rename(gene_list = v_gene)
      gene_list_levels <- IGKV_levels
      plot_title <- paste0(
        "IGKV usage - ",
        paste0(compare_groups, collapse = " vs ")
      )
    }

    # IGKV family usage
    if (compare_what == "IGKV_FAMILY") {
      input_file <- input_file %>%
        dplyr::mutate(
          gene_list = stringr::str_sub(string = v_gene, start = 1, end = 5)
        )
      gene_list_levels <- IGKV_family_levels
      plot_title <- paste0(
        "IGKV family usage - ",
        paste0(compare_groups, collapse = " vs ")
      )
    }

    # IGKJ usage
    if (compare_what == "IGKJ") {
      input_file <- input_file %>%
        dplyr::rename(gene_list = j_gene)
      gene_list_levels <- IGKJ_levels
      plot_title <- paste0(
        "IGKJ usage - ",
        paste0(compare_groups, collapse = " vs ")
      )
    }

    # <--- IGL ---> #
    # IGLV usage
    if (compare_what == "IGLV") {
      input_file <- input_file %>%
        dplyr::rename(gene_list = v_gene)
      gene_list_levels <- IGLV_levels
      plot_title <- paste0(
        "IGLV usage - ",
        paste0(compare_groups, collapse = " vs ")
      )
    }

    # IGLV family usage
    if (compare_what == "IGLV_FAMILY") {
      input_file <- input_file %>%
        dplyr::mutate(
          gene_list = stringr::str_sub(string = v_gene, start = 1, end = 5)
        )
      gene_list_levels <- IGLV_family_levels
      plot_title <- paste0(
        "IGLV family usage - ",
        paste0(compare_groups, collapse = " vs ")
      )
    }

    # IGLJ usage
    if (compare_what == "IGLJ") {
      input_file <- input_file %>%
        dplyr::rename(gene_list = j_gene)
      gene_list_levels <- IGLJ_levels
      plot_title <- paste0(
        "IGLJ usage - ",
        paste0(compare_groups, collapse = " vs ")
      )
    }

    # IGLC usage
    if (compare_what == "IGLC") {
      if (file_type == "IMGT") {
        stop(
          "IMGT does not call constant region genes, therefore it is missing in the data"
        )
      }

      input_file <- input_file %>%
        dplyr::rename(gene_list = c_gene)
      gene_list_levels <- IGLC_levels
      plot_title <- paste0(
        "IGLC usage - ",
        paste0(compare_groups, collapse = " vs ")
      )
    }

    # If plot title is supplied
    if (!is.null(custom_plot_title)) {
      plot_title <- custom_plot_title
    }

    # If custom gene_list_levels is supplied
    if (!is.null(custom_gene_levels)) {
      gene_list_levels <- custom_gene_levels
    }

    #############################
    ## Filter input data frame ##
    #############################

    # Retrieve genes present in both groups
    # Only validated for IGHV
    if (remove_genes_present_in_one_dataset) {
      present_genes <- input_file %>%
        dplyr::select(gene_list, dplyr::all_of(grouping_variable)) %>%
        dplyr::filter(.data[[grouping_variable]] %in% compare_groups) %>%
        dplyr::filter(!is.na(gene_list)) %>%
        dplyr::group_by(.data[[grouping_variable]]) %>%
        dplyr::count(gene_list) %>%
        dplyr::group_by(gene_list) %>%
        dplyr::summarise(
          n_groups = dplyr::n_distinct(.data[[grouping_variable]]),
          .groups = "drop"
        ) %>%
        dplyr::filter(n_groups == 2) %>%
        dplyr::pull(gene_list)

      # Filter data
      df_filtered <- input_file %>%
        dplyr::filter(.data[[grouping_variable]] %in% compare_groups) %>%
        dplyr::select(gene_list, dplyr::all_of(grouping_variable)) %>%
        dplyr::mutate(
          gene_list = factor(
            gene_list,
            levels = gene_list_levels[gene_list_levels %in% present_genes]
          ),
          grouping_variable = factor(
            .data[[grouping_variable]],
            levels = compare_groups
          )
        ) %>%
        dplyr::filter(
          !is.na(gene_list) & !is.na(.data[[grouping_variable]])
        ) %>%
        dplyr::group_by(grouping_variable) %>%
        dplyr::count(gene_list)
    } else {
      # Filter data
      df_filtered <- input_file %>%
        dplyr::filter(.data[[grouping_variable]] %in% compare_groups) %>% # Keep only cells from chosen groups
        dplyr::select(gene_list, dplyr::all_of(grouping_variable)) %>% # Remove excess columns
        dplyr::mutate(
          gene_list = factor(gene_list, levels = gene_list_levels),
          grouping_variable = factor(
            .data[[grouping_variable]],
            levels = compare_groups
          )
        ) %>% # Convert to factors and set levels
        dplyr::filter(!is.na(gene_list) & !is.na(.data[[grouping_variable]])) # Remove NA values in genes and groups
    }

    #######################
    ### Calculate stats ###
    #######################

    #
    # Note that you have to discard those genes which have 0 count in both samples
    # As this otherwise interferes with the p adjust!
    #
    # Stats
    if (return_only_data || add_significancies) {
      df_stats <- TS_compute_pairwise_fisher(
        df = df_filtered,
        gene_var = "gene_list",
        group_var = "grouping_variable",
        p_adj_method = p.adj.method
      )
    }
    # Return dataframe if specified
    if (return_only_data) {
      return(df_stats)
    }
    if (return_only_sig_data) {
      df_stats_sig <- df_stats %>%
        dplyr::filter(p_sig != "ns")
      return(df_stats_sig)
    }

    # Otherwise plot the resulting IG-usage chart with stats-results denoted
    if (!return_only_data && !return_only_sig_data) {
      ##################
      ### Plot usage ###
      ##################

      # Code to fill NAs with 0.1, so that lines are printed in the final plot (that are still present)
      if (plot_zeroes) {
        present_genes <- input_file %>%
          dplyr::pull(gene_list) %>%
          unique()
        gene_list_levels_plotting <- gene_list_levels[
          gene_list_levels %in% present_genes
        ]
      } else {
        gene_list_levels_plotting <- gene_list_levels
      }

      # If custom gene_list_levels is supplied
      if (!is.null(custom_gene_levels)) {
        gene_list_levels_plotting <- custom_gene_levels
      }

      if (length(gene_list_levels_plotting) == 0) {
        stop("Error: The gene list for plotting has 0 entries")
      }

      plot_data <- df_filtered %>%
        dplyr::count(gene_list, grouping_variable) %>%
        dplyr::mutate(
          gene_list = factor(gene_list, levels = gene_list_levels_plotting)
        ) %>%
        dplyr::group_by(grouping_variable, gene_list) %>%
        dplyr::summarise(count = n, .groups = "drop") %>%
        dplyr::mutate(count = as.double(count)) %>%
        tidyr::complete(
          gene_list,
          grouping_variable,
          fill = list(count = 0.01)
        ) %>%
        dplyr::group_by(grouping_variable) %>%
        dplyr::mutate(proportion = count / sum(count) * 100)
      if (return_usage_proportions) {
        usage_proportion_data <- plot_data %>%
          tidyr::pivot_wider(
            id_cols = gene_list,
            names_from = grouping_variable,
            values_from = proportion
          )
        return(usage_proportion_data)
      }
      if (!is.null(ylim_max)) {
        if (ylim_max < max(plot_data$proportion)) {
          message(
            "Warning: You have set your ylim_max to ",
            ylim_max,
            ". Increase the value of ylim_max to plot them.
  These values will NOT be plotted: "
          )
          print(plot_data[plot_data$proportion > ylim_max, ])
        }
      }
      if (!add_pattern) {
        # Implement something similar to this?
        # To add number of sequences information
        # ggplot2::scale_fill_manual(values = fill_colours,
        #                            name = "Sample",
        #                            labels = c(paste0(sample1, "\n(n = ", intermediate_export[1, 6], ")"), paste0(sample2, "\n(n = ", intermediate_export[1, 7], ")")))

        plot <- TS_plot_column_chart(
          data = plot_data,
          x = gene_list,
          y = proportion,
          fill = grouping_variable,
          fill_colours = fill_colours,
          name = legend_title, # Change this to set legend title - keep empty for now
          plot_n_sequences = plot_n_sequences,
          ggplot_theme = ggplot_theme,
          ylim_max = ylim_max,
          border_colour = border_colour,
          plot_title = plot_title,
          plot_subtitle = plot_subtitle
        )
      } else {
        plot <- TS_plot_column_pattern_chart(
          data = plot_data,
          x = gene_list,
          y = proportion,
          fill = grouping_variable,
          fill_colours = fill_colours,
          name = legend_title,
          plot_n_sequences = plot_n_sequences,
          ggplot_theme = ggplot_theme,
          ylim_max = ylim_max,
          border_colour = border_colour,
          plot_title = plot_title,
          plot_subtitle = plot_subtitle,
          pattern = pattern,
          pattern_colour = pattern_colour,
          pattern_fill = pattern_fill,
          pattern_density = pattern_density
        )
      }

      ######################
      ## Add stats labels ##
      ######################

      if (add_significancies) {
        if (plot_only_significant_comparisons) {
          df_stats <- df_stats %>%
            dplyr::filter(p_sig != "ns")
        }

        plot <- TS_add_sig_segments_and_labels(
          plot = plot,
          compare_groups = compare_groups,
          gene_list_levels_plotting = gene_list_levels_plotting,
          df_stats = df_stats,
          sig_distance_above_bar = sig_distance_above_bar,
          sig_distance_between_significancies = sig_distance_between_significancies,
          sig_line_thickness = sig_line_thickness,
          sig_label_distance_from_line = sig_label_distance_from_line,
          sig_label_text_size = sig_label_text_size
        )
      }

      # Fix export_pdf so that name includes all groups in compare_groups
      # just nest paste-statements (paste(compare_what, paste(compare_groups, collapse = "_")))
      # ezpz
      if (export_pdf) {
        if (!is.null(file_name_suffix)) {
          export_pdf_filename_suffix <- paste0("_", file_name_suffix)
        } else {
          export_pdf_filename_suffix <- NULL
        }
        compare_groups_filename <- gsub(
          x = paste(compare_groups, collapse = "_"),
          pattern = " ",
          replacement = "_"
        )
        export_pdf_filename <- paste0(
          compare_what,
          "_usage_",
          compare_groups_filename,
          export_pdf_filename_suffix,
          "_",
          file_type,
          ".pdf"
        )

        ggplot2::ggsave(
          filename = export_pdf_filename,
          plot = plot,
          units = "in",
          width = plot_width,
          height = plot_height
        )
        message(paste0(
          "\n  Plot exported to:\n  ",
          getwd(),
          "\n\n  With the filename:\n  ",
          export_pdf_filename
        ))
      } else {}
      return(plot)
      # return(plot_with_sig)
    }
  }
