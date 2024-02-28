#' Plots IGHV usage between two groups using Fisher's Exact Test
#'
#' @description
#' Plot the IGHV usage from two groups and perform Fisher's Exact Test to find significant differences.
#' Limited to comparing two groups.
#'
#' @param input_file Takes as input a dataframe from 10X Genomixs Cellranger VDJ pipeline that has been through QC with added sample/grouping information in the column "sample"
#' @param sample1 Grouping variable value for sample 1
#' @param sample2 Grouping variable value for sample 2
#' @param export_pdf Whether to automatically export the plot as PDF in your working directory. Defaults to 'FALSE'
#' @param return_only_data Whether to only return the the calculated frequencies and the test matrices. 'FALSE' returns the plot in standard out. Defaults to 'FALSE'
#' @param ylim_max Defaults to NULL. ggplot2 will adjust the y-axis scale to fit the data. Set this if you want to keep it consistent between plots.
#' @param p.adj.method Method for adjust for multiple tests. Defaults to 'Benjamini-Hochberg'. Alternative is 'FDR'.
#' @importFrom magrittr  %>%
#' @export

TS_IGHV_usage_prop_with_stats <-
  function(input_file,
           sample1,
           sample2,
           export_pdf = F,
           return_only_data = F,
           ylim_max = NULL,
           p.adj.method = "hochberg") {

    # Load required packages #
    library(ggplot2)
    library(ggrepel)
    library(dplyr)
    library(magrittr)
    library(purrr)
    library(stats)

    options(dplyr.summarise.inform = F) # TMI!

    # Set required variable names from user input #
    total_sample1 <- paste0("total_", sample1)
    total_sample2 <- paste0("total_", sample2)

    count_sample1 <- paste0("count_", sample1)
    count_sample2 <- paste0("count_", sample2)

    frequency_sample1 <- paste0("frequency_", sample1)
    frequency_sample2 <- paste0("frequency_", sample2)

    plot_title <- paste0("IGHV-usage - ", sample1, " vs ", sample2)
    plot_subtitle <- paste0("% of max")

    # Ordering of IGHV genes
    v_gene_levels <- c("IGHV(III)-82", "IGHV7-81", "IGHV4-80", "IGHV3-79", "IGHV(II)-78-1", "IGHV5-78", "IGHV7-77",
                       "IGHV(III)-76-1", "IGHV3-76", "IGHV3-75", "IGHV(II)-74-1", "IGHV3-74", "IGHV3-73", "IGHV3-72",
                       "IGHV3-71", "IGHV2-70", "IGHV1-69D", "IGHV1-68D", "IGHV(III)-67-4D", "IGHV(III)-67-3D", "IGHV1-69-2",
                       "IGHV3-69-1", "IGHV2-70D", "IGHV1-69", "IGHV1-68", "IGHV(III)-67-4", "IGHV(III)-67-3", "IGHV(III)-67-2",
                       "IGHV(II)-67-1", "IGHV1-67", "IGHV3-66", "IGHV(II)-65-1", "IGHV3-65", "IGHV3-64", "IGHV3-63",
                       "IGHV(II)-62-1", "IGHV3-62", "IGHV4-61", "IGHV(II)-60-1", "IGHV3-60", "IGHV4-59", "IGHV1-58",
                       "IGHV3-57", "IGHV7-56", "IGHV4-55", "IGHV3-54", "IGHV(II)-53-1", "IGHV3-53", "IGHV3-52", "IGHV(II)-51-2",
                       "IGHV3-51-1", "IGHV5-51", "IGHV3-50", "IGHV(II)-49-1", "IGHV3-49", "IGHV3-48", "IGHV(III)-47-1",
                       "IGHV3-47", "IGHV(II)-46-1", "IGHV1-46", "IGHV1-45", "IGHV(II)-44-2", "IGHV(IV)-44-1", "IGHV(III)-44",
                       "IGHV(II)-43-1", "IGHV3-43", "IGHV3-42", "IGHV3-41", "IGHV(II)-40-1", "IGHV7-40", "IGHV4-39",
                       "IGHV1-38-4", "IGHV(III)-38-1D", "IGHV3-38-3", "IGHV(III)-44D", "IGHV(II)-43-1D", "IGHV3-43D",
                       "IGHV3-42D", "IGHV7-40D", "IGHV4-38-2", "IGHV(III)-38-1", "IGHV3-38", "IGHV3-37", "IGHV3-36",
                       "IGHV3-35", "IGHV7-34-1", "IGHV4-34", "IGHV3-33-2", "IGHV(II)-33-1", "IGHV3-33", "IGHV3-32",
                       "IGHV(II)-31-1", "IGHV4-31", "IGHV3-30-52", "IGHV(II)-30-51", "IGHV3-30-5", "IGHV3-30-42",
                       "IGHV(II)-30-41", "IGHV4-30-4", "IGHV3-30-33", "IGHV(II)-30-32", "IGHV3-30-3", "IGHV3-30-22",
                       "IGHV(II)-30-21", "IGHV4-30-2", "IGHV4-30-1", "IGHV3-30-2", "IGHV(II)-30-1", "IGHV3-30", "IGHV3-29",
                       "IGHV(II)-28-1", "IGHV4-28", "IGHV7-27", "IGHV(II)-26-2", "IGHV(III)-26-1", "IGHV2-26",
                       "IGHV(III)-25-1", "IGHV3-25", "IGHV1-24", "IGHV3-23D", "IGHV(III)-22-2D", "IGHV(II)-22-1D",
                       "IGHV3-23", "IGHV(III)-22-2", "IGHV(II)-22-1", "IGHV3-22", "IGHV3-21", "IGHV(II)-20-1", "IGHV3-20",
                       "IGHV3-19", "IGHV1-18", "IGHV1-17", "IGHV(III)-16-1", "IGHV3-16", "IGHV(II)-15-1", "IGHV3-15",
                       "IGHV1-14", "IGHV(III)-13-1", "IGHV3-13", "IGHV1-12", "IGHV(III)-11-1", "IGHV3-11", "IGHV2-10",
                       "IGHV3-9", "IGHV1-8", "IGHV5-10-1", "IGHV3-64D", "IGHV3-7", "IGHV3-6", "IGHV(III)-5-2",
                       "IGHV(III)-5-1", "IGHV2-5", "IGHV7-4-1", "IGHV4-4", "IGHV1-3", "IGHV(III)-2-1", "IGHV1-2",
                       "IGHV(II)-1-1", "IGHV6-1")

    # Filter data and calculate stats #

    res <- input_file %>%
      dplyr::select(v_gene, sample) %>%
      dplyr::mutate(
        v_gene = factor(v_gene, levels = v_gene_levels),
        sample = factor(sample, levels = c(sample1, sample2))
      ) %>%
      dplyr::group_by(sample) %>%
      dplyr::count(v_gene) %>%
      dplyr::ungroup() %>%
      dplyr::group_by(v_gene, sample) %>%
      dplyr::summarise(count = n) %>%
      dplyr::ungroup() %>%
      dplyr::group_by(sample) %>%
      dplyr::mutate(total = sum(count),
                    frequency = count / sum(count)) %>%
      tidyr::pivot_wider(
        names_from = sample,
        values_from = c(count, frequency, total),
        values_fill = 0
      ) %>%
      dplyr::mutate(
        !!total_sample1 := get(total_sample1) %>% max(),
        !!total_sample2 := get(total_sample2) %>% max()
      ) %>%
      {. ->> intermediate_export} %>%
      dplyr::ungroup() %>%
      dplyr::group_by(v_gene) %>%
      dplyr::mutate(
        cm = list(# Not required for the analyses per se, but nice to have if you want to go back and have a look at a specific matrix.
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
    # Return dataframe if wanted
    if (return_only_data)
    {
      return(res)
    }

    # Otherwise plot the resulting IGHV-usage graph with stats denoted
    else{

      ################################
      ## Plot IGHV usage with stats ##
      ################################

      # Prepare data for plotting
      plot_data <- input_file %>%
        dplyr::mutate(
          v_gene = factor(v_gene, levels = v_gene_levels),
          sample = factor(sample, levels = c(sample1, sample2))
        ) %>%
        dplyr::group_by(sample) %>%
        dplyr::count(v_gene) %>%
        dplyr::ungroup() %>%
        dplyr::filter(!is.na(v_gene)) %>%
        dplyr::group_by(sample, v_gene) %>%
        dplyr::summarise(n) %>%
        dplyr::mutate(proportion = n / sum(n) * 100)

      if(!is.null(ylim_max)){
        if(ylim_max < max(plot_data$proportion)){
          message("Warning: You have set your ylim_max to ", ylim_max, ". Increase the value of ylim_max to plot them.
  These values will not be plotted: ")
          print(plot_data[plot_data$proportion > ylim_max,])
        }
      }
      plot <- ggplot2::ggplot() +

        ggplot2::geom_col(
          data = plot_data,
          mapping = aes(x = v_gene, y = proportion, fill = sample),
          width = 0.8,
          position = position_dodge(),
          color = "black"
        ) +
        ggplot2::scale_fill_manual(
          values = c("white", "black"),
          name = "Sample",
          labels = c(
            paste0(sample1, "\n(n = ", intermediate_export[1,6], ")"),
            paste0(sample2, "\n(n = ", intermediate_export[1,7], ")"))
        ) +
        ggplot2::theme_bw() +
        ggplot2::theme(axis.text.x = element_text(
          angle = 60,
          hjust = 1,
          vjust = 1,
          family = "Arial"
        )) +
        ggplot2::labs(title = plot_title, subtitle = plot_subtitle) +

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
            x = v_gene,
            y = purrr::map2_dbl((get(frequency_sample1) * 100), (get(frequency_sample2) * 100), max),
            label = p.sig
          )
        )
      if (export_pdf) {
        ggplot2::ggsave(
          filename = paste0("IGHV_usage_", sample1, "_vs_", sample2, ".pdf"),
          plot = plot,
          units = "in",
          width = 12,
          height = 5
        )
        cat(paste0("\n  Plot exported to:\n  ", getwd(), "\n\n  With the filename:\n  ", paste0("IGHV_usage_", sample1, "_vs_", sample2, ".pdf")))
      } else{}
      return(plot)
    }

  }
