#' RAM usage of objects.
#'
#' Lists RAM usage of n objects in the global environment.
#'
#' @param number_of_objects Numeric. Number of objects to list. Defaults to 5.
#'
#' @return List of the `number_of_objects` objects that uses the most RAM memory.
#'
#' @export
ram_usage <- function(number_of_objects = 5){

  # Retrieve ram usage for each object and slice top n objects
  top_n_ram_usage <- sapply(X = ls(envir = .GlobalEnv), FUN = function(x){object.size(get(x))}) %>%
    as.data.frame() %>%
    rownames_to_column(var = "object") %>%
    magrittr::set_colnames(c("object", "sizeMB")) %>%
    dplyr::mutate(sizeMB = round(sizeMB / 1e6, 1), # show size in megabytes
                  sizeGB = round(sizeMB / 1e3, 1)) %>% # show size in gigabytes
    dplyr::slice_max(order_by = sizeMB, n = number_of_objects)

  # Calculate total RAM usage for all objects
  total_size_all_objects <- sapply(X = ls(envir = .GlobalEnv), FUN = function(x){object.size(get(x))}) %>%
    as.data.frame() %>%
    rownames_to_column(var = "object") %>%
    magrittr::set_colnames(c("object", "sizeMB")) %>%
    dplyr::mutate(sizeMB = round(sizeMB / 1e6, 1), # show size in megabytes
                  sizeGB = round(sizeMB / 1e3, 1)) %>% # show size in gigabytes
    summarise(object = "Total", sizeMB = sum(sizeMB), sizeGB = sum(sizeGB))

  # Create divider data.frame
  divider <- data.frame(object = NA, sizeMB = "------", sizeGB = "----")

  # Combine dataframes
  output <- rbind(top_n_ram_usage, divider, total_size_all_objects)

  # Print output
  print(output, row.names = F, na.print = "----------------------------")

}
