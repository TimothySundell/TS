#' Filter list of markers from FindMarkers function, needs to be output into new df
#'
#' @param df The dataframe you wish to filter
#' @param markers_to_remove Vector of markers to remove
#' @export
TS_filter_marker_list <- function(df, markers_to_remove) {

  require(dplyr)

  for (i in 1:length(markers_to_remove)) {
    df <- dplyr::filter(df, gene != markers_to_remove[i])
  }
  return(df)
}
