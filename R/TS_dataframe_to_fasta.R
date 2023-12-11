#' Converts a dataframe and saves it in Fasta format
#'
#' @description
#' Requires a dataframe with two columns:
# - Column 1 == The ID for the cell (e.g. the 10X barcode)
# - Column 2 == The full nucleotide sequence of the contig (concatenated fwr and cdr sequences in nt-form)
#'
#' @param dataframe A 2-column dataframe where:
#' - Column 1 contains cell IDs
#' - Column 2 contains the full nucleotide sequence of the contig (concatenated fwr and cdr sequences in nt-form)
#' @param filename Filename for the output FILE.fasta. Defaults to 'system date and time.fasta'
#' @export
TS_dataframe_to_fasta <- function(dataframe, filename = gsub(pattern = " ", replacement = "_", x = gsub(pattern = ":", replacement = "-", x = paste0(Sys.time(), ".fasta")))){
  colnames(dataframe) <- c("a", "b")
  writeLines(paste(">", dataframe$a, "\n", dataframe$b, sep = ""), filename)
  cat(paste("Your file", filename, "has been saved to ", getwd(), "/", "\n"))
}
