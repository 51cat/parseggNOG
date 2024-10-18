#' Load and Process eggNOG Output
#'
#' This function loads an eggNOG output file and processes it by renaming the first column to 'Gene'.
#' It uses the readr package to read the delimited file and dplyr to rename the columns.
#'
#' @param file The path to the eggNOG output file.
#' @param col_names Logical indicating whether the file has column names. Defaults to TRUE.
#' @param delim The delimiter used in the file. Defaults to tab ('\t').
#' @param comment The character used to indicate comments in the file. Defaults to '##'.
#' @param ... Additional arguments passed to readr::read_delim.
#'
#' @return A data frame with the first column renamed to 'Gene'.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Load an eggNOG output file
#' eggno <- load.eggNOG("path/to/eggNOG_output.txt")
#'
#' # View the first few rows of the data frame
#' head(eggno)
#' }
load.eggNOG <- function(file,
                        col_names = T,
                        delim = "\t",
                        comment = "##",
                        ...){
  eggno <- readr::read_delim(file,
                             col_names = T, delim = "\t", comment = "##", show_col_types = FALSE, ...)

  return(eggno |>
           dplyr::rename(Gene = `#query`) )
}
