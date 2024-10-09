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
