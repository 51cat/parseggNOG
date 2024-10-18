#' Parse eggNOG CAZy Annotations
#'
#' This function takes an eggNOG output file and parse CAZy annotations,
#' filtering out genes without CAZy annotations, and returning a list of two data frames:
#' one mapping CAZy terms to genes and another mapping CAZy terms to their names.
#'
#' @param eggNOG.file The path to the eggNOG output file.
#' @param ... Additional arguments passed to load.eggNOG.
#'
#' @return A list containing two data frames:
#' \item{TERM2GENE}{A data frame with columns 'CAZy' and 'Gene', where each row represents a unique CAZy term associated with a gene.}
#' \item{TERM2NAME}{A data frame with columns 'CAZy' and 'name', where each row represents a unique CAZy term and its corresponding name.}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Parse an eggNOG output file for CAZy annotations
#' CAZy.res <- parse.eggNOG.CAZy("path/to/eggNOG_output.txt")
#'
#' # View the TERM2GENE data frame
#' print(CAZy.res$TERM2GENE)
#'
#' # View the TERM2NAME data frame
#' print(CAZy.res$TERM2NAME)
#' }

parse.eggNOG.CAZy <- function(eggNOG.file, ...) {
  eggnog <-load.eggNOG(eggNOG.file, ...)

  eggnog <- eggnog |>
    dplyr::filter(CAZy != "-") |>
    dplyr::select(Gene, CAZy) |>
    dplyr::group_by(Gene) |>
    dplyr::summarise(CAZy=stringr::str_split(CAZy, ",")) |>
    tidyr::unnest(CAZy) |>
    dplyr::distinct() |>
    dplyr::ungroup()

  return(
    list(
    TERM2GENE= eggnog |> dplyr::select(CAZy, Gene) |> dplyr::distinct(),

    TERM2NAME = eggnog |>
      dplyr::mutate(name = CAZy) |>
      dplyr::select(CAZy, name) |>
      dplyr::distinct()
  )
  )
}
