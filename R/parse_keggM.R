#' Parse eggNOG KEGG Module Annotations
#'
#' This function processes an eggNOG output file to parse KEGG Module annotations.
#' It filters out genes without KEGG module annotations, maps KEGG module identifiers to their corresponding genes,
#' and returns a list of data frames that map KEGG modules to genes and their names.
#'
#' @param eggNOG.file The path to the eggNOG output file.
#' @param ... Additional arguments passed to load.eggNOG.
#'
#' @return A list of two data frames:
#' \item{TERM2GENE}{A data frame mapping KEGG module identifiers to genes.}
#' \item{TERM2NAME}{A data frame mapping KEGG module identifiers to their names.}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Parse an eggNOG output file for KEGG module annotations
#' keggm_annotations <- parse.eggNOG.KEGGM("path/to/eggNOG_output.txt")
#'
#' # View the TERM2GENE data frame
#' print(keggm_annotations$TERM2GENE)
#'
#' # View the TERM2NAME data frame
#' print(keggm_annotations$TERM2NAME)
#' }
parse.eggNOG.KEGGM <- function(eggNOG.file , ...) {

  eggno <-load.eggNOG(eggNOG.file, ...)
  kegg.df <- get.KEGG.annotation(eggno, use = "Module")

  kegg.df <- kegg.df |>
    dplyr::mutate(KEGG_ko = gsub("ko:", "", KEGG_ko)) |>
    dplyr::group_by(Gene) |>
    dplyr::summarise(KO=stringr::str_split(KEGG_ko, ","),
                     Module=stringr::str_split(KEGG_Module, ",")) |>
    tidyr::unnest(KO) |>
    tidyr::unnest(Module) |>
    dplyr::distinct() |>
    dplyr::ungroup()

  moduel.all <- kegg_list('Module')
  colnames(moduel.all) <- c("Module", "name")


  kegg.df <- kegg.df |>
    dplyr::left_join(moduel.all, by = "Module") |>
    tidyr::drop_na()

  return(
    list(
      TERM2GENE= kegg.df |> dplyr::select(Module, Gene)  |> dplyr::distinct(),
      TERM2NAME = kegg.df |> dplyr::select(Module, name) |> dplyr::distinct()
    )
  )
}


