#' Parse eggNOG KEGG Reaction Annotations
#'
#' This function processes an eggNOG output file to parse KEGG reaction annotations.
#' It filters out genes without KEGG reaction annotations, maps KEGG reaction identifiers to their corresponding genes,
#' and returns a list of data frames that map KEGG reactions to genes and their names.
#' Optionally, it can filter reactions by a specified species.
#'
#' @param eggNOG.file The path to the eggNOG output file.
#' @param species A character string specifying the species for reaction filtering.
#'   Defaults to NULL, which means no species-specific filtering is applied.Supported organism listed in 'https://www.genome.jp/kegg/catalog/org_list.html'
#' @param ... Additional arguments passed to load.eggNOG.
#'
#' @return A list of two data frames:
#' \item{TERM2GENE}{A data frame mapping KEGG reaction identifiers to genes.}
#' \item{TERM2NAME}{A data frame mapping KEGG reaction identifiers to their names.}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Parse an eggNOG output file for KEGG reaction annotations without species filtering
#' keggr_annotations <- parse.eggNOG.KEGGR("path/to/eggNOG_output.txt")
#'
#' # Parse an eggNOG output file for KEGG reaction annotations with species filtering
#' keggr_annotations_species <- parse.eggNOG.KEGGR("path/to/eggNOG_output.txt", species = "hsa")
#'
#' # View the TERM2GENE data frame
#' print(keggr_annotations$TERM2GENE)
#'
#' # View the TERM2NAME data frame
#' print(keggr_annotations$TERM2NAME)
#' }
parse.eggNOG.KEGGR <- function(eggNOG.file, species = NULL, ...) {

  eggno <-load.eggNOG(eggNOG.file, ...)
  kegg.df <- get.KEGG.annotation(eggno, use = "Reaction")

  kegg.df <- kegg.df |>
    dplyr::mutate(KEGG_ko = gsub("ko:", "", KEGG_ko)) |>
    dplyr::group_by(Gene) |>
    dplyr::summarise(KO=stringr::str_split(KEGG_ko, ","),
                     Reaction=stringr::str_split(KEGG_Reaction, ",")) |>
    tidyr::unnest(KO) |>
    tidyr::unnest(Reaction) |>
    dplyr::distinct() |>
    dplyr::ungroup()

  reaction.all <- kegg_list('Reaction')
  colnames(reaction.all) <- c("Reaction", "name")


  kegg.df <- kegg.df |>
    dplyr::left_join(reaction.all, by = "Reaction") |>
    tidyr::drop_na()

  return(
    list(
      TERM2GENE= kegg.df |> dplyr::select(Reaction, Gene)  |> dplyr::distinct(),
      TERM2NAME = kegg.df |> dplyr::select(Reaction, name) |> dplyr::distinct()
    )
  )
}
