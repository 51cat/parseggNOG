#' Parse eggNOG KEGG Pathway Annotations
#'
#' This function processes an eggNOG output file to parse KEGG pathway annotations.
#' It filters out genes without KEGG pathway annotations, maps KEGG pathway identifiers to their corresponding genes,
#' and returns a list of data frames that map KEGG pathways to genes and their names.
#' Optionally, it can filter pathways by a specified species.
#'
#' @param eggNOG.file The path to the eggNOG output file.
#' @param species A character string specifying the species for pathway filtering.
#'   Defaults to NULL, which means no species-specific filtering is applied.supported organism listed in 'https://www.genome.jp/kegg/catalog/org_list.html'
#' @param ... Additional arguments passed to load.eggNOG.
#'
#' @return A list of two data frames:
#' \item{TERM2GENE}{A data frame mapping KEGG pathway identifiers to genes.}
#' \item{TERM2NAME}{A data frame mapping KEGG pathway identifiers to their names.}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Parse an eggNOG output file for KEGG pathway annotations without species filtering
#' keggp_annotations <- parse.eggNOG.KEGGP("path/to/eggNOG_output.txt")
#'
#' # Parse an eggNOG output file for KEGG pathway annotations with species filtering
#' keggp_annotations_species <- parse.eggNOG.KEGGP("path/to/eggNOG_output.txt", species = "hsa")
#'
#' # View the TERM2GENE data frame
#' print(keggp_annotations$TERM2GENE)
#'
#' # View the TERM2NAME data frame
#' print(keggp_annotations$TERM2NAME)
#' }

parse.eggNOG.KEGGP <- function(eggNOG.file, species = NULL, ...) {

  eggno <-load.eggNOG(eggNOG.file, ...)
  kegg.df <- get.KEGG.annotation(eggno, use = "Pathway")

  kegg.df <- kegg.df |>
    dplyr::mutate(KEGG_ko = gsub("ko:", "", KEGG_ko),
                  KEGG_Pathway = gsub("ko", "map", KEGG_Pathway)) |>
    dplyr::group_by(Gene) |>
    dplyr::summarise(KO=stringr::str_split(KEGG_ko, ","),
                     pathway=stringr::str_split(KEGG_Pathway, ",")) |>
    tidyr::unnest(KO) |>
    tidyr::unnest(pathway) |>
    dplyr::distinct() |>
    dplyr::ungroup()

    pathway.all <- kegg_list('pathway', species = species)
    colnames(pathway.all) <- c("pathway", "name")

    if (!is.null(species)) {
      kegg.df <- kegg.df |>
        dplyr::mutate(pathway = sub("^[a-zA-Z]+", species, pathway)) |>
        dplyr::filter(pathway %in% pathway.all[["pathway"]])
    }

    kegg.df <- kegg.df |>
      dplyr::left_join(pathway.all, by = "pathway") |>
    tidyr::drop_na()

    return(
      list(
        TERM2GENE= kegg.df |> dplyr::select(pathway, Gene)  |> dplyr::distinct(),
        TERM2NAME = kegg.df |> dplyr::select(pathway, name) |> dplyr::distinct()
      )
      )
}
