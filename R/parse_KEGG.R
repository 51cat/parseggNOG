library(magrittr)

parse.eggNOG.KEGG <- function(eggNOG.file, species = NULL, ...) {

  eggno <-load.eggNOG(eggNOG.file, ...)

  anno.KEGG <- eggno |>

    dplyr::filter(KEGG_Pathway != "-") |>
    dplyr::select(Gene, KEGG_ko, KEGG_Pathway) |>
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
      anno.KEGG <- anno.KEGG |>
        dplyr::mutate(pathway = sub("^[a-zA-Z]+", species, pathway)) |>
        dplyr::filter(pathway %in% pathway.all[["pathway"]])
    }

    anno.KEGG <- anno.KEGG |>
      dplyr::left_join(pathway.all, by = "pathway") |>
    tidyr::drop_na()

    return(
      list(
        TERM2GENE= anno.KEGG |> dplyr::select(pathway, Gene),
        TERM2NAME = anno.KEGG |> dplyr::select(pathway, name)
      )
      )
}
