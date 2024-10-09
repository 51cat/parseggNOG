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
        TERM2GENE= kegg.df |> dplyr::select(pathway, Gene),
        TERM2NAME = kegg.df |> dplyr::select(pathway, name) |> dplyr::distinct()
      )
      )
}
