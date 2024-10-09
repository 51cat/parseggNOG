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
