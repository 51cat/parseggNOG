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
      TERM2GENE= kegg.df |> dplyr::select(Module, Gene),
      TERM2NAME = kegg.df |> dplyr::select(Module, name) |> dplyr::distinct()
    )
  )
}


