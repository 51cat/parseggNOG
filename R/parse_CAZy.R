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
