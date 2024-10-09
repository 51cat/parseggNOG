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
    ungroup()

    pathway.all <- kegg_list('pathway', species = species)
    colnames(pathway.all) <- c("pathway", "name")

    if (!is.null(organism)) {
      anno.KEGG <- anno.KEGG |>
        dplyr::mutate(pathway = sub("^[a-zA-Z]+", organism, pathway)) |>
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


kegg_rest <- function(rest_url) {
  message('Reading KEGG annotation online: "', rest_url, '"...')
  content <- yulab.utils::yread(rest_url)

  content %<>% strsplit(., "\t") %>% do.call('rbind', .)
  res <- data.frame(from=content[,1],
                    to=content[,2])
  return(res)
}

kegg_list <- function(db, species = NULL) {
  if (db == "pathway") {
    url <- paste("https://rest.kegg.jp/list", db, species, sep="/")
  } else {
    ## module do not need species
    url <- paste("https://rest.kegg.jp/list", db, sep="/")
  }

  kegg_rest(url)
}




