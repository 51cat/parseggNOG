library(magrittr)

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

get.KEGG.annotation <- function(eggNOG_df,  use = "Pathway") {

  col <- stringr::str_glue("KEGG_{use}")
  cols <- c("Gene","KEGG_ko", col)

  return(
    eggNOG_df |>
    dplyr::filter(.data[[col]] != "-") |>
    dplyr::select(all_of(cols))
    )

}
