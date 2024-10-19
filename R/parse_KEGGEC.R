#' Parse eggNOG EC Number Annotations
#'
#' This function processes an eggNOG output file to parse EC number annotations.
#' It filters out genes without EC numbers, maps EC numbers to their corresponding genes,
#' and returns a list of data frames that map EC numbers to genes and their names.
#' The function can use either the KEGG database or the GotEnzymes database for EC number annoate.
#'
#' @param eggNOG.file The path to the eggNOG output file.
#' @param database A character string specifying the database to use for EC number lookup.
#'   Defaults to "KEGG". Can be set to "GotEnzymes" for slower, real-time lookup.
#'
#' @return A list of two data frames:
#' \item{TERM2GENE}{A data frame mapping EC numbers to genes.}
#' \item{TERM2NAME}{A data frame mapping EC numbers to their names.}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Parse an eggNOG output file for EC number annotations using the KEGG database
#' ec_annotations <- parse.eggNOG.EC("path/to/eggNOG_output.txt")
#'
#' # Parse an eggNOG output file for EC number annotations using the GotEnzymes database (very slowly)
#' ec_annotations_gote <- parse.eggNOG.EC("path/to/eggNOG_output.txt", database = "GotEnzymes")
#'
#' # View the TERM2GENE data frame
#' print(ec_annotations$TERM2GENE)
#'
#' # View the TERM2NAME data frame
#' print(ec_annotations$TERM2NAME)
#' }
parse.eggNOG.EC <- function(eggNOG.file, database = "KEGG"){
  eggno <- load.eggNOG(eggNOG.file)

  anno.EC <- eggno |>

    dplyr::filter(EC != "-") |>
    dplyr::group_by(Gene) |>

    dplyr::summarise(EC=stringr::str_split(EC, ",")) |>
    tidyr::unnest(EC)  |>
    dplyr::ungroup() |>
    dplyr::select(EC, Gene)

  if(database == "GotEnzymes") {
    EC.num <-  anno.EC$EC |> unique()
    EC.name <- lapply(EC.num, get.EC) |> unlist()

    EC.name.df <- dplyr::tibble(
      EC = EC.num,
      name = EC.name
    )

    drop <- EC.name.df |>
      dplyr::filter(name == "Not found from GotEnzymes") |>
      dplyr::pull(EC)
  }else{
    EC.name.df <- kegg_list("EC")
    colnames(EC.name.df) <- c("EC", "name")
    drop <- c()
  }


  TERM2GENE.df <- anno.EC |> dplyr::filter(!EC %in% drop)
  TERM2NAME.df <- EC.name.df |> dplyr::filter(!EC %in% drop)

  return(
    list(
      TERM2GENE = TERM2GENE.df,
      TERM2NAME = TERM2NAME.df
    )
  )
}


get.EC <- function(ec.number) {
  # databse: https://metabolicatlas.org/gotenzymes
  # API: https://metabolicatlas.org/api/v2/#/GotEnzymes/gotEnzymesECInfo
  #

  url <- stringr::str_glue("https://metabolicatlas.org/api/v2/gotenzymes/ecs/{ec.number}")
  ec.json <- RCurl::getURL(url)

  print(url)
  if (ec.json == 'Not Found') {
    return("Not found from GotEnzymes")
  }

  ec.info <- jsonlite::fromJSON(ec.json)$info$name

  if (!stringr::str_detect(ec.info, 'Transferred to', )) {
    return(ec.info)
  }else{
    ec.num.new <-stringr::str_replace_all(ec.info,"Transferred to ", "")

    if (stringr::str_detect(ec.num.new,  "and")) {
      ec.num.new <- stringr::str_split("5.6.2.3 and 5.6.2.4", " and ") |> unlist()
      ec.num.new <- ec.num.new[1]
    }

    warning(stringr::str_glue("{ec.number} --> {ec.info}"))
    return(get.EC(ec.num.new))
  }
}
