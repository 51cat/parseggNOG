# TO DO:
# 对于从kegg抓取的ecname需要将已经过时的名字进行修改
#  gotenzymes: 速度过慢,后期移除
parse.eggNOG.EC <- function(eggNOG.file, database = "KEGG"){
  eggno <- load.eggNOG(eggNOG.file)

  anno.EC <- eggno |>

    dplyr::filter(EC != "-") |>
    dplyr::group_by(Gene) |>

    dplyr::summarise(EC=stringr::str_split(EC, ",")) |>
    tidyr::unnest(EC)  |>
    dplyr::ungroup()

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

s <- parse.eggNOG.EC("test_data/out.emapper.annotations")
res <- s$TERM2NAME
