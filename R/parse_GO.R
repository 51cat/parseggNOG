#' Parse eggNOG GO Annotations
#'
#' This function processes an eggNOG output file to parse Gene Ontology (GO) annotations.
#' It filters out genes without GO annotations, maps GO terms to their corresponding genes and descriptions,
#' and returns a list of data frames that categorize genes by GO level (Cellular Component, Biological Process, Molecular Function).
#'
#' @param eggNOG.file The path to the eggNOG output file.
#' @param ... Additional arguments passed to load.eggNOG.
#'
#' @return A list of data frames:
#' \item{TERM2GENE.CC}{A data frame mapping GO terms (Cellular Component) to genes.}
#' \item{TERM2GENE.BP}{A data frame mapping GO terms (Biological Process) to genes.}
#' \item{TERM2GENE.MF}{A data frame mapping GO terms (Molecular Function) to genes.}
#' \item{TERM2GENE.ALL}{A data frame mapping all GO terms to genes.}
#' \item{TERM2NAME}{A data frame mapping GO terms to their names.}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Parse an eggNOG output file for GO annotations
#' go_annotations <- parse.eggNOG.GO("path/to/eggNOG_output.txt")
#'
#' # View the TERM2GENE.CC data frame
#' print(go_annotations$TERM2GENE.CC)
#'
#' # View the TERM2NAME data frame
#' print(go_annotations$TERM2NAME)
#' }
parse.eggNOG.GO <- function(eggNOG.file, ...) {

  eggno <-load.eggNOG(eggNOG.file, ...)

  # parse egoNOG
  go.info <- get.go.info()

  anno.GO <- eggno |>
    #select( `#query`, GOs, Description ) %>%
    dplyr::rename(GO = GOs)  |>

    dplyr::filter(GO != "-") |>
    dplyr::group_by(Gene, Description) |>

    dplyr::summarise(GO=stringr::str_split(GO, ",")) |>
    tidyr::unnest(GO)  |>
    dplyr::left_join(go.info, by = "GO")  |>
    dplyr::select(Gene, GO, level, Description)  |>
    dplyr::ungroup()

  TERM2GENE.df <- anno.GO[, c("GO", "Gene" ,"level")] |> tidyr::drop_na()
  TERM2NAME.df <- go.info[, c("GO", "name")] |> dplyr::distinct()

  return(
    list(
      TERM2GENE.CC = TERM2GENE.df |> get.level.df(level = 'CC'),
      TERM2GENE.BP = TERM2GENE.df |> get.level.df(level = 'BP'),
      TERM2GENE.MF = TERM2GENE.df |> get.level.df(level = 'MF'),
      TERM2GENE.ALL = TERM2GENE.df |> get.level.df(level = NULL),
      TERM2NAME = TERM2NAME.df
    )
  )
}

get.level.df <- function(df, level =NULL) {
  if (!is.null(level)) {
    df |>
      dplyr:: filter(level == level) |>
      dplyr::select(GO, Gene)
  }else{
    df |>
      dplyr::select(GO, Gene)
  }
}

get.go.info <- function(){
  goterms  <-  AnnotationDbi::Ontology(GO.db::GOTERM)
  termname  <-  AnnotationDbi::Term(GO.db::GOTERM)
  df.terms <- data.frame(GO = names(goterms),
                         level = goterms,
                         name = termname[names(goterms)]
  )
  return(df.terms)
}
