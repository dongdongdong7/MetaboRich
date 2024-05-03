#' @title filter_swiss_protein
#' @description
#' Filter swiss_protein based on probability
#'
#' @param swiss_protein swiss_protein tibble.
#' @param probability Probability threshold.
#'
#' @return A new swiss_protein tibble.
#' @export
#'
#' @examples
#' swiss_protein <- filter_swiss_protein(setsList$swiss_protein, probability = 0.9)
filter_swiss_protein <- function(swiss_protein, probability = 0.1){
  tmpList <- lapply(1:nrow(swiss_protein), function(i) {
    metabolites_tmp <- swiss_protein[i, ]$metabolites[[1]]
    probability_tmp <- swiss_protein[i, ]$probability[[1]]
    idx <- which(probability_tmp > probability)
    if(length(idx) == 0) return(NULL)
    metabolites_tmp <- metabolites_tmp[idx]
    probability_tmp <- probability_tmp[idx]
    tmp <- swiss_protein[i, ]
    tmp$metabolites[[1]] <- metabolites_tmp
    tmp$probability[[1]] <- probability_tmp
    return(tmp)
  })
  delete_vec <- sapply(tmpList, function(x) is.null(x))
  tmpList <- tmpList[!delete_vec]
  swiss_protein <- purrr::list_rbind(tmpList)
  return(swiss_protein)
}

#' @title UNIPROT2GO
#' @description
#' Enter a vector of uniprot and return these uniprot-related GO terms.
#'
#' @param uniprot_id A character of uniprot.
#'
#' @return A tibble.
#' @export
#'
#' @examples
#' swiss_protein <- filter_swiss_protein(setsList$swiss_protein, probability = 0.9)
#' uniprot_id <- swiss_protein$id
#' GO_tibble <- UNIPROT2GO(uniprot_id)
UNIPROT2GO <- function(uniprot_id){
  GO_tb <- AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db, keys = uniprot_id, columns = c("GO", "ONTOLOGY"), keytype = "UNIPROT")
  GO_tb <- dplyr::as_tibble(GO_tb) %>%
    dplyr::select(GO, UNIPROT, EVIDENCE, ONTOLOGY)
  return(GO_tb)
}

#' @title setsList_protein2setsList_GO
#' @description
#' Convert a setsList's protein-related metabolite set to a GO metabolite set.
#'
#' @param setsList_protein A setsList's protein-related metabolite set.
#'
#' @return A GO metabolite set
#' @export
#'
#' @examples
#' swiss_protein <- filter_swiss_protein(setsList$swiss_protein, probability = 0.9)
#' swiss_GO <- setsList_protein2setsList_GO(swiss_protein)
#' hmdb_GO <- setsList_protein2setsList_GO(setsList_protein = setsList$hmdb_protein)
setsList_protein2setsList_GO <- function(setsList_protein){
  GO_tb <- UNIPROT2GO(setsList_protein$id)

  setsList_GO <- GO_tb %>%
    dplyr::distinct(GO, .keep_all = TRUE) %>%
    dplyr::select(GO, ONTOLOGY, EVIDENCE)
  colnames(setsList_GO) <- c("id", "type", "evidence")

  setsList_GO$metabolites <- lapply(1:nrow(setsList_GO), function(i) {return(NULL)})
  pb <- utils::txtProgressBar(max = nrow(GO_tb), style = 3)
  for(i in 1:nrow(GO_tb)){
    utils::setTxtProgressBar(pb, i)
    GO_tb_tmp <- GO_tb[i, ]
    GO_term <- GO_tb_tmp$GO
    if(is.na(GO_term)) next
    uniprot_tmp <- GO_tb_tmp$UNIPROT
    metabolites_tmp <- setsList_protein[setsList_protein$id == uniprot_tmp, ]$metabolites[[1]]
    setsList_GO[which(setsList_GO$id == GO_term), ]$metabolites[[1]] <- unique(c(setsList_GO[which(setsList_GO$id == GO_term), ]$metabolites[[1]], metabolites_tmp))
  }
  goterms <- AnnotationDbi::Term(GO.db::GOTERM)
  goterms_tibble <- dplyr::tibble(term = names(goterms), description = as.vector(goterms))
  setsList_GO <- dplyr::left_join(setsList_GO, goterms_tibble, by = c("id" = "term"))
  setsList_GO <- setsList_GO %>%
    dplyr::select(id, description, type, metabolites, evidence)
  colnames(setsList_GO)[2] <- "name"
  return(setsList_GO)
}
