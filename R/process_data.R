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
