#' @title id_mapping
#' @description
#' To convert ids between databases,
#' you need to select appropriate metabolites_tibble based on requirements,
#' see metabolitesList.
#'
#' @param input Please edit the input as example_data.
#' @param from id from.
#' @param to id to.
#' @param metabolites_tibble See metabolitesList.
#'
#' @return An example_data format tibble.
#' @export
#'
#' @examples
#' data("example_data", package = "MetaboRich")
#' data("metabolitesList", package = "MetaboRich")
#' example_data <- id_mapping(input = example_data, from = "kegg_id", to = "hmdb_id", metabolites_tibble = metabolitesList$hmdb)
id_mapping <- function(input, from = "kegg_id", to = "hmdb_id", metabolites_tibble){
  input <- dplyr::as_tibble(input)
  n_start <- nrow(input)
  input <- dplyr::left_join(input, metabolites_tibble, by = c("id" = from))
  input <- input[, c(to, "logfc", "pvalue", "qvalue")]
  colnames(input)[1] <- "id"
  input <- input %>% dplyr::filter(!is.na(id))
  n_end <- nrow(input)
  mappingRatio <- paste0(round((n_end / n_start) * 100, 2), "%")
  message(paste0("mappingRatio: ", mappingRatio))
  return(input)
}

#' @title ORA
#' @description
#' The classic ORA algorithm uses
#' hypergeometric distribution for
#' pvalue calculation of function set.
#'
#' @param input
#' Please edit the input as example_data.
#' @param enrich_tibble
#' Please edit the enrich_tibble as setsList's function set.
#' @param enrich_type
#' "all": Enrichment analysis was performed using all differential metabolites;
#' "up": Enrichment analysis was performed using upregulated differential metabolites;
#' "down": Enrichment analysis was performed using downregulated differential metabolites.
#' @param N_type
#' "measure": The measured metabolites were used as background;
#' "database": Use the input metabolite function set as background.
#' @param adjust Method of P-value correction.See stats::p.adjust.
#' @param thread The number of parallel threads.
#'
#' @return An enrichmentRes tibble.
#' @export
#'
#' @examples
#' data("metabolitesList");data("setsList");data("example_data")
#' ora(input = example_data, enrich_tibble = setsList$kegg_pathway, enrich_type = "all", N_type = "database", adjust = "fdr", thread = 2)
ora <- function(input, enrich_tibble, enrich_type = "all", N_type = "measure", adjust = "fdr", thread = 1){
  t1 <- Sys.time()
  cal_p_ora <- function(N, n, M, m){
    #(n * M) / N 均值
    cal_P_hy <- function(N, n, M, m){
      (choose(M, m) * choose(N-M, n-m)) / choose(N, n)
    }
    sum(sapply(m:min(n, M), function(x) cal_P_hy(N = N, n = n, M = M, m = x)))
  }
  initialize_input <- function(input){
    input$direction <- "up"
    input$direction[which(input$logfc > 0)] <- "up"
    input$direction[which(input$logfc < 0)] <- "down"
    input$direction[which(input$pvalue > 0.05)] <- "no change"
    input$correlation <- -log(input$pvalue)
    return(input)
  }
  pcutils::dabiao("ORA")
  data <- initialize_input(input)
  message(paste0("enrich type: ", enrich_type))
  if(enrich_type == "all") differential_data <- data %>% dplyr::filter(direction != "no change")
  else if(enrich_type == "up") differential_data <- data %>% dplyr::filter(direction == "up")
  else if(enrich_type == "down") differential_data <- data %>% dplyr::filter(direction == "down")
  n <- nrow(differential_data)
  enrich_tibble <- enrich_tibble[which(purrr::map_int(enrich_tibble$metabolites, function(x) length(x)) != 0), ] # 去除metabolites长度为0的行
  sets_mets <- unique(purrr::list_c(enrich_tibble$metabolites))
  background_mets <- intersect(data$id, sets_mets)
  if(N_type == "measure") N <- length(background_mets)
  else if(N_type == "database") N <- length(sets_mets)
  else stop("N_type is wrong! {measure, database}")
  metabolitesRatio <- paste0(round((length(background_mets) / nrow(data)) * 100, 2), "%")
  message(paste0("metabolitesRatio: ", metabolitesRatio))
  sets_number <- nrow(enrich_tibble)
  loop <- function(i){
    enrich_tibble_tmp <- enrich_tibble[i, ]
    metabolites <- enrich_tibble_tmp$metabolites[[1]]
    if(N_type == "measure") M <- length(intersect(metabolites, background_mets))
    else M <- length(metabolites)
    target_id_in <- differential_data$id[differential_data$id %in% metabolites]
    m <- length(target_id_in)
    pvalue <- cal_p_ora(N, n, M, m)
    data_tmp <- differential_data[differential_data$id %in% target_id_in, ]
    if(enrich_type == "all") direaction <- sum(data_tmp$correlation[data_tmp$direction == "up"]) / sum(data_tmp$correlation)
    else if(enrich_type == "up") direaction <- 1
    else if(enrich_type == "down") direaction <- 0
    enrichment_res_tmp <- dplyr::tibble(id = enrich_tibble_tmp$id, name = enrich_tibble_tmp$name, type = enrich_tibble_tmp$type,
                                 set = M, inset = m, insetIDs = paste0(target_id_in, collapse = ";"),
                                 direaction = direaction, pvalue = pvalue, qvalue = pvalue)
    return(enrichment_res_tmp)
  }
  pb <- utils::txtProgressBar(max = sets_number, style = 3)
  if(thread == 1){
    enrichment_res_list <- lapply(1:sets_number, function(i) {
      utils::setTxtProgressBar(pb, i)
      loop(i)
    })
  }else if(thread > 1){
    cl <- snow::makeCluster(thread)
    doSNOW::registerDoSNOW(cl)
    opts <- list(progress = function(n) utils::setTxtProgressBar(pb,
                                                                 n))
    # enrichment_res_list <- foreach::foreach(i = 1:sets_number,
    #                                         .options.snow = opts,
    #                                         .packages = c("tidyverse")) %dopar% loop(i)
    enrichment_res_list <- foreach::`%dopar%`(foreach::foreach(i = 1:sets_number,
                                                               .options.snow = opts,
                                                               .packages = c("dplyr")),
                                              loop(i))
    snow::stopCluster(cl)
    gc()
  }else stop("thread is wrong!")
  enrichment_res <- purrr::list_rbind(enrichment_res_list)
  enrichment_res$qvalue <- stats::p.adjust(enrichment_res$pvalue, method = adjust)
  t2 <- Sys.time()
  stime <- sprintf("%.3f", t2 - t1)
  message(paste0("\nsets number: ", sets_number, "\n", "time used: ", stime, " ", attr(t2 - t1, "units"), "\n"))
  return(enrichment_res)
}

#' @title KSTESt
#' @description
#' The P-value of the function set of metabolites was calculated using kstest
#'
#' @param input Please edit the input as example_data.
#' @param enrich_tibble Please edit the enrich_tibble as setsList's function set.
#' @param adjust Method of P-value correction.See stats::p.adjust.
#' @param thread The number of parallel threads.
#'
#' @return An enrichmentRes tibble.
#' @export
#'
#' @examples
#' data("example_data", package = "MetaboRich")
#' data("metabolitesList", package = "MetaboRich")
#' data("setsList", package = "MetaboRich")
#' example_data <- id_mapping(example_data, from = "kegg_id", "hmdb_id", metabolites_tibble = metabolitesList$hmdb)
#' enrichmentRes <- kstest(input = example_data, enrich_tibble = setsList$smpdb_pathway, adjust = "fdr", thread = 8)
kstest <- function(input, enrich_tibble, adjust = "fdr", thread = 1){
  t1 <- Sys.time()
  initialize_input <- function(input){
    input$direction <- "up"
    input$direction[which(input$logfc > 0)] <- "up"
    input$direction[which(input$logfc < 0)] <- "down"
    input$direction[which(input$pvalue > 0.05)] <- "no change"
    input$correlation <- -log(input$pvalue)
    return(input)
  }
  pcutils::dabiao("KSTEST")
  data <- initialize_input(input)
  enrich_tibble <- enrich_tibble[which(purrr::map_int(enrich_tibble$metabolites, function(x) length(x)) != 0), ]
  sets_number <- nrow(enrich_tibble)
  sets_mets <- unique(purrr::list_c(enrich_tibble$metabolites))
  background_mets <- intersect(data$id, sets_mets)
  metabolitesRatio <- paste0(round((length(background_mets) / nrow(data)) * 100, 2), "%")
  message(paste0("metabolitesRatio: ", metabolitesRatio))
  loop <- function(i){
    enrich_tibble_tmp <- enrich_tibble[i, ]
    metabolites <- enrich_tibble_tmp$metabolites[[1]]
    target_id_in <- data$id[data$id %in% metabolites]
    data_tmp <- data[data$id %in% target_id_in, ]
    if(length(which(data_tmp$pvalue < 0.05)) > 0 & length(metabolites) > 2 & length(target_id_in) >= 2){
      pvalue_vec <- data_tmp$pvalue
      pvalue <- stats::ks.test(pvalue_vec, "punif", alternative = "greater")$p.value
    }else{
      pvalue <- 1
    }
    direaction <- sum(data_tmp$correlation[data_tmp$direction == "up"]) / sum(data_tmp$correlation[data_tmp$direction != "no change"])
    enrichment_res_tmp <- dplyr::tibble(id = enrich_tibble_tmp$id, name = enrich_tibble_tmp$name, type = enrich_tibble_tmp$type,
                                 set = length(metabolites), inset = length(target_id_in), insetIDs = paste0(target_id_in, collapse = ";"),
                                 direaction = direaction, pvalue = pvalue, qvalue = pvalue)
  }
  pb <- utils::txtProgressBar(max = sets_number, style = 3)
  if(thread == 1){
    enrichment_res_list <- lapply(1:sets_number, function(i) {
      utils::setTxtProgressBar(pb, i)
      loop(i)
    })
  }else if(thread > 1){
    cl <- snow::makeCluster(thread)
    doSNOW::registerDoSNOW(cl)
    opts <- list(progress = function(n) utils::setTxtProgressBar(pb,
                                                                 n))
    # enrichment_res_list <- foreach::foreach(i = 1:sets_number,
    #                                         .options.snow = opts,
    #                                         .packages = c("dplyr")) %dopar% loop(i)
    enrichment_res_list <- foreach::`%dopar%`(foreach::foreach(i = 1:sets_number,
                                                               .options.snow = opts,
                                                               .packages = c("dplyr")),
                                              loop(i))

    snow::stopCluster(cl)
    gc()
  }else stop("thread is wrong!")
  enrichment_res <- purrr::list_rbind(enrichment_res_list)
  enrichment_res$qvalue <- stats::p.adjust(enrichment_res$pvalue, method = adjust)
  t2 <- Sys.time()
  stime <- sprintf("%.3f", t2 - t1)
  message(paste0("\nsets number: ", sets_number, "\n", "time used: ", stime, " ", attr(t2 - t1, "units"), "\n"))
  return(enrichment_res)
}


#' @title GRSA
#' @description
#' The function set enrichment P-value is calculated using GRSA algorithm.
#'
#' @param input Please edit the input as example_data.
#' @param enrich_tibble Please edit the enrich_tibble as setsList's function set.
#' @param adjust Method of P-value correction.See stats::p.adjust.
#' @param thread The number of parallel threads.
#' @param perm The number of resampling times when calculating a random distribution.
#'
#' @return An enrichmentRes tibble.
#' @export
#'
#' @examples
#' data("example_data", package = "MetaboRich")
#' data("metabolitesList", package = "MetaboRich")
#' data("setsList", package = "MetaboRich")
#' example_data <- id_mapping(example_data, from = "kegg_id", "hmdb_id", metabolites_tibble = metabolitesList$hmdb)
#' enrichmentRes <- grsa(input = example_data, enrich_tibble = setsList$smpdb_pathway, adjust = "fdr", thread = 8)
grsa <- function(input, enrich_tibble, adjust = "fdr", thread = 1, perm = 499){
  transfor_data <- function(input){
    ko_pvalue_data <- as.data.frame(input)
    colnames(ko_pvalue_data) <- c("KO_id", "diff_mean", "p.value", "p.adjust")
    rownames(ko_pvalue_data) <- ko_pvalue_data$KO_id
    ko_stat_data <- ReporterScore::pvalue2zs(ko_pvalue_data, mode = "directed")
    return(ko_stat_data)
  }
  transfor_modulelist <- function(enrich_tibble){
    modulelist <- dplyr::tibble(id = enrich_tibble$id, K_num = purrr::map_vec(enrich_tibble$metabolites, function(x) length(x)),
                         KOs = purrr::map_vec(enrich_tibble$metabolites, function(x) paste0(x, collapse = ",")),
                         Description = enrich_tibble$type)
    return(modulelist)
  }
  initialize_input <- function(input){
    input$direction <- "up"
    input$direction[which(input$logfc > 0)] <- "up"
    input$direction[which(input$logfc < 0)] <- "down"
    input$direction[which(input$pvalue > 0.05)] <- "no change"
    input$correlation <- -log(input$pvalue)
    return(input)
  }
  data <-initialize_input(input)
  pcutils::dabiao("GRSA")
  data <- initialize_input(input)
  enrich_tibble <- enrich_tibble[which(purrr::map_int(enrich_tibble$metabolites, function(x) length(x)) != 0), ]
  ko_stat_data <- transfor_data(input)
  modulelist <- transfor_modulelist(enrich_tibble)
  reporter_data <- ReporterScore::get_reporter_score(ko_stat_data, perm = perm,
                                      feature = "compound",
                                      modulelist = modulelist,
                                      p.adjust.method2 = adjust,
                                      thread = thread)
  tmp <- dplyr::left_join(dplyr::as_tibble(reporter_data), enrich_tibble, by = c("ID" = "id"))
  tmp$direaction <- sapply(1:nrow(tmp), function(i){
    metabolites <- tmp[i, ]$metabolites[[1]]
    target_id_in <- data$id[data$id %in% metabolites]
    data_tmp <- data[data$id %in% target_id_in, ]
    direaction <- sum(data_tmp$correlation[data_tmp$direction == "up"]) / sum(data_tmp$correlation[data_tmp$direction != "no change"])
    return(direaction)
  })
  # v <- tmp$ReporterScore
  # direaction <- ifelse(is.na(v), NA,
  #                      ifelse(v > 0, 0.5 + (v / max(v[!is.na(v)], na.rm = TRUE) * 0.5),
  #                             ifelse(v < 0, 0.5 + (v / abs(min(v[!is.na(v)], na.rm = TRUE)) * 0.5), 0.5)))
  enrichment_res <- dplyr::tibble(id = tmp$ID, name = tmp$name, type = tmp$type, set = tmp$K_num, inset = tmp$Exist_K_num,
                           insetIDs = purrr::map_vec(tmp$metabolites, function(x) paste0(input$id[input$id %in% x], collapse = ";")),
                           direaction = tmp$direaction, pvalue = tmp$p.value, qvalue = tmp$p.adjust)
  return(enrichment_res)
}
