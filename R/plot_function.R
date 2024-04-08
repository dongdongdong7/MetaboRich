#' @title plot_enrichmentRes
#' @description
#' The results of functional set enrichment analysis were shown.
#'
#' @param enrichmentRes An enrichment result tibble.
#' @param top Select the top feature set before the p-value to display.
#' @param plot_type Display diagram types, 1,2,3.
#' @param legend.position The position of the diagram column
#'
#' @return A ggplot object.
#' @export
#'
#' @examples
#' data("metabolitesList");data("setsList");data("example_data")
#' enrichmentRes_ora <- ora(input = example_data, enrich_tibble = setsList$kegg_pathway, enrich_type = "all", N_type = "database", adjust = "fdr", thread = 2)
#' plot_enrichmentRes(enrichmentRes = enrichmentRes_ora, plot_type = 1)
#' plot_enrichmentRes(enrichmentRes = enrichmentRes_ora, plot_type = 2)
#' plot_enrichmentRes(enrichmentRes = enrichmentRes_ora, plot_type = 3)
plot_enrichmentRes <- function(enrichmentRes, top = 10, plot_type = 1, legend.position = "right"){
  enrichmentRes <- enrichmentRes %>%
    dplyr::arrange(pvalue)
  if(top > nrow(enrichmentRes)) top <- nrow(enrichmentRes)
  enrichmentRes <- enrichmentRes[1:top, ]
  enrichmentRes$id <- factor(enrichmentRes$id, levels = enrichmentRes$id)
  enrichmentRes$name <- factor(enrichmentRes$name, levels = enrichmentRes$name)
  if(plot_type == 1){
    p <- ggplot2::ggplot(enrichmentRes) +
      ggplot2::geom_point(ggplot2::aes(x = id, y = -log(pvalue), size = inset, color = direaction)) +
      ggplot2::scale_color_gradient2(low = "blue", mid = "gray", high = "red", midpoint = 0.5, limits = c(0,1)) +
      ggplot2::scale_size(range = c(5, 10)) +
      ggplot2::theme_bw() +
      ggrepel::geom_text_repel(ggplot2::aes(x = id, y = -log(pvalue), label = name), color = "gray20") +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, vjust = 0.6),
            legend.position = legend.position)
  }else if(plot_type == 2){
    enrichmentRes <- enrichmentRes %>% dplyr::mutate(richfc = inset / set)
    p <- ggplot2::ggplot(enrichmentRes) +
      ggplot2::geom_point(ggplot2::aes(x = richfc, y = -log(pvalue), size = inset, color = direaction)) +
      ggplot2::scale_color_gradient2(low = "blue", mid = "gray", high = "red", midpoint = 0.5, limits = c(0,1)) +
      #scale_color_gradientn(colours = (rev(RColorBrewer::brewer.pal(11,"RdBu")))) +
      ggplot2::scale_size(range = c(5, 10)) +
      ggplot2::theme_bw() +
      ggrepel::geom_text_repel(ggplot2::aes(x = richfc, y = -log(pvalue), label = name), color = "gray20") +
      ggplot2::theme(legend.position = legend.position)
  }else if(plot_type == 3){
    enrichmentRes <- enrichmentRes %>% dplyr::mutate(richfc = inset / set)
    enrichmentRes$name <- factor(enrichmentRes$name, levels = rev(levels(enrichmentRes$name)))
    p <- ggplot2::ggplot(enrichmentRes) +
      ggplot2::geom_point(ggplot2::aes(x = richfc, y = name, size = inset, color = pvalue)) +
      ggplot2::scale_color_gradient(low = "red", high = "blue") +
      ggplot2::scale_size(range = c(5, 10)) +
      ggplot2::theme_bw() +
      ggplot2::theme(legend.position = legend.position)
  }
  return(p)
}

#' @title plot_compareRes
#' @description
#' The results of various enrichment algorithms are compared by network diagram.
#'
#'
#' @param enrichmentResList An enrichment result tibble.
#' @param top Select the top feature set before the p-value to display.
#' @param plot_type Display diagram types, 1,2.
#'
#' @return visNetwork or igraph
#' @export
#'
#' @examples
#' data("metabolitesList");data("setsList");data("example_data")
#' enrichmentRes_ora <- ora(input = example_data, enrich_tibble = setsList$kegg_pathway, enrich_type = "all", N_type = "database", adjust = "fdr", thread = 2)
#' enrichmentRes_kstest <- kstest(input = example_data, enrich_tibble = setsList$kegg_pathway, adjust = "fdr", thread = 2)
#' enrichmentRes_grsa <- grsa(input = example_data, enrich_tibble = setsList$kegg_pathway, adjust = "fdr", thread = 2)
#' plot_compareRes(enrichmentResList = list(ora = enrichmentRes_ora, kstest = enrichmentRes_kstest, grsa = enrichmentRes_grsa), plot_type = 1)
plot_compareRes <- function(enrichmentResList, top = 15, plot_type = 1){
  resNum <- length(enrichmentResList)
  for(i in 1:resNum){
    if(nrow(enrichmentResList[[i]]) < top) top_tmp <- length(enrichmentResList[[i]])
    else top_tmp <- top
    enrichmentResList[[i]] <- enrichmentResList[[i]] %>% dplyr::arrange(pvalue)
    enrichmentResList[[i]] <- enrichmentResList[[i]][1:top_tmp,]
  }
  nodes <- dplyr::tibble(id2 = paste0("RES", 1:resNum), name = names(enrichmentResList), type = rep("res", resNum))
  enrichmentResAll <- purrr::list_rbind(enrichmentResList) %>% dplyr::distinct(id, .keep_all = TRUE)
  nodes_tmp <- enrichmentResAll %>% dplyr::select(id, name)
  nodes_tmp$type <- "sets";colnames(nodes_tmp)[1] <- "id2"
  nodes <- rbind(nodes, nodes_tmp)
  nodes <- tibble::rowid_to_column(nodes, var = "id")
  nodes_res <- nodes %>% dplyr::filter(type == "res")
  edges <- dplyr::tibble(start = character(0), end = character(0))
  for(i in 1:nrow(nodes_res)){
    enrichmentRes <- enrichmentResList[[i]]
    edges_tmp <- dplyr::tibble(start = rep(nodes_res[i, ]$id2, nrow(enrichmentRes)), end = enrichmentRes$id)
    edges <- rbind(edges, edges_tmp)
  }
  tmp <- dplyr::left_join(edges, nodes, by = c("start" = "id2")) %>% dplyr::select(start, end, id)
  colnames(tmp)[3] <- "from"
  tmp <- dplyr::left_join(tmp, nodes, by = c("end" = "id2")) %>% dplyr::select(start, end, from, id)
  colnames(tmp)[4] <- "to"
  edges <- tmp %>% dplyr::select(from, to, start, end)
  g <- igraph::graph_from_data_frame(edges, directed = FALSE, vertices = nodes)
  pal <- ggsci::pal_npg(palette = "nrc")(6)
  if(plot_type == 1){
    igraph::V(g)$color <- sapply(igraph::degree(g), function(x) {ifelse(x <= length(pal), pal[x], "gray")})
    igraph::V(g)$color[1:resNum] <- "gray"
    p <- plot(g)
  }
  else if(plot_type == 2){
    nodes$color <- sapply(igraph::degree(g), function(x) {ifelse(x <= length(pal), pal[x], "gray")})
    nodes$color[1:resNum] <- "gray"
    nodes$size <- ifelse(stringr::str_detect(nodes$id2, "RES"), 50, 50)
    nodes$font.size <- ifelse(stringr::str_detect(nodes$id2, "RES"), 30, 30)
    nodes$shape <- ifelse(stringr::str_detect(nodes$id2, "RES"), "diamond", "dot")
    nodes$title <- nodes$name
    nodes$label <- nodes$name
    p <- visNetwork::visNetwork(nodes, edges, height = "1000px", width = "100%") %>%
      visNetwork::visInteraction(hover = TRUE, hoverConnectedEdges = FALSE, selectConnectedEdges = FALSE) %>%
      visNetwork::visIgraphLayout()
  }
  return(p)
}
