data("setsList")
swiss_protein <- filter_swiss_protein(setsList$swiss_protein, probability = 0.9)
uniprot_id <- swiss_protein$id
GO_tb <- AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db, keys = uniprot_id, columns = c("GO", "ONTOLOGY"), keytype = "UNIPROT")
GO_tb <- dplyr::as_tibble(GO_tb) %>%
  dplyr::select(GO, UNIPROT, EVIDENCE, ONTOLOGY)

hmdb_GO <- GO_tb %>%
  dplyr::distinct(GO, .keep_all = TRUE) %>%
  dplyr::select(GO, ONTOLOGY, EVIDENCE)
colnames(hmdb_GO) <- c("id", "type", "evidence")

hmdb_GO$metabolites <- lapply(1:nrow(hmdb_GO), function(i) {return(NULL)})

for(i in 1:nrow(GO_tb)){
  print(i)
  GO_tb_tmp <- GO_tb[i, ]
  GO_term <- GO_tb_tmp$GO
  if(is.na(GO_term)) next
  uniprot_tmp <- GO_tb_tmp$UNIPROT
  metabolites_tmp <- swiss_protein[swiss_protein$id == uniprot_tmp, ]$metabolites[[1]]
  hmdb_GO[which(hmdb_GO$id == GO_term), ]$metabolites[[1]] <- c(hmdb_GO[which(hmdb_GO$id == GO_term), ]$metabolites[[1]], metabolites_tmp)
}

goterms <- AnnotationDbi::Term(GO.db::GOTERM)
goterms_tibble <- dplyr::tibble(term = names(goterms), description = as.vector(goterms))
hmdb_GO <- dplyr::left_join(hmdb_GO, goterms_tibble, by = c("id" = "term"))
hmdb_GO <- hmdb_GO %>%
  dplyr::select(id, description, type, metabolites, evidence)
colnames(hmdb_GO)[2] <- "name"

data("example_data")
data("metabolitesList")
example_data <- id_mapping(example_data, from = "kegg_id", "hmdb_id", metabolites_tibble = metabolitesList$hmdb)
enrichmentRes_kstest <- kstest(input = example_data, enrich_tibble = hmdb_GO, adjust = "fdr", thread = 8)
enrichmentRes_kstest <- enrichmentRes_kstest %>% dplyr::filter(pvalue < 0.05) %>% dplyr::arrange(pvalue)
enrichmentRes_grsa <- grsa(input = example_data, enrich_tibble = hmdb_GO, adjust = "fdr", thread = 8)
enrichmentRes_grsa <- enrichmentRes_grsa %>% dplyr::filter(pvalue < 0.05) %>% dplyr::arrange(pvalue)
plot_compareRes(list(kstest = enrichmentRes_kstest, grsa = enrichmentRes_grsa), top = 30, plot_type = 2)
