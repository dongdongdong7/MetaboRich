# Create swiss_GO and hmdb_GO
# data("setsList")
# swiss_protein <- filter_swiss_protein(setsList$swiss_protein, probability = 0.9)
# uniprot_id <- swiss_protein$id
# GO_tibble <- UNIPROT2GO(uniprot_id = uniprot_id)
# swiss_GO_0.9 <- setsList_protein2setsList_GO(swiss_protein)
# swiss_protein <- filter_swiss_protein(setsList$swiss_protein, probability = 0.1)
# uniprot_id <- swiss_protein$id
# GO_tibble <- UNIPROT2GO(uniprot_id = uniprot_id)
# swiss_GO_0.1 <- setsList_protein2setsList_GO(swiss_protein)
# hmdb_GO <- setsList_protein2setsList_GO(setsList_protein = setsList$hmdb_protein)
# hmdb_GO <- hmdb_GO %>%
#   dplyr::filter(!is.na(id))
# predicted_GO <- swiss_GO_0.9
# predcited_GO_0.1 <- swiss_GO_0.1
# predicted_GO <- predicted_GO %>%
#   dplyr::filter(!is.na(id))

data("setsList")
data("example_data")
data("metabolitesList")
data("predicted_GO")
data("idsl_GO")
data("hmdb_GO")
data("all_GO")
data("predicted_GO_0.1")
devtools::document()

# Tidy three GO library to one large GO library
# all_GO_id <- union(idsl_GO$id, hmdb_GO$id)
# all_GO_id <- union(all_GO_id, predicted_GO$id)
# all_GO_records <- rbind(idsl_GO, hmdb_GO, predicted_GO) %>%
#   dplyr::arrange(id)
# all_GO <- all_GO_records[0, ]
# for(i in 1:length(all_GO_id)){
#   print(paste0(i, "/", length(all_GO_id)))
#   GO_id_tmp <- all_GO_id[i]
#
#   all_GO_records_tmp <- all_GO_records[all_GO_records$id == GO_id_tmp, ]
#   name_tmp <- unique(all_GO_records_tmp$name)[1]
#   type_tmp <- unique(all_GO_records_tmp$type)[1]
#   metabolites_vec <- unique(purrr::list_c(all_GO_records_tmp$metabolites))
#   evidence_tmp <- unique(all_GO_records_tmp$evidence)[1]
#
#   all_GO_tmp <- dplyr::tibble(id = GO_id_tmp, name = name_tmp, type = type_tmp,
#                               metabolites = list(metabolites_vec), evidence = evidence_tmp)
#   all_GO <- rbind(all_GO, all_GO_tmp)
# }

example_data <- id_mapping(example_data, from = "kegg_id", "hmdb_id", metabolites_tibble = metabolitesList$hmdb)

enrichmentRes_kstest_swiss1_0.9 <- kstest(input = example_data, enrich_tibble = setsList$swiss_GO_0.9, adjust = "fdr", thread = 8)
enrichmentRes_kstest_swiss1_0.9 <- enrichmentRes_kstest_swiss1_0.9 %>% dplyr::filter(qvalue < 0.05) %>% dplyr::arrange(qvalue)
openxlsx::write.xlsx(enrichmentRes_kstest_swiss1_0.9, file = "D:/fudan/Projects/2024/meaWmrn/Progress/enrichGO/240503/enrichmentRes_kstest_swiss1_0.9.xlsx")
enrichmentRes_kstest_swiss1_0.1 <- kstest(input = example_data, enrich_tibble = setsList$swiss_GO_0.1, adjust = "fdr", thread = 8)
enrichmentRes_kstest_swiss1_0.1 <- enrichmentRes_kstest_swiss1_0.1 %>% dplyr::filter(qvalue < 0.05) %>% dplyr::arrange(qvalue)
openxlsx::write.xlsx(enrichmentRes_kstest_swiss1_0.1, file = "D:/fudan/Projects/2024/meaWmrn/Progress/enrichGO/240503/enrichmentRes_kstest_swiss1_0.1.xlsx")

enrichmentRes_kstest_hmdb1 <- kstest(input = example_data, enrich_tibble = setsList$hmdb_GO, adjust = "fdr", thread = 8)
enrichmentRes_kstest_hmdb1 <- enrichmentRes_kstest_hmdb1 %>% dplyr::filter(qvalue < 0.01) %>% dplyr::arrange(qvalue)
openxlsx::write.xlsx(enrichmentRes_kstest_hmdb1, file = "D:/fudan/Projects/2024/meaWmrn/Progress/enrichGO/240503/enrichmentRes_kstest_hmdb1.xlsx")

swiss_protein <- filter_swiss_protein(setsList$swiss_protein, probability = 0.1)
enrichmentRes_kstest_swiss_protein <- kstest(input = example_data, enrich_tibble = swiss_protein, adjust = "fdr", thread = 8)
enrichmentRes_kstest_swiss_protein <- enrichmentRes_kstest_swiss_protein %>% dplyr::filter(qvalue < 0.05) %>% dplyr::arrange(pvalue)

enrichmentRes_kstest_hmdb_protein <- kstest(input = example_data, enrich_tibble = setsList$hmdb_protein, adjust = "fdr", thread = 8)
enrichmentRes_kstest_hmdb_protein <- enrichmentRes_kstest_hmdb_protein %>% dplyr::filter(pvalue < 0.05) %>% dplyr::arrange(pvalue)

enrichGO_swiss_protein <- clusterProfiler::enrichGO(gene = enrichmentRes_kstest_swiss_protein$id,
                                                    OrgDb = org.Hs.eg.db::org.Hs.eg.db,
                                                    keyType = "UNIPROT", ont = "ALL")
write.csv(enrichGO_swiss_protein, file = "D:/fudan/Projects/2024/meaWmrn/Progress/enrichGO/240503/enrichGO_swiss_protein.csv")
enrichGO_hmdb_protein <- clusterProfiler::enrichGO(gene = enrichmentRes_kstest_hmdb_protein$id,
                                                   OrgDb = org.Hs.eg.db::org.Hs.eg.db,
                                                   keyType = "UNIPROT",ont = "ALL")
write.csv(enrichGO_hmdb_protein, file = "D:/fudan/Projects/2024/meaWmrn/Progress/enrichGO/240503/enrichGO_hmdb_protein.csv")

swiss_GO[1, ]
GO_tibble[which(GO_tibble$GO == "GO:0003331"), ]$UNIPROT
setsList$swiss_protein[setsList$swiss_protein$id == "Q96IY4",]

plot_enrichmentRes(enrichmentRes_kstest, top = 30, plot_type = 2)
enrichmentRes_grsa <- grsa(input = example_data, enrich_tibble = hmdb_GO, adjust = "fdr", thread = 8)
enrichmentRes_grsa <- enrichmentRes_grsa %>% dplyr::filter(pvalue < 0.05) %>% dplyr::arrange(pvalue)
plot_compareRes(list(kstest = enrichmentRes_kstest, grsa = enrichmentRes_grsa), top = 30, plot_type = 2)

enrichmentRes_ora <- ora(example_data, enrich_tibble = all_GO, N_type = "measure", thread = 8)
enrichmentRes_ora <- enrichmentRes_ora %>% dplyr::filter(pvalue < 0.05) %>% dplyr::arrange(pvalue)
plot_enrichmentRes(enrichmentRes_ora, plot_type = 2, top = 30)
