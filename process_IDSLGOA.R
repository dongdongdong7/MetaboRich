data("metabolitesList", package = "MetaboRich")
chemical_details <- readr::read_tsv("D:/fudan/Projects/2024/meaWmrn/Progress/database/IDSL.GOA/library/chemical_details.txt",
                                    col_names = FALSE)
colnames(chemical_details) <- c("inchikey", "name", "short_inchikey")
hmdb_metabolites <- metabolitesList$hmdb %>%
  dplyr::select(hmdb_id, name, inchikey)
hmdb_metabolites$my_inchikey <- sapply(hmdb_metabolites$inchikey, function(x) {
  tmp <- strsplit(x, "-")[[1]]
  return(paste0(tmp[1], "-", tmp[2]))
})
hmdb_metabolites$short_inchikey <- sapply(hmdb_metabolites$inchikey, function(x) {
  tmp <- strsplit(x, "-")[[1]]
  return(tmp[1])
})
chemical_details$my_inchikey <- sapply(chemical_details$inchikey, function(x) {
  tmp <- strsplit(x, "-")[[1]]
  return(paste0(tmp[1], "-", tmp[2]))
})
chemical_details <- chemical_details %>%
  dplyr::distinct(short_inchikey, .keep_all = TRUE)

overlap_chemical <- dplyr::left_join(chemical_details, hmdb_metabolites, by = c("short_inchikey" = "short_inchikey"))
overlap_chemical <- overlap_chemical %>%
  dplyr::select(hmdb_id, name.y, inchikey.x, inchikey.y, my_inchikey.x, my_inchikey.y, short_inchikey) %>%
  dplyr::filter(!is.na(hmdb_id))
overlap_chemical_single <- overlap_chemical[!(duplicated(overlap_chemical$short_inchikey, fromLast = TRUE) | duplicated(overlap_chemical$short_inchikey, fromLast = FALSE)),]
overlap_chemical_multi <- overlap_chemical[duplicated(overlap_chemical$short_inchikey, fromLast = TRUE) | duplicated(overlap_chemical$short_inchikey, fromLast = FALSE),]
overlap_chemical_multi <- overlap_chemical_multi %>%
  dplyr::filter(my_inchikey.x == my_inchikey.y)
overlap_chemical_multi <- overlap_chemical_multi %>%
  dplyr::distinct(short_inchikey, .keep_all = TRUE)
overlap_chemical <- rbind(overlap_chemical_single, overlap_chemical_multi)

go2inchikey_link <- readr::read_tsv("D:/fudan/Projects/2024/meaWmrn/Progress/database/IDSL.GOA/library/go_to_inchikey14_links.txt",
                                    col_names = FALSE)
colnames(go2inchikey_link) <- c("GO", "short_inchikey")
length(unique(go2inchikey_link$short_inchikey)) # 1856

length(which(unique(go2inchikey_link$short_inchikey) %in% overlap_chemical$short_inchikey)) # 1368

go2inchikey_link <- dplyr::left_join(go2inchikey_link, overlap_chemical, by = c("short_inchikey" = "short_inchikey"))
go2inchikey_link <- go2inchikey_link %>%
  dplyr::select(GO, hmdb_id) %>%
  dplyr::filter(!is.na(hmdb_id))

idsl_GO <- dplyr::tibble(id = unique(go2inchikey_link$GO))
idsl_GO$metabolites <- lapply(idsl_GO$id, function(x) {return(NULL)})
for(i in 1:nrow(idsl_GO)){
  print(i)
  id_tmp <- idsl_GO[i, ]$id
  idsl_GO[i, ]$metabolites[[1]] <- go2inchikey_link[go2inchikey_link$GO == id_tmp, ]$hmdb_id
}

goterms <- AnnotationDbi::Term(GO.db::GOTERM)
goterms_tibble <- dplyr::tibble(term = names(goterms), description = as.vector(goterms))
GO_all <- AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db,
                      keys = AnnotationDbi::keys(org.Hs.eg.db::org.Hs.eg.db, kyetype = "GO"),
                      columns = c("GO", "ONTOLOGY"))
GO_all <- dplyr::as_tibble(GO_all) %>%
  dplyr::select(GO, ONTOLOGY, EVIDENCE) %>%
  dplyr::filter(!is.na(GO))
GO_all <- dplyr::left_join(GO_all, goterms_tibble, by = c("GO" = "term")) %>%
  dplyr::distinct(GO, .keep_all = TRUE)

idsl_GO <- dplyr::left_join(idsl_GO, GO_all, by = c("id" = "GO")) %>%
  dplyr::select(id, description, ONTOLOGY, metabolites, EVIDENCE)
idsl_GO1 <- idsl_GO %>% dplyr::filter(!is.na(description))
idsl_GO2 <- idsl_GO %>% dplyr::filter(is.na(description))
idsl_GO2 <- dplyr::left_join(idsl_GO2, goterms_tibble, by = c("id" = "term")) %>%
  dplyr::select(id, description.y, ONTOLOGY, metabolites, EVIDENCE)
colnames(idsl_GO2)[2] <- "description"
idsl_GO <- rbind(idsl_GO1, idsl_GO2)
colnames(idsl_GO) <- c("id", "name", "type", "metabolites", "evidence")
idsl_GO <- idsl_GO %>% dplyr::arrange(id)

#usethis::use_data(idsl_GO)
