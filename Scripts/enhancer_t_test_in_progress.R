fragments_per_genes <- matchit_Blueprint_Pchic %>% split(., .$Blueprint_gene_names)

test <- fragments_per_genes[c(1:100)]

fragments_per_genes2 <- lapply(fragments_per_genes, function(Gene){
  connected_fragment_bait <- pchic %>% dplyr::filter(., IDbait %in% Gene$ID) %>% dplyr::select(., IDoe) %>% c(.$IDoe)
  connected_fragment_oe <- pchic %>% dplyr::filter(., IDoe %in% Gene$ID) %>% dplyr::select(., IDbait) %>% c(.$IDbait)
  tmp <- c(connected_fragment_bait, connected_fragment_oe) %>% unique(.)
  Gene <- unlist(tmp) %>% unique(.)
})

saveRDS(fragments_per_genes3, file = "GitHub/Epigenomic_integration/fragment_connected_to_genes.csv")

CpGs_per_fragments_EPIC <- matchit_CpGs_Pchic_EPIC %>% split(., .$ID)
  
saveRDS(CpGs_per_fragments_EPIC, file = "GitHub/Epigenomic_integration/CpGs_EPIC_per_fragment.csv")

CpGs_per_fragments_450k <- matchit_CpGs_Pchic_450 %>% split(., .$ID)

saveRDS(CpGs_per_fragments_450k, file = "GitHub/Epigenomic_integration/CpGs_450k_per_fragment.csv")

CpGs_per_fragments_450k_2 <- lapply(CpGs_per_fragments_450k, function(fragment){
  data.frame(fragment$CpG)
})

test <- CpGs_per_fragments_450k_2[c(1:100)]

CpG_in_BMIQ <- rownames(BMIQ_CD34_IDHm_WT)


test2 <- lapply(CpGs_per_fragments_450k_2, function(fragment){
  tmp <- na.omit(BMIQ_CD34_IDHm_WT[fragment$fragment.CpG,])
  if(length(rownames(tmp)) > 1){
    c(colMeans2(as.matrix(na.omit(BMIQ_CD34_IDHm_WT[fragment$fragment.CpG,]))))
  }else{
    c(na.omit(BMIQ_CD34_IDHm_WT[fragment$fragment.CpG,]))
  }
}) %>%
  list.filter(., is.double(.))

A <- Phenotype_BMIQ_CD34_IDHm_WT$Phenotype == "IDHm"
B <- Phenotype_BMIQ_CD34_IDHm_WT$Phenotype == "WT"

test3 <- lapply(test2, function(fragment){
  tmp <- t.test(fragment[A], fragment[B])
  data.frame("p.value" = tmp$p.value, "mean_of_x" = tmp$estimate[1], "mean_of_y" = tmp$estimate[2])
})

T_Test_on_enhancer_TCGA <- test3

saveRDS(T_Test_on_enhancer_TCGA, file = "GitHub/Epigenomic_integration/T_Test_on_enhancer_TCGA.rds")

CpGs_per_fragments_EPIC_2 <- lapply(CpGs_per_fragments_EPIC, function(fragment){
  data.frame(fragment$CpG)
})

CpG_in_BMIQ_whiele <- rownames(BMIQ_Whiele)


test2 <- lapply(CpGs_per_fragments_EPIC_2, function(fragment){
  tmp <- na.omit(BMIQ_Whiele[fragment$fragment.CpG,])
  if(length(rownames(tmp)) > 1){
    c(colMeans2(as.matrix(na.omit(BMIQ_Whiele[fragment$fragment.CpG,]))))
  }else{
    c(na.omit(BMIQ_Whiele[fragment$fragment.CpG,]))
  }
}) %>%
  list.filter(., is.double(.))

A <- Phenotype_BMIQ_CD34_IDHm_WT$Phenotype == "IDHm"
B <- Phenotype_BMIQ_CD34_IDHm_WT$Phenotype == "WT"

test3 <- lapply(test2, function(fragment){
  tmp <- t.test(fragment[A], fragment[B])
  data.frame("p.value" = tmp$p.value, "mean_of_x" = tmp$estimate[1], "mean_of_y" = tmp$estimate[2])
})