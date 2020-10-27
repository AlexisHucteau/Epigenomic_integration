Differential_analysis <- function(Focused_variable, DATA, type_of_data){
  design.pairs <- function(levels) {
    n <- length(levels)
    design <- matrix(0,n,choose(n,2))
    rownames(design) <- levels
    colnames(design) <- 1:choose(n,2)
    k <- 0
    for (i in 1:(n - 1))
      for (j in (i + 1):n) {
        k <- k + 1
        design[i,k] <- 1
        design[j,k] <- -1
        colnames(design)[k] <- paste(levels[i], "-", levels[j],sep = "")
      }
    design
  }
  design <- model.matrix(~0 + Focused_variable)
  if (type_of_data == "gene_expression"){
    DATA <- voom(DATA, design, plot = FALSE)
  }
  contr.matrix <- design.pairs(levels(factor(Focused_variable)))
  colnames(design) <- rownames(contr.matrix)   
  Fit <- lmFit(DATA, design) %>%
    contrasts.fit(., contr.matrix) %>%
    eBayes(., trend = TRUE)
  
  FitList <- list()
  for (i in 1:ncol(contr.matrix)) {
    FitList[[i]] <- topTable(Fit, coef = i, adjust.method = "BH", number = nrow(DATA)) %>%
      mutate(ID = rownames(.))
    
    message(paste0(i, " done"))
    
  }
  names(FitList) <- colnames(contr.matrix)
  return(FitList)
  
}

add_genes_coordinates <- function(vector_of_genes) {
  gene_name_annotation <- getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol', 'entrezgene_id'), filters = 'ensembl_gene_id', values = vector_of_genes$ID, mart = ensembl)
  genes_annoted <- merge(x = vector_of_genes, y = gene_name_annotation, by.x = "ID", by.y = "ensembl_gene_id", all.x = TRUE)
  return(genes_annoted)
}

prepare_pchic <- function(cell_lines = "all", minimum_interaction = 5){
  load("~/PCHIC/pchic.RData")
  if (cell_lines == "all") {
    cell_lines = c("Mon", "Mac0", "Mac1", "Mac2", "Neu", "MK", "EP", "Ery", "FoeT", "nCD4", "tCD4", "aCD4", "naCD4", "nCD8", "tCD8", "nB", "tB")
    }
  pchic <- data.frame(pchic[rowSums(pchic[,cell_lines] >= minimum_interaction) >= 1, 1:10]) %>% na.omit(.)
  colnames(pchic)[c(1:5, 6:10)] <- rep(c("chr", "start", "end", "ID", "Name"), 2)
  return(pchic)
}

Create_pchic_Grange <- function(pchic){
  PCHiC_bed <- unique(rbind(pchic[, c(1:3, 5)], pchic[, c(6:8, 10)]))
  PCHiC_GRange <- GRanges(
    seqnames = PCHiC_bed$chr,
    IRanges(start = PCHiC_bed$start, end = PCHiC_bed$end),
    Gene_Pchic = PCHiC_bed$Name,
    start_fragment = PCHiC_bed$start,
    end_fragment = PCHiC_bed$end
  )
  PCHiC_GRange$ID <- paste(PCHiC_bed$chr, PCHiC_bed$start, sep = "_")
  return(PCHiC_GRange)
}

Gene_name_finding <- function(gene){
  Blueprint_gene_name <- Blueprint$gene_names %>% unique(.)
  Illumina_gene_name <- anno$UCSC_RefGene_Name %>% strsplit(., ";") %>% unique(.)
  gene_referenced <- c(Blueprint_gene_name, Illumina_gene_name) %>% 
    unlist(.) %>% 
    unique(.) %>%
    .[str_detect(., gene)]
  return(gene_referenced)
}

Focus_Gene <- function(gene, DATA_met, Phenotype, matchit) {
  nb_sample <- length(colnames(DATA_met))
  match_hit_CpGs_Blueprint_focused <- matchit %>%
    dplyr::filter(., .$Blueprint_gene_names == gene | .$Illumina_Gene_name == gene)
  anno_focused <- anno %>% 
    dplyr::filter(., .$UCSC_RefGene_Name == gene) %>%
    dplyr::filter(., .$Name %ni% match_hit_CpGs_Blueprint_focused$CpG)
  df_length <- length(anno_focused$Name)
  to_add <- data.frame("Blueprint_gene_names" = rep(gene, df_length),
                       "type" = rep("Illumina_annotation", df_length),
                       "gene_type" = rep("", df_length),
                       "CpG" = anno_focused$Name,
                       "Illumina_Gene_name" = anno_focused$UCSC_RefGene_Name,
                       "position" = anno_focused$UCSC_RefGene_Group,
                       "Island" = anno_focused$Relation_to_UCSC_CpG_Island)
  final_matchit <- rbind(match_hit_CpGs_Blueprint_focused, to_add)
  Specific_CpGs_value <- merge(final_matchit, DATA_met, by.x = "CpG", by.y = 0, all.x = TRUE) 
  Methylation <- dplyr::filter(Specific_CpGs_value, Blueprint_gene_names == gene | Illumina_Gene_name == gene) %>%
    na.omit(.) %>%
    unique(.)
  Beta_value <- data.frame(c(Methylation[,1:nb_sample+7]) %>%
                             unlist(.) %>%
                             na.omit(.),
                           rep(Methylation$Island, nb_sample),
                           rep(Methylation$type, nb_sample),
                           rep(Methylation$position, nb_sample)) %>%
    unique(.)
  colnames(Beta_value) <- c("Beta_values", "cpg_localisation", "type", "Illumina_position_reference")
  Beta_value$Phenotype <- ifelse(str_detect(rownames(Beta_value), "GSM"), "CD34+", "IDHm")

  ggplot(Beta_value, aes(y=Beta_values, x=Phenotype, fill=cpg_localisation, colour= type))+
    geom_boxplot(alpha=0.25)+
    geom_jitter(width=0.25, alpha=0.25)+
    ggtitle(paste0(gene," Promoter Methylation"))+
    xlab("Position of CpGs")+
    ylab("Value of the methylation")+
    scale_shape_manual(values=c(0,1,2,5,6,7))+
    geom_signif(comparisons = list(c("CD34+", "IDHm")), map_signif_level=TRUE)

  ggsave(
    filename = paste0("../Results/", gene, "_", Phenotype, "_Promoter_state.png"),
    plot = last_plot(),
    device = NULL,
    path = NULL,
    scale = 1,
    width = NA,
    height = NA,
    units = c("in", "cm", "mm"),
    dpi = 300,
    limitsize = TRUE
  )
}

Focus_gene_neighborhood <- function(gene, DATA_met, Phenotype, matchit_blueprint, matchit_pchic){
  nb_sample <- length(colnames(DATA_met))
  match_hit_CpGs_Blueprint_focused <- matchit_blueprint %>%
    dplyr::filter(., Blueprint_gene_names == gene | Illumina_Gene_name == gene)
  anno_focused <- anno %>%
    dplyr::filter(., UCSC_RefGene_Name == gene) %>%
    dplyr::filter(., Name %ni% match_hit_CpGs_Blueprint_focused$CpG)
  df_length <- length(anno_focused$Name)
  to_add <- data.frame("Blueprint_gene_names" = rep(gene, df_length),
                       "type" = rep("Illumina_annotation", df_length),
                       "gene_type" = rep("", df_length),
                       "CpG" = anno_focused$Name,
                       "Illumina_Gene_name" = anno_focused$UCSC_RefGene_Name,
                       "position" = anno_focused$UCSC_RefGene_Group,
                       "Island" = anno_focused$Relation_to_UCSC_CpG_Island)
  final_matchit <- rbind(match_hit_CpGs_Blueprint_focused, to_add)
  CpGs_value <- merge(final_matchit, DATA_met, by.x = "CpG", by.y = 0, all.x = TRUE)
  Specific_CpGs_value <- CpGs_value %>%
    dplyr::filter(., type == "P" | position == "TSS1500" | position == "TSS200") %>%
    na.omit(.)
  focused_gene_fragment <- matchit_pchic %>%
    dplyr::filter(., CpG %in% Specific_CpGs_value$CpG)
  neighborhood <- pchic %>%
    dplyr::filter(., IDbait %in% focused_gene_fragment$ID | IDoe %in% focused_gene_fragment$ID) %>%
    .[,c(1:3,11, 5:8,12, 10)]
  colnames(colnames(neighborhood) <- rep(c("chr", "start", "end", "ID", "Gene_name"), 2))
  Network_nodes <- rbind(neighborhood[,c(1:5)], neighborhood[,c(6:10)]) %>%
    merge(., matchit_pchic, by.x = "ID", by.y = "ID") %>%
    dplyr::select(., c("ID", "CpG")) %>%
    merge(., matchit_blueprint, by.x = "CpG", by.y = "CpG") %>%
    dplyr::filter(., Blueprint_gene_names != gene & Illumina_Gene_name != gene) %>%
    merge(., DATA_met, by.x = "CpG", by.y = 0) %>%
    merge(., matchit_CpGs_Pchic)
  Beta_value_network <- data.frame(c(Network_nodes[,1:nb_sample+8]) %>% unlist(.) %>% na.omit(.),
                           rep(Network_nodes$Island, nb_sample),
                           rep(Network_nodes$type, nb_sample),
                           rep(Network_nodes$position, nb_sample),
                           rep(Network_nodes$CD34_enhancer, nb_sample),
                           rep(Network_nodes$AML_enhancer, nb_sample)) %>%
    unique(.)
  colnames(Beta_value_network) <- c("Beta_values", "cpg_localisation", "type", "Illumina_position_reference", "CD34_enhancer", "AML_enhancer")
  Beta_value_network$enhancer <- ifelse(str_detect(rownames(Beta_value_network), "GSM"), Beta_value_network$CD34_enhancer, Beta_value_network$AML_enhancer)
  Beta_value_network$Phenotype <- ifelse(str_detect(rownames(Beta_value_network), "GSM"), "CD34+", "IDHm")

  ggplot(Beta_value_network, aes(y=Beta_values, x=Phenotype, colour = CD34_enhancer))+
    geom_boxplot(alpha=0.25)+
    geom_jitter(width=0.25, alpha=0.25)+
    ggtitle(paste0(gene," Promoter Neighbor Methylation"))+
    xlab("Position of CpGs")+
    ylab("Value of the methylation")+
    scale_shape_manual(values=c(0,1,2,5,6,7))+
    geom_signif(comparisons = list(c("CD34+", "IDHm")), map_signif_level=TRUE)

  ggsave(
    filename = paste0("../Results/", gene, "_", Phenotype, "_Promoter_neighborhood_state.png"),
    plot = last_plot(),
    device = NULL,
    path = NULL,
    scale = 1,
    width = NA,
    height = NA,
    units = c("in", "cm", "mm"),
    dpi = 300,
    limitsize = TRUE
  )
}  

  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  