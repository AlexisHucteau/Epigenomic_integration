library(GenomicRanges)
library(ggplot2)

"%ni%" <- Negate("%in%")


BMIQ <- read.csv("~/GitHub/TCGA_Connections/BMIQ_met.csv")
rownames(BMIQ) <- BMIQ$X

Phenotype <- read.csv("~/GitHub/TCGA_Connections/Phenotype_met.csv")
rownames(Phenotype) <- Phenotype$X %>% str_replace_all(., "-", ".")
Phenotype_IDH <- Phenotype %>% 
  dplyr::filter(., DNMT3A == "DNMT3AWT", WT1 == "WT1WT", FLT3 == "FLT3AWT", TET2 == "TET2WT") 


BMIQ_focused <- BMIQ[,colnames(BMIQ) %in% rownames(Phenotype_IDH)]

Blueprint <- read.csv("~/GitHub/Epigenomic_integration/BLUEPRINT_fragments_good.tsv", sep = "\t") %>%
  dplyr::select(., "gene_names":"type", "gene_type") %>%
  separate_rows(., gene_names, sep = " ") %>%
  separate_rows(., gene_type, sep = " ") %>%
  unique(.)

anno <- read.csv("~/Illumina_Manifest/HumanMethylation450_15017482_v1-2.csv", skip = 7) %>% 
  dplyr::select(., "Name", "CHR", "MAPINFO", "UCSC_RefGene_Name", "UCSC_RefGene_Group", "Relation_to_UCSC_CpG_Island") %>%
  dplyr::filter(., CHR != "")

Gene_name_finding <- function(gene){
  Blueprint_gene_name <- Blueprint$gene_names %>% unique(.)
  Illumina_gene_name <- anno$UCSC_RefGene_Name %>% strsplit(., ";") %>% unique(.)
  gene_referenced <- c(Blueprint_gene_name, Illumina_gene_name) %>% 
    unlist(.) %>% 
    unique(.) %>%
    .[str_detect(., gene)]
  return(gene_referenced)
}

Focus_Gene <- function(gene) {
  nb_sample <- length(colnames(BMIQ_focused))
  Blueprint_focused <- dplyr::filter(Blueprint, str_detect(Blueprint$gene_names, gene))
  anno_focused <- anno %>% separate_rows(., UCSC_RefGene_Name, sep = ";") %>%
    dplyr::filter(., str_detect(.$UCSC_RefGene_Name, gene))
  Blueprint_Granges <- GRanges(
    seqnames = Blueprint_focused$chr,
    ranges = IRanges(start = Blueprint_focused$start, end = Blueprint_focused$end),
    Blueprint_gene_names = Blueprint_focused$gene_names,
    type = Blueprint_focused$type,
    gene_type = Blueprint_focused$gene_type
  )
  CpGs_Granges <- GRanges(
    seqnames = anno$CHR,
    ranges = IRanges(anno$MAPINFO, anno$MAPINFO +1),
    CpG = anno$Name,
    Illumina_Gene_name = anno$UCSC_RefGene_Name,
    position = anno$UCSC_RefGene_Group,
    Island = anno$Relation_to_UCSC_CpG_Island
  )
  overlaps_CpGs_Blueprint <- findOverlaps(Blueprint_Granges, CpGs_Granges)
  match_hit_CpGs_Blueprint <- data.frame(mcols(Blueprint_Granges[queryHits(overlaps_CpGs_Blueprint),]),
                                         data.frame(mcols(CpGs_Granges[subjectHits(overlaps_CpGs_Blueprint),])))
  anno_focused <- anno_focused %>% dplyr::filter(., .$Name %ni% match_hit_CpGs_Blueprint$CpG)
  df_length <- length(anno_focused$Name)
  to_add <- data.frame("Blueprint_gene_names" = rep(gene, df_length),
                       "type" = rep("Illumina_annotation", df_length),
                       "gene_type" = rep("", df_length),
                       "CpG" = anno_focused$Name,
                       "Illumina_Gene_name" = anno_focused$UCSC_RefGene_Name,
                       "position" = anno_focused$UCSC_RefGene_Group,
                       "Island" = anno_focused$Relation_to_UCSC_CpG_Island)
  final_matchit <- rbind(match_hit_CpGs_Blueprint, to_add)
  Specific_CpGs_value <- merge(final_matchit, BMIQ_focused, by.x = "CpG", by.y = 0, all.x = TRUE) %>%
    separate_rows(., Illumina_Gene_name, position, sep = ";")
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

  ggplot(Beta_value, aes(y=Beta_values, x=cpg_localisation, fill=cpg_localisation, colour= type))+
    geom_boxplot(alpha=0.25)+
    geom_jitter(width=0.25, alpha=0.25)+
    ggtitle(paste0(gene," Promoter Methylation"))+
    xlab("Position of CpGs")+
    ylab("Value of the methylation")+
  scale_shape_manual(values=c(0,1,2,5,6,7))
}

Gene_name_finding("AKT")

Focus_Gene("AKT2")
