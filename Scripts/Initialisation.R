source("~/GitHub/Epigenomic_integration/Scripts/functions.R")

gene_list <- c("CEBPA", "CPT1A", "CPT2", "SLC25A20", "AKT2", "PPARGC1A", "PTEN")

Blueprint_data <- list()
Blueprint_data[["Blueprint"]] <- read.csv("~/GitHub/Epigenomic_integration/DATA/BLUEPRINT_fragments_good.tsv", sep = "\t") %>%
  dplyr::select(., "gene_names":"type", "gene_type") %>%
  separate_rows(., gene_names, sep = " ") %>%
  separate_rows(., gene_type, sep = " ") %>%
  unique(.)

Pchic_data <- list()
Pchic_data[["pchic"]] <- prepare_pchic()
Pchic_data[["bed"]] <- unique(rbind(Pchic_data[["pchic"]][, c(1:3, 5)], Pchic_data[["pchic"]][, c(6:8, 10)]))
Pchic_data[["GRange"]] <- GRanges(
  seqnames = Pchic_data[["bed"]]$chr,
  IRanges(start = Pchic_data[["bed"]]$start, end = Pchic_data[["bed"]]$end),
  Gene_Pchic = Pchic_data[["bed"]]$Name,
  start_fragment = Pchic_data[["bed"]]$start,
  end_fragment = Pchic_data[["bed"]]$end
)
Pchic_data[["GRange"]]$ID <- paste(Pchic_data[["bed"]]$chr, Pchic_data[["bed"]]$start, sep = "_")
colnames(Pchic_data[["pchic"]]) <- c("chr_bait", "start_bait", "end_bait", "ID_bait", "Name_bait", "chr_oe", "start_oe", "end_oe", "ID_oe", "Name_oe")
Pchic_data[["pchic"]]$IDbait <- paste(Pchic_data[["pchic"]]$chr_bait, Pchic_data[["pchic"]]$start_bait, sep = "_")
Pchic_data[["pchic"]]$IDoe <- paste(Pchic_data[["pchic"]]$chr_oe, Pchic_data[["pchic"]]$start_oe, sep = "_")

Blueprint_data[["Blueprint"]] <- Blueprint
Blueprint_data[["GRanges"]] <- GRanges(
  seqnames = Blueprint_data[["Blueprint"]]$chr,
  ranges = IRanges(start = Blueprint_data[["Blueprint"]]$start, end = Blueprint_data[["Blueprint"]]$end),
  Blueprint_gene_names = Blueprint_data[["Blueprint"]]$gene_names,
  type = Blueprint_data[["Blueprint"]]$type,
  gene_type = Blueprint_data[["Blueprint"]]$gene_type
)

Illumina_annotations <- list()
## ILLUMINA 450K
Illumina_annotations[["450"]] <- read.csv("~/Illumina_Manifest/HumanMethylation450_15017482_v1-2.csv", skip = 7) %>% 
  dplyr::select(., "Name", "CHR", "MAPINFO", "UCSC_RefGene_Name", "UCSC_RefGene_Group", "Relation_to_UCSC_CpG_Island") %>%
  dplyr::filter(., CHR != "") %>%
  separate_rows(., UCSC_RefGene_Name, UCSC_RefGene_Group, sep = ";")

Illumina_annotations[["CpG_450_GRanges"]] <- GRanges(
  seqnames = Illumina_annotations[["450"]]$CHR,
  ranges = IRanges(Illumina_annotations[["450"]]$MAPINFO, Illumina_annotations[["450"]]$MAPINFO +1),
  CpG = Illumina_annotations[["450"]]$Name,
  Illumina_Gene_name = Illumina_annotations[["450"]]$UCSC_RefGene_Name,
  position = Illumina_annotations[["450"]]$UCSC_RefGene_Group,
  Island = Illumina_annotations[["450"]]$Relation_to_UCSC_CpG_Island
)
overlaps <- findOverlaps(Blueprint_data[["GRanges"]], Illumina_annotations[["CpG_450_GRanges"]])
Overlap_data <- list()
Overlap_data[["450_Blueprint"]] <- data.frame(mcols(Blueprint_data[["GRanges"]][queryHits(overlaps),]),
                                           data.frame(mcols(Illumina_annotations[["CpG_450_GRanges"]][subjectHits(overlaps),])))

Overlap_data[["450_Blueprint_promoter"]] <- dplyr::filter(Overlap_data[["450_Blueprint"]], type == "P")


####################################################################################################

anno_promoter_450 <- dplyr::filter(anno_450, UCSC_RefGene_Group == "TSS1500" | UCSC_RefGene_Group == "TSS200")

overlaps_CpGs_Pchic_450 <- findOverlaps(CpGs_Granges_450, PCHiC_GRange)
matchit_CpGs_Pchic_450 <- data.frame(mcols(CpGs_Granges_450[queryHits(overlaps_CpGs_Pchic_450),]),
                                     data.frame(mcols(PCHiC_GRange[subjectHits(overlaps_CpGs_Pchic_450),])))
gene_universe_450K <- unique(match_hit_CpGs_Blueprint_promoter_450$Blueprint_gene_names)
#########################
## ILLUMINA EPIC
anno_EPIC <- read.csv("~/Illumina_Manifest/MethylationEPIC_v-1-0_B4.csv", skip = 7) %>% 
  dplyr::select(., "Name", "CHR", "MAPINFO", "UCSC_RefGene_Name", "UCSC_RefGene_Group", "Relation_to_UCSC_CpG_Island") %>%
  dplyr::filter(., CHR != "") %>%
  separate_rows(., UCSC_RefGene_Name, UCSC_RefGene_Group, sep = ";")

CpGs_Granges_EPIC <- GRanges(
  seqnames = anno_EPIC$CHR,
  ranges = IRanges(anno_EPIC$MAPINFO, anno_EPIC$MAPINFO +1),
  CpG = anno_EPIC$Name,
  Illumina_Gene_name = anno_EPIC$UCSC_RefGene_Name,
  position = anno_EPIC$UCSC_RefGene_Group,
  Island = anno_EPIC$Relation_to_UCSC_CpG_Island
)
overlaps_CpGs_Blueprint_EPIC <- findOverlaps(Blueprint_Granges, CpGs_Granges_EPIC)
match_hit_CpGs_Blueprint_EPIC <- data.frame(mcols(Blueprint_Granges[queryHits(overlaps_CpGs_Blueprint_EPIC),]),
                                            data.frame(mcols(CpGs_Granges_EPIC[subjectHits(overlaps_CpGs_Blueprint_EPIC),])))

match_hit_CpGs_Blueprint_promoter_EPIC <- dplyr::filter(match_hit_CpGs_Blueprint_EPIC, type == "P")

anno_promoter_EPIC <- dplyr::filter(anno_EPIC, UCSC_RefGene_Group == "TSS1500" | UCSC_RefGene_Group == "TSS200")

overlaps_CpGs_Pchic_EPIC <- findOverlaps(CpGs_Granges_EPIC, PCHiC_GRange)
matchit_CpGs_Pchic_EPIC <- data.frame(mcols(CpGs_Granges_EPIC[queryHits(overlaps_CpGs_Pchic_EPIC),]),
                                      data.frame(mcols(PCHiC_GRange[subjectHits(overlaps_CpGs_Pchic_EPIC),])))
gene_universe_EPIC <- unique(match_hit_CpGs_Blueprint_promoter_EPIC$Blueprint_gene_names)
##########################

overlaps_Blueprint_Pchic <- findOverlaps(Blueprint_Granges, PCHiC_GRange)
matchit_Blueprint_Pchic <- data.frame(mcols(Blueprint_Granges[queryHits(overlaps_Blueprint_Pchic),]),
                                      data.frame(mcols(PCHiC_GRange[subjectHits(overlaps_Blueprint_Pchic),]))) %>%
  dplyr::filter(., type == "P")

signature <- read.csv("~/GitHub/Epigenomic_integration/DATA/BPmetCan.txt", sep = "\t")
rownames(signature) <- signature$CpGs
signature <- signature[,-1]

number_CpG_per_gene_450 <- table(match_hit_CpGs_Blueprint_promoter_450$Blueprint_gene_names) %>%
  as.data.frame(.) %>%
  .[-1,]

number_CpG_per_gene_EPIC <- table(match_hit_CpGs_Blueprint_promoter_EPIC$Blueprint_gene_names) %>%
  as.data.frame(.) %>%
  .[-1,]

Gene_lists_folder <- "~/GitHub/Epigenomic_integration/Results/Tables report/Gene lists/"
GO_folder <- "~/GitHub/Epigenomic_integration/Results/Tables report/GO/"

rm(overlaps_Blueprint_Pchic, overlaps_CpGs_Pchic_EPIC, overlaps_CpGs_Blueprint_EPIC, overlaps_CpGs_Pchic_450, overlaps_CpGs_Blueprint_450)
