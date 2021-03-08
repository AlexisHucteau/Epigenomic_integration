source("~/GitHub/Epigenomic_integration/Scripts/functions.R")

Gene_ontology <- list()

Blueprint_data <- list()
Blueprint_data[["Blueprint"]] <- read.csv("~/GitHub/Epigenomic_integration/DATA/BLUEPRINT_fragments_good.tsv", sep = "\t") %>%
  dplyr::select(., "gene_names":"type", "gene_type") %>%
  separate_rows(., gene_names, sep = " ") %>%
  separate_rows(., gene_type, sep = " ") %>%
  unique(.)

Pchic_data <- list()
Pchic_data[["pchic"]] <- prepare_pchic(cell_lines = c("Mon", "Mac1", "Mac0", "Mac2", "MK", "Ery", "EP"))
Pchic_data[["bed"]] <- unique(rbind(Pchic_data[["pchic"]][, c(1:3, 5)], Pchic_data[["pchic"]][, c(6:8, 10)]))
Pchic_data[["GRanges"]] <- GRanges(
  seqnames = Pchic_data[["bed"]]$chr,
  IRanges(start = Pchic_data[["bed"]]$start, end = Pchic_data[["bed"]]$end),
  Gene_Pchic = Pchic_data[["bed"]]$Name,
  start_fragment = Pchic_data[["bed"]]$start,
  end_fragment = Pchic_data[["bed"]]$end
)
Pchic_data[["GRanges"]]$ID <- paste(Pchic_data[["bed"]]$chr, Pchic_data[["bed"]]$start, sep = "_")
colnames(Pchic_data[["pchic"]]) <- c("chr_bait", "start_bait", "end_bait", "ID_bait", "Name_bait", "chr_oe", "start_oe", "end_oe", "ID_oe", "Name_oe")
Pchic_data[["pchic"]]$IDbait <- paste(Pchic_data[["pchic"]]$chr_bait, Pchic_data[["pchic"]]$start_bait, sep = "_")
Pchic_data[["pchic"]]$IDoe <- paste(Pchic_data[["pchic"]]$chr_oe, Pchic_data[["pchic"]]$start_oe, sep = "_")

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

Illumina_annotations[["450_promoter"]] <- dplyr::filter(Illumina_annotations[["450"]], UCSC_RefGene_Group == "TSS1500" | UCSC_RefGene_Group == "TSS200")

overlaps <- findOverlaps(Illumina_annotations[["CpG_450_GRanges"]], Pchic_data[["GRanges"]])
Overlap_data[["450_CpG_chromatin"]] <- data.frame(mcols(Illumina_annotations[["CpG_450_GRanges"]][queryHits(overlaps),]),
                                     data.frame(mcols(Pchic_data[["GRanges"]][subjectHits(overlaps),])))
Gene_lists <- list()
Gene_lists[["gene_universe_450K"]] <- unique(Overlap_data[["450_Blueprint_promoter"]]$Blueprint_gene_names)
#########################
## ILLUMINA EPIC
Illumina_annotations[["EPIC"]] <- read.csv("~/Illumina_Manifest/MethylationEPIC_v-1-0_B4.csv", skip = 7) %>% 
  dplyr::select(., "Name", "CHR", "MAPINFO", "UCSC_RefGene_Name", "UCSC_RefGene_Group", "Relation_to_UCSC_CpG_Island") %>%
  dplyr::filter(., CHR != "") %>%
  separate_rows(., UCSC_RefGene_Name, UCSC_RefGene_Group, sep = ";")

Illumina_annotations[["CpG_EPIC_GRanges"]] <- GRanges(
  seqnames = Illumina_annotations[["EPIC"]]$CHR,
  ranges = IRanges(Illumina_annotations[["EPIC"]]$MAPINFO, Illumina_annotations[["EPIC"]]$MAPINFO +1),
  CpG = Illumina_annotations[["EPIC"]]$Name,
  Illumina_Gene_name = Illumina_annotations[["EPIC"]]$UCSC_RefGene_Name,
  position = Illumina_annotations[["EPIC"]]$UCSC_RefGene_Group,
  Island = Illumina_annotations[["EPIC"]]$Relation_to_UCSC_CpG_Island
)
overlaps <- findOverlaps(Blueprint_data[["GRanges"]], Illumina_annotations[["CpG_EPIC_GRanges"]])
Overlap_data[["EPIC_Blueprint"]] <- data.frame(mcols(Blueprint_data[["GRanges"]][queryHits(overlaps),]),
                                            data.frame(mcols(Illumina_annotations[["CpG_EPIC_GRanges"]][subjectHits(overlaps),])))

Overlap_data[["EPIC_Blueprint_promoter"]] <- dplyr::filter(Overlap_data[["EPIC_Blueprint"]], type == "P")

Illumina_annotations[["EPIC_promoter"]] <- dplyr::filter(Illumina_annotations[["EPIC"]], UCSC_RefGene_Group == "TSS1500" | UCSC_RefGene_Group == "TSS200")

overlaps <- findOverlaps(Illumina_annotations[["CpG_EPIC_GRanges"]], Pchic_data[["GRanges"]])
Overlap_data[["EPIC_CpG_chromatin"]] <- data.frame(mcols(Illumina_annotations[["CpG_EPIC_GRanges"]][queryHits(overlaps),]),
                                      data.frame(mcols(Pchic_data[["GRanges"]][subjectHits(overlaps),])))
Gene_lists[["gene_universe_EPIC"]] <- unique(Overlap_data[["EPIC_CpG_chromatin"]]$Blueprint_gene_names)
##########################

overlaps <- findOverlaps(Blueprint_data[["GRanges"]], Pchic_data[["GRanges"]])
Overlap_data[["Blueprint_Chromatin"]] <- data.frame(mcols(Blueprint_data[["GRanges"]][queryHits(overlaps),]),
                                      data.frame(mcols(Pchic_data[["GRanges"]][subjectHits(overlaps),]))) %>%
  dplyr::filter(., type == "P")

Gene_lists[["BPmetCan"]] <- read.csv("~/GitHub/Epigenomic_integration/DATA/BPmetCan.txt", sep = "\t")
rownames(Gene_lists[["BPmetCan"]]) <- Gene_lists[["BPmetCan"]]$CpGs
Gene_lists[["BPmetCan"]] <- Gene_lists[["BPmetCan"]][,-1]

Gene_lists_folder <- "~/GitHub/Epigenomic_integration/Results/Tables report/Gene lists/"
GO_folder <- "~/GitHub/Epigenomic_integration/Results/Tables report/GO/"

tmp_data <- list()
tmp_data[["Promoter_fragments"]] <- Pchic_data[["pchic"]]$IDbait %>% unique(.)

Pchic_data[["Non_promoter_pchic"]] <- Pchic_data[["pchic"]] %>% 
  dplyr::filter(., IDoe %ni% tmp_data[["Promoter_fragments"]]) %>%
  dplyr::select(., IDbait, IDoe) %>%
  split(., .$IDbait)

Pchic_data[["Non_promoter_pchic"]] <- lapply(Pchic_data[["Non_promoter_pchic"]], function(x){
  x$IDoe %>% unique(.)
})

Overlap_data[["Fragment_connected_per_gene"]] <- split(Overlap_data[["Blueprint_Chromatin"]], Overlap_data[["Blueprint_Chromatin"]]$Blueprint_gene_names)
Overlap_data[["Fragment_connected_per_gene"]] <- Overlap_data[["Fragment_connected_per_gene"]][-1]

Overlap_data[["Fragment_connected_per_gene"]] <- lapply(Overlap_data[["Fragment_connected_per_gene"]], function(x){
  x[,c(1,7)] %>% unique(.)
})

Overlap_data[["Fragment_connected_per_gene"]] <- lapply(Overlap_data[["Fragment_connected_per_gene"]], function(x){
  Pchic_data[["Non_promoter_pchic"]][x$ID]
})

Overlap_data[["Fragment_connected_per_gene"]] <- lapply(Overlap_data[["Fragment_connected_per_gene"]], function(x){
  unlist(x) %>% unique(.)
})

rm(overlaps)


tmp_data[["Non_Promoter_fragments"]] <- Pchic_data[["Non_promoter_pchic"]] %>% unlist(.) %>% unique(.)

Overlap_data[["450_CpG_chromatin_enhancer"]] <- Overlap_data[["450_CpG_chromatin"]] %>%
  dplyr::filter(., ID %in% tmp_data[["Non_Promoter_fragments"]])

CpGs_per_fragment <- list()
CpGs_per_fragment[["450"]] <- Overlap_data[["450_CpG_chromatin_enhancer"]] %>%
  dplyr::select(., c("CpG", "ID")) %>%
  split(., Overlap_data[["450_CpG_chromatin_enhancer"]]$ID)

Overlap_data[["EPIC_CpG_chromatin_enhancer"]] <- Overlap_data[["EPIC_CpG_chromatin"]] %>%
  dplyr::filter(., ID %in% tmp_data[["Non_Promoter_fragments"]])

CpGs_per_fragment[["EPIC"]] <- Overlap_data[["EPIC_CpG_chromatin_enhancer"]] %>%
  dplyr::select(., c("CpG", "ID")) %>%
  split(., Overlap_data[["EPIC_CpG_chromatin_enhancer"]]$ID)
