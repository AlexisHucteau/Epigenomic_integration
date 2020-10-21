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
  load("../pchic.RData")
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
