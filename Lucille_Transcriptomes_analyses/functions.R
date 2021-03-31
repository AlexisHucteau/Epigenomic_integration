Prepare_Transcriptomic_microarray_data <- function(Raw_DATA_Dir, Phenotype_filename){
  Annot <- data.frame(ACCNUM = sapply(contents(hugene20sttranscriptclusterACCNUM), paste, collapse = ", "), 
                      SYMBOL = sapply(contents(hugene20sttranscriptclusterSYMBOL), paste, collapse = ", "), 
                      DESC = sapply(contents(hugene20sttranscriptclusterGENENAME), paste, collapse = ", "))
  SDRF <- read.csv(Phenotype_filename)
  celFiles <- list.celfiles(Raw_DATA_Dir, full.names = T)
  data <- read.celfiles(celFiles)
  eset <- rma(data)
  pheno <- SDRF
  exp_raw <- as.matrix(exprs(eset))
  my_frame <- exp_raw
  exprs(eset) <- exp_raw
  my_frame <- merge(Annot, my_frame, by.y = 0, by.x = 0, all = T)
  transcriptomes <- my_frame
  rownames(transcriptomes) <- transcriptomes[,1]
  transcriptomes <- transcriptomes[,c(5:29)]
  res <- list("Annot" = my_frame[,c(1:4)],
              "Transcriptomes" = transcriptomes,
              "Phenotype" = SDRF)
  return(res)
}

PCA_analysis <- function(DATA_for_PCA, Variable_to_analyse){
  res.pca <- PCA(DATA_for_PCA, graph = TRUE)
  fviz_eig(res.pca, addlabels = TRUE)
}


DGEs <- function(data, Phenotype, Annotations){
  require(limma)
  require(dplyr)
  #require(dendextend)
  
  
  message("[===========================]")
  message("[<<<<<<< DGEs START >>>>>>>>>]")
  message("[<<<< Pairwise analysis >>>>>]")
  message("-----------------------------")
  
  # This function creates the pairs for the pairwise matrices
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
  
  # This function creates the pairs for the pairwise matrices
  
  design <- model.matrix(~0 + Phenotype)
  
  #Removing heteroscedascity from data
  
  contr.matrix <- design.pairs(levels(factor(Phenotype)))
  colnames(design) <- rownames(contr.matrix)   
  
  # Fitting linear models for comparisons of interest
  Fit <- lmFit(data, design) %>%
    contrasts.fit(., contr.matrix) %>%
    eBayes(., trend = TRUE)
  
  FitList <- list()
  for (i in 1:ncol(contr.matrix)) {
    FitList[[i]] <- topTable(Fit, coef = i, adjust.method = "BH", number = nrow(data)) %>%
      mutate(ID = rownames(.))
    FitList[[i]] <- merge(FitList[[i]], Annotations, by.x = "ID", by.y = "Row.names", all.x = T)
    
    message(paste0(i, " done"))
    
  }
  names(FitList) <- colnames(contr.matrix)
  return(FitList)
}

volcanoplot_gene_expression <- function(DATA, title, xlim = 10, ylim = 5){
  Data <- dplyr::filter(DATA, SYMBOL != "NA")
  gene <- Data[,"SYMBOL"]
  
  enhancedvolcano_data <- data.frame(logFC = Data[, "logFC"],
                                      pvalue = Data[, "P.Value"],
                                      gene = gene)
  EnhancedVolcano(toptable = enhancedvolcano_data, 
                  lab = enhancedvolcano_data$gene, 
                  x = "logFC", 
                  y = "pvalue",
                  FCcutoff = 2.5,
                  pCutoff = 0.01,
                  title = title,
                  subtitle = NA,
                  legendPosition = "right",
                  subtitleLabSize = 0,
                  legendLabSize = 10,
                  ylim = c(0, ylim),
                  xlim = c(-xlim, xlim)
  )
  ggsave(paste("Volcano_plot_", title, ".png", sep=""))
  EnhancedVolcano(toptable = enhancedvolcano_data, 
                  lab = enhancedvolcano_data$gene, 
                  x = "logFC", 
                  y = "pvalue",
                  FCcutoff = 2.5,
                  pCutoff = 0.01,
                  title = title,
                  subtitle = NA,
                  legendPosition = "right",
                  subtitleLabSize = 0,
                  legendLabSize = 10,
                  ylim = c(0, ylim),
                  xlim = c(-xlim, xlim)
  )
  
}


Differential_TF_activities_t_test <- function(TF_activities_predictions, Phenotype_A, Phenotype_B, Variable){
  A <- Variable == Phenotype_A
  B <- Variable == Phenotype_B
  res <- sapply(colnames(TF_activities_predictions), function(x){
    test_tmp <- t.test(TF_activities_predictions[A, x], TF_activities_predictions[B, x])
    tmp <- c("p.value_activities" = test_tmp$p.value, 
      "A" = test_tmp$estimate[1], 
      "B" = test_tmp$estimate[2])
    tmp
  })
  
  rownames(res) <- c("p.value_activities", Phenotype_A, Phenotype_B)
  return(data.frame(t(res)))
}
