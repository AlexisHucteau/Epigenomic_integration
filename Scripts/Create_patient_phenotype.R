Patient_data <- data.frame("ID_patient" = Phenotype_Wang_Feng$Sample %>% unique(.))
res <- c()

for (i in Patient_data$ID_patient){
  tmp <- Phenotype_Wang_Feng[Phenotype_Wang_Feng$Sample == i,]
  if("Responder" %in% tmp$Phenotype){
    res <- c(res,"Responder")
  }else if("Non_Responder" %in% tmp$Phenotype){
    res <- c(res,"Non_Responder")
  }else{
    res <- c(res,NA)
  }
}

Patient_data$Baseline_phenotype <- res


res <- c()
for (i in Patient_data$ID_patient){
  tmp <- dplyr::filter(Phenotype_Wang_Feng, Sample == i & Phenotype == "Baseline")
  if(length(tmp$Sample_number) == 0){
    res <- c(res, NA)
  }else{
    res <- c(res, tmp$Sample_number)
  }
}

Patient_data$Baseline_Sample <- res

res <- c()
for (i in Patient_data$ID_patient){
  tmp <- dplyr::filter(Phenotype_Wang_Feng, Sample == i & Phenotype != "Baseline" & Phenotype != "Control")
  if(length(tmp$Sample_number) == 0){
    res <- c(res, NA)
  }else{
    res <- c(res, tmp$Sample_number)
  }
}
Patient_data$Post_treatment_sample <- res

Phenotype_Wang_Feng_RNAseq$Sample_name <- Phenotype_Wang_Feng_RNAseq$Sample %>% str_split(., "_") %>% lapply(., function(x){
  x[1]
}) %>%
  unlist(.)

Phenotype_Wang_Feng_RNAseq$Phenotype_RNAseq <- Phenotype_Wang_Feng_RNAseq$Phenotype

Phenotype_Wang_Feng_RNAseq <- Phenotype_Wang_Feng_RNAseq[,-2]

res <- c()
for (i in Patient_data$ID_patient){
  tmp <- dplyr::filter(Phenotype_Wang_Feng_RNAseq, Sample_name == i & Phenotype_RNAseq == "Baseline")
  if(length(tmp$Sample) == 0){
    res <- c(res, NA)
  }else{
    res <- c(res, tmp$Sample)
  }
}

Patient_data$Baseline_RNAseq_data <- res

res <- c()
for (i in Patient_data$ID_patient){
  tmp <- dplyr::filter(Phenotype_Wang_Feng_RNAseq, Sample_name == i & Phenotype_RNAseq == "Relapse")
  if(length(tmp$Sample) == 0){
    res <- c(res, NA)
  }else{
    res <- c(res, tmp$Sample)
  }
}

Patient_data$Relapse_RNAseq_data <- res

write.csv(Patient_data, "../DATA/Patient_phenotype_Koichi_cohort.csv", row.names = FALSE)