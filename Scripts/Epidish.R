library(EpiDISH)

library(biomaRt)
library(stringr)

signature <- read.delim("~/GitHub/Epigenomic_integration/DATA/BPRNACan_150.txt")



t1 <- read.table(snakemake@input[[1]], header = TRUE, row.names = 1)

res <- epidish(t1, signature, method = "RPC")
Fres <- as.data.frame(res$estF)
Fres <- cbind(Sample = rownames(Fres), Fres)
Fres[, -1] <- round(Fres[, -1], 3)

write.table(Fres, snakemake@output[[1]], sep = "\t", quote = F, row.names = F, col.names = T)



ensembl = useMart("ensembl", dataset = "hsapiens_gene_ensembl")

gene.annotations <- biomaRt::getBM(mart = ensembl, attributes=c("ensembl_gene_id", "external_gene_name", "start_position", "end_position"))
gene.annotations <- dplyr::transmute(gene.annotations, external_gene_name,  ensembl_gene_id, length = end_position - start_position)
RNAseq_Wang_Feng <- as.data.frame(RNAseq_Wang_Feng)

RNAseq_Wang_Feng$gene <- rownames(RNAseq_Wang_Feng)
RNAseq_Wang_Feng <- RNAseq_Wang_Feng %>%
  separate_rows(., gene, sep = "[|]")

RNAseq_Wang_Feng <- RNAseq_Wang_Feng %>%
  dplyr::filter(., gene %in% gene_formated)

gene_duplicated <- unique(RNAseq_Wang_Feng$gene[duplicated(RNAseq_Wang_Feng$gene)])

combined_dt <- RNAseq_Wang_Feng %>%
  dplyr::filter(., gene %ni% gene_duplicated)

for (gene_d in gene_duplicated) {
  tmp <- dplyr::filter(RNAseq_Wang_Feng, gene == gene_d) %>% 
    .[,c(1:51)] %>%
    as.matrix(.) %>%
    colSums2(.)
  tmp <- c(tmp, gene_d)
  combined_dt <- rbind(combined_dt, tmp)
}
 
rownames(combined_dt) <- combined_dt$gene
combined_dt <- combined_dt

# Filter and re-order gene.annotations to match the order in feature counts matrix
gene.annotations <- gene.annotations %>% dplyr::filter(external_gene_name %in% rownames(combined_dt))
gene.annotations <- gene.annotations[order(match(gene.annotations$external_gene_name, rownames(combined_dt))),]

# Assign feature lenghts into a numeric vector.
featureLength <- gene.annotations$length

combined_dt_with_fragment_length <- combined_dt[gene.annotations$external_gene_name,]

RNAseq_Wang_Feng_tpm <- counts_to_tpm(as.matrix(combined_dt_with_fragment_length), featureLength, rep(50,51))

counts_to_tpm <- function(counts, featureLength, meanFragmentLength)
{
  
  stopifnot(length(featureLength) == nrow(counts))
  stopifnot(length(meanFragmentLength) == ncol(counts))
  
  # Compute effective lengths of features in each library.
  effLen <- do.call(cbind, lapply(1:ncol(counts), function(i){featureLength - meanFragmentLength[i] + 1}))
  print(effLen)
  # Exclude genes with length less than the mean fragment length.
  idx <- apply(effLen, 1, function(x) min(x) > 1)
  counts <- counts[idx,]
  effLen <- effLen[idx,]
  featureLength <- featureLength[idx]
  
  # Process one column at a time.
  print("hello")
  tpm <- do.call(cbind, lapply(1:ncol(counts), function(i) 
  {
    rate = log(counts[,i]) - log(effLen[,i])
    denom = log(sum(exp(rate)))
    exp(rate - denom + log(1e6))
  }))
  colnames(tpm) <- colnames(counts)
  rownames(tpm) <- rownames(counts)
  return(tpm)
}
