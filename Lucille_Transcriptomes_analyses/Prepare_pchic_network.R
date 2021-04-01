library(dplyr)
library(dorothea)
library(tidyr)

load("~/PCHIC/pchic.RData")
pchic <- data.frame(pchic[rowSums(pchic[,c("Mon", "Mac1", "Mac0", "Mac2", "MK", "Ery", "EP")] >= 5) >= 1, 1:10]) %>% na.omit(.)

data(dorothea_hs, package = "dorothea")
regulons = dorothea_hs %>%
  dplyr::filter(confidence %in% c("A", "B"))

regulons_KLF4_targets <- dplyr::filter(regulons, tf == "KLF4") %>% .$target

KLF4_target_neighborhood <- dplyr::filter(pchic, baitName %in% regulons_KLF4_targets)

pchic_named_network <- pchic
pchic_named_network$IN <- pchic_named_network$baitName
pchic_named_network$OUT <- ifelse(pchic_named_network$oeName == ".", pchic_named_network$oeID, pchic_named_network$oeName)
pchic_named_network <- separate_rows(pchic_named_network, OUT, sep = ";")
pchic_named_network <- separate_rows(pchic_named_network, IN, sep = ";")
pchic_named_network <- dplyr::select(pchic_named_network, c("IN", "OUT"))

write.table(pchic_named_network, "GitHub/Epigenomic_integration/DATA/pchic_named_network.tsv", row.names = F, sep = '\t')
