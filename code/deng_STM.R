

#############  Deng et al Structured Topic Model analysis  #####################

library(singleCellRNASeqMouseDeng2014)
deng.counts <- Biobase::exprs(Deng2014MouseESC)
deng.meta_data <- Biobase::pData(Deng2014MouseESC)
deng.gene_names <- rownames(deng.counts)

counts_data <- t(deng.counts)

library(stm)
library(igraph)

documents <- vector(mode = "list", length = dim(counts_data)[1])

for(l in 1:dim(counts_data)[1]){
  idx <- which(counts_data[l,] != 0)
  counts_idx <- counts_data[l, idx]
  documents[[l]] <- as.matrix(rbind(idx, counts_idx))
}

genes_vocab <- colnames(counts_data)
samples_meta <- deng.meta_data

dengstmGoMFit <- stm(documents = documents, vocab = genes_vocab,
                       K = 2,
                       max.em.its = 30, data = samples_meta,
                       init.type = "Spectral")

