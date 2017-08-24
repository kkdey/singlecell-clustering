

################  PBMC data analysis  (+ sorted cells)  ############################

library(Seurat)
library(Matrix)

pbmc_data <- Read10X('../data/filtered_matrices_mex/hg19/')
idx <- which(!is.na(match(substring(rownames(pbmc_data), 1, 2), c("RP", "MT", "RN"))))
counts <- pbmc_data[-idx, ]

idx2 <- sample(1:ncol(counts), 5000, replace=FALSE)
counts2 <- counts[,idx2]

save(counts2, file = "../data/pbmc_68K_evan.rda")
