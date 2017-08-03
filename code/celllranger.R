

################  cellranger R kit   ################################

library(cellrangerRkit)

library(Seurat)
library(Matrix)

data <- Read10X("../data/GS_PBMC_8_2013_MCCB1050/filtered_gene_bc_matrices/hg38_merkerCellPolyomavirus_CellRanger_v2/")
data2 <- new("seurat",  raw.data = data)

data3 <- Setup(data2, min.cells = 3, min.genes = 200, project = "10X_PBMC", do.scale = F, do.center = F, names.field = 2,
      names.delim = "\\-")

mito.genes <- grep("^MT-", rownames(data3@data), value = T)
percent.mito <- colSums(expm1(data3@data[mito.genes, ])) / colSums(expm1(data3@data))
data3 <- AddMetaData(data3, percent.mito, "percent.mito")

colSums(data)

counts <- t(data)
idx <- which(!is.na(match(substring(colnames(counts), 1, 2), c("RP", "MT", "RN"))))
counts2 <- counts[, -idx]
counts3 <- counts2[,-(1:4)]
counts4 <- counts3[, - (grep("[.]", colnames(counts3)))]
counts5 <- counts4[, - (grep("or", colnames(counts4)))]
counts6 <- counts5[, - (grep("-", colnames(counts5)))]

save(counts6, file = "../output/GS_PBMC_8_2013_MCCB1050.rda")

counts6  <- get(load("../output/GS_PBMC_2_2015_GRN0304.rda"))
topic.clus.list <- list()
for(k in 2:4){
  topic.clus.list[[k]] <- maptpx::topics(counts6, K=k, tol = 10)
}
save(topic.clus.list, file = "../output/maptpx_GS_PBMC_2_2015_GRN0304.rda")
