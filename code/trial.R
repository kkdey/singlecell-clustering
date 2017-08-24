

library(maptpx)
library(Matrix)

counts <- get(load("../data/pbmc_68K_evan.rda"))
counts <- t(as.matrix(counts))
dim(counts)
library(Seurat)
rownames(counts) <- 1:dim(counts)[1]
seuratObj <- new('seurat', raw.data = t(counts))

seuratObj_TFH_global <- seuratObj


seuratObj_TFH_global <- FindVariableGenes(seuratObj_TFH_global, fxn.x = expMean, fxn.y = logVarDivMean, x.low.cutoff = 0.25, x.high.cutoff = 3, y.cutoff = 0.5, do.contour = FALSE)
length(seuratObj_TFH_global@var.genes)
seuratObj_TFH_global <- RunPCA(seuratObj_TFH_global, pc.genes = seuratObj_TFH_global@var.genes)
seuratObj_TFH_global <- ProjectPCA(seuratObj_TFH_global)
VizPCA(seuratObj_TFH_global, 1:2)

save(seuratObj_TFH_global, file = "../data/seurat_pbmc_68K_evan.rda")

genes1 <- rownames(counts)
out <- mygene::queryMany(genes1, scopes = "ensembl.gene", fields="symbol", species="human")
genes2 <- out$symbol
genes2_query <- out$query

idx <- match(genes2_query[!is.na(genes2_query)], rownames(counts))
counts2 <- counts[idx,]

counts3 <- counts2[-grep("[.]", genes2),]
genes3 <- genes2[-grep("[.]", genes2)]

idx2 <- which(!is.na(match(substring(genes3, 1, 2), c("RP", "MT", "RN"))))
counts4 <- counts3[-idx2, ]



