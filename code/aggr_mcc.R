

##########   aggregated data MCC analysis  ##############

counts61 <- get(load("../output/GS_PBMC_8_2013_MCCB1050.rda"))
counts61 <- as.matrix(counts61)
counts62 <- get(load("../output/GS_PBMC_2_2015_GRN0304.rda"))
counts62 <- as.matrix(counts62)
counts63 <- get(load("../output/GS_PBMC_2_2016_GRN0535.rda"))
counts63 <- as.matrix(counts63)
counts64 <- get(load("../output/GS_PBMC_10_2016_GRN0760.rda"))
counts64 <- as.matrix(counts64)
aggregated_data <- rbind(counts61, counts62, counts63, counts64)
dim(aggregated_data)

fac <- c(rep("Aug-2013", dim(counts61)[1]),
         rep("Feb-2015", dim(counts62)[1]),
         rep("Feb-2016", dim(counts63)[1]),
         rep("Oct-2016", dim(counts64)[1]))
fac <- factor(fac, levels = c("Aug-2013", "Feb-2015", "Feb-2016", "Oct-2016"))
length(fac)

rownames(aggregated_data) <- 1:dim(aggregated_data)[1]
seuratObj <- new('seurat', raw.data = t(aggregated_data))
seuratObj <- CreateSeuratObject(raw.data = t(aggregated_data), min.cells = 0, min.genes = 0, 
                                project = 'Pantaleo') # no need to filter - already done
#- Output
seuratObj_TFH_global <- seuratObj
seuratObj_TFH_global <- NormalizeData(object = seuratObj_TFH_global, normalization.method = "LogNormalize", 
                      scale.factor = 10000)
seuratObj_TFH_global <- FindVariableGenes(object = seuratObj_TFH_global, mean.function = ExpMean, dispersion.function = LogVMR, 
                          x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
length(x = seuratObj_TFH_global@var.genes)
mito.genes <- grep(pattern = "^MT-", x = rownames(x = seuratObj_TFH_global@data), value = TRUE)
percent.mito <- colSums(seuratObj_TFH_global@raw.data[mito.genes, ])/colSums(seuratObj_TFH_global@raw.data)

seuratObj_TFH_global <- AddMetaData(object = seuratObj_TFH_global, metadata = percent.mito, col.name = "percent.mito")
seuratObj_TFH_global <- ScaleData(object = seuratObj_TFH_global, vars.to.regress = c("nUMI", "percent.mito"))

seuratObj_TFH_global <- RunPCA(object = seuratObj_TFH_global, pc.genes = pbmc@var.genes, do.print = TRUE, pcs.print = 1:5, 
               genes.print = 5)
seuratObj_TFH_global <- RunTSNE(object = seuratObj_TFH_global, dims.use = 1:10, do.fast = TRUE)







seuratObj_TFH_global <- MeanVarPlot(seuratObj_TFH_global, fxn.x = expMean, fxn.y = logVarDivMean, x.low.cutoff = 0.25, x.high.cutoff = 3, y.cutoff = 0.5, do.contour = FALSE)
length(seuratObj_TFH_global@var.genes)
seuratObj_TFH_global <- PCA(seuratObj_TFH_global, pc.genes = seuratObj_TFH_global@var.genes)
seuratObj_TFH_global <- ProjectPCA(seuratObj_TFH_global)
VizPCA(seuratObj_TFH_global, 1:2)



