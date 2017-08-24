

################  Pooled sorted and unsorted data analysis   ##########################

pbmc_data <- get(load("../output/pbmc_68K_evan.rda"))
sorted_data <- get(load("../output/10X_genomics_evan_data.rda"))

counts <- sorted_data
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

genes4 <- genes3[-idx2]

rownames(counts4) <- genes4

common_genes <- intersect(rownames(counts4), rownames(pbmc_data))

counts5 <- counts4[match(common_genes, rownames(counts4)),]
pbmc_data_filtered <- as.matrix(pbmc_data[match(common_genes, rownames(pbmc_data)),])

pooled_data <- cbind(counts5, pbmc_data_filtered)
save(pooled_data, file = "../output/pooled_sorted_unsorted_68k.rda")
