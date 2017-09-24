
######################  check distrom on birds   ########################


library(ecostructure)
library(covtpx)
library(Biobase)
library(distrom)

data <- get(load(system.file("extdata", "HimalayanBirdsData.rda", package = "ecostructure")))
taxonomic_counts <- t(exprs(data))
grid_metadata <- pData(phenoData(data))


covars1 <- grid_metadata[,1:3]
covars2 <- model.matrix(~ factor(grid_metadata$WorE)-1)
covars <- cbind(covars1, covars2)

cl <- makeCluster(parallel::detectCores(),type=ifelse(.Platform$OS.type=="unix","FORK","PSOCK"))
print(cl)
system.time(fits <- dmr(cl, covars[,1:3], taxonomic_counts+1, verb=1))
stopCluster(cl)

coefs_birds <- as.matrix(coef(fits))

metadata4 <- cbind(rep(1, dim(covars)[1]), covars[, 1:3])
colnames(metadata4) <- c("Intercept", colnames(covars[, 1:3]))

tmp <- exp(as.matrix(metadata4) %*% coefs_birds)
tmp2 <- apply(tmp, 1, function(x) return(x/sum(x)))

tmp3 <- t(apply(taxonomic_counts+1, 2, function(x) return(x/sum(x))))

plot(tmp2[1,], col="red", pch=20, cex = 1)
plot(tmp3[1,], col="blue", pch=20, cex = 1)

plot(covars[,1], tmp2[which(rownames(tmp2) == "Phylloscopus_pulcher"),],
     xlab = "Elevation", col="red", pch=20, cex = 1, ylim= c(0,0.3))
points(covars[,1], tmp3[which(rownames(tmp2) == "Phylloscopus_pulcher"),],
       xlab = "Elevation", col="blue", pch=20, cex = 1, ylim=c(0,0.3))
legend("topright", legend =  c("fitted", "original"),
       fill = c("red", "blue"), col = c("red", "blue"), cex = 0.5)


plot(covars[,1], tmp2[which(rownames(tmp2) == "Phylloscopus_trochiloides"),],
     xlab = "Elevation", col="red", pch=20, cex = 1, ylim= c(0,0.3))
points(covars[,1], tmp3[which(rownames(tmp2) == "Phylloscopus_trochiloides"),],
       xlab = "Elevation", col="blue", pch=20, cex = 1, ylim=c(0,0.3))
legend("topright", legend =  c("fitted", "original"),
       fill = c("red", "blue"), col = c("red", "blue"), cex = 0.5)


which(rownames(tmp2) == "Phylloscopus_pulcher")
