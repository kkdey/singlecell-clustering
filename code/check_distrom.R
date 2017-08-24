
####################  check distrom (Matt Taddy)  ########################

install.packages("distrom")

library(MASS)
library(distrom)
data(fgl)
cl <- makeCluster(2,type=ifelse(.Platform$OS.type=="unix","FORK","PSOCK"))
print(cl)

fits <- dmr(cl, fgl[,1:9], fgl$type, verb=1)
stopCluster(cl)

par(mfrow=c(2,2))
for(j in 1:4){
  plot(fits[[j]])
  mtext(names(fits)[j],font=2,line=2) }


B <- coef(fits)

par(mfrow=c(1,1))
P <- predict(B, fgl[,1:9], type="response")
boxplot(P[cbind(1:214,fgl$type)]~fgl$type,
        ylab="fitted prob of true class")




library(singleCellRNASeqMouseDeng2014)
deng.counts <- Biobase::exprs(Deng2014MouseESC)
deng.meta_data <- Biobase::pData(Deng2014MouseESC)
deng.gene_names <- rownames(deng.counts)


cl <- makeCluster(4,type=ifelse(.Platform$OS.type=="unix","FORK","PSOCK"))
print(cl)


fits <- dmr(cl, deng.meta_data$cell_type, t(deng.counts), verb=1)


out <- collapse(deng.meta_data$cell_type,t(deng.counts),mu=NULL,bins=NULL)
fits <- dmr(cl=NULL, deng.meta_data$cell_type, t(deng.counts), verb=1)

