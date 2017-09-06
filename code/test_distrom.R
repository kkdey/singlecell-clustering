

#######################  Test distrom()   ######################

library(distrom)
library(MASS)
data(fgl)
## make your cluster
## FORK is faster but memory heavy, and doesn't work on windows.
cl <- makeCluster(2,type=ifelse(.Platform$OS.type=="unix","FORK","PSOCK"))
print(cl)

fits <- dmr(cl, fgl[,1:9], fgl$type, verb=1)
stopCluster(cl)

#################  tracing the skeleton of distrom  ##########################

covars <- do.call(rbind, replicate(1, fgl[,1:9], simplify=FALSE))
dim(covars)

mu=NULL
bins=NULL
verb=0
cv=FALSE

cl <- makeCluster(2,type=ifelse(.Platform$OS.type=="unix","FORK","PSOCK"))
print(cl)


counts <- fgl$type

argl <- list()
argl$family <- "poisson"
if(is.null(argl$nlambda))
  argl$nlambda <- formals(gamlr)$nlambda
argl$verb <- max(verb-1,0)
argl$cv <- cv
## collapse and clean
counts <- matrix(sample(1:1000, 214*10000, replace=TRUE), 214, 10000)
colnames(counts) <- paste("Gene:", 1:dim(counts)[2])
chk <- collapse(covars, counts, mu, bins)

if(verb)
  cat(sprintf("fitting %d observations on %d categories, %d covariates.\n",
              nrow(chk$v), ncol(chk$counts), ncol(chk$v)))
argl$x <- chk$v
argl$shift <- chk$mu
p <- ncol(chk$counts)
vars <- colnames(chk$counts)
## cleanup
rownames(argl$x) <- rownames(chk$counts) <- NULL
counts <- chk$counts
rm(covars,mu,chk)

C <- ifelse(is.null(cl),Inf,length(cl))
if(C < p/4){
  chunks <- round(seq(0,p,length.out=C+1))
  counts <- lapply(1:C,
                   function(i) counts[,(chunks[i]+1):chunks[i+1]])
  counts <- parLapply(cl,
                      counts,
                      function(x)
                        sapply(colnames(x),
                               function(j) x[,j,drop=FALSE]))
  counts <- unlist(counts,recursive=FALSE)
} else{
  counts <- sapply(vars,
                   function(j) counts[,j,drop=FALSE]) }



## lapply somehow, depending on cl and p
if(is.null(cl)){
  if(verb) cat("running in serial.\n ")
  mods <- lapply(counts,onerun,argl=argl)
}else{
  if(verb){
    cat("distributed run.\n")
    print(cl) }
  mods <- parLapply(cl,counts,onerun,argl=argl)
}

stopCluster(cl)

class(mods) <- "dmr"
attr(mods,"nobs") <- argl$nobs
attr(mods,"nlambda") <- argl$nlambda
attr(mods,"mu") <- argl$shift

coef.dmr(mods)
