
################    skeleton for gamlr    ###########################

n <- 1000
p <- 3
xvar <- matrix(0.9, nrow=p,ncol=p)
diag(xvar) <- 1
x <- matrix(rnorm(p*n), nrow=n)%*%chol(xvar)
y <- 4 + 3*x[,1] + -1*x[,2] + rnorm(n)

fitlasso <- gamlr(x, y, gamma=0, lambda.min.ratio=1e-3) # lasso
fitgl <- gamlr(x, y, gamma=2, lambda.min.ratio=1e-3) # small gamma


family="gaussian"
gamma=0
nlambda=100
lambda.start=Inf
lambda.min.ratio=0.01
free=NULL
standardize=TRUE
obsweight=NULL
varweight=NULL
prexx=(p<500)
tol=1e-7
maxit=1e5
verb=FALSE


famid = switch(family,
          "gaussian"=1, "binomial"=2, "poisson"=3)

n <- length(y)

################    eta represents the shift/offset in the model   ###################

eta <- rep(0.0,n)
if(!is.null(xtr$shift)){
  if(family=="gaussian") y = y-xtr$shift
  else eta <- xtr$shift   }
stopifnot(length(eta)==n)
eta <- as.double(eta)

###############    observation weights   #########################

if(!is.null(obsweight))
  if(family!="gaussian"){
    warning("non-null obsweight are ignored for family!=gaussian")
    obsweight <- NULL }
if(is.null(obsweight)) obsweight <- rep(1,n)
stopifnot(all(obsweight>0))
stopifnot(length(obsweight)==n)

# fixedcost (undocumented: additional fixed l1 penalty)

if(is.null(xtr$fixedcost))
  xtr$fixedcost <- 0
fixedcost = xtr$fixedcost
if(length(fixedcost)!=p){
  fixedcost <- rep(fixedcost[1],p) }


## variable (penalty) weights
if(is.null(varweight)) varweight <- rep(1,p)
stopifnot(all(varweight>=0))
stopifnot(length(varweight)==p)
varweight[free] <- 0

## check and clean all arguments
stopifnot(lambda.min.ratio<=1)
stopifnot(all(c(nlambda,lambda.min.ratio)>0))
stopifnot(all(c(lambda.start)>=0))
stopifnot(all(c(tol,maxit)>0))
if(lambda.start==0){
  nlambda <- 1
  standardize <- 0 }
lambda <- double(nlambda)
lambda[1] <- lambda.start

delta <- exp( log(lambda.min.ratio)/(nlambda-1) )

xbar <- double(p)
vxsum <- double(p)
vxy <- double(p)
vxx <- double(0)


###################   C code   ####################################


famid=as.integer(famid)
n=n
p=p
l=length(x@i)
xi=x@i
xp=x@p
xv=as.double(x@x)
y=y
prexx=as.integer(prexx)
xbar=xbar
xsum=vxsum
xx=vxx
xy=vxy
eta=eta
varweight=as.double(varweight)
obsweight=as.double(obsweight)
standardize=as.integer(standardize>0)
nlambda=as.integer(nlambda)
delta=as.double(delta)
gamma=gamvec
fixedcost=as.double(fixedcost)
tol=as.double(tol)
maxit=as.integer(rep(maxit,nlambda))
maxrw=as.integer(rep(maxrw,nlambda))
lambda=as.double(lambda)
deviance=double(nlambda)
df=double(nlambda)
alpha=as.double(rep(0,nlambda))
beta=as.double(rep(0,nlambda*p))
exits=integer(nlambda)
verb=as.integer(verb>0)
