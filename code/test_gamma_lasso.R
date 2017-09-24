
###############  gamma lasso testing mechanism   #######################

n <- 1000
p <- 3
xvar <- matrix(0.9, nrow=p,ncol=p)
diag(xvar) <- 1
x <- matrix(rnorm(p*n), nrow=n)%*%chol(xvar)

##########   Normal regression model  ###################

y <- 4 + 3*x[,1] + -1*x[,2] + rnorm(n)

## run models to extra small lambda 1e-3xlambda.start
fitlasso <- gamlr(x, y, gamma=0, lambda.min.ratio=1e-3) # lasso
fitgl <- gamlr(x, y, gamma=2, lambda.min.ratio=1e-3) # small gamma
fitglbv <- gamlr(x, y, gamma=10, lambda.min.ratio=1e-3) # big gamma


#############  Poisson regression model   ##########################


z <- sapply(y, function(x) return(rpois(1, exp(x))))

fitlasso <- gamlr(x, z, family = "poisson", gamma=0, lambda.min.ratio=1e-3) # lasso
fitgl <- gamlr(x, z, family = "poisson", gamma=2, lambda.min.ratio=1e-3) # small gamma
fitglbv <- gamlr(x, z, family = "poisson", gamma=10, lambda.min.ratio=1e-3) # big gamma

#############  Binomial regression model  #########################

w <- sapply(y, function(x) return (rbinom(1, 1, inv.logit(x))))

fitlasso <- gamlr(x, w, family = "binomial", gamma=0, lambda.min.ratio=1e-3) # lasso
fitgl <- gamlr(x, w, family = "binomial", gamma=2, lambda.min.ratio=1e-3) # small gamma
fitglbv <- gamlr(x, w, family = "binomial", gamma=10, lambda.min.ratio=1e-3) # big gamma


par(mfrow=c(1,3))
ylim = range(c(fitglbv$beta@x))
plot(fitlasso, ylim=ylim, col="navy")
plot(fitgl, ylim=ylim, col="maroon")
plot(fitglbv, ylim=ylim, col="darkorange")




data(hockey)
x <- cBind(config,team,player)
y <- goal$homegoal
## fit the plus-minus regression model
## (non-player effects are unpenalized)
fit <- gamlr(x, y, gamma=10, lambda.min.ratio=0.1,
             free=1:(ncol(config)+ncol(team)),
             standardize=FALSE, family="binomial")
plot(fit)
beta <- coef(fit)

mu <- beta[1,1] + x%*% beta[-1,]
inv_mu <- inv.logit(as.vector(mu))
w <- sapply(inv_mu, function(x) return (rbinom(1, 1, x)))
fitlasso <- gamlr(x, w, family = "binomial", gamma=0, lambda.min.ratio=1e-3) # lasso

mean(abs(coef(fitlasso) -  beta))

order(coef(fitlasso), decreasing = TRUE)[1:30]
order(beta, decreasing = TRUE)[1:30]



