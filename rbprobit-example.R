## ## rbprobitGibbs example ## 

require(bayesm)

if(nchar(Sys.getenv("LONG_TEST")) != 0) {R=2000} else {R=10}

set.seed(123)
simprobit = function(X, beta){
	y = ifelse((X%*%beta+rnorm(nrow(X)))<0, 0, 1)
	list(X=X, y=y, beta=beta)
}

nobs=200
X = cbind(rep(1,nobs), runif(nobs), runif(nobs))
beta = c(0,1,-1)
nvar = ncol(X)
simout = simprobit(X, beta)

Data1 = list(X=simout$X, y=simout$y)
Mcmc1 = list(R=R, keep=1)

out = rbprobitGibbs(Data=Data1, Mcmc=Mcmc1)

summary(out$betadraw, tvalues=beta)

plot(out$betadraw)

# see also:
# rhierBinLogit(Data, Prior, Mcmc)