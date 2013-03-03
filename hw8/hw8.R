# Stat 601
# HW8: Prediction for missing data
# Matt Dickenson


### Gibbs sampler 

library(MASS)

s = matrix(c(1, 0.8, 0.8, 1), nrow=2, byrow=T)
Y <- mvrnorm(150, mu=c(0,0), Sigma=s)

Y.oracle <- Y

Y[101:150,1] <- NA

# priors
n<-dim(Y)[1] ; p<-dim(Y)[2]
mu0 <- c(1,1)
sd0<-(mu0/2)
L0 <- matrix(c(1,.5,.5,1), nrow=2)
nu0<-p+2 ; S0<-L0

# starting vaues
Sigma<- S0
Y.full <- Y
O <- 1* (!is.na(Y))

# Gibbs sampler
THETA <- SIGMA <- Y.MISS <- NULL

for(j in 1:p){
	Y.full[ is.na(Y.full[, j]) , j] <- mean(Y.full[ , j ] ,na.rm=TRUE)
}

# From Hoff - http://www.stat.washington.edu/hoff/Book/Data/data/chapter7.r
# sample from the Wishart distribution
rwish<-function(n,nu0,S0)
{
  sS0 <- chol(S0)
  S<-array( dim=c( dim(S0),n ) )
  for(i in 1:n)
  {
     Z <- matrix(rnorm(nu0 * dim(S0)[1]), nu0, dim(S0)[1]) %*% sS0
     S[,,i]<- t(Z)%*%Z
  }
  S[,,1:n]
}

set.seed (1)

for(s in 1:1000)
{
	###update theta
	ybar<-apply(Y.full ,2 ,mean)
	Ln<-solve( solve (L0) + n*solve (Sigma) )
	mun<-Ln%*%( solve(L0)%*%mu0 + n*solve(Sigma)%*%ybar )
	theta<-mvrnorm(1 ,mun,Ln)
	###

	###update Sigma
	Sn<- S0 + ( t(Y.full)-c(theta) )%*%t( t(Y.full)-c(theta) )
	Sigma<-solve( rwish(1, nu0+n, solve(Sn)) )
	###

	###update missing data
	for(i in 101:n){
		b <- ( O[i,]==0 )
		a <- ( O[i,]==1 )
		iSa<- solve(Sigma[a,a])
		beta.j <- Sigma[b,a]%*%iSa 

		Sigma.j	<- Sigma[b,b] - Sigma[b,a] %*% iSa %*%Sigma[a,b]
		theta.j<- theta[b] + beta.j%*%(t(Y.full[i,a])-theta[a])
		Y.full[i,b]<-mvrnorm(1,theta.j,Sigma.j )
	}
	### save results
	THETA<-rbind(THETA,theta) ; SIGMA<-rbind(SIGMA,c(Sigma))
	Y.MISS<-rbind(Y.MISS, Y.full[O==0] )
	###
}


# obtain predictive interval
dim(Y.MISS)

Y.MISS.t <- t(Y.MISS)

yhat <- apply(Y.MISS.t,1, mean)
length(yhat)

quants <- function(x){
	quantile(x, probs=c(.025, .975))
}

yq <- apply(Y.MISS.t,1, quants)
ylb <- yq[1,]
yub <- yq[2,]

yest <- cbind(yhat, ylb, yub)

# estimate mean sq err
yobs <- Y.oracle[101:150,1]
y.sqerr <- (yest[,1] - yobs)^2
mean(y.sqerr) # 0.268
###


# estimate joint mle
y1mean <- mean(Y[1:100,1])
y2obs <- Y[101:150,2]
y2mean <- mean(y2obs)
vcv <- cov(Y[1:100,1:2])
vy1 <- vcv[1,1]
vy2 <- vcv[2,2]
cory1y2 <- vcv[1,2]

mle.mat <- matrix(NA, nrow=50, ncol=1000)
for(i in 1:1000){
	mle.mat[,i] <- cory1y2 * ( (y2obs - y2mean) * (vy1/vy2) * y1mean ) + rnorm(50, 0, sqrt(1-cory1y2))
}

mlest <- apply(mle.mat,1, mean)
mleq <- apply(mle.mat,1, quants)
mlelb <- mleq[1,]
mleub <- mleq[2,]

ml.sqerr <- (mlest - yobs)^2
mean(ml.sqerr) # 0.893

# estimate 'oracle model'

oracle.mat <- matrix(NA, nrow=50, ncol=1000)
for(i in 1:1000){
	oracle.mat[,i] <- 0.8 * ( y2obs  * (vy1/vy2)  ) + rnorm(50, 0, sqrt(0.2))
}

orest <- apply(oracle.mat, 1, mean)
oreq <- apply(oracle.mat, 1, quants)
orlb <- oreq[1,]
orub <- oreq[2,]

or.sqerr <- (orest - yobs)^2
mean(or.sqerr) # 0.262

# estimate regression mle
y.obs <- Y[1:100,]
names(y.obs)<-c('y1','y2')
y.obs <- as.data.frame(y.obs)
colnames(y.obs)<-c('y1','y2')
mod <- glm(y1 ~ y2, data=y.obs)

y2.obs <-Y.oracle[101:150,2]
y2.obs <- as.data.frame(y2.obs)
colnames(y2.obs)<-'y2'
y1.pred <- predict(mod, y2.obs)

y.mle.sqerr <- (y1.pred - yobs)^2
mean(y.mle.sqerr) #0.269

### Coverage

# Gibbs
covered <- ifelse(yobs>=yest[,2] & yobs<=yest[,3], 1, 0)
mean(covered) # 0.98

# Bivariate MLE
covmle <- ifelse(yobs>=mlelb & yobs<=mleub, 1, 0)
mean(covmle) # 0.74

# Oracle model
covoracle <- ifelse(yobs>=orlb & yobs<=orub, 1, 0)
mean(covoracle) # 0.94

# save data
setwd('~/desktop/stat601/hw8')
save.image('hw8.RData')
load('hw8.RData')