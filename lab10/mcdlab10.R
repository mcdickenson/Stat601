# Stat 601 - Lab 10
# Matt Dickenson

setwd('~/desktop/stat601/lab10')
library(tikzDevice)
data <- read.csv('data.csv')

# Truncated normal sampler:
rtnorm<-function(n,mu,sigma,lower,upper){ 
   lp<-pnorm(lower,mu,sigma) 
   up<-pnorm(upper,mu,sigma)  
   qnorm(runif(n,lp,up),mu,sigma) 
}

# Gibbs sampling function for probit model
myprobit<-function(y,x,sd.beta=10, iters=11000){
	n <- length(y)
  p <- ncol(x)
  low <- ifelse(y==1,0,-Inf)
  high <- ifelse(y==1,Inf,0)

  # Initial values
  ystar <- rep(0.5, n)
  beta <- rep(0,p)

  # Matrices for saving samples
  keep.beta<-matrix(0,iters,p)
  keep.ystar<-matrix(0,iters,10)

  # Matrix manipulations
  txx<-t(x)%*%x
  cov.beta<-solve(txx+diag(p)/sd.beta^2)
  P1<-cov.beta%*%t(x)
  P2<-t(chol(cov.beta))

  # Gibbs iterations
  for(i in 1:iters){

    # Update the latent probit variables, ystar:
    mn <- x%*%beta
    ystar <- rtnorm(n,mn,1,low,high)

    # Update beta:
    beta<-P1%*%ystar+P2%*%rnorm(p)

    keep.beta[i,]<-beta
    keep.ystar[i,]<-ystar[1:10]
   
  }

list(beta=keep.beta,
     ystar=keep.ystar)

}

myx <- rep(1, nrow(data))
myx <- cbind(myx, as.matrix(data[,1:3]))
colnames(myx)[1] <- 'intercept'

model <- myprobit(data$y, myx)

save.image('lab9.rda')

dim(model$beta)


# Posterior prediction
simx <- matrix(c(1, 26,0,1), nrow=1)
colnames(simx) <- colnames(myx)

xbeta <-  model$beta[1001:11000, ] %*% t(simx)
summary(xbeta)

y.hat.probs <- pnorm(xbeta)

tikz(file='postpred.tex', height=4, width=4, standAlone=F)
hist(y.hat.probs, breaks=20, main='Posterior Prediction', xlab='Pr(accident)')
dev.off()

y.hat <- ifelse(xbeta>0, 1, 0)
mean(y.hat) # 0.0246 

# Decision problem 
younger <- matrix(c(1, 17, 1, 0), nrow=1)
older <- matrix(c(1, 18, 1, 1), nrow=1)

youngerxb <-  model$beta[1001:11000, ] %*% t(younger)
olderxb <-  model$beta[1001:11000, ] %*% t(older)

youngerprob <- pnorm(youngerxb)
olderprob <- pnorm(olderxb)

mean(youngerprob < olderprob) # 67.1 percent chance the older driver is riskier
