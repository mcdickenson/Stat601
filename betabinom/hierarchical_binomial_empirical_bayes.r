# This implements an empirical Bayes analysis of the 
# hierarchical beta-binomial model of Section 5.1 of 
# Gelman et al (2004). Written by Marco A. R. Ferreira.
# This software is for instructional purposes only.

y <- c(rep(0,14),rep(1,8),rep(2,9),1,5,2,5,3,2,7,7,3,3,2,9,10,rep(4,7),10,4,4,4,5,11,12,5,5,6,5,6,6,6,6,16,15,15,9,4)

n <- c(rep(20,7),19,19,19,19,18,18,17,20,20,20,20,19,19,18,18,25,24,23,rep(20,6),10,49,19,46,27,17,49,47,20,20,13,48,50,rep(20,7),48,19,19,19,22,46,49,20,20,23,19,22,20,20,20,52,47,46,24,14)

plot(n,y)

plot(n,y/n)

logpygivenab <- function(n,y,alpha,beta)
{lgamma(n+1)+lgamma(alpha+beta)+lgamma(alpha+y)+lgamma(beta+n-y)-lgamma(y+1)-lgamma(n-y+1)-lgamma(alpha)-lgamma(beta)-lgamma(alpha+beta+n)}
	

loglikab <- function(n,y,alpha,beta)
{sum(logpygivenab(n,y,alpha,beta))}

alpha <- seq(from=0.1,to=10,by=0.1)
beta  <- seq(from=0.1,to=50,by=0.1)

loglik <- matrix(0,nrow=length(alpha),ncol=length(beta))
for(l in 1:length(alpha)) for(m in 1:length(beta))
  loglik[l,m] <- loglikab(n,y,alpha[l],beta[m])

contour(alpha,beta,loglik,xlab="a",ylab="b",font.lab=5,levels=c(-154.2,-155,-160,-165,-170,-175,-180))

contour(alpha,beta,exp(loglik),xlab="a",ylab="b",font.lab=5,main="Likelihood function")

#max(loglik) = 4184.582

index <- matrix(1:(length(alpha)*length(beta)),nrow=length(alpha))

index[1:10,1:10]

ind <- index[loglik==max(loglik)]

alphahat <- alpha[ind - length(alpha) * floor((ind-1)/length(alpha))]
alphahat

betahat <- beta[1+floor((ind-1)/length(alpha))]
betahat

th <- seq(0,1,0.01)

# Three different posterior densities:
# For theta_1:
plot(th,dbeta(th,alphahat+y[1],betahat+n[1]-y[1]),type="l")
# For theta_2:
lines(th,dbeta(th,alphahat+y[40],betahat+n[40]-y[40]),lty=2)
# For theta_3:
lines(th,dbeta(th,alphahat+y[70],betahat+n[70]-y[70]),lty=3)

# Now let's compare the posterior density for theta_1
# based on the EB analysis of the hierarchical model
# with the posterior density based on a beta-binomial model 
# with a uniform prior for theta_1 and considering only 
# sample 1.

# Density for theta_1 based on the hierarchical model:
plot(th,dbeta(th,alphahat+y[1],betahat+n[1]-y[1]),type="l",ylim=c(0,20),ylab="",xlab="q",font.lab=5)
# Density for theta_1 based on a model considering only sample 1:
lines(th,dbeta(th,1+y[1],1+n[1]-y[1]),lty=3)

# Density for theta_40 based on the hierarchical model:
plot(th,dbeta(th,alphahat+y[40],betahat+n[40]-y[40]),type="l",ylim=c(0,10),ylab="",xlab="q",font.lab=5)
# Density for theta_1 based on a model considering only sample 1:
lines(th,dbeta(th,1+y[40],1+n[40]-y[40]),lty=3)

# Density for theta_70 based on the hierarchical model:
plot(th,dbeta(th,alphahat+y[70],betahat+n[70]-y[70]),type="l",ylim=c(0,10),ylab="",xlab="q",font.lab=5)
# Density for theta_1 based on a model considering only sample 1:
lines(th,dbeta(th,1+y[70],1+n[70]-y[70]),lty=3)

plot(n,y/n)

postmeanth <- (alphahat+y)/(alphahat+betahat+n)

plot(y/n,postmeanth,xlab="Observed proportion",ylab="Posterior mean",main="Empirical Bayes")

postvarth <- (alphahat+y)*(alphahat+betahat+n)/((alphahat+betahat+n)^2 * (alphahat+betahat+n+1))

plot(y/n,postvarth,xlab="Observed proportion",ylab="Posterior variance",main="Empirical Bayes",ylim=c(0,0.008))

plot(n,postvarth,xlab="Sample size",ylab="Posterior variance",main="Empirical Bayes",ylim=c(0,0.008))

