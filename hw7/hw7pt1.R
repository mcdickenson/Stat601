# Stat 601
# Matt Dickenson
# HW 7: Multivariate Normal Distribution

library(MASS)
library(MCMCpack) # for Wishart
library(tikzDevice)

setwd('~/Desktop/Stat601/hw7/')

# 1. simulate data
s = matrix(c(1, 0.8, 0.8, 1), nrow=2, byrow=T)
samp <- mvrnorm(100, mu=c(0,0), Sigma=s)
hist(samp)

# 2a. estimate mle
mle.mu <- c(mean(samp[,1]), mean(samp[,2]))
mle.sigma <- cov(samp)

# > mle.mu
# [1] -0.1074701 -0.1235791
# > mle.sigma
#           [,1]      [,2]
# [1,] 1.0222084 0.7445078
# [2,] 0.7445078 0.9207407

# 2b. contour plot (code from Dr. Dunson)

# exact
x=seq(-3,3, length=100); y <- x; rho <- .8
bivnd <- function(x,y){ exp(-(x^2 - 2*rho *x *y +y^2)/(2*(1- rho^2))) }
z <- x %*% t(y)	# to set up z as a matrix of the right size 
for (i in 1:100){
     for (j in 1:100){ z[i,j] <- bivnd(x[i],y[j])
} }
tikz('contour1.tex', height=3, width=3, standAlone=F)
contour(x,y,z,xlab="y1",ylab="y2") 
title("bivariate normal", cex.main=0.9)
dev.off()

# sample
y1 <- sort(samp[,1]); y2 <- sort(samp[,2]); rho2 <- cor(samp)[1,2]
bivnd2 <- function(x,y){ exp(-(x^2 - 2*rho2 *x *y +y^2)/(2*(1- rho2^2))) }
z2 <- y1 %*% t(y2)	# to set up z as a matrix of the right size 
for (i in 1:100){
     for (j in 1:100){ z2[i,j] <- bivnd2(y1[i],y2[j])
} }
tikz('contour2.tex', height=3, width=3, standAlone=F)
contour(y1,y2,z2,xlab="y1",ylab="y2", xlim=c(-3,3), ylim=c(-2,3)) 
title("sample data", cex.main=0.9)
dev.off()

# 3. Gibbs sampler
mu0 <- c(1,1)
L0 <- matrix(c(1,0.5,0.5,1), nrow=2, byrow=T)

nu0 <- 4
S0 <- matrix(c(1,0.5,0.5,1), nrow=2, byrow=T)
n <- length(y1) 
ybar <- apply(samp, 2, mean)
Sigma <- cov(samp)
THETA <- SIGMA <- NULL



numIters <- 5000
for(s in 1:numIters){

	# sample theta
	Ls <- solve ( solve(L0) + n*solve(Sigma))
	mus <- Ls %*% (solve(L0)%*%mu0 + n*solve(Sigma)%*%ybar)
	theta <- mvrnorm(1,mus,Ls)

	# sample sigma
	Ss <- S0 + ( t(samp)-c(theta))%*%t( t(samp)-c(theta))
	#Sigma <- rwishart(nu0+n, solve(Ss))$IW
	#Sigma <- riwish(nu0+n, solve(Ss))
	Sigma <- solve(rwish(1, nu0+n, solve(Ss)))

	# record results
	THETA <- rbind(THETA, theta)
	SIGMA <- rbind(SIGMA, c(Sigma))

}


# bivariate traceplot of means
plot(THETA[,1], THETA[,2], type='o', col='grey5')

# 4. compare Bayes, MLE, truth
mean(THETA[10:5000,1]) 
mean(THETA[10:5000,2]) 
mean(SIGMA[10:5000,1])
mean(SIGMA[10:5000,2])
mean(SIGMA[10:5000,3])
mean(SIGMA[10:5000,4])


