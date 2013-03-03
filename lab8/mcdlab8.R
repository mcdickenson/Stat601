# Stat 601
# Lab 8: Mixture modeling
# Matt Dickenson

# setup
library(MASS)
library(gplots)
library(car)
library(MCMCpack) # for Dirichlet distribution

# data
Y <- c(3049.1702 , 2851.0109 , 3600.8916 , 2711.999 , 
	598.8806 , 3136.3589 , 2494.1068 , 653.5904 , 2929.5264 , 
	1948.5126 , 2435.2759 , 2430.166 , 2781.2692 , 3426.5447 , 
	1484.42 ,	1565.0709 , 1557.5712 , 3269.3378 , 2175.1284 , 
	2436.5775 , 2624.7031 , 3388.7005 , 2428.5287 , 4607.2945 , 
	2244.8497 , 5456.8745 , 2174.0161 , 3643.2164 , 3833.7689 , 
	1684.4714 , 3304.7266 , 2449.3057 , 2353.6666 , 3210.8027 , 
	3688.4642 , 1982.096 , 2143.1695 , 2795.4445 , 1211.8098 , 
	2119.1033 , 1622.4898 , 1944.0038, 3485.9889 , 2789.3708 , 
	3402.7669 , 876.2526 , 1057.3501 , 2809.6497 , 2005.777 , 
	2526.5537)

# priors
pivec <- c(.2, .2, .6)
priorsamp <- 5
alpha <- pivec * priorsamp
n1 <- n2 <- n3 <- 0 
n <- rep(0, 50)
z <- matrix(NA, nrow=length(Y), ncol=length(pivec))

# constants
lambda <- c(100, 200, 300)
gshape <- 10

M <- 1000 # num iterations
pimatrix <- matrix(NA, nrow=M, ncol=length(pivec))

for(m in 1:M){
	# step 1 - simulate 1 proportion vector from Dirichlet
	pisamp <- rdirichlet(1, c(alpha[1]+n1, alpha[2]+n2, alpha[3]+n3) )
	pimatrix[m, ] <- pisamp

	# step 2 - sample latent indicators
	for(i in 1:nrow(z)){
		for(j in 1:ncol(z)){
			z[i, j] = pisamp[1,j] * dgamma(Y[i], shape=gshape, scale=lambda[j])
		}
		z[i,] <- z[i,]/sum(z[i,])

		# step 3 - calculate new n's based on z draws 
		n[i] <- sample(c(1,2,3), size=1, prob = z[i,])
	}

	n1 = length(which(n==1))
	n2 = length(which(n==2))
	n3 = length(which(n==3))
}

# plot(density(pimatrix[-(1:100), 1]), xlim=c(0,1))
# lines(density(pimatrix[-(1:100), 2]), lty=2)
# lines(density(pimatrix[-(1:100), 3]), lty=3)

# describe joint distribution graphically
#plot(pimatrix[-(1:100), 2], pimatrix[-(1:100), 3])
# d <- kde2d(pimatrix[-(1:100), 2], pimatrix[-(1:100), 3])
# setwd('~/Desktop/Stat601/lab8')
# tiff(filename='credregion.tiff', width=600, height=600)
# image(d, main='Bivariate Density and Credible Region',
# 	xlab='pi2', ylab='pi3')
# contour(d, add=T) # indicates credible regions
# dev.off()
#persp(d)

setwd('~/Desktop/Stat601/lab8')
tiff(filename='credregion.tiff', width=600, height=600)
ci2d(pimatrix[-(1:100), 2], pimatrix[-(1:100), 3], show.points=T)
title(xlab="pi2", ylab="pi3", main="post burn-in points and ci")
dev.off()

pisum <- pimatrix[-(1:100), 1] + pimatrix[-(1:100), 2]
mean(pisum<.4) # 0.8455556
