# Stat 601
# HW9: Finite mixture models
# Matt Dickenson

# setup
library(MASS)
library(MCMCpack) # for Dirichlet distribution

# simulate data
Y <- 0.25*rnorm(100, 0, 0.5) + 0.75*rnorm(100, 2, 1)

### Gibbs sampler 

# first size of k
k <- 3 # or 10
M <- 5000
pis <- mus <- taus <- matrix(NA, nrow=M, ncol=k)
pis[1, ] <- rep(1/k, k)
mus[1,] <- seq(-1, 1, length.out=k)
taus[1,] <- rep(1, k)

S <- matrix(NA, nrow=length(Y), ncol=k)
S[1, ] <- 


for(m in 2:M){

	# update pr S_i = h from multinomial conditional posterior
	for(i in 1:length(Y)){
		for(h in 1:k){
			like <- prod( dnorm(Y ,mus[m-1, h], 1/taus[m-1, h]) )
			S[i, h] <- pivec[m-1, h] * like 
		}
		S[i, ] <- S[i, ]/sum(S[i, ]) # normalize
	}

	# update mu_h, tau^-1_h from conditional posterior

	for(h in 1:k){
		kappahat[h] <- (kappa^{-1} + n[h] )^{-1} # what is kappa? what is n[h]?
		muhat[h] <- kappahat[h] * ( kappa^{-1} * mu0 + n[h]*ybar[h] ) # mu0? ybar[h]? 
		ahat[h] <- a[h] + n[h]/2
		bhat[h] <- b[h] + 0.5 * { # BIG SUM HERE  }
		mus[m, h] <- rnorm(1, muhat[h], (kappahat*tau^-1) ) 
		taus[m, h] <- rgamma(1, ahat, bhat)
	}
	
	# update pi
	for(h in 1:k){
		pis[m, ] <- # Dirichlet
	}

}

# second size of k 
k <- 10

### estimate the density of the data

### examine and describe trace plots

### Increase the variance of the conjugate prior P0 greatly & repeat
