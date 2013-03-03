# Stat 601
# HW9: Finite mixture models
# Matt Dickenson

# setup
library(MCMCpack) # for Dirichlet distribution

# simulate data
Y <- 0.25*rnorm(100, 0, 0.5) + 0.75*rnorm(100, 2, 1)

# Run 1: k=3
k <- 3

# hyperparameters
M <- 5000
pis <- mus <- taus <- matrix(NA, nrow=M, ncol=k)
pis[1, ] <- rep(1/k, k)
alpha <- rep(1, k)
mus[1,] <- seq(-1, 1, length.out=k)
taus[1,] <- c(2, 4, 6)^-1
khat <- rep(0, k)
N <- rep(0, k)

S <- matrix(NA, nrow=length(Y), ncol=k)
S[1,] <- alpha/k

n <- rep(0, length(Y))
a <- b <- rep(0, k)

for(m in 2:M){

	# UPDATE PI

	# step 1 - simulate 1 proportion vector from Dirichlet
	pisamp <- rdirichlet(1, c(alpha[1]+N[1], alpha[2]+N[2], alpha[3]+N[3]) )
	pis[m, ] <- pisamp

	# step 2 - sample latent indicators
	for(i in 1:nrow(S)){
		for(j in 1:ncol(S)){
			S[i, j] = pisamp[1,j] * dnorm(Y[i], mus[m-1, j], taus[m-1, j])
		}
		S[i,] <- S[i,]/sum(S[i,])

		# step 3 - calculate new n's based on z draws 
		n[i] <- sample(c(1,2,3), size=1, prob = S[i,])
	}


	# UPDATE MU, TAU

	for(h in 1:k){
		N[h] = length(which(n==h))
		y.bar.h <- ifelse(N[h]==0, 1e-6, mean(Y[which(n==h)]) )
		khat[h] <- (1/k + N[h])^-1
		muhat <- khat[h]*(k^-1 * mus[m-1,h] + N[h] * y.bar.h)
		a[h] <- N[h]/2 
		a[h] <- ifelse(a[h]==0, 1e-6, a[h])
		b[h] <- 0.5*( sum( (Y-y.bar.h)^2 ) + (N[h]/ (1+ k*N[h]) ) *( y.bar.h-mus[m-1, h])^2 )
		b[h] <- ifelse(is.finite(b[h]), b[h], 1e-6)

		taus[m, h] <- rgamma(1, shape=a[h], rate=b[h]) 
		taus[m, h] <- ifelse(taus[m, h]<1e-6, 1e-6, taus[m, h])
		mus[m, h] <- rnorm(1, muhat, khat[h]*(1/taus[m,h]) )
	}
		
}

plot(mus[, 1], type='l') 
plot(mus[, 2], type='l') 
plot(mus[, 3], type='l') 

plot(taus[, 1], type='l')
plot(taus[, 2], type='l')
plot(taus[, 3], type='l')






# Run 2: k=10
k <- 10

# hyperparameters
M <- 5000
pis <- mus <- taus <- matrix(NA, nrow=M, ncol=k)
pis[1, ] <- rep(1/k, k)
alpha <- rep(1, k)
mus[1,] <- seq(-3, 3, length.out=k)
taus[1,] <- c(2, 3, 4, 2, 3, 4, 2, 3, 4, 2)^-1
khat <- rep(0, k)
N <- rep(0, k)

S <- matrix(NA, nrow=length(Y), ncol=k)
S[1,] <- alpha/k

n <- rep(0, length(Y))
a <- b <- rep(0, k)

for(m in 2:M){

	# UPDATE PI

	# step 1 - simulate 1 proportion vector from Dirichlet
	pisamp <- rdirichlet(1, alpha+N )
	pis[m, ] <- pisamp

	# step 2 - sample latent indicators
	for(i in 1:nrow(S)){
		for(j in 1:ncol(S)){
			S[i, j] = pisamp[1,j] * dnorm(Y[i], mus[m-1, j], taus[m-1, j])
		}
		S[i,] <- S[i,]/sum(S[i,])

		# step 3 - calculate new n's based on z draws 
		n[i] <- sample(1:k, size=1, prob = S[i,])
	}


	# UPDATE MU, TAU

	for(h in 1:k){
		N[h] = length(which(n==h))
		y.bar.h <- ifelse(N[h]==0, 1e-6, mean(Y[which(n==h)]) )
		khat[h] <- (1/k + N[h])^-1
		muhat <- khat[h]*(k^-1 * mus[m-1,h] + N[h] * y.bar.h)
		a[h] <- N[h]/2 
		a[h] <- ifelse(a[h]==0, 1e-6, a[h])
		b[h] <- 0.5*( sum( (Y-y.bar.h)^2 ) + (N[h]/ (1+ k*N[h]) ) *( y.bar.h-mus[m-1, h])^2 )
		b[h] <- ifelse(is.finite(b[h]), b[h], 1e-6)

		taus[m, h] <- rgamma(1, shape=a[h], rate=b[h]) 
		taus[m, h] <- ifelse(taus[m, h]<1e-6, 1e-6, taus[m, h])
		mus[m, h] <- rnorm(1, muhat, khat[h]*(1/taus[m,h]) )
	}
		
}

par(mfrow=c(5,2))
for(i in 1:k){
	plot(mus[, i], type='l')
}

# two of the clusters have reasonably non-zero mus

for(i in 1:k){
	plot(taus[, i], type='l')
}

# eighth cluster has a reasonably non-zero tau

setwd('~/Desktop/stat601/hw9')
save.image('hw9.rda')
