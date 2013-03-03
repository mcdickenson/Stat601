# Stat 601
# HW 6: MCMC with Gibbs sampling
# Matt Dickenson

lambda <- 5
nu <- 0.5
N <- 100

phis <- rgamma(N, shape=nu/2, scale=nu/2)

y <- rpois(N, phis*lambda)

lamba <- matrix(nrow=N, ncol=1)

set.seed(123)
for(i in 1:N){
	lam.a <- sum(y[1:i])
	lam.b <- sum(phis[1:i])
	lambda[i] <- rgamma(1, lam.a, lam.b)
}

plot(lambda, type='l')
abline(h=mean(lambda), col='red') # mean estimate
abline(h=quantile(lambda,.975), col='green')
abline(h=quantile(lambda,.025), col='green')
abline(h=5, col='blue') # true mean
