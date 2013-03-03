# Stat 601
# Lab 2 - Gamma-Poisson conjugacy
# Matt Dickenson

library(tikzDevice)

# 1. Gamma(shape=50, scale=.1) 
# 2. For theta in (3,4], p(theta)=0.07; in (4,5], p(theta)=0.45; in (5,6], p(theta)=0.39; in (6,7], p(theta)=0.09. 

# Plot both priors on the same plot. (I want to see this.)


prior2x <- c(3,4,4,5,5,6,6,7)
prior2y <- c(.07,.07,.45,.45,.39,.39,.09,.09)

tikz('priors.tex', height=3, width=4, standAlone=F)
plot(prior2x, prior2y, pch=c(16,1), 
	xlim=c(0,8), ylim=c(0,1),
	xlab=expression('$\\theta$'),
	ylab=expression('$p(\\theta)$'),
	main='Priors' )
points(3,0)	
points(7,0)
points(7,.09,pch=16)
for(ii in c(1,3,5,7)){
	lines(c(prior2x[ii], prior2x[ii+1]), rep(prior2y[ii], 2))
}
lines(density(rgamma(1000, shape=50, scale=.1)), col='blue', lwd=2)
dev.off()

# Now, suppose y|theta ~ Pois(theta). Suppose further that you observe the following data (n = 10): 
# y = (2, 1, 9, 4, 3, 3, 7, 7, 5, 7).
y = c(2,1,9,4,3,3,7,7,5,7)
n = length(y)

a.star = 50 + sum(y)
b.star = (n+ 1/.1)^-1

priors2p = c(0,0,.07,.45,.39,.09)
xvec = seq(3,7,by=.001)

# compute numerator for posteriors
calc.post <- function(xx){
	k = floor(xx)
	pk = priors2p[k]
	temp = pk * (exp(-n*xx)) * (xx^sum(y))
	return(temp)
}

post2p = sapply(xvec, FUN=calc.post)
post2p[1] <- 0
post2p[4001] <- 0 

# compute denominator for posteriors
scalar = 0 
for(ii in c(3:6)){
	temp = pgamma(ii+1,shape=a.star,scale=b.star) - pgamma(ii,shape=a.star,scale=b.star) 
	temp = temp * priors2p[ii]
	scalar = scalar + temp
}
scalar = scalar * ((b.star^a.star) * gamma(a.star))

post2p = post2p / scalar

# Plot the two posteriors on the same plot (but a different plot from the first). (I want to see this.)
tikz('posteriors.tex', height=3, width=4, standAlone=F)
plot(xvec, post2p*.25*(10^13), pch=19,
	xlim=c(0,8), ylim=c(0,1),
	xlab=expression('$\\theta$'),
	ylab=expression('$p(\\theta)$'),
	main='Posteriors' )
lines(density(rgamma(1000, shape=a.star, scale=b.star)), 
	lwd=2, col='blue')
dev.off()

# What is the 95% central credible interval (i.e., using 2.5th and 97.5th percentiles) for posterior 1 (the posterior derived from prior 1)? Describe in words how one would find the 95% central credible interval for posterior 2 (but you don't need to find it). 
qgamma(c(.025,.975), shape=a.star, scale=b.star)
# (3.98, 5.92)

# Find this interval, and you get a 10% bonus on this assignment (if Prof. Dunson lets me give one)! 

calc.post2 <- function(xx){
	k = floor(xx)
	pk = priors2p[k]
	temp = pk * (exp(-n*xx)) * (xx^sum(y))
	return(temp)
}

integrate(calc.post2, 3, 7)

# TODO: divid by normalizing constant

