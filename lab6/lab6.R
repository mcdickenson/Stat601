# Stat 601 Lab 6
# Metropolis-Hastings
# Matt Dickenson

# Set up workspace
library(tikzDevice)
setwd('~/desktop/stat601/lab6')

# 1. Find normalizing constant

integrand <- function(x){ exp(-0.5 * (x^2)) + 0.5*(exp(-0.5*(x-3)^2)) }
integrate(integrand, lower=-Inf, upper=Inf) # approx 3.76

# 2. Implement Metropolis-Hastings

theta0 <- 3
sig.cand <- 3
n.iters <- 5000

results <- matrix(NA, nrow=n.iters, ncol=1)
results[1] <- theta0

for(i in 2:n.iters){
	last.val <- results[i-1]
	p.last <- exp(-0.5 * (last.val^2)) + 0.5*(exp(-0.5*(last.val-3)^2))
	samp <- rnorm(1, mean=last.val, sd=sig.cand)
	p.samp <- exp(-0.5 * (samp^2)) + 0.5*(exp(-0.5*(samp-3)^2))
	acceptance.prob <- min(1, p.samp/p.last)

	result <- sample(c(last.val, samp), size=1, prob=c(1-acceptance.prob, acceptance.prob) )
	results[i] <- result 
}

save(results, file='results.rda')
load('results.rda')

# desired acceptance probability about .45

same = 0 
for(i in 2:nrow(results)){
	tmp <- ifelse(results[i,]==results[i-1,], 1, 0)
	same = same + tmp
}
same / nrow(results)

# 3. Plot analytic and sampled densitities

tikz('traceplot.tex', height=4, width=4, standAlone=F )
plot(results, type='l', xlab='Metropolis-Hastings Iteration')
dev.off()

x <- seq(-4,6,by=.01)
y <- density(results)

scalar <- 3*sqrt(pi/2)

tikz('density.tex', height=4, width=4, standAlone=F )
plot(x, integrand(x), type='l', ylab='Density')
lines(y$x, y$y*scalar, lty=2)
dev.off()

# 4. Questions...