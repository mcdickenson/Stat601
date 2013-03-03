theta0 <- 3
sig.cand <- 3
n.iters <- 5000

results <- matrix(NA, nrow=n.iters, ncol=1)
results[1] <- theta0

for(i in 2:n.iters){
	last.val <- results[i-1]
	p.last <- exp(-0.5 * (last.val^2)) + 
		0.5*(exp(-0.5*(last.val-3)^2))
	samp <- rnorm(1, mean=last.val, sd=sig.cand)
	p.samp <- exp(-0.5 * (samp^2)) + 0.5*(exp(-0.5*(samp-3)^2))
	acc.prob <- min(1, p.samp/p.last)

	result <- sample(c(last.val, samp), size=1, 
		prob=c(1-acc.prob, acc.prob) )
	results[i] <- result 
}