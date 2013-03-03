### Stat 601
### HW 1
### Matt Dickenson
### 1 Sept 2012

set.seed(123)
library(tikzDevice)

nBins = 200

# set up prior to have high density near 0
prior = rpois(nBins+1, .3)
prior=sort(prior, decreasing=T)
probs = seq(0,1,by=1/nBins)
names(prior) = probs
prior = prior/sum(prior)

#pdf(file="prior.pdf")
tikz("prior.tex", height=3, width=4, standAlone=F)
plot(probs,prior,,type="h",
	col="blue",lwd=2,
	xlab=paste('$','\\theta','$', sep=''), ylab="Density", main="Prior")
dev.off()
# not as good as beta, better than uniform


# we now observe data
trueTheta = 0.01
sampleSize = 1000
newData = rbinom(sampleSize, 1, trueTheta)


# compute likelihoods
likelihoods = rep(NA, nBins+1)
names(likelihoods) = probs
for(ii in c(1:nBins+1)){
	likelihoods[ii] = prod(dbinom(newData, 1, probs[ii]))
}

plot(probs,likelihoods, type="l",
	xlab=bquote(theta), ylab="L", main="Likelihoods")


# compute posteriors
posterior = (prior * likelihoods)/sum(prior*likelihoods, na.rm=T)

#pdf(file="posterior.pdf")
tikz("posterior.tex", height=3, width=4, standAlone=F)
plot(probs, posterior, type="h",
	col="blue",lwd=2,
	xlab=paste('$','\\theta','$', sep=''), ylab="Density", main="Posterior")
dev.off()

# Shahryar's method: 
# p = seq(.05, .25, .01)
# weight_prior <- c(2, 2, 2, 2, 3, 3, 3, 3, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1)
# posterior_binom <- function(successes, failures, p, weight_prior){
#		likelihood <- successes * log(p) + failures * log(1 - p)
# 	likelihood <- exp(likelihood - max(likelihood))
# 	like_prior <- likelihood * weight_prior
# 	posterior <- like_prior/sum(like_prior)
# }
# weight_prior <- weight_prior/sum(weight_prior)
# n <- 100
# successes <- sample(c(seq(0, 25, 1)), 1) # 12
# failures <- n - successes