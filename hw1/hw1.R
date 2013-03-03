nBins = 200

# set up prior to have high density near 0
prior = rpois(nBins+1, .3)
prior=sort(prior, decreasing=T)
probs = seq(0,1,by=1/nBins)
names(prior) = probs
prior = prior/sum(prior)

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

# compute posteriors
posterior = (prior * likelihoods)/sum(prior*likelihoods, na.rm=T)
