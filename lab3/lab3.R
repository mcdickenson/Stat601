# Stat 601
# Lab 3 - Normals + Gammas = T
# Matt Dickenson

library(tikzDevice)

# Suppose x|(tau^2) ~ N(0, 1/(tau^2)) and tau^2 ~ Gamma(shape=nu/2, rate=nu/2). Find the marginal distribution of x: p(x). Show ALL work.

# TODO

# Let nu=1. Get a sample of 10,000 from marginal distr. of x by drawing 10,000 tau^2's and then 10,000 x's given the tau^2's. Plot sample (either histogram or density is fine). Give two names for the actual marginal distribution p(x) when nu=1.

nu = 1
nsims = 10000

tausq = rgamma(nsims, shape=nu/2, rate=nu/2) # scale = 2/nu
x = rnorm(10000, 0, sd=sqrt(1/tausq))

tikz('density.tex', width=4, height=3, standAlone=F)
plot(density(x), main='Density of x given tau')
dev.off()

name1="t(df=1)"
name2="special case of generalized hyberbolic dist"
name3='normal-gamma dist'

# Use Kolmogorov-Smirnov test (ks.test in R) to test whether your observed distribution is equal to a t(df=1). Report p-value. What is the conclusion of the test?
ks.test(x, "pt", df=1)
# fail to reject null hypothesis that they're from the same distribution

pvals = c()
# Now, repeat ks.test 1000 times, using 100 draws from p(x) each time (instead of 10,000 draws as above). Record p-value at each iteration. (Do not report, but this will be used for the next step. The p-value can be grabbed using this R code: " ks.test(x,'pt',1)$p ".) Plot histogram of 1000 p-values and include this in report. What distribution should this be? Hint: it's a beta(a,b) for some a,b in {1,2,3,...,}.
for(ii in c(1:1000)){
	tmptausq = rgamma(100, shape=nu/2, rate=nu/2)
	tmpx = rnorm(100, 0, sd=sqrt(1/tmptausq))
	pvals[ii] = ks.test(tmpx, 'pt', df=1)$p
}

tikz('hist.tex', width=4, height=3, standAlone=F)
hist(pvals, main='Histogram of p values from KS test')
dev.off()


# Does the Central Limit Theorem hold for the mean of a sample from p(x) when nu=1? What about nu=2? nu=3? Why or why not? A quick explanation will do; an involved proof is NOT required.
nu2 = 2

pvals2 = c()
for(ii in c(1:1000)){
	tmptausq = rgamma(100, shape=nu2/2, rate=nu2/2)
	tmpx = rnorm(100, 0, sd=sqrt(1/tmptausq))
	pvals2[ii] = ks.test(tmpx, 'pt', df=nu2)$p
}
hist(pvals2)


nu3 = 3

pvals3 = c()
for(ii in c(1:1000)){
	tmptausq = rgamma(100, shape=nu3/2, rate=nu3/2)
	tmpx = rnorm(100, 0, sd=sqrt(1/tmptausq))
	pvals3[ii] = ks.test(tmpx, 'pt', df=nu3)$p
}
hist(pvals3)

# Yes, CLT continues to hold

tikz('nus.tex', width=4, height=3, standAlone=F)
par(mfrow=c(1,2))
hist(pvals2, main='Nu=2', xlab='p values of KS test')
hist(pvals3, main='Nu=3', xlab='p values of KS test')
dev.off()

