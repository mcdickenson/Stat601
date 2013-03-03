# Stat 601
# Lab 2 - Gamma-Poisson conjugacy
# Matt Dickenson

# Plot both priors on the same plot. 
prior2x <- c(3,4,4,5,5,6,6,7)
prior2y <- c(.07,.07,.45,.45,.39,.39,.09,.09)

# Now, suppose y|theta ~ Pois(theta). Observe:
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

# What is the 95% central credible interval?
qgamma(c(.025,.975), shape=a.star, scale=b.star)
# (3.98, 5.92)


