# Stat 601
# Homework 2
# Matt Dickenson

# set up workspace
rm(list=ls())
library(tikzDevice)

# set params
y = 2; n = 10; c = 1; d = 5; a = 0.1; b = 0.3
p = y/n
data = c(rep(1, y), rep(0, n-y))

# create prior distributions
probs = seq(0,1, length=1000)
regprior = dbeta(probs, c, d)
truncprior= (dbeta(probs, c, d) * as.numeric(probs>a & probs<b))/( 
  pbeta(b,c,d) - pbeta(a,c,d) )

# compute likelihood and posteriors
likelihood = prod(dbinom(data, 1, p))
regpost = (regprior * likelihood)/(sum(regprior * likelihood))
truncpost = (truncprior * likelihood)/(sum(truncprior * likelihood))

output = matrix(NA, nrow=2, ncol=4)
colnames(output) = c('mean', 'median', 'lower95', 'upper95')
row.names(output) = c('reg', 'trunc')

qtrunc <- function(p, spec, a = -Inf, b = Inf, ...)
{
  tt <- p 
  G <- get(paste("p", spec, sep = ""), mode = "function") 
  Gin <- get(paste("q", spec, sep = ""), mode = "function") 
  tt <- Gin(G(a, ...) + p*(G(b, ...) - G(a, ...)), ...) 
return(tt)
}

output[,'mean'] = c(sum(probs*regpost), sum(probs*truncpost))
output[,'median'] = c(qbeta(.5, a+y, b+n-y), 
  qtrunc(.5, spec="beta", a=, b=b, shape1=(c+y), shape2=(d+n-y)))
output[,'lower95'] = c(qbeta(.025, a+y, b+n-y), 
  qtrunc(.025, spec="beta", a=a, b=b, shape1=(c+y), shape2=(d+n-y)))
output[,'upper95'] = c(qbeta(.975, a+y, b+n-y), 
  qtrunc(.975, spec="beta", a=a, b=b, shape1=(c+y), shape2=(d+n-y)))

print(output)