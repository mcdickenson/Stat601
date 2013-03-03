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

# create distributions
probs = seq(0,1, length=1000)

regprior = dbeta(probs, c, d)
truncprior= (dbeta(probs, c, d) * as.numeric(probs>a & probs<b))/( pbeta(b,c,d) - pbeta(a,c,d) )

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
output[,'median'] = c(qbeta(.5, a+y, b+n-y), qtrunc(.5, spec="beta", a=, b=b, shape1=(c+y), shape2=(d+n-y)))
output[,'lower95'] = c(qbeta(.025, a+y, b+n-y), qtrunc(.025, spec="beta", a=a, b=b, shape1=(c+y), shape2=(d+n-y)))
output[,'upper95'] = c(qbeta(.975, a+y, b+n-y), qtrunc(.975, spec="beta", a=a, b=b, shape1=(c+y), shape2=(d+n-y)))

print(output)

xtrunc = seq(round(output['trunc','lower95'], 3), round(output['trunc','upper95'], 3), by=.001) 
ytrunc = truncpost[xtrunc*1000]

xreg = seq(round(output['reg','lower95'], 3), round(output['reg','upper95'], 3), by=.001) 
yreg = regpost[xreg*1000]

tikz('posteriors.tex', height=3, width=4, standAlone=F)
plot(probs, regpost, type='l', lwd=2, col='blue',
  main="Posterior Distributions", ylab="Density", 
  xlab=paste('$','\\theta','$', sep=''),
  ylim=c(0,.009))
lines(probs, truncpost, lwd=2, col='red')
lines(rep(output['reg','mean'],2), c(0,.006), col='blue',lty=2, lwd=2)
lines(rep(output['trunc','mean'],2), c(0,.006), col='red', lty=2, lwd=2)
polygon( c(xtrunc, rev(xtrunc)), c(ytrunc, rep(0,length(ytrunc))), col=rgb(255,0,0,100,maxColorValue=255), border=NA)
polygon( c(xreg, rev(xreg)), c(yreg, rep(0,length(yreg))), col=rgb(0,0,255,100,maxColorValue=255), border=NA)
legend('topright', legend=c('Truncated Prior', 'Unconstrained Prior'), fill=c('red', 'blue'), cex=.6, text.width=.33)
dev.off()
