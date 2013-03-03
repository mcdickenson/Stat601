# Stat 601
# Lab 1
# Matt Dickenson

### beta density example
xx <- seq(0,1,length.out=1000)
yy <- dbeta(xx, 4,7)
plot(xx, yy, type='l')

# alternately
curve(dbeta(x,4,7))

### Simulation
n = 30
p = .7
nSims = 10000

cicalc <- function(n, p, nSims){
	
  sims <- rbinom(nSims, n, p)
  #plot(density(sims))

  # matrix for confidence intervals
  freq.ci = matrix(NA, nrow=nSims, ncol=4)
  unibayes.ci = matrix(NA, nrow=nSims, ncol=4)
  betabayes.ci = matrix(NA, nrow=nSims, ncol=4)
  colnames(freq.ci) <- colnames(unibayes.ci) <- colnames(betabayes.ci) <- c('low', 'hi', 'length', 'capture')

  freq.ci[,'low'] <- mean(sims) - 1.96*sd(sims)
  freq.ci[,'hi'] <- mean(sims) + 1.96*sd(sims)
  freq.ci[,'length'] <- freq.ci[,'hi'] - freq.ci[,'low']
  freq.ci[,'capture'] <- ifelse(sims>=freq.ci[,'low'], 1, 0)*ifelse(sims<=freq.ci[,'hi'], 1, 0)

  unibayes.ci[,'low'] <- qbeta(.025, sims, n-sims)
  unibayes.ci[,'hi'] <- qbeta(.975, sims, n-sims)
  unibayes.ci[,'length'] <- unibayes.ci[,'hi'] - unibayes.ci[,'low']
  unibayes.ci[,'capture'] <- ifelse(.7>=unibayes.ci[,'low'], 1, 0)*ifelse(.7<=unibayes.ci[,'hi'], 1, 0)

  betabayes.ci[,'low'] <- qbeta(.025, sims+5-1, n-sims+3-1)
  betabayes.ci[,'hi'] <- qbeta(.975, sims+5-1, n-sims+3-1)
  betabayes.ci[,'length'] <- betabayes.ci[,'hi'] - unibayes.ci[,'low']
  betabayes.ci[,'capture'] <- ifelse(.7>=betabayes.ci[,'low'], 1, 0)*ifelse(.7<=betabayes.ci[,'hi'], 1, 0)

  lengths <- c(mean(freq.ci[,'length'])/(n*p), mean(unibayes.ci[,'length']),mean(betabayes.ci[,'length']))

  coverage <- c(mean(freq.ci[,'capture']), mean(unibayes.ci[,'capture']), mean(betabayes.ci[,'capture']))

  output <- matrix(c(lengths, coverage), byrow=T, nrow=2, ncol=3)
  colnames(output) <- c('freq','uniform','beta')
  row.names(output) <- c('length', 'coverage')

  return(output)

}

cicalc(n=30, p=.7, nSims=10000)

#             freq   uniform      beta
# length   0.469335 0.3126386 0.2982443
# coverage 0.930100 0.9525000 0.9739000

cicalc(n=5, p=.7, nSims=10000)

#             freq  uniform     beta
# length   1.145184 0.727276 0.688167
# coverage 0.969500 0.972400 0.997100
