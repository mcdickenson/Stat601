# Stat 601
# HW 5: Monte Carlo approximation
# Matt Dickenson


# Poisson example
a <- 2
b <- 1 
ycollege.sum <-66
ncollege <- 44

poisim <- function(a, b, ysum, n, numsims){
	theta.mc <- rgamma(numsims, a+ysum, b+n)
	y.mc <- rpois(numsims, theta.mc)
	lbl = paste('n=', numsims, sep='')
	hist(y.mc, main=lbl)
}

par(mfrow=c(2,2))
poisim(a=a, b=b, ysum=ycollege.sum, n=ncollege, numsims=10)
poisim(a=a, b=b, ysum=ycollege.sum, n=ncollege, numsims=100)
poisim(a=a, b=b, ysum=ycollege.sum, n=ncollege, numsims=10000)
poisim(a=a, b=b, ysum=ycollege.sum, n=ncollege, numsims=100000)

# Beta/Bernoulli example
# 1. Exact posterior is a beta(2+10,2+20)

curve(dbeta(x, 2+10, 2+20))

# 2. 
a=12
b=22

p1 <- rbeta(10, a, b)
summary(p1) 
# Min.    1st Qu.  Median  Mean    3rd Qu.    Max. 
#  0.2423  0.2919  0.3461  0.3556  0.4259  0.4672 

p2 <- rbeta(100, a, b)
summary(p2) 
# Min.    1st Qu.  Median    Mean  3rd Qu.    Max. 
#  0.1348  0.2962  0.3429  0.3519  0.4055  0.5825 

p3 <- rbeta(10000, a, b)
summary(p3) 
# Min.    1st Qu.  Median  Mean    3rd Qu.    Max. 
#  0.1032  0.2937  0.3487  0.3522  0.4074  0.6474 

# 3. 
y <- rbinom(10000, 30, p3)
plot(density(y))
summary(y)/30
# Min.     1st Qu.    Median     Mean     3rd Qu.   Max. 
# 0.0333333 0.2666667 0.3333333 0.3506667 0.4333333 0.8333333 

# in comparison to the beta, this density has:
# a lighter left tail
# a lower Median
# a heavier right tail 