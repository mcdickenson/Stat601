# Stat 601
# Exam 1 Prep
# Matt Dickenson

### Hoff Exercise 3.4

n <- 43; y <- 15

a <- 2; b <- 8

# A 

# prior
curve(dbeta(x, a, b), main='prior (2, 8)')

# likelihood
thetas <- seq(0,1, by=0.01)
y.given.theta <- (thetas^y) * (1-thetas)^(n-y)
plot(thetas, y.given.theta, type='l')

# posterior
curve(dbeta(x, a+y, b+n-y), main='posterior (2, 8)')

# B 

a <-8; b <- 2

curve(dbeta(x, a, b), main='prior (8, 2)')

# posterior
curve(dbeta(x, a+y, b+n-y), main='posterior (8, 2)')

# C 

p.theta <- 0.25 * ( gamma(10) / (gamma(2) * gamma(8))) * (3*thetas * (1-thetas)^7 + (thetas^7)*(1-thetas))
plot(thetas, p.theta, type='l', main='prior mixture')


# D 
library(Bolstad)
postmix <- binomixp(y, n, alpha0=c(2,8), alpha1=c(2,8), p=0.75, ret=T)


y.lims<-c(0,1.1*max(postmix$posterior,postmix$prior))
plot(postmix$pi,postmix$posterior,ylim=y.lims,type="l"
        ,xlab=expression(pi),ylab="Density",main="Posterior")




### Hoff exercise 4.2

# A: obtain pr(th_b < th_a | y_a, y_b) from prior in 3.3

y.a <- c(12,9,12,14,13,13,15,8,15,6)
y.b <- c(11,11,10,9,9,8,7,10,6,8,8,9,7)

s.a <- 120
r.a <- 10

s.b <- 12
r.b <- 1

theta.a <- rgamma(100000, s.a + sum(y.a), r.a + length(y.a))
theta.b <- rgamma(100000, s.b + sum(y.b), r.b + length(y.b))

mean(theta.b < theta.a )


# B 
for(n in c(1:12)){
	s.b <- 12*n
	r.b <- n
	theta.a <- rgamma(100000, s.a + sum(y.a), r.a + length(y.a))
	theta.b <- rgamma(100000, s.b + sum(y.b), r.b + length(y.b))
	out <- mean(theta.b < theta.a )
	print(c(n, out))
}

# The larger the prior sample size we put on theta.b,
# The less certain we are that the posterior is less than theta.a


# C 
s.a <- 120
r.a <- 10

s.b <- 12
r.b <- 1

theta.a <- rgamma(100000, s.a + sum(y.a), r.a + length(y.a))
theta.b <- rgamma(100000, s.b + sum(y.b), r.b + length(y.b))

samp.a <- rpois(100000, theta.a)
samp.b <- rpois(100000, theta.b)

mean(samp.b < samp.a)

for(n in c(1:12)){
	s.b <- 12*n
	r.b <- n
	theta.a <- rgamma(100000, s.a + sum(y.a), r.a + length(y.a))
	theta.b <- rgamma(100000, s.b + sum(y.b), r.b + length(y.b))
	samp.a <- rpois(100000, theta.a)
	samp.b <- rpois(100000, theta.b)
	out <- mean(samp.b < samp.a)
	print(c(n, out))
}








