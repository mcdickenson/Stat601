# Stat 601: Lab 9 
# Matt Dickenson

setwd('~/desktop/stat601/lab9')
library(tikzDevice)

### Pt 1: write your own Dirichlet function
mydirichlet <- function(a=c(1,1,1), k=1){
	DIRI <- matrix(rgamma(length(a)*k, a, 1), 
		nrow=k, ncol=length(a), byrow=T)
	DIRI/apply(DIRI,1,sum)
}

alpha <- c(6,1,6)
z <- mydirichlet(a=alpha, k=10)
print(z)

#            [,1]       [,2]      [,3]
#  [1,] 0.5241640 0.03148794 0.4443481
#  [2,] 0.5244181 0.07143341 0.4041485
#  [3,] 0.3085845 0.14148405 0.5499314
#  [4,] 0.3546736 0.02359379 0.6217326
#  [5,] 0.4817295 0.02394988 0.4943206
#  [6,] 0.3604205 0.10845767 0.5311218
#  [7,] 0.4758997 0.04773215 0.4763682
#  [8,] 0.4169551 0.05733824 0.5257067
#  [9,] 0.6367292 0.08229599 0.2809748
# [10,] 0.3731531 0.14817660 0.4786703

save(z, file='z.rda')


### Pt. 2: prior elicitation

alpha <- c(1,1,1)
a.over.b <- 0.007
a.over.b.squared <-  0.000007

b <- a.over.b/a.over.b.squared # 1000
a <- a.over.b*b # 7 

# prior for theta_j: gamma(7, 1000)
# prior probability:
theta1 <- rgamma(10000, shape=7, rate=1000)
lambda1 <- theta1^-1
ypred1 <- rgamma(10000, shape=10, scale=lambda1)
plot(density(ypred1))

tikz('priort1.tex', height=4, width=4, standAlone=F)
plot(density(ypred1), main='prior predictive distribution \n of batch 1 longevity')
dev.off()

mean(ypred1>5000) # 0.0105 -> one percent

save.image('lab9.rda')
load('lab9.rda')


### Pt. 3: posterior computation

# data
Y <- c(3049.1702 , 2851.0109 , 3600.8916 , 2711.999 , 
	598.8806 , 3136.3589 , 2494.1068 , 653.5904 , 2929.5264 , 
	1948.5126 , 2435.2759 , 2430.166 , 2781.2692 , 3426.5447 , 
	1484.42 ,	1565.0709 , 1557.5712 , 3269.3378 , 2175.1284 , 
	2436.5775 , 2624.7031 , 3388.7005 , 2428.5287 , 4607.2945 , 
	2244.8497 , 5456.8745 , 2174.0161 , 3643.2164 , 3833.7689 , 
	1684.4714 , 3304.7266 , 2449.3057 , 2353.6666 , 3210.8027 , 
	3688.4642 , 1982.096 , 2143.1695 , 2795.4445 , 1211.8098 , 
	2119.1033 , 1622.4898 , 1944.0038, 3485.9889 , 2789.3708 , 
	3402.7669 , 876.2526 , 1057.3501 , 2809.6497 , 2005.777 , 
	2526.5537)


k <- 3 # number of thetas 
M <- 10500 # num iterations 
pis <- thetas <- matrix(NA, nrow=M, ncol=k)

# hyperparameters
alpha <- c(6, 1, 6) # for dirichlet
pis[1, ] <- alpha/sum(alpha) 
thetas[1,] <- rep(0.007, k) # set to prior mean
tmp <- rep(0, k)
N <- rep(0, k) # how many datapoints belong to each component
n <- rep(0, length(Y))
gshape <- 10

z <- matrix(NA, nrow=length(Y), ncol=k) 
# prob that each datapoint belongs to given component

a <- rep(7, k) # gamma hyperparams to be updated
b <- rep(1000, k)

# Gibbs sampler
for(m in 2:M){
	# step 0 - set lambdas from thetas
	lambda <- 1/thetas[m-1, ]

	# step 1 - sample 1 proportion vector from Dirichlet
	pisamp <- mydirichlet(c(alpha[1]+N[1], alpha[2]+N[2], alpha[3]+N[3]), 1)
	pis[m, ] <- pisamp

	# step 2 - sample latent indicators
	for(i in 1:nrow(z)){
		for(j in 1:ncol(z)){
			z[i, j] = pisamp[1,j] * dgamma(Y[i], shape=gshape, scale=lambda[j])
		}
		z[i,] <- z[i,]/sum(z[i,])

		# step 3 - calculate new n's based on z draws 
		n[i] <- sample(c(1,2,3), size=1, prob = z[i,])
	}

	# step 4 - update thetas
	for(i in 1:k){
		N[i] = length(which(n==i))
		a[i] = a[i] + N[i]
		sum.y.i <- sum(Y[which(n==i)])
		b[i] = b[i] + sum.y.i 
		tmp[i] <- rgamma(1, a[i], b[i])
	}

	thetas[m, ] <- sort(tmp)

}

# Display Gibbs sampler results
tikz('thetas.tex', height=4, width=4, standAlone=F)
plot(thetas[500:1500,1], type='l', ylim=c(0,0.0008),
	main='Gibbs Iterations of Thetas',
	xlab='Index (minus 500 burn-in iterations)',
	ylab='thetas')
lines(thetas[500:1500,2], type='l', col='blue')
lines(thetas[500:1500,3], type='l', col='red')
legend("topright", c('theta1', 'theta2', 'theta3'), col=c('black', 'blue', 'red'), pch=16)
dev.off()

tikz('pis.tex', height=4, width=4, standAlone=F)
plot(pis[500:1500,1], type='l', ylim=c(0,1),
	main='Gibbs Iterations of Pis',
	xlab='Index (minus 500 burn-in iterations)',
	ylab='pis')
lines(pis[500:1500,2], type='l', col='blue')
lines(pis[500:1500,3], type='l', col='red')
legend("topright", c('pi1', 'pi2', 'pi3'), col=c('black', 'blue', 'red'), pch=16)
dev.off()

save.image('lab9all.rda')

# question about pis
mean(pis[,3]>=0.60) # 0.3774286
