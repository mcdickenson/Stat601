
lambda <- 25
jointbn <- matrix(NA, nrow=100000, ncol=2)
colnames(jointbn) <- c('b','n')
jointbn[1,] <- c(.05, 50) # set starting values

for(jj in 2:nrow(jointbn)){
	bprev <- jointbn[jj-1,'b']

	x <- rpois(1, lambda*(1-bprev)) # simulate Nhat using bprev
	Nhat <- x + 20

	bhat <- rbeta(1, 21, Nhat-19) # simulate bhat using Nhat

	jointbn[jj,] <- c(bhat, Nhat) # save in our data matrix
}

bvec <- jointbn[1:10,'b']
nvec <- jointbn[1:10,'n']

pdf(file='firstten.pdf')
plot(bvec, nvec, type='o', xlab='beta', ylab='N')
dev.off()

bfull <- jointbn[1001:nrow(jointbn),'b']
nfull <- jointbn[1001:nrow(jointbn),'n']

quantile(bfull, c(.05,.95)) # 0.55, 0.97
mean(bfull) # 0.775...

n20 <- ifelse(nfull==20, 1, 0)
mean(n20) # 0.073

n25 <- ifelse(nfull==25, 1, 0)
mean(n25) # 0.093