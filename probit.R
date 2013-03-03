#######  Code to fit the model    ########:
#y[i] =  Bern[Phi(x[i,]%*%beta)]
#x[i,] = covariates for sub i
# http://www4.stat.ncsu.edu/~reich/st740/probit.R

#Draw samples from a truncated normal:
rtnorm<-function(n,mu,sigma,lower,upper){ 
   lp<-pnorm(lower,mu,sigma) 
   up<-pnorm(upper,mu,sigma)  
   qnorm(runif(n,lp,up),mu,sigma) 
}


probit<-function(y,x,sd.beta=100,
            iters=10000,burn=1000,update=10){

  n<-length(y)
  p<-ncol(x)
  low<-ifelse(y==1,0,-Inf)
  high<-ifelse(y==1,Inf,0)

  #Initial values
  z<-y-.5
  beta<-rep(0,p)

  #store samples here
  keep.beta<-matrix(0,iters,p)
  keep.z<-matrix(0,iters,10)

  #Do some matrix manipulations offline
  txx<-t(x)%*%x
  cov.beta<-solve(txx+diag(p)/sd.beta^2)
  P1<-cov.beta%*%t(x)
  P2<-t(chol(cov.beta))

  #Let's go!
  for(i in 1:iters){

    #update the latent probit variables, z:
    mn<-x%*%beta
    z<-rtnorm(n,mn,1,low,high)

    #update beta:
    beta<-P1%*%z+P2%*%rnorm(p)

    keep.beta[i,]<-beta
    keep.z[i,]<-z[1:10]
   
  }


list(beta=keep.beta,
     z=keep.z)}

fit<-probit(y=dv,x=iv)
names(fit)
boxplot(fit$beta)
dim(fit$beta)
plot(fit$beta[,2], type='l')