# http://www4.stat.ncsu.edu/~reich/st740/Logistic.R

post<-function(y,x,beta,prior.mn,prior.sd){
    #full joint posterior beta|y for the logistic regression model:
    xbeta<-x%*%beta
    xbeta<-ifelse(xbeta>10,10,xbeta)
    like<-prod(dbinom(y,1,exp(xbeta)/(1+exp(xbeta))))
    prior<-prod(dnorm(beta,prior.mn,prior.sd))
like*prior}

Bayes.logistic<-function(y,x,prior.mn=0,prior.sd=10,n.samples=10000,can.sd=1){
    #MCMC code for the model:
    #y[i]~Bern(p[i])
    #Logit(p[i])=x%*%beta
    #beta[j]~N(prior.mn,prior.sd)

    p<-ncol(x)

    #Initial values:
    beta<-rnorm(p,0,1)

    keep.beta<-matrix(0,n.samples,p)
    acc<-rep(0,p)

    for(i in 1:n.samples){

     #Update beta using MH sampling:
     for(j in 1:p){
       canbeta<-beta
       canbeta[j]<-rnorm(1,beta[j],can.sd)      #Draw candidate:
       R<-post(y,x,canbeta,prior.mn,prior.sd)/  #Compute acceptance ratio:
          post(y,x,beta,prior.mn,prior.sd)  
       U<-runif(1)                          
       if(U<R){                                 #Accept the candidate w/ prob min(R,1)
         beta<-canbeta
         acc[j]<-acc[j]+1
       }
     }
     keep.beta[i,]<-beta

     #Plot the results thus far:
     if(i%%500==0){
       par(mfrow=c(p,1))
       for(j in 1:p){plot(keep.beta[1:i,j],type="l")}
     }

   }
   #Return the posterior sample beta and the MC acceptance rates acc.rate:
list(beta=keep.beta,acc.rate=acc/n.samples)}

#Generate data:
set.seed(2008)
n<-100
x<-cbind(1,rnorm(n))
true.beta<-c(0,0.5)
true.p<-exp(x%*%true.beta)/(1+exp(x%*%true.beta))
y<-rbinom(n,1,true.p)

#fit the model
fit<-Bayes.logistic(y,x,n.samples=50000)

#Display results:
par(mfrow=c(1,2))
hist(fit$beta[1000:10000,1],main="Intercept",xlab=expression(beta[0]),breaks=50)
hist(fit$beta[1000:10000,2],main="Slope",xlab=expression(beta[1]),breaks=50)

print("Posterior mean/sd")
print(round(apply(fit$beta[1000:10000,],2,mean),3))
print(round(apply(fit$beta[1000:10000,],2,sd),3))

#Compare with glm:
print(summary(glm(y~x[,2],family="binomial")))

