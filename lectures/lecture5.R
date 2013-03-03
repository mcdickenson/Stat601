# R code to generate stuff for lecture 5

# plot draws from a Poisson process on [0,1] with a couple different rate functions

tht = c(2,4,10)
y1 = rpois(1,tht[1])
y2 = rpois(1,tht[2])
y3 = rpois(1,tht[3])
t1 = runif(y1)
t2 = runif(y2)
t3 = runif(y3)

plot(t1,rep(0,y1),xlim=c(0,1),ylim=c(-0.1,0.1),xlab="event times", ylab="")
abline(h=0)

plot(t2,rep(0,y2),xlim=c(0,1),ylim=c(-0.1,0.1),xlab="event times", ylab="")
abline(h=0)
title("Event Times from Poisson Process (rate=4)")

# For this Poisson process the density of the gap times in between events is exponential(1/4)
y = seq(0,3,length=1000)
py = dexp(y,4)
plot(y,py,type="l",xlab="gap times",ylab="density")
abline(v=1/4)
title("Exponential density with rate 4 (mean 1/4)")

plot(t3,rep(0,y3),xlim=c(0,1),ylim=c(-0.1,0.1),xlab="event times", ylab="")
abline(h=0)
title("Event Times from Poisson Process (rate=10)")

y = seq(0,3,length=1000)
py = dexp(y,10)
plot(y,py,type="l",xlab="gap times",ylab="density")
abline(v=1/10)
title("Exponential density with rate 10 (mean 1/10)")

# plot the normal density 
y = seq(-3,3,length=1000)
plot(y,dnorm(y,1,.2),type="l",xlab="y",ylab="normal density")
lines(y,dnorm(y,0,1))
lines(y,dnorm(y,-1,.5))
title("normal densities")




