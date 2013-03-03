# R code to generate stuff for lecture 1

p = seq(0,1,length=100)

pd1 = dbeta(p,1,1)
pd2 = dbeta(p,1,20)
pd3 = dbeta(p,2,10)

plot(p,pd2,type="l",xlab="theta",ylab="prior density")
lines(p,pd1,lty=2)
lines(p,pd3,lty=3)
title("Examples of prior densities for a probability")
