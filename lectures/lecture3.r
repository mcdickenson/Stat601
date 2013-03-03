# R code for lecture 

p=seq(0,1,length=100)
fp = dbeta(p, a+y, b+n-y)
plot(p,fp,type="l",xlab="theta",ylab="posterior density")
title("beta posterior (y=2,n=10,a=b=1)")


