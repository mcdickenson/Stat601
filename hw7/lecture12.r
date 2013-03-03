# R code to show bivariate Gaussian contours 

# plot the multivariate normal density without sampling
# http://www.statslab.cam.ac.uk/~pat/misc.pdf - useful link to multivariate plotting & analysis in R
x=seq(-3,3, length=60); y <- x; rho <- .7
bivnd <- function(x,y){ exp(-(x^2 - 2*rho *x *y +y^2)/(2*(1- rho^2))) }
z <- x %*% t(y)	# to set up z as a matrix of the right size 
for (i in 1:60){
     for (j in 1:60){ z[i,j] <- bivnd(x[i],y[j])
} }
contour(x,y,z,xlab="y1",ylab="y2") 
title("contour plot of bivariate normal (rho=0.7, scale=1, mean=0)")

image(x,y,z,xlab="y1",ylab="y2") 
title("heat plot of bivariate normal (rho=0.7, scale=1, mean=0)")

persp(x,y,z,xlab="y1",ylab="y2",zlab="f(y1,y2)")
title("3d plot of bivariate normal (rho=0.7, scale=1, mean=0)")

# Above our moderately positively correlated Gaussians - repeat for a variety 

# Spherical Gaussians
rho <- 0
bivnd <- function(x,y){ exp(-(x^2 - 2*rho *x *y +y^2)/(2*(1- rho^2))) }
z <- x %*% t(y)	# to set up z as a matrix of the right size 
for (i in 1:60){
     for (j in 1:60){ z[i,j] <- bivnd(x[i],y[j])
} }
contour(x,y,z,xlab="y1",ylab="y2") 
title("contour plot of bivariate normal (rho=0, scale=1, mean=0)")

# Highly negatively correlated Gaussians
rho<- -.9
bivnd <- function(x,y){ exp(-(x^2 - 2*rho *x *y +y^2)/(2*(1- rho^2))) }
z <- x %*% t(y)	# to set up z as a matrix of the right size 
for (i in 1:60){
     for (j in 1:60){ z[i,j] <- bivnd(x[i],y[j])
} }
contour(x,y,z,xlab="y1",ylab="y2") 
title("contour plot of bivariate normal (rho=-0.9, scale=1, mean=0)")



