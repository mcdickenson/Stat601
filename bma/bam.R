# Exploring BMA

#install.packages('BMA')
library('BMA')

# Ex. 1
### logistic regression 
library("MASS") 
data(birthwt) 
y<- birthwt$lo
x<- data.frame(birthwt[,-1]) 
x$race<- as.factor(x$race) 
x$ht<- (x$ht>=1)+0 
x<- x[,-9]
x$smoke <- as.factor(x$smoke) 
x$ptl<- as.factor(x$ptl) 
x$ht <- as.factor(x$ht) 
x$ui <- as.factor(x$ui)

glm.out.FT<- bic.glm(x, y, strict = FALSE, 
	OR = 20, glm.family="binomial", factor.type=TRUE)
summary(glm.out.FT) 
imageplot.bma(glm.out.FT)
glm.out.FF<- bic.glm(x, y, strict = FALSE, 
	OR = 20, glm.family="binomial", factor.type=FALSE)
summary(glm.out.FF) 
imageplot.bma(glm.out.FF)
glm.out.TT<- bic.glm(x, y, strict = TRUE, 
	OR = 20, glm.family="binomial", factor.type=TRUE)
summary(glm.out.TT) 
imageplot.bma(glm.out.TT)
glm.out.TF<- bic.glm(x, y, strict = TRUE, 
	OR = 20, glm.family="binomial", factor.type=FALSE)
summary(glm.out.TF) 
imageplot.bma(glm.out.TF)