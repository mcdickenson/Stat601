# Lab 4 - Car Dealership Problem
# Matt Dickenson

# 1. What is probability that you exceed RA's performance?

# posterior is lambda~gamma(shape=5, rate=6)
lambdas = rgamma(100000, shape=5, rate=6)

x.week2 = rpois(100000, lambdas*5)

plot(density(x.week2), main="Posterior Predictive Distribution of x given lambda")
mean(x.week2>=6) # 0.27

# 2. What is probability that your bonus exceeds $100?

bonuses = (2*pi*x.week2)^sqrt(2)
mean(bonuses>100)
hist(bonuses, main='Predicted bonus') # 0.39

# 3. Why is it especially difficult to calculate the probability that you will sell 
# more cars in your second week than Ronald Aylmer sold in his?

# Because we are uncertain about:
# a) the true data generating process (assumed Poisson)
# b) the number of cars that RA will sell in the second week

# In other words, it's harder to calculate the probability that 
# our number of cars sold will be greater than his, since his
# is still a random variable (not a constant). 
