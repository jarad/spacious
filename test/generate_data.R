# generate data to use for fitting a block composite model

require(MASS)
require(fields)
require(geoR)

n <- 200

nugget <- 0.5
range <- 0.5   # range
smooth <- 0.5
tau2 <- 0.5   # partial sill
mu <- 20

set.seed(311)

# generate spatial locations S
S <- cbind(runif(n), runif(n))

# get distance matrix
D <- rdist(S)
D[row(D)==col(D)] <- 0

# get covariance matrix
#Sigma <- nugget * diag(n) + tau2 * matern(D, range, smooth)
Sigma <- nugget * diag(n) + tau2 * exp(-range * D)

# generate data
y <- mvrnorm(1, mu=rep(mu, n), Sigma=Sigma)

gd <- as.geodata(cbind(S,y))  # greate geodata object

#print(round(S,3))
#print(round(D,3))
#print(round(Sigma,3))
#print(summary(Y))
#print(summary(gd))

#plot(S[,1], S[,2], xlim=c(0,1), ylim=c(0,1))
#plot(gd)
#plot(variog(gd,max.dist=1))

#ml <- likfit(gd, ini=c(1,0.5))
#print(ml)
#reml <- likfit(gd, ini=c(1,0.5))
#print(ml)
