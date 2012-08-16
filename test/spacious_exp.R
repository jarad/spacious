# test out spacious

# unload the package
detach("package:spacious", unload=TRUE)

# load the package
require(spacious)

require(geoR)
set.seed(311)

# generate data to use for fitting a block composite model
n <- 200
np <- 1   # number to predict

# generate spatial locations S
S <- cbind(runif(n+np), runif(n+np))

# distance matrix
D <- rdist(S); D[row(D)==col(D)] <- 0

# matern covariance function parameters
nugget <- 0.5
tau2 <- 0.5
range <- 1.5
smooth <- 0.50  # 0.5 = exponential cov
mu <- 5
#Sigma <- nugget * diag(n) + tau2 * matern(D, range, smooth)
Sigma <- nugget * diag(n+np) + tau2 * exp(-range * D)

# generate data
y <- mvrnorm(1, mu=rep(mu, n+np), Sigma=Sigma)

# fit with spacious
X <- matrix(1, nrow=length(y), ncol=1)
x1 <- rnorm(n+np)
time.spacious <- proc.time()
#fit.spacious <- spacious(y, X, S, cov="exp", nblocks=1^2)
fit.spacious <- spacious(y~x2, data=data.frame(y=y[1:n], x2=x1[1:n]), S=S[1:n,], cov="exp", nblocks=2^2)
time.spacious <- proc.time() - time.spacious
beta.spacious <- fit.spacious$beta
theta.spacious <- fit.spacious$theta

cat("Spacious estimates:",beta.spacious,theta.spacious,"\n")
cat("Spacious execution time:\n")
print(time.spacious)

preds <- predict(fit.spacious, newdata=data.frame(x2=x1[(n+1):(n+np)]), newS=S[(n+1):(n+np),])
cat("Predictions:\n")
print(preds)
cat("Actual:\n")
print(y[(n+1):(n+np)])

if (0) {
# try likfit
gd <- as.geodata(cbind(S[1:n,],y[1:n]))
time.likfit <- proc.time()
fit.likfit <- likfit(gd, ini.cov.pars=c(1.22,1.22))
time.likfit <- proc.time() - time.likfit

beta.likfit <- fit.likfit$beta
theta.likfit <- c(fit.likfit$tausq, fit.likfit$sigmasq, 1/fit.likfit$phi)

cat("likfit estimates:",beta.likfit,theta.likfit,"\n")
cat("likfit execution time:\n")
print(time.likfit)
}
