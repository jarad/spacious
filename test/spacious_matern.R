# test out spacious

if ("package:spacious" %in% search()) {
	# unload the package
	detach("package:spacious", unload=TRUE)
}

# load the package
require(spacious)

require(geoR)
set.seed(311)

# generate data to use for fitting a block composite model
n <- 100
np <- 5   # number to predict

# generate spatial locations S
S <- cbind(runif(n+np), runif(n+np))

# distance matrix
D <- rdist(S); D[row(D)==col(D)] <- 0

# matern covariance function parameters
nugget <- 0.5
tau2 <- 0.5
range <- 0.25
smooth <- 3.50  # 0.5 = exponential cov
mu <- 5
Sigma <- nugget * diag(n+np) + tau2 * matern(D, range, smooth)
#Sigma <- nugget * diag(n+np) + tau2 * exp(-D/range)

if (0) {
nu <- smooth
    mid <- 2*D*range*sqrt(nu)
    rho <- mid^nu * besselK(mid, nu)/(2^(nu-1) * gamma(nu))
    rho[is.na(rho)] <- 1
    #theta[1] * diag(nrow(D)) + theta[2] * mid^nu * besselK(mid, nu)/(2^(nu-1) * gamma(nu))
Sigma <- nugget * diag(nrow(D)) + tau2 * rho
}

# generate data
y <- mvrnorm(1, mu=rep(mu, n+np), Sigma=Sigma)

# fit with spacious
X <- matrix(1, nrow=length(y), ncol=1)
x1 <- rnorm(n+np)
time.spacious <- proc.time()
#fit.spacious <- spacious(y, X, S, cov="exp", nblocks=1^2)
#fit.spacious <- spacious(y~x2, data=data.frame(y=y[1:n], x2=x1[1:n]), S=S[1:n,], cov="matern", smoothness=smooth, nblocks=2^2, verbose=TRUE)
fit.spacious <- spacious(y~x2, data=data.frame(y=y[1:n], x2=x1[1:n]), S=S[1:n,], cov="matern", fixed=list(smoothness=smooth), blocks=list(type="regular", nblocks=2^2), verbose=TRUE)
#fit.spacious <- spacious(y~x2, data=data.frame(y=y[1:n], x2=x1[1:n]), S=S[1:n,], cov="matern", fixed=NULL, nblocks=2^2, verbose=TRUE)
time.spacious <- proc.time() - time.spacious
beta.spacious <- fit.spacious$beta
theta.spacious <- fit.spacious$theta

print(summary(fit.spacious))

cat("Spacious estimates:",beta.spacious,theta.spacious,"\n")
cat("Spacious SEs:",fit.spacious$se.beta,fit.spacious$se.theta,"\n")
cat("Spacious execution time:\n")
print(time.spacious)

preds <- predict(fit.spacious, newdata=data.frame(x2=x1[(n+1):(n+np)]), newS=S[(n+1):(n+np),], interval="prediction", level=0.9)
cat("Predictions:\n")
print(preds)
cat("Actual:\n")
print(y[(n+1):(n+np)])

# consider using the rbenchmark package

if (0) {
# try likfit
gd <- as.geodata(cbind(S[1:n,],y[1:n]))
time.likfit <- proc.time()
fit.likfit <- likfit(gd, ini.cov.pars=c(0.5,0.5), fix.kappa=TRUE, kappa=smooth, cov.model="matern")
time.likfit <- proc.time() - time.likfit

beta.likfit <- fit.likfit$beta
theta.likfit <- c(fit.likfit$tausq, fit.likfit$sigmasq, 1/fit.likfit$phi)

cat("likfit estimates:",beta.likfit,theta.likfit,"\n")
cat("likfit execution time:\n")
print(time.likfit)
}
