
require(geoR)
set.seed(311) # TODO: should we keep this seed set?

# generate data to use for fitting a block composite model
n  <- 250
np <- 5   # number to predict

# generate spatial locations S
S <- cbind(runif(n+np), runif(n+np))

# distance matrix
D <- rdist(S); D[row(D)==col(D)] <- 0

# Covariance function parameters common to exponential and matern
nugget <- 0.5
tau2   <- 0.5
range  <- 1.5
mu     <- 5

# Explanatory variables
x <- rnorm(n+np)

####################### Exponential covariance function ###########################
Sigma  <- nugget * diag(n+np) + tau2 * exp(-range * D)

# generate data
y <- mvrnorm(1, mu=rep(mu, n+np), Sigma=Sigma) # TODO: should this include the explanatory variable?
# fit with spacious
sd <- data.frame(y=y[1:n], x=x[1:n])
time.spacious <- system.time(
                 fit.spacious <- spacious(y~x, sd, S=S[1:n,], cov="exp", nblocks=2^2))

# fit with geoR
gd <- as.geodata(cbind(S[1:n,],y[1:n])) # TODO: where is the explanatory variable?
time.likfit <- system.time(fit.likfit <- likfit(gd, ini.cov.pars=c(0.5,2), messages=F))


# Comparison
comp = data.frame(spacious=c(time.spacious[3], fit.spacious$beta[-2], fit.spacious$theta),
                  likfit  =c(time.likfit[3], fit.likfit$beta, fit.likfit$tausq, fit.likfit$sigmasq, 1/fit.likfit$phi))
rownames(comp) = c("time","beta","nugget","partial sill","range")

cat("\nComparison of spacious with geoR::likfit using exponential covariance function\n")
print(comp)

####################### Matern covariance function ###############################
smooth <- 1.50  
Sigma <- nugget * diag(n+np) + tau2 * matern(D, 1/range, smooth)
y <- mvrnorm(1, mu=rep(mu, n+np), Sigma=Sigma)

# fit with spacious
sd <- data.frame(y=y[1:n], x=x[1:n])
time.spacious <- system.time(
                 fit.spacious <- spacious(y~x, sd, S=S[1:n,], cov="exp", nblocks=2^2))

# fit with geoR
gd <- as.geodata(cbind(S[1:n,],y[1:n])) # TODO: where is the explanatory variable?
time.likfit <- system.time(fit.likfit <- likfit(gd, ini.cov.pars=c(0.5,2), messages=F))


# Comparison
comp = data.frame(spacious=c(time.spacious[3], fit.spacious$beta[-2], fit.spacious$theta),
                  likfit  =c(time.likfit[3], fit.likfit$beta, fit.likfit$tausq, fit.likfit$sigmasq, 1/fit.likfit$phi))
rownames(comp) = c("time","beta","nugget","partial sill","range")

cat("\nComparison of spacious with geoR::likfit using matern covariance function\n")
print(comp)



# TODO: what else should we be comparing? standard errors?

#preds <- predict(fit.spacious, newdata=data.frame(x2=x1[(n+1):(n+np)]), newS=S[(n+1):(n+np),], interval="prediction", level=0.9)
#cat("Predictions:\n")
#print(preds)
#cat("Actual:\n")
#print(y[(n+1):(n+np)])

