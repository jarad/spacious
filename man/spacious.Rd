\name{spacious}
\alias{spacious}
\title{Fit block composite model}
\description{
	Fits a block composite model.
}
\usage{
	spacious(y, X, S, cov = "exp", cov.inits = NULL, nblocks = 1, grid.type = "regular")
}
\arguments{
	\item{y}{response}
	\item{X}{model matrix}
	\item{S}{n by 2 matrix of spatial locations}
	\item{cov}{spatial covariance function}
	\item{cov.inits}{initial values for covariance function}
	\item{nblocks}{number of blocks to divide data into}
	\item{grid.type}{type of blocking structure to place points into}
}
\seealso{
	\code{\link{summary.spacious}}
}
\examples{

require(geoR)
set.seed(311)

# generate data to use for fitting a block composite model
n <- 200

# generate spatial locations S
S <- cbind(runif(n), runif(n))

# distance matrix
D <- rdist(S); D[row(D)==col(D)] <- 0

# matern covariance function parameters
nugget <- 0.5
tau2 <- 0.5
range <- 1.5
smooth <- 0.50  # 0.5 = exponential cov
mu <- 5
#Sigma <- nugget * diag(n) + tau2 * matern(D, range, smooth)
Sigma <- nugget * diag(n) + tau2 * exp(-range * D)

# generate data
y <- mvrnorm(1, mu=rep(mu, n), Sigma=Sigma)

# fit with spacious
X <- matrix(1, nrow=length(y), ncol=1)
time.spacious <- proc.time()
fit.spacious <- spacious(y, X, S, cov="exp", nblocks=1^2)
time.spacious <- proc.time() - time.spacious
beta.spacious <- fit.spacious$beta
theta.spacious <- fit.spacious$theta

cat("Spacious estimates:",beta.spacious,theta.spacious,"\n")
cat("Spacious execution time:\n")
print(time.spacious)

# try likfit
gd <- as.geodata(cbind(S,y))
time.likfit <- proc.time()
fit.likfit <- likfit(gd, ini.cov.pars=c(1.22,1.22))
time.likfit <- proc.time() - time.likfit

beta.likfit <- fit.likfit$beta
theta.likfit <- c(fit.likfit$tausq, fit.likfit$sigmasq, 1/fit.likfit$phi)

cat("likfit estimates:",beta.likfit,theta.likfit,"\n")
cat("likfit execution time:\n")
print(time.likfit)

}
\keyword{spacious}
