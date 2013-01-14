\name{spacious}
\alias{spacious}
\title{Fit block composite model}
\description{
Fits a block composite model.
}
\usage{
spacious(formula, data, S, cov = "exp", cov.inits, B, 
         neighbors, D, fixed = list(smoothness = 0.5),
         blocks = list(type = "cluster", nblocks = 4),
         verbose = FALSE, tol = 0.001, maxIter = 100)
}
\arguments{
	\item{formula}{model formula}
	\item{data}{data used in model fit}
	\item{S}{n by 2 matrix of spatial locations}
	\item{cov}{spatial covariance type}
	\item{cov.inits}{initial values for covariance function}
	\item{B}{block memberships for observations}
	\item{neighbors}{two column matrix specifying unique neighbor relationships}
	\item{D}{distance matrix of observations}
	\item{fixed}{parameters fixed and not estimated during fit}
	\item{blocks}{setup for grid of blocks}
	\item{verbose}{show informational messages during fit}
	\item{tol}{convergence tolerance}
	\item{maxIter}{maximum Fisher scoring iterations}
}
\seealso{
	\code{\link{mean_max_temps}}
}
\references{
Jo Eidsvik, Benjamin A. Shaby, Brian J. Reich, Matthew Wheeler, and Jarad Niemi.
``Estimation and prediction in spatial models with block composite likelihoods using
parallel computing.'' Journal of Computational and Graphical Statistics (accepted)
}
\examples{

\dontrun{
# generate data
n <- 100
S <- matrix(runif(n*2), n, 2)
D <- rdist(S); diag(D) <- 0
Sigma <- 0.5 * diag(n) + 0.5*exp(-1.5*D)
y <- mvrnorm(1, rep(5,n), Sigma)

# fit model
fit <- spacious(y~1, S=S, cov="exp", blocks=list(type="regular", nblocks=2^2))
summary(fit)

# make predictions
predict(fit, newS=cbind(0.5,0.5), interval="prediction", level=0.95)

}

\dontrun{
# fit data with cluster block type
data(mean_max_temps)

fit <- spacious(temp~1, data=mean_max_temps,
  S=cbind(mean_max_temps$lat,mean_max_temps$long),
  blocks=list(type="cluster", nblocks=250)
)

summary(fit)
}


}
\keyword{spacious}
