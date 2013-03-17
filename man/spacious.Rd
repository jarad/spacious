\encoding{UTF-8}
\name{spacious}
\alias{spacious}
\title{Fit spatial model with block composite likelihood}
\description{
Fits a spatial model using a block composite likelihood.
The block composite likelihood provides quick estimation of spatial
models at a cost of a small reduction in statistical efficiency.
}
\usage{
spacious(formula, data, S, cov = "exp", cov.inits, B, 
         neighbors, D, fixed = list(smoothness = 0.5),
         blocks = list(type = "cluster"),
         verbose = FALSE, tol = 0.001, maxIter = 100)
}
\arguments{
	\item{formula}{symbolic description of the mean structure of the model. See \code{\link{formula}}.}
	\item{data}{optional data frame containing the variables used in the model formula.}
	\item{S}{an n by 2 matrix of spatial locations for each observation}
	\item{cov}{spatial covariance type. Can be either one of \sQuote{exp} or \sQuote{matern}.}
	\item{cov.inits}{optional list of initial values for covariance function parameters. Specified as a named list(nugget, psill, range, smoothness): \sQuote{nugget} for the nugget, \sQuote{psill} for the partial sill, \sQuote{range} for the spatial range, and \sQuote{smoothness} for the \enc{Matérn}{Matern} smoothness.}
	\item{B}{an optional vector of block memberships for each observation.}
	\item{neighbors}{an optional two column matrix specifying the unique neighbor relationships of the blocks.}
	\item{D}{optional n by n distance matrix specifying distances between observations.}
	\item{fixed}{optional list of fixed covariance parameters that are not estimated. Specified as a named list(nugget, psill, range, smoothness): \sQuote{nugget} for the nugget, \sQuote{psill} for the partial sill, \sQuote{range} for the spatial range, and \sQuote{smoothness} for the \enc{Matérn}{Matern} smoothness (required for \sQuote{cov="matern"}).}
	\item{blocks}{optional list specifying the grid of blocks. Specified as a named list(type, nblocks): \sQuote{type} of blocks can be \sQuote{regular} for a square grid of blocks, \sQuote{cluster} (\bold{default}) for an irregular grid of blocks, \sQuote{pair} for the pairwise composite likelihood, or \sQuote{full} for the full likelihood; \sQuote{nblocks} specifying the number of blocks to place observations into, with the default being n/50, and \sQuote{type=regular} requiring \sQuote{nblocks} to be a perfect square. Also, \sQuote{nblocks <= 2} specifies the full likelihood.}
	\item{verbose}{boolean to show informational messages (\sQuote{verbose=TRUE}) during fit}
	\item{tol}{convergence tolerance of parameter estimates.}
	\item{maxIter}{maximum number of Fisher scoring iterations to compute.}
}
\seealso{
	\code{\link{predict.spacious}},
	\code{\link{plot.spacious}},
	\code{\link{mean_max_temps}}
}
\references{
Jo Eidsvik, Benjamin A. Shaby, Brian J. Reich, Matthew Wheeler, and Jarad Niemi.
\dQuote{Estimation and prediction in spatial models with block composite likelihoods using
parallel computing.} Journal of Computational and Graphical Statistics (accepted)
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
