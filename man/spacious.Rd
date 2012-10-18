\name{spacious}
\alias{spacious}
\title{Fit block composite model}
\description{
	Fits a block composite model.
}
\usage{
	spacious(formula, data, S, cov="exp", cov.inits=NULL, nblocks=1, grid.type="regular",
		verbose=FALSE, tol=0.001, maxIter=100)
}
\arguments{
	\item{formula}{model formula}
	\item{data}{data used in model fit}
	\item{S}{n by 2 matrix of spatial locations}
	\item{cov}{spatial covariance type}
	\item{cov.inits}{initial values for covariance function}
	\item{nblocks}{number of blocks to divide data into}
	\item{grid.type}{type of blocking structure to place points into}
	\item{verbose}{show messages during fit}
	\item{tol}{convergence tolerance}
	\item{maxIter}{maximum Fisher scoring iterations}
}
\seealso{
}
\examples{
  n = 100
  S = matrix(runif(n*2), n, 2)
  D = rdist(S); diag(D) = 0
  Sigma = 0.5 * diag(n) + 0.5*exp(-1.5*D)
  y = mvrnorm(1, rep(mu,n), Sigma)
  fit = spacious(y~1, S=S)
  summary(fit)
}
\keyword{spacious}
