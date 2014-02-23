\encoding{UTF-8}
\name{spacious}
\alias{spacious}
\title{Fit geostatistical spatial model with block composite likelihood}
\description{
Fits a geostatistical spatial model using a block composite likelihood.
The block composite likelihood provides quick estimation of spatial
models at a cost of a small reduction in statistical efficiency.
}
\usage{
spacious(formula, data, S, cov = "exp", cov.inits, B, neighbors,
         fixed = list(smoothness = 0.5),
         blocks = list(type = "cluster"),
         verbose = FALSE, tol = 1e-8, maxIter = 100,
         nthreads = 1, gpu = FALSE, engine = "C")
}
\arguments{
	\item{formula}{symbolic description of the mean structure of the model. See \code{\link{formula}}.}
	\item{data}{optional data frame containing the variables used in the model formula.}
	\item{S}{an n by 2 matrix of spatial locations for each observation}
	\item{cov}{spatial covariance type. Can be either one of \sQuote{exp} or \sQuote{matern}.}
	\item{cov.inits}{optional list of initial values for covariance function parameters. Specified as a named list(nugget, psill, range, smoothness): \sQuote{nugget} for the nugget effect, \sQuote{psill} for the partial sill, \sQuote{range} for the spatial range, and \sQuote{smoothness} for the \enc{Matérn}{Matern} smoothness.}
	\item{B}{an optional vector of block memberships for each observation.}
	\item{neighbors}{an optional two column matrix specifying the unique neighbor relationships of the blocks.}
	\item{fixed}{optional list of fixed covariance parameters that are not to be estimated. Specified as a named list(nugget, psill, range, smoothness) with the same meaning as in \sQuote{cov.inits}.}
	\item{blocks}{optional list specifying how the grid of blocks is created. Specified as a named list(type, nblocks): \sQuote{type} of blocks can be \sQuote{regular} for a square grid of blocks, \sQuote{cluster} (\bold{default}) for an irregular grid of blocks, or \sQuote{full} for the full likelihood; \sQuote{nblocks} specifies the number of blocks to place observations into, with the default being n/50, and \sQuote{type=regular} requires that \sQuote{nblocks} is a perfect square. Also, \sQuote{nblocks <= 2} specifies the full likelihood.}
	\item{verbose}{boolean to show informational messages (\sQuote{verbose=TRUE}) during fit}
	\item{tol}{convergence tolerance of log-likelihood.}
	\item{maxIter}{maximum number of Fisher scoring iterations before stopping.}
	\item{nthreads}{when \sQuote{pthreads} are available, the number of threads to use when fitting the block composite likelihood. In most cases you will want to set this to be the number of processor cores available to you. See the \sQuote{spacious} manual for help enabling \sQuote{pthreads} support.}
	\item{gpu}{boolean to use an available \sQuote{nVidia} GPU for fitting full likelihood models when \sQuote{CUDA} support is enabled. See the \sQuote{spacious} manual for help enabling \sQuote{CUDA} support.}
	\item{engine}{can be one of \sQuote{C} or \sQuote{R}. Use of \sQuote{C} is recommended, as \sQuote{R} does not have threading or GPU support is primarily for testing and prototyping.}
}
\value{
	\item{convergence}{boolean indicating if the algorithm converged.}
	\item{nIter}{number of Fisher scoring iterations performed.}
	\item{time}{timing of fit. Includes Fisher scoring and standard error calculations. Does not include initial value for grid creation timings.}

	\item{ll}{log-likelihood for estimated parameters.}
	\item{iters_ll}{log-likelihood at each iteration.}
	\item{beta}{mean parameter estimates.}
	\item{se_beta}{estimated standard errors of mean parameters.}
	\item{vcov_beta}{estimated covariance matrix of mean parameters.}
	\item{theta}{covariance paramter estimates.}
	\item{se_theta}{estimated standard errors of covariance parameters.}
	\item{vcov_theta}{estimated covariance matrix of covariance parameters.}
	\item{theta_fixed}{vector indicating which covariance parameters were held fixed.}
	\item{iters_theta}{matrix of covariance parameter values at each iteration.}

	\item{y}{response.}
	\item{X}{model matrix.}
	\item{S}{observation locations.}
	\item{grid}{grid of polygons that specify the region contained in each block.}
	\item{B}{block memberships.}
	\item{neighbors}{matrix of neighbor pairs.}
	\item{nblocks}{number of blocks.}
	\item{lik_form}{likelihood type.}
	\item{cov}{covariance type.}
	\item{fitted}{fitted values of data used in fit.}
	\item{resids}{residuals of fitted values.}
	\item{verbose}{was verbose used in fit?}
	\item{tol}{convergence tolerance.}
	\item{max_iter}{maximum number of iterations or fit.}
	\item{nthreads}{number of threads used in fit.}
	\item{gpu}{boolean indicating use of GPU.}
	\item{engine}{engine used to perform estimation.}
	\item{terms}{terms used in creating model information.}
}
\seealso{
	\code{\link{predict.spacious}},
	\code{\link{plot.spacious}},
	\code{\link{mean_max_temps}}
}
\references{
Jo Eidsvik, Benjamin A. Shaby, Brian J. Reich, Matthew Wheeler, and Jarad Niemi.
\dQuote{Estimation and prediction in spatial models with block composite likelihoods.}
Journal of Computational and Graphical Statistics (accepted)
}
\examples{

\dontrun{
# generate data
set.seed(1983)
n     <- 1000
S     <- matrix(runif(n*2), n, 2)
D     <- rdist(S); diag(D) <- 0
Sigma <- 0.5 * diag(n) + 0.5*exp(-D/0.10)
y     <- 5 + t(chol(Sigma)) \%*\% rnorm(n)

# fit full model
fit.full <- spacious(y~1, S=S, cov="exp",
	blocks=list(type="full"), verbose=TRUE)
summary(fit.full)
plot(fit.full)

# fit BCL model
fit.bcl <- spacious(y~1, S=S, cov="exp",
	blocks=list(type="regular", nblocks=4^2), verbose=TRUE, nthreads=4)
summary(fit.bcl)
plot(fit.bcl)

# make predictions
pred.full <- predict(fit.full, newS=cbind(0.5,0.5),
	interval="prediction", level=0.95)
print(pred.full)

pred.bcl <- predict(fit.bcl, newS=cbind(0.5,0.5),
	interval="prediction", level=0.95)
print(pred.bcl)
}

\dontrun{
# fit data with cluster block type
data(anom.2011)

fit.full <- spacious(anom~lon+lat+elev, data=anom.2011,
  S=cbind(anom.2011$lon, anom.2011$lat),
  blocks=list(type="full"), verbose=TRUE)
summary(fit.full)
plot(fit.full)

fit.bcl <- spacious(anom~lon+lat+elev, data=anom.2011,
  S=cbind(anom.2011$lon, anom.2011$lat),
  blocks=list(type="cluster", nblocks=4^2), verbose=TRUE, nthreads=4)
summary(fit.bcl)
plot(fit.bcl)
}


}
\keyword{spacious}
