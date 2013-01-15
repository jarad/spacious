\name{predict.spacious}
\alias{predict.spacious}
\title{Make predictions with spacious fit}
\description{
Predicts values at newS locations using a spacious fit.
}
\usage{
\method{predict}{spacious}(object, newdata, newS, newB, D,
        opts = list(type="block"),
        interval = "none", level = 0.95, ...)
}
\arguments{
	\item{object}{model fit with \code{\link{spacious}}.}
	\item{newdata}{data frame of variables used in the model fit. Not required with intercept only models.}
	\item{newS}{an n by 2 matrix of spatial locations to make predictions at.}
	\item{newB}{an optional vector of block memberships for prediction locations.}
	\item{D}{optional distance matrix for fitted and new locations specifying distances between observations.}
	\item{opts}{list of prediction options. Specified as a named list(type, num): \sQuote{type} specifies form of predictions, with \sQuote{block} (\bold{default}) kriging using the block composite likelihood, \sQuote{local} using the closest \sQuote{num} points for kriging, and \sQuote{all} kriging with the full likelihood.}
	\item{interval}{specifing if no interval (\sQuote{none}) or a prediction interval (\sQuote{prediction}) is computed.}
	\item{level}{confidence level of prediction interval.}
	\item{...}{further arguments passed to or from other methods.}
}
\seealso{
	\code{\link{spacious}}
}
