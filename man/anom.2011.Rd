\name{anom.2011}
\docType{data}
\alias{anom.2011}
\title{Tempearture Anomalies in Europe}
\description{
TODO
}
\usage{
data(anom.2011)
}
\format{
TODO
}
\source{TODO}

\references{
TODO
}

\seealso{
	\code{\link{spacious}}
}

\examples{

\dontrun{
data(anom.2011)

fit <- spacious(anom~lon+lat+elev, data=anom.2011,
  S=cbind(anom.2011$lon,anom.2011$lat)
)

summary(fit)

predict(fit, newdata=data.frame(elev=270), newS=cbind(50,10), interval="prediction")
}

}

\keyword{datasets}
