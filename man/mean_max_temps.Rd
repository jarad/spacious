\name{mean_max_temps}
\docType{data}
\alias{mean_max_temps}
\title{Mean Maximum Tempeartures in Europe}
\description{
This data has a response Y that is the mean of the annual
maximum temperatues from 1990-2010. Includes spatial coordinates
and elevation. Has 942 observations containing a large cluster
of observations in Germany.
}
\usage{
data(mean_max_temps)
}
\format{
Data frame of 942 observations with temperature (temp),
country (country), elevation (elev), latitude (lat), and longitude (long).
}
\source{European Climate Assessment (http://eca.knmi.nl/)}

\references{
Shaby BA, Reich BJ. Bayesian spatial extreme value analysis to
assess the changing risk of concurrent high temperatures across large
portions of European cropland. Accepted, Environmetrics.
}

\seealso{
	\code{\link{spacious}}
}

\examples{

\dontrun{
data(mean_max_temps)

fit <- spacious(temp~1, data=mean_max_temps,
  S=cbind(mean_max_temps$lat,mean_max_temps$long),
  blocks=list(type="cluster", nblocks=250)
)

summary(fit)
}

}

\keyword{datasets}
