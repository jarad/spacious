# Computes initial values for matern covariance matrix
"initial.theta" <- function(y, S, p.nugget=0.1, p.range=0.1) {
	# Calculate <p.range>th percentile of distance between 1000 random points
	n      <- length(y)
	points <- sample(n, min(n,1000))
	d      <- dist(S[points,])
	range  <- quantile(d,p.range)

	vg           <- variog(list(coords=S[points,], data=y[points]), messages=FALSE)
	partial.sill <- 0.5*vg$v[1]/exp(-vg$u[1]/range)
	nugget       <- vg$var.mark-partial.sill

	# make sure initial nugget and partial still aren't too small
	nugget       <- max(nugget, 0.05)
	partial.sill <- max(partial.sill, 0.05)

	smoothness <- 0.5

	return(as.vector(c(nugget, partial.sill, range, smoothness)))
}
