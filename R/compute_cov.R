# compute covariance matrix for specific types
"compute_cov" <- function(type="exp", theta, D) {
	# note that this function assumes there are no replicate sites
	# hence only the diagonal of D should be 0

	if (type == "exp") {
		theta[1] * diag(nrow(D)) + theta[2] * exp(-theta[3] * D)
	} else if (type == "matern") {
		theta[1] * diag(nrow(D)) + theta[2] * matern_rho(theta, D)
	} else {
		stop(paste("Unknown covariance type:",type))
	}
}

# compute correlation for matern covariance
"matern_rho" <- function(theta, D) {
	mid <- 2*D*theta[3]*sqrt(theta[4])
	rho <- mid^theta[4] * besselK(mid, theta[4])/(2^(theta[4]-1) * gamma(theta[4]))
	rho[is.na(rho)] <- 1

	rho
}
