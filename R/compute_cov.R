# compute covariance matrix for specific types
"compute_cov" <- function(type="exp", theta, D) {
	# note that this function assumes there are no replicate sites
	# hence only the diagonal of D should be 0

	if (type == "exp") {
		theta[1] * diag(nrow(D)) + theta[2] * exp(-theta[3] * D)
	} else if (type == "matern") {
#		mid <- 2*D*theta[3]*sqrt(theta[4])
#		rho <- mid^theta[4] * besselK(mid, theta[4])/(2^(theta[4]-1) * gamma(theta[4]))
#		rho[is.na(rho)] <- 1
		#theta[1] * diag(nrow(D)) + theta[2] * mid^theta[4] * besselK(mid, theta[4])/(2^(theta[4]-1) * gamma(theta[4]))
		#theta[1] * diag(nrow(D)) + theta[2] * rho
		theta[1] * diag(nrow(D)) + theta[2] * matern_rho(theta, D)
	} else {
		stop(paste("Unknown covariance type:",type))
	}
}

"matern_rho" <- function(theta, D) {
	mid <- 2*D*theta[3]*sqrt(theta[4])
	rho <- mid^theta[4] * besselK(mid, theta[4])/(2^(theta[4]-1) * gamma(theta[4]))
	rho[is.na(rho)] <- 1

	rho
}
