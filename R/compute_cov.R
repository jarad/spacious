# compute covariance matrix for specific types
"compute_cov" <- function(type="exp", theta, D, nu=0.5) {
	# note that this function assumes there are no replicate sites
	# hence only the diagonal of D should be 0

	if (type == "exp") {
		theta[1] * diag(nrow(D)) + theta[2] * exp(-theta[3] * D)
	} else if (type == "matern") {
		mid <- 2*D*theta[3]*sqrt(nu)
		rho <- mid^nu * besselK(mid, nu)/(2^(nu-1) * gamma(nu))
		rho[is.na(rho)] <- 1
		#theta[1] * diag(nrow(D)) + theta[2] * mid^nu * besselK(mid, nu)/(2^(nu-1) * gamma(nu))
		theta[1] * diag(nrow(D)) + theta[2] * rho
	} else {
		stop(paste("Unknown covariance type:",type))
	}
}
