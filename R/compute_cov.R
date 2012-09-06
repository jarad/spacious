# compute covariance matrix for specific types
"compute_cov" <- function(type="exp", theta, D) {
	# note that this function assumes there are no replicate sites
	# hence only the diagonal of D should be 0

	if (type == "exp") {
		theta[1] * diag(nrow(D)) + theta[2] * exp(-theta[3] * D)
	} else {
		stop(paste("Unknown covariance type:",type))
	}
}
