# compute covariance matrix for specific types
"compute_cov" <- function(type="exp", theta, D) {
	if (type == "exp") {
		theta[1] * diag(nrow(D)) + theta[2] * exp(-theta[3] * D)
	} else {
		stop(paste("Unknown covariance type:",type))
	}
}
