# compute covariance matrix for specific types
"compute_cov" <- function(type="exp", theta, D) {
	# note that this function assumes there are no replicate sites
	# hence only the diagonal of D should be 0

	if (type == "exp") {
		theta[1] * diag(nrow(D)) + theta[2] * exp(-D/theta[3])
	} else if (type == "matern") {
		theta[1] * diag(nrow(D)) + theta[2] * matern_rho(theta, D)
	} else {
		stop(paste("Unknown covariance type:",type))
	}
}

# compute correlation for matern covariance
"matern_rho" <- function(theta, D) {
#print(c(theta[3], theta[4]))
#print(c( (D/theta[3])^(theta[4]), besselK(D/theta[3], theta[4]), ( 2^(theta[4]-1) * gamma(theta[4])) ))
	rho <- (D/theta[3])^(theta[4]) * besselK(D/theta[3], theta[4])/( 2^(theta[4]-1) * gamma(theta[4]))
	rho[is.na(rho)] <- 1

	rho
}

"matern_rho_u" <- function(theta, D) {
#print(c(exp(-theta[3]), exp(theta[4])))
#print(c( (D / exp(-theta[3]))^(exp(theta[4])), besselK(D/exp(-theta[3]), exp(theta[4])), ( 2^(exp(theta[4])-1) * gamma(exp(theta[4])) ) ))
	rho <- (D*exp(theta[3]))^(exp(theta[4])) * besselK(D*exp(theta[3]), exp(theta[4]))/( 2^(exp(theta[4])-1) * gamma(exp(theta[4])) )
	rho[is.na(rho)] <- 1

	rho
}

# compute correlation for matern covariance
"p_matern_rho_u" <- function(param, theta, D) {
	if (param == 1) {
		return(0);
	} else if (param == 1) {
		return(matern_rho_u(theta, D));
	} else if (param == 2) {
		rho <- (D/theta[3])^(theta[4]) * besselK(D/theta[3], theta[4])/( 2^(theta[4]-1) * gamma(theta[4]))
		rho[is.na(rho)] <- 1
	}

	rho
}

if (FALSE) {
# compute correlation for matern covariance
"matern_rho" <- function(theta, D) {
	#mid <- 2*D*theta[3]*sqrt(theta[4])
	mid <- 2*D*sqrt(theta[4])/theta[3]
	rho <- mid^theta[4] * besselK(mid, theta[4])/(2^(theta[4]-1) * gamma(theta[4]))
	rho[is.na(rho)] <- 1

	rho
}
}
