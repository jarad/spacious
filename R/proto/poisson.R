# fit correlated Poisson data with copula

"corr_pois" <- function(y, V=0.5) {
	# assumes constant lambda and same correlation for all
	# Thus Sigma_ij = rho for all i != j and Sigma_ij = 1 otherwise

	invlogit <- function(x) { 1/(1+exp(-x)) }

	# to start let's just use nlm
	ull <- function(theta, y, n, V, nV) {
		lambda <- exp(theta[1])
		#rho <- 2*invlogit(theta[2]) - 1
		rho <- invlogit(theta[2])
		if (rho == 1) {
			return(1e4)
		}

		Sigma <- matrix(rho, nrow=n, ncol=n)
		diag(Sigma) <- 1

		cholSigma <- chol(Sigma)
		invSigma <- chol2inv(cholSigma)

		fy <- dpois(y, lambda=lambda, log=TRUE)
		efy <- exp(fy)

if (0) { # way to compute for a single V
		Fstar <- ppois(y-1, lambda=lambda) + (ystar-y) * exp(fy)
		z <- qnorm(Fstar)

		# -log lik
		ull <- sum(log(diag(cholSigma))) + 0.5 * (t(z) %*% ( chol2inv(cholSigma) - diag(n) ) %*% z) - sum(fy)
} else { # way to compute for multiple V

		# -log lik
		ull <- 0
		Fstar <- ppois(y-1, lambda=lambda)
		ull <- sum( unlist(
			mclapply(V, function(v) {
				z <- qnorm(Fstar + v * efy )
				0.5 * (t(z) %*% ( invSigma - diag(n) ) %*% z)
			}, mc.cores=4)
		) )
		ull <- ull/nV + sum(log(diag(cholSigma))) - sum(fy)
}

		if (is.na(ull)) { return(1e4) }

		ull
	}

	# initial values
	theta <- c(log(mean(y)), 0)

	#ystar <- y + V
	fit <- nlm(ull, theta, y=y, n=length(y), V=V, nV=length(V))
}

if (1) {

# simulate data
lambda <- 10
rho <- 0.75

n <- 250

Sigma <- matrix(rho, nrow=n, ncol=n)
diag(Sigma) <- 1

set.seed(311)
z <- mvrnorm(1, mu=rep(0, n), Sigma=Sigma)
y <- qpois(pnorm(z), lambda=lambda)

iy <- rpois(n, lambda=lambda)

#pdf("corr.pdf");hist(y,breaks=10);graphics.off()
#pdf("ind.pdf");hist(iy,breaks=10);graphics.off()


#fit <- corr_pois(y, V=0.5)
# TODO: change V!! it needs to be for each y; maybe call it U?
fit <- corr_pois(y, V=runif(10000))

}
