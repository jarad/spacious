# fit correlated Poisson data with copula

invlogit <- function(x) { 1/(1+exp(-x)) }

"corr_pois" <- function(y, nV=1) {
	# assumes constant lambda and same correlation for all
	# Thus Sigma_ij = rho for all i != j and Sigma_ij = 1 otherwise

	# to start let's just use nlm
	ull <- function(theta, y, ystar, n, V, nV) {
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

		Fstar <- ppois(y-1, lambda=lambda) + (ystar-y) * exp(fy)
		z <- qnorm(Fstar)

		# -log lik
		ull <- sum(log(diag(cholSigma))) + 0.5 * (t(z) %*% ( chol2inv(cholSigma) - diag(n) ) %*% z) - sum(fy)

		if (is.na(ull)) { return(1e4) }

		ull
	}

	# initial values
	theta <- c(log(mean(y)), 0)

	n <- length(y)

	if (nV == 1) {
		# use V = 0.5
		V <- matrix(0.5, nrow=1, ncol=n)
	} else {
		V <- matrix(runif(nV * n), nrow=nV, ncol=n)
	}
	fit <- nlm(ull, theta, y=y, ystar=y+V[1,], n=length(y), V=V, nV=nV)
}

if (1) {

require(MASS)

# simulate data
lambda <- 10
rho <- 0.75

n <- 500

Sigma <- matrix(rho, nrow=n, ncol=n)
diag(Sigma) <- 1

set.seed(311)
z <- mvrnorm(1, mu=rep(0, n), Sigma=Sigma)
y <- qpois(pnorm(z), lambda=lambda)

iy <- rpois(n, lambda=lambda)

#pdf("corr.pdf");hist(y,breaks=10);graphics.off()
#pdf("ind.pdf");hist(iy,breaks=10);graphics.off()

fit <- corr_pois(y, nV=1)
print(fit)

cat("Estimated lambda:",exp(fit$estimate[1]),"\n")
cat("Estimated rho:",invlogit(fit$estimate[2]),"\n")

}
