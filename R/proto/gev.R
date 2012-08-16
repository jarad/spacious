# fit correlated GEV data with copula

invlogit <- function(x) { 1/(1+exp(-x)) }

"corr_gev" <- function(y) {
	# assumes constant location mu in R, scale sigma > 0, shape alpha > 0 (Frechet), and same correlation for all
	# Thus Sigma_ij = rho for all i != j and Sigma_ij = 1 otherwise

	# to start let's just use nlm
	ull <- function(theta, y, n) {
		mu <- theta[1]
		sigma <- theta[2]
		alpha <- theta[3]
		rho <- theta[4]
#		sigma <- exp(theta[2])
#		alpha <- exp(theta[3])
#		rho <- invlogit(theta[4])

		if (rho >= 0.99) {
			return(1e4)
		}

		# must be within bounds
#		if ( sum( y < (mu - sigma/alpha) ) > 0 ) {
#			return(1e4)
#		}
print(round(c(mu,sigma,alpha,rho),3))

		Sigma <- matrix(rho, nrow=n, ncol=n)
		diag(Sigma) <- 1

		cholSigma <- chol(Sigma)
		invSigma <- chol2inv(cholSigma)

# TODO: compute pdf of GEV
		z <- qnorm( pgev(y, xi=alpha, mu=mu, beta=sigma) )
		gy <- 1 + (alpha/sigma)*(y-mu)
if ( sum(gy <= 0) > 1 ) {
	cat("gy <= 0!\n")
}

		# log like
		ll <- -sum(log(diag(cholSigma))) -0.5 * (t(z) %*% ( invSigma - diag(n) ) %*% z)
			-n * log(sigma) -(1+1/alpha) * sum(log(gy)) -sum( gy^(1/alpha) )

		if (is.na(ll)) { return(1e4) }

		# return negative log like
		-ll
	}

	# initial values
	theta <- c( min(y), 1, 1, 0.5)

	n <- length(y)

#	fit <- nlm(ull, theta, y=y, n=n)
	fit <- optim(par=theta, fn=ull, gr=NULL, y=y, n=n, method="L-BFGS-B",
		lower=c(0, 0.5, 0.5, 0), upper=c(Inf, Inf, Inf, 1)
	)
#optim(par, fn, gr = NULL, ...,
#method = c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN", "Brent"),
#lower = -Inf, upper = Inf, control = list(), hessian = FALSE)

}

if (1) {

require(MASS)
require(fExtremes)

# simulate data
mu <- 20
sigma <- 1
alpha <- 1
rho <- 0.50
#pdf("pdf/gev.pdf");plot(x<-seq(0,20,by=0.1),dgev(x,xi=alpha,mu=mu,beta=sigma),type="l");graphics.off()

n <- 500

Sigma <- matrix(rho, nrow=n, ncol=n)
diag(Sigma) <- 1

set.seed(311)
z <- mvrnorm(1, mu=rep(0, n), Sigma=Sigma)
y <- as.vector(qgev(pnorm(z), xi=alpha, mu=mu, beta=sigma))

fit <- corr_gev(y)
print(fit)

cat("Estimates:",fit$estimates,"\n")

}
