# estimate parameters for exponential model, where:
# - y ~ N(Xb, Sigma)
# - Sigma_ij = tau2 * I(i = j) + sigma2 * exp(-phi * d)
# - d = || s_i - s_j ||

"estimate_exp" <- function(y, X, D) {
	n <- nrow(X)
	p <- ncol(X)
	R <- 3

	# initial values
	beta <- rep(1, p)
	theta <- rep(0, R)  # tau2, sigma2, phi

	# partial derivative functions
	partials <- list(
		function(theta) {
			#diag(n)
			exp(theta[1])*diag(n)
		},
		function(theta) {
			#exp(-theta[3] * D)
			exp(theta[2])*exp(-exp(theta[3]) * D)
		},
		function(theta) {
			#-theta[2] * D * exp(-theta[3] * D)
			-exp(theta[2]+theta[3]) * D * exp(-exp(theta[3]) * D)
		}
	)

	# function to compute an iteration
	"iter" <- function(prev.beta, prev.theta) {
		#Sigma <- theta[1] * diag(n) + theta[2] * exp(-theta[3] * D)
		Sigma <- exp(prev.theta[1]) * diag(n) + exp(prev.theta[2]) * exp(-exp(prev.theta[3]) * D)
		cholSigma <- chol(Sigma)
		invSigma <- chol2inv(cholSigma)
#		cholInvSigma <- chol(invSigma)

		# obtain beta
		beta <- solve( t(X) %*% invSigma %*% X ) %*% t(X) %*% invSigma %*% y

		# obtain beta
		Xb <- X %*% beta
		ymXb <- y-Xb
		q <- invSigma %*% ymXb

		u <- rep(0, R)
		W <- vector("list", R)
		H <- rep(0, R*(R+1)/2)
		FI <- matrix(0, nrow=R, ncol=R)

		# compute the Ws
		for (r in 1:R) {
			partial <- partials[[r]](prev.theta)
			W[[r]] <- invSigma %*% partial
			u[r] <- -0.5 * sum( diag(W[[r]]) ) + 0.5 * t(q) %*% partial %*% q
		}

		index <- 1
		for (r in 1:R) {
			for (s in r:R) {
				H[index] <- 0.5 * sum(diag( W[[r]] %*% W[[s]] ))
				FI[r,s] <- H[index]
				if (r != s) {
					FI[s,r] <- H[index]
				}
				index <- index+1
			}
		}
#print(H)
#print(FI)
		invFI <- chol2inv(chol(FI))
#print(invFI)
#print(u)

		theta <- prev.theta + invFI %*% u

		list(beta=beta, theta=theta)
	}

	# estimate params
niter <- 20
	for (i in 1:niter) {
		prev.beta <- beta
		prev.theta <- theta

		step <- iter(beta, theta)
		beta <- step$beta
		theta <- step$theta
cat("iter",i,":");
print( c(beta, theta, exp(theta)) )
	}

	list(beta=beta, theta=exp(theta))
}

X <- matrix(1, nrow=length(y), ncol=1)
fit <- estimate_exp(y, X, D)

# print log likelihood
beta <- fit$beta
theta <- fit$theta
print(c(beta,theta))

"lik" <- function(x) {
	#Sigma <- theta[1] * diag(n) + theta[2] * exp(-theta[3] * D)
	Sigma <- theta[1] * diag(n) + theta[2] * exp(-x * D)
	invSigma <- chol2inv(chol(Sigma))
	-0.5 * log( det(Sigma) ) -0.5 * t(y-X%*%beta) %*% invSigma %*% (y-X%*%beta)
}

pdf("lik.pdf");
plot(x<-seq(0.2,2.2,by=0.01),sapply(x, function(x) { lik(x) }),type="l")
abline(v=theta[3])
abline(v=0.9031,lty=2)
graphics.off();

# try likfit
fit <- likfit(gd, ini=c(1,1))
#print(fit)
#print(summary(fit))

# print log likelihood
beta <- fit$beta
theta <- c(fit$tausq, fit$sigmasq, fit$phi)
print(c(beta,theta))
Sigma <- theta[1] * diag(n) + theta[2] * exp(-theta[3] * D)
invSigma <- chol2inv(chol(Sigma))
print( -0.5 * log( det(Sigma) ) -0.5 * t(y-X%*%beta) %*% invSigma %*% (y-X%*%beta))

