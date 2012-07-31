# function to run algorithm for fitting a block composite model

"spacious.fit" <- function(y, X, D, B, neighbors, cov, n, p, R, theta) {
	# y: response
	# X: model matrix
	# D: distance matrix
	# B: block memberships
	# neighbors: 2 column matrix listing unique neighbors of blocks
	# cov: covariance function type
	# n: number of observations
	# p: number of parameters for mean (columns of X)
	# R: number of covariance function parameters
	# theta: initial values for covariance parameters

	nH <- R*(R+1)/2

	seq.p <- 1:p
	seq.p2 <- 1:(p^2)
	seq.R <- 1:R
	seq.R2 <- 1:(R^2)
	seq.RH <- 1:nH

	# a list of functions to compute partial derivates
	# with respect to each covariance function param
	partials <- list()

	# setup fit based on covariance function
	if (cov == "exp") {
		partials <- list(
			function(theta, n.pair, in.pair) {
				exp(theta[1])*diag(n.pair)
			},
			function(theta, n.pair, in.pair) {
				exp(theta[2])*exp(-exp(theta[3]) * D[in.pair,in.pair])
			},
			function(theta, n.pair, in.pair) {
				-exp(theta[2]+theta[3]) * D[in.pair,in.pair] * exp(-exp(theta[3]) * D[in.pair,in.pair])
			}
		)
	} else {
		stop(paste("Unknown covariance type",cov))
	}

	# function to update beta based on theta
	A <- matrix(0, nrow=p, ncol=p)
	b <- matrix(0, nrow=p, ncol=1)
	"update_beta" <- function(theta) {
		# build A and b
		A[seq.p2] <<- 0
		b[seq.p] <<- 0
		apply(neighbors, 1, function(row) {
			in.pair <- B==row[1] | B==row[2]
			n.pair <- sum(in.pair)

# TODO: create a function to compute Sigma based on covariance type
			Sigma <- exp(theta[1])*diag(n.pair) + exp(theta[2])*exp(-exp(theta[3])*D[in.pair,in.pair])
			invSigma <- chol2inv(chol(Sigma))
			A <<- A + t(X[in.pair,]) %*% invSigma %*% X[in.pair,]
			b <<- b + t(X[in.pair,]) %*% invSigma %*% y[in.pair]
		})

		chol2inv(chol(A)) %*% b
	}

	# function to update theta with fisher scoring
	u <- rep(0, R)
	W <- vector("list", R)
	H <- rep(0, R*(R+1)/2)
	FI <- matrix(0, nrow=R, ncol=R)

	"update_theta" <- function(beta, theta) {
		u[seq.R] <- 0
		H[seq.RH] <- 0
		FI[seq.R2] <- 0

		apply(neighbors, 1, function(row) {
			in.pair <- B==row[1] | B==row[2]
			n.pair <- sum(in.pair)

# TODO: figure out how to do this once since it's done in beta update as well (maybe merge stuff?)
			Sigma <- exp(theta[1])*diag(n.pair) + exp(theta[2])*exp(-exp(theta[3])*D[in.pair,in.pair])
			invSigma <- chol2inv(chol(Sigma))

# TODO: see if any of this can be cleaned up
			Xb <- X[in.pair,] %*% beta
			ymXb <- y[in.pair]-Xb
			q <- invSigma %*% ymXb

			# compute the Ws
			for (r in 1:R) {
				partial <- partials[[r]](theta, n.pair, in.pair)
				W[[r]] <<- invSigma %*% partial
				u[r] <<- u[r] -0.5 * sum( diag(W[[r]]) ) + 0.5 * t(q) %*% partial %*% q
			}

			# compute the hessian H
			# this is stored in a compact form in H and then expanded into FI (fisher info)
			index <- 1
			sapply(seq.R, function(r) {
				sapply(r:R, function(s) {
					H[index] <<- H[index] + 0.5 * sum(diag( W[[r]] %*% W[[s]] ))
					index <<- index+1
				})
			})
		})

		index <- 1
		sapply(seq.R, function(r) {
			sapply(r:R, function(s) {
				FI[r,s] <<- H[index]
				if (r != s) {
					FI[s,r] <<- H[index]
				}
				index <<- index+1
			})
		})
#		invFI <- chol2inv(chol(FI))
#		invFI <- qr.solve(qr(FI))

		theta + chol2inv(chol(FI)) %*% u
	}

# TODO: allow user to set control parameters
	# estimate params

	# get initial beta from initial theta
	beta <- update_beta(theta)
	maxIter <- 100
	tol <- 1e-3
	for (i in 1:maxIter) {
		prev.beta <- beta
		prev.theta <- theta

		# update theta
		theta <- update_theta(prev.beta, prev.theta)

		# update beta
		beta <- update_beta(theta)
#cat("iter",i,":"); print( c(beta, exp(theta)) )

		if (i > 1) {
			# have we converged?
# TODO: figure out the best way to identify convergence
			if ( max( c(abs(prev.beta - beta),abs(prev.theta-theta)) ) <= tol ) {
cat("Converged at iteration",i,"\n")
				break
			}
		}
	}

	list(beta=beta, theta=exp(theta))
}
