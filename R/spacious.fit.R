# function to run algorithm for fitting a block composite model

"spacious.fit" <- function(y, X, D, nblocks, B, neighbors, cov, n, p, R, theta) {
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
		b[seq.p]  <<- 0

		apply(neighbors, 1, function(row) {
			in.pair <- B==row[1] | B==row[2]
			n.pair <- sum(in.pair)

			Sigma <- compute_cov(cov, exp(theta), D[in.pair,in.pair])
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
		u[seq.R]   <<- 0
		H[seq.RH]  <<- 0
		FI[seq.R2] <<- 0

		apply(neighbors, 1, function(row) {
			in.pair <- B==row[1] | B==row[2]
			n.pair <- sum(in.pair)

# TODO: figure out if it is posible to do this once since it's done in beta update as well (maybe merge stuff?)
			Sigma <- compute_cov(cov, exp(theta), D[in.pair,in.pair])
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

	# compute covariance matrix of parameters
	vcov.beta <- matrix(0, nrow=p, ncol=p)
	vcov.theta <- matrix(0, nrow=R, ncol=R)
	se.beta <- rep(0, p)
	se.theta <- rep(0, R)

	if (nblocks == 1) {
		# don't use the sandwich
		vcov.beta <- chol2inv(chol(A))
		vcov.theta <- chol2inv(chol(FI))
	} else {
		# use the sandwich

		J.beta  <- matrix(0, nrow=p, ncol=p)
		J.theta <- FI

		for (i in 1:nrow(neighbors)) {
			# which sites are in block i?
			in.i <- which(B==neighbors[i,1] | B==neighbors[i,2])
			n.i <- length(in.i)
			Sigma.i <- compute_cov(cov, exp(theta), D[in.i,in.i])
			invSigma.i <- chol2inv(chol(Sigma.i))

			for (j in 1:nrow(neighbors)) {
				# do pairs i and j have a common block?
				if (!any(neighbors[i,] == neighbors[j,1]) & !any(neighbors[i,] == neighbors[j,2])) {
					next;
				}

				# which sites are in block j?
				in.j <- which(B==neighbors[j,1] | B==neighbors[j,2])
				n.j <- length(in.j)
				Sigma.j <- compute_cov(cov, exp(theta), D[in.j,in.j])
				invSigma.j <- chol2inv(chol(Sigma.j))

				Sigma.ij <- compute_cov(cov, exp(theta), D[c(in.i,in.j),c(in.i,in.j)])

				J.beta <- J.beta + t(X[in.i,]) %*% invSigma.i %*% Sigma.ij[1:n.i,n.i+1:n.j] %*% invSigma.j %*% X[in.j,]

				if (j > i) {
					# update J.theta
# TODO: speed this up (maybe)
					sapply(seq.R, function(r) {
						sapply(r:R, function(s) {
							B.ir <- invSigma.i %*% partials[[r]](theta, n.i, in.i) %*% invSigma.i
							B.js <- invSigma.j %*% partials[[s]](theta, n.j, in.j) %*% invSigma.j
							add <- sum(diag( B.ir %*% Sigma.ij[1:n.i,n.i+1:n.j] %*% B.js %*% Sigma.ij[n.i+1:n.j,1:n.i] ))

							J.theta[r,s] <- J.theta[r,s] + add
							if (r != s) {
								J.theta[s,r] <- J.theta[s,r] + add
							}
						})
					})
				}
			}
		}

		vcov.beta <- chol2inv(chol(A %*% chol2inv(chol(J.beta)) %*% A))
		vcov.theta <- chol2inv(chol(FI %*% chol2inv(chol(J.theta)) %*% FI))
	}

	se.beta <- sqrt(diag(vcov.beta))
	for (i in 1:R) {
# TODO: is there a better way to do this?
		se.theta[i] <- exp(theta[i] + sqrt(vcov.theta[i,i])) - exp(theta[i])
	}

	# return estimates and standard errors
	list(beta=beta, theta=exp(theta),
		se.beta=se.beta, se.theta=se.theta
	)
}
