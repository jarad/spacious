# function to run algorithm for fitting a block composite model

"spacious.fit" <- function(y, X, D, nblocks, B, neighbors, cov, n, p, R, theta, theta.fixed,
	verbose, tol=1e-3, maxIter=100) {
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
	# theta.fixed: TRUE/FALSE indicating if parameter is fixed

	# verbose: print messages?
	# tol: error tolerance for identifying convergence
	# maxIter: maximum number of Fisher scoring iterations

	nH <- R*(R+1)/2

	seq.p <- 1:p
	seq.p2 <- 1:(p^2)
	seq.R <- 1:R
	seq.R2 <- 1:(R^2)
	seq.RH <- 1:nH

	# how many parameters are not fixed?
	Nnot_fixed <- length(theta.fixed)-sum(theta.fixed)
	which.fixed <- which(theta.fixed==TRUE)
	which.not_fixed <- which(theta.fixed==FALSE)

	# a list of functions to compute partial derivates
	# with respect to each covariance function param
	partials <- list()

	# transformation functions
	invlogit <- function(x) { 1/(1+exp(-x)) }
	#tsmooth <- function(x) { 0.5 + 3.5*invlogit(x) }
	tsmooth <- function(x) { exp(x) }
	t_theta <- function(x) { x }

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
		t_theta <- function(theta) { exp(theta) }
	} else if (cov == "matern") {
		partials <- list(
			function(theta, n.pair, in.pair) {
				exp(theta[1])*diag(n.pair)
			},
			function(theta, n.pair, in.pair) {
				theta4 <- tsmooth(theta[4])
				mid <- 2*D[in.pair,in.pair]*exp(theta[3])*sqrt(theta4)
				rho <- mid^theta4 * besselK(mid, theta4)/(2^(theta4-1) * gamma(theta4))
				rho[is.na(rho)] <- 1
				exp(theta[2])*rho
			},
			function(theta, n.pair, in.pair) {
				theta4 <- tsmooth(theta[4])
				mid <- 2*D[in.pair,in.pair]*exp(theta[3])*sqrt(theta4)
				rho <- -mid^(theta4+1)*besselK(mid, theta4-1)/(2^(theta4-1) * gamma(theta4))
				rho[is.na(rho)] <- 0
				exp(theta[2])*rho
			},
			function(theta, n.pair, in.pair) {
				e <- 1e-5
				p <- ( matern_rho( c(exp(theta[1:3]),tsmooth(theta[4]+e)), D[in.pair,in.pair] ) -
					matern_rho( c(exp(theta[1:3]),tsmooth(theta[4]-e)), D[in.pair,in.pair] ) )/(2*e)

				exp(theta[2]) * p
			}
		)
		t_theta <- function(theta) { c(exp(theta[1:3]),tsmooth(theta[4])) }
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

			Sigma <- compute_cov(cov, t_theta(theta), D[in.pair,in.pair])
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
		if (Nnot_fixed == 0) {   # don't do any updating
			return(theta)
		}

		u[seq.R]   <<- 0
		H[seq.RH]  <<- 0
		FI[seq.R2] <<- 0

		apply(neighbors, 1, function(row) {
			in.pair <- B==row[1] | B==row[2]
			n.pair <- sum(in.pair)

# TODO: figure out if it is posible to do this once since it's done in beta update as well (maybe merge stuff?)
			Sigma <- compute_cov(cov, t_theta(theta), D[in.pair,in.pair])
			invSigma <- chol2inv(chol(Sigma))

# TODO: see if any of this can be cleaned up
			Xb <- X[in.pair,] %*% beta
			ymXb <- y[in.pair]-Xb
			q <- invSigma %*% ymXb

			# compute the Ws
			for (r in 1:R) {
				if (theta.fixed[r]) { next; }

				partial <- partials[[r]](theta, n.pair, in.pair)
				W[[r]] <<- invSigma %*% partial
				u[r] <<- u[r] -0.5 * sum( diag(W[[r]]) ) + 0.5 * t(q) %*% partial %*% q
			}

			# compute the hessian H
			# this is stored in a compact form in H and then expanded into FI (fisher info)
			index <- 1
			sapply(seq.R, function(r) {
				sapply(r:R, function(s) {
					if (!theta.fixed[r] & !theta.fixed[s]) {
						H[index] <<- H[index] + 0.5 * sum(diag( W[[r]] %*% W[[s]] ))
					}
					index <<- index+1
				})
			})
		})

		index <- 1
		sapply(seq.R, function(r) {
			sapply(r:R, function(s) {
				if (!theta.fixed[r] & !theta.fixed[s]) {
					FI[r,s] <<- H[index]
					if (r != s) {
						FI[s,r] <<- H[index]
					}
				}
				index <<- index+1
			})
		})
#		invFI <- chol2inv(chol(FI))
#		invFI <- qr.solve(qr(FI))

#print(t_theta(theta))
#print( FI[which.not_fixed,which.not_fixed] )

		theta[which.not_fixed] <- theta[which.not_fixed] +
			chol2inv(chol(FI[which.not_fixed,which.not_fixed])) %*% u[which.not_fixed]

		theta
	}

	# compute log likelihood
	"loglik" <- function(beta, theta) {
		ll <- sum( apply(neighbors, 1, function(row) {
			in.pair <- B==row[1] | B==row[2]
			n.pair <- sum(in.pair)

			Sigma <- compute_cov(cov, t_theta(theta), D[in.pair,in.pair])
			cholSigma <- chol(Sigma)
			invSigma <- chol2inv(cholSigma)

			Xb <- X[in.pair,] %*% beta
			ymXb <- y[in.pair]-Xb

			-sum(log(diag(cholSigma))) -0.5 * t(ymXb) %*% invSigma %*% ymXb
		}) )
	}

# TODO: allow user to set control parameters
	# estimate params

	# get initial beta from initial theta
	beta <- update_beta(theta)

	ll <- -2 * loglik(beta,theta)

	# save theta and log lik at each iteration
	iters.theta <- t_theta(theta)
	iters.ll    <- ll

	for (iter in 1:maxIter) {
		prev.beta <- beta
		prev.theta <- theta

		# update theta
		theta <- update_theta(prev.beta, prev.theta)

		# update beta
		beta <- update_beta(theta)

		# get -2 * log likelihood
		ll <- -2 * loglik(beta, theta)

		# save values at each iteration
		iters.theta <- rbind( iters.theta, t_theta(theta) )
		iters.ll    <- c( iters.ll, ll )

		if (verbose) {
			cat("iter",iter,":"); print( c(beta, t_theta(theta), ll) )
		}

		if (iter > 1) {
# TODO: figure out the best way to identify convergence
			# have we converged?
			max_diff <- 0
			if (Nnot_fixed > 0) {
				max_diff <- max( c(
					abs(prev.beta - beta)/abs(beta),
					abs(prev.theta[which.not_fixed]-theta[which.not_fixed])/abs(theta[which.not_fixed])
				) )
			}

			if ( max_diff <= tol ) {
				if (verbose) {
					cat("Converged at iteration",iter,"\n")
				}

				break
			}
		}
	}

	convergence <- TRUE
	if (iter == maxIter) {
		warning("Possible issues with convergence: maximum number of iterations reached.")
		convergence <- FALSE
	}

	# compute covariance matrix of parameters
	vcov.beta <- matrix(0, nrow=p, ncol=p)
	vcov.theta <- matrix(0, nrow=R, ncol=R)
	se.beta <- rep(0, p)
	se.theta <- rep(0, R)

	# compute standard errors
	J.beta  <- matrix(0, nrow=p, ncol=p)
	J.theta <- FI

	for (i in 1:nrow(neighbors)) {
		# which sites are in block i?
		in.i       <- which(B==neighbors[i,1] | B==neighbors[i,2])
		n.i        <- length(in.i)
		Sigma.i    <- compute_cov(cov, t_theta(theta), D[in.i,in.i])
		invSigma.i <- chol2inv(chol(Sigma.i))

		for (j in 1:nrow(neighbors)) {
			# do pairs i and j have a common block?
			if (!any(neighbors[i,] == neighbors[j,1]) & !any(neighbors[i,] == neighbors[j,2])) {
				next;
			}

			if (i == j) {
				# this should only happen when we have one block
				J.beta <- J.beta + t(X[in.i,]) %*% invSigma.i %*% X[in.i,]
				next;
			}

			# which sites are in block j?
			in.j       <- which(B==neighbors[j,1] | B==neighbors[j,2])
			n.j        <- length(in.j)
			Sigma.j    <- compute_cov(cov, t_theta(theta), D[in.j,in.j])
			invSigma.j <- chol2inv(chol(Sigma.j))

			Sigma.ij <- compute_cov(cov, t_theta(theta), D[c(in.i,in.j),c(in.i,in.j)])

			J.beta <- J.beta + t(X[in.i,]) %*% invSigma.i %*% Sigma.ij[1:n.i,n.i+1:n.j] %*% invSigma.j %*% X[in.j,]

			if (j > i) {
				# update J.theta
# TODO: speed this up (maybe)
				sapply(seq.R, function(r) {
					sapply(r:R, function(s) {
						if (!theta.fixed[r] & !theta.fixed[s]) {
							B.ir <- invSigma.i %*% partials[[r]](theta, n.i, in.i) %*% invSigma.i
							B.js <- invSigma.j %*% partials[[s]](theta, n.j, in.j) %*% invSigma.j
							add <- sum(diag( B.ir %*% Sigma.ij[1:n.i,n.i+1:n.j] %*% B.js %*% Sigma.ij[n.i+1:n.j,1:n.i] ))

							J.theta[r,s] <- J.theta[r,s] + add
							if (r != s) {
								J.theta[s,r] <- J.theta[s,r] + add
							}
						}
					})
				})
			}
		}

		vcov.beta  <- chol2inv(chol(A %*% chol2inv(chol(J.beta)) %*% A))
		if (Nnot_fixed > 0) {
			vcov.theta[which.not_fixed,which.not_fixed] <- chol2inv(chol(FI[which.not_fixed,which.not_fixed] %*%
				chol2inv(chol(J.theta[which.not_fixed,which.not_fixed])) %*% FI[which.not_fixed,which.not_fixed]))
		}
	}

	# transform theta
	theta <- t_theta(theta)

	# get the standard errors
	se.beta  <- sqrt(diag(vcov.beta))
	se.theta <- as.vector(sqrt(diag(vcov.theta)) * theta)

	# transform theta[3] to range form
	theta[3]        <- 1/theta[3]
	se.theta[3]     <- se.theta[3] * theta[3]^2
	iters.theta[,3] <- 1/iters.theta[,3]

	# return estimates and standard errors
	list(
		convergence=convergence, nIter=iter, iters.theta=iters.theta, iters.ll=iters.ll,
		beta=beta, theta=theta, ll=ll,
		se.beta=se.beta, se.theta=se.theta,
		vcov.beta=vcov.beta, vcov.theta=vcov.theta
	)
}
