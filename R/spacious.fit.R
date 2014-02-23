# function to run algorithm for fitting a block composite model
"spacious.fit" <- function(y, X, S, nblocks, B, neighbors, cov, n, p, R, theta, theta_fixed,
	verbose, tol=1e-8, maxIter=100, compute_se=FALSE) {
	# y: response
	# X: model matrix
	# S: spatial locations
	# nblocks: number of blocks
	# B: block memberships
	# neighbors: 2 column matrix listing unique neighbors of blocks
	# cov: covariance function type
	# n: number of observations
	# p: number of parameters for mean (columns of X)
	# R: number of covariance function parameters
	# theta: initial values for covariance parameters
	# theta_fixed: TRUE/FALSE indicating if parameter is fixed

	# verbose: print messages?
	# tol: error tolerance for identifying convergence
	# maxIter: maximum number of Fisher scoring iterations

	# create distance matrix from spatial locations
	D <- rdist(S)
	diag(D) <- 0

	nH <- R*(R+1)/2

	seq.p <- 1:p
	seq.p2 <- 1:(p^2)
	seq.R <- 1:R
	seq.R2 <- 1:(R^2)
	seq.RH <- 1:nH

	# how many parameters are not fixed?
	Nnot_fixed <- length(theta_fixed)-sum(theta_fixed)
	which.fixed <- which(theta_fixed==TRUE)
	which.not_fixed <- which(theta_fixed==FALSE)

	# a list of functions to compute partial derivates
	# with respect to each covariance function param
	partials <- list()

	# transformation functions
	invlogit <- function(x) { 1/(1+exp(-x)) }
	#tsmooth <- function(x) { 0.5 + 3.5*invlogit(x) }
	tsmooth <- function(x) { exp(x) }
	t_theta <- function(theta) { theta }

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
		#t_theta <- function(theta) { exp(theta) }
		t_theta <- function(theta) { c(exp(theta[1]), exp(theta[2]), exp(-theta[3])) }
		u_theta <- function(theta) { c(log(theta[1]), log(theta[2]), -log(theta[3])) }
	} else if (cov == "matern") {
#rho <- (D/theta[3])^(theta[4]) * besselK(D/theta[3], theta[4])/( 2^(theta[4]-1) * gamma(theta[4]))
#rho <- (D exp(-theta[3]))^(exp(theta[4])) * besselK(D exp(-theta[3]), exp(theta[4]))/( 2^(exp(theta[4])-1) * gamma(exp(theta[4])))
		partials <- list(
			function(theta, n.pair, in.pair) {
				exp(theta[1])*diag(n.pair)
			},
			function(theta, n.pair, in.pair) {
#				theta4 <- tsmooth(theta[4])
#				mid <- 2*D[in.pair,in.pair]*exp(theta[3])*sqrt(theta4)
#				rho <- mid^theta4 * besselK(mid, theta4)/(2^(theta4-1) * gamma(theta4))
#				rho[is.na(rho)] <- 1
#				exp(theta[2])*rho

				exp(theta[2]) * matern_rho(t_theta(theta), D[in.pair,in.pair])
			},
			function(theta, n.pair, in.pair) {
#				theta4 <- tsmooth(theta[4])
#				mid <- 2*D[in.pair,in.pair]*exp(theta[3])*sqrt(theta4)
#				rho <- -mid^(theta4+1)*besselK(mid, theta4-1)/(2^(theta4-1) * gamma(theta4))
#				rho[is.na(rho)] <- 0
#				exp(theta[2])*rho
				mid <- D[in.pair,in.pair] * exp(theta[3])
				rho <- -(mid)^(exp(theta[4])+1) * besselK(mid, exp(theta[4])-1)/( 2^(exp(theta[4])-1) * gamma(exp(theta[4])))
				rho[is.na(rho)] <- 0
				exp(theta[2]) * rho
			},
			function(theta, n.pair, in.pair) {
				e <- 1e-5
				p <- ( matern_rho( c(exp(theta[1:3]),tsmooth(theta[4]+e)), D[in.pair,in.pair] ) -
					matern_rho( c(exp(theta[1:3]),tsmooth(theta[4]-e)), D[in.pair,in.pair] ) )/(2*e)

				exp(theta[2]) * p
			}
		)
		t_theta <- function(theta) { c(exp(theta[1]), exp(theta[2]), exp(-theta[3]), exp(theta[4])) }
		u_theta <- function(theta) { c(log(theta[1]), log(theta[2]), -log(theta[3]), log(theta[4])) }
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
			if (row[1] == row[2]) in.pair <- which(B==row[1]) # full likelihood
			else                  in.pair <- c(which(B==row[1]),  which(B==row[2]))

			Sigma    <- compute_cov(cov, t_theta(theta), D[in.pair,in.pair])
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

			Sigma <- compute_cov(cov, t_theta(theta), D[in.pair,in.pair])
			invSigma <- chol2inv(chol(Sigma))

			Xb <- X[in.pair,] %*% beta
			ymXb <- y[in.pair]-Xb
			q <- invSigma %*% ymXb

			# compute the Ws
			for (r in 1:R) {
				if (theta_fixed[r]) { next; }

				partial <- partials[[r]](theta, n.pair, in.pair)
				W[[r]] <<- invSigma %*% partial
				u[r] <<- u[r] -0.5 * sum( diag(W[[r]]) ) + 0.5 * t(q) %*% partial %*% q
			}

			# compute the hessian H
			# this is stored in a compact form in H and then expanded into FI (fisher info)
			index <- 1
			sapply(seq.R, function(r) {
				sapply(r:R, function(s) {
					if (!theta_fixed[r] & !theta_fixed[s]) {
						H[index] <<- H[index] + 0.5 * sum(diag( W[[r]] %*% W[[s]] ))
					}
					index <<- index+1
				})
			})

		})

		index <- 1
		sapply(seq.R, function(r) {
			sapply(r:R, function(s) {
				if (!theta_fixed[r] & !theta_fixed[s]) {
					FI[r,s] <<- H[index]
					if (r != s) {
						FI[s,r] <<- H[index]
					}
				}
				index <<- index+1
			})
		})

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
			#invSigma <- chol2inv(cholSigma)

			ymXb <- y[in.pair] - (X[in.pair,] %*% beta)
			L <- forwardsolve(t(cholSigma), ymXb)

			#-sum(log(diag(cholSigma))) -0.5 * t(ymXb) %*% invSigma %*% ymXb
			#-0.5*log(det(Sigma)) -0.5 * t(ymXb) %*% invSigma %*% ymXb
			-sum(log(diag(cholSigma))) -0.5 * t(L) %*% L
		}) )
	}

	# estimate params

	# incoming theta need to be transformed
	theta <- u_theta(theta)

	# get initial beta from initial theta
	beta   <- update_beta(theta)
	Nbeta  <- length(beta)
	Ntheta <- length(theta)

	names.show <- c("beta")
	if (Nbeta > 1) {
		names.show <- c(names.show, rep("", Nbeta-1))
	}
	names.show <- c(names.show, "theta", rep("", Ntheta-1))
	names.show <- c(names.show, "log lik")

	# initial log-likelihood
	ll <- loglik(beta,theta)

	# save theta and log lik at each iteration
	iters_theta <- t_theta(theta)
	iters_ll    <- ll

	for (iter in 1:maxIter) {
		prev.beta <- beta
		prev.theta <- theta

		# update theta
		theta <- update_theta(prev.beta, prev.theta)

		# update beta
		beta <- update_beta(theta)


		# get log likelihood
		ll <- loglik(beta, theta)

		# save values at each iteration
		iters_theta <- rbind( iters_theta, t_theta(theta) )
		iters_ll    <- c( iters_ll, ll )

		if (verbose) {
			cat(
				paste0("iter=",iter,
					": beta: ",paste(round(beta,2),collapse=" "),
					" ; theta: ",paste(round(t_theta(theta),2),collapse=" "),
					" ; ll=", round(ll,2)),
			"\n")
		}

		# have we converged?
		if (abs(ll - iters_ll[iter])/abs(0.1 + ll) <= tol) {
			if (verbose) {
				cat("Converged at iteration",iter,"\n")
			}

			break
		}
	}

	convergence <- TRUE
	if (iter == maxIter) {
		warning("Possible issues with convergence: maximum number of iterations reached.")
		convergence <- FALSE
	}

	# compute covariance matrix of parameters
	vcov_beta <- matrix(0, nrow=p, ncol=p)
	vcov_theta <- matrix(0, nrow=R, ncol=R)
	se_beta <- rep(0, p)
	se_theta <- rep(0, R)

if (TRUE) { # compute standard errors
	# update hessian
	update_theta(beta, theta)

	# compute standard errors
	if (nrow(neighbors) <= 1) { # full likelihood
		vcov_beta  <- chol2inv(chol(A))
		if (Nnot_fixed > 0) {
			vcov_theta[which.not_fixed,which.not_fixed] <- chol2inv(chol(FI[which.not_fixed,which.not_fixed]))
		}
	} else { # BCL
		J.beta  <- A #matrix(0, nrow=p, ncol=p)
		J.theta <- FI

		for (i in 1:nrow(neighbors)) {
			for (j in i:nrow(neighbors)) {
				if (j <= i) next;

				# do pairs i and j have a common block?
				if (!any(neighbors[i,] == neighbors[j,1]) & !any(neighbors[i,] == neighbors[j,2])) {
					# none in common; skip this pair
					next;
				}

				k1         <- neighbors[i,1]
				l1         <- neighbors[i,2]
				in.k1      <- which(B==k1)
				in.l1      <- which(B==l1)
				n.k1       <- length(in.k1)
				n.l1       <- length(in.l1)
				k2         <- neighbors[j,1]
				l2         <- neighbors[j,2]
				in.k2      <- which(B==k2)
				in.l2      <- which(B==l2)
				n.k2       <- length(in.k2)
				n.l2       <- length(in.l2)

				# which sites are in pair i?
				in.i       <- c(in.k1, in.l1)
				n.i        <- n.k1+n.l1
				Sigma.i    <- compute_cov(cov, t_theta(theta), D[in.i,in.i])
				invSigma.i <- chol2inv(chol(Sigma.i))

if (FALSE) {
				# zero out covariance when these pairs are not neighbors
				if (!any(neighbors[c(which(neighbors[,1]==k1),which(neighbors[,2]==k1)),] == k2)) {
					# pair (k1,k2) is zero
					Sigma.ij[1:n.k1,n.i+1:n.k2] <- 0
					Sigma.ij[n.i+1:n.k2,1:n.k1] <- 0
				}

				if (!any(neighbors[c(which(neighbors[,1]==l1),which(neighbors[,2]==l1)),] == k2)) {
					# pair (l1,k2) is zero
					Sigma.ij[n.k1+1:n.l1,n.i+1:n.k2] <- 0
					Sigma.ij[n.i+1:n.k2,n.k1+1:n.l1] <- 0
				}

				if (!any(neighbors[c(which(neighbors[,1]==k1),which(neighbors[,2]==k1)),] == l2)) {
					# pair (k1,l2) is zero
					Sigma.ij[1:n.k1,n.i+n.k2+1:n.l2] <- 0
					Sigma.ij[n.i+n.k2+1:n.l2,1:n.k1] <- 0
				}

				if (!any(neighbors[c(which(neighbors[,1]==l1),which(neighbors[,2]==l1)),] == l2)) {
					# pair (l1,l2) is zero
					Sigma.ij[n.k1+1:n.l1,n.i+n.k2+1:n.l2] <- 0
					Sigma.ij[n.i+n.k2+1:n.l2,n.k1+1:n.l1] <- 0
				}
}

				if (i == j) {
					# do not add contributions between neighbors in this case

if (FALSE) { # this already comes from A
					# add to beta
					J.beta <- J.beta + t(X[in.i,]) %*% invSigma.i %*% X[in.i,]
}

if (FALSE) { # this already comes from FI
					# add to theta
					sapply(seq.R, function(r) {
						sapply(r:R, function(s) {
							if (!theta_fixed[r] & !theta_fixed[s]) {
								add <- .5*sum(diag( invSigma.i %*% partials[[r]](theta, n.i, in.i) %*% invSigma.i %*% partials[[s]](theta, n.i, in.i) ))

								J.theta[r,s] <<- J.theta[r,s] + add
								if (r != s) {
									J.theta[s,r] <<- J.theta[s,r] + add
								}
							}
						})
					})
}
				} else {
					# which sites are in pair j?
					in.j       <- c(in.k2, in.l2)
					n.j        <- length(in.j)
					Sigma.j    <- compute_cov(cov, t_theta(theta), D[in.j,in.j])
					invSigma.j <- chol2inv(chol(Sigma.j))

					# compute covariance between pairs of blocks
					Sigma.ij   <- compute_cov(cov, t_theta(theta), D[c(in.i,in.j),c(in.i,in.j)])

					# add to beta
					#J.beta <- J.beta + t(X[in.i,]) %*% invSigma.i %*% Sigma.ij[1:n.i,n.i+1:n.j] %*% invSigma.j %*% X[in.j,]
					tmp <- t(X[in.i,]) %*% invSigma.i %*% Sigma.ij[1:n.i,n.i+1:n.j] %*% invSigma.j %*% X[in.j,]
					J.beta <- J.beta + tmp + t(tmp)

					# add to theta
					if (j > i) {
						sapply(seq.R, function(r) {
							sapply(r:R, function(s) {
								if (!theta_fixed[r] & !theta_fixed[s]) {
									B.ir <- invSigma.i %*% partials[[r]](theta, n.i, in.i) %*% invSigma.i
									B.js <- invSigma.j %*% partials[[s]](theta, n.j, in.j) %*% invSigma.j
									add  <- sum(diag( B.ir %*% Sigma.ij[1:n.i,n.i+1:n.j] %*% B.js %*% Sigma.ij[n.i+1:n.j,1:n.i] ))

									J.theta[r,s] <<- J.theta[r,s] + add
									if (r != s) {
										J.theta[s,r] <<- J.theta[s,r] + add
									}
								}
							})
						})
					}

				}

			}

		}


if (FALSE) {
		print(round(FI,3))
		print(round(J.theta,3))
		print(round(FI %*% solve(J.theta),3))
}

		vcov_beta  <- chol2inv(chol(A %*% chol2inv(chol(J.beta)) %*% A))
		if (Nnot_fixed > 0) {
			vcov_theta[which.not_fixed,which.not_fixed] <- chol2inv(chol(FI[which.not_fixed,which.not_fixed] %*%
				chol2inv(chol(J.theta[which.not_fixed,which.not_fixed])) %*% FI[which.not_fixed,which.not_fixed]))
		}

	} # end BCL
} # end std errors

	# transform theta
	theta <- t_theta(theta)

	# get the standard errors
	se_beta  <- sqrt(diag(vcov_beta))
	se_theta <- as.vector(sqrt(diag(vcov_theta)) * theta)

	# transform theta[3] to range form
	#theta[3]        <- 1/theta[3]
	#se_theta[3]     <- se_theta[3] * theta[3]^2
	#iters_theta[,3] <- 1/iters_theta[,3]

	# return estimates and standard errors
	list(
		convergence=convergence, nIter=iter,
		iters_theta=iters_theta, iters_ll=iters_ll,
		beta=beta, theta=theta, ll=ll,
		se_beta=se_beta, se_theta=se_theta,
		vcov_beta=vcov_beta, vcov_theta=vcov_theta
	)
}
