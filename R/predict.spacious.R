# predict values at new sites using fit from spacious block composite model
"predict.spacious" <- function(object, newdata=NULL, newS=NULL,
	interval = "none", level = 0.95, ...) {

# TODO: error check inputs

	# ensure this is a matrix
	newS <- matrix(newS, ncol=2)

	# number of blocks
	nB <- length(object$grid)
	nFit <- nrow(object$S)
	nNew <- nrow(newS)

	# build model matrix for the new data
	terms <- delete.response(terms(object))
	mf <- model.frame(terms, newdata)
	X <- model.matrix(terms, mf)

	# locations for all points
	S <- rbind(object$S, newS)

	# compute distance matrix
	D <- rdist(S)
	D[row(D)==col(D)] <- 0

	# vector to hold predictions
	y_0      <- rep(NA, nNew)

	if (interval == "prediction") {
		# vector to hold prediction std devs
		sd.pred <- rep(NA, nNew)
	}

	# figure out newB the location is in
	if (nB == 0) {  # no blocks
		# perform traditional kriging
		n <- nrow(S)

		# compute covariance matrix
		Sigma <- compute_cov(object$cov, object$theta, D)

		# get the predictions
		y_0[1:nNew] <- X %*% object$beta + Sigma[nFit+1:nNew,1:nFit] %*%
			chol2inv(chol(Sigma[1:nFit,1:nFit])) %*% object$resids

		if (interval == "prediction") {
			invSigma <- chol2inv(chol(Sigma))

			sd.pred <- sqrt(diag( chol2inv(chol( invSigma[nFit+1:nNew,nFit+1:nNew] )) ))
		}
	} else {
		# predict when we have blocks

		# figure out newBs contain these new locations
		newB <- rep(NA, nNew)

		if (is.null(grid)) {
			# we don't have polygons, so find the closest point and use that block
			warning("No polygons defining blocks. Using the closest point to find a block.")
# TODO: do this
		} else {
			for (b in 1:nB) {
				in_poly <- point.in.polygon(newS[,1], newS[,2],
					object$grid@polygons[[b]]@Polygons[[1]]@coords[,1],
					object$grid@polygons[[b]]@Polygons[[1]]@coords[,2]) == 1

				newB[in_poly] <- b
			}

			if (sum(is.na(newB))) {
				warning("Unable to find blocks for all prediction sites. Using closest point(s) to find block(s).")
# TODO: do this
			}

			# block memberships for all points
			B <- c(object$B, newB)

			# predict for each unique block
			y_0[1:nNew] <- X %*% object$beta

			# unique blocks where we have points to predict
			uB <- unique(newB)

# TODO: we could re-construct this so that we only invert each matrix once. worth considering if speed more important than memory.
			for (b in uB) {
				# neighbors of this block
				neighbors <- as.vector(object$neighbors[which( rowSums( object$neighbors==b ) == 1 ),])
				neighbors <- neighbors[neighbors != b]

				# what new points are in this block?
				in.new   <- which(newB == b)+nFit
				n.in.new <- length(in.new)

				# what observed points are in this block?
				in.obs   <- which(object$B == b)
				n.in.obs <- length(in.obs)

				# new and observed in this block
				in.b   <- c(in.new, in.obs)
				n.in.b <- n.in.new + n.in.obs

				A_0 <- matrix(0, nrow=n.in.new, ncol=n.in.new)
				b_0 <- rep(0, n.in.new)

				if (interval == "prediction") {
					J_0 <- matrix(0, nrow=n.in.new, ncol=n.in.new)
					B_0 <- matrix(0, nrow=n.in.new, ncol=n.in.new+n.in.obs)
				}

# TODO: identify overlap of n1 and n2 so that some of these matrix computations might be reduced
				for (n1 in neighbors) {
					in.n1   <- which(object$B == n1)
					n.in.n1 <- length(in.n1)

					# take points so that we have Y = (Y_{b,new}, Y_{b,obs}, Y_{n1})
					in.pair.n1 <- c(in.b, in.n1)

					# covariance for b and n1
					Sigma.n1    <- compute_cov(object$cov, object$theta, D[in.pair.n1,in.pair.n1])
					invSigma.n1 <- chol2inv(chol(Sigma.n1))

					A_0 <- A_0 + invSigma.n1[1:n.in.new,1:n.in.new]
					b_0 <- b_0 + invSigma.n1[1:n.in.new,n.in.new+1:n.in.obs] %*% object$resids[in.obs] +
						invSigma.n1[1:n.in.new,n.in.new+n.in.obs+1:n.in.n1] %*% object$resids[in.n1]

					if (interval == "prediction") {
						B_0 <- B_0 + cbind(
							matrix(invSigma.n1[1:n.in.new,1:n.in.new],nrow=n.in.new,ncol=n.in.new),
							matrix(invSigma.n1[1:n.in.new,n.in.new+1:n.in.obs],nrow=n.in.new,ncol=n.in.obs)
						)

						for (n2 in neighbors) {
							in.n2   <- which(object$B == n2)
							n.in.n2 <- length(in.n2)

							# take points so that we have Y = (Y_{b,new}, Y_{b,obs}, Y_{n2})
							in.pair.n2 <- c(in.b, in.n2)

							# points for Y = (Y_{n1}, Y_{n2})
							in.n1n2 <- c(in.n1, in.n2)

							# covariance for b and n2
							Sigma.n2    <- compute_cov(object$cov, object$theta, D[in.pair.n2,in.pair.n2])
							invSigma.n2 <- chol2inv(chol(Sigma.n2))

							# covariance for n1 and n2
							Sigma.n1n2 <- compute_cov(object$cov, object$theta, D[in.n1n2,in.n1n2])

							J_0 <- J_0 + invSigma.n1[1:n.in.new,n.in.new+n.in.obs+1:n.in.n1] %*%
								Sigma.n1n2[1:n.in.n1,n.in.n1+1:n.in.n2] %*% invSigma.n2[n.in.new+n.in.obs+1:n.in.n2,1:n.in.new]
						}
					}
				}

				# complete predictions
				b_0 <- -b_0
				y_0[in.new-nFit] <- y_0[in.new-nFit] + chol2inv(chol(A_0)) %*% b_0

				if (interval == "prediction") {
					# complete prediction variances
					H_0 <- -A_0

					Sigma <- compute_cov(object$cov, object$theta, D[in.b,in.b])
					J_0   <- J_0 + B_0 %*% Sigma %*% t(B_0)

					for (n1 in neighbors) {
						in.n1   <- which(object$B == n1)
						n.in.n1 <- length(in.n1)

						# take points so that we have Y = (Y_{b,new}, Y_{b,obs}, Y_{n1})
						in.pair.n1 <- c(in.b, in.n1)

						# covariance for b and n1
						Sigma.n1    <- compute_cov(object$cov, object$theta, D[in.pair.n1,in.pair.n1])
						invSigma.n1 <- chol2inv(chol(Sigma.n1))

						J_0 <- J_0 + 2 * B_0 %*% Sigma.n1[1:n.in.b,n.in.b+1:n.in.n1] %*% invSigma.n1[n.in.b+1:n.in.n1,1:n.in.new]
					}
					sd.pred[in.new-nFit] <- sqrt(diag( chol2inv(chol(H_0 %*% chol2inv(chol(J_0)) %*% H_0)) ))
				}
			}
		}
	}

	if (interval == "prediction") {
		tile <- qnorm((1-level)/2, lower.tail=FALSE)
		return(data.frame(y=y_0, lwr=y_0-tile*sd.pred, upr=y_0+tile*sd.pred))
	} else {
		return(data.frame(y=y_0))
	}
}
