# predict values at new sites using fit from spacious block composite model
"predict.spacious" <- function(object, newdata=NULL, newS=NULL,
	interval = c("none", "prediction"), level = 0.95, ...) {

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
	y_0 <- rep(NA, nNew)

	# figure out newB the location is in
	if (nB == 0) {  # no blocks
		# perform traditional kriging
		n <- nrow(S)

		# compute covariance matrix
		Sigma <- compute_cov(object$cov, object$theta, D)

		# get the predictions
		y_0[1:nNew] <- X %*% object$beta + Sigma[(nFit+1):n,1:nFit] %*%
			chol2inv(chol(Sigma[1:nFit,1:nFit])) %*% object$resids
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
					fit.spacious$grid@polygons[[b]]@Polygons[[1]]@coords[,1],
					fit.spacious$grid@polygons[[b]]@Polygons[[1]]@coords[,2]) == 1

				newB[in_poly] <- b
			}

			if (sum(is.na(newB))) {
				warning("Unable to find blocks for all prediction sites. Using closest point(s) to find block(s).")
# TODO: do this
			}

			# block memberships for all points
			B <- rbind(object$B, newB)

			# predict for each unique block
			y_0[1:nNew] <- X %*% object$beta

			# unique blocks we have data points for
			uB <- unique(newB)

# TODO: we could re-construct this so that we only invert each matrix once. worth considering if speed more important than memory.
			for (b in uB) {
				# neighbors of this block
				neighbors <- as.vector(object$neighbors[which( rowSums( object$neighbors==b ) == 1 ),])
				neighbors <- neighbors[neighbors != b]

				# what new points are in this block?
				in.new <- which(newB == b)+nFit
				n.in.new <- length(in.new)

				# what observed points are in this block?
				in.obs <- which(object$B == b)
				n.in.obs <- length(in.obs)

				A_0 <- matrix(0, nrow=n.in.new, ncol=n.in.new)
				b_0 <- rep(0, n.in.new)

				for (neighbor in neighbors) {
					in.neighbor <- which(object$B == neighbor)
					n.in.neighbor <- length(in.neighbor)

					# take points so that we have Y = (Y_{b,new}, Y_{b,obs}, Y_{neighbor})
					in.pair <- c(in.new, in.obs, in.neighbor)

					# compute covariance matrix
					Sigma <- compute_cov(object$cov, object$theta, D[in.pair,in.pair])
					invSigma <- chol2inv(chol(Sigma))

					A_0 <- A_0 + invSigma[1:n.in.new,1:n.in.new]
					b_0 <- b_0 + invSigma[1:n.in.new,n.in.new+1:n.in.obs] %*% object$resids[in.obs] +
						invSigma[1:n.in.new,n.in.new+n.in.obs+1:n.in.neighbor] %*% object$resids[in.neighbor]
				}

				b_0 <- -b_0

				y_0[in.new-nFit] <- y_0[in.new-nFit] + chol2inv(chol(A_0)) %*% b_0
			}
		}
	}

	return( as.vector(y_0) )
}
