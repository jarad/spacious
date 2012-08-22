# predict values at new sites using fit from spacious block composite model
"predict.spacious" <- function(object, newdata=NULL, newS=NULL,
	interval = c("none", "confidence", "prediction"), level = 0.95, ...) {

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
	y0 <- rep(NA, nNew)

	# figure out newB the location is in
	if (nB == 0) {  # no blocks
		# perform traditional kriging
		n <- nrow(S)

		# compute covariance matrix
		Sigma <- compute_cov(object$cov, object$theta, D)

		# get the predictions
		y0[1:nNew] <- X %*% object$beta + Sigma[(nFit+1):n,1:nFit] %*%
			chol2inv(chol(Sigma[1:nFit,1:nFit])) %*% (object$y - object$fitted)
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
			y0[1:nNew] <- X %*% object$beta
cat("Initial y0:",y0,"\n")

			# number of unique blocks we have data points for
			uB <- unique(newB)

			for (b in uB) {
				# neighbors of this block
				newNeighbors <- which( rowSums( object$neighbors==b ) == 1 )
cat("Block:",b,"\n")
cat("Neighbors:",newNeighbors,"\n")

				# what new points are in this block?
				newInBlock <- which(newB == b)
print(sum(newInBlock))
done

				for (pair in newNeighbors) {
					# compute covariance matrix
					#Sigma <- compute_cov(object$cov, object$theta, D[,])
				}

			}
		}
	}

	return( as.vector(y0) )
}
