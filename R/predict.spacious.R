# predict values at new sites using fit from spacious block composite model
"predict.spacious" <- function(object, newdata, newS, newB, D,
	opts = list(type="block"),
	interval = "none", level = 0.95,
	nthreads = 1, gpu=FALSE,
	engine, # use the C or R implementation?
	...) {

	if (missing(newS)) {
		stop("Must specify locations newS to make predictions at.")
	}

	if (!object$convergence) {
		warning("Predictions being made with fit that did not converge!")
	}

	nFit <- nrow(object$S)
	nNew <- nrow(newS)

	# build model matrix for the new data
	if (missing(newdata)) {
		if (length(object$beta) > 1) {
			# newdata missing, yet we have covariates...
			stop("newdata must be specified for this fit.")
		}

		# use intercept only
		X <- matrix(1, nrow=nNew, ncol=1)
	} else {
		if (is.data.frame(newdata)) {
			terms <- delete.response(terms(object))
			mf    <- model.frame(terms, newdata)
			X     <- model.matrix(terms, mf)
		} else if (is.matrix(newdata)) {
			# assume newdata is a matrix of covariates
			X <- cbind(1, newdata)
		} else {
			stop(paste0("Unable to handle newdata with class ", class(newdata)))
		}
	}

	# ensure this is a matrix
	newS <- matrix(newS, ncol=2)

	# locations for all points
	S <- rbind(object$S, newS)

	if (missing(engine)) {
		engine <- object$engine
	}

	if (object$lik_form == "block" & opts$type == "block") {
		# number of blocks
		nB   <- length(object$grid)

		# figure out newB the location is in
		if (!missing(newB)) {
			# make sure we have these blocks defined in the neighbors
			uniqueBlocks <- unique( as.vector( object$neighbors ) )
			sapply(newB, function(b) {
				if (sum(uniqueBlocks == b) == 0) {
					stop(paste0("No neighbors defined for block ",b,". Unable to make predictions."))
				}
			})
		} else {
			# user did not define blocks, so find some

			if (is.null(object$grid)) {
				# we don't have polygons, so find the closest point and use that block
				stop("No polygons defining blocks. Specify newB to make predictions.")
			}

			# figure out newBs contain these new locations
			newB <- rep(NA, nNew)

			for (b in 1:nB) {
				in_poly <- point.in.polygon(newS[,1], newS[,2],
					object$grid@polygons[[b]]@Polygons[[1]]@coords[,1],
					object$grid@polygons[[b]]@Polygons[[1]]@coords[,2]) >= 1

				newB[in_poly] <- b
			}

			if (sum(is.na(newB)) > 0) {
				warning("Unable to find blocks for all prediction locations. Assigning locations to blocks using block centroids.")

				which.naB <- which( is.na(newB) )   # which obs are missing blocks?
				N.naB     <- length(which.naB)      # how many of these are there?

				# get centers of blocks
				block.centers <- coordinates(object$grid)

				# compute distances between block centers and obs with missing blocks
				D.centers <- matrix( (rdist(rbind(newS[which.naB,], block.centers)))[1:N.naB, N.naB + 1:nrow(block.centers)], nrow=N.naB )

				# find block for each missing
				foundB <- apply(D.centers, 1, function(row) {
					which.min(row)
				})

				# updated new blocks
				newB[which.naB] <- foundB
			}
		}
	} else {
		nB <- 0
	}

	if (opts$type=="local") engine <- "R"  # C doesn't support local kriging (yet)

	if (engine == "C") {
		do_sd <- interval=="prediction"
		local <- opts$type=="local"

		if (missing(newB)) newB <- rep(1, nNew)
		preds <- .C(spacious_predict,
		          # new site information
		          n_0=as.integer(nNew), y_0=as.double(rep(NA, nNew)), newS=as.double(newS), newB=as.integer(newB-1),
		          newX=as.double(X), do_sd=as.logical(do_sd), sd=as.double(rep(NA, nNew)),
		          local=as.logical(local), Nlocal=as.integer(opts$num),
		          # data used in fit
		          y=as.double(object$y), X=as.double(object$X), S=as.double(object$S), B=as.integer(object$B-1),
		          neighbors=as.integer(object$neighbors-1), n=as.integer(length(object$y)), p=as.integer(ncol(object$X)),
		          nblocks=as.integer(object$nblocks), npairs=as.integer(nrow(object$neighbors)),
		          # type of fit
		          lik_form=as.character(object$lik_form), cov=as.character(object$cov),
		          # fitted values
		          beta=as.double(object$beta), theta=as.double(object$theta),
		          # parallelization options
		          nthreads=as.integer(nthreads), gpu=as.logical(gpu),
		          NAOK=TRUE
		)

		y_0 <- preds$y_0
		if (do_sd) {
 			sd.pred <- preds$sd
			tile <- qnorm((1-level)/2, lower.tail=FALSE)
			return(data.frame(y=y_0, sd=sd.pred, lwr=y_0-tile*sd.pred, upr=y_0+tile*sd.pred))
		} else {
			return(data.frame(y=y_0))
		}
	} else if (engine == "R") {
		if (missing(D)) {
			# compute distance matrix
			D <- rdist(S)
			diag(D) <- 0
		}

		# vector to hold predictions
		y_0      <- rep(NA, nNew)

		if (interval == "prediction") {
			# vector to hold prediction std devs
			sd.pred <- rep(NA, nNew)
		}

		# spatial params in terms of model structure
		theta <- object$theta

		if ( (nB == 0 && missing(newB)) || (opts$type == "all") ) {
			# use all data to perform predictions and perform traditional kriging
			n <- nrow(S)

			# compute covariance matrix
			Sigma <- compute_cov(object$cov, theta, D)

			# get the predictions
			y_0[1:nNew] <- X %*% object$beta + Sigma[nFit+1:nNew,1:nFit] %*%
				chol2inv(chol(Sigma[1:nFit,1:nFit])) %*% object$resids

			if (interval == "prediction") {
				invSigma <- chol2inv(chol(Sigma))

				sd.pred <- sqrt(diag( chol2inv(chol( invSigma[nFit+1:nNew,nFit+1:nNew] )) ))
			}
		} else {
			if (opts$type == "block") {
				# predict when we have blocks

				# block memberships for all points
				B <- c(object$B, newB)

				# predict for each unique block
				y_0[1:nNew] <- X %*% object$beta

				# unique blocks where we have points to predict
				uB <- unique(newB)

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

					for (n1 in neighbors) {
						in.n1   <- which(object$B == n1)
						n.in.n1 <- length(in.n1)

						# take points so that we have Y = (Y_{b,new}, Y_{b,obs}, Y_{n1})
						in.pair.n1 <- c(in.b, in.n1)

						# covariance for b and n1
						Sigma.n1    <- compute_cov(object$cov, theta, D[in.pair.n1,in.pair.n1])
						invSigma.n1 <- chol2inv(chol(Sigma.n1))

						A_0 <- A_0 + invSigma.n1[1:n.in.new,1:n.in.new]
						b_0 <- b_0 +
							matrix(invSigma.n1[1:n.in.new,n.in.new+1:n.in.obs],nrow=n.in.new,ncol=n.in.obs) %*% object$resids[in.obs] +
							matrix(invSigma.n1[1:n.in.new,n.in.new+n.in.obs+1:n.in.n1],nrow=n.in.new,ncol=n.in.n1) %*% object$resids[in.n1]

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
								Sigma.n2    <- compute_cov(object$cov, theta, D[in.pair.n2,in.pair.n2])
								invSigma.n2 <- chol2inv(chol(Sigma.n2))

								# covariance for n1 and n2
								Sigma.n1n2 <- compute_cov(object$cov, theta, D[in.n1n2,in.n1n2])

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

						Sigma <- compute_cov(object$cov, theta, D[in.b,in.b])
						J_0   <- J_0 + B_0 %*% Sigma %*% t(B_0)

						for (n1 in neighbors) {
							in.n1   <- which(object$B == n1)
							n.in.n1 <- length(in.n1)

							# take points so that we have Y = (Y_{b,new}, Y_{b,obs}, Y_{n1})
							in.pair.n1 <- c(in.b, in.n1)

							# covariance for b and n1
							Sigma.n1    <- compute_cov(object$cov, theta, D[in.pair.n1,in.pair.n1])
							invSigma.n1 <- chol2inv(chol(Sigma.n1))

							J_0 <- J_0 + 2 * B_0 %*% Sigma.n1[1:n.in.b,n.in.b+1:n.in.n1] %*% invSigma.n1[n.in.b+1:n.in.n1,1:n.in.new]
						}

						sd.pred[in.new-nFit] <- sqrt(diag( chol2inv(chol(H_0 %*% chol2inv(chol(J_0)) %*% H_0)) ))
					}
				}

			} else if (opts$type == "local") {
				nLocal <- opts$num

				if (is.null(opts$num)) {
					nLocal <- min(100, nFit)   # default number of local points to use
					warning(paste0("Number of points to use in local kriging not specified. Using ",nLocal," points"))
				}

				for (i in 1:nNew) {
					# use closest nLocal points for kriging
					pred.index <- i + nFit

					# the closest nLocal points to Snew[i,]
					which.local <- sort( D[pred.index, 1:nFit], decreasing=FALSE, index.return=TRUE )$ix[1:nLocal]

					# compute covariance matrix
					Sigma <- compute_cov(object$cov, theta, D[c(which.local, pred.index), c(which.local, pred.index)])

					# get the predictions
					y_0[i] <- X[i,] %*% object$beta + Sigma[nLocal+1,1:nLocal] %*%
						chol2inv(chol(Sigma[1:nLocal,1:nLocal])) %*% object$resids[which.local]

					if (interval == "prediction") {
						invSigma <- chol2inv(chol(Sigma))

						sd.pred[i] <- sqrt( 1/invSigma[nLocal+1,nLocal+1] )
					}
				}
			} else {
				stop(paste0("Unknown prediction type: ",opts$type))
			}
		} # end prediction

		if (interval == "prediction") {
			tile <- qnorm((1-level)/2, lower.tail=FALSE)
			return(data.frame(y=y_0, sd=sd.pred, lwr=y_0-tile*sd.pred, upr=y_0+tile*sd.pred))
		} else {
			return(data.frame(y=y_0))
		}

	} else {
		error(paste0("Unknown engine: ",engine,"\n"))
	}

}
