# function to fit block composite models
"spacious" <- function(
	formula, data=NULL, S,                  # input data
	cov="exp", cov.inits=NULL,              # covariance function
	nblocks=1, grid.type="regular",         # blocking style
	verbose=FALSE, tol=1e-3, maxIter=100    # algorithm control params
) {
# TODO: let user specify block memberships B
# TODO: let user specify neighbors

	# construct response y and model matrix X from formula
	mf <- model.frame(formula, data)
	y <- model.extract(mf, "response")
	X <- model.matrix(attr(mf, "terms"), mf, NULL)

	# information about model matrix
	n <- nrow(X)
	p <- ncol(X)

	# create distance matrix from spatial locations
	D <- rdist(S)
	D[row(D)==col(D)] <- 0

	# handle covariance types
	if (cov == "exp") {
		R <- 3
	} else {
		stop(paste("Unknown covariance type",cov))
	}

	# create grid
	B <- rep(NA, n)   # vector to hold block memberships
	neighbors <- c()
	grid <- c()
	if (nblocks == 1) {
		B[1:n] <- 1
		neighbors <- matrix( c(1, 1), nrow=1 )
	} else if (grid.type=="regular") {
		snb <- sqrt(nblocks)

		# ensure that nblocks is an integer squared
		if (snb != round(snb)) {
			stop("Number of blocks (nblocks) must be an integer squared")
		}

		# construct a bounding square
		s.x <- s.y <- c(floor(min(S)), ceiling(max(S[,1])))

		# construct a regular grid with snb rows and snb columns
		spacing <- (s.x[2]-s.x[1])/snb

		# create grid
		b <- 1
		for (i in seq(s.x[1],s.x[2]-spacing,by=spacing)) {
			for (j in seq(s.y[1],s.y[2]-spacing,by=spacing)) {
				poly.x <- c(i,i+spacing,i+spacing,i,i)
				poly.y <- c(j,j,j+spacing,j+spacing,j)
				in_poly <- point.in.polygon(S[,1], S[,2], poly.x, poly.y) == 1
				B[in_poly] <- b

				# save grid
				grid <- c(grid,list(Polygons(list(Polygon(cbind(
					poly.x,poly.y
				))),paste(b)) ))

				b <- b+1
			}
		}
		grid <- SpatialPolygons(grid)
#pdf("grid.pdf");plot(grid);points(S[,1],S[,2],pch='.');graphics.off();

		if (sum(is.na(B)) > 0) {
			stop("Some points not placed in grid.")
		}

		# get neighbors
		neighbor.mat <- nb2mat(poly2nb(grid), style='B', zero.policy=T)
		for (i in 1:nrow(neighbor.mat)) {
			for (j in i:ncol(neighbor.mat)) {
				if (neighbor.mat[i,j] == 1) {
					neighbors <- rbind(neighbors, c(i,j) )
				}
			}
		}
	}

# TODO: find a smart way to set these
	# initial values
	# NOTE: theta passed to spacious.fit are unconstrained via log()
	if (is.null(cov.inits)) {
		theta <- rep(0, R)
#theta[1] <- theta[2] <- log(var(y)/2)
#cat("Starting with:",exp(theta),"\n")
	} else {
		theta <- log(cov.inits$theta)
	}

# TODO: move computation of D to spacious.fit()
	fit <- spacious.fit(y, X, D, nblocks, B, neighbors, cov, n, p, R, theta,
		verbose, tol, maxIter)

	# construct output fit
	fit$beta  <- as.vector(fit$beta)
	fit$theta <- as.vector(fit$theta)

	fit$y         <- y
	fit$S         <- S
	fit$D         <- D
	fit$nblocks   <- nblocks
	fit$B         <- B
	fit$grid      <- grid
	fit$neighbors <- neighbors
	fit$cov       <- cov

	fit$terms  <- attr(mf, "terms")
	fit$fitted <- X %*% fit$beta
	fit$resids <- fit$y - fit$fitted

	class(fit) <- "spacious"

	fit
}
