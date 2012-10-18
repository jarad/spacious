# function to fit block composite models
"spacious" <- function(
	formula, data, S=NULL,                  # input data
	cov="exp", cov.inits=NULL,              # covariance function
	B=NULL, neighbors=NULL, D=NULL,
	fixed=list(smoothness=0.5),             # fixed parameters
	nblocks=1, grid.type="regular",         # blocking style
	verbose=FALSE, tol=1e-3, maxIter=100    # algorithm control params
) {
# TODO: error check inputs

	if (is.null(S)) {
		stop("No spatial locations S specified")
	}

	# construct response y and model matrix X from formula
	mf <- match.call(expand.dots = FALSE)
	m <- match(c("formula", "data"), names(mf), 0L)
	mf <- mf[c(1L, m)]
	mf$drop.unused.levels <- TRUE
	mf[[1L]] <- as.name("model.frame")
	mf <- eval(mf, parent.frame())
	mt <- attr(mf, "terms")

	y <- model.response(mf, "numeric")
	X <- model.matrix(mt, mf, NULL)

	# information about model matrix
	n <- nrow(X)
	p <- ncol(X)

	if (is.null(D)) {
		# create distance matrix from spatial locations
		D <- rdist(S)
		D[row(D)==col(D)] <- 0
	}

	# setup based on covariance types
	R <- NA    # number of covariance parameters
	theta <- c()
	theta.fixed <- c()

	# handle covariance types

# TODO: I suggest you use pmatch, e.g.
# cov = pmatch(cov, c("a","b","exp"))
# if (cov != 3) stop("Unknown covariance function")
# then replace R below with cov

	if (cov == "exp") {
		R <- 3

		theta <- rep(0, R)
		theta.fixed <- rep(FALSE, R)

		if (!is.null(fixed$nugget)) {
			theta[1] <- log(fixed$nugget)
			theta.fixed[1] <- TRUE
		}
		if (!is.null(fixed$psill)) {
			theta[2] <- log(fixed$psill)
			theta.fixed[2] <- TRUE
		}
		if (!is.null(fixed$range)) {
			theta[3] <- -log(fixed$range)
			theta.fixed[3] <- TRUE
		}
	} else if (cov == "matern") {
		R <- 4

		theta <- rep(0, R)
		theta.fixed <- rep(FALSE, R)

		if (!is.null(fixed$nugget)) {
			theta[1] <- log(fixed$nugget)
			theta.fixed[1] <- TRUE
		}
		if (!is.null(fixed$psill)) {
			theta[2] <- log(fixed$psill)
			theta.fixed[2] <- TRUE
		}
		if (!is.null(fixed$range)) {
			theta[3] <- -log(fixed$range)
			theta.fixed[3] <- TRUE
		}
		if (!is.null(fixed$smooth)) {
			theta[4] <- log(fixed$smooth)
			theta.fixed[4] <- TRUE
		}
	} else {
		stop(paste("Unknown covariance type",cov))
	}


if (0) {
# TODO: find a smart way to set initial values
	# initial values
	# NOTE: theta passed to spacious.fit are unconstrained via log()
	if (is.null(cov.inits)) {
		theta <- rep(0, R)
#theta[1] <- theta[2] <- log(var(y)/2)
#cat("Starting with:",exp(theta),"\n")
	} else {
		theta <- log(cov.inits)
	}
}

	# create grid
	grid <- c()
	if (!is.null(B) && !is.null(neighbors)) {
		# we have block memberships and neighbors, so don't do anything else
	} else if (nblocks == "pair") {
		# pairwise likelihood
		B <- 1:n
		neighbors <- do.call("rbind", sapply(1:(n-1), function(i) {
			cbind(i,(i+1):n)
		}))
	} else if (nblocks <= 1) {
		B <- rep(1,n)
		neighbors <- matrix( c(1, 1), nrow=1 )
	} else if (grid.type=="regular") {
		B <- rep(NA, n)   # vector to hold block memberships
		neighbors <- c()
		snb <- sqrt(nblocks)

		# ensure that nblocks is an integer squared
		if (snb != round(snb)) {
			stop("Number of blocks (nblocks) must be an integer squared")
		}

		# construct a bounding square
		s.x <- s.y <- c(floor(min(S)), ceiling(max(S)))

		# construct a regular grid with snb rows and snb columns
		spacing <- (s.x[2]-s.x[1])/snb

		# create grid
		b <- 1
		for (i in seq(s.x[1],s.x[2]-spacing,by=spacing)) {
			for (j in seq(s.y[1],s.y[2]-spacing,by=spacing)) {
				poly.x <- round(c(i,i+spacing,i+spacing,i,i), 4)
				poly.y <- round(c(j,j,j+spacing,j+spacing,j), 4)
				in_poly <- point.in.polygon(S[,1], S[,2], poly.x, poly.y) >= 1
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
# TODO: warn if too high a % of data is in a single block

# TODO: move computation of D (if needed) to spacious.fit()
	t1 <- proc.time()
	fit <- spacious.fit(y, X, D, nblocks, B, neighbors, cov, n, p, R, theta, theta.fixed,
		verbose, tol, maxIter)
	t <- proc.time()-t1

	# construct output fit
	fit$time  <- t
	fit$beta  <- as.vector(fit$beta)
	fit$theta <- as.vector(fit$theta)
	fit$theta.fixed <- theta.fixed

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
