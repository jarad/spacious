# function to fit block composite models
"spacious" <- function(
	formula, data, S,                       # input data
	cov="exp", cov.inits,                   # covariance function
	B, neighbors,
	fixed=list(smoothness=0.5),             # fixed parameters
	blocks=list(type="cluster"),            # blocking style
	verbose=FALSE, tol=1e-8, maxIter=25,    # algorithm control params
	compute_se=FALSE, compute_diag=FALSE,
	nthreads=1, gpu=FALSE,
	engine="C"                              # use C or R implementation?
) {

	if (missing(S)) {
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

	# setup based on covariance types
	R <- NA    # number of covariance parameters
	theta <- c()
	theta_fixed <- c()

	# handle covariance types
	if (cov == "exp") {
		R <- 3

		theta <- rep(0, R)
		theta_fixed <- rep(FALSE, R)

		if (!is.null(fixed$nugget)) {
			theta[1] <- fixed$nugget
			theta_fixed[1] <- TRUE
		}
		if (!is.null(fixed$psill)) {
			theta[2] <- fixed$psill
			theta_fixed[2] <- TRUE
		}
		if (!is.null(fixed$range)) {
			theta[3] <- fixed$range
			theta_fixed[3] <- TRUE
		}
	} else if (cov == "matern") {
		R <- 4

		theta <- rep(0, R)
		theta_fixed <- rep(FALSE, R)

		if (!is.null(fixed$nugget)) {
			theta[1] <- fixed$nugget
			theta_fixed[1] <- TRUE
		}
		if (!is.null(fixed$psill)) {
			theta[2] <- fixed$psill
			theta_fixed[2] <- TRUE
		}
		if (!is.null(fixed$range)) {
			theta[3] <- fixed$range
			theta_fixed[3] <- TRUE
		}
		if (!is.null(fixed$smooth)) {
			theta[4] <- fixed$smooth
			theta_fixed[4] <- TRUE
		} else {
			stop("Smoothness must be fixed")
		}
	} else {
		stop(paste("Unknown covariance type",cov))
	}

	if (sum(!theta_fixed) > 0) {
		# compute some initial values
		which.inits <- which(theta_fixed[1:3] == FALSE)
		if (length(which.inits) > 0) {
			theta[which.inits] <- (initial.theta(y, S))[which.inits]
		}
	}

	if (!missing(cov.inits)) {
		# set initial values from user
		if (!is.null(cov.inits$nugget)) {
			theta[1] <- cov.inits$nugget
		}
		if (!is.null(cov.inits$psill)) {
			theta[2] <- cov.inits$psill
		}
		if (!is.null(cov.inits$range)) {
			theta[3] <- cov.inits$range
		}
		if (cov == "matern" && !is.null(cov.inits$smooth)) {
			theta[4] <- cov.inits$smooth
		}
	}

	if (verbose) {
		cat("Initial values for theta:",theta,"\n")
	}

	# create grid
	grid <- c()
	lik_form <- "block"     # type of likelihood we want to use

	if (!missing(B) && !missing(neighbors)) {
		# we have block memberships and neighbors, so don't do anything else
	} else if (blocks$type == "full" || (!is.null(blocks$nblocks) && blocks$nblocks <= 2) ) {
		# full likelihood
		B         <- rep(1,n)
		neighbors <- matrix( c(1, 1), nrow=1 )
		lik_form <- "full"
	} else if (blocks$type == "cluster") {
		if (is.null(blocks$nblocks)) {
			# deafult to one block for each 50 obs
			blocks$nblocks <- round(n/50)
		}

		bc        <- blocks.cluster(S, blocks$nblocks)
		B         <- bc$B
		grid      <- bc$grid
		neighbors <- bc$neighbors
	} else if (blocks$type == "pair") {
		# pairwise likelihood
		B <- 1:n
		neighbors <- do.call("rbind", sapply(1:(n-1), function(i) {
			cbind(i,(i+1):n)
		}))
		lik_form <- "pair"
	} else if (blocks$type=="regular") {
		if (is.null(blocks$nblocks)) {
			# deafult to one block for each 50 obs
			blocks$nblocks <- round(sqrt(n/50))^2
		}

		B         <- rep(NA, n)   # vector to hold block memberships
		neighbors <- c()
		snb       <- sqrt(blocks$nblocks)

		# ensure that blocks$nblocks is an integer squared
		if (snb != round(snb)) {
			stop("Number of blocks (blocks$nblocks) must be an integer squared")
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

				if (sum(in_poly) == 0) { next; }

				B[in_poly] <- b

				# save grid
				grid <- c(grid,list(Polygons(list(Polygon(cbind(
					poly.x,poly.y
				))),paste(b)) ))

				b <- b+1
			}
		}
		grid <- SpatialPolygons(grid)

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

	# do fit
	t1 <- proc.time()
	if (engine == "C") {
		if (nthreads <= 0) {
			warning(paste0("nthreads specified as ",nthreads,". Using nthreads=1 instead."))
			nthreads <- 1
		}

		fit <- .C(spacious_fit,
		          # input data
		          y=as.double(y), X=as.double(X), S=as.double(S), B=as.integer(B-1), neighbors=as.integer(neighbors-1),
		          n=as.integer(n), p=as.integer(p), nblocks=as.integer(blocks$nblocks), npairs=as.integer(nrow(neighbors)),
		          # type of model to fit
		          lik_form=as.character(lik_form), cov=as.character(cov),
		          # parameter estimates and convergence info
		          theta=as.double(theta), theta_fixed=as.integer(theta_fixed), beta=as.double(rep(0, p)),
		          convergence=as.logical(FALSE), nIter=as.integer(0),
		          # standard errors
		          se_beta=as.double(rep(NA, p)), se_theta=as.double(rep(NA, R)),
		          vcov_beta=as.double(rep(NA, p*p)), vcov_theta=as.double(rep(NA, R*R)),
		          # fitted values and residuals
		          fitted=as.double(rep(0, n)), resids=as.double(rep(0, n)),
		          # values of theta and log likelihood at each iteration
		          iters_theta=as.double(rep(NA, maxIter)), iters_ll=as.double(rep(NA, maxIter)),
		          # fitting control parameters
		          verbose=as.logical(verbose), tol=as.double(tol), max_iter=as.integer(maxIter), compute_se=as.logical(compute_se),
		          # parallelization options
		          nthreads=as.integer(nthreads), gpu=as.logical(gpu),
		          NAOK=TRUE
		)

		fit$engine <- "C"

		# cleanup a few things
		fit$X <- NULL
		fit$n <- NULL
		fit$p <- NULL
		fit$npairs <- NULL
	} else if (engine == "R") {
		fit <- spacious.fit(y, X, S, blocks$nblocks, B, neighbors, cov, n, p, R, theta, theta_fixed,
			verbose, tol, maxIter)

		fit$engine      <- "R"
		fit$lik_form    <- lik_form
		fit$cov         <- cov
		fit$beta        <- as.vector(fit$beta)
		fit$theta       <- as.vector(fit$theta)
		fit$theta_fixed <- theta_fixed

		fit$nblocks <- blocks$nblocks
		fit$fitted  <- X %*% fit$beta
		fit$resids  <- y - fit$fitted
	} else {
		stop(paste0("Unknown engine ",engine,"\n"))
	}
	t <- proc.time()-t1

	# construct output fit
	fit$time      <- t

	fit$y         <- y
	fit$S         <- S
	fit$B         <- B
	fit$grid      <- grid
	fit$neighbors <- neighbors
	fit$cov       <- cov

	fit$terms  <- attr(mf, "terms")

	class(fit) <- "spacious"

	fit
}
