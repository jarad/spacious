# prototype for block composite fitting

# function to fit block composite models
"spacious" <- function(
	y, X, S,                                # input data
	cov="exp", cov.inits=NULL,              # covariance function
	nblocks=1, grid.type="regular"          # blocking style
) {
# TODO: move y/X to formula form
# TODO: let user specify block memberships B

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
		grid <- c()
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
pdf("grid.pdf");plot(grid);points(S[,1],S[,2],pch='.');graphics.off();

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
	} else {
		theta <- log(inits$theta)
	}

	fit <- spacious.fit(y, X, D, B, neighbors, cov, n, p, R, theta)

	class(fit) <- "spacious"

	fit
}

"spacious.fit" <- function(y, X, D, B, neighbors, cov, n, p, R, theta) {
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
		b[seq.p] <<- 0
		apply(neighbors, 1, function(row) {
			in.pair <- B==row[1] | B==row[2]
			n.pair <- sum(in.pair)

# TODO: create a function to compute Sigma based on covariance type
			Sigma <- exp(theta[1])*diag(n.pair) + exp(theta[2])*exp(-exp(theta[3])*D[in.pair,in.pair])
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
		u[seq.R] <- 0
		H[seq.RH] <- 0
		FI[seq.R2] <- 0

		apply(neighbors, 1, function(row) {
			in.pair <- B==row[1] | B==row[2]
			n.pair <- sum(in.pair)

# TODO: figure out how to do this once since it's done in beta update as well (maybe merge stuff?)
			Sigma <- exp(theta[1])*diag(n.pair) + exp(theta[2])*exp(-exp(theta[3])*D[in.pair,in.pair])
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
			for (r in 1:R) {
				for (s in r:R) {
					H[index] <<- H[index] + 0.5 * sum(diag( W[[r]] %*% W[[s]] ))
					FI[r,s] <<- H[index]
					if (r != s) {
						FI[s,r] <<- H[index]
					}
					index <- index+1
				}
			}
			invFI <<- chol2inv(chol(FI))
#			invFI <- qr.solve(qr(FI))
		})

		theta + invFI %*% u
	}

	# function to compute an iteration of fisher scoring
	"iter" <- function(prev.beta, prev.theta) {
		list(beta=beta, theta=theta)
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

	list(beta=beta, theta=exp(theta))
}

"summary.spacious" <- function(object, ...) {
	cat("Spacious summary!\n")
}

if (1) {

# test spacious()
require(MASS)
require(fields)
require(geoR)
require(sp)
require(spdep)

set.seed(311)

# generate data to use for fitting a block composite model
n <- 200

# generate spatial locations S
S <- cbind(runif(n), runif(n))

# distance matrix
D <- rdist(S); D[row(D)==col(D)] <- 0

# matern covariance function parameters
nugget <- 0.5
tau2 <- 0.5
range <- 1.5
smooth <- 0.50  # 0.5 = exponential cov
mu <- 5
#Sigma <- nugget * diag(n) + tau2 * matern(D, range, smooth)
Sigma <- nugget * diag(n) + tau2 * exp(-range * D)

# generate data
y <- mvrnorm(1, mu=rep(mu, n), Sigma=Sigma)

"ll" <- function(beta,theta) {
	Sigma <- theta[1] * diag(n) + theta[2] * exp(-theta[3] * D)
	invSigma <- chol2inv(chol(Sigma))
	#-0.5 * log( det(Sigma) ) -0.5 * t(y-X%*%beta) %*% invSigma %*% (y-X%*%beta)
	-0.5 * det(Sigma,log=TRUE) -0.5 * t(y-X%*%beta) %*% invSigma %*% (y-X%*%beta)
}

if (1) {
# fit with fisher scoring
X <- matrix(1, nrow=length(y), ncol=1)
time.fish <- proc.time()
fit.fish <- spacious(y, X, S, cov="exp", nblocks=1^2)
time.fish <- proc.time() - time.fish
beta.fish <- fit.fish$beta
theta.fish <- fit.fish$theta
ll.fish <- ll(beta.fish,theta.fish)

cat("fish estimates:",beta.fish,theta.fish,"\n")
cat("fish ll:",ll.fish,"\n")
cat("fish execution time:\n")
print(time.fish)

}

if (1) {
# try likfit
gd <- as.geodata(cbind(S,y))
time.likfit <- proc.time()
fit.likfit <- likfit(gd, ini.cov.pars=c(1.22,1.22))
time.likfit <- proc.time() - time.likfit

beta.likfit <- fit.likfit$beta
theta.likfit <- c(fit.likfit$tausq, fit.likfit$sigmasq, 1/fit.likfit$phi)
ll.likfit <- ll(beta.likfit,theta.likfit)

cat("likfit estimates:",beta.likfit,theta.likfit,"\n")
cat("likfit ll:",ll.likfit,"\n")
cat("likfit execution time:\n")
print(time.likfit)

}

if (0) {
pdf("lik.pdf");
plot(x<-seq(0,max(theta.likfit[3],theta.fish[3])+2.50,by=0.25),sapply(x, function(x) { ll(beta.likfit,c(theta.likfit[1],theta.likfit[2],x)) }),type="l")
abline(v=theta.fish[3])
abline(v=theta.likfit[3],lty=2)
graphics.off();
}


}
