# prototype for block composite fitting

# function to fit block composite models
"spacious" <- function(y, X, D,
	cov="exp", cov.inits=NULL) {

	n <- nrow(X)
	p <- ncol(X)

	if (cov == "exp") {
		R <- 3
	} else {
		stop(paste("Unknown covariance type",cov))
	}

# TODO: find a smart way to set these
	# initial values
	# NOTE: theta passed to spacious.fit are unconstrained via log()
	if (is.null(cov.inits)) {
		theta <- rep(0, R)
	} else {
		theta <- log(inits$theta)
	}

	fit <- spacious.fit(y, X, D, cov, n, p, R, theta)

	class(fit) <- "spacious"

	fit
}

"spacious.fit" <- function(y, X, D, cov, n, p, R, theta) {
	# a list of functions to compute partial derivates
	# with respect to each covariance function param
	partials <- list()

	# setup fit based on covariance function
	if (cov == "exp") {
		partials <- list(
			function(theta) {
				exp(theta[1])*diag(n)
			},
			function(theta) {
				exp(theta[2])*exp(-exp(theta[3]) * D)
			},
			function(theta) {
				-exp(theta[2]+theta[3]) * D * exp(-exp(theta[3]) * D)
			}
		)
	} else {
		stop(paste("Unknown covariance type",cov))
	}

	# function to compute an iteration of fisher scoring
	"iter" <- function(prev.beta, prev.theta) {
# TODO: create a compute Sigma function based on covariance type
		#Sigma <- theta[1] * diag(n) + theta[2] * exp(-theta[3] * D)
		Sigma <- exp(prev.theta[1]) * diag(n) + exp(prev.theta[2]) * exp(-exp(prev.theta[3]) * D)
		cholSigma <- chol(Sigma)
		invSigma <- chol2inv(cholSigma)
#		cholInvSigma <- chol(invSigma)

		# obtain beta
		beta <- solve( t(X) %*% invSigma %*% X ) %*% t(X) %*% invSigma %*% y

		Xb <- X %*% beta
		ymXb <- y-Xb
		q <- invSigma %*% ymXb

		u <- rep(0, R)
		W <- vector("list", R)
		H <- rep(0, R*(R+1)/2)
		FI <- matrix(0, nrow=R, ncol=R)

		# compute the Ws
		for (r in 1:R) {
			partial <- partials[[r]](prev.theta)
			W[[r]] <- invSigma %*% partial
			u[r] <- -0.5 * sum( diag(W[[r]]) ) + 0.5 * t(q) %*% partial %*% q
		}

		index <- 1
		for (r in 1:R) {
			for (s in r:R) {
				H[index] <- 0.5 * sum(diag( W[[r]] %*% W[[s]] ))
				FI[r,s] <- H[index]
				if (r != s) {
					FI[s,r] <- H[index]
				}
				index <- index+1
			}
		}
#		invFI <- chol2inv(chol(FI))
		invFI <- qr.solve(qr(FI))

		theta <- prev.theta + invFI %*% u

		list(beta=beta, theta=theta)
	}

# TODO: allow user to set control parameters
	# estimate params
	beta <- rep(0, p)
	maxIter <- 100
	tol <- 1e-3
	for (i in 1:maxIter) {
		prev.beta <- beta
		prev.theta <- theta

		step <- iter(beta, theta)
		beta <- step$beta
		theta <- step$theta
cat("iter",i,":"); print( c(beta, exp(theta)) )

		if (i > 1) {
			# have we converged?
# TODO: figure out the best way to identify convergence
			if ( sum( c(abs(prev.beta - beta),abs(prev.theta-theta)) ) <= tol ) {
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

set.seed(1983)

# generate data to use for fitting a block composite model
n <- 300

# generate spatial locations S
S <- cbind(runif(n), runif(n))

# distance matrix
D <- rdist(S); D[row(D)==col(D)] <- 0

# matern covariance function parameters
nugget <- 0.01
tau2 <- 0.25
range <- 0.25
smooth <- 0.50  # 0.5 = exponential cov
mu <- 0
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
fit.fish <- spacious(y, X, D, cov="exp")
time.fish <- proc.time() - time.fish
beta.fish <- fit.fish$beta
theta.fish <- fit.fish$theta
ll.fish <- ll(beta.fish,theta.fish)
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
}

cat("fish estimates:",beta.fish,theta.fish,"\n")
cat("fish ll:",ll.fish,"\n")
cat("fish execution time:\n")
print(time.fish)

cat("likfit estimates:",beta.likfit,theta.likfit,"\n")
cat("likfit ll:",ll.likfit,"\n")
cat("likfit execution time:\n")
print(time.likfit)

if (0) {
pdf("lik.pdf");
plot(x<-seq(0,max(theta.likfit[3],theta.fish[3])+2.50,by=0.25),sapply(x, function(x) { ll(beta.likfit,c(theta.likfit[1],theta.likfit[2],x)) }),type="l")
abline(v=theta.fish[3])
abline(v=theta.likfit[3],lty=2)
graphics.off();
}


}
