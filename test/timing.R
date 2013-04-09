# compute spacious timings

if ("package:spacious" %in% search()) {
	# unload the package
	detach("package:spacious", unload=TRUE)
	library.dynam.unload("spacious", libpath="~/Rlib/spacious")
}

# load the package
require(spacious)
require(multicore)

options(cores=3)

mclapply(seq(20, 60, by=10), function(i) {
	S <- as.matrix( expand.grid(seq(0,1,length=i), seq(0,1,length=i)) )
	D <- rdist(S)
	diag(D) <- 0
	n <- nrow(S)

	# covariance function parameters
	nugget <- 0.5
	tau2   <- 0.5
	range  <- 0.25
	mu <- 5
	Sigma <- nugget * diag(n) + tau2 * exp(-D/range)

	# generate data
	y <- round(mvrnorm(1, mu=rep(mu, n), Sigma=Sigma), 3)

	# 25 obs per block
	fit <- spacious(y~1, S=S, cov="exp", blocks=list(type="regular", nblocks=floor(sqrt(n/25))^2 ) )
	cat("n=", n, "; 25 obs per block: ", fit$time[3], "\n", sep="")
	# 100 obs per block
	fit <- spacious(y~1, S=S, cov="exp", blocks=list(type="regular", nblocks=floor(sqrt(n/100))^2 ) )
	cat("n=", n, "; 100 obs per block: ", fit$time[3], "\n", sep="")
	# full
	fit <- spacious(y~1, S=S, cov="exp", blocks=list(type="full"))
	cat("n=", n, "; full lik: ", fit$time[3], "\n", sep="")
})
