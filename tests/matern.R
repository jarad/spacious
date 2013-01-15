if ("package:spacious" %in% search()) {
  # unload the package
  detach("package:spacious", unload=TRUE)
}

require(spacious)
require(testthat)

# test matern fitting
test_that("simple matern", {

	set.seed(311)

	# generate data
	n <- 250
	np <- 5

	S <- cbind(runif(n+np), runif(n+np))
	D <- rdist(S); diag(D) <- 0

	# matern covariance function parameters
	nugget <- 0.5
	tau2 <- 0.5
	range <- 0.25
	smooth <- 3.50  # 0.5 = exponential cov
	mu <- 5
	Sigma <- nugget * diag(n+np) + tau2 * matern(D, range, smooth)

	y <- mvrnorm(1, mu=rep(mu, n+np), Sigma=Sigma)

	X <- matrix(1, nrow=length(y), ncol=1)
	x1 <- rnorm(n+np)

	# full
	fit <- spacious(y~x2, data=data.frame(y=y[1:n], x2=x1[1:n]), S=S[1:n,], cov="matern",
		fixed=list(smoothness=smooth), blocks=list(type="regular", nblocks=0))
	summary(fit)
	predict(fit, newdata=data.frame(x2=x1[(n+1):(n+np)]), newS=S[(n+1):(n+np),], interval="prediction")
	predict(fit, newdata=data.frame(x2=x1[(n+1):(n+np)]), newS=S[(n+1):(n+np),], interval="prediction", opts=list(type="local", num=5))
	# blocks
	fit <- spacious(y~x2, data=data.frame(y=y[1:n], x2=x1[1:n]), S=S[1:n,], cov="matern",
		fixed=list(smoothness=smooth), blocks=list(type="regular", nblocks=2^2))
	summary(fit)
	predict(fit, newdata=data.frame(x2=x1[(n+1):(n+np)]), newS=S[(n+1):(n+np),], interval="prediction")
	predict(fit, newdata=data.frame(x2=x1[(n+1):(n+np)]), newS=S[(n+1):(n+np),], interval="prediction", opts=list(type="local", num=5))

})
