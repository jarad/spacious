if ("package:spacious" %in% search()) {
  # unload the package
  detach("package:spacious", unload=TRUE)
}

require(spacious)
require(testthat)

# test exponential fitting
test_that("simple exponential", {
	set.seed(311)
	n <- 100
	S <- matrix(runif(n*2), n, 2)
	D <- rdist(S); diag(D) <- 0
	Sigma <- 0.5 * diag(n) + 0.5*exp(-D/0.75)
	y <- mvrnorm(1, rep(5,n), Sigma)
	# full likelihood
	fit <- spacious(y~1, S=S, cov="exp", blocks=list(type="regular", nblocks=0))
	summary(fit)
	predict(fit, newS=cbind(0.25,0.35), interval="prediction")
	predict(fit, newS=cbind(0.25,0.35), interval="prediction", opts=list(type="local", num=5))
	# block
	fit <- spacious(y~1, S=S, cov="exp", blocks=list(type="regular", nblocks=2^2))
	summary(fit)
	predict(fit, newS=cbind(0.25,0.35), interval="prediction")
	predict(fit, newS=cbind(0.25,0.35), interval="prediction", opts=list(type="local", num=5))
})
