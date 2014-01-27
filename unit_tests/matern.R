# test matern fitting
library(spacious)
library(testthat)
source("custom.R")

# simulate matern data
set.seed(311)
n <- 512
S <- matrix(runif(n*2), n, 2)
D <- rdist(S); diag(D) <- 0
smooth <- 3.50
Sigma <- 0.5 * diag(n) + 0.5*matern(D, 0.15, smooth)
y <- mvrnorm(1, rep(5,n), Sigma)

test_that("full matern converges", {
	fit <- spacious(y~1, S=S, cov="matern", fixed=list(smoothness=smooth), blocks=list(type="regular", nblocks=0))
	expect_equal(fit$conv, TRUE)
})

test_that("block matern converges", {
	fit <- spacious(y~1, S=S, cov="matern", fixed=list(smoothness=smooth), blocks=list(type="regular", nblocks=2^2))
	expect_equal(fit$conv, TRUE)
})

test_that("matern threads equals non-threads", {
	fit.N <- spacious(y~1, S=S, cov="matern", fixed=list(smoothness=smooth), blocks=list(type="regular", nblocks=2^2), nthreads=1)
	fit.T <- spacious(y~1, S=S, cov="matern", fixed=list(smoothness=smooth), blocks=list(type="regular", nblocks=2^2), nthreads=2)
	expect_equal(fit.N$ll, fit.T$ll)
})

test_that("matern gpu equals non-gpu", {
	fit.N <- spacious(y~1, S=S, cov="matern", fixed=list(smoothness=smooth), blocks=list(type="full"), gpu=FALSE)
	fit.G <- spacious(y~1, S=S, cov="matern", fixed=list(smoothness=smooth), blocks=list(type="full"), gpu=TRUE)
	expect_equal(fit.N$ll, fit.G$ll)
})

test_that("matern full engine is the same in R and C++", {
	fit.R <- spacious(y~1, S=S, cov="matern", fixed=list(smoothness=smooth), blocks=list(type="full"), engine="R")
	fit.C <- spacious(y~1, S=S, cov="matern", fixed=list(smoothness=smooth), blocks=list(type="full"), engine="C")
	expect_that(fit.R$ll, is_close(fit.C$ll))
})

test_that("matern blocks engine is the same in R and C++", {
	fit.R <- spacious(y~1, S=S, cov="matern", fixed=list(smoothness=smooth), blocks=list(type="regular", nblocks=2^2), engine="R")
	fit.C <- spacious(y~1, S=S, cov="matern", fixed=list(smoothness=smooth), blocks=list(type="regular", nblocks=2^2), engine="C")
	expect_that(fit.R$ll, is_close(fit.C$ll))
})

test_that("matern with smoothness=0.5 gives same as exponential", {
	fit.E <- spacious(y~1, S=S, cov="exp", blocks=list(type="regular", nblocks=2^2), engine="R")
	fit.M <- spacious(y~1, S=S, cov="matern", fixed=list(smoothness=0.5), blocks=list(type="regular", nblocks=2^2), engine="C")
	expect_that(fit.E$ll, is_close(fit.M$ll))
})
