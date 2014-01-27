# test exponential fitting
library(spacious)
library(testthat)
source("custom.R")

# simulate exponential data
set.seed(311)
n <- 512
S <- matrix(runif(n*2), n, 2)
D <- rdist(S); diag(D) <- 0
Sigma <- 0.5 * diag(n) + 0.5*exp(-D/0.15)
y <- mvrnorm(1, rep(5,n), Sigma)

test_that("full exponential converges", {
	fit <- spacious(y~1, S=S, cov="exp", blocks=list(type="regular", nblocks=0))
	expect_equal(fit$conv, TRUE)
})

test_that("block exponential converges", {
	fit <- spacious(y~1, S=S, cov="exp", blocks=list(type="regular", nblocks=2^2))
	expect_equal(fit$conv, TRUE)
})

test_that("exponential threads equals non-threads", {
	fit.N <- spacious(y~1, S=S, cov="exp", blocks=list(type="regular", nblocks=2^2), nthreads=1)
	fit.T <- spacious(y~1, S=S, cov="exp", blocks=list(type="regular", nblocks=2^2), nthreads=2)
	expect_equal(fit.N$ll, fit.T$ll)
})

test_that("exponential gpu equals non-gpu", {
	fit.N <- spacious(y~1, S=S, cov="exp", blocks=list(type="full"), gpu=FALSE)
	fit.G <- spacious(y~1, S=S, cov="exp", blocks=list(type="full"), gpu=TRUE)
	expect_equal(fit.N$ll, fit.G$ll)
})

test_that("exponential blocks engine is the same in R and C++", {
	fit.R <- spacious(y~1, S=S, cov="exp", blocks=list(type="regular", nblocks=2^2), engine="R")
	fit.C <- spacious(y~1, S=S, cov="exp", blocks=list(type="regular", nblocks=2^2), engine="C")
	expect_that(fit.R$ll, is_close(fit.C$ll))
})

test_that("exponential full engine is the same in R and C++", {
	fit.R <- spacious(y~1, S=S, cov="exp", blocks=list(type="full"), engine="R")
	fit.C <- spacious(y~1, S=S, cov="exp", blocks=list(type="full"), engine="C")
	expect_that(fit.R$ll, is_close(fit.C$ll))
})
