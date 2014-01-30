# test exponential fitting
library(spacious)
library(testthat)
source("custom.R")

# simulate exponential data
set.seed(311)
n     <- 512
np    <- 100
S     <- matrix(runif((n+np)*2), n+np, 2)
D     <- rdist(S); diag(D) <- 0
Sigma <- 0.5 * diag(n+np) + 0.5*exp(-D/0.15)
y     <- mvrnorm(1, rep(5,n+np), Sigma)

X <- matrix(1, nrow=length(y), ncol=1)
x1 <- round(rnorm(n+np),3) # noise predictor

data.fit  <- data.frame(y=y[1:n], x1=x1[1:n])
data.pred <- data.frame(y=y[n+1:np], x1=x1[n+1:np])

S.fit  <- S[1:n,]
S.pred <- S[n+1:np,]

# get initial values
inits <- spacious:::initial.theta(data.fit$y, S.fit)
cov.inits <- list(nugget=inits[1], psill=inits[2], range=inits[3])

test_that("full exponential converges", {
	fit <- spacious(y~x1, S=S.fit, cov="exp", blocks=list(type="regular", nblocks=0), data=data.fit, cov.inits=cov.inits)
	expect_equal(fit$conv, TRUE)
})

test_that("block exponential converges", {
	fit <- spacious(y~x1, S=S.fit, cov="exp", blocks=list(type="regular", nblocks=2^2), data=data.fit, cov.inits=cov.inits)
	expect_equal(fit$conv, TRUE)
})

test_that("exponential threads equals non-threads", {
	fit.N <- spacious(y~x1, S=S.fit, cov="exp", blocks=list(type="regular", nblocks=2^2), data=data.fit, cov.inits=cov.inits, nthreads=1)
	fit.T <- spacious(y~x1, S=S.fit, cov="exp", blocks=list(type="regular", nblocks=2^2), data=data.fit, cov.inits=cov.inits, nthreads=2)
	expect_equal(fit.N$ll, fit.T$ll)
})

test_that("exponential gpu equals non-gpu", {
	fit.N <- spacious(y~x1, S=S.fit, cov="exp", blocks=list(type="full"), data=data.fit, cov.inits=cov.inits, gpu=FALSE)
	fit.G <- spacious(y~x1, S=S.fit, cov="exp", blocks=list(type="full"), data=data.fit, cov.inits=cov.inits, gpu=TRUE)
	expect_equal(fit.N$ll, fit.G$ll)
})

test_that("exponential full engine is the same in R and C++", {
	fit.R <- spacious(y~x1, S=S.fit, cov="exp", blocks=list(type="full"), data=data.fit, cov.inits=cov.inits, engine="R")
	fit.C <- spacious(y~x1, S=S.fit, cov="exp", blocks=list(type="full"), data=data.fit, cov.inits=cov.inits, engine="C")
	expect_equal(fit.R$ll, fit.C$ll)

	preds.R <- predict(fit.C, newdata=data.pred, newS=S.pred, interval="prediction", level=0.9, engine="R")
	preds.C <- predict(fit.C, newdata=data.pred, newS=S.pred, interval="prediction", level=0.9, engine="C")

	for (i in 1:np) {
		expect_equal(preds.R$y[i], preds.C$y[i], label=paste0("y for ",i))
		expect_equal(preds.R$sd[i], preds.C$sd[i], label=paste0("sd for ",i))
	}
})

test_that("exponential blocks engine is the same in R and C++", {
	fit.R <- spacious(y~x1, S=S.fit, cov="exp", blocks=list(type="regular", nblocks=2^2), data=data.fit, cov.inits=cov.inits, engine="R")
	fit.C <- spacious(y~x1, S=S.fit, cov="exp", blocks=list(type="regular", nblocks=2^2), data=data.fit, cov.inits=cov.inits, engine="C")
	expect_equal(fit.R$ll, fit.C$ll)

	preds.R <- predict(fit.C, newdata=data.pred, newS=S.pred, interval="prediction", level=0.9, engine="R")
	preds.C <- predict(fit.C, newdata=data.pred, newS=S.pred, interval="prediction", level=0.9, engine="C")

	for (i in 1:np) {
		expect_equal(preds.R$y[i], preds.C$y[i], label=paste0("y for ",i))
		expect_equal(preds.R$sd[i], preds.C$sd[i], label=paste0("sd for ",i))
	}
})
