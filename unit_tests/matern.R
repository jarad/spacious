# test matern fitting
library(spacious)
library(testthat)

# simulate matern data
set.seed(311)
n     <- 512
np    <- 100
S     <- matrix(runif((n+np)*2), n+np, 2)
D     <- rdist(S); diag(D) <- 0
smooth<- 3.50
Sigma <- 0.5 * diag(n+np) + 0.5*matern(D, 0.15, smooth)
y     <- mvrnorm(1, rep(5,n+np), Sigma)

X <- matrix(1, nrow=length(y), ncol=1)
x1 <- round(rnorm(n+np),3) # noise predictor

data.fit  <- data.frame(y=y[1:n], x1=x1[1:n])
data.pred <- data.frame(y=y[n+1:np], x1=x1[n+1:np])

S.fit  <- S[1:n,]
S.pred <- S[n+1:np,]

# get initial values
inits <- spacious:::initial.theta(data.fit$y, S.fit)
cov.inits <- list(nugget=inits[1], psill=inits[2], range=inits[3], smoothness=smooth)

source("custom.R")

test_that("full matern converges", {
	fit <- spacious(y~x1, S=S.fit, cov="matern", fixed=list(smoothness=smooth), blocks=list(type="regular", nblocks=0), data=data.fit, cov.inits=cov.inits)
	expect_equal(fit$conv, TRUE)
})

test_that("block matern converges", {
	fit <- spacious(y~x1, S=S.fit, cov="matern", fixed=list(smoothness=smooth), blocks=list(type="regular", nblocks=2^2), data=data.fit, cov.inits=cov.inits)
	expect_equal(fit$conv, TRUE)
})

test_that("matern threads equals non-threads", {
	fit.N <- spacious(y~x1, S=S.fit, cov="matern", fixed=list(smoothness=smooth), blocks=list(type="regular", nblocks=2^2), nthreads=1, data=data.fit, cov.inits=cov.inits)
	fit.T <- spacious(y~x1, S=S.fit, cov="matern", fixed=list(smoothness=smooth), blocks=list(type="regular", nblocks=2^2), nthreads=2, data=data.fit, cov.inits=cov.inits)
	expect_equal(fit.N$ll, fit.T$ll)
})

test_that("matern gpu equals non-gpu", {
	fit.N <- spacious(y~x1, S=S.fit, cov="matern", fixed=list(smoothness=smooth), blocks=list(type="full"), gpu=FALSE, data=data.fit, cov.inits=cov.inits)
	fit.G <- spacious(y~x1, S=S.fit, cov="matern", fixed=list(smoothness=smooth), blocks=list(type="full"), gpu=TRUE, data=data.fit, cov.inits=cov.inits)
	expect_equal(fit.N$ll, fit.G$ll)
})

test_that("matern full engine is the same in R and C++", {
	fit.R <- spacious(y~x1, S=S.fit, cov="matern", fixed=list(smoothness=smooth), blocks=list(type="full"), engine="R", data=data.fit, cov.inits=cov.inits)
	fit.C <- spacious(y~x1, S=S.fit, cov="matern", fixed=list(smoothness=smooth), blocks=list(type="full"), engine="C", data=data.fit, cov.inits=cov.inits)
	expect_equal(fit.R$ll, fit.C$ll)

	preds.R <- predict(fit.C, newdata=data.pred, newS=S.pred, interval="prediction", level=0.9, engine="R")
	preds.C <- predict(fit.C, newdata=data.pred, newS=S.pred, interval="prediction", level=0.9, engine="C")

	for (i in 1:np) {
		expect_equal(preds.R$y[i], preds.C$y[i], label=paste0("y for ",i))
		expect_equal(preds.R$sd[i], preds.C$sd[i], label=paste0("sd for ",i))
	}
})

test_that("matern blocks engine is the same in R and C++", {
	fit.R <- spacious(y~x1, S=S.fit, cov="matern", fixed=list(smoothness=smooth), blocks=list(type="regular", nblocks=2^2), engine="R", data=data.fit, cov.inits=cov.inits)
	fit.C <- spacious(y~x1, S=S.fit, cov="matern", fixed=list(smoothness=smooth), blocks=list(type="regular", nblocks=2^2), engine="C", data=data.fit, cov.inits=cov.inits)
	expect_equal(fit.R$ll, fit.C$ll)

	preds.R <- predict(fit.C, newdata=data.pred, newS=S.pred, interval="prediction", level=0.9, engine="R")
	preds.C <- predict(fit.C, newdata=data.pred, newS=S.pred, interval="prediction", level=0.9, engine="C")

	for (i in 1:np) {
		expect_equal(preds.R$y[i], preds.C$y[i], label=paste0("y for ",i))
		expect_equal(preds.R$sd[i], preds.C$sd[i], label=paste0("sd for ",i))
	}
})

test_that("matern with smoothness=0.5 gives same as exponential", {
	fit.E <- spacious(y~x1, S=S.fit, cov="exp", blocks=list(type="regular", nblocks=2^2), engine="R", data=data.fit, cov.inits=cov.inits)
	fit.M <- spacious(y~x1, S=S.fit, cov="matern", fixed=list(smoothness=0.5), blocks=list(type="regular", nblocks=2^2), engine="C", data=data.fit, cov.inits=cov.inits)
	expect_equal(fit.E$ll, fit.M$ll)
})
