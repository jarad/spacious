# fit models to generated data
require(geoR)

# S holds the spatial locations, D is the distance matrix, Sigma is the covariance matrix,
# Y is the simulated data, and gd is a geodata object of S and Y

# fit with geoR's likfit
#likfit.ml <- likfit(gd, lik.method="REML", cov.model="matern", ini=c(.05,.05), fix.kappa=FALSE, kappa=0.2)
#cat("Nugget:", likfit.ml$nugget,"\n")
#cat("Covariance params:", likfit.ml$cov.pars, "\n")
#cat("Smoothness:", likfit.ml$kappa,"\n")
#cat("Mean:", likfit.ml$beta,"\n")

#likfit.reml <- likfit(gd, ini=c(1,0.5), lik.method="REML")
#print(likfit.reml)

"loglik.matern" <- function(beta, lnugget, lsill, lrange, lsmoothness,   # params
	X, Y, D   # data
) {
	nugget <- exp(lnugget)
	sill <- exp(lsill)
	range <- exp(lrange)
	smoothness <- exp(lsmoothness)

	cholSigma <- chol(Sigma<-nugget * diag(n) + sill * matern(D, range, 0.5))

	loglik <- -sum(log(diag(cholSigma))) -0.5*crossprod(backsolve(cholSigma, Y-X %*% beta))
	-loglik
}

"mf" <- function(params, X, Y, D) {
	loglik.matern(params[1], params[2], params[3], params[4], params[5], X, Y, D)
}

X <- matrix(1, nrow=length(Y), ncol=1)
beta <- matrix(0, nrow=1, ncol=1)

print(loglik.matern(mu, nugget, tau2, range, smooth, X, Y, D))
print(loglik.matern(mu+0.5, nugget, tau2, range, smooth, X, Y, D))
print(loglik.matern(mu-0.5, nugget, tau2, range, smooth, X, Y, D))

#res <- optim( par=c(mu, nugget, tau2, range, smooth), fn=mf, X=X, Y=Y, D=D, method="BFGS" )
res <- nlm( f=mf, p=c(mu, nugget, tau2, range, smooth), X=X, Y=Y, D=D )
print(res)
