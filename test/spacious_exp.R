# test out spacious

if ("package:spacious" %in% search()) {
	# unload the package
	detach("package:spacious", unload=TRUE)
	library.dynam.unload("spacious", libpath="~/Rlib/spacious")
}

# load the package
require(spacious)

if (TRUE) {
set.seed(311)

# generate data to use for fitting a block composite model
n <- 1*1024
np <- 100   # number to predict

# generate spatial locations S
S <- round(cbind(runif(n+np), runif(n+np)),3)

# distance matrix
D <- rdist(S); D[row(D)==col(D)] <- 0

# matern covariance function parameters
nugget <- 0.25
tau2   <- 0.25
range  <- 0.25
smooth <- 0.50  # 0.5 = exponential cov
mu <- 5
#Sigma <- nugget * diag(n) + tau2 * matern(D, range, smooth)
Sigma <- nugget * diag(n+np) + tau2 * exp(-D/range)

# generate data
y <- round(mvrnorm(1, mu=rep(mu, n+np), Sigma=Sigma),3)

# fit with spacious
X <- matrix(1, nrow=length(y), ncol=1)
x1 <- round(rnorm(n+np),3)

y.fit <- y[1:n]
y.pred <- y[n+1:np]
X.fit <- matrix(x1[1:n], nrow=n, ncol=1)
X.pred <- matrix(x1[n+1:np], nrow=np, ncol=1)
S.fit <- S[1:n,]
S.pred <- S[n+1:np,]
}

if (FALSE) { # test with lm()
	fit <- lm(y.fit~X.fit)
	preds <- predict(fit, newdata=as.data.frame(X.pred))
print(fit)
print(preds)
done
}

if (FALSE) {
# try this 500 times
for (i in 1:500) {
	fit.spacious <- spacious(y~x, data=data.frame(y=y[1:n], x=x1[1:n]), S=S[1:n,], cov="exp", blocks=list(type="full"), nthreads=1, gpu=TRUE, verbose=TRUE)
	print(fit.spacious$time)
}
done
}

if (FALSE) { # GPU vs not
fit.spacious <- spacious(y~x, data=data.frame(y=y[1:n], x=x1[1:n]), S=S[1:n,], cov="exp", blocks=list(type="full"), nthreads=1, gpu=TRUE, verbose=TRUE)
print(fit.spacious$time)

fit.spacious <- spacious(y~x, data=data.frame(y=y[1:n], x=x1[1:n]), S=S[1:n,], cov="exp", blocks=list(type="full"), nthreads=1, gpu=FALSE, verbose=TRUE)
print(fit.spacious$time)
done
}

if (FALSE) { # threads
cat("Fitting...\n")
fit.spacious <- spacious(y~x, data=data.frame(y=y[1:n], x=x1[1:n]), S=S[1:n,], cov="exp", blocks=list(type="regular", nblocks=6^2), nthreads=1, gpu=FALSE, verbose=TRUE)
print(fit.spacious$time)
fit.spacious <- spacious(y~x, data=data.frame(y=y[1:n], x=x1[1:n]), S=S[1:n,], cov="exp", blocks=list(type="regular", nblocks=6^2), nthreads=4, gpu=FALSE, verbose=TRUE)
print(fit.spacious$time)
#fit.spacious <- spacious(y~x, data=data.frame(y=y[1:n], x=x1[1:n]), S=S[1:n,], cov="exp", blocks=list(type="regular", nblocks=6^2), nthreads=3, gpu=FALSE)
#print(fit.spacious$time)
#fit.spacious <- spacious(y~x, data=data.frame(y=y[1:n], x=x1[1:n]), S=S[1:n,], cov="exp", blocks=list(type="regular", nblocks=6^2), nthreads=2, gpu=FALSE)
#print(fit.spacious$time)
#fit.spacious <- spacious(y~x, data=data.frame(y=y[1:n], x=x1[1:n]), S=S[1:n,], cov="exp", blocks=list(type="regular", nblocks=6^2), nthreads=1, gpu=FALSE)
#print(fit.spacious$time)
}

if (TRUE) { # compare R vs C
	#fit.R <- spacious(y~x2, data=data.frame(y=y[1:n], x2=x1[1:n]), S=S[1:n,], cov="exp", blocks=list(type="full"), verbose=TRUE, engine="R")
	#fit.C <- spacious(y~x2, data=data.frame(y=y[1:n], x2=x1[1:n]), S=S[1:n,], cov="exp", blocks=list(type="full"), verbose=TRUE, engine="C")
	#fit.CG <- spacious(y~x2, data=data.frame(y=y[1:n], x2=x1[1:n]), S=S[1:n,], cov="exp", blocks=list(type="full"), verbose=TRUE, engine="C", gpu=TRUE)
#preds.C <- predict(fit.CG, newdata=data.frame(x2=x1[(n+1):(n+np)]), newS=S[(n+1):(n+np),], interval="prediction", level=0.9, engine="C")
#preds.R <- predict(fit.CG, newdata=data.frame(x2=x1[(n+1):(n+np)]), newS=S[(n+1):(n+np),], interval="prediction", level=0.9, engine="R")
#print(preds.C)
#print(preds.R)
#done
	#fit.R <- spacious(y~x2, data=data.frame(y=y[1:n], x2=x1[1:n]), S=S[1:n,], cov="exp", blocks=list(type="regular", nblocks=6^2), verbose=TRUE, engine="R")
#preds.R <- predict(fit.R, newdata=data.frame(x2=x1[(n+1):(n+np)]), newS=S[(n+1):(n+np),], interval="prediction", level=0.9)
	fit.C <- spacious(y~x2, data=data.frame(y=y[1:n], x2=x1[1:n]), S=S[1:n,], cov="exp", blocks=list(type="regular", nblocks=6^2), verbose=TRUE, nthreads=4, engine="C")
t1 <- proc.time()
preds.R <- predict(fit.C, newdata=data.frame(x2=x1[(n+1):(n+np)]), newS=S[(n+1):(n+np),], interval="prediction", level=0.9, engine="R")
print(proc.time()-t1)
t1 <- proc.time()
preds.C <- predict(fit.C, newdata=data.frame(x2=x1[(n+1):(n+np)]), newS=S[(n+1):(n+np),], interval="prediction", level=0.9, engine="C")
print(proc.time()-t1)
print(head(preds.R))
print(head(preds.C))
done
	fit.CT <- spacious(y~x2, data=data.frame(y=y[1:n], x2=x1[1:n]), S=S[1:n,], cov="exp", blocks=list(type="regular", nblocks=6^2), verbose=TRUE, nthreads=4, engine="C")
done
	#fit.R <- spacious(y~x2, data=data.frame(y=y[1:n], x2=x1[1:n]), S=S[1:n,], cov="matern", fixed=list(smoothness=1/2), blocks=list(type="regular", nblocks=6^2), verbose=TRUE, engine="R")
	#fit.R <- spacious(y~x2, data=data.frame(y=y[1:n], x2=x1[1:n]), S=S[1:n,], cov="matern", fixed=list(smoothness=3/2), blocks=list(type="regular", nblocks=6^2), verbose=TRUE, engine="R")
	#fit.R <- spacious(y~x2, data=data.frame(y=y[1:n], x2=x1[1:n]), S=S[1:n,], cov="matern", fixed=list(smoothness=5/2), blocks=list(type="regular", nblocks=6^2), verbose=TRUE, engine="R")
	#fit.R <- spacious(y~x2, data=data.frame(y=y[1:n], x2=x1[1:n]), S=S[1:n,], cov="matern", fixed=list(smoothness=7/2), blocks=list(type="regular", nblocks=6^2), verbose=TRUE, engine="R")
	#fit.R <- spacious(y~x2, data=data.frame(y=y[1:n], x2=x1[1:n]), S=S[1:n,], cov="matern", fixed=list(smoothness=9/2), blocks=list(type="regular", nblocks=6^2), verbose=TRUE, engine="R")
	#fit.R <- spacious(y~x2, data=data.frame(y=y[1:n], x2=x1[1:n]), S=S[1:n,], cov="matern", fixed=list(smoothness=11/2), blocks=list(type="regular", nblocks=6^2), verbose=TRUE, engine="R")
	#fit.R <- spacious(y~x2, data=data.frame(y=y[1:n], x2=x1[1:n]), S=S[1:n,], cov="matern", fixed=list(smoothness=13/2), blocks=list(type="regular", nblocks=6^2), verbose=TRUE, engine="R")
#done
	#fit.CE <- spacious(y~x2, data=data.frame(y=y[1:n], x2=x1[1:n]), S=S[1:n,], cov="exp", blocks=list(type="regular", nblocks=6^2), verbose=TRUE, engine="C", nthreads=4)
	#fit.CM <- spacious(y~x2, data=data.frame(y=y[1:n], x2=x1[1:n]), S=S[1:n,], cov="matern", blocks=list(type="regular", nblocks=6^2), verbose=TRUE, engine="C", nthreads=4)
	#fit.RE <- spacious(y~x2, data=data.frame(y=y[1:n], x2=x1[1:n]), S=S[1:n,], cov="exp", blocks=list(type="regular", nblocks=6^2), verbose=TRUE, engine="R")
	#fit.RM <- spacious(y~x2, data=data.frame(y=y[1:n], x2=x1[1:n]), S=S[1:n,], cov="matern", blocks=list(type="regular", nblocks=6^2), verbose=TRUE, engine="R")
done
	preds.R <- predict(fit.R, newdata=data.frame(x2=x1[(n+1):(n+np)]), newS=S[(n+1):(n+np),], interval="prediction", level=0.9)
	preds.C <- predict(fit.C, newdata=data.frame(x2=x1[(n+1):(n+np)]), newS=S[(n+1):(n+np),], interval="prediction", level=0.9)
}

done

time.spacious <- proc.time()
#fit.spacious <- spacious(y, X, S, cov="exp", nblocks=1^2)
#fit.spacious <- spacious(y~x2, data=data.frame(y=y[1:n], x2=x1[1:n]), S=S[1:n,], cov="exp", nblocks=2^2, verbose=TRUE)
fit.spacious <- spacious(y~x2, data=data.frame(y=y[1:n], x2=x1[1:n]), S=S[1:n,], cov="exp", blocks=list(type="regular", nblocks=6^2), verbose=TRUE, nthreads=4, gpu=FALSE)
#fit.spacious <- spacious(y.fit~X.fit, S=S.fit, cov="exp", blocks=list(type="regular", nblocks=2^2), verbose=TRUE)
#fit.spacious <- spacious(y~x, data=data.frame(y=y[1:n], x=x1[1:n]), S=S[1:n,], cov="exp", blocks=list(type="full"), verbose=TRUE)
time.spacious <- proc.time() - time.spacious
beta.spacious <- fit.spacious$beta
theta.spacious <- fit.spacious$theta

print(summary(fit.spacious))

cat("Spacious estimates:",beta.spacious,theta.spacious,"\n")
cat("Spacious SEs:",fit.spacious$se.beta,fit.spacious$se.theta,"\n")
cat("Spacious execution time:\n")
print(time.spacious)

preds <- predict(fit.spacious, newdata=data.frame(x2=x1[(n+1):(n+np)]), newS=S[(n+1):(n+np),], interval="prediction", level=0.9)
#preds <- predict(fit.spacious, newS=S[(n+1):(n+np),], interval="prediction", level=0.9)
#preds <- predict(fit.spacious, newdata=X.pred, newS=S.pred, interval="prediction", level=0.9)
cat("Predictions:\n")
print(preds)
cat("Actual:\n")
print(y.pred)

done


# predict at a far out location
preds <- predict(fit.spacious, newdata=data.frame(x2=0), newS=matrix(rbind( c(10,10), c(-5, -5) ),nrow=2,ncol=2), interval="prediction", level=0.9)
print(preds)

# consider using the rbenchmark package

if (0) {
# try likfit

require(geoR) # for likfit

gd <- as.geodata(cbind(S[1:n,],y[1:n]))
time.likfit <- proc.time()
fit.likfit <- likfit(gd, ini.cov.pars=c(1.22,1.22))
time.likfit <- proc.time() - time.likfit

beta.likfit <- fit.likfit$beta
theta.likfit <- c(fit.likfit$tausq, fit.likfit$sigmasq, 1/fit.likfit$phi)

cat("likfit estimates:",beta.likfit,theta.likfit,"\n")
cat("likfit execution time:\n")
print(time.likfit)
}
