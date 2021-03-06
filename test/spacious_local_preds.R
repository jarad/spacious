# make local predictions

if ("package:spacious" %in% search()) {
  # unload the package
  detach("package:spacious", unload=TRUE)
}

require(spacious)

set.seed(311)

nrow<-20
range<-0.1
nug<-0.1
ps<-1

S<-as.matrix(expand.grid(1:nrow/nrow,1:nrow/nrow))
n<-nrow(S)
d<-rdist(S)
C<-exp(-d/range)
#Y<-10+S[,1]+sqrt(ps)*t(chol(C))%*%rnorm(n)+rnorm(n,0,sqrt(nug))
Sigma <- nug * diag(n) + ps * exp(-d/range)
Y<-mvrnorm(1, mu=rep(10, n), Sigma=Sigma)

# fit model
cat("Fitting model...\n")
X<-S
#fit1 <- spacious(Y~X, S=S, cov="exp", blocks=list(type="regular", nblocks=10^2), verbose=TRUE)
#fit1 <- spacious(Y~1, S=S, cov="exp", blocks=list(type="regular", nblocks=10^2), verbose=TRUE)

nrowp<-nrow
Sp<-as.matrix(expand.grid(1:nrowp/nrowp,1:nrowp/nrowp))
Xp<-Sp

# make predictions
cat("Making predictions...\n")
#pred1 <- predict(fit1, newdata=Xp, newS=Sp)
pred1 <- predict(fit1, newS=Sp[1:10,], opts=list(type="local", num=100), interval="prediction")
print(pred1)
