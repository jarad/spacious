if ("package:spacious" %in% search()) {
	# unload the package
	detach("package:spacious", unload=TRUE)
	library.dynam.unload("spacious", libpath="~/Rlib/spacious")
}

require(spacious)

#Testing out the spatial prediction functions:

#prediction grid
sss<-seq(0,1,.05)
sp<-expand.grid(sss,sss)
sp<-as.matrix(sp)

#Data locations
s<-cbind(runif(500),runif(500))

#Generate data
sall<-rbind(s,sp)
d<-rdist(sall)
diag(d)<-0
n<-nrow(d)
C<-0.1*diag(n)+1.2*exp(-d/.1)
yall<-t(chol(C))%*%rnorm(n)
y<-yall[1:500]
yp<-yall[-c(1:500)]


#Fit model with different block sizes
fit0 <- spacious(y~1, S=s,blocks=list(type="cluster", nblocks=1^2))
fit1 <- spacious(y~1, S=s,blocks=list(type="cluster", nblocks=5^2))
fit2 <- spacious(y~1, S=s,blocks=list(type="cluster", nblocks=7^2))

