# compare to Jo's matlab results
if ("package:spacious" %in% search()) {
	# unload the package
	detach("package:spacious", unload=TRUE)
	library.dynam.unload("spacious", libpath="~/Rlib/spacious")
}

require(spacious)

S <- cbind(read.table("test/jo/synt_data_pos1.txt")$V1,read.table("test/jo/synt_data_pos2.txt")$V1)
covar <- S[,1]
y <- read.table("test/jo/synt_data.txt")$V1

#fit.FL_R <- spacious(y~covar, S=S, blocks=list(type="full"), engine="R")
#fit.FL_C <- spacious(y~covar, S=S, blocks=list(type="full"), engine="C")
#cat("Full R:\n"); print(summary(fit.FL_R))
#cat("Full C:\n"); print(summary(fit.FL_C))
fit.CL_R <- spacious(y~covar, S=S, blocks=list(type="regular", nblocks=5^2), engine="R", verbose=TRUE)
fit.CL_C1 <- spacious(y~covar, S=S, blocks=list(type="regular", nblocks=5^2), engine="C", verbose=TRUE)
fit.CL_C4 <- spacious(y~covar, S=S, blocks=list(type="regular", nblocks=5^2), engine="C", nthreads=4, verbose=TRUE)
cat("Comp R:\n"); print(summary(fit.CL_R))
cat("Comp C1:\n"); print(summary(fit.CL_C1))
cat("Comp C4:\n"); print(summary(fit.CL_C4))
done

# predict at new sites
newS <- rbind(
	c(0.43,0.18),
	c(0.23,0.37),
	c(0.49,0.51),
	c(0.30,0.70)
)

preds.FL <- predict(fit.FL, newS=newS, newdata=data.frame(covar=newS[,1]), interval="prediction", opts=list(type="all"))
preds.CL <- predict(fit.CL, newS=newS, newdata=data.frame(covar=newS[,1]), interval="prediction")
print(preds.FL)
print(preds.CL)
