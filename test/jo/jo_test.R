# compare to Jo's matlab results
if ("package:spacious" %in% search()) {
	# unload the package
	detach("package:spacious", unload=TRUE)
}

require(spacious)

S <- cbind(read.table("test/jo/synt_data_pos1.txt")$V1,read.table("test/jo/synt_data_pos2.txt")$V1)
covar <- S[,1]
y <- read.table("test/jo/synt_data.txt")$V1

#fit.FL <- spacious(y~covar, S=S, blocks=list(type="regular", nblocks=1^2))
fit.CL <- spacious(y~covar, S=S, blocks=list(type="regular", nblocks=5^2), verbose=TRUE)
cat("Full:\n"); print(summary(fit.FL))
cat("Comp:\n"); print(summary(fit.CL))

# predict at new sites
newS <- rbind(
	c(0.43,0.18),
	c(0.23,0.37),
	c(0.49,0.51),
	c(0.30,0.70)
)

#preds.FL <- predict(fit.FL, newS=newS, newdata=data.frame(covar=newS[,1]), interval="prediction", opts=list(type="all"))
#preds.CL <- predict(fit.CL, newS=newS, newdata=data.frame(covar=newS[,1]), interval="prediction")
#print(preds.FL)
#print(preds.CL)
