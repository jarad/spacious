# test out spacious

if ("package:spacious" %in% search()) {
	# unload the package
	detach("package:spacious", unload=TRUE)
}

# load the package
require(spacious)

data(mean_max_temps)

#fit.full <- spacious(temp~1, data=mean_max_temps, S=cbind(mean_max_temps$lat,mean_max_temps$long), blocks=NULL, verbose=TRUE)
#fit.re <- spacious(temp~1, data=mean_max_temps, S=cbind(mean_max_temps$lat,mean_max_temps$long), blocks=list(type="regular", nblocks=100), verbose=TRUE)
#fit <- spacious(temp~1, data=mean_max_temps, S=cbind(mean_max_temps$lat,mean_max_temps$long), nblocks=9, verbose=TRUE)

#print(summary(fit.full))
#print(summary(fit.clus))

# fit model
fit.clus <- spacious(temp~1, data=mean_max_temps, S=cbind(mean_max_temps$lat,mean_max_temps$long), blocks=list(type="cluster", nblocks=250), verbose=TRUE)

if (FALSE) {
# plot data
pdf("pdf/mmt.pdf")
	quilt.plot(mean_max_temps$lat, mean_max_temps$lon, mean_max_temps$temp,
		nx=250, ny=250, main="Observation Locations", add.legend=TRUE)
	world(add=TRUE, lwd=1, col="grey")
graphics.off()

pdf("pdf/mmt_grid.pdf")
	quilt.plot(mean_max_temps$lat, mean_max_temps$lon, mean_max_temps$temp,
		nx=250, ny=250, main="Observation Locations", add.legend=TRUE)
	world(add=TRUE, lwd=1, col="grey")
	plot(fit.clus$grid,lty=2,add=TRUE)
graphics.off()

}

if (FALSE) {
# make predictions

newS <- as.matrix(expand.grid(seq(25,30,length=25), seq(47,52,length=25)))
#preds <- predict(fit.clus, newdata=data.frame( rep(1, nrow(newS)) ), newS=newS, newB=rep(5, nrow(newS)))

preds <- rep(NA, nrow(newS))

sapply(1:nrow(newS), function(i) {
	try({
		preds[i] <<- predict(fit.clus, newdata=data.frame( rep(1, 1) ), newS=newS[i,])$y
	})
})

print(preds)

}

if (FALSE) { 
# plot data + predictions
pdf("pdf/mmt_pred.pdf")
	quilt.plot(c(mean_max_temps$lat,newS[,1]), c(mean_max_temps$lon,newS[,2]), c(mean_max_temps$temp,preds),
		nx=250, ny=250, main="Observations + Predictions", add.legend=TRUE)
	world(add=TRUE, lwd=1, col="grey")
graphics.off()
}
