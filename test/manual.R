# run code for generating user manual
options(digits=2)

if ("package:spacious" %in% search()) {
  # unload the package
  detach("package:spacious", unload=TRUE)
	library.dynam.unload("spacious", libpath="~/Rlib/spacious")
}

# load the package
library(spacious)
data(anom.2011)

library(ggmap)
if (!("euro" %in% search())) {
	euro <- get_map(location=c(mean(anom.2011$lon),mean(anom.2011$lat)), zoom=4, maptype="satellite")
}

S <- cbind(anom.2011$lon, anom.2011$lat)

if (FALSE) {
# plot data

fig  <- ggmap(euro) +
        geom_point(aes(lon, lat, color=anom), shape=15, size=1.5, data=anom.2011) +
        scale_colour_gradient(low="green", high="red")
pdf("pdf/manual_plot.pdf")
	print(fig)
graphics.off()
}

if (FALSE) {
# basic call
fit <- spacious(anom ~ lon + lat + elev, S=S, data=anom.2011)
pdf("pdf/manual_fit.pdf")
	plot(fit)
graphics.off()
print(summary(fit))
}

if (FALSE) {
# fix range=5
fit <- spacious(anom ~ lon + lat + elev, S=S, data=anom.2011, fixed=list(range=5))
print(summary(fit))
}

if (FALSE) {
# matern covariance with smoothness=1
fit <- spacious(anom ~ lon + lat + elev, S=S, data=anom.2011, cov="matern", fixed=list(smoothness=1))
print(summary(fit))
}

if (TRUE) {
# blocking structures
library(maps)
#fit_c <- spacious(anom ~ lon + lat + elev, S=S, data=anom.2011, blocks=list(type="cluster"), gpu=TRUE)
#fit_r <- spacious(anom ~ lon + lat + elev, S=S, data=anom.2011, blocks=list(type="regular"), gpu=TRUE)

pdf("pdf/manual_blocks_c.pdf") #, height=7/1.5, width=7/1.5)
	plot(S, pch=4, xlab="lon", ylab="lat", main="Cluster blocks", cex=0.1)
	map("world", add=TRUE, col="darkgreen")
	plot(fit_c$grid,lty=1,lwd=4,border="blue",cex=0.25,add=TRUE)
graphics.off()

pdf("pdf/manual_blocks_r.pdf") #, height=7/1.5, width=7/1.5)
	plot(S, pch=4, xlab="lon", ylab="lat", main="Regular blocks", cex=0.1)
	map("world", add=TRUE, col="darkgreen")
	plot(fit_r$grid,lty=1,lwd=4,border="blue",cex=0.25,add=TRUE)
graphics.off()
}

if (FALSE) {
# predictions
#fit     <- spacious(anom~lon + lat + elev, S=S, data=anom.2011)
S.p     <- cbind(anom.pred.grid$lon, anom.pred.grid$lat)
preds   <- predict(fit, newdata=anom.pred.grid, newS=S.p, interval="prediction")
pred    <- preds$y
pred.sd <- preds$sd

fig <- ggmap(euro) + ggtitle("Predictions") +
       geom_point(aes(lon,lat,color=pred), shape=15, size=1.25, data=data.frame(anom.pred.grid, pred=pred))+
       scale_colour_gradient(low="green", high="red")
pdf("pdf/manual_pred.pdf")
	print(fig)
graphics.off()

fig <- ggmap(euro) + ggtitle("Prediction SD") +
       geom_point(aes(lon,lat,color=pred.sd), shape=15, size=1.25, data=data.frame(anom.pred.grid, pred.sd=pred.sd))+
       scale_colour_gradient(low="green", high="red")
pdf("pdf/manual_pred_sd.pdf")
	print(fig)
graphics.off()
}
