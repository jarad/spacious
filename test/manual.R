# run code for generating user manual
options(digits=2)

if ("package:spacious" %in% search()) {
  # unload the package
  detach("package:spacious", unload=TRUE)
}

# load the package
library(spacious)
data(anom.2011)

if (FALSE) {
# plot data
library(ggmap)

euro <- get_map(location=c(mean(anom.2011$lon),mean(anom.2011$lat)), zoom=4, maptype="satellite")
fig  <- ggmap(euro) +
        geom_point(aes(lon, lat, color=anom), shape=15, size=1.5, data=anom.2011) +
        scale_colour_gradient(low="green", high="red")
pdf("pdf/manual_plot.pdf")
	print(fig)
graphics.off()
}

if (FALSE) {
# basic call
S   <- cbind(anom.2011$lon, anom.2011$lat)
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

if (FALSE) {
# blocking structures
library(maps)
fit_c <- spacious(anom ~ lon + lat + elev, S=S, data=anom.2011, blocks=list(type="cluster"))
fit_r <- spacious(anom ~ lon + lat + elev, S=S, data=anom.2011, blocks=list(type="regular"))

pdf("pdf/manual_blocks_c.pdf")
	plot(S, col=fit_c$B, pch=15, xlab="lon", ylab="lat", main="Cluster blocks")
	map("world", add=TRUE)
	plot(fit_c$grid,lty=1,border="gray",cex=0.25,add=TRUE)
graphics.off()

pdf("pdf/manual_blocks_r.pdf")
	plot(S, col=fit_r$B, pch=15, xlab="lon", ylab="lat", main="Regular blocks")
	map("world", add=TRUE)
	plot(fit_r$grid,lty=1,border="gray",cex=0.25,add=TRUE)
graphics.off()
}

if (FALSE) {
# predictions
fit <- spacious(anom~lon + lat + elev, S=S, data=anom.2011)
S.p <- cbind(anom.pred.grid$lon, anom.pred.grid$lat)
pred <- predict(fit, newdata=anom.pred.grid, newS=S.p)$y

fig <- ggmap(euro) +
       geom_point(aes(lon,lat,color=pred), shape=15, size=1.25, data=data.frame(anom.pred.grid, pred=pred))+
       scale_colour_gradient(low="green", high="red")
pdf("pdf/manual_pred.pdf")
	print(fig)
graphics.off()
}
