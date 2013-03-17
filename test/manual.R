# run code for generating user manual

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

if (TRUE) {
# basic call
S   <- cbind(anom.2011$lon, anom.2011$lat)
fit <- spacious(anom ~ lon + lat + elev, S=S, data=anom.2011, verbose=TRUE)
pdf("pdf/manual_fit.pdf")
	plot(fit)
graphics.off()
}
