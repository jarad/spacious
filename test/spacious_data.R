# test out spacious

if ("package:spacious" %in% search()) {
	# unload the package
	detach("package:spacious", unload=TRUE)
}

# load the package
require(spacious)

data(mean_max_temps)

fit.full <- spacious(temp~1, data=mean_max_temps, S=cbind(mean_max_temps$lat,mean_max_temps$long), blocks=NULL, verbose=TRUE)
fit.clus <- spacious(temp~1, data=mean_max_temps, S=cbind(mean_max_temps$lat,mean_max_temps$long), blocks=list(type="cluster", nblocks=25), verbose=TRUE)
#fit <- spacious(temp~1, data=mean_max_temps, S=cbind(mean_max_temps$lat,mean_max_temps$long), nblocks=9, verbose=TRUE)

print(summary(fit.full))
print(summary(fit.clus))

