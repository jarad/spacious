# test data set

if ("package:spacious" %in% search()) {
  # unload the package
  detach("package:spacious", unload=TRUE)
}

require(spacious)
require(testthat)

# test fitting to mean_max_temps data
test_that("mean max temps fit", {
	# cluster blocks
	fit <- spacious(temp~elev, data=mean_max_temps, S=cbind(mean_max_temps$lat,mean_max_temps$long), blocks=list(type="cluster", nblocks=250))
	summary(fit)
	p <- predict(fit, newdata=data.frame(elev=270), newS=cbind(10,50), interval="prediction")
	print(p)
	p <- predict(fit, newdata=data.frame(elev=270), newS=cbind(10,50), interval="prediction", opts=list(type="local", num=5))
	print(p)
})
