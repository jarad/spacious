# plot dianostics of model fit
"plot.spacious" <- function(x, ...) {

	# plot spatial parameter estimates vs iteration
	par(mfrow=c(2,2))
	plot(0:x$nIter, x$iters_theta[,1], type="b", xlab="Iteration", ylab="Nugget")
	plot(0:x$nIter, x$iters_theta[,2], type="b", xlab="Iteration", ylab="Partial Sill")
	plot(0:x$nIter, x$iters_theta[,3], type="b", xlab="Iteration", ylab="Range")
	plot(0:x$nIter, x$iters_ll, type="b", xlab="Iteration", ylab="Log-Likelihood")

}
