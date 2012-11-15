# plot dianostics of model fit
"plot.spacious" <- function(object) {

	# plot spatial parameter estimates vs iteration
	par(mfrow=c(2,2))
	plot(0:object$nIter, object$iters.theta[,1], type="b", xlab="Iteration", ylab="Nugget")
	plot(0:object$nIter, object$iters.theta[,2], type="b", xlab="Iteration", ylab="Partial Sill")
	plot(0:object$nIter, object$iters.theta[,3], type="b", xlab="Iteration", ylab="Range")
	plot(0:object$nIter, object$iters.ll, type="b", xlab="Iteration", ylab="-2 * log likelihood")

}
