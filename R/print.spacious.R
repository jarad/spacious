# print the summary
"print.spacious" <- function(object, ...) {
	print(summary(object))
	cat("Convergence:", object$conv,"\n")
	cat("Iterations:", object$nIter,"\n")
	cat("-2 x Log Likelihood:", object$ll,"\n")
}
