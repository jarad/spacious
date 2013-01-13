# print the summary
"print.spacious" <- function(x, ...) {
	print(summary(x))
	cat("Convergence:", x$conv,"\n")
	cat("Iterations:", x$nIter,"\n")
	cat("-2 x Log Likelihood:", x$ll,"\n")
}
