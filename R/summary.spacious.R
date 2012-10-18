# summarizes fit of block composite model from spacious
"summary.spacious" <- function(object, ...) {
# TODO: figure out how to make nice tables of results

	df.beta <- data.frame(beta=object$beta, se=object$se.beta)
	df.theta <- data.frame(theta=object$theta, se=object$se.theta)

	ss <- list(
		df.beta=df.beta,
		df.theta=df.theta
	)
	class(ss) <- "summary.spacious"

	ss
}

# print the summary
"print.summary.spacious" <- function(x, ...) {
	cat("Spacious fit:\n")
	print(x$df.beta)
	print(x$df.theta)
}
