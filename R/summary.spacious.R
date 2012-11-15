# summarizes fit of block composite model from spacious
"summary.spacious" <- function(object, ...) {
# TODO: figure out how to make nice tables of results

	p.beta <- sapply(object$beta/object$se.beta,
		function(z) {
			2*( 1-pnorm(abs(z)) )
		})

	df.beta <- data.frame(beta=object$beta, se=object$se.beta, p=p.beta)
	df.theta <- data.frame(theta=object$theta, se=object$se.theta)

	colnames(df.beta) <- c("Estimate","Std Err","P-value")
	rownames(df.beta) <- paste0("b",0:(nrow(df.beta)-1) )

	colnames(df.theta) <- c("Estimate","Std Err")
	rownames(df.theta) <- c("Nugget","Partial Sill","Range")

	ss <- list(
		df.beta=df.beta,
		df.theta=df.theta
	)
	class(ss) <- "summary.spacious"

	ss
}

# print the summary
"print.summary.spacious" <- function(x, ...) {
	cat("Coefficients:\n")
	print(x$df.beta)

	cat("Spatial parameters\n")
	print(x$df.theta)
}
