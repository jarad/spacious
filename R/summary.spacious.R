# summarizes fit of block composite model from spacious
"summary.spacious" <- function(object, ...) {
# TODO: figure out how to make nice tables of results
	cat("Spacious fit:\n")

	coefs <- data.frame(beta=object$beta, se=object$se.beta)
	covp <- data.frame(theta=object$theta, se=object$se.theta)

	print(coefs)
	print(covp)

	NULL
}
