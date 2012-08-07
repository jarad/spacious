# summarizes fit of block composite model from spacious
"summary.spacious" <- function(object, ...) {
	cat("Spacious summary!\n")
	cat("Coefficients:",object$beta,"\n")
	cat("Covariance parameters:",object$theta,"\n")
}
