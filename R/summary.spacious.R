# summarizes fit of block composite model from spacious
"summary.spacious" <- function(object, ...) {
	p_beta <- rep(NA, length(object$beta))
	if (!is.null(object$se_beta)) {
		p_beta <- sapply(object$beta/object$se_beta,
			function(z) {
				2*( 1-pnorm(abs(z)) )
			})
	}

	df.beta <- data.frame(beta=object$beta, se=object$se_beta, p=p_beta)
	df.theta <- data.frame(theta=object$theta, se=object$se_theta)

	colnames(df.beta) <- c("Estimate","Std Err","P-value")
	rownames(df.beta) <- paste0("b",0:(nrow(df.beta)-1) )

	colnames(df.theta) <- c("Estimate","Std Err")
	if (object$cov == "exp") {
		rownames(df.theta) <- c("Nugget","Partial Sill","Range")
	} else if (object$cov == "matern") {
		rownames(df.theta) <- c("Nugget","Partial Sill","Range","Smoothness")
	} else {
		stop(paste0("Unknown covariance type ",object$cov))
	}

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
