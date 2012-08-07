# predict values at new sites using fit from spacious block composite model
"predict.spacious" <- function(object, newdata=NULL, newS=NULL,
	interval = c("none", "confidence", "prediction"), level = 0.95, ...) {

	cat("Spacious predict!\n")

	# figure out which block the location is in
print(newS)
	nB <- length(object$grid)
	if (nB == 0) {  # no blocks
	} else {
		which.block <- NULL
		for (i in 1:nB) {
		}
	}
}
