# compare estimation time for GPU vs non-GPU for full likelihood

if ("package:spacious" %in% search()) {
	# unload the package
	detach("package:spacious", unload=TRUE)
	library.dynam.unload("spacious", libpath="~/Rlib/spacious")
}

# load the package
library(spacious)
library(RandomFields)
library(multicore)

"get_timings" <- function(gpu=TRUE, obs_per_block) {

	if (gpu) {
		fn.lapply <- lapply
	} else {
		fn.lapply <- mclapply
		options(cores=3)
	}

	if (missing(obs_per_block)) {
		full <- TRUE
	} else {
		full <- FALSE
	}

	timings <- fn.lapply(sample.sizes, function(i) {
		# generate data
		model <- "exponenti"
		mean <- 5
		variance <- 1
		nugget <- 0.5
		scale <- 1
		x1 <- seq(0, 1, length=round(sqrt(i)))
		x2 <- seq(0, 1, length=round(sqrt(i)))
		S <- as.matrix(expand.grid(rev(x1), x2))
		S <- cbind(S[,2],S[,1])

		set.seed(311)  # make sure GPU vs CPU are estimating the same data
		y <- as.vector( GaussRF(x=x1, y=x2, model=model, grid=TRUE, param=c(mean, variance, nugget, scale)) )
		n <- nrow(S)

		if (full) {
			fit.spacious <- spacious(y~1, S=S, cov="exp", blocks=list(type="full"), nthreads=1, gpu=gpu, verbose=TRUE)
			cat("n=",n,", elapsed=",fit.spacious$time[3],", nIter=",fit.spacious$nIter, ", tpi=", fit.spacious$time[3]/fit.spacious$nIter, "\n", sep="")
			list(n=n, time=fit.spacious$time[3], nIter=fit.spacious$nIter, tpi=fit.spacious$time[3]/fit.spacious$nIter)
		} else {
			nblocks <- round(sqrt((n / obs_per_block)))^2
			fit.spacious <- spacious(y~1, S=S, cov="exp", blocks=list(type="regular", nblocks=nblocks), nthreads=1, verbose=TRUE)
			cat("n=",n,", elapsed=",fit.spacious$time[3],", nIter=",fit.spacious$nIter, ", tpi=", fit.spacious$time[3]/fit.spacious$nIter, ", nblocks=", nblocks, "\n", sep="")
			list(n=n, time=fit.spacious$time[3], nIter=fit.spacious$nIter, tpi=fit.spacious$time[3]/fit.spacious$nIter, nblocks=nblocks)
		}
	})

}

if (FALSE) { # GPU vs CPU (full likelihood)
	sample.sizes <- rev( c(1000, 2000, 3000, 4000, 5000, 6000) )
	timings.gpu <- do.call(rbind, get_timings(TRUE))
	timings.cpu <- do.call(rbind, get_timings(FALSE))
}

if (TRUE) { # compute timings for various obs per block
	sample.sizes <- 25 * rev( c(6,9,11,13:20)^2 )
	timings.b25 <- do.call(rbind, get_timings(FALSE, 25)); write.csv(timings.b25, file="test/timings/b25.csv", row.names=FALSE)
	#sample.sizes <- 100 * rev( seq(3,10)^2 )
	#timings.b100 <- do.call(rbind, get_timings(FALSE, 100)); write.csv(timings.b100, file="test/timings/b100.csv", row.names=FALSE)
}
