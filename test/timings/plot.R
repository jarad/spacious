# plots timings

if (FALSE) {
	# plot timing of full likelihood for GPU vs CPU

	t.cpu <- read.csv("test/timings/full_cpu.csv")
	t.gpu <- read.csv("test/timings/full_gpu.csv")

	nIter <- round(mean(t.cpu$nIter))   # use same number of iterations for each

	timings.cpu <- (nIter*t.cpu$tpi/60/60)^(1/3)
	timings.gpu <- (nIter*t.gpu$tpi/60/60)^(1/3)

	y.lims <- range(c(timings.cpu,timings.gpu))

	axis.at     <- (c(1,15,60,5*60,15*60,30*60,60*60,2*60*60,3*60*60,4*60*60)/60/60)^(1/3)
	axis.labels <- c("1s","15s","1m","5m","15m","30m","1h","2h","3h","4h")

	pdf("pdf/timing_full_gpu_vs_cpu.pdf");
		plot(t.cpu$n, timings.cpu, type="b", lty=1, xlim=c(1000, 6000), ylim=y.lims,
			main="Execution Time vs Sample Size (Cube Root Scale)",
			xlab="Sample Size", ylab="Execution Time", xaxt="n", yaxt="n", lwd=3, col="black");

		lines(t.gpu$n, timings.gpu, lty=1, lwd=3, col="red", type="b");

		legend("topleft",c("CPU","GPU"), lty=1, ncol=1, inset=0.05, col=c("black","red"), lwd=3)

		axis(1, at = seq(1000,6000,by=1000), labels=c(paste0(1:6,"k")), cex.axis=1 )
		axis(2, at = axis.at, labels = axis.labels, las = 2, cex.axis=1)
	graphics.off();
}

if (TRUE) {
	# plot timing of full likelihood for CPU vs 25 and 100 obs per block

	t.full <- read.csv("test/timings/full_cpu.csv")
	t.gpu  <- read.csv("test/timings/full_gpu.csv")
	t.b100 <- read.csv("test/timings/b100.csv")
	t.b25  <- read.csv("test/timings/b25.csv")

	nIter.full  <- round(mean(t.full$nIter))   # use same number of iterations for each
	nIter.block <- round(mean(c(t.b100$nIter,t.b25$nIter)))

	timings.full <- (nIter.full*t.full$tpi/60/60)^(1/3)
	timings.gpu  <- (nIter.full*t.gpu$tpi/60/60)^(1/3)
	timings.b100 <- (nIter.block*t.b100$tpi/60/60)^(1/3)
	timings.b25  <- (nIter.block*t.b25$tpi/60/60)^(1/3)

	y.lims <- range(c(timings.full,timings.gpu,timings.b100,timings.b25))

	axis.at     <- (c(1,15,60,5*60,15*60,30*60,60*60,2*60*60,3*60*60,4*60*60)/60/60)^(1/3)
	axis.labels <- c("1s","15s","1m","5m","15m","30m","1h","2h","3h","4h")

	pdf("pdf/timing_full_vs_blocks.pdf");
		plot(t.full$n, timings.full, type="b", lty=1, xlim=c(900, 10000), ylim=y.lims,
			main="Execution Time vs Sample Size (Cube Root Scale)",
			xlab="Sample Size", ylab="Execution Time", xaxt="n", yaxt="n", lwd=3, col="black");
		#lines(t.gpu$n, timings.gpu, lty=2, lwd=3, col="black", type="b");

		lines(t.b100$n, timings.b100, lty=1, lwd=3, col="blue", type="b");
		lines(t.b25$n, timings.b25, lty=1, lwd=3, col="red", type="b");

		legend("topright",c("Full","100 Obs per Block","25 Obs per Block"),lty=c(1,1,1), ncol=1, inset=0.05, col=c("black","blue","red"), lwd=3)
		#legend("topright",c("Full (CPU)","Full (GPU)","100 Obs per Block","25 Obs per Block"),lty=c(1,2,1,1), ncol=1, inset=0.05, col=c("black","black","blue","red"), lwd=3)

		axis(1, at = seq(1000,10000,by=1000), labels=c(paste0(1:10,"k")), cex.axis=1 )
		axis(2, at = axis.at, labels = axis.labels, las = 2, cex.axis=1)
	graphics.off();
}
