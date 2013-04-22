if (FALSE) { # compares full to various block sizes

timing <- rbind(
        c(3.652,72.628,356.954,2098.471,3771.413,14370.3,NA,NA,NA),        # full
        c(2.724,7.428,12.482,27.121,33.251,47.728,64.492,102.969,103.83),  # 100 per block
        c(0.359,0.793,1.51,2.466,2.926,4.085,5.402,6.895,8.639)            # 25 per block
)

sizes <- 100*(2:10)^2
blocks <- c(1,100,25)
sb <- c("Full Likelihood","100 Obs per Block","25 Obs per Block")
cols <- c("black","red","blue")

# cube root these
timing <- (timing/60/60)^(1/3)
#at <- c(0.05,0.10,0.25,0.5,.1,seq(1,4,by=1))^(1/3)
#tl <- at^3
at <- (c(1,15,60,5*60,15*60,30*60,60*60,2*60*60,3*60*60,4*60*60)/60/60)^(1/3)
tl <- c("1s","15s","1m","5m","15m","30m","1h","2h","3h","4h")

# plot timings
pdf("figures/timing.pdf");
	#plot(sizes, timing[1,]/60/60,type="l",lty=1, ylim=c(min(timing),max(timing))/60/60,
	plot(sizes, timing[1,],type="b",lty=1, xlim=c(0, 10000), ylim=c(min(timing,na.rm=TRUE),max(timing,na.rm=TRUE)),
		main="Execution Time vs Sample Size (Cube Root Scale)",
		xlab="Sample Size",ylab="Execution Time", xaxt="n", yaxt="n", lwd=3, col=cols[1]);
	for (j in 2:nrow(timing)) {
		#lines(sizes, timing[j,]/60/60,lty=j);
		lines(sizes, timing[j,],lty=1, lwd=3, col=cols[j],type="b");
	}
	legend("topleft",sb,lty=1,ncol=1,inset=0.05,col=cols,lwd=3)
	axis(1, at = c(0,seq(1000,10000,by=1000)), labels=c("0",paste0(1:10,"k")), cex.axis=1 )
	axis(2, at = at, labels = tl, las = 2, cex.axis=1)
graphics.off();
}

if (TRUE) { # compares thread timings

times <- read.csv("test/times_1k_to_1m.csv")
n <- unique(times$n) #c(1024,5041,10000,19881,50176,99856,250000)
nthreads <- unique(times$nthreads) #c(24,16,8,4,2,1)
cols <- rgb(0, seq(0.25,1,length=length(nthreads)), 0)
symbs <- c("T", "8", "4", "2", "1")

timing <- matrix(NA, nrow=length(nthreads), ncol=length(n))
sapply(1:length(nthreads), function(i) {
	timing[i,] <<- times[times$nthreads==nthreads[i],]$time
})

at <- c(10*60,30*60,1*60*60,1.5*60*60,2*60*60,2.5*60*60,3*60*60,3.5*60*60,4*60*60,4.5*60*60,5*60*60)
tl <- c("10m","30m","1h","1.5h","2h","2.5h","3h","3.5h","4h","4.5h","5h")

pdf("pdf/timing_threads.pdf");
	plot(n, timing[1,],type="b",lty=1, xlim=c(0, max(n)), ylim=c(min(timing,na.rm=TRUE),max(timing,na.rm=TRUE)),
		main="Execution Time vs Sample Size (100 Obs per Block)",
		#xlab="Sample Size",ylab="Execution Time",
		xlab="",ylab="",
		xaxt="n", yaxt="n", lwd=3, col=cols[1], pch=symbs[1], cex=1.25);
	for (j in 2:nrow(timing)) {
		lines(n, timing[j,],lty=1, lwd=3, col=cols[j], type="b", pch=symbs[j], cex=1.25);
	}
	#legend("topleft",paste(nthreads,"threads"),lty=1,ncol=1,inset=0.05,col=cols,pch=symbs,lwd=2,merge=TRUE)
	legend("topleft",paste(nthreads,"threads"),ncol=1,inset=0.05,col=cols,pch=symbs,cex=1.5)
	axis(1, at = c(1000,100000,250000,500000,1000000), labels=c("1k","100k","250k","500k","1m"), cex.axis=1.5)
	axis(2, at = at, labels = tl, las = 2, cex.axis=1.5)
graphics.off();
}
