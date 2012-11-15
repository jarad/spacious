# create blocks by clustering locations
"blocks.cluster" <- function(S, nblocks) {
	km <- kmeans(S, centers=nblocks)

	# create polygons from clustering
	r.x <- c(floor(min(S[,1])),ceiling(max(S[,1])))
	r.y <- c(floor(min(S[,2])),ceiling(max(S[,2])))

	scale <- 100
	sp  <- expand.grid(seq(r.x[1],r.x[2],length=scale), seq(r.y[1],r.y[2],length=scale))
	ddd <- rdist(sp, km$centers)
	ks  <- apply(ddd, 1, which.min)

	grid <- c()
	for (b in 1:nblocks) {
		if (length(which(ks==b)) == 0) {
			# skip this cluster
			next;
		}

		bpts <- sp[which(ks==b),]
		hull <- chull(bpts)
		poly <- rbind(bpts[hull,], bpts[hull[1],])

		grid <- c(grid,list(Polygons(list(Polygon(cbind(
          poly[,1],poly[,2]
        ))),paste(b)) ))
	}
	grid <- SpatialPolygons(grid)

	# get neighbors
	neighbors <- c()
	sep <- max(c(abs(r.x[2]-r.x[1])/scale,abs(r.y[2]-r.y[1])/scale))
	neighbor.mat <- nb2mat(poly2nb(grid, snap=sep), style='B', zero.policy=T)
	for (i in 1:nrow(neighbor.mat)) {
		for (j in i:ncol(neighbor.mat)) {
			if (neighbor.mat[i,j] == 1) {
				neighbors <- rbind(neighbors, c(i,j) )
			}
		}
	}

	list(B=km$cluster, neighbors=neighbors, grid=grid)
}

if (0) {
bc <- blocks.cluster(cbind(mean_max_temps$lat,mean_max_temps$long), 10)
if (0) {
print(neighbors)
#print(dim(sp));print(length(ks));print(sp[which(ks==1),])
pdf("polys.pdf")
par(mfrow=c(2,2))
plot(S,col=km$cluster,xlim=range(sp[,1]),ylim=range(sp[,2]),main="data and cluster ID")
plot(km$centers,col=1:nblocks,pch=19,cex=2,xlim=range(sp[,1]),ylim=range(sp[,2]),main="cluster centers")
plot(table(km$cluster),main="Number points in each cluster")
plot(sp,col=ks,pch=19,main="domain of each cluster")
	text(km$centers,labels=paste(1:nblocks),col="orange",cex=0.5)
graphics.off()
}

}
