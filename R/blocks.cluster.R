# create blocks by clustering locations
"blocks.cluster" <- function(S, nblocks) {
	km <- kmeans(S, centers=nblocks)

	# create polygons from clustering
	r.x <- c(floor(min(S[,1])),ceiling(max(S[,1])))
	r.y <- c(floor(min(S[,2])),ceiling(max(S[,2])))

	scale <- 100 #max(100,nblocks)
	sp  <- expand.grid(seq(r.x[1],r.x[2],length=scale), seq(r.y[1],r.y[2],length=scale))
	ddd <- rdist(sp, km$centers)
	ks  <- apply(ddd, 1, which.min)

	grid <- c()
	for (b in 1:nblocks) {
		if (length(which(ks==b)) == 0) {
			# no points in grid for this block; use actual points for polygon
			bpts <- matrix(S[which(km$cluster==b),], nrow=sum(km$cluster==b), ncol=2)
			hull <- chull(bpts)
			poly <- rbind(bpts[hull,], bpts[hull[1],])

			grid <- c(grid,list(Polygons(list(Polygon(cbind(
				poly[,1],poly[,2]
			))),paste(b)) ))

			next;
		} else {
			# use points from grid for polygon
			bpts <- sp[which(ks==b),]
			hull <- chull(bpts)
			poly <- rbind(bpts[hull,], bpts[hull[1],])

			grid <- c(grid,list(Polygons(list(Polygon(cbind(
				poly[,1],poly[,2]
			))),paste(b)) ))
		}
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
