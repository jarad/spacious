# create blocks by clustering locations
"blocks.cluster" <- function(S, nblocks) {
	km <- kmeans(S, centers=S[round(seq(1,nrow(S),length=nblocks)),])

	# create polygons from clustering
	r.x <- c(floor(min(S[,1])),ceiling(max(S[,1])))
	r.y <- c(floor(min(S[,2])),ceiling(max(S[,2])))

	scale <- 100 #max(100,nblocks)
	sp  <- expand.grid(seq(r.x[1],r.x[2],length=scale), seq(r.y[1],r.y[2],length=scale))
	ddd <- rdist(sp, km$centers)
	ks  <- apply(ddd, 1, which.min)

	# create tessellation from centers
	tv  <- deldir(km$centers[,1], km$centers[,2], list(ndx=0, ndy=0), c(r.x[1],r.x[2],r.y[1],r.y[2]))
	tvl <- tile.list(tv)

	B    <- km$cluster
	grid <- c()
	for (b in 1:nblocks) {
		poly.x <- c(tvl[[b]]$x, tvl[[b]]$x[1])
		poly.y <- c(tvl[[b]]$y, tvl[[b]]$y[1])

		grid <- c(grid,list(Polygons(list(Polygon(cbind( poly.x, poly.y ))),paste(b)) ))

		in_poly <- point.in.polygon(S[,1], S[,2], poly.x, poly.y) >= 1

		if (sum(in_poly) == 0) { next; }

		B[in_poly] <- b
	}
	grid <- SpatialPolygons(grid)

	# get neighbors
	neighbors <- c()
	sep <- max(c(abs(r.x[2]-r.x[1])/scale,abs(r.y[2]-r.y[1])/scale))
	neighbor.mat <- nb2mat(poly2nb(grid), style='B', zero.policy=TRUE)
	for (i in 1:nrow(neighbor.mat)) {
		for (j in i:ncol(neighbor.mat)) {
			if (neighbor.mat[i,j] == 1) {
				neighbors <- rbind(neighbors, c(i,j) )
			}
		}
	}

	list(B=B, neighbors=neighbors, grid=grid)
}
