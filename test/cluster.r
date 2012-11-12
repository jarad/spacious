

n<-1000
M<-25
group<-sample(1:5,n,replace=T,prob=(1:5)^2)

x<-c(rnorm(n-10,group,.3),rnorm(10,3,3))
y<-c(rnorm(n-10,group,.3),rnorm(10,3,3))
s<-cbind(x,y)
d<-dist(s)

plot(s)

km<-kmeans(s,centers=M)


sp<-expand.grid(seq(-5,8,length=100),
                seq(-5,8,length=100))

library(fields)
ddd<-rdist(sp,km$centers)
ks<-apply(ddd,1,which.min)

pdf("test.pdf")
par(mfrow=c(2,2))
plot(s,col=km$cluster,xlim=range(sp[,1]),ylim=range(sp[,2]),main="data and cluster ID")
plot(km$centers,col=1:M,pch=19,cex=2,xlim=range(sp[,1]),ylim=range(sp[,2]),main="cluster centers")
plot(table(km$cluster),main="Number points in each cluster")
plot(sp,col=ks,pch=19,main="domain of each cluster")
graphics.off()


