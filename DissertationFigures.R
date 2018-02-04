set.seed(45)

# generate points
x = runif(100)
y = runif(100)

# Draw them
plot(x,y,pch=16,xlim=c(0,1),ylim=c(0,1),
     bty='n',xaxt='n',yaxt='n',ann=F)

# Add survey area box horizontal lines
segments(x0=c(0,0,0,1),
         y0=c(0,0,1,0),
         x1=c(0,1,1,1),
         y1=c(1,0,1,1))

# Add strip lines and vertical box lines
segments(x0=seq(0,1,length=6),
         y0=rep(0,6),
         x1=seq(0,1,length=6),
         y1=rep(1,6))