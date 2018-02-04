#################### FIGURE 1
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



#################### FIGURE 2
set.seed(46)

hn <- function(x,sigma){
  # half normal detection function
  return(exp(-(x**2)/(2*sigma**2)))
}

AnimalDistances <- runif(1000,min=0,max=100)

DetectionProbs <- hn(AnimalDistances, 9)
acceptance <- vector()

for(i in 1:length(DetectionProbs)){
  if (rbinom(1,1,DetectionProbs[i])){
    acceptance <- append(acceptance,
                         AnimalDistances[i])
  }
}

hist(acceptance,main='',
     xlab='Perpendicular Distance')

################## FIGURE 3


hist(acceptance,main='',
     xlab='Perpendicular Distance')

x <- seq(0,30,length=500)
y1 <- rep(41,500)
y2 <- hn(x,9)*41

# polygon(c(x,rev(x)),
#         c(y2,rev(y1)),col=rgb(1,0,0,0.5))

polygon(c(x,rev(x)),
        c(y2,rep(0,500)),col=rgb(0,0,1,0.4))

curve(hn(x,sigma=9)*41,add=TRUE)
      
segments(x0=c(0,30),
         y0=c(41,41),
         x1=c(30,30),
         y1=c(41,0))

