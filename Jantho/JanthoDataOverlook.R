```{r, warning=F}


library('LT2D')
Data = read.csv('JanthoData.csv',header=T)

DataBySpecies = list()                             # break up data by species
species.number=1
while (dim(subset(Data,Data$species==species.number))[1]!=0){
  DataBySpecies[[species.number]]=subset(Data,Data$species==species.number)
  species.number = species.number + 1
}


par(mfrow=c(2,2))
n = dim(Data)[1]
hist(Data$Forward.Distance,
     main='Forward Distance',
     xlab='Distance',col=rgb(1,0,0,0.5))
legend(0.022, 30,paste('n =',n),bty='n')
hist(Data$PP.Distance,
     main='Perpendicular Distance',
     xlab='Distance',col=rgb(0,0,1,0.5))
hist(Data$Radial.Distance,
     main='Radial Distance',
     xlab='Distance',col=rgb(0,1,0,0.5))
hist(Data$Bearing..0.180..,
     main='Bearing 0-180',
     xlab='Distance',col=rgb(1,1,0,0.5))

for (i in DataBySpecies){
  spec.name = i$species[1]
  n = dim(i)[1]
  dat = c(i$Forward.Distance, i$PP.Distance, i$Radial.Distance)
  dat = na.omit(dat)
  maxval = max(dat)
  ylim = c(0,10)
  
  hist(i$Forward.Distance,
       main = paste('Species',as.character(spec.name),'Frwd Dist'),
       xlim = c(0,maxval), xlab = 'Distance',
       col = rgb(red=1,green=0,blue=0,alpha=0.5)
  )
  legend(0.016, 8,paste('n =',n),bty='n')
  
  hist(i$PP.Distance,
       main = paste('Species',as.character(spec.name),'Perp Dist'),
       xlab='Distance',
       col = rgb(red=0,green=0,blue=1,alpha=0.5)
  )
  
  #   legend.coordinates = list()
  #   legend.coordinates[[1]] = c(0.035,10)
  #   legend.coordinates[[2]] = c(0.025,11)
  #   legend.coordinates[[3]] = c(0.022,10)
  #   legend.coordinates[[4]] = c(0.025,10)
  #   legend(x=legend.coordinates[[spec.name]][1],
  #          y=legend.coordinates[[spec.name]][2],
  #          legend=c('Forward','Perpendicular'),
  #          fill=c(rgb(1,0,0,0.5),rgb(0,0,1,0.5)),
  #          cex = 0.8
  #         )
  
  hist(i$Radial.Distance,
       main = paste('Species',spec.name,'Radial Distance'),
       col = rgb(red=0,green=1,blue=0,alpha=0.5),
       xlab = 'Distance'
  )
  plot(i$PP.Distance,i$Forward.Distance,xlab='Perpendicular Distance',
       ylab='Forward Distance', main='')
}


# Analysis assuming constant detection between species...
Data$x = Data$PP.Distance ; Data$y = Data$Forward.Distance
Data$transect = Data$Sample.Label ; Data$Sample.Label = NULL
Data$stratum = Data$Region.Label ; Data$Region.Label = NULL
Data$Bearing..Rel..to.North. = NULL
Data$Bearing..0.180.. = NULL
Data$PP.Distance = NULL ; Data$Forward.Distance = NULL 
Data$Area = NULL

GetDataSetEffort = function(Data){
  TransectTotal = 0
  for (i in c((1:20)[-c(11,14)])){
    A = subset(Data, Data$transect==i)$Effort[1]
    TransectTotal = TransectTotal + A
  }
  return(TransectTotal)
}

addObjectNumbers = function(Dataset){
  length = length(Dataset$x)
  Dataset$object = rep(NA, length)
  Dataset$size = rep(NA, length)
  
  counter = 1 
  for (i in (1:length)){ # loop through rows
    if (!is.na(Dataset$x[i])){
      Dataset$object[i] = counter 
      Dataset$size[i] = 1
      counter = counter + 1 
    }
  }
  return(Dataset)
}

w = 0.05 ; ystart = max(na.omit(Data$y))
Data = addObjectNumbers(Data)
Data$L = Data$Effort
Data$area = rep(2*w*GetDataSetEffort(Data))

## fitting H1 detection functions
b=c(-2.3,0.4) ; logphi=c(0,-4.4)

Data2 = subset(Data, !is.na(Data$x))
fityx(Data2$y, Data2$x, b, 'h1',ystart,'pi.norm', logphi, w, hessian = F)
fit.h1.norm = LT2D.fit(DataFrameInput=Data,hr='h1',b=b,ystart=ystart,
                       pi.x='pi.norm',logphi=logphi,w=w, hessian=TRUE)

fit.h1.unif = LT2D.fit(DataFrameInput=Data,hr='h1',b=b,ystart=ystart,
                       pi.x='pi.const',logphi=NULL,w=w, hessian=TRUE)

fit.h1.unif$fit$par 

b=c(-9.5,-0.7) ; logphi=c(0.5,-3.9) # HESSIAN PROBLEM
fit.h1.chnorm = LT2D.fit(DataFrameInput=Data,hr='h1',b=b,ystart=ystart,
                         pi.x='pi.chnorm',logphi=logphi,w=w, hessian=F)

b=c(-7,1) ; logphi = c(-5,-0.1)
fit.h1.chnorm.i = LT2D.fit(DataFrameInput=Data,hr='h1',b=b,ystart=ystart,
                           pi.x='pi.chnorm.i',logphi=logphi,w=w, hessian=T)
# par(mfrow=c(2,2))
# plot(fit.h1.chnorm.i)
# gof.LT2D(fit.h1.chnorm.i)
#fit.h1.chnorm.i$ests ; fit.h1.chnorm.i$fit$AIC

## Analyses with truncation 0.03 and H1 detection 
wtrunc = 0.03 ; logphi = c(-5,-0.1)
fit.h1.chnorm.i.03trunc = LT2D.fit(DataFrameInput=Data,hr='h1',b=b,
                                   ystart=ystart, pi.x='pi.chnorm.i',
                                   logphi=logphi,w=wtrunc, hessian=T)
plot(fit.h1.chnorm.i.03trunc)
gof.LT2D(fit.h1.chnorm.i.03trunc, plot=T)
b=c(-2.3,0.4) ; logphi=c(0,-4.4) # BEST AIC for 0.03 H1 models
fit.h1.norm.i.03trunc = LT2D.fit(DataFrameInput=Data,hr='h1',b=b,
                                 ystart=ystart, pi.x='pi.norm',
                                 logphi=logphi,w=wtrunc, hessian=T)
par(mfrow=c(1,1))
fit.h1.const.i.03trunc = LT2D.fit(DataFrameInput=Data,hr='h1',b=b,
                                  ystart=ystart, pi.x='pi.const',
                                  logphi=NULL,w=wtrunc, hessian=T)
#plot(fit.hi.chnorm.i.03trunc)
#gof.LT2D(fit.hi.chnorm.i.03trunc,plot=T)

b=c(-7.8,0.3,8.8) ; logphi=c(0.1,-4.6)
fit.ip1.norm = LT2D.fit(DataFrameInput=Data,hr='ip1',b=b,ystart=ystart,
                        pi.x='pi.norm',logphi=logphi,w=wtrunc, hessian=T)

fit.ip1.unif = LT2D.fit(DataFrameInput=Data,hr='ip1',b=b,ystart=ystart,
                        pi.x='pi.const',logphi=NULL,w=wtrunc, hessian=F)

b=c(-4.046,0.1295,8.641)
fit.ip1.chnorm = LT2D.fit(DataFrameInput=Data,hr='ip1',b=b,ystart=ystart,
                          pi.x='pi.chnorm',logphi=logphi,w=wtrunc, hessian=F)

logphi = c(-5,-0.1)
fit.ip1.chnorm.i = LT2D.fit(DataFrameInput=Data,hr='ip1',b=b,ystart=ystart,
                      pi.x='pi.chnorm.i',logphi=logphi,w=wtrunc, hessian=F)

summarise.LT2D.models()
```