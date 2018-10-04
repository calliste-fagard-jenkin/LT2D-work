library(plotrix) # load the library needed to draw the arc on our graph

library('DEVDEV5') # load the development version of the LT2D package
Data = read.csv('GibbonSiamang.csv',header=T) # load the x and y data
Data = subset(Data, Data$Species == 1)
x = Data$PP.Distance ; y = Data$Forward.Distance

par(mfrow=c(1,3))

# we want to add a small amount of noise to data points at (0,0) to make it 
# clearer in the graph that there is a cluster here:
x_to_plot = x                                      # create vectors to 
y_to_plot = y                                      # fill with values

for (i in (1:length(x))){                          # loop through and add noise
  if (x[i]==0 & y[i]==0){                          # where it's required
    x_to_plot[i] = x[i] + runif(1,0,0.002)
    y_to_plot[i] = y[i] + runif(1,0,0.002)
  }         
}

pdlab="Perpendicular distance"
fdlab="Forward distance"
plot(x_to_plot,y_to_plot,pch="+",ylab=fdlab,xlab=pdlab)
rmin = sort(unique(Data$Radial.Distance))[2] # Calculate Rmin
draw.arc(0,0,radius=rmin,deg1=0,deg2=360,col="red",lwd=1)
legend(0.0113,0.065,"Radius Rmin",lty=1,col="red")
hist(y,breaks=seq(0,max(na.omit(y)),length=9),xlab=fdlab,main="")
hist(x,breaks=seq(0,max(na.omit(x)),length=13),xlab=pdlab,main="")

par(mfrow=c(1,1))
R = Data$Radial.Distance
hist(R, breaks = 22, main="Radial distances of detections")

w=0.05;ystart=0.065 # set these to use in all our fits

stratum = rep(1, length(x))
transect = rep(1, length(x))
object = seq(1, length(x))
size = rep(1, length(x))
L = rep(4*52, length(x))    # do we times this by the 52 visits? 
area = rep(4*52*2*w, length(x))

df = data.frame(y,x,stratum, transect, object, size, L, area)
df.zeros = df # so that we can use the data set with the zeros later, after 
# we've changed what df is

# H1 in the LT2D package is the Haynes and Buckland Hazard rate. 
# Normal bump with hazard h1:
b=c(-6.82664896,1.09392846) ; logphi=c(-0.03503625,0.51044374)  # start values

fit.h1 = LT2D.fit(DataFrameInput=df,hr='h1',b=b,ystart=ystart,pi.x='pi.norm',
                  logphi=logphi,w=w, rmin=rmin, hessian=TRUE)

fit.h1$ests
# par(mfrow=c(1,1))
# gof.LT2D(fit.n,plot=T)
# plot(fit.n,xbins=20,ybins=32,smooth.fy=TRUE,addrug=TRUE) 
# phatInterval(fit.n)*w #EHSW
# fit.n$p0 # p(0)

# Normal dip: --  bad ks gof in x --  problem with hessian inversion 
b=c(-6,1) ; logphi=c(1,-4)
fit.h1.chn = LT2D.fit(DataFrameInput=df,hr='h1',b=b,ystart=ystart,pi.x='pi.chnorm',
                      logphi=logphi,w=w, rmin=rmin, hessian=TRUE)
fit.h1.chn$ests

# Uniform with h1 -- bad ks gof in x -- 
b=c(-7.0744154,0.9876447)                 
fit.h1.unif = LT2D.fit(DataFrameInput=df,hr='h1',b=b,ystart=ystart,pi.x='pi.const',
                       logphi=NULL,w=w, rmin=rmin, hessian=TRUE)
fit.h1.unif$ests

df = data.frame(x,y)                  # fast way to remove all the 
df = subset(df, sqrt(x**2+y**2)>rmin) # (0,0) values
x = df$x ; y = df$y 

stratum = rep(1, length(x))
transect = rep(1, length(x))
object = seq(1, length(x))
size = rep(1, length(x))
L = rep(4*52, length(x))    # do we times this by the 52 visits? 
area = rep(4*52*2*w, length(x))

df = data.frame(y,x,stratum, transect, object, size, L, area)

# Normal bump with hazard h1:
b=c(-6.82664896,1.09392846) ; logphi=c(-0.03503625,0.51044374)  # start values

fit.h1.nz = LT2D.fit(DataFrameInput=df,hr='h1',b=b,ystart=ystart,pi.x='pi.norm',
                     logphi=logphi,w=w, hessian=TRUE)
fit.h1.nz$ests

# Normal dip: --- doesn't fit ---
b=c(-6,1) ; logphi=c(1,-3)
fit.h1.chn.nz = LT2D.fit(DataFrameInput=df,hr='h1',b=b,ystart=ystart,pi.x='pi.chnorm',
                         logphi=logphi,w=w, hessian=TRUE)
fit.h1.chn.nz$ests

# Uniform with h1 -- bad ks gof in x -- 
b = c(-7,1)
fit.h1.unif.nz = LT2D.fit(DataFrameInput=df,hr='h1',b=b,ystart=ystart,pi.x='pi.const',
                          logphi=NULL,w=w, hessian=TRUE)
fit.h1.unif.nz$ests

# H1 AICs
fit.h1.nz$fit$AIC ; fit.h1.chn.nz$fit$AIC ; fit.h1.unif.nz$fit$AIC

library(Distance)
#df = df.zeros - we want to use the dataset with the zeros excluded 
# We add some columns to have the correct names
df$Area = df$area ; df$Sample.Label = df$transect ; df$Region.Label = df$stratum
df$distance = df$x ; df$Effort = df$L
df$x = NULL ; df$y = NULL


# fitting models:
mod.hn <- ds(data=df, transect = "line",
             key = "hn", adjustment = NULL, truncation = w, quiet = TRUE)
mod.hz <- ds(data=df, transect = "line",
             key = "hr", adjustment = NULL,truncation = w, quiet = TRUE)
mod.hn.cos <- ds(data=df, transect = "line",
                 key = "hn", adjustment = "cos",truncation = w, quiet = TRUE)
mod.hr.cos <- ds(data=df, transect = "line",
                 key = "hr", adjustment = "cos",truncation = w, quiet = TRUE)

summarize_ds_models(mod.hn, mod.hz, mod.hn.cos, mod.hr.cos)

GOF = ds.gof(mod.hz, main = 'CDS') # CDS fit seems very reasonable... 
gof.LT2D(fit.h1.chn.nz, plot = TRUE)
plot(mod.hz, xlab="Distance / Km", main='CDS')
plot(fit.h1.unif.nz, xbins = 7)

h1.to.HB(fit.h1.unif.nz$fit$par)
h1.to.HB(exp(fit.h1.unif.nz$fit$par))

mod.hz$ddf$par