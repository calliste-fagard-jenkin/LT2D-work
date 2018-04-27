library('LT2D')

# set perp trunc, forward trunc
# total transect length and survey
# area:
w = 0.15 ; ystart = 0.55
L = 10 ; A = 2*w*L

# set value of 'true' parameters
# for simulated data:
b=c(-7.3287948, 0.9945317)
logphi <- c(0.02,-4.42)

# produce simulated data:
set.seed(3)
simDat = simXY(50, 'pi.norm',
               logphi, 'h1', 
               b, w, 
               ystart)$locs

# create the data.frame:
all.1s <- rep(1,length(simDat$x))
obj <- 1:length(simDat$x)
sim.df <- data.frame(x = simDat$x,
                    y = simDat$y,
                    stratum = all.1s,
                    transect = all.1s,
                    L = L,
                    area = A,
                    object = obj,
                    size = all.1s)

# fit an LT2D model
fit <- LT2D.fit(DataFrameInput = sim.df,
                hr = 'h1',
                # start values for b:
                b = b,
                ystart = ystart,
                pi.x = 'pi.norm',
                # start values for logphi:
                logphi = logphi,
                w = w,
                hessian = TRUE)

gof.LT2D(fit)
boot <- LT2D.bootstrap(fit,r=999,alpha = 0.05)
boot$ci

#used fort he figure in the bootstrap section
hist(boot$Ns,main='',xlab='Estimate of Abundance')
plot(fit)
fit$fit$AIC


set.seed(3)
n <- length(simDat$x)

sim.df$observer.id <- factor(sample(
  1:3,n,replace=T))

sim.df$altitude <-  rnorm(n,2,1)


iform <- formula(i~sim.df$observer.id)
xform <- formula(x~sim.df$altitude)
yform <- formula(y~sim.df$altitude)
ParamNumRequired(sim.df, iform)

# covariates in intercept only
fit2 <- LT2D.fit(DataFrameInput = sim.df,
                 hr = 'h1',
                 # start values for b:
                 b = c(-3.72,0.74),
                 ystart = ystart,
                 pi.x = 'pi.norm',
                 # start values for logphi:
                 logphi = c(0.02, -4.42),
                 w = w,
                 hessian = TRUE,
                 formulas = list(iform),
                 ipars=rep(0,2))

# covariates in the x and y directions
# (which must have the same effect in 
# an h1 detection function)
fit3 <- LT2D.fit(DataFrameInput = sim.df,
                 hr = 'h1',
                 # start values for b:
                 b = c(-3.72,0.74),
                 ystart = ystart,
                 pi.x = 'pi.norm',
                 # start values for logphi:
                 logphi = c(0.02, -4.42),
                 w = w,
                 hessian = TRUE,
                 formulas=list(xform,yform),
                 xpars=0,
                 ypars=0)

# covariates in all of the h1 slots
fit4 <- LT2D.fit(DataFrameInput = sim.df,
                 hr = 'h1',
                 # start values for b:
                 b = c(-3.72,0.74),
                 ystart = ystart,
                 pi.x = 'pi.norm',
                 # start values for logphi:
                 logphi = c(0.02, -4.42),
                 w = w,
                 hessian = TRUE,
                 formulas=list(iform,xform,yform),
                 ipars=rep(0,2),
                 xpars=0,
                 ypars=0)

plot(fit4, covar.row=1,smooth.fy=T)
boot <- LT2D.bootstrap(fit4,r=999,alpha=0.05)
boot$ci
############### SIMULATION STUDY #############################
set.seed=(0)
S <- 250
N <- 50
nocovars <- rep(NA, S)
covars <- rep(NA, S)
Lsim <- L
Asim <- A
n1s <- rep(NA, S)
n2s <- rep(NA, S)

for (i in 1:S){
  simDat1 <- simXY(N, 'pi.norm', logphi, 'h1', c(b[1],b[2]), w, ystart)$locs
  simDat2 <- simXY(N, 'pi.norm', logphi, 'h1', c(b[1]+2,b[2]), w, ystart)$locs
  
  n1 <- length(simDat1$x)
  n2 <- length(simDat2$x)
  
  # We can consider fakeFactor to be a factor covariate representing species, 
  # where the second species has a far higher detectability.
  sim.df2 <- data.frame(x = c(simDat1$x, simDat2$x),
                       y = c(simDat1$y, simDat2$y),
                       stratum = rep(1,n1+n2),
                       transect = rep(1, n1+n2),
                       L = Lsim,
                       area = Asim,
                       object = 1:(n1+n2),
                       size = rep(1,n1+n2),
                       fakeFactor = factor(rep(c(1,2),c(n1,n2))))
  
  # with covariates:
  T6 <- try(LT2D.fit(DataFrameInput = sim.df2, hr = 'h1', b=c(1,-0.55),
                ystart=ystart, pi.x='pi.norm', logphi=logphi, w=w,
                formulas = list(formula(i~fakeFactor)), ipars = c(0)),
           silent=T)
  
  # without taking covariates into account:
  T7 <- try(LT2D.fit(DataFrameInput = sim.df2, hr = 'h1', b=c(1,-0.55),
                ystart=ystart, pi.x='pi.norm', logphi=logphi, w=w),
           silent=T)
  
  n1s[i] <- n1
  n2s[i] <- n2
  
  if(class(T6)!='try-error' & class(T6)!='try-error'){
    covars[i] <- T6$ests$N[2]
    nocovars[i] <- T7$ests$N[2]
  }
}

indicesC <- which(is.na(covars))
indicesNC <- which(is.na(nocovars))
covars2 <- covars[-indicesC]
nocovars2 <- nocovars[-indicesNC]

u.ind <- rep(0, length=length(covars2))
for (i in 1:length(u.ind)){
  if (covars2[i]>10*N | nocovars2[i]>10*N) u.ind[i] <- 1
}

((mean(covars2)-2*N)/(2*N))*100
((mean(nocovars2)-2*N)/(2*N))*100

covars3 <- covars2[-which(u.ind==1)]
nocovars3 <- nocovars2[-which(u.ind==1)]

((mean(covars3)-2*N)/(2*N))*100
((mean(nocovars3)-2*N)/(2*N))*100

closer <- rep(NA,length(covars3))

for(i in 1:length(closer)){
  if (abs(2*N-covars3[i])<abs(2*N-nocovars3[i])) closer[i] <- TRUE
  else if (abs(2*N-covars3[i])==abs(2*N-nocovars3[i])) print('equal')
  else closer[i] <- FALSE
}

closer <- na.omit(closer)
print(sum(closer)/length(closer))
##############################################################

# Plots in presentation:
x = seq(0,1,length=100)
y = pi.hnorm(x,-0.8,1)
plot(x,y,type='l',lty=2,lwd=2,xlab='perpendicular distance',
     ylab='density')

x = seq(0,1,length=100)
y = pi.norm(x,c(7,0.7),1)
lines(x,y,type='l',lty=3,lwd=2)
lines(x,rep(0.94,100),lwd=2)

n=100
x = seq(0.01,0.03,length=n)
y = seq(0.15,0.60,length=n)
z = matrix(nrow=n,ncol=n)
for(i in 1:n){
  for(j in 1:n){
    #z[i,j] = h1(y[j],x[i],list(-7.3287948, 0.9945317))
    z[i,j] = ep2(y[j],x[i],
                 list(10, -3, 1, 0)
                 )
  }
}
hist3D(x,y,z)

################################################################

#######################################################################
set.seed(1)
simDat1 <- simXY(500, 'pi.norm', logphi, 'ep1', c(20,10,10), w, ystart)$locs
simDat2 <- simXY(3*50, 'pi.norm', logphi, 'ep1', c(20,10,1), w, ystart)$locs

n1 <- length(simDat1$x)
n2 <- length(simDat2$x)

# We can consider fakeFactor to be a factor covariate representing species, 
# where the second species has a far higher detectability.
sim.df2 <- data.frame(x = c(simDat1$x, simDat2$x),
                      y = c(simDat1$y, simDat2$y),
                      stratum = rep(1,n1+n2),
                      transect = rep(1, n1+n2),
                      L = Lsim,
                      area = Asim,
                      object = 1:(n1+n2),
                      size = rep(1,n1+n2),
                      fakeFactor = factor(rep(c(1,2),c(n1,n2))))

# with covariates:
T6 <- try(LT2D.fit(DataFrameInput = sim.df2, hr = 'h1', b=c(1,-0.55),
                   ystart=ystart, pi.x='pi.norm', logphi=logphi, w=w,
                   formulas = list(formula(i~fakeFactor)), ipars = c(-2)),
          silent=T)