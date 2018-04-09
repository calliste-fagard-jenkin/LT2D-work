# This document acts as proof of functionality for the LT2D package covariate
# functionality additions

library(LT2D)
source('Covar Stuff.R')
# First we create a couple plain data sets (no covariates) and test that 
# abundance estimation is reasonable, and that poorer models work less well.
# After this we will examine the effect of covariates on the data:

ystart = 0.1 ; w = 0.05 # These 2 parameters will remain constant for all data

# We consistently pick parameters that give us an average detection probability
# of around 0.8, to provide good but realistic data for testing.
set.seed(1)
simDat = simXY(500, 'pi.const', NULL, 'h1', c(1,-0.55), w, ystart)$locs
Lsim = 20 ; Asim = Lsim*w*2

sim.df = data.frame(x = simDat$x, y = simDat$y, stratum=rep(1,length(simDat$x)),
                    transect = rep(1,length(simDat$x)), L = Lsim, area = Asim,
                    object = 1:length(simDat$x),size = rep(1, length(simDat$x)))

logphi = c(0.01322874, -4.80826844)
# we fit a model with the wrong perpendicular density:
T1 = LT2D.fit(DataFrameInput = sim.df, hr = 'h1', b=c(1,-0.55),
              ystart=ystart,pi.x='pi.norm', logphi=logphi, w=w, hessian=T)

# we fit a model with the wrong detection function:
T2 = LT2D.fit(DataFrameInput = sim.df, hr = 'ep1',
              b=c(12.25,4.41,-24),ystart=ystart,
              pi.x='pi.const', logphi=NULL, w=w, hessian=T)

# we fit a model with the correct perpendicular density and hazard:
T3 = LT2D.fit(DataFrameInput = sim.df, hr = 'h1', b=c(1,-0.55), ystart=ystart,
              pi.x='pi.const', logphi=NULL, w=w, hessian=T)


# As we can see by looking at the abundance estimates, the estimate of abundance
# is far better with the true model T3 (sample has been simulated from
# a population of N=500)
T1$ests ; T2$ests ; T3$ests

# A quick test of the bootstrap functionality to prove it works:
LT2D.bootstrap(T3, r=499)


# Now we create a data frame with covariates artificially introduced and see
# how the system performs. We also change the perpendicular density function
# from flat to a normal bump, so that our sample may be more representative of
# the responsivement movement data that the package is designed to deal with:
set.seed(314)
simDat1 = simXY(100, 'pi.norm', logphi, 'h1', c(1,-0.55), w, ystart)$locs
simDat2 = simXY(100, 'pi.norm', logphi, 'h1', c(0.7,-0.55), w, ystart)$locs

n1 = length(simDat1$x)
n2 = length(simDat2$x)

# We can consider fakeFactor to be a factor covariate representing species, 
# where the second species has a far higher detectability.
sim.df2 = data.frame(x = c(simDat1$x, simDat2$x),
                     y = c(simDat1$y, simDat2$y),
                     stratum = rep(1,n1+n2),
                     transect = rep(1, n1+n2),
                     L = Lsim,
                     area = Asim,
                     object = 1:(n1+n2),
                     size = rep(1,n1+n2),
                     fakeFactor = factor(rep(c(1,2),c(n1,n2))))

sim.df2.1 = data.frame(x = simDat1$x,
                       y = simDat1$y,
                       stratum = rep(1,n1),
                       transect = rep(1, n1),
                       L = Lsim,
                       area = Asim,
                       object = 1:n1,
                       size = rep(1,n1))

sim.df2.2 = data.frame(x = simDat2$x,
                       y = simDat2$y,
                       stratum = rep(1,n2),
                       transect = rep(1, n2),
                       L = Lsim,
                       area = Asim,
                       object = 1:n2,
                       size = rep(1,n2))

T4 = LT2D.fit(DataFrameInput = sim.df2.1, hr = 'h1', b=c(1,-0.55), ystart=ystart,
               pi.x='pi.norm', logphi=logphi, w=w)

T5 = LT2D.fit(DataFrameInput = sim.df2.2, hr = 'h1', b=c(1,-0.55), ystart=ystart,
               pi.x='pi.norm', logphi=logphi, w=w)

# with covariates:
T6 = LT2D.fit(DataFrameInput = sim.df2, hr = 'h1', b=c(1,-0.55), ystart=ystart,
               pi.x='pi.norm', logphi=logphi, w=w,
               formulas = list(formula(i~fakeFactor)), 
               ipars = c(0))

T7 = LT2D.fit(DataFrameInput = sim.df2, hr = 'h1', b=c(1,-0.55), ystart=ystart,
              pi.x='pi.norm', logphi=logphi, w=w)

# trying to get abundance without taking into account the added covar:
B11 = LT2D.fit(DataFrameInput = sim.df2, hr = 'h1', b=c(1,-0.55), ystart=ystart,
               pi.x='pi.norm', logphi=logphi, w=w)