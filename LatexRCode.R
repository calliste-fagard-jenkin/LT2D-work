library('LT2D')

# set perp trunc, forward trunc
# total transect length and survey
# area:
w = 0.03 ; ystart = 0.05
L = 10 ; A = 2*w*L

# set value of 'true' parameters
# for simulated data:
b <- c(-4.04,0.79)
logphi <- c(0.02,-4.42)

# produce simulated data:
set.seed(3)
simDat = simXY(500, 'pi.norm',
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
                b = c(-3.72,0.74),
                ystart = ystart,
                pi.x = 'pi.norm',
                # start values for logphi:
                logphi = c(0.02, -4.42),
                w = w,
                hessian = TRUE)

gof.LT2D(fit)
LT2D.bootstrap(fit,999)$ci
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
                 ipars=0,
                 xpars=0,
                 ypars=0)