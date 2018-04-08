# An attempt to replicate results obtained in the original LT2D paper, to prove
# that the software still works as before, with covariate capabilities
# included. 

library(LT2D)              # We load the original package as of August 2017
source('Covar stuff.R')    # and add in the changes made for the project

# NOTE : Loading the LT2D package also loads primate.dat and dolphin.dat into
#        the global environment, since these data sets are included within the 
#        package.

# Set w and ystart to paper values:
w = 0.03 ; ystart = 0.05
L = 10 ; A = 2*w*L
n.primate = length(primate.dat$x)

# Primate Data Selected Model
primate.df = data.frame(x = primate.dat$x,
                        y = primate.dat$y,
                        stratum = rep(1,n.primate),
                        transect = rep(1,n.primate),
                        L = L,
                        area = A,
                        object = 1:n.primate,
                        size = rep(1,n.primate))


b=c(5.2919208, 8.4701307)
logphi=c(0.01784102, -4.42209067)

Primate.Fit = LT2D.fit(DataFrameInput =  primate.df,
                       hr = 'ip0',
                       b = b,
                       ystart = ystart,
                       pi.x = 'pi.norm',
                       logphi = logphi,
                       w = w,
                       hessian = T)

gof.LT2D(Primate.Fit) # We confirm that these are the same as the paper's

# out of curiosity, we now also look at the estimate of abundance using the
# method implemented in July 2017:
Primate.Fit$ests


# Change w and ystart global variables to the appropriate values for the 
# dolphin data:
w = 0.15 ; ystart = 0.55
n.dolphin = length(dolphin.dat$x)

# Create a data frame for the dolphin data:
dolphin.df = data.frame(x = dolphin.dat$x,
                        y = dolphin.dat$y,
                        stratum = rep(1,n.dolphin),
                        transect = rep(1,n.dolphin),
                        L = L,
                        area = A,
                        object = 1:n.dolphin,
                        size = rep(1,n.dolphin))

# Dolphin Data Selected Model
b=c(-7.3287948, 0.9945317)
logphi=-0.4811025

Dolphin.Fit = LT2D.fit(DataFrameInput =  dolphin.df,
                       hr = 'h1',
                       b = b,
                       ystart = ystart,
                       pi.x = 'pi.hnorm',
                       logphi = logphi,
                       w = w,
                       hessian = , rmin=-0.1)

gof.LT2D(Dolphin.Fit)
#EHSW:
phatInterval(dfit.hn)
phatInterval(dfit.hn)*wd
# p(0):
p0.hn=1-Sy(0,0,ystartd,dfit.hn$b,h1);p0.hn
plotfit.smoothfy(dfit.hn,xmax=0.01)
# Density estimate:
n=length(dfit.hn$dat$x)
L=1672.77 # from Canadas et al. (in nm)
Dhat=(n/phatInterval(dfit.hn)$phat)/(2*wd*L)
Dhat
