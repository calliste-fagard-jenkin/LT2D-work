library('LT2D')

# SET WORKING DIRECTORY TO SOURCE FILE LOCATION

Data = read.csv('JanthoData.csv',header=T)
Data = subset(Data, !is.na(Data$PP.Distance))   # Remove NAs for testing

x = Data$PP.Distance ; y = Data$Forward.Distance
species = as.factor(Data$species)

DF = data.frame(x,y,species) ; DF=subset(DF, DF$x<0.03)

fi = i~species # set the formula 
DM = DesignMatrix(DF, fi)

x = DF$x ; y = DF$y

# b = list(c(-10.2,-2.7,-9.0,0.0), -0.14) ; logphi = c(0,-4)
# pars = list(beta=b,logphi=logphi)
# vectorp = unlist(pars)
# 
# negloglik.yx(vectorp,y,x,'h1',0.05,'pi.norm',0.03, rounded.points=0,
#              DesignMatrices = list(DM),skeleton=pars)
# 
# 
# A = optim(vectorp, fn=negloglik.yx, hessian=F, y=y,x=x,hr='h1',ystart=0.05,
#       pi.x='pi.norm',w=0.03,rounded.points=0,
#       DesignMatrices=list(DM),skeleton=pars)
# 
# A$par

b = c(1.0122061, -0.6613409)
i.parameters = c(-11.14648, -13.34820, -14.56664)
logphi = c(0.01322874, -4.80826844)
w = 0.03 ; ystart = 0.03

A = fityx(DataFrameInput = DF, b=b, hr='h1', ystart=ystart,
          pi.x='pi.norm', logphi = logphi, w=w, formulas=list(fi),
          covarPars = list(i=i.parameters), hessian=F)
A$par
A$value
A$convergence
A$AIC

b = c(A$par[1],A$par[5]) ; names(b) = NULL
i.parameters = A$par[2:4] ; names(i.parameters) = NULL
logphi = A$par[6:7]; names(logphi) = NULL

# print(b)
# print(i.parameters)
# print(logphi)

A = fityx(DataFrameInput = DF, b=b, hr='h1', ystart=ystart,
          pi.x='pi.norm', logphi = logphi, w=w, hessian=F)

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

### Add some columns to DF so that it is compatible with LT2D.fit:
Data = read.csv('JanthoData.csv',header=T)
Data <- subset(Data, Data$PP.Distance<0.03 | is.na(Data$PP.Distance) )

x = Data$PP.Distance ; y = Data$Forward.Distance
species = as.factor(Data$species)
DF <- data.frame(x,y,species)
DF$stratum <- rep(1, length(DF$x))
DF$transect <- Data$Sample.Label
DF$L <- Data$Effort
DF$area <- rep(101*2*0.03, length(DF$x))
DF <- addObjectNumbers(DF)



# Testing covariate inclusion on h1
B2 = LT2D.fit(DataFrameInput = DF, hr = 'h1', b=b, ystart = ystart, 
             pi.x = 'pi.norm', logphi = logphi, w=w, formulas = list(fi), 
             ipars = i.parameters)

# Testing covariate inclusion on h1, with wrong start parameters
B2.1 = LT2D.fit(DataFrameInput = DF, hr = 'h1', b=b, ystart = ystart, 
              pi.x = 'pi.norm', logphi = logphi, w=w, formulas = list(fi), 
              ipars = x.parameters)

# Making sure the model still works without the covariates included
B3 = LT2D.fit(DataFrameInput = DF, hr = 'h1', b=b, ystart = ystart, 
              pi.x = 'pi.norm', logphi = logphi, w=w)

# and with a different detection function:
B4 = LT2D.fit(DataFrameInput = DF, hr = 'ep1',
              b=c(10.23357070,2.38792702,-20.23029177),
              ystart = ystart, pi.x = 'pi.norm', logphi = logphi, w=w)

# the same with attempt at covariates:
B5 = LT2D.fit(DataFrameInput = DF, hr = 'ep1',
              b=c(10.23357070,2.38792702,-20.23029177),
              ystart = ystart, pi.x = 'pi.norm', logphi = logphi, w=w,
              formulas = list(fi), ipars = i.parameters)

# with the covariates in the x dimension instead:
B6 = LT2D.fit(DataFrameInput = DF, hr = 'ep1',
              b=c(10.23357070,2.38792702,-20.23029177),
              ystart = ystart, pi.x = 'pi.norm', logphi = logphi, w=w,
              formulas = list(formula(x~species)), xpars = i.parameters,
              ypars=i.parameters)

# covariates in x and i dimension at the same time:
B7 = LT2D.fit(DataFrameInput = DF, hr = 'ep1',
              b=c(10.23357070,2.38792702,-20.23029177),
              ystart = ystart, pi.x = 'pi.norm', logphi = logphi, w=w,
              formulas = list(formula(x~species), formula(i~species)),
              xpars = i.parameters, ipars=i.parameters,
              ypars = i.parameters)

set.seed(0)
simDat = simXY(500, 'pi.const', NULL, 'h1', b, w, ystart)$locs

Lsim = 20 ; Asim = Lsim*w*2

sim.df = data.frame(x = simDat$x, y = simDat$y, stratum=rep(1,length(simDat$x)),
                   transect = rep(1,length(simDat$x)), L = Lsim, area = Asim,
                   object = 1:length(simDat$x), size = rep(1, length(simDat$x)))



B8 = LT2D.fit(DataFrameInput = sim.df, hr = 'h1', b=b, ystart=ystart,
              pi.x='pi.norm', logphi=logphi, w=w)
B8$fit$par


# Try using a different perp density:
B9 = LT2D.fit(DataFrameInput = sim.df, hr = 'h1', b=b, ystart=ystart,
              pi.x='pi.const', logphi=NULL, w=w)


# Create another sample, but with a detection function with a different
# intercept parameter, so that we may create an artificial covariate
simDat2 = simXY(500, 'pi.const', NULL, 'h1', c(b[1]-1,b[2]), w, ystart)$locs

n1 = length(simDat$x)
n2 = length(simDat2$x)

sim.df2 = data.frame(x = c(simDat$x, simDat2$x),
                     y = c(simDat$y, simDat2$y),
                     stratum = rep(1,n1+n2),
                     transect = rep(1, n1+n2),
                     L = Lsim,
                     area = Asim,
                     object = 1:(n1+n2),
                     size = rep(1,n1+n2),
                     fakeFactor = factor(rep(c(1,2),c(n1,n2))))

# And with covariates:
B10 = LT2D.fit(DataFrameInput = sim.df2, hr = 'h1', b=b, ystart=ystart,
              pi.x='pi.const', logphi=NULL, w=w,
              formulas = list(formula(i~fakeFactor)), 
              ipars = c(-1))

# trying to get abundance without taking into account the added covar:
B11 = LT2D.fit(DataFrameInput = sim.df2, hr = 'h1', b=b, ystart=ystart,
               pi.x='pi.const', logphi=NULL, w=w)

# remove(Sy, px, phat, ParamNumRequired, p.pi.x, negloglik.yx, NDest, LT2D.fit,
#        LinPredictor, invp1_replacement, HazardBCheck, HazardCovarsAllowed,
#        HazardCovarSlots, HazardNumberLookup, h1, fyx, FormulaChecking, fityx,
#        ep1, DesignMatrix, DensityNumberLookup,data.with.b.conversion)


gof.LT2D(B10, plot=T)
Betas = data.with.b.conversion(B10$fit)$beta
Sy(sim.df2$x[2],0,ystart,as.list(Betas[[2]]),'h1')

# Try with an ep1 model: 

# Create another sample, but with a detection function with a different
# intercept parameter, so that we may create an artificial covariate
simDat3 = simXY(500, 'pi.const', NULL, 'ep1',
               b=c(1000,1000,1000), w, ystart)$locs

n3 = length(simDat3$x);n3 

number.seen = function(par) -length(simXY(1000, 'pi.const',NULL,'ep1',
                                         par,w,ystart)$locs)

optim(c(1.5,1.7,-20.22), number.seen, hessian=F)

sim.df3 = data.frame(x = simDat3$x,
                     y = simDat3$y,
                     stratum = rep(1,n3),
                     transect = rep(1, n3),
                     L = Lsim,
                     area = Asim,
                     object = 1:(n3),
                     size = rep(1,n3))

