library('LT2D')

setwd("C:/Users/Cal/OneDrive/covariate functions")

Data = read.csv('JanthoData.csv',header=T)
Data = subset(Data, !is.na(Data$PP.Distance))   # Remove NAs for testing

x = Data$PP.Distance ; y = Data$Forward.Distance
species = as.factor(Data$species)

DF = data.frame(x,y,species) ; DF=subset(DF, DF$x<0.03)

fi = i~species # set the formula 
DM = DesignMatrix(DF, fi)

x = DF$x ; y = DF$y

b = list(c(-10.2,-2.7,-9.0,0.0), -0.14) ; logphi = c(0,-4)
pars = list(beta=b,logphi=logphi)
vectorp = unlist(pars)

negloglik.yx(vectorp,y,x,'h1',0.05,'pi.norm',0.03, rounded.points=0,
             DesignMatrices = list(DM),skeleton=pars)


A = optim(vectorp, fn=negloglik.yx, hessian=F, y=y,x=x,hr='h1',ystart=0.05,
      pi.x='pi.norm',w=0.03,rounded.points=0,
      DesignMatrices=list(DM),skeleton=pars)

A$par

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

# b = c(A$par[1],A$par[5]) ; names(b) = NULL
# i.parameters = A$par[2:4] ; names(i.parameters) = NULL
# logphi = A$par[6:7]; names(logphi) = NULL

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
              formulas = list(formula(x~species)), xpars = i.parameters)

# covariates in x and i dimension at the same time:
B7 = LT2D.fit(DataFrameInput = DF, hr = 'ep1',
              b=c(10.23357070,2.38792702,-20.23029177),
              ystart = ystart, pi.x = 'pi.norm', logphi = logphi, w=w,
              formulas = list(formula(x~species),formula(i~species)),
              xpars = i.parameters) ### Hasn't attempted to do the i formula??
