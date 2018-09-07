# This document serves as a simulation study on the use of 2 part mixture models
# for the perpendicular densities of LT2D analyses

# Load the LT2D package, this also loads parallelization libraries:
library(LT2D)

produce.mixt.df <- function(N, pi.x, logphi1, logphi2, hr, b, w, ystart, lambda,
                            Lsim, Asim){
  # purpose : Produces a smiluated data frame from a population with chosen 
  #           parameters. Currently not designed to work in conjuction with
  #           covariate inclusion.
  #
  # output  : A data.frame in the style required by the LT2D software to 
  #           estimate MLEs and produce an estimate of abundance
  
  # produce the simulated data from each component:
  comp1 <- simXY(floor(N*lambda), pi.x, logphi1, hr, b, w, ystart)$locs
  comp2 <- simXY(floor(N*(1-lambda)), pi.x, logphi2, hr, b, w, ystart)$locs
  
  # get the number of observations:
  n <- (length(c(comp1$x,comp2$x)))
  
  # produce the LT2D format data.frame:
  mixt.df <- data.frame(x = c(comp1$x, comp2$x),
                        y = c(comp1$y, comp2$y),
                        stratum = 1,
                        transect = 1,
                        area = Asim,
                        L = Lsim,
                        object = 1:n,
                        size = 1)
  return(mixt.df)
}

simulation <- function(R, N, pi.x, logphi1, logphi2, hr, b, w, ystart, lambda,
                       Lsim, Asim, produce.mixt.df){
  # purpose : Produces R simulated data sets with the chosen parameters,
  #           then fits the data using a mixture and regular LT2D model, to 
  #           produce estimates of bias
  #
  # output  : A list, with elements mixt.bias, reg.bias, mixt.ests and reg.ests
  
  # Create the clusters and initiate paralellisation:
  cl<-makeCluster(detectCores())
  registerDoParallel(cl)

  # Do the parallel loop:
  ls <- foreach(i=1:R, .packages = 'LT2D') %dopar% {

    # 1. Get the data:
    df <- produce.mixt.df(N, pi.x, logphi1, logphi2, hr, b, w, ystart, lambda,
                          Lsim, Asim)
    
    # 2. Fit a normal model:
    F1 <- try(LT2D.fit(DataFrameInput = df, hr = hr, b=b, ystart=ystart,
                       pi.x=pi.x, logphi=logphi1, w=w), silent = T)

    # 3. Fit a mixture model:
    F2 <- try(LT2D.mixture(df, hr, b, ystart, pi.x, logphi1, logphi2, w,
                           qlogis(lambda)), silent = T)
    
    # 4. If either model failed, we can't use the information:
    if (class(F1)=='try-error' | class(F2)=='try-error'){
      to.ls <- c(NA, NA)
    }
    
    # 5. Otherwise, update the information
    else{
      N.mixt <- F1$ests$N ; N.mixt <- N.mixt[length(N.mixt)]
      N.reg <- F2$ests$N ; N.reg <- N.reg[length(N.reg)]
      
      to.ls <- c(N.mixt, N.reg)
    }
  }
  
  # Silently close the ports:
  try(stopCluster(cl), silent = T)
  
  # 6. Calculate bias and return info
  ls <- matrix(unlist(ls),nrow=2)
  mixt.Ns <- na.omit(ls[1,])
  regu.Ns <- na.omit(ls[2,])
  N <- floor(N*lambda) + floor(N*(1-lambda))
  
  mixt.bias <- 100*(mean(mixt.Ns)-N)/N
  regu.bias <- 100*(mean(regu.Ns)-N)/N
  
  output <- list(mixt.bias=mixt.bias,
                 regu.bias=regu.bias,
                 mixt.ests=mixt.Ns,
                 regu.ests=regu.Ns)
  return(output)
}


# Select simulation parameters:
N <- 1000
pi.x <- 'pi.norm'
hr <- 'h1'
lambda <- 0.6
logphi1 <- c(0,-3)
logphi2 <- c(0,-4)
b <- c(1,-0.5)
w <- 0.05
ystart <- 0.1
Lsim <- 20
Asim <- Lsim*w*2


set.seed(962018)
s1 <- simulation(500, N, pi.x, logphi1, logphi2, hr, b, w, ystart, lambda,
                 Lsim, Asim, produce.mixt.df)

set.seed(972018)
s2 <- simulation(500, N, pi.x, logphi1, c(0.5,-4), hr, b, w, ystart, 0.8,
                 Lsim, Asim, produce.mixt.df)

sum(abs(s2$mixt.ests-1000)<abs(s2$regu.ests-1000))
sum(abs(s1$mixt.ests-1000)<abs(s1$regu.ests-1000))

