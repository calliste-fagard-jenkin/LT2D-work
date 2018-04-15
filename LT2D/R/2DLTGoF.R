#' @title Kolmogarov-Smirnov goodness-of-fit p-value.
#'
#' @description
#' Kolmogarov-Smirnov goodness-of-fit p-value calculation.
# 
#' @param x value of Kolmogarov-distributed random variable at which to evaluate.
#' @param inf approximation to infinity (a large number).
#' @param dp approximation convergence criterion.
#' @return 1 - p-value for Kolmogarov distribution
#' @details
#' Calculates p-value for Kolmogarov distribution at x, approximating infite sum
#' \code{sqrt(2*pi)/x*sum{i=1}^infty exp(-(2i-1)^2*pi^2/(8x^2)))}
#' by a finite sum to inf (default 1000) if sum to inf and inf+1 differ by less
#' than dp (default 1e-4), else sum until difference is less than dp.
#' @export
p.kolmogarov=function(x,inf=1000,dp=1e-4)
{
  infsum=rep(0,inf)
  i=1:inf
  K=sqrt(2*pi)/x
  p=p1=K*sum(exp(-(2*i-1)^2*pi^2/(8*x^2)))
  dp=1
  while(dp>1e-4) {
    inf=inf+1
    p=p1+K*exp(-(2*inf-1)^2*pi^2/(8*x^2))
    dp=abs(p-p1)
  }
  return(p=p)
}


#' @title Goodness-of-fit for LT2D models.
#'
#' @description
#' Calculates goodness-of-fit in perpendicular and forward dimensions,
#' plots fit, and returns p-values. Returns two p-values for each dimension:
#' the Kolmogorov-Smirnov, and the Cramer von Mises. 
# 
#'@param LT2D.fit.modek, as output by \code{\link{fit.yx}}
#'@param plot If TRUE, does Q-Q plot. Point corresponding to largest difference
#'between empirical and theoretical cdf (on which the Kolmogarov-Smirnov test is
#'based) is circled in red.
#'@return
#'list containing vectors 'x' and 'y', each containing their 2 p-values 
#'@seealso \code{\link{fityx}} \code{\link{p.kolmogarov}}
#'@export
#' @title Goodness-of-fit for LT2D models.
#'
#' @description
#' Calculates goodness-of-fit in perpendicular and forward dimensions,
#' plots fit, and returns p-values. Returns two p-values for each dimension:
#' the Kolmogorov-Smirnov, and the Cramer von Mises. 
# 
#'@param LT2D.fit.modek, as output by \code{\link{fit.yx}}
#'@param plot If TRUE, does Q-Q plot. Point corresponding to largest difference
#'between empirical and theoretical cdf (on which the Kolmogarov-Smirnov test is
#'based) is circled in red.
#'@return
#'list containing vectors 'x' and 'y', each containing their 2 p-values 
#'@seealso \code{\link{fityx}} \code{\link{p.kolmogarov}}
#'@export
gof.LT2D = function(fit, plot=FALSE){
  require('goftest')      # used for the CvM test
  
  # First, we verify that the input data is an LT2D model:
  CLASS = class(fit)=="LT2D.fit.function.object"
  message = 'Input data type invalid, expecting an LT2D.model.fit object'
  if (class(fit)!="LT2D.fit.object" & CLASS==FALSE){
    stop(message)}
  
  if (CLASS==TRUE){fit = fit$fit} # Extract the correct 
  # fit object from the LT2D.fit function output 
  
  # Now we extract the information in this model which we require for the GoF:
  ystart = fit$ystart
  w = fit$w  
  hrname = fit$hr
  piname = fit$pi.x
  logphi=fit$logphi
  FittedData = data.with.b.conversion(fit)
  x=FittedData$x
  y=FittedData$y
  B=FittedData$beta
  n=length(x)
  
  # GoF tests for Y direction: 
  
  # Calculating the CDF and empirical CDF for the Y data:
  
  if (fit$covariates==T){
    # If we have covariates, be sure to calculate the value which is 
    # appropriate for the fitted values of the detection function for each 
    # coordinate: 
    Fy = rep(NA, n)
    F0 = rep(NA, n)
    for (i in 1:n){
      Fy[i]=(1-Sy(x=x[i],y=y[i],ymax=ystart,b=as.list(B[[i]]),hr=hrname))
      F0[i]=(1-Sy(x=x[i],y=0,ymax=ystart,b=as.list(B[[i]]),hr=hrname))
    }
  }
  
  else{
    # If there are no covariates (and hence the beta values are the same for
    # each coordinate), we can save a little computation time by doing things
    # this way instead:
    Fy=(1-Sy(x=x,y=y,ymax=ystart,b=as.list(B[[1]]),hr=hrname))
    F0=(1-Sy(x=x,y=rep(0,n),ymax=ystart,b=as.list(B[[1]]),hr=hrname))
  }
  
  Fy0=Fy/F0
  Fy0.order=order(Fy0)
  
  yy=y[Fy0.order]
  xx=x[Fy0.order]
  
  cdf.y=Fy0[Fy0.order]
  e.cdf.y=order(cdf.y)/n
  
  # Calculating the statistics:
  dF.y=cdf.y-e.cdf.y
  Dn=max(abs(dF.y))*sqrt(n)            # K-S test stat
  
  p.cvm.y=goftest::cvm.test(Fy0)$p.value
  p.ks.y=1-p.kolmogarov(Dn)  # Working out p-value from K-S test stat
  
  
  # GoF tests for X direction: 
  edf=(1:n)/n
  
  # We need a for loop to calculate F0v with or without covariate inclusion,
  # so we do this now regardless:
  f0V=vector(mode = 'numeric',length=length(x))
  for(i in 1:n){
    f0V[i]=integrate(p.pi.x,0,x[i],as.list(B[[i]]),hrname,
                     ystart,piname,logphi,w)$value
  }
  
  if (fit$covariates==T){
    Af0=rep(NA, n)
    for (i in 1:n){
      Af0[i]=integrate(p.pi.x,0,w,as.list(B[[i]]),hrname,
                       ystart,piname,logphi,w)$value
    }
  }
  
  else{
    Af0=integrate(p.pi.x,0,w,as.list(B[[1]]),hrname,
                  ystart,piname,logphi,w)$value
  }
  
  
  cdf.x=sort(f0V/Af0)          # CDF
  #cdf.x.order=order(cdf.x)    # here, we order the cdf to be able
  #cdf.x=cdf.x[cdf.order.x]    # to make comparisons with the edf
  
  # K-S statistic
  dF.x=cdf.x-edf              # vector of differences between EDF and CDF
  Dn.x=max(abs(dF.x))*sqrt(n) # K-S stats == largest diff
  p.ks.x=1-p.kolmogarov(Dn.x) # p val from Kolmogorov distribution
  
  # Under model, cdf values are from uniform; default for cvm.test is "punif"
  p.cvm.x=goftest::cvm.test(cdf.x)$p.value
  # There is a nortest::cvm.test which tests against normality instead,
  # we want to ensure we use the correct one.
  
  
  
  if(plot){ # Plotting QQs, if the user requested we do so: 
    
    # Y fit plot:
    worst.y=which(abs(dF.y)==max(abs(dF.y)))
    plot(1-e.cdf.y,cdf.y,xlab="Empirical Distribution Function",
         ylab="Cumulative Distribution Function",
         main="Forward Dist. Q-Q Plot",
         xlim=c(0,1),ylim=c(0,1),pch="+")
    lines(c(0,1),c(1,0))
    
    mtext(paste('p-values: Cramer-von Mises=',round(p.cvm.y,2),' ; kolmogarov=',
                round(p.ks.y,2)))
    
    points(1-e.cdf.y[worst.y],cdf.y[worst.y],
           col="red") # red circle round K-S point
    
    # X fit plot:
    worst.x=which(abs(dF.x)==max(abs(dF.x)))
    plot(edf,cdf.x,pch="+",xlim=c(0,1),ylim=c(0,1),
         xlab="Empirical Distribution Function",
         ylab="Cumulative Distribution Function",main="Perp. Dist. Q-Q Plot")
    lines(c(0,1),c(0,1))
    mtext(paste('p-values: Cramer-von Mises=',round(p.cvm.x,2),' ; kolmogarov=',
                round(p.ks.x,2)))
    points(edf[worst.x],cdf.x[worst.x],col="red")
  }
  
  # Making an object which contains the four calculated p values:
  NAMES = c("K-S","CvM")
  Xpvals = c(p.ks.x, p.cvm.x) ; names(Xpvals) = NAMES
  Ypvals = c(p.ks.y, p.cvm.y) ; names(Ypvals) = NAMES
  output = list(Xpvals, Ypvals) ; names(output) = c("X","Y")
  
  return(output)
}

