#'@title Detection hazard function \code{h1} prob(detect | available at x,y)
#'
#'@description  This hazard function has the form k(r,y)=a*r^(-b) from Hayes and Buckland (1983)
#' p36. Note: This function uses x for perp. dist., they use y.
#'
#'@references Hayes, R. J., and S. T. Buckland. "Radial-distance models for the line-transect method." Biometrics (1983): 29-42.
#'@param y Forward distance
#'@param x perpendicular distance
#'@param b parameter vector, where \code{b[2]} is log(theta), and the function returns
#'theta[1]*(y^2+x^2)^(-theta[2]/2).
#'@return probability of detection given that an animal is availabe at location x,y
#'@examples
#'h1(0.5,0.5,b=log(c(0.001,1)))
#'@seealso \code{\link{h2}}, \code{\link{ghy}}, \code{\link{ghy2}}
#'@export
h1=function(y,x,b)
{

  # Check parameters have been passed as list
  #if (class(b)!='list'){stop(
  #  'paramaters must be given as a list')}

  # Check correct number of parameters have been passed
  if(length(b)!=2) {
    cat(b,"\n")
    stop("b must be vector of length 2.")
  }

  # Check x and y dimensions, and that parameter vectors
  # are of the same dimension
  theta1 = exp(b[1]) ; theta2 = exp(b[2])
  l1 = length(theta1) ; l2 = length(theta2)
  lx = length(x) ; ly = length(y)

  #if (lx!=ly){stop('x and y lengths differ')}

  ErrorMessage = 'Covariate inclusion dimension error'

  if (l1>1){            # If covariates in intercept...
    if(l1!=lx){stop(ErrorMessage)}
  }

  if (l2>1){            # If covariates in theta2...
    if(l2!=lx){stop(ErrorMessage)}
  }

  return(theta1*(y^2+x^2)^(-theta2/2))
}

#'@title Detection hazard function \code{ghy} prob(detect | available at x,y)
#'
#'@description  This hazard function is a generalization of the form k(r,y)=a*r^(-b) from
#'Hayes and Buckland (1983) p36, the generalization being that a parameter to be estimated
#'is added to y. When this parameter is zero you get the form of Hayes and Buckland (1983) p36.
#'Note: This function uses x for perp. dist., they use y.
#'
#'@references Hayes, R. J., and S. T. Buckland. "Radial-distance models for the line-transect method." Biometrics (1983): 29-42.
#'@param y Forward distance
#'@param x perpendicular distance
#'@param b parameter vector, where \code{b[2]} is log(theta); theta[1] is as per \link{\code{h1}}
#'parameter theta[1]; theta[2] is as per \link{\code{h1}} theta[2]; theta[3] is parameter that
#'is added to y to shift forward distance origin and allow p(0)<1.
#'@return probability of detection given that an animal is availabe at location x,y
#'@examples
#'h1(0.5,0.5,b=log(c(0.001,1)))
#'ghy(0.5,0.5,b=log(c(0.001,1,0.0)))
#'ghy(0.5,0.5,b=log(c(0.001,1,0.01)))
#'@seealso \code{\link{h2}}, \code{\link{ghy2}}
#'@export
ghy=function(y,x,b)
{
  if(length(b)!=3) {
    cat(b,"\n")
    stop("b must be vector of length 3.")
  }
  theta=exp(b)
  #theta=c(exp(b[1:2]),b[3])
  theta1=theta[1]^(1/theta[2])
  return(((x/theta1)^2+((y+theta[3])/theta1)^2)^(-theta[2]/2))
}

#'@title Detection hazard function \code{ghy2} prob(detect | available at x,y)
#'
#'@description  This hazard function is a generalization of the form k(r,y)=a*r^(-b) from
#'Hayes and Buckland (1983) p36, the generalization being that (1) a parameter to be estimated
#'is added to y, and (2) x and y have separate scale parameters. It is a generalization of
#'\link{\code{ghy}} to allow x and y to have separate scale parameters.
#'Note: This function uses x for perp. dist., they use y.
#'
#'@references Hayes, R. J., and S. T. Buckland. "Radial-distance models for the line-transect method." Biometrics (1983): 29-42.
#'@param y Forward distance
#'@param x perpendicular distance
#'@param b parameter vector, where \code{b[2]} is log(theta); theta[1] is as per \link{\code{h1}}
#'parameter theta[1]; theta[2] is as per \link{\code{h1}} theta[2]; theta[3] is parameter that
#'is added to y to shift forward distance origin and allow p(0)<1; theta[4] is the equivalent of
#'theta[1], but specific to y, whereas theta[1] is specific to x in this function.
#'@return probability of detection given that an animal is availabe at location x,y
#'@examples
#'h1(0.5,0.5,b=log(c(0.001,1)))
#'ghy(0.5,0.5,b=log(c(0.001,1,0.0)))
#'ghy(0.5,0.5,b=log(c(0.001,1,0.01)))
#'ghy2(0.5,0.5,b=log(c(0.001,1,0.01,0.001)))
#'ghy2(0.5,0.5,b=log(c(0.001,1,0.01,0.005)))
#'@seealso \code{\link{h2}}, \code{\link{ghy2}}
#'@export
ghy2=function(y,x,b)
{
  if(length(b)!=4) {
    cat(b,"\n")
    stop("b must be vector of length 4.")
  }
  theta=exp(b)
  thetax=theta[1]^(1/theta[2])
  thetay=theta[4]^(1/theta[2])
  return(((x/thetax)^2+((y+theta[3])/thetay)^2)^(-theta[2]/2))
}

#'@title Detection hazard function \code{h2} prob(detect | available at x,y)
#'
#'@description  This hazard function has the form k(r,y)=a*sqrt(r^2-y^2)/r^(b+1); b>2 from Hayes and Buckland (1983) p37.
#'
#'@references Hayes, R. J., and S. T. Buckland. "Radial-distance models for the line-transect method." Biometrics (1983): 29-42.
#'@param y Forward distance
#'@param x perpendicular distance
#'@param b parameter vector, where \code{b[2]} is log(theta)
#'@return probability of detection given that an animal is availabe at location x,y
#'@examples
#'h2(0.5,0.5,b=log(c(0.75,1)))
#'@export
#'@seealso \code{\link{h1}}
h2=function(y,x,b)
  #-------------------------------------------------------------------------------
# Detection hazard function prob(detect | available at x,y),
# Corresponding to Hayes and Buckland (1983) k(r,y)=a*sqrt(r^2-y^2)/r^(b+1); b>2
# on p37.
# Note: I use x for perp. dist., they use y.
# Inputs:
#  b: log(theta), where theta is vector of hazard rate parameters
#-------------------------------------------------------------------------------
{
  if(length(b)!=2) {
    cat(b,"\n")
    stop("b must be vector of length 2.")
  }
  theta=exp(b)
  return(theta[1]*y*(y^2+x^2)^(-(theta[2]+3)/2))
}

#' @title Hazard detection function of form \code{h2} with g(0)<1
#'
#'@description  This hazard function has the form \eqn{k(r,y)=c*[a*\sqrt(r^2-y^2)/r^{(b+1)}]; b>2} modified
#'from Hayes and Buckland (1983) p37, where c is imperfect detectability i.e. \eqn{g(0)<1}.
#'
#'@references Hayes, R. J., and S. T. Buckland. "Radial-distance models for the line-transect method." Biometrics (1983): 29-42.
#'@param y Forward distance
#'@param x perpendicular distance
#'@param b parameter vector, where \code{b[2]} is log(theta) and \code{b[3]} is qlogis(b[3]) detectability at the observer i.e. \eqn{g(0)=c}
#'@return probability of detection given that an animal is availabe at location x,y
#'#'@examples
#'h21(0.5,0.5,b=c(log(c(0.75,1)),qlogis(0.9)))
#' @export
h21=function(y,x,b)
  #-------------------------------------------------------------------------------
# Detection hazard function prob(detect | available at x,y),
# Corresponding to Hayes and Buckland (1983) k(r,y)=a*sqrt(r^2-y^2)/r^(b+1); b>2
# on p37.
# Note: I use x for perp. dist., they use y.
# Inputs:
#  b: log(theta), where theta is vector of hazard rate parameters
#-------------------------------------------------------------------------------
{
  if(length(b)!=3) {
    cat(b,"\n")
    stop("b must be vector of length 3.")
  }
  theta=exp(b[1:2])
  dF=function(y,theta) theta[1]*y*(y^2+x^2)^(-(theta[2]+3)/2)
  g0=plogis(b[3])
  return(dF(y,theta)*g0)
}

#' @title Hayes+Buckland hazard rate model.
#' @description
#' Evaluates the Hayes and Buckland (1983) hazard rate model (HB model) as parameterised on
#' page 36 of that paper.
#' @param x perpendicular distance(s).
#' @param theta parameter vector with theta[1] being HB model parameter a1 and theta[2] being b.
#' @return
#' g(x)=1-exp(-a1*x^(-(b-1))).
#' @seealso \code{\link{h1.to.HB}}
#' @examples
#' xx=seq(0,0.03,length=100)
#' p=HBhr(xx,h1.to.HB(c(-7, 0.85)))
#' plot(xx,p,type="l",xlab="Perpendicular distance (x)",ylab="p(x)",ylim=c(0,1))
HBhr=function(x,theta) 1-exp(-theta[1]*x^(-(theta[2]-1)))

#' @title Parameter conversion for Hayes+Buckland hazard rate model.
#'
#' @description
#' Converts parameters from the 2-dimensional a*r^{-b} form of the Hayes and Buckland (1983)
#' hazard rate model used in \link{\code{h1}}, to parameters of the perpendicular distance form
#' of the model used in conventional distance sampling.
#
#' @param b parameter vector of the \link{\code{h1}}.
#' @return
#' parameter vector c(a1,b1) for Hayes and Buckland (1983) hazard rate model
#' g(x)=1-exp(-a1*x^(-(b1-1))).
#' @seealso \code{\link{h1}}
#' @examples
#' h1.to.HB(c(-7, 1))
#' @export
h1.to.HB=function(b){
  if(length(b)!=2) {
    cat(b,"\n")
    stop("b must be vector of length 2.")
  }
  b1=exp(b[2])
  a=exp(b[1])
  a1=a*gamma((b1-1)/2)*gamma(0.5)/(2*gamma(b1/2))
  return(c(a1,b1))
}

#'@title Three-parameter inverse power hazard detection function
#'
#'@description  Inverse power hazard function, as per Borchers and Langrock (in press):
#'Has form h(y,x)=theta[1]*[theta[2]/(sqrt{theta[2]^2+x^2+y^2})]^(theta[3]+1).
#'
#'@references Borchers, D.L and Langrock, R."Double-observer line transect surveys with Markov-
#'modulated Poisson process models for animal availability" Biometrics (in press).
#'@param y Forward distance
#'@param x perpendicular distance
#'@param b parameter vector, where \code{b[1]} is plogis(theta[1])  \code{b[2]} is
#'log(theta[2]) and \code{b[3]} is log(b[3]).
#'@return probability of detection given that an animal is availabe at location x,y
#'@examples
#'b=c(-23.725809, -3.136638,2.122910)
#'ip1(0.5,0.5,b=b)
#'yy=seq(0,0.03,length=100);xx=rep(0,100)
#'hh=ip1(yy,xx,b=b)
#'plot(yy,hh,type="l")
#' @export
ip1=function(y,x,b)
{
  if(length(b)!=3) {
    cat(b,"\n")
    stop("b must be vector of length 3.")
  }
#  theta=exp(b[2:3])
#  dF=function(y,x,theta) (theta[1]/(theta[1]^2+x^2+y^2))^(theta[2]+1)
#  g0=plogis(b[1])
  theta=exp(b)
  p=theta[1]*(1/sqrt(1+(x/theta[2])^2+(y/theta[2])^2))^(theta[3]+1)
#  p=theta[1]*(theta[2]/sqrt(theta[2]^2+x^2+y^2))^(theta[3]+1)
  return(p)
}


#'@title Inverse power hazard detection function
#'
#'@description  Inverse power hazard function, as per Borchers and Langrock (in press):
#'Has form h(y,x)=theta[1]*(1/sqrt(1+(x)^2+(y)^2))^(theta[2]+1).
#'
#'@references Borchers, D.L and Langrock, R."Double-observer line transect surveys with Markov-
#'modulated Poisson process models for animal availability" Biometrics (in press).
#'@param y Forward distance
#'@param x perpendicular distance
#'@param b 2-parameter vector, where \code{b[1]} is log(theta[1]) and  \code{b[2]} is
#'log(theta[2]).
#'@return probability of detection given that an animal is availabe at location x,y
#'@examples
#'b=c(5.2919208, 8.4701307)
#'ip0(0.05,0.05,b=b)
#'yy=seq(0,0.03,length=100);xx=rep(0,100)
#'hh=ip0(yy,xx,b=b)
#'plot(yy,hh,type="l")
#' @export
ip0=function(y,x,b)
{
  if(length(b)!=2) {
    cat(b,"\n")
    stop("b must be vector of length 2.")
  }
  theta=exp(b)
  p=theta[1]*(1/sqrt(1+(x)^2+(y)^2))^(theta[2]+1)
  return(p)
}


#'@title Three-parameter exponential power hazard detection function
#'
#'@description  Inverse power hazard function, as per Borchers and Langrock (in press):
#'Has form h(y,x)=theta[1]*exp(-(x^theta[3]+y^theta[3])/(theta[2]^theta[3])).
#'
#'@references Borchers, D.L and Langrock, R."Double-observer line transect surveys with Markov-
#'modulated Poisson process models for animal availability" Biometrics (in press).
#'@param y Forward distance
#'@param x perpendicular distance
#'@param b parameter vector, where \code{b[1]} is plogis(theta[1])  \code{b[2]} is
#'log(theta[2]) and \code{b[3]} is log(b[3]).
#'@return probability of detection given that an animal is availabe at location x,y
#'#'@examples
#'b=c(1, -4, 1)
#'ep1(0.5,0.5,b=b)
#'yy=seq(0,0.03,length=100);xx=rep(0,100)
#'hh=ep1(yy,xx,b=b)
#'plot(yy,hh,type="l")
#' @export
ep1=function(y,x,b)
{
  if(length(b)!=3) {
    cat(b,"\n")
    stop("b must be vector of length 3.")
  }
  theta=exp(b[2:3])
  dF=function(y,x,theta) exp(-(x^theta[2]+y^theta[2])/(theta[1]^theta[2]))
  g0=plogis(b[1])
  return(dF(y,x,theta)*g0)
}


#'@title Four-parameter inverse power hazard detection function
#'
#'@description  Inverse power hazard function, as per Borchers and Langrock (in press):
#'Has form h(y,x)=theta[1]*(1/(1+(x/theta[2])^2+(y/theta[4])^2))^(theta[3]+1).
#'
#'@references Borchers, D.L and Langrock, R."Double-observer line transect surveys with Markov-
#'modulated Poisson process models for animal availability" Biometrics (in press).
#'@param y Forward distance
#'@param x perpendicular distance
#'@param b parameter vector, where \code{b[1]} is plogis(theta[1]) \code{b[2]} is
#'log(theta[2]), where \code{theta[2]} is the scale parameter for x, \code{b[3]} is log(theta[3])
#'and \code{b[2]} is log(theta[4]), where \code{theta[4]} is the scale parameter for y.
#'@return probability of detection given that an animal is availabe at location x,y
#'#'@examples
#'b=c(-23.725809,-3.136638,2.122910,-3.136638)
#'ip2(0.5,0.5,b=b)
#'yy=seq(0,0.03,length=100);xx=rep(0,100)
#'hh=ip2(yy,xx,b=b)
#'plot(yy,hh,type="l")
#' @export
ip2=function(y,x,b)
{
  if(length(b)!=4) {
    cat(b,"\n")
    stop("b must be vector of length 3.")
  }
#  theta=exp(b[2:4])
#  dF=function(y,x,theta) (1/(1+(x/theta[1])^2+(y/theta[3])^2))^(theta[2]+1)
#  g0=plogis(b[1])
#  return(dF(y,x,theta)*g0)
  theta=exp(b)
  p=theta[1]*(1/sqrt(1+(x/theta[2])^2+(y/theta[4])^2))^(theta[3]+1)
  return(p)
}



#'@title Four-parameter exponential power hazard detection function
#'
#'@description  Inverse power hazard function, as per Borchers and Langrock (in press):
#'Has form h(y,x)=theta[1]*exp(-(x^theta[3]+y^theta[3])/(theta[2]^theta[3])).
#'
#'@references Borchers, D.L and Langrock, R."Double-observer line transect surveys with Markov-
#'modulated Poisson process models for animal availability" Biometrics (in press).
#'@param y Forward distance
#'@param x perpendicular distance
#'@param b parameter vector, where \code{b[1]} is plogis(theta[1]) \code{b[2]} is
#'log(theta[2]), where \code{theta[2]} is the scale parameter for x, \code{b[3]} is log(theta[3])
#'and \code{b[2]} is log(theta[4]), where \code{theta[4]} is the scale parameter for y.
#'@return probability of detection given that an animal is availabe at location x,y
#'#'@examples
#'b=c(1, -4, 1)
#'ep1(0.5,0.5,b=b)
#'yy=seq(0,0.03,length=100);xx=rep(0,100)
#'hh=ep1(yy,xx,b=b)
#'plot(yy,hh,type="l")
#' @export
ep2=function(y,x,b)
{
  if(length(b)!=4) {
    cat(b,"\n")
    stop("b must be vector of length 3.")
  }
  theta=exp(b[2:4])
  dF=function(y,x,theta) exp(-((x/theta[1])^theta[2]+(y/theta[3])^theta[2]))
  g0=plogis(b[1])
  return(dF(y,x,theta)*g0)
}


#'@title Detection hazard function \code{h.exp2} prob(detect | available at x,y)
#'
#'@description  2-paramter Exponential power hazard model of Skaug & Schweder 1999. The gamma parameter is fixed at 2.

#'@references Skaug, Hans J., and Tore Schweder. "Hazard models for line transect surveys with independent observers." Biometrics 55.1 (1999): 29-36.
#'@param y Forward distance
#'@param x perpendicular distance
#'@param b parameter vector, may be logged parameter values
#'@return probability of detection given that an animal is availabe at location x,y
#'@examples
#'h.exp2(0.5,0.5,b=log(c(0.75,0.9)))
#'@seealso \code{\link{h1}} \code{\link{h2}}
#'@export
h.exp2=function(y,x,b=c(0,0))
  #----------------------------------------------------------
# Detection hazard function prob(detect | available at x,y),
# 2-paramter Exponential power hazard model of Skaug & Schweder 1999.
# (gama fixed equal to 2).
#----------------------------------------------------------
{
  mu=exp(b[1])
  sigma=exp(b[2])
  gama=2
  hr=mu*exp(-(x^gama + y^gama)/sigma^gama)
  return(hr)
}

#'@title Detection hazard function of Hiroshi et al. (2003)
#'#'
#'@description  2-parameter hazard model of Hiroshi, et al. 2003.
#'
#'@param y Forward distance
#'@param x perpendicular distance
#'@param b parameter vector, may be logged parameter values
#'@references Okamura, Hiroshi et al. (2003) Abundance Estimation of Diving Animals by the Double-Platoform Line Transect Method, Biometrics 59(3):512-520.
#'@return probability of detection given that an animal is available at location x,y
#'@examples
#'h.okamura(0.5,0.5,b=log(c(0.5,0.9)))
#'@seealso \code{\link{h1}} \code{\link{h2}} \code{\link{h.exp2}}
#'@export
h.okamura=function(y,x,b=c(0,0))
  #----------------------------------------------------------
# Detection hazard function prob(detect | available at x,y).
# From Okamura's paper
#----------------------------------------------------------
{
  if(length(b)!=2) {
    cat(b,"\n")
    stop("b must be vector of length 2.")
  }
  sigma.x=exp(b[1])
  sigma.y=exp(b[2])
  return(exp(-(x/sigma.x + y/sigma.y)))
}


#'@title Detection hazard function \code{h.const} prob(detect | available at x,y)
#'
#'@description  This function is for constant detectability throughout
#'
#'@param y Forward distance
#'@param x perpendicular distance
#'@param b single value parameter vector (giving probability of detection).
#'@return probability of detection given that an animal is available at location x,y
#'@examples
#'h.const(0.5,0.5,b=1)
#'@export
#'@seealso \code{\link{h1}} \code{\link{h2}} \code{\link{h.exp2}} \code{\link{h.okamura}}
h.const=function(y,x,b=1) {fName='h.const'
  return(rep(b[1],length(y)))}

#'@title Half-normal form for perpendicular animal density function
#'
#'@description Half-normal distribution of perpendicular animal density.
#'Truncation occurs at x=w (the perpendicular truncation distance)
#'
#'@param x prependicular trackline distance
#'@param logphi numeric vector; parameters, some of which may be logged
#'@param w perpendicular truncation distance
#'@return \eqn{\pi(x)} animal density at distance x
#'@examples
#'plot(seq(0,1,length=100),pi.hnorm(x=seq(0,1,length=100),logphi=0.5,w=1),
#'type='l',xlab='Perp. distance, x',ylab=expression(pi(x)))
#'@export
pi.hnorm=function(x,logphi,w){
  hnF=function(x,logphi) exp(-x^2/(2*exp(logphi[1])^2))
  return(hnF(x,logphi)/integrate(hnF,0,w,logphi)$value)
}

#'@title Complementary half-normal form for perpendicular animal density function
#'
#'@description Complementary half-normal distribution of perpendicular animal density.
#'Truncation occurs at x=w (the perpendicular truncation distance)
#'
#'@param x prependicular trackline distance
#'@param logphi numeric vector; parameters, some of which may be logged
#'@param w perpendicular truncation distance
#'@return \eqn{\pi(x)} animal density at distance x
#'@examples
#'plot(seq(0,1,length=100),pi.hnorm(x=seq(0,1,length=100),logphi=0.5,w=1),
#'type='l',xlab='Perp. distance, x',ylab=expression(pi(x)))
#'@export
pi.chnorm=function(x,logphi,w){
  chnF=function(x,logphi) 1-exp(-(x-logphi[1])^2/(2*exp(logphi[2])^2))
  return(chnF(x,logphi)/integrate(chnF,0,w,logphi)$value)
}


#'@title Truncated normal form for perpendicular animal density function
#'
#'@description Truncated normal distribution of perpendicular animal density.
#'Truncation occurs at x=0 (on the track line) and x=w (perpendicular truncation distance)
#'
#'@param x prependicular trackline distance
#'@param logphi numeric vector; parameters, some of which may be logged
#'@param w perpendicular truncation distance
#'@return \eqn{\pi(x)} animal density at distance x
#'@examples
#'plot(seq(0,1,length=100),pi.norm(x=seq(0,1,length=100),logphi=c(0.5,log(0.3)),w=1))
#'@export
pi.norm=function(x,logphi,w)
  #-------------------------------------------------------------------------------
# Animal density function (with respect to perp dist, x).
# Inputs:
#  logtphi: theta is vector of parameters, some of which may be logged
#           (see below).
#  w      : perp. truncation dist.
#-------------------------------------------------------------------------------
{
  if(length(logphi)!=2) {
    cat(logphi,"\n")
    stop("logphi must be vector of length 2.")
  }
  if(any(x>w)) stop("x can't be greater than w")
  mu=logphi[1]
  sigma=exp(logphi[2])
  f=dnorm(x,mean=mu,sd=sigma)
  denom=(pnorm(w,mean=mu,sd=sigma)-pnorm(0,mean=mu,sd=sigma))
  if(denom>0) f=f/denom else f=0
  return(f)
}

#'@title Uniform perpendicular animal density function
#'
#'@description Uniform distribution for perpendicular animal density.
#'
#'@param x prependicular trackline distance
#'@param logphi Not used
#'@param w perpendicular truncation distance
#'@return \eqn{\pi(x)} constant animal density of 1/w
#'@examples
#'plot(seq(0,1,length=100),pi.const(x=seq(0,1,length=100),w=1))
#'@export
pi.const=function(x,logphi=NULL,w){
  return(rep(1/w,length(x)))}



#' @title Calculate HR perp dist function from hazard
#'
#' @description Implements the g(y) on page 36 or that on page 37 of Hayes and Buckland (1983)
#' from the hazard function given on page 37 or 38, respectively, of that paper (and implemented
#' in function \code{\link{h2}}).
#'
#' @param x perpendicular distance.
#' @param b vector of parameters of \code{\link{h2}} (on log scale).
#' @param hr hazard rate function name (must be character): only "h1" and "h2" valid
#' @return value of hazard rate detection function at x
#' @examples
#' b=log(c(0.75,0.1))
#' x=seq(0,1,length=50)
#' plot(x,hr.to.p(x,b,hr='h2'),type="l",ylim=c(0,1),xlab="Perpendicular distance",ylab="P(detect)")
#'@export
hr.to.p=function(x,b,hr){
  if(length(b)!=2) {
    cat(b,"\n")
    stop("b must be vector of length 2.")
  }
  theta=exp(b)
  if(theta[2]<=1) {
    warning("exp(b[2])<=1 so setting it equal to 1+1e-10")
    theta[2]=1+1e-10 # gamma can't deal with zero, so get close to it
  }
  a1=(theta[1]*gamma((theta[2]-1)/2)*gamma(0.5))/(2*gamma(theta[2]/2))
  return(1-exp(-a1*x^(-(theta[2]-1))))
}


#' @title Perpendicular animal density function calulated from hazard rate \code{h1}
#'
#' @description Calculates perpendicular animal density \eqn{\pi(x)} using the prependicular
#' distance function of the hazard rate function \eqn{k(x,y)=a*r^{-b}} on page 36  of Hayes and Buckland (1983)
#' from the hazard function given on page 37, of that paper (and implemented
#' in function \code{\link{hr1.to.p}}).
#'
#' @param x perpendicular distance.
#' @param logphi vector of parameters of \code{\link{h1}} (on log scale).
#' @param w perpendicular truncation distance
#' @return Animal density at x calculated from hazard rate \code{\link{h1}}
#' @examples
#' logphi=log(c(0.01,1.01))
#' x=seq(0,1,length=50)
#' plot(x,pi.hr1(x,logphi,w=1),type="l",xlab="Perpendicular distance",
#' ylab=expression(pi(x)))
#' @seealso \code{\link{h1}} \code{\link{hr1.to.p}}
#' @export
pi.hr1=function(x,logphi,w)
  #-------------------------------------------------------------------------------
# Animal density function (with respect to perp dist, x).
# Inputs:
#  logtphi: theta is vector of parameters, some of which may be logged
#           (see below).
#  w      : perp. truncation dist.
#-------------------------------------------------------------------------------
{
  if(length(logphi)!=2) {
    cat(logphi,"\n")
    stop("logphi must be vector of length 2.")
  }
  if(any(x>w)) stop("x can't be greater than w")
  f=hr1.to.p(x,b=logphi)/integrate(hr1.to.p,lower=0,upper=w,b=logphi)$value
  return(f)
}

#' @title Calculate harzard rate perpendicular distance function from hazard rate \code{h1}
#'
#' @description Implements the a prependicular distance function of the hazard rate function
#' k(x,y)=a*r^(-b) on page 36  of Hayes and Buckland (1983)
#' from the hazard function given on page 37 of that paper (and implemented
#' in function \code{\link{h1}}).
#'
#' @param x perpendicular distance.
#' @param b vector of parameters of \code{\link{h1}} (on log scale).
#' @param w perpendicular truncation distance
#' @return value of hazard rate detection function at x
#' @examples
#' b=log(c(0.01,1.01))
#' x=seq(0,1,length=50)
#' plot(x,hr1.to.p(x,b),type="l",ylim=c(0,1),xlab="Perpendicular distance",ylab="P(detect)")
#' @seealso \code{\link{h1}} \code{\link{pi.hr1}}
#' @export
hr1.to.p=function(x,b,w){
  if(length(b)!=2) {
    cat(b,"\n")
    stop("b must be vector of length 2.")
  }
  theta=exp(b)
  if(theta[2]<=1) {
    warning("exp(b[2])<=1 so setting it equal to 1+1e-10")
    theta[2]=1+1e-10 # gamma can't deal with zero, so get close to it
  }
  a1=(theta[1]*gamma((theta[2]-1)/2)*gamma(0.5))/(2*gamma(theta[2]/2))
  return(1-exp(-a1*x^(-(theta[2]-1))))
}


#' @title Perpendicular animal density function calulated from hazard rate \code{hr2}
#'
#' @description Calculates perpendicular animal density \eqn{\pi(x)} using the prependicular
#' distance function of the hazard rate function \eqn{k(r,y)=a \sqrt{(r^2-y^2)}/r^{(b+1)}; b>2}
#' on page 36  of Hayes and Buckland (1983) #' from the hazard function given on page 37,
#' of that paper (and implemented in function \code{\link{hr1.to.p}}).
#'
#' @param x perpendicular distance.
#' @param logphi vector of parameters of \code{\link{h1}} (on log scale).
#' @param w perpendicular truncation distance
#' @return Animal density at x calculated from hazard rate \code{\link{h1}}

#' @examples
#' x=seq(0,1,length=50)
#'plot(x,pi.hr2(x,logphi=c(-0.2876821, -2.3025851),w=1),
#'xlab='x',ylab=expression(pi(x)),type='l')
#' @seealso \code{\link{h1}} \code{\link{h2}} \code{\link{pi.hr1}}
#' @export
pi.hr2=function(x,logphi,w)
#-------------------------------------------------------------------------------
# Animal density function (with respect to perp dist, x).
# Inputs:
#  logtphi: theta is vector of parameters, some of which may be logged
#           (see below).
#  w      : perp. truncation dist.
#-------------------------------------------------------------------------------
{
  if(length(logphi)!=2) {
    cat(logphi,"\n")
    stop("logphi must be vector of length 2.")
  }
  if(any(x>w)) stop("x can't be greater than w")
  f=hr2.to.p(x,b=logphi)/integrate(hr2.to.p,lower=0,upper=w,b=logphi)$value
  return(f)
}

#' @title Calculate hazard rate perpendicular distance function from hazard rate \code{h2}
#'
#' @description Implements the a prependicular distance function of the hazard rate function
#' \eqn{k(r,y)=a \sqrt(r^2-y^2)/r^{(b+1)}; b>2} on page 37  of Hayes and Buckland (1983)
#' from the hazard function given on page 38 of that paper (and implemented
#' in function \code{\link{h2}}).
#'
#' @param x perpendicular distance.
#' @param b vector of parameters of \code{\link{h2}} (on log scale).
#' @param w perpendicular truncation distance
#' @return value of hazard rate detection function at x
#' @examples
#' x=seq(0,1,length=50)
#' plot(x,hr2.to.p(x,b=c(-0.2876821, -2.3025851),w=1),ylim=c(0,1),
#' xlab='Perp. distance, x', ylab='P(detect)',type='l')
#' @export
hr2.to.p=function(x,b,w){
  if(length(b)!=2) {
    cat(b,"\n")
    stop("b must be vector of length 2.")
  }
  theta=exp(b)
  return(1-exp(-theta[1]/(theta[2]+1)*x^(-(theta[2]+1))))
}

#'@title Waiting distance pdf
#'
#'@description Calculates the pdf of the 'waiting distance' \eqn{f(y,x)=h(y,x)*\exp(-\int_y^{ystart}) h(t,x) dt)}.
#'
#'@param y scalar or vector; forward distance
#'@param x scale or vector; perp. distance
#'@param b two-element vector of hazard rate parameters, some of which may be logged
#'@param hr hazard rate function
#'@param ystart max forward distance at which could possibly detect animal (see details).
#'@param nint number of intervals in numerical integration.
#'@return pdf of waiting distance at x,y
#'@details Need to ensure the hazard function has decayed to (very close to) zero by \code{ystart}.
#'@examples
#'w=1; ystart=4
#'gridx=seq(0,w,length=50); gridy=seq(0,ystart,length=50)
#'f=outer(gridy,gridx,FUN=fyx,b=c(-0.2876821, -2.3025851),hr=h2,ystart)
#'persp(gridx,gridy,t(f),theta=45,phi=35,zlab="f(y|x)")
#'@export
fyx=function(y,x,b,hr,ystart,nint=100)
{
  if(length(y)!=length(x)) stop("Lengths of x and y must be the same.")
  n=length(x)
  f=intval=rep(NA,n)
  if (!class(hr)=='character'){stop('message from fyx: hr must be passed as
                                    a character')}
  hr=match.fun(hr)

  for(i in 1:n) {
    dy=(ystart-y[i])/nint/2                           # for crude integration
    yy=seq(y[i],ystart,length=(nint+1))[-(nint+1)]+dy # for crude integration
    h=hr(yy,rep(x[i],nint),b)
    int=sum(hr(yy,rep(x[i],nint),b)*dy*2)  # crude integration
    intval[i]=exp(-int)
  }
  hrval=hr(y,x,b)
  bads=which(hrval>=.Machine$double.xmax) # identify infinite hazards
  if(length(bads)>0) { # infinite hazard so p(detect)=0
    f[bads]=.Machine$double.xmax
    f[-bads]=hr(y[-bads],x[-bads],b)*intval[-bads]
  }else{
    f=hr(y,x,b)*intval
  }
return(f)
}

#'@title Numerical calculation of perpendicular detection function from a hazard
#'
#'@description Calculates the perpendicular detection function, \eqn{p(x)}, for a given hazard.
#'
#'@param x scale or vector; perp. distance
#'@param b two-element vector of hazard rate parameters, some of whihc may be logged
#'@param hr hazard rate function
#'@param ystart max forward distance at which could possibly detect animal (see details).
#'@param nint number of intervals in numerical integration.
#'@return probability of detection at x
#'@examples
#'gridx=seq(0,1,length=50)
#'p.x=px(gridx,b=c(-0.2876821, -2.3025851),
#'  hr=h2,ystart=4,nint=100)
#'plot(gridx,p.x,type="l",ylim=c(0,max(p.x)),
#' xlab="prep. distance, x",ylab="p(x)")
#'@export
px=function(x,b,hrname,ystart,nint=100){
  if (!class(hrname)=='character'){stop('message from px: hr must be supplied as
                                    a character')}
  return(1-Sy(x,rep(0.0001,length(x)),ystart,b,hrname))
}

#'@title Product of p(x) and pi(x)
#'
#'@description Returns product of perp dist det prob px() and animal dbn pi.x
#'
#'@param y scalar or vector; forward distance observations
#'@param x scale or vector; perp. distance observations
#'@param pars c(b,logphi); hazard rate and density log-parameters in a vector (see details).
#'@param hr hazard rate function
#'@param ystart max forward distance at which could possibly detect animal (see details).
#'@param pi.x perpendicular distance density distribution
#'@param w perpendicular truncation distance.
#'@param length.b length of the hazard rate parameter vector
#'@return negative log likelihood for forward distance, \code{y} and perpendicular distance \code{x}.
#'@export
#'@details
#'Must to ensure the hazard function has decayed to (very close to) zero by \code{ystart}.
#'The parameter vector, \code{pars}, must be passed in with two parameters for the hazard rate first,
#'then two parameters for perpendicular density gradient \eqn{\pi(x)} i.e. \code{c(b,logphi)}.
#'@examples
#'p.pi.x(x,b,hr,ystart,pi.x,logphi,w)
p.pi.x=function(x,b,hr,ystart,pi.x,logphi,w){
  if (!class(hr)=='character'){stop('message from p.pi.x: hr must be supplied
                                    as character')}
  if (!class(pi.x)=='character'){stop('message from p.pi.x: pi.x must be supplied
                                    as character')}
  pi.x = match.fun(pi.x) # So that we can evaluate it as a function below

  return(px(x,b,hr,ystart)*pi.x(x,logphi,w))
}

# A few functions added to deal with rounding of radial distances to 0
# when they are below a value rmin. - Cal

#'@title y(x) for Sy(y(x)|x).
#'
#'@description retrieves the forward distance y at which an animal meets the
#' circumference of the truncation zone, given its perpendicular distance x,
# and the radius of the truncation zone, rmin.
#'
#'@param x numeric ; perpendicular distance of animal from the transect
#'@param rmin numeric ;  radius of the truncation zone. The distance below which
#' all perpendicular distances are rounded down to 0 - if this was done in the
#' the field.
#'@return y ; forward distance
#'@export
get.y.from.x = function(x,rmin){
  # we mst be able to deal with a vector of x inputs
  if (class(rmin)!='numeric'){stop('get.y.from.x: rmin must be numeric')}

  for (i in x){
    if (i>rmin){stop('get.y.from.x: x must be smaller than rmin')}
  }

  rmin = rep(rmin, length(x)) # vectorise rmin
  return(sqrt(rmin**2-x**2))
}

#'@title Sy(y(x)|x), when radial distances close to 0 have been rounded
#'
#'@description evaluates the Survival function at ymin, for a given x.
#' ymin is the forward distance at which an animal enters the truncation
#' zone of radius rmin.
#'
#'@param x numeric ; perpendicular distance of animal from the transect
#'@param rmin numeric ;  radius of the truncation zone. The distance below which
#' all perpendicular distances are rounded down to 0 - if this was done in the
#' the field.
#'@param ymax numeric ; the maximum forward distance at which which we can
#' detect animals.
#'@param b vector of numerics ; parameters for the detection function
#'@param hr character ; name of the chosen hazard rate
#'@return numeric ; survival function evaluated with relevant parameters
#'@export
Sy.of.x.rounded.data = function(x,rmin,ymax,hr,b){
  # integrate passes a vector of inputs, rather than calling
  # the function many times, we must be able to deal with this
  if (class('hr')!='character'){stop('hr must be a character')}
  y = get.y.from.x(x,rmin)

  return(Sy(x,y,ymax,b,hr))
}

#'@title density function with vecotr inputs
#'
#'@description The integrate function requires the ability to call
#' function it integrates over with a vector of values. The
#' density distribution functions do not allow this, this function
#' acts as a wrapper, accepting a vector of x value inputs.'
#'
#'@param xVector numeric vector; perpendicular distances of animals from the transect
#'@param pi.x.Function function ; the density function we wish to call
#'@param logphi numeric vector; parameters for the density function
#'@param w numeric ; perpendicular truncation distance
#'@return vector ; pi.x evaluated at each of the values in xVector
#'@export
pi.x.vector.input = function(xVector, pi.x.Function, logphi, w){
  n = length(xVector)
  outputVals = rep(NA, n)
  for (i in (1:n)){
    outputVals[i] = pi.x.Function(xVector[i], logphi=logphi, w=w)
  }

  if (anyNA(outputVals)){stop('Unexpected error evaluating density function')}

  return(outputVals)
}



#'@title intermediary function for rounded likelihood
#'
#'@description the function which we integrate to obtain the numerator
#' of the likelihood for rounded data
#'
#'@param x numeric ; perpendicular distance of animal from the transect
#'@param pi.x character ; user-chosen density model
#'@param logphi numeric vector ; parameters for pi.x
#'@param rmin numeric ; value below which radial distances have been
#' rounded to 0
#'@param ymax numeric ; maximum forward distance at which we can detect animals
#'@param hr character ; user-chosen hazard rate
#'@param b numeric vector ; parameters for hr
#'@param w numeric ; perpendicular truncation distance
#'@return numeric
#'@export
round.survival.prod.pi = function(x,pi.x,logphi,rmin,ymax,hr,b,w){
  if (length(x)==1){stop('expecting vector input of x values')}
  pi.x = match.fun(pi.x)

  Symin = Sy.of.x.rounded.data(x,rmin,ymax,hr,b) # returns a vector
  # Sy can only deal with both x and y being vectors, so we vecotrise
  # 0s to obtain the correct output:
  zeros = rep(0, length(x)) ; S0 = Sy(x,zeros,ymax,b,hr)

  PiVals = pi.x.vector.input(x, pi.x,logphi,w)

  if (length(Symin)==1|length(S0)==1|length(PiVals)==1){
    stop('Incorrect handling of vector inputs')
  }

  return( (Symin - S0)*PiVals)
}

#'@title numerator for rounded data likelihood
#'
#'@description numerator for rounded data likelihood
#'
#'@param pi.x character ; user-chosen density model
#'@param logphi numeric vector ; parameters for pi.x
#'@param rmin numeric ; value below which radial distances have been
#' rounded to 0
#'@param ymax numeric ; maximum forward distance at which we can detect animals
#'@param hr character ; user-chosen hazard rate
#'@param b numeric vector ; parameters for hr
#'@param w numeric ; perpendicular truncation distance
#'@return numeric
#'@export
round.lik.num = function(pi.x,logphi,rmin,ymax,hr,b,w){
  int=integrate(f=round.survival.prod.pi,
    lower=0,upper=rmin,pi.x=pi.x,logphi=logphi,rmin=rmin,
    ymax=ymax,b=b,hr=hr,w=w)
  return(int$value)
}

# left in for context but has been removed from the package
# after simulations confirmed suspisions that the likelihood
# did not require this term.
#'@title denominator for rounded data likelihood
#'
#'@description denominator for rounded data likelihood
#'
#'@param pi.x character ; user-chosen density model
#'@param logphi numeric vector ; parameters for pi.x
#'@param rmin numeric ; value below which radial distances have been
#' rounded to 0
#'@return numeric
#'@export
round.lik.denom = function(pi.x,logphi,rmin,w){
  pi.x = match.fun(pi.x)

  # We need pi.x to be able to accept a vector of x values:
  output = integrate(f=pi.x.vector.input, lower=0, upper=rmin,
    pi.x.Function = pi.x,logphi=logphi,w=w)$value

  return(output)
}


#'@title negative log likelihood for rounded data
#'
#'@description log likelihood for rounded data
#'
#'@param x numeric ; perpendicular distance of animal from the transect
#'@param pi.x character ; user-chosen density model
#'@param logphi numeric vector ; parameters for pi.x
#'@param rmin numeric ; value below which radial distances have been
#' rounded to 0
#'@param ymax numeric ; maximum forward distance at which we can detect animals
#'@param hr character ; user-chosen hazard rate
#'@param b numeric vector ; parameters for hr
#'@param w numeric ; perpendicular truncation distance
#'@param undrounded numeric ; the number of observations in the data set with
#' radial distance <= rmin
#'@return numeric ; negative log likelihood of the rounded data
#'@export
round.lik = function(rounded,pi.x,logphi,rmin,ymax,hr,b,w){
  if (class(hr)!='character'){stop('round.lik: hr must be a character')}
  if (class(pi.x)!='character'){stop('round.lik: pi.x must be a character')}

  int=integrate(f=p.pi.x,lower=0,   # We divide by the same
    upper=w,b=b,hr=hr,              # integral as the original
    ystart=ystart,pi.x=pi.x,        # likelihood. This is the 1/
    logphi=logphi,w=w)$value         # term in the paper

  TLN = round.lik.num(pi.x,logphi,rmin,ymax,hr,b,w)   # num
  #TLD = round.lik.denom(pi.x,logphi,rmin,w)           # denom

  #if(DENOM==TRUE) frac=TLN/TLD else
  frac=TLN  # The above line used to be needed when we
  # included an extra term in the likelihood. It is
  # uneccessary and has been removed, however the code
  # has been left in the file for completeness.
  return(-log(((frac)/int)**rounded))
}

#F.x=function(x,b,hr,ystart,pi.x,logphi,w) return((1-px(x,b,hr,ystart))*pi.x(x,logphi,w))

#'@title Negative log-likelihood for forward distance and perpendicular distance
#'
#'@description Calculates the negative log-likelihood for forward distance, \code{y}, and
#'perpendicular distance, \code{x}, for a given hazard and perpendicular density distribution.
#'
#'@param y scalar or vector; forward distance observations
#'@param x scale or vector; perp. distance observations
#'@param pars c(b,logphi); hazard rate and density log-parameters in a vector (see details).
#'@param hr hazard rate function
#'@param ystart max forward distance at which could possibly detect animal (see details).
#'@param pi.x perpendicular distance density distribution
#'@param w perpendicular truncation distance.
#'@param length.b length of the hazard rate parameter vector
#'@return negative log likelihood for forward distance, \code{y} and perpendicular distance \code{x}.
#'@details
#'Must to ensure the hazard function has decayed to (very close to) zero by \code{ystart}.
#'The parameter vector, \code{pars}, must be passed in with two parameters for the hazard rate first,
#'then two parameters for perpendicular density gradient \eqn{\pi(x)} i.e. \code{c(b,logphi)}.
#'@examples
#'ystart=4;w=1
#'hr=h2; b=log(c(0.75,1))
#'pi.x=pi.norm; logphi=c(0.5,log(0.2))
#'N=50 #true number of animals
#'#generate some observations
#'simDat=simXY(N=N,pi.x=pi.x,logphi=logphi,
#'hr=hr,b=b,w=w,ystart=ystart)
#'x=simDat$locs$x; y=simDat$locs$y
#'pars=c(b,logphi)
#'negloglik.yx(y,x,pars,hr,ystart,pi.x,w)
#'@seealso \code{\link{simXY}}
#'@export
negloglik.yx=function(pars,y,x,hr,ystart,pi.x,w,rmin,length.b=2,debug=FALSE,DENOM)
{
  if(length(y)!=length(x)) stop("Lengths of x and y must be the same.")
  if(debug) print(pars)
  if (!class(hr)=='character'){stop('message from negloglik: hr must
                                    be passed as a character')}
  if (!class(pi.x)=='character'){stop('message from negloglik: pi.x must
                                    be passed as a character')}

  # determine which values are rounded or not:
  x.rounded = x[sqrt(y**2 + x**2)<rmin]
  y.rounded = y[sqrt(y**2 + x**2)<rmin]
  new.x = x[sqrt(y**2 + x**2)>=rmin]
  new.y = y[sqrt(y**2 + x**2)>=rmin]
  x = new.x ; y = new.y # avoids writing over vars as we reassign
  rounded = length(x.rounded)

  if (length(x)!=length(y)){stop('X + Y')} # Safety check to make sure
  # rounding didn't make a mistake with data dimensions

  n=length(y) ; hrname = hr ; piname = pi.x
  # unpack parameters *** need to change if hr and pi.x don't have 2 pars each
  b=pars[1:length.b]
  if(piname=="pi.const") logphi=NULL else logphi=pars[(1+length.b):length(pars)]

  hr=match.fun(hr) ; pi.x=match.fun(pi.x)
  llik=rep(NA,n)
  # calculate numerator:
  num=sum(log(fyx(y,x,b,hrname,ystart)) + log(pi.x(x,logphi,w)))
  # calculate denominator:
  int=integrate(f=p.pi.x,lower=0,upper=w,b=b,hr=hrname,
                ystart=ystart,pi.x=piname,logphi=logphi,w=w)
  denom=log(int$value)

  # likelihood:
  llik=-(num-n*denom) # 2016 paper likelihood, for the un-rounded data points

  if(rounded>0){
    negllik.rounded = round.lik(rounded,pi.x=piname,
      logphi,rmin,ymax=ystart,hr=hrname,b,w,DENOM=DENOM)
  } # we calculate the rounded data points part of the likelihood

  else{negllik.rounded=0} # No addition to normal likelihood needed


  return(llik + negllik.rounded) #round.lik returns a neg.log.lik, so we add.
}

#' Simulate sightings given a known perpendicular density distribution and hazard function
#' @description Simulates sightings from a known population given a density distribution
#' and hazard rate.  This function has been replaced by \code{\link{sim.n}}.
#' @param N animal population
#' @param pi.x function describing the perpendicular distance density distribution
#' @param logphi parameters for pi.x (some maybe logged)
#' @param hr function describing the hazard rate
#' @param b hazard rate parameter vector
#' @param w truncation distance
#' @param ystart max forward distance at which could possibly detect animal (see details).
#' @param xSampL length of x-dimension vector to sample perpendicular distances from.
#' @param discardNotSeen boolean; discard individuals not detected.  See details
#' @param ... arguments to be passed into \code{\link{simnhPP}}
#' @details if \code{discardNotSeen=FALSE} individuals that are not detected are
#' assigned y-dimension distances = -999, otherwise \code{discardNotSeen=TRUE}
#' invididuals are removed and not returned
#' @return list of \code{$locs} x and y coordinates for simulated sightings and
#' \code{$settings} simulation settings.
#' @export
#'@examples
#'ystart=4;w=1
#'hr=h2; b=log(c(0.75,1))
#'pi.x=pi.norm; logphi=c(0.5,log(0.2))
#'N=50 #true number of animals
#'#generate some observations
#'simDat=simXY(N=N,pi.x=pi.x,logphi=logphi,
#'hr=hr,b=b,w=w,ystart=ystart)
#'x=simDat$locs$x; y=simDat$locs$y
#'
simXY=function(N,pi.x,logphi,hr,b,w,ystart,xSampL=5*N,discardNotSeen=TRUE,...)
{
  if (class(hr)!='character'|class(pi.x)!='character'){
    stop('SimXY: hr and pi.x must be supplied as characters')}

  xV=seq(0,w,length=xSampL)

  # It seems that for even for reasonable values, the call to simnhPP produces
  # an error roughly 2% of the time. To avoid this becoming an inconvenience
  # to the user more often than it should, I put in place a loop to retry the
  # call some reasonable number of times, before deciding that the user chosen
  # start values consistently lead to errors. - Cal

  lt2d.tryCounter = 0                         # odd name to avoid naming
  lt2d.Error = TRUE                           # conflict with user's Global Env

  while (lt2d.tryCounter<10 & lt2d.Error==TRUE){

    x=sample(x=xV,size=N,replace=TRUE,        # resample from the Xs, since its
           prob=match.fun(pi.x)               # these values which lead to
           (x=xV,logphi=logphi,w=w))          # the errors during integration

    y = try({simnhPP(x,b,ystart,hr,...)},     # we try to call the function,
      silent = T )                            # supressing the error to the user

    if (class(y)=="try-error"){
        lt2d.tryCounter = lt2d.tryCounter + 1 #+1 to count of unsuccessful tries
        lt2d.simXY.error.message = y[[1]]     # remember the error message
      }
      else{lt2d.Error = FALSE}                # succesful call;  stop loop
  }

  if (lt2d.Error){                            # If after 10 tries we still have
      stop(lt2d.simXY.error.message)          # an error, pass it to the user
  }

  remove(lt2d.tryCounter,lt2d.Error)          # remove floating global vars


  if(discardNotSeen){
    keep=which(y>=0)
    n=length(keep)
    x=x[keep]
    y=y[keep]}

  output = list(locs=cbind.data.frame(x,y),settings=
                list(N=N,pi.x=pi.x,logphi=logphi,
                     hr=hr,b=b,w=w,ystart=ystart,
                     discardNotSeen=discardNotSeen))
  class(output) = 'LT2D.simulated.data'
  return(output)
}

#'@title Plot the simulated positions
#'
#' @description Plots the simulated data.
#'
#' @param simDat object from a call of \code{\link{simXY}}
#' @examples
#'ystart=4;w=1
#'hr=h2; b=log(c(0.75,1))
#'pi.x=pi.norm; logphi=c(0.5,log(0.2))
#'N=500 #true number of animals
#'#generate some observations
#'simDat=simXY(N=N,pi.x=pi.x,logphi=logphi,
#'hr=hr,b=b,w=w,ystart=ystart)
#' plotSimNotUsed(simDat)
#' @export
plot.LT2D.simulated.data=function(simDat){
  x=simDat$locs$x; y=simDat$locs$y
  pi.x=simDat$settings$pi.x
  par(mfrow=c(2,2),mar=c(3,3,3,3),mgp=c(1.5,0.5,0))
  gridx=seq(0,simDat$settings$w,length=1000)

  adbn=match.fun(pi.x)(gridx,logphi=simDat$settings$logphi,
    w=simDat$settings$w)

  plot(gridx,adbn,type="l",
       xlim=c(0,simDat$settings$w),
       xlab="perp. dist (x)",ylab="pi(x)",
       main='Perp. density function')
  rug(x,ticksize=0.1)
  plot(x,y,
       xlim=c(0,simDat$settings$w),
       ylim=c(0,simDat$settings$ystart),
       cex=0.6,
       main='Sighting locations',
       xlab='perp. distance (x)',
       ylab='forward dist. (y)')
  mtext(paste('N=',simDat$settings$N,'; n=',nrow(simDat$locs)))
  hist(x,freq=FALSE,xlab="perp. dist. (x)",
       main="perp. dist. (x)",
       xlim=c(0,simDat$settings$w))
  hist(y,freq=FALSE,xlab="forward dist. (y)",
       main='Forward dist. (y)',
       xlim=c(0,simDat$settings$ystart))
}

histline <-
  function(height,breaks,lineonly=FALSE,outline=FALSE,fill=FALSE,xlim=range(breaks),
           ylim=range(height),xlab="x",ylab="y",
           transpose=FALSE,...){
    #-------------------------------------------------------------------------------------
    # Takes bar heights (height) and cutbpoints (breaks), and constructs a line-only
    # histogram from them using the function plot() (if lineonly==FALSE) or lines()
    # (if lineonly==TRUE).
    # If fill==TRUE, uses polygon() to fill bars
    # If fill==TRUE, valid arguments to plot() or lines() are passed via argument(s) "..."
    # If outline==TRUE, only outline of histogram is plotted
    # If fill!=TRUE, valid arguments to polygon() are passed via argument(s) "..."
    #
    # DLB 2009
    #-------------------------------------------------------------------------------------

    n=length(height)
    if(length(breaks)!=(n+1)) stop("breaks must be 1 longer than height")
    if(outline) {
      y=c(0,rep(height,times=rep(2,n)),0)
      x=rep(breaks,times=rep(2,(n+1)))
    }   else {
      y=rep(0,4*n)
      x=rep(0,4*n+2)
      for(i in 1:n) {
        y[((i-1)*4+1):(i*4)]=c(0,rep(height[i],2),0)
        x[((i-1)*4+1):(i*4)]=c(rep(breaks[i],2),rep(breaks[i+1],2))
      }
      x=x[1:(4*n)]
    }
    if(transpose)
    {
      xstore=x
      x=y
      y=xstore
      xlimstore=xlim
    }
    if(lineonly) {
      if(!fill) lines(x,y,...)
      else polygon(x,y,...)
    } else {
      if(!fill) plot(x,y,type="l",xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,...)
      else {
        plot(x,y,type="n",xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab)
        polygon(x,y,...)
      }
    }
  }


#'@title Plot the simulated positions
#'
#' @description Plots the simulated data.
#'
#' @param simDat object from a call of \code{\link{simXY}}
#' @param nclass number of bins in x and y histograms
#' @param image Boolean \code{TRUE} image background on sightings scatter plot.
#' @param xlab x-axis label for sightings scatter plot
#' @param ylab y-axis label for sightings scatter plot
#' @param ... other arguments to be passed into \code{\link{plotfit.y}}
#' @seealso \code{\link{simXY}} \code{\link{plotSim}} \code{\link{plotfit.y}}
#'@export
#'@examples
#'n=100;ymin=0.01;ymax=5;W=1
#' b=log(c(0.75,1));logphi=c(0.5,log(0.2))
#' simDat=sim.n(n,ymin,ymax,W,h2,b,pi.norm,logphi)
#'plotSim(simDat=simDat,nclass=12)
plotSim = function(simDat, nclass=10,xlab="", ylab="",image=FALSE,...){
  b=simDat$settings$b; hr=simDat$settings$hr; ystart=simDat$settings$ystart
  pi.x=simDat$settings$pi.x;logphi=simDat$settings$logphi; w=simDat$settings$w
  x=simDat$locs$x; y=simDat$locs$y
  pi.x=simDat$settings$pi.x

  zones=matrix(c(2,0,1,3), ncol=2, byrow=TRUE)
  layout(zones, widths=c(4/5,1/5), heights=c(1/5,4/5))
  xhist = hist(x, plot=FALSE,breaks=nclass)
  yhist = hist(y, plot=FALSE,breaks=nclass)
  top = max(c(xhist$counts, yhist$counts))
  par(mar=c(3,3,1,1))
  gridx=seq(0,w,length=50); gridy=seq(0,ystart,length=50)
  #f=outer(gridy,gridx,FUN=fyx,b=b,hr=hr,ystart=ystart) * pi.x(gridx,logphi=logphi,w=w)
  #gridx=seq(0,w,length=100)#c(0+dx,xhist$mids,w-dx)
  adbn=pi.x(gridx,logphi=logphi,w=w)

  est=list(b=b,hr=hr,ystart=ystart,pi.x=pi.x,logphi=logphi,w=w)
  ly=plotfit.y(y,x,est=est,nclass=nclass,plot=FALSE,lineonly=FALSE,nint=50,...)
  yd=ly$fy.[floor(seq(1,length(ly$fy.),length.out=50))]
  if(image){
    image(x=gridx,y=gridy,
          t(matrix(rep(adbn,length(gridy)),nrow=length(gridy),byrow=T)*yd),
          xlab=xlab,ylab=ylab)
    points(x,y,cex=0.4,pch=19)}
  if(!image)
    plot(x,y,xlim=c(0,w),ylim=c(0,ystart),xlab=xlab,ylab=ylab,pch=19,cex=0.4)
  par(mar=c(0,3,1,1))
  histline(height=xhist$density,breaks=xhist$breaks,
           xlim=c(0,w),xaxt='n',ylab='Density')
  #dx=mean(diff(xhist$mids))/2
  #adbn=adbn/sum(adbn)
  lines(gridx,adbn,lwd=2)#*length(x),lwd=2)

  par(mar=c(3,0,1,1))

  histline(height=yhist$density,breaks=yhist$breaks,transpose=TRUE,xlab='Density',
           xlim=c(0,max(ly$fy.)),ylim=c(0,ystart),yaxt='n',ylab='')
  #deny<<-ly$fy./sum(ly$fy.)
  #lines(ly$fy.,ly$gridy,lwd=2)
  lines(yd,gridy,lwd=2)
  par(oma=c(3,3,0,0))
  mtext(xlab, side=1, line=1, outer=TRUE, adj=0,
        at=.8 * (mean(x) - min(x))/(max(x)-min(x)))
  mtext(ylab, side=2, line=1, outer=TRUE, adj=0,
        at=(.8 * (mean(y) - min(y))/(max(y) - min(y))))
}

negloglik.yx2=function(y,x,ps,hr,b,ys,pi.x,logphi,w)
  #-------------------------------------------------------------------------------
# Returns negative log likelihood for forward dist y and perp. dist. x.
# Inputs:
#  y       : forward distances (scalar or vector)
#  x       : perp. distances (scalar or vector)
#  pars    : c(b,logphi); haz rate and density log-parameters in a vector
#  hr      : name of hazard rate function to use.
#  ystart  : max forward distance at which could possibly detect animal.
#            NB: need to ensure hazard has decayed to (very close to) zero by
#            this distance
#  pi.x    : name of animal density function to use.
#  w       : perp. truncation dist.
# ------
#  NOTE: code at *** must be changed if hr and pi.x don't have 2 pars each
# ------
#-------------------------------------------------------------------------------
{
  if(length(y)!=length(x)) stop("Lengths of x and y must be the same.")
  ystart=ys
  hr=match.fun(hr)
  pi.x=match.fun(pi.x)
  n=length(y)
  # unpack parameters *** need to change if hr and pi.x don't have 2 pars each
  #b=pars[1:2]
  #logphi=pars[3:4]
  llik=rep(NA,n)
  # caluclate numerator:
  num=sum(log(fyx(y,x,b,hr,ystart)) + log(pi.x(x,logphi,w)))
  # calculate denominator:
  int=integrate(f=p.pi.x,lower=0,upper=w,b=b,hr=hr,ystart=ystart,pi.x=pi.x,logphi=logphi,w=w)
  #  F.x=function(x,b,hr,ystart,pi.x,logphi,w) return((1-px(x,b,hr,ystart))*pi.x(x,logphi,w))
  #  int=integrate(f=F.x,lower=0,upper=w,b=b,hr=hr,ystart=ystart,pi.x=pi.x,logphi=logphi,w=w)
  denom=log(int$value)
  # likelihood:
  llik=num-n*denom
  #message(-llik)
  return(-llik)
}

# first stage in my effort to the S3 the code -C
#'@title LT2D fit object maker
#'
#'@description Used by \code{\link{fityx}} to produce a custom object which
#'contains relevant information pertaining to the produced fit, and which could
#'be needed by fitting and GOF functions. It type-checks all of the values to
#'ensure the user is never presented with an object with unexpected content.
#'
#'@param par numeric scalar or vector; values of fitted MLE parameters.
#'@param value numeric; value of the negative log likelihood evaluated at the MLEs.
#'@param counts integer, two-element vector which represents the number of
#'likelihood (and gradient function) evaluations performed by the optimisation
#'routine (\code{\link{optim}})
#'@param convergence integer; optim's system for describing fit behaviour. see
#'\code{\link{optim}}
#'@param message NULL or character; optim's output message, if any.
#'@param hessian matrix or NULL; solved hessian matrix of MLEs
#'@param error logical; boolean; has an error occurred or were any flags raised?
#'@param hr character; name of the user-chosen hazard function
#'@param pi.x character; name of the user-chosen density
#'@param ystart numeric; value of user-chosen ystart value
#'@param w numeric; user-chosen perpendicular truncation distance
#'@param b numeric; vector of user-chosen start values for hazard rate MLE
#'@param logphi numeric; scalar or vector of user-chosen start values
#'for density rate MLE
#'@param AICval numeric; Akaike's Information Criterion of the model (uncorrected)
#'@param dat data.frame; data set used to fit the model.
#'@param MLEvcov matrix; variance covariance matrix of MLEs
#'@param CVpar numeric vector; Coefficient of variation of the MLEs
#'@param corr matrix; matrix of correlations for the MLEs
#'@param p0 numeric; estimated probability of detection at distance (0,0)
#'@return LT2Dfit object
#'@seealso \code{\link{fityx}}
#'@export
LT2D.FitObjectMaker = function(par,value,counts,convergence,message,hessian,
                               error,hr,pi.x,ystart,w,b,logphi,AICval,dat,
                               MLEvcov,CVpar,corr, p0, rmin)
# largely useless function, type-checks and sanitises, but likely nothing
# incorrect ever gets passed to it anyway
{
  # The section that type-checks all of our inputs for consistency:
  if (!class(par)=='numeric' & class(value)=='numeric' &
      class(counts)=='integer' & class(convergence)=='integer' &
      (class(message)=='integer' || class(message)=='NULL') &
      (class(hessian)=='matrix'  || class(hessian)=='NULL') &
      class(error)=='logical' & class(hr)=='character' &
      class(pi.x)=='character' & class(ystart)=='numeric' &
      class(w)=='numeric' & class(b)=='numeric' & class(logphi)=='numeric' &
      class(AICval)=='numeric' & class(dat)=='data.frame' &
      class(MLEvcov)=='matrix' & class(CVpar)=='numeric' &
      class(corr)=='matrix' & class(p0)=='numeric' & class(rmin)=='numeric')

  {stop('One or more incorrect data types supplied to LT2D.FitObjectMaker ')}

  # The section that creates the output:
  output = list()
  output$par = par ; output$value = value ; output$counts = counts
  output$convergence = convergence ; output$message = message
  output$hessian = hessian ; output$error = error ; output$hr = hr
  output$pi.x = pi.x ; output$ystart = ystart ; output$w = w ; output$b = b
  output$logphi = logphi ; output$AIC = AICval ; output$dat = dat
  output$vcov = MLEvcov ; output$CVpar = CVpar ; output$corr = corr
  output$p0 = p0 ; output$rmin = rmin

  class(output) = 'LT2D.fit.object'

  return(output)
}



#'@title Maximum likelihood estimation for unknown hazard and perpendicular distance distribution
#'
#'@description Uses \code{\link{optim}} to obtain a MLE for the hazard function and animal
#'perpendicular distance distribution.  Functional forms for the hazard and perpendicular distance
#'distribution must be specified.
#'
#'@param y scalar or vector; forward distance observations
#'@param x scale or vector; perp. distance observations
#'@param b two-element vector of hazard rate parameters, some of whihc may be logged
#'@param hr hazard rate function
#'@param ystart max forward distance at which could possibly detect animal (see details).
#'@param pi.x perpendicular distance density distribution
#'@param logphi parameters for pi.x (some maybe logged)
#'@param w perpendicular truncation distance.
#'@param control see \code{\link{optim}} control
#'@param hessian return hessian.  See also \code{\link{optim}}.
#'@param corrFlag=0.7 Absolute parameter correlation value above which a warning is issued.
#'@param ... arguments to be passed into \code{\link{optim}}
#'@return
#'\code{\link{optim}} fit object and \cr
#'\code{$hr} = hazard rate function used.\cr
#'\code{$pi.x} = perpendicular distance function used.\cr
#'\code{$ystart} = ystart max forward distance detection used.\cr
#'\code{$w} = perpendicular truncation distance used.\cr
#'\code{$b} = estimated hazard parameters\cr
#'\code{$dat} = data frame with data (\code{$x} and \code{$y})\cr
#'\code{$logphi} \cr
#'\code{AIC} AIC value\cr
#'And if \code{hessian=TRUE}:\cr
#'\code{vcov} variance covariance matrix.  Will warn if there is a problem inverting
#'the hessian.\cr
#'\code{CVpar} Coefficient of variation for each paramter estimate. \cr
#'\code{error} Boolean, \code{TRUE} if convergence!=0 or problem inverting the hessian,
#'or parameter correlation is exceeded.\cr
#'@details Must to ensure the hazard function has decayed to (very close to) zero by
#'\code{ystart}.
#'
#'@examples
#'\dontrun{
#'ystart=4;w=1
#'hr=h2; b=log(c(0.75,1))
#'pi.x=pi.norm; logphi=c(0.5,log(0.2))
#'N=50 #true number of animals
#'#generate some observations
#'simDat=simXY(N=N,pi.x=pi.x,logphi=logphi,
#'hr=hr,b=b,w=w,ystart=ystart)
#'x=simDat$locs$x; y=simDat$locs$y
#'fityx(y,x,b,hr,ystart,pi.x,logphi,w)
#'}
#'@seealso \code{\link{negloglik.yx}}
#'@export
fityx = function(y,x,b,hr,ystart,pi.x,logphi,w,rmin=0,control = list(),hessian =
                   FALSE,corrFlag = 0.7,debug = FALSE,DENOM=FALSE,...)
{
  # The methodology used to retrieve function names made life very difficult
  # when functions were passed, as opposed to their character name. The
  # solution used unpassed variables within Hazard Rate and Density Distribution
  # functions, and lead to the use of eval and parse statements, whose use are
  # generally considered questionable.

  # I felt like this functionality could be expended, and that giving the user
  # the ability to choose between passing either the function or its name was
  # not worth this complicated overhead (which didn't work at all unless the
  # functions were in the Global Environment, as having them loaded by the
  # package only, lead to an error whereby this function could not find 'fName').

  # These fName variables have been removed from the hr and pi functions, and
  # replaced with a check which reminds the user to pass the appropriate string
  # (character) name into the fit function.  I have also removed the fNameFinder
  # function, since it was the eval / parse culprit, and is no longer needed
  # with these changes. - Cal
  
  print ('in use 101')
  HazardWarning = 'Incorrect option supplied for Hazard Rate. Please supply
  the name of desired function as a character. Use SeeHazardOptions()
  to see a list of available choices'

  DensityWarning = 'Incorrect option supplied for Density. Please supply
  the name of desired function as a character. Use SeeDensityOptions()
  to see a list of available choices'
  if (!length(y)==length(x)){stop('x must have the same length as y')}

  # Truncate the data to perp <= w. It's important to calculate y first,
  # or else nothing will get truncated:

  if (w){
    # we don't want to truncate points with radial distance smaller than Rmin

    y = y[x<=w] ; x = x[x<=w]
    warning('data truncated according to user\'s chosen perpendicular truncation
  distance. Data in model object may be a different dimension to the
  supplied data')
  }

  if (!class(hr)=='character'){           # We check that the passed functions
    stop(HazardWarning)                   # are strings, and hence a name -C
  }

  if (!class(pi.x)=='character'){
    stop(DensityWarning)
  }

  # We're going to exception handle the function names passed, to see that they
  # do indeed exist. I impose a further restriction, which is that the functions
  # exist within the package itself, to avoid match.fun() false positives
  # with functions in the users' GlobalEnv.


  hrname = hr                             # These lines keep the variable names
  piname = pi.x                           # consistent with the 'old' code below
  length.b = length(b)

  #change when name of package changes...
  if (!(piname %in% ls("package:LT2D"))&(hrname %in% ls("package:LT2D"))){
    stop('Invalid choice for hr or pi.x')
  }

  if (piname == "pi.const") pars = b else pars = c(b,logphi)
  # logphi has no default. The above line implies that if pi.const is chosen,
  # logphi is not needed... Hence the user will have no reason to pass this to
  # the fit function, will this cause an issue - Verify - Cal

  fit = optim(
    par = pars,fn = negloglik.yx,y = y,x = x,hr = hrname,ystart = ystart,
    pi.x = piname,w = w, length.b = length.b,
    hessian = hessian, debug = debug,
    control = control, rmin=rmin,DENOM = DENOM,...
  )

  error = FALSE
  if (fit$convergence != 0) {
    warning('Convergence issue (code = ',
            fit$convergence,') . Check optim() help.')
    error = TRUE
  }

  # ***
  if (length.b != length(pars)) {
    logphi = fit$par[(1 + length.b):length(pars)]
  }else{
    logphi = NA
  }
  # ***

  if (hessian) {
    mNames = paste('b',1:length.b,sep = '')
    if (!all(is.na(logphi)))
      mNames = c(mNames,paste('logphi',1:length(logphi),sep = ''))
    if (any(diag(fit$vcov) <= 0)) {
      warning('Failed to invert hessian.  Model covergance problem in fityx?')
      error = TRUE
      CVpar = rep(NA,length(fit$par))
    } else {
      CVpar = sqrt(diag(solve(fit$hessian))) / abs(fit$par)
    }

    vcov=solve(fit$hessian)
    corr = cov2cor(vcov)
    row.names(corr) = mNames
    colnames(corr) = mNames
    corr[upper.tri(corr,diag = TRUE)] = NA

    corrIND = which(abs(corr) > corrFlag,arr.ind = T)
    if (nrow(corrIND)) {
      warning(
        'absolute correlation exceeds ',corrFlag,' in parameter estimates: ',
        paste(paste(mNames[corrIND[,1]],mNames[corrIND[,2]],sep =
                      ' to '),collapse = '; ')
      )
      error = TRUE
    }
  }

  AICval = 2 * fit$value + 2 * length(fit$par)
  b = fit$par[1:length.b]
  dat = data.frame(x = x,y = y)

  p0 = 1 - Sy(0,0,ystart,b,hrname)
  
  # REMOVE AFTER DEBUG:
  
  # finalfit = LT2D.FitObjectMaker(par = fit$par, value = fit$value,
  #                     counts = fit$counts, convergence = fit$convergence,
  #                     message = fit$message, hessian = fit$hessian,
  #                     error = error, hr = hrname, pi.x = piname, w = w,
  #                     b = b, logphi = logphi, AICval = AICval, dat = dat,
  #                     MLEvcov = vcov, CVpar = CVpar, corr = corr,
  #                     ystart = ystart, p0=p0, rmin=rmin)
  
  finalfit = list()
  finalfit$counts = fit$counts
  finalfit$convergence = fit$convergence
  finalfit$message = fit$message
  finalfit$hessian = fit$hessian
  finalfit$hr = hrname
  finalfit$pi.x = piname
  finalfit$w = w
  finalfit$b = b
  finalfit$logphi = logphi
  finalfit$AICval = AICval
  finalfit$dat = dat
  finalfit$ystart = ystart
  finalfit$p0 = p0
  finalfit$rmin = rmin
  
  if (hessian==TRUE){
    finalfit$error = error
    finalfit$MLEvcov = vcov
    finalfit$CVpar = CVpar
    finalfit$corr = corr
  }
  
  class(finalfit) = 'LT2D.fit.object'
  return(finalfit)
}

#' @title Estimates density and abundance.
#'
#' @description
#' Horvitz-Thompson like estimation of density and abundance of groups and of individuals, as well as
#' of group size (estimated as the ratio of individual density and group density estimates). Produces
#' estimates by stratum and over all strata.
#'
#' @param dat DS data frame. Must have cols "stratum","area","transect","L","size","object","x","y"
#' (and possibly others).
#' @param hmltm.fit output from \code{\link{fit.hmltm}}.
#' @param W perpendicular truncation distance for estimation.
#'
#' @export
NDest <- function(dat,hmltm.fit){

  W = hmltm.fit$w # hmltm.fit should be an LT2D.model.fit object, which has attribute w

  # remove the smaller than w observations from the data frame which won't have been
  # used to fit the likelihood:
  dat = subset(dat, dat$x<W | is.na(dat$x)) # We musn't remove the NAs

  # Add 1/p column
  dat$invp <- rep(NA,dim(dat)[1])
  invp <- invp1_replacement(dat,hmltm.fit)

  for(i in 1:length(invp$object)) {
    row <- which(dat$stratum==invp$stratum[i]
      & dat$transect==invp$transect[i]
      & dat$object==invp$object[i])

    if(length(row)>1) {
      cat("Target stratum:",invp$stratum[i],"\n")
      cat("Target transect:",invp$transect[i],"\n")
      cat("Target sighting:",invp$object[i],"\n")
      cat("Found >1: at rows",row,"\n")
      stop("")
    }

    dat$invp[row] <- invp$invp[i]
  }

  # Calculate density and abundance by stratum

  strat <- unique(dat$stratum)
  nstrat <- length(strat)
  n <- L <- a <- A <- Dg <- D <- Ng <- N <- sbar <- rep(0,nstrat+1)
  stratname <- rep("",nstrat+1)

  for(i in 1:nstrat){
    stratname[i] <- as.character(strat[i])
    vdat <- dat[dat$stratum==strat[i],]
    trans <- unique(vdat$transect)
    L.tr <- 0

    for(tr in 1:length(trans)){
      L.tr <- L.tr+vdat$L[min(which(vdat$transect==trans[tr]))]
    }

    L[i] <- L.tr
    a[i] <- L[i]*2*W
    A[i] <- vdat$area[1]
    svdat <- vdat[!is.na(vdat$object),]
    n[i] <- length(svdat$invp)
    Dg[i] <- sum(svdat$invp)/a[i]
    D[i] <- sum(svdat$size*svdat$invp)/a[i]
    sbar[i] <- D[i]/Dg[i]
    Ng[i] <- Dg[i]*A[i]
    N[i] <- D[i]*A[i]
  }

  stratname[nstrat+1] <- "Total"
  Ng[nstrat+1] <- sum(Ng[1:nstrat])
  N[nstrat+1] <- sum(N[1:nstrat])
  A[nstrat+1] <- sum(A[1:nstrat])
  Dg[nstrat+1] <- Ng[nstrat+1]/sum(A[1:nstrat])
  D[nstrat+1] <- N[nstrat+1]/sum(A[1:nstrat])
  n[nstrat+1] <- sum(n[1:nstrat])
  L[nstrat+1] <- sum(L[1:nstrat])
  a[nstrat+1] <- sum(a[1:nstrat])
  sbar[nstrat+1] <- D[nstrat+1]/Dg[nstrat+1]

  # add transect frequency:
  tfreq <- apply(table(dat$stratum,dat$transect)>0,1,sum)
  k <- c(tfreq,sum(tfreq))

  return(list(invp=invp,
              ests=data.frame(stratum=stratname,n=n,k=k,L=L,covered.area=a,stratum.Area=A,
                              Dgroups=signif(Dg,3),Ngroups=signif(Ng,3),mean.size=round(sbar,1),
                              D=signif(D,5),N=round(N,1))
  )
  )
}

#' @title Adds inverse detection probability column to dataset
#'
#' @description
#' Adds inverse detection probability column to dataset, so that the NDest function
#' has the input it requires to produce its estimations.
#'
#' @param LT2D.df; an input data frame of the same format as required by NDest.
#' @param LT2D.fit output from \code{\link{fit.yx}}.
#' @export
invp1_replacement = function(LT2D.df,LT2D.fit){
  # From the fitted model we extract the invp value:
  p = phat(LT2D.fit)
  inversep = 1/p
  invp.vector = rep(inversep, length(LT2D.df$x))
  # And we add it to data frame:
  LT2D.df.new = LT2D.df              # copy the input fitted data frame
  LT2D.df.new$invp = invp.vector     # add the invp column
  return(LT2D.df.new)
}

# A few functions to deal with adding in the influence of
# covariates into the shape paramaters of detection hazards,
# followed by the top level fitting function available to
# the package user:


#' @title Top level fitting and abundance estimation function for LT2D user
#' @description This function is a wrapper for the fityx function. It takes
#' in a 'Distance' style data frame, type-checks it, finds the MLEs using
#' \code{\link{fityx}} and then calls an abundance function \code{\link{NDest}}
#' to obtain density and abundance estimates by stratum.
#' @param DataFrameInput, data.frame ; data frame with required columns
#'  stratum (stratum number of observation), transect (transect number of
#' observation, object(object number of detection, NA if transect with no
#' detections), size (group size of detection), area (stratum area of observation)
#' and L (transect length of transect).
#' @param hr character; name of hazard rate to fit
#' @param b numeric; vector of start parameters for hr
#' @param ystart numeric; furthest possible forwards distance at which we can
#' detect animals.
#' @param pi.x character; name of perpendicular density to fit
#' @param logphi numeric; vector of start parameters for pi.x
#' @param w numeric; perpendicular truncation distance
#' @param rmin numeric; radial distance below which all values were rounded down to 0
#' @param hessian boolean
#' @param corrFlag numeric; value above which correlation flag is raised
#' @export
LT2D.fit = function(DataFrameInput,hr,b,ystart,pi.x,logphi,w,rmin=0,
                control = list(),hessian=TRUE,corrFlag = 0.7,
                debug = FALSE,DENOM=FALSE){
  
  print ('in use 102')
  # Basic type-checking of inputs:
  if (class(DataFrameInput)!='data.frame'){
    stop('First arg must be data.frame')
  }

  OnlyCallFityx = TRUE # Set flag saying we can't work out N or D
  if ( is.null(DataFrameInput$stratum)
      | is.null(DataFrameInput$transect)
      | is.null(DataFrameInput$object)
      | is.null(DataFrameInput$size)
      | is.null(DataFrameInput$area)
      | is.null(DataFrameInput$L)){
     warning('Insufficient information for abundance calculation')
  }
  else{OnlyCallFityx=FALSE}

  if (is.null(DataFrameInput$x) | is.null(DataFrameInput$y)){
    stop('LT2D: Perpendicular and forward distances must be included in data.frame')
  }

  # We remove the NAs (non-detections) to fit the likelihood, but we
  # retain the rows so that we can add them back in when it comes to
  # calling the abundance function.
  NoNAs = subset(DataFrameInput, !is.na(DataFrameInput$object))

  X = NoNAs$x ; Y = NoNAs$y
  # run the fityx call:
    # We don't exception handle the call, to allow its own errors to be flagged
  fitted.model = fityx(y=Y,x=X,b=b,hr=hr,ystart=ystart,pi.x=pi.x,logphi=logphi,
                       w=w,rmin=rmin,control=control,
                       hessian,corrFlag)
  if (OnlyCallFityx==TRUE){return(fitted.model)}

  # Now we pass the data frame and fitted model objects to the abundance
  # estimation function and return the output
  output = NDest(DataFrameInput, fitted.model)
  output$fit = fitted.model # We attach the fityx model to the output
  class(output) = 'LT2D.fit.function.object'
  return(output)
}

negloglik.yx.w=function(y,x,pars,hr,ystart,pi.x,logphi,w)
#-------------------------------------------------------------------------------
# Returns negative log likelihood for forward dist y and perp. dist. x, taking
# pi(x) as known.
# Corresponding to Hayes and Buckland (1983) k(r,y)=a*sqrt(r^2-y^2)/r^(b+1); b>2
# on p37.
# Note: I use x for perp. dist., they use y.
# Inputs:
#  y       : forward distances (scalar or vector)
#  x       : perp. distances (scalar or vector)
#  pars    : b; haz rate log-parameters in a vector
#  hr      : name of hazard rate function to use.
#  ystart  : max forward distance at which could possibly detect animal.
#            NB: need to ensure hazard has decayed to (very close to) zero by
#            this distance
#  pi.x    : name of animal density function to use.
#  logphi  : vector of animal density function parameters (some logged)
#  w       : perp. truncation dist.
# ------
#  NOTE: code at *** must be changed if hr and pi.x don't have 2 pars each
# ------
#-------------------------------------------------------------------------------
{
  if(length(y)!=length(x)) stop("Lengths of x and y must be the same.")

  hr=match.fun(hr)
  pi.x=match.fun(pi.x)
  #print(hr);print(pi.x)
  n=length(y)
  llik=rep(NA,n)
  # caluclate numerator:
  num=sum(log(fyx(y,x,b=pars,hr,ystart)) + log(pi.x(x,logphi,w)))
  # calculate denominator:
  int=integrate(f=p.pi.x,lower=0,upper=w,b=pars,hr=hr,ystart=ystart,pi.x=pi.x,logphi=logphi,w=w)
  #  F.x=function(x,b,hr,ystart,pi.x,logphi,w) return((1-px(x,b,hr,ystart))*pi.x(x,logphi,w))
  #  int=integrate(f=F.x,lower=0,upper=w,b=b,hr=hr,ystart=ystart,pi.x=pi.x,logphi=logphi,w=w)
  denom=log(int$value)
  # likelihood:
  llik=num-n*denom

  return(-llik)
}

fityx.w=function(y,x,b,hr,ystart,pi.x,logphi,w,control)
{
  pars=b
  fit=optim(par=pars,fn=negloglik.yx.w,y=y,x=x,hr=hr,ystart=ystart,pi.x=pi.x,logphi=logphi,w=w,hessian=FALSE,control=control)
  fit$hr=hr
  fit$pi.x=pi.x
  fit$ystart=ystart
  fit$w=w
  fit$b=fit$par
  fit$logphi=logphi
  return(fit)
}




negloglik.x=function(x,pars,hr,ystart,pi.x,logphi,w,nint=100)
  #-------------------------------------------------------------------------------
# Returns negative log likelihood for perp. dist. x. given dbn pi.x and logphi
# Corresponding to Hayes and Buckland (1983) k(r,y)=a*sqrt(r^2-y^2)/r^(b+1); b>2
# on p37.
# Note: I use x for perp. dist., they use y.
# Inputs:
#  x       : perp. distances (scalar or vector)
#  pars    : b; haz rate log-parameters in a vector
#  hr      : name of hazard rate function to use.
#  ystart  : max forward distance at which could possibly detect animal.
#            NB: need to ensure hazard has decayed to (very close to) zero by
#            this distance
#  pi.x    : name of animal density function to use.
#  logphi  : vector of animal density function parameters (some logged)
#  w       : perp. truncation dist.
#-------------------------------------------------------------------------------
{
  hr=match.fun(hr)
  pi.x=match.fun(pi.x)
  n=length(x)
  # unpack parameters *** need to change if hr and pi.x don't have 2 pars each
  b=pars
  llik=rep(NA,n)
  # caluclate numerator:
  num=sum(log(px(x,b,hr,ystart,nint)) + log(pi.x(x,logphi,w)))
  # calculate denominator:
  int=integrate(f=p.pi.x,lower=0,upper=w,b=b,hr=hr,ystart=ystart,pi.x=pi.x,logphi=logphi,w=w)
  #  F.x=function(x,b,hr,ystart,pi.x,logphi,w) return((1-px(x,b,hr,ystart))*pi.x(x,logphi,w))
  #  int=integrate(f=F.x,lower=0,upper=w,b=b,hr=hr,ystart=ystart,pi.x=pi.x,logphi=logphi,w=w)
  denom=log(int$value)
  # likelihood:
  llik=num-n*denom

  return(-llik)
}


fitx=function(x,b,hr,ystart,pi.x,logphi,w,control)
{
  pars=b
  fit=optim(par=pars,fn=negloglik.x,x=x,hr=hr,ystart=ystart,pi.x=pi.x,logphi=logphi,w=w,hessian=FALSE,control=list(trace=5))
  fit$hr=hr
  fit$pi.x=pi.x
  fit$ystart=ystart
  fit$w=w
  fit$b=fit$par           # ***
  fit$logphi=logphi
  return(fit)
}




Eyx=function(y,x,b,hr,ystart)
  #Eyx=function(y,x,b,hr,ystart,nint=500)
  #-------------------------------------------------------------------------------
# Returns Expectation \int_y^ystart h(t,x) dt.
# Inputs:
#  y       : forward dist. (scalar or vector)
#  x       : perp. dist. (scalar or vector)
#  b: log(theta), where theta is vector of hazard rate parameters
#  hr      : name of hazard rate function to use.
#  ystart  : max forward distance at which could possibly detect animal.
#            NB: need to ensure hazard has decayed to (very close to) zero by
#            this distance
#-------------------------------------------------------------------------------
{
  if(length(y)!=length(x)) stop("Lengths of x and y must be the same.")
  if (class(hr)!='character'){stop('Eyx: hr must be supplied as character')}
  n=length(x)
  int=rep(NA,n)
  hr=match.fun(hr)
  ylo=1e-5  # set to avoid evaluating hr at y=0, which gives Inf
  for(i in 1:n) {
    y0=max(y[i],ylo)
    #    dy=(ystart-y0)/nint/2                           # for crude integration
    #    yy=seq(y0,ystart,length=(nint+1))[-(nint+1)]+dy # for crude integration
    #    int[i]=sum(hr(yy,rep(x[i],nint),b)*dy*2)  # crude integration
    int[i]=integrate(f=hr,lower=max(y[i],ylo),upper=ystart,x=x[i],b=b)$value
  }
  return(int)
}



#' @title Simulate forward distances
#'
#' @description Simulate non-homogeneous Poisson Process data by solving inverse CDF.
#' This function has been replaced by \code{\link{sim.n}}.
#'
#'@param x scale or vector; perpendicular distance observations
#'@param b two-element vector of hazard rate parameters, some of whihc may be logged
#'@param ystart max forward distance at which could possibly detect animal (see details).
#'@param hr hazard rate function
#'@param miss If TRUE, allows animals not detected by y=0 to have no y.
#            Misses are indicated by y==-999
#'@param ylo minimum forward distance.
#'@return vector of y distances if \code{miss=TRUE} missed animals represented by -999, otherwise
#'no animals are missed.
#'@details \code{miss} allows a fixed number of animal to be simulated i.e. known $N$ in strip,
#'or $n$ where animals may have been missed.
#' @examples
#'ystart=4
#'hr=h2; b=log(c(0.75,1))
#'X=runif(50,0,1)
#'Y=simnhPP(x=X,b=b,ystart=ystart,
#'hr=h2,miss=TRUE)
#'length(Y[Y!=-999])
#'Y=simnhPP(x=X,b=b,ystart=ystart,
#'hr=h2,miss=FALSE)
#'length(Y[Y!=-999])
#' @seealso \code{\link{simXY}}
simnhPP=function(x,b,ystart,hr,miss=TRUE,ylo=1e-5)
#-------------------------------------------------------------------------------
# Simulate non-homogeneous Poisson Process data by solving inverse CDF.
# Inputs:
#  x       : perp. dist. (scalar or vector)
#  b: log(theta), where theta is vector of hazard rate parameters
#  hr      : name of hazard rate function to use.
#  ystart  : max forward distance at which could possibly detect animal.
#            NB: need to ensure hazard has decayed to (very close to) zero by
#            this distance.
#  miss    : If TRUE, allows animals not detected by y=0 to have no y.
#            Misses are indicated by y==-999.
#-------------------------------------------------------------------------------
{
  if (class(hr)!='character'){stop('simnhPP: hr must be supplied as character')}
  obj=function(y,x,b,hr,ystart,u) return(((1-exp(-Eyx(y,x,b,hr,ystart)))-u)^2)
  n=length(x)
  u=runif(n)
  y=rep(-999,n)
  for(i in 1:n) {
    if(!miss) {
      while(y[i]<0) {
        if(u[i]>(1-exp(-Eyx(ylo,x[i],b,hr,ystart)))) u[i]=runif(1) # y<0, so missed: try again
        else {
          ymin=optimize(f=obj,interval=c(ylo,ystart),x[i],b,hr,ystart,u[i])
          y[i]=ymin$minimum
        }
      }
    } else if(u[i]<=(1-exp(-Eyx(ylo,x[i],b,hr,ystart)))) {
      ymin=optimize(f=obj,interval=c(ylo,ystart),x[i],b,hr,ystart,u[i])
      y[i]=ymin$minimum
    }
  }
  return(y)
}



#'@title Plot fitted hazard and perpendicular denisty distribution
#'
#'@description Plot fitted hazard and perpendicular denisty distribution
#'resulting from a call of \code{\link{fityx}}.  When simulating, true hazard functions
#'and perpendicular density distribution can also be added to the plot.

#'@param x perpendicular distance observations
#'@param fit return from a call of \code{\link{fityx}}
#'@param nclass number of histogram classes
#'@param nint number of intervals in numerical integration
#'@param plot boolean, plot results
#'@param addTruth boolean add true hazard and perp. density when simulating
#'@param true.pi.x true perpendicular density distribution function used when simulating
#'@param true.logphi true perpendicular density distribution function
#'parameters used when simulating
#'@param true.hr true hazard rate function used when simulating.
#'@param true.b true hazard rate function parameters used when simulating.
#'@param N true number of animals in simulated distribution, used to calculate bias (see details).
#'@param true.legend If true (and addTruth) plots legend for true functions in bottom left.
#'
#'@details When \code{N} is specified, bias in estimated number of animals ,\eqn{\hat N},
#' is calculated.
#'
#'@return
#'list with:
#'\code{$gridx} = x values used in plotting
#'\code{$p.xpifit} = product of perpendicular distance det probability p(x)
#'and perpendicular animal  distribution \eqn{\pi(x)}.
#'\code{$mufit} = effective strip width \eqn{\hat p}
#'\code{$f.xfit} =
#'\code{$p.xfit} =
#'\code{$ptot} =
#'\code{$p.xfit.std} =
#'\code{$adbn} =
#'\code{$N} = true number of animals in population
#'\code{$n} = number of seen animals
#'\code{$Nhat} = estimated number of animals.
#'\code{$bias} = \eqn{\hat N} bias.

#'@examples
#'\dontrun{
#'ystart=4;w=1
#'hr=h2; b=log(c(0.75,1))
#'pi.x=pi.norm; logphi=c(0.5,log(0.2))
#'N=50 #true number of animals
#'#generate some observations
#'simDat=simXY(N=N,pi.x=pi.x,logphi=logphi,
#'hr=hr,b=b,w=w,ystart=ystart)
#'x=simDat$locs$x; y=simDat$locs$y
#'est.yx=fityx(y,x,b,hr,ystart,pi.x,logphi,w)
#'plotdat.yx=plotfit.x(x,est.yx,addTruth=TRUE,true.logphi=logphi,true.b=b,N=N)
#'}
#'@seealso \code{\link{fityx}}
#'@export
plotfit.x=function(est,nclass=10,nint=100,    # est is a fitted LT2D model
                   plot=TRUE,dotitle="FALSE",
                   addTruth=FALSE,
                   true.pi.x=NULL,
                   true.logphi=NULL,
                   true.hr=NULL,
                   true.b=NULL,
                   N=NULL,...)
{
  # Some type-checking to ensure we can extract the data in the usual way
  if (class(est)!='LT2D.fit.object'){stop('Can only plot LT2D objects')}

  x = est$dat$x     # We exract the perpendicular distances from the fit
  Nhat.yx=bias=NULL
  b=est$b; hrname=est$hr; ystart=est$ystart; piname=est$pi.x
  logphi=est$logphi; w=est$w
  f.x=p.x.std=adbnTRUE=0
  # calculate stuff to plot:
  gridx=seq(1e-10,w,length=100)
  # first do f(x)
  p.xpifit=p.pi.x(gridx,b,hr=hrname,ystart,pi.x=piname,logphi,w)
  mufit=integrate(f=p.pi.x,lower=0,upper=w,b=b,hr=hrname
                  ,ystart=ystart,pi.x=piname,logphi=logphi,w=w)$value
  f.xfit=p.xpifit/mufit
  p.xfit=px(gridx,b,hr=hrname,ystart,nint=nint)
  ptot=integrate(f=px,lower=0,upper=w,b=b,hr=hrname,ystart=ystart)$value
  p.xfit.std=p.xfit/ptot
  adbn=match.fun(piname)(gridx,logphi,w) # changed to use piname instead of pi.x

  if(addTruth) {   # Haven't really checked this yet - Cal
    if(!is.null(true.pi.x)) pi.x=true.pi.x
    if(!is.null(true.logphi)) logphi=true.logphi
    if(!is.null(true.hr)) hr=true.hr
    if(!is.null(true.b)) b=true.b
    p.xpi=p.pi.x(gridx,b,hr,ystart,pi.x,logphi,w)
    mu=integrate(f=p.pi.x,lower=0,upper=w,b=b,hr=hr,ystart=ystart,pi.x=pi.x,logphi=logphi,w=w)$value
    f.x=p.xpi/mu
    adbnTRUE=pi.x(gridx,logphi,w)
  }
  if(plot){
    breaks=seq(0,w,length=(nclass+1))
    hx=hist(x,breaks=breaks,plot=FALSE) # get hist bar heights
    ymax=max(f.xfit,p.xfit.std,adbn,f.x,p.x.std,adbnTRUE,hx$density)
    main=""
    if(dotitle) main="Fitted curves"
    if(addTruth) main="Fitted curves (grey=true)"
    hx=hist(x,breaks=breaks,freq=FALSE,ylim=c(0,ymax),
            main=main,xlab="perpendicular distance (x)",ylab="pdf")
    lines(gridx,f.xfit,lwd=1)
    # overlay p(x), scaled to have area=1
    lines(gridx,p.xfit.std,lty=2,col="black",lwd=2)
    # overlay animal pdf:
    lines(gridx,adbn,lty=3,col="black",lwd=2)

    if(addTruth) legend("topright",title="Estimated",legend=c("f(x)","p(x)",expression(pi(x)),
           col=c("black","black","black"),lwd=c(2,2,2),lty=c(1,2,3)))
    else legend("topright",legend=c("f(x)","p(x)",expression(pi(x))),
                col=c("black","black","black"),lwd=c(2,2,2),lty=c(1,2,3))
    if(addTruth){
      lines(gridx,f.x,col="grey",lwd=2)
      p.x=px(gridx,b,hr,ystart,nint=nint)
      ptot=integrate(f=px,lower=0,upper=w,b=b,hr=hr,ystart=ystart)$value
      #      p.x.std=p.xfit.std=p.xfit/ptot
      #      lines(gridx,p.x*p.x.std,col="grey",lty=2,lwd=2)
      p.x.std=p.x/ptot
      lines(gridx,p.x.std,col="grey",lty=2,lwd=2)

      lines(gridx,adbnTRUE,col="grey",lty=3,lwd=2)
    }
  }
  if(!is.null(N)){
    n=length(x)
    Nhat.yx=n/mufit
    bias=(Nhat.yx/N-1)*100
    msg=paste("N=",N,"; n=",n,"; Nhat.yx=",signif(Nhat.yx,3),";``bias''=",signif(bias,3),"%\n",sep="")
    message(msg)
    if(plot) {
      mtext(msg,cex=0.8)
    }
  }else{N=N;n=NULL;Nhat=NULL;bias=NULL}
  invisible(list(gridx=gridx,p.xpifit=p.xpifit,mufit=mufit,
                 f.xfit=f.xfit,p.xfit=p.xfit,ptot=ptot,p.xfit.std=p.xfit.std,adbn=adbn,
                 N=N,n=n,Nhat=Nhat.yx,bias=bias))
}


#'@title Plot fitted f(y) and forward distance distribution
#'
#'@description Plot f(y) and forward distance distribution
#'resulting from a call of \code{\link{fityx}}.
#'
#'@param y forward distance observations (if NULL, uses est$dat$y)
#'@param x perpendicular distance observations (if NULL, uses est$dat$x)
#'@param est return from a call of \code{\link{fityx}}
#'@param nclass number of histogram classes
#'@param breaks break points passed to hist (overrides nclass if not NULL)
#'@param plot boolean, plot results
#'@param lineonly if TRUE plots only f(y), else plots histogram of forward distances too
#'@param nint number of intervals to use in calculating f(y)
#'@param max.obs If TRUE, plots only up to maximum observed forward distance, else plots
#'up to est$ystart (the forward distance beyond which detection is assumed impossible).
#'@param add if TRUE, adds line to existing plot, else creates new plot. Only applicable
#'if lineonly==TRUE.
#'@param ... other parameters passed to \code{hist} and \code{plot} (need to separate these two!)
#'
#'@details Plot f(y) and forward distance distribution resulting from a call of
#'\code{\link{fityx}}, optionall with overlaid histogram of foward distances.
#'Invisibly returns various f(y)'s, as detailed below.
#'
#'@return
#'Invisibly returns a list with these elements
#'\code{$gridx} = x values used in plotting
#'\code{$fy.x} = Unscaled f(y|x) for all the xs observed
#'\code{$fy.} = Unscaled mean of f(y|x) for all the xs observed
#'\code{$scaled. fy.} = Mean of f(y|x), scaled to integrate to 1
#'@examples
#'\dontrun{
#'ystart=4;w=1
#'hr=h2; b=log(c(0.75,1))
#'pi.x=pi.norm; logphi=c(0.5,log(0.2))
#'N=50 #true number of animals
#'#generate some observations
#'simDat=simXY(N=N,pi.x=pi.x,logphi=logphi,
#'hr=hr,b=b,w=w,ystart=ystart)
#'x=simDat$locs$x; y=simDat$locs$y
#'est.yx=fityx(y,x,b,hr,ystart,pi.x,logphi,w)
#'plotfit.y(y,x,est.yx,nclass=10)
#'}
#'@seealso \code{\link{fityx}}
#'@export
plotfit.y=function(est,nclass=10,breaks=NULL,plot=TRUE,dotitle=FALSE,
                   lineonly=FALSE,nint=100,max.obs=TRUE,add=FALSE,
                   y=NULL,x=NULL,...)
{
  # Some type-checking to ensure we can extract the data in the usual way:
  if (class(est)!='LT2D.fit.object')            # We want the fit object,
    #& y==NULL & x==NULL)                       # or data passed by plot.smoothfy
  {stop('Can only plot LT2D objects')}

  # barely made changes to this function (just type checking), worked fine - Cal

  b=est$b; hr=est$hr; ystart=est$ystart; pi.x=est$pi.x
  logphi=est$logphi; w=est$w
  if(is.null(y)) y=est$dat$y
  if(is.null(x)) x=est$dat$x
  # calculate stuff to plot:
  n=length(y)
  res=100
  gridy=seq(1e-10,ystart,length=res)
  fy.x=matrix(rep(NA,n*res),nrow=n)
  for(i in 1:n) {
    fy.x[i,]=fyx(gridy,rep(x[i],res),b,hr,ystart,nint=nint)
  }

  fy.=apply(fy.x,2,mean)
  fy.area=sum((fy.[-1]+fy.[-length(fy.)])/2*diff(gridy))
  scaled.fy.=fy./fy.area
  if(plot){
    ymax=ystart
    if(max.obs) ymax=max(y)
    if(is.null(breaks)) breaks=seq(min(y,1e-10),ymax,length=(nclass+1))
    fy.area=sum((fy.[-1]+fy.[-length(fy.)])/2*diff(gridy))
    scaled.fy.=fy./fy.area
    if(lineonly) {
      if(add) lines(gridy,scaled.fy.,...)
      else plot(gridy,scaled.fy.,ylim=c(0,max(scaled.fy.)),type="l",
                xlab="forward distance (y)",ylab="f(y)",...)
    }
    else {
      # hst=hist(y,plot=FALSE)
      # hist(y,freq=FALSE,xlab="forward distance (y)",nclass=nclass,ylim=c(0,max(hst$intensities,fy.)))
      hst=hist(y,breaks=breaks,plot=FALSE,...)
      # cat("hist area=",hst$desity*diff(hst$breaks),"\n")
      hmax=max(scaled.fy.,hst$density)
      if(dotitle) hist(y,freq=FALSE,xlab="forward distance (y)",breaks=breaks,ylim=c(0,hmax),...)
      else hist(y,freq=FALSE,xlab="forward distance (y)",breaks=breaks,ylim=c(0,hmax),main="",...)
      lines(gridy,scaled.fy.,...)
      # cat("fy area=",sum((fy.[-1]+fy.[-length(fy.)])/2*diff(gridy)),"\n")
    }}

  invisible(list(gridy=gridy,fy.x=fy.x,fy.=fy.,scaled.fy.=scaled.fy.))
}


#'@title Plot smooth fitted f(y) and forward distance distribution for small xs
#'
#'@description Plot spline smooth of f(y) and forward distance distribution
#'resulting from a call of \code{\link{plotfit.y}}.
#'
#'@param fit fitted object output by \code{\link{fityx}}
#'@param nclass number of histogram bins to use
#'@param nfys number of points to use in plotting smooth
#'@param xmax maxumum perp. dist. to use
#'
#'@details Plot f(y) and forward distance distribution resulting from a call of
#'\code{\link{fityx}}. This is a post-hoc fix of \link{\code{plotfit.y}}, which
#'produces f(y) with some sharp and implausible bends.
#'
#'@return
#'Invisibly returns a list with these elements
#'hst=hst,y=ys,smfy=smfy
#'\code{$hst} = histogram object from call to \link{\code{hist}}, containing data for
#'histogram of the detections within perp. dist \code{xmax}.
#'\code{$y} = y values for plot of smooth f(y)
#'\code{$smfy} = smooth f(y)
#'
#'@seealso \code{\link{fityx}}
#'@export
plotfit.smoothfy=function(fit,nclass=12,nfys=200,xmax=max(fit$dat$x),main="",
  plot=TRUE) {

  if (class(fit)!='LT2D.fit.object'){stop('Can only plot LT2D objects')}

  near0=which(fit$dat$x<=xmax)
  ys=fit$dat$y[near0]
  ymax=max(ys)
  breaks=seq(0,ymax,length=(nclass+1))
  fy=plotfit.y(fit,nclass=nclass,nint=100,plot=FALSE,y=ys,x=fit$dat$x[near0])
  sm=splinefun(fy$gridy,fy$scaled.fy.,method="monoH.FC")
  hst=hist(ys,breaks=breaks,plot=FALSE)
  xs=seq(0,ymax,length=nfys)
  smfy=sm(xs) # smooth of curve
  ymax=max(smfy,hst$density)

  hist(ys,xlab="Forward distance (y)",ylab="Density",breaks=breaks,main=main,
    freq=FALSE,ylim=c(0,ymax))

  lines(xs,smfy)
  invisible(list(hst=hst,y=ys,smfy=smfy))
}

#'@title Plot fitted curves of an LT2D.fit.object
#'
#'@description Produces plots of fitted functions, pi, p and f for perpendicular
#'  distance x, as well as a plot of fitted fy for forward distance y. Calls
#'  \link{\code{plotfit.y}}, \link{\code{plotfit.x}} and
#'  \link{\code{plotfit.smothfy}} to produce the figures..
#'
#'@param fit fitted object output by \code{\link{fityx}}
#'@param xbins int;number of histogram bins to use for perpendicular distance
#'@param ybins int;number of histogram bins to use for forward distance

#'@return invisibly returns the outputs of the three functions to which it
#'delagates tasks
#'
#'@seealso \code{\link{fityx}}
#'@export
plot.LT2D.fit.object = function(fit,
  xbins=20,
  ybins=20,
  smooth.fy=FALSE,
  addrug=FALSE){

  # Wrapper function which calls the appropriate functions to plot LT2D fits.

  # type-check input:
  if (class(fit)!='LT2D.fit.object'){stop('Can only plot LT2D objects')}

  X = plotfit.x(fit,nclass=xbins)
  if (addrug){rug(x[x<=w])}                        # perpendicular distance plot

  if (smooth.fy){
    # plotfit.smoothfy truncates troublesome data
    # before calling the appropriate function,
    # to produce a smoother plot
    Y = plotfit.smoothfy(fit,nclass = ybins)       # forward distance plot
    if (addrug){rug(x=y[x<=w])}
  }
  else{ # We plot without calling .smoothfy
    plotfit.y(fit,nclass=ybins)
    Y = if (addrug){rug(x=y[x<=w])}
  }
}

#'@title Plot LT2D fit
#'
#'@description acts as a wrapper for the appropriate
#' plotting function.
#'
#'@param fit fitted object output by \code{\link{LT2D.fit}}
#'@param ..., parameters to be passed to \code{\link{plot.LT2D.fit.object}}
#'@seealso \code{\link{LT2D.fit}}
#'@export
plot.LT2D.fit.function.object = function(fit, ...){
  plot(fit$fit, ...) # this rather cryptic line extracts
  # the fit object created by fit.yx from the one
  # created by LT2D.fit, and calls plot, which knows
  # how to deal with it correctly, using the above
  # function.
}


#'@title Calculate coverage probabilities of \eqn{\hat p} for simulated data
#'
#'@description Calculate coverage probabilities of \eqn{\hat p} using the delta
#'method and assuming a log-normal error distribution.
#'
#'@param fit object resulting from a call of \code{\link{fityx}} (see details)
#'@param interval the interval used to determine coverage probability
#'@param true.hr true form of hazard rate function
#'@param true.b true values of hazard rate parameters
#'@param true.pi.x true form of perpendicular density distribution pi(x)
#'@param true.logphi true values for pi(x) perpendicular density distribution parameters
#'@param verbose boolean.  FALSE only covered returned.  TRUE a data frame of covered and associated calculations returned (see returns).
#'@param type \code{LOGNORM} log-normal confidence intervals; \code{NORM} normal confidence intervals.
#'@details In the call of \code{\link{fityx}} must have \code{hessian=TRUE}
#'@return boolean; \code{TRUE} if within log-normal confidence interval, otherwise \code{FALSE}.
#'\code{verbose=TRUE}   a data frame of p phat, and CV[phat], interval, lower bound returned.
#'@examples
#'\dontrun{
#'ystart=4;w=1
#'hr=h2; b=log(c(0.75,1))
#'pi.x=pi.norm; logphi=c(0.5,log(0.2))
#'N=50 #true number of animals
#'#generate some observations
#'simDat=simXY(N=N,pi.x=pi.x,logphi=logphi,
#'hr=hr,b=b,w=w,ystart=ystart)
#'x=simDat$locs$x; y=simDat$locs$y
#'est.yx=fityx(y,x,b,hr,ystart,pi.x,logphi,w)
#'coveragep(fit=est.yx,true.hr=hr,true.b=b,interval=0.95,
#'  true.pi.x=pi.x,true.logphi=logphi,verbose=TRUE)
#'}
#'@seealso \code{\link{phat}} \code{\link{fityx}}
#'@export
coveragep=function(fit,true.hr,true.b,true.pi.x,true.logphi,type='LOGNORM',
                   interval=0.95,verbose=FALSE){
  lnci.nmin=function(stat,cv,stat.min=0,interval=interval){
    q=Mod(qnorm((1-interval)/2,0,1))
    varNhat=(stat*cv)^2
    cfactor=exp(q*sqrt(log(1+varNhat/(stat-stat.min)^2)))
    lower=stat.min+(stat-stat.min)/cfactor
    upper=stat.min+(stat-stat.min)*cfactor
    return(list(lower=lower,upper=upper))
  }

  if(!'hessian' %in% names(fit))
    stop('fit ARG must include a hessian matrix')
  pars=fit$par;
  hr=match.fun(fit$hr); b=fit$b;
  ystart=fit$ystart; w=fit$w
  pi.x=fit$pi.x; logphi=fit$logphi
  #true p
  p=phat(w=w,hr=true.hr,b=true.b,ystart=ystart,pi.x=true.pi.x,logphi=true.logphi)
  #estimated p
  p.hat=phat(w=w,hr=hr,b=b,ystart=ystart,pi.x=pi.x,logphi=logphi)
  #calc. variance-covariance
  vcov=solve(fit$hessian)
  #Implement the delta method:
  #numerical differentiation
  dbyd=numericDeriv(quote(phat(w=w,b=b,hr=hr,ystart=ystart,pi.x=pi.x,logphi=logphi)), c("b","logphi"))
  dbyd=as.vector(slot(dbyd,'gradient'))
  var.p.hat=as.vector(t(dbyd)%*%vcov%*%dbyd) #$var[hat{p}(\hat{\Beta})$
  if(type=='LOGNORM')
    bounds=lnci.nmin(stat=p.hat,cv=sqrt(var.p.hat)/p.hat,interval=interval)
  if(type=='NORM'){
    bounds=list(lower=qnorm((1-interval)/2,p.hat,sqrt(var.p.hat)),
                upper=qnorm(interval+(1-interval)/2,p.hat,sqrt(var.p.hat)))
  }
  if(any(sapply(bounds,is.nan))){
    warning('One or both p.hat bounds NaN')
    return(data.frame(covered=NA,p=p,phat=p.hat,CV.phat=NA,
                      interval=interval,lower.bound=bounds$lower,upper.bound=bounds$upper))
  }
  covered=TRUE
  if(p < bounds$lower | p > bounds$upper) covered=FALSE
  if(!verbose)
    return(covered)
  return(data.frame(covered=covered,p=p,phat=p.hat,CV.phat=sqrt(var.p.hat)/p.hat,
                    interval=interval,lower.bound=bounds$lower,upper.bound=bounds$upper))
}

#'@title Calculate effective strip width
#'
#'@description Calculate effective strip width, \eqn{\hat p} for a given hazard rate function and
#'perpendicular density distribution.
#'
#' @param w truncation distance
#' @param hr function describing the hazard rate
#' @param b hazard rate parameter vector
#' @param ystart maximum forward distance
#' @param pi.x function describing the perpendicular distance density distribution
#' @param logphi parameters for pi.x (some maybe logged)
#' @param fit=NULL alternatively just pass in an object resulting from a call of \link{fityx}.
#'@return Effective strip widith \eqn{\hat p}
#'@examples phat(w=1,hr=h2,b=log(c(0.75,1)),ystart=4,pi.x=pi.norm,logphi=c(0.5,log(0.2)))
#'@export
phat=function(fit=NULL,w=NULL,hr=NULL,b=NULL,ystart=NULL,pi.x=NULL,logphi=NULL)
{
  if (!is.null(hr)){
    if (class(hr)!='character'){
      stop('phat: hr must be a character')}
    hrname = hr
  }

  if (!is.null(pi.x)){
    if (class(pi.x)!='character'){
      stop('phat: pi.x must be a character')}
    piname = pi.x
  }
  if(!is.null(fit)){
    #f=fit$p.pi.x;
    upper=fit$w   ; b=fit$b ; hrname=fit$hr ; ystart=fit$ystart
    piname=fit$pi.x ; logphi=fit$logphi ; w=fit$w
  }

  int=integrate(f=p.pi.x,lower=0,upper=w,b=b,hr=hrname,
                ystart=ystart,pi.x=piname,logphi=logphi,w=w)$value
  return(int)
}

#--------------- Functions added by DLB 24/7/14 ------------------------------

#' @title Simulate n sightings from NHPP
#'
#' @description Simulates n sighting locations (x,y) given a perp.dist
#' distribution model and detection hazard model
#'
#' @param n sample size
#' @param ymin smallest forward distance
#' @param ymax largest forward distance
#' @param W perpendicular truncation distance
#' @param hfun detection hazard function
#' @param b vector of detection hazard function parameters
#' @param pi.x perpendicular distance distribution function
#' @param logphi vector with log of pi.x parameters
#' @param fix.n if TRUE sample size of exactly n is generated, else sample is
#'        generated from model with expected sample size n
#' @param intscale amount by which to multiply detection location pdf f(x,y)
#'        in order to get required sample size. Either an object of class
#'        "ppscale" output by \code{\link{calc.lpars}} or NULL (in which case
#'        \code{\link{calc.lpars}} is called inside \code{sim.n}.
#' @param nbuffer amount by which to multiply the expected n by (given all model
#'        parameters and intscale) to reduce probability that generated n is less
#'        that n on first call to NHPP generating funciton rpoispp. If NULL,
#'        it is set to 1.25 inside \code{sim.n}.
#' @details Uses the \code{spatstat} function \code{rpoispp} to
#' generate detections from a NHPP, with intensity parameter such that the
#' expected (if fix.n is FALSE) or actual (if fix.n is TRUE) sample size is n.
#'
#' @return a list object with element 1 a data.frame of simulated \code{x} and \code{y} sightings locations; element 2 \code{spatstat} object of class "ppp" with x- and y-coordinates
#' of detections and element 2 simulation settings, comprising of \code{n}
#'@examples
#' \dontrun{
#' # simulate with fixed n:
#' n=100;ymin=0.01;ymax=5;W=2
#' b=log(c(0.75,1));logphi=c(0.5,log(0.3))
#' dat=sim.n(n,ymin,ymax,W,h2,b,pi.norm,logphi)
#' dat$locs$n
#'
#' plot(density(dat))
#' contour(density(dat),add=TRUE)
#' plot(dat,pch="+",cex=0.75,add=TRUE)
#' hist(abs(dat$y),xlab="Perpendicular distance",main="")
#' hist(dat$x,xlab="Forward distance",main="")
#' # do same with random n:
#' dat=sim.n(n,ymin,ymax,W,h2,b,pi.norm,logphi,fix.n=FALSE)
#'
#' # compare time taken if calculate intscale on the run vs pass it:
#' # first calculate each time:
#' system.time(for(i in 1:20) dat<-sim.n(n,ymin,ymax,W,h2,b,pi.norm,logphi))
#' # then calculate once and pass:
#' intscale=calc.lpars(n,ymin,ymax,W,h2,b,pi.norm,logphi)
#' system.time(for(i in 1:20) dat<-sim.n(n,ymin,ymax,W,h2,b,pi.norm,logphi,intscale=intscale))
#' }
#' @export
sim.n=function(n,ymin,ystart,w,hr,b,pi.x,logphi,fix.n=TRUE,intscale=NULL,nbuffer=NULL){
  simDat=list()
  simDat[[3]]=list(n=n,ymin=ymin,ystart=ystart,b=b,hr=hr,pi.x=pi.x,logphi=logphi,w=w,
                   fix.n=fix.n,intscale=intscale,nbuffer=nbuffer)
  names(simDat)[3]='settings'
  #I'd like to keep the function ARGS to keep them in line with the other functions-
  #- but I don't want to mess with the args inside the function.  So:
  # change args to David's
  ymax=ystart; W=w; hfun=hr

  # calculate scaling needed to give E[n]=n
  if(is.null(intscale))
    intscale=calc.lpars(n,ymin,ymax,W,hfun,b,pi.x,logphi)
  if(class(intscale)!="ppscale") stop("intscale must be class `ppscale' (output of calc.lpars).")
  window=owin(c(0,ymax),c(-W,W)) # create observation window
  lmax=intscale$lmax # maximum value of intensity function on grid used by calc.lpars
  lscale=intscale$lscale # multiplier required to get intensity with expected sampls size n
  if(is.null(nbuffer)){
    # increase E[n] by 25% to reduce prob that sample size is < n:
    nbuffer=ifelse(fix.n,1.25,1)
  }
  pp=rpoispp(poisint,lmax,window,ymin=ymin,ymax=ymax,b=b,hfun=hfun,
             pi.x=pi.x,logphi=logphi,W=W,lscale=lscale*nbuffer)
  if(fix.n) {
    while(pp$n<n){ # crude way of generating big enough sample size:

      pp2=rpoispp(poisint,lmax,window,ymin=ymin,ymax=ymax,b=b,hfun=h2,
                  pi.x=pi.norm,logphi=logphi,W=W,lscale=lscale*nbuffer)
      print('RPOIS PP SUCCESS')
      pp$n=pp$n+pp2$n
      pp$x=c(pp$x,pp2$x)
      pp$y=c(pp$y,pp2$y)
    }
    pp$n=n
    pp$x=pp$x[1:n]
    pp$y=pp$y[1:n]
  }
  x=abs(pp$y);y=pp$x
  simDat[[1]]=data.frame(x=x,y=y);names(simDat)[1]='locs'
  simDat[[2]]=pp;names(simDat)[2]='pp'
  return(simDat)
}




#' @title Calculates scaling parameters required by \code{sim.n}
#'
#' @description Calculates scaling parameters required by \code{sim.n} in
#'              order to generate sample of given size.
#'
#' @param n expected sample size reqired from NHPP
#' @param ymin smallest forward distance
#' @param ymax largest forward distance
#' @param W perpendicular truncation distance
#' @param hfun detection hazard function
#' @param b vector of detection hazard function parameters
#' @param pi.x perpendicular distance distribution function
#' @param logphi vector with log of pi.x parameters
#' @param nx number of intervals on x-axis at which to calculate intensity
#' @param ny number of intervals on y-axis at which to calculate intensity
#' @param inflate: multiplier by which to increase max intensity in region; this is
#'        just a safeguard against the max at some point not calculated being
#'        greater than the points at which it was calculated. (The sample
#'        generator in \code{rpoispp} uses rejection sampling so needs the
#'        global max.)
#'
#' @details Calls \code{\link{poisint}} on grid of (x,y) points, calculates total
#' intensity in area and then calculates scaling needed to make this equal to
#' n, and the maximum intensity in the region.
#'
#' Output from this function can be passed as the \code{lscale} argument to
#' \code{\link{sim.n}}.
#'
#' @return Object of class "ppscale", being alist with two parameters:
#' \code{$lscale} is the multiplier needed to make the total intensity equal
#' to n, while \code{$lmax} is maximum lintensity in the region.
#'
#'@examples
#' \dontrun{
#' n=100;ymin=0.01;ymax=5;W=2
#' b=log(c(0.75,1));logphi=c(0.5,log(0.3))
#' intscale=calc.lpars(n,ymin,ymax,W,h2,b,pi.norm,logphi)
#' intscale
#' }
#' @export
calc.lpars=function(n,ymin,ymax,W,hfun,b,pi.x,logphi,nx=100,ny=100,inflate=1.05){
  x=seq(-W,W,length=nx)
  y=seq(ymin,ymax,length=ny)
  a=diff(x)[1]*diff(y)[1] # grid cell area
  lambda=outer(y,x,poisint,ymin=ymin,ymax=ymax,hfun=hfun,b=b,pi.x=pi.x,logphi=logphi,W=W)
  En=sum(lambda*a)
  lscale=n/En
  lmax=max(lambda*lscale)*inflate # bigger than observed in case higher between grid cells
  outlist=list(lscale=lscale,lmax=lmax)
  class(outlist)="ppscale"
  return(outlist)
}


#' @title NHPP intensity calculation in region
#'
#' @description Calculates NHPP intensity in given region, using
#'              given detection hazard and perp. dist. distribution.
#' @param y forward distances at which to calculate intensities (vector)
#' @param x perpendicular distances at which to calculate intensities (must
#'          be same length as y).
#' @param ymin smallest forward distance
#' @param ymax largest forward distance
#' @param W perpendicular truncation distance
#' @param hfun detection hazard function
#' @param b vector of detection hazard function parameters
#' @param pi.x perpendicular distance distribution function
#' @param logphi vector with log of pi.x parameters
#' @param lscale output of \code{\link{calc.lpars}} (object of class ``ppscale'')
#'
#' @details Calculates survival model pdf of forward distance \code{y}, given
#' perpendicular distance \code{x} (\code{f(y|x).}), multiplies this by the
#' perpendicular distance pdf \code{pi.x} and then scales it by multiplying
#' by lscale (in order to get some total intensity - typically that to
#' generate some expected sample size). See \code{\link{calc.lpars}} for
#' details of the scaling.
#'
#' This function is called by \code{rpoispp} inside \code{sim.n} to generate
#' samples from NHPPs.
#'
#' @return NHPP intensities at all \code{(x,y)}s input.
#'
#'@examples
#' \dontrun{
#' n=100;ymin=0.01;ymax=5;W=2
#' b=log(c(0.75,1));logphi=c(0.5,log(0.3))
#' nf=100
#' ys=seq(ymin,ymax,length=nf)
#' intscale=calc.lpars(n,ymin,ymax,W,h2,b,pi.norm,logphi)
#' f=poisint(ys,rep(0,nf),ymin=ymin,ymax=ymax,hfun=h2,b=b,pi.x=pi.norm,logphi=logphi,W=W,lscale=intscale$lscale)
#' plot(ys,f,type="l",xlab="Forward distance (y)",ylab="f(y)")
#' }
#' @export
poisint=function(y,x,ymin,ymax,hfun,b,pi.x,logphi,W,lscale=1){
  if (!class(hfun)=="character"){stop('Message from poisint:
                                      hfun must be supplied as character')}
  h=match.fun(hfun)
  pix=match.fun(pi.x)
  nx=length(x)
  ax=abs(x)
  if(length(y)!=nx) stop("Lengths of x and y must be the same.")
  pxx=pix(ax,logphi,W)
  f=p=p0=rep(NA,nx)
  for(i in 1:nx){
    #    p0[i]=1-Sy(ax[i],ymin,ymax,b,hfun)
    p[i]=1-Sy(ax[i],y[i],ymax,b,hfun)
    f[i]=h(y[i],ax[i],b)*(1-p[i])
  }
  #  return(list(f=f,p=p,p0=p0))
  return(f*pxx*lscale)
}


#' @title Calculates survivor function
#'
#' @description Calculates survivor function to forward distance \code{y},
#' for given perpendicular distance \code{x} and given forward distance range.
#'
#' @param x perpendicular distance (scalar)
#' @param y forward distance  (scalar)
#' @param ymax largest forward distance
#' @param hfun detection hazard function
#' @param b vector of detection hazard function parameters
#' @details Calculates probability of making it from forward distance
#' \code{ymax} to \code{y} without being detected.
#'
#' @return Probability of making it from forward distance
#' \code{ymax} to \code{y} without being detected.
#'
#'@examples
#' \dontrun{
#' ymax=5
#  b=log(c(0.75,1))
#' Sy(0,0.1,ymax,b,h2)
#' }
#' @export
Sy=function(x,y,ymax,b,hr) {

  if (!class(hr)=='character'){stop('message from Sy: hr must be passed as a
                                    character')}
  n=length(x)
  if(length(y)!=n) stop("Lengths of x and y must be the same.")
  pS=rep(NA,n)

  if(hr=="h1") { # Hayes & Buckland hazard rate model, so can do analytically
    hmax=h1(y,x,b)
    for(i in 1:n){
      if(y[i]==0 | hmax[i]>1e10){ # computer will think integral divergent for large hmax
        pS[i]=1-HBhr(x[i],h1.to.HB(b))
      } else {
        pS[i]=exp(-integrate(match.fun(hr),y[i],ymax,x=x[i],b=b,subdivisions = 1000L)$value)
        pS<<-pS[i]
      }
    }
  } else { # Not Hayes & Buckland hazard rate model, so can't do analytically
    for(i in 1:n){
      pS[i]=exp(-integrate(match.fun(hr),y[i],ymax,x=x[i],b=b,subdivisions = 1000L)$value)
      pS2<<-pS[i]
    }
  }
  return(pS)
}


#' AIC-based model selection for models fitted using \link{fitxy}
#'
#'AIC-based model selection of a list of models fitted using \link{fitxy}.  Parameter
#'estimates for the hazard and perpendicular distribution functions along with
#'associated coefficients of variation are also provided.
#'Models are ranked in order of increasing dAIC.  The model list can also be
#'returned in a data frame suitable making a tex-type table.
#'@param modList a list of \link{fitxy}-type models.
#'@param modNames a vector of model names.
#'@param tab boolean - return a data frame suitable for generating a table for report or paper.
#'@param digits - number of digits of rounding (see \link{round}) in the table
#'@param maxblength=NULL; maximum number of parameters for the hazard function.  If NULL
#'the maximum number of parameters for the models in modList is used.
#'@param maxlogphilength=NULL; maximum number of parameters for the perpendicular density function.  If NULL
#'the maximum number of parameters for the models in modList is used.
#'@export
#'@return tab=FALSE list; element 1 =data frame of AIC-ranked models, with parameter estimates,
#'coefficients of variation and AIC; element 2= vector of the order of AIC-ranked models. tab=TRUE, a three element list with AIC-ranked model order and two tables, the first as tab=FALSE, the
#'second, a table for use in reports or manuscripts.
#'@seealso \link{fitxy}
#'@examples
#'\dontrun{
#'ystart=4;w=1
#'hr=h2; b=log(c(0.75,1))
#'pi.x=pi.norm; logphi=c(0.5,log(0.2))
#'N=50 #true number of animals
#'#generate some observations
#'simDat=simXY(N=N,pi.x=pi.x,logphi=logphi,
#'hr=hr,b=b,w=w,ystart=ystart)
#'x=simDat$locs$x; y=simDat$locs$y
#'fit1=fityx(y,x,b,hr,ystart,pi.x,logphi,w)
#'fit2=fityx(y,x,b=log(c(0.001,1)),hr=h1,ystart,pi.x,logphi,w)
#'modSelect(modList=list(fit1,fit2),modNames=c('h1','h2'),tab=TRUE)
#'}
modSelect=function(modList,modNames=NULL,tab=FALSE,digits=2,
                   maxblength=NULL,maxlogphilength=NULL)
{
  if(is.null(maxblength))
    maxblength=max(sapply(modList,function(x) length(x$b)))
  if(is.null(maxlogphilength))
    maxlogphilength=max(sapply(modList,function(x) length(x$logphi)))

  fitVal=function(fit,maxblength=maxblength,
                  maxlogphilength=maxlogphilength)
  {
    bhat=CVbhat=rep(NA,maxblength)
    logphi=CVlogphi=rep(NA,maxlogphilength)
    bhat[1:length(fit$b)]=fit$b[1:length(fit$b)]
    CVbhat[1:length(fit$b)]=fit$CVpar[1:length(fit$b)]
    if(all(!is.na(fit$logphi)))
    {
      valLoc=which(!is.na(fit$logphi))
      logphi[valLoc]=fit$logphi[valLoc]
      CVlogphi[valLoc]=fit$CVpar[length(fit$b)+valLoc]
    }
    out=c(bhat,CVbhat,logphi,CVlogphi,length(fit$par),fit$value,fit$AIC)
    names(out)=c(paste('b',1:maxblength,sep=''),paste('CVb',1:maxblength,sep=''),
                 paste('logphi',1:maxlogphilength,sep=''),paste('CVlogphi',1:maxlogphilength,sep=''),
                 'n','logLik','AIC')
    return(out)
  }
  tab1=t(sapply(modList,fitVal,maxblength=maxblength,
                maxlogphilength=maxlogphilength))
  out1=data.frame(tab1)
  names(out1)=colnames(tab1)
  rownames(out1)=modNames
  minAICLoc=which.min(out1$AIC)
  out1$dAIC=out1$AIC-out1$AIC[minAICLoc]
  AICorder=order(out1$dAIC)
  out1=out1[AICorder,]
  ww=exp( -0.5 * out1$dAIC)
  out1$w=ww/sum(ww)
  if(tab)
  {
    bhat=logphihat=vector(mode='character',length=nrow(out1))
    for(i in 1:nrow(out1)){
      bhat[i]=paste(paste(round(tab1[i,1:maxblength],digits),'(',round(tab1[i,(maxblength+1):(2*maxblength)],digits),')',sep=''),collapse='; ')
      logphihat[i]=paste(paste(round(tab1[i,(2*maxblength+1):(2*maxblength+maxlogphilength)],digits),
                               '(',round(tab1[i,(2*maxblength+maxlogphilength+1):(2*maxblength*maxlogphilength)],digits),')',sep=''),collapse='; ')
    }
    out2=data.frame(bhat=bhat,logphihat=logphihat)
    out2=out2[AICorder,]
    out2$n=out1$n; out2$logLik=out1$logLik; out2$AIC=out1$AIC
    rownames(out2)=rownames(out1)
    out2$dAIC=out1$dAIC; out2$w=out1$w
    out2[,c('logLik','AIC','dAIC','w')]=round(out2[,c('logLik','AIC','dAIC','w')],digits)
    return(list(AICorder=AICorder,res=out1,tab=out2))
  }
  return(list(AICorder=AICorder,res=out1))
}

#'@title Calculates coverage probabilities of \eqn{\hat p}
#'
#'@description Calculate coverage probabilities of \eqn{\hat p} using the delta method and assuming a
#'log-normal error distribution.
#'
#'@param fit object resulting from a call of \code{\link{fityx}} (see details)
#'@param interval the interval used to determine coverage probability
#'@param type \code{LOGNORM} log-normal confidence intervals; \code{NORM} normal confidence intervals.
#'@details In the call of \code{\link{fityx}} in the \code{fit} argument must have \code{hessian=TRUE}
#'@return a data frame of p phat, and CV[phat], interval, lower bound returned.
#'@examples
#'\dontrun{
#'ystart=4;w=1
#'hr=h2; b=log(c(0.75,1))
#'pi.x=pi.norm; logphi=c(0.5,log(0.2))
#'N=50 #true number of animals
#'#generate some observations
#'simDat=simXY(N=N,pi.x=pi.x,logphi=logphi,
#'hr=hr,b=b,w=w,ystart=ystart)
#'x=simDat$locs$x; y=simDat$locs$y
#'est.yx=fityx(y,x,b,hr,ystart,pi.x,logphi,w)
#'phatInterval(fit=est.yx,interval=0.95)
#'}
#'@seealso \code{\link{phat}} \code{\link{fityx}}
#'@export
phatInterval=function(fit,type='LOGNORM',
                      interval=0.95){
  lnci.nmin=function(stat,cv,stat.min=0,interval=interval){
    q=Mod(qnorm((1-interval)/2,0,1))
    varNhat=(stat*cv)^2
    cfactor=exp(q*sqrt(log(1+varNhat/(stat-stat.min)^2)))
    lower=stat.min+(stat-stat.min)/cfactor
    upper=stat.min+(stat-stat.min)*cfactor
    return(list(lower=lower,upper=upper))
  }

  test = 'test'

  if(!'hessian' %in% names(fit))
    stop('fit ARG must include a hessian matrix')
  pars=fit$par;
  hr=fit$hr; b=fit$b;
  ystart=fit$ystart; w=fit$w
  pi.x=fit$pi.x; logphi=fit$logphi
  #estimated p
  p.hat=phat(fit=fit)
  #variance-covariance matrix
  vcov=fit$vcov
  #Implement the delta method:

  #numerical differentiation
  if(!is.numeric(fit$logphi))
  {
    dbyd=numericDeriv(quote(phat(w=w,b=b,hr=hr,ystart=ystart,pi.x=pi.x,
                                logphi=logphi)), c("b"))  }
  else {
    dbyd=numericDeriv(quote(phat(w=w,b=b,hr=hr,ystart=ystart,pi.x=pi.x,
                                 logphi=logphi)), c("b","logphi")) }

  dbyd=as.vector(slot(dbyd,'gradient'))
  var.p.hat=as.vector(t(dbyd)%*%vcov%*%dbyd) #$var[hat{p}(\hat{\Beta})$
  if(type=='LOGNORM')
    bounds=lnci.nmin(stat=p.hat,cv=sqrt(var.p.hat)/p.hat,interval=interval)
  if(type=='NORM'){
    bounds=list(lower=qnorm((1-interval)/2,p.hat,sqrt(var.p.hat)),
                upper=qnorm(interval+(1-interval)/2,p.hat,sqrt(var.p.hat)))
  }
  if(any(sapply(bounds,is.nan))){
    warning('One or both p.hat bounds NaN')
    return(data.frame(covered=NA,p=p,phat=p.hat,CV.phat=NA, interval=interval,
                      lower.bound=bounds$lower,upper.bound=bounds$upper))
  }

  return(data.frame(phat=p.hat,CV.phat=sqrt(var.p.hat)/p.hat,interval=interval,
                    lower.bound=bounds$lower,upper.bound=bounds$upper))
}

#'Calculate \eqn{\hat p} and optionally \eqn{\hat N} for a list of models
#'
#'Calculate \eqn{\hat p} along with variance, \eqn{Var[\hat p]}, using the delta method.  Optionally \eqn{\hat N} can also be calculated.
#'@param modList list object of models created by \link{fityx}.
#'@param n=NULL number of animals detected. If n!=NULL \eqn{\hat N} is calculated
#'@param tab boolean - return a data frame suitable for generating a table for report or paper.
#'@param digits - number of digits of rounding (see \link{round}) in the table
#'@param ... arguments to be passed into \link{phatInterval}
#'@export
#'@return data frame with:
#'\code{phat} estimate of \eqn{\hat p}
#'\code{CV.phat} estimate of \eqn{CV[\hat p]}
#'\code{interval} confidence interval specified in the \code{interval} argument in \link{phatInterval}
#'\code{lower.bound} lower bound of \eqn{\hat p}
#'\code{upper.bound} upper bound of \eqn{\hat p}
#'and optionally if n!=NULL
#'\code{n} number of detected animals
#'\code{Nhat} estimated number of animals in covered region.
#'\code{NhatLower} lower bound of \eqn{\hat N}
#'\code{NhatUpper} upper boudn of \eqn{\hat N}
#'@export
#'@seealso \link{fityx} \link{phatInterval}

phatModels=function(modList,n=NULL,tab=FALSE,digits=2,...)
{
  phatTab=sapply(modList,function(x) phatInterval(fit=x,...),...)
  colName=row.names(phatTab)
  phatTab=as.data.frame(t(matrix(as.numeric(phatTab),nrow(phatTab),ncol(phatTab))),
                        row.names=colnames(phatTab))
  names(phatTab)=colName
  if(!is.null(n))
  {
    phatTab$n=n
    phatTab$Nhat=n/phatTab$phat
    phatTab$NhatLower=n/phatTab$upper.bound
    phatTab$NhatUpper=n/phatTab$lower.bound
  }
  if(!tab){
    return(phatTab)}else{
      phatV=vector(length=length(modList))
      for(i in 1:length(modList))
        phatV[i]=paste(paste(round(phatTab$phat[i],digits),
                             '(',round(phatTab$CV.phat[i],digits),')',sep=''),collapse='; ')
      tab=data.frame(phat=phatV,row.names=row.names(phatTab))
      if(!is.null(n)){
        NhatV=vector(length=length(modList))
        for(i in 1:length(modList))
          NhatV[i]=paste(paste(round(phatTab$Nhat[i],0),
                               '(',round(phatTab$NhatLower[i],0),',',
                               round(phatTab$NhatUpper[i],0),
                               ')',sep=''),collapse='; ')
        tab=cbind.data.frame(tab,Nhat=NhatV)
      }
      return(list(res=phatTab,tab=tab))
    }
}

#' Plot the 2D fit of a model
#'
#' Plot the 2D fit of a model resulting from a call of \link{fityx}
#' @param fit object resulting from a call of \link{fityx}
#' @param ... other parameters passed into \link{plotSim}
#' @details This function is a wrapper for \link{plotSim}.
#' @seealso \link{plotSim} \link{fityx}
plotFit=function(fit,...){
  obj=list(locs=data.frame(x=fit$dat$x,y=fit$dat$y),
           settings=list(pi.x=match.fun(fit$pi.x),
                         logphi=fit$logphi,
                         hr=fit$hr,
                         b=fit$b,
                         w=fit$w,
                         ystart=fit$ystart))
  plotSim(simDat=obj, nclass=10,xlab="perpendicular distance (x)",
          ylab="forward distance (y)",image=TRUE,...)
}

