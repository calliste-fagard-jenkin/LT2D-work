ParamNumRequired = function(DataFrameInput, Formula){
  # This function acts as a way for the user to verify how many parameters 
  # the formula they specified will require. It returns the number of rows
  # of the generated design matrix minus one, since the intercept is already 
  # specified in the vector b of start parameters.
  DM = DesignMatrix(DataFrameInput, Formula) # Handles formula LHS terms
  return(dim(DM)[2]-1)
}

HazardBCheck = function(y,x,b, HazardName){
  if (class(b)!='list'){stop('b must be a list')}
  
  if (length(b)!=HazardNumberLookup(HazardName)){stop('incorrect b length')}
  lx = length(x) ; ly = length(y)
  
  if (lx==1 & ly==1){anchor=1} # If they're both scalars, set anchor to 1
  
  # x and y can either be vectors or scalars, however, if they are both 
  # vectors, they should be the same length:
  else if (lx>1 & ly>1 & lx!=ly){
    stop('x and y vectors must be the same length')}
  
  # If only one of x and y is a vector, we need to use its length as the anchor
  # for the required length of beta entries, if they are also vectors:
  else if (lx>1){anchor=lx}
  else if (ly>1){anchor=ly}
  else {anchor = 1} # x and y have length one - single data point, so LPs 
  # must be of length 1
  
  # Now we need to check that only slots which are allowed to have linear 
  # predictors are vectors...
  SlotsAllowed = HazardCovarSlots(HazardName)
  
  for (index in (1:length(b))){           # For each slot
    beta = b[[index]]
    lb = length(beta)
    # cat(beta)
    # cat('lb:',lb)
    # cat(' Slots Allowed:',SlotsAllowed)
    # cat(' index',index,'\n\n')
    if (lb>1 &                            # if its entry is a vector and the 
        !(index %in% SlotsAllowed)){      # slot isn't meant to be an LP...
      # If it's meant to be a slot with no covars, ensure that any vectorisation
      # has preserved the uniqueness of the values in that column 
      if (!length(unique(beta))==1){ 
      stop('Linear predictor present in illegal slot for specified hazard')}
    }
    
    if (lb>1 & lb!=anchor){ # if it's a vector, it must match the data set dim
      stop('incorrect linear predictor length for supplied data')}
  }
}

HazardNumberLookup = function(HazardName){
  # Not super elegant, but it seems R doesn't have dictionaries?
  
  if (class(HazardName)!='character'){
    stop('HazardName must be a character')
  }
  if      (HazardName %in% c('h1','h2','ip0','h.okamura')){return(2)}
  else if (HazardName %in% c('ghy','h21','ip1','ep1' )){return(3)}
  else if (HazardName %in% c('gh2','ip2','ep2')){return(4)}
  else if (HazardName=='h.const'){return(1)}
  else{stop('Hazard function not found')}
}

DensityNumberLookup = function(DensityName){
  if (class(DensityName)!='character'){
    stop('DensityName must be a character')
  }
  
  if      (DensityName=='pi.norm' | DensityName=='pi.chnorm'){return(2)} 
  else if (DensityName=='pi.const' | DensityName=='pi.hnorm'){return(1)}
  else{stop('Density function not found')}
}

HazardCovarsAllowed = function(HazardName){
  # Only the hazards in the 2016 paper have been made 
  # to allow covariate inclusion, for now. 
  
  # i - represents covars affecting intercept terms
  # x - represents covars affecting shape in x direction
  # y - represents covars affecting shape in y direction
  # xy - indicates the formulas for x and y must be the same 
  
  if (HazardName=='h1'){return(c('i'))}
  else if (HazardName=='ip1' | HazardName=='ep1'){return(c('i','x','y','xy'))}
  else if (HazardName=='ip2' | HazardName=='ep2'){return(c('i','x','y'))}
  else {return(FALSE)}  # Covariates not supported yet
}

HazardCovarSlots = function(HazardName){
  if (HazardName=='h1'){return(c(1))}
  else if (HazardName=='ip1' | HazardName=='ep1'){return(c(1,2))}
  else if (HazardName=='ip2' | HazardName=='ep2'){return(c(1,2,4))}
}

FormulaChecking = function(HazardName, Formulas){
  # We expect that the DataFrameInput has already been 
  # checked and sanitised by top level functions which call 
  # HazardWrapper.
  
  formulas.to.output = list(NULL) # Blank entries for now
  # We tag the object using a class name so that later, when it is passed to 
  # a different function, we can be sure that the object came from here
  class(formulas.to.output)='LT2D.verified.formulas' 
  
  hazard.function = match.fun(HazardName) # Check the hazard exists
  
  AllowedCovars = HazardCovarsAllowed(HazardName)
  # So the next function knows where to place the linear predictors:
  formulas.to.output[[4]]=AllowedCovars 
  
  if (class(AllowedCovars)!='character'){
    stop('Hazard function does not allow covariates')}
  
  # We want to make sure that the formulas the user has
  # specified lead to inclusion of covariates in a way which
  # is allowed by their chosen detection function:
  for (formula in Formulas){
    # formula[[2]] is the LHS, formula[[3]] is the RHS.
    # The formula LHS is saved as a name, so we convert it to a character
    # to make our comparison with the characters stored in AllowedCovars
    if (! as.character(formula[[2]]) %in% AllowedCovars){
      s1 = 'LHS of one or more specified formulae is not appropriate for'
      s2 = 'the chosen hazard function'
      stop(paste(s1,s2))
    }
    
    if (formula[[2]]=='i'){ # Save objects for later use
      i.formula = formula
      formulas.to.output[[1]]=i.formula # No checks needed, we can add it in now
      } 
    else if (formula[[2]]=='x'){x.formula = formula} # We need to check x and y 
    else if (formula[[2]]=='y'){y.formula = formula} # for consistency first
  } 
  
  if ('xy' %in% AllowedCovars){   # if x and y have to be the same:
    
    # if both have been specified, check they're valid:
    if (exists('x.formula') & exists('y.formula') & x[[3]]!=y[[3]]){
      y[[3]]=x[[3]]          # If the RHSs are different, we have a problem...
      s1 = 'x and y covariate formulas must be the same for the'
      s2 = 'specificed hazard. User chosen formulas were'
      s3 = 'different, or one was missinand so the x formula was used for both'
      s4 = 'x and y'
      warning(paste(s1,s2,s3,s4))
    }
    
    else if (exists('x.formula')|exists('y.formula')){ # only one of them exists
      
      warningStart = 'x and y formulas must be the same,'# To avoid typing again
      
      # if x is the defined one, set y equal to it:
      if (exists('x.formula')){
        y.formula=x.formula
        y.formula[[2]]=as.name('y') # change the LHS for consistency
        warning(paste(warningStart, 'x formula has been used for y'))
        }
      else{
        x.formula=y.formula
        x.formula[[2]]=as.name('x') # change the LHS for consistency
        warning(paste(warningStart, 'y formula has been used for x'))
      }
    }
    
    # This means no x or y formulas have been specified, so we return now
    else{return(formulas.to.output)} 
    
    # skipping the above else statement means we have passed through
    # all of the xy checks, and so we can add x and y formulas to the output
    # object and return it 
    formulas.to.output[[2]]=x.formula
    formulas.to.output[[3]]=y.formula
    return(formulas.to.output)
    
  }
   return(formulas.to.output)
}

DesignMatrix = function(DataFrameInput, formula){
  # We must remove the LHS of the formula, if it exists (it usually does):
  RHS = try(formula[[3]], silent=T)
  if (class(RHS)!='try-error'){ # This means there is a LHS as well as a RHS
    formula[[2]] = NULL         # This deletes the LHS and fixes the indices
  }
  # Make the design matrix for the whole data set:
    DM = model.matrix(object = formula, DataFrameInput)
    class(DM) = 'LT2D.design.matrix'
    return(DM)
}

LinPredictor = function(parameters, DM){
  # Multiply DM by the column matrix of covariate parameters to obtain 
  # the linear predictor for the theta parameter of interest:
  parameters = matrix(parameters, ncol=1)
  
  # Exception handle the matrix multiplication to ensure the correct number 
  # of parameters has been supplied for the relationship specified by the 
  # passed formula object:
  LP = try({DM%*%parameters},silent=T)
  if (class(LP)=='try-error'){
    stop('Incorrect number of parameters supplied to LinPredictor')
  }
  return(LP)
}

h1=function(y,x,b)
{ 
  HazardBCheck(y,x,b,'h1')
  theta1 = as.vector(exp(b[[1]])) ; theta2 = as.vector(exp(b[[2]]))
  return(theta1*(y^2+x^2)^(-theta2/2))
}

fityx = function(y=NULL,x=NULL,b,hr,ystart,pi.x,logphi,w,rmin=0,formulas=NULL,
                 covarPars=NULL,control = list(),hessian = FALSE,corrFlag = 0.7,
                 debug = FALSE, DataFrameInput=NULL,...){
  
  hrname = hr                             # These lines keep the variable names 
  piname = pi.x                           # consistent with previous code
  
  output = list() # We add to this object as the function progresses
  
  HazardWarning = 'Incorrect option supplied for Hazard Rate. Please supply
  the name of desired function as a character.' 
  
  DensityWarning = 'Incorrect option supplied for Density. Please supply
  the name of desired function as a character.'
  
  if ( (is.null(x) | is.null(y)) & is.null(DataFrameInput)){ 
    # If either x or y is missing, and there is no data frame to get both from:
   stop('missing data values')
  }

  if (!length(y)==length(x)){stop('x must have the same length as y')}
  # Truncate the data to perp <= w. It's important to calculate y first,
  # or else nothing will get truncated:
  if (w){

    if (is.null(y)==FALSE & is.null(x)==FALSE){
      y = y[x<=w] ; x = x[x<=w]                       # Truncate the data
    }
    
    else if (!is.null(DataFrameInput)){
      DataFrameInput = subset(DataFrameInput,         # Truncate covar values
          DataFrameInput$x < w)
      
      x = DataFrameInput$x  # get x and y from supplied data.frame
      y = DataFrameInput$y
    }
    else{stop('No data supplied')} # Shouldn't ever reach this line 
    
    warning('data truncated according to user\'s chosen perpendicular truncation
            distance. Data in model object may be a different dimension to the
            supplied data')
  }
  else{stop('w must be supplied to fityx')}

  if (!class(hr)=='character'){           # We check that the passed functions
    stop(HazardWarning)                   # are strings, and hence a name
  }

  if (!class(pi.x)=='character'){          
    stop(DensityWarning)                    
  }
  
  rounded.points.count = length(x[sqrt(x**2 + y**2)<=rmin])# save for likelihood
  new.x = x[sqrt(x**2 + y**2)>rmin]  # We will Only pass unrounded data and 
  new.y = y[sqrt(x**2 + y**2)>rmin]  # the number of rounded points to the nll
  
  if (length(new.x)!=length(new.y)){stop('x and y dimension mismatch')}
  
  DesignMatrices = NULL # If there are no formulas, this variable will not get
  # overwritten, and so DesignMatrices = NULL will get passed to the negative
  # log-likelihood, indicating that no covariates have been included
  
  if (!is.null(formulas)){ # This massive 'if' block deals with producing design 
    # matrices and recreating the vector b for the hazard function with the
    # right structure, so that the neg.log.lik can create the linear predictors
    
    if (is.null(DataFrameInput) | is.null(covarPars)){
      stop('When formulas have been specified,
           DataFrameInput and covarPars are required')}
    
    # only the unrounded data points will be passed to the likelihood, so we 
    # must ensure that the dimension of the data frame we use to create our 
    # design matrices is the same as the x and y vectors
    unrounded.points.data.frame = subset(DataFrameInput,      
        sqrt(DataFrameInput$x**2 + DataFrameInput$y**2)>rmin)
    #
    #
    #
    #
    # When adding covars to RMIN stuff. Remember to remove the right 
    # rows from the design matrices 
    #
    #
    #
    #
    
    slotsAllowed = HazardCovarSlots(hr)
    # Calling the below function ensures the user hasn't defined any formulas
    # where they aren't supposed to, and checks that other types of consistency
    # requirements are also satisfied. Errors are raised if they are not.
    checkedFormulas = FormulaChecking(hr, formulas)
    i.formula = checkedFormulas[[1]]
    x.formula = checkedFormulas[[2]]
    y.formula = checkedFormulas[[3]]
    # Knowing if x and y formulas are the same tells us if the y formula 
    # LP goes in slot 2 or 4
    xy = FALSE ; if ('xy' %in% checkedFormulas[[4]]){xy=TRUE}

    DesignMatrices = list() # Make a list of all design matrices to pass to the 
    # negative log likelihood through optim. We insert the DM at the same index
    # as the parameters for the theta whose formula it describes, for simplicity
    
    b.for.optim = as.list(b) # To preserve the intercepts when no DM for entry
    # covar pars must be a list with names i, x and/or y, where appropriate
    if(!is.null(i.formula)){ # we have an i formula - always slot 1
      DesignMatrices[[1]] = DesignMatrix(unrounded.points.data.frame, i.formula)
      b.for.optim[[1]] = c(b[1],covarPars$i)
    }
    
    if (xy & (!is.null(x.formula))){ # if x and y are the same, one formula
      # describes both of them, and it goes in slot 2:
      DesignMatrices[[2]] = DesignMatrix(unrounded.points.data.frame, x.formula)
      b.for.optim[[2]] = c(b[2], covarPars$x)
    }
   
    else if (!xy){ # x and y formulas can be different
      if (!is.null(x.formula)){ 
        DesignMatrices[[2]]=DesignMatrix(unrounded.points.data.frame, x.formula)
        b.for.optim[[2]]=c(b[2], covarPars$x)# if x is there, it goes in slot 2
      }
      
      if (!is.null(y.formula)){ 
        DesignMatrices[[4]]=DesignMatrix(unrounded.points.data.frame, x.formula)
        b.for.optim[[4]]=c(b[4], covarPars$y)# if y is there, it goes in slot 4
      }
    }
  }
  
  else{b.for.optim = as.list(b)} # No covariates, so b is simply the one we have
  
  # We pack the parameters as a vector:
  
  if (piname == "pi.const"){pars = list(beta = b.for.optim)} # construct a list
  else{pars = list(beta=b.for.optim,logphi=logphi)} # of the right form, then 
  u.pars = unlist(pars)                             # vectorise it for optim
  
  fit = optim(
    par = u.pars,fn = negloglik.yx,y = new.y,x = new.x,hr = hrname,
    ystart = ystart, pi.x = piname,w = w, hessian = hessian,
    debug = debug, rounded.points = rounded.points.count,
    control = control,DesignMatrices=DesignMatrices, skeleton=pars,...
  )
  
  
  # Extremely stupid because it calls the same function twice for no reason,
  # needed a very quick fix - will make this better if there's time:
  # Gets the betas for the fitted parameter values, for each row of the data
  B.per.line <-  negloglik.yx(fit$par,
                              y=new.y, x=new.x, hr=hrname, ystart=ystart,
                              pi.x=piname, w=w,
                              DesignMatrices=DesignMatrices,
                              skeleton=pars, 
                              returnB = T)[[1]]
  
  
  unrounded.points.with.betas <- list(data.frame(x=new.x,y=new.y),b=B.per.line)
  
  error = FALSE
  if (fit$convergence != 0) {
    warning('Convergence issue (code = ',
            fit$convergence,') . Check optim() help.')
    error = TRUE
  }

  # ***
  par.as.list = relist(fit$par, skeleton = pars)
  b = par.as.list[[1]]
  b.as.vector = unlist(b)
  logphi = try(par.as.list[[2]])
  if (class(logphi)=='try-error'){logphi = NA}
  # ***

  if (hessian){
    
    mNames = names(fit$par)
    
    if (any(diag(fit$vcov) <= 0)) {
      warning('Failed to invert hessian.  Model convergance problem in fityx?')
      error = TRUE
      CVpar = rep(NA,length(fit$par))
    }
    
    else {CVpar = sqrt(diag(solve(fit$hessian))) / abs(fit$par)}
    
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
    output$hessian = hessian
    output$CVpar = CVpar
    output$corr = corr 
  }
  
  AICval = 2 * fit$value + 2 * length(fit$par)
  
  dat = data.frame(x = x,y = y)
  #p0 = 1 - Sy(0,0,ystart,b,hrname)   # reinclude once we decide how to deal 
  
  output$par = fit$par ; output$value = fit$value
  output$counts = fit$counts ; output$convergence = fit$convergence
  output$message = fit$message
  output$error = error ; output$hr = hr
  output$pi.x = pi.x ; output$ystart = ystart ; output$w = w ; output$b = b
  output$logphi = logphi ; output$AIC = AICval ; output$dat = dat
  output$unrounded.points.with.betas = unrounded.points.with.betas
  #output$p0 = p0
  output$rmin = rmin
  
  class(output) = 'LT2D.fit.object'
  return(output)
}

negloglik.yx=function(pars,y,x,hr,ystart,pi.x,w,rounded.points=0,
                      DesignMatrices=NULL, skeleton=NULL, debug=FALSE,
                      returnB = FALSE){
  # The only work we want to have to redo is make the linear predictors. 
  # so we allow design matrices to be passed as arguments, and we make the 
  # linear predictors and remake b appropriately
  unpacked = relist(pars, skeleton=skeleton)    # Unpack the parameters.
  b = unpacked[[1]]                             # See fityx to see how they
  logphi = try (unpacked[[2]],silent=T)         # they were packed.
  
  if (class(logphi)=='try-error'){logphi = NA}
  if(length(y)!=length(x)) stop("Lengths of x and y must be the same.")
  
  if(debug) print(pars)
  
  if (!class(hr)=='character'){stop('message from negloglik: hr must 
                                    be passed as a character')}
  if (!class(pi.x)=='character'){stop('message from negloglik: pi.x must 
                                      be passed as a character')}
  
  n=length(y) ; hrname = hr ; piname = pi.x 
  hr=match.fun(hr) ; pi.x=match.fun(pi.x)
  
  # if any b entry has non-scalar entries, replace these entries with the 
  # appropriate linear predictor...
  if (!is.null(DesignMatrices)){  # We have been given design matrices
    # The DM for entry 'i' of b, is also at entry 'i' of DesignMatrices
    for (index in (1:length(b))){
      if (length(b[[index]])>1){  # We need to create the LP:
        LP = LinPredictor(b[[index]], DesignMatrices[[index]]) # Create LP
        b[[index]] = as.vector(LP)# Replace the start parameters with the LP
      }
    }
    
    if (rounded.points!=0){
      S1 = 'Covariate inclusion has not yet been implemented for Rmin > 0.'
      S2 = 'The truncated part of the likelihood has assumed common detection'
      S3 = 'throughout all covariate levels'
      warning(paste(S1,S2,S3))}
  }
  # This loop turns entries in b of length 1 into entries of length n
  j = 1 # At the end of this loop, j-1 will be the number of hr parameters
  while (class(try(b[[j]],silent=T))!='try-error'){
    # if there's only one element in a column, make a column of that 
    # element repeated the correct number of times to match the number of 
    # (x,y) coordinates we have:
    if (length(b[[j]])==1){b[[j]]=as.numeric(rep(b[[j]],n))}
    j = j + 1
  }
  j = j - 1 # Remove 1 from j so that it now represents the dimension of b
  likelihood.level.b = list(b,j)
  class(likelihood.level.b) = 'likelihood.level.b'
  if (returnB) return(likelihood.level.b)

  # calculate numerator:
  num=sum(log(fyx(y,x,likelihood.level.b,hrname,ystart))
                                + log(pi.x(x,logphi,w)))
  
  # calculate denominator:
  if (is.null(DesignMatrices)){                 # Easy with no covariates:
    b.normal <- as.list(sapply(b, '[[', 1))     # Extract the 'normal' b
    
    int=integrate(f=p.pi.x,lower=0,upper=w,b=b.normal,hr=hrname,
                  ystart=ystart,pi.x=piname,logphi=logphi,w=w)
    
    denom=n*log(int$value)
  }
  
  else{ # We must work out the integral individually for each (x,y) point:
    integrals = rep(NA, n)
    for (i in (1:n)){
      # We know j exists, and we know b will have been turned into a list
      # of all vectors, so we use b and j from the (not directly) above 
      # while loop
      integration.b = list()
      for (column.num in (1:j)){
        integration.b[[column.num]] = b[[column.num]][i]
      }
      integrals[i] = integrate(f=p.pi.x,lower=0,upper=w,b=integration.b,
                               hr=hrname, ystart=ystart, pi.x=piname, 
                               logphi=logphi, w=w)$value
    }
    denom=sum(log(integrals))
  }
  
  # likelihood:
  llik=-(num-denom) # 2016 paper likelihood, for the un-rounded data points
  
  if(rounded.points>0){
    negllik.rounded = round.lik(rounded.points,pi.x=piname,
                                logphi,rmin,ymax=ystart,hr=hrname,b,w)
  } # we calculate the rounded data points part of the likelihood
  
  else{negllik.rounded=0} # No addition to normal likelihood needed
  
  
  return(llik + negllik.rounded) #round.lik returns a neg.log.lik, so we add. 
}

LT2D.fit = function(DataFrameInput,hr,b,ystart,pi.x,logphi,w,formulas=NULL,
                    ipars=NULL,xpars=NULL,ypars=NULL,rmin=0,
                    control = list(),hessian=FALSE,corrFlag = 0.7,
                    debug = FALSE){
  
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
  
  covarPars = NULL      # if any formulas were specified by the user, 
  covarPars$x = xpars   # we add them to the covar pars object in the
  covarPars$y = ypars   # correct slot. Otherwise, if none were specified, 
  covarPars$i = ipars   # we set covarPars to NULL before the call to fityx
  
  if (length(covarPars)==0) covarPars = NULL
  
  # run the fityx call:
  # We don't exception handle the call, to allow its own errors to be flagged
  fitted.model = fityx(x = NULL, y=NULL,
                       DataFrameInput = NoNAs,                  # Data
                       b=b,hr=hr,                               # Hazard stuff
                       pi.x=pi.x,logphi=logphi,                 # Density stuff
                       ystart=ystart,w=w,rmin=rmin,             # Settings
                       formulas=formulas,covarPars = covarPars, # Covar stuff
                       control=control, hessian=hessian,        # Optim stuff
                       corrFlag=corrFlag)       
  
  if (OnlyCallFityx==TRUE){return(fitted.model)}
  return(fitted.model)
  # Now we pass the data frame and fitted model objects to the abundance
  # estimation function and return the output
  output = NDest(DataFrameInput, fitted.model)
  print(output)
  output$fit = fitted.model # We attach the fityx model to the output
  class(output) = 'LT2D.fit.function.object'
  return(output)
}

fyx=function(y,x,b,hr,ystart,nint=100)
{
  if(length(y)!=length(x)) stop("Lengths of x and y must be the same.")
  n=length(x)
  f=intval=rep(NA,n)
  if (!class(hr)=='character'){stop('message from fyx: hr must be passed as 
                                    a character')}
  hr=match.fun(hr)

  # The below loop was added to generalise the function enough to deal with 
  # covariate inclusion.
  
  if (class(b)=='likelihood.level.b'){
    j = b[[2]] # keep track of the number of parameters
    b = b[[1]] # the b list with all columns of equal dimension is found here
  }
  
  else {
    # This loop is largely redundant. It is expected this function will
    # only get called by the likelihood function (negloglik.yx), which performs
    # all the necessary checks and produces the 'likelihood.level.b' object.
    # Leaving it in means that if this fyx function is called manually, or by
    # a different function, it is able to handle a list b with a mixture 
    # of dimension 1 and linear predictor parameters. 

    j = 1 # At the end of this loop, j-1 will be the number of hr parameters
    while (class(try(b[[j]],silent=T))!='try-error'){
      
      # if there's only one element in a column, make a column of that 
      # element repeated the correct number of times to match the number of 
      # (x,y) coordinates we have:
      
     if (length(b[[j]])==1){b[[j]]=rep(b[[j]],n)}
     j = j + 1
    }
    j = j - 1 # Remove 1 from j so that it now represents the dimension of b
  }
  
  for(i in (1:n)) {
    # create an empty list to insert b values for this individual (x,y):
    integration.b = list()

    for (column.number in (1:j)){
      # Take the i th element of each beta column to obtain the correct 
      # parameter values for the point (x,y) in position i
      integration.b[[column.number]] = b[[column.number]][i]
    }
    
    dy=(ystart-y[i])/nint/2                           # for crude integration
    yy=seq(y[i],ystart,length=(nint+1))[-(nint+1)]+dy # for crude integration
    
    int=sum(hr(yy,rep(x[i],nint),integration.b)*dy*2) # crude integration
    intval[i]=exp(-int)
  }
  
  hrval=hr(y,x,b) # value of hr at each (x,y)
  bads=which(hrval>=.Machine$double.xmax)  # identify infinite hazards
  if(length(bads)>0) { # infinite hazard so p(detect)=0
    
    # We can't simply do b[-bads], so we are forced to redo a for loop:
    for (col.num in (1:j)){b.bads = b[[col.num]][-bads]}
    
    f = rep(NA,length(x))
    f[bads]=.Machine$double.xmax
    f[-bads]=hr(y[-bads],x[-bads],b.bads)*intval[-bads]
  }else{
    f=hr(y,x,b)*intval
  }
  return(f)
}

p.pi.x=function(x,b,hr,ystart,pi.x,logphi,w){
  
  if (!class(hr)=='character'){stop('hr must be supplied as character')}
  if (!class(pi.x)=='character'){stop('pi.x must be supplied as character')}
  
  pi.x = match.fun(pi.x) # So that we can evaluate it as a function below
  
  return(px(x,b,hr,ystart)*pi.x(x,logphi,w))
}

px=function(x,b,hrname,ystart){
  if (!class(hrname)=='character'){stop('message from px: hr must be supplied as 
                                    a character')}
  return(1-Sy(x,rep(0.0001,length(x)),ystart,b,hrname))
}

Sy=function(x,y,ymax,b,hr) {
  if (!class(hr)=='character'){stop('hr must be passed as a 
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

NDest <- function(dat, hmltm.fit){
  
  W = hmltm.fit$w # hmltm.fit should be an LT2D.model.fit object,
                  # which has attribute w 
  
  # remove the smaller than w observations from the data frame which won't 
  # have been used to fit the likelihood, but we musn't remove the NAs:
  dat = subset(dat, dat$x<W | is.na(dat$x))
  
  # Add 1/p column
  dat$invp <- rep(NA,dim(dat)[1])
  invp <- invp1_replacement(dat, hmltm.fit) 
  
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
              ests=data.frame(stratum=stratname,
                              n=n,
                              k=k,
                              L=L,
                              covered.area=a,
                              stratum.Area=A,
                              Dgroups=signif(Dg,3),
                              Ngroups=signif(Ng,3),
                              mean.size=round(sbar,1),
                              D=signif(D,5),N=round(N,1))
  )
  )
}

invp1_replacement = function(LT2D.df, LT2D.fit){
  # From the fitted model we extract the invp value:
  # p = phat(LT2D.fit)
  # inversep = 1/p
  # invp.vector = rep(inversep, length(LT2D.df$x))
  # # And we add it to the data.frame:
  # LT2D.df.new = LT2D.df              # copy the input fitted data frame
  # LT2D.df.new$invp = invp.vector     # add the invp column
  w = LT2D.fit$w
  ystart = LT2D.fit$ystart
  hr = LT2D.fit$hr
  pi.x = LT2D.fit$pi.x
  logphi = LT2D.fit$logphi
  
  unrounded.points.with.betas = LT2D.fit$unrounded.points.with.betas
  
  LT2D.df$invp <- rep(NA, dim(LT2D.df)[1])
  
  # For each row of the data.frame:
  
  # check if covariates were used here, return invp everywhere if there were not
  # and proceed to the for loop if they were...
  
  for (i in (1:dim(LT2D.df)[1])){
    data.row <- LT2D.df[i,]
    x = data.row$x
    y = data.row$y
    
    if (is.na(data.row$object)){
      # if a transect has no detections, it will not be used to calculate the 
      # Horvitz-Thompson estimate, and so we may just set the invp value to 
      # zero:
      
      LT2D.df$invp[i] <- 0
    }
    
    else{ # if the data point is a detection
      
      # if we had covariates included:
      if(!is.null(LT2D.fit$COVARIATESINCLUDED)){}
    }
  }
  return(LT2D.df)
}