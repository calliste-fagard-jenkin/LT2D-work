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
    
    row <- which(dat$stratum==invp$stratum[i]     # ensures data is consistent
                 & dat$transect==invp$transect[i] # with invp version, and that
                 & dat$object==invp$object[i])    # only detection rows count
    
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

invp1_replacement = function(LT2D.df, LT2D.fit.obj){
  # purpose : Adds a column of inverse detection probability to a data frame 
  #           used as the data.frame argument to a call of LT2D.fit
  # inputs  : LT2D.df      - The data.frame containing all the transect
  #                          detections (and non-detections), as it was passed
  #                          to LT2D.fit
  #           LT2D.fit.obj - The fitted model object produced by fityx
  # output  : A data.drame, LT2D.df, with an extra column, invp, for the inverse
  #           of estimated detection probability 
  w = LT2D.fit.obj$w            # we extract all the values needed
  ystart = LT2D.fit.obj$ystart  # to calculate inverse p from 
  hr = LT2D.fit.obj$hr          # the fit object
  pi.x = LT2D.fit.obj$pi.x
  logphi = LT2D.fit.obj$logphi
  
  unrounded.points.with.betas <- LT2D.fit.obj$unrounded.points.with.betas
  converted.betas <- data.with.b.conversion(fityx.output.object = LT2D.fit.obj)
  betas <- converted.betas[[3]]
  x <- converted.betas[[1]]
  y <- converted.betas[[2]]
  
  LT2D.df$invp <- rep(NA, dim(LT2D.df)[1])
  
  # check if covariates were used, return invp everywhere if there were not
  # and proceed to the for loop if they were...
  if (LT2D.fit.obj$covariates!= TRUE){
    
    # since no covariates were included, the value of the betas should be the
    # the same everywhere:
    B.no.covar = as.list(betas[[1]])
    LT2D.df$invp = rep(1/phat(w = w, hr = hr, b = B.no.covar, ystart = ystart, 
                              pi.x = pi.x, logphi = logphi))
    return(LT2D.df)
  }
  
  # We want to ensure that the data for which we have detection function values
  # is of the same dimension as the number of rows of data which represent a 
  # detection:
  if (dim(subset(LT2D.df, !is.na(LT2D.df$object)))[1] != length(x)){
    stop('Input data frame and fitted data dimensions do not match')
  }
  
  NAcounter <- 0  # Keep track of how many non-detection rows we have considered
  for (i in (1:dim(LT2D.df)[1])){
    data.row <- LT2D.df[i,]
    
    if (!is.na(data.row$object)){
      # if a transect has a detection, we must calculate its value of invp
      # and add it to the appropriate row of the invp column. 
      j <- i - NAcounter
      
      # Check that the x and y values in the data.frame match up the ones which
      # give us the detection function values
      if(data.row$x != x[j] | data.row$y!= y[j]) stop('mismatched data entry')
      
      B <- as.list(betas[[j]])
      LT2D.df$invp[i] <- 1/phat(w = w, hr = hr, b = B, ystart = ystart, 
                                pi.x = pi.x, logphi = logphi)
    }
    
    else{NAcounter <- NAcounter + 1}
  }
  return(LT2D.df)
}