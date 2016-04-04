# IDENTIFY ABSORBING STATES IN MULTI-STATE MODEL -------------------------------
absorbing <- function(trans){
  # Determines which states in a transition matrix are absorbing
  #
  # Args:
  #   trans: Square transition matrix whose r, s entry is i if the ith 
  #          transition type is r, s
  # Returns:
  #   Integer i of abosrbing transition types
  which(apply(trans, 1, function(x) all(is.na(x))))
}

transient <- function(trans){
  # Determines which states in a transition matrix are not absorbing
  #
  # Args:
  #   trans: Square transition matrix whose r, s entry is i if the ith 
  #          transition type is r, s
  # Returns:
  #   Integer i of non-absorbing transition types
  which(apply(trans, 1, function(x) any(!is.na(x))))
}

# EXTRACT DISTRIBUTION INFO AND PARAMETER ESTIMATES FROM FITTED MODEL ----------
SimPrep <- function(x){
  # Prepares a model fit with flexsurv object for simulation
  #
  # Args:
  #   x: A model fit with flexsurvreg
  #   dlist: list
  #
  # Returns:
  #   List containing information about the parametric distribution and the
  #   estimated parameters.
  if(class(x) == "flexsurvreg"){
    d <- x$dlist
    d$r <- x$dfns$r
    int.indx <- which(d$ßpars == d$location)
    d$parvals <- x$res.t[x$basepars, "est"]
    d$parvals[int.indx] <- NA
    d$beta <- c(x$res.t[int.indx, "est"], x$res.t[x$covpars,"est"])
  } else{
    print(paste0("SimPrep does not work for an object of class ", class(z)))
  }
  return(d)
}

# RETURN SIMULATION PARAMETERS FROM MODEL(S) -----------------------------------
SimPars <- function(simdist, newdata=NULL) {
  # Returns model parameters needed to simulate a transition for
  # a particular indiviudal
  #
  # Args:
  #   simdist: Information about the parametric distribution and the
  #            estimated parameters for a given transition.
  #   newdata: A n by k dataframe with n individuals and k covariates specifying 
  #            the values of the covariates in the fitted model. 
  #
  # Returns:
  #   Model parameters corresponding to the parametric distribution of choice
  X <- if (length(newdata) == 0) matrix(0) else as.matrix(newdata)
  rpars <- vector(length(simdist$pars), mode = "list")
  for (i in seq_along(simdist$pars)){
    if (simdist$pars[i] == simdist$location){
      rpars[[i]] <- X %*% simdist$beta
    } else{
      rpars[[i]] <- simdist$parvals[i]
    }
    rpars[[i]] <- simdist$inv.transform[[i]](rpars[[i]])
  }
  names(rpars) <- simdist$pars
  return(rpars)
}

# SIMULATE MULTI-STATE MODEL ---------------------------------------------------
simMS <- function(x, trans, t, newdata=NULL, start=1, tcovs=NULL, 
                  rtreat=NULL, ptreat=NULL,
                  r=.03, qol=NULL){
  # Simulate a fitted multi-state model
  #
  # Args:
  #   x: List containing information about the parametric distribution and the
  #      estimated parameters.
  #   trans: Matrix indicating allowed transitions.
  #   newdata: A n by k dataframe with n individuals and k covariates specifying 
  #            the values of the covariates in the fitted model. 
  #   start: Starting state of the simulation.
  #   tcovs: Names of time-dependent covariates whose values change at the
  #          same rate as time. During each transition, the current time will be
  #          added to the value given in newdata. Age is a typical example.
  #   r: Continuous time discount rate. default is 3%.
  #   qol: Vector containing quality of life weights for each state. Default is 
  #        weight of 1 for every state.
  #          
  # Returns:
  #   Dataframe with columns for model state, simulation time, and discounted qalys
  #   at each transition.
  if(length(qol)==0){
    qol <- rep(1, nrow(trans))
  }
  N <- nrow(newdata)
  if (length(t)==1) {
    t <- rep(t, N) 
  } else if (length(t)!= N) 
    stop("length of t should be 1 or number of rows in newdata")
  if (length(start)==1) {
    start <- rep(start, N)
  } else if (length(start)!=N) {
    stop("length of start should be 1 or number of rows in newdata")
  }
  
  ### INITIALIZE LOOP
  nst <- nrow(trans)
  res.st <- cur.st <- start
  res.t <- cur.t <- rep(0, N)
  res.dqaly <- cur.dqaly <- rep(0, N)
  todo <- seq_len(N)
  
  ### MAIN LOOP: CONTINUE UNTIL ALL IN ABOSRIBING STATE OR SURVIVED PAST TIME t
  while (any(todo)){ # note that todo is individual's who have not died or lived past time t
    cur.st.out <- cur.st[todo]
    cur.t.out <- cur.t[todo]
    cur.dqaly.out <- cur.dqaly[todo]
    # update newdata to account for aging
    if (length(tcovs)>0){
      newdata[todo, tcovs] <- newdata[todo, tcovs] + cur.t.out
    } 
    done <- numeric()
    
    ### SECOND LOOP: SIMULATE NEXT TIME AND STATES FOR PEOPLE WHOSE CURRENT STATE IS i
    for (i in unique(cur.st[todo])){            
      if (i %in% transient(trans)) { # don't do simulation for people whose current state is aborsbing
        transi <- na.omit(trans[i,])
        ni <- sum(cur.st[todo]==i)
        t.trans1 <- matrix(0, nrow=ni, ncol=length(transi))
        
        ## treatment effect during state 1
        # if (length(rtreat) > 0 & i == 1){
        #   todo1 <- todo[which(cur.st[todo]==1)]
        #   newdata[todo1, rtreat] <- rbinom(ni, 1, ptreat)
        # }
        
        ### THIRD LOOP: SIMULATE TIMES FROM STATE i TO ALL POTENTIAL DESTINATION STATES j
        for (j in seq_along(transi)) { 
          x.trans <- x[[transi[[j]]]]
          simpars <- SimPars(x.trans, newdata[todo, ])
          t.trans1[, j] <- do.call(x.trans$r, c(list(n=ni), simpars))
          basepars <- BaseparsFlexsurv(model[[1]], newdata[todo, ])
          #t.trans1[,j] <- do.call(modbase$dfns$r, c(list(n=ni, basepars[, 1], basepars[, 2])))
          if (is.null(x.trans$r)) stop("No random sampling function found for this model")
        } 
        ### END THIRD LOOP
        ## simulated next state is the one with minimum simulated time
        mc <- max.col(-t.trans1)
        next.state <- match(transi[mc], trans[i,])
        next.time <- t.trans1[cbind(seq_along(next.state), mc)]
        inds <- which(cur.st[todo]==i)
        next.t <- cur.t.out[inds] + next.time
        ## if final simulated state is greater than target time, censor at target time
        cens <- next.t > t[inds]
        next.t[cens] <- t[inds][cens]    
        cur.dqaly.out[inds] <- qol[i] * ((exp(-r*cur.t.out[inds]) -
                                            exp(-r * next.t))/r) + cur.dqaly.out[inds]
        cur.t.out <- next.t
        cur.st.out[!cens] <- next.state[!cens]
        done <- todo[inds][cens]
      }
    } 
    ### END SECOND LOOP
    cur.st[todo] <- cur.st.out
    cur.t[todo] <- cur.t.out
    cur.dqaly[todo] <- cur.dqaly.out
    res.st <- cbind(res.st, cur.st)
    res.t <- cbind(res.t, cur.t)
    res.dqaly <- cbind(res.dqaly, cur.dqaly)
    done <- union(done, which(cur.st %in% absorbing(trans)))
    todo <- setdiff(todo, done)
  }
  ### END MAIN LOOP
  list(st=unname(res.st), t=unname(res.t), dqaly=unname(res.dqaly))
} 

# LIST OF DISTRIBUTIONS --------------------------------------------------------
sim.dists <- list(
  genf = list(
    name="genf",
    pars=c("mu","sigma","Q","P"),
    location="mu",
    transforms=c(identity, log, identity, log),
    inv.transforms=c(identity, exp, identity, exp)
  ),
  genf.orig = list(
    name="genf.orig",
    pars=c("mu","sigma","s1","s2"),
    location="mu",
    transforms=c(identity, log, log, log),
    inv.transforms=c(identity, exp, exp, exp)
  ),
  gengamma = list(
    name="gengamma",
    pars=c("mu","sigma","Q"),
    location="mu",
    transforms=c(identity, log, identity),
    inv.transforms=c(identity, exp, identity)
  ),
  gengamma.orig = list(
    name="gengamma.orig",
    pars=c("shape","scale","k"),
    location="scale",
    transforms=c(log, log, log),
    inv.transforms=c(exp, exp, exp)
  ),
  exp = list(
    name="exp",
    pars=c("rate"),
    location="rate",
    transforms=c(log),
    inv.transforms=c(exp)
  ),
  weibull = list(
    name = "weibull.quiet",
    pars=c("shape","scale"),
    location="scale",
    transforms=c(log, log),
    inv.transforms=c(exp, exp)
  ),
  weibullPH = list(
    name="weibullPH", 
    pars=c("shape","scale"),
    location="scale",
    transforms=c(log, log),
    inv.transforms=c(exp, exp)
  ),
  lnorm = list(
    name="lnorm",
    pars=c("meanlog","sdlog"),
    location="meanlog",
    transforms=c(identity, log),
    inv.transforms=c(identity, exp)
  ),
  gamma = list(
    name="gamma",
    pars=c("shape","rate"),
    location="rate",
    transforms=c(log, log),
    inv.transforms=c(exp, exp)
  ),
  gompertz = list(
    name="gompertz",
    pars=c("shape","rate"),
    location="rate",
    transforms=c(identity, log),
    inv.transforms=c(identity, exp)
  ),
  llogis = list(
    name="llogis",
    pars=c("shape","scale"),
    location="scale",
    transforms=c(log, log),
    inv.transforms=c(exp, exp)
  )
)