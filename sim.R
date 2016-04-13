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
    int.indx <- which(d$pars == d$location)
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
  # a particular indiviudal. Used internally in simMS.
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
                  r=.03, qol=NULL, output = "matrix.list"
#                  rtreat=NULL, ptreat=NULL,
                  ){
  # Simulate a fitted parametric semi-Markov multi-state model
  #
  # Args:
  #   x: List containing information about the parametric distribution and the
  #      estimated parameters.
  #   trans: Matrix indicating allowed transitions.
  #   t: time to end simulation
  #   newdata: A n by k dataframe with n individuals and k covariates specifying 
  #            the values of the covariates in the fitted model. 
  #   start: Starting state of the simulation.
  #   tcovs: Names of time-dependent covariates whose values change at the
  #          same rate as time. During each transition, the current time will be
  #          added to the value given in newdata. Age is a typical example.
  #   r: Continuous time discount rate. default is 3%.
  #   qol: Vector containing quality of life weights for each state. Default is 
  #        weight of 1 for every state.
  #   output: how should simulation present output? default is matrix.list which is
  #           a list of 3 matrices st, t, and dqaly with rows for each 
  #           individual. Columns of t contain times when individual changes states
  #           which are shown in st. Columns of dqaly contain cumulative discounted QALYs 
  #           at each transition. Option data.table provides same output in an
  #           enhanced data.frame from the data.tabe package. Simlation is slightly slower 
  #           using this option but output is perhaps easier to manipulate.
  #          
  # Returns:
  #   See output argument above
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
  if (output == "matrix.list"){
      list(st=unname(res.st), t=unname(res.t), dqaly=unname(res.dqaly))
  } else if (output == "data.table"){
      nevents <- ncol(res.t)
      res <- data.table(st = c(t(res.st)),
                            t = c(t(res.t)),
                            dqaly = c(t(res.dqaly)),
                            event = rep(1:nevents, N),
                            id = rep(seq(N), each = nevents))
  }
} 

# SIMULATION RESULTS -----------------------------------------------------------
simLOS <- function(sim, trans){
  # Calculated length of stay and discounted QALYs in each non-absorbing state
  # from a semi-Markov multi-state model.
  #
  # Args:
  #   sim: Simulation output from simMS. May be either a list or a data.table.
  #            estimated parameters for a given transition.
  #   trans: Matrix indicating allowed transitions.
  #
  # Returns:
  #   Dataframe with columns for state, length of stay, and discounted QALYs
  if (class(sim)[1] == "data.table"){
      x = copy(sim)
      N <- length(unique(x$id))
      x[, diff_t := c(diff(t), 0), by = id]
      x[, diff_dqaly := c(diff(dqaly), 0), by = id]
      los <- x[, .(los = sum(diff_t), dqalys = sum(diff_dqaly)),
                     by = st][order(st)]
      los$los <- los$los/N
      los$dqalys <- los$dqalys/N
      los <- los[-absorb, ]
      los$st <- factor(rownames(trans)[-absorb])
      setnames(los, c("state", "los", "dqalys"))
  } else if (class(sim) == "list"){
      n.states <- nrow(trans)
      n <- nrow(sim$t)
      absorb <- absorbing(trans)
      d.t <- diff(t(cbind(sim$t, 0)))
      d.dqaly <- diff(t(cbind(sim$dqaly, 0)))
      st <- factor(t(sim$st), levels = 1:n.states)
      los <- matrix(NA, nrow = n.states, ncol = 2)
      colnames(los) <- c("los", "dqalys")
      los[, 1] <- c(tapply(d.t, st, sum) / n)
      los[, 2] <- c(tapply(d.dqaly, st, sum) / n)
      los <- los[-absorb, ]
      los <- data.frame(state = rownames(trans)[-absorb], los)
  }
  return(data.frame(los))
}

# RANDOM SAMPLE OF PARAMETER DISTRIBUTIONS -------------------------------------
rmvnormPars <- function(x, B){
  # Randomly sample parameters from flexsurv model from their multivariate 
  # asymptotic distribution. This is a Quasi-Bayesian approach because 
  # (for large samples), the posterior distribution of a parameter is close
  # to the distribution of the MLE around the truth. Similar to 
  # normboot.flexsurvreg but does not return predicted values. Used internally
  # in simPSA
  #
  # Args:
  #   x: A fitted model from flexsurvreg
  #   B: Number of samples.
  #
  # Returns:
  #   List of matrixes with one matrix of sampled parameters not a function of 
  #   covariates (parvals) and one matrix of covariate values for parameter 
  #   that is a function of covariates (beta)
  sim <- matrix(nrow = B, ncol = nrow(x$res))
  colnames(sim) <- rownames(x$res)
  sim[, x$optpars] <- rmvnorm(B, x$opt$par, x$cov)
  sim[, x$fixedpars] <- rep(x$res.t[x$fixedpars, "est"], each = B)
  d <- x$dlist
  par.indx <- which(x$dlist$pars != d$location)
  sim <- list(parvals = sim[, par.indx, drop = FALSE],
              beta = sim[, -par.indx, drop = FALSE])
  return(sim)
}

# PROBABILISTIC SENSITIVITY ANALYSIS -------------------------------------------
simPSA <- function(simdist, x, B, trans, t, newdata){
  # Conduct probabilistic sensitivity analysis for semi-Markov multi-state 
  # simulation
  #
  # Args:
  #   simdist: Information about the parametric distribution and the
  #            estimated parameters for a given transition. 
  #   x: List of models (one for each transition) fit with flexsurvreg
  #   B: Number of random samples
  #   trans: Matrix indicating allowed transitions
  #   t: time to end simulation
  #   newdata: A n by k dataframe with n individuals and k covariates specifying 
  #            the values of the covariates in the fitted model. 
  #
  # Returns:
  #   Vector of discounted QALYs for each random sample of parameters
  ntrans <- length(simdist)
  rpars <- vector(4, mode = "list")
  for (i in 1:ntrans){
    rpars[[i]] <- rmvnormPars(x[[i]], B)
  }
  dqaly <- rep(NA, B)
  for (i in 1:B){
    for (j in 1:ntrans){
      simdist[[j]]$parvals <- rpars[[j]]$parvals[i, ]
      simdist[[j]]$beta <- rpars[[j]]$beta[i, ]
    }
    sim <- simMS(x = simdist, trans = trans, t = 30, newdata = newdata)
    dqaly[i] <- mean(sim$dqaly[, ncol(sim$dqaly)])
  }
  return(dqaly)
}

# LIST OF DISTRIBUTIONS --------------------------------------------------------
sim.dists <- list(
  # This is a list of the distributions used in the flexsurv package. It does 
  # not currently include the rng functions.
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

