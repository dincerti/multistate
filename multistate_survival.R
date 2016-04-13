rm(list = ls())
setwd("/Users/devinincerti/Dropbox/Projects/Multi-State Survival")

## @knitr setup
source("sim.R")
source("msplots.R")
source("misc.R")
packages <- c("mstate", "flexsurv", "mvtnorm",  
              "ggplot2", "data.table", "tables")
lapply(packages, library, character.only = TRUE)
theme_set(theme_bw())
set.seed(101)

## ---- MULTI-STATE DATA -------------------------------------------------------
## @knitr create_data
library("mstate")
data(prothr)
dat <- prothr
rm(prothr)
dat <- rename(dat, ftreat = treat)
dat$years <- (dat$Tstop - dat$Tstart)/365.25
dat <- dat[dat$years != 0 & dat$Tstop < 3500 & dat$years < 9, ]
dat$trans <- as.factor(dat$trans)
dat$Tstart <- dat$Tstart/365.25
dat$Tstop <- dat$Tstop/365.25

## @knitr show_data
tmp <- data.frame(dat[, c("id", "from", "to", "trans",
                          "Tstart", "Tstop", "years", 
                          "status", "ftreat")])
tmp$Tstart <- round(tmp$Tstart, 2)
tmp$Tstop <- round(tmp$Tstop, 2)
tmp$years <- round(tmp$years, 2)
print(tmp[tmp$id %in% c(33), ], row.names = FALSE)

## @knitr trans
tmat <- attributes(dat)$trans
events(dat)$Freq

## ---- FITTING MULTI-STATE MODELS ---------------------------------------------
## @knitr fit_cox
dat$treat <- ifelse(dat$ftreat == "Placebo", 0, 1)
dat <- expand.covs(dat, "treat")
crcox <- coxph(Surv(years, status) ~ strata(trans) + treat.1 + treat.2 + 
                 treat.3 + treat.4, data  = dat)
cfcox <- coxph(Surv(Tstart, Tstop, status) ~ strata(trans) + treat.1 + treat.2 + 
                 treat.3 + treat.4, data  = dat)

## @knitr fit_weibull
crjwei <- flexsurvreg(Surv(years, status) ~  trans + shape(trans) +
                         treat.1 + treat.2 + treat.3 + treat.4,
                       data = dat, dist = "weibull")
crwei <- vector(4, mode = "list")
for (i in 1:4){
  crwei[[i]] <- flexsurvreg(Surv(years, status) ~ treat,
                            subset = (trans==i),  data = dat,
                            dist = "weibull")
}

## @knitr fit_gompertz
crgomp <- vector(4, mode = "list")
for (i in 1:4){
  crgomp[[i]] <- flexsurvreg(Surv(years, status) ~ treat,
                            subset = (trans==i),  data = dat,
                            dist = "gompertz")
}

## ---- MODEL CHECKING ---------------------------------------------------------
## @knitr cumhaz_cox
mspat.t <- data.frame(strata = 1:4, diag(1, 4))
mspat.c <- data.frame(strata = 1:4, matrix(0, 4, 4))
colnames(mspat.t) <- colnames(mspat.c) <- c("strata", paste0("treat.", seq(1, 4)))
crcox.msfit.t <- msfit(crcox, newdata = mspat.t, variance = FALSE, 
                     trans = tmat)
crcox.msfit.c <- msfit(crcox, newdata = mspat.c, variance = FALSE, 
                       trans = tmat)
cfcox.msfit.t <- msfit(cfcox, newdata = mspat.t, variance = FALSE, 
                     trans = tmat)
cfcox.msfit.c <- msfit(cfcox, newdata = mspat.c, variance = FALSE, 
                       trans = tmat)

## @knitr cumhaz_parametric
tgrid <- seq(0, 30, by = 0.1)
flexpat.t <- data.frame(treat = 1)
crwei.msfit.t <- msfit.flexsurvreg(crwei, t = tgrid, newdat = flexpat.t, 
                                 trans = tmat, variance = FALSE)
crgomp.msfit.t <- msfit.flexsurvreg(crgomp, t = tgrid, newdat = flexpat.t, 
                                 trans = tmat, variance = FALSE)

## @knitr cumhaz_plot
cumhazPlot(cumhaz = list(crcox.msfit.t$Haz, crwei.msfit.t$Haz, crgomp.msfit.t$Haz),
           model.names = c("Cox", "Weibull", "Gompertz"),
           trans = tmat)

## ---- TRANSITION PROBABILITIES -----------------------------------------------
## @knitr cf_cox_tp
cfcox.tp.t <- probtrans(cfcox.msfit.t, predt = 0, direction = "forward",
                 variance = FALSE)[[1]]
cfcox.tp.c <- probtrans(cfcox.msfit.c, predt = 0, direction = "forward",
                        variance = FALSE)[[1]]
tpPlot(probs = list(cfcox.tp.t, cfcox.tp.c),
       treat.names = c("Prednisone", "Placebo"),
       state.names = c("Normal", "Low", "Death"))

## @knitr cr_cox_tp
tv <- unique(crcox.msfit.t$Haz$time)
crcox.sim.t <- mssample(crcox.msfit.t$Haz, trans = tmat, clock = "reset",
                        M = 1000, tvec = tv, do.trace = NULL)
crcox.sim.c <- mssample(crcox.msfit.c$Haz, trans = tmat, clock = "reset",
                      M = 1000, tvec = tv, do.trace = NULL)
tpPlot(probs = list(crcox.sim.t, crcox.sim.c),
       treat.names = c("Prednisone", "Placebo"),
       state.names = c("Normal", "Low", "Death"))

## ---- LIFETIME SIMULATION ----------------------------------------------------
## @knitr cr_parametric_sim
# specify distribution info for each transition 
cr.dist <- vector(4, mode = "list")
for (i in 1:4){
  cr.dist[[i]] <- SimPrep(crgomp[[i]])
}

# individuals to simulate
sim.dat.t <- data.frame(int = 1, treat = 1)
sim.dat.t <- sim.dat.t[rep(seq_len(nrow(sim.dat.t)), 
                             each = 100000), ]
sim.dat.c <- data.frame(int = 1, treat = 0)
sim.dat.c <- sim.dat.c[rep(seq_len(nrow(sim.dat.c)), 
                           each = 100000), ]

# simulate
sim.t <- simMS(x = cr.dist, trans = tmat, t = 30, newdata = sim.dat.t)
sim.c <- simMS(x = cr.dist, trans = tmat, t = 30, newdata = sim.dat.c)
los.t <- simLOS(sim.t, tmat)
los.c <- simLOS(sim.c, tmat)

# make pretty table
los <- rbind(los.t, los.c)
colnames(los) <- c("state", "LOS", "dQALYs")
los$treat <- rep(c("Prednisone", "Placebo"), each = 2)
los$treat <- factor(los$treat)
los.tab <- tabular(Heading() * state ~ Format(digits=2)*Heading()*treat*(LOS + dQALYs) * 
                     Heading() * (identity), data = los)
latex(los.tab)

## @knitr cr_parametric_sim_flexsurv
# compare to sim.fmsm in flexsurv. note that it is only the same because
# we are using the same distribution for each transition and individuals
# in newdata are all the same!
sim.flexsurv.t <- sim.fmsm(x = crgomp, trans = tmat, newdata = data.frame(treat = 1),
                         t = 30, M = 100000)
mean(sim.flexsurv.t$t[, ncol(sim.flexsurv$t)])
mean(sim.t$t[, ncol(sim.t$t)])

## @knitr psa
sim.psa <- simPSA(simdist = cr.dist, x = crgomp, B = 1000, trans = tmat, t = 30, 
       newdata = sim.dat.t[1:1000, ])
dqalys <- c(mean = mean(sim.psa), quantile(sim.psa, probs = c(.025, .5, .975)))
dqalys <- t(data.frame(dqalys))
rownames(dqalys) <- "dQALYs using Prednisone"
print(round(dqalys, 2))
