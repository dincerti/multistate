rm(list = ls())
setwd("/Users/devinincerti/Dropbox/Projects/Multi-State Survival")

## ---- SLIDE 1 ----------------------------------------------------------------
## @knitr simdata

## @knitr ms_data
library("flexsurv")
tmat <- rbind(c(NA, 1, 2), c(NA, NA, 3), c(NA, NA, NA))
dat <- bosms3
age <- data.frame(id = unique(dat$id), age = rnorm(204, 50, 5)) 
dat <- merge(dat, age, by = "id")
tmat <- rbind(c(NA, 1, 2), c(NA, NA, 3), c(NA, NA, NA))

## @knitr flexsurvreg
cr.jreg <- flexsurvreg(Surv(years, status) ~ trans + shape(trans),
                     data = bosms3, dist = "weibull")
cr.reg <- vector(3, mode = "list")
for (i in 1:3){
  cr.reg[[i]] <- flexsurvreg(Surv(years, status) ~ 1 + age,
                                 subset = (trans==i),  data = dat,
                                 dist = "weibull")
}

## @knitr flexsurv simprep
cr.dist <- vector(3, mode = "list")
for (i in 1:3){
  cr.dist[[i]] <- SimPrep(cr.reg[[i]])
}

## @knitr flexsurv sim
source("sim.R")
sim.dat <- data.frame(int = 1, age = 51.2)
sim.dat <- sim.dat[rep(seq_len(nrow(sim.dat)), 
                             each = 100000), ]

sim <- simMS(x = cr.dist, trans = tmat, t = 30, newdata = sim.dat)
sim2 <- sim.fmsm(x = cr.reg, trans = tmat, newdata = data.frame(age = 51.2),
                 t = 30, M = 100000)
mean(sim$t[, ncol(sim$t)])
mean(sim2$t[, ncol(sim2$t)])