
rm(list = ls())
library(dplyr)
library(gtools)
library(magrittr)
library(plyr)
library(reshape2)

source('./lsl_functions.R')
source('./lsl_sem.R')

beta_vf <- matrix(NA, 9, 3)
beta_vf[c(1, 2, 3), 1] <-
 beta_vf[c(4, 5, 6), 2] <- beta_vf[c(7, 8, 9), 3] <- 1

pattern = list(beta_vf = beta_vf)

dta <- lavaan::HolzingerSwineford1939

rc_sem <- lslSEM()
rc_sem$input(raw_obs = dta, var_subset = c(7:15))
rc_sem$specify(pattern)
rc_sem$learn(penalty = "scad",
 gamma = seq(0.01,0.2,0.01),
 delta = 2.5)

rc_sem$knowledge


