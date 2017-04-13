rm(list=ls())
library(dplyr);library(gtools);library(magrittr);library(plyr);library(reshape2)

source('./lsl_functions.R')
source('./lsl_sem.R')

beta_vf <- matrix(NA, 9, 3)
beta_vf[c(1,2,3), 1] <- beta_vf[c(4,5,6), 2] <- beta_vf[c(7,8,9), 3] <- 1

pattern=list(beta_vf=beta_vf)

dta       <- lavaan::HolzingerSwineford1939

rc_sem <- lslSEM()
rc_sem$input(raw_obs=dta,var_subset = c(7:15))
rc_sem$specify(pattern)
rc_sem$learn(penalty="scad",gamma = seq(0.1),delta=2.5)

rc_sem$knowledge

model<-rc_sem$model
data<-rc_sem$data

pl<-"l1"
gamma  <-seq(0.025,0.1,0.025)
delta  <-c(2.5)
control<-list(max_iter=500,rel_tol=10^(-5))
mat       <-attributes(model)$mat
allpen<-expand.grid(pl=pl,delta=delta,gamma=gamma)
x=1
penalize<-list(pl=pl,delta=allpen[x,2],gamma=allpen[x,3])
maxit=control[[1]]
cri=control[[2]]


