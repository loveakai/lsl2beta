rm(list=ls())
set.seed(4869)
library(dplyr);library(gtools);library(magrittr);library(plyr);library(reshape2)

source('./lsl_tool.R')
source('./lsl_functions.R')

model.cfa<-'
F1=~0.8*x1+0.8*x2+0.8*x3
F2=~0.8*x4+0.8*x5+0.8*x6
F3=~0.8*x7+0.8*x8+0.8*x9
x1~~(1-0.8^2)*x1
x2~~(1-0.8^2)*x2
x3~~(1-0.8^2)*x3
x4~~(1-0.8^2)*x4
x5~~(1-0.8^2)*x5
x6~~(1-0.8^2)*x6
x7~~(1-0.8^2)*x7
x8~~(1-0.8^2)*x8
x9~~(1-0.8^2)*x9
F2~~1*F2
F3~~1*F3
F1~~1*F1
F1~~0*F2
F1~~0*F3
F2~~0*F3
'

model.cfa2<-'
F1=~0.8*x1+0.6*x2+0.6*x3
F2=~0.8*x4+0.6*x5+0.6*x6
F3=~0.8*x7+0.6*x8+0.6*x9
x1~~(1-0.8^2)*x1
x2~~(1-0.6^2)*x2
x3~~(1-0.6^2)*x3
x4~~(1-0.8^2)*x4
x5~~(1-0.6^2)*x5
x6~~(1-0.6^2)*x6
x7~~(1-0.8^2)*x7
x8~~(1-0.6^2)*x8
x9~~(1-0.6^2)*x9
F2~~1*F2
F3~~1*F3
F1~~1*F1
F1~~0*F2
F1~~0*F3
F2~~0*F3
'



model.cfa3<-'
F1=~0.8*x1+0.4*x2+0.4*x3
F2=~0.8*x4+0.4*x5+0.4*x6
F3=~0.8*x7+0.4*x8+0.4*x9
x1~~(1-0.8^2)*x1
x2~~(1-0.4^2)*x2
x3~~(1-0.4^2)*x3
x4~~(1-0.8^2)*x4
x5~~(1-0.4^2)*x5
x6~~(1-0.4^2)*x6
x7~~(1-0.8^2)*x7
x8~~(1-0.4^2)*x8
x9~~(1-0.4^2)*x9
F2~~1*F2
F3~~1*F3
F1~~1*F1
F1~~0*F2
F1~~0*F3
F2~~0*F3
'


 dta       <-list()
 dta[[1]]  <- lavaan::simulateData(model.cfa,sample.nobs = 10000L)  %>% cbind(group="g1")#%>% cbind(.,sample(c(1,2),size=nrow(.),rep=T))
 dta[[2]]  <- lavaan::simulateData(model.cfa2,sample.nobs = 10000L) %>% cbind(group="g2")
 dta[[3]]  <- lavaan::simulateData(model.cfa3,sample.nobs = 10000L) %>% cbind(group="g3")
 dta       <- do.call(rbind,dta)



#mat       <- matgen(lambda=beta_vf)


#dta       <- lavaan::HolzingerSwineford1939[7:15]
data      <- import(dta,var_group = 10)

beta_vf <- matrix(NA, 9, 3)
beta_vf[c(1,2,3), 1] <- beta_vf[c(4,5,6), 2] <- beta_vf[c(7,8,9), 3] <- 1
pattern <-list()
pattern$beta_vf<-beta_vf
beta_vf <- matrix(0,9,3)
beta_vf[c(2,3),1]    <- beta_vf[c(5,6),2]    <- beta_vf[c(8,9),3]    <- 1
beta_vf[1,1]         <- beta_vf[4,2]         <- beta_vf[7,3]         <- 0.8
value <-list()
value$beta_vf<-beta_vf
model     <- specify(pattern,value)
#model <- specify(pattern)

learn     <- function(penalty,gamma,delta,control=list(max_iter,rel_tol)){
  if (missing(penalty)) {pl<-"l1"} else {pl<-penalty}
  if (missing(gamma))   gamma  <-seq(0.025,0.1,0.025)
  if (missing(delta))   delta  <-2.5
  if (missing(control)) {control<-list(max_iter=500,rel_tol=10^(-5))} else {
    if (is.null(control$max_iter)) {control$max_iter<-500}
    if (is.null(control$rel_tol))  {control$rel_tol <-10^(-5)}
    }
  mat       <-attributes(model)$mat
  eta_label <-c(attributes(model)$labels[[1]],attributes(model)$labels[[2]])
  n_groups  <-attributes(data)$n_groups
  n_eta     <-attributes(mat)$n_eta
  n_v       <-attributes(mat)$n_v
  n_f       <-attributes(mat)$n_f
  raw_obs   <-data$raw_obs
  
  ide       <- diag(1, ncol = n_eta, nrow = n_eta)  %>% `colnames<-`(eta_label) %>% `rownames<-`(eta_label) 
  G_eta     <- c(rep(T,n_v),rep(F,n_f)) %>%`names<-`(eta_label)
  sigma     <- lapply(1:n_groups, function(i_groups) { data$raw_cov[[i_groups]] + data$raw_mean[[i_groups]] %>% tcrossprod })
  e_v       <- lapply(1:n_groups, function(i_groups) { data$raw_mean[[i_groups]] })
  
  

  
  allpen<-expand.grid(pl=pl,delta=delta,gamma=gamma)
  
  for (p in (1:nrow(allpen))){
  penalize  <- list(pl=allpen[p,1],delta=allpen[p,2],gamma=allpen[p,3])
  
  #ECM
  
  ecmm<-ecm(mat=mat,ide=ide,G_eta=G_eta,maxit=control[[1]],cri=control[[2]],penalize=penalize)
  
}

  ecmm$theta$mat$value$beta_r+ecmm$theta$mat$value$beta_i[[1]]
  ecmm$theta$mat$value$beta_r+ecmm$theta$mat$value$beta_i[[2]]
  ecmm$theta$mat$value$beta_r+ecmm$theta$mat$value$beta_i[[3]]
  