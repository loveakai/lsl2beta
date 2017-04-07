rm(list=ls())
set.seed(4869)
library(dplyr);library(gtools);library(magrittr);library(plyr);library(reshape2)

source('./lsl_functions.R')
source('./lsl_sem.R')
pic<-readLines("yeah.txt",warn=F)

beta_vf <- matrix(0, 9, 3)
#beta_vf <- matrix(0,9,3)
beta_vf[c(1,2,3), 1] <- beta_vf[c(4,5,6), 2] <- beta_vf[c(7,8,9), 3] <- 1
#alpha_v<-rep(0,9)
#alpha_f<-rep(1,3)
#pattern=list(beta_vf=beta_vf,alpha_v=alpha_v,alpha_f=alpha_f)
pattern=list(beta_vf=beta_vf)

# alpha_v<-rep(0,9)
# alpha_f<-rep(0,3)
# beta_vf <- matrix(0,9,3)
# beta_vf[c(1,2,3), 1] <- beta_vf[c(4,5,6), 2] <- beta_vf[c(7,8,9), 3] <- 1
# beta_vf[c(2,3),1]<-beta_vf[c(5,6),2]<-beta_vf[c(8,7),3]<-1
# value = list(alpha_v=alpha_v,beta_vf=beta_vf)
dta       <- lavaan::HolzingerSwineford1939

rc_sem <- lslSEM()
rc_sem$import(raw_obs=dta,var_subset = c(7:15))
#rc_sem$specify(pattern = list(beta_vf = beta_vf),value = list(alpha_v=alpha_v))
#rc_sem$specify(pattern,value)
rc_sem$specify(pattern)
rc_sem$learn(penalty="scad",gamma = 0.1,delta=2.5)


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


rc_sem2 <- lslSEM()
rc_sem2$import(raw_obs=dta,var_group=10)
rc_sem2$specify(pattern = list(beta_vf = beta_vf),ref_group = "g2")
rc_sem2$learn(penalty="mcp")



model<-rc_sem2$model
data<-rc_sem2$data

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



