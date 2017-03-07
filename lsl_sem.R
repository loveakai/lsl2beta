rm(list=ls())
set.seed=4869
library(dplyr);library(gtools)

source('/Volumes/phaksie/Dropbox/lsl2_beta/lsl_tool.R')
source('/Volumes/phaksie/Dropbox/lsl2_beta/lsl_functions.R')

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
F1~~0.4*F2
F1~~0.4*F3
F2~~0.4*F3
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
F1~~0.4*F2
F1~~0.4*F3
F2~~0.4*F3
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
F1~~0.4*F2
F1~~0.4*F3
F2~~0.4*F3
'

dta       <-list()
dta[[1]]  <- lavaan::simulateData(model.cfa,sample.nobs = 10000L) #%>% cbind(.,sample(c(1,2),size=nrow(.),rep=T))
dta[[2]]  <- lavaan::simulateData(model.cfa2,sample.nobs = 10000L)
dta[[3]]  <- lavaan::simulateData(model.cfa3,sample.nobs = 10000L)
#dta       <- lavaan::HolzingerSwineford1939[7:15]

n_groups  <- length(dta)
n_v       <- ncol(dta[[1]])
n_f       <- 3
n_eta     <- n_v + n_f
nm        <- c(paste0("v",1:n_v),paste0("f",1:n_f))
sigma     <- lapply(1:n_groups, function(i_groups) {dta[[i_groups]] %>% as.matrix %>% crossprod %>% "/"(nrow(dta[[i_groups]])) %>%
    `colnames<-`(nm[1:n_v]) %>% `rownames<-`(colnames(.))})
e_v       <- lapply(1:n_groups, function(i_groups) sapply(dta[[i_groups]],mean)[1:n_v] %>% `names<-`(nm[1:n_v]))


lambda <- matrix(NA, 9, 3)
lambda[c(1,2,3), 1] <- lambda[c(4,5,6), 2] <- lambda[c(7,8,9), 3] <- 1

# beta_p    <- matrix(0, ncol = n_eta, nrow = n_eta) %>% `colnames<-`(nm) %>% `rownames<-`(nm)
# beta_p[c(1,2,3), 10] <- beta_p[c(4,5,6), 11] <- beta_p[c(7,8,9), 12] <- 1  
# Beta      <- Beta <- 0.8*.is_one(beta_p) #starting value of Beta
# Beta[c(2,3), 10] <- Beta[c(5,6), 11] <- Beta[c(8,9), 12] <- 1  

#phi_p     <- matrix(0,n_eta,n_eta)
#phi_p[c(11,12),10] <- phi_p[c(10,12),11] <- phi_p[c(10,11),12] <- NA
# Phi       <- matrix(0,n_eta,n_eta)
# Phi[(n_v+1):n_eta,(n_v+1):n_eta] <- 0.4
# diag(Phi) <- 1-0.8^2
# Phi[10,10]<-Phi[11,11]<-Phi[12,12]<-1

#mat       <- matgen(lambda=lambda,Beta = Beta,scale=T)

mat       <- matgen(lambda=lambda)

eta       <- vector(mode = "numeric",n_eta)   %>%`names<-`(nm)
zeta      <- vector(mode = "numeric",n_eta)   %>%`names<-`(nm)
ide       <- diag(1, ncol = n_eta, nrow = n_eta)  %>% `colnames<-`(nm) %>% `rownames<-`(nm) 
G_eta     <- c(rep(T,n_v),rep(F,n_f)) %>%`names<-`(nm)
v         <- subset(eta,G_eta)
penalize  <- list(type="L1",delta=0.025,gamma=2.5)

pl <-penalize

#ECM

ecmm<-ecm(mat=mat,ide=ide,G_eta=G_eta,maxit=500,cri=10^(-5),penalize=pl)
