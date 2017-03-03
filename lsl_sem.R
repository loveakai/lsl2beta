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
#dta[[2]]  <- lavaan::simulateData(model.cfa2,sample.nobs = 10000L)
#dta[[3]]  <- lavaan::simulateData(model.cfa3,sample.nobs = 10000L)
#dta       <- lavaan::HolzingerSwineford1939[7:15]

n_gps     <- length(dta)
n_obs     <- ncol(dta[[1]])
n_lat     <- 3
M         <- n_obs + n_lat
nm        <- c(paste0("v",1:n_obs),paste0("f",1:n_lat))
Sigma     <- lapply(1:n_gps, function(x) {dta[[x]] %>% as.matrix %>% t %>% as.data.frame %>% lapply(tcrossprod) %>% simplify2array %>% apply(1:2,mean) %>%
    `colnames<-`(nm[1:n_obs]) %>% `rownames<-`(colnames(.))})
e_v       <- lapply(1:n_gps, function(x) sapply(dta[[x]],mean)[1:n_obs] %>% `names<-`(nm[1:n_obs]))


lambda <- matrix(NA, 9, 3)
lambda[c(1,2,3), 1] <- lambda[c(4,5,6), 2] <- lambda[c(7,8,9), 3] <- 1

# Beta_p    <- matrix(0, ncol = M, nrow = M) %>% `colnames<-`(nm) %>% `rownames<-`(nm)
# Beta_p[c(1,2,3), 10] <- Beta_p[c(4,5,6), 11] <- Beta_p[c(7,8,9), 12] <- 1  
# Beta      <- Beta <- 0.8*.is_one(Beta_p) #starting value of Beta
# Beta[c(2,3), 10] <- Beta[c(5,6), 11] <- Beta[c(8,9), 12] <- 1  

Phi_p     <- matrix(0,M,M)
Phi_p[c(11,12),10] <- Phi_p[c(10,12),11] <- Phi_p[c(10,11),12] <- NA
# Phi       <- matrix(0,M,M)
# Phi[(n_obs+1):M,(n_obs+1):M] <- 0.4
# diag(Phi) <- 1-0.8^2
# Phi[10,10]<-Phi[11,11]<-Phi[12,12]<-1

#mat       <- matgen(lambda=lambda,Beta = Beta,scale=T)

mat       <- matgen(lambda=lambda,Phi_p = Phi_p)

eta       <- vector(mode = "numeric",M)   %>%`names<-`(nm)
zeta      <- vector(mode = "numeric",M)   %>%`names<-`(nm)
ide       <- diag(1, ncol = M, nrow = M)  %>% `colnames<-`(nm) %>% `rownames<-`(nm) 
G_obs     <- c(rep(T,n_obs),rep(F,n_lat)) %>%`names<-`(nm)
v         <- subset(eta,G_obs)

#ECM

#ecm       <- function(mat=mat,ide=ide,G_obs=G_obs){
            alpha_p   <- mat$pattern$alpha_p
            Beta_p    <- mat$pattern$Beta_p
            Phi_p     <- mat$pattern$Phi_p
            alpha_u   <- mat$value$alpha_u
            Beta_u    <- mat$value$Beta_u
            Phi_u     <- mat$value$Phi_u
            alpha_g   <- mat$value$alpha_g
            Beta_g    <- mat$value$Beta_g
            Phi_g     <- mat$value$Phi_g
            
            #model-implied matrix
            
            IBinv     <- lapply(1:n_gps, function(x) solve(ide-(Beta_u+Beta_g[[x]])))
            mu_eta    <- lapply(1:n_gps, function(x) IBinv[[x]]%*%(alpha_u+alpha_g[[x]]))
            Sigma_etaeta<-lapply(1:n_gps, function(x) IBinv[[x]]%*%(Phi_u+Phi_g[[x]])%*%t(IBinv[[x]]))
            
            w_g       <- rep(list(0.5),n_gps)
            #alpha_g   <- rep(list(rep(0,M)),n_gps)
            JK        <- expand.grid(1:M,1:M)[2:1]
            JLK       <- expand.grid(1:(M-1),1:M)[2:1]
            #Beta_g    <- rep(list(matrix(0, M, M)),n_gps)
            #Phi_g     <- matrix(0, M, M) %>% `diag<-`(0.1) %>%  list() %>%  rep(n_gps)
            
            ini       <- list(IBinv=IBinv,mu_eta=mu_eta,Sigma_etaeta=Sigma_etaeta,Sigma=Sigma,G_obs=G_obs,e_v=e_v,mat=mat)
 
            for (it in 1:500){
              print(it)
              e_step    <- estep(ini)
              cm_step   <- cmstep(w_g=w_g,JK=JK,JLK=JLK,mat=ini$mat,e_step=e_step,type="L1")
              ini$mat$value$Phi_u   <-Phi_u     <- cm_step$Phi_u
              ini$mat$value$Beta_u  <-Beta_u    <- cm_step$Beta_u
              ini$mat$value$alpha_u <-alpha_u   <- cm_step$alpha_u
              ini$mat$value$Phi_g   <-Phi_g     <- cm_step$Phi_g
              ini$mat$value$Beta_g  <-Beta_g    <- cm_step$Beta_g
              ini$mat$value$alpha_g <-alpha_g   <- cm_step$alpha_g
              ini$IBinv          <- lapply(1:n_gps, function(x) solve(ide-(Beta_u+Beta_g[[x]])))
              ini$mu_eta         <- lapply(1:n_gps, function(x) ini$IBinv[[x]]%*%(alpha_u+alpha_g[[x]]))
              ini$Sigma_etaeta   <- lapply(1:n_gps, function(x) ini$IBinv[[x]]%*%(Phi_u+Phi_g[[x]])%*%t(ini$IBinv[[x]]))
            }
            
            theta     <- c(cm_step$alpha[.is_est(ini$mat$pattern$alpha_p)],cm_step$Beta[.is_est(ini$mat$pattern$Beta_p)],cm_step$Phi[.is_est(ini$mat$pattern$Phi_p)])
            length(theta)
            
            dml_cal(Sigma=Sigma,e_v=e_v,Sigma_vv=subset(ini$Sigma_etaeta,G_obs,G_obs),mu_v=subset(ini$mu_eta,G_obs))
            
          
#}



