library(dplyr);library(gtools)

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

#dta<-lavaan::HolzingerSwineford1939[,-c(1:6)]
dta <-lavaan::simulateData(model.cfa,sample.nobs = 500L)
n_obs<- ncol(dta)
n_lat<- 3
M    <- n_obs + n_lat

#matrices generation

eta       <- vector(mode = "numeric",M)
alpha     <- vector(mode = "numeric",M)
Beta      <- matrix(0, ncol = M, nrow = M)
eta       <- c(rep(0.5,9),rep(0.5,3))
Beta[1:9,10:12]<-0.5
zeta      <- vector(mode = "numeric",M)
Phi       <- matrix(diag(0.1,M,M), ncol = M, nrow = M)
ide       <- diag(1, ncol = M, nrow = M)
G_obs     <- c(rep(T,n_obs),rep(F,n_lat))
v         <- subset(eta,G_obs)
e_v    <- sapply(dta,mean)[1:n_obs]

#E-step

estep     <- function(alpha=alpha,Beta=Beta){

IBinv     <- solve(ide-Beta)
mu_eta    <- IBinv%*%alpha
mu_v      <- subset(mu_eta,G_obs)
Sigma_etaeta<-IBinv%*%Phi%*%t(IBinv)
Sigma_veta<- subset(Sigma_etaeta,G_obs)
Sigma_etav<- t(Sigma_veta)
Sigma_vv  <- subset(Sigma_etaeta,G_obs,G_obs)
J         <- mu_eta - Sigma_etav %*% solve(Sigma_vv) %*% mu_v
K         <- Sigma_etav %*% solve(Sigma_vv)
C_vv      <- cov(dta)
e_eta     <- J+K%*%e_v
C_etaeta  <- Sigma_etaeta - Sigma_etav %*% solve(Sigma_vv) %*% Sigma_veta +
  J %*% t(J) + J %*% t(e_v) %*% t(K) + K %*% e_v %*% t(J) + K %*% C_vv %*% t(K)
C_zetazeta<- C_etaeta - e_eta%*%alpha -  C_etaeta %*% t(Beta) - alpha %*% t(e_eta) +
  alpha %*% t(alpha) + alpha %*% t(e_eta) %*% t(Beta) - Beta %*% t(C_etaeta) + Beta %*% e_eta %*% t(alpha) +
  Beta %*% C_etaeta %*% t(Beta)
#zeta_tilda<- sapply(c(1:M), function(j) {solve(Phi[-j,-j]) %*% zeta[-j]})
C_zetatildazeta<- lapply(c(1:M), function(j) {solve(Phi[-j,-j])%*%C_zetazeta[-j,]})
C_zetatildazetatilda<- lapply(c(1:M), function(j) {solve(Phi[-j,-j])%*%C_zetazeta[-j,-j]%*%solve(Phi[-j,-j])})
return(list(e_eta=e_eta,
            C_etaeta=C_etaeta,
            C_zetazeta=C_zetazeta,
            C_zetatildazeta=C_zetatildazeta,
            C_zetatildazetatilda=C_zetatildazetatilda))
}


#CM

w_g       <- 1
alpha_u   <- vector(mode = "numeric",M)
alpha_hat <- vector(mode = "numeric",M)
JK        <- expand.grid(1:M,1:M)[2:1]
JLK       <- expand.grid(1:(M-1),1:M)[2:1]
Beta_u    <- matrix(0, M, M)
beta_hat  <- matrix(0, M, M)
Phi_u     <- matrix(0, M, M)
Phi_hat   <- matrix(0, M, M)
varphi_hat<- vector(mode = "numeric",M)

cmstep    <- function(w_g=w_g,
                      alpha_u=alpha_u,
                      Beta_u=Beta_u,
                      Phi=Phi,
                      e_eta=e_eta,C_etaeta=C_etaeta,
                      C_zetazeta=C_zetazeta,
                      C_zetatildazeta=C_zetatildazeta,
                      C_zetatildazetatilda=C_zetatildazetatilda){
  
  phi       <- solve(Phi)
  varphi    <- sapply(c(1:M), function(m) {diag(Phi)[m]-Phi[m,-m]%*%solve(Phi[-m,-m])%*%Phi[-m,m]})
  w_alpha   <- sapply(c(1:M), function(j) {1/(w_g*phi[j,j])})
  w_beta    <- mapply(function(j,k) 1/(w_g*solve(Phi[j,j])*C_etaeta[k,k]), j=JK[,1], k=JK[,2] ,SIMPLIFY = T) %>% matrix(nrow=M,byrow=T)
  diag(w_beta)<-0
  w_phi     <- matrix(0, M, M)
  w_phiq    <- mapply(function(j,lk) 1/((w_g/varphi[j])*C_zetatildazetatilda[[j]][lk,lk]),j=JLK[,1],lk=JLK[,2],SIMPLIFY = "matrix") %>% matrix(nrow=M,byrow=T)
  w_phi[upper.tri(w_phi)]<-w_phiq[upper.tri(w_phiq,diag=T)]
  w_phi[lower.tri(w_phi)]<-w_phiq[lower.tri(w_phiq)]
  
  for (j in 1:M){
    alpha_hat[j]   <- w_alpha[j] * ( w_g * phi[j,j] * (e_eta[j] - alpha_u[j] - Beta[j,] %*% e_eta) + 
                                    w_g * phi[j,-j] %*% (e_eta - alpha - alpha_u - Beta %*% e_eta)[-j])
  }
  
  for (j in 1:M){
    for (k in 1:M){
      beta_hat[j,k]<- w_beta[j,k] * (w_g * phi[j,j] * (C_etaeta[j,k] - e_eta[k] * alpha[j] - Beta[j,-k] %*% C_etaeta[k,-k] - Beta_u[j,] %*% C_etaeta[,k])+
                       w_g *  t(phi[j,-j]) %*% (C_etaeta[-j,k] - e_eta[k] * alpha[-j] - (Beta[-j,] + Beta_u[-j,]) %*% C_etaeta[,k]))
    }
  }
  
  for (j in 1:M){
    for ( k in 1:M){
      if (j==k) (Phi_hat[j,k]<-Phi[j,k])
      else {
       if (k<j) (lk<-k) else (lk<-c(k-1))
        Phi_hat[j,k]<- w_phi[j,k] * (w_g / phi[j,j]) * (C_zetatildazeta[[j]][lk,j] - Phi[j,-c(j,k)] %*% matrix(C_zetatildazetatilda[[j]][-lk,lk],c(M-2),1) - 
                                                       Phi_u[j,-j] %*% C_zetatildazetatilda[[j]][,lk])
      }
    }
  }

  for (j in 1:M){
    varphi_hat[j] <- w_g * (C_zetazeta[j,j] - 2 * Phi[j,-j] %*% C_zetatildazeta[[j]][,j] + Phi[j,-j] %*% C_zetatildazetatilda[[j]] %*% Phi[-j,j])
  }
  
  Beta     <- beta_hat
  alpha    <- alpha_hat
  Phi      <- Phi_hat
  varphi   <- varphi_hat

return(list(Beta=Beta,alpha=alpha,Phi=Phi,varphi=varphi))
}

# e_eta     <-estep(alpha=alpha,Beta=Beta)$e_eta
# C_etaeta  <-estep(alpha=alpha,Beta=Beta)$C_etaeta
# C_zetazeta<-estep(alpha=alpha,Beta=Beta)$C_zetazeta
# C_zetatildazeta<-estep(alpha=alpha,Beta=Beta)$C_zetatildazeta
# C_zetatildazetatilda<-estep(alpha=alpha,Beta=Beta)$C_zetatildazetatilda
# Beta      <-cmstep(w_alpha=w_alpha,w_g=w_g,e_eta=e_eta,alpha_u=alpha_u,Beta=Beta,alpha=alpha,C_etaeta=C_etaeta)$Beta
# alpha     <-cmstep(w_alpha=w_alpha,w_g=w_g,e_eta=e_eta,alpha_u=alpha_u,Beta=Beta,alpha=alpha,C_etaeta=C_etaeta)$alpha
# 
# tryy<-function(){
#   e_eta     <-estep(J=J,K=K,e_v=e_v,Sigma_etaeta=Sigma_etaeta,Sigma_vv=Sigma_vv,Sigma_veta=Sigma_veta,C_vv=C_vv,e_eta=e_eta,alpha=alpha,Beta=Beta)$e_eta
#   C_etaeta  <-estep(J=J,K=K,e_v=e_v,Sigma_etaeta=Sigma_etaeta,Sigma_vv=Sigma_vv,Sigma_veta=Sigma_veta,C_vv=C_vv,e_eta=e_eta,alpha=alpha,Beta=Beta)$C_etaeta
#   C_zetazeta<-estep(J=J,K=K,e_v=e_v,Sigma_etaeta=Sigma_etaeta,Sigma_vv=Sigma_vv,Sigma_veta=Sigma_veta,C_vv=C_vv,e_eta=e_eta,alpha=alpha,Beta=Beta)$C_zetazeta
#   Beta      <-cmstep(w_alpha=w_alpha,w_g=w_g,e_eta=e_eta,alpha_u=alpha_u,Beta=Beta,alpha=alpha,C_etaeta=C_etaeta)$Beta
#   alpha     <-cmstep(w_alpha=w_alpha,w_g=w_g,e_eta=e_eta,alpha_u=alpha_u,Beta=Beta,alpha=alpha,C_etaeta=C_etaeta)$alpha
#   list(C_zetazeta=C_zetazeta,Beta=Beta,alpha=alpha)
#   }
# 
# list(C_zetazeta=C_zetazeta,Beta=Beta,alpha=alpha)
