library(dplyr);library(gtools)

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
#varphi    <- sapply(c(1:M),function(i) {phi[i]-Phi[i,-i]%*%solve(Phi[-i,-i])%*%Phi[-i,i]} )
ide       <- diag(1, ncol = M, nrow = M)


#mu_v      <- mu_eta[1:n_obs]


G_obs     <- c(rep(T,n_obs),rep(F,n_lat))
v         <- subset(eta,G_obs)
e_v    <- sapply(dta,mean)[1:n_obs]

#E-step

#

estep     <- function(alpha_alpha,Beta=Beta){

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
phi       <- solve(Phi)
alpha_u   <- vector(mode = "numeric",M)
alpha_hat <- vector(mode = "numeric",M)
w_alpha   <- sapply(c(1:M), function(j) {1/(w_g*phi[j,j])})


JK        <- permutations(n = M, r = 2, v = 1:M , repeats.allowed = T)
w_beta    <- mapply(function(j,k) 1/(w_g*solve(Phi[j,j])*C_etaeta[k,k]), j=JK[,1], k=JK[,2] ,SIMPLIFY = T) %>% matrix(nrow=M,byrow=T)
Beta_u    <- matrix(0, ncol = M, nrow = M)
beta_hat  <- matrix(0, ncol = M, nrow = M)

cmstep    <- function(w_alpha=w_alpha,w_g=w_g,e_eta=e_eta,alpha_u=alpha_u,Beta=Beta,alpha=alpha,C_etaeta=C_etaeta){
  for (j in 1:M){
    alpha_hat[j]<- w_alpha[j] * ( w_g * phi[j,j] * (e_eta[j] - alpha_u[j] - Beta[j,] %*% e_eta) + 
                                    w_g * phi[j,-j] %*% (e_eta - alpha - alpha_u - Beta %*% e_eta)[-j])
  }
  
  for (j in 1:M){
    for (k in 1:M){
      beta_hat[j,k]<- w_beta[j,k] * (w_g * phi[j,j] * (C_etaeta[j,k] - e_eta[k] * alpha[j] - Beta[j,-k] %*% C_etaeta[k,-k] - Beta_u[j,] %*% C_etaeta[,k])+
                       w_g *  t(phi[j,-j]) %*% (C_etaeta[-j,k] - e_eta[k] * alpha[-j] - (Beta[-j,] + Beta_u[-j,]) %*% C_etaeta[,k]))
    }
  }
  
  Beta     <- beta_hat
  alpha    <- alpha_hat

return(list(Beta=Beta,alpha=alpha))
}

e_eta     <-estep(J=J,K=K,e_v=e_v,Sigma_etaeta=Sigma_etaeta,Sigma_vv=Sigma_vv,Sigma_veta=Sigma_veta,C_vv=C_vv,e_eta=e_eta,alpha=alpha,Beta=Beta)$e_eta
C_etaeta  <-estep(J=J,K=K,e_v=e_v,Sigma_etaeta=Sigma_etaeta,Sigma_vv=Sigma_vv,Sigma_veta=Sigma_veta,C_vv=C_vv,e_eta=e_eta,alpha=alpha,Beta=Beta)$C_etaeta
C_zetazeta<-estep(J=J,K=K,e_v=e_v,Sigma_etaeta=Sigma_etaeta,Sigma_vv=Sigma_vv,Sigma_veta=Sigma_veta,C_vv=C_vv,e_eta=e_eta,alpha=alpha,Beta=Beta)$C_zetazeta
Beta      <-cmstep(w_alpha=w_alpha,w_g=w_g,e_eta=e_eta,alpha_u=alpha_u,Beta=Beta,alpha=alpha,C_etaeta=C_etaeta)$Beta
alpha     <-cmstep(w_alpha=w_alpha,w_g=w_g,e_eta=e_eta,alpha_u=alpha_u,Beta=Beta,alpha=alpha,C_etaeta=C_etaeta)$alpha

tryy<-function(){
  e_eta     <-estep(J=J,K=K,e_v=e_v,Sigma_etaeta=Sigma_etaeta,Sigma_vv=Sigma_vv,Sigma_veta=Sigma_veta,C_vv=C_vv,e_eta=e_eta,alpha=alpha,Beta=Beta)$e_eta
  C_etaeta  <-estep(J=J,K=K,e_v=e_v,Sigma_etaeta=Sigma_etaeta,Sigma_vv=Sigma_vv,Sigma_veta=Sigma_veta,C_vv=C_vv,e_eta=e_eta,alpha=alpha,Beta=Beta)$C_etaeta
  C_zetazeta<-estep(J=J,K=K,e_v=e_v,Sigma_etaeta=Sigma_etaeta,Sigma_vv=Sigma_vv,Sigma_veta=Sigma_veta,C_vv=C_vv,e_eta=e_eta,alpha=alpha,Beta=Beta)$C_zetazeta
  Beta      <-cmstep(w_alpha=w_alpha,w_g=w_g,e_eta=e_eta,alpha_u=alpha_u,Beta=Beta,alpha=alpha,C_etaeta=C_etaeta)$Beta
  alpha     <-cmstep(w_alpha=w_alpha,w_g=w_g,e_eta=e_eta,alpha_u=alpha_u,Beta=Beta,alpha=alpha,C_etaeta=C_etaeta)$alpha
  list(C_zetazeta=C_zetazeta,Beta=Beta,alpha=alpha)
  }

list(C_zetazeta=C_zetazeta,Beta=Beta,alpha=alpha)
