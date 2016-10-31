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
dta       <- lavaan::simulateData(model.cfa,sample.nobs = 1000L)
n_obs     <- ncol(dta)
n_lat     <- 3
M         <- n_obs + n_lat
Sigma     <- cov(dta)
e_v       <- sapply(dta,mean)[1:n_obs]

matgen    <- function(alpha_p,Beta_p,Phi_p,alpha,Beta,Phi,scale=T){

              #pattern matarices generation
              #only 3 matrices have pattern matrix, including alpha, beta, Phi
              
                if (missing(alpha_p)) {alpha_p <- c(rep(0,n_obs),rep(1,n_lat))}
                if (missing(Beta_p)) {
                  Beta_p              <- matrix(0, ncol = M, nrow = M) 
                  Beta_p[1:n_obs,(n_obs+1):M]   <-1
                  Beta_p[(n_obs+1):M,(n_obs+1):M] <-diag(1,n_lat,n_lat)
                }
                if (missing(Phi_p)) {Phi_p <- matrix(diag(1,M,M), ncol = M, nrow = M)}
                
              #matrices generation
              
                if (missing(alpha)) {alpha <- c(sapply(dta,mean)[1:n_obs],rep(0,n_lat))}
                if (missing(Beta)) {Beta <- 0.1*.is_one(Beta_p)}
                if (missing(Phi)) {Phi <- matrix(diag(c(rep(0.1,n_obs,n_obs),rep(1,n_lat,n_lat))), ncol = M, nrow = M)}
                if (scale) {Beta_p[apply(Beta_p[(1:n_obs),(n_obs+1):M] == 1, 2, function(x) min(which(x))) %>% cbind((n_obs+1):M)] <- 0}　
                
                return(list(pattern=list(alpha_p=alpha_p,Beta_p=Beta_p,Phi_p=Phi_p),value=list(alpha=alpha,Beta=Beta,Phi=Phi)))
                
}

estep     <- function(ini){
              mu_v      <- subset(ini$mu_eta,ini$G_obs)
              Sigma_veta<- subset(ini$Sigma_etaeta,ini$G_obs)
              Sigma_etav<- t(Sigma_veta)
              Sigma_vv  <- subset(ini$Sigma_etaeta,ini$G_obs,ini$G_obs)
              J         <- ini$mu_eta - Sigma_etav %*% solve(Sigma_vv) %*% mu_v
              K         <- Sigma_etav %*% solve(Sigma_vv)
              C_vv      <- ini$Sigma
              e_eta     <- J+K%*%ini$e_v
              C_etaeta  <- ini$Sigma_etaeta - Sigma_etav %*% solve(Sigma_vv) %*% Sigma_veta +
                J %*% t(J) + J %*% t(ini$e_v) %*% t(K) + K %*% ini$e_v %*% t(J) + K %*% C_vv %*% t(K)
              alpha     <- ini$mat$value$alpha
              Beta      <- ini$mat$value$Beta
              Phi       <- ini$mat$value$Phi
              M         <- length(ini$G_obs)
              C_zetazeta<- C_etaeta - e_eta%*%alpha -  C_etaeta %*% t(Beta) - alpha %*% t(e_eta) +
                alpha %*% t(alpha) + alpha %*% t(e_eta) %*% t(Beta) - Beta %*% t(C_etaeta) + 
                Beta %*% e_eta %*% t(alpha) + Beta %*% C_etaeta %*% t(Beta)
              C_zetatildazeta<- lapply(c(1:M), function(j) {solve(Phi[-j,-j])%*%C_zetazeta[-j,]})
              C_zetatildazetatilda<- lapply(c(1:M), function(j) {solve(Phi[-j,-j])%*%C_zetazeta[-j,-j]%*%solve(Phi[-j,-j])})
              return(list(e_eta=e_eta,
                          C_etaeta=C_etaeta,
                          C_zetazeta=C_zetazeta,
                          C_zetatildazeta=C_zetatildazeta,
                          C_zetatildazetatilda=C_zetatildazetatilda))
}

cmstep    <- function(w_g=w_g,alpha_u=alpha_u,Beta_u=Beta_u,mat=mat,e_step=e_step){
  
  e_eta     <-e_step$e_eta
  C_etaeta  <-e_step$C_etaeta
  C_zetazeta<-e_step$C_zetazeta
  C_zetatildazeta<-e_step$C_zetatildazeta
  C_zetatildazetatilda<-e_step$C_zetatildazetatilda
  phi       <- solve(mat$value$Phi)
  varphi    <- sapply(c(1:M), function(m) {diag(Phi)[m]-Phi[m,-m]%*%solve(Phi[-m,-m])%*%Phi[-m,m]})
  w_alpha   <- sapply(c(1:M), function(j) {1/(w_g*phi[j,j])})
  w_beta    <- mapply(function(j,k) 1/(w_g*solve(Phi[j,j])*C_etaeta[k,k]), j=JK[,1], k=JK[,2] ,SIMPLIFY = T) %>% matrix(nrow=M,byrow=T)
  diag(w_beta)<-0
  w_phi     <- matrix(0, M, M)
  w_phiq    <- mapply(function(j,lk) 1/((w_g/varphi[j])*C_zetatildazetatilda[[j]][lk,lk]),j=JLK[,1],lk=JLK[,2],SIMPLIFY = "matrix") %>% matrix(nrow=M,byrow=T)
  w_phi[upper.tri(w_phi)]<-w_phiq[upper.tri(w_phiq,diag=T)]
  w_phi[lower.tri(w_phi)]<-w_phiq[lower.tri(w_phiq)]
  
  for (j in which(.is_est(alpha_p))){
    alpha_hat[j]   <- w_alpha[j] * ( w_g * phi[j,j] * (e_eta[j] - alpha_u[j] - Beta[j,] %*% e_eta) + 
                                       w_g * phi[j,-j] %*% (e_eta - alpha - alpha_u - Beta %*% e_eta)[-j])
  }
  
  for (j in which(.is_est(Beta_p),arr.ind = T)[,1]){
    for (k in which(.is_est(Beta_p),arr.ind = T)[,2]){
      beta_hat[j,k]<- w_beta[j,k] * (w_g * phi[j,j] * (C_etaeta[j,k] - e_eta[k] * alpha[j] - Beta[j,-k] %*% C_etaeta[k,-k] - Beta_u[j,] %*% C_etaeta[,k])+
                                       w_g *  t(phi[j,-j]) %*% (C_etaeta[-j,k] - e_eta[k] * alpha[-j] - (Beta[-j,] + Beta_u[-j,]) %*% C_etaeta[,k]))
    }
  }
  
  if (any(!.is_est(Phi_p))) {Phi_hat<-Phi}
  else {
    for (j in which(.is_est(Phi_p),arr.ind = T)[,1]){
      for ( k in which(.is_est(Phi_p),arr.ind = T)[,2]){
        if (j==k) (Phi_hat[j,k]<-Phi[j,k])
        else {
          if (k<j) (lk<-k) else (lk<-c(k-1))
          Phi_hat[j,k]<- w_phi[j,k] * (w_g / phi[j,j]) * (C_zetatildazeta[[j]][lk,j] - Phi[j,-c(j,k)] %*% matrix(C_zetatildazetatilda[[j]][-lk,lk],c(M-2),1) - 
                                                            Phi_u[j,-j] %*% C_zetatildazetatilda[[j]][,lk])
        }
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

Beta_p    <- matrix(0, ncol = M, nrow = M)
Beta_p[c(1,2,3), 10] <- Beta_p[c(4,5,6), 11] <- Beta_p[c(7,8,9), 12] <- 1


mat       <- matgen(Beta_p = Beta_p)

eta       <- vector(mode = "numeric",M)
eta       <- c(rep(0.5,9),rep(0.5,3))
zeta      <- vector(mode = "numeric",M)
ide       <- diag(1, ncol = M, nrow = M)
G_obs     <- c(rep(T,n_obs),rep(F,n_lat))
v         <- subset(eta,G_obs)

#ECM

ecm       <- function(mat=mat,ide=ide,G_obs=G_obs){
            alpha_p   <- mat$pattern$alpha_p
            Beta_p    <- mat$pattern$Beta_p
            Phi_p     <- mat$pattern$Phi_p
            alpha     <- mat$value$alpha
            Beta      <- mat$value$Beta
            Phi       <- mat$pattern$Phi_p
            #model-implied matrix
            
            IBinv     <- solve(ide-Beta)
            mu_eta    <- IBinv%*%alpha
            Sigma_etaeta<-IBinv%*%Phi%*%t(IBinv)
            
            Sigma     <- Sigma
            e_v       <- e_v
            
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
            
            ini       <- list(IBinv=IBinv,mu_eta=mu_eta,Sigma_etaeta=Sigma_etaeta,G_obs=G_obs,Sigma=Sigma,e_v=e_v,mat=mat)
            
            for (it in 1:1000){
            e_step    <- estep(ini)
            cm_step   <- cmstep(w_g=w_g,alpha_u=alpha_u,Beta_u=Beta_u,mat=mat,e_step=e_step)
            ini$mat$value$Beta<- cm_step$Beta
            ini$mat$value$alpha<-cm_step$alpha
            ini$mat$value$Phi <- cm_step$Phi
            }
            
}

ls(e_step)

 
#  tryy<-function(){
#    
#    e_eta     <-estep(alpha=alpha,Beta=Beta)$e_eta
#    C_etaeta  <-estep(alpha=alpha,Beta=Beta)$C_etaeta
#    C_zetazeta<-estep(alpha=alpha,Beta=Beta)$C_zetazeta
#    C_zetatildazeta<-estep(alpha=alpha,Beta=Beta)$C_zetatildazeta
#    C_zetatildazetatilda<-estep(alpha=alpha,Beta=Beta)$C_zetatildazetatilda
#    Beta      <-cmstep(w_g=w_g,alpha_u=alpha_u,Beta_u=Beta_u,Phi=Phi,e_eta=e_eta,C_etaeta=C_etaeta,C_zetazeta=C_zetazeta,C_zetatildazeta=C_zetatildazeta,
#                       C_zetatildazetatilda=C_zetatildazetatilda)$Beta
#    alpha     <-cmstep(w_g=w_g,alpha_u=alpha_u,Beta_u=Beta_u,Phi=Phi,e_eta=e_eta,C_etaeta=C_etaeta,C_zetazeta=C_zetazeta,C_zetatildazeta=C_zetatildazeta,
#                       C_zetatildazetatilda=C_zetatildazetatilda)$alpha
#    Phi       <-cmstep(w_g=w_g,alpha_u=alpha_u,Beta_u=Beta_u,Phi=Phi,e_eta=e_eta,C_etaeta=C_etaeta,C_zetazeta=C_zetazeta,C_zetatildazeta=C_zetatildazeta,
#                       C_zetatildazetatilda=C_zetatildazetatilda)$Phi
#    varphi    <-cmstep(w_g=w_g,alpha_u=alpha_u,Beta_u=Beta_u,Phi=Phi,e_eta=e_eta,C_etaeta=C_etaeta,C_zetazeta=C_zetazeta,C_zetatildazeta=C_zetatildazeta,
#                       C_zetatildazetatilda=C_zetatildazetatilda)$varphi
#    
#    list(e_eta=e_eta,C_etaeta=C_etaeta,C_zetazeta=C_zetazeta,C_zetatildazeta=C_zetatildazeta,
#         C_zetatildazetatilda=C_zetatildazetatilda,Beta=Beta,alpha=alpha,Phi=Phi,varphi=varphi)
#    }
#  
# # list(C_zetazeta=C_zetazeta,Beta=Beta,alpha=alpha)
#  tryy()$Beta
 