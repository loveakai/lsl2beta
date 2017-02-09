betagen   <- function(lambda){
  
  Beta_p                          <- matrix(0, ncol = M, nrow = M)
  Beta_p[(1:n_obs),(n_obs+1):M]   <- lambda
  Beta_p[(n_obs+1):M,(n_obs+1):M] <- diag(1,n_lat,n_lat)
  return(Beta_p)
  
}

matgen    <- function(alpha_p,Beta_p,Phi_p,alpha,Beta,Phi,lambda,scale=T){
  
  #pattern matarices generation
  #only 3 matrices have pattern matrix, including alpha, beta, Phi
    
  if (missing(alpha_p)) {   alpha_p <- c(rep(1,n_obs),rep(0,n_lat)) %>% `names<-`(nm) }
  if (missing(Beta_p))  {
    if (missing(lambda)) { 
      stop("lambda matrix is not specified")
    } else { 
      Beta_p                        <- betagen(lambda) 
    }
  }
  if (missing(Phi_p)) {
    Phi_p                           <- matrix(diag(1,M,M), ncol = M, nrow = M)
    Phi_p[(n_obs+1):M,(n_obs+1):M]  <- 1
    colnames(Phi_p)                 <- nm
    rownames(Phi_p)                 <- nm
  }
  
  #matrices generation
  
  if (missing(alpha)) {
    alpha <- c(sapply(dta,mean)[1:n_obs],rep(0,n_lat)) %>% `names<-`(nm)}
  if (missing(Beta))  {
    Beta  <- 0.1*.is_one(Beta_p)
    colnames(Beta)  <- nm
    rownames(Beta)  <- nm
    }
  if (missing(Phi))   {
    Phi   <- matrix(diag(c(rep(0.1,n_obs,n_obs),rep(1,n_lat,n_lat))), ncol = M, nrow = M)
    colnames(Phi)   <- nm
    rownames(Phi)   <- nm
    }
  if (scale) {Beta_p[apply(Beta_p[(1:n_obs),(n_obs+1):M] == 1, 2, function(x) min(which(x))) %>% cbind((n_obs+1):M)] <- 0}　
  
  
  return(list(pattern=list(alpha_p=alpha_p,Beta_p=Beta_p,Phi_p=Phi_p),value=list(alpha=alpha,Beta=Beta,Phi=Phi)))
  
}

threshold <- function(theta,gamma){
   sign(theta)*max(abs(theta)-gamma,0) 
}

varphi    <- function(m,Phi) {diag(Phi)[m]-Phi[m,-m]%*%solve(Phi[-m,-m])%*%Phi[-m,m]}

penalty   <- function(theta,gamma,cw,delta,type){
  if (type == "l1") {
    theta <- threshold(theta, gamma * cw)
  } else if (type == "SCAD") {
    if (abs(theta) <= gamma * (1 + cw)) {
      theta <- threshold(theta, cw * gamma)
    } else if (gamma * (1 + cw) < theta & theta <= gamma * delta) {
      theta <-
        thereshold(theta, (cw * gamma * delta) / (delta - 1)) * solve(1 - cw / (delta - 1))
    } else if (gamma * delta < theta) {
      theta <- theta
    }
  } else if (type == "MCP") {
    if (theta <= gamma * delta) {
      theta <- threshold(theta, cw * gamma) * solve(1 - cw / delta)
    } else if (gamma * delta < theta) {
      theta <- theta
    }
  }
  return(theta)
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
  
  return(list(e_eta=e_eta,
              C_etaeta=C_etaeta))
}

cmstep    <- function(w_g=w_g,JK=JK,JLK=JLK,alpha_u=alpha_u,Beta_u=Beta_u,Phi_u=Phi_u,mat=ini$mat,e_step=e_step){
  
  e_eta     <- e_step$e_eta
  C_etaeta  <- e_step$C_etaeta
  Phi       <- mat$value$Phi
  Beta      <- mat$value$Beta
  alpha     <- mat$value$alpha
  phi       <- solve(Phi) # for alpha and Beta
  M         <- length(mat$pattern$alpha_p)

  ## reference components updating

  # alpha
  
  w_alpha   <- sapply(c(1:M), function(j) {1/(w_g*phi[j,j])})
  for (j in which(.is_est(mat$pattern$alpha_p))){
    alpha[j]   <- w_alpha[j] * ( w_g * phi[j,j] * (e_eta[j] - alpha_u[j] - Beta[j,] %*% e_eta) +
                                   w_g * phi[j,-j] %*% (e_eta - alpha - alpha_u - Beta %*% e_eta)[-j])
  }

  # Beta
  
  w_beta    <- mapply(function(j,k) 1/(w_g*phi[j,j]*C_etaeta[k,k]), j=JK[,1], k=JK[,2] ,SIMPLIFY = T) %>% matrix(nrow=M,byrow=T)
  diag(w_beta)<-0
  for (i in which(.is_est(mat$pattern$Beta_p))){
    k      <- ceiling(i/M)
    j      <- i-(k-1)*M
    Beta[j,k]<- w_beta[j,k] * (w_g * phi[j,j] * (C_etaeta[j,k] - e_eta[k] * alpha[j] - Beta[j,-k] %*% C_etaeta[k,-k] - Beta_u[j,] %*% C_etaeta[,k])+
                                 w_g * phi[j,-j] %*% (C_etaeta[-j,k] - e_eta[k] * alpha[-j] - (Beta[-j,] + Beta_u[-j,]) %*% C_etaeta[,k]))
    
  }
  
  # Phi
  
  C_zetazeta<- C_etaeta - e_eta%*%alpha -  C_etaeta %*% t(Beta) - alpha %*% t(e_eta) +
    alpha %*% t(alpha) + alpha %*% t(e_eta) %*% t(Beta) - Beta %*% t(C_etaeta) + 
    Beta %*% e_eta %*% t(alpha) + Beta %*% C_etaeta %*% t(Beta)
  
  for (i in which(.is_est(mat$pattern$Phi_p))){
    k      <- ceiling(i/M)
    j      <- i-(k-1)*M
    if (j<=k) {} else {
      lk<-k
      
      C_zetatildazeta     <- solve(Phi[-j,-j])%*%C_zetazeta[-j,]
      C_zetatildazetatilda<- solve(Phi[-j,-j])%*%C_zetazeta[-j,-j]%*%solve(Phi[-j,-j])
      var_phi   <- sapply(c(1:M), varphi, Phi)
      w_phi     <- 1/((w_g/var_phi[j])*C_zetatildazetatilda[lk,lk])
      
      Phi[j,k]<- w_phi * (w_g / var_phi[j]) * (C_zetatildazeta[lk,j] - Phi[j,-c(j,k)] %*% matrix(C_zetatildazetatilda[-lk,lk],c(M-2),1) - 
                                                 Phi_u[j,-j] %*% C_zetatildazetatilda[,lk])
    }
  }
  Phi[upper.tri(Phi)]<-t(Phi)[upper.tri(Phi)]


  for (j in 1:M){
    C_zetatildazeta<- solve(Phi[-j,-j])%*%C_zetazeta[-j,]
    C_zetatildazetatilda<- solve(Phi[-j,-j])%*%C_zetazeta[-j,-j]%*%solve(Phi[-j,-j])
    Phi[j,j] <- w_g * (C_zetazeta[j,j] - 2 * Phi[j,-j] %*% C_zetatildazeta[,j] + Phi[j,-j] %*% C_zetatildazetatilda %*% Phi[-j,j]) +
      Phi[j,-j] %*% solve(Phi[-j,-j]) %*% Phi[-j,j]
  }

  ## increment components weights
  w_alpha_u <- sapply(c(1:M), function(j) {1/w_g*phi[j,j]})
  
  
  
  return(list(alpha=alpha,
              Beta=Beta,
              Phi=Phi))
} 

dml_cal   <- function(Sigma=Sigma,e_v=e_v,Sigma_vv=subset(ini$Sigma_etaeta,G_obs,G_obs),mu_v=subset(ini$mu_eta,G_obs)){
  Sigma_vv_iv <- solve(Sigma_vv)
  dml <- -log(det(Sigma_vv_iv %*% Sigma)) + sum(diag(Sigma_vv_iv%*%Sigma)) - dim(Sigma_vv)[1] + t(e_v-mu_v) %*% Sigma_vv_iv %*% (e_v-mu_v)
  return(dml)
  }
