betagen   <- function(lambda){
  
  Beta_p                          <- matrix(0, ncol = M, nrow = M)
  Beta_p[(1:n_obs),(n_obs+1):M]   <- lambda
  #Beta_p[(n_obs+1):M,(n_obs+1):M] <- diag(1,n_lat,n_lat)
  return(Beta_p)
  
}

matgen    <- function(alpha_p,Beta_p,Phi_p,alpha,Beta,Phi,lambda,scale=T){ #**diminfo could be modified

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
  ## reference component
  if (missing(alpha)) {
    alpha_u <- c(sapply(dta,mean)[1:n_obs],rep(0,n_lat))
  } else {
    alpha_u<-alpha }
  
  if (missing(Beta))  {
    Beta_u  <- 1 *.is_one(Beta_p)
  } else {
    Beta_u <-Beta} 
  
  if (missing(Phi))   {
    #Phi_u   <- matrix(diag(c(rep(0.1,n_obs,n_obs),rep(1,n_lat,n_lat))), ncol = M, nrow = M)
    Phi_u <- diag(0,M,M)
    } else {
    Phi_u <- Phi
  }
  
  rownames(Phi_u)   <- colnames(Phi_u)   <- rownames(Beta_u)  <- colnames(Beta_u)  <-names(alpha_u)      <- nm

  if (scale) {Beta_p[apply(Beta_p[(1:n_obs),(n_obs+1):M] == 1, 2, function(x) min(which(x))) %>% cbind((n_obs+1):M)] <- 0}　
  
  ## increment component
  alpha_g<-rep(list(0*alpha_u),n_gps)
  Beta_g <-rep(list(0*Beta_u),n_gps)
  
  g<-0.1*diag(ncol(Phi_u))
  g[(n_obs+1):M,(n_obs+1):M]<-(diag(n_lat))
  Phi_g  <-rep(list(g %>%`colnames<-`(nm) %>% `rownames<-`(nm)),n_gps) 
  
  return(list(pattern=list(alpha_p=alpha_p,Beta_p=Beta_p,Phi_p=Phi_p),value=list(alpha_u=alpha_u,Beta_u=Beta_u,Phi_u=Phi_u,alpha_g=alpha_g,Beta_g=Beta_g,Phi_g=Phi_g)))
  
}

threshold <- function(theta,gma){
   sign(theta)*max(abs(theta)-gma,0) 
}

varphi    <- function(m,Phi) {diag(Phi)[m]-Phi[m,-m]%*%solve(Phi[-m,-m])%*%Phi[-m,m]}

penalty   <- function(theta,gamma,cth,w,delta,type){
  if (type == "L1") {
    theta <- threshold(theta, gamma * cth * w)
  } else if (type == "SCAD") {
    if (abs(theta) <= gamma * (1 + cth * w)) {
      theta <- threshold(theta, cth * w * gamma)
    } else if (gamma * (1 + cth * w) < abs(theta) & abs(theta) <= gamma * delta) {
      theta <- threshold(theta, (cth * w * gamma * delta) / (delta - 1)) / (1 - ((cth * w) / (delta - 1)))
    } else { }
  } else if (type == "MCP") {
    if (abs(theta) <= gamma * delta) {
      theta <- threshold(theta, cth * w * gamma) / (1 - ((cth * w) / delta))
    } else { }
  }
  return(theta)
}

estep     <- function(ini){
  mu_v      <- lapply(1:n_gps, function(x) subset(ini$mu_eta[[x]],ini$G_obs))
  Sigma_veta<- lapply(1:n_gps, function(x) subset(ini$Sigma_etaeta[[x]],ini$G_obs))
  Sigma_etav<- lapply(1:n_gps, function(x) t(Sigma_veta[[x]]))
  Sigma_vv  <- lapply(1:n_gps, function(x) subset(ini$Sigma_etaeta[[x]],ini$G_obs,ini$G_obs))
  J         <- lapply(1:n_gps, function(x) ini$mu_eta[[x]] - Sigma_etav[[x]] %*% solve(Sigma_vv[[x]]) %*% mu_v[[x]])
  K         <- lapply(1:n_gps, function(x) Sigma_etav[[x]] %*% solve(Sigma_vv[[x]]))
  C_vv      <- lapply(1:n_gps, function(x) ini$Sigma[[x]])
  e_eta     <- lapply(1:n_gps, function(x) J[[x]]+K[[x]]%*%ini$e_v[[x]])
  C_etaeta  <- lapply(1:n_gps, function(x) {ini$Sigma_etaeta[[x]] - Sigma_etav[[x]] %*% solve(Sigma_vv[[x]]) %*% Sigma_veta[[x]] +
    J[[x]] %*% t(J[[x]]) + J[[x]] %*% t(ini$e_v[[x]]) %*% t(K[[x]]) + K[[x]] %*% ini$e_v[[x]] %*% t(J[[x]]) + K[[x]] %*% C_vv[[x]] %*% t(K[[x]])})
  
  return(list(e_eta=e_eta, C_etaeta=C_etaeta))
}

cmstep    <- function(w_g=w_g,JK=JK,JLK=JLK,mat=ini$mat,e_step=e_step,type=type){
  if (missing(type)) {type<-"L1"}
  
  e_eta     <- e_step$e_eta
  C_etaeta  <- e_step$C_etaeta
  Phi_u     <- mat$value$Phi_u
  Beta_u    <- mat$value$Beta_u
  alpha_u   <- mat$value$alpha_u
  Phi_g     <- mat$value$Phi_g
  Beta_g    <- mat$value$Beta_g
  alpha_g   <- mat$value$alpha_g
  #phi       <- solve(Phi_u) # for alpha and Beta
  M         <- length(mat$pattern$alpha_p)

  ## reference components updating

  # alpha
  for (j in which(.is_est(mat$pattern$alpha_p))){
    w_alpha_u    <- lapply(1:n_gps, function(x) {1/(w_g[[x]]*solve(Phi_u+Phi_g[[x]])[j,j])}) %>% unlist %>% sum
    alpha_u[j]   <- 
      w_alpha_u * (
        sum(unlist(lapply(1:n_gps, function(x) {w_g[[x]] * solve(Phi_u+Phi_g[[x]])[j,j] * (e_eta[[x]][j] - alpha_g[[x]][j] - (Beta_u[j,]+Beta_g[[x]][j,]) %*% e_eta[[x]])}))) +
        sum(unlist(lapply(1:n_gps, function(x) {w_g[[x]] * solve(Phi_u+Phi_g[[x]])[j,-j]  %*% (e_eta[[x]] - (alpha_u + alpha_g[[x]]) - (Beta_u+Beta_g[[x]]) %*% e_eta[[x]])[-j]})))
        )
        }

  # Beta
  
  for (i in which(.is_est(mat$pattern$Beta_p))){
    k        <- ceiling(i/M)
    j        <- i-(k-1)*M
    w_beta_u <- lapply(1:n_gps, function(x) {1/(w_g[[x]]*solve(Phi_u+Phi_g[[x]])[j,j]*C_etaeta[[x]][k,k])}) %>% unlist %>% sum
    bet_u    <- 
      w_beta_u * (
        sum(unlist(lapply(1:n_gps, function(x) {w_g[[x]] * solve(Phi_u+Phi_g[[x]])[j,j] *
            (C_etaeta[[x]][j,k] - e_eta[[x]][k] * (alpha_u+alpha_g[[x]])[j] - Beta_u[j,-k] %*% C_etaeta[[x]][-k,k] - Beta_g[[x]][j,] %*% C_etaeta[[x]][,k])})))+
        sum(unlist(lapply(1:n_gps, function(x) {w_g[[x]] * solve(Phi_u+Phi_g[[x]])[j,-j] %*% 
            (C_etaeta[[x]][-j,k] - e_eta[[x]][k] * (alpha_u+alpha_g[[x]])[-j] - (Beta_u[-j,] + Beta_g[[x]][-j,]) %*% C_etaeta[[x]][,k])})))
      )
    cth<-is.na(mat$pattern$Beta_p[i])
    Beta_u[j,k]<-penalty(theta=bet_u,gamma=0.025,cth=cth,w=w_beta_u,delta=2.5,type=type)
  }
  
  # Phi
  
  C_zetazeta<- lapply(1:n_gps, function(x) {C_etaeta[[x]] - e_eta[[x]]%*%(alpha_u+alpha_g[[x]]) -  C_etaeta[[x]] %*% t(Beta_u+Beta_g[[x]]) - (alpha_u+alpha_g[[x]]) %*% t(e_eta[[x]]) +
    (alpha_u+alpha_g[[x]]) %*% t(alpha_u+alpha_g[[x]]) + (alpha_u+alpha_g[[x]]) %*% t(e_eta[[x]]) %*% t(Beta_u+Beta_g[[x]]) - (Beta_u+Beta_g[[x]]) %*% t(C_etaeta[[x]]) + 
    (Beta_u+Beta_g[[x]]) %*% e_eta[[x]] %*% t(alpha_u+alpha_g[[x]]) + (Beta_u+Beta_g[[x]]) %*% C_etaeta[[x]] %*% t(Beta_u+Beta_g[[x]])})
  
  for (i in which(.is_est(mat$pattern$Phi_p))){
    k      <- ceiling(i/M)
    j      <- i-(k-1)*M
    if (j<=k) {} else {
      lk<-k
      C_zetatildazeta     <- lapply(1:n_gps, function(x) {solve(Phi_u+Phi_g[[x]])[-j,-j]%*%C_zetazeta[[x]][-j,]})
      C_zetatildazetatilda<- lapply(1:n_gps, function(x) {solve(Phi_u+Phi_g[[x]])[-j,-j]%*%C_zetazeta[[x]][-j,-j]%*%solve(Phi_u+Phi_g[[x]])[-j,-j]})
      var_phi   <- lapply(1:n_gps, function(x) sapply(c(1:M), varphi, Phi_u+Phi_g[[x]]))
      w_phi_u   <- lapply(1:n_gps, function(x) 1/((w_g[[x]]/var_phi[[x]][j])*C_zetatildazetatilda[[x]][lk,lk])) %>% unlist %>% sum
      Phi_u[j,k]<- 
        w_phi_u * (
          lapply(1:n_gps, function(x) {(w_g[[x]] / var_phi[[x]][j]) * (C_zetatildazeta[[x]][lk,j] - Phi_u[j,-c(j,k)] %*% matrix(C_zetatildazetatilda[[x]][-lk,lk]) - 
                                                 Phi_g[[x]][j,-j] %*% C_zetatildazetatilda[[x]][,lk])}) %>% unlist %>% sum )
    }
  }
  Phi_u[upper.tri(Phi_u)]<-t(Phi_u)[upper.tri(Phi_u)]

  # increment components updating
  for (x in 1:n_gps) {
  #alpha
    if (x==1){} else {
      for (j in which(.is_est(mat$pattern$alpha_p))){
        w_alpha_g         <- 1/(w_g[[x]]*solve(Phi_g[[x]])[j,j])
        alpha_g[[x]][j]   <- w_alpha_g * ( w_g[[x]] * solve(Phi_u+Phi_g[[x]])[j,j] * (e_eta[[x]][j] - alpha_u[j] - (Beta_u+Beta_g[[x]])[j,] %*% e_eta[[x]]) +
                                       w_g[[x]] * ((Phi_u+Phi_g[[x]])[j,-j] %*% (e_eta[[x]] - (alpha_u + alpha_g[[x]]) - (Beta_u+Beta_g[[x]]) %*% e_eta[[x]])[-j]))
      }
  
      for (i in which(.is_est(mat$pattern$Beta_p))){
        k        <- ceiling(i/M)
        j        <- i-(k-1)*M
        w_beta_g <- 1/(w_g[[x]]*solve(Phi_u+Phi_g[[x]])[j,j]*C_etaeta[[x]][k,k])
        bet_u    <- w_beta_g * (
            w_g[[x]] * solve(Phi_u+Phi_g[[x]])[j,j] * (C_etaeta[[x]][j,k] - e_eta[[x]][k] * (alpha_u+alpha_g[[x]])[j] - Beta_u[j,] %*% C_etaeta[[x]][,k] - Beta_g[[x]][j,-k] %*% C_etaeta[[x]][-k,k])+
            w_g[[x]] * solve(Phi_u+Phi_g[[x]])[j,-j] %*% (C_etaeta[[x]][-j,k] - e_eta[[x]][k] * (alpha_u+alpha_g[[x]])[-j] - (Beta_u[-j,] + Beta_g[[x]][-j,]) %*% C_etaeta[[x]][,k]))
        cth<-is.na(mat$pattern$Beta_p[i])
        Beta_g[[x]][j,k]<-penalty(theta=bet_u,gamma=0.075,cth=cth,w=w_beta_g,delta=2.5,type=type)
      }
  
      
      C_zetazeta<- lapply(1:n_gps, function(x) {C_etaeta[[x]] - e_eta[[x]]%*%(alpha_u+alpha_g[[x]]) -  C_etaeta[[x]] %*% t(Beta_u+Beta_g[[x]]) - (alpha_u+alpha_g[[x]]) %*% t(e_eta[[x]]) +
          (alpha_u+alpha_g[[x]]) %*% t(alpha_u+alpha_g[[x]]) + (alpha_u+alpha_g[[x]]) %*% t(e_eta[[x]]) %*% t(Beta_u+Beta_g[[x]]) - (Beta_u+Beta_g[[x]]) %*% t(C_etaeta[[x]]) + 
          (Beta_u+Beta_g[[x]]) %*% e_eta[[x]] %*% t(alpha_u+alpha_g[[x]]) + (Beta_u+Beta_g[[x]]) %*% C_etaeta[[x]] %*% t(Beta_u+Beta_g[[x]])})
  
      
      for (i in which(.is_est(mat$pattern$Phi_p))){
        k      <- ceiling(i/M)
        j      <- i-(k-1)*M
        if (j<=k) {} else {
          lk<-k
          C_zetatildazeta     <- solve(Phi_u+Phi_g[[x]])[-j,-j]%*%C_zetazeta[[x]][-j,]
          C_zetatildazetatilda<- solve(Phi_u+Phi_g[[x]])[-j,-j]%*%C_zetazeta[[x]][-j,-j]%*%solve(Phi_u+Phi_g[[x]])[-j,-j]
          var_phi   <- varphi(j,(Phi_u+Phi_g[[x]]))
          w_phi_g   <- 1/(w_g[[x]]*var_phi*C_zetatildazetatilda[lk,lk]) %>% unlist
          Phi_g[[x]][j,k]<- w_phi_g * (w_g[[x]] / var_phi) * (C_zetatildazeta[lk,j] - Phi_u[j,-j] %*% matrix(C_zetatildazetatilda[,lk]) - Phi_g[[x]][j,-c(j,k)] %*% matrix(C_zetatildazetatilda[-lk,lk]))
        }
       }
      Phi_g[[x]][upper.tri(Phi_g[[x]])]<-t(Phi_g[[x]])[upper.tri(Phi_g[[x]])]
    
    }
      
    
      C_zetazeta<- lapply(1:n_gps, function(x) {C_etaeta[[x]] - e_eta[[x]]%*%(alpha_u+alpha_g[[x]]) -  C_etaeta[[x]] %*% t(Beta_u+Beta_g[[x]]) - (alpha_u+alpha_g[[x]]) %*% t(e_eta[[x]]) +
          (alpha_u+alpha_g[[x]]) %*% t(alpha_u+alpha_g[[x]]) + (alpha_u+alpha_g[[x]]) %*% t(e_eta[[x]]) %*% t(Beta_u+Beta_g[[x]]) - (Beta_u+Beta_g[[x]]) %*% t(C_etaeta[[x]]) + 
          (Beta_u+Beta_g[[x]]) %*% e_eta[[x]] %*% t(alpha_u+alpha_g[[x]]) + (Beta_u+Beta_g[[x]]) %*% C_etaeta[[x]] %*% t(Beta_u+Beta_g[[x]])})
      
    
      for (j in 1:M){
        C_zetatildazeta      <- solve(Phi_u+Phi_g[[x]])[-j,-j]%*%C_zetazeta[[x]][-j,]
        C_zetatildazetatilda <- solve(Phi_u+Phi_g[[x]])[-j,-j]%*%C_zetazeta[[x]][-j,-j]%*%solve(Phi_u+Phi_g[[x]])[-j,-j]
        Phi_g[[x]][j,j] <- C_zetazeta[[x]][j,j] - 2 * (Phi_u[j,-j]+Phi_g[[x]][j,-j]) %*% C_zetatildazeta[,j] + (Phi_u[j,-j]+Phi_g[[x]][j,-j]) %*% C_zetatildazetatilda %*% (Phi_u[-j,j]+Phi_g[[x]][-j,j]) +
          (Phi_u[j,-j]+Phi_g[[x]][j,-j]) %*% solve(Phi_u[-j,-j]+Phi_g[[x]][-j,-j]) %*% (Phi_u[-j,j]+Phi_g[[x]][-j,j])
      }
  }
 
  return(list(alpha_u=alpha_u,
              Beta_u=Beta_u,
              Phi_u=Phi_u,
              alpha_g=alpha_g,
              Beta_g=Beta_g,
              Phi_g=Phi_g))
} 

dml_cal   <- function(Sigma=Sigma,e_v=e_v,Sigma_vv=subset(ini$Sigma_etaeta,G_obs,G_obs),mu_v=subset(ini$mu_eta,G_obs)){
  Sigma_vv_iv <- solve(Sigma_vv)
  dml <- -log(det(Sigma_vv_iv %*% Sigma)) + sum(diag(Sigma_vv_iv%*%Sigma)) - dim(Sigma_vv)[1] + t(e_v-mu_v) %*% Sigma_vv_iv %*% (e_v-mu_v)
  return(dml)
  }
