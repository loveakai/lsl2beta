import <-
  function(raw_obs,var_subset,var_group,obs_subset,obs_weight,raw_cov,raw_mean,obs_size) {
    if (missing(raw_obs)) {
      output<-list(raw_obs = raw_obs, raw_mean = raw_mean)
      if (exists(obs_size)) {attr(output,"obs_size")<-obs_size}
    }  else {
      if (!is.data.frame(raw_obs)) {
        as.data.frame(raw_obs)
        }
      if (missing(obs_subset)) {
        obs_subset <- 1:nrow(raw_obs)
      }
      if (missing(var_group)) {
        if (missing(var_subset)) {
          var_subset <- 1:ncol(raw_obs)
        }
        raw_obs %<>% .[obs_subset, ] 
        if (is.null(colnames(raw_obs))) colnames(raw_obs)<-paste0("v",1:ncol(raw_obs))
        raw_obs %<>% cbind(group=1)
        var_group <-ncol(raw_obs)
      } else {
        if (is.character(var_group)) {
          var_group <- which(colnames(raw_obs)%in%(var_group))
        }
        if (missing(var_subset)) {
          var_subset <- (1:ncol(raw_obs)) %>% .[!. %in% var_group]
        }
        raw_o   <- raw_obs[obs_subset,var_subset]
        if (is.null(colnames(raw_o))) colnames(raw_o)<-paste0("v",1:ncol(raw_o))
        raw_obs <- cbind(raw_o,raw_obs[,var_group])
      }
      output<-list(
        raw_obs = raw_obs,
        raw_cov = split(raw_obs, raw_obs[,ncol(raw_obs)]) %>% lapply(function(x) {
          cov(x[, -ncol(x)])
        }),
        raw_mean = split(raw_obs, raw_obs[,ncol(raw_obs)]) %>% lapply(function(x) {
          apply(x[, -ncol(x)], 2, mean)
        })
      )
      attr(output,"obs_size")  <- plyr::count(raw_obs[,ncol(raw_obs)]) %>% .[,2]
      
    }
    attr(output,"var_group") <- output$raw_cov %>% length
    if (!is.null(colnames(output$raw_cov[[1]]))) {
      attr(output,"v_label")<-colnames(output$raw_cov[[1]])
      }
    return(output)
  }


betagen   <- function(lambda){
  
  beta_p                          <- matrix(0, ncol = n_eta, nrow = n_eta)
  beta_p[(1:n_v),(n_v+1):n_eta]   <- lambda
  #beta_p[(n_v+1):n_eta,(n_v+1):n_eta] <- diag(1,n_f,n_f)
  return(beta_p)
  
}

matgen    <- function(alpha_p,beta_p,phi_p,alpha_r,beta_r,phi_r,lambda,n_groups,scale=T){ #**diminfo could be modified

  #pattern matarices generation
  #only 3 matrices have pattern matrix, including alpha, beta, Phi
  
  n_v       <- nrow(lambda)
  n_f       <- ncol(lambda)
  n_eta     <- n_v + n_f
  nm <- c(v_label,f_label)
  
  if (missing(beta_p))  {
    if (missing(lambda)) { 
      stop("lambda matrix is not specified")
    } else { 
      beta_p                        <- betagen(lambda) 
    }
  }
  
  if (missing(alpha_p)) {   alpha_p <- c(rep(1,n_v),rep(0,n_f)) }
  
  if (missing(phi_p)) {
    phi_p                           <- matrix(diag(1,n_eta,n_eta), ncol = n_eta, nrow = n_eta)
    phi_p[(n_v+1):n_eta,(n_v+1):n_eta]  <- 1
  }
  

  #matrices generation
  ## reference component
  if (missing(alpha_r)) {
    alpha_r <- c(rep(0.01,n_v),rep(0,n_f))
  } else {
    alpha_r<-alpha_r }
  
  if (missing(beta_r))  {
    beta_r  <- 0.8 *.is_one(beta_p)
  } else {
    beta_r <-beta_r} 
  
  if (missing(phi_r))   {
    #phi_r   <- matrix(diag(c(rep(0.1,n_v,n_v),rep(1,n_f,n_f))), ncol = n_eta, nrow = n_eta)
    phi_r <- diag(0,n_eta,n_eta)
    } else {
    phi_r <- phi_r
    }
  
  names(alpha_p)    <-colnames(beta_p)   <- rownames(beta_p)  <- colnames(phi_p)  <-rownames(phi_p)      <- nm
  rownames(phi_r)   <- colnames(phi_r)   <- rownames(beta_r)  <- colnames(beta_r) <-names(alpha_r)      <- nm

  if (scale) {beta_p[apply(beta_p[(1:n_v),(n_v+1):n_eta] == 1, 2, function(x) min(which(x))) %>% cbind((n_v+1):n_eta)] <- 0}　
  
  ## increment component
  alpha_i<-rep(list(0*alpha_r),n_groups)
  beta_i <-rep(list(0*beta_r),n_groups)
  
  g<-0.1*diag(ncol(phi_r))
  g[(n_v+1):n_eta,(n_v+1):n_eta]<-(diag(n_f))
  phi_i  <-rep(list(g %>%`colnames<-`(nm) %>% `rownames<-`(nm)),n_groups) 
  
  return(list(pattern=list(alpha_p=alpha_p,beta_p=beta_p,phi_p=phi_p),value=list(alpha_r=alpha_r,beta_r=beta_r,phi_r=phi_r,alpha_i=alpha_i,beta_i=beta_i,phi_i=phi_i)))
  
}

threshold <- function(theta,gma){
   sign(theta)*max(abs(theta)-gma,0) 
}

varphi    <- function(x,Phi) {diag(Phi)[x]-Phi[x,-x]%*%solve(Phi[-x,-x])%*%Phi[-x,x]}

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
  mu_v      <- lapply(1:n_groups, function(i_groups) subset(ini$mu_eta[[i_groups]],ini$G_eta))
  sigma_vf  <- lapply(1:n_groups, function(i_groups) subset(ini$sigma_eta[[i_groups]],ini$G_eta))
  sigma_fv  <- lapply(1:n_groups, function(i_groups) t(sigma_vf[[i_groups]]))
  sigma_v   <- lapply(1:n_groups, function(i_groups) subset(ini$sigma_eta[[i_groups]],ini$G_eta,ini$G_eta))
  J         <- lapply(1:n_groups, function(i_groups) ini$mu_eta[[i_groups]] - sigma_fv[[i_groups]] %*% solve(sigma_v[[i_groups]]) %*% mu_v[[i_groups]])
  K         <- lapply(1:n_groups, function(i_groups) sigma_fv[[i_groups]] %*% solve(sigma_v[[i_groups]]))
  c_v       <- lapply(1:n_groups, function(i_groups) ini$sigma[[i_groups]])
  e_eta     <- lapply(1:n_groups, function(i_groups) J[[i_groups]]+K[[i_groups]]%*%ini$e_v[[i_groups]])
  c_eta  <- lapply(1:n_groups, function(i_groups) {ini$sigma_eta[[i_groups]] - sigma_fv[[i_groups]] %*% solve(sigma_v[[i_groups]]) %*% sigma_vf[[i_groups]] +
    J[[i_groups]] %*% t(J[[i_groups]]) + J[[i_groups]] %*% t(ini$e_v[[i_groups]]) %*% t(K[[i_groups]]) + K[[i_groups]] %*% ini$e_v[[i_groups]] %*% t(J[[i_groups]]) +
    K[[i_groups]] %*% c_v[[i_groups]] %*% t(K[[i_groups]])})
  
  return(list(e_eta=e_eta, c_eta=c_eta))
}

cmstep    <- function(w_g=w_g,JK=JK,JLK=JLK,mat=ini$mat,e_step=e_step,type=type,gamma=gamma,delta=delta){
  
  e_eta     <- e_step$e_eta
  c_eta     <- e_step$c_eta
  phi_r     <- mat$value$phi_r
  beta_r    <- mat$value$beta_r
  alpha_r   <- mat$value$alpha_r
  phi_i     <- mat$value$phi_i
  beta_i    <- mat$value$beta_i
  alpha_i   <- mat$value$alpha_i
  n_eta     <- length(mat$pattern$alpha_p)

  ## reference components updating

  # alpha
  for (j in which(.is_est(mat$pattern$alpha_p))){
    w_alpha_r    <- 1/(sapply(1:n_groups, function(i_groups) {w_g[[i_groups]]*solve(phi_r+phi_i[[i_groups]])[j,j]}) %>% sum)
    alpha_r[j]   <- 
      w_alpha_r * (
        sum(unlist(lapply(1:n_groups, function(i_groups) {
          w_g[[i_groups]] * solve(phi_r+phi_i[[i_groups]])[j,j] * (e_eta[[i_groups]][j] - alpha_i[[i_groups]][j] - (beta_r[j,]+beta_i[[i_groups]][j,]) %*% e_eta[[i_groups]])}))) +
        sum(unlist(lapply(1:n_groups, function(i_groups) {
          w_g[[i_groups]] * solve(phi_r+phi_i[[i_groups]])[j,-j]  %*% (e_eta[[i_groups]] - (alpha_r + alpha_i[[i_groups]]) - (beta_r+beta_i[[i_groups]]) %*% e_eta[[i_groups]])[-j]})))
        )
        }

  # Beta
  
  for (i in which(.is_est(mat$pattern$beta_p))){
    k        <- ceiling(i/n_eta)
    j        <- i-(k-1)*n_eta
    w_beta_r <- 1/(sapply(1:n_groups, function(i_groups) {w_g[[i_groups]]*solve(phi_r+phi_i[[i_groups]])[j,j]*c_eta[[i_groups]][k,k]}) %>% sum)
    beta    <- 
      w_beta_r * (
        sum(sapply(1:n_groups, function(i_groups) {w_g[[i_groups]] * solve(phi_r+phi_i[[i_groups]])[j,j] *
            (c_eta[[i_groups]][j,k] - e_eta[[i_groups]][k] * (alpha_r+alpha_i[[i_groups]])[j] - beta_r[j,-k] %*% c_eta[[i_groups]][-k,k] - beta_i[[i_groups]][j,] %*% c_eta[[i_groups]][,k])}))+
        sum(sapply(1:n_groups, function(i_groups) {w_g[[i_groups]] * solve(phi_r+phi_i[[i_groups]])[j,-j] %*% 
            (c_eta[[i_groups]][-j,k] - e_eta[[i_groups]][k] * (alpha_r+alpha_i[[i_groups]])[-j] - (beta_r[-j,] + beta_i[[i_groups]][-j,]) %*% c_eta[[i_groups]][,k])}))
      )
    cth<-is.na(mat$pattern$beta_p[i])
    beta_r[j,k]<-penalty(theta=beta,gamma=gamma,cth=cth,w=w_beta_r,delta=delta,type=type)
  }
  
  # Phi
  
  c_zeta<- lapply(1:n_groups, function(i_groups) {c_eta[[i_groups]] - e_eta[[i_groups]]%*%(alpha_r+alpha_i[[i_groups]]) -  c_eta[[i_groups]] %*% t(beta_r+beta_i[[i_groups]]) - (alpha_r+alpha_i[[i_groups]]) %*% t(e_eta[[i_groups]]) +
    (alpha_r+alpha_i[[i_groups]]) %*% t(alpha_r+alpha_i[[i_groups]]) + (alpha_r+alpha_i[[i_groups]]) %*% t(e_eta[[i_groups]]) %*% t(beta_r+beta_i[[i_groups]]) - (beta_r+beta_i[[i_groups]]) %*% t(c_eta[[i_groups]]) + 
    (beta_r+beta_i[[i_groups]]) %*% e_eta[[i_groups]] %*% t(alpha_r+alpha_i[[i_groups]]) + (beta_r+beta_i[[i_groups]]) %*% c_eta[[i_groups]] %*% t(beta_r+beta_i[[i_groups]])})
  
  for (i in which(.is_est(mat$pattern$phi_p))){
    k      <- ceiling(i/n_eta)
    j      <- i-(k-1)*n_eta
    if (j<=k) {} else {
      lk<-k
      c_zetatzeta <- lapply(1:n_groups, function(i_groups) {solve(phi_r[-j,-j]+phi_i[[i_groups]][-j,-j])%*%c_zeta[[i_groups]][-j,]})
      c_zetat     <- lapply(1:n_groups, function(i_groups) {solve(phi_r[-j,-j]+phi_i[[i_groups]][-j,-j])%*%c_zeta[[i_groups]][-j,-j]%*%solve(phi_r[-j,-j]+phi_i[[i_groups]][-j,-j])})
      var_phi     <- lapply(1:n_groups, function(i_groups) sapply(c(1:n_eta), varphi, phi_r+phi_i[[i_groups]]))
      w_phi_r     <- 1/(sapply(1:n_groups, function(i_groups) {(w_g[[i_groups]]/var_phi[[i_groups]][j])*c_zetat[[i_groups]][lk,lk]}) %>% sum)
      phi<- 
        w_phi_r * (
          sapply(1:n_groups, function(i_groups) {(w_g[[i_groups]] / var_phi[[i_groups]][j]) * (c_zetatzeta[[i_groups]][lk,j] - phi_r[j,-c(j,k)] %*% matrix(c_zetat[[i_groups]][-lk,lk]) - 
                                                 phi_i[[i_groups]][j,-j] %*% c_zetat[[i_groups]][,lk])})  %>% sum )
      cth<-is.na(mat$pattern$phi_p[i])
      phi_r[j,k]<-penalty(theta=phi,gamma=gamma,cth=cth,w=w_phi_r,delta=delta,type=type)
    }
  }
  phi_r[upper.tri(phi_r)]<-t(phi_r)[upper.tri(phi_r)]

  ## increment components updating
  for (i_groups in 1:n_groups) {
  #alpha
    if (i_groups==1){} else {
      for (j in which(.is_est(mat$pattern$alpha_p))){
        w_alpha_i         <- 1/(w_g[[i_groups]]*solve(phi_i[[i_groups]])[j,j])
        alpha_i[[i_groups]][j]   <- w_alpha_i * ( w_g[[i_groups]] * solve(phi_r+phi_i[[i_groups]])[j,j] * (e_eta[[i_groups]][j] - alpha_r[j] - (beta_r+beta_i[[i_groups]])[j,] %*% e_eta[[i_groups]]) +
                                       w_g[[i_groups]] * ((phi_r+phi_i[[i_groups]])[j,-j] %*% (e_eta[[i_groups]] - (alpha_r + alpha_i[[i_groups]]) - (beta_r+beta_i[[i_groups]]) %*% e_eta[[i_groups]])[-j]))
      }
  
      for (i in which(.is_est(mat$pattern$beta_p))){
        k        <- ceiling(i/n_eta)
        j        <- i-(k-1)*n_eta
        w_beta_i <- 1/(w_g[[i_groups]]*solve(phi_r+phi_i[[i_groups]])[j,j]*c_eta[[i_groups]][k,k])
        beta    <- w_beta_i * (
            w_g[[i_groups]] * solve(phi_r+phi_i[[i_groups]])[j,j] * (c_eta[[i_groups]][j,k] - e_eta[[i_groups]][k] * (alpha_r+alpha_i[[i_groups]])[j] - beta_r[j,] %*% c_eta[[i_groups]][,k] - beta_i[[i_groups]][j,-k] %*% c_eta[[i_groups]][-k,k])+
            w_g[[i_groups]] * solve(phi_r+phi_i[[i_groups]])[j,-j] %*% (c_eta[[i_groups]][-j,k] - e_eta[[i_groups]][k] * (alpha_r+alpha_i[[i_groups]])[-j] - (beta_r[-j,] + beta_i[[i_groups]][-j,]) %*% c_eta[[i_groups]][,k]))
        cth<-is.na(mat$pattern$beta_p[i])
        beta_i[[i_groups]][j,k]<-penalty(theta=beta,gamma=gamma,cth=cth,w=w_beta_i,delta=delta,type=type)
      }
  
      
      c_zeta<- lapply(1:n_groups, function(i_groups) {c_eta[[i_groups]] - e_eta[[i_groups]]%*%(alpha_r+alpha_i[[i_groups]]) -  c_eta[[i_groups]] %*% t(beta_r+beta_i[[i_groups]]) - (alpha_r+alpha_i[[i_groups]]) %*% t(e_eta[[i_groups]]) +
          (alpha_r+alpha_i[[i_groups]]) %*% t(alpha_r+alpha_i[[i_groups]]) + (alpha_r+alpha_i[[i_groups]]) %*% t(e_eta[[i_groups]]) %*% t(beta_r+beta_i[[i_groups]]) - (beta_r+beta_i[[i_groups]]) %*% t(c_eta[[i_groups]]) + 
          (beta_r+beta_i[[i_groups]]) %*% e_eta[[i_groups]] %*% t(alpha_r+alpha_i[[i_groups]]) + (beta_r+beta_i[[i_groups]]) %*% c_eta[[i_groups]] %*% t(beta_r+beta_i[[i_groups]])})
  
      
      for (i in which(.is_est(mat$pattern$phi_p))){
        k      <- ceiling(i/n_eta)
        j      <- i-(k-1)*n_eta
        if (j<=k) {} else {
          lk<-k
          c_zetatzeta <- solve(phi_r[-j,-j]+phi_i[[i_groups]][-j,-j])%*%c_zeta[[i_groups]][-j,]
          c_zetat     <- solve(phi_r[-j,-j]+phi_i[[i_groups]][-j,-j])%*%c_zeta[[i_groups]][-j,-j]%*%solve(phi_r[-j,-j]+phi_i[[i_groups]][-j,-j])
          var_phi     <- varphi(j,(phi_r+phi_i[[i_groups]]))
          w_phi_i     <- 1/(w_g[[i_groups]]/var_phi*c_zetat[lk,lk])
          phi       <- w_phi_i * (w_g[[i_groups]] / var_phi) * (c_zetatzeta[lk,j] - phi_r[j,-j] %*% matrix(c_zetat[,lk]) - phi_i[[i_groups]][j,-c(j,k)] %*% matrix(c_zetat[-lk,lk]))
          cth<-is.na(mat$pattern$phi_p[i])
          phi_i[[i_groups]][j,k]<-penalty(theta=phi,gamma=gamma,cth=cth,w=w_phi_i,delta=delta,type=type)
        }
       }
      phi_i[[i_groups]][upper.tri(phi_i[[i_groups]])]<-t(phi_i[[i_groups]])[upper.tri(phi_i[[i_groups]])]
    
    }
      
    
      for (j in 1:n_eta){
        c_zetatzeta  <- solve(phi_r[-j,-j]+phi_i[[i_groups]][-j,-j])%*%c_zeta[[i_groups]][-j,]
        c_zetat      <- solve(phi_r[-j,-j]+phi_i[[i_groups]][-j,-j])%*%c_zeta[[i_groups]][-j,-j]%*%solve(phi_r[-j,-j]+phi_i[[i_groups]][-j,-j])
        phi          <- c_zeta[[i_groups]][j,j] - 2 * (phi_r[j,-j]+phi_i[[i_groups]][j,-j]) %*% c_zetatzeta[,j] + (phi_r[j,-j]+phi_i[[i_groups]][j,-j]) %*% c_zetat %*% (phi_r[-j,j]+phi_i[[i_groups]][-j,j]) +
          (phi_r[j,-j]+phi_i[[i_groups]][j,-j]) %*% solve(phi_r[-j,-j]+phi_i[[i_groups]][-j,-j]) %*% (phi_r[-j,j]+phi_i[[i_groups]][-j,j])
        cth<-is.na(mat$pattern$phi_p[i])
        phi_i[[i_groups]][j,j]<-penalty(theta=phi,gamma=gamma,cth=cth,w=1,delta=delta,type=type)
      }
    
  }
 
  return(list(alpha_r=alpha_r,
              beta_r=beta_r,
              phi_r=phi_r,
              alpha_i=alpha_i,
              beta_i=beta_i,
              phi_i=phi_i))
} 

ecm       <- function(mat=mat,ide=ide,G_eta=G_eta,maxit=500,cri=10^(-5),penalize=pl){
  
  if (missing(penalize)) {
    message("penalize information is not specified, using default")
    type<-"L1"
    gamma<-0.025
    delta<-2.5
  } else{
    if (exists("type",penalize))  {type<-penalize$type}     else {message("type is not specified,using default:'L1'")  ;type <-"L1"}
    if (exists("gamma",penalize)) {gamma<-penalize$gamma}   else {message("gamma is not specified,using default:0.025");gamma<-0.025} 
    if (exists("delta",penalize)) {delta<-penalize$delta}   else {message("delta is not specified,using default:2.5")  ;delta<-2.5}  
  }
  

  alpha_p   <- mat$pattern$alpha_p
  beta_p    <- mat$pattern$beta_p
  phi_p     <- mat$pattern$phi_p
  alpha_r   <- mat$value$alpha_r
  beta_r    <- mat$value$beta_r
  phi_r     <- mat$value$phi_r
  alpha_i   <- mat$value$alpha_i
  beta_i    <- mat$value$beta_i
  phi_i     <- mat$value$phi_i
  
  #initialization
  
  IBinv     <- lapply(1:n_groups, function(i_groups) solve(ide-(beta_r+beta_i[[i_groups]])))
  mu_eta    <- lapply(1:n_groups, function(i_groups) IBinv[[i_groups]]%*%(alpha_r+alpha_i[[i_groups]]))
  sigma_eta <- lapply(1:n_groups, function(i_groups) IBinv[[i_groups]]%*%(phi_r+phi_i[[i_groups]])%*%t(IBinv[[i_groups]]))
  
  w_g       <- rep(list(1/3),n_groups)
  JK        <- expand.grid(1:n_eta,1:n_eta)[2:1]
  JLK       <- expand.grid(1:(n_eta-1),1:n_eta)[2:1]
  
  ini       <- list(IBinv=IBinv,mu_eta=mu_eta,sigma_eta=sigma_eta,sigma=sigma,G_eta=G_eta,e_v=e_v,mat=mat)

  
  for (it in 1:maxit){
    cat(it, "...")
    e_step    <- estep(ini)
    cm_step   <- cmstep(w_g=w_g,JK=JK,JLK=JLK,mat=ini$mat,e_step=e_step,type=type,gamma=gamma,delta=delta)
    ini$mat$value$phi_r   <- phi_r     <- cm_step$phi_r
    ini$mat$value$beta_r  <- beta_r    <- cm_step$beta_r
    ini$mat$value$alpha_r <- alpha_r   <- cm_step$alpha_r
    ini$mat$value$phi_i   <- phi_i     <- cm_step$phi_i
    ini$mat$value$beta_i  <- beta_i    <- cm_step$beta_i
    ini$mat$value$alpha_i <- alpha_i   <- cm_step$alpha_i
    ini$IBinv             <- lapply(1:n_groups, function(i_groups) solve(ide-(beta_r+beta_i[[i_groups]])))
    ini$mu_eta            <- lapply(1:n_groups, function(i_groups) ini$IBinv[[i_groups]]%*%(alpha_r+alpha_i[[i_groups]]))
    ini$sigma_eta         <- lapply(1:n_groups, function(i_groups) ini$IBinv[[i_groups]]%*%(phi_r+phi_i[[i_groups]])%*%t(ini$IBinv[[i_groups]]))
  }
  
  #theta     <- c(alpha_r[.is_est(ini$mat$pattern$alpha_p)],beta_r[.is_est(ini$mat$pattern$beta_p)],phi_r[.is_est(ini$mat$pattern$phi_p)])
  #length(theta)
  
  dml<-dml_cal(sigma=sigma,e_v=e_v,ini=ini,G_eta=G_eta,n_groups=n_groups,w_g=w_g)
  
  return(list(theta=ini,dml=dml,gamma=gamma,delta=delta,type=type,iteration=it))
  
}

dml_cal   <- function(sigma=sigma,e_v=e_v,ini=ini,G_eta=G_eta,n_groups=n_groups,w_g=w_g){
  mu_v        <- lapply(1:n_groups, function(i_groups) subset(ini$mu_eta[[i_groups]],ini$G_eta))
  sigma_v     <- lapply(1:n_groups, function(i_groups) subset(ini$sigma_eta[[i_groups]],ini$G_eta,ini$G_eta))
  sigma_v_iv  <- lapply(1:n_groups, function(i_groups) solve(sigma_v[[i_groups]]))
  dml <-
    (sapply(1:n_groups, function(i_groups) {
      w_g[[i_groups]] * (sum(diag(sigma[[i_groups]] %*% sigma_v_iv[[i_groups]])) -
                           log(det(sigma[[i_groups]] %*% sigma_v_iv[[i_groups]]))  -
                           dim(sigma_v[[i_groups]])[1])
    })  +
      sapply(1:n_groups, function(i_groups) {
        w_g[[i_groups]] * (t(e_v[[i_groups]] - mu_v[[i_groups]]) %*% sigma_v_iv[[i_groups]] %*% (e_v[[i_groups]] -
                                                                                                    mu_v[[i_groups]]))
      })) %>% sum
  return(dml)
}


  
getpar    <-function(pattern,value,v_label,f_label,mat_label){
 mapply(function(val,p,v,nm) {
    val[p|v] %>% `names<-`(nm[p|v])},
    val=value,
    p=lapply(pattern,.is_est),
    v=lapply(value,function(x) {x!=0}),
    nm=list(c(v_label,f_label),mat_label,mat_label)) %>% `names<-`(c("alpha","beta","gamma"))
}





