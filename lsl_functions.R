.matgen    <- function(pattern,value,n_groups,scale,labels,ref_group,data=data){ #**diminfo could be modified
  
  #pattern matarices generation
  #only 3 matrices have pattern matrix, including alpha, beta, Phi
  
  n_v       <- nrow(pattern$beta_vf)
  n_f       <- ncol(pattern$beta_vf)
  n_eta     <- n_v + n_f
  nm        <- c(labels$v_label,labels$f_label)
  nm_g      <- attributes(data)$g_label[c(ref_group,(1:n_groups)[-ref_group])]
  
  beta_p                                 <- matrix(0, ncol = n_eta, nrow = n_eta)
  beta_p[(1:n_v),(n_v+1):n_eta]          <- pattern$beta_vf
  if(exists("beta_ff",pattern)){
    beta_p[(n_v+1):n_eta,(n_v+1):n_eta]  <- pattern$beta_ff
  }
  if(exists("beta_vv",pattern)){
    beta_p[(1:n_v),(1:n_v)]              <- pattern$beta_vv
  }
  if(exists("beta_fv",pattern)){
    beta_p[(n_v+1):n_eta,(1:n_v)]        <- pattern$beta_fv
  }
  
  alpha_p                                <- c(rep(1,n_v),rep(0,n_f))
  if(exists("alpha_v",pattern)){
    alpha_p[1:n_v]                       <- pattern$alpha_v
  }
  if(exists("alpha_f",pattern)){
    alpha_p[(n_v+1):n_eta]               <- pattern$alpha_f
  }
  
  phi_p                                  <- matrix(diag(1,n_eta,n_eta), ncol = n_eta, nrow = n_eta)
  phi_p[(n_v+1):n_eta,(n_v+1):n_eta]     <- 1
  if(exists("phi_ff",pattern)){
    phi_p[(n_v+1):n_eta,(n_v+1):n_eta]   <- pattern$phi_ff
  }
  if(exists("phi_vv",pattern)){
    phi_p[(1:n_v),(1:n_v)]               <- pattern$phi_vv
  }
  
  
  #matrices generation
  ## reference component
  if(missing(value)){
    alpha_r <- c(data$raw_mean[[nm_g[[1]]]],rep(0,n_f))
    beta_r  <- 1 *.is_one(beta_p)
    phi_r   <- diag(0,n_eta,n_eta)
  } else {
    
    alpha_r                                <- c(data$raw_mean[[nm_g[[1]]]],rep(0,n_f))
    if(exists("alpha_v",value)){
      alpha_r[1:n_v]                       <- value$alpha_v
    }
    if(exists("alpha_f",value)){
      alpha_r[(n_v+1):n_eta]               <- value$alpha_f
    }
    
    beta_r  <- 1 *.is_one(beta_p)
    if(exists("beta_vf",value)){
      beta_r[(1:n_v),(n_v+1):n_eta]        <- value$beta_vf
    }
    if(exists("beta_ff",value)){
      beta_r[(n_v+1):n_eta,(n_v+1):n_eta]  <- value$beta_ff
    }
    if(exists("beta_vv",value)){
      beta_r[(1:n_v),(1:n_v)]              <- value$beta_vv
    }
    if(exists("beta_fv",value)){
      beta_r[(n_v+1):n_eta,(1:n_v)]        <- value$beta_fv
    }
    
    phi_r                                  <- diag(0,n_eta,n_eta)
    if(exists("phi_ff",value)){
      phi_r[(n_v+1):n_eta,(n_v+1):n_eta]   <- value$phi_ff
    }
    if(exists("phi_vv",value)){
      phi_r[(1:n_v),(1:n_v)]               <- value$phi_vv
    }
  }
  
  names(alpha_p)    <-colnames(beta_p)   <- rownames(beta_p)  <- colnames(phi_p)  <-rownames(phi_p)      <- nm
  rownames(phi_r)   <-colnames(phi_r)    <- rownames(beta_r)  <- colnames(beta_r) <-names(alpha_r)       <- nm
  
  if (scale) {beta_p[apply(beta_p[(1:n_v),(n_v+1):n_eta] == 1, 2, function(x) min(which(x))) %>% cbind((n_v+1):n_eta)] <- 0}　
  
  ## increment component
  alpha_i<-lapply(data$raw_mean, function(x) {c(x,rep(0,n_f))-alpha_r})[c(ref_group,(1:n_groups)[-ref_group])]
  #alpha_i<-rep(list(rep(0,n_eta)),n_groups)
  beta_i <-rep(list(0*beta_r),n_groups)
  
  g<-0.1*diag(ncol(phi_r))
  g[(n_v+1):n_eta,(n_v+1):n_eta]<-(diag(n_f))
  phi_i  <-rep(list(g %>%`colnames<-`(nm) %>% `rownames<-`(nm)),n_groups) 
  
  names(alpha_i)    <-names(beta_i)      <-names(phi_i)       <-nm_g
  
  output<-list(pattern=list(alpha_p=alpha_p,beta_p=beta_p,phi_p=phi_p),value=list(alpha_r=alpha_r,beta_r=beta_r,phi_r=phi_r,alpha_i=alpha_i,beta_i=beta_i,phi_i=phi_i))
  attr(output,"n_v") <- n_v
  attr(output,"n_f") <- n_f
  attr(output,"n_eta") <- n_eta
  return(output)
  
}

.getpar    <- function(pattern,value,v_label,f_label,mat_label,group){
  pars<-mapply(function(ma,val,pat,p,v,co,ro,nm) {
    rbind(nm[p|v],ma,group,ro[p|v],co[p|v],pat[p|v],val[p|v]) %>% 
      `rownames<-`(c("name","matrix","group","row","col","type","value"))},
    val=lapply(value,as.matrix),
    pat=lapply(pattern,as.matrix),
    p=lapply(pattern,function(x) {as.matrix(x) %>% .is_est}),
    v=lapply(value,function(x) {as.matrix(x)!=0}),
    co=lapply(pattern,function(x) {as.matrix(x) %>% col}),
    ro=lapply(pattern,function(x) {as.matrix(x) %>% row}),
    nm=list(c(v_label,f_label),mat_label,mat_label) %>% lapply(as.matrix),
    ma=list("alpha","beta","phi")) %>% 
    `names<-`(c("alpha","beta","phi")) %>% lapply(.,t) %>% do.call(rbind,.)
   within(as.data.frame(pars),{
     col    <-as.numeric(levels(col)[col])
     row    <-as.numeric(levels(row)[row])
     value  <-as.numeric(levels(value)[value])
     type   <-as.numeric(levels(type)[type])
   }) %>% return
  }

.varphi    <- function(x,Phi) {diag(Phi)[x]-Phi[x,-x]%*%solve(Phi[-x,-x])%*%Phi[-x,x]}

.penalty   <- function(theta,gamma,cth,w,delta,type){
  if (type == "l1") {
    theta <- sign(theta)*max(abs(theta) - gamma * cth * w,0)
  } else if (type == "scad") {
    if (abs(theta) <= gamma * (1 + cth * w)) {
      theta <- sign(theta)*max(abs(theta) - gamma * cth * w,0)
    } else if (gamma * (1 + cth * w) < abs(theta) & abs(theta) <= gamma * delta) {
      theta <-  sign(theta)*max(abs(theta) - ((cth * w * gamma * delta) / (delta - 1)),0) / (1 - ((cth * w) / (delta - 1)))
    } else { }
  } else if (type == "mcp") {
    if (abs(theta) <= gamma * delta) {
      theta <- sign(theta)*max(abs(theta) - gamma * cth * w,0) / (1 - ((cth * w) / delta))
    } else { }
  }
  return(theta)
}

.estep     <- function(ini,data=data){
  n_groups  <- attributes(data)$n_groups
  mu_eta    <- ini$mu_eta
  mu_v      <- lapply(1:n_groups, function(i_groups) subset(ini$mu_eta[[i_groups]],ini$G_eta))
  sigma_veta<- lapply(1:n_groups, function(i_groups) subset(ini$sigma_eta[[i_groups]],ini$G_eta))
  sigma_etav<- lapply(1:n_groups, function(i_groups) t(sigma_veta[[i_groups]]))
  sigma_v   <- lapply(1:n_groups, function(i_groups) subset(ini$sigma_eta[[i_groups]],ini$G_eta,ini$G_eta))
  J         <- lapply(1:n_groups, function(i_groups) mu_eta[[i_groups]] - sigma_etav[[i_groups]] %*% solve(sigma_v[[i_groups]]) %*% mu_v[[i_groups]])
  K         <- lapply(1:n_groups, function(i_groups) sigma_etav[[i_groups]] %*% solve(sigma_v[[i_groups]]))
  c_v       <- ini$c_v
  e_eta     <- lapply(1:n_groups, function(i_groups) J[[i_groups]]+K[[i_groups]]%*%ini$e_v[[i_groups]])
  c_eta     <- lapply(1:n_groups, function(i_groups) {ini$sigma_eta[[i_groups]] - sigma_etav[[i_groups]] %*% solve(sigma_v[[i_groups]]) %*% sigma_veta[[i_groups]] +
      J[[i_groups]] %*% t(J[[i_groups]]) + J[[i_groups]] %*% t(ini$e_v[[i_groups]]) %*% t(K[[i_groups]]) + K[[i_groups]] %*% ini$e_v[[i_groups]] %*% t(J[[i_groups]]) +
      K[[i_groups]] %*% c_v[[i_groups]] %*% t(K[[i_groups]])})
  return(list(e_eta=e_eta, c_eta=c_eta))
}


.cmstep    <- function(w_g=w_g,JK=JK,JLK=JLK,mat=ini$mat,e_step=e_step,type=type,gamma=gamma,delta=delta,data=data){
  
  e_eta     <- e_step$e_eta
  c_eta     <- e_step$c_eta
  phi_r     <- mat$value$phi_r
  beta_r    <- mat$value$beta_r
  alpha_r   <- mat$value$alpha_r
  phi_i     <- mat$value$phi_i
  beta_i    <- mat$value$beta_i
  alpha_i   <- mat$value$alpha_i
  n_eta     <- length(mat$pattern$alpha_p)
  n_groups  <- attributes(data)$n_groups

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
    beta_r[j,k]<-.penalty(theta=beta,gamma=gamma,cth=cth,w=w_beta_r,delta=delta,type=type)
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
      var_phi     <- lapply(1:n_groups, function(i_groups) sapply(c(1:n_eta), .varphi, phi_r+phi_i[[i_groups]]))
      w_phi_r     <- 1/(sapply(1:n_groups, function(i_groups) {(w_g[[i_groups]]/var_phi[[i_groups]][j])*c_zetat[[i_groups]][lk,lk]}) %>% sum)
      phi<- 
        w_phi_r * (
          sapply(1:n_groups, function(i_groups) {(w_g[[i_groups]] / var_phi[[i_groups]][j]) * (c_zetatzeta[[i_groups]][lk,j] - phi_r[j,-c(j,k)] %*% matrix(c_zetat[[i_groups]][-lk,lk]) - 
                                                                                                 phi_i[[i_groups]][j,-j] %*% c_zetat[[i_groups]][,lk])})  %>% sum )
      cth<-is.na(mat$pattern$phi_p[i])
      phi_r[j,k]<-.penalty(theta=phi,gamma=gamma,cth=cth,w=w_phi_r,delta=delta,type=type)
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
        beta_i[[i_groups]][j,k]<-.penalty(theta=beta,gamma=gamma,cth=cth,w=w_beta_i,delta=delta,type=type)
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
          var_phi     <- .varphi(j,(phi_r+phi_i[[i_groups]]))
          w_phi_i     <- 1/(w_g[[i_groups]]/var_phi*c_zetat[lk,lk])
          phi       <- w_phi_i * (w_g[[i_groups]] / var_phi) * (c_zetatzeta[lk,j] - phi_r[j,-j] %*% matrix(c_zetat[,lk]) - phi_i[[i_groups]][j,-c(j,k)] %*% matrix(c_zetat[-lk,lk]))
          cth<-is.na(mat$pattern$phi_p[i])
          phi_i[[i_groups]][j,k]<-.penalty(theta=phi,gamma=gamma,cth=cth,w=w_phi_i,delta=delta,type=type)
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
      phi_i[[i_groups]][j,j]<-.penalty(theta=phi,gamma=gamma,cth=cth,w=1,delta=delta,type=type)

    }
    
  }
  
  return(list(alpha_r=alpha_r,
              beta_r=beta_r,
              phi_r=phi_r,
              alpha_i=alpha_i,
              beta_i=beta_i,
              phi_i=phi_i))
} 

.ecm       <- function(mat=mat,maxit,cri,penalize,model=model,data=data){

  alpha_p   <- mat$pattern$alpha_p 

  beta_p    <- mat$pattern$beta_p
  phi_p     <- mat$pattern$phi_p
  alpha_r   <- mat$value$alpha_r
  beta_r    <- mat$value$beta_r
  phi_r     <- mat$value$phi_r
  alpha_i   <- mat$value$alpha_i
  beta_i    <- mat$value$beta_i
  phi_i     <- mat$value$phi_i
  
  type      <- penalize$pl
  delta     <- penalize$delta
  gamma     <- penalize$gamma
  
  n_groups  <-attributes(data)$n_groups
  n_eta     <-attributes(mat)$n_eta
  n_v       <-attributes(mat)$n_v
  n_f       <-attributes(mat)$n_f
  
  v_label   <-attributes(model)$labels$v_label
  f_label   <-attributes(model)$labels$f_label
  eta_label <-c(v_label,f_label)
  mat_label <-attributes(model)$labels$mat_label
  
  
  ide       <- diag(1, ncol = n_eta, nrow = n_eta)  %>% `colnames<-`(eta_label) %>% `rownames<-`(eta_label) 
  G_eta     <- c(rep(T,n_v),rep(F,n_f)) %>%`names<-`(eta_label)
  c_v       <- lapply(1:n_groups, function(i_groups) { data$raw_cov[[i_groups]] + data$raw_mean[[i_groups]] %>% tcrossprod })
  e_v       <- lapply(1:n_groups, function(i_groups) { data$raw_mean[[i_groups]] })
  #initialization
  
  IBinv     <- lapply(1:n_groups, function(i_groups) solve(ide-(beta_r+beta_i[[i_groups]])))
  mu_eta    <- lapply(1:n_groups, function(i_groups) IBinv[[i_groups]]%*%(alpha_r+alpha_i[[i_groups]]))
  sigma_eta <- lapply(1:n_groups, function(i_groups) IBinv[[i_groups]]%*%(phi_r+phi_i[[i_groups]])%*%t(IBinv[[i_groups]]))
  
  w_g       <- attributes(data)$obs_size/sum(attributes(data)$obs_size)
  JK        <- expand.grid(1:n_eta,1:n_eta)[2:1]
  JLK       <- expand.grid(1:(n_eta-1),1:n_eta)[2:1]
  
  ini       <- list(IBinv=IBinv,mu_eta=mu_eta,sigma_eta=sigma_eta,c_v=c_v,G_eta=G_eta,e_v=e_v,mat=mat)
  ref      <- .getpar(pattern = ini$mat$pattern,value = list(ini$mat$value$alpha_r,ini$mat$value$beta_r,ini$mat$value$phi_r),v_label,f_label,mat_label,group="r")
  output      <- lapply((1:n_groups),function(x){
    .getpar(pattern = ini$mat$pattern, value = list(ini$mat$value$alpha_i[[x]],ini$mat$value$beta_i[[x]],ini$mat$value$phi_i[[x]]),v_label,f_label,mat_label,group=names(ini$mat$value$alpha_i[x]))
  } ) %>% do.call(rbind,.) %>% rbind(ref,.)
  output   <-output[!(output$matrix=="phi"&(output$row>output$col)),]
  inival   <-output$value
  
  for (it in 1:maxit){
    
    #cat("...",it)
    #mod10<-ifelse(mod(it,75)==0,75,mod(it,75))
    #print(noquote(pic[mod10]))
    
    e_step    <- .estep(ini,data=data)
    cm_step   <- .cmstep(w_g=w_g,JK=JK,JLK=JLK,mat=ini$mat,e_step=e_step,type=type,gamma=gamma,delta=delta,data=data)
    ini$mat$value$phi_r   <- phi_r     <- cm_step$phi_r
    ini$mat$value$beta_r  <- beta_r    <- cm_step$beta_r
    ini$mat$value$alpha_r <- alpha_r   <- cm_step$alpha_r
    ini$mat$value$phi_i   <- phi_i     <- cm_step$phi_i
    ini$mat$value$beta_i  <- beta_i    <- cm_step$beta_i
    ini$mat$value$alpha_i <- alpha_i   <- cm_step$alpha_i
    ini$IBinv             <- lapply(1:n_groups, function(i_groups) solve(ide-(beta_r+beta_i[[i_groups]])))
    ini$mu_eta            <- lapply(1:n_groups, function(i_groups) ini$IBinv[[i_groups]]%*%(alpha_r+alpha_i[[i_groups]]))
    ini$sigma_eta         <- lapply(1:n_groups, function(i_groups) ini$IBinv[[i_groups]]%*%(phi_r+phi_i[[i_groups]])%*%t(ini$IBinv[[i_groups]]))
    ref      <- .getpar(pattern = ini$mat$pattern,value = list(ini$mat$value$alpha_r,ini$mat$value$beta_r,ini$mat$value$phi_r),v_label,f_label,mat_label,group="r")
    inc      <- lapply((1:n_groups),function(x){
      .getpar(pattern = ini$mat$pattern, value = list(ini$mat$value$alpha_i[[x]],ini$mat$value$beta_i[[x]],ini$mat$value$phi_i[[x]]),v_label,f_label,mat_label,group=names(ini$mat$value$alpha_i[x]))
    } ) %>% do.call(rbind,.) %>% rbind(ref,.) 
    output<-inc
    output   <-output[!(output$matrix=="phi"&(output$row>output$col)),]
    newval   <-output$value
    if (sum((newval-inival)^2)<cri) {break} else {inival<-newval} 
  }
  
  output_par<-output[.is_one(output$type)|(output$type==0&output$value!=0)|(is.na(output$type)&output$value!=0),]
  q<-output_par[.is_est(output_par$type),]
  ref_group<-attributes(model)$ref_group
  
  #n_par counting
  n_ref<-nrow(q[(q$group=="r")&!((q$matrix=="phi")&(q$row==q$col)),])
  n_inc1<-nrow(q[(q$group!="r")&(q$group!=ref_group),])
  n_inc2<-nrow(q[(q$group==ref_group)&((q$matrix=="phi")&(q$col==q$row)),])
  n_par<-n_ref+n_inc1+n_inc2
  
  sigma_v<-lapply((1:n_groups),function(i_groups) subset(ini$sigma_eta[[i_groups]],G_eta,G_eta))
  mu_v   <-lapply(1:n_groups, function(i_groups) subset(ini$mu_eta[[i_groups]],G_eta))
  df     <- n_v*(n_v+3)/2-n_par
  df_b   <- n_v*(n_v+3)/2-2*n_v
  n_obs  <-attributes(data)$obs_size %>% sum

  
  dml<-.dml_cal(sigma=data$raw_cov,
                e_v=e_v,
                mu_v=mu_v,
                sigma_v=sigma_v,
                n_groups=n_groups,
                w_g=w_g,
                n_v=n_v)
  dml_b<-.dml_cal(sigma=data$raw_cov,
                    e_v=e_v,
                    mu_v=lapply((1:n_groups),function(i_groups) as.matrix(e_v[[i_groups]])),
                    sigma_v=lapply((1:n_groups),function(i_groups) diag(diag(data$raw_cov[[i_groups]]))),
                    n_groups=n_groups,
                    w_g=w_g,
                    n_v=n_v)
  rpl  <-sapply(output$value[is.na(output$type)],function(x) .rpl_cal(x, gamma=gamma, type=type, delta=delta)) %>% sum
  dpl  <-dml+rpl
  lrt  <-n_obs*dml
  lrt_b<-n_obs*dml_b
  srmr <-sapply((1:n_groups),function(i_groups) {
    sqrt(
      sum(.ltri((sigma_v[[i_groups]]-data$raw_cov[[i_groups]])^2/tcrossprod(diag(sigma_v[[i_groups]]))))/(n_v*(n_v+1)/2)
      ) 
    }) %>% sum
  
  rmsea<-sqrt(max((lrt-df)/(df*n_obs),0))
  mc   <-exp(-0.5*(lrt-df)/n_obs)
  ghat <-n_v/(n_v+2*(lrt-df)/n_obs)
  cfi  <-1 - max(lrt - df, 0) / max(lrt- df, lrt_b - df_b, 0)
  nnfi <-(lrt_b / df_b - lrt / df)/(lrt_b/ df_b - 1)
  bl89 <-(lrt_b - lrt) / (lrt_b- df)
  rni  <-((lrt_b - df_b) - (lrt - df))/(lrt_b - df_b)
  aic  <-dml + (2 / n_obs) * n_par
  bic  <-dml + (log(n_obs) / n_obs) * n_par

  
  return(list(parameter=output$value,
              optimization=list(pl=type, gamma=gamma, delta=delta, iteration=it, n_par=n_par, dml=dml),
              goodness=list(df=df,rpl=rpl,dpl=dpl,lrt=lrt,srmr=srmr,rmsea=rmsea,mc=mc,ghat=ghat,cfi=cfi,nnfi=nnfi,bl89=bl89,rni=rni,aic=aic,bic=bic)))
  
}

.dml_cal   <- function(sigma=data$raw_cov,e_v=e_v,mu_v,sigma_v=sigma_v,n_groups=n_groups,w_g=w_g,n_v=n_v){
  sigma_v_iv  <- lapply(1:n_groups, function(i_groups) solve(sigma_v[[i_groups]]))
  dml <-
    (sapply(1:n_groups, function(i_groups) {
      w_g[[i_groups]] * (sum(diag(sigma[[i_groups]] %*% sigma_v_iv[[i_groups]])) -
                           log(det(sigma[[i_groups]] %*% sigma_v_iv[[i_groups]]))  -
                           n_v)
    })  +
      sapply(1:n_groups, function(i_groups) {
        w_g[[i_groups]] * (t(e_v[[i_groups]] - mu_v[[i_groups]]) %*% sigma_v_iv[[i_groups]] %*% (e_v[[i_groups]] -
                                                                                                   mu_v[[i_groups]]))
      })) %>% sum
  return(dml)
}

.rpl_cal    <- function(theta, gamma, type, delta) {
  if (type == "l1") {
    rho <- abs(theta) * gamma
  } else if (type == "scad") {
    if (abs(theta) <= gamma) {
      rho <- abs(theta) * gamma
    } else if (gamma < abs(theta) & abs(theta) <= gamma * delta) {
      rho <-
        -(abs(theta) ^ 2 + gamma ^ 2 - 2 * gamma * delta * abs(theta)) / (2 * (delta - 1))
    } else {
      rho <- (gamma ^ 2 * (delta ^ 2 - 1)) / (2 * (delta - 1))
    }
  } else if (type == "mcp") {
    if (abs(theta) <= gamma * delta) {
      rho <- gamma * delta - abs(theta) ^ 2 / (2 * delta)
    } else{
      rho <- 1 / 2 * gamma ^ 2 * delta
    }
  }
  return(rho)
}

invspecify <- function(model, value) {
  split(model, model$group) %>% lapply(., function(w) {
    split(w, w$matrix) %>% lapply(., function(x) {
      if (any(x$col  !=  1)) {
        y  <-  diag(0, attributes(model)$mat %>% attributes %$% n_eta)
      } else {
        y <- matrix(0, attributes(model)$mat %>% attributes %$% n_eta)
      }
      if(value=="type"){
        y[cbind(x$row, x$col)] <- x$type
      } else if(value=="initial"){
        y[cbind(x$row, x$col)] <- x$initial
      } else if (value=="current") {
        y[cbind(x$row, x$col)] <- x$current
      }
      return(y)
    })
  })
}

.is_est    <- function(x){
  x <- is.na(x) | x == 1
  return(x)
} # 判定是否需要估計，如果X為1為需要估計，X為NA為懲罰項，均返回"TRUE"。X為0為不需要估計，返回"FALSE"。

.is_one    <- function(x){
  x <- !is.na(x) & x == 1
  return(x)
}

.ltri <- function(x) {
  x[upper.tri(x == 1)] = 0
  return(x)
}

