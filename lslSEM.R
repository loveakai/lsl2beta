
#' A Reference Class for Learning a SEM model via penalized likelihood.
#' @field $data
#' A list of stored data.

#' @field $model
#' A list of pattern matrices and starting values for model parameters.

#' @field $structure 
#' A list of learned structure based on the inputed data and specified model.

#' @export lslSEM
#' @exportClass lslSEM
#' @import methods reshape2 ggplot2 lavaan
#' @examples
#' #Example 1: Holzinger and Swineford (1939) Data#
#' lambda <- matrix(NA, 9, 3)
#' lambda[c(1,2,3), 1] <- lambda[c(4,5,6), 2] <- lambda[c(7,8,9), 3] <- 1
#' 
#' rc_sem <- lslSEM()
#' rc_sem$input(raw = lavaan::HolzingerSwineford1939)
#' rc_sem$specify(pattern = list(lambda = lambda))
#' rc_sem$learn(penalty = list(type = "scad"), variable = 7:15)
#' rc_sem$draw(type = "overall")
#' rc_sem$draw(type = "individual", object = "lambda")
#' rc_sem$summarize(type = "overall")
#' rc_sem$summarize(type = "individual")
#'
#' #Example 2: Industrialization and Political Democracy Data#
#' lambda <- matrix(0, 11, 3)
#' lambda[1:4, 1] <- lambda[5:8, 2] <- lambda[9:11, 3] <- 1
#' psi <- diag(1, 11)
#' psi[5, 1] <- psi[4, 2] <- psi[6, 2] <- psi[7, 3] <- psi[8, 4] <- psi[8, 6] <- NA
#' beta <- matrix(0, 3, 3)
#' beta[1, 3] <- beta[2, 1] <- beta[2, 3] <- 1
#' phi <- diag(1, 3)
#' 
#' rc_sem <- lslSEM()
#' rc_sem$input(raw = lavaan::PoliticalDemocracy)
#' rc_sem$specify(pattern = list(lambda = lambda, psi = psi, beta = beta, phi = phi))
#' rc_sem$learn(penalty = list(type = "scad"))
#' rc_sem$draw(type = "overall")
#' rc_sem$draw(type = "individual", object = "psi")
#' rc_sem$summarize(type = "overall")
#' rc_sem$summarize(type = "individual")
#'
#'  
lslSEM <- methods::setRefClass(Class = "lslSEM", 
                               fields = list(
                                 data = "list", 
                                 model = "list", 
                                 structure = "list"),
                               
                               methods = list(
                                 input = function(raw) {
                                   "Method input() is used to input a data set for further statistical analysis. \n
                                   Argument raw is a raw data matrix that will be used for analysis."
                                   data$raw <<- as.data.frame(raw)
                                 },
                                 
                                 specify = function(pattern, value, scale) {
                                   "Method specify() is used to specify a SEM model for fitting.\n
                                   Argument pattern is a list of pattern matrices for the six model parameter matrices in SEM. 
                                   For each pattern matrix, the element can be set as 1, 0, or NA.
                                   1 means that the corresponding parameter should be freely estimated, 
                                   0 indicates that the parameter sould be fixed, 
                                   and NA makes the parameter to be estimated with penalization.
                                   These pattern matrices for the six model parameter matrices include \n
                                   (1) lambda: a P by M pattern matrix for loading matrix; \n
                                   (2) psi: a P by P pattern matrix covariance matrix of measurement errors. Its default value is diag(1, P). 
                                       Note that the diagonal element of psi must be set as 1; \n
                                   (3) beta: a M by M pattern matrix for path coefficient matrix. Its default value is matrix(0, M, M); \n
                                   (4) phi: a M by M pattern matrix for covariance matrix of residuals. Its default value is diag(1, M). 
                                       Note that the diagonal element of phi must be set as 1; \n
                                   (5) nu: a P by 1 pattern matrix for intercepts of observed variables. Its default value is matrix(1, P, 1);\n 
                                   (6) alpha: a M by 1 pattern matrix for intercepts of latent factors. Its default value is matrix(0, M, 1).\n
                                   Argument value is a list of starting value matrices for the six model parameter matrices in SEM.
                                   These starting value matrices for the six model parameter matrices include \n
                                   (1) lambda: a P by M starting value matrix for the loading matrix. The default value for its element is 
                                       (a) 1 for the freely estimated parameters and (2) 0 for the fixed or penealized parameter.  \n
                                   (2) psi: a P by P starting value matrix for the covariance matrix of measurement errors. Its default value is diag(0.1, P). \n
                                   (3) beta: a M by M starting value matrix for the path coefficient matrix. Its default value is matrix(0, M, M); \n
                                   (4) phi: a M by M starting value matrix for the covariance matrix of residuals. Its default value is diag(1, M). \n
                                   (5) nu: a P by 1 starting value matrix for the intercepts of observed variables. Its default value is sample means of observed variables;\n 
                                   (6) alpha: a M by 1 starting value matrix for the intercepts of latent factors. Its default value is matrix(0, M, 1).\n
                                   Argument scale is a logical indicator for specifying whether the scale of latent variable should be determined automatically.
                                   If scale = TRUE, the first freely estimated loading of each latent factor will be set as one for scale setting.
                                   The default value of scale is 1.
                                   "
                                   if (missing(pattern) | is.null(pattern$lambda)) {
                                     stop("lambda matrix in pattern is not specified")
                                   }
                                   
                                   if (is.null(pattern$psi)) {
                                     pattern$psi <- diag(1, dim(pattern$lambda)[1])
                                   }
                                   if (is.null(pattern$beta)) {
                                     pattern$beta <- matrix(0, dim(pattern$lambda)[2], dim(pattern$lambda)[2])
                                   }
                                   if (is.null(pattern$phi)) {
                                     pattern$phi <- matrix(1, dim(pattern$lambda)[2], dim(pattern$lambda)[2])
                                   }
                                   if (is.null(pattern$nu)) {
                                     pattern$nu <- matrix(1, dim(pattern$lambda)[1], 1)
                                   }
                                   if (is.null(pattern$alpha)) {
                                     pattern$alpha <- matrix(0, dim(pattern$lambda)[2], 1)
                                   }
                                   if (missing(value)) {
                                     value <- list()
                                     value$lambda <- 1*(.is_one(pattern$lambda))
                                   }
                                   if (is.null(value$lambda)) {
                                     value$lambda <- 1*(.is_one(pattern$lambda))
                                   }
                                   
                                   if (is.null(value$psi)) {
                                     value$psi <- diag(.1, dim(pattern$lambda)[1])
                                   }
                                   if (is.null(value$beta)) {
                                     value$beta <- matrix(0, dim(pattern$lambda)[2], dim(pattern$lambda)[2])
                                   }
                                   if (is.null(value$phi)) {
                                     value$phi <- diag(1, dim(pattern$lambda)[2])
                                   }
                                   if (is.null(value$nu)) {
                                     value$nu <- matrix(0, dim(pattern$lambda)[1], 1)
                                   }
                                   if (is.null(value$alpha)) {
                                     value$alpha <- matrix(0, dim(pattern$lambda)[2], 1)
                                   }
                                   if (missing(scale)) {
                                     scale <- TRUE
                                   }
                                   if (scale) {
                                     scale_idx <- apply(pattern$lambda == 1, 2, function(x) min(which(x)))
                                     pattern$lambda[cbind(scale_idx, 1:ncol(pattern$lambda))] <- 0
                                   }
                                   pattern$psi[upper.tri(pattern$psi)] <- pattern$psi[lower.tri(pattern$psi)]
                                   pattern$phi[upper.tri(pattern$phi)] <- pattern$phi[lower.tri(pattern$phi)]
                                   
                                   value$psi[upper.tri(value$psi)] <- value$psi[lower.tri(value$psi)]
                                   value$phi[upper.tri(value$phi)] <- value$phi[lower.tri(value$phi)]
                                   
                                   model$pattern <<- pattern  
                                   model$value <<- value
                                 },
                                 
                                 learn = function(penalty, control, variable) {
                                   "Method learn() is used to calculate penalized likelihood estimates based on the inputed data and specified model.\n 
                                   Argument penalty is a list for the information of penalty function and regularization parameters.
                                   Argument penalty includes three elements \n
                                   (1) type = c('l1', 'scad', 'mcp'): the penalty function to be implemented. The default value is 'l1'; \n 
                                   (2) gamma = seq(0.025, .10, .025): the values of regularization parameter gamma to be considered; \n 
                                   (3) delta = 2.5: the values of shape parameter delta to be considered.\n
                                   Argument control is a list for controlling the optimization process. 
                                   Argument control includes two elements: \n       
                                   (1) maxit = 500, the maximum number of ECM iterations; \n 
                                   (2) epsilon = 10^-5, the convergence criterion of ECM algorithm. \n
                                   Argument variable is a vector of index of variable names to specify which variables in data should be used for analysis.
                                   The default value is 1:nrow(data$raw), which select all of the variables."
                                   if (missing(penalty)) {
                                     penalty <- list(type = "l1", gamma = seq(0.025, .10, .025), delta = 2.5)  
                                   } else {
                                     if (is.null(penalty$type)) {
                                       penalty$type <- "l1" 
                                     }
                                     if (is.null(penalty$gamma)) {
                                       penalty$gamma <- seq(0.025, .10, .025)
                                     }
                                     if (is.null(penalty$delta)) {
                                       penalty$delta <- 2.5 
                                     }
                                   }
                                   
                                   if (missing(control)) {
                                     control <- list(maxit = 500, epsilon = 10^-5)
                                   } else {
                                     if (is.null(control$maxit)) {
                                       control$maxit <- 500
                                     }
                                     if (is.null(control$epsilon)) {
                                       control$epsilon <- 10^-5
                                     }
                                   }
                                   
                                   if (missing(variable)) {
                                     variable <- seq(1, ncol(data$raw), 1)
                                   }
                                   structure$obs_moment$Sg <<- cov(data$raw[,variable], use = "pairwise.complete.obs")
                                   structure$obs_moment$mu <<- as.matrix(colMeans(data$raw[,variable]))
                                   
                                   structure$info$N <<- nrow(data$raw)
                                   structure$info$P <<- nrow(model$pattern$lambda)
                                   structure$info$M <<- ncol(model$pattern$lambda)
                                   theta <- .sem_theta_cal(Ldp = model$pattern$lambda, Psp = model$pattern$psi, 
                                                           Btp = model$pattern$beta, Php = model$pattern$phi,
                                                           nup = model$pattern$nu, app = model$pattern$alpha,
                                                           Ld = model$value$lambda, Ps = model$value$psi,
                                                           Bt = model$value$beta, Ph = model$value$phi, 
                                                           nu = model$value$nu, ap = model$value$alpha)
                                   structure$info$Qall <<- length(theta)
                                                                      
                                   A <- length(penalty$gamma)
                                   B <- length(penalty$delta)
                                   penalty_fit <- list(type = penalty$type)
                                   theta_names <- .sem_theta_names(model$pattern$lambda, model$pattern$psi, 
                                                                   model$pattern$beta, model$pattern$phi, 
                                                                   model$pattern$nu, model$pattern$alpha)
                                   gm_name_all <- paste("gm=", as.character(penalty$gamma), sep="")
                                   dt_name_all <- paste("dt=", as.character(penalty$delta), sep="")
                                   structure$overall <<- array(NA, c(18, A, B), dimnames = list(c("gm", "dt", "dpl", "dml", 
                                                                                               "Q", "df", "it", 
                                                                                               "lrt", "srmr", "rmsea", "mc", "ghat", 
                                                                                               "cfi", "nnfi", "bl89", "rni", "aic", "bic"), gm_name_all, dt_name_all))
                                   structure$individual <<- array(NA, c(structure$info$Qall, A, B), dimnames = list(theta_names, gm_name_all, dt_name_all))
                                   for (a in A:1) {
                                     penalty_fit$gm <- penalty$gamma[a]
                                     for (b in B:1) {
                                       penalty_fit$dt <- penalty$delta[b]
                                       fit_summary <- rep(NA, 18)
                                       names(fit_summary) <- c("gm", "dt", "dpl", "dml", 
                                                               "Q", "df", "it", 
                                                               "lrt", "srmr", "rmsea", "mc", "ghat", 
                                                               "cfi", "nnfi", "bl89", "rni", "aic", "bic")
                                       rst_ecm <- .sem_ecm(model$pattern, model$value, penalty_fit, control, structure$info, structure$obs_moment)
                                       Sg <- rst_ecm$implied_moment$Sg
                                       mu <- rst_ecm$implied_moment$mu
                                       N <- structure$info$N
                                       P <- structure$info$P
                                       dml <- rst_ecm$dml
                                       Q <- rst_ecm$Q
                                       lrt <- N * dml
                                       df <- rst_ecm$df
                                       dml_b <- .sem_dml_cal(structure$obs_moment$Sg, structure$obs_moment$mu, 
                                                             diag(diag(structure$obs_moment$Sg)), structure$obs_moment$mu) 
                                       df_b <- P *(P + 3)/2 - 2 * P
                                       lrt_b <- N * dml_b
                                       
                                       fit_summary["gm"] <- penalty_fit$gm
                                       fit_summary["dt"] <- penalty_fit$dt
                                       fit_summary["dpl"] <- rst_ecm$dpl
                                       fit_summary["dml"] <- rst_ecm$dml
                                       fit_summary["Q"] <- rst_ecm$Q
                                       fit_summary["df"] <- rst_ecm$df
                                       fit_summary["it"] <- rst_ecm$it
                                       fit_summary["lrt"] <- lrt
                                       fit_summary["srmr"] <- sqrt(sum(.ltri(((Sg - Sg)^2) / tcrossprod(diag(Sg))))/(P * (P + 1) / 2))
                                       fit_summary["rmsea"] <- sqrt(max((lrt - df) / (df * N), 0))
                                       fit_summary["mc"] <- exp(-0.5 * (lrt - df) / N)
                                       fit_summary["ghat"] <- P / (P + 2 * (lrt - df) / N)
                                       fit_summary["cfi"] <- 1 - max(lrt - df, 0) / max(lrt- df, lrt_b - df_b, 0)
                                       fit_summary["nnfi"] <- (lrt_b / df_b - lrt / df)/(lrt_b/ df_b - 1)
                                       fit_summary["bl89"] <- (lrt_b - lrt) / (lrt_b- df)
                                       fit_summary["rni"] <- ((lrt_b - df_b) - (lrt - df))/(lrt_b - df_b)
                                       fit_summary["aic"] <- dml + (2 / N) * Q
                                       fit_summary["bic"] <- dml + (log(N) / N) * Q
                                       
                                       structure$overall[, a, b] <<- fit_summary
                                       structure$individual[, a, b] <<- rst_ecm$theta
                                     }    
                                   }
                                 },
                                 
                                 summarize = function(type) {
                                   "Method summarize() is used to obtained a summary for the learned structure.
                                     Argument type specify which type of summary should be made. \n
                                   If type = 'overall', summarize() will give the overall model information (including goodness-of-fit indices) under best AIC and BIC models. \n
                                   If type = 'individual', summarize() will give the parameter estimates under best AIC and BIC model."
                                   if (type == "overall") {
                                     select_idx <- which(as.matrix(structure$overall["aic", , ] == min(structure$overall["aic", , ])), arr.ind = T) 
                                     select_idx <- select_idx[dim(select_idx)[1], , drop = F]
                                     aic <- structure$overall[, select_idx[1], select_idx[2]]
                                     select_idx <- which(as.matrix(structure$overall["bic", , ] == min(structure$overall["bic", , ])), arr.ind = T) 
                                     select_idx <- select_idx[dim(select_idx)[1], , drop = F]
                                     bic <- structure$overall[, select_idx[1], select_idx[2]]
                                     rst <- data.frame(aic, bic)
                                     colnames(rst) <- c("aic optimal", "bic optimal")
                                     return(rst)                                     
                                   }
                                   if (type == "individual") {
                                     select_idx <- which(as.matrix(structure$overall["aic", , ] == min(structure$overall["aic", , ])), arr.ind = T) 
                                     select_idx <- select_idx[dim(select_idx)[1], , drop = F]
                                     aic <- structure$individual[, select_idx[1], select_idx[2]]
                                     select_idx <- which(as.matrix(structure$overall["bic", , ] == min(structure$overall["bic", , ])), arr.ind = T) 
                                     select_idx <- select_idx[dim(select_idx)[1], , drop = F]
                                     bic <- structure$individual[, select_idx[1], select_idx[2]]
                                     rst <- data.frame(aic, bic)
                                     colnames(rst) <- c("aic optimal", "bic optimal")
                                     return(rst)    
                                   }

                                 },
                                 draw = function(type, object) {
                                   "Method draw() is used to draw a plot for the learned structure.
                                   Argument type specify which type of plot should be drawn. \n 
                                   If type = 'overall', draw() will give a plot for the goodness-of-fit indices across different values of gamma and delta.
                                   The argument object can be used to specify which goodness-of-fit indices should be plot. 
                                   Its value can be any combination of 'dml', 'srmr', 'rmsea', 'mc', 'ghat', 'cfi', 'nnfi', 'bl89', 'rni', 'aic', and 'bic'.
                                   The default value is c('dml', 'aic', 'bic'). \n
                                   If type = 'individual', draw() will give a plot for the solution paths of parameter estimates under across different values of gamma and delta.
                                   The argument object can be used to specify parameters in which matrix should be plot.
                                   Its value can be 'lambda', 'psi', 'beta', 'phi', 'nu', and 'alpha'.
                                   The default value is 'lambda'."
                                   if (type == "overall") {
                                     if (missing(object)) {
                                       object = c("dml","aic","bic")
                                     }
                                     validation <- reshape2::melt(structure$overall[object, , , drop = F])
                                     validation[, 2] <- as.numeric(substring(validation[, 2], 4))
                                     validation[, 3] <- as.numeric(substring(validation[, 3], 4))
                                     colnames(validation) <- c("criterion", "gm", "dt", "value")
                                     validation_plot <- ggplot2::ggplot(validation, ggplot2::aes(gm, value, colour = criterion)) + ggplot2::geom_line(size = 1) + ggplot2::geom_point(size = 3) + 
                                       ggplot2::facet_grid(.~ dt, labeller = label_both) + 
                                       ggplot2::labs(title = "Plot of Goodness-of-Fit", x = "gamma", y = "value")
                                     print(validation_plot) 
                                   } 
                                   if (type == "individual") {
                                     if (missing(object)) {
                                       object = "lambda"
                                     }
                                     theta_names_mat <- matrix(unlist(strsplit(row.names(structure$individual), split="[][,]")), dim(structure$individual)[1], 3, byrow=T)
                                     matrix_name_idc <- theta_names_mat[, 1] == object
                                     thetap <- .sem_theta_cal(Ldp = model$pattern$lambda, Psp = model$pattern$psi, 
                                                              Btp = model$pattern$beta, Php = model$pattern$phi, 
                                                              nup = model$pattern$nu, app = model$pattern$alpha,
                                                              Ld = model$pattern$lambda, Ps = model$pattern$psi, 
                                                              Bt = model$pattern$beta, Ph = model$pattern$phi, 
                                                              nu = model$pattern$nu, ap = model$pattern$alpha)        
                                     theta_names_penalized <- row.names(structure$individual)[is.na(thetap)]
                                     path <- reshape2::melt(structure$individual[matrix_name_idc, , , drop = F])
                                     path[, 1] <- as.character(path[,1])
                                     path[, 2] <- as.numeric(substring(path[,2],4))
                                     path[, 3] <- as.numeric(substring(path[,3],4))
                                     penalized <- is.element(path[, 1], theta_names_penalized)
                                     path[, 5] <- penalized
                                     colnames(path) <- c("parameter", "gm", "dt", "value", "penalized")
                                     path_plot <- ggplot2::ggplot(path, ggplot2::aes(gm, value, colour = parameter, linetype = penalized)) + ggplot2::geom_line(size = 1) + 
                                       ggplot2::scale_linetype_manual(values = c("FALSE" = "solid", "TRUE" = "dashed")) +
                                       ggplot2::facet_grid(.~ dt, labeller = label_both) +
                                       ggplot2::labs(title = paste("Solution Paths of Parameters in ", object, sep = ""), x = "gamma", y = "value")
                                     print(path_plot)
                                   }
                                 }
                               )
)





