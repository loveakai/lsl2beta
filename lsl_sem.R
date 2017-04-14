lslSEM <- methods::setRefClass(
  Class = "lslSEM",
  fields = list(
    data = "list",
    model = "data.frame",
    knowledge = "list"
  ),
  
  methods = list(
    input = function(raw_obs,
                     var_subset,
                     var_group,
                     obs_subset,
                     obs_weight,
                     raw_cov,
                     raw_mean,
                     obs_size) {
      if (missing(raw_obs)) {
        if (is.list(raw_cov)) {
          output <- list(raw_cov = raw_cov, raw_mean = raw_mean)
        } else {
          output <- list(raw_cov = list(raw_cov),
                         raw_mean = list(raw_mean))
        }
        if (!missing(obs_size)) {
          attr(output, "obs_size") <- obs_size
        }
      } else {
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
          raw_obs %<>% .[obs_subset, var_subset]
          if (is.null(colnames(raw_obs)))
            colnames(raw_obs) <- paste0("v", 1:ncol(raw_obs))
          raw_obs %<>% cbind(group = 1)
          var_group <- ncol(raw_obs)
        } else {
          if (is.character(var_group)) {
            var_group <- which(colnames(raw_obs) %in% (var_group))
          }
          if (missing(var_subset)) {
            var_subset <- (1:ncol(raw_obs)) %>% .[!. %in% var_group]
          }
          raw_o <-
            raw_obs[obs_subset, var_subset]
          if (is.null(colnames(raw_o)))
            colnames(raw_o) <- paste0("v", 1:ncol(raw_o))
          raw_obs <-
            cbind(raw_o, group = raw_obs[, var_group])
        }
        output <- list(
          raw_obs = raw_obs,
          raw_cov = split(raw_obs, raw_obs[, ncol(raw_obs)]) %>% lapply(function(x) {
            cov(x[, -ncol(x)])
          }),
          raw_mean = split(raw_obs, raw_obs[, ncol(raw_obs)]) %>% lapply(function(x) {
            apply(x[, -ncol(x)], 2, mean)
          })
        )
        attr(output, "obs_size") <-
          plyr::count(raw_obs[, ncol(raw_obs)]) %>% .[, 2]
        
      }
      attr(output, "n_groups") <-
        output$raw_cov %>% length
      if (!is.null(colnames(output$raw_cov[[1]]))) {
        attr(output, "v_label") <- colnames(output$raw_cov[[1]])
      }
      attr(output, "g_label") <-
        unique(raw_obs[, ncol(raw_obs)]) %>% as.character
      data <<- output
    },
    
    specify = function(pattern,
                       value,
                       difference,
                       ref_group,
                       auto_scale = T,
                       v_label,
                       f_label) {
      if (!exists("beta_vf", pattern))
        stop("beta_vf must be specified")
      if (missing(v_label)) {
        v_label <- attributes(data)$v_label
      }
      if (missing(f_label)) {
        f_label <- paste0("f", 1:ncol(pattern$beta_vf))
      }
      if (missing(ref_group)) {
        ref_group <- 1L
      }
      if (is.character(ref_group)) {
        ref_group <- which(attributes(data)$g_label == ref_group)
      }
      vf_label <- paste0(v_label, "<-", rep(f_label, each = length(v_label)))
      fv_label <- paste0(f_label, "<-", rep(v_label, each = length(f_label)))
      par_mat_label <-
        sapply(c(v_label, f_label), function(x)
          paste0(c(v_label, f_label), "<-", x)) %>% `rownames<-`(c(v_label, f_label))
      labels <-
        list(
          v_label = v_label,
          f_label = f_label,
          vf_label = vf_label,
          fv_label = fv_label,
          par_mat_label = par_mat_label
        )
      n_groups <-
        attributes(data)$n_groups
      
      par_mat <-
        .par_mat_cal(
          pattern,
          value,
          n_groups = n_groups,
          labels = labels,
          scale = auto_scale,
          ref_group = ref_group,
          data = data
        )
      ref <-
        .par_tab_cal(
          pattern = par_mat$pattern,
          value = list(
            par_mat$value$alpha_r,
            par_mat$value$beta_r,
            par_mat$value$phi_r
          ),
          v_label,
          f_label,
          par_mat_label,
          group = "r"
        )
      inc <-
        lapply((1:n_groups), function(x) {
          .par_tab_cal(
            pattern = par_mat$pattern,
            value = list(
              par_mat$value$alpha_i[[x]],
              par_mat$value$beta_i[[x]],
              par_mat$value$phi_i[[x]]
            ),
            v_label,
            f_label,
            par_mat_label,
            group = names(par_mat$value$alpha_i[x])
          )
        }) %>% do.call(rbind, .) %>% rbind(ref, .)
      output <- within(inc, {
        ini <- value
        rm(value)
      })
      
      output <- output[!(output$matrix == "phi" & (output$row > output$col)),]
      attr(output, "par_mat") <- par_mat
      attr(output, "labels") <- labels
      attr(output, "ref_group") <- attributes(data)$g_label[ref_group]
      model <<- output
    },
    
    learn = function(penalty,
                     gamma,
                     delta,
                     control = list(max_iter, rel_tol)) {
      if (missing(penalty)) {
        pl <- "l1"
      } else {
        pl <- penalty
      }
      if (missing(gamma))
        gamma <- seq(0.025, 0.1, 0.025)
      if (missing(delta))
        delta <- c(2.5)
      if (missing(control)) {
        control <- list(max_iter = 500, rel_tol = 10 ^ (-5))
      } else {
        if (is.null(control$max_iter)) {
          control$max_iter <- 500
        }
        if (is.null(control$rel_tol)) {
          control$rel_tol <- 10 ^ (-5)
        }
      }
      
      
      par_mat <- attributes(model)$par_mat
      allpen <- expand.grid(pl = pl,
                            delta = delta,
                            gamma = gamma)
      
      knowledge <<-
        lapply(1:nrow(allpen), function(x) {
          penalize <- list(pl = pl,
                           delta = allpen[x, 2],
                           gamma = allpen[x, 3])
          .ecm_alg_cal(
            par_mat = par_mat,
            maxit = control[[1]],
            cri = control[[2]],
            penalize = penalize,
            model = model,
            data = data
          )
        })
    }
  )
)
