lslSEM2 <- methods::setRefClass(Class = "lslSEM2",
                                fields = list(
                                  data = "list",
                                  model= "data.frame",
                                  knowledge = "array"),
                                methods= list (
                                  input = list(raw_obs,var_subset,var_group,obs_subset,obs_weight,raw_cov,raw_mean,obs_size){
                                    if (missing(raw_obs)) {
                                      if (is.list(raw_cov)) {
                                        output<-list(raw_cov = raw_cov, raw_mean = raw_mean)
                                      } else {
                                        output<-list(raw_cov = list(raw_cov),raw_mean=list(raw_mean))
                                      }
                                      if (!missing(obs_size)) {attr(output,"obs_size")<-obs_size}
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
                                ))