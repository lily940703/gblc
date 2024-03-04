gb_motality_pi<-function(matrix_motality, matrix_W, 
                      lambda=1, h=1,ln,
                      state_include = FALSE, states,
                      alpha=c(0.025,0.975),J=100){
  T = dim(matrix_motality)[1]
  if(state_include==TRUE){
    lca = leecarter_state_pi(matrix_motality, h=h, ln=ln, states=states,
                             alpha =alpha, J=J)
  }else{
    lca = leecarter_pi(matrix_motality,h=h,ln=ln,
                       alpha =alpha, J=J)
  }
  fitted_motality = lca$point$fitted
  forec_motality = lca$point$mean
  paths = lca$paths
  opt = optim(par=1, obj_motality, method = "BFGS",
              matrix_motality=matrix_motality,
              fitted_motality=fitted_motality,
              matrix_W=matrix_W,
              lambda=lambda, h=h)
  gamma_opt = opt$par
  value_opt = opt$value
  # T*n matrix
  gradient = c()
  for (t in 1:T) {
    gradient_row = ((matrix_motality[t,]-
                       gamma_opt*fitted_motality[t,])
                    -lambda*(matrix_W+t(matrix_W))%*%fitted_motality[t,])
    gradient=rbind(gradient,t(gradient_row))
  }
  return(list(gradient = gradient,
              forec_motality=forec_motality,
              fitted_motality = fitted_motality,
              paths = paths,
              gamma_opt = gamma_opt,
              value_opt = value_opt))
}

gb_motality_pi_ss<-function(matrix_motality, matrix_W_1, matrix_W_2, 
                      lambda_1=1, lambda_2=1,
                      h=1,ln,state_include = FALSE, states,
                      alpha=c(0.025,0.975),J=100){
  T = dim(matrix_motality)[1]
  if(state_include==TRUE){
    lca = leecarter_state_pi(matrix_motality, h=h, ln=ln, states=states,
                             alpha =alpha, J=J)
  }else{
    lca = leecarter_pi(matrix_motality,h=h,ln=ln,
                       alpha =alpha, J=J)
  }
  fitted_motality = lca$point$fitted
  forec_motality = lca$point$mean
  quantiles = lca$quantile
  paths = lca$paths
  opt = optim(par=1, obj_motality_0, method = "BFGS",
              matrix_motality=matrix_motality,
              fitted_motality=fitted_motality,
              matrix_W_1=matrix_W_1,
              matrix_W_2=matrix_W_2,
              lambda_1=lambda_1,
              lambda_2=lambda_2,
              h=h)
  gamma_opt = opt$par
  value_opt = opt$value
  
  gradient = c()
  for (t in 1:T) {
    gradient_row = ((matrix_motality[t,]-gamma_opt*fitted_motality[t,])
                    -lambda_1*(matrix_W_1+t(matrix_W_1))%*%fitted_motality[t,]
                    -lambda_2*(matrix_W_2+t(matrix_W_2))%*%fitted_motality[t,])
    gradient=rbind(gradient,t(gradient_row))
  }
  return(list(gradient = gradient,
              forec_motality=forec_motality,
              quantiles = quantiles,
              paths= paths,
              gamma_opt = gamma_opt,
              value_opt = value_opt))
}

gb_forecast_pi <- function(matrix_input, matrix_W, 
                        lambda=1, h=1, max.interation=50,
                        scale=FALSE, threshold = 1e-07,
                        state_include = FALSE, states,
                        alpha=c(0.025,0.975),J=100){
  forec_sum = 0
  fitted_sum = 0
  forec_inter = list()
  fitted_inter=list()
  # quant_inter = list()
  # quant_sum = list(0,0)
  process_gamma = c()
  process_value = c()
  
  path_sum=list()
  for (j in 1:J) {
    path_sum[[j]] = matrix(0,nrow =h,ncol=dim(matrix_input)[2])
  }    
  for (n in 1:max.interation) {
    if(n==1){
      ln=TRUE
    }else{ln=FALSE}
    gb = gb_motality_pi(matrix_motality=matrix_input,
                     matrix_W = matrix_W,
                     lambda=lambda, h=h,ln=ln,
                     state_include = state_include, states,
                     alpha=alpha,J=J)
    matrix_input = gb$gradient
    fitted = gb$fitted_motality
    forec = gb$forec_motality
    paths = gb$paths
    gamma = gb$gamma_opt
    process_value = c(process_value,gb$value_opt)
    process_gamma = c(process_gamma,gamma)
    forec_sum = forec_sum + gamma*forec
    fitted_sum = fitted_sum+gamma*fitted
    for (j in 1:J) {
      path_sum[[j]] = path_sum[[j]]+paths[[j]]
    }
    
    forec_inter[[n]] = forec_sum
    fitted_inter[[n]] = fitted_sum
    if(n>1){
      if(abs(process_value[n-1]-process_value[n])<threshold) break
    }  
    quantiles_list=list()
    for (a in 1:length(alpha)) {
      quantiles_list[[a]] = matrix(nrow =h,ncol=dim(matrix_input)[2])
    }    
    for (time in 1:h) {
      for (age in 1:86) {
        quantiles = quantile(sapply(path_sum, function(x)x[time,age]),
                             alpha)
        for (a in 1:length(alpha)) {
          quantiles_list[[a]][time,age] = quantiles[a]
        }
      }
    }
    names(quantiles_list) = alpha 
  }
  return(list(final_forecast = forec_sum,
              final_quant = quantiles_list,
              process_fitted = fitted_inter,
              process_forecast = forec_inter,
              process_gamma = process_gamma,
              process_value = process_value))
}

gb_forecast_pi_ss <- function(matrix_input, matrix_W_1, matrix_W_2, 
                          lambda_1=1, lambda_2=1,
                          h=1, max.interation=50,
                        threshold = 1e-07,
                        state_include = TRUE, states,
                        alpha=c(0.025,0.975),J=100){
  forec_sum = 0
  forec_inter = list()
  process_gamma = c()
  process_value = c()
  # quant_sum = list(0,0)
  # quant_inter = list()
  path_sum=list()
  for (j in 1:J) {
    path_sum[[j]] = matrix(0,nrow =h,ncol=dim(matrix_input)[2])
  }
  for (n in 1:max.interation) {
    if(n==1){
      ln=TRUE
    }else{ln=FALSE}
    gb = gb_motality_pi_ss(matrix_motality=matrix_input,
                       matrix_W_1=matrix_W_1,
                       matrix_W_2=matrix_W_2,
                       lambda_1=lambda_1,
                       lambda_2=lambda_2,
                       h=h,ln=ln,
                     state_include = state_include, states,
                     alpha=alpha,J=J)
    matrix_input = gb$gradient
    forec = gb$forec_motality
    quant = gb$quantiles
    gamma = gb$gamma_opt
    paths = gb$paths
    process_value = c(process_value,gb$value_opt)
    process_gamma = c(process_gamma,gamma)
    forec_sum = forec_sum + gamma*forec
    forec_inter[[n]] = forec_sum
    for (j in 1:J) {
      path_sum[[j]] = path_sum[[j]]+paths[[j]]
    }
    if(n>1){
      if(abs(process_value[n-1]-process_value[n])<threshold) break
    } 
    
  }
  quantiles_list=list()
  for (a in 1:length(alpha)) {
    quantiles_list[[a]] = matrix(nrow =h,ncol=dim(matrix_input)[2])
  }    
  for (time in 1:h) {
    for (age in 1:dim(matrix_input)[2]) {
      quantiles = quantile(sapply(path_sum, function(x)x[time,age]),
                           alpha)
      for (a in 1:length(alpha)) {
        quantiles_list[[a]][time,age] = quantiles[a]
      }
    }
  }
  names(quantiles_list) = alpha 
  
  return(list(final_forecast = forec_sum,
              final_quant = quantiles_list,
              process_forecast = forec_inter,
              # process_quant = quant_inter,
              process_gamma = process_gamma,
              process_value = process_value))
}

