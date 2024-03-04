obj_motality<-function(matrix_motality,
                       fitted_motality, 
                       matrix_W, 
                       lambda=1, 
                       gamma, h=1){
  T = dim(matrix_motality)[1]
  obj_T = 0
  for (t in 1:T) {
    forec = fitted_motality[t,]
    
      error = matrix_motality[t,] - gamma*(forec)
     # distance = lambda*(gamma^2)*t(forec)%*%matrix_W%*%forec
      distance = lambda*t(forec)%*%matrix_W%*%forec
      
    sq_error = t(error)%*%error
    obj = sq_error+distance
    obj_T = obj_T+obj
  }
  return(obj_T)
}

obj_motality_state<-function(matrix_motality,fitted_motality, 
                         matrix_W_1, matrix_W_2, 
                         lambda_1=1, lambda_2=1, 
                         gamma, h=1){
  T = dim(matrix_motality)[1]
  obj_T = 0
  for (t in 1:T) {
    forec = fitted_motality[t,]
    error = matrix_motality[t,] - gamma*(forec)
    distance = lambda_1*t(forec)%*%matrix_W_1%*%forec+
      lambda_2*t(forec)%*%matrix_W_2%*%forec
    sq_error = t(error)%*%error
    obj = sq_error+distance
    obj_T = obj_T+obj
  }
  return(obj_T)
}

gb_motality<-function(matrix_motality, matrix_W, 
                      lambda=1, h=1,ln,
                      state_include = FALSE, states){
  T = dim(matrix_motality)[1]
  if(state_include==TRUE){
    lca = leecarter_state(matrix_motality, h=h, ln=ln, states=states)
  }else{
    lca = leecarter(matrix_motality,h=h,ln=ln)
  }
  fitted_motality = lca$fitted_motality
  forec_motality = lca$forec_motality
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
              gamma_opt = gamma_opt,
              value_opt = value_opt))
}


gb_motality_state<-function(matrix_motality, matrix_W_1, matrix_W_2, 
                      lambda_1=1, lambda_2=1,
                      h=1,ln,state_include = FALSE, states){
  T = dim(matrix_motality)[1]
  if(state_include==TRUE){
    lca = leecarter_state(matrix_motality, h=h, ln=ln, states=states)
  }else{
    lca = leecarter(matrix_motality,h=h,ln=ln)
  }
  fitted_motality = lca$fitted_motality
  forec_motality = lca$forec_motality
  opt = optim(par=1, obj_motality_state, method = "BFGS",
              matrix_motality=matrix_motality,
              fitted_motality=fitted_motality,
              matrix_W_1=matrix_W_1,
              matrix_W_2=matrix_W_2,
              lambda_1=lambda_1,
              lambda_2=lambda_2,
              h=h)
  gamma_opt = opt$par
  value_opt = opt$value
  # T*n matrix
  gradient = c()
  for (t in 1:T) {
    gradient_row = ((matrix_motality[t,]-gamma_opt*fitted_motality[t,])
                    -lambda_1*(matrix_W_1+t(matrix_W_1))%*%fitted_motality[t,]
                    -lambda_2*(matrix_W_2+t(matrix_W_2))%*%fitted_motality[t,])
    gradient=rbind(gradient,t(gradient_row))
  }
  return(list(gradient = gradient,
              forec_motality=forec_motality,
              gamma_opt = gamma_opt,
              value_opt = value_opt))
}

gb_forecast <- function(matrix_input, matrix_W, 
                         lambda=1, h=1, max.interation=50,
                         scale=FALSE, threshold = 1e-07,
                        state_include = FALSE, states){
  forec_sum = 0
  fitted_sum = 0
  forec_inter = list()
  fitted_inter=list()
  process_gamma = c()
  process_value = c()
  for (n in 1:max.interation) {
    if(n==1){
      ln=TRUE
    }else{ln=FALSE}
    gb = gb_motality(matrix_motality=matrix_input,
                     matrix_W = matrix_W,
                     lambda=lambda, h=h,ln=ln,
                     state_include = state_include, states)
    matrix_input = gb$gradient
    fitted = gb$fitted_motality
    forec = gb$forec_motality
    gamma = gb$gamma_opt
    process_value = c(process_value,gb$value_opt)
    process_gamma = c(process_gamma,gamma)
    forec_sum = forec_sum + gamma*forec
    fitted_sum = fitted_sum+gamma*fitted
    forec_inter[[n]] = forec_sum
    fitted_inter[[n]] = fitted_sum    if(n>1){
      if(abs(process_value[n-1]-process_value[n])<threshold) break
    }  
  }
  return(list(final_forecast = forec_sum,
              process_fitted = fitted_inter,
              process_forecast = forec_inter,
              process_gamma = process_gamma,
              process_value = process_value))
}

gb_forecast_state <- function(matrix_input, matrix_W_1, matrix_W_2, 
                          lambda_1=1, lambda_2=1,
                          h=1, max.interation=50,
                        threshold = 1e-07,
                        state_include = FALSE, states){
  forec_sum = 0
  forec_inter = list()
  process_gamma = c()
  process_value = c()
  for (n in 1:max.interation) {
    if(n==1){
      ln=TRUE
    }else{ln=FALSE}
    gb = gb_motality_state(matrix_motality=matrix_input,
                       matrix_W_1=matrix_W_1,
                       matrix_W_2=matrix_W_2,
                       lambda_1=lambda_1,
                       lambda_2=lambda_2,
                       h=h,ln=ln,
                     state_include = state_include, states)
    matrix_input = gb$gradient
    forec = gb$forec_motality
    gamma = gb$gamma_opt
    process_value = c(process_value,gb$value_opt)
    process_gamma = c(process_gamma,gamma)
    forec_sum = forec_sum + gamma*forec
    forec_inter[[n]] = forec_sum
    if(n>1){
      if(abs(process_value[n-1]-process_value[n])<threshold) break
    }  
  }
  return(list(final_forecast = forec_sum,
              process_forecast = forec_inter,
              process_gamma = process_gamma,
              process_value = process_value))
}

