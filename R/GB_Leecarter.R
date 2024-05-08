leecarter <- function(matrix,h=1, ln=FALSE){
  if(ln==TRUE){
    matrix=log(matrix)
  }
  T = dim(matrix)[1]
  ax = apply(matrix,2,mean, na.rm=TRUE) 
  matrix_sd = sweep(matrix,2,ax)
  U = svd(matrix_sd)$u
  V = svd(matrix_sd)$v
  d = svd(matrix_sd)$d
  kapa = d[1]*U[,1]*sum(V[,1])
  
  matrix_sd_fitted = as.matrix(kapa)%*%t(V[,1]/sum(V[,1]))
  
  # arima_kapa = forecast::auto.arima(kapa)
  # fitted_kapa = arima_kapa$fitted
  # forec_kapa = forecast::forecast(arima_kapa,h=h)$mean
  # matrix_sd_fitted = as.matrix(fitted_kapa)%*%t(V[,1]/sum(V[,1]))
  
  matrix_fitted = sweep(matrix_sd_fitted,2,-ax)
  
  # forec_kb = as.numeric(forec_kapa)%*%t(V[,1]/sum(V[,1]))
  # forec_motality = sweep(forec_kb,2,-ax)
  
  rwf.model = forecast::rwf(kapa, drift=TRUE,h=h)
  kt.forecast = rwf.model$mean
  forec_kb = as.numeric(kt.forecast)%*%t(V[,1]/sum(V[,1]))
  forec_motality = sweep(forec_kb,2,-ax)
  
  if(ln==TRUE){
    matrix_fitted=exp(matrix_fitted)
    forec_motality = exp(forec_motality)
  }
  return(list(fitted_motality = matrix_fitted,
              forec_motality = forec_motality))
}

leecarter_seasonal <- function(matrix, h=1, ln=FALSE,
                               sea_len){
  if(ln==TRUE){
    matrix=log(matrix)
  }
  T = dim(matrix)[1]
  ax = apply(matrix,2,mean, na.rm=TRUE) 
  matrix_sd = sweep(matrix,2,ax)
  U = svd(matrix_sd)$u
  V = svd(matrix_sd)$v
  d = svd(matrix_sd)$d
  kapa = d[1]*U[,1]*sum(V[,1])
  
  kapa = ts(kapa,frequency =sea_len)
  arima_kapa = forecast::auto.arima(kapa,allowdrift = TRUE)
  fitted_kapa = arima_kapa$fitted
  forec_kapa = forecast::forecast(arima_kapa,h=h)$mean
  
  matrix_sd_fitted = as.matrix(fitted_kapa)%*%t(V[,1]/sum(V[,1]))
  matrix_fitted = sweep(matrix_sd_fitted,2,-ax)
  
  forec_kb = as.numeric(forec_kapa)%*%t(V[,1]/sum(V[,1]))
  forec_motality = sweep(forec_kb,2,-ax)
  
  if(ln==TRUE){
    matrix_fitted=exp(matrix_fitted)
    forec_motality = exp(forec_motality)
  }
  return(list(fitted_motality = matrix_fitted,
              forec_motality = forec_motality))
}

leecarter_state <- function(data, h=1, ln=FALSE
                            , states){
  # data = matrix_motality_train
  # states=states_abb
  data_states = list()
  for (i in 1:length(states)) {
  data_states_0 = data[,stringr::str_detect(colnames(data),states[i])]
  data_states[[i]] = matrix(as.numeric(data_states_0),nrow=dim(data)[1])
  }
  
  res_states = lapply(data_states, leecarter, h=h, ln=ln)
  fitted_motality_states=c()
  forec_motality_states=c()
  for (i in 1:length(res_states)) {
    fitted_motality_states=cbind(fitted_motality_states,res_states[[i]]$fitted_motality)
    forec_motality_states = cbind(forec_motality_states,res_states[[i]]$forec_motality)
  }
  return(list(fitted_motality = fitted_motality_states,
              forec_motality = forec_motality_states))
} 

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

obj_motality_0<-function(matrix_motality,fitted_motality, 
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

obj_motality_se<-function(matrix_motality,fitted_motality, 
                         gamma, h=1){
  T = dim(matrix_motality)[1]
  obj_T = 0
  for (t in 1:T) {
    forec = fitted_motality[t,]
    error = matrix_motality[t,] - gamma*(forec)
    sq_error = t(error)%*%error
    obj_T = obj_T+sq_error
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

gb_motality_seasonal<-function(matrix_motality, matrix_W, 
                      lambda=1, h=1,ln,iteration=1,sea_len = 1){
  T = dim(matrix_motality)[1]
  if(iteration==1){
    lca = leecarter_seasonal(matrix_motality,h=h,ln=ln,
                             sea_len = sea_len)
    }else{
    lca = leecarter_seasonal(matrix_motality,h=h,ln=ln,
                             sea_len = 1)
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

gb_motality_0<-function(matrix_motality, matrix_W_1, matrix_W_2, 
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
    # gb = gb_motality_seasonal(matrix_motality=matrix_input,
    #                  matrix_W = matrix_W, 
    #                  lambda=lambda, h=h,ln=ln,iteration = n,sea_len = 1)
    matrix_input = gb$gradient
    fitted = gb$fitted_motality
    forec = gb$forec_motality
    gamma = gb$gamma_opt
    process_value = c(process_value,gb$value_opt)
    process_gamma = c(process_gamma,gamma)
    forec_sum = forec_sum + gamma*forec
    fitted_sum = fitted_sum+gamma*fitted
    forec_inter[[n]] = forec_sum
    fitted_inter[[n]] = fitted_sum
   # if(gb$value_opt<1e-05) break
    if(n>1){
      # if(abs(process_value[n-1]-process_value[n])<threshold) break
      if(process_value[n]<threshold) break
    }  
  }
  return(list(final_forecast = forec_sum,
              process_fitted = fitted_inter,
              process_forecast = forec_inter,
              process_gamma = process_gamma,
              process_value = process_value))
}

gb_forecast_seasonal <- function(matrix_input, matrix_W, 
                        lambda=1, h=1, max.interation=50,
                        scale=FALSE, threshold = 1e-07,
                        sea_len = 1){
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
    gb = gb_motality_seasonal(matrix_motality=matrix_input,
                     matrix_W = matrix_W,
                     lambda=lambda, h=h,ln=ln,
                     iteration = n,sea_len = sea_len)
    matrix_input = gb$gradient
    fitted = gb$fitted_motality
    forec = gb$forec_motality
    gamma = gb$gamma_opt
    process_value = c(process_value,gb$value_opt)
    process_gamma = c(process_gamma,gamma)
    forec_sum = forec_sum + gamma*forec
    fitted_sum = fitted_sum+gamma*fitted
    forec_inter[[n]] = forec_sum
    fitted_inter[[n]] = fitted_sum
    # if(gb$value_opt<1e-05) break
    if(n>1){
      if(abs(process_value[n-1]-process_value[n])<threshold) break
      # if(process_value[n]<threshold) break
    }  
  }
  return(list(final_forecast = forec_sum,
              process_fitted = fitted_inter,
              process_forecast = forec_inter,
              process_gamma = process_gamma,
              process_value = process_value))
}

gb_forecast_0 <- function(matrix_input, matrix_W_1, matrix_W_2, 
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
    gb = gb_motality_0(matrix_motality=matrix_input,
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
    # if(gb$value_opt<1e-05) break
    if(n>1){
      if(abs(process_value[n-1]-process_value[n])<threshold) break
      # if(process_value[n]<threshold) break
    }  
  }
  return(list(final_forecast = forec_sum,
              process_forecast = forec_inter,
              process_gamma = process_gamma,
              process_value = process_value))
}


get_mse <- function(forecast_mortality, real_mortality, iteration=F) {
  h=dim(real_mortality)[1]
  if(iteration==F){
    mse = (1/h)*sum(rowMeans((forecast_mortality-real_mortality)^2))
  }else{
    mse = sapply(forecast_mortality, 1, 
           function(x){(1/h)*sum(rowMeans((x-real_mortality)^2))})
  }
  return(mse)
}

get_mape <- function(forecast_mortality, real_mortality, iteration=F) {
  h=dim(real_mortality)[1]
  if(iteration==F){
    mape = (1/h)*sum(rowMeans(abs(forecast_mortality-real_mortality)/
                              abs(real_mortality)))
  }else{
    mape = sapply(forecast_mortality, 1, 
                 function(x){(1/h)*sum(rowMeans(abs(x-real_mortality)/
                                                  abs(real_mortality)))})
  }
  return(mape)
}

get_mase <- function(forecast_mortality, real_mortality, 
                     hist_mortality){
  h=dim(real_mortality)[1]
  scale = colMeans(abs(utils::head(hist_mortality, -1) - utils::tail(hist_mortality, -1)))
  scale_matrix = t(matrix(rep(scale,h), ncol = h))
  mase = (1/h)*sum(rowMeans(abs(forecast_mortality-real_mortality)/
                                scale_matrix))
  return(mase)
}


get_mase_states <- function(forecast_mortality, real_mortality, 
                            hist_mortality, states){
  h=dim(real_mortality)[1]
  mase = c()
  for (s in 1:length(states)) {
    states_col = which(stringr::str_detect(colnames(hist_mortality),states[s]))
    forecast_mortality_s = matrix(forecast_mortality[,states_col],nrow=h)
    real_mortality_s = matrix(real_mortality[,states_col],nrow=h)
    hist_mortality_s = hist_mortality[,states_col]
    mase_s = get_mase(forecast_mortality_s, real_mortality_s,hist_mortality_s)
    mase = c(mase,mase_s)
  }
  return(mase)
}

get_mase_age <- function(forecast_mortality, real_mortality, 
                            hist_mortality){
  h=dim(real_mortality)[1]
  mase = c()
  for (a in 1:86) {
    age_col = seq(from = a,by=86,length=51)
    forecast_mortality_a = matrix(forecast_mortality[,age_col],nrow=h)
    real_mortality_a = matrix(real_mortality[,age_col],nrow=h)
    hist_mortality_a = hist_mortality[,age_col]
    mase_a = get_mase(forecast_mortality_a, real_mortality_a,hist_mortality_a)
    mase = c(mase,mase_a)
  }
  return(mase)
}

get_mase_states_age <- function(forecast_mortality, real_mortality, 
                            hist_mortality, states){
  h=dim(real_mortality)[1]
  mase = c()
  for (s in 1:length(states)) {
    states_col = which(stringr::str_detect(colnames(hist_mortality),states[s]))
    forecast_mortality_s = matrix(forecast_mortality[,states_col],nrow=h)
    real_mortality_s = matrix(real_mortality[,states_col],nrow=h)
    hist_mortality_s = hist_mortality[,states_col]
    scale = colMeans(abs(utils::head(hist_mortality_s, -1) - 
                           utils::tail(hist_mortality_s, -1)))
    scale_matrix = t(matrix(rep(scale,h), ncol = h))
    mase_s = colMeans(abs(forecast_mortality_s-real_mortality_s)/
                                scale_matrix)
    mase = cbind(mase,mase_s)
  }
  colnames(mase) = states
  return(mase)
}
