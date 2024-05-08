leecarter_pi <- function(matrix, h=10, ln=FALSE,
                         alpha=c(0.025,0.975),J=100){
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
  matrix_fitted = sweep(matrix_sd_fitted,2,-ax)
  epsilon_hat = matrix-matrix_fitted
  
  rwf.model = forecast::rwf(kapa, drift=TRUE,h=h)
  kt.fitted = rwf.model$fitted
  eta = rwf.model$residuals
  kt.forecast = rwf.model$mean
  gamma = diff(kt.forecast)[1]
  
  tt = c(2:(T-h+1))
  
  y_hat = list()
  for (i in 1:J) {
    j=sample(tt,1)
    eta_j = eta[c(j:(j+9))]
    epsilon_j = epsilon_hat[c(j:(j+9)),]
    
    kapa0 = tail(kapa,1)
    kapa_hat_h = c()
    for (year in 1:h) {
      kapa_hat = kapa0+gamma+eta_j[year]
      kapa0 =kapa_hat
      kapa_hat_h = c(kapa_hat_h,kapa_hat)
    }
    kb = as.matrix(kapa_hat_h)%*%t(V[,1]/sum(V[,1]))
    y_hat[[i]] =  sweep(kb,2,-ax)+epsilon_j
    if(ln==TRUE){
      y_hat[[i]] = exp(y_hat[[i]])
    }
    }
  #percentile
  quantiles_list=list()
  for (a in 1:length(alpha)) {
    quantiles_list[[a]] = matrix(nrow =h,ncol=dim(matrix)[2])
  }

  for (time in 1:h) {
      for (age in 1:86) {
      quantiles = quantile(sapply(y_hat, function(x)x[time,age]),
                           alpha)
      for (a in 1:length(alpha)) {
        quantiles_list[[a]][time,age] = quantiles[a]
      }
      }
    }
  names(quantiles_list) = alpha 
  
  forec_kb = as.numeric(kt.forecast)%*%t(V[,1]/sum(V[,1]))
  forec_motality = sweep(forec_kb,2,-ax)
  if(ln==TRUE){
    matrix_fitted=exp(matrix_fitted)
    forec_motality = exp(forec_motality)
  }
  output=list(quantile = quantiles_list,
              paths = y_hat,
              point = list(fitted = matrix_fitted,
                           mean = forec_motality))
  # for each alpha, a h*86 matrix
  return(output)
}

winkler_score <- function(lt, ut, actual, level) {
  alpha <- 1 - level / 100
  score <- ifelse(
    actual < lt,
    (ut - lt) + (2 / alpha) * (lt - actual),
    ifelse(
      actual > ut,
      (ut - lt) + (2 / alpha) * (actual - ut),
      ut - lt
    )
  )
  # Age-year-specific
  #return(score)
  # Age-specific
  #return(rowMeans(score, na.rm=TRUE))
  # Year-specific
  #return(colMeans(score, na.rm=TRUE))
  # Overall
  return(mean(score, na.rm = TRUE))
}

quantile_score <- function(qt, actual, p) {
  score <- ifelse(
    actual < qt,
    2 * (1 - p) * (qt - actual),
    2 * (p) * (actual - qt)
  )
  # Age-year-specific
  #return(score)
  # Age-specific
  #return(rowMeans(score, na.rm=TRUE))
  # Year-specific
  #return(colMeans(score, na.rm=TRUE))
  # Overall
  return(mean(score, na.rm = TRUE))
}

coverage <- function(lt, ut, actual) {
  row = dim(actual)[1]
  col = dim(actual)[2]
  cover = matrix(nrow = row, ncol = col)
  for(i in 1:row){
    for(j in 1:col){
      if(ut[i,j] >= actual[i,j] & lt[i,j] <= actual[i,j]){
        cover[i,j] = TRUE
      }else{
        cover[i,j] = FALSE
      }
    }
  }
  return(sum(cover)/(row*col))
}


leecarter_state_pi <- function(data,
                               h = 1,
                               ln = FALSE,
                               states,
                               alpha = c(0.025, 0.975),
                               J = 100) {
  # data = matrix_motality_train
  # states=states_abb
  data_states = list()
  for (i in 1:length(states)) {
    data_states_0 = data[, stringr::str_detect(colnames(data), states[i])]
    data_states[[i]] = matrix(as.numeric(data_states_0), nrow = dim(data)[1])
  }
  
  res_states = lapply(
    data_states,
    leecarter_pi,
    h = h,
    ln = ln,
    alpha = alpha,
    J = J
  )
  fitted_motality_states = c()
  forec_motality_states = c()
  
  for (i in 1:length(res_states)) {
    fitted_motality_states = cbind(fitted_motality_states,
                                   res_states[[i]]$point$fitted)
    forec_motality_states = cbind(forec_motality_states,
                                  res_states[[i]]$point$mean)
  }
  paths = list()
  for (j in 1:J) {
    paths_motality_states = c()
    for (i in 1:length(res_states)) {
      paths_motality_states = cbind(paths_motality_states,
                                    res_states[[i]]$paths[[j]])
    }
    paths[[j]] = paths_motality_states
  }
  
  quantiles_list = list()
  for (a in 1:length(alpha)) {
    quantiles = c()
    for (i in 1:length(res_states)) {
      quantiles = cbind(quantiles, res_states[[i]]$quantile[[a]])
    }
    quantiles_list[[a]] = quantiles
  }
  names(quantiles_list) = alpha
  output = list(
    quantile = quantiles_list,
    paths = paths,
    point = list(fitted = fitted_motality_states,
                 mean = forec_motality_states)
  )
  return(output)
} 

obj_motality <- function(matrix_motality,
                         fitted_motality,
                         matrix_W,
                         lambda = 1,
                         gamma,
                         h = 1) {
  T = dim(matrix_motality)[1]
  obj_T = 0
  for (t in 1:T) {
    forec = fitted_motality[t, ]
    
    error = matrix_motality[t, ] - gamma * (forec)
    # distance = lambda*(gamma^2)*t(forec)%*%matrix_W%*%forec
    distance = lambda * t(forec) %*% matrix_W %*% forec
    
    sq_error = t(error) %*% error
    obj = sq_error + distance
    obj_T = obj_T + obj
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
  # quantiles = lca$quantile
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
              # quantiles = quantiles,
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
  
  # set gamma =1
  # gamma_opt = 1
  # value_opt = obj_motality_0(matrix_motality=matrix_motality,
  #                            fitted_motality=fitted_motality,
  #                            matrix_W_1=matrix_W_1,
  #                             matrix_W_2=matrix_W_2,
  #                             lambda_1=lambda_1,
  #                             lambda_2=lambda_2,
  #                            gamma=1,
  #                             h=h)
  
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
    #quantile
    # for (a in 1:length(alpha)) {
    #   quant_sum[[a]] = quant_sum[[a]]+gamma*quant[[a]]
    # }
    # quant_inter[[n]] = quant_sum
    for (j in 1:J) {
      path_sum[[j]] = path_sum[[j]]+paths[[j]]
    }
    
    forec_inter[[n]] = forec_sum
    fitted_inter[[n]] = fitted_sum
    # if(gb$value_opt<1e-05) break
    if(n>1){
      if(abs(process_value[n-1]-process_value[n])<threshold) break
      # if(process_value[n]<threshold) break
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
              # final_quant = quant_sum,
              # process_quant = quant_inter,
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
    #quantile
    # for (a in 1:length(alpha)) {
    #   quant_sum[[a]] = quant_sum[[a]]+gamma*quant[[a]]
    # }
    # quant_inter[[n]] = quant_sum
    for (j in 1:J) {
      path_sum[[j]] = path_sum[[j]]+paths[[j]]
    }
    if(n>1){
      if(abs(process_value[n-1]-process_value[n])<threshold) break
      # if(process_value[n]<threshold) break
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


