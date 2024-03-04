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

leecarter_state_pi <- function(data, h=1, ln=FALSE
                            , states,
                            alpha=c(0.025,0.975),J=100){
  data_states = list()
  for (i in 1:length(states)) {
  data_states_0 = data[,stringr::str_detect(colnames(data),states[i])]
  data_states[[i]] = matrix(as.numeric(data_states_0),nrow=dim(data)[1])
  }
  
  res_states = lapply(data_states, leecarter_pi, 
                      h=h, ln=ln,alpha=alpha,J=J)
  fitted_motality_states=c()
  forec_motality_states=c()
  
  for (i in 1:length(res_states)) {
    fitted_motality_states=cbind(fitted_motality_states,
                                 res_states[[i]]$point$fitted)
    forec_motality_states = cbind(forec_motality_states,
                                  res_states[[i]]$point$mean)
  }
  paths = list()
  for (j in 1:J) {
    paths_motality_states = c()
    for (i in 1:length(res_states)){
      paths_motality_states = cbind(paths_motality_states,
                                    res_states[[i]]$paths[[j]])
    }
    paths[[j]] = paths_motality_states
  }
  
  quantiles_list =list()
  for (a in 1:length(alpha)) {
    quantiles=c()
    for (i in 1:length(res_states)) {
      quantiles=cbind(quantiles,res_states[[i]]$quantile[[a]])
    }
    quantiles_list[[a]] = quantiles
  }
  names(quantiles_list) = alpha 
  output=list(quantile = quantiles_list,
              paths = paths,
              point = list(fitted = fitted_motality_states,
                           mean = forec_motality_states))
  return(output)
} 

leecarter_state_pi <- function(data,
                               h = 1,
                               ln = FALSE,
                               states,
                               alpha = c(0.025, 0.975),
                               J = 100) {
 
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
