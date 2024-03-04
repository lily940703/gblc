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
  
  matrix_fitted = sweep(matrix_sd_fitted,2,-ax)
  
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


leecarter_state <- function(data, h=1, ln=FALSE
                            , states){
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
