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
