calculate_mase_expanding <- function(idx, matrix_W,data_use, data_use_0,
                                     lambda_interval=seq(0,0.1,0.01),
                                     lambda=rep(0,idx),
                                     num.interation =50){
  year = dim(data_use)[1]
  j=idx
  matrix_motality_train = data_use_0[c(1:(year-10-j+1)),]
  matrix_motality_test = data_use[c((year-10-j+2):(year-j+1)),]
  
  h=1
  x=dim(matrix_motality_train)[1]
  
  if(lambda== "selected"){
    mase_lambda = c()
    for (l in lambda_interval) {
      mase = c()
      for (i in c(1:3)) {
        matrix_motality_train_y = matrix_motality_train[c(1:(x-i+1)),]
        y=dim(matrix_motality_train_y)[1]
        matrix_motality_train_lambda = matrix_motality_train_y[-c((y-h+1):y),]
        matrix_motality_train_test_lambda = matrix_motality_train_y[c((y-h+1):y),]
        fore_GB_lca = gb_forecast(matrix_input=matrix_motality_train_lambda,
                                  matrix_W = matrix_W,
                                  lambda=l, h=h,
                                  max.interation=5,scale=F,
                                  threshold=1e-06,
                                  state_include = TRUE, states = states_abb)
        mase0 = get_mase(fore_GB_lca$final_forecast,
                         matrix(matrix_motality_train_test_lambda,nrow=h),
                         matrix_motality_train_lambda)
        mase = c(mase,mase0)
      }
      mase_lambda = c(mase_lambda,mean(mase))
    }
    lambda_opt = lambda_interval[which.min(mase_lambda)]
    cat("lambda for the expanding window", idx, ":", lambda_opt)
  }else{
    lambda_opt=lambda[idx]
  }
  
  fore_GB_lca = gb_forecast(matrix_input=matrix_motality_train, 
                            matrix_W = matrix_W, 
                            lambda=lambda_opt, h=10,
                            max.interation=num.interation,
                            scale=F,threshold=1e-06,
                            state_include = TRUE, 
                            states = states_abb)
  
  mase_s_all_new = c()
  for (i in 1:10) {
    mase_s = get_mase_states(matrix(fore_GB_lca$final_forecast[c(1:i),],nrow=i), 
                             matrix(matrix_motality_test[c(1:i),],nrow=i), 
                             matrix_motality_train, states=states_abb)
    mase_s_all_new = rbind(mase_s_all_new, mase_s)
  }
  colnames(mase_s_all_new)=states_abb
  mase_a_all_new = c()
  for (i in 1:10) {
    mase_a = get_mase_age(matrix(fore_GB_lca$final_forecast[c(1:i),],nrow=i), 
                          matrix(matrix_motality_test[c(1:i),],nrow=i), 
                          matrix_motality_train)
    mase_a_all_new = rbind(mase_a_all_new, mase_a)
  }
  colnames(mase_a_all_new)=c(0:85)
  return(list(mase_state = mase_s_all_new,
              mase_age = mase_a_all_new,
              mase_all = round(rowMeans(mase_s_all_new),4),
              lambda = lambda_opt))
}


calculate_mase_expanding_0 <- function(idx, matrix_W_1, matrix_W_2, 
                                       data_use, data_use_0,
                                     lambda_interval_1=seq(0,0.1,0.01),
                                     lambda_interval_2=seq(0,0.1,0.01),
                                     lambda_1=rep(0,idx),
                                     lambda_2=rep(0,idx),
                                     max.interation = 50,threshold=1e-06){
  year = dim(data_use)[1]
  j=idx
  matrix_motality_train = data_use_0[c(1:(year-10-j+1)),]
  matrix_motality_test = data_use[c((year-10-j+2):(year-j+1)),]
  
  h=1
  x=dim(matrix_motality_train)[1]
  
  lambda_1_opt=lambda_1[idx]
  lambda_2_opt=lambda_2[idx]
  
  fore_GB_lca = gb_forecast_0(matrix_input=matrix_motality_train, 
                              matrix_W_1=matrix_W_1, 
                              matrix_W_2=matrix_W_2, 
                              lambda_1=lambda_1_opt, 
                              lambda_2=lambda_2_opt, h=10,
                            max.interation=max.interation,threshold=threshold,
                            state_include = TRUE, states = states_abb)
  
  mase_s_all_new = c()
  for (i in 1:10) {
    mase_s = get_mase_states(matrix(fore_GB_lca$final_forecast[c(1:i),],nrow=i), 
                             matrix(matrix_motality_test[c(1:i),],nrow=i), 
                             matrix_motality_train, states=states_abb)
    mase_s_all_new = rbind(mase_s_all_new, mase_s)
  }
  colnames(mase_s_all_new)=states_abb
  mase_a_all_new = c()
  for (i in 1:10) {
    mase_a = get_mase_age(matrix(fore_GB_lca$final_forecast[c(1:i),],nrow=i), 
                          matrix(matrix_motality_test[c(1:i),],nrow=i), 
                          matrix_motality_train)
    mase_a_all_new = rbind(mase_a_all_new, mase_a)
  }
  colnames(mase_a_all_new)=c(0:85)
  return(list(mase_state = mase_s_all_new,
              mase_age = mase_a_all_new,
              mase_all = round(rowMeans(mase_s_all_new),4)))
}




calculate_mase_expanding_mcs <- function(idx, matrix_W_1, matrix_W_2, 
                                       data_use, data_use_0,
                                       lambda_interval_1=seq(0,0.1,0.01),
                                       lambda_interval_2=seq(0,0.1,0.01),
                                       lambda_1=rep(0,idx),
                                       lambda_2=rep(0,idx)){
  year = dim(data_use)[1]
  j=idx
  matrix_motality_train = data_use_0[c(1:(year-10-j+1)),]
  matrix_motality_test = data_use[c((year-10-j+2):(year-j+1)),]
  
  h=1
  x=dim(matrix_motality_train)[1]
  
  lambda_1_opt=lambda_1[idx]
  lambda_2_opt=lambda_2[idx]
  
  fore_GB_lca = gb_forecast_0(matrix_input=matrix_motality_train, 
                              matrix_W_1=matrix_W_1, 
                              matrix_W_2=matrix_W_2, 
                              lambda_1=lambda_1_opt, 
                              lambda_2=lambda_2_opt, h=10,
                              max.interation=50,threshold=1e-07,
                              state_include = TRUE, states = states_abb)
  
  mase_all = list()
  for (i in 1:10) {
    mase_s = get_mase_states_age(matrix(fore_GB_lca$final_forecast[c(1:i),],nrow=i), 
                             matrix(matrix_motality_test[c(1:i),],nrow=i), 
                             matrix_motality_train, states=states_abb)
    mase_all[[i]] = mase_s
  }
  
  return(mase_all)
}


calculate_mase_expanding_compare_mcs <- function(idx, data_use, data_use_0,
                                                 data_use_pop, method="L-C",
                                                 states_abb = states_abb){
  year = dim(data_use)[1]
  j=idx
  matrix_motality_train = data_use_0[c(1:(year-10-j+1)),]
  matrix_motality_test = data_use[c((year-10-j+2):(year-j+1)),]
  pop_test = data_use_pop[year-10-j+1,]
  pop_train = data_use_pop[c(1:(year-10-j+1)),]
  
  # #cfdm#
  if(method=="H-B-Y"){
    demogdata_motality = list(rate = list(),
                              pop=list(),
                              ages=c(0:85),
                              years=c(1969:(1969+dim(matrix_motality_train)[1]-1)),
                              type="mortality", label="US",
                              lambda=0)
    for (s in 1:length(states_abb)) {
      states_col = which(stringr::str_detect(colnames(matrix_motality_train),states_abb[s]))
      matrix_motality_train_s = matrix_motality_train[,states_col]
      pop_train_s = pop_train[,states_col]
      demogdata_motality$rate[[s]] = t(matrix_motality_train_s)
      demogdata_motality$pop[[s]] = t(pop_train_s)
    }
    names(demogdata_motality$rate) =  states_abb
    names(demogdata_motality$pop) =  states_abb
    fit.fdm = demography::coherentfdm(demogdata_motality)
    forec.fdm = forecast::forecast(fit.fdm,h=10)
    fore_s_list = lapply(forec.fdm[1:length(states_abb)], function(x){t(x$rate[[1]])})
    fore_s = c()
    for (s in 1:length(states_abb)) {
      fore_s = cbind(fore_s,fore_s_list[[s]])
    }
    colnames(fore_s) = colnames(matrix_motality_train)
    
  }else{
    fore_s = c()
    for (s in 1:length(states_abb)) {
      states_col = which(stringr::str_detect(colnames(matrix_motality_train),states_abb[s]))
      matrix_motality_train_s = matrix_motality_train[,states_col]
      matrix_motality_test_s = matrix_motality_test[,states_col]
      pop_train_s = pop_train[,states_col]
      demogdata_motality = demography::demogdata(data =t(matrix_motality_train_s),
                                     pop=t(pop_train_s),
                                     ages=c(0:85),
                                     years=c(1969:(1968+dim(matrix_motality_train_s)[1])),
                                     type="mortality", name="male",label="US")
      # lee-carter
      if(method=="L-C"){
        lca_motality = demography::lca(demogdata_motality, max.age = 85,
                           adjust="none", chooseperiod = FALSE)
        fore = t(forecast::forecast(lca_motality,h=10)$rate$male)
      }else if(method=="H-U"){
        fdm_motality3 = demography::fdm(demogdata_motality, max.age = 85,method="classical")
        fore = t(forecast::forecast(fdm_motality3,h=10)$rate$male)
      }
      
      fore_s = cbind(fore_s,fore)
    }
    colnames(fore_s) = colnames(matrix_motality_train)
  }
 
  
  
  mase_all = list()
  for (i in 1:10) {
    mase_s = get_mase_states_age(matrix(fore_s[c(1:i),],nrow=i), 
                                 matrix(matrix_motality_test[c(1:i),],nrow=i), 
                                 matrix_motality_train, states=states_abb)
    mase_all[[i]] = mase_s
  }
  
  return(mase_all)
}
