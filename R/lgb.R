mortality_data_generate_lgb <- function(ts,num_lag=10){
  data_train = ts$x
  for (i in 1:num_lag) {
    data_train = cbind(data_train,stats::lag(ts$x,-i))
  }
  colnames(data_train)=c("Value",paste0("L",c(1:num_lag)) )
  data_train_0 = data_train[(num_lag+1):(length(ts$x)),]
  data_forec = t(as.matrix(data_train[(length(ts$x)+1),]))
  return(list(data_train = data_train_0,
              data_forec = data_forec))
}

mortality_lgb_func <- function(param, dataset.train,
                     nround, dataset.forec, h ){
  lgb.dtrain <- lgb.Dataset(as.matrix(dataset.train[,-1]), 
                            label = dataset.train[,1])
  lgb.model <- lgb.train(param, lgb.dtrain, nrounds = nround,verbose = -1L)
  
  point_forec = predict(lgb.model, (matrix(dataset.forec[,-1],
                                              nrow = dim(dataset.forec)[1])))
  
  num_lag = dim(dataset.train)[2]-1
  if(h>1){
    xreg = cbind(point_forec,(matrix(dataset.forec[,-1],
                                     nrow = dim(dataset.forec)[1])))
    point_forec_h = point_forec
    for (i in 2:h) {
      point_forec = predict(lgb.model, matrix(xreg[,c(1:num_lag)],
                                              nrow = dim(dataset.forec)[1]))
      xreg = cbind(point_forec,xreg)
      point_forec_h = rbind(point_forec_h,point_forec)
    }
    return(point_forec_h)
  }else
    return(point_forec)
}
