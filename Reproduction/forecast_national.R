# ----------------------------------------------------#
# Data preparation: male mortality for 1969-2019
data = read.table("D:/gblc-main/data/Mortality_national_male.txt",head=T)
data_y = data[data[,"Year"]%in%c(1969:2019),]
data_pop = read.table("D:/gblc-main/data/Population_national.txt",head=T)
data_pop_y = data_pop[data_pop[,"Year"]%in%c(1969:2019),]

matrix_motality_male = c()
matrix_pop_male = c()
for (i in 1969:2019) {
  row =  data_y[(data_y[,"Year"]==i & data_y[,"Age"]%in%c(0:84)),"Male"]
  row_pop = data_pop_y[(data_pop_y[,"Year"]==i & data_y[,"Age"]%in%c(0:84)),"Male"]
  motality_comb = data_y[(data_y[,"Year"]==i & !data_y[,"Age"]%in%c(0:84)),"Male"]
  pop_comb = data_pop_y[(data_pop_y[,"Year"]==i & !data_pop_y[,"Age"]%in%c(0:84)),"Male"]
  age_comb = sum(motality_comb*pop_comb)/sum(pop_comb)
  row_comb = c(row, age_comb)
  row_pop_comb = c(row_pop, sum(pop_comb))
  matrix_motality_male = rbind(matrix_motality_male, row_comb)
  matrix_pop_male = rbind(matrix_pop_male, row_pop_comb)
}
colnames(matrix_motality_male) = c(0:84,"85+")
rownames(matrix_motality_male) = c(1969:2019)
colnames(matrix_pop_male) = c(0:84,"85+")
rownames(matrix_pop_male) = c(1969:2019)

#---------------------------------------------------------------
# We implement 13 expanding windows to
# generate forecasts for horizons of h = 1, 2, . . . , 10. 
# More specifically, our training sample for each
# expanding window is from 1969 to 1996, 1997, . . . , 2009. 
# The subsequent ten years of mortality data
# after each training sample are used for testing purposes. 

#--- Table 4.1: MASE of out-of-sample forecasts for national-level male mortality----
source("D:/gblc-main/R/GB_Leecarter.R")
library(demography)
# L–C and H–U   
mase_lca_male = list()
mase_fdm_male = list()
for (j in 1:13) {
  year = dim(matrix_motality_male)[1]
  matrix_motality_train = matrix_motality_male[c(1:(year-10-j+1)),]
  matrix_motality_test = matrix_motality_male[c((year-10-j+2):(year-j+1)),]
  matrix_pop_train = matrix_pop_male[c(1:(year-10-j+1)),] 
  
  # Lca(adjust="none")
  demogdata_motality = demogdata(data =t(matrix_motality_train),
                                 pop=t(matrix_pop_train),
                                 ages=c(0:85), 
                                 years=as.numeric(rownames(matrix_motality_train)),
                                 type="mortality", name="male",label="US")
  

  lca_motality = lca(demogdata_motality, max.age = 85, 
                     adjust="none", chooseperiod = F)
  forec_motality = forecast(lca_motality,h=10)$rate$male
  
  mase_lca = c()
  for (h in 1:10) {
    forecast_mortality=t(forec_motality)[c(1:h),]
    real_mortality=matrix(matrix_motality_test[c(1:h),],nrow=h)
    hist_mortality=matrix_motality_train
    scale = colMeans(abs(utils::head(hist_mortality, -1) - utils::tail(hist_mortality, -1)))
    scale_matrix = t(matrix(rep(scale,h), ncol = h))
    mase = colMeans(abs(forecast_mortality-real_mortality)/
                      scale_matrix) 
    mase_lca = rbind(mase_lca,mase)
  }
  mase_lca_male[[j]] = mase_lca
  
  # fdm() 
  fdm_motality = fdm(demogdata_motality, max.age = 85,method="classical")
  fdm_forec = forecast(fdm_motality,h=10)$rate$male
  
  mase_fdm = c()
  for (h in 1:10) {
    forecast_mortality=t(fdm_forec)[c(1:h),]
    real_mortality=matrix(matrix_motality_test[c(1:h),],nrow=h)
    hist_mortality=matrix_motality_train
    scale = colMeans(abs(utils::head(hist_mortality, -1) - 
                           utils::tail(hist_mortality, -1)))
    scale_matrix = t(matrix(rep(scale,h), ncol = h))
    mase = colMeans(abs(forecast_mortality-real_mortality)/
                      scale_matrix) 
    mase_fdm = rbind(mase_fdm,mase)
  }
  mase_fdm_male[[j]] = mase_fdm
}

mase = c()
for (h in 1:10) {
  mase0 = mean(rowMeans(sapply(mase_lca_male, 
                               function(x){x[h,]})))
  mase = c(mase,mase0)
}
round(mase,3)
# [1] 1.468 1.539 1.606 1.672 1.741 1.814 1.895 1.980 2.066 2.151
mase = c()
for (h in 1:10) {
  mase0 = mean(rowMeans(sapply(mase_fdm_male, 
                               function(x){x[h,]})))
  mase = c(mase,mase0)
}
round(mase,3)
# [1] 0.695 0.835 0.981 1.121 1.262 1.402 1.545 1.681 1.806 1.919

# LightGBM
library(lightgbm)
source("D:/gblc-main/R/lgb.R")

motality_male_list = list()
for (j in 1:13) {
  year = dim(matrix_motality_male)[1]
  age = dim(matrix_motality_male)[2]
  matrix_motality_train = matrix_motality_male[c(1:(year-10-j+1)),]
  matrix_motality_test = matrix_motality_male[c((year-10-j+2):(year-j+1)),]
  for (a in 1:age) {
    ts =list()
    ts$h=10
    ts$age=a-1
    ts$windowID=j
    ts$x = ts(matrix_motality_train[,a],
              start  = names(matrix_motality_train[,a])[1] )
    ts$xx = ts(matrix_motality_test[,a],
               start  = names(matrix_motality_test[,a])[1] )
    motality_male_list[[(j-1)*age+a]]=ts
  }
}

param_best = list(learning_rate=0.05,
                  lambda_l2=0.03,
                  bagging_fraction=0.8,
                  feature_fraction=0.8,
                  num_leaves=30,
                  max_depth=4,
                  min_data_in_leaf=10,
                  objective = "regression")

mase_lgb_male = list()
for (j in 1:13) {
  windowID = sapply(motality_male_list,function(x)x$windowID==j)
  motality_male_list_use = motality_male_list[windowID]
  
  data_train = c()
  data_forec = c()
  motality_global_lgb = list()
  for (a in 1:86) {
    data_local = mortality_data_generate_lgb(ts = motality_male_list_use[[a]],
                                             num_lag=3)
    data_train = rbind(data_train,data_local$data_train)
    data_forec = rbind(data_forec,data_local$data_forec)
    #log(mortality)
    motality_global_lgb$data_train = log(data_train)
    motality_global_lgb$data_forec = log(data_forec)
  }
  
  
  # params = hyperparameter_search (dataset=motality_global_lgb$data_train, 
  #                                 num_model=1,
  #                                 baye_n_iter=10, lgb_nround = 10)
  # Best_par = params[[1]]$Best_Par
  # param_best = list(learning_rate=Best_par["learning_rate"],
  #                   lambda_l2=Best_par["lambda_l2"],
  #                   bagging_fraction=Best_par["bagging_fraction"],
  #                   feature_fraction=Best_par["feature_fraction"],
  #                   num_leaves=Best_par["num_leaves"],
  #                   max_depth=Best_par["max_depth"],
  #                   min_data_in_leaf=Best_par["min_data_in_leaf"],
  #                   objective = "regression")
  
  ff = mortality_lgb_func(param = param_best, nround = 100,
                          dataset.train = motality_global_lgb$data_train,
                          dataset.forec = motality_global_lgb$data_forec, 
                          h=10 )
  ff = exp(ff)
  matrix_motality_train = matrix_motality_male[c(1:(year-10-j+1)),]
  matrix_motality_test = matrix_motality_male[c((year-10-j+2):(year-j+1)),]
  
  mase_us = c()
  for (h in 1:10) {
    forecast_mortality=ff[c(1:h),]
    real_mortality=matrix(matrix_motality_test[c(1:h),],nrow=h)
    hist_mortality=matrix_motality_train
    scale = colMeans(abs(utils::head(hist_mortality, -1) - utils::tail(hist_mortality, -1)))
    scale_matrix = t(matrix(rep(scale,h), ncol = h))
    mase = colMeans(abs(forecast_mortality-real_mortality)/
                      scale_matrix) 
    mase_us = rbind(mase_us,mase)
  }
  mase_lgb_male[[j]] = mase_us
}

mase = c()
for (h in 1:10) {
  mase0 = mean(rowMeans(sapply(mase_lgb_male, 
                               function(x){x[h,]})))
  mase = c(mase,mase0)
}
round(mase,3)
# 0.824 0.870 0.956 1.043 1.129 1.219 1.312 1.410 1.507 1.602

# Boosted Lee-Carter
# W matrix that captures the age structure
matrix_W=diag(c(0,1,rep(2,82),1,0))
for (i in 2:85) {
  for (j in 2:85) {
    if(abs(i-j)==1){
      matrix_W[i,j]=-1
    }
  }
}

# GBLC

mase_nos_male=list()
for (j in 1:13) {
  year = dim(matrix_motality_male)[1]
  matrix_motality_train = matrix_motality_male[c(1:(year-10-j+1)),]
  matrix_motality_test = matrix_motality_male[c((year-10-j+2):(year-j+1)),]
  # The parameters max.interation and threshold together control the number of iterations. 
  fore_GB_lca = gb_forecast(matrix_input=matrix_motality_train, 
                            matrix_W = matrix_W, 
                            lambda=0, h=10,
                            max.interation=100,scale=F,threshold=1e-06)
  mase_us = c()
  for (h in 1:10) {
    forecast_mortality=(fore_GB_lca$final_forecast)[c(1:h),]
    real_mortality=matrix(matrix_motality_test[c(1:h),],nrow=h)
    hist_mortality=matrix_motality_train
    scale = colMeans(abs(utils::head(hist_mortality, -1) - utils::tail(hist_mortality, -1)))
    scale_matrix = t(matrix(rep(scale,h), ncol = h))
    mase = colMeans(abs(forecast_mortality-real_mortality)/
                      scale_matrix) 
    mase_us = rbind(mase_us,mase)
  }
  mase_nos_male[[j]] = mase_us
}

mase = c()
for (h in 1:10) {
  mase0 = mean(rowMeans(sapply(mase_nos_male, 
                       function(x){x[h,]})))
  mase = c(mase,mase0)
}
round(mase,3)
# [1] 0.591 0.683 0.785 0.881 0.973 1.072 1.181 1.294 1.404 1.512

# GBLC-age
mase_agesk_male=list()
for (j in 1:13) {
  year = dim(matrix_motality_male)[1]
  matrix_motality_train = matrix_motality_male[c(1:(year-10-j+1)),]
  matrix_motality_test = matrix_motality_male[c((year-10-j+2):(year-j+1)),]
  
  # Select lambda (The process is time-consuming, more than 10 minutes)
  # We provide optimized parameter lambda

  # lambda_interval = seq(0,0.1,0.005)
  # x=dim(matrix_motality_train)[1]
  # 
  # mase_lambda = c()
  # for (l in lambda_interval) {
  #   mase = c()
  #   for (i in c(1:3)) {
  #     matrix_motality_train_y = matrix_motality_train[c(1:(x-i+1)),]
  #     y=dim(matrix_motality_train_y)[1]
  #     matrix_motality_train_lambda = matrix_motality_train_y[-c((y-h+1):y),]
  #     matrix_motality_train_test_lambda = matrix_motality_train_y[c((y-h+1):y),]
  #     fore_GB_lca = gb_forecast(matrix_input=matrix_motality_train_lambda,
  #                               matrix_W = matrix_W,
  #                               lambda=l, h=h,
  #                               max.interation=30,scale=F,threshold=1e-06)
  #     mase0 = get_mase(fore_GB_lca$final_forecast,
  #                      matrix(matrix_motality_train_test_lambda,nrow=h),
  #                      matrix_motality_train_lambda)
  #     mase = c(mase,mase0)
  #   }
  #   mase_lambda = c(mase_lambda,mean(mase))
  # }
  # lambda_opt = lambda_interval[which.min(mase_lambda)]
  # cat("lambda:",lambda_opt)
  
  lambda_a_m = c(0.025,0.03,0.03,0.015,0.05,0.035,
                 0.02,0.02,0.02,0.02,0.015,0.015,0.015)
  
  fore_GB_lca = gb_forecast(matrix_input=matrix_motality_train, 
                            matrix_W = matrix_W,
                            lambda=lambda_a_m[j], h=10, max.interation=100,
                            threshold=1e-06)
  mase_us = c()
  for (h in 1:10) {
    forecast_mortality=(fore_GB_lca$final_forecast)[c(1:h),]
    real_mortality=matrix(matrix_motality_test[c(1:h),],nrow=h)
    hist_mortality=matrix_motality_train
    scale = colMeans(abs(utils::head(hist_mortality, -1) - utils::tail(hist_mortality, -1)))
    scale_matrix = t(matrix(rep(scale,h), ncol = h))
    mase = colMeans(abs(forecast_mortality-real_mortality)/
                      scale_matrix) 
    mase_us = rbind(mase_us,mase)
  }
  mase_agesk_male[[j]] = mase_us
}

mase = c()
for (h in 1:10) {
  mase0 = mean(rowMeans(sapply(mase_agesk_male, 
                       function(x){x[h,]})))
  mase = c(mase,mase0)
}
round(mase,3)
# [1] 0.580 0.678 0.781 0.876 0.968 1.067 1.176 1.290 1.402 1.511


#---Table 4.2: MCS results for national-level male mortality----
# Need a few minutes!
library(MCS)
ssm = list()
for (h in c(1:10)) {
  MCS_matrix = cbind(
    as.vector(sapply(mase_lca_male,function(x){as.vector(x[h,])})),
    as.vector(sapply(mase_fdm_male,function(x){as.vector(x[h,])})),
    as.vector(sapply(mase_nos_male,function(x){as.vector(x[h,])})),
    as.vector(sapply(mase_agesk_male,function(x){as.vector(x[h,])}))
  )
  colnames(MCS_matrix) = c("LC","HU","GBLC","GBLC-age")
  MCS = MCSprocedure(MCS_matrix,alpha = 0.05,verbose =F)
  ssm[[h]] = attributes(MCS)$Info$model.names
}
ssm

#--- Figure 4.1: Improvement in MASE for national-level male mortality ---
# Figure 4.2: Improvement in MASE for national-level male mortality, from H–U to GBLC-age.
age_data_h=c()
for (h in c(1:10)) {
  fdm=rowMeans(sapply(mase_fdm_male, 
                       function(x){x[h,]}))
  fdm.mean = c()
  for (i in c(1:4)) {
    fdm.mean=c(fdm.mean,mean(fdm[((i-1)*20+1):(i*20)]))
  }
  fdm.mean = c(fdm.mean, mean(fdm[81:86]))
  
  gblc=rowMeans(sapply(mase_nos_male, 
                       function(x){x[h,]}))
  gblc.mean = c()
  for (i in c(1:4)) {
    gblc.mean=c(gblc.mean,mean(gblc[((i-1)*20+1):(i*20)]))
  }
  gblc.mean = c(gblc.mean, mean(gblc[81:86]))
  
  gblc.a=rowMeans(sapply(mase_agesk_male, 
                         function(x){x[h,]}))
  gblc.a.mean = c()
  for (i in c(1:4)) {
    gblc.a.mean=c(gblc.a.mean,mean(gblc.a[((i-1)*20+1):(i*20)]))
  }
  gblc.a.mean = c(gblc.a.mean, mean(gblc.a[81:86]))
  
  age_data=data.frame(cbind(
    c(1:5),
    fdm.mean,
    gblc.mean, gblc.a.mean,
    rep(h,5)
  ))
  age_data_h = rbind(age_data_h, age_data)
}

colnames(age_data_h) = c("Age","HU","GBLC","GBLC-age","Horizon")
age_data_h = data.frame(age_data_h)
age_data_h$improv = age_data_h$GBLC-age_data_h$GBLC.age
age_data_h$improv_1 = age_data_h$HU-age_data_h$GBLC
age_data_h$improv_2 = age_data_h$HU-age_data_h$GBLC.age


library (ggplot2)
library (reshape2)
p1=ggplot(age_data_h,aes(x=Age,y=Horizon,fill=improv)) +
  geom_raster()+ 
  scale_fill_gradient2(low="#D6604D", high="#0571B0", mid="#DEEBF7",
                       name = "Improv",
                       limits = c(-1,1)*max(abs(age_data_h$improv))
  )+
  scale_y_continuous(breaks = c(0, 2,4,6,8,10))+
  scale_x_continuous(breaks = c(1:5), 
                     labels = c("0-19","20-39","40-59"
                                ,"60-79","80-85+"))+
  labs(x="Age",y="Horizon")

p2=ggplot(age_data_h,aes(x=Age,y=Horizon,fill=improv_1)) +
  geom_raster()+ 
  scale_fill_gradient2(low="#D6604D", high="#0571B0", mid="#DEEBF7",
                       name = "Improv",
                       limits = c(-1,1)*max(abs(age_data_h$improv_2))
  )+
  scale_y_continuous(breaks = c(0, 2,4,6,8,10))+
  scale_x_continuous(breaks = c(1:5), 
                     labels = c("0-19","20-39","40-59"
                                ,"60-79","80-85+"))+
  labs(x="Age",y="Horizon")

library(gridExtra)
library(grid)

titles = c("(a) H-U to GBLC","(b) GBLC to GBLC-age")
plots = mapply(arrangeGrob, list(p2,p1),
               bottom = titles, SIMPLIFY=FALSE)
grid.arrange(grobs = plots, ncol=2)


ggplot(age_data_h,aes(x=Age,y=Horizon,fill=improv_2)) +
  geom_raster()+ 
  scale_fill_gradient2(low="#D6604D", high="#0571B0", mid="#DEEBF7",
                       name = "Improv",
                       limits = c(-1,1)*max(abs(age_data_h$improv_2))
                       )+
  scale_y_continuous(breaks = c(0, 2,4,6,8,10))+
  scale_x_continuous(breaks = c(1:5), 
                     labels = c("0-19","20-39","40-59"
                                ,"60-79","80-85+"))+
  labs(x="Age",y="Horizon")

#-------------------------------------------------------------------
# Figure 4.3: Intermediate results of GBLC in the first two iterations 
# based on national-level male mortality (1969–2009). 

matrix_motality_train = matrix_motality_male[c(1:41),]
matrix_motality_test = matrix_motality_male[c(42:51),]
forec_sum = 0
fitted_sum = 0
process_gamma = c()
process_value = c()
kapa_all = c()
b_all = c()
a_all = c()
process_fitted=list()
process_actual=list()
process_fitted_kb=list()
for (n in 1:5) {
  if(n==1){
    ln=TRUE
  }else{ln=FALSE}
  T = dim(matrix_motality_train)[1]
  if(ln==TRUE){
    matrix=log(matrix_motality_train)
  }else{
    matrix = matrix_motality_train
  }
  ax = apply(matrix,2,mean, na.rm=TRUE) 
  matrix_sd = sweep(matrix,2,ax)
  U = svd(matrix_sd)$u
  V = svd(matrix_sd)$v
  d = svd(matrix_sd)$d
  kapa = d[1]*U[,1]*sum(V[,1])
  
  #rwd  
  matrix_sd_fitted = as.matrix(kapa)%*%t(V[,1]/sum(V[,1]))
  kapa_all = rbind(kapa_all,kapa)
  b_all = rbind(b_all,V[,1]/sum(V[,1]))
  a_all = rbind(a_all,ax)
  fitted_motality = sweep(matrix_sd_fitted,2,-ax)
  rwf.model = forecast::rwf(kapa, drift=TRUE,h=10)
  kt.forecast = rwf.model$mean
  forec_kb = as.numeric(kt.forecast)%*%t(V[,1]/sum(V[,1]))
  forec_motality = sweep(forec_kb,2,-ax)
  
  if(ln==TRUE){
    fitted_motality=exp(fitted_motality)
    forec_motality = exp(forec_motality)}
  
  process_fitted[[n]] = fitted_motality
  process_actual[[n]] = matrix_motality_train
  process_fitted_kb[[n]] = matrix_sd_fitted
  
  opt = optim(par=1, obj_motality, method = "BFGS",
              matrix_motality=matrix_motality_train,
              fitted_motality=fitted_motality,
              matrix_W=matrix_W,
              lambda=0, h=10)
  gamma_opt = opt$par
  value_opt = opt$value
  # T*n matrix
  gradient = c()
  for (t in 1:T) {
    gradient_row = ((matrix_motality_train[t,]-
                       gamma_opt*fitted_motality[t,])
                    -0*(matrix_W+t(matrix_W))%*%fitted_motality[t,])
    gradient=rbind(gradient,t(gradient_row))
  }
  matrix_motality_train = gradient
  fitted = fitted_motality
  forec = forec_motality
  gamma = gamma_opt
  process_value = c(process_value,value_opt)
  forec_sum = forec_sum + gamma*forec
  fitted_sum = fitted_sum+gamma*fitted
  if(n>1){
    if(abs(process_value[n-1]-process_value[n])<1e-07) break
  }  
}

a.data=c()
for (n in c(1:5)) {
  data0 = cbind(a_all[n,], rep(n-1,86),c(0:85))
  a.data=rbind(a.data, data0)
}
colnames(a.data) = c("a", "Iteration", "Age")
a.data = data.frame(a.data)

library(latex2exp)
library(scales)
a1 = ggplot()+
  geom_line( data = a.data[a.data$Iteration ==0,], aes(x = Age,y = exp(a)), size=0.6)+
  labs(x="Age",y=TeX('$a_x$'),title="")+
  scale_y_continuous(trans = log_trans(),
                     breaks = trans_breaks("log", function(x) exp(x)),
                     labels = function(x){(round(x,4))})
a2=ggplot()+
  geom_line( data = a.data[a.data$Iteration ==1,], 
             aes(x = Age,y = (a)), size=0.6)+
  labs(x="Age",y=TeX('$a_x$'),title="")

b.data=c()
for (n in c(1:5)) {
  data0 = cbind(b_all[n,], rep(n-1,86),c(0:85))
  b.data=rbind(b.data, data0)
}
colnames(b.data) = c("b", "Iteration", "Age")
b.data = data.frame(b.data)
b1 = ggplot()+
  geom_line( data = b.data[b.data$Iteration ==0,], aes(x = Age,y = (b)), size=0.6)+
  labs(x="Age",y=TeX('$b_x$'),title="")
b2=ggplot()+
  geom_line( data = b.data[b.data$Iteration ==1,], aes(x = Age,y = (b)), size=0.6)+
  labs(x="Age",y=TeX('$b_x$'),title="")


kapa.data=c()
for (n in c(1:5)) {
  data0 = cbind(kapa_all[n,], rep(n-1,41),c(1969:2009))
  kapa.data=rbind(kapa.data, data0)
}
colnames(kapa.data) = c("kapa", "Iteration", "Year")
kapa.data = data.frame(kapa.data)
k1 = ggplot()+
  geom_line( data = kapa.data[kapa.data$Iteration ==0,], 
             aes(x = Year,y = (kapa)), size=0.6)+
  labs(x="Year",y=TeX('$kappa_t$'),title="")
k2=ggplot()+
  geom_line( data = kapa.data[kapa.data$Iteration ==1,], 
             aes(x = Year,y = (kapa)), size=0.6)+
  ylim(-0.1,0.1)+
  labs(x="Year",y=TeX('$kappa_t$'),title="")

# Plot residual
process_resid=list()
for (n in 1:length(process_fitted)) {
  process_resid[[n]]=process_actual[[n]]-process_fitted[[n]]
}
Age0 = c()
for (x in 1:86) {
  Age0 = c(Age0,rep(x-1,41))
}
resid_plot1 = data.frame(Residual= as.vector(process_resid[[1]]),
                         Year = rep(c(1969:2009),86),
                         Age = Age0)
resid_plot2 = data.frame(Residual= as.vector(process_resid[[2]]),
                         Year = rep(c(1969:2009),86),
                         Age = Age0)
resid_plot3 = data.frame(Residual= as.vector(process_resid[[3]]),
                         Year = rep(c(1969:2009),86),
                         Age = Age0)
r1=ggplot()+
  geom_line( data = resid_plot1, aes(x = Age,y = ((Residual)),
                                     group = Year, color =Year), size=0.6)+
  labs(x="Age",y="Residual",title=TeX('$l=1$'))+
  #ylim(-0.01,0.01)+
  scale_color_continuous(type = "viridis",breaks = seq(1969,2019,10))+
  theme(plot.title = element_text(hjust=-0.5))
r2=ggplot()+
  geom_line( data = resid_plot2, aes(x = Age,y = ((Residual)),
                                     group = Year, color =Year), size=0.6)+
  labs(x="Age",y="Residual",title=TeX('$l=2$'))+
  ylim(-0.01,0.01)+
  scale_color_continuous(type = "viridis",breaks = seq(1969,2019,10))+
  theme(plot.title = element_text(hjust=0.5))

mortality_plot = data.frame(Mortality= as.vector((process_actual[[1]])),
                            Year = rep(c(1969:2009),86),
                            Age = Age0)
m1=ggplot()+
  geom_line( data = mortality_plot, aes(x = Age,y = ((Mortality)),
                                        group = Year, color =Year), size=0.6)+
  labs(x="Age",y="Mortality rate",title=TeX('$l=0$'))+
  #ylim(-0.01,0.01)+
  scale_color_continuous(type = "viridis",breaks = seq(1969,2019,10))+
  scale_y_continuous(trans = log_trans(),
                     breaks = trans_breaks("log", function(x) exp(x)),
                     labels = function(x){(round(x,4))})+
  theme(plot.title = element_text( hjust=-0.5))


kb_plot1 = data.frame(kb= as.vector((process_fitted_kb[[1]])),
                      Year = rep(c(1969:2009),86),
                      Age = Age0)
kb1=ggplot()+
  geom_line( data = kb_plot1, aes(x = Age,y = ((kb)),
                                  group = Year, color =Year), size=0.6)+
  labs(x="Age",y="kt*bx",title="l=0")+
  #ylim(-0.01,0.01)+
  scale_color_continuous(breaks = seq(1969,2019,10))

kb_plot2 = data.frame(kb= as.vector((process_fitted_kb[[2]])),
                      Year = rep(c(1969:2009),86),
                      Age = Age0)
kb2=ggplot()+
  geom_line( data = kb_plot2, aes(x = Age,y = ((kb)),
                                  group = Year, color =Year), size=0.6)+
  labs(x="Age",y="kt*bx",title="l=1")+
  ylim(-0.01,0.01)+
  scale_color_continuous(breaks = seq(1969,2019,10))


gridExtra::grid.arrange(m1,r1,
                        a1,a2,
                        k1,k2,
                        b1,b2,
                        ncol = 2)
#------------------------------------------------------------------------
# Appendix
# Table B.1: Winkler score (×100) of out-of-sample forecasts 
# for national-level male mortality

source("D:/gblc-main/R/GB_Leecarter_PI.R")


# L–C
ws_lca=c()
for (j in 1:13) {
  year = dim(matrix_motality_male)[1]
  matrix_motality_train = matrix_motality_male[c(1:(year-10-j+1)),]
  matrix_motality_test = matrix_motality_male[c((year-10-j+2):(year-j+1)),]
  matrix_pop_train = matrix_pop_male[c(1:(year-10-j+1)),] 
  
  demogdata_motality = demogdata(data =t(matrix_motality_train),
                                 pop=t(matrix_pop_train),
                                 ages=c(0:85), 
                                 years=as.numeric(rownames(matrix_motality_train)),
                                 type="mortality", name="male",label="US")
  
  lca_motality = lca(demogdata_motality, max.age = 85,
                     adjust="none", chooseperiod = F)
  fore = forecast(lca_motality,h=10,level = 95)
  forel = t(fore$rate$lower)
  foreu = t(fore$rate$upper)
  
  mase_s_all = c()
  for (i in 1:10) {
    mase_s = winkler_score (lt = matrix(forel[c(1:i),],nrow=i), 
                            ut = matrix(foreu[c(1:i),],nrow=i), 
                            actual = matrix(matrix_motality_test[c(1:i),],nrow=i), 
                            level = 95)
    mase_s_all = c(mase_s_all, mase_s)
  }
  
  ws_lca = cbind(ws_lca,mase_s_all)
  
}

round(rowMeans(ws_lca)*100,3)
# [1] 1.980 2.014 2.066 2.127 2.198 2.268 2.347 2.421 2.495 2.564


# H–U
ws_lca=c()
for (j in 1:13) {
  year = dim(matrix_motality_male)[1]
  matrix_motality_train = matrix_motality_male[c(1:(year-10-j+1)),]
  matrix_motality_test = matrix_motality_male[c((year-10-j+2):(year-j+1)),]
  matrix_pop_train = matrix_pop_male[c(1:(year-10-j+1)),] 
  
  demogdata_motality = demogdata(data =t(matrix_motality_train),
                                 pop=t(matrix_pop_train),
                                 ages=c(0:85), 
                                 years=as.numeric(rownames(matrix_motality_train)),
                                 type="mortality", name="male",label="US")
  
  fdm_motality3 = fdm(demogdata_motality, max.age = 85,method="classical")
  fdm_motality3_forec = forecast(fdm_motality3,h=10,level = 95)
  forel = t(fdm_motality3_forec$rate$lower)
  foreu = t(fdm_motality3_forec$rate$upper)
  
  mase_s_all = c()
  for (i in 1:10) {
    mase_s = winkler_score (lt = matrix(forel[c(1:i),],nrow=i), 
                            ut = matrix(foreu[c(1:i),],nrow=i), 
                            actual = matrix(matrix_motality_test[c(1:i),],nrow=i), 
                            level = 95)
    mase_s_all = c(mase_s_all, mase_s)
  }
  
  ws_lca = cbind(ws_lca,mase_s_all)
  
}
round(rowMeans(ws_lca)*100,3)
# [1] 0.304 0.450 0.612 0.806 1.006 1.201 1.393 1.566 1.724 1.865


#LightGBM
param_best = list(learning_rate=0.05,
                  lambda_l2=0.03,
                  bagging_fraction=0.8,
                  feature_fraction=0.8,
                  num_leaves=30,
                  max_depth=4,
                  min_data_in_leaf=10,
                  objective = "quantile",
                  alpha = 0.025)

motality_male_list = list()
for (j in 1:13) {
  year = dim(matrix_motality_male)[1]
  age = dim(matrix_motality_male)[2]
  matrix_motality_train = matrix_motality_male[c(1:(year-10-j+1)),]
  matrix_motality_test = matrix_motality_male[c((year-10-j+2):(year-j+1)),]
  for (a in 1:age) {
    ts =list()
    ts$h=10
    ts$age=a-1
    ts$windowID=j
    ts$x = ts(matrix_motality_train[,a],
              start  = names(matrix_motality_train[,a])[1] )
    ts$xx = ts(matrix_motality_test[,a],
               start  = names(matrix_motality_test[,a])[1] )
    motality_male_list[[(j-1)*age+a]]=ts
  }
}

lgb_lower = list()
lgb_upper = list()
for (j in 1:13) {
  windowID = sapply(motality_male_list,function(x)x$windowID==j)
  motality_male_list_use = motality_male_list[windowID]
  
  # motality_global_lgb: $data_train 
  # $data_forec : lags used to forecast, true Value is NA
  #  Value L1 L2 L3 L4 L5 L6 L7 L8 L9 L10 (for all ages)
  data_train = c()
  data_forec = c()
  motality_global_lgb = list()
  for (a in 1:86) {
    data_local = mortality_data_generate_lgb(ts = motality_male_list_use[[a]],
                                             num_lag=3)
    data_train = rbind(data_train,data_local$data_train)
    data_forec = rbind(data_forec,data_local$data_forec)
    #log(mortality)
    motality_global_lgb$data_train = log(data_train)
    motality_global_lgb$data_forec = log(data_forec)
  }
  
  param_best$alpha = 0.025
  ff_l = mortality_lgb_func(param = param_best, nround = 100,
                            dataset.train = motality_global_lgb$data_train,
                            dataset.forec = motality_global_lgb$data_forec, 
                            h=10 )
  ff_l = exp(ff_l)
  param_best$alpha = 0.975
  ff_u = mortality_lgb_func(param = param_best, nround = 100,
                            dataset.train = motality_global_lgb$data_train,
                            dataset.forec = motality_global_lgb$data_forec, 
                            h=10 )
  ff_u = exp(ff_u)
  lgb_lower[[j]] =ff_l
  lgb_upper[[j]] =ff_u
  
}


ws_lgb =c()
for (j in c(1:13)) {
  year = dim(matrix_motality_male)[1]
  matrix_motality_train = matrix_motality_male[c(1:(year-10-j+1)),]
  matrix_motality_test = matrix_motality_male[c((year-10-j+2):(year-j+1)),]
  fore_l = lgb_lower[[j]]
  fore_u = lgb_upper[[j]]
  mase_s_all = c()
  for (i in 1:10) {
    mase_s = winkler_score (lt = matrix(fore_l[c(1:i),],nrow=i), 
                            ut = matrix(fore_u[c(1:i),],nrow=i), 
                            actual = matrix(matrix_motality_test[c(1:i),],nrow=i), 
                            level = 95)
    mase_s_all = c(mase_s_all, mase_s)
  }
  ws_lgb = cbind(ws_lgb,mase_s_all)
}
round(rowMeans(ws_lgb)*100,3)
# [1] 4.891 4.887 4.885 4.889 4.892 4.895 4.899 4.902 4.905 4.909

# GBLC (a few minutes!)
ws_lca=c()
for (j in 1:13) {
  year = dim(matrix_motality_male)[1]
  matrix_motality_train = matrix_motality_male[c(1:(year-10-j+1)),]
  matrix_motality_test = matrix_motality_male[c((year-10-j+2):(year-j+1)),]
  
  set.seed(224)
  pi = gb_forecast_pi(matrix_input=matrix_motality_train, 
                      matrix_W,
                      lambda=0, h=10, max.interation=5,
                      scale=FALSE, threshold = 1e-06,
                      state_include = FALSE, states,
                      alpha=c(0.025,0.975),J=100)
  forel = pi$final_quant[[1]]
  foreu = pi$final_quant[[2]]
  
  ws_s_all = c()
  for (i in 1:10) {
    ws_s = winkler_score (lt = matrix(forel[c(1:i),],nrow=i), 
                          ut = matrix(foreu[c(1:i),],nrow=i), 
                          actual = matrix(matrix_motality_test[c(1:i),],
                                          nrow=i), 
                          level = 95)
    ws_s_all = c(ws_s_all, ws_s)
  }
  ws_lca = cbind(ws_lca,ws_s_all)
}

round(rowMeans(ws_lca)*100,3)
# [1] 0.361 0.381 0.396 0.417 0.441 0.471 0.509 0.564 0.619 0.685


# GBLC-age (a few minutes!)
ws_lca=c()
for (j in 1:13) {
  year = dim(matrix_motality_male)[1]
  matrix_motality_train = matrix_motality_male[c(1:(year-10-j+1)),]
  matrix_motality_test = matrix_motality_male[c((year-10-j+2):(year-j+1)),]
  
  set.seed(224)
  #lambda=0.02
  pi = gb_forecast_pi(matrix_input=matrix_motality_train, 
                      matrix_W,
                      lambda=0.02, h=10, max.interation=5,
                      scale=FALSE, threshold = 1e-06,
                      state_include = FALSE, states,
                      alpha=c(0.025,0.975),J=100)
  forel = pi$final_quant[[1]]
  foreu = pi$final_quant[[2]]
  ws_s_all = c()
  for (i in 1:10) {
    ws_s = winkler_score (lt = matrix(forel[c(1:i),],nrow=i), 
                          ut = matrix(foreu[c(1:i),],nrow=i), 
                          actual = matrix(matrix_motality_test[c(1:i),],
                                          nrow=i), 
                          level = 95)
    ws_s_all = c(ws_s_all, ws_s)
  }
  ws_lca = cbind(ws_lca,ws_s_all)
}
round(rowMeans(ws_lca)*100,3)
# [1] 0.362 0.381 0.396 0.416 0.440 0.470 0.507 0.561 0.615 0.680

# Table C.1: MASE of out-of-sample forecasts for all age groups, h=10,
lc = rowMeans(sapply(mase_lca_male, function(x)x[10,]))
hu = rowMeans(sapply(mase_fdm_male, function(x)x[10,]))
lgb = rowMeans(sapply(mase_lgb_male, function(x)x[10,]))
gblc = rowMeans(sapply(mase_nos_male, function(x)x[10,]))
gblcage = rowMeans(sapply(mase_agesk_male, function(x)x[10,]))

national_male_ages_mase = data.frame(lc,hu,lgb,
                                     gblc,gblcage)
head(national_male_ages_mase)
#          lc        hu       lgb      gblc   gblcage
# 0 2.8028063 2.1877107 0.8021056 1.5298579 1.5389705
# 1 1.0861330 1.0405549 0.6899928 1.0563945 1.2013801
# 2 0.7540019 0.7362433 0.4905470 0.6996250 0.6263446
# 3 0.6484622 0.6313851 0.4947108 0.5321738 0.4671170
# 4 0.7303939 0.7371057 0.5581981 0.7679878 0.7654318
# 5 0.5390781 0.4932304 0.6498877 0.5721789 0.5978181