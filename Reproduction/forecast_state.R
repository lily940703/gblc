data_male = readRDS("D:/gblc-main/data/Mortality_state_male.rds")

states_abb = data_male$State_abb[!duplicated(data_male$State_abb)]

# preprocessing for the zero mortality
data_male_0 = data.frame(data_male)
data_male_0[data_male_0$count==0, "count"]=0.1
data_male_0[data_male_0$count==0.1, "mortality"]=0.1/data_male_0[data_male_0$count==0.1, "total"]

#------------------------------------------------------------------------------#
# The weight matrix of age and state shrinkage
library(splm)
data(usaww)

usaww_matrix = as.matrix(usaww)
usaww_matrix1 = cbind(usaww_matrix[,1], rep(0,48),
                      usaww_matrix[,c(2:7)],rep(0,48),
                      usaww_matrix[,c(8:9)],rep(0,48),
                      usaww_matrix[,c(10:48)])
usaww_matrix2 = rbind(matrix(usaww_matrix1[1,], nrow = 1), matrix(rep(0,51),nrow = 1),
                      usaww_matrix1[c(2:7),],matrix(rep(0,51),nrow = 1),
                      usaww_matrix1[c(8:9),],matrix(rep(0,51),nrow = 1),
                      usaww_matrix1[c(10:48),])
D_s = diag(apply(usaww_matrix2, 1, function(x){sum(x!=0)})) 
A_s = matrix(0, nrow = 51,ncol=51)
for (i in 1:51) {
  for (j in 1:51) {
    if(usaww_matrix2[i,j]!=0){
      A_s[i,j]=1
    }
  }
}
rownames(D_s) = colnames(D_s) = states_abb
rownames(A_s) = colnames(A_s) = states_abb
# DC has two neighbours, Maryland (MD) and Virginia (VA).
A_s["DC","MD"]=1
A_s["DC","VA"]=1
A_s["MD","DC"]=1
A_s["VA","DC"]=1
W_s = D_s-A_s

ages = c(0:85)
D_a = matrix_W=diag(c(0,1,rep(2,82),1,0))
A_a = matrix(0, nrow = 86, ncol=86)
for (i in 2:85) {
  for (j in 2:85) {
    if(abs(i-j)==1){
      A_a[i,j]=1
    }
  }
}
rownames(D_a) = colnames(D_a) = ages
rownames(A_a) = colnames(A_a) = ages
W_a = D_a-A_a

I_s = diag(rep(1,51))
I_a = diag(rep(1,86))

A_comb = kronecker(I_s,A_a) + kronecker(A_s, I_a)
names_comb = c()
for (i in 1:51) {
  name = paste0(states_abb[i],"-",ages)
  names_comb = c(names_comb,name)
}
rownames(A_comb) = colnames(A_comb) = names_comb
D_comb = kronecker(I_s,D_a) + kronecker(D_s, I_a)
rownames(D_comb) = colnames(D_comb) = names_comb
W_comb_as = D_comb - A_comb

# The weight matrix for only state shrinkage
D_a_state = A_a_state = matrix(0, nrow = 86, ncol=86)
D_comb_state = kronecker(I_s,D_a_state) + kronecker(D_s, I_a)
A_comb_state = kronecker(I_s,A_a_state) + kronecker(A_s, I_a)
W_comb_s = D_comb_state - A_comb_state
rownames(W_comb_s) = colnames(W_comb_s) = names_comb

# The weight matrix for only age shrinkage
D_s_age = A_s_age = matrix(0, nrow = 51, ncol=51)
D_comb_age = kronecker(I_s,D_a) + kronecker(D_s_age, I_a)
A_comb_age = kronecker(I_s,A_a) + kronecker(A_s_age, I_a)
W_comb_a = D_comb_age - A_comb_age
rownames(W_comb_a) = colnames(W_comb_a) = names_comb

#-------------------------------------------------------------------------------#

# Boosted Lee Carter
source("D:/gblc-main/R/calculate_expanding.R")
source("D:/gblc-main/R/GB_Leecarter.R")


# Train and test
data = as.matrix(data_male)
data_use = matrix(nrow = length(c(1969:2019)), ncol = dim(W_comb_as)[1])
colnames(data_use) = names_comb
for (s in 1:length(states_abb)) {
  for (a in 1:length(ages)) {
    data_s = data[which(data[,"State_abb"] == states_abb[s]),] 
    data_a = data_s[which(as.numeric(data_s[,"Age"]) == ages[a]),]
    col = paste0(states_abb[s],"-",ages[a])
    data_use[,col]  = as.numeric(data_a[,"mortality"])    
  }
}
data_0 = as.matrix(data_male_0)
data_use_0 = matrix(nrow = length(c(1969:2019)), ncol = dim(W_comb_as)[1])
colnames(data_use_0) = names_comb
for (s in 1:length(states_abb)) {
  for (a in 1:length(ages)) {
    data_s = data_0[which(data_0[,"State_abb"] == states_abb[s]),] 
    data_a = data_s[which(as.numeric(data_s[,"Age"]) == ages[a]),]
    col = paste0(states_abb[s],"-",ages[a])
    data_use_0[,col]  = as.numeric(data_a[,"mortality"])    
  }
}


# GBLC methods are time-consuming for state data (for hours), 
# especially when considering shrinkage.
# Parallel computing is provided below.
# As this part of the calculation takes a lot of time, 
# we provide the final calculation results 
# for the convenience of subsequent analysis.

# If select the optimal lambda, set `lambda="selected"` (add a few extra hours!)
# If lambda is fixed, set lambda to be a vector, e.g.,`lambda=rep(0,13)`.
# We also provide optimized parameters.
lambda_s_m = c(0.004,0.004,0.004,0.000,0.000,0.002,
               0.004,0.004,0.004,0.002,0.002,0.002,0.004)
lambda_a_m = c(0.03,0.03,0.04,0.03,0.02,0.08,0.08,
               0.08,0.05,0.08,0.09,0.10,0.10)

# To check the usability of the code, 
# the reader can make a prediction only once (e.g., idx=10),
# A smaller maximum number of iterations can be set 
# to save running time (e.g., max.interation = 3).

# A quick demo (a few minutes)
# GBLC
gblc_res = calculate_mase_expanding_0(
  idx = 10,
  matrix_W_1 = W_comb_a,
  matrix_W_2 = W_comb_s,
  data_use = data_use,
  data_use_0 = data_use_0,
  lambda_1 = rep(0,13),
  lambda_2 = rep(0,13),
  max.interation = 3,
  threshold = 1e-06
)
# MASE for horizons 1-10
gblc_res$mase_all


# Parallel computing for reproducibility purpose
library(parallel)
library(doParallel)
library(foreach)

# GBLC
cl <- makeCluster(2)
registerDoParallel(cl)
male_mase_expending_gblc <- foreach(i = 1:13) %dopar%
  calculate_mase_expanding_0(idx=i, matrix_W_1=W_comb_a,
                             matrix_W_2=W_comb_s,
                             data_use=data_use, 
                             data_use_0=data_use_0,
                             lambda_1=rep(0,13),
                             lambda_2=rep(0,13),
                             max.interation = 50,threshold = 1e-06)
stopCluster(cl)

# GBLC-state
cl <- makeCluster(2)
registerDoParallel(cl)
male_mase_expending_s <- foreach(i = 1:13) %dopar%
  calculate_mase_expanding_0(idx=i, matrix_W_1=W_comb_a,
                             matrix_W_2=W_comb_s,
                             data_use=data_use, 
                             data_use_0=data_use_0,
                             lambda_1=rep(0,13),
                             lambda_2=lambda_s_m,
                             max.interation = 50,threshold = 1e-06)
stopCluster(cl)

# GBLC-age
cl <- makeCluster(2)
registerDoParallel(cl)
male_mase_expending_a <- foreach(i = 1:13) %dopar%
  calculate_mase_expanding_0(idx=i, matrix_W_1=W_comb_a,
                             matrix_W_2=W_comb_s,
                             data_use=data_use, 
                             data_use_0=data_use_0,
                             lambda_1=lambda_a_m,
                             lambda_2=rep(0,13),
                             max.interation = 50,threshold = 1e-06)
stopCluster(cl)

# GBLC-age-state
cl <- makeCluster(2)
registerDoParallel(cl)
male_mase_expending_ss <- foreach(i = 1:13) %dopar%
  calculate_mase_expanding_0(idx=i, matrix_W_1=W_comb_a,
                             matrix_W_2=W_comb_s,
                             data_use=data_use, 
                             data_use_0=data_use_0,
                             lambda_1=lambda_a_m,
                             lambda_2=lambda_s_m,
                             max.interation = 50,threshold = 1e-06)
stopCluster(cl)
# To facilitate subsequent comparative analysis,
# the aforementioned data are saved in RDS files:
# male_mase_expending_gblc
# male_mase_expending_a
# male_mase_expending_s
# male_mase_expending_ss

#----------------------------------------------------------------------------#
# Compare with mortality forecasting methods 
# The default calculation is for the forecast results of L-C. 
# If you need to calculate H-U and H-B-Y, run the code now commented out. 
# H-U uses the fdm function, H–B–Y uses the coherentfdm function (time-consuming).
  
# population
data = as.matrix(data_male)
data_use_pop = matrix(nrow = length(c(1969:2019)), ncol = dim(W_comb_as)[1])
colnames(data_use_pop) = names_comb
for (s in 1:length(states_abb)) {
  for (a in 1:length(ages)) {
    data_s = data[which(data[,"State_abb"] == states_abb[s]),] 
    data_a = data_s[which(as.numeric(data_s[,"Age"]) == ages[a]),]
    col = paste0(states_abb[s],"-",ages[a])
    data_use_pop[,col]  = as.numeric(data_a[,"total"])    
  }
}


library(demography)
mase_cfdm_all_state_age=list()
mase_fdm_all_state_age=list()
mase_lca_all_state_age=list()
for (j in 1:13) {
  year = dim(data_use)[1]
  matrix_motality_train = data_use_0[c(1:(year-10-j+1)),]
  matrix_motality_test = data_use[c((year-10-j+2):(year-j+1)),]
  pop_test = data_use_pop[year-10-j+1,]
  pop_train = data_use_pop[c(1:(year-10-j+1)),]
  #--------------------------------------------------------#
  # lca fdm #
  fore_s = c()

  for (s in 1:length(states_abb)) {
    states_col = which(stringr::str_detect(colnames(matrix_motality_train),states_abb[s]))
    matrix_motality_train_s = matrix_motality_train[,states_col]
    matrix_motality_test_s = matrix_motality_test[,states_col]
    pop_train_s = pop_train[,states_col]
    demogdata_motality = demogdata(data =t(matrix_motality_train_s),
                                   pop=t(pop_train_s),
                                   ages=c(0:85),
                                   years=c(1969:(1968+dim(matrix_motality_train_s)[1])),
                                   type="mortality", name="male",label="US")
    # lee-carter
    lca_motality = lca(demogdata_motality, max.age = 85,
                       adjust="none", chooseperiod = FALSE)
    fore = t(forecast(lca_motality,h=10)$rate$male)

    # # fdm
    # fdm_motality3 = fdm(demogdata_motality, max.age = 85,method="classical")
    # fore = t(forecast(fdm_motality3,h=10)$rate$male)

    fore_s = cbind(fore_s,fore)
  }
  colnames(fore_s) = colnames(matrix_motality_train)
  #--------------------------------------------------------#
  # #cfdm#
  # demogdata_motality = list(rate = list(),
  #                           pop=list(),
  #                           ages=c(0:85),
  #                           years=c(1969:(1969+dim(matrix_motality_train)[1]-1)),
  #                           type="mortality", label="US",
  #                           lambda=0)
  # for (s in 1:length(states_abb)) {
  #   states_col = which(stringr::str_detect(colnames(matrix_motality_train),states_abb[s]))
  #   matrix_motality_train_s = matrix_motality_train[,states_col]
  #   pop_train_s = pop_train[,states_col]
  #   demogdata_motality$rate[[s]] = t(matrix_motality_train_s)
  #   demogdata_motality$pop[[s]] = t(pop_train_s)
  # }
  # names(demogdata_motality$rate) =  states_abb
  # names(demogdata_motality$pop) =  states_abb
  # fit.fdm = coherentfdm(demogdata_motality)
  # forec.fdm = forecast(fit.fdm,h=10)
  # fore_s_list = lapply(forec.fdm[1:length(states_abb)], function(x){t(x$rate[[1]])})
  # fore_s = c()
  # for (s in 1:length(states_abb)) {
  #   fore_s = cbind(fore_s,fore_s_list[[s]])
  # }
  # colnames(fore_s) = colnames(matrix_motality_train)
  #--------------------------------------------------------#
  
  mase_s_all = c()
  for (i in 1:10) {
    mase_s = get_mase_states(matrix(fore_s[c(1:i),],nrow=i), 
                             matrix(matrix_motality_test[c(1:i),],nrow=i), 
                             matrix_motality_train, states=states_abb)
    mase_s_all = rbind(mase_s_all, mase_s)
  }
  colnames(mase_s_all)=states_abb
  mase_a_all = c()
  for (i in 1:10) {
    mase_a = get_mase_age(matrix(fore_s[c(1:i),],nrow=i), 
                          matrix(matrix_motality_test[c(1:i),],nrow=i), 
                          matrix_motality_train)
    mase_a_all = rbind(mase_a_all, mase_a)
  }
  colnames(mase_a_all)=c(0:85)
  mase_lca_all_state_age[[j]] = list(mase_state = mase_s_all,
                                     mase_age = mase_a_all,
                                     mase_all = round(rowMeans(mase_s_all),4))
}
# To facilitate subsequent comparative analysis
# mase_lca_all_state_age,
# mase_fdm_all_state_age,
# mase_cfdm_all_state_age
# are saved in male_mase_lca_fdm_cfdm.RData

#  LightGBM (time-consuming!)
param_best = list(learning_rate=0.05,
                  lambda_l2=0.03,
                  bagging_fraction=0.8,
                  feature_fraction=0.8,
                  num_leaves=30,
                  max_depth=4,
                  min_data_in_leaf=10,
                  objective = "regression")
mase_lgb_all_state_age = list()
for (j in 1:13) {
  year = dim(data_use)[1]
  matrix_motality_train = data_use_0[c(1:(year-10-j+1)),]
  matrix_motality_test = data_use[c((year-10-j+2):(year-j+1)),]
  
  # motality_global_lgb: $data_train 
  # $data_forec : lags used to forecast, true Value is NA
  #  Value L1 L2 L3 L4 L5 L6 L7 L8 L9 L10 (for all ages)
  
  ff = c()
  for (s in 1:length(states_abb)) {
    states_col = which(stringr::str_detect(colnames(matrix_motality_train),
                                           states_abb[s]))
    matrix_motality_train_s = matrix_motality_train[,states_col]
    data_train = c()
    data_forec = c()
    motality_global_lgb = list()
    for (a in 1:dim(matrix_motality_train_s)[2]) {
      ts_sa = list(x = ts(matrix_motality_train_s[,a]))
      data_local = mortality_data_generate_lgb(ts = ts_sa,
                                               num_lag=2)
      data_train = rbind(data_train,data_local$data_train)
      data_forec = rbind(data_forec,data_local$data_forec)
    }
    motality_global_lgb$data_train = log(data_train)
    motality_global_lgb$data_forec = log(data_forec)
    ff_s = mortality_lgb_func(param = param_best, nround = 100,
                              dataset.train = motality_global_lgb$data_train,
                              dataset.forec = motality_global_lgb$data_forec, 
                              h=10 )
    ff_s = exp(ff_s)
    ff = cbind(ff,ff_s)
  }
  
  mase_s_all = c()
  for (i in 1:10) {
    mase_s = get_mase_states(matrix(ff[c(1:i),],nrow=i), 
                             matrix(matrix_motality_test[c(1:i),],nrow=i), 
                             matrix_motality_train, states=states_abb)
    mase_s_all = rbind(mase_s_all, mase_s)
  }
  colnames(mase_s_all)=states_abb
  mase_a_all = c()
  for (i in 1:10) {
    mase_a = get_mase_age(matrix(ff[c(1:i),],nrow=i), 
                          matrix(matrix_motality_test[c(1:i),],nrow=i), 
                          matrix_motality_train)
    mase_a_all = rbind(mase_a_all, mase_a)
  }
  colnames(mase_a_all)=c(0:85)
  mase_lgb_all_state_age[[j]] = list(mase_state = mase_s_all,
                                     mase_age = mase_a_all,
                                     mase_all = round(rowMeans(mase_s_all),3))
}

# To facilitate subsequent comparative analysis
# mase_lgb_all_state_age
# is saved in male_mase_lgb.RData

#-----------------------------------------------------------------------------
# Table 4.3: MASE of out-of-sample forecasts for state-level male mortality
load("D:/gblc-main/data/male_mase_lca_fdm_cfdm.RData")
cat("L–C")
round(rowMeans(sapply(mase_lca_all_state_age, function(x)x$mase_all)),3)
# 0.700  0.714  0.729  0.743  0.757  0.771  0.786  0.801  0.817  0.832 

cat("H–U")
round(rowMeans(sapply(mase_fdm_all_state_age, function(x)x$mase_all)),3)
# 0.637  0.655  0.672  0.688  0.705  0.722  0.739  0.756  0.773  0.790 

cat("H–B–Y")
round(rowMeans(sapply(mase_cfdm_all_state_age, function(x)x$mase_all)),3)
# 0.622  0.635  0.649  0.664  0.679  0.694  0.710  0.725  0.740  0.755 


load("D:/gblc-main/data/male_mase_lgb.RData")
cat("LightGBM")
round(rowMeans(sapply(mase_lgb_all_state_age, function(x)x$mase_all)),3)
# 0.602  0.628  0.656  0.685  0.714  0.742  0.771  0.800  0.828  0.855 


male_mase_expending_NOshrinkage = readRDS("D:/gblc-main/data/male_mase_expending_gblc")
male_mase_expending_age_optlamb = readRDS("D:/gblc-main/data/male_mase_expending_a")
male_mase_expending_state_optlamb = readRDS("D:/gblc-main/data/male_mase_expending_s")
male_mase_expending_ss = readRDS("D:/gblc-main/data/male_mase_expending_ss")
  
cat("GBLC")
round(rowMeans(sapply(male_mase_expending_NOshrinkage, function(x)x$mase_all)),3)
# 0.614  0.628  0.644  0.661  0.678  0.696  0.717  0.738  0.760  0.781 

cat("GBLC-age")
round(rowMeans(sapply(male_mase_expending_age_optlamb, function(x)x$mase_all)),3)
# 0.608  0.623  0.639  0.657  0.675  0.694  0.706  0.718  0.731  0.743 

cat("GBLC-state")
round(rowMeans(sapply(male_mase_expending_state_optlamb, function(x)x$mase_all)),3)
# 0.612  0.624  0.641  0.657  0.676  0.694  0.714  0.735  0.755  0.776 

cat("GBLC-age-state")
round(rowMeans(sapply(male_mase_expending_ss, function(x)x$mase_all)),3)
# 0.605  0.620  0.636  0.654  0.671  0.691  0.702  0.715  0.727  0.739 



#----------------------------------------------------------------------------#
# Table 4.4: MCS results for state-level male mortality
#  (time-consuming!)
library(MCS)
for (h in c(1:10)) {
    MCS_matrix = cbind(
      as.vector(sapply(mase_lca_all_state_age,function(x){x$mase_state[h,]})),
      as.vector(sapply(mase_fdm_all_state_age,function(x){x$mase_state[h,]})),
      as.vector(sapply(mase_lgb_all_state_age,function(x){x$mase_state[h,]})),
      as.vector(sapply(male_mase_expending_NOshrinkage,function(x){x$mase_state[h,]})),
      as.vector(sapply(male_mase_expending_age_optlamb,function(x){x$mase_state[h,]})),
      as.vector(sapply(male_mase_expending_state_optlamb,function(x){x$mase_state[h,]})),
      as.vector(sapply(male_mase_expending_ss,function(x){x$mase_state[h,]}))
    )
    colnames(MCS_matrix) = c("LC","HU","LightGBM",
                             "GBLC","GBLC-age",
                             "GBLC-state", "GBLC-age-state")
    MCS = MCSprocedure(MCS_matrix,alpha = 0.05,B=1000,verbose =F)
    ssm = attributes(MCS)$Info$model.names
    cat(ssm,"\n")
}

#---------------------------------------------------------
# Figure 4.4: Improvement in MASE for state-level male mortality,
library (ggplot2)
library (reshape2)
library(gridExtra)


age_data_h=c()
for (h in c(1:10)) {
  cfdm=rowMeans(sapply(mase_cfdm_all_state_age, 
                       function(x){x$mase_age[h,]}))
  cfdm.mean = c()
  for (i in c(1:4)) {
    cfdm.mean=c(cfdm.mean,mean(cfdm[((i-1)*20+1):(i*20)]))
  }
  cfdm.mean = c(cfdm.mean, mean(cfdm[81:86]))
  
  gblc=rowMeans(sapply(male_mase_expending_NOshrinkage, 
                       function(x){x$mase_age[h,]}))
  gblc.mean = c()
  for (i in c(1:4)) {
    gblc.mean=c(gblc.mean,mean(gblc[((i-1)*20+1):(i*20)]))
  }
  gblc.mean = c(gblc.mean, mean(gblc[81:86]))
  
  gblc.a=rowMeans(sapply(male_mase_expending_age_optlamb, 
                         function(x){x$mase_age[h,]}))
  gblc.a.mean = c()
  for (i in c(1:4)) {
    gblc.a.mean=c(gblc.a.mean,mean(gblc.a[((i-1)*20+1):(i*20)]))
  }
  gblc.a.mean = c(gblc.a.mean, mean(gblc.a[81:86]))
  
  gblc.s=rowMeans(sapply(male_mase_expending_state_optlamb, 
                         function(x){x$mase_age[h,]}))
  gblc.s.mean = c()
  for (i in c(1:4)) {
    gblc.s.mean=c(gblc.s.mean,mean(gblc.s[((i-1)*20+1):(i*20)]))
  }
  gblc.s.mean = c(gblc.s.mean, mean(gblc.s[81:86]))
  
  gblc.as=rowMeans(sapply(male_mase_expending_ss, 
                          function(x){x$mase_age[h,]}))
  gblc.as.mean = c()
  for (i in c(1:4)) {
    gblc.as.mean=c(gblc.as.mean,mean(gblc.as[((i-1)*20+1):(i*20)]))
  }
  gblc.as.mean = c(gblc.as.mean, mean(gblc.as[81:86]))
  
  age_data=data.frame(cbind(
    c(1:5),
    cfdm.mean, gblc.mean, gblc.a.mean,
    gblc.s.mean, gblc.as.mean,
    rep(h,5)
  ))
  age_data_h = rbind(age_data_h, age_data)
}


colnames(age_data_h) = c("Age","HBY","GBLC","GBLC-age",
                         "GBLC-state","GBLC-age-state","Horizon")
age_data_h = data.frame(age_data_h)
age_data_h$improv_1 = age_data_h$HBY-age_data_h$GBLC
age_data_h$improv_2 = age_data_h$GBLC-age_data_h$GBLC.age
age_data_h$improv_3 = age_data_h$GBLC.age-age_data_h$GBLC.age.state
age_data_h$improv_4 = age_data_h$HBY-age_data_h$GBLC.age.state


p1=ggplot(age_data_h,aes(x=Age,y=Horizon,fill=improv_1)) +
  geom_raster()+ 
  scale_fill_gradient2(low="#D6604D", high="#0571B0", mid="#DEEBF7",
                       name = "Improv",
                       limits = c(-1,1)*max(abs(age_data_h$improv_4))
  )+
  scale_y_continuous(breaks = c(0, 2,4,6,8,10))+
  scale_x_continuous(breaks = c(1:5), 
                     labels = c("0-19","20-39","40-59"
                                ,"60-79","80-85+"))+
  labs(x="Age",y="Horizon")


p2=ggplot(age_data_h,aes(x=Age,y=Horizon,fill=improv_2)) +
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


p3=ggplot(age_data_h,aes(x=Age,y=Horizon,fill=improv_3)) +
  geom_raster()+ 
  scale_fill_gradient2(low="#D6604D", high="#0571B0", mid="#DEEBF7",
                       name = "Improv",
                       limits = c(-1,1)*max(abs(age_data_h$improv_3))
  )+
  scale_y_continuous(breaks = c(0, 2,4,6,8,10))+
  scale_x_continuous(breaks = c(1:5), 
                     labels = c("0-19","20-39","40-59"
                                ,"60-79","80-85+"))+
  labs(x="Age",y="Horizon")


titles = c("(a) H-B-Y to GBLC","(b) GBLC to GBLC-age",
           "(c) GBLC-age to GBLC-age-state")
plots = mapply(arrangeGrob, list(p1,p2,p3),
               bottom = titles, SIMPLIFY=FALSE)
grid.arrange(grobs = plots, ncol=3)


ggplot(age_data_h,aes(x=Age,y=Horizon,fill=improv_4)) +
  geom_raster()+ 
  scale_fill_gradient2(low="#D6604D", high="#0571B0", mid="#DEEBF7",
                       name = "Improv",
                       limits = c(-1,1)*max(abs(age_data_h$improv_4))
  )+
  scale_y_continuous(breaks = c(0, 2,4,6,8,10))+
  scale_x_continuous(breaks = c(1:5), 
                     labels = c("0-19","20-39","40-59"
                                ,"60-79","80-85+"))+
  labs(x="Age",y="Horizon")

#--------------------------------------------------------
# Table 4.5: MASE of out-of-sample forecasts for selected states, h=10, males
state_rank=data.frame(cbind(
  rowMeans(sapply(mase_lca_all_state_age, 
                  function(x){x$mase_state[10,]})),
  rowMeans(sapply(mase_fdm_all_state_age, 
                  function(x){x$mase_state[10,]})),
  rowMeans(sapply(mase_cfdm_all_state_age, 
                  function(x){x$mase_state[10,]})),
  rowMeans(sapply(mase_lgb_all_state_age, 
                  function(x){x$mase_state[10,]})),
  rowMeans(sapply(male_mase_expending_NOshrinkage, 
                  function(x){x$mase_state[10,]})),
  rowMeans(sapply(male_mase_expending_age_optlamb, 
                  function(x){x$mase_state[10,]})),
  rowMeans(sapply(male_mase_expending_state_optlamb, 
                  function(x){x$mase_state[10,]})),
  rowMeans(sapply(male_mase_expending_ss, 
                  function(x){x$mase_state[10,]}))
))
colnames(state_rank) = c("L-C","H-U","H-B-Y","LGB","GBLC","GBLC-age",
                         "GBLC-state","GBLC-age-state") 

state_rank_imprv = data.frame(state = rownames(state_rank),
                              GBLC = state_rank$GBLC,
                              GBLCage = state_rank$`GBLC-age`,
                              imprv1 = state_rank$GBLC-state_rank$`GBLC-state`,
                              GBLCstate = state_rank$`GBLC-state`,
                              GBLCagestate = state_rank$`GBLC-age-state`,
                              imprv2 = state_rank$`GBLC-age`-state_rank$`GBLC-age-state`)
head(state_rank_imprv[order(state_rank_imprv$imprv1,decreasing = TRUE),])
#    state      GBLC   GBLCage      imprv1 GBLCstate GBLCagestate      imprv2
# 36    OH 0.9021320 0.8826705 0.010919450 0.8912125    0.8687324 0.013938051
# 13    ID 0.7370361 0.6934329 0.010763052 0.7262731    0.6845977 0.008835141
# 6     CO 0.6524722 0.6307633 0.010715713 0.6417565    0.6118072 0.018956140
# 28    NE 0.7993532 0.7549011 0.010569321 0.7887839    0.7419147 0.012986418
# 45    UT 0.6827320 0.6383565 0.010440160 0.6722919    0.6290616 0.009294963
# 51    WY 1.0390874 0.9698567 0.009602884 1.0294845    0.9607547 0.009101998
head(state_rank_imprv[order(state_rank_imprv$imprv2,decreasing = TRUE),])
#    state      GBLC   GBLCage      imprv1 GBLCstate GBLCagestate      imprv2
# 6     CO 0.6524722 0.6307633 0.010715713 0.6417565    0.6118072 0.018956140
# 16    IA 0.7572824 0.7236528 0.008510858 0.7487716    0.7078507 0.015802152
# 36    OH 0.9021320 0.8826705 0.010919450 0.8912125    0.8687324 0.013938051
# 28    NE 0.7993532 0.7549011 0.010569321 0.7887839    0.7419147 0.012986418
# 49    WV 0.9759804 0.9399623 0.004224633 0.9717557    0.9292608 0.010701503
# 17    KS 0.7561514 0.7155624 0.007151211 0.7490002    0.7057592 0.009803154

#--------------------------------------------------------
# Figure 4.6: Improvement in MASE for state-level male mortality from GBLC to GBLC-state, h=10.
# Figure 4.7: Improvement in MASE for state-level male mortality from GBLC-age to GBLC-age-state, h=10.


library(usmap) #import the package
library(ggplot2) #use ggplot2 to add layer for visualization
improve = data.frame(state =state_rank_imprv$state
                     , improvement1 = (state_rank_imprv$imprv1),
                     improvement2 = (state_rank_imprv$imprv2))

plot_usmap(data = improve, values = "improvement1", 
           include = state.abb, 
           color = "black") + 
  scale_fill_distiller(
    type = "div",
    palette = "RdBu",
    name = "The improvement \n in MASE",
    direction=1,
    limits = c(-1,1)*max(abs(improve$improvement1))
  )+
  theme(legend.position = "right")

plot_usmap(data = improve, values = "improvement2", 
           include = state.abb, 
           color = "black") + 
  scale_fill_distiller(
    type = "div",
    palette = "RdBu",
    name = "The improvement \n in MASE",
    direction=1,
    limits = c(-1,1)*max(abs(improve$improvement2))
  )+
  theme(legend.position = "right")


#-------------------Appendix--------------------------------
# Table A.2: MASE of out-of-sample forecasts for all states and District of Columbia, h=10, males

state_rank=data.frame(cbind(
  rowMeans(sapply(mase_lca_all_state_age, 
                  function(x){x$mase_state[10,]})),
  rowMeans(sapply(mase_fdm_all_state_age, 
                  function(x){x$mase_state[10,]})),
  rowMeans(sapply(mase_cfdm_all_state_age, 
                  function(x){x$mase_state[10,]})),
  rowMeans(sapply(mase_lgb_all_state_age, 
                  function(x){x$mase_state[10,]})),
  rowMeans(sapply(male_mase_expending_NOshrinkage, 
                  function(x){x$mase_state[10,]})),
  rowMeans(sapply(male_mase_expending_age_optlamb, 
                  function(x){x$mase_state[10,]})),
  rowMeans(sapply(male_mase_expending_state_optlamb, 
                  function(x){x$mase_state[10,]})),
  rowMeans(sapply(male_mase_expending_ss, 
                  function(x){x$mase_state[10,]}))
))
colnames(state_rank) = c("L-C","H-U","H-B-Y","LGB","GBLC","GBLC-age",
                         "GBLC-state","GBLC-age-state") 
head(state_rank)
#          L-C       H-U     H-B-Y       LGB      GBLC  GBLC-age GBLC-state GBLC-age-state
# AL 0.7850912 0.7432296 0.7646414 0.8907137 0.7666506 0.7201835  0.7627662      0.7208659
# AK 0.8661689 0.8395078 0.7391491 0.9110229 0.8793709 0.8149384  0.8783709      0.8149384
# AZ 0.5803605 0.5693170 0.5668184 0.6916979 0.6401783 0.6140098  0.6360093      0.6105017
# AR 0.7031229 0.7173536 0.7258102 0.8488125 0.7659479 0.7172188  0.7596606      0.7126269
# CA 0.7059773 0.6457623 0.8588497 0.7459960 0.5675464 0.5497963  0.5645826      0.5489028
# CO 0.6524175 0.6390560 0.6119814 0.6590579 0.6524722 0.6307633  0.6417565      0.6118072 

# Table B.2: Winkler score (×100) of out-of-sample forecasts 
# for state-level male mortality

#L-c
library(demography)
ws_lca = c()
for (j in 1:13) {
  year = dim(data_use)[1]
  matrix_motality_train = data_use_0[c(1:(year-10-j+1)),]
  matrix_motality_test = data_use[c((year-10-j+2):(year-j+1)),]
  pop_test = data_use_pop[year-10-j+1,]
  pop_train = data_use_pop[c(1:(year-10-j+1)),]

  fore_l = c()
  fore_u = c()
  for (s in 1:length(states_abb)) {
    states_col = which(stringr::str_detect(colnames(matrix_motality_train),states_abb[s]))
    matrix_motality_train_s = matrix_motality_train[,states_col]
    matrix_motality_test_s = matrix_motality_test[,states_col]
    pop_train_s = pop_train[,states_col]
    
    demogdata_motality = demogdata(data =t(matrix_motality_train_s),
                                   pop=t(pop_train_s),
                                   ages=c(0:85),
                                   years=c(1969:(1968+dim(matrix_motality_train_s)[1])),
                                   type="mortality", name="male",label="US")
    lca_motality = lca(demogdata_motality, max.age = 85,
                       adjust="none", chooseperiod = FALSE)
    fore = forecast(lca_motality,h=10,level = 95)
    forel = t(fore$rate$lower)
    foreu = t(fore$rate$upper)
    
    fore_l = cbind(fore_l,forel)
    fore_u = cbind(fore_u,foreu)
    
  }
  colnames(fore_l) = colnames(matrix_motality_train)
  colnames(fore_u) = colnames(matrix_motality_train)

  mase_s_all = c()
  for (i in 1:10) {
    mase_s = winkler_score (lt = matrix(fore_l[c(1:i),],nrow=i),
                            ut = matrix(fore_u[c(1:i),],nrow=i),
                            actual = matrix(matrix_motality_test[c(1:i),],nrow=i)
                            , level = 95
    )
    mase_s_all = c(mase_s_all, mase_s)
  }
  ws_lca = cbind(ws_lca,mase_s_all)
}
round(rowMeans(ws_lca)*100,3)
# [1] 2.546 2.444 2.407 2.408 2.427 2.457 2.495 2.539 2.592 2.651

#H-U (time-consuming!)
ws_lca = c()
for (j in 1:13) {
  year = dim(data_use)[1]
  matrix_motality_train = data_use_0[c(1:(year-10-j+1)),]
  matrix_motality_test = data_use[c((year-10-j+2):(year-j+1)),]
  pop_test = data_use_pop[year-10-j+1,]
  pop_train = data_use_pop[c(1:(year-10-j+1)),]
  fore_l = c()
  fore_u = c()
  for (s in 1:length(states_abb)) {
    states_col = which(stringr::str_detect(colnames(matrix_motality_train),states_abb[s]))
    matrix_motality_train_s = matrix_motality_train[,states_col]
    matrix_motality_test_s = matrix_motality_test[,states_col]
    pop_train_s = pop_train[,states_col]
    
    demogdata_motality = demogdata(data =t(matrix_motality_train_s),
                                   pop=t(pop_train_s),
                                   ages=c(0:85),
                                   years=c(1969:(1968+dim(matrix_motality_train_s)[1])),
                                   type="mortality", name="male",label="US")

    fdm_motality3 = fdm(demogdata_motality, max.age = 85,method="classical")
    fdm_motality3_forec = forecast(fdm_motality3,h=10,level = 95)
    forel = t(fdm_motality3_forec$rate$lower)
    foreu = t(fdm_motality3_forec$rate$upper)
    
    
    fore_l = cbind(fore_l,forel)
    fore_u = cbind(fore_u,foreu)
    
  }
  colnames(fore_l) = colnames(matrix_motality_train)
  colnames(fore_u) = colnames(matrix_motality_train)

  mase_s_all = c()
  for (i in 1:10) {
    mase_s = winkler_score (lt = matrix(fore_l[c(1:i),],nrow=i),
                            ut = matrix(fore_u[c(1:i),],nrow=i),
                            actual = matrix(matrix_motality_test[c(1:i),],nrow=i)
                            , level = 95
    )
    mase_s_all = c(mase_s_all, mase_s)
  }
  ws_lca = cbind(ws_lca,mase_s_all)
}
round(rowMeans(ws_lca)*100,3)

# H-B-Y (time-consuming!)
ws_lca = c()
for (j in 1:13) {
  year = dim(data_use)[1]
  matrix_motality_train = data_use_0[c(1:(year-10-j+1)),]
  matrix_motality_test = data_use[c((year-10-j+2):(year-j+1)),]
  pop_test = data_use_pop[year-10-j+1,]
  pop_train = data_use_pop[c(1:(year-10-j+1)),]

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
  fit.fdm = coherentfdm(demogdata_motality)
  forec.fdm = forecast(fit.fdm,h=10,level = 95)
  fore_l_list = lapply(forec.fdm[1:length(states_abb)],
                       function(x){t(x$rate$lower)})
  fore_u_list = lapply(forec.fdm[1:length(states_abb)],
                       function(x){t(x$rate$upper)})

  fore_l = c()
  for (s in 1:length(states_abb)) {
    fore_l = cbind(fore_l,fore_l_list[[s]])
  }
  colnames(fore_l) = colnames(matrix_motality_train)
  fore_u = c()
  for (s in 1:length(states_abb)) {
    fore_u = cbind(fore_u,fore_u_list[[s]])
  }
  colnames(fore_u) = colnames(matrix_motality_train)
  
  mase_s_all = c()
  for (i in 1:10) {
    mase_s = winkler_score (lt = matrix(fore_l[c(1:i),],nrow=i),
                            ut = matrix(fore_u[c(1:i),],nrow=i),
                            actual = matrix(matrix_motality_test[c(1:i),],nrow=i)
                            , level = 95
    )
    mase_s_all = c(mase_s_all, mase_s)
  }
  ws_lca = cbind(ws_lca,mase_s_all)
}
round(rowMeans(ws_lca)*100,3)



#LightGBM  (time-consuming, hours!)
param_best = list(learning_rate=0.05,
                  lambda_l2=0.03,
                  bagging_fraction=0.8,
                  feature_fraction=0.8,
                  num_leaves=30,
                  max_depth=4,
                  min_data_in_leaf=10,
                  objective = "quantile",
                  alpha = 0.025)

lgb_lower = list()
lgb_upper = list()
for (j in 1:13) {
  year = dim(data_use)[1]
  matrix_motality_train = data_use_0[c(1:(year-10-j+1)),]
  matrix_motality_test = data_use[c((year-10-j+2):(year-j+1)),]
  
  # motality_global_lgb: $data_train 
  # $data_forec : lags used to forecast, true Value is NA
  #  Value L1 L2 L3 L4 L5 L6 L7 L8 L9 L10 (for all ages)
  
  ff_lower = c()
  ff_upper = c()
  for (s in 1:length(states_abb)) {
    states_col = which(stringr::str_detect(colnames(matrix_motality_train),
                                           states_abb[s]))
    matrix_motality_train_s = matrix_motality_train[,states_col]
    data_train = c()
    data_forec = c()
    motality_global_lgb = list()
    for (a in 1:dim(matrix_motality_train_s)[2]) {
      ts_sa = list(x = ts(matrix_motality_train_s[,a]))
      data_local = mortality_data_generate_lgb(ts = ts_sa,
                                               num_lag=9)
      data_train = rbind(data_train,data_local$data_train)
      data_forec = rbind(data_forec,data_local$data_forec)
    }
    motality_global_lgb$data_train = log(data_train)
    motality_global_lgb$data_forec = log(data_forec)
    param_best$alpha = 0.025
    ff_l = mortality_lgb_func(param = param_best, nround = 100,
                              dataset.train = motality_global_lgb$data_train,
                              dataset.forec = motality_global_lgb$data_forec, 
                              h=10 )
    ff_l = exp(ff_l)
    ff_lower = cbind(ff_lower,ff_l)
    param_best$alpha = 0.975
    ff_u = mortality_lgb_func(param = param_best, nround = 100,
                              dataset.train = motality_global_lgb$data_train,
                              dataset.forec = motality_global_lgb$data_forec, 
                              h=10 )
    ff_u = exp(ff_u)
    ff_upper = cbind(ff_upper,ff_u)
  }
  lgb_lower[[j]] =ff_lower
  lgb_upper[[j]] =ff_upper
}


ws_lgb =c()
for (j in c(1:13)) {
  year = dim(data_use)[1]
  matrix_motality_train = data_use_0[c(1:(year-10-j+1)),]
  matrix_motality_test = data_use[c((year-10-j+2):(year-j+1)),]
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

#GBLC (time-consuming, hours!)
ws_lca = c()
for (j in 1:13) {
  year = dim(data_use)[1]
  matrix_motality_train = data_use_0[c(1:(year-10-j+1)),]
  matrix_motality_test = data_use[c((year-10-j+2):(year-j+1)),]
  set.seed(225)
  quantile = gb_forecast_pi_ss (matrix_input=matrix_motality_train,
                                matrix_W_1=W_comb_a,
                                matrix_W_2=W_comb_s,
                                lambda_1=0,
                                lambda_2=0,
                                h=10,
                                max.interation=3,
                                threshold = 1e-07,
                                state_include = TRUE,
                                states=states_abb,
                                alpha=c(0.025,0.975),
                                J=100)
  fore_u = quantile$final_quant[[2]]
  fore_l = quantile$final_quant[[1]]
  
  colnames(fore_l) = colnames(matrix_motality_train)
  colnames(fore_u) = colnames(matrix_motality_train)
  
  mase_s_all = c()
  for (i in 1:10) {
    mase_s = winkler_score (lt = matrix(fore_l[c(1:i),],nrow=i),
                            ut = matrix(fore_u[c(1:i),],nrow=i),
                            actual = matrix(matrix_motality_test[c(1:i),],nrow=i)
                            , level = 95
    )
    mase_s_all = c(mase_s_all, mase_s)
  }
  ws_lca = cbind(ws_lca,mase_s_all)
}
round(rowMeans(ws_lca)*100,3)



# Table C.2: MASE of out-of-sample forecasts for all age groups, h=10, state-level males

age_rank=data.frame(cbind(
  rowMeans(sapply(mase_lca_all_state_age, 
                  function(x){x$mase_age[10,]})),
  rowMeans(sapply(mase_fdm_all_state_age, 
                  function(x){x$mase_age[10,]})),
  rowMeans(sapply(mase_cfdm_all_state_age, 
                  function(x){x$mase_age[10,]})),
  rowMeans(sapply(mase_lgb_all_state_age, 
                  function(x){x$mase_age[10,]})),
  rowMeans(sapply(male_mase_expending_NOshrinkage, 
                  function(x){x$mase_age[10,]})),
  rowMeans(sapply(male_mase_expending_age_optlamb, 
                  function(x){x$mase_age[10,]})),
  rowMeans(sapply(male_mase_expending_state_optlamb, 
                  function(x){x$mase_age[10,]})),
  rowMeans(sapply(male_mase_expending_ss, 
                  function(x){x$mase_age[10,]}))
)
)
colnames(age_rank) = c("L-C","H-U","H-B-Y","LightGBM","GBLC","GBLC-age",
                       "GBLC-state","GBLC-age-state") 
head(age_rank)
#         L-C       H-U     H-B-Y  LightGBM      GBLC  GBLC-age GBLC-state GBLC-age-state
# 0 0.9124422 0.7008736 0.6919630 0.6047085 0.7877963 0.7707243  0.7830493      0.7642860
# 1 0.6507482 0.6490994 0.6626906 0.8016105 0.7561388 0.7599935  0.7510382      0.7525488
# 2 0.7556295 0.7547692 0.7643534 0.7399812 0.8762575 0.8288598  0.8704025      0.8229246
# 3 0.8162862 0.8287962 0.8440232 0.8041800 0.9133002 0.8800138  0.9082634      0.8661177
# 4 0.8124087 0.8169980 0.8520376 0.8541099 0.8840611 0.8612331  0.8781560      0.8515777
# 5 0.7749375 0.7824588 0.8269271 0.8878179 0.8406141 0.8340131  0.8357962      0.8274833





