library(rmgarch)
library(vars)
library(dplyr)
library(aTSA)
library(tseries)
library(forecast)
library(FinTS)

gfc_from<-"2007-06-29"
post_gfc_from<-"2009-02-06"
low_oil_from<-"2014-11-30"
fed_rates_from<-"2017-08-10"
covid_from<-"2020-01-01"
gfc_to<-"2009-02-05"
post_gfc_to<-"2014-11-29"
low_oil_to<-"2017-08-09"
fed_rates_to<-"2019-12-31"
covid_to<-"2022-01-01"
ranges<-data.frame(cbind(c(gfc_from,post_gfc_from,low_oil_from,fed_rates_from,covid_from),
                         c(gfc_to,post_gfc_to,low_oil_to,fed_rates_to,covid_to)))
colnames(ranges)<-c("from","to")

best_garch<-function(fit){
  p_max <- 5
  q_max <- 5
  aic_min <- Inf
  best_p <- 0
  best_q <- 0
  x <- fit$x %>%as.vector()
  ar <- ifelse(length(fit_auto$coef)>0,fit_auto$arma[2],0)
  ma<- ifelse(length(fit_auto$coef)>0,fit_auto$arma[1],0)
  for (i1 in 1:p_max) {
    for (i2 in 1:q_max) {
      model_specification <- ugarchspec(mean.model = list(armaOrder = c(ar,ma), include.mean = TRUE),
                                        variance.model = list(garchOrder = c(i1, i2)))
      fit <- ugarchfit(spec = model_specification, data = x)
      inf_crit <- infocriteria(fit)[1]
      aic_min <- ifelse(inf_crit < aic_min, inf_crit, aic_min)
      best_p <- ifelse(inf_crit == aic_min, i1, best_p)
      best_q <- ifelse(inf_crit == aic_min, i2, best_q)
    }
  }

  return(c(best_p,best_q))
}

adcc_output<-function(mat_list,data=data,
                      statistics = c("mean","median","sd","min","max"),
                      type = "weight"){
  list_subset<-mat_list
  z <- array(unlist(list_subset), c(dim(list_subset[[1]]), length(list_subset)))
  colnames(z)<-colnames(data)
  rownames(z)<-colnames(data)
  temp<-data.frame(matrix(0,ncol=length(statistics)))
  for(stat in statistics){
    fun_mat<-apply(z, 1:2, stat)
    fun_arr<-mat_to_long(fun_mat,type=type,stat=stat)
    if(dim(temp)[1]==1){
      temp<-fun_arr
    } else{
      temp<-left_join(temp,fun_arr,by="name")
    }
  }
  final <- temp[,-1]
  rownames(final) <- temp[,1]
  return(final)
}
mat_to_long<-function(mat,type="weight",stat="mean"){
  vals<-as.vector(mat)
  df<-data.frame(matrix(ncol=3))
  counter<-1
  for(j in 1:ncol(mat)){
    for (i in 1:nrow(mat)){
      nm1<-rownames(mat)[j]
      nm2<-colnames(mat)[i]
      if(i != j){
        df[counter,]<-c(nm1,nm2,mat[i,j])
        counter=counter+1
      }
    }
  }
  clean_df<-df
  if (type=="weight"){
    cols = c(1,2)
    zz = df[,cols]
    for (i in 1:nrow(df)){
      zz[i, ] = sort(df[i,cols])
    }
    clean_df<-df[!duplicated(zz),]
    rownames(clean_df) <- NULL
  }
  clean_df<-cbind(paste0(clean_df[,1], " ", clean_df[,2]),clean_df[,3])
  clean_df<-as.data.frame(clean_df)
  clean_df[,2]<-as.numeric(clean_df[,2])
  colnames(clean_df)<-c("name",stat)
  return(clean_df)
}
weightMat<-function(table,dt){
  n<-length(dt)
  nms<-colnames(dt)
  weight_mat<-data.frame(matrix(ncol = n, nrow = n))
  colnames(weight_mat)<-nms
  rownames(weight_mat)<-nms

  for (i in 1:nrow(weight_mat)){
    for (j in 1:length(weight_mat)){
      h_jj<-table[j,j]
      h_ij<-table[i,j]
      h_ii<-table[i,i]
      w_ij<-(h_jj-h_ij)/(h_ii-2*h_ij+h_jj)
      weight_mat[i,j]<-w_ij
    }
  }
  return(weight_mat)
}

hedgeMat<-function(table,dt){
  n<-length(dt)
  nms<-colnames(dt)
  hedge_mat<-data.frame(matrix(ncol = n, nrow = n))
  colnames(hedge_mat)<-nms
  rownames(hedge_mat)<-nms

  for (i in 1:nrow(hedge_mat)){
    for (j in 1:length(hedge_mat)){
      h_jj<-table[j,j]
      h_ij<-table[i,j]
      h_ii<-table[i,i]
      hedge_ij<-h_ij/h_jj
      hedge_mat[i,j]<-hedge_ij
    }
  }
  return(hedge_mat)
}


log_ret_df<-readRDS("data/log_ret.rds")
dummy_thursday<-readRDS("data/dummy_log_ret.rds")
dummy_thursday <-dummy_thursday %>% dplyr::filter(rownames(dummy_thursday)>=covid_from,rownames(dummy_thursday)<=covid_to) %>%
  select(Weekday_Thursday) %>% as.matrix()
data<-log_ret_df[-1] %>%as.data.frame() %>% select(Crude_Oil_WTI,Diesel,Heating_Oil,Natural_Gas,NY_Gasoline)
dates<-log_ret_df[1]%>%pull()
rownames(data)<-dates

garch_0011<-ugarchspec(mean.model = list(armaOrder = c(0,0)),
                           variance.model = list(garchOrder=c(1,1)))
garch_1011<-ugarchspec(mean.model = list(armaOrder = c(1,0)),
                       variance.model = list(garchOrder=c(1,1)))
garch_0021<-ugarchspec(mean.model = list(armaOrder = c(0,0)),
                       variance.model = list(garchOrder=c(2,1)))
garch_0031<-ugarchspec(mean.model = list(armaOrder = c(0,0)),
                       variance.model = list(garchOrder=c(3,1)))
garch_0041<-ugarchspec(mean.model = list(armaOrder = c(0,0)),
                       variance.model = list(garchOrder=c(4,1)))
garch_0022<-ugarchspec(mean.model = list(armaOrder = c(0,0)),
                       variance.model = list(garchOrder=c(2,2)))
garch_0011_ged<-ugarchspec(mean.model = list(armaOrder = c(0,0)),
                  variance.model = list(garchOrder=c(1,1)),
                  distribution.model="ged")
garch_2021<-ugarchspec(mean.model = list(armaOrder = c(2,0)),
                       variance.model = list(garchOrder=c(2,1)))
garch_0051_thursday<-ugarchspec(
  mean.model = list(armaOrder = c(0,0),external.regressors=dummy_thursday),
  variance.model = list(garchOrder=c(5,1))
  )
garch_1041<-ugarchspec(mean.model = list(armaOrder = c(1,0)),
                       variance.model = list(garchOrder=c(4,1)))
garch_1021<-ugarchspec(mean.model = list(armaOrder = c(1,0)),
                       variance.model = list(garchOrder=c(4,1)))

## POST-GFC


### AR and ARCH test
test <- arima.sim(model = list(ar = c(0.75)), n = 2000)
test <- arima.sim(n = 2000, list(ma=0.8), innov=rnorm(100))
test<-gfc_data[1] %>% as.matrix()
test<-low_oil_data[5] %>% as.matrix()
test<-fed_rates_data[4]  %>% as.matrix()
fit_auto <- auto.arima(test, ic = c("bic"), stationary = TRUE)
fit_auto<-arima(test, order =c(5,0,2))
checkresiduals(fit_auto)
k<-Box.test(fit_auto$residuals, lag = 20, type = c("Ljung-Box"))
l<-ArchTest(fit_auto$residuals, lags = 12, demean = FALSE)



## FIN KRIZE
gfc_data <-data %>% dplyr::filter(rownames(data)>=gfc_from,rownames(data)<=gfc_to)
gfc_spec<-multispec(c(garch_0011,garch_0011_ged,garch_0011_ged,garch_0011,garch_0011))
dcc.11mn = dccspec(uspec=gfc_spec, VAR = FALSE, lag = 1,
                   dccOrder = c(1, 1), model="aDCC",
                   distribution = 'mvnorm')

oil<-ugarchfit(garch_0011,gfc_data[1])
diesel<-ugarchfit(garch_0011_ged,gfc_data[2])
heating_oil<-ugarchfit(garch_0011,gfc_data[3])
natural_gas<-ugarchfit(garch_0011,gfc_data[4])
ny_gas<-ugarchfit(garch_0011,gfc_data[5])
# gfc_3<-multispec(c(garch_0011,garch_0011_ged,garch_0011_ged))
# gfc_12<-multispec(c(garch_0011,garch_0011_ged))
# gfc_23<-multispec(c(garch_0011_ged,garch_0011_ged))
# gfc_13<-multispec(c(garch_0011,garch_0011_ged))
#
# dcc3<-dccspec(uspec=gfc_3, model="aDCC",distribution = 'mvnorm')
# dcc12<-dccspec(uspec=gfc_12, model="aDCC",distribution = 'mvnorm')
# dcc23<-dccspec(uspec=gfc_23, model="aDCC",distribution = 'mvnorm')
# dcc13<-dccspec(uspec=gfc_13, model="aDCC",distribution = 'mvnorm')
#
# gfc_3_fit<-dccfit(dcc.11mn,data=gfc_data[c(1,2,3)])
# gfc_12_fit<-dccfit(dcc12,data=gfc_data[c(1,2)])
# gfc_23_fit<-dccfit(dcc23,data=gfc_data[c(2,3)])
# gfc_13_fit<-dccfit(dcc13,data=gfc_data[c(1,3)])
#
#
# gfc_3_fit@model$pars[1:3,1]+
# gfc_12_fit@model$pars[1:3,1]+
# gfc_23_fit@model$pars[1:3,1]
# gfc_13_fit@model$pars[1:3,1]


gfc_fit<-dccfit(dcc.11mn,data=gfc_data)
gfc_fit@mfit$
dcc_cov<-rcov(gfc_fit)
dcc_cov_list<-lapply(seq(dim(dcc_cov)[3]), function(x) dcc_cov[ , , x])
hedge_list<-lapply(dcc_cov_list,hedgeMat,dt=gfc_data)
weight_list<-lapply(dcc_cov_list,weightMat,dt=gfc_data)
name_dates<-as.Date(dimnames(dcc_cov)[3][[1]])
names(hedge_list)<-name_dates
names(weight_list)<-name_dates
gfc_hedge<-adcc_output(hedge_list,data=gfc_data,type="hedge") %>% cbind(period = "gfc")
gfc_weight<-adcc_output(weight_list,data=gfc_data,type="weight") %>% cbind(period = "gfc")


## Post GFC
post_gfc_data <-data %>% dplyr::filter(rownames(data)>=post_gfc_from,rownames(data)<=post_gfc_to)
post_gfc_spec<-multispec(c(garch_0011,garch_0011,garch_0022,garch_2021,garch_0031))
dcc.11mn = dccspec(uspec=post_gfc_spec,
                   VAR = FALSE, lag = 1,
                   dccOrder = c(1, 1), model="aDCC",
                   distribution = 'mvnorm')
post_gfc_fit<-dccfit(dcc.11mn,data=post_gfc_data)
dcc_cov<-rcov(post_gfc_fit)
dcc_cov_list<-lapply(seq(dim(dcc_cov)[3]), function(x) dcc_cov[ , , x])
hedge_list<-lapply(dcc_cov_list,hedgeMat,dt=post_gfc_data)
weight_list<-lapply(dcc_cov_list,weightMat,dt=post_gfc_data)
name_dates<-as.Date(dimnames(dcc_cov)[3][[1]])
names(hedge_list)<-name_dates
names(weight_list)<-name_dates
post_gfc_hedge<-adcc_output(hedge_list,data=post_gfc_data,type="hedge") %>% cbind(period = "post_gfc")
post_gfc_weight<-adcc_output(weight_list,data=post_gfc_data,type="weight") %>% cbind(period = "post_gfc")
g<-multifit(post_gfc_spec,gfc_data)
## Low_oil
low_oil_data <-data %>% dplyr::filter(rownames(data)>=low_oil_from,rownames(data)<=low_oil_to)
low_oil_spec<-multispec(c(garch_0011,garch_0011,garch_0041,garch_0011,garch_0021))
dcc.11mn = dccspec(uspec=low_oil_spec,
                   VAR = FALSE, lag = 1,
                   dccOrder = c(1, 1), model="aDCC",
                   distribution = 'mvnorm')
low_oil_fit<-dccfit(dcc.11mn,data=low_oil_data)
dcc_cov<-rcov(low_oil_fit)
dcc_cov_list<-lapply(seq(dim(dcc_cov)[3]), function(x) dcc_cov[ , , x])
hedge_list<-lapply(dcc_cov_list,hedgeMat,dt=low_oil_data)
weight_list<-lapply(dcc_cov_list,weightMat,dt=low_oil_data)
name_dates<-as.Date(dimnames(dcc_cov)[3][[1]])
names(hedge_list)<-name_dates
names(weight_list)<-name_dates
low_oil_hedge<-adcc_output(hedge_list,data=low_oil_data,type="hedge") %>% cbind(period = "low_oil")
low_oil_weight<-adcc_output(weight_list,data=low_oil_data,type="weight") %>% cbind(period = "low_oil")


## fed_rates
fed_rates_data <-data %>% dplyr::filter(rownames(data)>=fed_rates_from,rownames(data)<=fed_rates_to)
fed_rates_spec<-multispec(c(garch_0011,garch_0011,garch_0011,garch_0011,garch_0021))
dcc.11mn = dccspec(uspec=fed_rates_spec,
                   VAR = FALSE, lag = 1,
                   dccOrder = c(1, 1), model="aDCC",
                   distribution = 'mvnorm')
fed_rates_fit<-dccfit(dcc.11mn,data=fed_rates_data)
dcc_cov<-rcov(fed_rates_fit)
dcc_cov_list<-lapply(seq(dim(dcc_cov)[3]), function(x) dcc_cov[ , , x])
hedge_list<-lapply(dcc_cov_list,hedgeMat,dt=fed_rates_data)
weight_list<-lapply(dcc_cov_list,weightMat,dt=fed_rates_data)
name_dates<-as.Date(dimnames(dcc_cov)[3][[1]])
names(hedge_list)<-name_dates
names(weight_list)<-name_dates
fed_rates_hedge<-adcc_output(hedge_list,data=fed_rates_data,type="hedge") %>% cbind(period = "fed_rates")
fed_rates_weight<-adcc_output(weight_list,data=fed_rates_data,type="weight") %>% cbind(period = "fed_rates")


## covid
covid_data <-data %>% dplyr::filter(rownames(data)>=covid_from,rownames(data)<=covid_to)
covid_spec<-multispec(c(garch_0051_thursday,garch_0021,garch_0021,garch_1041,garch_1021))
covid_dates<-rownames(covid_data)

p_max <- 7
q_max <- 7
aic_min <- Inf
best_p <- 0
best_q <- 0

for (i1 in 1:p_max) {
  for (i2 in 1:q_max) {
    model_specification <- ugarchspec(mean.model = list(armaOrder = c(0, 0), include.mean = TRUE),
                                      variance.model = list(garchOrder = c(i1, i2)))
    fit <- ugarchfit(spec = model_specification, data = covid_data[1])
    inf_crit <- infocriteria(fit)[1]
    aic_min <- ifelse(inf_crit < aic_min, inf_crit, aic_min)

    best_p <- ifelse(inf_crit == aic_min, i1, best_p)
    best_q <- ifelse(inf_crit == aic_min, i2, best_q)
  }
}

oil<-ugarchfit(garch_0051_thursday,covid_data[1])
diesel<-ugarchfit(garch_0021,covid_data[2])
heating_oil<-ugarchfit(garch_0021,covid_data[3])
natural_gas<-ugarchfit(garch_1041,covid_data[4])
ny_gas<-ugarchfit(garch_1021,covid_data[5])
# library(garchmodels)
# garch_reg(arch_order = 5,garch_order = 1,ar_order = 0,ma_order = 0) %>%
#   set_engine("rugarch", mean.model = list(external.regressors = dummy_thursday)) %>%
#   fit(Crude_Oil_WTI~covid_dates+Weekday_Thursday,data = covid_dt)
# library(fGarch)
# spec = garchSpec(model = list(alpha = c(0.1, 0.1,0.1,0.1,0.1), beta = 0.08))
# garch_0051<-ugarchspec(
#   mean.model = list(armaOrder = c(1,1)),
#   variance.model = list(garchOrder=c(5,1))
# )
# k<-ugarchfit(garch_0051,covid_data[1])
# garch(covid_data[1],order = c(5,1))
# fgarch_vec<-garchSpec(model = list(alpha = c(0.01,0.01,0.01, 0.01,0.01), beta = 0.08))
# z<-garchFit(~garch(5,1),data=covid_data[1])
dcc.11mn = dccspec(uspec=covid_spec,
                   VAR = FALSE, lag = 1,
                   dccOrder = c(1, 1), model="aDCC",
                   distribution = 'mvnorm')
covid_fit<-dccfit(dcc.11mn,data=covid_data)
covid_fit
dcc_cov<-rcov(covid_fit)
dcc_cov_list<-lapply(seq(dim(dcc_cov)[3]), function(x) dcc_cov[ , , x])
hedge_list<-lapply(dcc_cov_list,hedgeMat,dt=covid_data)
weight_list<-lapply(dcc_cov_list,weightMat,dt=covid_data)
name_dates<-as.Date(dimnames(dcc_cov)[3][[1]])
names(hedge_list)<-name_dates
names(weight_list)<-name_dates
covid_hedge<-adcc_output(hedge_list,data=covid_data,type="hedge") %>% cbind(period = "covid")
covid_weight<-adcc_output(weight_list,data=covid_data,type="weight") %>% cbind(period = "covid")

## summary
fits<-list(gfc_fit,post_gfc_fit,low_oil_fit,fed_rates_fit,covid_fit)
names(fits)<-c("gfc_fit","post_gfc_fit","low_oil_fit","fed_rates_fit","covid_fit")
hedge<-rbind(gfc_hedge,post_gfc_hedge,low_oil_hedge,fed_rates_hedge,covid_hedge)
weight<-rbind(gfc_weight,post_gfc_weight,low_oil_weight,fed_rates_weight,covid_weight)
saveRDS(fits,"data/adcc_fits.rds")
write.csv(hedge,"data/hedge.csv")
write.csv(weight,"data/weight.csv")
final <- matrix(nrow=0,ncol=6)
for(i in 1:length(fits)){
  z<-apply(fits[[i]]@mfit$H,1:2,mean)
  colnames(z)<-colnames(data)
  rownames(z)<-colnames(data)
  z<-z%>%as.data.frame()
  z[1:5,6]<-names(fits)[i]
  final<-rbind(final,z)
}

write.csv(final,"data/H_matrices.csv")

