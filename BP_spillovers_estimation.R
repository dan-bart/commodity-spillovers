## ----packages, echo = F----------------------------------------------------------------------------------------------------------------------------------
library(ggplot2)
library(readxl)
library(ggpubr)
library(kableExtra)
library(stargazer)
library(tseries)
library(lmtest)
library(stringr)
library(frequencyConnectedness)
library(imputeTS)
library(ggcorrplot)
library(lubridate)
library(dplyr)
library(reshape2)
library(strucchange)
library(parallel)
library(BigVAR)
library(scales)
library(gridExtra)


## ----load_data, include =F-------------------------------------------------------------------------------------------------------------------------------
log_ret<-readRDS("data/log_ret.rds")
df<-readRDS("data/df.rds")
data<-log_ret[-1] %>%as.data.frame()
dates<-log_ret[1]%>%pull()
rownames(data)<-dates


gfc_from<-"2007-06-29"
post_gfc_from<-"2009-02-06"
low_oil_from<-"2014-11-30"
fed_rates_from<-"2017-08-10"
covid_from<-"2020-01-01"
gfc_to<-"2009-02-05"
post_gfc_to<-"2011-11-29"
low_oil_to<-"2017-08-09"
fed_rates_to<-"2019-12-31"
covid_to<-"2022-01-01"
#----struct_change, eval = F, include = F----------------------------------------------------------------------------------------------------------------
chow_list<-list()
for(i in 5:(length(df))){
  col<-colnames(df)[i]
  print(col)
  val<-df %>% dplyr::select(i)%>%pull()
  ##Chow Test: H0: no structural break in the data
  chow<-sctest(I(val)~1,type="Chow")
  chow_list<-c(chow_list,list(chow))
  names(chow_list)[[i-3]]<-col

  # ## Bai Perron
  bp_name<-paste("data/struc_change/bp_prices_v2","_", col,".rds",sep="")
  bp <- breakpoints(I(val) ~ I(lag(val)),h=600)
  saveRDS(bp,file=bp_name)
}

saveRDS(chow_list,file="data/struc_change/chow.rds")
input<-colnames(log_ret_df)[-c(1)]
#save Plots and Breakpoints
bp_plots<-list()
bp_breaks<-list()
for (i in 1:length(input)){
  bp_name<-paste("data/struc_change/bp_prices_v2","_", input[i],".rds",sep="")
  print(bp_name)
  bp<-readRDS(bp_name)

  bp_plot<-plot(bp)
  sum_bp<-summary(bp)
  sum_bp_df<-sum_bp$breakpoints%>%as.data.frame()
  bp_break_dates<-sapply(sum_bp_df,function(x) char_dates[x]) %>% as.data.frame()

  bp_plots<-c(bp_plots,list(bp_plot))
  bp_breaks<-c(bp_breaks,list(bp_break_dates))

  names(bp_plots)[[i]]<-input[i]
  names(bp_breaks)[[i]]<-input[i]
}

saveRDS(bp_plots,file="data/struc_change/bp_plots_v2.rds")
saveRDS(bp_breaks,file="data/struc_change/break_dates_v2.rds")



##STATIC

est_lag1<-vars::VAR(log_ret_df[-1], p = 1, type = "const")
est_lag2<-vars::VAR(log_ret_df[-1], p = 2, type = "const")
est_lag4<-vars::VAR(log_ret_df[-1], p = 4, type = "const")
sp_lag1<-spilloverDY12(est_lag1, n.ahead = 100, no.corr = F)
sp_lag2<-spilloverDY12(est_lag2, n.ahead = 100, no.corr = F)
sp_lag4<-spilloverDY12(est_lag4, n.ahead = 100, no.corr = F)
saveRDS(sp_lag1,"data/sp_data/total/sp_lag1.rds")
saveRDS(sp_lag2,"data/sp_data/total/sp_lag2.rds")
saveRDS(sp_lag4,"data/sp_data/total/sp_lag4.rds")


o <- big_var_est(volatilities)
sp_lasso_lag1<-spilloverDY12(o, n.ahead = 100, no.corr = F)
saveRDS(sp_lasso_lag1,"data/sp_data/total/sp_lasso_lag1.rds")

##ROLLING
dates<-log_ret_df %>% dplyr::select(1)%>%pull()
## ----DY_est, include = F, eval = F-----------------------------------------------------------------------------------------------------------------------
cl <- makeCluster(8)
sp <- spilloverRollingDY12(data, n.ahead = 100, no.corr = F, func_est = "VAR",
                           params_est = list(p = 1,type = "const"),
                           window = 100,cluster  = cl)
plotOverall(sp)
nm<-paste0("sp_w",100,"_lag1")
fullnm<-paste0("data/sp_data/",nm,".rds")
saveRDS(sp,fullnm)
stopCluster(cl)

sp_ls<-list.files("./data/sp_data")
sp_ls<-sp_ls[-1]
sp_ls<-sp_ls[str_detect(sp_ls,"sp")]
sp_ls<-sp_ls[!str_detect(sp_ls,"lasso")]
sp_ls<-sp_ls[!str_detect(sp_ls,"freq")]
for (i in sp_ls){
  print(i)
  sp<-readRDS(paste0("data/sp_data/",i))
  space<-length(dates)-length(sp[[1]])+1
  print(space)
  dates_temp<-dates[space:length(dates)]
  for (j in 1:length(sp[[1]])){
    sp[[1]][[j]]$date<-as.POSIXct(dates_temp[j])
  }

  saveRDS(sp,paste0("data/sp_data/",i))
  sp_nm<-substr(i,3,nchar(i))
  plot_nm<-paste0("data/sp_data/plots/sp_plot",sp_nm)
  print(plot_nm)
  ov_plot<-plotOverall(sp)
  saveRDS(ov_plot,plot_nm)
}


## ----rolling_freq_est, eval=F,include=F------------------------------------------------------------------------------------------------------------------
freq_con<-function(from="1996-01-01",to="2022-01-01",name="ALL"){
  dates<-log_ret%>% filter(Date>=from,Date<=to) %>% dplyr::select(1)%>%pull()
  volatilities<-log_ret%>% filter(Date>=gfc_from,Date<=gfc_to) %>% select(-Date)
  #data("volatilities")
  #volatilities<-volatilities %>% imputeTS::na_ma(k = 1, weighting = "simple")
  cl <- makeCluster(8)
  sp_roll_BK12 <- spilloverRollingBK12(volatilities, n.ahead = 100, no.corr = F, func_est = "VAR",
                                       params_est= list(p = 1, type = "const"), window = 100,
                                       partition = c(pi+0.00001, pi/5, pi/20, pi/200, 0),cluster  = cl)
  stopCluster(cl)
  space<-length(dates)-length(sp_roll_BK12[[1]])+1
  print(space)

  dates_temp<-dates[space:length(dates)]
  for (j in 1:length(sp_roll_BK12[[1]])){
    sp_roll_BK12[[1]][[j]]$date<-as.POSIXct(dates_temp[j])
  }
  plotOverall(sp_roll_BK12)
  saveRDS(sp_roll_BK12,paste0("data/sp_data/sp_roll_freq_lag1_",name,".rds"))
}

freq_con(gfc_from,gfc_to,"gfc")

## ----big_var_est, eval=F,include=F-----------------------------------------------------------------------------------------------------------------------
volatilities<-log_ret_df[-1]
big_var_est <- function(data) {
    Model1 = BigVAR::constructModel(as.matrix(data), p = 1, struct = "Basic", gran = c(50, 50), VARX = list(), verbose = F)
    Model1Results = BigVAR::cv.BigVAR(Model1)
}


# Perform the estimation

cl <- makeCluster(8)
sp <- spilloverRollingDY12(volatilities, n.ahead = 100, no.corr = F, func_est = "big_var_est",
                           params_est = list(),
                           window = 100,cluster  = cl)
plotOverall(sp)
nm<-paste0("sp_lasso_w",100,"_lag1")
fullnm<-paste0("data/sp_data/",nm,".rds")
saveRDS(sp,fullnm)
stopCluster(cl)

sp_ls<-list.files("./data/sp_data")
sp_ls<-sp_ls[str_detect(sp_ls,"sp")]
for (i in sp_ls){
  sp<-readRDS(paste0("data/sp_data/",i))
  space<-length(dates)-length(sp[[1]])+1
  print(space)
  dates_temp<-dates[space:length(dates)]
  for (j in 1:length(sp[[1]])){
    sp[[1]][[j]]$date<-as.POSIXct(dates_temp[j])
  }

  saveRDS(sp,paste0("data/sp_data/",i))
  sp_nm<-substr(i,3,nchar(i))
  plot_nm<-paste0("data/sp_data/plots/sp_plot",sp_nm)
  ov_plot<-plotOverall(sp)
  saveRDS(ov_plot,plot_nm)
}

sp<-readRDS(paste0("data/sp_data/sp_w100_lag1.rds"))
space<-length(dates)-length(sp[[1]])+1
