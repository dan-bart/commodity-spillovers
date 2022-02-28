library(rmgarch)
library(vars)
library(dplyr)
library(aTSA)
library(tseries)
library(forecast)


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

log_ret_df<-readRDS("data/log_ret.rds")
dummy_thursday<-readRDS("data/dummy_log_ret.rds")
dummy_thursday <-dummy_thursday %>% dplyr::filter(rownames(dummy_thursday)>=covid_from,rownames(dummy_thursday)<=covid_to) %>%
  select(Weekday_Thursday) %>% as.matrix()
data<-log_ret_df[-1] %>%as.data.frame() %>% select(Crude_Oil_WTI,Diesel,Heating_Oil,Natural_Gas,NY_Gasoline)
dates<-log_ret_df[1]%>%pull()
rownames(data)<-dates


## Applying proven GARCH settings from Chan (2016)

# POST-GFC
post_gfc_data <-data %>% dplyr::filter(rownames(data)>=post_gfc_from,rownames(data)<=post_gfc_to)

## NY Gasoline - GARCH-GJR and GARCH-M
garch_gjr<-ugarchspec(mean.model = list(armaOrder = c(0,0)),
                           variance.model = list(model = "gjrGARCH",garchOrder=c(1,1)))
garch_11<-ugarchspec(mean.model = list(armaOrder = c(0,0)),
                      variance.model = list(model = "sGARCH",garchOrder=c(1,1)))
garch_ar2<-ugarchspec(mean.model = list(armaOrder = c(2,0)),
                      variance.model = list(model = "sGARCH",garchOrder=c(1,1)))
garch_m<-ugarchspec(mean.model = list(armaOrder = c(0,0),archm=TRUE,archpow=2),
                      variance.model = list(model = "sGARCH",garchOrder=c(1,1)))
garch_25<-ugarchspec(mean.model = list(armaOrder = c(0,0)),
                     variance.model = list(model = "sGARCH",garchOrder=c(2,5)))
garch_2155<-ugarchspec(mean.model = list(armaOrder = c(2,1)),
                     variance.model = list(model = "sGARCH",garchOrder=c(5,5)))
garch_13<-ugarchspec(mean.model = list(armaOrder = c(0,0)),
                     variance.model = list(model = "sGARCH",garchOrder=c(1,3)))

ny_gasoline_own<-ugarchfit(garch_25,post_gfc_data[5])
resid_norm <- residuals(ny_gasoline_own, standardize = TRUE)
ArchTest(resid_norm, lags = 12, demean = FALSE)
resid_sqrd_norm <- resid_norm^2 # create squared residuals
par(mfrow = c(1,2))
forecast::Acf(resid_sqrd_norm, lag.max = 50, main = "ACF")
forecast::Pacf(resid_sqrd_norm, lag.max = 50, main = "PACF")

ny_gasoline_gjr<-ugarchfit(garch_gjr,post_gfc_data[5])
resid_norm <- residuals(ny_gasoline_gjr, standardize = TRUE)
ArchTest(resid_norm, lags = 12, demean = FALSE)
resid_sqrd_norm <- resid_norm^2 # create squared residuals
par(mfrow = c(1,2))
forecast::Acf(resid_sqrd_norm, lag.max = 50, main = "ACF")
forecast::Pacf(resid_sqrd_norm, lag.max = 50, main = "PACF")

ny_gasoline_m<-ugarchfit(garch_m,post_gfc_data[5])
resid_norm <- residuals(ny_gasoline_m, standardize = TRUE)
ArchTest(resid_norm, lags = 12, demean = FALSE)
resid_sqrd_norm <- resid_norm^2 # create squared residuals
par(mfrow = c(1,2))
forecast::Acf(resid_sqrd_norm, lag.max = 50, main = "ACF")
forecast::Pacf(resid_sqrd_norm, lag.max = 50, main = "PACF")

infocriteria(ny_gasoline_own)
infocriteria(ny_gasoline_gjr)
infocriteria(ny_gasoline_m)

## HEATING OIL- GARCH-GJR and GARCH-M
heating_oil_own<-ugarchfit(garch_13,post_gfc_data[3])
resid_norm <- residuals(heating_oil_own, standardize = TRUE)
ArchTest(resid_norm, lags = 12, demean = FALSE)
resid_sqrd_norm <- resid_norm^2 # create squared residuals
par(mfrow = c(1,2))
forecast::Acf(resid_sqrd_norm, lag.max = 50, main = "ACF")
forecast::Pacf(resid_sqrd_norm, lag.max = 50, main = "PACF")

heating_oil_gjr<-ugarchfit(garch_gjr,post_gfc_data[3])
resid_norm <- residuals(heating_oil_gjr, standardize = TRUE)
ArchTest(resid_norm, lags = 12, demean = FALSE)
resid_sqrd_norm <- resid_norm^2 # create squared residuals
par(mfrow = c(1,2))
forecast::Acf(resid_sqrd_norm, lag.max = 50, main = "ACF")
forecast::Pacf(resid_sqrd_norm, lag.max = 50, main = "PACF")

heating_oil_m<-ugarchfit(garch_m,post_gfc_data[3])
resid_norm <- residuals(heating_oil_m, standardize = TRUE)
ArchTest(resid_norm, lags = 12, demean = FALSE)
resid_sqrd_norm <- resid_norm^2 # create squared residuals
par(mfrow = c(1,2))
forecast::Acf(resid_sqrd_norm, lag.max = 50, main = "ACF")
forecast::Pacf(resid_sqrd_norm, lag.max = 50, main = "PACF")

infocriteria(heating_oil_own)
infocriteria(heating_oil_gjr)
infocriteria(heating_oil_m)

## NATURAL GAS - garch normal, garch ar2
nat_gas_own<-ugarchfit(garch_2155,post_gfc_data[4])
resid_norm <- residuals(nat_gas_own, standardize = TRUE)
ArchTest(resid_norm, lags = 12, demean = FALSE)
resid_sqrd_norm <- resid_norm^2 # create squared residuals
par(mfrow = c(1,2))
forecast::Acf(resid_sqrd_norm, lag.max = 50, main = "ACF")
forecast::Pacf(resid_sqrd_norm, lag.max = 50, main = "PACF")

nat_gas_11<-ugarchfit(garch_11,post_gfc_data[4])
resid_norm <- residuals(nat_gas_11, standardize = TRUE)
ArchTest(resid_norm, lags = 12, demean = FALSE)
resid_sqrd_norm <- resid_norm^2 # create squared residuals
par(mfrow = c(1,2))
forecast::Acf(resid_sqrd_norm, lag.max = 50, main = "ACF")
forecast::Pacf(resid_sqrd_norm, lag.max = 50, main = "PACF")

nat_gas_ar2<-ugarchfit(garch_ar2,post_gfc_data[4])
resid_norm <- residuals(nat_gas_ar2, standardize = TRUE)
ArchTest(resid_norm, lags = 12, demean = FALSE)
resid_sqrd_norm <- resid_norm^2 # create squared residuals
par(mfrow = c(1,2))
forecast::Acf(resid_sqrd_norm, lag.max = 50, main = "ACF")
forecast::Pacf(resid_sqrd_norm, lag.max = 50, main = "PACF")

infocriteria(nat_gas_own)
infocriteria(nat_gas_11)
infocriteria(nat_gas_ar2)

