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
library(tidyverse)
library(sjmisc)



# Reasons for NAs in data:
#
#   * In the month of March 2020, Natural_Gas values have been reported on Sunday
# * The first column "Natural_Gas" is not an american index. Therefore, its values have been reported on days, when there were holidays in the US, resulting in approximately 25 NAs in years 2015 onward.
# * NG has not been reported for two consecutive weeks in September 2005. Since the futures are priced similarly for this period, we will substitute the missing data with the prices of futures

# #### Spikes in Nat_Gas
# Natural gas spot prices at several trading hubs approached their record highs briefly during the week of February 14 amid significantly colder-than-normal weather that affected most of the Lower 48 states. The cold weather led to natural gas supply and demand imbalances. Natural gas production declined because of freeze-offs (temporary interruptions in production caused by cold weather) amid high demand for heating and power. At the benchmark Henry Hub, natural gas prices reached $23.86 per million British thermal units (MMBtu) on February 17, the highest real (inflation-adjusted) price since an Arctic blast on February 25, 2003. Henry Hub prices averaged $5.49/MMBtu in February, the highest monthly average since February 2014.
#
# #### Spikes in NY_Gas
# NEW YORK (CNN/Money) - Consumers can expect retail gas prices to rise to $4 a gallon soon, but whether they stay there depends on the long-term damage to oil facilities from Hurricane Katrina, oil and gas analysts said Wednesday.

 ## ----load_data, include =F-------------------------------------------------------------------------------------------------------------------------------
stylized_names <-c("Natural Gas","Crude Oil WTI","NY Gasoline","Heating Oil","Diesel")

holidays_df<-read_csv("data/US Holiday Dates (2004-2021).csv") %>% dplyr::select(1,2)

nat_gas_spot <- read_excel("data/NG_PRI_FUT_S1_D.xls", sheet = "Data 1", skip = 2) %>%
  rename("Natural_Gas" = 2)
crude_oil_spot <- read_excel("data/PET_PRI_SPT_S1_D.xls", sheet = "Data 1", skip = 2) %>%
  rename("Crude_Oil_WTI" = 2, "Crude_Oil_Europe_Brent" = 3)
gasoline <- read_excel("data/PET_PRI_SPT_S1_D.xls", sheet = "Data 2", skip = 2) %>%
  rename("NY_Gasoline" = 2, "Gulf_Gasoline" = 3)
heating_oil <- read_excel("data/PET_PRI_SPT_S1_D.xls", sheet = "Data 4", skip = 2) %>%
  rename("Heating_Oil" = 2)
diesel <- read_excel("data/PET_PRI_SPT_S1_D.xls", sheet = "Data 5", skip = 2) %>%
  dplyr::select(-2, -3) %>%
  rename("Diesel" = 2)

df <- nat_gas_spot %>%
  full_join(crude_oil_spot, by = "Date") %>%
  full_join(gasoline, by = "Date") %>%
  full_join(heating_oil, by = "Date") %>%
  full_join(diesel, by = "Date") %>%
  left_join(holidays_df,by="Date") %>%
  arrange(Date) %>%
  dplyr::filter(Date >= "1997-01-07") %>%
  mutate(Date=as.Date(Date),
         Weekday = weekdays(Date)) %>%
  dplyr::select(Date,Weekday,Holiday,Natural_Gas,Crude_Oil_WTI,NY_Gasoline,Heating_Oil,Diesel)
df_orig<-df


## ----nat_gas_fut, include = F----------------------------------------------------------------------------------------------------------------------------
nat_gas_fut<- read_excel("data/NG_PRI_FUT_S1_D.xls", sheet = "Data 2", skip = 2) %>%
  rename("Natural_Gas_Fut" = 2) %>% dplyr::select(1,2)

nat_gas_df<-nat_gas_spot%>% full_join(nat_gas_fut) %>% filter(Date >= "2005-09-16",Date<="2005-10-17") %>% arrange(Date)

## ----NA before, include = F------------------------------------------------------------------------------------------------------------------------------
df$Month<-month(df$Date)
df$Day<-day(df$Date)
df<-df%>%
  filter(!(Day %in%c(24,25,26,31) & Month == 12)) %>%
  filter(!(Day %in%c(1,2) & Month == 1)) %>% select(-Day,-Month) %>%
  dplyr::filter(!Weekday %in% c("Sunday","Saturday"))


NA_dates <- df %>%
  filter_at(vars(-Date,-Weekday,-Holiday), all_vars(is.na(.))) %>%
  pull(Date)
df <- df %>% dplyr::filter(!Date %in% NA_dates)

## ----NA after, include = F-------------------------------------------------------------------------------------------------------------------------------
##replace missing Nat_Gas spot prices with futures
na_2005<-nat_gas_df[!complete.cases(nat_gas_df),] %>% mutate(Date=as.Date(Date)) %>% dplyr::select(-"Natural_Gas")
df<-df%>%left_join(na_2005) %>% mutate(Natural_Gas = coalesce(Natural_Gas,Natural_Gas_Fut)) %>%dplyr::select(-"Natural_Gas_Fut")

#Remove data where there was a holiday and the metrics were not reported
df<-df %>% dplyr::filter(!(!is.na(Holiday) &
                             (is.na(Crude_Oil_WTI) | is.na(Natural_Gas) |
                                is.na(NY_Gasoline) | is.na(Heating_Oil) | is.na(Diesel))
                           ))
df_NA<-df %>% filter_at(vars(-Date,-Weekday,-Holiday), any_vars(is.na(.)))
# rest is substituted with the previous days price
for(i in 4:length(df)){
  df[,i]<-imputeTS::na_ma(df[,i], k = 1, weighting = "simple")
}

## ----crude_oil_sub, include=F----------------------------------------------------------------------------------------------------------------------------
#We see that this **crude oil WTI metric** had negative prices at one point. We substitute it with the average of its previous and following days value.
for (i in 1:nrow(df)){
  val<-df[i,"Crude_Oil_WTI"]
  if(val<0){
    df[i,"Crude_Oil_WTI"] <-(df[i-1,"Crude_Oil_WTI"]+df[i+1,"Crude_Oil_WTI"])/2
  }
}
df[df$Date=="2020-04-20","Crude_Oil_WTI"]


## ----log_ret, include=F----------------------------------------------------------------------------------------------------------------------------------
log_ret_df <- tibble(df[2:nrow(df), "Date"])
for (i in 4:(ncol(df))) {
  var_name <- paste(colnames(df)[i])
  vals <- log(df[, i]) %>% pull()
  lag_vals <- base::diff(vals, lag = 1)
  log_ret_df[var_name] <- lag_vals
}
weekdays<-df[2:nrow(df), "Weekday"] %>% as.vector()
weekdays<-sjmisc::to_dummy(weekdays,suffix="label")
log_ret_df_dum<-cbind(log_ret_df,weekdays) %>% select(-Weekday)
log_ret_df_dum$Covid<-ifelse(log_ret_df_dum$Date>="2020-01-01",1,0)
head(log_ret_df_dum)
write_csv(log_ret_df_dum,"data/dummy_log_ret.csv")
saveRDS(log_ret_df,"data/log_ret.rds")
saveRDS(df,"data/df.rds")
