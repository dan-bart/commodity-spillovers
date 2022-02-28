
# Commodity - SPGSCITR
# Energy - GSPE
# Technology - NASDAQ
# Global - S&P
# VIX, NEWS
library(tibble)
library(xts)
library(quantmod)
library(readxl)
library(dplyr)
library(readr)
library(tseries)
library(ggcorrplot)
library(psych)
library(pls)


log_ret<-readRDS("data/log_ret.rds")


GSCI_from <- read_csv("data/pca_indexes/GSCI_from.csv") %>%
  mutate(Date = as.Date(Date,format = "%b %d, %Y")) %>% arrange(Date) %>%
  dplyr::select(Date,Price) %>% rename("GSCI" = "Price")

GSCI_to <- read_csv("data/pca_indexes/GSCI_to.csv") %>%
  mutate(Date = as.Date(Date,format = "%b %d, %Y")) %>% arrange(Date) %>%
  filter(Date<="2021-11-09") %>% dplyr::select(Date,Price) %>% rename("GSCI" = "Price")
GSCI<-rbind(GSCI_from,GSCI_to)


GSPE<-getSymbols.yahoo(c("^GSPE"),from = as.Date("1997-01-07"), to = as.Date("2021-11-09"),
                       auto.assign = F,periodicity = "daily") %>%
  as.data.frame() %>% rownames_to_column("Date") %>% rename("GSPE" = "GSPE.Open") %>%
  dplyr::select(Date,GSPE) %>% mutate(Date = as.Date(Date))

NASDAQ<-getSymbols.yahoo(c("^IXIC"),from = as.Date("1997-01-07"), to = as.Date("2021-11-09"),
                         auto.assign = F,periodicity = "daily") %>%
  as.data.frame() %>% rownames_to_column("Date") %>% rename("NASDAQ" = "IXIC.Open") %>%
  dplyr::select(Date,NASDAQ) %>% mutate(Date = as.Date(Date))
SP<-getSymbols.yahoo(c("^GSPC"),from = as.Date("1997-01-07"), to = as.Date("2021-11-09"),
                     auto.assign = F,periodicity = "daily") %>%
  as.data.frame() %>% rownames_to_column("Date") %>% rename("SP" = "GSPC.Open") %>%
  dplyr::select(Date,SP) %>% mutate(Date = as.Date(Date))

VIX<-getSymbols.yahoo(c("^VIX"),from = as.Date("1997-01-07"), to = as.Date("2021-11-09"),
                auto.assign = F,periodicity = "daily") %>%
  as.data.frame() %>% rownames_to_column("Date") %>% rename("VIX" = "VIX.Open") %>%
  dplyr::select(Date,VIX) %>% mutate(Date = as.Date(Date))

NEWS<-read_excel("data/pca_indexes/news_sentiment_data.xlsx", sheet = "Data") %>% rename("Date" = "date","Sentiment" = 2) %>%
  mutate(Date = as.Date(Date))  %>% filter(Date>="1997-01-07",Date<="2021-11-09")

df<- GSCI%>%
  left_join(NEWS,by="Date")%>%
  left_join(GSPE,by="Date")%>%
  left_join(NASDAQ,by="Date")%>%
  left_join(VIX,by="Date")%>%
  left_join(SP,by="Date")

df[!complete.cases(df),]
df<-df %>% imputeTS::na_ma(k = 1, weighting = "simple")
adf.test(df[-1,"Sentiment"]%>%pull()) # reject nn-stationarity with pval 0.01

log_df <-df  %>% dplyr::select(-Date,-Sentiment) %>% apply(2,function(x) diff(log(x), lag=1)) %>% as.data.frame() %>%cbind("Sentiment"=df$Sentiment[-1],"Date"=as.Date(df$Date[-1]) )
norm_df<-df  %>% dplyr::select(-Date) %>% apply(2,function(x) (x-mean(x))/sd(x))%>% as.data.frame() %>%cbind("Date"=as.Date(df$Date))

# final_df<-log_ret%>%left_join(norm_df,by="Date")
# final_df[!complete.cases(final_df),]
# final_df<-final_df %>% imputeTS::na_ma(k = 1, weighting = "simple")
final_df<-norm_df
cor_df<-cor(final_df[-c(length(final_df))])
cor_ret_df<-cor(log_df[-c(length(log_df))])
#cor_df<-cor(final_df[,2:6],final_df[,7:12])
ggcorrplot(cor_df,lab = TRUE)
fa.parallel(cor_df, n.obs=nrow(norm_df)-3, fa="pc", n.iter=100,
            main="Screen plot with parallel analysis")

PC<-principal(norm_df[-length(norm_df)], nfactors=3, rotate="none",score=TRUE)
z<-pcr(GSCI~Sentiment+GSPE+NASDAQ+VIX+SP,data=df,scale = TRUE, validation = "CV",ncomp=3)
summary(z)
z$coefficients

saveRDS(cor_ret_df,"data/pca/cor_ret_df.rds")
saveRDS(cor_df,"data/pca/pca_corr.rds")
saveRDS(PC$scores,"data/pca/pca_vals.rds")


