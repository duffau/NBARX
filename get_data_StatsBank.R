rm(list=ls())
library(httr)
library(zoo)
library(seasonal)
source("init.R")
setwd(data_path)
f_m <- function(x) as.yearmon(format(x, nsmall = 2), "%YM%m")
f_q <- function(x) as.yearqtr(format(x, nsmall = 2), "%YQ%q")
f_a <- function(x) ceiling(as.numeric(as.yearmon(x)))

functions <- list("m"=f_m,"q"=f_q,"a"=f_a)

jsons <- c("aus01","aus07","aulaar","tvang1","bev3c","deta2","deta21","ejen5","konk9","dnsekt1","dnpin","dnpud")
freq_str <- c("m","m","a","m","m","m","m","q","m","m","m","m")
col_class <- list(c("NULL","NULL","character","numeric"),
                  c("NULL","NULL","character","numeric"),
                  NA,
                  NA,
                  NA,
                  c("NULL","character","numeric"),
                  c("NULL","character","numeric"),
                  NA,
                  NA,
                  c("character","NULL","character","character","numeric"),
                  c("NULL", "NULL","character","NULL","character","numeric"),
                  c("NULL","character","NULL","character","numeric"))
split_var <- list(NULL, NULL,"PERPCT", "TYPE","BEVAEGELSEV",NULL,NULL,"EJENDOMSKATE","SAESON",c("BALPOSTNAT1","SEKTORNAT"),"SEKTOR2","SEKTORNAT")

df_list <- list()

for(i in 1:length(jsons)){
  freq <- freq_str[i]
  r <- POST("http://api.statbank.dk/v1/data", body = upload_file(paste0("./POST_json/",jsons[i],".json")))
  str <- content(r,as="text",encoding="UTF-8")
  str <- gsub("Æ","AE",str)   # reported bug in read.zoo. Cannot split by UTF-8 column names.
  df_list[[jsons[i]]] <- read.zoo(
                          text = str,
                          split=split_var[[i]],
                          index.column ="TID",
                          regular=T,
                          FUN=functions[[freq]],
                          sep=";",
                          na.string=c(".."),
                          header=T,
                          colClasses=col_class[[i]],
                          encoding="UTF-8")
}

data <- do.call("cbind", df_list)
head(data)
plot(data,type="p")



# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# "Jonining" series (stiching together in the time dimension)  
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# "Join" function
join <- function(x) if (is.na(x[1])) x[2] else x[1]

# Rescaling and joinining retail two sales index series into a single series
scale <- mean(data$deta2.deta2/data$deta21.deta21,na.rm=T)
data$deta21.deta21 <- data$deta21.deta21*scale
data$retail_sale <- rollapply(data[,c('deta21.deta21','deta2.deta2')], 1, join, by.column = FALSE) 
plot(data[,c('deta2.deta2','deta21.deta21','retail_sale')],plot.type="single",col=c(1,2,4),lwd=3)
# Deleting original variables
data$deta21.deta21 <- NULL
data$deta2.deta2 <- NULL
head(data)

# Joining unemployment
data$aus07.aus07 <- data$aus07.aus07/1000
plot(data[,c('aus01.aus01','aus07.aus07','Unemployed (thousands)')],plot.type="single",type="p",col=c(1,2,4))
data$unemployed <- rollapply(data[,c('aus07.aus07','aus01.aus01')], 1, join, by.column = FALSE) 
plot(data$unemployed)
data$unemployed <- rollapply(data[,c('unemployed','Unemployed (thousands)')], 1, join, by.column = FALSE) 
plot(data$unemployed,type="p")
# Deleting original variables
data$aus01.aus01 <- NULL
data$aus07.aus07 <- NULL
data$`Unemployed (thousands)` <- NULL
data$`Per cent of the labour force` <- NULL

# Joinning loan and deposit variables
join_vars1 <- c("Deposits.1000: All sectors","X000: All sectors domestic and foreign.dnpin")
data$deposit_all <- rollapply(data[,join_vars1], 1, join, by.column = FALSE) 

join_vars2 <- c("Deposits.1100: Non-financial corporations","- X100: Non-financial corporations.dnpin")
data$deposit_corp <- rollapply(data[,join_vars2], 1, join, by.column = FALSE) 

join_vars3 <- c("Deposits.1415: Households","- X400: Households.dnpin")
data$deposit_household <- rollapply(data[,join_vars3], 1, join, by.column = FALSE) 

join_vars4 <- c("Lending.1000: All sectors" ,"X000: All sectors domestic and foreign.dnpud")
data$loan_all <- rollapply(data[,join_vars4], 1, join, by.column = FALSE)

join_vars5 <- c("Lending.1100: Non-financial corporations","- X100: Non-financial corporations.dnpud")
data$loan_corp <- rollapply(data[,join_vars5], 1, join, by.column = FALSE) 

join_vars6 <- c("Lending.1415: Households","- X400: Households.dnpud")
data$loan_household <- rollapply(data[,join_vars6], 1, join, by.column = FALSE) 
data <- data[,! (names(data) %in% c(join_vars1,join_vars2,join_vars3,join_vars4,join_vars5,join_vars6))]

# Renaming
names(data)
names(data)[names(data)=="Non-seasonally adjusted"] <- "bankrupt"
names(data)[names(data)=="Seasonally adjusted"] <- "bankrupt_sa"
names(data)[names(data)=="Forced sales of real property total"] <- "foreclose"
names(data)[names(data)=="Forced sales of real property total (seasonally adjusted)"] <- "foreclose_sa"
names(data)[names(data)=="Divorces"] <- "divorces"
names(data)[names(data)=="Marriages"] <- "marriages"
names(data)[names(data)=="One-family houses"] <- "house_price"
names(data)[names(data)=="Owner-occupied flats, total"] <- "flat_price"

# Seasonal adjustment
plot(data$divorces)
out <- final(seas(x = as.ts(data$divorces)))
divorces_sa <- zoo(round(out),freq=12)
lines(divorces_sa,col=2)

plot(data$marriages)
out <- final(seas(x = as.ts(data$marriages)))
marriages_sa <- zoo(round(out),freq=12)
lines(marriages_sa,col=2)

data <- merge(data,divorces_sa,marriages_sa)

head(data)
plot(data,type="p",pch=20)
write.zoo(data,file="data.csv")
saveRDS(data,file="data.RDS")
