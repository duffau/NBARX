rm(list= ls())
cat("\014")
setwd("~/Dropbox/Speciale")
source("./r_scripts/init.R")

# ++++++++++++
# Loading data
# ++++++++++++
data <- readRDS("./data/data.RDS")
PAR.filter <- readRDS("./EstOutput/foreclose_sa.PAR.output.RDS")$filter.out
f <- function(x) as.yearmon(format(x, nsmall = 2), "%YM%m")
df <- read_excel(path = "./data/dors_data/konj_data_monthly.xlsx", sheet = 1, skip=1)
industry <- read.zoo(df,FUN=f)$GIKI.M

f <- function(x) as.yearmon(as.Date(x, origin="1899-12-30"))
cpi <- read.zoo(read_excel(path = "./data/dors_data/mb_data.xlsx", sheet = 2),FUN=f)
yield <- read.zoo(read_excel(path = "./data/dors_data/mb_data.xlsx", sheet = 1),FUN=f)

# Quarterly data
# --------------
f <- function(x) as.yearqtr(format(x, nsmall = 2), "%YQ%q")
df <- read_excel(path = "./data/dors_data/konj_data_quarterly.xlsx", sheet = 1, skip=1)
quarter <-read.zoo(df,FUN=f)

# ++++++++++++++++++++++++++++++++++++++++++
# Tranforming variables and merging datasets
# ++++++++++++++++++++++++++++++++++++++++++
data <- na.approx(data) # Filling gaps with linear interpolation

# Transforming and combining variables
gdp <- pchya(quarter$GFY.Q)
data$ltd_all <- data$loan_all/data$deposit_all # Loan to deposit ratio whole danish economy
data$ltd_household <- data$loan_household/data$deposit_household # Loan to deposit ratio danish housholds 
data$ltd_corp <- data$loan_corp/data$deposit_corp # Loan to deposit ratio non-financial sector

d_data <- diff(data, na.pad = TRUE) # Calculating first difference for all vars
pch_data <- pchya.m(data) # Calculating yea-to-year growhth rate for all vars

# Renaming fist diff and growth rate data
names(d_data) <- paste0(names(data),".d") 
names(pch_data) <- paste0(names(data),".pch")

# Merging first diff and growth rate to original level data
data <- cbind(data,d_data)
data <- cbind(data,pch_data)
# head(data)

# Calculating negative and positive part of time series
retail_neg <- -pmin(pchya.m(data$retail_sale),0)
retail_pos <- pmax(pchya.m(data$retail_sale),0) 
industry_neg <- -pmin(industry,0)
industry_pos <- pmax(industry,0)
gdp_neg <- -pmin(gdp,0)
gdp_pos <- pmax(gdp,0)
unemployed.d_neg <- -pmin(data$unemployed.d,0)
unemployed.d_pos <- pmax(data$unemployed.d,0)

# 30 year danish bond yield deflated with consummer price index
real_interest <- pmax(yield-pchya.m(cpi)*100,0)
real_interest.d <- diff(real_interest)

# Choose data for model
data_est <- merge(
              "foreclose_sa"=data$foreclose_sa,
              "divorce.pch"=lag(data$divorces.pch,-6),
              "house_price.d"=lag(data$house_price.d,-1), #-9
              "unemployed.d"=lag(data$unemployed.d,-4),
              #"ltd_household.d"=lag(data$ltd_household.d,-1),
              #"real_interest"=lag(real_interest.d,-3),
              #"retail_sale.pch"=lag(data$retail_sale.pch,-9),
              #retail_pos,
              #gdp_neg,
              #gdp_pos,
              "constant"=1,
              fill=NA,all=F)

data_est <- na.trim(data_est,side="both")
data_est$y <- data_est$foreclose_sa
# plot(data_est)
# saveRDS(data_est,"./EstOutput/foreclose_sa.data_est.RDS")


# +++++++++++++++++++++++++++++++++++++++++
# Preliminary linear model on PAR-residuals
# +++++++++++++++++++++++++++++++++++++++++
foo  <- merge(y=data_est$y,lambda.hat=PAR.filter$lambda.hat)
res  <- foo$y - foo$lambda.hat
lm_y <- coredata(window(res,start=start(data_est),end=end(data_est)))
lm_X <- coredata(data_est[,!names(data_est) %in% c("foreclose_sa","y","constant")])
summary_lm <- summary(lm(lm_y ~ lm_X))
summary_lm

y <- coredata(data_est$y)
N <- length(y)
print(N)
x_t <- data_est[,!names(data_est) %in% c("y","foreclose_sa"),drop=F]
x <- coredata(x_t)
x_names <- colnames(x)
# x <- matrix(0,N,1)

# ++++++++++++++++
# Model estimation
# ++++++++++++++++

# Lags
# ----
p <- 1:3
q <- 1
p_len <- length(p)
q_len <- length(q)
fun_x <- fun_x_pos
k_x <- dim(x)[2]
PARX.par.names <- c("omega",paste0("alpha",1:p_len),paste0("beta",1:q_len),x_names)
NBARX.par.names <- c("omega",paste0("alpha",1:p_len),paste0("beta",1:q_len),x_names,"nu")
dates <- c(start(data_est$y),end(data_est$y))
conf_vec <- c(0.99,0.95,0.90)

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Estimating AR(1)-models for each x-variable for later stationarity check
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
A <- rep(NA,k_x-1)
for(i in 1:(k_x-1)){
  foo_xt   <- x[-1,i]
  foo_xtm1 <- cbind(1,x[-N,i])
  A[i] <- solve(t(foo_xtm1)%*%foo_xtm1,t(foo_xtm1)%*%foo_xt)[2]
}

# theta_sim <- c(1,0.45,0.52,rep(0,),5)
# out <- sim.PARX(N=N,theta=theta_sim,p=p,q=q,x=x,fun_x=fun_x,seed=NULL,all.out=T)
# plot(out$fx)
# plot(out$y,ylim=c(0,max(out$y,out$fx)),type="l",col=2,lwd=3)
# lines(out$fx,col=4,lwd=3)
# lines(out$lambda,col=1,lwd=3)
# gamma_sim <- repam.inv.PARX(theta_sim,p,q,k_x)
# loglike.PARX(gamma_sim,p=p,q=q,y,x,fun_x=fun_x)
# 
# y <- out$y
# data_est$y <- out$y

# +++++++++++++++++++++++++++++++++++++
# Estimate PARX model and make forecast
# +++++++++++++++++++++++++++++++++++++
# Initializing parameter vector
omega0 <- mean(y)*(1-0.45-0.52)
theta0 <- c(omega0,rep(0.4/p_len,p_len),rep(0.5/q_len,q_len),rep(0,k_x))+runif(1+p_len+q_len+k_x,-0.05,0.05)
gamma0 <- repam.inv.PARX(theta0,p,q,k_x)
loglike.PARX(gamma0,p=p,q=q,y,x,fun_x=fun_x)

# Minimizing negative of the log likelihood
out <- optim(par=gamma0,loglike.PARX,p=p,q=q,y=y,x=x,fun_x=fun_x,print.par=F,hessian=T,method="BFGS",control=list(trace=1,maxit=1500,reltol=1e-8,REPORT=15))

# Reparametrizing estimated parameters and saving in list object
theta.hat <- repam.PARX(out$par,p,q,k_x)
par.list <- PARX.par(theta.hat,p,q,k_x)

# Calculating standard errors and p-values from invers Hessian
summary.out <- summary.avg_mle(N,out,repam.PARX,par.names=PARX.par.names,p=p,q=q,k_x=k_x)
summary.out

# Checking stationarity
s <- 1
pow <- 4
rho <- (sum(abs(A)**pow)**1/pow)**s
eps <- max(0,(1-sum(par.list$beta))*(rho-par.list$alpha[1]))
sum(par.list$beta) + sum(par.list$alpha)
1-eps
summary_lm

# Calculating one-step ahead filter
filter.out <- filter.PARX(theta.hat,p,q,data_est$y,x,fun_x,conf=conf_vec,zoo=T)

# Kupiec Test for unconditional coverage
delta <- 1-conf_vec
N_forecast <- N-max(p,q)
hit.count.out <- hit.count(filter.out,data_est$y)
kupiec.out <- kupiec(N_forecast,hit.count.out,delta)
print(kupiec.out)

# Saving output
out_rds <- list(data=data_est,
                par.list=par.list,
                filter.out=filter.out,
                summary.out=summary.out,
                kupiec.out=kupiec.out)
saveRDS(out_rds,paste0("./EstOutput/foreclose_sa.PARX.output.RDS"))

# ++++++++++++++++++++
# Estimate NBARX model
# ++++++++++++++++++++
omega0 <- mean(y)*(1-0.45-0.52)
theta0 <- c(omega0,rep(0.4/p_len,p_len),rep(0.5/q_len,q_len),rep(0,k_x),10)+runif(1+p_len+q_len+k_x+1,-0.05,0.05)
gamma0 <- repam.inv.NBARX(theta0,p,q,k_x)
loglike.NBARX(gamma0,p=p,q=q,y,x,fun_x=fun_x)

out <- optim(par=gamma0,loglike.NBARX,p=p,q=q,y=y,x=x,fun_x=fun_x,print.par=F,hessian=T,method="BFGS",control=list(trace=1,maxit=1500,reltol=1e-8,REPORT=15))
theta.hat <- repam.NBARX(out$par,p,q,k_x)
par.list <- NBARX.par(theta.hat,p,q,k_x)

summary.out <- summary.avg_mle(N,out,repam.NBARX,par.names=NBARX.par.names,p=p,q=q,k_x=k_x)
summary.out

s <- 1
pow <- 4
rho <- (sum(abs(A)**pow)**1/pow)**s
eps <- max(0,(1-sum(par.list$beta))*(rho-par.list$alpha[1]))
sum(par.list$beta) + sum(par.list$alpha)
1-eps
summary_lm

filter.out <- filter.NBARX(theta.hat,p,q,data_est$y,x,fun_x,conf=conf_vec,zoo=T)

# Kupiec Test
delta <- 1-conf_vec
N_forecast <- N-max(p,q)
hit.count.out <- hit.count(filter.out,data_est$y)
kupiec.out <- kupiec(N_forecast,hit.count.out,delta)
print(kupiec.out)

# Saving output
out_rds <- list(data=data_est,
                par.list=par.list,
                filter.out=filter.out,
                summary.out=summary.out,
                kupiec.out=kupiec.out)
saveRDS(out_rds,paste0("./EstOutput/foreclose_sa.NBARX.output.RDS"))
