rm(list=ls())
source("init.R")


# Load data from StatBank Denmark (Statistics Denmark)
# ----------------------------------------------------
y_str_vec <- c(
  #"marriages_sa", 
  "divorces_sa",
  #"foreclose_sa",
  #"bankrupt_sa",
  #"marriages", 
  "divorces"
  #"foreclose",
  #"bankrupt"
  )


neg_pars <- c(
  "marriages" = TRUE,
  "divorces"  = TRUE,
  "foreclose" = TRUE,
  "bankrupt"  = TRUE,
  "marriages_sa" = FALSE,
  "divorces_sa"  = FALSE,
  "foreclose_sa" = FALSE,
  "bankrupt_sa"  = FALSE
)

lags <- list(
  "marriages_sa" = list("p.PAR"=4,"q.PAR"=2,"p.NBAR"=2,"q.NBAR"=2),
  "divorces_sa"  = list("p.PAR"=1:4,"q.PAR"=1:1,"p.NBAR"=1:4,"q.NBAR"=1:1),
  "foreclose_sa" = list("p.PAR"=4,"q.PAR"=1,"p.NBAR"=4,"q.NBAR"=1),
  "bankrupt_sa"  = list("p.PAR"=4,"q.PAR"=1,"p.NBAR"=4,"q.NBAR"=1),
  
  "marriages" = list("p.PAR"=13,"q.PAR"=2,"p.NBAR"=13,"q.NBAR"=2),
  "divorces"  = list("p.PAR"=c(1,12,13),"q.PAR"=1:2,"p.NBAR"=c(1,12,13),"q.NBAR"=1:2),
  "foreclose" = list("p.PAR"=13,"q.PAR"=2,"p.NBAR"=13,"q.NBAR"=2),
  "bankrupt"  = list("p.PAR"=13,"q.PAR"=2,"p.NBAR"=13,"q.NBAR"=2)
  )


for(i in 1:length(y_str_vec)){
  y_str <- y_str_vec[i]
  data <- readRDS(paste0(data_path,"data.RDS"))
  data <-  data[!is.na(data[,y_str]),]
  data$y <- data[,y_str]
  y <- coredata(data$y)
  N <- length(y)
  neg_par <- neg_pars[y_str]
  
  plot.time(data$y)
  
  # Mean by month
  plot(aggregate(data$y, by=format(time(data$y), "%m"), FUN=mean),type="h",lwd=10,lend=1)
  
  # ++++++++++++++++++
  # Estimate PAR model
  # ++++++++++++++++++
  p <- lags[[y_str]][["p.PAR"]]
  q <- lags[[y_str]][["q.PAR"]]
  p_len <- length(p)
  q_len <- length(q)
  par.names <- c("omega",paste0("alpha_",p),paste0("beta_",q))
  
  # Initializing parameters
  omega0 <- mean(y)*(1-0.3-0.5)
  theta0 <- c(omega0,rep(0.3/p_len,p_len),rep(0.5/q_len,q_len))+runif(1+p_len+q_len,-0.05/max(p_len,q_len),0.05/max(p_len,q_len))
  gamma0 <- repam.inv.PAR(theta0,neg_par=neg_par)
  
  # Evaluating log-likelihood at intial parameter point
  print(loglike.PAR(gamma0,p=p,q=q,y))
  
  # Minimizing negative log-likelihood
  out <- optim(par=gamma0,loglike.PAR,p=p,q=q,y=y,neg_par=neg_par,method="BFGS",hessian=T,control=list(trace=2,maxit=500,reltol=1e-8,REPORT=5))
  summary.out <- summary.avg_mle(N,out,repam.PAR,par.names=par.names)
  print("---- Summary PAR ----")
  print(summary.out)
  
  # Transforming parameters from reparametrized version to "regular" version
  theta.hat <- repam.PAR(out$par)
  par.list <- PAR.par(theta.hat,p,q)
  print(sum(par.list$alpha)+sum(par.list$beta))
  
  filter.out <- filter.PAR(theta.hat,p,q,data$y,conf=c(0.99,0.95,.90),zoo=T)
  plot.time(filter.out$lambda.hat,data$y)
  
  # Kupiec Test
  delta <- 1-as.numeric(names(filter.out$l.conf))
  N_forecast <- N-max(p,q)
  hit.count.out <- hit.count(filter.out,data$y)
  kupiec.out <- kupiec(N_forecast,hit.count.out,delta)
  print(kupiec.out)
  
  # Plot HP-filtered residuals and PAR-filtered residuals
  # res_hp <- hpfilt(data$y)-data$y
  # res_PAR <- PAR.filter.out$lambda.hat-data$y
  # plot.time(res_hp,res_PAR,diff(data$y),l.pos=1,lwd=2,startplot="2000-01-01",endplot="2006-12-31")
  
  # Export parameter estimates, filtered values, summary etc.
  out_rds <- list(data=data,
                  par.list=par.list,
                  filter.out=filter.out,
                  summary.out=summary.out,
                  kupiec.out=kupiec.out)
  saveRDS(out_rds,paste0(est_output_path,y_str,".PAR.output.RDS"))
  
  # Plot empirical and theoretical ACF
  h <- 50
  acf.PAR(h,data=y,plot="new",pch=20)
  acf.PAR(h,theta.hat,p,q,plot="add",col=col_vec[1],pch=20)
  
  # +++++++++++++++++++
  # Estimate NBAR model
  # +++++++++++++++++++
  p <- lags[[y_str]][["p.NBAR"]]
  q <- lags[[y_str]][["q.NBAR"]]
  p_len <- length(p)
  q_len <- length(q)
  par.names <- c("omega",paste0("alpha",p),paste0("beta",q),"nu")
  
  omega0 <- mean(y)*(1-0.3-0.5)
  theta0 <- c(omega0,rep(0.3/p_len,p_len),rep(0.5/q_len,q_len),1)+runif(2+p_len+q_len,-0.05/max(p_len,q_len),0.05/max(p_len,q_len))
  gamma0 <- repam.inv.NBAR(theta0,neg_par=neg_par)
  loglike.NBAR(gamma0,p=p,q=q,y)
  
  out <- optim(par=gamma0,loglike.NBAR,p=p,q=q,y=y,neg_par=neg_par,method="BFGS",hessian=T,control=list(trace=2,maxit=500,reltol=1e-8,REPORT=5))
  summary.out <- summary.avg_mle(N,out,repam.NBAR,par.names=par.names)
  print("---- Summary NBAR ----")
  print(summary.out)
  
  theta.hat <- repam.NBAR(out$par)
  par.list <- NBAR.par(theta.hat,p,q)
  print(sum(par.list$alpha)+sum(par.list$beta))
  
  filter.out <- filter.NBAR(theta.hat,p,q,data$y,conf=c(0.99,0.95,0.90),zoo=T)
  plot.time(filter.out$lambda.hat,data$y)
  
  # Kupiec Test
  delta <- 1-as.numeric(names(filter.out$l.conf))
  N_forecast <- N-max(p,q)
  hit.count.out <- hit.count(filter.out,data$y)
  kupiec.out <- kupiec(N_forecast,hit.count.out,delta)
  print(kupiec.out)
  
  h <- 50
  acf.PAR(h,data=y,plot="new",pch=20)
  acf.PAR(h,theta.hat,p,q,plot="add",col=col_vec[1],pch=20)
  
  
  # Export parameter estimates, filtered values etc.
  out_rds <- list(data=data,
                  par.list=par.list,
                  filter.out=filter.out,
                  summary.out=summary.out,
                  kupiec.out=kupiec.out)
  saveRDS(out_rds,paste0(est_output_path,y_str,".NBAR.output.RDS"))
  
}