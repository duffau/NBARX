rm(list=ls())
source("init.R")

y_str_vec <- c(
  #"foreclose_sa",
  "marriages_sa" 
  #"divorces_sa",
  #"bankrupt_sa"
  )

lgd_pos_vec <- c("foreclose_sa"=2,"marriages_sa"=1,"divorces_sa"=1,"bankrupt_sa"=1)

i <- 1
for(i in 1:length(y_str_vec)){
  y_str  <- y_str_vec[i]
  lgd_pos <- lgd_pos_vec[y_str]
  
  # =========
  # PAR plots
  # =========
  out <- readRDS(paste0(est_output_path,y_str,".PAR.output.RDS"))
  data <- out$data
  dates <- c(start(data$y),end(data$y))
  N <- length(data$y)
  y <- coredata(data$y)
  
  PAR.par.list <- out$par.list
  PAR.filter.out <- out$filter.out
  p <- PAR.par.list$p
  q <- PAR.par.list$q
  print(sum(PAR.par.list$alpha)+sum(PAR.par.list$beta))
  
  # +++++++++++++++++++++++++++++++++++++++++++++  
  # Plot of 1-step ahead forecast (lambda-filter)
  # +++++++++++++++++++++++++++++++++++++++++++++
  
  # Hit rate
  delta <- 1-as.numeric(names(PAR.filter.out$l.conf))
  N_forecast <- N-max(p,q)
  PAR.hit.count <- hit.count(PAR.filter.out,data$y)
  print(kupiec(N_forecast,PAR.hit.count,delta))
  PAR.hit.rate <- PAR.hit.count/N_forecast
  
  # Time plot
  pdf(paste0(plots_path,"/emp/",y_str,"_PAR_forecast.pdf"),w=w_pdf,h=h_pdf)
  par(par_plot)
  plot.time(data$y,
            PAR.filter.out$lambda.hat,
            PAR.filter.out$l.conf[,"0.95"],
            PAR.filter.out$h.conf[,"0.95"],
            data$y[coredata(data$y<PAR.filter.out$l.conf[,"0.95"])],
            data$y[coredata(data$y>PAR.filter.out$h.conf[,"0.95"])],
            lwd=c(2,2,1,1,1,1),lty=c(1,1,2,2,0,0),col=c(col_vec[1:2],greys[6],greys[6],1,1),
            ylab=paste("Counts"),type=c("S","l","S","S","p","p"),pch=c(32,32,32,32,20,20),ylim=range(data$y))
  l.text <- c("Counts","lambda_t+1|t","95% Conf. bands","Outside conf. bands")
  y.txt <- min(data$y) + diff(range(data$y))*0.825
  if(lgd_pos == 1){
    legend(as.Date(dates[1]),max(data$y),l.text,ncol=2,col=c(col_vec[1:2],greys[6],1),pch=c(32,32,32,19),lty=c(1,1,2,0),lwd=c(2,2,1,0),xjust=0,yjust=1,cex=0.9)
    text(as.Date(dates[1]),y.txt,paste0("Unconditional coverage = ",round(PAR.hit.rate["0.95"]*100,1),"%"),adj=c(-0.02,1),cex=0.9)
  } else if(lgd_pos==2){
    legend(as.Date(dates[2]),max(data$y),l.text,ncol=2,col=c(col_vec[1:2],greys[6],1),pch=c(32,32,32,19),lty=c(1,1,2,0),lwd=c(2,2,1,0),xjust=1,yjust=1,cex=0.9)
    text(as.Date(dates[2]),y.txt,paste0("Unconditional coverage = ",round(PAR.hit.rate["0.95"]*100,1),"%"),adj=c(0.98,1),cex=0.9)
  }
  dev.off()
  
  # ==========
  # NBAR plots
  # ==========
  out <- readRDS(paste0(est_output_path,y_str,".NBAR.output.RDS"))
  data <- out$data
  dates <- c(start(data$y),end(data$y))
  N <- length(data$y)
  y <- coredata(data$y)
  
  NBAR.par.list <- out$par.list
  NBAR.filter.out <- out$filter.out
  p <- NBAR.par.list$p
  q <- NBAR.par.list$q
  print(sum(NBAR.par.list$alpha)+sum(NBAR.par.list$beta))
  
  # +++++++++++++++++++++++++++++++++++++++++++++  
  # Plot of 1-step ahead forecast (lambda-filter)
  # +++++++++++++++++++++++++++++++++++++++++++++
  # Hit rate (unconditional coverage)
  NBAR.hit.rate <- hit.rate(NBAR.filter.out,data$y,N=N-max(p,q))
  print(NBAR.hit.rate)
  
  pdf(paste0(plots_path,"emp/",y_str,"_NBAR_forecast.pdf"),w=w_pdf,h=h_pdf)
  par(par_plot)
  plot.time(data$y,
            NBAR.filter.out$lambda.hat,
            NBAR.filter.out$l.conf[,"0.95"],
            NBAR.filter.out$h.conf[,"0.95"],
            data$y[coredata(data$y<NBAR.filter.out$l.conf[,"0.95"])],
            data$y[coredata(data$y>NBAR.filter.out$h.conf[,"0.95"])],
            lwd=c(2,2,1,1,1,1),lty=c(1,1,2,2,0,0),col=c(col_vec[1:2],greys[6],greys[6],1,1),
            ylab=paste("Counts"),type=c("S","l","S","S","p","p"),pch=c(32,32,32,32,20,20),ylim=range(data$y))
  l.text <- c("Counts","lambda_t+1|t","95% Conf. bands","Outside conf. bands")
  y.txt <- min(data$y) + diff(range(data$y))*0.825
  if(lgd_pos == 1){
    legend(as.Date(dates[1]),max(data$y),l.text,ncol=2,col=c(col_vec[1:2],greys[6],1),pch=c(32,32,32,19),lty=c(1,1,2,0),lwd=c(2,2,1,0),xjust=0,yjust=1,cex=0.9)
    text(as.Date(dates[1]),y.txt,paste0("Unconditional coverage = ",round(NBAR.hit.rate["0.95"]*100,1),"%"),adj=c(0,1),cex=0.9)
  } else if(lgd_pos==2){
    legend(as.Date(dates[2]),max(data$y),l.text,ncol=2,col=c(col_vec[1:2],greys[6],1),pch=c(32,32,32,19),lty=c(1,1,2,0),lwd=c(2,2,1,0),xjust=1,yjust=1,cex=0.9)
    text(as.Date(dates[2]),y.txt,paste0("Unconditional coverage = ",round(NBAR.hit.rate["0.95"]*100,1),"%"),adj=c(1,1),cex=0.9)
  }
  dev.off()
  
  # ++++++++++++++++++++++++++++++++++++++++++++
  # q-q Plots of pseudo-residuals of both models
  # ++++++++++++++++++++++++++++++++++++++++++++
  PAR.res.out <- pseudo_residuals.P(y,PAR.filter.out$lambda.hat)
  NBAR.res.out <- pseudo_residuals.NB(y,lambda.hat=NBAR.filter.out$lambda.hat,nu.hat=NBAR.par.list$nu)
  
  pdf(paste0(plots_path,"emp/",y_str,"_y_qq.pdf"),w=w_pdf,h=h_pdf)
  qq_plot(x=PAR.res.out$res.m,type="p",pch=20,cex=0.5,conf=c("KS","Sim"),col=col_vec[1],xlim=c(-4,4),ylim=c(-4,4))
  # PAR.q.l <- qq_plot(x=PAR.res.out$res.l,plot=F)
  # PAR.q.h <- qq_plot(x=PAR.res.out$res.h,plot=F)
  # segments(PAR.q.l$q_teo, PAR.q.l$q_emp, PAR.q.h$q_teo, PAR.q.h$q_emp,col=col_vec3[1],lwd=3,lend=1)
  
  qq_plot(x=NBAR.res.out$res.m, plot="add",pch=20,cex=0.5,col=col_vec[2])
  # NBAR.q.l <- qq_plot(x=NBAR.res.out$res.l,plot=F)
  # NBAR.q.h <- qq_plot(x=NBAR.res.out$res.h,plot=F)
  # segments(NBAR.q.l$q_teo, NBAR.q.l$q_emp, NBAR.q.h$q_teo,NBAR.q.h$q_emp,col=col_vec3[2],lwd=3,lend=1)
  
  legend(-4,4,c("PAR","NBAR"),col=col_vec,pch=20)
  dev.off()
  
  # ++++++++
  # ACF of y
  # ++++++++
  PAR.theta.hat <- par.list2vec.PAR(PAR.par.list)
  NBAR.theta.hat <- par.list2vec.NBAR(NBAR.par.list)
  
  h <- 24
  lwd <- 2
  ylim <- NULL
  # ylim <- c(0,1)
  pdf(paste0(plots_path,"emp/",y_str,"_acf.pdf"),w=w_pdf,h=h_pdf)
  par(par_plot_acf)
  acf.PAR(data=y,h=h,plot="new",type="b",lend=1,pch=20,lwd=lwd,xlab="Lags",ylab="ACF",ylim=ylim)
  acf.PAR(h,PAR.theta.hat,off_set=-0.2,type="b",lend=1,pch=20,lwd=lwd,p=PAR.par.list$p,q=PAR.par.list$q,col=col_vec[1],plot="add")
  acf.NBAR(h,NBAR.theta.hat,off_set=0.2,type="b",lend=1,pch=20,lwd=lwd,p=NBAR.par.list$p,q=NBAR.par.list$q,col=col_vec[2],plot="add",conf.alpha=0.05)
  legend(h,1,c("Empirical","PAR","NBAR"),col=c(1,col_vec),lwd=lwd,xjust=1)
  dev.off()
  
  # ++++++++++++++++++
  # Predicted variance
  # ++++++++++++++++++
  pdf(paste0("./plots/emp/",y_str,"_y_dispertion.pdf"),w=w_pdf,h=h_pdf)
  par(par_plot)
  plot.time(sqrt(PAR.filter.out$lambda.hat),sqrt(NBAR.filter.out$lambda.hat+NBAR.filter.out$delta.hat),lwd=2,
            ylab="Conditional std. deviation",l.pos=lgd_pos+1,
            l.text=c(expression(paste("PAR:     ", sqrt(lambda[t]))),
                     expression(paste("NBAR: ", sqrt(lambda[t]+lambda[t]^2/nu)))))
  dev.off()
  
  # ++++++++++++++++++++++++++++++++++++
  # Autocorrelation in Pearson-residuals
  # ++++++++++++++++++++++++++++++++++++
  PAR.pearson.res <- (y - PAR.filter.out$lambda.hat)/sqrt(PAR.filter.out$lambda.hat)
  NBAR.pearson.res <- (y - NBAR.filter.out$lambda.hat)/sqrt(NBAR.filter.out$lambda.hat + NBAR.filter.out$delta.hat)
  plot.time(PAR.pearson.res,NBAR.pearson.res,lwd=2,l.pos=1)
  hist(PAR.pearson.res)
  hist(NBAR.pearson.res)
  
  pdf(paste0("./plots/emp/",y_str,"_y_acf_pearson.pdf"),w=w_pdf,h=h_pdf)
  par(par_plot_acf)
  acf.PAR(1:24,data=na.omit(coredata(PAR.pearson.res)),off_set=-0.2,plot="new",type="h",lend=1,lwd=10,ylim=c(-0.2,0.2),col=col_vec[1],xlab="Lags",ylab="ACF")
  acf.PAR(1:24,data=na.omit(coredata(NBAR.pearson.res)),off_set=0.2,plot="add",type="h",lend=1,lwd=10,ylim=c(-0.2,0.2),col=col_vec[2],conf.alpha = c(0.01,0.05,0.1))
  legend(24,-0.2,c("PAR","NBAR"),pch=15,col=col_vec,xjust=1,yjust=0)
  dev.off()
}