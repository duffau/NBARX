rm(list=ls())
setwd("/home/christian/Dropbox/Speciale")
source("./r_scripts/init.R")

y_str <- "foreclose_sa"
out <- readRDS(paste0("./EstOutput/",y_str,".PARX.output.RDS"))
data <- out$data
dates <- c(start(data$y),end(data$y))
N <- length(data$y)
y <- coredata(data$y)
x_t <- data_est[,!names(data_est) %in% c("y",y_str),drop=F]
x <- coredata(x_t)
x_names <- colnames(x)

# PARX parammeters and filtered values
PARX.par.list <- out$par.list
PARX.filter.out <- out$filter.out
PARX.summary.out <- out$summary.out
PARX.filter.out <-  out$filter.out

out <- readRDS(paste0("./EstOutput/",y_str,".NBARX.output.RDS"))
NBARX.par.list <- out$par.list
NBARX.summary <- out$summary.out
NBARX.filter.out <- out$filter.out

# ++++++++++++++++++++++++
# Plot exogenous variables
# ++++++++++++++++++++++++
for(name in setdiff(names(data_est),c("y",y_str,"constant"))){
  # png(paste0("./plots/emp/",y_str,"_",name,".png"),w=700,h=500)
  pdf(paste0("./plots/emp/",y_str,"_",name,".pdf"),w=w_pdf,h=h_pdf)
  par(par_plot)
  plot.time(data_est[,name],main=name)
  dev.off()
}

# ++++++++++
# PARX plots
# ++++++++++
p <- PARX.par.list$p
q <- PARX.par.list$q
k_x <- PARX.par.list$k_x
print(PARX.summary)


# Plot of 1-step ahead pmf-forecast
# ---------------------------------

# Hit rate
delta <- 1-as.numeric(names(PARX.filter.out$l.conf))
N_forecast <- N-max(p,q)
PARX.hit.count <- hit.count(PARX.filter.out,data_est$y)
kupiec(N_forecast,PARX.hit.count,delta)
PARX.hit.rate <- PARX.hit.count/N_forecast

# Forecast plot
pdf(paste0("./plots/emp/",y_str,"_PARX_forecast.pdf"),w=w_pdf,h=h_pdf)
par(par_plot)
plot.time(data_est$y,
          PARX.filter.out$lambda.hat,
          PARX.filter.out$l.conf[,"0.95"],
          PARX.filter.out$h.conf[,"0.95"],
          data_est$y[coredata(data_est$y<PARX.filter.out$l.conf[,"0.95"])],
          data_est$y[coredata(data_est$y>PARX.filter.out$h.conf[,"0.95"])],
          lwd=c(2,2,1,1,1,1),lty=c(1,1,2,2,0,0),col=c(col_vec[1:2],greys[6],greys[6],1,1),
          ylab=paste("Counts"),type=c("S","l","S","S","p","p"),pch=c(32,32,32,32,20,20),ylim=range(y))
l.text <- c("Counts",expression(hat(lambda)[paste(t+1,"|",t)]),"95% Conf. bands","Outside conf. bands")
legend(as.Date(dates[2]),max(y),l.text,col=c(col_vec[1:2],greys[6],1),pch=c(32,32,32,19),lty=c(1,1,2,0),lwd=c(2,2,1,0),xjust=1,yjust=1)
text(as.Date(dates[2]),mean(range(y)),paste0("Unconditional coverage = ",round(PARX.hit.rate["0.95"]*100,1),"%"),adj=c(1,0.5))
dev.off()

# Bar plot of decomposed effect from exogeneous variables 
# -------------------------------------------------------
x_times_gamma <- x_t*matrix(PARX.par.list$gam_x,N,k_x,byrow=T)
m <- apply(x_times_gamma,2,mean)
ord.m <- order(m,decreasing=F)
names_ordered <- setdiff(names(m)[ord.m],"constant")
exo.comp.pos <-exo.comp.neg <- matrix(NA,N,k_x)
time <- as.Date(time(x_times_gamma))
const <- PARX.par.list$gam_x[k_x]

exo.comp.pos[,1] <- 0
exo.comp.neg[,1] <- 0
for(i in 2:(k_x)){
  exo.comp.pos[,i] <- exo.comp.pos[,i-1] + pmax(x_times_gamma[,names_ordered[i-1]],0)
  exo.comp.neg[,i] <- exo.comp.neg[,i-1] + pmin(x_times_gamma[,names_ordered[i-1]],0)
}
exo.comp.neg <- pmax(exo.comp.neg,-const)

par(mar=c(2,4,1,1))
foo_zoo <- zoo(PARX.filter.out$fx-const,order.by=time)
rng <- c(-const,max(c(exo.comp.pos,exo.comp.neg)))

pdf(paste0("./plots/emp/",y_str,"_PARX_xeffect_decomposed.pdf"),w=w_pdf,h=h_pdf)
par(par_plot)
plot.time(foo_zoo,ylab="Foreclosure counts",type="n",yaxt="n",ylim=rng)
at <- seq(rng[1],rng[2],10)
axis(2,at=at,lab=seq(0,(length(at)-1)*10,10))
for(i in (k_x):2){
  lines(time,exo.comp.pos[,i],type="h",lend=1,lwd=1.3,col=col_vec5[i-1])
  lines(time,exo.comp.neg[,i],type="h",lend=1,lwd=1.3,col=col_vec5[i-1])
}
lines(time+15,PARX.filter.out$fx-const,lwd=1.5,type="s")
legend(time[1],rng[2],ncol=1,c("Total exogenous effect",names_ordered),xjust=0,
       lwd=c(1.2,1,1,1),
       lty=c(1,NA,NA,NA),
       pch=c(NA,15,15,15),
       pt.cex = 2,
       col=c(1,col_vec5))
dev.off()

# +++++++++++
# NBARX plots
# +++++++++++
p <- NBARX.par.list$p
q <- NBARX.par.list$q
k_x <- NBARX.par.list$k_x
print(NBARX.summary)

# Hit rate
delta <- 1-as.numeric(names(NBARX.filter.out$l.conf))
N_forecast <- N-max(p,q)
NBARX.hit.count <- hit.count(NBARX.filter.out,data_est$y)
kupiec(N_forecast,NBARX.hit.count,delta)
NBARX.hit.rate <- NBARX.hit.count/N_forecast

# Plot of 1-step ahead pmf-forecast
pdf(paste0("./plots/emp/",y_str,"_NBARX_forecast.pdf"),w=w_pdf,h=h_pdf)
# png(paste0("./plots/emp/",y_str,"_NBARX_forecast.png"),w=700,h=500)
par(par_plot)
plot.time(data_est$y,
          NBARX.filter.out$lambda.hat,
          NBARX.filter.out$l.conf[,"0.95"],
          NBARX.filter.out$h.conf[,"0.95"],
          data_est$y[coredata(data_est$y<NBARX.filter.out$l.conf[,"0.95"])],
          data_est$y[coredata(data_est$y>NBARX.filter.out$h.conf[,"0.95"])],
          lwd=c(2,2,1,1,1,1),lty=c(1,1,2,2,0,0),col=c(col_vec[1:2],greys[6],greys[6],1,1),
          ylab=paste("Counts"),type=c("S","l","S","S","p","p"),pch=c(32,32,32,32,20,20),ylim=range(y))
l.text <- c("Counts",expression(hat(lambda)[paste(t+1,"|",t)]),"95% Conf. bands","Outside conf. bands")
legend(as.Date(dates[2]),max(y),l.text,col=c(col_vec[1:2],greys[6],1),pch=c(32,32,32,19),lty=c(1,1,2,0),lwd=c(2,2,1,0),xjust=1,yjust=1)
text(as.Date(dates[2]),mean(range(y)),paste0("Unconditional coverage = ",round(NBARX.hit.rate["0.95"]*100,1),"%"),adj=c(1,0.5))
dev.off()

# Bar plot of decomposed effect from exogeneous variables 
# -------------------------------------------------------
x_times_gamma <- x_t*matrix(NBARX.par.list$gam_x,N,k_x,byrow=T)
m <- apply(x_times_gamma,2,mean)
ord.m <- order(m,decreasing=F)
names_ordered <- setdiff(names(m)[ord.m],"constant")
exo.comp.pos <-exo.comp.neg <- matrix(NA,N,k_x)
time <- as.Date(time(x_times_gamma))
const <- NBARX.par.list$gam_x[k_x]

exo.comp.pos[,1] <- 0
exo.comp.neg[,1] <- 0
for(i in 2:(k_x)){
  exo.comp.pos[,i] <- exo.comp.pos[,i-1] + pmax(x_times_gamma[,names_ordered[i-1]],0)
  exo.comp.neg[,i] <- exo.comp.neg[,i-1] + pmin(x_times_gamma[,names_ordered[i-1]],0)
}
exo.comp.neg <- pmax(exo.comp.neg,-const)

par(mar=c(2,4,1,1))
foo_zoo <- zoo(NBARX.filter.out$fx-const,order.by=time)
rng <- c(-const,max(c(exo.comp.pos,exo.comp.neg)))

# png(paste0("./plots/emp/",y_str,"_NBARX_xeffect_decomposed.png"),w=700,h=500)
pdf(paste0("./plots/emp/",y_str,"_NBARX_xeffect_decomposed.pdf"),w=w_pdf,h=h_pdf)
par(par_plot)
plot.time(foo_zoo,ylab="Foreclosure counts",type="n",yaxt="n",ylim=rng) #,startplot="2000-01-01")
at <- seq(rng[1],rng[2],10)
axis(2,at=at,lab=seq(0,(length(at)-1)*10,10))
for(i in (k_x):2){
  lines(time,exo.comp.pos[,i],type="h",lend=1,lwd=1.3,col=col_vec5[i-1])
  lines(time,exo.comp.neg[,i],type="h",lend=1,lwd=1.3,col=col_vec5[i-1])
}
lines(time+15,NBARX.filter.out$fx-const,lwd=1.1,type="S")
legend(time[1],rng[2],ncol=1,c("Total exogenous effect",names_ordered),xjust=0,
       lwd=c(1.2,1,1,1),
       lty=c(1,NA,NA,NA),
       pch=c(NA,15,15,15),
       pt.cex = 2,
       col=c(1,col_vec5))
dev.off()

# ++++++++++++++++++++++++++++++++++++
# Comparinson of PARX and NBARX models
# ++++++++++++++++++++++++++++++++++++

# q-q Plots both models
# ---------------------
PARX.res.out <- pseudo_residuals.P(y,PARX.filter.out$lambda.hat)
NBARX.res.out <- pseudo_residuals.NB(y,NBARX.filter.out$lambda.hat,NBARX.par.list$nu)

pdf(paste0("./plots/emp/",y_str,"_X_qq.pdf"),w=w_pdf,h=h_pdf)

qq_plot(x=PARX.res.out$res.m,type="n",pch=20,cex=0.75,conf=c("KS","Sim"),col=col_vec[1],ylim=c(-4,4))
PARX.q.l <- qq_plot(x=PARX.res.out$res.l,plot="n")
PARX.q.h <- qq_plot(x=PARX.res.out$res.h,plot="n")
segments(PARX.q.l$q_teo,PARX.q.l$q_emp,PARX.q.l$q_teo,PARX.q.h$q_emp,col=col_vec[1],lwd=2,lend=1)

qq_plot(x=NBARX.res.out$res.m,plot="n",pch=20,cex=0.75,col=col_vec[2])
NBARX.q.l <- qq_plot(x=NBARX.res.out$res.l,plot="n")
NBARX.q.h <- qq_plot(x=NBARX.res.out$res.h,plot="n")
segments(NBARX.q.l$q_teo,NBARX.q.l$q_emp,NBARX.q.l$q_teo,NBARX.q.h$q_emp,col=col_vec[2],lwd=2,lend=1)
legend(min(NBARX.q.l$q_teo),4,c("PARX","NBARX"),lwd=2,col=col_vec)

dev.off()
# Conditional variance
# --------------------
pdf(paste0("./plots/emp/",y_str,"_X_dispertion.pdf"),w=w_pdf,h=h_pdf)
par(par_plot)
par(mar=c(3,4,1,1))
plot.time(sqrt(PARX.filter.out$lambda.hat),sqrt(NBARX.filter.out$lambda.hat+NBARX.filter.out$delta.hat),lwd=2,
          ylab="Conditional std. deviation",l.pos=3,
          l.text=c(expression(paste("PARX: ", sqrt(hat(lambda)[t]))),
                   expression(paste("NBARX: ", sqrt(hat(lambda)[t]+hat(lambda)[t]^2/nu)))))
dev.off()
