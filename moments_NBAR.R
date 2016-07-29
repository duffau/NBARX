rm(list=ls())
setwd("/home/christian/Dropbox/Speciale")
source("r_scripts/init.R")

beta.moment2 <- function(alpha,nu){
  sqrt((nu-alpha**2)/nu)-alpha
}

moment4 <- function(alpha,beta,nu){
  (alpha+beta)^4+(6*alpha^2*(alpha+beta)^2)/nu+(alpha^3*(11*alpha+8*beta))/nu^2+(6*alpha^4)/nu^3
}

# Choosing parameters
alpha <- beta <- seq(0,1,by=0.01)
nu <- 10
y_text <- 0.2
plot_text <- F

# Calculating parameter boundaries for 2nd order moments
beta_m2 <- beta.moment2(alpha,nu)
id <- beta_m2>=0
alpha_plot <- alpha[id]
beta_plot <- beta_m2[id]

pdf(paste("./plots/theo/NBAR_finte_moments_nu_",nu,".pdf"),w=w_pdf,h=h_pdf)
par(par_plot_acf)
plot(alpha_plot,beta_plot,type="l",xlim=c(0,1),ylim=c(0,1),col=col_vec[1],lwd=2,xlab=~alpha,ylab=~beta)
max_alpha_m2 <- max(alpha_plot[beta_plot>y_text],na.rm=T)
x_text <- mean(c(1-y_text,max_alpha_m2))
if(plot_text) text(x_text,0.2,~E*group("|",X,"|")<infinity,adj=c(0.5,0.5))

# Calculating parameter boundaries for 4th order moments
val <- outer(alpha,beta,function(x,y) moment4(x,y,nu))
n <- length(beta)
alpha.n <- rep(alpha,n)
n <- length(alpha)
beta.n <- rep(beta,each=n)

d <- 0.01
alpha_plot <- alpha.n[val>1-d & val<1+d]
beta_plot <- beta.n[val>1-d & val<1+d]
lines(alpha_plot,beta_plot,col=col_vec[2],lwd=2)
max_alpha_m4 <- max(alpha_plot[beta_plot>y_text],na.rm=T)
x_text <- mean(c(max_alpha_m2,max_alpha_m4))
if(plot_text) text(x_text,0.2,~E*group("|",X,"|")^2<infinity,adj=c(0.5,0.5))
x_text <- mean(c(0,max_alpha_m4))
if(plot_text) text(x_text,0.2,~E*group("|",X,"|")^4<infinity,adj=c(0.5,0.5))


segments(c(0,0),c(0,0),c(0,1),c(1,0))
lines(0:1,1:0)
dev.off()
