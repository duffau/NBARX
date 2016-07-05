# Forskellige nyttige funktioner til tidsserier i zoo-format
# ==========================================================
library(RColorBrewer)

# Utilities
# ---------
round.up <- function(x,to=10) to*(x%/%to + as.logical(x%%to))
round.do <- function(x,to=10) to*(x%/%to)

# Vstrate frstedifferens og HP-filter
# --------------------------------------
pchya <- pchya.q <- function(x) x/lag(x,-4)-1
pchya.m <- function(x) x/lag(x,-12)-1

dlog <- function(y){
  cl <- class(y)
  if(cl=="zoo"){
    log(y)-log(lag(y,-1))  
  } else {
    N <- length(y)
    log(y[2:N])-log(y[1:(N-1)])  
  }
}

pch <- function(y){
  cl <- class(y)
  if(cl=="zoo"){
    y/lag(y,-1)-1 
  } else {
    N <- length(y)
    y[2:N]/y[1:(N-1)]-1  
  }
}


hpfilt  <-  function(y,lambda=1600){
  x <- coredata(y)
  eye <- diag(length(x))
  hpout <- solve(eye+lambda*crossprod(diff(eye,d=2)),x)
  zoo(hpout,order.by=time(y))
}

# Plot of zoo-series
# ------------------
plot.time <- function(x1,x2=NULL,x3=NULL,x4=NULL,x5=NULL,x6=NULL,lwd=1,type="l",lty=1,pch=32,
                     startplot=NULL,endplot=NULL,
                     majortick="1 years",minortick="3 months",starttick=NULL,endtick=NULL,
                     xlab="",
                     zero=F, ylim=NULL, ylimR=NULL,ylab=NULL,yRlab=NULL,
                     col=1,
                     l.pos=0,l.text=NULL,l.cex=1,RL=c("L","L","L","L","L","L"),...){
  
  rapcol <- rgb(c(9,164,19,204,159,0),c(51,133,145,0,159,173),c(83,68,35,0,159,204),maxColorValue=255)
  x.indik <- logical(0)
  for(i in 1:6) x.indik[i] <- !is.null(get(paste0("x",i))) 
  x.list <- which(x.indik)
  # Repeating graphical attributes if only a single one is provided
  if(length(type)!=max(x.list)) type <- rep(type[1],max(x.list))
  if(length(lwd)!=max(x.list)) lwd <- rep(lwd[1],max(x.list))
  if(length(lty)!=max(x.list)) lty <- rep(lty[1],max(x.list))
  if(length(col)!=max(x.list)) col.vec <- brewer.pal(9,"Set1") else col.vec <- col
  if(length(pch)!=max(x.list)) pch <- rep(pch[1],max(x.list))
  
  L.varset <- which(RL=="L" & x.indik)
  R.varset <- which(RL=="R" & x.indik)
  
  L.var <- numeric(0)
  for(i in L.varset){
     assign("foo",get(paste0("x",i)))
     L.var <- c(L.var,foo) 
  } 
  
  R.var <- numeric(0)
  for(i in R.varset){
    assign("foo",get(paste0("x",i)))
    R.var <- c(R.var,foo) 
  }
  
  ymin <- if(is.null(ylim)) min(L.var,na.rm=T) else ylim[1]
  ymax <- if(is.null(ylim)) max(L.var,na.rm=T) else ylim[2]
  if(length(R.varset)>0){
  yminR <- if(is.null(ylimR)) min(R.var,na.rm=T) else ylimR[1]
  ymaxR <- if(is.null(ylimR)) max(R.var,na.rm=T) else ylimR[2]
  }
  if("zoo" %in% class(x1)){
    startplot <- if(!is.null(startplot)) {as.Date(startplot)} else as.Date(start(x1))
    endplot <- if(!is.null(endplot)) as.Date(endplot) else as.Date(end(x1))
    start <-if(!is.null(starttick)) as.Date(starttick) else if(!is.null(startplot)) as.Date(startplot) else as.Date(start(x1))
    slut <- if(!is.null(endtick)) as.Date(endtick) else if(!is.null(endplot)) as.Date(endplot) else as.Date(end(x1))
    
    for(i in setdiff(x.list,1)){
      assign("foo",get(paste0("x",i)))
      if(!"zoo" %in% class(x1)){assign(paste0("x",i),zoo(foo,order.by=time(x1)))}  
    }
      
  } else {
    startplot <- if(!is.null(startplot)) {startplot} else 0
    endplot <- if(!is.null(endplot)) {endplot} else length(x1)
    start <-if(!is.null(starttick)) starttick else if(!is.null(startplot)) startplot else 0
    slut <- if(!is.null(endtick)) endtick else if(!is.null(endplot)) endplot else length(x1)
  }
  
  plot(NULL,xaxt="n",xlim=c(startplot,endplot),ylim=c(ymin,ymax),ylab=ylab,xlab=xlab,...)
  minorticks <- seq(start, slut, by = minortick)
  majorticks <- seq(start, slut, by = majortick)
  tlab <- format(majorticks,"%Y")
  legd <- l.pos!=0 | !is.null(l.text)
  
  axis(1, at = minorticks, labels =F,tcl = -0.3)
  axis(1, at = majorticks, labels =tlab, tcl = -0.6)
  
  abline(v=majorticks,col="lightgrey",lty=3)
  grid(NA,NULL)
  if(zero) abline(h=0)
  
  for(i in L.varset){
     assign("foo",get(paste0("x",i)))  
     lines(as.Date(time(foo)),foo,col=col.vec[i],type=type[i],lty=lty[i],lwd=lwd[i],pch=pch[i],...)
  }
  
  if(length(R.varset)>0){
    par(new=T)
    plot(NULL,xaxt="n",yaxt="n",xlim=c(startplot,endplot),ylim=c(yminR,ymaxR),ann=F)
    
    for(i in R.varset){
      assign("foo",get(paste0("x",i)))  
      lines(as.Date(time(foo)),foo,col=col.vec[i],type=type[i],lty=lty[i],lwd=lwd[i],pch=pch[i],...)
    }  
    axis(4,ylab=yRlab)
  }
  
  
  par(new=T)
  plot(NULL,xaxt="n",yaxt="n",xlim=c(startplot,endplot),ylim=c(ymin,ymax),ann=F)
  posx <- startplot; posy <- ymin; xjust<-0; yjust<-0
  if(l.pos==1){posx <- startplot; posy <- ymin; xjust<-0; yjust<-0}
  if(l.pos==2){posx <- startplot; posy <- ymax; xjust<-0; yjust<-1}
  if(l.pos==3){posx <- endplot ; posy <- ymax; xjust<-1; yjust<-1}
  if(l.pos==4){posx <- endplot ; posy <- ymin; xjust<-1; yjust<-0}
  v1 <- if(!is.null(x1)) deparse(substitute(x1)) else NULL
  v2 <- if(!is.null(x2)) deparse(substitute(x2)) else NULL
  v3 <- if(!is.null(x3)) deparse(substitute(x3)) else NULL
  v4 <- if(!is.null(x4)) deparse(substitute(x4)) else NULL
  v5 <- if(!is.null(x5)) deparse(substitute(x5)) else NULL
  v6 <- if(!is.null(x6)) deparse(substitute(x6)) else NULL
  if(legd & is.null(l.text)) l.text <- c(v1,v2,v3,v4,v5,v6)
  if(legd){legend(posx,posy,l.text,lwd=lwd,lty=lty,pch=pch,col=col.vec,xjust=xjust,yjust=yjust,cex=l.cex,bg="white")}   
}

# Data fectchers
# --------------

# Daily and  monthly data
get.yahoo <- function(symb="^DJI", dates=c("2000-01-01","2016-05-30"),freq="m"){
  dates_split <- strsplit(dates, "-")
  read.zoo(paste0("http://real-chart.finance.yahoo.com/table.csv?",
                  "s=",symb,
                  "&a=",as.numeric(dates_split[[1]][2])-1,
                  "&b=",dates_split[[1]][3],
                  "&c=",dates_split[[1]][1],
                  "&d=",as.numeric(dates_split[[2]][2])-1,
                  "&e=",dates_split[[2]][3],
                  "&f=",dates_split[[2]][1],
                  "&g=",freq,
                  "&ignore=.csv"),sep=",",index.column = "Date",header=T)
}



