library(sfsmisc) # Kolmogorov-Smirnov's D critical values.
library(numDeriv) # For numerical derivatives. Used in Jacobian for standard erros in summary-function
library(zoo) # Package for irregular time series objects

# Fractional integrated process
sim.frac.x <- function(N, d=1/4, seed=123, burn_in=100){
  x <- rep(NA,N)
  set.seed(seed)
  z <- rnorm(N+burn_in)
  x[1] <- z[1]
  for(t in 2:(N+burn_in)){
    lag <- min(t-1, 100)
    x[t] <- z[t] - sum((-1)^(1:lag)*choose(d,1:lag)*x[(t-1):(t-lag)])
  }
  x[(burn_in+1):(burn_in+N)]
}

# AR process
sim.AR.x <- function(N,rho=0.8,seed=123,burn_in=100){
  set.seed(seed)
  z <- rnorm(N+burn_in)
  x[1] <- z[1] + 1/(1-rho^2)
  for(t in 2:(N+burn_in)){
    x[t]  <- rho*x[t-1] + z[t]
  }
  x[(burn_in+1):(burn_in+N)]
}

# Link functions for exogenous varaibles
# --------------------------------------
fun_x_id <- function(x,gam){
  x
}

fun_x_lin <- function(x,gam){
  x%*%gam
}

fun_x_abs <- function(x,gam){
  abs(x)%*%gam
}

fun_x_square <- function(x,gam){
  (x^2)%*%gam
}

fun_x_exp <- function(x,gam){
  exp(x%*%gam)
}

fun_x_pos <- function(x,gam){
  pmax(x%*%gam,0)
}

IHS <- function(x,theta) log(theta*x+sqrt((theta*x)^2+1))/theta


# ===================================
# Estimation and filtration functions
# ===================================

# PAR model - Fokianos, Rahbek & Tjøstheim (2009)
# -----------------------------------------------
PAR.par <- function(theta,p,q){
  list(
    "omega"=theta[1],
    "alpha"=theta[2:(1+p)],
    "beta"=theta[(2+p):(1+p+q)],
    "p"=p,
    "q"=q
    )
}

par.list2vec.PAR <- function(par.list){
  c(
    par.list$omega,
    par.list$alpha,
    par.list$beta
  )
}


sim.PAR <- function(N,theta,p,q,seed=123,lambda.out=F){
  y <- lambda <- rep(NA,N)
  par <- PAR.par(theta,p,q)
  omega <- par$omega
  alpha <- par$alpha
  beta <- par$beta
  lag_max <- max(p,q)
  
  lambda[1:lag_max] <- omega/(1-sum(alpha)-sum(beta))
  y[1:lag_max] <- rpois(lag_max,lambda[1:lag_max])
  
  if(is.null(seed)){
    for(t in (1+lag_max):N){
      lambda[t] <- max(omega + alpha%*%y[(t-1):(t-p)] + beta%*%lambda[(t-1):(t-q)],0)
      y[t] <- rpois(1,lambda[t])
    }
    
  }else{
    for(t in (1+lag_max):N){
      lambda[t] <- max(omega + alpha%*%y[(t-1):(t-p)] + beta%*%lambda[(t-1):(t-q)],0)
      set.seed(seed=seed+t)
      y[t] <- rpois(1,lambda[t])
    }
  }
  if(lambda.out){list(y=y,lambda=lambda)} else {y}
}


repam.PAR <- function(gamma,neg_par=F) {
  if(neg_par){
    c(exp(gamma[1]),(1/(1+exp(-gamma[-1]))-1/2)*2)
  } else {
    exp(gamma)
  }
}

repam.inv.PAR <- function(theta,neg_par=F) {
  if(neg_par){  
    c(log(theta[1]),-log(1/(theta[-1]*0.5 + 0.5) - 1))
  } else {
    log(theta)
  }
}

loglike.PAR <- function(gamma,p,q,y,neg_par=F){
  theta <- repam.PAR(gamma,neg_par)
  par <- PAR.par(theta,p,q)
  omega <- par$omega
  alpha <- par$alpha
  beta <- par$beta
  lag_max <- max(p,q)
  
  N <- length(y)

  lambda <- c(y[1:lag_max],rep(NA,N-lag_max))
  llikecontr <- c(rep(0,lag_max),rep(NA,N-lag_max))
  
  for(t in (1+lag_max):N){
    # lambda[t] <- omega + alpha%*%y[(t-1):(t-p)] + beta%*%lambda[(t-1):(t-q)]
    lambda[t] <- max(omega + alpha%*%y[(t-1):(t-p)] + beta%*%lambda[(t-1):(t-q)],0)
    llikecontr[t] <- log(max(dpois(y[t],lambda=lambda[t]),10^-8))
    #if(llikecontr[t]==-Inf) {print(c(lambda[t],llikecontr[t],t))}
  }
  #if(-sum(llikecontr)==-Inf){0} else {if(-sum(llikecontr)==Inf){50000} else {-sum(llikecontr)} }
  -sum(llikecontr)/N
}

filter.PAR <- function(theta,p,q,y,conf=.90,zoo=F){   
  par <- PAR.par(theta,p,q)
  omega <- par$omega
  alpha <- par$alpha
  beta <- par$beta
  lag_max <- max(p,q)
  
  N <- length(y)
  
  if(zoo){
    time <- time(y)
    y <- coredata(y)
  }
  lambda.hat <- rep(NA,N)
  l.conf <- h.conf <- matrix(NA,N,length(conf))
  colnames(l.conf) <- colnames(h.conf) <- conf
  lambda.hat[1:lag_max] <- y[1:lag_max]
  
  for(t in (1+lag_max):N){
    lambda.hat[t] <- omega + alpha%*%y[(t-1):(t-p)] + beta%*%lambda.hat[(t-1):(t-q)]
    l.conf[t,] <- qpois((1-conf)/2,lambda.hat[t])
    h.conf[t,] <- qpois(1-(1-conf)/2,lambda.hat[t])
  }
  if(zoo){
    lambda.hat <- zoo(lambda.hat,order.by=time)
    l.conf <- zoo(l.conf,order.by=time)
    h.conf <- zoo(h.conf,order.by=time)
  }
  lambda.hat[1:lag_max] <- NA
  list("lambda.hat"=lambda.hat,"l.conf"=l.conf,"h.conf"=h.conf)
}


# PARX model - Rahbek (forthcoming)
# ---------------------------------
PARX.par <- function(theta,p,q,k_x){
  list(
  "omega" = theta[1],
  "alpha" = theta[2:(1+p)],
  "beta"  = theta[(2+p):(1+p+q)],
  "gam_x" = theta[(2+p+q):(1+p+q+k_x)],
  "p" = p,
  "q" = q,
  "k_x" =k_x
  )
}

sim.PARX <- function(N,theta,p,q,x,fun_x,seed=123,all.out=F){
  k_x <- if(is.vector(x)) 1 else dim(x)[2]
  par <- PARX.par(theta,p,q,k_x)
  omega <- par$omega
  alpha <- par$alpha
  beta <- par$beta
  gam_x <- par$gam_x
  lag_max <- max(p,q)
  
  y <- lambda <- rep(NA,N)
  fx <- fun_x(x,gam_x)
  
  lambda[1:lag_max] <- (omega+mean(fx,na.rm=T))/(1-sum(alpha)-sum(beta))
  y[1:lag_max] <- rpois(lag_max,lambda[1:lag_max])

  if(is.null(seed)){
    for(t in (1+lag_max):N){
      par <- omega + alpha%*%y[(t-1):(t-p)] + beta%*%lambda[(t-1):(t-q)]
      print("----------------")
      print(paste("t =",t))
      print(paste("par =",par))
      print(paste("fx[t] =",fx[t]))
      lambda[t] <- par + fx[t]
      y[t] <- rpois(1,lambda[t])
    }
  }else{
    for(t in (1+lag_max):N){
      lambda[t] <- max(omega + alpha%*%y[(t-1):(t-p)] + beta%*%lambda[(t-1):(t-q)] + fx[t],0)
      set.seed(seed=seed+t)
      y[t] <- rpois(1,lambda[t])
    }
  }
  if(all.out){list(y=y,lambda=lambda,fx=fx)} else {y}
}

repam.PARX <- function(gamma,p,q,k_x,neg_par=F){
  if(neg_par){
    theta_par <- c(exp(gamma[1]),(1/(1+exp(-gamma[2:(1+p+q)]))-1/2)*2)
    gam_x <- gamma[(2+p+q):(1+p+q+k_x)]
    c(theta_par, gam_x)
  } else {
    theta_par <- exp(gamma[1:(1+p+q)])
    gam_x <- gamma[(2+p+q):(1+p+q+k_x)]
    c(theta_par, gam_x)
  }
} 

repam.inv.PARX <- function(theta,p,q,k_x,neg_par=F){
  if(neg_par){
    gamma_par <- c(log(theta[1]),-log(1/(theta[2:(1+p+q)]*0.5 + 0.5) - 1))
    gam_x <- theta[(2+p+q):(1+p+q+k_x)]
    c(gamma_par, gam_x)
  } else {
    gamma_par <- log(theta[1:(1+p+q)])
    gam_x <- theta[(2+p+q):(1+p+q+k_x)]
    c(gamma_par, gam_x)
  } 
}

loglike.PARX <- function(gamma,p,q,y,x,fun_x,neg_par=F,print.par=F){
  k_x <- if(is.vector(x)) 1 else dim(x)[2]
  theta <- repam.PARX(gamma,p,q,k_x,neg_par)
  if(print.par) print(theta)
  par <- PARX.par(theta,p,q,k_x)
  omega <- par$omega
  alpha <- par$alpha
  beta <- par$beta
  gam_x <- par$gam_x
  fx <- fun_x(x,gam_x)
  lag_max <- max(p,q)
  
  N <- length(y)
  
  lambda <- c(y[1:lag_max],rep(NA,N-lag_max))
  llikecontr <- c(rep(0,lag_max),rep(NA,N-lag_max))
  
  for(t in (lag_max+1):N){
    lambda[t] <- omega + alpha%*%y[(t-1):(t-p)] + beta%*%lambda[(t-1):(t-q)] + fx[t]
    # print(t)
    # print(lambda[t])
    llikecontr[t] <- log(max(dpois(y[t],lambda=lambda[t]),10^-8))
    #if(llikecontr[t]==-Inf) {print(c(lambda[t],llikecontr[t],t))}
  }
  #if(-sum(llikecontr)==-Inf){0} else {if(-sum(llikecontr)==Inf){50000} else {-sum(llikecontr)} }
  -sum(llikecontr,na.rm=T)/N
}

filter.PARX <- function(theta,p,q,y,x,fun_x,conf=.90,zoo=F){   
  k_x <- if(is.vector(x)) 1 else dim(x)[2]
  par <- PARX.par(theta,p,q,k_x)
  omega <- par$omega
  alpha <- par$alpha
  beta <- par$beta
  gam_x <- par$gam_x
  fx <- fun_x(x,gam_x)
  lag_max <- max(p,q)
  
  N <- length(y)
  if(zoo){
    time <- time(y)
    y <- coredata(y)
  }
  lambda.hat <- rep(NA,N)
  l.conf <- h.conf <- matrix(NA,N,length(conf))
  colnames(l.conf) <- colnames(h.conf) <- conf
  lambda.hat[1:lag_max]  <- y[1:lag_max]

  for(t in (1+lag_max):N){
    lambda.hat[t] <- omega + alpha%*%y[(t-1):(t-p)] + beta%*%lambda.hat[(t-1):(t-q)] + fx[t] 
    l.conf[t,] <- qpois((1-conf)/2,lambda.hat[t])
    h.conf[t,] <- qpois(1-(1-conf)/2,lambda.hat[t])
  }
  if(zoo){
    lambda.hat <- zoo(lambda.hat,order.by=time)
    fx <- zoo(fx,order.by=time)
    l.conf <- zoo(l.conf,order.by=time)
    h.conf <- zoo(h.conf,order.by=time)
  }
  lambda.hat[1:lag_max] <- NA
  list("lambda.hat"=lambda.hat,"fx"=fx,"l.conf"=l.conf,"h.conf"=h.conf)
}

# NBAR model - Christou & Fokianos (2014)
# ---------------------------------------
NBAR.par <- par.vec2list.NBAR <- function(theta,p,q){
  list(
  "omega" = theta[1],
  "alpha" = theta[2:(1+p)],
  "beta"  = theta[(2+p):(1+p+q)],
  "nu" = theta[2+p+q],
  "p" = p,
  "q" = q
  )
}

par.list2vec.NBAR <- function(par.list){
  c(
    par.list$omega,
    par.list$alpha,
    par.list$beta,
    par.list$nu
  )
}


sim.NBAR <- function(N,theta,p,q,seed=123,lambda.out=F){
  y <- lambda <- rep(NA,N)
  par <- NBAR.par(theta,p,q)
  omega <- par$omega
  alpha <- par$alpha
  beta <- par$beta
  nu <- par$nu
  lag_max <- max(p,q)
  
  lambda[1:lag_max] <- omega/(1-sum(alpha)-sum(beta))
  y[1:lag_max] <- rnbinom(1,mu=lambda[1],size=nu)
  
  if(is.null(seed)){
    for(t in (1+lag_max):N){
      lambda[t] <- max(omega + alpha%*%y[(t-1):(t-p)] + beta%*%lambda[(t-1):(t-q)],0)
      y[t] <- rnbinom(1, mu=lambda[t], size=nu)
    }
  }else{
    for(t in (1+lag_max):N){
      lambda[t] <- max(omega + alpha%*%y[(t-1):(t-p)] + beta%*%lambda[(t-1):(t-q)],0)
      set.seed(seed=seed+t)
      y[t] <- rnbinom(1, mu=lambda[t], size=nu)
    }
  }
  if(lambda.out){cbind(y,lambda)} else {y}
}

repam.NBAR <- function(gamma,neg_par=F){
  if(neg_par){
    n <- length(gamma)
    c(exp(gamma[1]),(1/(1+exp(-gamma[c(-1,-n)]))-1/2)*2,exp(gamma[n]))
  } else {
    exp(gamma)
  }
}

repam.inv.NBAR <- function(theta,neg_par=F) {
  if(neg_par){
    n <- length(theta)
    c(log(theta[1]),-log(1/(theta[c(-1,-n)]*0.5 + 0.5) - 1),log(theta[n]))
  } else {
    log(theta)
  }
}

loglike.NBAR <- function(gamma,p,q,y,neg_par=F){
  theta <- repam.NBAR(gamma,neg_par)
  par <- NBAR.par(theta,p,q)
  omega <- par$omega
  alpha <- par$alpha
  beta <- par$beta
  nu <- par$nu
  lag_max <- max(p,q)
  
  N <- length(y)
  
  lambda <- c(y[1:lag_max],rep(NA,N-lag_max))
  llikecontr <- c(rep(0,lag_max),rep(NA,N-lag_max))
  
  for(t in (lag_max+1):N){
    lambda[t] <- omega + alpha%*%y[(t-1):(t-p)] + beta%*%lambda[(t-1):(t-q)]
    llikecontr[t] <- log(max(dnbinom(y[t],mu=lambda[t],size=nu),10^-8))
    #if(llikecontr[t]==-Inf) {print(c(lambda[t],llikecontr[t],t))}
  }
  #if(-sum(llikecontr)==-Inf){0} else {if(-sum(llikecontr)==Inf){50000} else {-sum(llikecontr)} }
  -sum(llikecontr)/N
}

filter.NBAR <- function(theta,p,q,y,conf=.90,zoo=F){   
  par <- NBAR.par(theta,p,q)
  omega <- par$omega
  alpha <- par$alpha
  beta <- par$beta
  nu <- par$nu
  lag_max <- max(p,q)
  
  N <- length(y)
  if(zoo){
    time <- time(y)
    y <- coredata(y)
  }
  lambda.hat <- delta.hat <- rep(NA,N)
  l.conf <- h.conf <- matrix(NA,N,length(conf))
  colnames(l.conf) <- colnames(h.conf) <- conf
  lambda.hat[1:lag_max] <- y[1:lag_max]

  for(t in (lag_max+1):N){
    lambda.hat[t] <-  omega + alpha%*%y[(t-1):(t-p)] + beta%*%lambda.hat[(t-1):(t-q)]
    delta.hat[t] <- lambda.hat[t]^2/nu 
    l.conf[t,] <- qnbinom((1-conf)/2,mu=lambda.hat[t],size=nu)
    h.conf[t,] <- qnbinom(1-(1-conf)/2,mu=lambda.hat[t],size=nu)
  }
  if(zoo){
    lambda.hat <- zoo(lambda.hat,order.by=time)
    delta.hat <- zoo(delta.hat,order.by=time)
    l.conf <- zoo(l.conf,order.by=time)
    h.conf <- zoo(h.conf,order.by=time)
  }
  
  delta.hat[1:lag_max] <- lambda.hat[1:lag_max] <- NA  
  list("lambda.hat"=lambda.hat,"delta.hat"=delta.hat,"l.conf"=l.conf,"h.conf"=h.conf)
}




# NBARX model
# -----------
NBARX.par <- function(theta,p,q,k_x){
  list(
    "omega" = theta[1],
    "alpha" = theta[2:(1+p)],
    "beta"  = theta[(2+p):(1+p+q)],
    "gam_x" = theta[(2+p+q):(1+p+q+k_x)],
    "nu" = theta[2+p+q+k_x],
    "p" = p,
    "q" = q,
    "k_x" = k_x
    )
}

sim.NBARX <- function(N,theta,p,q,x,fun_x,seed=123,lambda.out=F){
  k_x <- if(is.vector(x)) 1 else dim(x)[2]
  par <- NBARX.par(theta,p,q,k_x)
  omega <- par$omega
  alpha <- par$alpha
  beta  <- par$beta
  gam_x <- par$gam_x
  nu <- par$nu
  lag_max <- max(p,q)
  
  y <- lambda <- rep(NA,N)
  fx <- fun_x(x,gam_x)
  
  lambda[1:lag_max] <- (omega+mean(fx,na.rm=T))/(1-sum(alpha)-sum(beta))
  y[1:lag_max] <- rnbinom(lag_max,mu=lambda[1:lag_max],size=nu)
  
  if(is.null(seed)){
    for(t in (1+lag_max):N){
      lambda[t] <- omega + alpha%*%y[(t-1):(t-p)] + beta%*%lambda[(t-1):(t-q)] + fx[t] 
      y[t] <- rnbinom(1,mu=lambda[t],size=gam)
    }
    
  }else{
    for(t in (1+lag_max):N){
      lambda[t] <- omega + alpha%*%y[(t-1):(t-p)] + beta%*%lambda[(t-1):(t-q)] + fx[t] 
      set.seed(seed=seed+t)
      y[t] <- rnbinom(1,mu=lambda[t],size=gam)
    }
  }
  if(lambda.out==T){cbind(y,lambda)} else {y}
}


repam.NBARX <- function(gamma,p,q,k_x,neg_par=F){
  if(neg_par){
    theta_par <- c(exp(gamma[1]),(1/(1+exp(-gamma[2:(1+p+q)]))-1/2)*2)
    gam_x <- gamma[(2+p+q):(1+p+q+k_x)]
    nu <- exp(gamma[2+p+q+k_x])
    c(theta_par, gam_x, nu)
  } else {
    theta_par <- exp(gamma[1:(1+p+q)])
    gam_x <- gamma[(2+p+q):(1+p+q+k_x)]
    nu <- exp(gamma[2+p+q+k_x])
    c(theta_par, gam_x, nu)
  }
} 

repam.inv.NBARX <- function(theta,p,q,k_x,neg_par=F){
  if(neg_par){
    gamma_par <- c(log(theta[1]),-log(1/(theta[2:(1+p+q)]*0.5 + 0.5) - 1))
    gam_x <- theta[(2+p+q):(1+p+q+k_x)]
    nu_inv <- log(theta[2+p+q+k_x])
    c(gamma_par, gam_x, nu_inv)
  } else {
    gamma_par <- log(theta[1:(1+p+q)])
    gam_x <- theta[(2+p+q):(1+p+q+k_x)]
    nu_inv <- log(theta[2+p+q+k_x])
    c(gamma_par, gam_x, nu_inv)
  }
} 

loglike.NBARX <- function(gamma,p,q,y,x,fun_x,neg_par=F,print.par=F){
  k_x <- if(is.vector(x)) 1 else dim(x)[2]
  theta <- repam.NBARX(gamma,p,q,k_x,neg_par)
  if(print.par) print(theta)
  par <- NBARX.par(theta,p,q,k_x)
  omega <- par$omega
  alpha <- par$alpha
  beta  <- par$beta
  gam_x <- par$gam_x
  nu <- par$nu
  lag_max <- max(p,q)
  
  N <- length(y)
  
  lambda <- c(y[1:lag_max],rep(NA,N-lag_max))
  fx <- fun_x(x,gam_x)
  llikecontr <- c(rep(0,lag_max),rep(NA,N-lag_max))
  
  for(t in (lag_max+1):N){
    lambda[t] <- omega + alpha%*%y[(t-1):(t-p)] + beta%*%lambda[(t-1):(t-q)] + fx[t] 
    llikecontr[t] <- log(max(dnbinom(y[t],mu=lambda[t],size=nu),10^-8))
    #if(llikecontr[t]==-Inf) {print(c(lambda[t],llikecontr[t],t))}
  }
  #if(-sum(llikecontr)==-Inf){0} else {if(-sum(llikecontr)==Inf){50000} else {-sum(llikecontr)} }
  -sum(llikecontr)/N
}

filter.NBARX <- function(theta,p,q,y,x,fun_x,conf=.90,zoo=F){   
  k_x <- if(is.vector(x)) 1 else dim(x)[2]
  par <- NBARX.par(theta,p,q,k_x)
  omega <- par$omega
  alpha <- par$alpha
  beta  <- par$beta
  gam_x <- par$gam_x
  nu <- par$nu
  lag_max <- max(p,q)
  
  N <- length(y)
  
  if(zoo){
    time <- time(y)
    y <- coredata(y)
  }
  lambda.hat <- delta.hat <- rep(NA,N)
  l.conf <- h.conf <- matrix(NA,N,length(conf))
  colnames(l.conf) <- colnames(h.conf) <- conf
  lambda.hat[1:lag_max] <- y[1:lag_max]
  fx <- fun_x(x,gam_x)
  
  for(t in (lag_max+1):N){
    lambda.hat[t] <- omega + alpha%*%y[(t-1):(t-p)] + beta%*%lambda.hat[(t-1):(t-q)] + fx[t]
    delta.hat[t] <- lambda.hat[t]^2/nu
    l.conf[t,] <- qnbinom((1-conf)/2,mu=lambda.hat[t],size=nu)
    h.conf[t,] <- qnbinom(1-(1-conf)/2,mu=lambda.hat[t],size=nu)
  }
  if(zoo){
    lambda.hat <- zoo(lambda.hat,order.by=time)
    delta.hat <- zoo(delta.hat,order.by=time)
    fx <- zoo(fx,order.by=time)
    l.conf <- zoo(l.conf,order.by=time)
    h.conf <- zoo(h.conf,order.by=time)
  }
  delta.hat[1:lag_max] <- lambda.hat[1:lag_max] <- NA
  
  list("lambda.hat"=lambda.hat,"fx"=fx,"delta.hat"=delta.hat,"l.conf"=l.conf,"h.conf"=h.conf)
}


# =========
# Utilities
# =========
summary.avg_mle <- function(N,optim.out,repam=NULL,par.names=NULL,...){
  par <- optim.out$par
  hessian <- optim.out$hessian*N
  if(!is.null(repam)){
    out.par <- repam(par,...)
    J <- jacobian(repam,par,...)
    inv_hessian <- chol2inv(chol(hessian))
    cov <- t(J)%*%inv_hessian%*%J
  } else {
    out.par <- par
    cov <- solve(hessian)
  }    
  std <- diag(cov)^0.5
  t.values <- out.par/std 
  p.values <- (1-pnorm(abs(t.values)))*2
  
  #   if(!is.null(outer.score)){
  #     std.robust <- diag(cov%*%outer.score%*%cov)^0.5
  #     t.robust <- out.par/std.robust 
  #     p.robust <- (1-pnorm(abs(t.robust)))*2
  #     out <- rbind("Parameter estimates"=out.par,std,t.values,p.values,std.robust,t.robust,p.robust)
  #   }
  #   else out <- rbind("Parameter estimates"=out.par,std,t.values,p.values)
  par.table <- rbind("Est."=out.par,std,t.values,p.values)
  k <- length(par)
  log_like <- -optim.out$value*N
  aic <- 2*k - 2*log_like
  bic <- k*log(N) - 2*log_like
  hq <- 2*k*log(log(N)) - 2*log_like
  info_crit.table <- rbind("N"=N,"k"=k,"Log-like"=log_like,"AIC"=aic,"BIC"=bic,"Hannan-Quinn"=hq)
  colnames(par.table) <- par.names
  list("par.table"=par.table,"info_crit.table"=info_crit.table)
}


pseudo_residuals.P <- function(y,lambda.hat){
  p.l <- ppois(pmax(y-1,0),lambda=lambda.hat)
  p.h <- ppois(y,lambda=lambda.hat)
  p.m <- (p.l+p.h)/2
  res.l <- qnorm(p.l)
  res.h <- qnorm(p.h)
  res.m <- qnorm(p.m)
  list("res.l"=res.l,"res.h"=res.h,"res.m"=res.m)
}

pseudo_residuals.NB <- function(y,lambda.hat,nu.hat){
  p.l <- pnbinom(pmax(y-1,0),mu=lambda.hat,size=nu.hat)
  p.h <- pnbinom(y,mu=lambda.hat,size=nu.hat)
  p.m <- (p.l+p.h)/2
  res.l <- qnorm(p.l)
  res.h <- qnorm(p.h)
  res.m <- qnorm(p.m)
  list("res.l"=res.l,"res.h"=res.h,"res.m"=res.m)
}

acf.PAR <- acf.NBAR <- acf.PARX <- acf.NBARX <- function(h,theta,p,q,data=NULL,plot=NULL,off_set=0,conf.alpha=NULL,...){
  if(length(h)==1){
    h_vec <- 0:h
  } else {
    h_vec <- h
    h <- max(h)
  }
  if(is.null(data)){
    alpha <- theta[2:(p+1)]
    beta  <- theta[(p+2):(p+q+1)]
    lag_max <- max(p,q)
    a <- c(alpha,rep(0,lag_max-p))
    b <- c(beta,rep(0,lag_max-q))
    acf <- ARMAacf(ar = a+b, ma = -beta, lag.max = h)[h_vec+1]
  } else {
    acf <- acf(data,h,plot=F)$acf[h_vec+1]
  }
  if(!is.null(plot)){
    if(plot=="new"){
      plot(h_vec+off_set,acf,...)
    } else if(plot=="add"){
      points(h_vec+off_set,acf,...)
    } else stop(paste("plot =",plot,"is not a valid value."))
    if(!is.null(conf.alpha)){
      conf <- qnorm(1-conf.alpha/2)/sqrt(length(y))
      abline(h=c(conf,-conf),lty=2,col=rep(greys[4:6],each=2))
      text(max(h_vec),conf,paste(conf.alpha*100,"%"),col=greys[4:6],cex=0.8)
    }
  }else{
   acf  
  }
}



hit.rate <- function(filter.out,y,N) apply(coredata(y<filter.out$l.conf)+coredata(y>filter.out$h.conf),2,function(x) sum(x,na.rm=T)/N)
hit.count <- function(filter.out,y) apply(coredata(y<filter.out$l.conf)+coredata(y>filter.out$h.conf),2,function(x) sum(x,na.rm=T))

# Kupiec (1995) - Kupiec LR-test for uncondidtional coverage
kupiec <- function(N,n,delta){
  uncon_cov <- n/N
  LR <- 2*(n*log(n/N) + (N-n)*log(1-n/N) - n*log(delta) - (N-n)*log(1-delta))
  p <- 1-pchisq(LR,1)
  out <- rbind(paste0(formatC(uncon_cov*100,1,format="f"),"%"),formatC(LR,2,format="f"),formatC(p,4,format="f"))
  colnames(out) <- paste0(formatC(delta*100,0,format="f"),"%")
  out
}


# ========================
# Q-Q plot confidens bands
# ========================

# Kolmogorov-Smirnov
KS.conf <- function(p,n,dist="normal",df=NULL){
  d <- KSd(n) # 95% confidence level
  
  if(dist=="normal"){
    conf.u <- qnorm((p+d)*(p+d<1)+1*(p+d>=1))
    conf.l <- qnorm((p-d)*(p-d>0))
  } else if(dist=="t"){
    conf.u <- qt((p+d)*(p+d<1)+1*(p+d>=1),df)*sqrt((df-2)/df)
    conf.l <- qt((p-d)*(p-d>0),df)*sqrt((df-2)/df)  
  }
  cbind(conf.l,conf.u)
}

# Aldor-Noiman et al. (2013)
sim.conf <- function(alpha=0.05,eps=0.0001,M=1000,n=200,dist="normal",df=12){
  U.sim <- matrix(runif(n*M),n,M)
  U.sort <- apply(U.sim,2,sort)
  P <- alpha-1
  a <- alpha
  while(P<1-alpha){
    a <- a-eps
    conf.u <- qbeta(1-a/2,1:n,n-1:n+1)
    conf.l <- qbeta(a/2,1:n,n-1:n+1)
    U.indik <- U.sort<matrix(conf.u,n,M) & U.sort>matrix(conf.l,n,M)
    P <- sum(U.indik)/(M*n)
  }
  
  if(dist=="normal"){
    cbind(qnorm(conf.l),qnorm(conf.u))
  } else if(dist=="t"){
    cbind(qt(conf.l,df)*sqrt((df-2)/df),qt(conf.u,df)*sqrt((df-2)/df))
  }
}

# Q-Q plot function
qq_plot <- function(x,dist="normal",grid=F,conf=NULL,mark=FALSE,alpha=0.05,plot="new",ylim=NULL,xlim=NULL,asp=NA,...){
  x <- coredata(x)
  x <- x[!is.na(x)]
  x <- x[x>-Inf & x<Inf]
  n <- length(x)
  x <- sort(x)
  p <- ppoints(n)

  if(dist=="normal"){
    q <- qnorm(p)
    xlab <-"Normal quantiles"
  } else if(dist =="t"){
    q <- qt(p)*(df-2)/df
    xlab <- "t quantiles"
  }
  
  if(is.null(xlim)){
    xlim <- range(q)
  }
  if(is.null(ylim)){
    ylim <- range(x)
  }
  
  if(plot == "new"){
    par(mar=c(4,4,1,1))
    plot(q,x,type="n",xlab=xlab,ylab="Empirical quantiles",xlim=xlim,ylim=ylim,asp=asp)
    if(grid) grid()
    lines(q,q,col="#969696")
    points(q,x,...)
  } else if(plot == "add"){
    points(q,x,...)
  }
  
  bands <- list()
  if(!is.null(conf)){
    if("KS" %in% conf){
      bands[["KS"]] <- KS.conf(p,n,dist,df)
    }
    if("Sim" %in% conf){
      bands[["Sim"]] <- sim.conf(alpha=0.05,eps=0.0001,M=1000,n,dist,df)
    }
    i <- 1
    if(plot=="new"){
      for(band in conf){
        lines(q,bands[[band]][,1],lty=i,col="#969696")
        lines(q,bands[[band]][,2],lty=i,col="#969696")
        if(mark){
          points(q[x<bands[[band]][,1]],x[x<bands[[band]][,1]],col="red")
          points(q[x>bands[[band]][,2]],x[x>bands[[band]][,2]],col="red")
        }
        i <- i + 1
      }
    legend(xlim[2],ylim[1],conf,lty=1:i,xjust=1,yjust=0,col="#969696")
    }
  }
  return(list("q_teo"=q,"q_emp"=x,"bands"=bands))
}
