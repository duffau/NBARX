library(sfsmisc) # Kolmogorov-Smirnov's D critical values.
library(numDeriv) # For numerical derivatives. Used in Jacobian for standard erros in summary-function
library(zoo) # Package for irregular time series objects

# Fractional integrated process
sim.frac.x <- function(N, d=1/4,sigma2=1, seed=123, burn_in=100){
  x <- rep(NA,N+burn_in)
  set.seed(seed)
  z <- rnorm(N+burn_in,0,sd=sqrt(sigma2))
  psi <- c(1,choose(1:(N+burn_in-1)+d-1,1:(N+burn_in-1)))
  for(t in 1:(N+burn_in)){
    x[t] <- sum(z[t:1]*psi[1:t])
  }
  x[(burn_in+1):(burn_in+N)]
}

# AR process
sim.AR.x <- function(N,rho=0.8,sigma2=1,seed=123,burn_in=100){
  x <- rep(NA,N+burn_in)
  set.seed(seed)
  z <- rnorm(N+burn_in,0,sd=sqrt(sigma2))
  x[1] <- z[1]
  for(t in 2:(N+burn_in)){
    x[t]  <- rho*x[t-1] + z[t]
  }
  x[(burn_in+1):(burn_in+N)]
}



# Link functions for exogenous variables
# --------------------------------------
fun_x_id <- function(x,gam) x

fun_x_lin <- function(x,gam) x%*%gam

fun_x_abs <- function(x,gam) abs(x)%*%gam

fun_x_square <- function(x,gam) (x^2)%*%gam

fun_x_exp <- function(x,gam) exp(x%*%gam)

fun_x_pos <- function(x,gam) pmax(x%*%gam,0)

IHS <- function(x,theta) log(theta*x+sqrt((theta*x)^2+1))/theta


# ===================================
# Estimation and filtration functions
# ===================================

# PAR model - Fokianos, Rahbek & Tjøstheim (2009)
# -----------------------------------------------
PAR.par <- function(theta,p,q){
  p_len <- length(p)
  q_len <- length(q)
  list(
    "omega"=theta[1],
    "alpha"=theta[2:(1+p_len)],
    "beta"=theta[(2+p_len):(1+p_len+q_len)],
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
      lambda[t] <- max(omega + alpha%*%y[t-p] + beta%*%lambda[t-q],0)
      y[t] <- rpois(1,lambda[t])
    }
    
  }else{
    for(t in (1+lag_max):N){
      lambda[t] <- max(omega + alpha%*%y[t-p] + beta%*%lambda[t-q],0)
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
    lambda[t] <- max(omega + alpha%*%y[t-p] + beta%*%lambda[t-q],0)
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
    lambda.hat[t] <- omega + alpha%*%y[t-p] + beta%*%lambda.hat[t-q]
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
  p_len <- length(p)
  q_len <- length(q)
  list(
  "omega" = theta[1],
  "alpha" = theta[2:(1+p_len)],
  "beta"  = theta[(2+p_len):(1+p_len+q_len)],
  "gam_x" = theta[(2+p_len+q_len):(1+p_len+q_len+k_x)],
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
      par <- omega + alpha%*%y[t-p] + beta%*%lambda[t-q]
      print("----------------")
      print(paste("t =",t))
      print(paste("par =",par))
      print(paste("fx[t] =",fx[t]))
      lambda[t] <- par + fx[t]
      y[t] <- rpois(1,lambda[t])
    }
  }else{
    for(t in (1+lag_max):N){
      lambda[t] <- max(omega + alpha%*%y[t-p] + beta%*%lambda[t-q] + fx[t],0)
      set.seed(seed=seed+t)
      y[t] <- rpois(1,lambda[t])
    }
  }
  if(all.out){list(y=y,lambda=lambda,fx=fx)} else {y}
}

repam.PARX <- function(gamma,p,q,k_x,neg_par=F){
  p_len <- length(p)
  q_len <- length(q)
  if(neg_par){
    theta_par <- c(exp(gamma[1]),(1/(1+exp(-gamma[2:(1+p_len+q_len)]))-1/2)*2)
    gam_x <- gamma[(2+p_len+q_len):(1+p_len+q_len+k_x)]
    c(theta_par, gam_x)
  } else {
    theta_par <- exp(gamma[1:(1+p_len+q_len)])
    gam_x <- gamma[(2+p_len+q_len):(1+p_len+q_len+k_x)]
    c(theta_par, gam_x)
  }
} 

repam.inv.PARX <- function(theta,p,q,k_x,neg_par=F){
  p_len <- length(p)
  q_len <- length(q)
  if(neg_par){
    gamma_par <- c(log(theta[1]),-log(1/(theta[2:(1+p_len+q_len)]*0.5 + 0.5) - 1))
    gam_x <- theta[(2+p_len+q_len):(1+p_len+q_len+k_x)]
    c(gamma_par, gam_x)
  } else {
    gamma_par <- log(theta[1:(1+p_len+q_len)])
    gam_x <- theta[(2+p_len+q_len):(1+p_len+q_len+k_x)]
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
    lambda[t] <- omega + alpha%*%y[t-p] + beta%*%lambda[t-q] + fx[t]
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
    lambda.hat[t] <- omega + alpha%*%y[t-p] + beta%*%lambda.hat[t-q] + fx[t] 
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
  p_len <- length(p)
  q_len <- length(q)
  list(
  "omega" = theta[1],
  "alpha" = theta[2:(1+p_len)],
  "beta"  = theta[(2+p_len):(1+p_len+q_len)],
  "nu" = theta[2+p_len+q_len],
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
      lambda[t] <- max(omega + alpha%*%y[t-p] + beta%*%lambda[t-q],0)
      y[t] <- rnbinom(1, mu=lambda[t], size=nu)
    }
  }else{
    for(t in (1+lag_max):N){
      lambda[t] <- max(omega + alpha%*%y[t-p] + beta%*%lambda[t-q],0)
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
    lambda[t] <- omega + alpha%*%y[t-p] + beta%*%lambda[t-q]
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
    lambda.hat[t] <-  omega + alpha%*%y[t-p] + beta%*%lambda.hat[t-q]
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
  p_len <- length(p)
  q_len <- length(q)
  list(
    "omega" = theta[1],
    "alpha" = theta[2:(1+p_len)],
    "beta"  = theta[(2+p_len):(1+p_len+q_len)],
    "gam_x" = theta[(2+p_len+q_len):(1+p_len+q_len+k_x)],
    "nu" = theta[2+p_len+q_len+k_x],
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
      lambda[t] <- omega + alpha%*%y[t-p] + beta%*%lambda[t-q] + fx[t] 
      y[t] <- rnbinom(1,mu=lambda[t],size=nu)
    }
    
  }else{
    for(t in (1+lag_max):N){
      lambda[t] <- omega + alpha%*%y[t-p] + beta%*%lambda[t-q] + fx[t] 
      set.seed(seed=seed+t)
      y[t] <- rnbinom(1,mu=lambda[t],size=nu)
    }
  }
  if(lambda.out==T){cbind(y,lambda)} else {y}
}


repam.NBARX <- function(gamma,p,q,k_x,neg_par=F){
  p_len <- length(p)
  q_len <- length(q)
  if(neg_par){
    theta_par <- c(exp(gamma[1]),(1/(1+exp(-gamma[2:(1+p_len+q_len)]))-1/2)*2)
    gam_x <- gamma[(2+p_len+q_len):(1+p_len+q_len+k_x)]
    nu <- exp(gamma[2+p_len+q_len+k_x])
    c(theta_par, gam_x, nu)
  } else {
    theta_par <- exp(gamma[1:(1+p_len+q_len)])
    gam_x <- gamma[(2+p_len+q_len):(1+p_len+q_len+k_x)]
    nu <- exp(gamma[2+p_len+q_len+k_x])
    c(theta_par, gam_x, nu)
  }
} 

repam.inv.NBARX <- function(theta,p,q,k_x,neg_par=F){
  p_len <- length(p)
  q_len <- length(q)
  if(neg_par){
    gamma_par <- c(log(theta[1]),-log(1/(theta[2:(1+p_len+q_len)]*0.5 + 0.5) - 1))
    gam_x <- theta[(2+p_len+q_len):(1+p_len+q_len+k_x)]
    nu_inv <- log(theta[2+p_len+q_len+k_x])
    c(gamma_par, gam_x, nu_inv)
  } else {
    gamma_par <- log(theta[1:(1+p_len+q_len)])
    gam_x <- theta[(2+p_len+q_len):(1+p_len+q_len+k_x)]
    nu_inv <- log(theta[2+p_len+q_len+k_x])
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
    lambda[t] <- omega + alpha%*%y[t-p] + beta%*%lambda[t-q] + fx[t] 
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
    lambda.hat[t] <- omega + alpha%*%y[t-p] + beta%*%lambda.hat[t-q] + fx[t]
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
    # inv_hessian <- chol2inv(chol(hessian))
    inv_hessian <- solve(hessian)
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

# Autocorrelation and autocoavriance functions
# --------------------------------------------

# Autocovariance of stationary AR(1) 
acovf.ar1 <- function(h,rho,sigma2=1){
  rho^h/(1-rho^2)*sigma2
}

# Autocovariance of fractionally integrated process with -1/2 < d < 1/2.
acovf.frac <- function(h,d,sigma2=1){
  (-1)^h*gamma(1-2*d)/(gamma(1+h-d)*gamma(1-h-d))*sigma2
}

# Populate A-matrix for ARMA acf
pop_A_matrix <- function(ar_expand, ma_expand, m){
  
  # Maximum number of lags plus one
  m_1 <- m + 1
  
  # Empty matrix of same dimension as A. 
  # Used for generating upper and lower triangular-matrices
  M <- matrix(NA,m_1,m_1)
  
  # Two (m+1) x (m+1) matrices of zeroes to be populated.
  A1 <- A2 <- matrix(0,m_1,m_1) 
  
  # Making index for populating lower an upper part of A
  sym_count <- matrix(1:(m_1),m_1,m_1,byrow=T) + 0:m
  id_seq_upper <- sym_count[upper.tri(M, diag = T)]
  id_seq_lower <- sym_count[lower.tri(M, diag = F)]
  
  # The upper tringular part of A (including the diagonal) 
  A1[upper.tri(M, diag = T)]  <- -ar_expand[id_seq_upper] 
  A1 <- A1 + diag(m+1)
  A1[1,1] <- 1
  
  # The lower tringular part of A (excluding the diagonal)
  A2[lower.tri(M,diag = F)] <- -ar_expand[id_seq_lower]
  A2[,1] <- 0
  A2[lower.tri(M,diag = F)] <- A2[lower.tri(M,diag = F)] -ar_expand[sequence(m:1)+1] 
  
  # Combining the upper and lower trinagular matrices into A
  A <- A1 + A2
  return(A)
}

# Autocovaraince for ARMA where alpha and beta are parameter vectors
# from a GARCH-type recusion like the PAR, NBAR, PARX or NBARX models. 
acovf.ARMA <- function(h_vec, alpha, beta, sigma2=1){
  
  # Length of input vectors
  p <- length(alpha)
  q <- length(beta)
  m <- max(p,q+1)  # Max. number of lags
  m_1 <- m + 1
  
  # Expanding alpha and beta vectors with zeroes, to have length m.
  a <- b <- rep(0,m)
  a[1:p] <- alpha
  b[1:q] <- beta
  
  # Calculating AR and MA parameters in ARMA-representation
  ar <- a+b
  ma <- -b
  # Defining ar_0 and ma_0 as one and padding with 
  # m+1 zeroes, to have length 1+2*(m+1)
  ar_expand <- c(1,ar,rep(0,m_1))
  ma_expand <- c(1,ma,rep(0,m_1)) 
  
  # MA-infinity coefficients with built-in R-function
  psi <- ARMAtoMA(ar=a+b, ma=-beta, q)
  # Defining psi_0 as one and psi_k for k>q as zero
  psi <- c(1,psi,rep(0,m-q))
  
  # making right hand side vector b
  b <- rep(NA,m_1)
  for(i in 1:m_1) b[i] <- sum(ma_expand[i:m_1]*psi[1:(m_1-i+1)])
  
  # Populating A-matrix according to Brockwell et al. (1991)
  A <- pop_A_matrix(ar_expand,ma_expand,m)
  
  # Solving the system of m+1 linear equations
  auto_cov_m1 <- solve(A,b*sigma2)
  
  # Recursively calculating the autocovariance if h > m+1
  if(max(h_vec) > m){
    h_max <- max(h_vec)+1
    auto_cov_h <- rep(NA,h_max)
    auto_cov_h[1:m_1] <- auto_cov_m1
    for(k in (m_1+1):h_max){
      auto_cov_h[k] <- auto_cov_h[k-1:p] %*% ar[1:p]
    }
  } else {
    auto_cov_h[h_vec+1] <- auto_cov_m1
  }
  return(auto_cov_h)
}

# Wrapper for calculating and plotting acf of the count time series models
acf.INGAR <- acf.NBAR <- acf.PARX <- acf.NBARX <- function(h,theta,p,q,sigma2=1,acf_type="correlation",model="PAR",data=NULL,plot=NULL,off_set=0,conf.alpha=NULL,...){
  if(length(h)==1){
    h_vec <- 0:h
  } else {
    h_vec <- h
    h <- max(h)
  }
  if(is.null(data)){
    alpha <- rep(0,max(p))
    beta  <- rep(0,max(q))
    if(model=="PAR"){
      alpha[p] <- PAR.par(theta,p,q)$alpha
      beta[q] <- PAR.par(theta,p,q)$beta
    } else if(model=="NBAR") {
      alpha[p] <- NBAR.par(theta,p,q)$alpha
      beta[q] <- NBAR.par(theta,p,q)$beta
    } else {stop("Model type not implemented") }
    acovf <- acovf.ARMA(h_vec,alpha, beta,sigma2=sigma2)
    if(acf_type=="correlation"){
      acf <- acovf/acovf[1]
    } else if(acf_type=="covariance"){
      acf <- acovf
      print("length acf")
      print(length(acf))
      print("legnth h_vec")
      print(length(h_vec))
    } else {
      stop(paste("type =",acf_type,"is not a valid value. Valid values are 'correlation' or 'covariance'."))
    }
  } else {
    acf <- acf(data,h,plot=F,type=acf_type)$acf[h_vec+1]
  }
  if(!is.null(plot)){
    if(plot=="new"){
      plot(h_vec+off_set,acf,...)
    } else if(plot=="add"){
      points(h_vec+off_set,acf,...)
    } else stop(paste("plot =",plot,"is not a valid value. Valid values are 'new' or 'add'."))
    if(!is.null(conf.alpha)){
      conf <- qnorm(1-conf.alpha/2)/sqrt(length(y))
      abline(h=c(conf,-conf),lty=2,col=rep(greys[4:6],each=2))
      text(max(h_vec),conf,paste(conf.alpha*100,"%"),col=greys[4:6],cex=0.8)
    }
  }else{
   return(acf)  
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
