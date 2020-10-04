# Following Thoning 1989
# A time series is decomposed into: a trend: n polynomial + interannual variability (long time variability)
#                                   a seasonal cycle: m harmonic + short time variability
#                                   residual term
# n polynomial + m harmonic is first fitted to the original data, the residual is then filtered by a bandpass filter to compute the short 
# time variability and filtered by a lowpass filter to compute the interannual variability (longtime variability)
# Coded by Yuming Jin (y2jin@ucsd.edu)
# 2020/07/02


thoning_fit <- function(dt,value,poly = NA,harm = NA, intercept = T, long.var.pass = 1, short.var.pass=c(1/6,1),interp=T){
  
  library(lubridate)
  library(dplR)
  library(signal)

  if(!exists("dt")){
    stop("Fitting time step missing")
  }
  if(!exists("value")){
    stop("Fitting value is missing")
  }
  
  #construct polynomial fit
  dtfit <- decimal_date(dt)* 2 * pi
  
  interp_dt <- seq(dt[1],dt[length(dt)],by="day")
  interp_date <- decimal_date(interp_dt)*2*pi
  
  if(poly < 0 | poly > 3){
    stop("0 to 3 Polynomial is required")
  }else if(poly == 1){
    polynomial <- cbind(dtfit)
    polynomial_interp <- cbind(interp_date)
  }else if(poly == 2){
    polynomial <- cbind(dtfit,dtfit^2)
    polynomial_interp <- cbind(interp_date,interp_date^2)
  }else if(poly == 3){
    polynomial <- cbind(dtfit,dtfit^2,dtfit^3)
    polynomial_interp <- cbind(interp_date,interp_date^2,interp_date^3)
  }
  
  #construct harmonic fit
  if(harm < 0 | harm > 4){
    stop("1 to 4 harmonic is required")
  }else if(harm == 1){
    harmonic <- cbind(cos(dtfit),sin(dtfit))
    harmonic_interp <- cbind(cos(interp_date),sin(interp_date))
  }else if(harm == 2){
    harmonic <- cbind(cos(dtfit),sin(dtfit),cos(2*dtfit),sin(2*dtfit))
    harmonic_interp <- cbind(cos(interp_date),sin(interp_date),cos(2*interp_date),sin(2*interp_date))
  }else if(harm == 3){
    harmonic <- cbind(cos(dtfit),sin(dtfit),cos(2*dtfit),sin(2*dtfit),cos(3*dtfit),sin(3*dtfit))
    harmonic_interp <- cbind(cos(interp_date),sin(interp_date),cos(2*interp_date),sin(2*interp_date),cos(3*interp_date),sin(3*interp_date))
  }else if(harm == 4){
    harmonic <- cbind(cos(dtfit),sin(dtfit),cos(2*dtfit),sin(2*dtfit),cos(3*dtfit),sin(3*dtfit),cos(4*dtfit),sin(4*dtfit))
    harmonic_interp <- cbind(cos(interp_date),sin(interp_date),cos(2*interp_date),sin(2*interp_date),cos(3*interp_date),sin(3*interp_date),cos(4*interp_date),sin(4*interp_date))
  }
  
  #merge polynomial harmonic and offset
  if(!exists("intercept")){
    stop('Need to specify whether to have intercept')
  }
  
  if(poly == 0){
    if(intercept){
      x <- cbind(1,harmonic)
      x_interp <- cbind(1,harmonic_interp)
    }else{
      x <- harmonic
      x_interp <- harmonic_interp
    }
  }else{
    if(intercept){
      x <- cbind(1,polynomial,harmonic)
      x_interp <- cbind(1,polynomial_interp,harmonic_interp)
    }else{
      x <- cbind(polynomial,harmonic)
      x_interp <- cbind(polynomial_interp,harmonic_interp)
    }
  }
  
  #Fit
  fit <- lm(value~x+0)
  
  if(any(is.na(fit$coefficients))){
    stop('Has NA in fit coefficient. Try Using a lower polynomial or harmonic')
  }
  #Reconstruct trend and season
  if(intercept){
    if(poly != 0){
      trend <- x[,seq(1,(poly+1))] %*% fit$coefficients[seq(1,(poly+1))]
      trend_interp <- x_interp[,seq(1,(poly+1))] %*% fit$coefficients[seq(1,(poly+1))]
    }else{
      trend <- NA
      trend_interp <- NA
    }
  }else{
    if(poly != 0){
      trend <- x[,seq(1,(poly))] %*% fit$coefficients[seq(1,(poly))] 
      trend_interp <- x_interp[,seq(1,(poly))] %*% fit$coefficients[seq(1,(poly))] 
    }else{
      trend <- NA
      trend_interp <- NA
    }
  }
  
  season <- x[,seq((ncol(x)-harm*2 + 1),ncol(x))] %*% fit$coefficients[seq((ncol(x)-harm*2 + 1),ncol(x))]
  season_interp <- x_interp[,seq((ncol(x)-harm*2 + 1),ncol(x))] %*% fit$coefficients[seq((ncol(x)-harm*2 + 1),ncol(x))]
  
  decomp <- data.frame('date'=dt,'year'=year(dt),'month'=month(dt),'DOY'=yday(dt),'value'=value,'trend'=trend,'season'=season)
  decomp_interp <- data.frame('date'=interp_dt,'year'=year(interp_dt),'month'=month(interp_dt),'DOY'=yday(interp_dt),
                              'trend'=trend_interp,'season'=season_interp)
  # compute the residual
  if(any(!is.na(trend))){
    res <- value - trend - season 
  }else{
    res <- value - season
  }
  
  if(!exists("long.var.pass")){
    stop('A cutoff frequency is required to filter Interannual Variability (Recommend 1 year)')
  }
  if(!exists("short.var.pass")){
    stop("A cutoff frequency range is required to filter short time variability (Recommend 2 - 12 months)")
  }
  
  step <- 365
  # res_interp <- interp1(x=dtfit,y=as.numeric(res),xi=interp_date,method = 'linear')
  res_interp <- approx(x=dt,y=as.numeric(res),xout=interp_dt,method = 'linear')$y
  short <- pass.filt(res_interp,n=2,W=1/step/short.var.pass, type="pass", method="Butterworth") 
  long <- pass.filt(res_interp,n=2,W=1/step/long.var.pass, type="low", method="Butterworth")
  
  decomp_interp$residual <- res_interp
  decomp_interp$short_var <- short
  decomp_interp$long_var <- long
  decomp_interp$trend_fit <- decomp_interp$trend+decomp_interp$long_var
  decomp_interp$season_fit <- decomp_interp$season + decomp_interp$short_var 
  
  decomp$short_var <- approx(x=interp_date,y=decomp_interp$short_var,xout=dtfit,rule=1,method = 'linear')$y
  decomp$long_var <- approx(x=interp_date,y=decomp_interp$long_var,xout=dtfit,rule=1,method='linear')$y
  decomp$trend_fit <- decomp$trend + decomp$long_var
  decomp$season_fit <- decomp$season+decomp$short_var
  decomp$res <- decomp$value - decomp$trend_fit - decomp$season_fit
  if(interp){
    print('Return decomposed time series with daily resolution')
    return(decomp_interp)
  }else{
    print('Return decomposed time series at original timestep')
    return(decomp)
  }
}