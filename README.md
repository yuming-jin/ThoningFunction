# ThoningFunction
Curve Fitting Methods Applied to Time Series
(Adpated from Thoning et al. 1989 https://esrl.noaa.gov/gmd/ccgg/mbl/crvfit/crvfit.html)

------------------------------

**Required R packages**

*lubridate*

*dplR*

*signal*

------------------------------

**Defination of input variable**

dt: time 

value: data

poly: number of polynomials

harm: number of harmonics

long.var.pass: cutoff frequency to low pass filter for decomposing long-term variability (unit: year)

shot.var.pass: frequency range to band pass filter for decomposing short-term variability (unit: year)

interp: whether output interpolated data (daily resolution), default as True

------------------------------


**Example**

  mlo <- read.csv('/monthly_in_situ_co2_mlo.csv',skip=59,header = F)
  
  mlo[which(mlo[,5] < 0),5] <- NA
  
  mlo$date <- as.Date(paste(mlo[,1],mlo[,2],'16',sep='-'))
  
  colnames(mlo)[5] <- 'co2'
  
  mlo <- mlo[min(which(!is.na(mlo$co2))):max(which(!is.na(mlo$co2))),]
  
  mlo <- mlo[,c('date','co2')]
  
  mlo$date <- as.Date(mlo$date)
  
  mlo <- thoning_fit(dt=mlo$date,value=mlo$co2,poly = 3,harm = 2, intercept = T, short.var.pass=c(1/4,3/2),long.var.pass = 2,interp=T)
  
