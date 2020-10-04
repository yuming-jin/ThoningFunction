# ThoningFunction
Curve Fitting Methods Applied to Time Series
(Adpated from Thoning et al. 1989 https://esrl.noaa.gov/gmd/ccgg/mbl/crvfit/crvfit.html)

------------------------------

Defination of input variable

dt: time 

value: data

poly: number of polynomials

harm: number of harmonics

long.var.pass: cutoff frequency to low pass filter for decomposing long-term variability (unit: year)

shot.var.pass: frequency range to band pass filter for decomposing short-term variability (unit: year)

interp: whether output interpolated data (daily resolution), default as True

------------------------------


Example:

load()
