# Experiment with Generalized Pareto
# SRC 2016-12-25

rm(list = ls())
graphics.off()

library(extRemes)
library(quantreg)

# Functions

# Determine if a year is a leap year and set the day number for 30 Sept
isleap = function(yr) {
  rem = (2048-yr)%%4
  leapyr = ifelse(rem==0,1,0)
  Sept30 = 273+leapyr
  return(Sept30)
}

# Make a matrix with columns of water year, row DOY, and some variable X between years Y1 and YN
# Water year is year t-1 day 274 (1 Oct) to year t day 273 (30 Sept)
WatYear = function(year,doy,X,Y1,YN) {
  NY = YN-Y1+1 # of years
  # extract first year
  sep30.0 = isleap(Y1-1)
  sep30.1 = isleap(Y1)
  # make x 
  x0 = subset(X,subset=(year==Y1-1 & doy>sep30.0))
  x1 = subset(X,subset=(year==Y1 & doy<=sep30.1))
  x.wyr = c(x0,x1)
  # make water year
  n.wyr = length(x.wyr)
  y.wyr = rep(Y1,n.wyr)
  # make water doy
  d.wyr = (1:n.wyr)
  # save calendar years
  y0 = subset(year,subset=(year==Y1-1 & doy>sep30.0))
  y1 = subset(year,subset=(year==Y1 & doy<=sep30.1))
  y.cal = c(y0,y1)
  # save calendar doys
  d0 = subset(doy,subset=(year==Y1-1 & doy>sep30.0))
  d1 = subset(doy,subset=(year==Y1 & doy<=sep30.1))
  d.cal = c(d0,d1)
  #print(c(sep30.0,sep30.1))
  #print(c(NY,length(x.wyr),length(y.wyr),length(d.wyr),length(y.cal),length(d.cal)))
  # loop over subsequent years
  for(iy in (Y1+1):YN) {
    # extract first year
    sep30.0 = isleap(iy-1)
    sep30.1 = isleap(iy)
    # make x 
    x0 = subset(X,subset=(year==iy-1 & doy>sep30.0))
    x1 = subset(X,subset=(year==iy & doy<=sep30.1))
    x.wyr = c(x.wyr,x0,x1)
    # make water year
    n.wyr = length(x0) + length(x1)
    y.iy = rep(iy,n.wyr)
    y.wyr = c(y.wyr,y.iy)
    # make water doy
    d.iy = (1:n.wyr)
    d.wyr = c(d.wyr,d.iy)
    # save calendar years
    y0 = subset(year,subset=(year==iy-1 & doy>sep30.0))
    y1 = subset(year,subset=(year==iy & doy<=sep30.1))
    y.cal = c(y.cal,y0,y1)
    # save calendar DOYs
    d0 = subset(doy,subset=(year==iy-1 & doy>sep30.0))
    d1 = subset(doy,subset=(year==iy & doy<=sep30.1))
    d.cal = c(d.cal,d0,d1)
    #print(iy)
    #print(c(sep30.0,sep30.1))
    #print(c(length(x.wyr),length(y.wyr),length(d.wyr),length(y.cal),length(d.cal)))
  }
  print('columns: water year, water doy, x, calendar year, calendar doy',quote=F)
  mat.wyr = matrix(c(y.wyr,d.wyr,x.wyr,y.cal,d.cal),nr=length(y.wyr),nc=5,byrow=F)
  return(mat.wyr)
}

# INPUT ----------------------------------------------------------------------------------------------------
# read data
#dat0 = read.csv("Madison_DCRA_precip_1940-2015.csv", stringsAsFactors=FALSE)
load('DCRA_Precip_1940-2021.Rdata')
dat0=datall  # rename for compatibility with existing code
print('Madison DCRA precip dataset',quote=F)
print(dat0[1,]) # show first line with column names

# read data
#dat1 = read.csv("Yahara_Windsor_Flow+P.csv", stringsAsFactors=FALSE)
load('AnnLoads_PB+YP_1995-2021.Rdata')
dat1 = Pload  # rename for compatibility with existing code
dat1$TotPload = dat1$PB.kg.d + dat1$YW.kg.d   # PB + YW combined loads
rm(Pload)
print('Yahara-Windsor P load dataset',quote=F)
print(dat1[1,]) # show first line with column names


# Make water year matrices and test plot ------------------------------------------------------------------
print('Water year matrices for PPT and P load',quote=F)

# make water year matrix for daily precip
PPT = WatYear(dat0$YEAR,dat0$DOY,dat0$PRCP,1995,2021)

# test plot
doy.dec = PPT[,1] + (PPT[,2]/366)
windows(width=15,height=3)
par(mfrow=c(1,1),mar=c(4, 4.2, 1, 2) + 0.1, cex.axis=1.5,cex.lab=1.5)
plot(doy.dec,PPT[,3],type='l',col='blue',xlab='year',ylab='precip')

# make water year matrix for load
Pload = WatYear(dat1$year,dat1$doy,dat1$TotPload,1995,2021)

# test plot
doy.dec = Pload[,1] + (Pload[,2]/366)
windows(width=15,height=3)
par(mfrow=c(1,1),mar=c(4, 4.2, 1, 2) + 0.1, cex.axis=1.5,cex.lab=1.5)
plot(doy.dec,Pload[,3],type='l',col='red',xlab='year',ylab='P load')

# Quantiles
print('',quote=F)
print('Quantiles of Precipitation',quote=F)
Q.ppt = quantile(PPT[,3],probs=c(0.75,0.8,0.9,0.95,0.99,0.995),na.rm=T)
print(Q.ppt)
print('Quantiles of P Load',quote=F)
Q.load = quantile(Pload[,3],probs=c(0.75,0.8,0.9,0.95,0.99,0.995),na.rm=T)
print(Q.load)

# Merge data
dPPT = as.data.frame(PPT)
colnames(dPPT) = c('WatYr','WatDoY','ppt','CalYr','CalDoY')
dPL = as.data.frame(Pload)
colnames(dPL) = c('WatYr','WatDoY','Pload','CalYr','CalDoY')
#datP0 = merge(dPPT,dPL,by=c('WatYr','WatDoY'))
datP0 = merge(dPPT,dPL,by=c('CalYr','CalDoY'))
datP1 = na.omit(datP0)
datP = subset(datP1,select=c('CalYr','CalDoY','ppt','Pload'))
datP$Ydec = datP$CalYr + (datP$CalDoY/366) # year index for plotting

# Make a data frame
#colnames(Pload) = c('WatYr','WatDoY','Pload','CalYr','CalDoY')
#colnames(PPT) = c('WatYr','WatDoY','ppt','CalYr','CalDoY')
#datmat = cbind(doy.dec,Pload,PPT[,3])
#datP = as.data.frame(datmat)
#colnames(datP) = c('doy.dec','WatYr','WatDoY','Pload','CalYr','CalDoY','ppt')

# Stationary Gen. Pareto
THR = 50  # specify threshold; 25 and 50 are reasonable, 50 has lower AIC
fit1 = fevd(x=Pload,data=datP,threshold=THR,type='GP')
print('',quote=F)
print('Stationary GP model --------------------------------------------------------',quote=F)
print(c('Threshold = ',THR),quote=F)
print(fit1)

# MODELS DEPENDING ON PRECIP FIT BETTER THAN THOSE DEPENDING ON TIME ++++++++++++++++++++++++++++++++++++++++++
print('==========================================================',quote=F)

THR = 50  # specify threshold; 25 or 50 kg/day are reasonable, 50 has lower AIC

fit6 = fevd(x=Pload,data=datP,threshold=THR,type='GP',
            scale.fun=~ppt)
print('',quote=F)
print('GP with scale depending on precip-----------------------------------------',quote=F)
print(c('Threshold = ',THR),quote=F)
print(fit6)

# save parameters
parvec = fit6$results$par

# Calculate mean exceedances per year
Y0 = min(datP$CalYr)
YN = max(datP$CalYr)
nWY = YN - Y0 + 1
exc.per.y = rep(0,nWY)
for(iy in Y0:YN) {
  yvalue = subset(datP$Pload,subset=(datP$CalYr == iy))
  xc = rep(0,length(yvalue))
  xc = ifelse(yvalue>THR,1,0)
  exc.per.y[iy-Y0+1] = sum(xc)
}

print('',quote=F)
print('Stats on Exceedances per year',quote=F)
print(summary(exc.per.y))
print(c('SD = ',sd(exc.per.y)),quote=F)
xc.mean = mean(exc.per.y)

# Compute event sizes for a series of return times versus precip
Precip.vec = (10:50)
sigmas = parvec[1] + Precip.vec*parvec[2]
Z.Nyr = function(Nyr,sigmas) 
  { THR + (sigmas/parvec[3])*( ((xc.mean*Nyr)^parvec[3]) -1) }
Z.2yr = Z.Nyr(2,sigmas)
Z.5yr = Z.Nyr(5,sigmas)
Z.10yr = Z.Nyr(10,sigmas)
Z.20yr = Z.Nyr(20,sigmas)

windows()
par(mfrow=c(1,1),mar=c(4, 4.2, 1, 2) + 0.1, cex.axis=1.5,cex.lab=1.5)
yrange = range(Z.2yr,Z.20yr)
plot(Precip.vec,Z.2yr,type='l',lwd=2,col='blue',xlab='Daily Precip',ylim=yrange,ylab='Return Level, kg/day',
     main='P Load Return Levels vs Precip')
points(Precip.vec,Z.5yr,type='l',lwd=2,col='forestgreen')
points(Precip.vec,Z.10yr,type='l',lwd=2,col='orange')
points(Precip.vec,Z.20yr,type='l',lwd=2,col='red')
legend('topleft',legend=c('20 yr','10 yr','5 yr','2 yr'),title='Return Times',
       lwd=c(2,2,2,2),col=c('red','orange','forestgreen','blue'),bty='n',cex=1.3)

# Compute return level versus return time for contrasting precip
Nyrvec = (2:20)
sigmas = parvec[1] + 10*parvec[2]
Z10mm = Z.Nyr(Nyrvec,sigmas)
sigmas = parvec[1] + 25*parvec[2]
Z25mm = Z.Nyr(Nyrvec,sigmas)
sigmas = parvec[1] + 50*parvec[2]
Z50mm = Z.Nyr(Nyrvec,sigmas)

windows()
par(mfrow=c(1,1),mar=c(4, 4.2, 1, 2) + 0.1, cex.axis=1.5,cex.lab=1.5)
yrange = range(Z10mm,Z50mm)
plot(Nyrvec,Z10mm,type='l',lwd=2,col='blue',ylim=yrange,
     xlab='Return Interval, y',ylab='Return Level, kg/day',
     main='P Load Return Level vs Return Interval')
points(Nyrvec,Z25mm,type='l',lwd=2,col='forestgreen')
points(Nyrvec,Z50mm,type='l',lwd=2,col='red')
legend('topleft',legend=c('50 mm','25 mm','10 mm'),title='Precipitation',
       lwd=c(2,2,2,2),col=c('red','forestgreen','blue'),bty='n',cex=1.3)

windows(width=12,height=6)
par(mfrow=c(1,2),mar=c(4, 4.2, 1, 2) + 0.1, cex.axis=1.5,cex.lab=1.5)
yrange = range(Z.2yr,Z.20yr)
plot(Precip.vec,Z.2yr,type='l',lwd=2,col='blue',xlab='Daily Precip',ylim=yrange,ylab='Return Level, kg/day',
     main='A. P Load Return Levels vs Precip')
points(Precip.vec,Z.5yr,type='l',lwd=2,col='forestgreen')
points(Precip.vec,Z.10yr,type='l',lwd=2,col='orange')
points(Precip.vec,Z.20yr,type='l',lwd=2,col='red')
legend('topleft',legend=c('20 yr','10 yr','5 yr','2 yr'),title='Return Times',
       lwd=c(2,2,2,2),col=c('red','orange','forestgreen','blue'),bty='n',cex=1.3)
#
yrange = range(Z10mm,Z50mm)
plot(Nyrvec,Z10mm,type='l',lwd=2,col='blue',ylim=yrange,
     xlab='Return Interval, y',ylab='Return Level, kg/day',
     main='B. P Load Return Level vs Return Interval')
points(Nyrvec,Z25mm,type='l',lwd=2,col='forestgreen')
points(Nyrvec,Z50mm,type='l',lwd=2,col='red')
legend('topleft',legend=c('50 mm','25 mm','10 mm'),title='Precipitation',
       lwd=c(2,2,2,2),col=c('red','forestgreen','blue'),bty='n',cex=1.3)
