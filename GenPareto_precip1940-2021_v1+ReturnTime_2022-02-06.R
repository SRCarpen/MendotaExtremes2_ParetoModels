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
# final data has 1 NA value for PRCP
# variates:  c('YEAR','DOY','TMAX','TMIN','TAVE','PRCP','SNWD')
#save(datall,file='DCRA_precip_1940-2021.Rdata')
load(file='DCRA_precip_1940-2021.Rdata')
dat0 = datall # rename to fit the code below
print('Madison DCRA precip dataset',quote=F)
print(dat0[1,]) # show first line with column names

# Make water year matrices and test plot ------------------------------------------------------------------
print('Water year matrices for PPT and P load',quote=F)

# make water year matrix for daily precip
PPT = WatYear(dat0$YEAR,dat0$DOY,dat0$PRCP,1940,2021)

# test plot
doy.dec = PPT[,1] + (PPT[,2]/366)
windows(width=15,height=3)
par(mfrow=c(1,1),mar=c(4, 4.2, 1, 2) + 0.1, cex.axis=1.5,cex.lab=1.5)
plot(doy.dec,PPT[,3],type='l',col='blue',xlab='year',ylab='precip')

# Make a frequency plot
rppt = rank(PPT[,3])/length(PPT[,3])
windows()
par(mfrow=c(1,1),mar=c(4, 4.2, 1, 2) + 0.1, cex.axis=1.5,cex.lab=1.5)
plot(rppt,PPT[,3]+0.1,type='p',pch=19,col='blue',cex=0.5,xlab='rank prop.',ylab='Precip',log='xy')
grid(ny=15)

# Make a data frame
colnames(PPT) = c('WatYr','WatDoY','ppt','CalYr','CalDoY')
datmat = cbind(doy.dec,PPT)
datF0 = as.data.frame(datmat)
datF = na.omit(datF0)

# Proportion of rain days < 0.2 mm
nrain = length(datF$ppt)
r02 = ifelse(datF$ppt <= 0.2,0,1) 
n02 = sum(r02)
print(c('total days = ',nrain),quote=F)
print(c('days with <= 0.2 mm = ',n02,', proportion of total = ',n02/nrain),quote=F)

# Make a density plot
XD = density(datF$ppt,kernel='epanechnikov',n=256)
windows()
par(mfrow=c(1,1),mar=c(4, 4.2, 1, 2) + 0.1, cex.axis=1.5,cex.lab=1.5)
plot(XD$x,XD$y,type='l',lwd=2,col='blue',xlab='Precip mm/day',ylab='density')

# Make a density plot omitting days with precip < 0.2 mm
dat02 = subset(datF,subset=(ppt > 0.2))
XD2 = density(dat02$ppt,kernel='epanechnikov',n=256)
windows()
par(mfrow=c(1,1),mar=c(4, 4.2, 1, 2) + 0.1, cex.axis=1.5,cex.lab=1.5)
plot(XD2$x,XD2$y,type='l',lwd=2,col='blue',xlab='Precip mm/day',ylab='density',
     main='Days with ppt <= 0.2 mm are removed')

# Time-dependent Gen. Pareto
THR = 20  # specify threshold
fit3 = fevd(x=ppt,data=datF,threshold=THR,type='GP',
            scale.fun=~doy.dec)
            #shape.fun=~doy.dec)
print('',quote=F)
print('GP with scale time-dependent-----------------------------------------',quote=F)
print(c('Threshold = ',THR),quote=F)
print(fit3)
# save parameters
parvec = fit3$results$par

# Calculate mean exceedances per year
Y0 = min(datF$WatYr)
YN = max(datF$WatYr)
nWY = YN - Y0 + 1
exc.per.y = rep(0,nWY)
for(iy in Y0:YN) {
  yppt = subset(datF$ppt,subset=(datF$WatYr == iy))
  xc = rep(0,length(yppt))
  xc = ifelse(yppt>THR,1,0)
  exc.per.y[iy-Y0+1] = sum(xc)
}

print('',quote=F)
print('Stats on Exceedances per year',quote=F)
print(summary(exc.per.y))
print(c('SD = ',sd(exc.per.y)),quote=F)
xc.mean = mean(exc.per.y)

# Compute event sizes for a series of return times versus year
Yrs = (Y0:YN)
sigmas = parvec[1] + Yrs*parvec[2]
Z.Nyr = function(Nyr) { THR + (sigmas/parvec[3])*( ((xc.mean*Nyr)^parvec[3]) -1) }
Z.1yr = Z.Nyr(1)
Z.2yr = Z.Nyr(2)
Z.5yr = Z.Nyr(5)
Z.10yr = Z.Nyr(10)

windows()
par(mfrow=c(1,1),mar=c(4, 4.2, 1, 2) + 0.1, cex.axis=1.5,cex.lab=1.5)
yrange = range(Z.1yr,Z.10yr)
plot(Yrs,Z.1yr,type='l',lwd=2,col='blue',xlab='Year',ylim=yrange,ylab='Return Level, mm',
     main='Tippett et al. function')
points(Yrs,Z.2yr,type='l',lwd=2,col='forestgreen')
points(Yrs,Z.5yr,type='l',lwd=2,col='orange')
points(Yrs,Z.10yr,type='l',lwd=2,col='red')

RLs = return.level(fit3, return.period=c(2,5,10))
NT = length(doy.dec)
windows()
par(mfrow=c(1,1),mar=c(4, 4.2, 1, 2) + 0.1, cex.axis=1.5,cex.lab=1.5)
yrange = range(RLs)
plot(doy.dec[2:NT],RLs[,1],type='l',lwd=2,col='forestgreen',xlab='Year',ylim=yrange,ylab='Return Level, mm',
     main='Rainfall in 2, 5, 10 Year Storms')
points(doy.dec[2:NT],RLs[,2],type='l',lwd=2,col='orange')
points(doy.dec[2:NT],RLs[,3],type='l',lwd=2,col='red')

# Extract first-of-year values from extRemes method
RL.vec = subset(RLs,subset=(datF$WatDoY == 1))
RL.ann = matrix(RL.vec,nr=length(Yrs),nc=3,byrow=F)

windows(width=15,height=5)
par(mfrow=c(1,3),mar=c(4, 4.2, 1, 2) + 0.1, cex.axis=1.5,cex.lab=1.5)
plot(Z.2yr,RL.ann[,1],type='p',pch=19,cex=0.5,col='forestgreen',xlab='Tippett Level',
     ylab='extRemes level',main='Tippett vs extRemes')
plot(Z.5yr,RL.ann[,2],type='p',pch=19,cex=0.5,col='orange',xlab='Tippett Level',
     ylab='extRemes level',main='Tippett vs extRemes')
plot(Z.10yr,RL.ann[,3],type='p',pch=19,cex=0.5,col='red',xlab='Tippett Level',
     ylab='extRemes level',main='Tippett vs extRemes')

# Play with return levels function
#RL2 = return.level(fit3, return.period=5, do.ci=T, alpha=0.95)
#RL3 = rlevd(period=5,loc=0,scale=sigmas[30],shape=parvec[3],threshold=THR,type='GP',npy=1)
#RL4 = ci.fevd(fit3, alpha=0.05, type='return.level', return.period=5 )

# Calculate return period (return time) with Coles eq 4.12
xm = 100  # size of extreme event 
Yrs = (Y0:YN)
sigmas = parvec[1] + Yrs*parvec[2]
shape = parvec[3]
part1 = 1 + shape*((xm-THR)/sigmas)
part2 = part1^(-1/shape)
part3 = shape*part2
RT = (1/part3)/365  # divide by 365 to convert to years

windows()
par(mfrow=c(1,1),mar=c(4, 4.2, 1, 2) + 0.1, cex.axis=1.5,cex.lab=1.5)
plot(Yrs,RT,type='b',pch=19,lwd=2,col='blue',
     xlab='Year',ylab='Return Time, Years',
     main='100 mm Rain Day')
grid()

windows(width=10,height=5)
par(mfrow=c(1,2),mar=c(4, 4.2, 1, 2) + 0.1, cex.axis=1.8,cex.lab=1.8)
yrange = range(RLs)
plot(doy.dec[2:NT],RLs[,1],type='l',lwd=3,col='forestgreen',xlab='Year',ylim=yrange,ylab='Return Level, mm',
     main='A. 2, 5, & 10 Year Events')
points(doy.dec[2:NT],RLs[,2],type='l',lwd=3,col='orange')
points(doy.dec[2:NT],RLs[,3],type='l',lwd=3,col='red')
grid()
legend('topleft',legend=c('10 yr','5 yr','2 yr'),#title='Return Times',
       lwd=c(3,3,3),col=c('red','orange','forestgreen'),bty='n',cex=1.3)
#
plot(Yrs,RT,type='l',lwd=3,col='blue',
     xlab='Year',ylab='Return Time, Years',
     main='B. Days with 100 mm Precipitation')
grid()

#