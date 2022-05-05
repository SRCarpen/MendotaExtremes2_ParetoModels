# Experiment with Generalized Pareto for Phycocyanin daily
# SRC 2016-12-25

rm(list = ls())
graphics.off()

library(extRemes)

# load merged daily precip, P load, & buoys 2008-2021
# Note same data used for CCFs
load(file='PPT_Pload_BGAdark_2008-2021.Rdata')
print('Daily data summer stratification 2008-2021',quote=F)
print(dat3[1,])

# Summarize data
print('summary of zlBGA',quote=F)
print(summary(dat3$zlBGA))
pvec = c(0.75,0.8,0.9,0.95,0.975,0.99,0.995)
Qphyco = quantile(dat3$zlBGA,probs=pvec)
print(Qphyco)

# stationary Pareto for z-transform log BGA (z removes differences among sensors)
THR = 1  # specify threshold; 
fit0 = fevd(x=zlBGA,data=dat3,threshold=THR,type='GP')
print('',quote=F)
print('Stationary GP model --------------------------------------------------------',quote=F)
print(c('Threshold = ',THR),quote=F)
print(fit0)

# Stationary Gen. Pareto for a range of thresholds
THRvec = c(0.5,0.6,0.75,0.85,1,1.1,1.6,2.2,2.9)
NTH = length(THRvec)
parmat = matrix(0,nr=NTH,nc=2)
SEmat = matrix(0,nr=NTH,nc=2)

for(i in 1:NTH) {
  THR = THRvec[i]
  fit1 = fevd(x=zlBGA,data=dat3,threshold=THR,type='GP')
  parvec = fit1$results$par
  covpar = fit1$results$hessian
  # transform sigma; see Coles p 83
  parmat[i,1]=parvec[1] - parvec[2]*THR
  derivs = matrix(c(1,-THR),nr=2,nc=1)
  varsigma = t(derivs)%*%covpar%*%derivs
  SEmat[i,1] = sqrt(varsigma) 
  # xi
  parmat[i,2]=parvec[2]
  SEmat[i,2]=sqrt(covpar[2,2])
}

parplus = parmat+SEmat
parminus = parmat-SEmat

windows(height=9,width=5)
par(mfrow=c(2,1),mar=c(4, 4.2, 1, 2) + 0.1, cex.axis=1.5,cex.lab=1.5)
yrange = range(parminus[,1],parplus[,1])
plot(THRvec,parmat[,1],type='b',lwd=2,pch=19,col='blue',ylim=yrange,
     xlab='Threshold, kg/d',ylab='sigma',main='ZlBGA')
points(THRvec,parplus[,1],type='l',lty=2,col='blue')
points(THRvec,parminus[,1],type='l',lty=2,col='blue')
yrange = range(parminus[,2],parplus[,2])
plot(THRvec,parmat[,2],type='b',lwd=2,pch=19,col='red',ylim=yrange,
     xlab='Threshold, kg/d',ylab='xi',main='ZlBGA')
points(THRvec,parplus[,2],type='l',lty=2,col='red')
points(THRvec,parminus[,2],type='l',lty=2,col='red')

windows(height=9,width=5)
par(mfrow=c(2,1),mar=c(4, 4.2, 1, 2) + 0.1, cex.axis=1.5,cex.lab=1.5)
plot(THRvec,parmat[,1],type='b',lwd=2,pch=19,col='blue',log='x',
     xlab='Threshold, Z',ylab='sigma',main='ZlBGA')
grid()
plot(THRvec,parmat[,2],type='b',lwd=2,pch=19,col='red',log='x',
     xlab='Threshold, Z',ylab='xi',main='ZlBGA')
grid()
