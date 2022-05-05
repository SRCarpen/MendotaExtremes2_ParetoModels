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
pvec = c(0.5,0.75,0.9,0.95,0.99)
quantile(dat3$zlBGA,probs=pvec)

# stationary Pareto for z-transform log BGA (z removes differences among sensors)
THR = 1  # specify threshold; 
fit1 = fevd(x=zlBGA,data=dat3,threshold=THR,type='GP')
print('',quote=F)
print('Stationary GP model --------------------------------------------------------',quote=F)
print(c('Threshold = ',THR),quote=F)
print(fit1)

# Pareto with zlBGA parameters related to P load
fit2 = fevd(x=zlBGA,data=dat3,threshold=THR,type='GP',
            scale.fun = ~ TotPload)
print('',quote=F)
print('GP model with scale or shape related to P load --------------------------------------------------------',quote=F)
print(c('Threshold = ',THR),quote=F)
print(fit2)
print('! Neither dependency on P load fits better than stationary GP',quote=F)

# Pareto with zlBGA parameters related to precip
fit3 = fevd(x=zlBGA,data=dat3,threshold=THR,type='GP',
            scale.fun = ~ PRCP)
print('',quote=F)
print('GP model with scale or shape related to Precipitation --------------------------------------------------------',quote=F)
print(c('Threshold = ',THR),quote=F)
print(fit3)
print('! Neither dependency on P load fits better than stationary GP',quote=F)
print('Model with scale dependent on PRCP has the same AIC as stationary',quote=F)
print(' but slope of PRCP effect is not significant',quote=F)

