library(readr)
dat.infl <- read_csv("D:/ѧϰ/MSc/Practicals/CompStats-Week1-TT21/infl.csv")

inflFR<-dat.infl$inflFR
year<-dat.infl$year

################################
############# EDA ##############
################################

par(mfrow=c(1,1))

plot(inflFR~year, type='l',
     xlab='Year', ylab='Rate')
points(year,inflFR, pch=16,cex=0.7)
library(zoo)
ma3 = rollmean(inflFR,3,fill = list(NA,NULL,NA))
ma7 = rollmean(inflFR,7,fill = list(NA,NA,NA,NULL,NA,NA,NA))
abline(lm(inflFR~year),lty=2)
lines(year,ma3, col='red')
lines(year,ma7, col='blue')
legend(1996,15,legend=c("Inflation Rate", "3-year Moving Average", "7-year Moving Average"),
       col=c("black","red", "blue"),lty=1,cex=0.8)
title('Real Consumer Price Annual Inflation Rate')

par(mfrow=c(1,2))
acf(inflFR, ylab='Autocorrelation',xlab='Lag',main="ACF")
pacf(inflFR,ylab='Partial Autocorrelation',xlab='Lag',main="PACF")

inflFR <- inflFR - mean(inflFR)

####################################
############## AR(1) ###############
####################################

AR <- arima(inflFR, order = c(1,0,0))
print(AR)
# intercept is the value of \mu. 
# c would be \mu*(1-\phi_1)
residuals <- residuals(AR)
AR.fitted <- inflFR - residuals

plot(dat.infl$year,inflFR, type='l',col='blue',xlab='Year',ylab='Inflation Rate')
points(year[2:65],AR.fitted[2:65], type = "l", col = 'orange', lty = 2)
legend(2000,11,legend=c("Inflation Rate", "Fitted Value"),
       col=c("blue","orange"),lty=c(1,2),cex=0.8)

par(mfrow=c(1,2))
plot(array(AR.fitted),array(residuals),type='p',
     xlab='Fitted Value',ylab='Residual')
abline(a=2,b=0,lty=2)
abline(a=-2,b=0,lty=2)

qqnorm(residuals)
qqline(residuals)
############################################################
##################### Kalman Filter ########################
############################################################


kalman = function(y, F, G, H, Q, R, mu0, Sigma0){

  T = ncol(y)
  ## INITIALIZATION ##
  
  mu.p =  array(0, T)
  Sigma.p = array(0, T)
  mu.f =  array(0, T)
  Sigma.f = array(0, T)
  
  ## FORWARD RECURSION 
  ## Time 1
  mu.p[1] = F * mu0
  Sigma.p[1] = F**2*Sigma0 + G**2 * Q
  nu1 = y[1] - H[1]*mu.p[1]
  S1 = H[1]**2 * Sigma.p[1] + R
  K1 = Sigma.p[1] * H[1] / S1
  mu.f[1] = mu.p[1] + K1  * nu1
  Sigma.f[1] = (1 - K1*H[1]) * Sigma.p[1]
  log.lik = dnorm(y[1],H[1]*mu.p[1],sqrt(S1),log=TRUE)
  # Time 2:T
  for (t in (2:T)) {
    # Prediction
    mu.p[t] = F * mu.f[t-1]
    Sigma.p[t] = F**2 * Sigma.f[t-1]+ G**2*Q
    
    # Update
    nut = y[t] - H[t] * mu.p[t]
    St = H[t]**2 * Sigma.p[t]+ R
    Kt = Sigma.p[t]*H[t]/St
    mu.f[t] = mu.p[t] + Kt * nut
    Sigma.f[t] = (1- Kt*H[t])*Sigma.p[t]
    log.lik = log.lik + dnorm(y[t],H[t]*mu.p[t],sqrt(St),log=TRUE)
  }
  
  return(list(mu.f = mu.f, Sigma.f = Sigma.f, 
              mu.p = mu.p, Sigma.p = Sigma.p,
              log.lik = log.lik))
  
}

F = 1
G = 1
Qa = 0.01
Ra = 4
Sigma0 = 1
mu0 = 0
T0 = length(inflFR)
y = t(inflFR[2:T0])
H = t(inflFR[1:(T0-1)])
  #matrix(1,nrow=1,ncol=T0-1)


results.KF = kalman(y, F, G, H, Qa, Ra, mu0, Sigma0)

mu.f = results.KF$mu.f
Sigma.f = results.KF$Sigma.f
se.f = sqrt(Sigma.f)
mu.f0 = c(mu0,mu.f)
se.f0 = c(sqrt(Sigma0),se.f)
alpha = 0.05
cv95 = qnorm(1-alpha/2)
CIupper = mu.f0+cv95*se.f0
CIlower = mu.f0-cv95*se.f0

mu.p = results.KF$mu.p
y.p = array(H)*mu.f

residual.KF <- y.p - array(y)
par(mfrow=c(1,2))
plot(y.p,residual.KF,type='p',
     xlab='Fitted Value',ylab='Residual')
abline(a=2,b=0,lty=2)
abline(a=-2,b=0,lty=2)

qqnorm(residual.KF)
qqline(residual.KF)

install.packages('latex2exp')
library(latex2exp)

plot(year[-1],y/H, col='darkgreen',pch=15,
     xlab='Year',ylab='Value',ylim=c(-2.5,3.5))
points(year[-1],mu.f,col='red',pch=16)
points(dat.infl$year,CIupper,col='blue',type="l",lty=2,lwd=2)
points(dat.infl$year,CIlower,col='blue',type="l",lty=2,lwd=2)
legend('bottomright',pch=c(15,16),col=c('darkgreen','red'),
       legend=c(expression(frac('y'[t],'y'[t-1]) ),expression(mu[t~"|"~t])) )

plot(year[-1],y,cex=1,col='darkgreen',pch=15, 
     xlab='Year',ylab='Inflation Rate')
points(year[-1],y.p,cex=1, col='red',pch=16)
points(year[-1],array(H)*mu.p,col='blue')
points(year[-1],AR.fitted[-1],col='pink')

Q.range <- c(1:100)/1000
R.range <- c(380:420)/100
lis <- matrix(nrow=4100,ncol=3)
i=1
for (q in Q.range){
  for (r in R.range){
    results.KF.qr <- kalman(y, F, G, H, q, r, mu0, Sigma0)
    log.lik.qr <- results.KF.qr$log.lik
    lis[i,] = c(q,r,log.lik.qr)
    i = i+1
    }
}
n = which.max(lis[,3])
lis[n,]


