## ---- echo=FALSE-------------------------------------------------------------------------------------------------------------------
library(tsdl)


## ----echo=FALSE--------------------------------------------------------------------------------------------------------------------
library(MASS)
k=18
femp.ts <- tsdl[[k]]
length(femp.ts)
attr(femp.ts, "subject")
attr(femp.ts, "source")
attr(femp.ts, "description")


## ----echo=FALSE--------------------------------------------------------------------------------------------------------------------
ts.plot(femp.ts, col="blue",xlab="Years From 1948-1981", ylab="Unemployment Figures (thousands)",main="Monthly US Female (16-19 Years) Unemployment Figures from 1948-1981")


## ----echo=FALSE--------------------------------------------------------------------------------------------------------------------
acf(femp.ts, lag.max = 100, main = "ACF of Original Time Series")
pacf(femp.ts, lag.max = 100, main = "PACF of Original Time Series")


## ----echo=FALSE--------------------------------------------------------------------------------------------------------------------
plot(decompose(femp.ts))


## ----echo=FALSE--------------------------------------------------------------------------------------------------------------------
#Box-Cox Transformation
library(MASS)
t = 1:length(femp.ts)
fit = lm(femp.ts ~ t)
bcTransform = boxcox(femp.ts ~ t,plotit = TRUE)
lambda = bcTransform$x[which(bcTransform$y == max(bcTransform$y))]
lambda
femp.bc = (1/lambda)*(femp.ts^lambda-1)
op <- par(mfrow = c(1,2))
#plot box-cox transformation
ts.plot(femp.ts,main = "Original data",ylab = expression(X[t]))
ts.plot(femp.bc,main = "Box-Cox transformed data", ylab = expression(Y[t]))
#log transformation
femp.tr <- log(femp.ts)
femp.log <- femp.ts^(1/3)
#plot log transformation
ts.plot(femp.ts,main = "Original data", ylab = expression(X[t]))
ts.plot(femp.tr, xlab = "Months Since Jan 1948", ylab = "log(Monthly Sales)", main = "Log Transformed Data")


## ----------------------------------------------------------------------------------------------------------------------------------
#variance of original data compared to boxcox transformed
var(femp.ts)
var(femp.bc)


## ----echo=FALSE--------------------------------------------------------------------------------------------------------------------
#ACF/PACF of boxcox transformed data
op = par(mfrow = c(1,2))
acf(femp.bc,lag.max = 60,main = "")
pacf(femp.bc,lag.max = 60,main = "")
title("Box-Cox Transformed Time Series", line = -1, outer=TRUE)


## ----echo=FALSE--------------------------------------------------------------------------------------------------------------------
# Difference at lag = 12 (cycle determined by the ACF) to remove seasonal component
y1 = diff(femp.bc, 12)
ts.plot(y1,main = "De-seasonalized Time Series",ylab = expression(nabla^{12}~nabla~Y[t]))
abline(h = 0,lty = 2)


## ----echo=FALSE--------------------------------------------------------------------------------------------------------------------
#Differenced seasonality   
op = par(mfrow = c(1,2))
acf(y1,lag.max = 60,main = "")
pacf(y1,lag.max = 60,main = "")
title("De-seasonalized Time Series", line = -1, outer=TRUE)


## ----echo=FALSE--------------------------------------------------------------------------------------------------------------------
#Variance of differenced seasonality   
var(femp.bc)
var(y1)


## ----echo=FALSE--------------------------------------------------------------------------------------------------------------------
#differenced once
y2 = diff(y1, 1)
ts.plot(y2,main = "De-trended/seasonalized Time Series",ylab = expression(nabla^{1}~nabla~Y[t]))
abline(h = 0,lty = 2)
abline(h=mean(y2), col="red")


## ----echo=FALSE--------------------------------------------------------------------------------------------------------------------
#Differenced once acf/pacf 
op = par(mfrow = c(1,2))
acf(y2,lag.max = 100,main = "")
pacf(y2,lag.max = 100,main = "")
title("De-trended Time Series", line = -1, outer=TRUE)

## ----echo=FALSE--------------------------------------------------------------------------------------------------------------------
var(femp.ts)
var(y1)
var(y2)


## ----echo=FALSE--------------------------------------------------------------------------------------------------------------------
library(qpcR)
mod1 <- arima(femp.ts, order = c(0, 1, 3), seasonal = list(order = c(0, 1, 1), period = 12), method="ML")
mod1
AICc(mod1)
#confidence interval
confint(mod1, level=0.95)
#new fit
fit1 <- arima(femp.bc, order = c(0, 1, 3), seasonal = list(order = c(0, 1, 1), period = 12), fixed=c(NA,0,NA,NA), method="ML")
fit1
AICc(fit1)


## ----echo=FALSE--------------------------------------------------------------------------------------------------------------------
mod2 <- arima(femp.ts, order = c(0, 1, 3), seasonal = list(order = c(0, 1, 2), period = 12), method="ML")
mod2
AICc(mod2)
#confidence interval
confint(mod2, level=0.95)
#new fit
fit2 <- arima(femp.bc, order = c(0, 1, 3), seasonal = list(order = c(0, 1, 2), period = 12), fixed=c(NA,0,NA,NA,0), method="ML")
fit2
AICc(fit2)


## ----echo=FALSE--------------------------------------------------------------------------------------------------------------------
mod3 <- arima(femp.ts, order = c(0, 1, 3), seasonal = list(order = c(5, 1, 1), period = 12), method="ML")
mod3
AICc(mod3)
#confidence interval
confint(mod3, level=0.95) #all sig
#new fit
fit3 <- arima(femp.bc, order = c(0, 1, 3), seasonal = list(order = c(5, 1, 1), period = 12), fixed=c(NA,0,NA,0,0,0,0,0,0), method="ML")
fit3
AICc(fit3)


## ----echo=FALSE--------------------------------------------------------------------------------------------------------------------
library("astsa")
fit1 <- arima(femp.bc, order = c(0, 1, 3), seasonal = list(order = c(0, 1, 1), period = 12), fixed=c(NA,0,NA,NA), method="ML")
sarima(femp.ts,0,1,2,0,1,1,12)
#test for independence of residuals
Box.test(residuals(fit1), type="Ljung") #H0: the model does not exhibit lack of fit
#test for
Box.test(residuals(fit1), type ="Box-Pierce")
#test for normality of residuals
shapiro.test(residuals(fit1))
#plot fitted residuals
op <- par(mfrow=c(2,2))
ts.plot(residuals(fit1),main = "Fitted Residuals SARIMA(0,1,2)(0,1,1)[12]")
par(mfrow=c(1,2),oma=c(0,0,2,0))
#Histogram
hist(residuals(fit1),main = "SARIMA(0,1,2)(0,1,1)[12] Residuals")


## ----------------------------------------------------------------------------------------------------------------------------------
#Fit to AR(0)
ar(residuals(fit1), aic = TRUE, order.max = NULL, method = c("yule-walker"))
tsdiag(arima(residuals(fit1), order=c(0,0,0)))


## ----echo=FALSE--------------------------------------------------------------------------------------------------------------------
plot.roots <- function(ar.roots=NULL, ma.roots=NULL, size=2, angles=FALSE, special=NULL, sqecial=NULL,my.pch=1,first.col="blue",second.col="red",main=NULL)
{xylims <- c(-size,size)
      omegas <- seq(0,2*pi,pi/500)
      temp <- exp(complex(real=rep(0,length(omegas)),imag=omegas))
      plot(Re(temp),Im(temp),typ="l",xlab="x",ylab="y",xlim=xylims,ylim=xylims,main=main)
      abline(v=0,lty="dotted")
      abline(h=0,lty="dotted")
      if(!is.null(ar.roots))
        {
          points(Re(1/ar.roots),Im(1/ar.roots),col=first.col,pch=my.pch)
          points(Re(ar.roots),Im(ar.roots),col=second.col,pch=my.pch)
        }
      if(!is.null(ma.roots))
        {
          points(Re(1/ma.roots),Im(1/ma.roots),pch="*",cex=1.5,col=first.col)
          points(Re(ma.roots),Im(ma.roots),pch="*",cex=1.5,col=second.col)
        }
      if(angles)
        {
          if(!is.null(ar.roots))
            {
              abline(a=0,b=Im(ar.roots[1])/Re(ar.roots[1]),lty="dotted")
              abline(a=0,b=Im(ar.roots[2])/Re(ar.roots[2]),lty="dotted")
            }
          if(!is.null(ma.roots))
            {
              sapply(1:length(ma.roots), function(j) abline(a=0,b=Im(ma.roots[j])/Re(ma.roots[j]),lty="dotted"))
            }
        }
      if(!is.null(special))
        {
          lines(Re(special),Im(special),lwd=2)
        }
      if(!is.null(special))
        {
          lines(Re(special),Im(special),lwd=2)
        }
}
#Check causality and invertibility
plot.roots(NULL,polyroot(c(1, -0.4809, -0.0937)), main="roots of ma part")
plot.roots(NULL,polyroot(c(1, -0.6889)), main="roots of sma part")


## ----echo=FALSE--------------------------------------------------------------------------------------------------------------------
#We will forecast the next 12 observations of the original data. We expect to see its confidence interval include the removed values 799 832 846 824 678.
newfemp <- femp.ts[1:403]
rfemp <- femp.ts[404:408]
femp.test=ts(rfemp, start = c(404,1))
sarima.for(newfemp, 11,0,1,2,0,1,1,12)
points(femp.test, col="blue")

