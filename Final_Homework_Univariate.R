# UNIVARIATE ANALYSIS

#############################
####### Load packages #######
#############################
par(ask = FALSE)
library(readxl)
library(CADFtest)
library(forecast)
library(ggplot2)
library(ggfortify)
library(dplyr)

############################
####### Data loading #######
############################

# Load data
data <- data_energy <- read_excel("C:/Users/leonw/OneDrive - KU Leuven/3rd semester/Advanced Time Series/Final Homework/Data_Energy.xlsx", sheet = "Tabelle1")
data["gdp_capita"] <- data$gdp/data$population
data_emission <- read_excel("C:/Users/leonw/OneDrive - KU Leuven/3rd semester/Advanced Time Series/Final Homework/Data_Emission.xlsx")

# subset German data
data_DEU <- subset(data, country == "Germany")

###############################
####### Plot Country TS #######
###############################

Germany_renewable <- ts(data_DEU$renewables_consumption[6:59], frequency = 1, start = 1970)
Germany_emission <- ts(data_emission$CO2_DEU, frequency = 1, start = 1970)
Germany_GDP <- ts(data_DEU$gdp_capita[6:59], frequency = 1, start = 1970)

# Figure 1 in slides
par(mfrow=c(3,1))
ts.plot(Germany_renewable, gpars = list(col = c("darkblue"), lwd = 2,ylab= "Terrawatt-hours", xaxt = "n",lty=c(1)))
#legend("topleft", legend = c("Renewable energy consumption"), col = c("darkblue"), lty = c(1))

ts.plot(Germany_emission, gpars=list(col = c("firebrick"), lwd = 2, ylab= "Tons",xaxt = "n", lty=c(2)))
#legend("topright", legend = c("CO2 emission per capita"), col = c("firebrick"), lty = c(1))

ts.plot(Germany_GDP, gpars=list(col = c("forestgreen"), lwd = 2,ylab= "2011 int. $", lty=c(3)))
legend("bottomright", legend = c("Renewable cons.","CO2 emission","GDP/capita"), col = c("darkblue","firebrick","forestgreen"), lty = c(1:3))

par(mfrow = c(1,1))



#######################################################
####### TS Germany Renewable Energy Consumption #######
#######################################################

#####################################
# time series renewable energy
ts_renewable <- ts(data_DEU$renewables_consumption[6:59], frequency = 1, start = c(1970))
ts.plot(ts_renewable)
# There is a trend in the data

#####################################
# time series renewable energy
ts_logrenewable <- log(ts_renewable)
ts.plot(ts_logrenewable)
# There is a trend in the data

#####################################
# Perform test for stochastic vs. deterministic trend
max.lag <- round(sqrt(54))
CADFtest(ts_logrenewable, type = "trend", criterion = "BIC", max.lag.y = max.lag)
# I do not reject, that there is a stochastic trend.

#####################################
# Take difference
dts_logrenewable <- diff(ts_logrenewable)
ts.plot(dts_logrenewable)
# Unit root test for stationarity
CADFtest(dts_logrenewable, type = "drift", criterion = "BIC", max.lag.y = max.lag)
# Significant result thus there is stationarity

#####################################
# Correlogram and ARIMA model
par(mfrow=c(1,1))
pacf(ts_logrenewable)
acf(ts_logrenewable)

# Fit four possible ARIMA model
fit_ar_DE1 <- arima(ts_logrenewable, order = c(2,1,0),seasonal=c(0,0,0))
fit_ar_DE2 <- arima(ts_logrenewable, order = c(1,1,1),seasonal=c(0,0,0))
fit_ar_DE3 <- arima(ts_logrenewable, order = c(1,1,0),seasonal=c(0,0,0))
fit_ar_DE4 <- arima(ts_logrenewable, order = c(1,1,2),seasonal=c(0,0,0))


#####################################
# Significant?

abs(fit_ar_DE1$coef/sqrt(diag(fit_ar_DE1$var.coef)))
# Not all terms are significant. Ar2 is not significant
abs(fit_ar_DE2$coef/sqrt(diag(fit_ar_DE2$var.coef)))
# All terms are significant.
abs(fit_ar_DE3$coef/sqrt(diag(fit_ar_DE3$var.coef)))
# All terms are significant.
abs(fit_ar_DE4$coef/sqrt(diag(fit_ar_DE4$var.coef)))
# Not all terms are significant. MA2 is not significant

#####################################
# Validation

acf(fit_ar_DE1$res,plot=T,main="residual correlogram")
Box.test(fit_ar_DE1$res, lag=max.lag, type = "Ljung-Box")
# Residuals of the model behave like white noise.
# There are two borderline significant residuals at 5 and 16 
# This model is valid

acf(fit_ar_DE2$res,plot=T,main="residual correlogram")
Box.test(fit_ar_DE2$res, lag=max.lag, type = "Ljung-Box")
# Residuals of the model behave like white noise.
# There are two borderline significant residuals at 5 and 16 
# This model is valid
# check for heteroscedasticity
acf(fit_ar_DE2$residuals^2, main="squared residual correlogram")


acf(fit_ar_DE3$res,plot=T,main="residual correlogram")
Box.test(fit_ar_DE3$res, lag=max.lag, type = "Ljung-Box")
# Residuals of the model behave like white noise.
# There are two borderline significant residuals at 5 and 16 
# This model is valid but barely

acf(fit_ar_DE4$res,plot=T,main="residual correlogram")
Box.test(fit_ar_DE4$res, lag=max.lag, type = "Ljung-Box")
# Residuals of the model behave like white noise.
# There are two borderline significant residuals at 5 and 16 
# This model is valid


######################################
# Model comparisson

# AIC and BIC
AIC(fit_ar_DE1)
AIC(fit_ar_DE2)
AIC(fit_ar_DE3)
AIC(fit_ar_DE4)
AIC(fit_ar_DE1,k = log(64))
AIC(fit_ar_DE2,k = log(64))
AIC(fit_ar_DE3,k = log(64))
AIC(fit_ar_DE4,k = log(64))

# For both the AIC and BIC the the second model performs the best.
# Thus, ARMA(1,1,1) seems to be a good match.

# Expanding-window approach
# one lag forcasted
y <-ts_logrenewable
S=round(0.75*length(y))
h=1
error1.h<-c()
for (i in S:(length(y)-h))
{
  mymodel.sub<-arima(y[1:i], order = c(1,1,1),seasonal=c(0,0,0))
  predict.h<-predict(mymodel.sub,n.ahead=h)$pred[h]
  error1.h<-c(error1.h,y[i+h]-predict.h)
}
error2.h<-c()
for (i in S:(length(y)-h))
{
  mymodel.sub<-arima(y[1:i], order = c(1,1,0),seasonal=c(0,0,0))
  predict.h<-predict(mymodel.sub,n.ahead=h)$pred[h]
  error2.h<-c(error2.h,y[i+h]-predict.h)
}
error3.h<-c()
for (i in S:(length(y)-h))
{
  mymodel.sub<-arima(y[1:i], order = c(2,1,0),seasonal=c(0,0,0))
  predict.h<-predict(mymodel.sub,n.ahead=h)$pred[h]
  error3.h<-c(error3.h,y[i+h]-predict.h)
}
error4.h<-c()
for (i in S:(length(y)-h))
{
  mymodel.sub<-arima(y[1:i], order = c(1,1,2),seasonal=c(0,0,0))
  predict.h<-predict(mymodel.sub,n.ahead=h)$pred[h]
  error4.h<-c(error4.h,y[i+h]-predict.h)
}


cbind(error1.h, error2.h, error3.h, error4.h)

# Mean Absolut error for the different predictions
MAE1<-mean(abs(error1.h))
MAE2<-mean(abs(error2.h))
MAE3<-mean(abs(error3.h))
MAE4<-mean(abs(error4.h))

paste("Mean absolute error ARIMA(1,1,1): ",MAE1,
      ", ARIMA(1,1,0): ",MAE2, 
      ", ARIMA(2,1,0): ", MAE3,  
      ", ARIMA(1,1,2): ", MAE4)

# Diebold-Marino test
dm.test(error2.h, error1.h, h=h,power=1)
# Values are not significantly different
dm.test(error2.h, error3.h, h=h,power=1)
dm.test(error2.h, error4.h, h=h,power=1)



# 4 lags forcasted
y <-ts_logrenewable
S=round(0.75*length(y))
h=4
error1.h<-c()
for (i in S:(length(y)-h))
{
  mymodel.sub<-arima(y[1:i], order = c(1,1,1),seasonal=c(0,0,0))
  predict.h<-predict(mymodel.sub,n.ahead=h)$pred[h]
  error1.h<-c(error1.h,y[i+h]-predict.h)
}
error2.h<-c()
for (i in S:(length(y)-h))
{
  mymodel.sub<-arima(y[1:i], order = c(1,1,0),seasonal=c(0,0,0))
  predict.h<-predict(mymodel.sub,n.ahead=h)$pred[h]
  error2.h<-c(error2.h,y[i+h]-predict.h)
}
error3.h<-c()
for (i in S:(length(y)-h))
{
  mymodel.sub<-arima(y[1:i], order = c(2,1,0),seasonal=c(0,0,0))
  predict.h<-predict(mymodel.sub,n.ahead=h)$pred[h]
  error3.h<-c(error3.h,y[i+h]-predict.h)
}
error4.h<-c()
for (i in S:(length(y)-h))
{
  mymodel.sub<-arima(y[1:i], order = c(1,1,2),seasonal=c(0,0,0))
  predict.h<-predict(mymodel.sub,n.ahead=h)$pred[h]
  error4.h<-c(error4.h,y[i+h]-predict.h)
}


cbind(error1.h, error2.h, error3.h, error4.h)

# Mean Absolut error for the different predictions
MAE1<-mean(abs(error1.h))
MAE2<-mean(abs(error2.h))
MAE3<-mean(abs(error3.h))
MAE4<-mean(abs(error4.h))

paste("Mean absolute error ARIMA(1,1,1): ",MAE1,
      ", ARIMA(1,1,0): ",MAE2, 
      ", ARIMA(2,1,0): ", MAE3,  
      ", ARIMA(1,1,2): ", MAE4)

# Diebold-Marino test
dm.test(error2.h, error1.h, h=h,power=1)
# Values are not significantly different
dm.test(error2.h, error3.h, h=h,power=1)
dm.test(error2.h, error4.h, h=h,power=1)


############################
######## Predictions #######
############################

# prediction of the best model
myforecastAR <- predict(fit_ar_DE2, n.ahead = 8)
expectedAR <- myforecastAR$pred

# confidence bounds
lower<-myforecastAR$pred-qnorm(0.975)*myforecastAR$se
upper<-myforecastAR$pred+qnorm(0.975)*myforecastAR$se

# plot the forecasted values with 95% prediction interval
par(mfrow=c(1,1))
plot.ts(ts_logrenewable, xlim=c(2002,2035), ylim=c(4.8,8.5))
lines(expectedAR, col="red")
lines(lower, col="blue")
lines(upper, col="blue")

# Nicer plot
# Differs from given code. But as discussed in the lecture i hope this is fine.
# Source: https://cran.r-project.org/web/packages/ggfortify/vignettes/plot_ts.html
fc1 <- forecast(fit_ar_DE2, level = c(95), h = 10)
autoplot(ts_logrenewable, size = 1.5) +
  xlim(2002,2035) +
  autolayer(fc1, series="Prediction ARMA(1,1,1)", size=1.5) +
  ggtitle("Forecast from ARMA(1,1,1)") +
  xlab("Year") + ylab("log(Terrawatt-hours consumed)") +
  theme(axis.text=element_text(size=15), axis.title = element_text(size=15))+
  guides(colour=guide_legend(title="Forecast"))




