


#############################
####### Load packages #######
#############################

library(readxl)
library(CADFtest)
library(forecast)
library(ggplot2)
library(ggfortify)
library(dplyr)
library(vars)
library(urca)
library(BETS) # For normalization in the plot
par(ask = FALSE)
par(mfrow=c(1,1))
############################
####### Data loading #######
############################

# Load data
data_energy <- read_excel("C:/Users/leonw/OneDrive - KU Leuven/3rd semester/Advanced Time Series/Final Homework/Data_Energy.xlsx", sheet = "Tabelle1")
data_energy["gdp_capita"] <- data_energy$gdp/data_energy$population
data_emission <- read_excel("C:/Users/leonw/OneDrive - KU Leuven/3rd semester/Advanced Time Series/Final Homework/Data_Emission.xlsx")

# subset German data
data_DEU <- subset(data_energy, country == "Germany")

###########################
## Ts for Germany

## TS Germany GDP per capita
ts_GDP_cap_DEU <- ts(data_DEU$gdp_capita[6:58], frequency = 1, start = 1970)
ts_logGDP_cap_DEU <- log(ts_GDP_cap_DEU)
# Perform test for stochastic vs. deterministic trend
max.lag <- round(sqrt(length(ts_logGDP_cap_DEU)))
CADFtest(ts_logGDP_cap_DEU, type = "trend", criterion = "BIC", max.lag.y = max.lag)
# I do not reject, that there is a stochastic trend.
## Take difference
dts_logGDP_cap_DEU <- diff(ts_logGDP_cap_DEU)
# Unit root test for stationarity
CADFtest(dts_logGDP_cap_DEU, type = "drift", criterion = "BIC", max.lag.y = max.lag)
# Significant result thus there is stationarity

## TS Germany renewable energy
ts_renewable <- ts(data_DEU$renewables_share_energy[6:58], frequency = 1, start = c(1970))
ts.plot(ts_renewable)
ts_logrenewable <- log(ts_renewable)
ts.plot(ts_logrenewable)
CADFtest(ts_logrenewable, type = "trend", criterion = "BIC", max.lag.y = max.lag)
# I do not reject, that there is a stochastic trend.
## Take difference
dts_logrenewable <- diff(ts_logrenewable)
CADFtest(dts_logrenewable, type = "drift", criterion = "BIC", max.lag.y = max.lag)
# Significant result thus there is stationarity

## TS Germany CO2 emission per capita
ts_emission <- ts(data_emission$CO2_DEU, frequency = 1, start = c(1970))
ts.plot(ts_emission)
ts_logemission <- log(ts_emission)
ts.plot(ts_logemission)
CADFtest(ts_logemission, type = "trend", criterion = "BIC", max.lag.y = max.lag)
# I do not reject, that there is a stochastic trend.
## Take difference
dts_logemission <- diff(ts_logemission)
CADFtest(dts_logemission, type = "drift", criterion = "BIC", max.lag.y = max.lag)
# Significant result thus there is stationarity

##############################################
####### Plot GDP vs Energy Consumption #######
##############################################
ts_energy <- ts(data_DEU$primary_energy_consumption[6:58], frequency = 1, start = c(1970))
ts_logenergy <- log(ts_energy)

# Normalization done with a different package not included in the course
GDP_normalized <- normalize(ts_GDP_cap_DEU)
Renewable_normalized <- normalize(ts_renewable)
Energy_normalized <- normalize(ts_energy)
ts.plot(GDP_normalized, Renewable_normalized, Energy_normalized, gpars=list(col = c("darkblue", "forestgreen", "orange"), lwd = 3, xlab="Year", lty=c(1:3)))
legend("topleft", legend = c("GDP", "Renewable", "Energy"), col = c("darkblue", "forestgreen", "orange"), lty = c(1:3))
unloadNamespace("BETS")
#######################################
####### Looking for correlation #######
#######################################
ccf(ts_logrenewable,ts_logGDP_cap_DEU)
ccf(ts_logemission,ts_logGDP_cap_DEU)
ccf(ts_logrenewable,ts_logemission)

###########################
####### VAR Modell ########
###########################
#DL<N-L with L=lags 
dlogdata<-data.frame(dts_logemission, dts_logrenewable)
plot.ts(dlogdata)
names(dlogdata)<-c("dlogEmission","dlogRenew")
attach(dlogdata)
# Arbitrarily fit VAR(1)
fit_var1<-VAR(dlogdata,type="const",p=1)
summary(fit_var1)
# Only constant is significant. R-squared is only 0.01852
## Model validation
var1_residuals<-resid(fit_var1)
par(mfrow=c(2,2))
acf(var1_residuals[,1])
acf(var1_residuals[,2])
ccf(var1_residuals[,1],var1_residuals[,2])
#The model seems to be not valid as the cross-correlogram still
# shows signficant results
par(mfrow=c(1,1))

# Finidng the perfect lag
VARselect(dlogdata,lag.max=15,type="const")
fit_varautom<-VAR(dlogdata,type="const",p=15)
summary(fit_varautom)
varautom_residuals<-resid(fit_varautom)

# win.graph()
par(mfrow=c(2,2))
acf(varautom_residuals[,1])
acf(varautom_residuals[,2])
ccf(varautom_residuals[,1],varautom_residuals[,2])

# Does not look valid.

#######################################
####### Check for Cointegration #######
#######################################
# Save log transformed time series data 
logdata<-data.frame(ts_logemission, ts_logrenewable)
plot.ts(logdata)
names(logdata)<-c("logemission","logRenew")
attach(logdata)
VARselect(logdata,lag.max=15,type="const")

# See if Cointegration exists
trace_test<-ca.jo(logdata,type="trace",K=14,ecdet="const",spec="transitory")
summary(trace_test)
maxeigen_test<-ca.jo(logdata,type="eigen",K=14,ecdet="const",spec="transitory")
summary(maxeigen_test)
# both tests show cointegration!!!
# -> Use the VECM model

##########################
####### VECM Model #######
##########################

fit_vecm1<-cajorls(trace_test,r=1)
summary(fit_vecm1$rlm)

##################
# Predictions
# transform the VECM to a VAR
fit_var<-vec2var(trace_test,r=1)

myforecast<-predict(fit_var, n.ahead = 10)
par(mfrow=c(1,1))
logemission_forecast<-ts(myforecast$fcst$logemission[,1],frequency=1,start=c(2022))
logemission_lower<-ts(myforecast$fcst$logemission[,2],frequency=1,start=c(2022))
logemission_upper<-ts(myforecast$fcst$logemission[,3],frequency=1,start=c(2022))
ts.plot(logemission_forecast,logemission_lower,logemission_upper,col=c("black","red","red"))
title(main = "10-step-ahead forecast of log(Co2 Emission per capita)")

logRenew_forecast<-ts(myforecast$fcst$logRenew[,1],frequency=1,start=c(2022))
logRenew_lower<-ts(myforecast$fcst$logRenew[,2],frequency=1,start=c(2022))
logRenew_upper<-ts(myforecast$fcst$logRenew[,3],frequency=1,start=c(2022))
ts.plot(logRenew_forecast,logRenew_lower,logRenew_upper,col=c("black","red","red"))
title(main = "10-step-ahead forecast of log(Renewable Energy Consumption)")





#-----------------------------------------------------------------------------------------
###################################
######### ROBUSTNESS CHECK ########
###################################


##################################
####### Adding GDP to VECM #######
##################################

dlogdata<-data.frame(dts_logemission, dts_logrenewable, dts_logGDP_cap_DEU)
names(dlogdata)<-c("dlogEmission","dlogRenew", "dlogGDP")
attach(dlogdata)
# Arbitrarily fit VAR(1)
fit_var1<-VAR(dlogdata,type="const",p=1)
summary(fit_var1)
# Only dlogRenew.l1 is significant. R-squared is only 0.1536
## Model validation
var1_residuals<-resid(fit_var1)
par(mfrow=c(2,2))
acf(var1_residuals[,1])
acf(var1_residuals[,2])
acf(var1_residuals[,3])
ccf(var1_residuals[,1],var1_residuals[,2])
ccf(var1_residuals[,1],var1_residuals[,3])
#The model seems to be not valid as some of the cross-correlogram still
# show signficant results as well as residuals of logrenewable
par(mfrow=c(1,1))

# Select right amount of lags
VARselect(dlogdata,lag.max=11,type="const")
# Most criteria indicate a lag of 5 having the best fit
fit_varautom<-VAR(dlogdata,type="const",p=11)
summary(fit_varautom)
varautom_residuals<-resid(fit_varautom)
# Might have a long time GDP effect

# win.graph()
par(mfrow=c(2,2))
acf(varautom_residuals[,1])
acf(varautom_residuals[,2])
acf(varautom_residuals[,3])
ccf(varautom_residuals[,1],varautom_residuals[,2])
ccf(varautom_residuals[,1],varautom_residuals[,3])
# Does not look perfectly valid.

#######################################
####### Check for Cointegration #######
#######################################

logdata<-data.frame(ts_logemission, ts_logrenewable, ts_logGDP_cap_DEU)
names(logdata)<-c("logemission","logRenew", "logGDP")
attach(logdata)
VARselect(logdata,lag.max=11,type="const")

# See if Cointegration exists
trace_test<-ca.jo(logdata,type="trace",K=10,ecdet="const",spec="transitory")
summary(trace_test)
maxeigen_test<-ca.jo(logdata,type="eigen",K=10,ecdet="const",spec="transitory")
summary(maxeigen_test)
# both tests show cointegration!!! There are 2 cointegration relation
# -> Use the VECM model
fit_vecm1<-cajorls(trace_test, r=2)
summary(fit_vecm1$rlm)


##########################################
####### Effects in other countries #######
##########################################
# 3 High income countries
data_KOR <- subset(data_energy, country == "South Korea")
data_BEL <- subset(data_energy, country == "Belgium")
data_USA <- subset(data_energy, country == "United States")

# 3 upper-middle income countries
data_CHI <- subset(data_energy, country == "China")
data_COL <- subset(data_energy, country == "Colombia")
data_ZAF <- subset(data_energy, country == "South Africa")

# 3 lower-middle income countries
data_EGY <- subset(data_energy, country == "Egypt")
data_MAR <- subset(data_energy, country == "Morocco")
data_PHL <- subset(data_energy, country == "Philippines")
data_LKA <- subset(data_energy, country == "Sri Lanka")
data_VNM <- subset(data_energy, country == "Vietnam")

#----------------------------------
###############################
####### Time Series KOR #######
###############################
data_KOR["renewables_share"] <- (data_KOR$renewables_consumption/data_LKA$primary_energy_consumption)*100
ts_renewable_KOR <- ts(data_KOR$renewables_share[6:58], frequency = 1, start = c(1970))
ts_logrenewable_KOR <- log(ts_renewable_KOR)
CADFtest(ts_logrenewable_KOR, type = "trend", criterion = "BIC", max.lag.y = max.lag)
dts_logrenewable_KOR <- diff(ts_logrenewable_KOR)
CADFtest(dts_logrenewable_KOR, type = "drift", criterion = "BIC", max.lag.y = max.lag)
# Significant result thus there is stationarity

## Emission per capita
ts_emission_KOR <- ts(data_emission$CO2_KOR, frequency = 1, start = c(1970))
ts_logemission_KOR <- log(ts_emission_KOR)
CADFtest(ts_logemission_KOR, type = "trend", criterion = "BIC", max.lag.y = max.lag)
dts_logemission_KOR <- diff(ts_logemission_KOR)
CADFtest(dts_logemission_KOR, type = "drift", criterion = "BIC", max.lag.y = max.lag)
# Significant result thus there is stationarity

#######################################
####### Check for Cointegration #######
#######################################

logdata<-data.frame(ts_logemission_KOR, ts_logrenewable_KOR)
names(logdata)<-c("logemission","logRenew")
attach(logdata)
VARselect(logdata,lag.max=15,type="const")
# AIC indicates lag 14
# See if Cointegration exists
trace_test<-ca.jo(logdata,type="trace",K=14,ecdet="const",spec="transitory")
summary(trace_test)
maxeigen_test<-ca.jo(logdata,type="eigen",K=14,ecdet="const",spec="transitory")
summary(maxeigen_test)
# both tests show cointegration!!!
# -> Use the VECM model
fit_vecm1<-cajorls(trace_test, r=1)
summary(fit_vecm1$rlm)
# Significant long term effect
# significant negative renewable energy effect
# R-squared high, Adjusted R-squared a bit lower
# F-statistic significant
# Low Residual standarad error

# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.04500 -0.01353  0.00191  0.01032  0.04895 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)   
# ect1             -0.11176    0.04298  -2.600  0.02321 * 
#   logemission.dl1   0.22857    0.26847   0.851  0.41122   
# logRenew.dl1     -0.09143    0.05339  -1.713  0.11248   
# logemission.dl2   0.17275    0.31427   0.550  0.59261   
# logRenew.dl2     -0.09771    0.05474  -1.785  0.09955 . 
# logemission.dl3  -0.30200    0.24982  -1.209  0.24999   
# logRenew.dl3     -0.10988    0.05609  -1.959  0.07377 . 
# logemission.dl4   0.24109    0.18726   1.287  0.22221   
# logRenew.dl4     -0.06719    0.05720  -1.175  0.26294   
# logemission.dl5   0.11208    0.18281   0.613  0.55127   
# logRenew.dl5     -0.13006    0.05229  -2.487  0.02858 * 
#   logemission.dl6  -0.21175    0.18468  -1.147  0.27391   
# logRenew.dl6     -0.12998    0.06511  -1.996  0.06912 . 
# logemission.dl7   0.10646    0.16286   0.654  0.52563   
# logRenew.dl7     -0.07066    0.05830  -1.212  0.24889   
# logemission.dl8  -0.22962    0.16476  -1.394  0.18870   
# logRenew.dl8     -0.21932    0.05390  -4.069  0.00156 **
#   logemission.dl9  -0.31038    0.21355  -1.453  0.17175   
# logRenew.dl9     -0.19092    0.07241  -2.637  0.02171 * 
#   logemission.dl10 -0.26875    0.17762  -1.513  0.15615   
# logRenew.dl10    -0.04876    0.05174  -0.942  0.36458   
# logemission.dl11  0.03310    0.18075   0.183  0.85774   
# logRenew.dl11    -0.12185    0.04064  -2.998  0.01110 * 
#   logemission.dl12 -0.37070    0.16089  -2.304  0.03990 * 
#   logRenew.dl12    -0.08496    0.04760  -1.785  0.09955 . 
# logemission.dl13 -0.38693    0.19225  -2.013  0.06715 . 
# logRenew.dl13    -0.01601    0.03933  -0.407  0.69114   
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.03614 on 12 degrees of freedom
# Multiple R-squared:  0.8809,	Adjusted R-squared:  0.6129 
# F-statistic: 3.287 on 27 and 12 DF,  p-value: 0.01707

###############################
####### Time Series BEL #######
###############################

ts_renewable_BEL <- ts(data_BEL$renewables_share_energy[6:58], frequency = 1, start = c(1970))
ts_logrenewable_BEL <- log(ts_renewable_BEL)
CADFtest(ts_logrenewable_BEL, type = "trend", criterion = "BIC", max.lag.y = max.lag)
dts_logrenewable_BEL <- diff(ts_logrenewable_BEL)
CADFtest(dts_logrenewable_BEL, type = "drift", criterion = "BIC", max.lag.y = max.lag)
# Significant result thus there is stationarity

## Emission per capita
ts_emission_BEL <- ts(data_emission$CO2_BEL, frequency = 1, start = c(1970))
ts_logemission_BEL <- log(ts_emission_BEL)
CADFtest(ts_logemission_BEL, type = "trend", criterion = "BIC", max.lag.y = max.lag)
dts_logemission_BEL <- diff(ts_logemission_BEL)
CADFtest(dts_logemission_BEL, type = "drift", criterion = "BIC", max.lag.y = max.lag)
# Significant result thus there is stationarity

#######################################
####### Check for Cointegration #######
#######################################

logdata<-data.frame(ts_logemission_BEL, ts_logrenewable_BEL)
names(logdata)<-c("logemission","logRenew")
attach(logdata)
VARselect(logdata,lag.max=15,type="const")
# best model with K=1 but cannot specify that. --> next best model according to AIC is lag = 4
# See if Cointegration exists
trace_test<-ca.jo(logdata,type="trace",K=4,ecdet="const",spec="transitory")
summary(trace_test)
maxeigen_test<-ca.jo(logdata,type="eigen",K=4,ecdet="const",spec="transitory")
summary(maxeigen_test)
# both tests show cointegration!!!
# -> Use the VECM model
fit_vecm1<-cajorls(trace_test, r=1)
summary(fit_vecm1$rlm)
# None significant but positive effect of renewable energy consumption
# significatn long-term effect
# No significant emission lag
# low Rsquared
# F-statistic Significant

# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.09832 -0.02142  0.01220  0.02970  0.11152 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# ect1            -0.18478    0.04817  -3.836 0.000364 ***
#   logemission.dl1 -0.25478    0.12389  -2.056 0.045196 *  
#   logRenew.dl1     0.02966    0.02605   1.139 0.260521    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.04331 on 48 degrees of freedom
# Multiple R-squared:  0.2759,	Adjusted R-squared:  0.2307 
# F-statistic: 6.097 on 3 and 48 DF,  p-value: 0.001336



###########################
####### Time Series USA #######
###########################

ts_renewable_USA <- ts(data_USA$renewables_share_energy[6:58], frequency = 1, start = c(1970))
ts_logrenewable_USA <- log(ts_renewable_USA)
CADFtest(ts_logrenewable_USA, type = "trend", criterion = "BIC", max.lag.y = max.lag)
dts_logrenewable_USA <- diff(ts_logrenewable_USA)
CADFtest(dts_logrenewable_USA, type = "drift", criterion = "BIC", max.lag.y = max.lag)
# Significant result thus there is stationarity

## Emission per capita
ts_emission_USA <- ts(data_emission$CO2_USA, frequency = 1, start = c(1970))
ts_logemission_USA <- log(ts_emission_USA)
CADFtest(ts_logemission_USA, type = "trend", criterion = "BIC", max.lag.y = max.lag)
dts_logemission_USA <- diff(ts_logemission_USA)
CADFtest(dts_logemission_USA, type = "drift", criterion = "BIC", max.lag.y = max.lag)
# Significant result thus there is stationarity

#######################################
####### Check for Cointegration #######
#######################################

logdata<-data.frame(ts_logemission_USA, ts_logrenewable_USA)
names(logdata)<-c("logemission","logRenew")
attach(logdata)
VARselect(logdata,lag.max=15,type="const")
# See if Cointegration exists
trace_test<-ca.jo(logdata,type="trace",K=14,ecdet="const",spec="transitory")
summary(trace_test)
maxeigen_test<-ca.jo(logdata,type="eigen",K=14,ecdet="const",spec="transitory")
summary(maxeigen_test)
# both tests show cointegration!!!
# -> Use the VECM model
fit_vecm1<-cajorls(trace_test, r=1)
summary(fit_vecm1$rlm)
# log renew dl9 significant and negative
# log emission dl 11 significant and positiv
# High R squared, much lower adjusted R-squared
# F-statstic not significant

# Residuals:
#   Min        1Q    Median        3Q       Max 
# -0.042747 -0.011279  0.001050  0.007087  0.028176 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)  
# ect1              0.54618    0.29017   1.882   0.0843 .
# logemission.dl1  -0.37950    0.41440  -0.916   0.3778  
# logRenew.dl1     -0.13577    0.09908  -1.370   0.1957  
# logemission.dl2  -0.51973    0.33152  -1.568   0.1429  
# logRenew.dl2     -0.03814    0.09207  -0.414   0.6860  
# logemission.dl3  -0.10982    0.32545  -0.337   0.7416  
# logRenew.dl3     -0.10426    0.09327  -1.118   0.2855  
# logemission.dl4  -0.06991    0.26020  -0.269   0.7928  
# logRenew.dl4     -0.10535    0.09266  -1.137   0.2778  
# logemission.dl5  -0.03110    0.25850  -0.120   0.9062  
# logRenew.dl5     -0.04033    0.07984  -0.505   0.6226  
# logemission.dl6  -0.39919    0.25862  -1.544   0.1487  
# logRenew.dl6     -0.08773    0.08370  -1.048   0.3152  
# logemission.dl7   0.08480    0.28057   0.302   0.7676  
# logRenew.dl7      0.08753    0.07909   1.107   0.2901  
# logemission.dl8   0.09249    0.29382   0.315   0.7583  
# logRenew.dl8      0.05734    0.08227   0.697   0.4991  
# logemission.dl9  -0.11409    0.26884  -0.424   0.6788  
# logRenew.dl9     -0.17945    0.07449  -2.409   0.0330 *
#   logemission.dl10  0.02691    0.26636   0.101   0.9212  
# logRenew.dl10     0.05970    0.08158   0.732   0.4783  
# logemission.dl11  0.62582    0.26467   2.365   0.0358 *
#   logRenew.dl11     0.08865    0.07610   1.165   0.2667  
# logemission.dl12 -0.42740    0.31650  -1.350   0.2018  
# logRenew.dl12    -0.10745    0.06777  -1.585   0.1388  
# logemission.dl13 -0.43751    0.33950  -1.289   0.2218  
# logRenew.dl13    -0.02274    0.07450  -0.305   0.7654  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.02763 on 12 degrees of freedom
# Multiple R-squared:  0.8028,	Adjusted R-squared:  0.3592 
# F-statistic:  1.81 on 27 and 12 DF,  p-value: 0.1402



###########################
####### Time Series CHI #######
###########################

ts_renewable_CHI <- ts(data_CHI$renewables_share_energy[6:58], frequency = 1, start = c(1970))
ts_logrenewable_CHI <- log(ts_renewable_CHI)
CADFtest(ts_logrenewable_CHI, type = "trend", criterion = "BIC", max.lag.y = max.lag)
dts_logrenewable_CHI <- diff(ts_logrenewable_CHI)
CADFtest(dts_logrenewable_CHI, type = "drift", criterion = "BIC", max.lag.y = max.lag)
# Significant result thus there is stationarity

## Emission per capita
ts_emission_CHI <- ts(data_emission$CO2_CHN, frequency = 1, start = c(1970))
ts_logemission_CHI <- log(ts_emission_CHI)
CADFtest(ts_logemission_CHI, type = "trend", criterion = "BIC", max.lag.y = max.lag)
dts_logemission_CHI <- diff(ts_logemission_CHI)
CADFtest(dts_logemission_CHI, type = "drift", criterion = "BIC", max.lag.y = max.lag)
# Significant result thus there is stationarity

#######################################
####### Check for Cointegration #######
#######################################

logdata<-data.frame(ts_logemission_CHI, ts_logrenewable_CHI)
names(logdata)<-c("logemission","logRenew")
attach(logdata)
VARselect(logdata,lag.max=15,type="const")
# See if Cointegration exists
trace_test<-ca.jo(logdata,type="trace",K=2,ecdet="const",spec="transitory")
summary(trace_test)
maxeigen_test<-ca.jo(logdata,type="eigen",K=2,ecdet="const",spec="transitory")
summary(maxeigen_test)
# both tests show cointegration!!!
# -> Use the VECM model
fit_vecm1<-cajorls(trace_test, r=1)
summary(fit_vecm1$rlm)
# renew not significant and very small negative effect
# emission significant and positive effect
# small and significant long-term effect
# F-statistic significant

# Residuals:
#   Min        1Q    Median        3Q       Max 
# -0.071847 -0.023745  0.001352  0.032907  0.112615 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)  
# ect1             0.030539   0.013664   2.235   0.0301 *
#   logemission.dl1  0.388031   0.166369   2.332   0.0239 *
#   logRenew.dl1    -0.002399   0.078180  -0.031   0.9757  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.04183 on 48 degrees of freedom
# Multiple R-squared:  0.5331,	Adjusted R-squared:  0.504 
# F-statistic: 18.27 on 3 and 48 DF,  p-value: 4.795e-08



###########################
####### Time Series COL #######
###########################

ts_renewable_COL <- ts(data_COL$renewables_share_energy[6:58], frequency = 1, start = c(1970))
ts_logrenewable_COL <- log(ts_renewable_COL)
CADFtest(ts_logrenewable_COL, type = "trend", criterion = "BIC", max.lag.y = max.lag)
dts_logrenewable_COL <- diff(ts_logrenewable_COL)
CADFtest(dts_logrenewable_COL, type = "drift", criterion = "BIC", max.lag.y = max.lag)
# Significant result thus there is stationarity

## Emission per capita
ts_emission_COL <- ts(data_emission$CO2_COL, frequency = 1, start = c(1970))
ts_logemission_COL <- log(ts_emission_COL)
CADFtest(ts_logemission_COL, type = "trend", criterion = "BIC", max.lag.y = max.lag)
dts_logemission_COL <- diff(ts_logemission_COL)
CADFtest(dts_logemission_COL, type = "drift", criterion = "BIC", max.lag.y = max.lag)
# Significant result thus there is stationarity

#######################################
####### Check for Cointegration #######
#######################################

logdata<-data.frame(ts_logemission_COL, ts_logrenewable_COL)
names(logdata)<-c("logemission","logRenew")
attach(logdata)
VARselect(logdata,lag.max=15,type="const")
# See if Cointegration exists
# I cant specify K=1 --> next best model according to AIC, K=14 
trace_test<-ca.jo(logdata,type="trace",K=14,ecdet="const",spec="transitory")
summary(trace_test)
maxeigen_test<-ca.jo(logdata,type="eigen",K=2,ecdet="const",spec="transitory")
summary(maxeigen_test)
# both tests show cointegration!!!
# -> Use the VECM model
fit_vecm1<-cajorls(trace_test, r=1)
summary(fit_vecm1$rlm)
# High R-squared, slightly smaller Adjusted R-squared
# Significant long term effects. Multiple positive and significant effects of emissions
# Multiple negative and signifiant effects of Renewable energy consumption
# F-statistic significant

# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.04070 -0.01082  0.00187  0.01031  0.03561 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)   
# ect1             -1.36686    0.45470  -3.006  0.01094 * 
#   logemission.dl1   0.83177    0.34887   2.384  0.03450 * 
#   logRenew.dl1     -0.69644    0.24022  -2.899  0.01335 * 
#   logemission.dl2   1.04492    0.27202   3.841  0.00235 **
#   logRenew.dl2     -0.36211    0.23288  -1.555  0.14593   
# logemission.dl3   1.06919    0.31750   3.368  0.00560 **
#   logRenew.dl3     -0.36146    0.17313  -2.088  0.05881 . 
# logemission.dl4   0.82914    0.37222   2.228  0.04581 * 
#   logRenew.dl4      0.03150    0.11413   0.276  0.78727   
# logemission.dl5   0.85860    0.40889   2.100  0.05756 . 
# logRenew.dl5     -0.08861    0.09675  -0.916  0.37778   
# logemission.dl6   0.88108    0.48104   1.832  0.09193 . 
# logRenew.dl6     -0.05151    0.08276  -0.622  0.54529   
# logemission.dl7   1.50994    0.46028   3.280  0.00657 **
#   logRenew.dl7      0.26187    0.08679   3.017  0.01071 * 
#   logemission.dl8   1.00834    0.39361   2.562  0.02492 * 
#   logRenew.dl8      0.22588    0.10184   2.218  0.04661 * 
#   logemission.dl9   0.88232    0.33906   2.602  0.02313 * 
#   logRenew.dl9      0.17239    0.09928   1.736  0.10808   
# logemission.dl10  0.43252    0.34187   1.265  0.22983   
# logRenew.dl10     0.17142    0.11311   1.515  0.15554   
# logemission.dl11  0.85047    0.34881   2.438  0.03126 * 
#   logRenew.dl11     0.33780    0.13686   2.468  0.02959 * 
#   logemission.dl12  0.49655    0.35346   1.405  0.18543   
# logRenew.dl12     0.29280    0.14769   1.983  0.07079 . 
# logemission.dl13  0.53550    0.27623   1.939  0.07643 . 
# logRenew.dl13     0.02182    0.12078   0.181  0.85967   
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.02817 on 12 degrees of freedom
# Multiple R-squared:  0.8659,	Adjusted R-squared:  0.564 
# F-statistic: 2.869 on 27 and 12 DF,  p-value: 0.02929



###########################
####### Time Series ZAF #######
###########################

ts_renewable_ZAF <- ts(data_ZAF$renewables_share_energy[6:58], frequency = 1, start = c(1970))
ts_logrenewable_ZAF <- log(ts_renewable_ZAF)
CADFtest(ts_logrenewable_ZAF, type = "trend", criterion = "BIC", max.lag.y = max.lag)
dts_logrenewable_ZAF <- diff(ts_logrenewable_ZAF)
CADFtest(dts_logrenewable_ZAF, type = "drift", criterion = "BIC", max.lag.y = max.lag)
# Significant result thus there is stationarity

## Emission per capita
ts_emission_ZAF <- ts(data_emission$CO2_ZAF, frequency = 1, start = c(1970))
ts_logemission_ZAF <- log(ts_emission_ZAF)
CADFtest(ts_logemission_ZAF, type = "trend", criterion = "BIC", max.lag.y = max.lag)
dts_logemission_ZAF <- diff(ts_logemission_ZAF)
CADFtest(dts_logemission_ZAF, type = "drift", criterion = "BIC", max.lag.y = max.lag)
# Significant result thus there is stationarity

#######################################
####### Check for Cointegration #######
#######################################

logdata<-data.frame(ts_logemission_ZAF, ts_logrenewable_ZAF)
names(logdata)<-c("logemission","logRenew")
attach(logdata)
VARselect(logdata,lag.max=15,type="const")
# Best model has K=14 acording to AIC
# See if Cointegration exists
trace_test<-ca.jo(logdata,type="trace",K=14,ecdet="const",spec="transitory")
summary(trace_test)
maxeigen_test<-ca.jo(logdata,type="eigen",K=14,ecdet="const",spec="transitory")
summary(maxeigen_test)
# both tests show cointegration!!!
# -> Use the VECM model
fit_vecm1<-cajorls(trace_test, r=1)
summary(fit_vecm1$rlm)

## Very interesting way stronger relations than Germany.
# High R squared, much smaller adjusted R-squared --> overfitting?
# F-statistic not significant

# Residuals:
#   Min        1Q    Median        3Q       Max 
# -0.054580 -0.010177  0.000456  0.014916  0.035747 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)  
# ect1              0.879788   0.331585   2.653   0.0211 *
#   logemission.dl1  -1.232577   0.486010  -2.536   0.0261 *
#   logRenew.dl1      0.216664   0.078321   2.766   0.0171 *
#   logemission.dl2  -0.763428   0.430968  -1.771   0.1019  
# logRenew.dl2      0.197336   0.076287   2.587   0.0238 *
#   logemission.dl3  -0.304947   0.388758  -0.784   0.4480  
# logRenew.dl3      0.187749   0.073472   2.555   0.0252 *
#   logemission.dl4  -0.341487   0.349322  -0.978   0.3476  
# logRenew.dl4      0.154637   0.072316   2.138   0.0537 .
# logemission.dl5  -0.523623   0.329157  -1.591   0.1376  
# logRenew.dl5      0.119841   0.065444   1.831   0.0920 .
# logemission.dl6  -0.160887   0.269838  -0.596   0.5621  
# logRenew.dl6      0.148972   0.058134   2.563   0.0249 *
#   logemission.dl7  -0.756254   0.281901  -2.683   0.0199 *
#   logRenew.dl7      0.150675   0.054130   2.784   0.0165 *
#   logemission.dl8  -0.756378   0.334402  -2.262   0.0431 *
#   logRenew.dl8      0.145638   0.051122   2.849   0.0147 *
#   logemission.dl9  -0.359725   0.288472  -1.247   0.2362  
# logRenew.dl9      0.093594   0.047185   1.984   0.0707 .
# logemission.dl10 -0.441649   0.268189  -1.647   0.1255  
# logRenew.dl10     0.058248   0.037400   1.557   0.1453  
# logemission.dl11 -0.043032   0.259329  -0.166   0.8710  
# logRenew.dl11     0.038623   0.031886   1.211   0.2491  
# logemission.dl12 -0.484715   0.228536  -2.121   0.0554 .
# logRenew.dl12     0.018276   0.020725   0.882   0.3952  
# logemission.dl13 -0.512372   0.227433  -2.253   0.0438 *
#   logRenew.dl13     0.005529   0.012985   0.426   0.6778  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.03372 on 12 degrees of freedom
# Multiple R-squared:  0.7938,	Adjusted R-squared:  0.3298 
# F-statistic: 1.711 on 27 and 12 DF,  p-value: 0.1647

###########################
####### Time Series EGY #######
###########################

ts_renewable_EGY <- ts(data_EGY$renewables_share_energy[6:58], frequency = 1, start = c(1970))
ts_logrenewable_EGY <- log(ts_renewable_EGY)
CADFtest(ts_logrenewable_EGY, type = "trend", criterion = "BIC", max.lag.y = max.lag)
dts_logrenewable_EGY <- diff(ts_logrenewable_EGY)
CADFtest(dts_logrenewable_EGY, type = "drift", criterion = "BIC", max.lag.y = max.lag)
# Significant result thus there is stationarity

## Emission per capita
ts_emission_EGY <- ts(data_emission$CO2_EGY, frequency = 1, start = c(1970))
ts_logemission_EGY <- log(ts_emission_EGY)
CADFtest(ts_logemission_EGY, type = "trend", criterion = "BIC", max.lag.y = max.lag)
dts_logemission_EGY <- diff(ts_logemission_EGY)
CADFtest(dts_logemission_EGY, type = "drift", criterion = "BIC", max.lag.y = max.lag)
# Significant result thus there is stationarity

#######################################
####### Check for Cointegration #######
#######################################

logdata<-data.frame(ts_logemission_EGY, ts_logrenewable_EGY)
names(logdata)<-c("logemission","logRenew")
attach(logdata)
VARselect(logdata,lag.max=15,type="const")
# AIC and HQ show 14 lags to be the best.
# See if Cointegration exists
trace_test<-ca.jo(logdata,type="trace",K=14,ecdet="const",spec="transitory")
summary(trace_test)
maxeigen_test<-ca.jo(logdata,type="eigen",K=14,ecdet="const",spec="transitory")
summary(maxeigen_test)
# both tests show cointegration but only at the 1% level and not 5% level!!!
# -> Use the VECM model
fit_vecm1<-cajorls(trace_test, r=1)
summary(fit_vecm1$rlm)


###########################
####### Time Series LKA #######
###########################
data_LKA["renewables_share"] <- (data_LKA$renewables_consumption/data_LKA$primary_energy_consumption)*100
ts_renewable_LKA <- ts(data_LKA$renewables_share[6:58], frequency = 1, start = c(1970))
ts_logrenewable_LKA <- log(ts_renewable_LKA)
ts.plot(ts_logrenewable_LKA)
CADFtest(ts_logrenewable_LKA, type = "trend", criterion = "BIC", max.lag.y = max.lag)
# Borederline no trend
dts_logrenewable_LKA <- diff(ts_logrenewable_LKA)
CADFtest(dts_logrenewable_LKA, type = "drift", criterion = "BIC", max.lag.y = max.lag)
# Significant result thus there is stationarity

## Emission per capita
ts_emission_LKA <- ts(data_emission$CO2_LKA, frequency = 1, start = c(1970))
ts_logemission_LKA <- log(ts_emission_LKA)
CADFtest(ts_logemission_LKA, type = "trend", criterion = "BIC", max.lag.y = max.lag)
dts_logemission_LKA <- diff(ts_logemission_LKA)
CADFtest(dts_logemission_LKA, type = "drift", criterion = "BIC", max.lag.y = max.lag)
# Significant result thus there is stationarity

#######################################
####### Check for Cointegration #######
#######################################

logdata<-data.frame(ts_logemission_LKA, ts_logrenewable_LKA)
names(logdata)<-c("logemission","logRenew")
attach(logdata)
VARselect(logdata,lag.max=15,type="const")
# AIC suggests lag 14
# See if Cointegration exists
trace_test<-ca.jo(logdata,type="trace",K=14,ecdet="const",spec="transitory")
summary(trace_test)
maxeigen_test<-ca.jo(logdata,type="eigen",K=14,ecdet="const",spec="transitory")
summary(maxeigen_test)
# both tests show cointegration but only at the 1% level and not 5% level!!!
# both tests show cointegration!!!
# VER Model
dlogdata<-data.frame(dts_logemission_LKA, dts_logrenewable_LKA)
names(dlogdata)<-c("dlogemission","dlogRenew")
attach(dlogdata)
VARselect(dlogdata,lag.max=15,type="const")
# AIC indicates 15
fit_varautom<-VAR(dlogdata,type="const",p=1)
summary(fit_varautom)
varautom_residuals<-resid(fit_varautom)

# win.graph()
par(mfrow=c(2,2))
acf(varautom_residuals[,1])
acf(varautom_residuals[,2])
ccf(varautom_residuals[,1],varautom_residuals[,2])
# Model does not seem valid

# Also try the VECM model
fit_vecm1<-cajorls(trace_test, r=1)
summary(fit_vecm1$rlm)
# No significant results
# Bad discrepancy between R-squared and Adjusted R-squared
# None significant F-test

# Residuals:
# Min       1Q   Median       3Q      Max 
# -0.12349 -0.05172 -0.01871  0.05130  0.17245 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)
# ect1             -0.02436    0.02316  -1.052    0.314
# logemission.dl1   0.11515    0.42704   0.270    0.792
# logRenew.dl1     -0.01720    0.30556  -0.056    0.956
# logemission.dl2   0.16567    0.49480   0.335    0.744
# logRenew.dl2     -0.07628    0.39483  -0.193    0.850
# logemission.dl3  -0.49303    0.57109  -0.863    0.405
# logRenew.dl3     -0.40271    0.44893  -0.897    0.387
# logemission.dl4   0.22817    0.56567   0.403    0.694
# logRenew.dl4     -0.26744    0.47418  -0.564    0.583
# logemission.dl5  -0.35483    0.61789  -0.574    0.576
# logRenew.dl5     -0.41054    0.50959  -0.806    0.436
# logemission.dl6  -0.61523    0.51667  -1.191    0.257
# logRenew.dl6     -0.36176    0.46909  -0.771    0.456
# logemission.dl7  -0.50020    0.52412  -0.954    0.359
# logRenew.dl7     -0.30104    0.43301  -0.695    0.500
# logemission.dl8  -0.40899    0.53914  -0.759    0.463
# logRenew.dl8     -0.24883    0.42136  -0.591    0.566
# logemission.dl9  -0.35823    0.56153  -0.638    0.535
# logRenew.dl9     -0.41127    0.37071  -1.109    0.289
# logemission.dl10 -0.58018    0.57250  -1.013    0.331
# logRenew.dl10    -0.29684    0.29126  -1.019    0.328
# logemission.dl11 -0.04765    0.47608  -0.100    0.922
# logRenew.dl11    -0.09194    0.24202  -0.380    0.711
# logemission.dl12 -0.17395    0.44032  -0.395    0.700
# logRenew.dl12    -0.09816    0.22559  -0.435    0.671
# logemission.dl13 -0.51640    0.43002  -1.201    0.253
# logRenew.dl13    -0.02692    0.22673  -0.119    0.907
# 
# Residual standard error: 0.1217 on 12 degrees of freedom
# Multiple R-squared:  0.6283,	Adjusted R-squared:  -0.2079 
# F-statistic: 0.7513 on 27 and 12 DF,  p-value: 0.7415

###########################
####### Time Series MAR #######
###########################
data_MAR["renewables_share"] <- (data_MAR$renewables_consumption/data_MAR$primary_energy_consumption)*100
ts_renewable_MAR <- ts(data_MAR$renewables_share[6:58], frequency = 1, start = c(1970))
ts_logrenewable_MAR <- log(ts_renewable_MAR)
ts.plot(ts_logrenewable_MAR)
CADFtest(ts_logrenewable_MAR, type = "trend", criterion = "BIC", max.lag.y = max.lag)
dts_logrenewable_MAR <- diff(ts_logrenewable_MAR)
CADFtest(dts_logrenewable_MAR, type = "drift", criterion = "BIC", max.lag.y = max.lag)
# Significant result thus there is stationarity

## Emission per capita
ts_emission_MAR <- ts(data_emission$CO2_MAR, frequency = 1, start = c(1970))
ts_logemission_MAR <- log(ts_emission_MAR)
CADFtest(ts_logemission_MAR, type = "trend", criterion = "BIC", max.lag.y = max.lag)
dts_logemission_MAR <- diff(ts_logemission_MAR)
CADFtest(dts_logemission_MAR, type = "drift", criterion = "BIC", max.lag.y = max.lag)
# Significant result thus there is stationarity

#######################################
####### Check for Cointegration #######
#######################################

logdata<-data.frame(ts_logemission_MAR, ts_logrenewable_MAR)
names(logdata)<-c("logemission","logRenew")
attach(logdata)
VARselect(logdata,lag.max=15,type="const")
# Lag 14 shows the best AIC
# See if Cointegration exists
trace_test<-ca.jo(logdata,type="trace",K=14,ecdet="const",spec="transitory")
summary(trace_test)
maxeigen_test<-ca.jo(logdata,type="eigen",K=14,ecdet="const",spec="transitory")
summary(maxeigen_test)
# both tests show cointegration!!!
# -> Use the VECM model
fit_vecm1<-cajorls(trace_test, r=1)
summary(fit_vecm1$rlm)
# Significant negative renew effects
# Also one negative significant effect of log emission
# Borderline none significant long-term effect
# High R-squared, lower adjusted R-squared
# Borderline significant F-statistic

# Residuals:
#   Min        1Q    Median        3Q       Max 
# -0.027071 -0.012160 -0.001213  0.011557  0.036870 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)  
# ect1              0.026825   0.013410   2.000   0.0686 .
# logemission.dl1  -0.364220   0.266927  -1.364   0.1975  
# logRenew.dl1      0.061930   0.024587   2.519   0.0270 *
#   logemission.dl2   0.034712   0.292111   0.119   0.9074  
# logRenew.dl2     -0.016670   0.022837  -0.730   0.4794  
# logemission.dl3  -0.380701   0.283919  -1.341   0.2048  
# logRenew.dl3      0.060331   0.020519   2.940   0.0124 *
#   logemission.dl4   0.450869   0.343564   1.312   0.2140  
# logRenew.dl4      0.012727   0.025638   0.496   0.6286  
# logemission.dl5  -0.901802   0.351881  -2.563   0.0249 *
#   logRenew.dl5     -0.005556   0.020591  -0.270   0.7919  
# logemission.dl6  -0.191546   0.280453  -0.683   0.5076  
# logRenew.dl6      0.031647   0.021829   1.450   0.1727  
# logemission.dl7   0.389236   0.290786   1.339   0.2055  
# logRenew.dl7     -0.054699   0.023450  -2.333   0.0379 *
#   logemission.dl8  -0.803547   0.308170  -2.607   0.0229 *
#   logRenew.dl8      0.018222   0.021377   0.852   0.4107  
# logemission.dl9   0.414091   0.320712   1.291   0.2210  
# logRenew.dl9     -0.013581   0.019203  -0.707   0.4929  
# logemission.dl10 -0.331751   0.272815  -1.216   0.2474  
# logRenew.dl10    -0.017731   0.023923  -0.741   0.4729  
# logemission.dl11  0.294087   0.224470   1.310   0.2147  
# logRenew.dl11     0.013479   0.022415   0.601   0.5588  
# logemission.dl12 -0.055091   0.239811  -0.230   0.8222  
# logRenew.dl12    -0.070173   0.023245  -3.019   0.0107 *
#   logemission.dl13 -0.392632   0.229157  -1.713   0.1123  
# logRenew.dl13     0.017756   0.027880   0.637   0.5362  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.03024 on 12 degrees of freedom
# Multiple R-squared:  0.8488,	Adjusted R-squared:  0.5086 
# F-statistic: 2.495 on 27 and 12 DF,  p-value: 0.04922
# 
