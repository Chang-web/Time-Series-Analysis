
# Time Series Project

# 研究目標: 以時間序列分析美國政府升息政策是否造成核心通膨率下降?

# 研究期間: 1990年1月 ~ 2017年12月

# 資料來源: FRED

setwd("D:/NCCU/double major and minor/Statistics/time series analysis/report/data")
library(readr)
corecpi <- read_csv("CPILFENS_1990-2017.csv", show_col_types = F)
colnames(corecpi) <- c("date", "cpi")
# 原始資料時間序列圖
cpi <- ts(corecpi$cpi, start=c(1990,1),end=c(2017,12),frequency=12)
plot(cpi, xlab="year", ylab="cpi index", main = "Monthly CPI index", lwd = 2) 

# 除了讓時間數列更為平穩，因為是要觀察核心通膨率，所以對 CPI 取 log，再開始處理並進行分析。
corelogcpi <- data.frame("date" = corecpi$date, "logcpi" = log(corecpi$cpi))
logcpi <- ts(corelogcpi$logcpi, start=c(1990,1),end=c(2017,12),frequency=12)
# 取對數後的時間序列圖，有時間趨勢
plot(logcpi, xlab="year", ylab="logcpi index", main = "Monthly CPI index", lwd = 2);abline(reg = lm(logcpi ~ time(logcpi)), col = "blue", lwd = 2) # time trend

# 切割訓練和測試集: 切割標準為聯準會在2006年後的首度升息2015年12月，因為接近年底，以2016年1月為分割點。
## 訓練集1990年1月 ~ 2015年12月(共312筆資料)
## 測試集2016年1月 ~ 2017年12月(共24筆資料)
#### 註記1: 因為2019年後的疫情，經濟數據變動可能受到其他外生因素干擾，會影響預測結果，因此取較早的升息時間點做本次報告的分析。
#### 註記2: 測試集原先只需要取5 ~ 10筆，但因為政策對經濟的影響有時間延滯，所以增加測試集區間為2年(24筆)。

corecpi_train <- corelogcpi[-c(313:336), ]
corecpi_test <- corelogcpi[313:336, ]
str(corecpi_train)

logcpi <- ts(corecpi_train$logcpi, start=c(1990,1),end=c(2015,12), frequency = 12)
library(astsa)
logcpiacf2 <- acf2(logcpi, max.lag = 100) # acf: slow decay


# 分析流程:

## 方法一、差分處理:
### a. 有時間趨勢，取一階差分
d1logcpi <- diff(logcpi, 1) 
ts.plot(d1logcpi, main = "first difference of logcpi")
### 此序列雖然沒有考慮季節性差分，因為已經是定態，所以我們還是對其配適模型
library(tseries)
adf.test(d1logcpi) # stationary
### b. 一階差分後，季節性部分acf緩慢遞減，再取一階lag = 12之季節性差分
d1logcpiacf <- acf2(d1logcpi,  max.lag = 180)
d12d1logcpi <- diff(d1logcpi, 12) 
ts.plot(d12d1logcpi, main = "first and seasonal difference of logcpi")
d12d1logcpiacf <- acf2(d12d1logcpi,  max.lag = 150) # 沒有明顯趨勢，不須再做差分
adf.test(d12d1logcpi) # stationary


## 方法二、去除趨勢前處理: 先去除CPI增長的趨勢，使用時間t作為迴歸變數，得到的fitted values計算殘差，觀察殘差時間序列圖，再取適當的差分。
t <- 1:312
reg1 <- lm(logcpi ~ 1 + t, data = corecpi_train)
res1 <- reg1$residuals
ts.plot(res1, main = "residual time series of regression on time")
res1acf <- acf2(res1, max.lag = 100)
### res1的acf有緩慢遞減，對其取一階差分
d1res1 <- diff(res1, 1) 
ts.plot(d1res1, main = c(paste("residual time series of regression on time"),
                         paste("(first order difference)")))
d1res1acf <- acf2(d1res1, max.lag = 180)
### 一階差分後，季節性部分緩慢遞減，再取一階lag = 12之季節性差分
d12d1res1 <- diff(d1res1, 12) 
ts.plot(d12d1res1, main = c(paste("residual time series of regression on time"),
                         paste("(first and seasonal difference)")))
d12d1res1acf <- acf2(d12d1res1, max.lag = 150) # 沒有明顯趨勢，不須再做差分
adf.test(d12d1res1) # stationary
### 我們發現 d12d1res1 和 d12d1logcpi 的 acf, pacf 相同，但因為是從不同方式得到的序列，且後續計算PMSE的方式也相異，因此雖然配適模型order選取相同，仍視為兩種情形。


## 方法三、迴歸前處理: 使用 IPI 和 UE 作為迴歸變數，觀察得到殘差的時間序列圖，再取適當的差分。
### 次要變數一: 工業生產指數 industry production index
ipi <- read_csv("IPI_1990-2017.csv", show_col_types = F)
## train and test set
ipi_train <- ipi[-c(313:336), ]
### 次要變數二: 失業率 unemployment rate
ue <-  read_csv("UE_1990-2017.csv", show_col_types = F)
### train and test set
ue_train <- ue[-c(313:336), ]
traindata <- data.frame(corecpi_train, IPI = ipi_train$IPB50001N, UE = ue_train$UNRATENSA) 
str(traindata)  # 全部都是數值型變數
### regression
### a. 此迴歸殘差序列非定態，但因為ACF沒有緩慢遞減，我們嘗試對其建模。
reg2 <- lm(logcpi ~ IPI + UE, data = traindata)
res2 <- reg2$residuals
ts.plot(res2, main = "residual time series of multiple linear regression")
res2acf <- acf2(res2, max.lag = 100) 
adf.test(res2) # not stationary

### 迴歸殘差序列非定態，且時間序列圖還有些上升的趨勢，考慮一階差分
d1res2 <- diff(res2, 1) 
ts.plot(d1res2, main = c(paste("residual time series of multiple linear regression"),
                         paste("(first order difference)")))
d1res2acf <- acf2(d1res2, max.lag = 180) 
adf.test(d1res2)
### 一階差分為定態序列，但是acf在seasonal有緩慢遞減的現象，考慮再一次季節性差分
d12d1res2 <- diff(d1res2, 12) 
ts.plot(d12d1res2, main = c(paste("residual time series of multiple linear regression"),
                         paste("(first and seasonal difference)")))
### b. 此迴歸殘差序列沒有明顯的趨勢，不須再差分
d12d1res2acf <- acf2(d12d1res2, max.lag = 150) 
adf.test(d12d1res2) # stationary


# 分別對四種序列(方法一:a.b.; 方法三:a.b.)建立 SARIMA 模型: 畫出acf和pacf，決定order

## 方法一、二之時間數列圖都在 0 附近震盪，所以配適模型皆不考慮估計常數項
## 方法一:
### a.
d1logcpiacf2 <- acf2(d1logcpi, max.lag = 150)
### seasonal part: acf and pacf tail off => seasonal ARMA(1,1)
### non-seasonal part: acf tails off; pacf 在前兩期 seasonal lag 左側有一期顯著 => AR(1)
fit00 <- stats::arima(logcpi, order=c(1,1,0), seasonal=list(order=c(1,0,1), period=12), method = "CSS-ML")
fit00 # 係數都有顯著
### seasonal part: acf and pacf tail off => seasonal ARMA(1,1)
### non-seasonal part: acf 在seasonal 左右各有一期顯著 => MA(1)
fit01 <- stats::arima(logcpi, order=c(0,1,1), seasonal=list(order=c(1,0,1), period=12), method = "CSS-ML")
fit01 # 係數都有顯著
### seasonal part: acf and pacf tail off => seasonal ARMA(1,1)
### non-seasonal part: acf and pacf tail off => ARMA(1,1)
fit02 <- stats::arima(logcpi, order=c(1,1,1), seasonal=list(order=c(1,0,1), period=12), method = "CSS-ML")
fit02 # 係數都有顯著
### seasonal part: acf tails off; pacf lag = 2 => seasonal AR(2)
### non-seasonal part: acf tails off;  pacf 在前兩期 seasonal lag 左側有一期顯著 => AR(1)
fit10 <- stats::arima(logcpi, order=c(1,1,0), seasonal=list(order=c(2,0,0), period=12), method = "CSS-ML")
fit10 # 係數都有顯著
### seasonal part: acf tails off; pacf lag = 2 => seasonal AR(2)
### non-seasonal part: acf 在seasonal 左右各有一期顯著 => MA(1)
fit11 <- stats::arima(logcpi, order=c(0,1,1), seasonal=list(order=c(2,0,0), period=12), method = "CSS-ML")
fit11 # 係數都有顯著
### seasonal part: acf tails off; pacf lag = 2 => seasonal AR(2)
### non-seasonal part:  acf and pacf tail off => ARMA(1,1)
fit12 <- stats::arima(logcpi, order=c(1,1,1), seasonal=list(order=c(2,0,0), period=12), method = "CSS-ML")
fit12 # 係數都有顯著

### b.
d12d1logcpiacf2 <- acf2(d12d1logcpi, max.lag = 150)
### seasonal part: acf cuts off at lag = 1; pacf tails off => seasonal MA(1)
### non-seasonal part: acf 在seasonal lag = 1 左側兩期顯著; pacf tails off => MA(2)
fit20 <- stats::arima(logcpi, order=c(0,1,2), seasonal=list(order=c(0,1,1), period=12), method = "CSS-ML")
fit20 # ma2 係數不顯著，需要修正
fit21 <- stats::arima(logcpi, order=c(0,1,1), seasonal=list(order=c(0,1,1), period=12), method = "CSS-ML")
fit21 # 係數都有顯著
### seasonal part: acf cuts off at lag = 1; pacf tails off => seasonal MA(1)
### non-seasonal part: acf and pacf tail off => ARMA(1,1)
fit22 <- stats::arima(logcpi, order=c(1,1,1), seasonal=list(order=c(0,1,1), period=12), method = "CSS-ML")
fit22 # 係數都有顯著
### seasonal part: acf and pacf tail off => seasonal ARMA(1,1)
### non-seasonal part: acf 在seasonal lag = 1 左側兩期顯著; pacf tails off => MA(2)
fit30 <- stats::arima(logcpi, order=c(0,1,2), seasonal=list(order=c(1,1,1), period=12), method = "CSS-ML")
fit30 # ma2, sar1 係數不顯著，需要修正，修正後結果為 fit21
### seasonal part: acf and pacf tail off => seasonal ARMA(1,1)
### non-seasonal part: acf and pacf tail off => ARMA(1,1)
fit31 <- stats::arima(logcpi, order=c(1,1,1), seasonal=list(order=c(1,1,1), period=12), method = "CSS-ML")
fit31 # sar1 係數不顯著，需要修正, 修正後為 fit22


## 方法二: (選取order的理由與方法一a.相同，差別在於使用資料為迴歸後的殘差)
d12d1res1acf2 <- acf2(d12d1res1, max.lag = 150)
### seasonal part: acf cuts off at lag = 1; pacf tails off => seasonal MA(1)
### non-seasonal part: acf 在seasonal lag = 1 左側兩期顯著; pacf tails off => MA(2)
fit40 <- stats::arima(res1, order=c(0,1,2), seasonal=list(order=c(0,1,1), period=12), method = "CSS-ML")
fit40 # ma2 係數不顯著，需要修正
fit41 <- stats::arima(res1, order=c(0,1,1), seasonal=list(order=c(0,1,1), period=12), method = "CSS-ML")
fit41 # 係數都有顯著
### seasonal part: acf cuts off at lag = 1; pacf tails off => seasonal MA(1)
### non-seasonal part: acf and pacf tail off => ARMA(1,1)
fit42 <- stats::arima(res1, order=c(1,1,1), seasonal=list(order=c(0,1,1), period=12), method = "CSS-ML")
fit42 # 係數都有顯著
### seasonal part: acf and pacf tail off => seasonal ARMA(1,1)
### non-seasonal part: acf 在seasonal lag = 1 左側兩期顯著; pacf tails off => MA(2)
fit50 <- stats::arima(res1, order=c(0,1,2), seasonal=list(order=c(1,1,1), period=12), method = "CSS-ML")
fit50 # ma2, sar1 係數不顯著，需要修正，修正後結果為 fit41
### seasonal part: acf and pacf tail off => seasonal ARMA(1,1)
### non-seasonal part: acf and pacf tail off => ARMA(1,1)
fit51 <- stats::arima(res1, order=c(1,1,1), seasonal=list(order=c(1,1,1), period=12), method = "CSS-ML")
fit51 # sar1 係數不顯著，需要修正, 修正後為 fit42



## 方法三:
## 方法三a.之時間數列圖不在 0 附近震盪，所以配適模型考慮估計常數項
### a.
res2acf2 <- acf2(res2, max.lag = 150)
### seasonal part: acf tails off; pacf cut off at lag = 1 => seasonal AR(1)
### nonseasonal part: acf and pacf tail off => ARMA(1,1)
### 因為沒有在 0 附近震盪，所以估計常數項
fit60 <- stats::arima(res2, order=c(1,0,1), seasonal=list(order=c(0,0,1), period=12), method = "CSS-ML")
fit60 # 係數都有顯著

## 方法三b.之時間數列圖在 0 附近震盪，所以配適模型不考慮估計常數項
### b.
d12d1res2acf <- acf2(d12d1res2, max.lag = 150) 
### seasonal part: acf cut off at lag = 1; pacf tails off => seasonal MA(1)
### nonseasonal part:  acf 在seasonal 右或左側有一個顯著的lag => MA(1)
fit70 <- stats::arima(res2, order=c(0,1,1), seasonal=list(order=c(0,1,1), period=12), method = "CSS-ML")
fit70 # 係數都有顯著
### seasonal part: acf cut off at lag = 1; pacf tails off => seasonal MA(1)
### nonseasonal part:  acf tails off; pacf 在seasonal lag = 0,1 右側顯著 => AR(1)
fit71 <- stats::arima(res2, order=c(1,1,0), seasonal=list(order=c(0,1,1), period=12), method = "CSS-ML")
fit71 # 係數都有顯著
### seasonal part: pacf cut off at lag = 2 => seasonal AR(2)
### nonseasonal part:  acf 在seasonal 右或左側有一個顯著的lag => MA(1)
fit80 <- stats::arima(res2, order=c(0,1,1), seasonal=list(order=c(2,1,0), period=12), method = "CSS-ML")
fit80 # 係數都有顯著
### seasonal part: pacf cut off at lag = 2 => seasonal AR(2)
### nonseasonal part:  acf tails off; pacf 在seasonal lag = 0,1 右側顯著 => AR(1)
fit81 <- stats::arima(res2, order=c(1,1,0), seasonal=list(order=c(2,1,0), period=12), method = "CSS-ML")
fit81 # 係數都有顯著

# 分別列出四種序列下的候選模型，以AIC為標準，在三種方法下各自選出最佳模型: fit12; fit22; fit51
## 方法一: fit12
### a.
fit00$aic # -3373.177
fit01$aic # -3369
fit02$aic # -3354.22
fit10$aic # -3351.049
fit11$aic # -3348.36
fit12$aic # -3373.339
### b.
fit21$aic # -3259.66
fit22$aic # -3285.293


## 方法二: fit42
fit41$aic # -3260.736
fit42$aic # -3285.47

## 方法三: fit71
### a.
fit60$aic # -1655.65
### b.
fit70$aic # -1951.541
fit71$aic # -1952.119
fit80$aic # -1950.68
fit81$aic # -1951.195


# 確認三種模型都有通過殘差檢定。
## 註記1: 殘差檢定要用老師上課的 Ljung--Box test statistic 得到 p-value，用迴圈畫出圖形觀察是否累積lag結果都有超過0.05。
## 註記2: 可以再用 SARIMA 的輸出結果 QQplot, residuals acf 確認模型通過殘差檢定。

## 方法一: fit12 沒有通過殘差檢定
library(TSA)
stdresfit12 <- rstandard(fit12)  
ts.plot(stdresfit12, ylab="", main="standardized residuals(fit12)")
stdresfit12acf2 <- acf2(stdresfit12, max.lag = 150) 
B_text_p_value = c(0,0)
for(hh in 1:30){
  B_text_p_value[hh] = Box.test(stdresfit12, lag=hh, type="Ljung-Box")$p.value
}
plot(1:30, B_text_p_value[1:30], type="p", 
     main="p values for Ljung-Box statistic (fit12)", 
     xlab="lag", ylab="p value", ylim=c(0,1));abline(h=0.05, lty=2, col=4) 

### 修正模型: 
#### 觀察 fit12 標準化殘差的acf & pacf
#### seasonal part: acf tail off and pacf cut off at lag = 2. => add seasonal AR(2)
mfit120 <- stats::arima(logcpi, order=c(1,1,1), seasonal=list(order=c(4,0,0), period=12), method = "CSS-ML")
mfit120 # sar3 係數不顯著，需要修正
# 修正模型為 mfit121
mfit121 <- stats::arima(logcpi, order=c(1,1,1), seasonal=list(order=c(3,0,0), period=12), method = "CSS-ML")
mfit121 # 係數都有顯著
mfit121$aic # -3378.061 更小的 AIC

### mfit121 有通過殘差檢定
stdresmfit121 <- rstandard(mfit121)  
ts.plot(stdresmfit121, ylab="", main="standardized residuals(mfit121)")
stdresmfit121acf2 <- acf2(stdresmfit121, max.lag = 150) # close to white noise
B_text_p_value = c(0,0)
for(hh in 1:30){
  B_text_p_value[hh] = Box.test(stdresmfit121, lag=hh, type="Ljung-Box")$p.value
}
plot(1:30, B_text_p_value[1:30], type="p", 
     main="p values for Ljung-Box statistic (mfit121)", 
     xlab="lag", ylab="p value", ylim=c(0,1));abline(h=0.05, lty=2, col=4) 

## check
sarima(logcpi,1,1,1,3,0,0,12)

## 方法二: fit42 有通過殘差檢定
stdresfit42 <- rstandard(fit42)  
ts.plot(stdresfit42, ylab="", main="standardized residuals(fit42)")
stdresfit42acf2 <- acf2(stdresfit42, max.lag = 100) # close to white noise
B_text_p_value = c(0,0)
for(hh in 1:100){
  B_text_p_value[hh] = Box.test(stdresfit42, lag=hh, type="Ljung-Box")$p.value
}
plot(1:100, B_text_p_value[1:100], type="p", 
     main="p values for Ljung-Box statistic (fit42)", 
     xlab="lag", ylab="p value", ylim=c(0,1));abline(h=0.05, lty=2, col=4)

## check
sarima(res1, 1,1,1,0,1,1,12) 

## 方法三: fit71 沒有通過殘差檢定
stdresfit71 <- rstandard(fit71)  
ts.plot(stdresfit71, ylab="", main="standardized residuals(fit71)")
stdresfit71acf2 <- acf2(stdresfit71, max.lag = 100) # close to white noise
B_text_p_value = c(0,0)
for(hh in 1:100){
  B_text_p_value[hh] = Box.test(stdresfit71, lag=hh, type="Ljung-Box")$p.value
}
plot(1:100, B_text_p_value[1:100], type="p", 
     main="p values for Ljung-Box statistic (fit71)", 
     xlab="lag", ylab="p value", ylim=c(0,1));abline(h=0.05, lty=2, col=4) 

### 修正模型: 
#### 觀察 fit71 標準化殘差的acf & pacf
#### non-seasonal part: acf and pacf tail off => ARMA(1,1)
mfit71 <- stats::arima(res2, order=c(2,1,1), seasonal=list(order=c(0,1,1), period=12), method = "CSS-ML")
mfit71 # 係數都有顯著
mfit71$aic # -1951.004 AIC 有稍微大一些

### mfit71 有通過殘差檢定
stdresmfit71 <- rstandard(mfit71)  
ts.plot(stdresmfit71, ylab="", main="standardized residuals(mfit71)")
stdresmfit71acf2 <- acf2(stdresmfit71, max.lag = 100) # close to white noise
B_text_p_value = c(0,0)
for(hh in 1:100){
  B_text_p_value[hh] = Box.test(stdresmfit71, lag=hh, type="Ljung-Box")$p.value
}
plot(1:100, B_text_p_value[1:100], type="p", 
     main="p values for Ljung-Box statistic (mfit71)", 
     xlab="lag", ylab="p value", ylim=c(0,1));abline(h=0.05, lty=2, col=4) 

## check
sarima(res2, 2,1,1,0,1,1,12) 

## 最後在三種方式有通過殘差檢定的模型分別為: mfit121, fit42, mfit71


# 預測
library(forecast)
## 24-steps ahead predictions and confidence intervals
forcastm121 <- forecast(mfit121, level=c(95), h=2*12)
forcast42 <- forecast(fit42, level=c(95), h=2*12)
forcastm71 <- forecast(mfit71, level=c(95), h=2*12)
forcastm121
forcast42
forcastm71

## 視覺化呈現預測結果
## Fan charts
plot(forcastm121, main = "Froecasts from SARIMA(1,1,1)*(3,0,0)[12]")

tscorelogcpi <- ts(corelogcpi$logcpi, start=c(1990,1),end=c(2017,12),frequency=12)
plot(tscorelogcpi, xlab="year", ylab="logcpi index", 
     main = "Froecasts from SARIMA(1,1,1)*(3,0,0)[12]", 
     lwd = 2, lty = 1, ylim = c(4.85, 5.7))
lines(forcastm121$mean, type = "l", col = "green", lwd = 2)
lines(forcastm121$lower, col = "blue", lwd = 2, lty = 1)
lines(forcastm121$upper, col = "red", lwd = 2, lty = 1)
legend("topleft", 
       c("upperbound", "prediction", "lowerbound", "core_cpi% in dataset"),
       lty = 1, lwd = 2, 
       col = c("red", "green", "blue", "black"))


### fan chart of residuals
plot(forcast42, main = "Froecast of residuals from SARIMA(1,1,1)*(0,1,1)[12]")
### 從迴歸我們有係數的估計值，將這些係數和所對應測試集的解釋變數相乘並加總，再將這些數值加上預測殘差，結果為測試集的配適值(fitted value)。
beta0_hat <- reg1$coefficients[[1]]
beta1_hat <- reg1$coefficients[[2]]
prederror1 <- forcast42$mean
t <- 313:336
fittedreg1 <- beta0_hat + t*beta1_hat + prederror1
lbd_fit42 <- beta0_hat + t*beta1_hat + forcast42$lower
ubd_fit42 <- beta0_hat + t*beta1_hat + forcast42$upper 
fittedreg1 <- ts(fittedreg1, start = c(2016,1), end = c(2017,12), frequency = 12)
lbd_fit42 <- ts(lbd_fit42, start = c(2016,1), end = c(2017,12), frequency = 12)
ubd_fit42 <- ts(ubd_fit42, start = c(2016,1), end = c(2017,12), frequency = 12)
plot(tscorelogcpi, xlab="year", ylab="logcpi index", 
     main = "Froecasts from SARIMA(1,1,1)*(0,1,1)[12]", 
     lwd = 2, lty = 1, ylim = c(4.85, 5.7))
lines(fittedreg1, type = "l", col = "green", lwd = 2)
lines(lbd_fit42, col = "blue", lwd = 2, lty = 1)
lines(ubd_fit42, col = "red", lwd = 2, lty = 1)
legend("topleft", 
       c("upperbound", "prediction", "lowerbound", "core_cpi%  in dataset"),
       lty = 1, lwd = 2, 
       col = c("red", "green", "blue", "black"))


### fan chart of residuals
plot(forcastm71, main = "Froecast of residuals from SARIMA(2,1,1)*(0,1,1)[12]")

beta0_hat <- reg2$coefficients[[1]]
beta1_hat <- reg2$coefficients[[2]]
beta2_hat <- reg2$coefficients[[3]]
ipi_test <- ipi[313:336, 2]
ue_test <- ue[313:336, 2]
prederror2 <- forcastm71$mean
fittedreg2 <- c() 
lbd_fitm71 <- c()
ubd_fitm71 <- c()
for(i in 1:24){
  
  fittedreg2[i] <- beta0_hat + beta1_hat * ipi_test$IPB50001N[i] + beta2_hat * ue_test$UNRATENSA[i] + prederror2[i]
  lbd_fitm71[i] <- beta0_hat + beta1_hat * ipi_test$IPB50001N[i] + beta2_hat * ue_test$UNRATENSA[i] + forcastm71$lower[i]
  ubd_fitm71[i] <- beta0_hat + beta1_hat * ipi_test$IPB50001N[i] + beta2_hat * ue_test$UNRATENSA[i] + forcastm71$upper[i]
  
}
fittedreg2 <- ts(fittedreg2, start = c(2016,1), end = c(2017,12), frequency = 12)
lbd_fitm71 <- ts(lbd_fitm71, start = c(2016,1), end = c(2017,12), frequency = 12)
ubd_fitm71 <- ts(ubd_fitm71, start = c(2016,1), end = c(2017,12), frequency = 12)
plot(tscorelogcpi, xlab="year", ylab="logcpi index", 
     main = "Froecasts from SARIMA(2,1,1)*(0,1,1)[12]",
     lwd = 2, lty = 1, ylim = c(4.85, 5.7))
lines(fittedreg2, type = "l", col = "green", lwd = 2)
lines(lbd_fitm71, col = "blue", lwd = 2, lty = 1)
lines(ubd_fitm71, col = "red", lwd = 2, lty = 1)
legend("topleft", 
       c("upperbound", "prediction", "lowerbound", "core_cpi% in dataset"),
       lty = 1, lwd = 2, 
       col = c("red", "green", "blue", "black"))


## 觀察測試集的部分
#### 由下圖， fit42(SARIMA(1,1,1)*(0,1,1)[12]) 折線與實際值折線最為接近，mfit70(SARIMA(2,1,1)*(0,1,1)[12]) 的預測情形最不理想。
#### 另外，可以觀察到在前15期，fit42的預測較為準確，然而在20期後，反而是 mfit121(SARIMA(1,1,1)*(3,0,0)[12]) 較接近實際值。
#### 從圖可看出 fit42, fitm71 的預測結果，都在年底(第12期和24期)有下降的現象，此結果和測試集的實際值一致。
plot(forcastm121$mean,
     ylim = c(5.495, 5.56), xlim = c(2016, 2018), xaxt="n",
     type = "l", col = "darkorange", lwd = 2, xlab = "",
     main = "Forcasts from models and observed data in 2016-2017", ylab = "core_cpi%")
axis(1, at=seq(2016, 2018, by=0.5), labels = c("2016 Jan", "2016 June", "2017 Jan", "2017 June", "2018 Jan"))
lines(fittedreg1, type = "l", col = "green", lwd = 2)
lines(fittedreg2, type = "l", col = "blue", lwd = 2)
observed <- ts(corecpi_test$logcpi, start = c(2016,1), end = c(2017,12), frequency = 12)
lines(observed,type = "l", col = "black", lwd = 2)
points(forcastm121$mean, col = "darkorange", pch = 16)
points(fittedreg1, col = "green", pch = 16)
points(fittedreg2, col = "blue", pch = 16)
points(observed, col = "black", pch = 16)
legend("topleft", 
       c("SARIMA(1,1,1)*(3,0,0)[12]", "SARIMA(1,1,1)*(0,1,1)[12]", "SARIMA(2,1,1)*(0,1,1)[12]", "core_cpi% in test set"),
       lty = 1, lwd = 2,pch = 16, 
       col = c("darkorange", "green", "blue", "black"))




## 透過比較PMSE大小選出最佳模型
### 比較測試集實際和預測的差距
#### 方法一
pred_corecpi <- data.frame(forcastm121$mean)
error_mfit121 <- c()
errorsq_mfit121 <- c()
for(i in 1:24){
  
  error_mfit121[i] <- corecpi_test[i,2] - pred_corecpi[i,1]
  errorsq_mfit121[i] <- error_mfit121[i]^2
  PMSE_mfit121 <- mean(errorsq_mfit121)
  
}


### 我們將測試集的實際資料與其配適值相減，計算平方和再取平均，可以得到 PMSE。方法二、三的PMSE皆是透過此流程計算得出的結果。
#### 方法二
error_fit42 <- c()
errorsq_fit42 <- c()
for(i in 1:24){
  
  error_fit42[i] <- corecpi_test[i,2] - fittedreg1[i]
  errorsq_fit42[i] <- error_fit42[i]^2
  PMSE_fit42 <- mean(errorsq_fit42)
}

#### 方法三
error_mfit71 <- c()
errorsq_mfit71 <- c()
for(i in 1:24){
  
  error_mfit71[i] <- corecpi_test[i,2] - fittedreg2[i]
  errorsq_mfit71[i] <- error_mfit71[i]^2
  PMSE_mfit71 <- mean(errorsq_mfit71)
  
}

## 以表格呈現
result <- round(data.frame(errorsq_mfit121, errorsq_fit42, errorsq_mfit71), 6)
k <- c()
for(i in 1:24){
  k[i] <- paste(i,"-step")
}
rownames(result) <- k
colnames(result) <- c("mfit121", "fit42", "mfit71")

## fit42 有最小的 PMSE，其次是 mfit121，最後則是 mfit71。
## 因此 fit42(以時間t為解釋變數的迴歸之配飾模型) 有較良好的預測品質，以此模型作為本次報告通過殘差檢定的最佳配適模型。
## PMSE
PMSE <- round(c(PMSE_mfit121, PMSE_fit42, PMSE_mfit71), 6)
result <- data.frame(rbind(result, PMSE))
rownames(result)[25] <- "PMSE"
result


# 結論
## 討論 fit42(SARIMA(1,1,1)*(0,1,1)[12]) 的預測情形：
### 1. 
#### 將實際值(observed_data)、配適值(fitted values)、方差(square_error)、預測區間(prediction interval) 以表格呈現。
#### prediction intervals
#### 實際值都有在prediction interval內。
data.frame("observed_data" = corecpi_test$logcpi, 
           "fitted_value" = fittedreg1,
           "square_error" = round(errorsq_fit42,6),
           "95%_lower_bound" = lbd_fit42,
           "95%_upper_bound" = ubd_fit42
           )

### 2. 
#### 透過計算前後期的差距，可以發現在2016年與2017年底，都有負向的變化量，主要原因是美國聯準會在2015(0.25% –> 0.50%)和2016(0.50% –> 0.75%)年底都有宣布升息。
#### 因此選定的最佳配適模型(fit42)有確實反映出升息政策實施後，核心通膨率在未來下降的現象，表示升息政策有抑制通膨的效果，符合本次報告的預期結果。
diff <- c()
for(i in 2:24){
  diff[i] <- round(fittedreg1[i] - fittedreg1[i-1], 6)
}
diff

## 附錄
# 資料來源:
# https://fred.stlouisfed.org/series/CPILFENS
# https://fred.stlouisfed.org/series/UNRATENSA
# https://fred.stlouisfed.org/series/IPB50001N

# 參考資料:
# http://www.math.chalmers.se/Stat/Grundutb/GU/MSA220/S16/Lecture9-2015.pdf
# https://www.investopedia.com/ask/answers/12/inflation-interest-rate-relationship.asp
# https://money.cnn.com/2015/12/16/news/economy/federal-reserve-interest-rate-hike/index.html
# https://en.wikipedia.org/wiki/History_of_Federal_Open_Market_Committee_actions#December_2015_historic_interest_rate_hike

