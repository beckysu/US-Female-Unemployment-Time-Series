# US Female Unemployment
## Abstract
The purpose of this project is to accurately forecast the unemployment figures of females aged 16-19 based on previous rates. As the unemployment rate is dropping, unemployment figures are increasing exponentially, and we will examine this in my analysis. I will use RStudio for this time series analysis.

## Contributor
Becky Su

## Reproduction
Monthly U.S. female (16-19 years) unemployment figures (thousands) 1948-1981, data available on Prof. Rob Hyndmans Time Series Data Library (TSDL), accessed in RStudio.

## Statistical Analysis
In my analysis of this time series, I used various techniques, including Box-Cox transformations and differencing to remove trend and seasonality to achieve stationarity, thus allowing me to identify potential models using ACF and PACF plots. I then performed diagnostic checking on potential models to identify the best model to use for forecasting. During my analysis, I was able to come up with various candidate models, however, only one was most suitable for forecasting: SARIMA(0,1,3)x(0,1,1)[12] model. The SARIMA(0,1,3)x(0,1,1)[12] model passed the Box-Ljung and Box-Pierce tests, possesses the lowest AICc, and is viable for forecasting compared to our other candidates. It does not pass the Shapiro-Wilk test for normality of residuals but that can be attributed to its heavy tailed distribution. Through forecasting I was able to plot a potential trajectory with 95% confidence for 6 months in the future. Despite certain validation points being outside the confidence interval, the forecasted values still remain valid.

## Conclusion
My goal for this project was to predict the unemployment figures of females age 16-19 for the next year after 1981. I forecasted the unemployment rate using the model $SARIMA(0,1,2)(0,1,1)_{12}$. This is a viable model because none of my prediction intervals contain 0. Though not all of my validation points lie within the 95% CI, the majority of them are and, furthermore,the validation points converge towards our predictions as time goes on. Looking at the original time series, we can perhaps explain the 406th data point not lying within the 95% confidence interval with the seasonality, since the trend oscillates between the months of every year. However, with time, we can see that our predictions are indeed valid and accurate, as our prediction for December 1981 and our validation point are nearly the same. In conclusion, we can expect the unemployment figures to continue to fluctuate drastically between the months of every year but the rate of fluctuation will remain stagnant as it has since 1978. I've accomplished my goal of understanding the future unemployment rates of US female teens of ages 16-19, arriving at my final model of $SARIMA(0,1,2)(0,1,1)_{12}$.
Model: $X_t-X_{t-1}-X_{t-12}+X_{t-13}=Z_t-0.48Z_{t-1}-0.09Z_{t-2}-0.69Z_{t-12}+0.33Z_{t-13}+0.06Z_{t-14}$