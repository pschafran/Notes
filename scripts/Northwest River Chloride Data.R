setwd("~/Desktop")
NWRchlorideData <- read.csv("NorthwestRiverChlorideData.csv", header = TRUE, sep = ",")
View(NWRchlorideData)
NWRchlorideData[[1]] <- as.Date(NWRchlorideData[[1]])

chlorideLM <- lm(NWRchlorideData$Chloride_PPM ~ NWRchlorideData$Date)
colorLM <- lm(NWRchlorideData$Color ~ NWRchlorideData$Date)
turbidityLM <- lm(NWRchlorideData$Turbidity_NTU ~ NWRchlorideData$Date)

summary(chlorideLM)
summary(colorLM)
summary(turbidityLM)
par(mfrow=c(2,2))
plot(NWRchlorideData$Chloride_PPM ~ NWRchlorideData$Date, pch = 20, cex = 0.5, main = "Northwest River Chloride", ylab = "Chloride (PPM)", xlab = "")
abline(lm(NWRchlorideData$Chloride_PPM ~ NWRchlorideData$Date), col = red, lwd = 5)
text(16000, 1300, "Multiple R-squared: 0.08707\nAdjusted R-squared: 0.08687\nF-statistic: 428.5 on 1 and 4493 DF\np-value: < 2.2e-16", cex = 0.5)

plot(NWRchlorideData$Color ~ NWRchlorideData$Date, pch = 20, cex = 0.5, main = "Northwest River Color", xlab = "", ylab = "Color")
abline(lm(NWRchlorideData$Color ~ NWRchlorideData$Date), col = red, lwd = 5)
text(15000, 950, "Multiple R-squared: 0.02354\nAdjusted R-squared: 0.02332\nF-statistic: 108.3 on 1 and 4493 DF\np-value: < 2.2e-16", cex = 0.5)

plot(NWRchlorideData$Turbidity_NTU ~ NWRchlorideData$Date, pch = 20, cex = 0.5, main = "Northwest River Turbidity", xlab = "", ylab = "Turbidity (NTU)")
abline(lm(NWRchlorideData$Turbidity_NTU ~ NWRchlorideData$Date), col = red, lwd = 5)
text(16000, 60, "Multiple R-squared: 0.001052\nAdjusted R-squared: 0.0008293\nF-statistic: 4.73 on 1 and 4493 DF\np-value: 0.02969", cex = 0.5)


NWRchlorideDaysAbove250PPM <- read.csv("NorthwestRiverChlorideDaysAbove250PPM.csv", header = TRUE, sep = ",")
barplot(NWRchlorideDaysAbove250PPM$Number.of.Days.Chloride...250.PPM, names.arg = NWRchlorideDaysAbove250PPM$YEAR, ylab = "Num. of Days >250 PPM", main = "Northwest River Chloride")

library(imputeTS)
library(tseries)
timeseries <- ts(NWRchlorideData$Chloride_PPM, frequency = 365, start = 2006)
decomposedTS <- decompose(na_seadec(timeseries), type = "mult")
plot(decomposedTS)
adf.test(timeseries)
kpss.test(timeseries)


