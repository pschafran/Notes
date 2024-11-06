# Adapted from https://www.polarmicrobes.org/analyzing-flow-cytometry-data-with-r/

library("flowCore")
library("ggplot2")
library("ggpubr")
library("dplyr")
setwd("~/your/path/here/")  # directory your .fcs files are in
f.name <- 'your_file_here.fcs'  # name of the file you want to analyze, file must have extension ".fcs"

fcm <- read.FCS(f.name)
fcm <- as.data.frame((exprs(fcm)))

# View whole range of data incl. noise
ggplot(fcm, aes(x=fcm$`YL1-A`, y=fcm$`YL1-H`) ) +
  geom_bin2d(bins = 100) +
  scale_fill_continuous(type = "viridis") +
  scale_x_continuous(breaks = scales::pretty_breaks(n=15)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n=15)) +
  xlab("YL1-A") +
  ylab("YL1-H") +
  theme_bw()

## eliminate noise that is below or equal to PI fluorescence thresholds
pia.ll <- 30000 # YL1-A lower limit
pia.ul <- 200000 # YL1-A upper limit
pih.ll <- 20000 # YL1-H lower limit 
pih.ul <- 200000 # YL1-H upper limit

fcm$`YL1-A`[fcm$`YL1-A` <= pia.ll| fcm$`YL1-A` >= pia.ul | fcm$`YL1-H` <= pih.ll | fcm$`YL1-H` >= pih.ul] <- NA
fcm <- na.omit(fcm)

# 2D Density Plot of YL1-A (PI area) vs YL1-H (PI height)
ggplot(fcm, aes(x=fcm$`YL1-A`, y=fcm$`YL1-H`) ) +
  geom_bin2d(bins = 100) +
  #scale_fill_gradient(low = "white", high = "black") +
  scale_fill_continuous(type = "viridis") +
  scale_x_continuous(breaks = scales::pretty_breaks(n=10)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n=10)) +
  xlab("YL1-A") +
  ylab("YL1-H") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# Histogram of YL1-A (PI area)
histA <- ggplot(fcm, aes(x=fcm$`YL1-A`) ) +
  geom_histogram(bins = 1000) +
  scale_x_continuous(breaks = scales::pretty_breaks(n=10)) +
  xlab("YL1-A") +
  ylab("Count") +
  theme_bw()

# Histogram of YL1-H (PI height)
histH <- ggplot(fcm, aes(x=fcm$`YL1-H`) ) +
  geom_histogram(bins = 1000) +
  scale_x_continuous(breaks = scales::pretty_breaks(n=10)) +
  xlab("YL1-H") +
  ylab("Count") +
  theme_bw()

#Plot histograms
ggarrange(histA, histH, ncol = 1)

# Calculate Median of Each Peak (based on 2D density plot or histogram)
pia.standard.ll <- 40000  # YL1-A lower limit of peak
pia.standard.ul <- 60000 # YL1-A upper limit of peak
pih.standard.ll <- 40000  # YL1-H lower limit of peak (if unsure, use same as pia.ll)
pih.standard.ul <- 60000  # YL1-H upper limit of peak (if unsure, use same as pia.ul)
median(fcm$`YL1-A`[fcm$`YL1-A` >= pia.standard.ll & fcm$`YL1-A` <= pia.standard.ul & fcm$`YL1-H` >= pih.standard.ll & fcm$`YL1-H` <= pih.standard.ul])

pia.sample.ll <- 140000  # YL1-A lower limit of peak
pia.sample.ul <- 160000 # YL1-A upper limit of peak
pih.sample.ll <- 120000  # YL1-H lower limit of peak (if unsure, use same as pia.ll)
pih.sample.ul <- 150000  # YL1-H upper limit of peak (if unsure, use same as pia.ul)
median(fcm$`YL1-A`[fcm$`YL1-A` >= pia.sample.ll & fcm$`YL1-A` <= pia.sample.ul & fcm$`YL1-H` >= pih.sample.ll & fcm$`YL1-H` <= pih.sample.ul])



# Count number of nuclei in the target range
length(fcm$`YL1-A`[fcm$`YL1-A` >= pia.standard.ll & fcm$`YL1-A` <= pia.standard.ul & fcm$`YL1-H` >= pih.standard.ll & fcm$`YL1-H` <= pih.standard.ul])
length(fcm$`YL1-A`[fcm$`YL1-A` >= pia.sample.ll & fcm$`YL1-A` <= pia.sample.ul & fcm$`YL1-H` >= pih.sample.ll & fcm$`YL1-H` <= pih.sample.ul])

