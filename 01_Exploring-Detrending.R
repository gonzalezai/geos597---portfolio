#GEOS597 Dendrochronology Methods Workshop
#spring semester 2023
#script01

#author: Isabel Gonzalez

#Script description: Throughout this script I explored the different detrending methods and applied them to Pinus hartwegii data. 

################################################

###     EXPLORING DETRENDING         ###
##        PINUS HARTWEGII            ##

################################################

#load libraries 

library(dplR)
library(graphics)
library(utils)
library(ggplot2)
library(cowplot)
library(patchwork)
library(gridExtra)
library(tidyverse)

#set working directory 
#change directory
setwd('/Users/potatoknish/Documents/Arizona/spring2023/davidfworkshop/scripts')


#read data
data <- read.rwl('SQL_TWH_9.rwl')

#plot raw data: can choose also "seg"
plot(data, plot.type="spag")

#build chronology from raw data
datacrn <- chron(data)
plot(datacrn, add.spline=TRUE, nyrs=50)

### Explore raw data into time series object ###

ph <- read.rwl("SQL_TWH_9.rwl")

#convert to a time-series object
ph.ts <- ts(ph,start=as.numeric(rownames(scand)[1]))

#plot the raw data
ts.plot(ph.ts,ylab="Tree Ring Width (TRW)",xlab="Years")

#calculate a tree-ring chronology of the raw data, i.e., mean value for every year
ph.raw.crn <- ts(rowMeans(ph.ts,na.rm=T),start=1720)

#add the chronology to the plot with the raw data
lines(ph.raw.crn,col="orange")

#Look closely at the next periods: 1900-2000, 1900-1950, 1950-2000. 
ts.plot(scand.ts,ylab="TRW",xlab="Years",xlim=c(1950,2000));lines(scand.raw.crn,col="orange")

#Let's calculate the number of samples over time
ph.raw.repl <- ts(rowSums(ph.ts*0+1,na.rm=T),start=1720)
ts.plot(ph.raw.repl, ylab="Number of Samples Over Time",xlab="Years")

####################################

### DETRENDING STARTS HERE  ###
##      PINUS HARTWEGII     ##

####################################

####################################

# STANDARD CHRONOLOGY

####################################


## Detrending using power transforming the data ## 


#30 year spline: flexible curve to retain interannual variability 
data01.rwi <- detrend(rwl = powt(data), method = "Spline", nyrs = 30, f = 0.5)
                    

#50 year spline: flexible curve to retain interannual variability 
data02.rwi <- detrend(rwl = powt(data), method = "Spline", nyrs = 50, f = 0.5)

#Negative exponential: deterministic method trying to retain longer-term climate. 
data03.rwi <- detrend(rwl = powt(data), method = "ModNegExp")

# get the mean interseries correlation #

stats_data01 <- interseries.cor(data01.rwi, prewhiten=TRUE,
                                method="pearson")
stats_data02 <- interseries.cor(data02.rwi, prewhiten=TRUE,
                                method="pearson")
stats_data03 <- interseries.cor(data03.rwi, prewhiten=TRUE,
                                method="pearson")

## Build standard chronology from different applied methods ##
#Note: the three applied methods have been calculate by subtraction.

#The chron function builds the chronology by averaging the
#rows of the rwi using Tukey's biweight robust mean

crn01 <- chron(data01.rwi)
str(crn01)

crn02 <- chron(data02.rwi)
str(crn02)

crn03 <- chron(data03.rwi)
str(crn03)

#plot the standard chronology using "ggplot" and "gridExtra" libraries to plot three plots on the same panel. 
#plot(crn02, add.spline=TRUE, nyrs=30) this is the conventional dlpR function. 

p1 <- ggplot()+
  geom_hline(yintercept = 1,linetype="dashed") +
  geom_line(aes(x=time(crn01),y=crn01$std)) +
  geom_line(aes(x=time(crn01),
                y=caps(crn01$std,nyrs = 30)),
            color="turquoise4",size=1.5) +
  labs(x="Year",y="RWI") + theme_minimal()+ ggtitle("30 year spline")

p2 <- ggplot()+
  geom_hline(yintercept = 1,linetype="dashed") +
  geom_line(aes(x=time(crn02),y=crn02$std)) +
  geom_line(aes(x=time(crn02),
                y=caps(crn02$std,nyrs = 30)),
            color="lightslateblue",size=1.5) +
  labs(x="Year",y="RWI") + theme_minimal()+ ggtitle("50 year spline")

p3 <- ggplot()+
  geom_hline(yintercept = 1,linetype="dashed") +
  geom_line(aes(x=time(crn03),y=crn03$std)) +
  geom_line(aes(x=time(crn03),
                y=caps(crn03$std,nyrs = 30)),
            color="tomato1",size=1.5) +
  labs(x="Year",y="RWI") + theme_minimal()+ ggtitle("Negative exponential")

grid.arrange(p2, ncol = 1)

### Detrending by ratios ###

#Detrending by division (default) has the benefit of mitigating the commonly occuring heteroscedastic
#The heteroscedastic refers to the often observed relationship between the local variance and the mean
#(spread vs level)

data001.rwi <- detrend(rwl = data, method = "ModNegExp")
data002.rwi <- detrend(rwl = data, method = "Spline", nyrs = 30, f = 0.5)
data003.rwi <- detrend(rwl = data, method = "Spline", nyrs = 50, f = 0.5)

###build standard chronology from different applied methods###

chron001 <- chron(data001.rwi)
str(chron001)

chron002 <- chron(data002.rwi)
str(chron002)

chron003 <- chron(data003.rwi)
str(chron003)

#plot chronology using "ggplot" and "gridExtra" libraries to plot three plots -- 
#on the same panel. 

#plot(crn02, add.spline=TRUE, nyrs=30) this is the conventional dlpR function. 

p01 <- ggplot()+
  geom_hline(yintercept = 1,linetype="dashed") +
  geom_line(aes(x=time(chron001),y=chron001$std)) +
  geom_line(aes(x=time(chron001),
                y=caps(chron001$std,nyrs = 30)),
            color="turquoise4", size=1.5) +
  labs(x="Year",y="RWI") + theme_minimal()+ ggtitle("Negative Exponential")

p02 <- ggplot()+
  geom_hline(yintercept = 1,linetype="dashed") +
  geom_line(aes(x=time(chron002),y=chron002$std)) +
  geom_line(aes(x=time(chron002),
                y=caps(chron002$std,nyrs = 30)),
            color="tomato1", size=1.5) +
  labs(x="Year",y="RWI") + theme_minimal()+ ggtitle("30 year spline")

p03 <- ggplot()+
  geom_hline(yintercept = 1,linetype="dashed") +
  geom_line(aes(x=time(chron003),y=chron003$std)) +
  geom_line(aes(x=time(chron003),
                y=caps(chron003$std, nyrs = 30)),
            color="lightslateblue", size=1.5) +
  labs(x="Year",y="RWI") + theme_minimal()+ ggtitle("50 year spline")

grid.arrange(p02,ncol = 1)

#Now, I'll be truncating the data and re build the standard chronology. 
#The data will be truncating to a sample depth of > 10. 
#We can also be more strict and truncate the data to >10-15
#Note that here I'll be using the "ratio by division". 


crntrunc01 <- chron(detrend(data[chron001$samp.depth > 15,], method='ModNegExp'))
crntrunc02 <- chron(detrend(data[chron002$samp.depth > 15,], method = "Spline", nyrs = 30, f = 0.5))
crntrunc03 <- chron(detrend(data[chron003$samp.depth > 15,],method = "Spline", nyrs = 50, f = 0.5))   

#Now, we re plot the truncate chronology. 

p01 <- ggplot()+
  geom_hline(yintercept = 1,linetype="dashed") +
  geom_line(aes(x=time(crntrunc01),y=crntrunc01$std)) +
  geom_line(aes(x=time(crntrunc01),
                y=caps(crntrunc01$std,nyrs = 30)),
            color="turquoise4", size=1.5) +
  labs(x="Year",y="RWI") + theme_minimal()+ ggtitle("Negative exponential")

p02 <- ggplot()+
  geom_hline(yintercept = 1,linetype="dashed") +
  geom_line(aes(x=time(crntrunc02),y=crntrunc02$std)) +
  geom_line(aes(x=time(crntrunc02),
                y=caps(crntrunc02$std,nyrs = 30), size=3.5),
            color="tomato1", size=1.5) +
  labs(x="Year",y="RWI") + theme_minimal()+ ggtitle("30 year spline")

p03 <- ggplot()+
  geom_hline(yintercept = 1,linetype="dashed") +
  geom_line(aes(x=time(crntrunc03),y=crntrunc03$std)) +
  geom_line(aes(x=time(crntrunc03),
                y=caps(crntrunc03$std,nyrs = 30)),
            color="lightslateblue", size=1.5) +
  labs(x="Year",y="RWI") + theme_minimal()+ ggtitle("50 year spline")

#grid.arrange(p01, p02, p03, ncol = 3)
grid.arrange(p03, ncol = 1)


#####################

# RESIDUAL CHRONOLOGY

#####################

#here we explore the residual chronology, refering to the process were we detrend the raw data, 
#then remove the autocorrelation structure of the TRW and then we calculate the mean. 
#This is also known as the "pre whitening" method. 

#Let's re write our variables. 
#From now on, I'll be working only with the ratios by division data. 

res01.rwi <- detrend(rwl = data, method = "ModNegExp")
res02.rwi <- detrend(rwl = data, method = "Spline", nyrs = 30, f = 0.5)
res03.rwi <- detrend(rwl = data, method = "Spline", nyrs = 50, f = 0.5)

###build residual chronology from different applied methods###
#prewhiten = removing of autocorrelation

res01 <- chron(res01.rwi, prewhiten = TRUE)
str(res01)

res02 <- chron(res02.rwi, prewhiten = TRUE)
str(res02)

res03 <- chron(res03.rwi, prewhiten = TRUE)
str(res03)

#plot chronology using "ggplot" and "gridExtra" libraries to plot three plots -- 
#on the same panel. 

#plot(crn02, add.spline=TRUE, nyrs=30) this is the conventional dlpR function. 

px1 <- ggplot()+
  geom_hline(yintercept = 1,linetype="dashed") +
  geom_line(aes(x=time(res01),y=res01$res)) +
  geom_line(aes(x=time(res01),
                y=caps(res01$res,nyrs = 30)),
            color="turquoise4", size=1.5) +
  labs(x="Year",y="RWI") + theme_minimal()+ ggtitle("Negative Exponential")

px2 <- ggplot()+
  geom_hline(yintercept = 1,linetype="dashed") +
  geom_line(aes(x=time(res02),y=res02$res)) +
  geom_line(aes(x=time(res02),
                y=caps(res02$res,nyrs = 30)),
            color="tomato1", size=1.5) +
  labs(x="Year",y="RWI") + theme_minimal()+ ggtitle("30 year spline")

px3 <- ggplot()+
  geom_hline(yintercept = 1,linetype="dashed") +
  geom_line(aes(x=time(res03),y=res03$res)) +
  geom_line(aes(x=time(res03),
                y=caps(res03$res,nyrs = 30)),
            color="lightslateblue", size=1.5) +
  labs(x="Year",y="RWI") + theme_minimal()+ ggtitle("50 year spline")

grid.arrange(px1,ncol = 1)

#Now, I'll be truncating the data and re build the residual chronology. 
#The data will be truncating to a sample depth of > 10. 
#We can also be more strict and truncate the data to >10-15
#Note that here I'll be using the "ratio by division". 


restrunc01 <- chron(detrend(data[res01$samp.depth > 15,], method='ModNegExp'))
restrunc02 <- chron(detrend(data[res02$samp.depth > 15,], method = "Spline", nyrs = 30, f = 0.5))
restrunc03 <- chron(detrend(data[res03$samp.depth > 15,],method = "Spline", nyrs = 50, f = 0.5))   

#Now, we re plot the truncate chronology. 

py1 <- ggplot()+
  geom_hline(yintercept = 1,linetype="dashed") +
  geom_line(aes(x=time(restrunc01),y=restrunc01$res)) +
  geom_line(aes(x=time(restrunc01),
                y=caps(restrunc01$res,nyrs = 30)),
            color="turquoise4", size=1.5) +
  labs(x="Year",y="RWI") + theme_minimal()+ ggtitle("Negative Exponential")

py2 <- ggplot()+
  geom_hline(yintercept = 1,linetype="dashed") +
  geom_line(aes(x=time(restrunc02),y=restrunc02$res)) +
  geom_line(aes(x=time(restrunc02),
                y=caps(restrunc02$res,nyrs = 30), size=3.5),
            color="tomato1", size=1.5) +
  labs(x="Year",y="RWI") + theme_minimal()+ ggtitle("30 year spline")

py3 <- ggplot()+
  geom_hline(yintercept = 1,linetype="dashed") +
  geom_line(aes(x=time(restrunc03),y=restrunc03$res)) +
  geom_line(aes(x=time(restrunc03),
                y=caps(restrunc03$res,nyrs = 30)),
            color="lightslateblue", size=1.5) +
  labs(x="Year",y="RWI") + theme_minimal()+ ggtitle("50 year spline")

grid.arrange(py3, ncol = 1)


### USING SUBSAMPLE SIGNAL STRENGHT (SSS) ###

#here we'll try to truncate the data using the subsample signal strenght (sss) 
#we'll keep on working only with the 50 year spline "res03" data frame, which is a data frame with ratios by division. 

#first we define our variables 
dataIds <- autoread.ids(data)
sssThresh <- 0.85
dataSSS <- sss(res03.rwi, dataIds)
yrs <- time(data)
yrCutoff <- max(yrs[dataSSS < sssThresh])

ggplot() +
  geom_rect(aes(ymin=-Inf, ymax=Inf,xmin=-Inf,xmax=yrCutoff),
            fill="orange",alpha=0.3) +
  annotate(geom = "text",y=1.5,x=1725,label="SSS < 0.85")+
  geom_hline(yintercept = 1,linetype="dashed") +
  geom_line(aes(x=yrs,y=res03$res)) +
  labs(x="Year",y="RWI") + theme_minimal()+ ggtitle("50 year spline - residual")

#Now we truncate the res02.rwi data frame and then we plot. 

dataRwlSSS <- data[dataSSS > sssThresh,]
dataRwiSSS <- detrend(dataRwlSSS, method = "Spline", nyrs = 50, f = 0.5)

#we build residual chronology. 
#now our chronology files will be "crnresid" for the residual chronology. 
crnresid <- chron(dataRwiSSS,  prewhiten = TRUE)
str(crnresid)

#plot two figures on the same panel. 

po1 <- ggplot() +
  geom_hline(yintercept = 1,linetype="dashed") +
  geom_line(aes(x=time(crnresid),y=crnresid$res)) +
  geom_line(aes(x=time(crnresid),
                y=caps(crnresid$res,nyrs = 30)),
            color="tomato1", size=2) +
  labs(x="Year",y="RWI") + theme_minimal_grid() + ggtitle("50 year spline - residual chronology")

grid.arrange(po1, ncol = 1)


###########################

# CHRONOLOGY UNCERTAINTY

###########################

#Make chronology and get the mean 
#plus two standard errors of the yearly growth quite simply

AvgCrn01 <- apply(dataRwiSSS,1,mean,na.rm=TRUE)
se <- function(x){
  x2 <- na.omit(x)
  n <- length(x2)
  sd(x2)/sqrt(n)
}

AvgCrnSE01 <- apply(dataRwiSSS,1,se)

#plot

dat <- data.frame(yrs =as.numeric(names(AvgCrn01)), 
                  std = AvgCrn01,
                  lwr = AvgCrn01 - AvgCrnSE01*2,
                  upr = AvgCrn01 + AvgCrnSE01*2)

ggplot(dat,aes(x=yrs)) +
  geom_hline(yintercept = 1,linetype="dashed") +
  geom_ribbon(aes(ymin=lwr,ymax=upr),
              alpha=0.5,fill="blue") +
  geom_line(aes(y=std),col="grey30") +
  labs(x="Year",y="RWI") + theme_minimal()


#panels configurations
p123 <- po1 | (po2 / po3)
p123 + plot_annotation(tag_levels = "I") 

#write csv file with 50 year spline chronology for the residual chronology. 
#standard chronology = "crnresid.csv"
#residual chronology = 

write.csv(crnresid, file = "crnresid.csv")
#### THIS IS THE END  ####




