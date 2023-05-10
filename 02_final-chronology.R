#GEOS597 Dendrochronology Methods Workshop
#spring semester 2023
#script02

#author: Isabel Gonzalez

#Script description: This script contains the chronology building from Pinus hartwegii samples collected in 2009 in the western
#highlands of Guatemala. Note that after exploring various detrending methods (see script01), the one considered to be the best fit was a 50 year spline.
#two chronologies are presented here: standard chronology and residual chronology. 
#Here I also explore the variance stabilization in the chronologies. 

################################################

### DETRENDING AND CHRONOLOGY BUILING ###
 ##        PINUS HARTWEGII            ##

################################################


#load libraries 

library(dplR)
library(graphics)
library(utils)
library(ggplot2)

#set working directory 
path = setwd('/Users/potatoknish/Documents/Arizona/spring2023/davidfworkshop/scripts')

#read data
data <- read.rwl('SQL_TWH_9.rwl')

#plot raw data: can choose also "seg"
plot(data, plot.type="spag")

#chronology from raw data
datacrn <- chron(data)
plot(datacrn, add.spline=TRUE, nyrs=30)

#we can also plot with ggplot## 
ggplot()+
  geom_hline(yintercept = 1,linetype="dashed") +
  geom_line(aes(x=time(datacrn),y=datacrn$std)) +
  geom_line(aes(x=time(datacrn),
                y=caps(datacrn$std,nyrs = 50)),
            color="darkred") +
  labs(x="Year",y="RWI") + theme_minimal()+ ggtitle("CHRONOLOGY FROM RAW DATA")

#plot only sample depth#
ggplot()+
  #geom_hline(yintercept = 1,linetype="dashed") +
  #geom_line(aes(x=time(datacrn),y=datacrn$std)) +
  geom_line(aes(x=time(datacrn),y=datacrn$samp.depth)) 
  #geom_line(aes(x=time(datacrn), +
               # y=caps(datacrn$std,nyrs = 50))
            #color="darkred") 
  labs(x="Year",y="Sample depth") + theme_minimal()+ ggtitle("SAMPLE DEPTH")
  
## Get some stadistics to explore our raw (withouth detrending) chronology ##
  
# Mean interseries correlation
  
#prewhitening the data = removes autocorrelation from raw tree ring width 
inter <- interseries.cor(data, n = NULL, prewhiten = TRUE, biweight = TRUE, method = "pearson")
write.csv(inter, file = "interseries_corrs.csv") #write to a file to save data. One can explore the file in the R Environment. 

################################################
### DETRENDING AND STANDARIZATION STARTS HERE ##
################################################

## I'll be using a 50 year spline to preserve the interannual variability. 
## I'll be using the default "ratios by division" feature of the dplR function. 

#50 year spline
data.rwi <- detrend(rwl = data, method = "Spline", nyrs = 50, f = 0.5)

#First, I'll be using the dplR function to build the standard chronology calculating the Tukey's Biweight Robust Mean of the detrending data. 
sdcrn <- chron(data)
str(sdcrn )

#plot standard chronology using ggplot without truncating data. 
ggplot()+
  geom_hline(yintercept = 1,linetype="dashed") +
  geom_line(aes(x=time(sdcrn),y=sdcrn$std)) +
  geom_line(aes(x=time(sdcrn),
                y=caps(sdcrn$std,nyrs = 50)),
            color="darkred") +
  labs(x="Year",y="RWI") + theme_minimal()+ ggtitle("50 YEAR SPLINE STANDARD CHRONOLOGY - PINUS HARTWEGII")

# Looking at the first plot, we can see that the sample depth is not representative since around 1860 - 1870. 
# Now, I'll be truncating the data and re build the standard chronology. 
# The data will be truncating to a sample depth of > 10. We can also be more strict and truncate the data to >10-15

sdcrn <- subset(sdcrn, samp.depth > 10)

#plot truncated standard chronology. 

ggplot()+
  geom_hline(yintercept = 1,linetype="dashed") +
  geom_line(aes(x=time(sdcrn),y=sdcrn$std)) +
  geom_line(aes(x=time(sdcrn),
                y=caps(sdcrn$std,nyrs = 30)),
            color="darkred") +
  labs(x="Year",y="RWI") + theme_minimal()+ ggtitle("50 YEAR SPLINE STANDARD CHRONOLOGY - PINUS HARTWEGII")


##residual chronology#

#50 year spline
chronres <- chron(data.rwi, prewhiten = TRUE)
str(chronres)

#truncat the residual chronology#

#50 year spline
phcrn <- subset(chronres, samp.depth > 10)

#plot residual chronology

ggplot()+
  geom_hline(yintercept = 1,linetype="dashed") +
  geom_line(aes(x=time(crnres),y=crnres$std)) +
  geom_line(aes(x=time(crnres),
                y=caps(crnres$std,nyrs = 50)),
            color="darkred") +
  labs(x="Year",y="RWI") + theme_minimal()+ ggtitle("50 YEAR SPLINE RESIDUAL CHRONOLOGY - PINUS HARTWEGII")

ggplot()+
  geom_hline(yintercept = 1,linetype="dashed") +
  geom_line(aes(x=time(crnres01),y=crnres01$std)) +
  geom_line(aes(x=time(crnres01),
                y=caps(crnres01$std,nyrs = 30)),
            color="darkred") +
  labs(x="Year",y="RWI") + theme_minimal()+ ggtitle("30 YEAR SPLINE RESIDUAL CHRONOLOGY - PINUS HARTWEGII")

#chronology uncertainty 

#Make chronology and get the mean 
#plus two standard errors of the yearly growth 

#50 YEAR SPLINE#
avgchron <- apply(data.rwi,1,mean,na.rm=TRUE)
se <- function(x){
  x2 <- na.omit(x)
  n <- length(x2)
  sd(x2)/sqrt(n)
}

avgron <- apply(data.rwi,1,se)

#plot

dat <- data.frame(yrs =as.numeric(names(avgchron)), 
                  std = avgchron,
                  lwr = avgchron - avgron*2,
                  upr = avgchron + avgron*2)
ggplot(dat,aes(x=yrs)) +
  geom_hline(yintercept = 1,linetype="dashed") +
  geom_ribbon(aes(ymin=lwr,ymax=upr),
              alpha=0.5,fill="blue") +
  geom_line(aes(y=std),col="grey30") +
  labs(x="YEAR",y="INDEX") + theme_minimal()+ ggtitle("50 YEAR SPLINE CHRONOLOGY UNCERTAINTY")

#write csv file with 50 year spline chronology 

write.csv(crnres, file = "phchron.csv")






