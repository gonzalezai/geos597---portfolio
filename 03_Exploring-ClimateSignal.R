#GEOS597 Dendrochronology Methods Workshop
#spring semester 2023
#script03

#author: Isabel Gonzalez

#Script description: This script was written to explore the climate signal with the
#tree chronology built on script 01 and 03. 

#I'll be performing some correlations with instrumental precipitation and temperature data. 
#I'll be using the residual chronology. 

#The analysis will be perfomed in three temporal resolutions:  

#- monthly correlations
#- seasonal correlations
#- moving correlations 

#Datasets in use: 
#precipitation: "ppt_gpcc" monthly gridded climate data. 
#temperature: "temp_cru2" monthly gridded climate data. 
#bothdatasets: "ppt_temp"

#LET'S RE BUILD CHRONOLOGY TO RE DO STEPS AND BE EXTRA CAUTIOUS# 

#load libraries
library(dplR)
library(treeclim)
library(TRADER)
library(graphics)
library(utils)
library(ggplot2)

#set working directory 
setwd('/Users/potatoknish/Documents/Arizona/spring2023/davidfworkshop/scripts')

#read data
data <- read.rwl('SQL_TWH_9.rwl')

#rebuild chronology
data.rwi <- detrend(rwl = data, method = "Spline", nyrs = 50, f = 0.5)
chron <- chron(data.rwi, prewhiten = TRUE)
str(chron)

#Truncate data using SSS
dataIds <- autoread.ids(data)
sssThresh <- 0.85
dataSSS <- sss(data.rwi, dataIds)
yrs <- time(data)
yrCutoff <- max(yrs[dataSSS < sssThresh])

ggplot() +
  geom_rect(aes(ymin=-Inf, ymax=Inf,xmin=-Inf,xmax=yrCutoff),
            fill="orange",alpha=0.3) +
  annotate(geom = "text",y=1.5,x=1725,label="SSS < 0.85")+
  geom_hline(yintercept = 1,linetype="dashed") +
  geom_line(aes(x=yrs,y=chron$res)) +
  labs(x="Year",y="RWI") + theme_minimal()+ ggtitle("50 year spline - residual")

#Now we truncate the data.rwi data frame and then we plot. 

dataRwlSSS <- data[dataSSS > sssThresh,]
dataRwiSSS <- detrend(dataRwlSSS, method = "Spline", nyrs = 50, f = 0.5)

#we build residual chronology. 
#now our chronology files will be "crnresid" for the residual chronology. 
phcrn <- chron(dataRwiSSS,  prewhiten = TRUE)
str(phcrn)

##################################

### CLIMATE  ###
##  ANALYSIS ##

##################################

#recall: "crnresid" chronology, "ppt_gpcc", "temp_cru2", "ppt_temp"

precip <- read.csv("ppt_gpcc.csv")
temp <- read.csv("temp_cru2.csv")

#this file contains precip and temp data
clim <- read.csv("ppt_temp.csv")

#######################

#MONTHLY CORRELATIONS

#######################

#Can use = "static", "moving", "evolving"
#I'll be using (2:12) based on past results on Abies being sensitive to early growing season moisture. 
#Then I'll try using past growing season. From past January to the present November.  

precipcorr <- dcc(chrono = phcrn, climate = precip, selection = -2:12, 
                method = "correlation", dynamic = "static", win_size = 35, win_offset = 1, start_last = FALSE,
                timespan = NULL, var_names = "Precipitation", ci = 0.05, boot = "std", sb = FALSE) 

coef <- coef(precipcorr)
plot(precipcorr) 

precipres <- dcc(chrono = phcrn, climate = precip, selection = -2:12, 
                  method = "response", dynamic = "static", win_size = 35, win_offset = 1, start_last = FALSE,
                  timespan = NULL, var_names = "Precipitation", ci = 0.05, boot = "std", sb = FALSE)

coefr <- coef(precipres)
plot(precipres) 

#observations: correlation from the present February to December show "strong" signal in february, with a weaker signal in March. 
#showing a negative correlation with precipitation in May.
#With the response function, we also see a negative relationship with May. 
#Correlation taking the past growing season is showing a strong relationship again with the current February and March and a negative
#correlation with current May (onset of the rainy season over the region). 

#Let's try with temperature. Again same seasons. (2:11), (-02:12)

tempcorr <- dcc(chrono = phcrn, climate = temp, selection = -02:11, 
                  method = "correlation", dynamic = "static", win_size = 35, win_offset = 1, start_last = FALSE,
                  timespan = NULL, var_names = "Temperature", ci = 0.05, boot = "std", sb = FALSE)

coef <- coef(tempcorr)
plot(tempcorr) 

#Again, #let's change the method function to "response"

tempres <- dcc(chrono = phcrn, climate = temp, selection = 2:11, 
                method = "response", dynamic = "static", win_size = 35, win_offset = 1, start_last = FALSE,
                timespan = NULL, var_names = "Temperature", ci = 0.05, boot = "std", sb = FALSE) 

coef <- coef(tempres)
plot(tempres) 

#observations: I don't see any significant correlation with temperature and tree ring width. 
#this results imply sensitivity specifically to precipitation in the hartwegii chronology. 

#######################

#MOVING CORRELATIONS

#######################

#we will try 35, 25 and 15 year moving correlation. [2:12]

ppmov <- dcc(chrono = phcrn, climate = precip, selection = -2:12, 
                  method = "correlation", dynamic = "moving", win_size = 25, win_offset = 1, start_last = FALSE,
                  timespan = NULL, var_names = "Precip", ci = 0.05, boot = "std", sb = FALSE)

coef <- coef(ppmov) 
plot(ppmov)

##Observations##

#first thing that call my attention, specifically for the 35 year correlation is the strong correlation
#with February precipitation and the negative correlation with May. 
#correlation with precipitation exceed the 0.4 values with a positive behavior. 

#Changing the "win_size" to a 25 year moving correlation, we notice some strong negavtive
#correlations with the previous June and August precipitation and then again, a strong 
#positive with the current February and March. 

#let's see what temperature is doing. 

tempmov <- dcc(chrono = phcrn, climate = temp, selection = 02:12, 
             method = "correlation", dynamic = "moving", win_size = 25, win_offset = 1, start_last = FALSE,
             timespan = NULL, var_names = "Temp", ci = 0.05, boot = "std", sb = FALSE) #this is the main function in treeclim

coef <- coef(tempmov) 
plot(tempmov) 

##Observations##

#Current December precipitation signal really stands out. What's happening there?
#Looking at the previous growing season, we see a strong signal with December. 
#Changing the "win_size" to a 25 year moving correlation, I noticed a positive correlation
#with the previous February and previous March but the consistent relationship persists for current December. 


#######################

#SEASONAL CORRELATIONS

#######################

## seasonal correlations with seascorr ##
#here I'll use the "clim" variable again. 

sc <- seascorr(chrono = phcrn, climate = clim,var_names = NULL, timespan = NULL, 
         complete = 12, season_lengths = c(1, 2, 3),primary = 2,
         secondary = 1,ci = 0.05)
plot(sc)

##Observations##

#Strong seasonal correlationtion with precipitation with FMA. This results are similar to the ones computed from Abies Guatemalensis in the highlands. 
#Next step: performing PCA. 