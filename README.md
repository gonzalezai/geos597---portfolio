# Geos597-portfolio

This repository contains three R codes as part of my portfolio for GEOS597 Dendrochronology Methods Workshop. 

Research goals: 

- Explore detrending methods on P. hartwegii chronology.

- Build residual and standard chronology. 

- Investigate climate signal in the chronology . 

Scripts descriptions: 

01_Exploring-Detrending.R: Throughout this script I explored the different detrending methods and applied them to Pinus hartwegii data. 

02_final-chronology.R: This script contains the chronology building from Pinus hartwegii samples collected in 2009 in the western highlands of Guatemala. Note that after exploring various detrending methods (see script01), the one considered to be the best fit was a 50 year spline. two chronologies are presented here: standard chronology and residual chronology. 

03_Exploring-ClimateSignal.R: This script was written to explore the climate signal with the tree chronology built on script 01 and 02. I'll be performing some correlations with instrumental precipitation and temperature data. I'll be using the residual chronology. 

The analysis was perfomed in three temporal resolutions:  

- monthly correlations
- seasonal correlations
- moving correlations 

Datasets in use: 

- precipitation: "ppt_gpcc" monthly gridded climate data. 
- temperature: "temp_cru2" monthly gridded climate data. 
- bothdatasets: "ppt_temp"

