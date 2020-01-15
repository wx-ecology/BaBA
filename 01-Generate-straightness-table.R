################ Barrier Behavior Analysis (BaBA): straightness table  ############
## Input ##
# Fence shp, movement points (with consistent time intervals since the behavior type assumption 

## output ##
## straightness table 

## Note ##
# this script sometimes take a long time to run. Might consider using parallel looping if your data size is big

## last updated ##
# Dec 19, 2019

###############################################################################


# Set up the R session ------------------------------------------------------------------

## clear environment ####
rm(list = ls())

## load libraries ####

### general
library(dplyr)

### spatial analysis
library(rgdal)
library(rgeos)
library(sp)
library(raster)

#### trajectory analysis
library(adehabitatLT)

## set your working directory ####
setwd(".")


# Set up Parameters Parameters --------------------------------------------------------------

## define input and ouput paths ####

input_file_path <- "data/PRONG_TEST.csv" # path to your movement data (relative to your working directory)

ouput_file_path <- "results/straigthness_output.csv" # path you want your output data to be saved, including .csv extension (relative to your working directory)

## define the CRS of your data ####
target.crs <- "+proj=utm +zone=12 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"

## define BaBA parameters ####

# define time interval (in hours) of the movement data
interval <- 2 

## ** parameters after this should be determined by local management requirements and characteristics of the target species. ** ##

# define, based on the data time interval and animal ecology, the maximum encounter duration that you'd call it a "bounce" or "Quick Cross".Aka, what is quick for you? note: b.hours + interval is the minimum duration that straightness table will calculate
b.hours <- 4
b <- b.hours/interval

# define the minimum encounter duration in the burst that you'd call the encounter event a "trapped" condition. Must be divisible by interval. Aka, what is too long for you? p.hours is also the maximum duration that straightness table will calculate.
p.hours <- 36
p <- p.hours/interval


# Function ----------------------------------------------------------------


# create the function that calculating straightness. Input must be a dataframe with Easting and Northing. 

strtns <- function(mov.seg) {
  
  pts <- cbind(mov.seg$Easting, mov.seg$Northing)
  pts.sp <- SpatialPointsDataFrame(pts, mov.seg, proj4string = CRS(target.crs))
  traj <- as.ltraj(xy =  pts, date = mov.seg$date, id = as.character(mov.seg$Location.ID))
  
  #moving distance from first pt to last pt in the burst
  traj.dist <- sqrt(
    (traj[[1]]$x[1] - traj[[1]]$x[nrow(traj[[1]])]) * (traj[[1]]$x[1] - traj[[1]]$x[nrow(traj[[1]])]) +
      (traj[[1]]$y[1] - traj[[1]]$y[nrow(traj[[1]])]) * (traj[[1]]$y[1] - traj[[1]]$y[nrow(traj[[1]])]) 
  )
  
  #sum of all step lengths
  traj.lgth <- sum(traj[[1]]$dist, na.rm = TRUE)
  
  #straightness ranges from 0 to 1. More close to 0 more sinuous it is.
  straightness <- traj.dist/traj.lgth
  
  return(straightness)
}

# Data --------------------------------------------------------------------


# read in movement data
#ideally, the movement data should not have missing point. This trial file does have missing points.

movement.df.all <- read.csv(input_file_path)  # read in your movement data. 
movement.df.all$date <- as.POSIXct(strptime(as.character(movement.df.all$date),"%m/%d/%Y %H:%M")) 
movement.df.all <- movement.df.all[(!is.na(movement.df.all$date))&(!is.na(movement.df.all$Easting)),]


# Analysis ----------------------------------------------------------------

# create an empty dataframe 
animal.stn.df <- data.frame(AnimalID = integer(), window.size = numeric(), Date = character(), Straightness = numeric())

# run a for loop for each individual
system.time(for (i in unique(movement.df.all$Location.ID)) {
  
  movement.df.i <- movement.df.all[movement.df.all$Location.ID == i,]
  
  # e.g. interval = 2, b = 2, p = 36. Calculating moving window with size 3 - 24 (including 24)
  for (ii in (b+1):min(p, nrow(movement.df.i))) { # ii is the window size (duration of trajectories that straightness will be calculated)
    straightness.ii <- vector()
    date.ii <- character()
    
    for (iii in seq(1, (nrow(movement.df.i)-ii), by = 2)) {  # can change "by" for different sampling rate to calculate strightness to save some time
      mov.seg.iii <- movement.df.i[iii:(iii+ii),]
      # save the starting time of the trajectory
      date.ii <-  c(date.ii, as.character(strftime(mov.seg.iii$date[1], "%Y-%m-%d %H:%M:%S")))
      # calculate and save straightness of this traijectory
      straightness.ii <- c(straightness.ii, strtns(mov.seg.iii))
    }
    
    n <- length(straightness.ii)
    rows.i <- data.frame(AnimalID = rep(i, n), window.size = rep(ii, n), Date = date.ii, Straightness = straightness.ii)
    animal.stn.df <- rbind(animal.stn.df, rows.i)
  }
})


# Save output -------------------------------------------------------------

write.csv(animal.stn.df, paste(getwd(), ouput_file_path, sep = "/"), row.names = F)
