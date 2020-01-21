################ Barrier Behavior Analysis (BaBA): classification  ##########################
## Input ##
# Fence shp, movement points (with the same time intervals since the behavior type assumption
# Straighness table generated from step 01

## output ##
## 1. encounter event table; 2. step 1 classification table (bounce/trapped/quick cross); 3. Final classification table

## data requirement ##
# 1. Movement coordinates - projected coordinations (m). X/Y coordinate column name: Easting/Northing. Ideally, the movement data should not the less missing points the better.
# 2. Fence polylines - as accurate as possible

## last updated ##
# Dec 19, 2019

#############################################################################################

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

### trajectory analysis
library(adehabitatLT)

### for parallel multi-core calculation
library(foreach)
library(doParallel)
library(doSNOW)

# Set up Parameters --------------------------------------------------------------


## define input and ouput paths ####

movement_data_input_file_path <- "data/PRONG_TEST.csv" # path to your movement data (relative to your working directory)

fence_data_input_file_path <- "data/Fence_AOI_FINAL_Dec2019.shp" # path to your linear structure shapefile (relative to your working directory)

straightness_reference_file_path <- "results/straigthness_output.csv"

ouput_encounter_events_file_path <- "results/Encounter_Events" # path you want your encounter event data frame to to be saved, NOTE: an extension indicating the BaBA parameters will be attached to that

ouput_classification_events_file_path <- "results/Classification_Events" # path you want the classification of your events to be saved, NOTE: an extension indicating the BaBA parameters will be attached to that


## define the CRS of your data ####
target.crs <- "+proj=utm +zone=12 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"


## define BaBA parameters ####

### Fence buffer distance in meters
#### this might be the most important parameter to set the BaBA.
#### the buffer distance will affect numbers and durations of trajectory. For different species, barrier effect distance might be different.
#### you can decide the distance either by experience, or testing a range of distance and compare the results.
FB.dist <- 50

### define time interval (in hours) of the movement data
interval <- 2


### ** parameters after this should be determined by local management requirements and characteristics of the target species. ** ##

#### define, based on the data time interval and animal ecology, the maximum encounter duration that you'd call it a "bounce" or "Quick Cross".Aka, what is quick for you? note: b.hours + interval is the minimum duration that straightness table will calculate
b.hours <- 4
b <- b.hours / interval

#### define the minimum encounter duration in the burst that you'd call the encounter event a "trapped" condition. Must be divisible by interval. Aka, what is too long for you? p.hours is also the maximum duration that straightness table will calculate.
p.hours <- 36
p <- p.hours / interval


#### When differenciating trace/back-n-forth from "normal movement",
# we compare the straightness of the encounter event to the average straightness around the time the encounter event happens.
# how long would you set this window? Current default is 7 days.
ave.window <- 7
half.window <- ave.window / 2


#### Sometimes when barriers are densely distributed, a trace behavior might actually cross some barrier lines with two points are connected together
# here set the maximum number of crosses allowed for considering a trajectory being tracing behavior
max.cross <- 4

##### for high temporal interval movement data
## tolerance parameter. If "a" is point in the buffer, "b" is a point outside of the buffer
# x is the number of b that is allowed in between of a to allow the point series be considered as a continuous encounter event
tolerance <- 0

# Function ----------------------------------------------------------------



## create function that extracts movement segment between the time of m and predetermined time lap b*interval.
movement.segment.b <- function(pt1, pt2) {
  segments <- movement.df[which(movement.df$ptsID >= pt1 - 1 &
                                  movement.df$ptsID <= pt2 + 1), ]
  seg.line <- Lines(Line(cbind(segments$coords.x1, segments$coords.x2)), 
                    ID = segments$date[1])
  
  segments.sp <- SpatialLines(list(seg.line), proj4string = CRS(target.crs))
  
  return(segments.sp)
}

## create the function that calculating straightness. Input must be a dataframe with Easting and Northing.
strtns <- function(mov.seg) {
  
  pts <- cbind(mov.seg$Easting, mov.seg$Northing)
  pts.sp <- SpatialPointsDataFrame(pts, mov.seg, proj4string = CRS(target.crs))
  traj <- as.ltraj(xy =  pts,
                   date = mov.seg$date,
                   id = as.character(mov.seg$Location.ID))
  
  # moving distance from first pt to last pt in the burst
  traj.dist <- sqrt(as.numeric((traj[[1]]$x[1] - traj[[1]]$x[nrow(traj[[1]])])) *  as.numeric((traj[[1]]$x[1] - traj[[1]]$x[nrow(traj[[1]])])) + as.numeric((traj[[1]]$y[1] - traj[[1]]$y[nrow(traj[[1]])])) * as.numeric((traj[[1]]$y[1] - traj[[1]]$y[nrow(traj[[1]])])))
  
  #sum of all step lengths
  traj.lgth <- sum(traj[[1]]$dist, na.rm = TRUE)
  
  #straightness ranges from 0 to 1. More close to 0 more sinuous it is.
  straightness <- traj.dist / traj.lgth
  
  return(straightness)
}

# pick the minimum non-negative number
min.nonneg <- function(x)  min(x[x >= 0])



# Data --------------------------------------------------------------------
## read in data ####
## read in fence shapefile
fence.sp <- readOGR(fence_data_input_file_path)


## read in movement data
movement.df.all <- read.csv(movement_data_input_file_path)


## prepare data ####

### create a buffer around de fences
fence.buffer <- raster::buffer(fence.sp, width = FB.dist)

### format date in movement data
movement.df.all$date <- as.POSIXct(strptime(as.character(movement.df.all$date),"%m/%d/%Y %H:%M")) #change the format based on the data

### remove fixes with missing dates and missing easting
movement.df.all <- movement.df.all[(!is.na(movement.df.all$date)) & (!is.na(movement.df.all$Easting)), ]


### add point ID by individual

movement.df.all$ptsID <- numeric(nrow(movement.df.all))

for (i in unique(movement.df.all$Location.ID)) {
  mov.seg.i <- movement.df.all[movement.df.all$Location.ID == i, ]
  movement.df.all[movement.df.all$Location.ID == i, ]$ptsID <-
    seq(nrow(mov.seg.i))
}


### make movement pts into spatial point dataframe
xy <- cbind(movement.df.all$Easting, movement.df.all$Northing)

movement.sp.all <- SpatialPointsDataFrame (coords = xy,
                                           data = movement.df.all,
                                           proj4string = CRS(target.crs))

### read in straightness table from last step
animal.stn.df <- read.csv(straightness_reference_file_path)


# Analysis ----------------------------------------------------------------

## ---- classification step 1: generate encountering event dataframe ---- ####

### set up parallel looping
cores <- detectCores()
cl <- makeSOCKcluster(cores[1] - 1) #to not overload your computer
registerDoSNOW(cl)

### set up progress bar
pb <- txtProgressBar(max = 100, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

### generate the event dataframes

encounter.df = foreach (
  i = unique(movement.df.all$Location.ID),
  .combine = rbind,
  .options.snow = opts,
  .packages = c('raster', 'sp', 'rgdal', 'rgeos', 'adehabitatLT', 'sp', 'dplyr')
) %dopar% {
  # ---------------------  build dataframe, burst ID = first timestamp of the traj -------------------------
  movement.sp <- movement.sp.all[movement.sp.all$Location.ID == i,]
  movement.df <- as.data.frame(movement.sp)
  #extract points that fall inside of the buffer
  encounter.sp <- raster::intersect(movement.sp, fence.buffer)
  #Creat a data frame with encounter event marked as burst ID
  encounter.df <- as.data.frame(encounter.sp)
  encounter.df <-
    encounter.df[which(!is.na(encounter.df$coords.x1)), ]
  encounter.df$date <-
    as.POSIXct(strptime(as.character(encounter.df$date), "%Y-%m-%d %H:%M"))
  
  encounter.df$burst <- vector(length = nrow(encounter.df))
  for (ii in 2:nrow(encounter.df)) {
    # first group events based on time
    #if the time intervals between two points are within the set tolerance time, they are likely belog to the same encounter event, marked as "m"
    if (difftime(encounter.df$date[ii], encounter.df$date[ii - 1], units = "hours") <= (tolerance +
                                                                                        1) * interval) {
      encounter.df$burst[ii - 1] <- "m"  #mark the cell for the burst
      if (ii == nrow(encounter.df)) {
        encounter.df$burst[ii] <- "m"
      }
    }
    else {
      #all points before this point that are marked as "m" are in the same burst
      encounter.df$burst[ii - 1] <- "m"
      burst.list <- encounter.df[which(encounter.df$burst == "m"), ]
      encounter.df[which(encounter.df$burst == "m"), ]$burst <-
        format(burst.list$date[1], "%Y-%m-%d %H:%M") #name each burst using the starting time stamp
      if (ii == nrow(encounter.df)) {
        encounter.df$burst[ii] <- "m"
      }
    }
    if (ii == nrow(encounter.df)) {
      burst.list <- encounter.df[which(encounter.df$burst == "m"), ]
      if (nrow(burst.list) > 0) {
        encounter.df[which(encounter.df$burst == "m"), ]$burst <-
          format(burst.list$date[1], "%Y-%m-%d %H:%M")
      }
    }
  }
  # clean the data frame
  encounter.df <-
    encounter.df[which(!is.na(encounter.df$burst)), ]  #all points that are in buffer in one dataframe
  encounter.df
}

### save encouter dataframe
write.csv(encounter.df, paste0(ouput_encounter_events_file_path, "_I", interval,
                               "_FB", FB.dist, "_B",b.hours, "_P", p.hours, ".csv"))

### close progress bar
close(pb)

### stop cluster
stopCluster(cl)

## ---- classification step 2: classify bounce, quick cross, and trap ---- ####
###  bounce, quick cross, and trap are all based on duration 

### set up parallel looping
cores <- detectCores()
cl <- makeSOCKcluster(cores[1] - 1) #to not overload your computer
#registerDoParallel(cl)
registerDoSNOW(cl)

### set up progress bar
pb <- txtProgressBar(max = 100, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

### start the loop
event.df = foreach (
  i = unique(encounter.df$Location.ID),
  .combine = rbind,
  .options.snow = opts,
  .packages = c('raster', 'sp', 'rgdal', 'rgeos', 'adehabitatLT', 'sp', 'dplyr')
) %dopar% {
  # focus on the i animal
  movement.df <-
    as.data.frame(movement.sp.all[movement.sp.all$Location.ID == i, ])
  encounter.df.i <- encounter.df[encounter.df$Location.ID == i,]
  # empty data frame
  l <- length(unique(encounter.df.i$burst))
  event.df <- data.frame(
    AnimalID = character(l),
    burstID = character(l),
    easting = numeric(l),
    northing = numeric(l),
    duration = numeric(l),
    cross = numeric(l),
    straightness = numeric(l),
    eventTYPE = character(l),
    stringsAsFactors = FALSE
  )
  event.df$burstID <- unique(encounter.df.i$burst)
  event.df$AnimalID <- rep(i, l)
  # looping through each encounter event
  for (ii in unique(encounter.df.i$burst)) {
    burst.i <- encounter.df.i[encounter.df.i$burst == ii, ]
    
    start.time <- burst.i[1,]$date
    end.time <- burst.i[nrow(burst.i),]$date
    
    event.df[which(event.df$burst == ii),]$duration <- difftime (end.time, start.time, units = "hours")
    event.df[which(event.df$burst == ii),]$easting <- burst.i$coords.x1[1]
    event.df[which(event.df$burst == ii),]$northing <- burst.i$coords.x2[1]
    
    #first for short encounter
    #no more than b*interval H, only spend small amount of time in this burst
    if (difftime (end.time, start.time, units = "hours") <= b * interval) {
      pt.first <- burst.i[1,] #first point in the burst
      pt.last <- burst.i[nrow(burst.i),]
      #extract movement segment with one point before and one point after the segmentation
      mov.seg.i <-
        movement.segment.b(pt.first$ptsID, pt.last$ptsID)
      if (nrow(coordinates(mov.seg.i)[[1]][[1]]) <= b) {
        #which means this event is at the end of the individual's traj. Not enough points to tell catagories.
        event.df[which(event.df$burst == ii),]$eventTYPE <-
          "unknown"
        next
      }
      int.num <- length(gIntersection(mov.seg.i, fence.sp))
      ############# Need to acknowledge possibility for mistakes here since the links between points are not real movement routes
      ############# so the intersecting points might not actually mean crossing.
      if (int.num == 0) {
        event.df[which(event.df$burst == ii),]$eventTYPE <- "Bounce"
        event.df[which(event.df$burst == ii),]$cross <- 0
      }
      else {
        event.df[which(event.df$burst == ii),]$eventTYPE <- "Quick Cross"
        event.df[which(event.df$burst == ii),]$cross <- int.num
      }
    }
    
    #next for longer encounter
    else {
      # all points in the burst into one movement trajectory, calculate intersections with fence lines
      mov.seg.i <- encounter.df.i[which(encounter.df.i$burst == ii),]
      seg.line.i <- Lines(Line(cbind(mov.seg.i$coords.x1, mov.seg.i$coords.x2)), 
                          ID = mov.seg.i$date[1])
      seg.sp.i <- SpatialLines(list(seg.line.i), proj4string = CRS(target.crs))
      int.num <- length(gIntersection(seg.sp.i, fence.sp))
      
      if (difftime(end.time, start.time, units = "hours") > p * interval) {
        event.df[which(event.df$burst == ii),]$eventTYPE <- "Trapped"
        event.df[which(event.df$burst == ii),]$cross <- int.num
      }
      else {
        # calculating straightness of the encounter event
        straightness <- strtns(mov.seg.i)
        event.df[which(event.df$burst == ii),]$straightness <-
          straightness
        event.df[which(event.df$burst == ii),]$eventTYPE <-
          "TBD" # these will be further classified in the next loop
        event.df[which(event.df$burst == ii),]$cross <- int.num
      }
    }
  }
  event.df #this will be attached to the event.df in the "for each" loop
}

### close progress bar
close(pb)

### stop cluster
stopCluster(cl)

## ---- classification step 3: classify back-n-forth and trace ---- ####
### back-n-forth and trace are based on comparing average straightness around the encounter event

for (i in 1:nrow(event.df)) {
  if (event.df[i, ]$eventTYPE == "TBD") {
    Animal.str.i <-
      animal.stn.df[(
        animal.stn.df$AnimalID == event.df[i, ]$AnimalID &
          animal.stn.df$window.size == event.df[i, ]$duration / interval
      ), ]
    date.i <- event.df[i, ]$burstID
    straightness.i <- event.df[i, ]$straightness
    Animal.str.i$diff.time <-
      as.numeric(difftime(Animal.str.i$Date, date.i, units = "days"))
    Animal.str.i <- Animal.str.i[!is.na(Animal.str.i$diff.time), ]
    #need to make sure there are enough data to calculate average monthly straightness before and after the encounter event
    if ((sum(Animal.str.i$diff.time[!is.na(Animal.str.i$diff.time)] <= (-half.window)) == 0)
        |
        (sum(Animal.str.i$diff.time[!is.na(Animal.str.i$diff.time)] >=  (half.window)) == 0)) {
      event.df[i, ]$eventTYPE = "unknown"
    }
    else {
      ##the row that is the closest to starting date of the encounter event
      #row.n <- which(Animal.str.i$diff.time == min.nonneg(Animal.str.i$diff.time))[1]
      #the row half.window days before and half.window days after the starting of the encountering event
      row.max <-
        which((Animal.str.i$diff.time - half.window) == min.nonneg (Animal.str.i$diff.time - half.window)
        )
      row.min <-
        which((Animal.str.i$diff.time + half.window) == min.nonneg ((
          Animal.str.i$diff.time + half.window
        )))
      Animal.str.i.subset <- Animal.str.i[row.min:row.max, ]
      upper <-
        mean(Animal.str.i.subset$Straightness) + sd(Animal.str.i.subset$Straightness)
      lower <-
        mean(Animal.str.i.subset$Straightness) - sd(Animal.str.i.subset$Straightness)
      if (event.df[i, ]$straightness < lower) {
        if (event.df[i, ]$cross  < max.cross) {
        event.df[i, ]$eventTYPE <- "Back-n-forth"
        }
        else {
          event.df[i, ]$eventTYPE <- "unknown"
        }
      } else if (event.df[i, ]$straightness > upper) {
        if (event.df[i, ]$cross  < max.cross) {
          event.df[i, ]$eventTYPE <- "Trace"
        }
        else {
          event.df[i, ]$eventTYPE <- "unknown"
        }
      } else {
        event.df[i, ]$eventTYPE <- "Average Movement"
      }
    }
  }
}

# Save output -------------------------------------------------------------

write.csv(event.df, paste0(ouput_classification_events_file_path, "_I", interval,
                           "_FB", FB.dist, "_B",b.hours, "_P", p.hours, ".csv"))
