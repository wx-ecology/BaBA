################ Barrier Behavior Analysis (BaBA): classification  ##########################
# Time: 11272019
# Description: This is the classification step of BaBA 
## Input ##
# Fence shp, movement points (with the same time intervals since the behavior type assumption 
# Straighness table 

## output ##
## 1. encounter event table; 2. step 1 classification table (bounce/trapped/quick cross); 3. Final classification table

## data requirement ##
# Movement coordinates - projected coordinations (m). In this analysis we use 2-hour interval.
# # ideally, the movement data should not have missing point
# Fence polylines - as accurate as possible

# Update from v1: seperated step out of a single loop 

#############################################################################################
# ----- libraries -------
library(dplyr)
# spatial analysis
library(rgdal)
library(rgeos)
library(sp)
library(raster)
#trajectory analysis
library(adehabitatLT)
# for parallel multi-core calculation 
library(foreach)
library(doParallel)
library(doSNOW)

#############################
#########Parameters##########
#############################s
target.crs <- "+proj=utm +zone=12 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"
interval <- 2 #define time interval (in hours) of the movement data

## Fence buffer distance in meters 
# advised fence buffer distance is the 1st qu. of all points distance to fences
# there is a seperate script that can calculate distance from all points to fences (Dist2FecneAnalysis_Official)
FB.dist <- 50

## tolerance parameter. If "a" is point in the buffer, "b" is a buffer outside of the buffer
# x is the number of b that is allowed in between of a to allow the point series be considered as a continuous encounter event
# it is more useful for high temporal resolution data
tolerance <- 0

## parameters after this should be determined by local management requirements and characteristics of the target species.
# based on the data time interval and animal ecology, maximum encounter duration that you'd call it a "bounce" or "Quick Cross".
# aka, what is quick for you?
b.hours <- 4
b <- b.hours/interval
# minimum encounter duration in the burst that you'd call the animal might be trapped, thus no classification will be followed.
p.hours <- 36
p <- p.hours/interval

# When differenciating trace/back-n-forth from "normal movement", 
# we compare the straightness of the encounter event to the average straightness around the time the encounter event happens.
# how long would you set this window? Current default is 7 days. 
ave.window <- 7
half.window <- ave.window/2

# maximum crosses allowed for tracing behavior 
max.cross <- 4 
# 

#############################
#########Functions###########
#############################
# extract movement segment between the time of m and predetermined time lap b*interval.  
movement.segment.b <- function(pt1, pt2) {
  segments <- movement.df[which(movement.df$ptsID >= pt1-1 & movement.df$ptsID <= pt2 + 1),]
  seg.line <- Lines(Line(cbind(segments$coords.x1, segments$coords.x2)), ID = segments$date[1])
  segments.sp <- SpatialLines(list(seg.line), proj4string = CRS(target.crs))
  return(segments.sp)
}

# calculating straightness. Input a dataframe with Easting and Northing. 
strtns <- function(mov.seg) {
  pts <- cbind(mov.seg$Easting, mov.seg$Northing)
  pts.sp <- SpatialPointsDataFrame(pts, mov.seg, proj4string = CRS(target.crs))
  traj <- as.ltraj(xy =  pts, date = mov.seg$date, id = mov.seg$Location.ID)
  #moving distance from first pt to last pt in the burst
  traj.dist <- sqrt(
    as.numeric((traj[[1]]$x[1]-traj[[1]]$x[nrow(traj[[1]])]))*  as.numeric((traj[[1]]$x[1]-traj[[1]]$x[nrow(traj[[1]])])) +
      as.numeric((traj[[1]]$y[1]-traj[[1]]$y[nrow(traj[[1]])]))*as.numeric((traj[[1]]$y[1]-traj[[1]]$y[nrow(traj[[1]])])) 
  )
  #sum of all step lengths
  traj.lgth <- sum(traj[[1]]$dist, na.rm = TRUE)
  #straightness ranges from 0 to 1. More close to 0 more sinuous it is.
  straightness <- traj.dist/traj.lgth
  return(straightness)
}

# pick the minimum non-negative number
min.nonneg <- function(x) min(x[x >= 0])

#############################
######### Set-up ###########
#############################
setwd("C:\\Users\\wenjing.xu\\Google Drive\\RESEARCH\\Pronghorn\\Analysis\\FenceBehavior_Official")
# prepare spatial dataframe -----------------------------------
#read in fence data
fence.filename <- 'Fence_AOI_FINAL_Dec2019'
fence.sp <- readOGR(".", fence.filename)
#fence.sp <- spTransform(fence.sp,target.crs)
fence.buffer <- raster::buffer(fence.sp, width=FB.dist)

#read in movement data
#ideally, the movement data should not have missing point. This trial file does have missing points.
movement.df.all <- read.csv("Int2_MULE_Raw_Final.csv") 
movement.df.all$date <- as.POSIXct(strptime(as.character(movement.df.all$date),"%Y-%m-%d %H:%M")) #change the format based on the data
movement.df.all <- movement.df.all <- movement.df.all[(!is.na(movement.df.all$date))&(!is.na(movement.df.all$Easting)),]

# add point ID by individual
movement.df.all$ptsID <- numeric(nrow(movement.df.all))
for (i in unique(movement.df.all$Location.ID)) {
  mov.seg.i <- movement.df.all[movement.df.all$Location.ID==i,]
  movement.df.all[movement.df.all$Location.ID==i,]$ptsID <- seq(nrow(mov.seg.i))
}
# make movement pts into spatial point dataframe
xy <- cbind(movement.df.all$Easting, movement.df.all$Northing)
movement.sp.all <- SpatialPointsDataFrame (coords = xy, data = movement.df.all, proj4string = CRS(target.crs))

# read in straightness table 
animal.stn.df <- read.csv("I2_MULE_Straightness.csv")
################################################################################
# classification step 1: generate encountering event dataframe --------

# for parallel looping 
cores <- detectCores()
cl <- makeSOCKcluster(cores[1]-1) #to not overload your computer
#registerDoParallel(cl)
registerDoSNOW(cl)

pb <- txtProgressBar(max=100, style=3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress=progress)


encounter.df = foreach (i = unique(movement.df.all$Location.ID), 
                    .combine=rbind,
                    .options.snow=opts,
                    .packages=c('raster', 'sp', 'rgdal', 'rgeos', 'adehabitatLT', 'sp', 'dplyr')) %dopar% { 
                      
  # ---------------------  build dataframe, burst ID = first timestamp -----------------------------
  movement.sp <- movement.sp.all[movement.sp.all$Location.ID == i, ]
  movement.df <- as.data.frame(movement.sp)
  #extract points that fall inside of the buffer 
  encounter.sp <- raster::intersect(movement.sp, fence.buffer)
  #Creat a data frame with encounter event marked as burst ID
  encounter.df <- as.data.frame(encounter.sp) 
  encounter.df <- encounter.df[which(!is.na(encounter.df$coords.x1)),]
  encounter.df$date <- as.POSIXct(strptime(as.character(encounter.df$date),"%Y-%m-%d %H:%M"))
  
  encounter.df$burst <- vector(length = nrow(encounter.df))
  for (ii in 2:nrow(encounter.df)){
    # first group events based on time
    #if the time intervals between two points are within the set tolerance time, they are likely belog to the same encounter event, marked as "m"
    if (difftime(encounter.df$date[ii], encounter.df$date[ii-1],units = "hours") <= (tolerance+1)*interval) {
      encounter.df$burst[ii-1] <- "m"  #mark the cell for the burst
      if (ii == nrow(encounter.df)) {
        encounter.df$burst[ii] <- "m"
      }
    }
    else {
      #all points before this point that are marked as "m" are in the same burst
      encounter.df$burst[ii-1] <- "m"
      burst.list <- encounter.df[which(encounter.df$burst == "m"),]
      encounter.df[which(encounter.df$burst == "m"),]$burst <- format(burst.list$date[1], "%Y-%m-%d %H:%M") #name each burst using the starting time stamp
      if (ii == nrow(encounter.df)) {
        encounter.df$burst[ii] <- "m"
      }
    }
    if (ii == nrow(encounter.df)) {
      burst.list <- encounter.df[which(encounter.df$burst == "m"),]
      if (nrow(burst.list) > 0) {
        encounter.df[which(encounter.df$burst == "m"),]$burst <- format(burst.list$date[1], "%Y-%m-%d %H:%M")
      }
    }
  }
  # clean the data frame
  encounter.df <- encounter.df[which(!is.na(encounter.df$burst)),]  #all points that are in buffer in one dataframe
  encounter.df
}
write.csv(encounter.df, paste0("I2_PRON_FB", FB.dist, "_B4_P36_EncounterEvents.csv"))

close(pb)
#stop cluster
stopCluster(cl)

################################################################################
#classification step 2: classify bounce, quick cross, and trap ---------------

# for parallel looping 
cores <- detectCores()
cl <- makeSOCKcluster(cores[1]-1) #to not overload your computer
#registerDoParallel(cl)
registerDoSNOW(cl)

pb <- txtProgressBar(max=100, style=3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress=progress)

#start the loop
event.df = foreach (
  i = unique(encounter.df$Location.ID),
  .combine = rbind,
  .options.snow = opts,
  .packages = c('raster', 'sp', 'rgdal', 'rgeos', 'adehabitatLT', 'sp', 'dplyr')) %dopar% {
    
  movement.df <- as.data.frame(movement.sp.all[movement.sp.all$Location.ID == i,])
  encounter.df.i <- encounter.df[encounter.df$Location.ID == i, ]
  
  #------------------------  identify bounce, trapped, and cross, calculating duration and straightess  -------------------------------------
  # make a dataframe for this animal for layer be joined to the event.df later
  
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
  
  for (ii in unique(encounter.df.i$burst)) {
    burst.i <- encounter.df.i[encounter.df.i$burst == ii,]
    start.time <- burst.i[1, ]$date
    end.time <- burst.i[nrow(burst.i), ]$date
    event.df[which(event.df$burst == ii), ]$duration <- difftime (end.time, start.time, units = "hours")
    event.df[which(event.df$burst == ii), ]$easting <- burst.i$coords.x1[1]
    event.df[which(event.df$burst == ii), ]$northing <- burst.i$coords.x2[1]
    
    #first for short encounter
    if (difftime (end.time, start.time, units = "hours") <= b * interval) {  #no more than b*interval H, only spend small amount of time in this burst
      pt.first <- burst.i[1, ] #first point in the burst
      pt.last <- burst.i[nrow(burst.i), ]
      mov.seg.i <- movement.segment.b(pt.first$ptsID, pt.last$ptsID) #extract movement segment with one point before and one point after the segmentation
      #here is another measurement that prefer smaller time interval movement data. Big interval data between point distance can misrepresent real trajectory
      if (nrow(coordinates(mov.seg.i)[[1]][[1]]) <= b) {
        #which means this event is at the end of the individual's traj. Not enough points to tell catagories.
        event.df[which(event.df$burst == ii), ]$eventTYPE <- "unknown"
        next
      }
      int.num <- length(gIntersection(mov.seg.i, fence.sp))
      ############# Need to acknowledge possibility for mistakes here since the links between points are not real movement routes.
      if (int.num == 0) {
        event.df[which(event.df$burst == ii), ]$eventTYPE <- "Bounce"
        event.df[which(event.df$burst == ii), ]$cross <- 0
      }
      else {
        event.df[which(event.df$burst == ii), ]$eventTYPE <- "Quick Cross"
        event.df[which(event.df$burst == ii), ]$cross <- int.num
      }
    }
    
    #next for longer encounter
    else {
      # all points in the burst into one movement trajectory, calculate intersect with fence lines
      mov.seg.i <- encounter.df.i[which(encounter.df.i$burst == ii), ]
      seg.line.i <- Lines(Line(cbind(mov.seg.i$coords.x1, mov.seg.i$coords.x2)), ID = mov.seg.i$date[1])
      seg.sp.i <- SpatialLines(list(seg.line.i), proj4string = CRS(target.crs))
      int.num <- length(gIntersection(seg.sp.i, fence.sp))
      if (difftime(end.time, start.time, units = "hours") > p * interval) {
        event.df[which(event.df$burst == ii), ]$eventTYPE <- "Trapped" #sometimes could also be a home ranging behavior. Can be differenciate by time
        event.df[which(event.df$burst == ii), ]$cross <- int.num
      }
      else {
        # calculating straightness of the encounter event
        straightness <- strtns(mov.seg.i)
        event.df[which(event.df$burst == ii), ]$straightness <- straightness
        event.df[which(event.df$burst == ii), ]$eventTYPE <- "TBD"
        event.df[which(event.df$burst == ii), ]$cross <- int.num
      }
    }
  }
  event.df #this will be attached to the event.df in the "for each" loop
}
#event.df.temp <- event.df
#event.df <- event.df.temp
write.csv(event.df, paste0("I2_PRON_FB", FB.dist, "_B4_P36_Step1Cls.csv"))

close(pb)
#stop cluster
stopCluster(cl)


################################################################################
# classification step 3: classify back-n-forth and trace -----------------------
# based on comparing average straightness around the encounter event -----------

#animal.stn.df <- read.csv(paste0(getwd(), "I2_all_Straightness.csv"))
#animal.stn.df$Date <- as.POSIXct(strptime(as.character(animal.stn.df$Date),"%Y-%m-%d %H:%M"))

#event.df <- read.csv(paste0(getwd(), "I2_All_FB180_B4_P36_Step1Cls.csv"))
#event.df$Date <- as.POSIXct(strptime(as.character(event.df$burstID),"%Y-%m-%d %H:%M"))
#event.df$eventTYPE <- as.character(event.df$eventTYPE)
#event.df[is.na(event.df$eventTYPE),]$eventTYPE <- "unknown"

for (i in 1:nrow(event.df)) {
  if (event.df[i,]$eventTYPE == "TBD") {
    Animal.str.i <- animal.stn.df[(animal.stn.df$AnimalID == event.df[i,]$AnimalID & animal.stn.df$window.size == event.df[i,]$duration/interval),]
    date.i <- event.df[i,]$burstID
    straightness.i <- event.df[i,]$straightness
    Animal.str.i$diff.time <- as.numeric(difftime(Animal.str.i$Date, date.i, units = "days"))
    Animal.str.i <- Animal.str.i[!is.na(Animal.str.i$diff.time),]
    #need to make sure there are enough data to calculate average monthly straightness before and after the encounter event
    if ((sum(Animal.str.i$diff.time[!is.na(Animal.str.i$diff.time)] <= (-half.window)) == 0)
        | (sum(Animal.str.i$diff.time[!is.na(Animal.str.i$diff.time)] >=  (half.window)) == 0)) {
      event.df[i,]$eventTYPE = "unknown"
    }
    else {
      ##the row that is the closest to starting date of the encounter event
      #row.n <- which(Animal.str.i$diff.time == min.nonneg(Animal.str.i$diff.time))[1]
      #the row half.window days before and half.window days after the starting of the encountering event
      row.max <- which((Animal.str.i$diff.time - half.window) == min.nonneg (Animal.str.i$diff.time - half.window))
      row.min <- which((Animal.str.i$diff.time + half.window) == min.nonneg ((Animal.str.i$diff.time + half.window)))
      Animal.str.i.subset<- Animal.str.i[row.min:row.max,]
      upper <- mean(Animal.str.i.subset$Straightness) + sd(Animal.str.i.subset$Straightness) 
      lower <-  mean(Animal.str.i.subset$Straightness) - sd(Animal.str.i.subset$Straightness) 
      if (event.df[i,]$straightness < lower) {
        event.df[i,]$eventTYPE <- "Back-n-forth"
      } else if (event.df[i,]$straightness > upper) {
          if (event.df[i,]$cross  < max.cross) {
            event.df[i,]$eventTYPE <- "Trace"
          }
          else {
            event.df[i,]$eventTYPE <- "unknown"
          }
      } else {
        event.df[i,]$eventTYPE <- "Average Movement"
      }
    }
  }
}
event.df.1 <- event.df
write.csv(event.df.1, paste0("I2_PRON_FB", FB.dist, "_B4_P36_FinalCls.csv"))

