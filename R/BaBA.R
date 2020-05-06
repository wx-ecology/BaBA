BaBA <- function(animal, barrier, d = 50, interval = NULL, b_hours = 4, p_hours = 36, w = 7, max_cross = 4, tolerance = 0, sd_multiplier = 1, exclude_buffer = F, export_images = T,img_path = "event_imgs", img_suffix = NULL) {

  if(export_images) {
    if(!dir.exists(img_path)) dir.create(img_path)
  }
  # prepare parameters ####
  if(is.null(interval)) { # figure out interval (as the most frequent difference in timestamp) if not provided but give an error if not the same for all individuals
    interval_per_individual <- tapply(animal$date, animal$Location.ID, function(x) names(which.max(table(diff(x)))))
    
    if(all(interval_per_individual == interval_per_individual[1])) interval <- as.numeric(interval_per_individual[1]) else stop("You did not provide a time interval and not all individuals have been sampled at the same frequency.")
  }
  
  b <- b_hours / interval
  p <- p_hours / interval

  # create point ID by individual ----

  animal$ptsID <- NA

  for (i in unique(animal$Location.ID)) {
    mov.seg.i <- animal[animal$Location.ID == i, ]
    animal@data$ptsID[animal$Location.ID == i] <-
      seq(nrow(mov.seg.i))
  }
  
  # create buffer around barrier ----
  barrier_buffer <- raster::buffer(barrier, width = d)

  # ---- classification step 1: generate encountering event dataframe ---- ####

    ## extract points that fall inside the buffer ----
    encounter <- raster::intersect(animal, barrier_buffer)

    ## create a burstID ----
 
    for(i in unique(encounter$Location.ID)){

      encounter_i <- encounter[encounter$Location.ID == i,]
   ## first get time difference between all points in the buffer
      encounter_i$timediff <- c(interval, diff(encounter_i$date ))
    ## then remove the interval from that
    encounter_i$timediff2 <-  encounter_i$timediff  - interval
    
    
    ## then, if any timediff2 is >0 but <= tolerance, we need to bring in the missing points from outside the buffer.
    if(any(encounter_i$timediff2 > 0 & encounter_i$timediff2 <= tolerance)) {
      
      idx_pts_of_interest <- which(encounter_i$timediff2 > 0 & encounter_i$timediff2 <= tolerance)
      
      for(pt in idx_pts_of_interest) {
        # find out what pts to fetch
        ptsID_of_interest_B <- encounter_i$ptsID[pt]
        ptsID_of_interest_A <- encounter_i$ptsID[pt-1]
        
        # fetch the points and add timediff as NA and timediff2 as 0
        fetched_pt <- animal[animal$Location.ID == i & animal$ptsID > ptsID_of_interest_A & animal$ptsID < ptsID_of_interest_B, ]
        fetched_pt$timediff <- NA
        fetched_pt$timediff2 <- 0
        
        # replace timediff2 of pts_of_interests to 0
        encounter_i$timediff2[pt] <- 0
        
        # append fetched points to each other 
        if(pt == idx_pts_of_interest[1]) fetched_pts <- fetched_pt else fetched_pts <- rbind(fetched_pts, fetched_pt)
        
      }
     
     # append pts and reorder
      encounter_i <- rbind(encounter_i, fetched_pts)
      encounter_i <- encounter_i[order(encounter_i$ptsID), ]
      
    }
    
  
    ## then do the cum sum of that, and that is the burst ID (with animalID) (since now any same number is from the same burst)
    
    encounter_i$burstID <- paste(i, cumsum( encounter_i$timediff2), sep = "_")

    # save into encounter_complete ####
    if(i == unique(encounter$Location.ID[1])) encounter_complete <- encounter_i else encounter_complete <- rbind(encounter_complete, encounter_i)



  }

  encounter <- encounter_complete # save back as encounter (encoutner_complete is bigger as it includes extra points that  are within tolerance)
   
  
    
  
  ## ---- classification step 2: classify bounce, quick cross, and trap based on duration---- ####
  ### open progress bar ----
  pb <- txtProgressBar( style = 3)

  ### create empty object that will hold results ----
  event_df <- NULL

  ## run classification procedure for each encounter ####
  for(i in unique(encounter$burstID)) {

    # update progressbar
    setTxtProgressBar(pb, which(unique(encounter$burstID) == i)/length(unique(encounter$burstID)))

    # get what we need from the encounter ####
    encounter_i <- encounter[encounter$burstID == i, ]
    animal_i <- animal[animal$Location.ID == encounter_i$Location.ID[1],]

    start_time <- encounter_i$date[1]
    end_time <- encounter_i$date[nrow(encounter_i)]
    duration <-  difftime (end_time, start_time, units = "hours")
    # calculating straightness of the encounter event ###
    ## this will be used for median duration event but is output for reference for other event ####
    straightness_i <- strtns(encounter_i)

    # classify short encounters (bounce and quick cross) ####
    #no more than b*interval H, only spend small amount of time in this burst
    if (duration <= b * interval) {
      pt.first <- encounter_i$ptsID[1]#first point in the burst
      pt.last <- encounter_i$ptsID[nrow(encounter_i)]

      # extract movement segment with one point before and one point after the segmentation ####
      mov_seg_i <- movement.segment.b(animal_i, pt.first, pt.last)

      # count the number of crossing ####
      int.num <- length(rgeos::gIntersection(mov_seg_i, barrier))

      # if no crossing and we didn't have both points (before and after), then we can't tell if it crossed
      if (int.num == 0 & nrow(coordinates(mov_seg_i)[[1]][[1]]) != (nrow(encounter_i)+2)) {
        # means that no points were before or after the encounter and we can't tell if the animal crossed
        classification <- "unknown"

      } else {
        classification <- ifelse(int.num == 0, "Bounce", "Quick Cross")

      }

# plot these examples to check later
      if(export_images) {
        png(paste0(img_path, "/", classification, "_", i, "_", img_suffix, ".png"), width = 6, height = 6, units = "in", res = 90)
        plot(mov_seg_i, main = classification, sub = paste("cross =", int.num, "- duration =", duration))
        plot(barrier_buffer, border = scales::alpha("red", 0.5), add = T)
        lines(barrier, col = "red")
        points(encounter_i, pch = 16, col = "blue")
        points(animal_i[abs(difftime(animal_i$date, start_time, units = "days")) <= 0.5,], type = "o") # one day of data around encounter
        dev.off()
      }

    }

    # classify longer encounters (only trap for now) ########
    if (duration > b * interval) {

      ## first calculate number of crossings (without looking at extra points like we did for short encounter)
      mov_seg_i <- SpatialLines(list(Lines(Line(coordinates(encounter_i)),                            ID = encounter_i$date[1])), proj4string = CRS(proj4string(animal)))
      int.num <- length(gIntersection(mov_seg_i, barrier))
 ## then check if duration is smaller of bigger than p and classify accordingly
      if(duration > p * interval) {

        classification <- "Trapped"

      } else {

        classification <- "TBD" # these will be further classified in the next loop

      }


      # plot these examples to check later
      if(export_images & !classification %in% "TBD") {
        png(paste0(img_path, "/", classification, "_", i, "_", img_suffix, ".png"), width = 6, height = 6, units = "in", res = 90)
        plot(mov_seg_i, main = classification, sub = paste("cross =", int.num, "- duration =", duration))
        plot(barrier_buffer, border = scales::alpha("red", 0.5), add = T)
        lines(barrier, col = "red")
        points(encounter_i, pch = 16, col = "blue")
        points(animal_i[abs(difftime(animal_i$date, start_time, units = "days")) <= 0.5,], type = "o") # one day of data around encounter
        dev.off()
      }
    }

    # save output ####

    event_df <- rbind(event_df, data.frame(
      AnimalID = encounter_i$Location.ID[1],
      burstID = i,
      easting = coordinates(encounter_i)[1, 1],
      northing = coordinates(encounter_i)[1, 2],
      start_time,
      end_time,
      duration,
      cross = int.num,
      straightness = straightness_i,
      eventTYPE = classification,
      stringsAsFactors = F

    ))
  }


  ### close progress bar ----
  close(pb)

  ## ---- classification step 3: classify back-n-forth and trace ---- ####
  ### back-n-forth and trace are based on comparing average straightness around the encounter event

  for (i in 1:nrow(event_df)) {
    if (event_df[i, ]$eventTYPE == "TBD") {
      event_i <- event_df[i, ]
      duration_i <- event_i$duration
      straightness_i <- event_i$straightness

      # get movement data of the individual of interest
      animal_i <- animal[animal$Location.ID == event_i$AnimalID, ]
      encounter_i <- encounter[encounter$burstID %in% event_i$burstID,]

      # remove points that are inside the buffer if used said so, if not, at least remove the points of the event
      if(exclude_buffer) {
        animal_i <- animal_i[!animal_i$ptsID %in% encounter$ptsID[encounter$Location.ID == event_i$AnimalID], ]
        } 
      # else {
      #   animal_i <- animal_i[!animal_i$ptsID %in% encounter$ptsID[encounter$Location.ID == event_i$AnimalID] & encounter$burstID %in% event_i$burstID, ]
      # }

      # keep only data X days before and X after event
      animal_i <- animal_i[animal_i$date >= event_i$start_time - as.difftime(w/2, units = "days") & animal_i$date <= event_i$end_time +  as.difftime(w/2, units = "days"), ]

      # identify continuous sections in the remaining movement data
      animal_i$continuousID <- cumsum(c(interval, diff(animal_i$date)) - interval)

      # for each continuous sections, calculate straigness of all movements lasting the duration of our event (moving window of the size of the encounter)

      straightnesses_i <- NULL
      for(ii in unique(animal_i$continuousID)) {
        animal_ii <- animal_i[animal_i$continuousID == ii, ]
        #duration of period
        duration_ii <- difftime(animal_ii$date[nrow(animal_ii)], animal_ii$date[1])

        # calculate straigness only if at least as long as encounter event
        if(duration_ii >= duration_i) {
          for(iii in c(1: (which(animal_ii$date > animal_ii$date[nrow(animal_ii)] - as.difftime(duration_i, units = "hours"))[1] -1))) {
            mov_seg <- animal_ii[iii:(iii+duration_i/interval), ]

            straightnesses_i <- c(straightnesses_i, strtns(mov_seg))

          }

        }
      }


      #need to make sure there are enough data to calculate average monthly straightness before and after the encounter event
      if(length(straightnesses_i) >= ((w*24/interval+1)/2)) {
        # minimum max number possible/2 to calculate sd
        upper <- mean(straightnesses_i) + sd_multiplier * sd(straightnesses_i)
        lower <- mean(straightnesses_i) - sd_multiplier * sd(straightnesses_i)
        if(straightness_i < lower) event_df[i, ]$eventTYPE <- ifelse(event_i$cross < max_cross, "Back-n-forth", "unknown")
        if (straightness_i > upper) event_df[i, ]$eventTYPE <- ifelse(event_i$cross < max_cross, "Trace", "unknown")
        if(straightness_i >= lower & event_i$straightness <= upper) event_df[i, ]$eventTYPE <- "Average Movement"
      }
      else {
        event_df[i, ]$eventTYPE = "unknown"
      }


      # plot to check later ####
      if(export_images) {
        png(paste0(img_path, "/",  event_df[i, ]$eventTYPE, "_", event_i$burstID, "_", img_suffix, ".png"), width = 6, height = 6, res = 96, units  = "in")
      A = by(as.data.frame(coordinates(animal_i)), animal_i$continuousID, Line, simplify = T)

      A = SpatialLines(mapply(sp::Lines, A, ID = names(A), SIMPLIFY = F))

        plot(A,
           main = event_df[i, ]$eventTYPE, sub = paste("cross = ", event_df[i, ]$cross, ", duration =", event_df[i, ]$duration, ", stri =", round(straightness_i, 2), ", str_mean = ",  round(mean(straightnesses_i), 2), ", str_sd = ",  round(sd(straightnesses_i), 2)))
        plot(barrier_buffer, border = scales::alpha("red", 0.5), add = T)
      lines(barrier, col = "red")
      points(encounter_i, pch = 16, col = "blue")
      points(animal_i[abs(difftime(animal_i$date, event_i$start_time, units = "days")) <= 0.5,], type = "o") # one day of data around encounter
      dev.off()
}
      }
  }
  ## return output as a lits ####
  return(list(encounters = encounter,
              classification = event_df))
  }


#increase movement segment by one points before and one point after the focused encounter ####
movement.segment.b <- function(animal, pt1, pt2) {
  segments <- animal[animal$ptsID >= pt1 - 1 &
                       animal$ptsID <= pt2 + 1, ]
  seg.line <- Lines(Line(coordinates(segments)),
                    ID = segments$date[1])

  segments.sp <- SpatialLines(list(seg.line), proj4string = CRS(proj4string(animal)))

  return(segments.sp)
}


# calculate straigness of movement segment ####
strtns <- function(mov_seg) {



  # calculate trajectory
  traj <- adehabitatLT::as.ltraj(xy = coordinates(mov_seg), date = mov_seg$date, id = as.character(mov_seg$Location.ID))

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



