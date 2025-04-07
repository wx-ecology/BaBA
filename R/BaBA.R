# # # for debugging
# animal = pronghorn
# barrier = fences
# d = 110
# interval = NULL
# b_time = 4
# p_time = 36
# w = 168
# units = "hours"
# max_cross = 0
# tolerance = 0
# sd_multiplier = 1
# round_fixes = T
# exclude_buffer = F
# export_images = T


BaBA <-
  function(animal, barrier, d, interval = NULL, b_time = 4, p_time = 36, w = 168,
           tolerance = 0, units = "hours", max_cross = 0,  sd_multiplier = 1,
           exclude_buffer = F, round_fixes = F, export_images = F,
           img_path = "event_imgs", img_suffix = NULL) {
    
    # initial checks ----------------------------------------------------------
    if(export_images) {
      if(!dir.exists(img_path)) dir.create(img_path)
    }
    
    ## prepare parameters and check input
    if (class(animal)[1] != "sf") stop("animal needs to be an sf object")
    if (sf::st_geometry_type(animal)[1] != 'POINT') stop("animal needs to have a POINT geometry")
    if (class(barrier)[1] != "sf") stop("barrier needs to be an sf object")
    if (!(sf::st_geometry_type(barrier)[1] %in% c('MULTILINESTRING', 'LINESTRING'))) stop("barrier needs to have either a LINESTRING or MULTILINESTRING geometry")
    if (!"date" %in% names(animal)) stop("please rename the date column to 'date'")
    if (!"Animal.ID" %in% names(animal)) stop("please rename the individual ID column to 'Animal.ID'")
    if (!(inherits(animal$date, "POSIXct"))) stop("date needs to be 'POSIXct' format")
    if (sum(is.na(animal$date)) > 0) stop("please exclude rows where date is NA")
    
    if(round_fixes){
      interval_per_individual <- tapply(animal$date, animal$Animal.ID, function(x) names(which.max(table(round(as.numeric(diff(x), units = units),0)))))
    } else {
      interval_per_individual <- tapply(animal$date, animal$Animal.ID, function(x) names(which.max(table(as.numeric(diff(x), units = units)))))
    }
    if(is.null(interval)) { ## figure out interval (as the most frequent difference in timestamp) if not provided but give an error if not the same for all individuals
      if(all(interval_per_individual == interval_per_individual[1])) interval <- as.numeric(interval_per_individual[1]) else stop("Not all individuals have been sampled at the same frequency. Run individuals with different intervals seperately, or double-check whether your date column is cleaned.")
    } else {
      if (any(as.numeric(interval_per_individual) > interval, na.rm = T)) stop("BaBA interval needs to be no smaller than the actual data interval. Also double-check whether your date column is cleaned.") 
    }
    
    b <- b_time / interval
    if(b < 1) stop("interval needs to be set no bigger than b_time")
    if (round(b) != b) stop("b_time must be divisible by interval")
    p <- p_time / interval
    if (round(p) != p) stop("p_time must be divisible by interval")
    
    
    # classification step 1: generate encounter event data.frame -------------------------------------------------
    
    ## create point ID by individual
    animal <-
      animal %>% 
      dplyr::arrange(Animal.ID, date) %>% 
      dplyr::group_by(Animal.ID) %>% 
      dplyr::mutate(ptsID = 1:dplyr::n()) %>% 
      dplyr::ungroup()
    
    ## explicitly suppress constant geometry assumption warning by confirming attribute is constant throughout the geometry. See https://github.com/r-spatial/sf/issues/406 for details.
    sf::st_agr(animal) <- 'constant'   
    sf::st_agr(barrier) <- 'constant'
    
    ## create buffer around barrier
    print("locating encounter events...")
    barrier_buffer <- 
      barrier %>% 
      sf::st_buffer(dist = d, nQuadSegs = 5) %>%   ## Note that nQuadSegs is set to 5 as this was the default value for rgeos::gBuffer in previous versions of BaBA
      sf::st_union()
    
    ## extract points that fall inside the buffer
    encounter <- sf::st_intersection(animal, barrier_buffer)
    
    if (nrow(encounter) == 0) stop("no barrier encounter detected.")
    
    ## create unique burstIDs
    for(i in unique(encounter$Animal.ID)){
      if (nrow(encounter %>% dplyr::filter(Animal.ID == i)) == 0) {
        warning(paste0 ("Individual ", i, " has no locations overlapped with the barrier buffer and is eliminated from analysis." ))
        next()
      }
      encounter_i <-
        encounter %>% 
        dplyr::filter(Animal.ID == i) %>% 
        ## add time difference
        dplyr::mutate(
          ## time difference between all points in the buffer
          timediff = c(interval, as.numeric(diff(date), units = units)),
          ## remove the interval from that so when there is no missing point, timediff2 should be 0. If <0, replicated timestamp; if >0, missing timestamp
          timediff2 = round(timediff - interval, digits = 1))
      
      ## if any timediff2 is >= interval but <= tolerance, bring in the missing points from outside the buffer
      if(any(encounter_i$timediff2 >= interval & encounter_i$timediff2 <= tolerance, na.rm = T )) {
        idx_pts_of_interest <- which(encounter_i$timediff2 >= interval & encounter_i$timediff2 <= tolerance)
        
        for(pt in idx_pts_of_interest) {
          ## find out what pts to fetch
          ptsID_of_interest_B <- encounter_i$ptsID[pt]
          ptsID_of_interest_A <- encounter_i$ptsID[pt-1]
          
          ## fetch the points outside of the buffer and placehold timediff as NA and timediff2 as 0
          fetched_pt <- 
            animal %>% 
            dplyr::filter(Animal.ID == i & 
                            ptsID > ptsID_of_interest_A & 
                            ptsID < ptsID_of_interest_B)
          
          if (nrow(fetched_pt) == 0) {  ## if there's no point outside of the buffer between the timestamp that means there's missing data
            ## since the missing data is still within the tolerance, we consider timediff2=0 so the points before and after will be in the same event
            encounter_i$timediff2[pt] <- 0
            next() } 
          else {
            fetched_pt$timediff <- NA
            fetched_pt$timediff2 <- 0 
            ## reset timediff2 of pts_of_interests to 0
            encounter_i$timediff2[pt] <- 0 
            ## append fetched points to each other 
            if(pt == idx_pts_of_interest[1]) {fetched_pts <- fetched_pt} else if (exists("fetched_pts")) { fetched_pts <- rbind(fetched_pts, fetched_pt) } else {fetched_pts <- fetched_pt}
          }
        }
        
        ## append fetched pts
        encounter_i <- rbind(encounter_i, fetched_pts)
        ## recorder animal i's encounter event data.frame
        encounter_i <- encounter_i[order(encounter_i$ptsID), ]
        
        # reset fetched_pts for the next loop
        fetched_pts <- fetched_pts[0,]
      }
      
      ## do the cumulative sum of the new data.frame based on timediff2, using that as the updated unique burst ID (with animalID) 
      encounter_i$burstID <- paste(i, cumsum(encounter_i$timediff2), sep = "_")
      
      ## save into encounter_complete
      if(i == unique(encounter$Animal.ID[1])) encounter_complete <- encounter_i else encounter_complete <- rbind(encounter_complete, encounter_i)
    }
    
    encounter <- encounter_complete ## save back as encounter (encounter_complete is bigger as it includes extra points that are within tolerance)
    
    
    # classification step 2: classify events -------------------------------------------------
    
    print("classifying behaviors...") 
    ## open progress bar
    pb <- utils::txtProgressBar(style = 3)
    
    ## create empty object that will hold results
    event_df <- NULL
    
    ## run classification procedure for each encounter
    for(i in unique(encounter$burstID)) {
      ## update progressbar
      utils::setTxtProgressBar(pb, which(unique(encounter$burstID) == i)/length(unique(encounter$burstID)))
      
      ## get what we need from the encounter
      encounter_i <- encounter[encounter$burstID == i, ]
      animal_i <- animal[animal$Animal.ID == encounter_i$Animal.ID[1],]
      start_time <- encounter_i$date[1]
      end_time <- encounter_i$date[nrow(encounter_i)]
      duration <-  difftime (end_time, start_time, units = units)
      
      ## calculating straightness of the encounter event
      ## this will be used for median duration events but is output for reference for other events
      straightness_i <- strtns(encounter_i)
      
      
      ### classify short events (bounce and quick cross) ---------------------------------------------------
      
      ## if no more than b*interval, only spend small amount of time in this burst
      if (duration <= b_time) {
        pt.first <- encounter_i$ptsID[1] ## first point in the burst
        pt.last <- encounter_i$ptsID[nrow(encounter_i)]
        
        ## extract movement segment with one point before and one point after the segmentation
        mov_seg_i <- movement.segment.b(animal_i, pt.first, pt.last)
        
        ## count the number of crossings
        int.num <-
          mov_seg_i %>% 
          sf::st_intersection(barrier) %>% 
          sf::st_cast(to = 'MULTIPOINT') %>% 
          sf::st_coordinates() %>% 
          nrow()
        
        ## if no crossing is indicated and both before and after points were missing then we cannot tell if the animal crossed
        if (int.num == 0 & nrow(sf::st_coordinates(mov_seg_i)) != (nrow(encounter_i)+2)) {
          classification <- "unknown"
        } else {
          ## if there was not a crossing, classify as bounce, otherwise quick cross
          classification <- ifelse(int.num == 0, "Bounce", "Quick_Cross")
        }
      }
      
      ## dummy variable to ensure desired plotting and output
      tbd.plot <- 0
      
      if (duration > b_time) {
        
        ### classify trapped events -------------------------------------------------
        
        ## first calculate number of crossings (without looking at extra points like we did for short encounter)
        mov_seg_i <- 
          encounter_i %>% 
          dplyr::summarize(do_union = FALSE) %>% 
          sf::st_cast(to = 'LINESTRING')
        
        int.num <-
          mov_seg_i %>% 
          sf::st_intersection(barrier) %>% 
          sf::st_cast(to = 'MULTIPOINT') %>% 
          sf::st_coordinates() %>% 
          nrow()
        
        ## check if duration is smaller of bigger than p and classify accordingly
        if(duration > p_time) {
          classification <- "Trapped"
        } else {
          classification <- "TBD" ## these will be further classified in the next loop
        }
        
        
        ### classify back-n-forth, trace, and average movement -----------------------------------------
        
        ## process the "TBD" event types as back-n-forth, trace, or average movement
        ## back-n-forth and trace are based on comparing average straightness around the encounter event
        if(classification == 'TBD'){
          
          tbd.plot <- 1
          
          ## remove points that are inside the buffer if user said so
          if (exclude_buffer) {
            animal_i <- animal_i[!animal_i$ptsID %in% encounter$ptsID[encounter$Animal.ID == animal_i$Animal.ID[1]], ]
          }
          
          ## keep only data w/2 units before and w/2 after event
          animal_i <- animal_i[animal_i$date >= start_time - as.difftime(w/2, units = units) & animal_i$date <= end_time +  as.difftime(w/2, units = units), ]
          
          ## identify continuous sections in the remaining movement data
          animal_i$continuousID <- cumsum(abs(c(interval, round(as.numeric(diff(animal_i$date), units = units), digits = 1) - interval))) # abs() is to accommodate potential data points with smaller time intervals
          ## for each continuous sections, calculate straightness of all movements lasting the duration of our event (moving window of the size of the encounter)
          straightnesses_i <- NULL
          for(ii in unique(animal_i$continuousID)) {
            animal_ii <- animal_i[animal_i$continuousID == ii, ]
            
            ## duration of period
            duration_ii <- difftime(animal_ii$date[nrow(animal_ii)], animal_ii$date[1], units = units)
            
            ## calculate straightness only if at least as long as encounter event
            if(duration_ii >= duration) {
              for(iii in 1:(which(animal_ii$date > (animal_ii$date[nrow(animal_ii)] - as.difftime(duration, units = units)))[1] -1)) {
                mov_seg <- animal_ii[iii:(iii + duration/interval), ]
                straightnesses_i <- c(straightnesses_i, strtns(mov_seg))
              }
            }
          }
          
          ## make sure there are enough data to calculate average straightness before and after the encounter event
          ## (w/interval + 1) is the total possible segments if all data are present. 
          ## We define "enough" as at least 1/4 of the total possible segments are present to calculate average straightness.
          if (length(straightnesses_i) >= (w/interval + 1)/4) {
            ## minimum max number possible/2 to calculate sd
            upper <- mean(straightnesses_i) + sd_multiplier * stats::sd(straightnesses_i)
            lower <- mean(straightnesses_i) - sd_multiplier * stats::sd(straightnesses_i)
            if(straightness_i < lower) classification <- ifelse(int.num <= max_cross, "Back_n_forth", "unknown")
            if (straightness_i > upper) classification <- ifelse(int.num <= max_cross, "Trace", "unknown")
            if(straightness_i >= lower & straightness_i <= upper) classification <- "Average_Movement"
          } else {
            classification <- "unknown"
            if(is.null(straightnesses_i)) {straightnesses_i <- NA} ## add this to avoid warning message when plotting
          }
        }
      }
      
      
      ### Consolidate outputs -----------------------------------------------------
      
      ## plot the encounters to check later, if desired
      if (export_images) {
        grDevices::png(paste0(img_path, "/", classification, "_", i, "_", img_suffix, ".png"), width = 6, height = 6, units = "in", res = 90)
        if(tbd.plot == 0){
          bbox <- sf::st_bbox(mov_seg_i)
          expanded_bbox <- bbox
          expanded_bbox["xmin"] <- bbox["xmin"] - d*1.1
          expanded_bbox["xmax"] <- bbox["xmax"] + d*1.1
          expanded_bbox["ymin"] <- bbox["ymin"] - d*1.1
          expanded_bbox["ymax"] <- bbox["ymax"] + d*1.1
          
          plot(sf::st_as_sfc(expanded_bbox), col = "white", border = "white", main = classification, sub = paste("cross =", int.num, ", duration =", round(duration, 0), units))
          plot(mov_seg_i,add = T)
          plot(barrier_buffer, border = scales::alpha("red", 0.5), lty = "dashed", add = T)
          plot(sf::st_geometry(barrier), col = 'red', lwd = 2, add = TRUE)
          plot(sf::st_geometry(encounter_i), pch = 20, col = "cyan3", type = "o", lwd = 2, add = TRUE)
        } else{
          A <-
            animal_i %>% 
            dplyr::group_by(continuousID) %>% 
            ## Check sample size per group
            dplyr::mutate(n = dplyr::n()) %>%
            ## Exclude groups with only a single point to avoid errors
            dplyr::filter(n > 1) %>% 
            ## Convert to a line
            dplyr::summarize(do_union = FALSE) %>% 
            sf::st_cast(to = 'LINESTRING')
          plot(sf::st_geometry(A), main = classification,
               sub = paste0("cross = ", int.num, ", duration =", round(duration, 0), ", stri =", round(straightness_i, 2), ", str_mean = ",  round(mean(straightnesses_i), 2), ", str_sd = ",  round(stats::sd(straightnesses_i), 2)))
          plot(sf::st_geometry(barrier_buffer), border = scales::alpha("red", 0.5), lty = "dashed", add = TRUE)
          plot(sf::st_geometry(barrier), col = "red", lwd = 2, add = TRUE)
          plot(sf::st_geometry(encounter_i), pch = 20, col = "cyan3", type = "o", lwd = 2, add = TRUE)
        }
        grDevices::dev.off()
      }
      
      ## combine output
      event_df <- rbind(event_df, data.frame(
        AnimalID = encounter_i$Animal.ID[1],
        burstID = i,
        easting = sf::st_coordinates(encounter_i)[1, 1],
        northing = sf::st_coordinates(encounter_i)[1, 2],
        start_time,
        end_time,
        duration,
        cross = int.num,
        str_i = straightness_i,
        str_mean = ifelse(tbd.plot == 0, NA, mean(straightnesses_i)),
        str_sd = ifelse(tbd.plot == 0, NA, stats::sd(straightnesses_i)),
        eventTYPE = classification,
        stringsAsFactors = F
      ))
    }
    
    ## close progress bar
    close(pb)
    
    
    # finalize data -----------------------------------------------------------
    
    print("creating dataframe...")
    ## clean the encounter data
    encounter_final <- 
      encounter %>% 
      dplyr::filter(!duplicated(burstID)) %>% 
      dplyr::left_join(event_df %>% 
                         dplyr::select(burstID, eventTYPE),
                       by = 'burstID') %>% 
      dplyr::select(Animal.ID, burstID, date, eventTYPE)
    
    ## return output as a list
    return(list(encounters = encounter_final,
                classification = event_df))
  }


## increase movement segment by one points before and one point after the focused encounter
movement.segment.b <- function(animal, pt1, pt2) {
  pts_tmp <- animal[animal$ptsID >= pt1 - 1 & animal$ptsID <= pt2 + 1, ]
  pts_comb <- dplyr::summarize(pts_tmp, do_union = FALSE)
  segments_out <- sf::st_cast(pts_comb, to = 'LINESTRING')
  return(segments_out)
}

## helper function on for calculating Euclidean distance
calc_dist <- function(x.start, x.end, y.start, y.end){
  sqrt((x.end - x.start)^2 + (y.end - y.start)^2)
}

## calculate straightness of movement segment
strtns <- function(mov_seg) {
  
  locs_tmp <- sf::st_coordinates(mov_seg)
  
  if (sum(duplicated(mov_seg$date)) > 0 ) {
    straightness <- NA
    # warning("There are duplicated timestamps")
  } else if(nrow(locs_tmp) == 1){
    straightness <- NA
  } else {
    ## calculate Euclidean distance between the first and last point
    euc_dist <- 
      as.numeric(
        calc_dist(x.start = locs_tmp[1,1], x.end = locs_tmp[nrow(locs_tmp),1],
                  y.start = locs_tmp[1,2], y.end = locs_tmp[nrow(locs_tmp),2]))
    
    ## calculate path distance as the sum of all step lengths
    mov_seg$dist <- NA
    for(j in 2:nrow(mov_seg)){
      mov_seg$dist[j] <- calc_dist(x.start = locs_tmp[j - 1, 1],
                                   x.end = locs_tmp[j, 1],
                                   y.start = locs_tmp[j - 1, 2],
                                   y.end = locs_tmp[j, 2])
    }
    path_dist <- sum(mov_seg$dist, na.rm = TRUE)
    
    ## calculate straightness as the ratio of Euclidean to path distance.
    ## straightness ranges from 0 to 1 with values closer to 0 being more
    ## sinuous.
    straightness <- euc_dist/path_dist
  }
  
  return(straightness)
}
