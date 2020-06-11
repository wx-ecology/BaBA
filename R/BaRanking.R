

range01 <- function(x){(x-min(x))/(max(x)-min(x))}

# default index_fun is  to calculate the index we used in Scenario 2
BaRanking <- function(classification = results_pron$classification, barrier = fences, d = 110, min_total_enc = 10, index_fun = expression((alt_enc/total_enc)*unique_ind), show_plot = T) {
 
  barrier_sf <- as(barrier, "sf")
  
  # recode eventTYPE for later column calculation
  eventTYPE <- as.character(classification$eventTYPE)
  eventTYPE <- factor(recode(eventTYPE, "Quick Cross" = "Quick_Cross", "Average Movement" = "Average_Movement", "Back-n-forth" = "Back_n_forth"), levels = c("Quick_Cross", "Bounce", "unknown", "Average_Movement", "Trace", "Back_n_forth", "Trapped"))
  classification$eventTYPE <- eventTYPE
  
  # read as sf object
  classification_sf <- st_as_sf(classification, coords = c("easting", "northing"), crs = st_crs(barrier))
  
  # spatial join by the fence buffer distance used in BaBA
  barrier_sf_joined <- st_join(barrier_sf, classification_sf, join = st_is_within_distance, dist = d)
  
  # calculate # of each encounter envent types
  by_type <- as.data.frame(barrier_sf_joined)[!is.na(barrier_sf_joined$AnimalID),]  %>% group_by(FID_Fence0, eventTYPE, .drop = F) %>% 
    summarise(count = n())
  
  by_type_df <- as.data.frame(by_type) %>% select(FID_Fence0, eventTYPE, count) %>% pivot_wider(names_from = eventTYPE, values_from = count)
  
  # calculate # of unique invidiuals interacted with each fence
  by_ID <- as.data.frame(barrier_sf_joined)[!is.na(barrier_sf_joined$AnimalID),]  %>% 
    group_by(FID_Fence0, .drop = F) %>%
    summarise(unique_ind = length(unique(AnimalID))) 
  
  # combine the two dataframes and calculate total encounters and total bad encounters
  barrier_sf_joined <- by_ID %>% 
    left_join(by_type_df, by = "FID_Fence0") %>% 
    replace(is.na(.), 0) %>%
    mutate(total_enc = Bounce + Quick_Cross + Average_Movement + Back_n_forth + Trace + unknown + Trapped, 
           alt_enc = Bounce + Back_n_forth + Trace + Trapped) 
  
  # change back the col names for display
  # names(barrier_sf_joined) <- recode(names(barrier_sf_joined),  "Quick_Cross" = "Quick Cross", "Average_Movement" = "Average Movement" , "Back_n_Forth" = "Back-n-forth")
  
  
 

  barrier_sf_joined <- barrier_sf_joined %>% filter(total_enc >= min_total_enc) %>% mutate(index = range01(eval(index_fun)))
  
  # put backk into spatial
  barrier_sf_joined <- merge(barrier_sf, barrier_sf_joined, by = "FID_Fence0", all.x = T)
  barrier_sf_joined$index[is.na(barrier_sf_joined$index)] <- 0
  
  
  if(show_plot) plot(barrier_sf_joined['index'])
  
  
  return(barrier_sf_joined)
  
  
}
