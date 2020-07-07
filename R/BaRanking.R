

range01 <- function(x){(x-min(x))/(max(x)-min(x))}

# default index_fun is  to calculate the index we used in the manuscript
BaRanking <- function(classification, barrier, d, min_total_enc = 0, 
                      index_fun = expression(((Bounce + Back_n_forth + Trace + Trapped)/total_enc)*unique_ind), 
                      show_plot = T) {
 
  barrier_sf <- as(barrier, "sf")
  
  # read as sf object
  classification_sf <- st_as_sf(classification, coords = c("easting", "northing"), crs = st_crs(barrier))
  
  # spatial join by the fence buffer distance used in BaBA
  barrier_sf_joined <- st_join(barrier_sf, classification_sf, join = st_is_within_distance, dist = d)
  
  # calculate # of each encounter envent types
  by_type <- as.data.frame(barrier_sf_joined)[!is.na(barrier_sf_joined$AnimalID),]  %>% group_by(FID_Fence0, eventTYPE, .drop = F) %>% 
    summarise(count = n())
  
  by_type_df <- as.data.frame(by_type) %>% select(FID_Fence0, eventTYPE, count) %>% pivot_wider(names_from = eventTYPE, values_from = count) %>% replace(is.na(.), 0)
  
  # find out whether all behavioral types are listed and added the column into the dataframe and fill with 0
  if (length(names(by_type_df)) < 8) {
    ba.names <- names(by_type_df)[2:length(names(by_type_df))]
    ba.all <- c("Quick_Cross", "Bounce", "Average_Movement", "Back_n_forth", "Trace", "Trapped", "unknown") # a full list of names
    ba.miss <- ba.all[!(ba.all %in% ba.names)]
    n.miss <- length(ba.miss)
    add.df <- as.data.frame(matrix(0, nrow = nrow(by_type_df), ncol = n.miss))
    colnames(add.df) <- ba.miss
    by_type_df <- cbind(by_type_df, add.df)
  }
  
  # calculate # of unique invidiuals interacted with each fence
  by_ID <- as.data.frame(barrier_sf_joined)[!is.na(barrier_sf_joined$AnimalID),]  %>% 
    group_by(FID_Fence0, .drop = F) %>%
    summarise(unique_ind = length(unique(AnimalID))) 
  
  # combine the two dataframes and calculate total encounters and total bad encounters
  barrier_sf_joined <- by_type_df %>% 
    left_join(by_ID, by = "FID_Fence0") %>% 
    replace(is.na(.), 0) %>%
    mutate(total_enc = Bounce + Quick_Cross + Average_Movement + Back_n_forth + Trace + unknown + Trapped) 
  
  # calculate impermeability index based on user set expression
  barrier_sf_joined <- barrier_sf_joined %>% filter(total_enc >= min_total_enc) %>% mutate(index = range01(eval(index_fun)))
  
  # put back into spatial
  barrier_sf_joined <- merge(barrier_sf, barrier_sf_joined, by = "FID_Fence0", all.x = T)
  
  if(show_plot) {
    plot(st_geometry(barrier_sf), col = "grey")
    plot(barrier_sf_joined['index'], add = T)
  }
  
  return(barrier_sf_joined)
  
}
