
range01 <- function(x){(x-min(x))/(max(x)-min(x))}

# default index_fun is  to calculate the index we used in the manuscript
BaRanking <- function(classification, barrier, d, Barrier_ID, min_total_enc = 0, 
                      index_fun = expression(((Bounce + Back_n_forth + Trace + Trapped)/total_enc)*unique_ind), 
                      show_plot = F) {
 
  barrier_sf <- as(barrier, "sf")
  Barrier_ID <- sym(Barrier_ID)
  
  # read as sf object
  classification_sf <- st_as_sf(classification, coords = c("easting", "northing"), crs = st_crs(barrier))
  
  # spatial join by the fence buffer distance used in BaBA
  barrier_sf_joined <- st_join(barrier_sf, classification_sf, join = st_is_within_distance, dist = d)
  
  # calculate # of each encounter envent types
  
  by_type <- as.data.frame(barrier_sf_joined)[!is.na(barrier_sf_joined$AnimalID),]  %>% 
    group_by(!!Barrier_ID, eventTYPE, .drop = F) %>%   
    summarise(count = n(), .groups = 'drop') # The bang-bang operator !! forces a single object. One common case for !! is to substitute an environment-variable (created with <-) with a data-variable (inside a data frame).
  
  by_type_df <- as.data.frame(by_type) %>% dplyr::select(!!Barrier_ID, eventTYPE, count) %>% pivot_wider(names_from = eventTYPE, values_from = count) %>% replace(is.na(.), 0)
  
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
    group_by(!!Barrier_ID, .drop = F) %>%
    summarise(unique_ind = length(unique(AnimalID)), .groups = 'drop') 
  
  # combine the two dataframes and calculate total encounters and total bad encounters
  barrier_sf_joined <- by_type_df %>% 
    left_join(by_ID, by = quo_name(Barrier_ID)) %>% 
    replace(is.na(.), 0) %>%
    mutate(total_enc = Bounce + Quick_Cross + Average_Movement + Back_n_forth + Trace + unknown + Trapped) 
  
  # calculate impermeability index based on user set expression
  barrier_sf_joined <- barrier_sf_joined %>% filter(total_enc >= min_total_enc) %>% mutate(index = range01(eval(index_fun)))
  
  # put back into spatial
  barrier_sf_joined <- merge(barrier_sf, barrier_sf_joined, by = quo_name(Barrier_ID), all.x = T)
  
  if(show_plot) {
    plot(st_geometry(barrier_sf), col = "grey")
    plot(barrier_sf_joined['index'], add = T)
  }
  
  return(barrier_sf_joined)
  
}
