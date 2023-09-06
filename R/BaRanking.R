
range01 <- function(x){(x-min(x))/(max(x)-min(x))}

## default index_fun is to calculate the index used in Xu et al. 2021
BaRanking <- function(classification, barrier, d, Barrier_ID, min_total_enc = 0, 
                      index_fun = expression(((Bounce + Back_n_forth + Trace + Trapped)/total_enc)*unique_ind), 
                      show_plot = F) {
 
  Barrier_ID <- rlang::sym(Barrier_ID)
  
  ## make classification spatial
  classification_sf <- sf::st_as_sf(classification, coords = c("easting", "northing"), crs = sf::st_crs(barrier))
  
  ## spatial join by the fence buffer distance used in BaBA()
  barrier_sf_joined <- sf::st_join(barrier_sf, classification_sf, join = st_is_within_distance, dist = d)
  
  ## calculate # of each encounter event types
  by_type <- 
    barrier_sf_joined %>%
    ## only keep those with animal encounters
    dplyr::filter(!is.na(barrier_sf_joined$AnimalID)) %>% 
    ## The bang-bang operator !! forces a single object. One common case for !! is to substitute an environment-variable (created with <-) with a data-variable (inside a data.frame).
    dplyr::group_by(!!Barrier_ID, eventTYPE, .drop = F) %>%    
    dplyr::summarise(count = dplyr::n(), .groups = 'drop') %>% 
    ## No longer want a spatial object
    sf::st_drop_geometry() %>% 
    ## change to wide format with one column per fence segment x eventTYPE
    tidyr::pivot_wider(names_from = 'eventTYPE', values_from = 'count') %>%
    replace(is.na(.), 0)
  
  ## ensure all behaviore types are listed in by_type, filling in 0s for any that are missing
  ba.all <- c("Quick_Cross", "Bounce", "Average_Movement", "Back_n_forth", "Trace", "Trapped", "unknown") # a full list of names
  ba.miss <- ba.all[!(ba.all %in% names(by_type))]
  add.df <- as.data.frame(matrix(0, nrow = nrow(by_type), ncol = length(ba.miss)))
  colnames(add.df) <- ba.miss
  by_type <- cbind(by_type, add.df)
  
  ## calculate # of unique individuals encountering each fence
  by_ID <- 
    barrier_sf_joined %>% 
    dplyr::filter(!is.na(barrier_sf_joined$AnimalID))  %>% 
    dplyr::group_by(!!Barrier_ID, .drop = F) %>%
    dplyr::summarise(unique_ind = length(unique(AnimalID)), .groups = 'drop') %>% 
    sf::st_drop_geometry()
  
  ## combine the two tibbles
  barrier_encounters <- 
    by_type %>% 
    left_join(by_ID, by = rlang::as_name(Barrier_ID)) %>% 
    replace(is.na(.), 0) %>%
    dplyr::mutate(
      ## calculate total encounters
      total_enc = Bounce + Quick_Cross + Average_Movement + Back_n_forth + Trace + unknown + Trapped,
      ## calculate the impermeability index based on the user-set expression for all fence segments with sufficient encounters (total_enc >= min_total_enc).
      ## this must be split in two steps so that rescaling the index only considers values where total_enc >= min_total_enc.
      calc_expr = dplyr::if_else(total_enc >= min_total_enc, eval(index_fun), NA),
      index = range01(calc_expr)) %>% 
    ## calc_expr is no longer needed  
    select(-calc_expr)
  
  ## put back into spatial format
  barrier_encounters_sf <- merge(barrier_sf, barrier_encounters, by = rlang::as_name(Barrier_ID), all.x = TRUE)
  
  if(show_plot) {
    plot(st_geometry(barrier_encounters_sf), col = "grey")
    plot(barrier_encounters_sf['index'], add = TRUE)
  }
  
  return(barrier_encounters_sf)
  
}
