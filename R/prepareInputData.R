#' Prepare line transect, known fate CMR, and rodent covariate data for integrated analysis
#'
#' @param d_trans tibble containing information on transects (events). Output of
#' wrangleData_LineTrans(). 
#' @param d_obs tibble containing information on observations made along transects 
#' (distance to transect line, numbers of birds in each age/sex class observed,
#' etc.). Output of wrangleData_LineTrans(). 
#' @param d_cmr list with 2 elements. Surv1 and Surv2 are matrices of individuals 
#' released (column 1) and known to have survived (column 2) in each year (row)
#' for season 1 and season 2, respectively. Output of wrangleData_CMR().
#' @param d_rodent list containing the matrix with average number of transects 
#' with rodent observations per area and year and the mean and standard deviation
#' of the original covariate. 
#' @param d_pca dataframe containing grid cell specific environmental habitat information (PCA1 and PCA2)
#' @param grid_df dataframe containing grid cell centroids and area information
#' @param grid_area_lookup dataframe containing info about which transect comes from which area and lies within which grid cell
#' @param localities vector of strings listing localities to consider. Either localities or areas must be provided. 
#' @param areas vector of strings listing areas to consider. Either localities or areas must be provided. 
#' @param areaAggregation logical. If TRUE, areas are used as smallest spatial unit. If FALSE, locations (within areas) are used as smallest spatial unit.
#' @param excl_neverObs logical. If TRUE (default), transects on which ptarmigans were never observed are excluded. If FALSE, all transects are included.
#' @param R_perF logical. If TRUE, treats recruitment rate as juvenile per adult female.
#' If FALSE, treats recruitment rate as juvenile per adult (sum of both sexes).
#' @param R_parent_drop0 logical. If TRUE, removes observations of juveniles without adults
#' from recruitment data. If FALSE, sets 1 as the number of adults/adults females when none
#' are observed. 
#' @param sumR.Level character string. Default ("group") summarises reproduction/recruitment
#' data at the group/observation level. Setting to "line" summarises data at the 
#' transect line level instead. 
#' @param dataVSconstants logical. If TRUE (default) returns a list of 2 lists
#' containing data and constants for analysis with Nimble. If FALSE, returns a
#' list containing all data and constants. 
#' @param addDummyDim logical. If TRUE (default) adds a dummy dimension when only
#' one area is included to allow running the multi-area setup. If FALSE, drops
#' the dummy dimension for running with the single area setup. 
#' @param save logical. If TRUE (default) saves prepared data in working 
#' directory as .rds file.
#'
#' @return A list or list of lists, depending on argument `dataVSconstants` 
#' (see above).
#' @export
#'
#' @examples



prepareInputData <- function(d_trans, d_obs, d_rodent, d_pca, 
                             # basis_df, 
                             grid_df, 
                             transect_grid_lookup, obs_segment_lookup,
                             localities = NULL, areas = NULL, areaAggregation, 
                             excl_neverObs = TRUE, R_perF, R_parent_drop0, 
                             sumR.Level = "group", dataVSconstants = TRUE, 
                             addDummyDim = TRUE, save = TRUE){
  
  
  # Multi-area setup #
  #------------------#
  
  ## Assignment of spatial units
  if(areaAggregation){
    sUnits <- areas
  }else{
    sUnits <- localities
  }
  N_sUnits <- length(sUnits)
  
  ## Rename appropriate column in line transect data to reflect level of spatial aggregation
  if(areaAggregation){
    colnames(d_trans)[which(colnames(d_trans) == "verbatimLocality")] <- "spatialUnit"
    colnames(d_obs)[which(colnames(d_obs) == "verbatimLocality")] <- "spatialUnit"
  }else{
    colnames(d_trans)[which(colnames(d_trans) == "locality")] <- "spatialUnit"
    colnames(d_obs)[which(colnames(d_obs) == "locality")] <- "spatialUnit"
  }
  
  ## If desired: remove all transects on which willow ptarmigans were never encountered
  if(excl_neverObs){
    d_trans <- d_trans %>%
      dplyr::filter(locationID %in% d_obs$locationID)
  }
  
  ## Build observation-grid lookup so each observation is linked to a grid cell through coordinates
  # Should this be replaced by observation-segment lookup? Or added?
  d_obs <- d_obs %>%
    left_join(obs_segment_lookup %>%
        select(eventID, grid_id, segmentID), by = "eventID")
  
  ## Variables shared across areas
  # Number of age classes
  N_ageC <- 2
  
  # (Number of) years
  range_yearsTot <- min(d_trans$Year):max(d_trans$Year)     
  idx_yearsTot <- range_yearsTot - min(d_trans$Year) + 1  
  N_yearsTot <- length(range_yearsTot)
  
  # Truncation distance
  W <- 200
  
  ## Count sites per area (for defining correct array dimensions)
  site_count <- d_trans %>%                             
    dplyr::group_by(spatialUnit) %>%
    dplyr::summarise(count = n_distinct(locationID))
  
  if(length(which(site_count$count == 0)) > 0){
    stop("One or more specified areas contain no willow ptarmigan observations")
  }
  
  ## Count grid cells per area and set up indexing structure for grid_df 
  grid_count <- grid_df %>%
    dplyr::group_by(spatialUnit) %>%
    dplyr::summarise(count = n())
  # grid_df <- grid_df %>%
  #   mutate(rowX = row_number())
  
  ## Count observations per area
  obs_count <- d_obs %>% 
    dplyr::filter(between(DistanceToTransectLine, -0.1, W)) %>% 
    dplyr::group_by(spatialUnit) %>%
    dplyr::summarise(count = n())
    
  ## Set up arrays for area-specific data
  N_sites <- N_years <- min_years <- max_years <- N_obs <- N_sumR_obs <- rep(NA, N_sUnits) 
  N_grid <- rep(NA, N_sUnits) # Number of grid cells within an area
  N_segments <- rep(NA, N_sUnits) # Number of transect segments within an area. Like N_sites but for parts of the transect line
  max_segments <- max(transect_grid_lookup %>%
      group_by(spatialUnit) %>%
      summarise(n = n()) %>%
      pull(n))
  
  transect_grid <- matrix(NA, nrow = N_sUnits,ncol = max(site_count$count))
  
  segment_grid <- matrix(NA, nrow = N_sUnits, ncol = max_segments)
  segment_prop <- matrix(NA, nrow = N_sUnits, ncol = max_segments)
  segment_length <- matrix(NA, nrow = N_sUnits, ncol = max_segments)
  
  A <- matrix(NA, nrow = N_sUnits, ncol = N_yearsTot)
  y <- Year_obs <- zeros_dist <- matrix(NA, nrow = N_sUnits, ncol = max(obs_count$count))
  sumR_obs <- sumR_obs_year <- sumAd_obs <- matrix(NA, nrow = N_sUnits, ncol = max(obs_count$count))
  sumR_obs_grid <- matrix(NA, nrow = N_sUnits, ncol = max(obs_count$count))
  
  L <- N_line_year <- N_J_line_year <- N_A_line_year <- array(0, dim = c(N_sUnits, max(site_count$count), N_yearsTot))
  N_segment_year <- N_J_segment_year <- N_A_segment_year <- array(0, dim = c(N_sUnits, max_segments, N_yearsTot))
  L_segment <- array(0, dim = c(N_sUnits, max_segments, N_yearsTot))
  L_grid <- array(0, dim = c(N_sUnits, max(grid_count$count), N_yearsTot))
  
  N_a_line_year <- array(0, dim = c(N_sUnits, N_ageC, max(site_count$count), N_yearsTot))
  N_a_segment_year <- array(0, dim = c(N_sUnits, N_ageC, max_segments, N_yearsTot))
  
  # X_spatial <- array(NA, dim = c(N_sUnits, max(site_count$count), K)) # Spatial spline structure
  # X_grid <- array(0, dim = c(N_sUnits, max(grid_count$count), K)) # Spatial spline structure
  
  # basis_lookup <- basis_df %>%
  #   mutate(rowX = 1:n()) # Not sure why we are doing this
  # grid_lookup <- grid_df %>%
  #   mutate(rowX = 1:n())
  
  transect_grid_segments <- vector("list", N_sUnits)
  
  for(x in 1:N_sUnits){
    
    ## Subset data (specific area)
    # if(areaAggregation){
    #   d_trans_sub <- subset(d_trans, spatialUnit == areas[x])
    #   d_obs_sub <- subset(d_obs, spatialUnit == areas[x])
    # }else{
    #   d_trans_sub <- subset(d_trans, spatialUnit == localities[x])
    #   d_obs_sub <- subset(d_obs, spatialUnit == localities[x])
    # }
    if(!(sUnits[x] %in% d_trans$spatialUnit)){
      stop(paste0("Spatial unit ", sUnits[x], " (index ", x, ") is not in the data."))
    }
    
    d_trans_sub <- subset(d_trans, spatialUnit == sUnits[x])
    d_obs_sub <- subset(d_obs, spatialUnit == sUnits[x])
    
    # Constants #
    #-----------#
    
    ## Numbers of sites
    N_sites[x] <- n_distinct(d_trans_sub$locationID)
    
    ## (Number of) years with data
    years_mon <- sort(unique(d_trans_sub$Year)) - min(range_yearsTot) + 1 
    min_years[x] <- min(years_mon)
    max_years[x] <- max(years_mon)
    
    # Spatial stuff #
    #---------------#
    
    # # Simplest spatial version: only spline on initial density, no indexing changes
    # sites_area <- d_trans_sub %>%
    #   distinct(locationID) %>%
    #   arrange(locationID)
    # 
    # sites_area <- sites_area %>%
    #   left_join(basis_lookup,
    #             by = "locationID")
    # 
    # X_area <- X[sites_area$rowX, , drop = FALSE]
    # X_spatial[x, 1:N_sites[x], 1:K] <- X_area
    
    # # New version of spatial structure: density on grid level rather than site level
    # sites_area <- d_trans_sub %>%
    #   distinct(locationID) %>%
    #   arrange(locationID)
    # 
    # sites_area <- sites_area %>%
    #   left_join(grid_area_lookup %>%
    #       select(locationID, grid_id), by = "locationID")
    # 
    # grid_area <- grid_df[grid_df$spatialUnit == sUnits[x], ]
    # 
    # grid_area <- dplyr::arrange(grid_area, grid_id)
    #
    # N_grid[x] <- nrow(grid_area)
    # 
    # local_grid_index <- match(sites_area$grid_id, grid_area$grid_id)
    # 
    # # X_area <- X[grid_area$rowX, , drop = FALSE]
    # # 
    # # X_grid[x, 1:N_grid[x], 1:K] <- X_area
    # transect_grid[x, 1:N_sites[x]] <- local_grid_index
    
    ## Once again a new version, now taking into account that transects can belong to multiple grid cells 
    ## Observations are split between grid cells
    
    grid_area <- grid_df[grid_df$spatialUnit == sUnits[x], ]
    grid_area <- dplyr::arrange(grid_area, grid_id)
    
    N_grid[x] <- nrow(grid_area)
    
    transect_segments_area <- transect_grid_lookup %>%
      dplyr::filter(spatialUnit == sUnits[x])
    transect_segments_area <- transect_segments_area %>%
      arrange(segmentID)
    
    transect_segments_area$grid_local <- match(transect_segments_area$grid_id, grid_area$grid_id)
    
    N_segments[x] <- nrow(transect_segments_area)
    
    segment_grid[x, 1:N_segments[x]] <- transect_segments_area$grid_local
    
    segment_prop[x, 1:N_segments[x]] <- transect_segments_area$propLength
    
    segment_length[x, 1:N_segments[x]] <- transect_segments_area$segmentLength
    
    transect_grid_segments[[x]] <- transect_segments_area
    
    
    # Transect characteristics #
    #--------------------------#
    
    ## Transect lengths
    TransLen <- d_trans_sub %>% 
      dplyr::select(locationID, Year, sampleSizeValue) %>%
      reshape2::dcast(locationID~Year, value.var = "sampleSizeValue", sum) %>%
      arrange(locationID)
    
    TransLen_mat <- TransLen %>% select(-locationID) %>% as.matrix()
    colnames(TransLen_mat) <- as.numeric(colnames(TransLen_mat)) - min(range_yearsTot) + 1
    
    segment_effort_lookup <- transect_segments_area %>%
      dplyr::select(segmentID, locationID, propLength)
    
    for(t in 1:N_yearsTot){
      
      ## ORiginal transect effort
      
      if(t %in% years_mon){
        
        L[x, 1:N_sites[x], t] <- TransLen_mat[
          1:N_sites[x],
          which(years_mon == t)
        ]
      }
      
     ## New transect effort
      
      yr <- range_yearsTot[t]
      
      effort_year <- d_trans_sub %>%
        dplyr::filter(Year == yr) %>%
        dplyr::select(locationID, sampleSizeValue)
      
      segment_effort_year <- segment_effort_lookup %>%
        dplyr::left_join(effort_year, by = "locationID") %>%
        dplyr::mutate(sampleSizeValue = ifelse(is.na(sampleSizeValue), 0, sampleSizeValue), segmentEffort = sampleSizeValue * propLength)
      
      L_segment[x, 1:N_segments[x], t] <- segment_effort_year$segmentEffort
    }
    
    for(t in 1:N_yearsTot){
      for(g in 1:N_grid[x]){
        seg_idx <- which(segment_grid[x, 1:N_segments[x]] == g)
        if(length(seg_idx) > 0){
          L_grid[x, g, t] <- sum(L_segment[x, seg_idx, t], na.rm = TRUE)
        }
      }
    }
    
    ## Total covered area 
    A[x, 1:N_yearsTot] <- colSums(L[x,,])*W*2
    
    
    # Observation distance from transect #
    #------------------------------------#
    
    temp_dist <- d_obs_sub %>% 
      dplyr::filter(between(DistanceToTransectLine, -0.1, W)) %>% 
      dplyr::select(Year, DistanceToTransectLine) %>%
      dplyr::mutate(Year2 = Year - (minYear) + 1)
    
    ## Number of observations
    N_obs[x] <- length(temp_dist$DistanceToTransectLine)
    
    ## Distance to transect line
    y[x, 1:N_obs[x]] <- temp_dist$DistanceToTransectLine
    
    ## Observation year
    Year_obs[x, 1:N_obs[x]] <- temp_dist$Year2
    
    ## Vector of 0's of same length as y
    zeros_dist[x, 1:N_obs[x]] <- rep(0, N_obs[x])
    
    
    # Number of birds/line (pooled age classes) #
    #-------------------------------------------#
    
    temp <- TransLen %>% select(locationID)

    TaksObs <- d_obs_sub %>% 
      filter(between(DistanceToTransectLine, -0.1, W)) %>% # 
      dplyr::mutate(cs = unknownJuvenile+unknownunknown+FemaleAdult+MaleAdult) %>%
      reshape2::dcast(locationID~Year, value.var = "cs", sum) %>%
      dplyr::right_join(., temp, by = c("locationID" = "locationID")) %>%
      replace(., is.na(.), 0) %>%
      dplyr::arrange(locationID)
  
    TaksObs_mat <- TaksObs %>% dplyr::select(-locationID) %>% as.matrix() 
    colnames(TaksObs_mat) <- as.numeric(colnames(TaksObs_mat)) - min(range_yearsTot) + 1
    
    years_mon <- as.numeric(colnames(TaksObs_mat)) # filtering observations on DistanceToTransectLine will in some cases remove all observations for a transect line within 1 year. This recalculates number of years with samples
    
    for(t in 1:N_yearsTot){
      
      if(t %in% years_mon){
        N_line_year[x, 1:N_sites[x], t] <- TaksObs_mat[1:N_sites[x], which(years_mon == t)]
      }
    }
    
    # Number of birds/line (by age class) #
    #-------------------------------------#
    
    ## Juveniles (& unknowns)
    TaksObs_J <- d_obs_sub %>% 
      filter(between(DistanceToTransectLine, -0.1, W)) %>%
      dplyr::mutate(cs = unknownJuvenile + unknownunknown) %>%
      reshape2::dcast(locationID~Year, value.var="cs", sum) %>%
      dplyr::right_join(., temp, by=c("locationID"="locationID")) %>%
      replace(., is.na(.), 0) %>%
      dplyr::arrange(locationID)
    
    TaksObs_J_mat <- TaksObs_J %>% dplyr::select(-locationID) %>% as.matrix() 
    colnames(TaksObs_J_mat) <- as.numeric(colnames(TaksObs_J_mat)) - min(range_yearsTot) + 1
    
    for(t in 1:N_yearsTot){
      
      if(t %in% years_mon){
        N_J_line_year[x, 1:N_sites[x], t] <- TaksObs_J_mat[1:N_sites[x], which(years_mon == t)]
      }
    }
   
    ## Adults 
    TaksObs_A <- d_obs_sub %>% 
      filter(between(DistanceToTransectLine, -0.1, W)) %>%
      dplyr::mutate(cs = FemaleAdult + MaleAdult) %>%
      reshape2::dcast(locationID~Year, value.var="cs", sum) %>%
      dplyr::right_join(., temp, by=c("locationID"="locationID")) %>%
      replace(., is.na(.), 0) %>%
      dplyr::arrange(locationID)
    
    TaksObs_A_mat <- TaksObs_A %>% dplyr::select(-locationID) %>% as.matrix() 
    colnames(TaksObs_A_mat) <- as.numeric(colnames(TaksObs_A_mat)) - min(range_yearsTot) + 1
    
    for(t in 1:N_yearsTot){
      
      if(t %in% years_mon){
        N_A_line_year[x, 1:N_sites[x], t] <- TaksObs_A_mat[1:N_sites[x], which(years_mon == t)]
      }
    }
    
    ## Check and combine in array
    if(!(all(N_J_line_year + N_A_line_year == N_line_year))){
      warning("Number of observed adults and juveniles does not add up correctly. Double-check data.")
    }
    
    N_a_line_year[x, 1,,] <- N_J_line_year[x,,]
    N_a_line_year[x, 2,,] <- N_A_line_year[x,,]
    
    
    # Number of birds/segment (by age class) #
    #----------------------------------------#
   
    ## Essentially the same as above, but now aggregated in segments rather than full transect line
    
    ## Pooled age classes
    segment_temp <- transect_segments_area %>%
      select(segmentID)
    
    TaksObs_seg <- d_obs_sub %>%
      filter(between(DistanceToTransectLine, -0.1, W)) %>%
      filter(!is.na(segmentID)) %>%
      mutate(cs = unknownJuvenile + unknownunknown + FemaleAdult + MaleAdult) %>%
      reshape2::dcast(segmentID ~ Year, value.var = "cs", sum) %>%
      right_join(segment_temp, by = "segmentID") %>%
      replace(is.na(.), 0) %>%
      arrange(segmentID)
    
    TaksObs_seg_mat <- TaksObs_seg %>% dplyr::select(-segmentID) %>% as.matrix() 
    colnames(TaksObs_seg_mat) <- as.numeric(colnames(TaksObs_seg_mat)) - min(range_yearsTot) + 1
    
    for(t in 1:N_yearsTot){
      
      if(t %in% years_mon){
        N_segment_year[x, 1:N_segments[x], t] <- TaksObs_seg_mat[1:N_segments[x], which(years_mon == t)]
      }
    }
    
    ## Juveniles & unknowns
    TaksObs_J_seg <- d_obs_sub %>% 
      filter(between(DistanceToTransectLine, -0.1, W)) %>%
      dplyr::mutate(cs = unknownJuvenile + unknownunknown) %>%
      filter(!is.na(segmentID)) %>%
      reshape2::dcast(segmentID~Year, value.var="cs", sum) %>%
      dplyr::right_join(., segment_temp, by=c("segmentID"="segmentID")) %>%
      replace(., is.na(.), 0) %>%
      dplyr::arrange(segmentID)
    
    TaksObs_J_seg_mat <- TaksObs_J_seg %>% dplyr::select(-segmentID) %>% as.matrix() 
    colnames(TaksObs_J_seg_mat) <- as.numeric(colnames(TaksObs_J_seg_mat)) - min(range_yearsTot) + 1
    
    for(t in 1:N_yearsTot){
      
      if(t %in% years_mon){
        N_J_segment_year[x, 1:N_segments[x], t] <- TaksObs_J_seg_mat[1:N_segments[x], which(years_mon == t)]
      }
    }
    
    ## Adults
    TaksObs_A_seg <- d_obs_sub %>% 
      filter(between(DistanceToTransectLine, -0.1, W)) %>%
      dplyr::mutate(cs = FemaleAdult + MaleAdult) %>%
      filter(!is.na(segmentID)) %>%
      reshape2::dcast(segmentID~Year, value.var="cs", sum) %>%
      dplyr::right_join(., segment_temp, by=c("segmentID"="segmentID")) %>%
      replace(., is.na(.), 0) %>%
      dplyr::arrange(segmentID)
    
    TaksObs_A_seg_mat <- TaksObs_A_seg %>% dplyr::select(-segmentID) %>% as.matrix() 
    colnames(TaksObs_A_seg_mat) <- as.numeric(colnames(TaksObs_A_seg_mat)) - min(range_yearsTot) + 1
    
    for(t in 1:N_yearsTot){
      
      if(t %in% years_mon){
        N_A_segment_year[x, 1:N_segments[x], t] <- TaksObs_A_seg_mat[1:N_segments[x], which(years_mon == t)]
      }
    }
    
    N_a_segment_year[x,1,,] <- N_J_segment_year[x,,]
    N_a_segment_year[x,2,,] <- N_A_segment_year[x,,]
    
    # Recruitment #
    #-------------#
    
    # NOTE: For consistency with the population model, we need to define R as the number of recruits
    # per adult (female)
    
    # Reformat data
    temp_Rec <- d_obs_sub %>% filter(between(DistanceToTransectLine, -0.1, W)) %>%
      dplyr::mutate(sumR = unknownJuvenile + unknownunknown,
                    sumAd = MaleAdult + FemaleAdult,
                    sumAdF = FemaleAdult) %>%
      dplyr::mutate(Year2 = Year - (min(Year)) + 1)
  
    # # Optional: summarise data at line level (per year)
    # if(sumR.Level == "line"){
    #   temp_Rec <- temp_Rec %>%
    #     dplyr::group_by(locationID, Year, Year2) %>%
    #     dplyr::summarise(sumR = sum(sumR),
    #                      sumAd = sum(sumAd),
    #                      sumAdF = sum(sumAdF), .groups = "keep")
    # }
    
    # # Attach grid information BEFORE aggregation
    # temp_Rec <- temp_Rec %>%
    #   left_join(obs_grid_lookup %>%
    #       select(eventID, grid_id), by = "eventID")
    # 
    
    # Aggregate by grid-year (earlier this was transect-year level, shoudl eventually change the wording in this toggle from "line" to "grid")
    if(sumR.Level == "line"){
      temp_Rec <- temp_Rec %>%
        dplyr::group_by(grid_id, Year, Year2) %>%
        dplyr::summarise(sumR = sum(sumR), sumAd = sum(sumAd), sumAdF = sum(sumAdF), .groups = "drop")
    }
    
    
    # Set which adult measure to use
    if(R_perF){
      temp_Rec$sumAd_use <- temp_Rec$sumAdF
    }else{
      temp_Rec$sumAd_use <- temp_Rec$sumAd
    }
    
    # Deal with instances of 0 adults observed
    if(R_parent_drop0){ # --> Drop all cases of 0 adults observed
      temp_Rec <- temp_Rec %>%
        dplyr::filter(sumAd_use != 0)
      
    }else{ # --> Drop only cases when neither adult (females) nor juveniles were observed and add +1 to adults otherwise
      
      temp_Rec <- temp_Rec %>%
        dplyr::filter(sumAd_use + sumR != 0) %>%
        dplyr::mutate(sumAd_use = ifelse(sumAd_use == 0, 1, sumAd_use))
    }
    
    # Add column with grid ID information
    # temp_Rec <- temp_Rec %>%
    #   left_join(grid_area_lookup %>%
    #       select(locationID, grid_id),
    #     by = "locationID")
    
    temp_Rec$grid_local <- match(temp_Rec$grid_id, grid_area$grid_id)
    
    # Extract relevant data vectors
    N_sumR_obs[x] <- length(temp_Rec$sumR)
    sumR_obs[x, 1:N_sumR_obs[x]] <- temp_Rec$sumR
    sumAd_obs[x, 1:N_sumR_obs[x]] <- temp_Rec$sumAd_use
    sumR_obs_year[x, 1:N_sumR_obs[x]] <- temp_Rec$Year2
    sumR_obs_grid[x, 1:N_sumR_obs[x]] <- temp_Rec$grid_local
  }
  
  ## Drop excess columns in recruitment data
  sumR_obs <- sumR_obs[1:N_sUnits, 1:max(N_sumR_obs), drop = FALSE]
  sumAd_obs <- sumAd_obs[1:N_sUnits, 1:max(N_sumR_obs), drop = FALSE]
  sumR_obs_year <- sumR_obs_year[1:N_sUnits, 1:max(N_sumR_obs), drop = FALSE]
  sumR_obs_grid <- sumR_obs_grid[1:N_sUnits, 1:max(N_sumR_obs),drop = FALSE]
  
  ## Covariate toggles
  PCA1_on <- ifelse(fitPCA1, 1, 0)
  PCA2_on <- ifelse(fitPCA2, 1, 0)
  MTemp_on <- ifelse(fitMeanTemp, 1, 0)
  TempA_on <- ifelse(fitTempAnom, 1, 0)
  
  
  # Data assembly #
  #---------------#
  
  ## Add dummy dimensions if running for only one spatial unit
  if(N_sUnits == 1 & addDummyDim){
    N_sites <- c(N_sites, NA)
  }
  
  ## Assembling all data in a list
  input.data <- list(
    sumR_obs = sumR_obs, # Observed numbers of recruits
    sumAd_obs = sumAd_obs, # Observed numbers of adults/adult females
    sumR_obs_year = sumR_obs_year, # Year of observed numbers of recruits
    sumR_obs_grid = sumR_obs_grid, # Grid ID of observed number of recruits
    N_sumR_obs = N_sumR_obs, # Total number of observations of numbers of recruits
    
    ## Transect segment information
    segment_grid = segment_grid, # which grid cell does each segment belong to
    segment_prop = segment_prop, # segment proportion of the entire transect line
    segment_length = segment_length, # length of the segments
    N_segments = N_segments, # number of segments per ... grid cell?
    L_segment = L_segment, # segment 'length' corrected for effort
    L_grid = L_grid, # segment length corrected for effort, summed for all segments in a grid cell
    
    y = y, # Distance to transect line for each individual observation
    zeros_dist = zeros_dist, # Vector of 0's of same length as y
    Year_obs = Year_obs, # Year of each observation
    N_obs = N_obs, # Total number of observations
    
    N_line_year = N_line_year, # Number of birds observed per site per year
    N_segment_year = N_segment_year, # Number of birds observed per transect segment per year
    N_a_line_year = N_a_line_year, # Number of birds observed per ageclass per site per year
    N_a_segment_year = N_a_segment_year, # Number of birds observed per age class per segment per year
    L = L, # Transect length per site and year
    
    N_years = N_yearsTot, # Max number of years with data
    min_years = min_years, # Earliest year (index) of monitoring per area
    max_years = max_years, # Last year (index) of monitoring per area
    N_sites = N_sites, # Total number of monitored sites per area
    
    A = A, # Total covered area per year
    W = W, # Truncation distance
    N_ageC = N_ageC, # Number of age classes
    
    ## Rodent covariate
    RodentOcc = d_rodent$rodentAvg,
    RodentOcc_meanCov = d_rodent$meanCov,
    RodentOcc_sdCov = d_rodent$sdCov,
    
    ## Spring temperature covariate
    MeanTemp = d_temp$MeanTemp,
    # SpringTemp_meanCov = d_temp$mean,
    # SpringTemp_sdCov = d_temp$sd,
    TempAnom = d_temp$TempAnom,
    
    ## PCA habitat information 
    PCA1 = d_pca$PCA1,
    PCA2 = d_pca$PCA2,
    
    N_areas = N_sUnits,
    area_names = sUnits,
    
    ## Covariate toggles
    PCA1_on = PCA1_on,
    PCA2_on = PCA2_on,
    MTemp_on = MTemp_on,
    TempA_on = TempA_on
  )
  
  ## Drop dummy dimension when only one area is included
  if(N_sUnits == 1 & !addDummyDim){
    
    for(i in 1:length(input.data)){
      input.data[[i]] <- drop(input.data[[i]])
    }
    
  }
  
  ## Assembling Nimble data
  nim.data <- list(sumR_obs = input.data$sumR_obs, sumAd_obs = input.data$sumAd_obs,
                   y = input.data$y, 
                   zeros.dist = input.data$zeros_dist, 
                   L = input.data$L,
                   L_segment = input.data$L_segment,
                   L_grid = input.data$L_grid,
                   N_line_year = input.data$N_line_year,
                   N_segment_year = input.data$N_segment_year,
                   N_a_line_year = input.data$N_a_line_year, 
                   N_a_segment_year = input.data$N_a_segment,
                   A = input.data$A,
                   RodentOcc = input.data$RodentOcc,
                   MeanTemp = input.data$MeanTemp,
                   TempAnom = input.data$TempAnom
                   # X_spatial = X_spatial
                   # X_grid = X_grid,
                   )
  
  ## Assembling Nimble constants
  nim.constants <- list(N_years = input.data$N_years, min_years = input.data$min_years, max_years = input.data$max_years,
                        W = input.data$W,
                        N_obs = input.data$N_obs, Year_obs = input.data$Year_obs,
                        N_sites = input.data$N_sites, 
                        N_ageC = N_ageC,
                        N_areas = input.data$N_areas, area_names = input.data$area_names,
                        sumR_obs_year = input.data$sumR_obs_year, N_sumR_obs = input.data$N_sumR_obs,
                        sumR_obs_grid = input.data$sumR_obs_grid,
                        N_ageC = N_ageC,
                        # transect_grid_segments = input.data$transect_grid_segments,
                        N_segments = input.data$N_segments,
                        segment_grid = input.data$segment_grid,
                        segment_prop = input.data$segment_prop,
                        segment_length = input.data$segment_length,
                        RodentOcc_meanCov = input.data$RodentOcc_meanCov,
                        RodentOcc_sdCov = input.data$RodentOcc_sdCov,
                        PCA1 = input.data$PCA1,
                        PCA2 = input.data$PCA2,
                        # SpringTemp_meanCov = input.data$SpringTemp_meanCov,
                        # SpringTemp_sdCov = input.data$SpringTemp_sdCov,
                        N_grid = N_grid,
                        # transect_grid = transect_grid)
                        PCA1_on = input.data$PCA1_on,
                        PCA2_on = input.data$PCA2_on,
                        MTemp_on = input.data$MTemp_on,
                        TempA_on = input.data$TempA_on
                        )
  
  ## Make final data list to return
  if(dataVSconstants){
    rype.data <- list(nim.data = nim.data,
                      nim.constants = nim.constants)
  }else{
    rype.data <- input.data
  }
  
  ## Optional: save data as .rds
  if(save){
    saveRDS(rype.data, file = "RypeData_forIM.rds")
  }
  
  ## Return data
  return(rype.data)
}

