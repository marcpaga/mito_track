require(ijtiff)

video_files_names <- sapply(strsplit(list.files('data/rawdata'), '.', fixed = T),'[', 1)


process_mito_track <- function(video_name) {

  
  print(paste0('Processing: ', video_name))
  data_dir <- paste0('data/processed/', video_name)
  tiff_file <- ijtiff::read_tif(path = paste0('data/rawdata/', video_name, '.tif'))
  
  x_max <- dim(tiff_file)[2]
  y_max <- dim(tiff_file)[1]
  n_frames <- dim(tiff_file)[4]
  memb_threshold <- 1.5
  
  print(paste0('Number of frames: ' , n_frames))
  
  # cut file to the number of frames for which we have tracking
  
  ### DATA PROCESSING ### ======================================================
  
  
  # retrieve coordiantes of the membrane
  membrane_coords <- function(mat) {
    
    coor_list <- list()
    u <- 0
    for (i in seq_len(nrow(mat))) {
      for (j in seq_len(ncol(mat))) {
        if (is.na(mat[i, j])) {
          next
        }
        if (i == 1 | i == nrow(mat)) {
          next
        }
        if (j == 1 | j == ncol(mat)) {
          next
        }
        
        top_left <- is.na(mat[i-1,j+1])
        top <- is.na(mat[i,j+1])
        top_right <- is.na(mat[i+1,j+1])
        left <- is.na(mat[i-1,j])
        right <- is.na(mat[i+1,j])
        bot_left <- is.na(mat[i-1,j-1])
        bot <- is.na(mat[i,j-1])
        bot_right <- is.na(mat[i+1,j-1])
        surr <- c(top_left,top,top_right,left,right,bot_left,bot_right)
        
        if (any(surr)) {
          u <- u + 1
          coor_list[[u]] <- c(j,i)
        }
        
      }
    }
    return(coor_list)
  }
  
  print('Calculating membrane coordinates')
  
  # membrane coords of all the frames
  # frames as a list of matrices
  frames_mat <- list()
  for (i in seq_len(dim(tiff_file)[4])) {
    
    mat <- tiff_file[,,,i]
    mat[log10(mat) < memb_threshold] <- NA
    frames_mat[[i]] <- mat
    
  }
  
  mem_coords <- lapply(frames_mat, membrane_coords)
  
  temp_list <- list()
  for (i in seq_along(mem_coords)) {
    
    tdf <- as.data.frame(do.call('rbind', mem_coords[[i]]))
    tdf$t <- i
    temp_list[[i]] <- tdf
    
  }
  mem_df <- do.call('rbind', temp_list)
  write.csv(mem_df, file = paste0('data/processed/membrane_', video_name, '.csv'))
  
  print(paste0('Number of membrane coordinates: ' , length(mem_coords[[1]])))
  
  
  
  
  print('Reading tracking files')
  
  insertRow2 <- function(existingDF, newrow, r) {
    existingDF <- rbind(existingDF,newrow)
    existingDF <- existingDF[order(c(1:(nrow(existingDF)-1),r-0.5)),]
    row.names(existingDF) <- 1:nrow(existingDF)
    return(existingDF)  
  }
  
  # read the files
  for (i in seq_along(list.files(data_dir))) {
    if (i == 1) {
      temp_list <- list()
      file_name <- paste0(data_dir, '/', 'points_', i, '.csv')
      new <- read.csv(file = file_name, 
                      header = FALSE, 
                      col.names = c('x', 'y'))
      new$keep <- 1
      ghost_ids <- c()
    } else {
      file_name <- paste0(data_dir, '/', 'points_', i, '.csv')
      new <- read.csv(file = file_name, 
                      header = FALSE, 
                      col.names = c('x', 'y', 'keep'))
    }
    
    if (length(ghost_ids) != 0) {
      for (j in sort(ghost_ids, decreasing = FALSE)) {
        new_row <- c(NA, NA, NA)
        new <- insertRow2(new, new_row, j)
      }
    }
    
    if (0 %in% new$keep) {
      ghost_ids <- c(ghost_ids, which(new$keep == 0))
    }
    
    new$t <- i
    new$id <- seq_len(nrow(new))
    temp_list[[i]] <- new
    
  }
  trace_df <- do.call('rbind', temp_list)
  trace_df <- subset(trace_df, !is.na(x))
  
  print('Processing tracking data')
  
  # center coordinates 
  trace_df$x_centered <- trace_df$x - x_max/2
  trace_df$y_centered <- trace_df$y - y_max/2
  
  # calculate the speed
  for (i in seq_along(unique(trace_df$id))) {
    if (i == 1) {
      temp_list <- list()
    }
    tempDf <- subset(trace_df, id == unique(trace_df$id)[i])
    tempDf$v <- 0
    for (j in seq_len(nrow(tempDf))) {
      if (j == 1) {
        next
      }
      tempDf$v[j] <- sqrt((tempDf$x[j] - tempDf$x[j-1])^2 + 
                            (tempDf$y[j] - tempDf$y[j-1])^2)
    }
    temp_list[[i]] <- tempDf
  }
  trace_df <- do.call('rbind', temp_list)
  
  # calculate the movement direction
  for (i in seq_along(unique(trace_df$id))) {
    if (i == 1) {
      temp_list <- list()
    }
    tempDf <- subset(trace_df, id == unique(trace_df$id)[i])
    tempDf$angle <- NA
    for (j in seq_len(nrow(tempDf))) {
      if (j == 1) {
        next
      }
      v2 <- c(tempDf$x_centered[j], tempDf$y_centered[j])
      v1 <- c(tempDf$x_centered[j-1], tempDf$y_centered[j-1])
      angle <- (atan2(y = (v2[2] - v1[2]), 
                      x = (v2[1] - v1[1])) * 180)/pi
      if (angle < 0) {
        angle <- angle + 360
      }
      tempDf$angle[j] <- angle
      
    }
    temp_list[[i]] <- tempDf
  }
  trace_df <- do.call('rbind', temp_list)
  
  # relativize the coordinates as t=0 is (0,0)
  temp_list <- list()
  for (i in seq_along(unique(trace_df$id))) {
    tempDf <- subset(trace_df, id == unique(trace_df$id)[i])
    tempDf$x_sca <- tempDf$x - tempDf$x[which(tempDf$t == 1)]
    tempDf$y_sca <- tempDf$y - tempDf$y[which(tempDf$t == 1)]
    temp_list[[i]] <- tempDf
  }
  trace_df <- do.call('rbind', temp_list)
  
  print('Calculating distance to membrane')
  
  # calculate the minimum distance of each tracking thing to the membrane in 
  # each frame
  trace_df$dist_memb <- NA
  for (i in seq_len(nrow(trace_df))) {
    
    time <- trace_df$t[i] 
    temp_coords <- mem_coords[[time]]
    tdf <- do.call('rbind', temp_coords)
    dot_pos <- c(trace_df$x[i], trace_df$y[i])
    
    trace_df$dist_memb[i] <- min(sqrt((dot_pos[1] - tdf[,1])^2 + (dot_pos[2] - tdf[,2])^2))
    
  }
  
  write.csv(trace_df, file = paste0('data/processed/', video_name, '.csv'))
  
}

process_mito_track(video_files_names)

