#####################################
############DEPENDENCIES#############
#####################################

library(foveafit)

#####################################
##############FUNCTIONS##############
#####################################


# Gets a smoothened version of the ellipse by first making a mean of the two rotations. Then the mean ellipse is split into 4 and each segment is smoothed before it is stuck together. 
# Currently thee is a slight error in that there are too many points at the top of the ellipse and ~2 points missing from the bootom. Needs to be fixed but still works pretty well. 
# Currently generally running using a smoothing factor of 30, but this is random so can be changed

smooth_ellipse_max_slope <- function(e, smoothf){
  mean_ellipse_table <- NULL
  mean_ellipse_table <- cbind(mean_ellipse_table, e$fitted_center$max_slope[1:360,])
  mean_ellipse_table <- rbind(mean_ellipse_table, e$fitted_center$max_slope[541:720,])
  mean_ellipse_table <- rbind(mean_ellipse_table, e$fitted_center$max_slope[361:540,])
  mean_ellipse <- NULL
  mean_ellipse <- cbind(mean_ellipse, (mean_ellipse_table[1:360,1]+mean_ellipse_table[361:720,1])/2)
  mean_ellipse <- cbind(mean_ellipse, (mean_ellipse_table[1:360,2]+mean_ellipse_table[361:720,2])/2)
  smoothed_ellipse <-cbind(runmed(mean_ellipse[,1], smoothf), runmed(mean_ellipse[,2], smoothf))
  plot(smoothed_ellipse)
  return(smoothed_ellipse)
}

# Finds the rotation from the y axis in degrees. Have to remove 90o because calculates as degree from x axis going clockwise. 
find_rotation <- function(smoothed_ellipse){
  rotation_data <- NULL
  rotation_data <- cbind(rotation_data, smoothed_ellipse[1:180,])
  rotation_data <- cbind(rotation_data, smoothed_ellipse[181:360,])
  distance_value <- sqrt(((rotation_data[,1]-rotation_data[,3])^2)+(rotation_data[,2]-rotation_data[,4])^2)
  rotation_data <- cbind(rotation_data, distance_value)
  rotation_data_order <- rotation_data[order(rotation_data[,5]),]
  degree<- match(rotation_data_order[1,], rotation_data)
  degree <- degree[1]
  degree <- degree - 90
  return(degree)
}


#######################################
###############ANALYSIS################
#######################################
rotates <- NULL
for (x in all_eye_paths){
  image_file <- x
  e <- registar(image_file, debug=FALSE)
  recur <-(is.recursive(e))
  recur2 <- (is.recursive(e$fitted_center))
  print(recur)
  print(recur2)
  smooth_ellipse <- smooth_ellipse_max_slope(e, 30)
  rotate_value <- find_rotation(smooth_ellipse)
  rotates <- rbind(rotates, rotate_value)
}