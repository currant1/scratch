###################################
############DEPENDECIES############
###################################

library(foveafit)
library(colorspace)

library(NISTunits)
library(mgcv)

library(gridExtra)
library(cowplot)

library(dplyr)


###################################
############FUNCTIONS##############
###################################

convert_coords <- function(v360,divs){
  f.xy1 = NULL
  f.xy2 = NULL
  for (ang in 1:nrow(divs)) {
    ang2 = abs(ang + 179)
    if (ang2 > 360) {
      ang2 = ang2 - 360
    }
    inds1 = unlist(abs(floor(v360$coors[ang, ])))
    inds2 = unlist(abs(floor(v360$coors[ang2, ])))
    inds1[inds1[1] == 0][1] = 1
    inds1[inds1[2] == 0][2] = 1
    inds2[inds2[1] == 0][1] = 1
    inds2[inds2[2] == 0][2] = 1
    s1 = seq(inds1[1], inds2[1], length = v360$n)
    s2 = seq(inds1[2], inds2[2], length = v360$n)
    f.xy1 = rbind(f.xy1, s1[divs[ang, ]])
    f.xy2 = rbind(f.xy2, s2[divs[ang, ]])
  }
  return(list(f.xy1 = f.xy1, f.xy2 = f.xy2))
}

# creates a reference matrix that can be scaled based on the size of the fovea and centred so that the fovea is in the middle of the grid
scale_and_centre <- function(amodel, fovea_height, fovea_coords, n){
  fovea_height <- fovea_height/dim(amodel)[1]
  fovea_coords[1] <- fovea_coords[1]/dim(amodel)[2]
  fovea_coords[2] <- fovea_coords[2]/dim(amodel)[1]
  tmp_mat <- create_matrix(0, ((dim(amodel)[2])-1), 0, ((dim(amodel)[1])-1))
  rownames(tmp_mat) <- paste(tmp_mat[,1], ".", tmp_mat[,2], sep="")
  id_mapping <- data.frame(id=paste(tmp_mat[,1], ".", tmp_mat[,2], sep=""), x=tmp_mat[,1], y=tmp_mat[,2], stringsAsFactors = FALSE)
  d1 <- dim(amodel)[1]*(1/3)
  fov_height_scale <- fovea_height*dim(amodel)[1]
  fovea_coords_big <- c(fovea_coords[1]*dim(amodel)[2], fovea_coords[2]*dim(amodel)[1])
  scale_co <- d1/fov_height_scale
  fovea_coords_scale <- fovea_coords_big*scale_co
  tmp_mat_scale <- tmp_mat*scale_co
  centre_mat <- c((dim(amodel)[1]*scale_co)/2, (dim(amodel)[2]*scale_co)/2)
  foveal_difference <- c(centre_mat[1]-fovea_coords_scale[1], centre_mat[2]-fovea_coords_scale[2])
  scale_and_shift <- tmp_mat_scale
  scale_and_shift[,1] <- scale_and_shift[,1] + foveal_difference[2]
  scale_and_shift[,2] <- scale_and_shift[,2] + foveal_difference[1]
  scale_shift_small <- scale_and_shift
  scale_shift_small[,1] <- scale_shift_small[,1]/dim(amodel)[2]
  scale_shift_small[,2] <- scale_shift_small[,2]/dim(amodel)[1]
  scale_shift_small <- scale_shift_small[complete.cases(scale_shift_small),]
  return(list(ref_table = scale_shift_small, id_mapping=id_mapping))
  
}

# Makes the standard grid of triangles. takes in amodel and n, the number of triangles desired across the y axis (currently using 20)
make_tri_grid <- function(amodel, n){
  negative_range <- seq(from=0, to=1 + 1.732051, by=1/n)
  positive_range <- seq(from=1, to=-1.732051, by=-1/n)
  horizontal <- seq(from=0, to=1, by=1/(n*2))
  image(amodel)
  s = sapply(positive_range, function(n) abline(a=n, b=1.732051))
  s = sapply(negative_range, function(n) abline(a=n, b=-1.732051))
  s = sapply(horizontal, function(n) abline(a=n, b=0))
  indicies=NULL
  for (l in 1:length(positive_range)){
    for (j in 1:length(negative_range)){
      coords <- get_intersect(1.732051, positive_range[l], -1.732051, negative_range[j])
      indicies <- rbind(indicies, coords)
      if(sum(coords>=0)==2) {
        matplot(coords$x, coords$y, col="blue", pch=20, add=T)
      }
    }
  }
  return(indicies)
}

#removes all indices that lie outside of the 1x1 boundary
remove_to_square <- function(inds){
  inds_square <- inds[inds [,1]>=0,]
  inds_square <- inds_square[inds_square[,1]<=1,]
  inds_square <- inds_square[inds_square[,2]>=0,]
  inds_square <- inds_square[inds_square[,2]<=1,]
  plot(inds_square)
  return(inds_square)
}

#get three vertices of each triangle and scale to values of model
get_triangle_vertices <- function(inds_square, amodel){
  sorted_square <- inds_square[order( unlist(inds_square[,2]), unlist(inds_square[,1])),]
  sorted_square <- lapply(sorted_square, round, 8)
  sorted_square <- matrix(unlist(sorted_square), ncol=2)
  sorted_square <- unique(sorted_square[,1:2])
  yunique <- unique(round(sorted_square[,2], 8))
  twolines_table_order <- NULL
  for (x in 1:(length(yunique)-1)) {
    print(x)
    twolines <- NULL
    twolines <- rbind(twolines, sorted_square[which(sorted_square[,2]==yunique[[x]]),])
    twolines <- rbind(twolines, sorted_square[which(sorted_square[,2]==yunique[[x+1]]),])  
    twolines_order <- twolines[order(unlist(twolines[,1])),]
    twolines_order <- unique(round(twolines_order, 6))
    twolines_order <- unique(twolines_order)
    twolines_table_order <- cbind2(twolines_table_order, twolines_order)
    #twolines_table <- cbind2(twolines_table, twolines)
    #twolines_table_order <- twolines_table[order(unlist(twolines_table[,1])),]
  }
  vertices<- NULL
  for (y in seq(from=1, to = dim(twolines_table_order)[2], by=2)){
    for (z in 1:(dim(twolines_table_order)[1]-2)){
      trianglelong <- NULL
      triangle1 <- twolines_table_order[z:(z+2),y:(y+1)]
      #trianglelong <- cbind(trianglelong, unlist(triangle1[1,]), unlist(triangle1[2,]), unlist(triangle1[3,]))
      vertices <- rbind(vertices, triangle1)
    }
  }
  #vertices[,1] <- (unlist(vertices[,1]))*dim(amodel)[1]
  #vertices[,2] <- (unlist(vertices[,2]))*dim(amodel)[2]
  verticesdf <- matrix(c(unlist(vertices)), nrow=dim(vertices), ncol=2)
  return(verticesdf)
}

#get grid saying which triangle each coordiante of the image is part of
get_triangle_ref_coords <- function(m, verticesdf){
  all_refs <- NULL
  triangle_groups <- seq(from=1, to=(dim(verticesdf)[1]), by=3)
  for (x in triangle_groups){
    small_matrix <- verticesdf[x:(x+2),]
    inoutvalues <- in.out(small_matrix, m)
    TRUTHS <- which(inoutvalues)
    
    refcoords <- NULL
    for (y in 1:length(TRUTHS)){
      tmp <- NULL
      tmp <- rbind(tmp,(m[TRUTHS[y],]))
      tmp <- cbind(tmp, as.integer((x+2)/3))
      refcoords <- rbind(refcoords, tmp)
    }
    rownames(refcoords) <- rownames(m)[TRUTHS]
    #refcoors <- cbind(rownames(m)[TRUTHS], refcoords)
    all_refs<-rbind(all_refs,refcoords)
    
  }
  return(all_refs)
}

# find the means across the triangles
tesselate_mean <- function(all_refs, amodel){
  triangleids <- unique(all_refs[,3])
  all_triangle_means <- NULL
  for (x in 1:length(triangleids)){
    xth_triangle <- subset(all_refs, all_refs[,3]==triangleids[x])
    triangle_values=NULL
    for (y in 1:dim(xth_triangle)[1]){
      triangle_values <- rbind(triangle_values, amodel[(xth_triangle[y,])[1], (xth_triangle[y,])[2]])
    }
    triangle_mean <- colMeans(triangle_values)
    triangle_mean <- cbind(triangle_mean, triangleids[x])
    all_triangle_means <- rbind(all_triangle_means, triangle_mean)
  }
  return(all_triangle_means)
}

# Finds the pixel value at each of the coords
Ids2coordinates <- function(id, original) {
  splitid <- as.numeric(unlist(strsplit(id[1], "[.]")))
  xcoord <- splitid[1]+1
  ycoord <- splitid[2]+1
  data.frame(id[1], original[ycoord,xcoord], stringsAsFactors=FALSE)
}

# Finds the pixel value at each of the coords
Ids2coordinates_rotate <- function(id, original) {
  splitid <- as.numeric(unlist(strsplit(id[1], "[.]")))
  xcoord <- splitid[2]+1
  ycoord <- splitid[1]+1
  data.frame(id[1], original[xcoord,ycoord], stringsAsFactors=FALSE)
}

#make the mean for each triangle and assign to triangle id
tri_means <- function(id, merged_all_refs_tmp2){
  small_matrix <- merged_all_refs_tmp2[which(merged_all_refs_tmp2[,9]==id),]
  tri_mean <- mean(small_matrix[,10])
  data.frame(id, tri_mean)
}

tri_means_just_scale <- function(id, merged_all_refs_tmp2){
  small_matrix <- merged_all_refs_tmp2[which(merged_all_refs_tmp2[,6]==id),]
  tri_mean <- mean(small_matrix[,7])
  data.frame(id, tri_mean)
}

#make the mean for each triangle and assign to triangle id but for triangles that have been rotated
tri_means_rotate <- function(id){
  small_matrix <- merged_all_refs_tmp2[which(merged_all_refs_tmp2[6]==id),]
  tri_mean <- mean(small_matrix[,7])
  data.frame(id, tri_mean)
}

# create a final table that is used for plotting triangles
get_final_table <- function(the_row, merged_all_refs_tmp3){
  data.frame(the_row[3],the_row[2], merged_all_refs_tmp3[which(the_row[1]==merged_all_refs_tmp3$coord_id),11])
}

# create a final table that is used for plotting triangles after rotating
get_final_table_rotate <- function(the_row, merged_all_refs_tmp3){
  data.frame(the_row[4],the_row[3], merged_all_refs_tmp3[which(the_row[1]==merged_all_refs_tmp3$tri_id),11])
}

# Creates a 2 column matrix used for reference
create_matrix <- function(xmin, xmax, ymin, ymax){
  tmp <-lapply(ymin:ymax, function(y) cbind(xmin:xmax, rep(y,length(xmin:xmax))))
  final <- do.call(rbind,tmp)
  return(final)
}

# find the intersect between the lines having given them the coefficients of the linear equation
get_intersect <- function(slope1, intercept1, slope2, intercept2) { 
  x_cor = (intercept2-intercept1)/ifelse(slope2<0, slope1-slope2, slope1+slope2)
  y_cor = (slope1*abs(x_cor))+intercept1
  return(list("x"= x_cor, "y"= y_cor))
}

#Smoothen the values of the max slope ellipse using a running mean on the ellipse split into 4
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

smooth_ellipse_plateu <- function(e, smoothf){
  mean_ellipse_table <- NULL
  mean_ellipse_table <- cbind(mean_ellipse_table, e$fitted_center$plateu[1:360,])
  mean_ellipse_table <- rbind(mean_ellipse_table, e$fitted_center$plateu[541:720,])
  mean_ellipse_table <- rbind(mean_ellipse_table, e$fitted_center$plateu[361:540,])
  mean_ellipse <- NULL
  mean_ellipse <- cbind(mean_ellipse, (mean_ellipse_table[1:360,1]+mean_ellipse_table[361:720,1])/2)
  mean_ellipse <- cbind(mean_ellipse, (mean_ellipse_table[1:360,2]+mean_ellipse_table[361:720,2])/2)
  smoothed_ellipse <- cbind(runmed(mean_ellipse[,1], smoothf), runmed(mean_ellipse[,2], smoothf))
  plot(smoothed_ellipse)
  return(smoothed_ellipse)
}

#find the rotation value assigning the shortest diameter to the y axis
new_find_rotate <- function(smoothed_ellipse){
  rotation_data <- NULL
  rotation_data <- cbind(rotation_data, smoothed_ellipse[1:180,])
  rotation_data <- cbind(rotation_data, smoothed_ellipse[181:360,])
  distance_value <- sqrt(((rotation_data[,1]-rotation_data[,3])^2) + (rotation_data[,2]- rotation_data[,4])^2)
  rotation_data <- cbind(rotation_data, distance_value)
  rotation_data_order <- rotation_data[order(rotation_data[,5]),]
  degree <- match(rotation_data_order[1,], rotation_data)
  degree <- degree[1]
  degree <- degree - 90
}

#get a list of all of the dicom files within a file
get_eye_paths <- function(pathname){
  filenames <- list.files(path=pathname)
  all_eye_paths <- list()
  for (x in filenames){
    pathname1 <- paste(pathname, x, sep="")
    pathname1 <- paste(pathname1, "/OCT/", sep="")
    eye_files <- list.files(pathname1, pattern=".dcm")
    for(y in eye_files){
      full_eye_path <- paste(pathname1, y, sep="")
      all_eye_paths <- append(all_eye_paths, full_eye_path)
    }
  }
  return(all_eye_paths)
}
#get either the first scan or the duplicate by specifying val which, either 1 or 2 which represents whether you want the first or second repear
get_eye_pathsrepeats <- function(pathname, val){
  filenames <- list.files(path=pathname)
  all_eye_paths <- list()
  for (x in filenames){
    pathname1 <- paste(pathname, x, sep="")
    pathname1 <- paste(pathname1, "/OCT/", sep="")
    eye_files <- list.files(pathname1, pattern=".dcm")
    full_eye_path <- paste(pathname1, eye_files[val], sep="")
    all_eye_paths <- append(all_eye_paths, full_eye_path)
  }
  return(all_eye_paths)
}

#collect all of the paramaters including the name, fovea coordiantes, fovea height, left or right scan and sex of paient
get_needed_params <- function(all_eye_paths){
  all_eyes <- NULL
  for (x in all_eye_paths){
    #image_file <- all_eye_paths[[27]]
    image_file <- x
    e = registar(image_file, debug=FALSE)
    ordered_ring <- e$e2[order(e$e2[,1]),]
    fovea_centre_x <- e$zcenter/(dim(e$model$top_band)[1])
    fovea_centre_y <- e$ycenter/(dim(e$model$top_band)[2])
    fovea_coords <- c(fovea_centre_x, fovea_centre_y)
    fovea_height <- abs((ordered_ring[1,1]) - (ordered_ring[dim(ordered_ring)[1],1]))
    fovea_height <- fovea_height/dim(e$centered_band)[1]
    eye_name <- strsplit(x, "[/]")
    eye_name <- eye_name[length(eye_name)]
    info <- dicomInfo(image_file)
    lorr <- info$hdr$value[info$hdr$name=="ImageLaterality"]
    sex <- info$hdr$value[info$hdr$name=="PatientsSex"]
    eye_data <- NULL
    eye_data <- cbind(x, eye_name, fovea_centre_x,fovea_centre_y, fovea_height, lorr, sex)
    all_eyes <- rbind( all_eyes, eye_data)
  }
  return(all_eyes)
}

#creates a reference matrix that can be scaled based on the size of the fovea and centred so that the fovea is in the middle of the grid, this can then be fed to rotate matrix
scale_and_centre_rotate <- function(amodel, fovea_height, fovea_coords, n){
  fovea_height <- fovea_height/dim(amodel)[1]
  fovea_coords[1] <- fovea_coords[1]/dim(amodel)[2]
  fovea_coords[2] <- fovea_coords[2]/dim(amodel)[1]
  tmp_mat <- create_matrix(0, ((dim(amodel)[2])-1), 0, ((dim(amodel)[1])-1))
  rownames(tmp_mat) <- paste(tmp_mat[,1], ".", tmp_mat[,2], sep="")
  id_mapping <- data.frame(id=paste(tmp_mat[,1], ".", tmp_mat[,2], sep=""), x=tmp_mat[,1], y=tmp_mat[,2], stringsAsFactors = FALSE)
  d1 <- dim(amodel)[1]*(1/3)
  fov_height_scale <- fovea_height*dim(amodel)[1]
  fovea_coords_big <- c(fovea_coords[1]*dim(amodel)[2], fovea_coords[2]*dim(amodel)[1])
  scale_co <- d1/fov_height_scale
  fovea_coords_scale <- fovea_coords_big*scale_co
  tmp_mat_scale <- tmp_mat*scale_co
  centre_mat <- c((dim(amodel)[1]*scale_co)/2, (dim(amodel)[2]*scale_co)/2)
  foveal_difference <- c(centre_mat[1]-((dim(amodel)[1])/2), centre_mat[2]-((dim(amodel)[2])/2))
  scale_and_shift <- tmp_mat_scale
  scale_and_shift[,1] <- scale_and_shift[,1] - foveal_difference[2]
  scale_and_shift[,2] <- scale_and_shift[,2] - foveal_difference[1]
  scale_shift_small <- scale_and_shift
  #scale_shift_small[,1] <- scale_shift_small[,1]/dim(amodel)[2]
  #scale_shift_small[,2] <- scale_shift_small[,2]/dim(amodel)[1]
  scale_shift_small <- scale_shift_small[complete.cases(scale_shift_small),]
  return(list(ref_table = scale_shift_small, id_mapping=id_mapping))
  
}

#technique to rotate the matrix 
#rotate_matrix <- function(amodel, rotate_value){
rotate_matrix <- function(tst_matrix, rotate_value, amodel){
  #tst_matrix <- create_matrix(0, (dim(amodel)[2])-1, 0, (dim(amodel)[1])-1)
  #rownames(tst_matrix) <- paste( tst_matrix[,1], tst_matrix[,2], sep=".")
  id_mapping <- data.frame(id=paste(tst_matrix[,1], ".", tst_matrix[,2], sep=""), x=tst_matrix[,1], y=tst_matrix[,2], stringsAsFactors = FALSE)
  tst_matrix[,1] <- tst_matrix[,1]-(dim(amodel)[2])%/%2
  tst_matrix[,2] <- tst_matrix[,2]-(dim(amodel)[1])%/%2
  rotate_value_clock <- rotate_value
  rotate_value_clock_rad <- deg2rad(rotate_value_clock)
  rotation_matrix_tst <- matrix(c(0,1,-1,0), nrow=2)
  rotation_matrix <- matrix(c(cos(rotate_value_clock_rad), sin(rotate_value_clock_rad), -sin(rotate_value_clock_rad), cos(rotate_value_clock_rad)), nrow=2, ncol=2)
  rotated_amodel <- tst_matrix %*% rotation_matrix
  rotated_amodel[,1] <- rotated_amodel[,1]+(dim(amodel)[2])%/%2
  rotated_amodel[,2] <- rotated_amodel[,2]+(dim(amodel)[1])%/%2
  rotated_amodel[,1] <- rotated_amodel[,1]/(dim(amodel)[2])
  rotated_amodel[,2] <- rotated_amodel[,2]/(dim(amodel)[1])
  #rotated_amodeldf <- as.data.frame(rotated_amodel)
  #rotate_and_id <- merge(rotated_amodeldf, ref_ids, by=0, all=TRUE, stringsAsFactors=FALSE )
  return(list(rotated_amodel=rotated_amodel, id_mapping=id_mapping))
}

#Scale, centre and tesselate an image without rotation
tesselate_and_plot <- function(image_file, centre_f, fovea_height, fovea_coords, n, triangle_num = 5){
  model = run_generate_model_and_find_fovea_dip_poisition(image_file)
  e <- registar(image_file, debug=FALSE)
  
  amodel = center_align(model, center_f)
  scale_centre_ids <- scale_and_centre(amodel, fovea_height, fovea_coords, n)
  
  indicies <- make_tri_grid(amodel,triangle_num)
  indicies[,2]<- round(as.numeric(indicies[,2]),3)
  sqr_inds <- remove_to_square(indicies)
  verticesdf <- get_triangle_vertices(sqr_inds,amodel)
  all_refs <- get_triangle_ref_coords(scale_centre_ids$ref_table, verticesdf)
  all_refs <- all_refs[complete.cases(all_refs),]
  merged_all_refs <- merge(scale_centre_ids, all_refs, by=0, all=TRUE, stringsAsFactors=FALSE)
  
  test <- do.call(rbind, apply(merged_all_refs, 1, Ids2coordinates, original=amodel))
  merged_all_refs_tmp2 <- merge(merged_all_refs, test, by=1, all=TRUE, stringsAsFactors=FALSE)
  colnames(merged_all_refs_tmp2) <- c("coord_id", "x", "y","id_mapping", "x1", "y1", "x2", "x3", "tri_id", "value")
  
  tri_ids <- unique(merged_all_refs_tmp2$tri_id)
  all_tri_means <- do.call(rbind, lapply(tri_ids, tri_means))
  colnames(all_tri_means) <- c("tri_id", "mean")
  
  merged_all_refs_tmp3 <- merge(merged_all_refs_tmp2, all_tri_means, by="tri_id", all=TRUE, stringsAsFactors=FALSE)
  tri_id_table <- data.frame(merged_all_refs_tmp3$y1, merged_all_refs_tmp3$x1, merged_all_refs_tmp3$tri_id)
  end_info <- data.frame(merged_all_refs_tmp3$y1, merged_all_refs_tmp3$x1, merged_all_refs_tmp3$mean)
  verticesdf2 <- verticesdf
  verticesdf2[,1] <- verticesdf2[,1]*dim(amodel)[2]
  verticesdf2[,2] <- verticesdf2[,2]*dim(amodel)[1]
  
  p <- ggplot(data=end_info)
  p + 
    geom_point(aes(x=merged_all_refs_tmp3.x1, y=merged_all_refs_tmp3.y1, col=merged_all_refs_tmp3.mean), na.rm=TRUE) + 
    geom_point(data=as.data.frame(verticesdf2), aes(x=V1, y=V2), color='red')  +
    theme_bw()
}

#Scale, centres, rotates and tesselates an image ready for lapply across angles
rotate_and_plot_byangle <- function(rotate_value, image_file, centre_f, triangle_num=5, fovea_height, fovea_coords, n){
  model <- run_generate_model_and_find_fovea_dip_poisition(image_file)
  e <- registar(image_file, debug=FALSE)
  amodel <- center_align(model, center_f)
  smoothed_ellipse <- smooth_ellipse_max_slope(e, 31)
  #rotate_value <- new_find_rotate(smoothed_ellipse)
  scale_centre_ids <- scale_and_centre_rotate(amodel, fovea_height, fovea_coords, n)
  #rotated_matrix <- rotate_matrix(amodel, rotate_value)
  #rotate_value <- 90
  rotated_matrix <- rotate_matrix(scale_centre_ids$ref_table, rotate_value)
  
  indicies <- make_tri_grid(amodel,triangle_num)
  indicies[,2]<- round(as.numeric(indicies[,2]),3)
  sqr_inds <- remove_to_square(indicies)
  verticesdf <- get_triangle_vertices(sqr_inds,amodel)
  all_refs <- get_triangle_ref_coords(rotated_matrix$rotated_amodel, verticesdf)
  all_refs <- all_refs[complete.cases(all_refs),]
  merged_all_refs <- merge(rotated_matrix, all_refs, by=0, all=TRUE, stringsAsFactors=FALSE)
  
  test <- do.call(rbind, apply(merged_all_refs, 1, Ids2coordinates_rotate, original=amodel))
  merged_all_refs_tmp2 <- merge(merged_all_refs, test, by=1, all=TRUE, stringsAsFactors=FALSE)
  colnames(merged_all_refs_tmp2) <- c("coord_id", "x", "y","id_mapping", "x1", "y1", "x2", "x3", "tri_id", "value")
  
  tri_ids <- unique(merged_all_refs_tmp2$tri_id)
  all_tri_means <- do.call(rbind, lapply(tri_ids, tri_means))
  colnames(all_tri_means) <- c("tri_id", "mean")
  
  merged_all_refs_tmp3 <- merge(merged_all_refs_tmp2, all_tri_means, by="tri_id", all=TRUE, stringsAsFactors=FALSE)
  end_info2 <-  data.frame(merged_all_refs_tmp3$y, merged_all_refs_tmp3$x, merged_all_refs_tmp3$mean)
  
  p <- ggplot(data=end_info2)
  p + 
    geom_point(aes(x=merged_all_refs_tmp3.x, y=merged_all_refs_tmp3.y, col=merged_all_refs_tmp3.mean), na.rm=TRUE) + 
    geom_point(data=as.data.frame(verticesdf), aes(x=V1, y=V2), color='red')  +
    theme_bw()
  return(list(end_info2=end_info2, verticesdf=verticesdf, all_refs= all_refs, p=p))
}

create_rotate_and_scale_table <- function(x, param_table, centre_f){
  image_file <- param_table[[x,1]]
  fovea_height <- param_table[[x,5]]
  model <- run_generate_model_and_find_fovea_dip_poisition(image_file)
  e <- registar(image_file, debug=FALSE)
  amodel <- center_align(model, centre_f)
  smoothed_ellipse <- smooth_ellipse_max_slope(e, 31)
  rotate_val <- new_find_rotate(smoothed_ellipse)
  d1 <- dim(amodel)[1]*(1/3)
  fovea_height_scale <- fovea_height*dim(amodel)[1]
  scale_co <- d1/fovea_height_scale
  return(c(image_file=image_file, rotate_val=rotate_val, scale_co=scale_co))
}


#Scale, centres, rotates and tesselates an image ready for lapply across images
rotate_and_plot <- function(image_file, center_f, triangle_num=5, fovea_height, fovea_coords, n){
  model <- run_generate_model_and_find_fovea_dip_poisition(image_file)
  e <- registar(image_file, debug=FALSE)
  amodel <- center_align(model, center_f)
  #################
  #image_file <- all_params_odd[[x,1]]
  #fovea_height <- all_params_odd[[x,5]]
  fovea_height_scale <- fovea_height*(dim(amodel)[1])
  #fovea_coordx <- all_params_odd[[x,3]]
  fovea_coordx_scale <- fovea_coords[2]*(dim(amodel)[2])
  #fovea_coordy <- all_params_odd[[x,4]]
  fovea_coordy_scale <- fovea_coords[1]*(dim(amodel)[1])
  fovea_coords <- c(fovea_coordy_scale, fovea_coordx_scale)
  #################
  smoothed_ellipse <- smooth_ellipse_max_slope(e, 31)
  rotate_value <- new_find_rotate(smoothed_ellipse)
  scale_centre_ids <- scale_and_centre_rotate(amodel, fovea_height_scale, fovea_coords, n)
  #rotated_matrix <- rotate_matrix(amodel, rotate_value)
  #rotate_value <- 90
  rotated_matrix <- rotate_matrix(scale_centre_ids$ref_table, rotate_value, amodel)
  
  indicies <- make_tri_grid(amodel,triangle_num)
  indicies[,2]<- round(as.numeric(indicies[,2]),3)
  sqr_inds <- remove_to_square(indicies)
  verticesdf <- get_triangle_vertices(sqr_inds,amodel)
  all_refs <- get_triangle_ref_coords(rotated_matrix$rotated_amodel, verticesdf)
  all_refs <- all_refs[complete.cases(all_refs),]
  merged_all_refs <- merge(rotated_matrix, all_refs, by=0, all=TRUE, stringsAsFactors=FALSE)
  
  test <- do.call(rbind, apply(merged_all_refs, 1, Ids2coordinates_rotate, original=amodel))
  merged_all_refs_tmp2 <- merge(merged_all_refs, test, by=1, all=TRUE, stringsAsFactors=FALSE)
  colnames(merged_all_refs_tmp2) <- c("coord_id", "x", "y","id_mapping", "x1", "y1", "x2", "x3", "tri_id", "value")
  
  tri_ids <- unique(merged_all_refs_tmp2$tri_id)
  all_tri_means <- do.call(rbind, lapply(tri_ids, tri_means, merged_all_refs_tmp2=merged_all_refs_tmp2))
  colnames(all_tri_means) <- c("tri_id", "mean")
  
  merged_all_refs_tmp3 <- merge(merged_all_refs_tmp2, all_tri_means, by="tri_id", all=TRUE, stringsAsFactors=FALSE)
  end_info2 <-  data.frame(merged_all_refs_tmp3$y, merged_all_refs_tmp3$x, merged_all_refs_tmp3$mean)
  
  p <- ggplot(data=end_info2)
  p + 
    geom_point(aes(x=merged_all_refs_tmp3.x, y=merged_all_refs_tmp3.y, col=merged_all_refs_tmp3.mean), na.rm=TRUE) + 
    geom_point(data=as.data.frame(verticesdf), aes(x=V1, y=V2), color='red')  +
    theme_bw()
  return(list(end_info2=end_info2, verticesdf=verticesdf, all_refs= all_refs, p=p))
}

#Scale, centres, rotates and tesselates an image ready for lapply across images
rotate_and_plot_mean_rotate <- function(image_file, center_f, triangle_num=5, fovea_height, fovea_coords, n){
  model <- run_generate_model_and_find_fovea_dip_poisition(image_file)
  e <- registar(image_file, debug=FALSE)
  amodel <- center_align(model, center_f)
  #################
  #image_file <- all_params_odd[[x,1]]
  #fovea_height <- all_params_odd[[x,5]]
  fovea_height_scale <- fovea_height*(dim(amodel)[1])
  #fovea_coordx <- all_params_odd[[x,3]]
  fovea_coordx_scale <- fovea_coords[2]*(dim(amodel)[2])
  #fovea_coordy <- all_params_odd[[x,4]]
  fovea_coordy_scale <- fovea_coords[1]*(dim(amodel)[1])
  fovea_coords <- c(fovea_coordy_scale, fovea_coordx_scale)
  #################
  smoothed_ellipse1 <- smooth_ellipse_max_slope(e, 31)
  smoothed_ellipse2 <- smooth_ellipse_plateu(e,31)
  rotate_value1 <- new_find_rotate(smoothed_ellipse1)
  rotate_value2 <- new_find_rotate(smoothed_ellipse2)
  rotate_value <- mean(rotate_value1, rotate_value2)
  scale_centre_ids <- scale_and_centre_rotate(amodel, fovea_height_scale, fovea_coords, n)
  #rotated_matrix <- rotate_matrix(amodel, rotate_value)
  #rotate_value <- 90
  rotated_matrix <- rotate_matrix(scale_centre_ids$ref_table, rotate_value, amodel)
  
  indicies <- make_tri_grid(amodel,triangle_num)
  indicies[,2]<- round(as.numeric(indicies[,2]),3)
  sqr_inds <- remove_to_square(indicies)
  verticesdf <- get_triangle_vertices(sqr_inds,amodel)
  all_refs <- get_triangle_ref_coords(rotated_matrix$rotated_amodel, verticesdf)
  all_refs <- all_refs[complete.cases(all_refs),]
  merged_all_refs <- merge(rotated_matrix, all_refs, by=0, all=TRUE, stringsAsFactors=FALSE)
  
  test <- do.call(rbind, apply(merged_all_refs, 1, Ids2coordinates_rotate, original=amodel))
  merged_all_refs_tmp2 <- merge(merged_all_refs, test, by=1, all=TRUE, stringsAsFactors=FALSE)
  colnames(merged_all_refs_tmp2) <- c("coord_id", "x", "y","id_mapping", "x1", "y1", "x2", "x3", "tri_id", "value")
  
  tri_ids <- unique(merged_all_refs_tmp2$tri_id)
  all_tri_means <- do.call(rbind, lapply(tri_ids, tri_means, merged_all_refs_tmp2=merged_all_refs_tmp2))
  colnames(all_tri_means) <- c("tri_id", "mean")
  
  merged_all_refs_tmp3 <- merge(merged_all_refs_tmp2, all_tri_means, by="tri_id", all=TRUE, stringsAsFactors=FALSE)
  end_info2 <-  data.frame(merged_all_refs_tmp3$y, merged_all_refs_tmp3$x, merged_all_refs_tmp3$mean)
  
  p <- ggplot(data=end_info2)
  p + 
    geom_point(aes(x=merged_all_refs_tmp3.x, y=merged_all_refs_tmp3.y, col=merged_all_refs_tmp3.mean), na.rm=TRUE) + 
    geom_point(data=as.data.frame(verticesdf), aes(x=V1, y=V2), color='red')  +
    theme_bw()
  return(list(end_info2=end_info2, verticesdf=verticesdf, all_refs= all_refs, p=p))
}

rotate_and_plot_no_scale <-function(image_file, center_f, triangle_num=5, fovea_height, fovea_coords, n){
  model <- run_generate_model_and_find_fovea_dip_poisition(image_file)
  e <- registar(image_file, debug=FALSE)
  amodel <- center_align(model, center_f)
  #################
  #image_file <- all_params_odd[[x,1]]
  #fovea_height <- all_params_odd[[x,5]]
  fovea_height_scale <- fovea_height*(dim(amodel)[1])
  #fovea_coordx <- all_params_odd[[x,3]]
  fovea_coordx_scale <- fovea_coords[2]*(dim(amodel)[2])
  #fovea_coordy <- all_params_odd[[x,4]]
  fovea_coordy_scale <- fovea_coords[1]*(dim(amodel)[1])
  fovea_coords <- c(fovea_coordy_scale, fovea_coordx_scale)
  #################
  tmp_mat <- create_matrix(0, ((dim(amodel)[2])-1), 0, ((dim(amodel)[1])-1))
  rownames(tmp_mat) <- paste(tmp_mat[,1], ".", tmp_mat[,2], sep="")
  smoothed_ellipse <- smooth_ellipse_max_slope(e, 31)
  rotate_value <- new_find_rotate(smoothed_ellipse)
  #scale_centre_ids <- scale_and_centre_rotate(amodel, fovea_height_scale, fovea_coords, n)
  #rotated_matrix <- rotate_matrix(amodel, rotate_value)
  #rotate_value <- 90
  rotated_matrix <- rotate_matrix(tmp_mat, rotate_value, amodel)
  
  
  indicies <- make_tri_grid(amodel,triangle_num)
  indicies[,2]<- round(as.numeric(indicies[,2]),3)
  sqr_inds <- remove_to_square(indicies)
  verticesdf <- get_triangle_vertices(sqr_inds,amodel)
  all_refs <- get_triangle_ref_coords(rotated_matrix$rotated_amodel, verticesdf)
  all_refs <- all_refs[complete.cases(all_refs),]
  merged_all_refs <- merge(rotated_matrix, all_refs, by=0, all=TRUE, stringsAsFactors=FALSE)
  
  test <- do.call(rbind, apply(merged_all_refs, 1, Ids2coordinates_rotate, original=amodel))
  merged_all_refs_tmp2 <- merge(merged_all_refs, test, by=1, all=TRUE, stringsAsFactors=FALSE)
  colnames(merged_all_refs_tmp2) <- c("coord_id", "x", "y","id_mapping", "x1", "y1", "x2", "x3", "tri_id", "value")
  
  tri_ids <- unique(merged_all_refs_tmp2$tri_id)
  all_tri_means <- do.call(rbind, lapply(tri_ids, tri_means, merged_all_refs_tmp2=merged_all_refs_tmp2))
  colnames(all_tri_means) <- c("tri_id", "mean")
  
  merged_all_refs_tmp3 <- merge(merged_all_refs_tmp2, all_tri_means, by="tri_id", all=TRUE, stringsAsFactors=FALSE)
  end_info2 <-  data.frame(merged_all_refs_tmp3$y, merged_all_refs_tmp3$x, merged_all_refs_tmp3$mean)
  
  p <- ggplot(data=end_info2)
  p + 
    geom_point(aes(x=merged_all_refs_tmp3.x, y=merged_all_refs_tmp3.y, col=merged_all_refs_tmp3.mean), na.rm=TRUE) + 
    geom_point(data=as.data.frame(verticesdf), aes(x=V1, y=V2), color='red')  +
    theme_bw()
  return(list(end_info2=end_info2, verticesdf=verticesdf, all_refs= all_refs, p=p))
}

rotate_and_plot_no_scale_mean_rotate <-function(image_file, center_f, triangle_num=5, fovea_height, fovea_coords, n){
  model <- run_generate_model_and_find_fovea_dip_poisition(image_file)
  e <- registar(image_file, debug=FALSE)
  amodel <- center_align(model, center_f)
  #################
  #image_file <- all_params_odd[[x,1]]
  #fovea_height <- all_params_odd[[x,5]]
  fovea_height_scale <- fovea_height*(dim(amodel)[1])
  #fovea_coordx <- all_params_odd[[x,3]]
  fovea_coordx_scale <- fovea_coords[2]*(dim(amodel)[2])
  #fovea_coordy <- all_params_odd[[x,4]]
  fovea_coordy_scale <- fovea_coords[1]*(dim(amodel)[1])
  fovea_coords <- c(fovea_coordy_scale, fovea_coordx_scale)
  #################
  tmp_mat <- create_matrix(0, ((dim(amodel)[2])-1), 0, ((dim(amodel)[1])-1))
  rownames(tmp_mat) <- paste(tmp_mat[,1], ".", tmp_mat[,2], sep="")
  smoothed_ellipse1 <- smooth_ellipse_max_slope(e, 31)
  smoothed_ellipse2 <- smooth_ellipse_plateu(e, 31)
  rotate_value1 <- new_find_rotate(smoothed_ellipse1)
  rotate_value2 <- new_find_rotate(smoothed_ellipse2)
  rotate_value <- mean(rotate_value1, rotate_value2)
  #scale_centre_ids <- scale_and_centre_rotate(amodel, fovea_height_scale, fovea_coords, n)
  #rotated_matrix <- rotate_matrix(amodel, rotate_value)
  #rotate_value <- 90
  rotated_matrix <- rotate_matrix(tmp_mat, rotate_value, amodel)
  
  
  indicies <- make_tri_grid(amodel,triangle_num)
  indicies[,2]<- round(as.numeric(indicies[,2]),3)
  sqr_inds <- remove_to_square(indicies)
  verticesdf <- get_triangle_vertices(sqr_inds,amodel)
  all_refs <- get_triangle_ref_coords(rotated_matrix$rotated_amodel, verticesdf)
  all_refs <- all_refs[complete.cases(all_refs),]
  merged_all_refs <- merge(rotated_matrix, all_refs, by=0, all=TRUE, stringsAsFactors=FALSE)
  
  test <- do.call(rbind, apply(merged_all_refs, 1, Ids2coordinates_rotate, original=amodel))
  merged_all_refs_tmp2 <- merge(merged_all_refs, test, by=1, all=TRUE, stringsAsFactors=FALSE)
  colnames(merged_all_refs_tmp2) <- c("coord_id", "x", "y","id_mapping", "x1", "y1", "x2", "x3", "tri_id", "value")
  
  tri_ids <- unique(merged_all_refs_tmp2$tri_id)
  all_tri_means <- do.call(rbind, lapply(tri_ids, tri_means, merged_all_refs_tmp2=merged_all_refs_tmp2))
  colnames(all_tri_means) <- c("tri_id", "mean")
  
  merged_all_refs_tmp3 <- merge(merged_all_refs_tmp2, all_tri_means, by="tri_id", all=TRUE, stringsAsFactors=FALSE)
  end_info2 <-  data.frame(merged_all_refs_tmp3$y, merged_all_refs_tmp3$x, merged_all_refs_tmp3$mean)
  
  p <- ggplot(data=end_info2)
  p + 
    geom_point(aes(x=merged_all_refs_tmp3.x, y=merged_all_refs_tmp3.y, col=merged_all_refs_tmp3.mean), na.rm=TRUE) + 
    geom_point(data=as.data.frame(verticesdf), aes(x=V1, y=V2), color='red')  +
    theme_bw()
  return(list(end_info2=end_info2, verticesdf=verticesdf, all_refs= all_refs, p=p))
}


just_scale_and_plot <- function(image_file, center_f, triangle_num=5, fovea_height, fovea_coords, n){
  model <- run_generate_model_and_find_fovea_dip_poisition(image_file)
  e <- registar(image_file, debug=FALSE)
  amodel <- center_align(model, center_f)
  #################
  #image_file <- all_params_odd[[x,1]]
  #fovea_height <- all_params_odd[[x,5]]
  fovea_height_scale <- fovea_height*(dim(amodel)[1])
  #fovea_coordx <- all_params_odd[[x,3]]
  fovea_coordx_scale <- fovea_coords[2]*(dim(amodel)[2])
  #fovea_coordy <- all_params_odd[[x,4]]
  fovea_coordy_scale <- fovea_coords[1]*(dim(amodel)[1])
  fovea_coords <- c(fovea_coordy_scale, fovea_coordx_scale)
  #################
  smoothed_ellipse <- smooth_ellipse_max_slope(e, 31)
  #rotate_value <- new_find_rotate(smoothed_ellipse)
  scale_centre_ids <- scale_and_centre_rotate(amodel, fovea_height_scale, fovea_coords, n)
  scale_centre_ids$ref_table[,1] <- scale_centre_ids$ref_table[,1]/(dim(amodel)[2])
  scale_centre_ids$ref_table[,2] <- scale_centre_ids$ref_table[,2]/(dim(amodel)[1])
  #rotated_matrix <- rotate_matrix(amodel, rotate_value)
  #rotate_value <- 90
  #rotated_matrix <- rotate_matrix(scale_centre_ids$ref_table, rotate_value, amodel)
  
  indicies <- make_tri_grid(amodel,triangle_num)
  indicies[,2]<- round(as.numeric(indicies[,2]),3)
  sqr_inds <- remove_to_square(indicies)
  verticesdf <- get_triangle_vertices(sqr_inds,amodel)
  all_refs <- get_triangle_ref_coords(scale_centre_ids$ref_table, verticesdf)
  all_refs <- all_refs[complete.cases(all_refs),]
  merged_all_refs <- merge(scale_centre_ids$ref_table, all_refs, by=0, all=TRUE, stringsAsFactors=FALSE)
  
  test <- do.call(rbind, apply(merged_all_refs, 1, Ids2coordinates_rotate, original=amodel))
  merged_all_refs_tmp2 <- merge(merged_all_refs, test, by=1, all=TRUE, stringsAsFactors=FALSE)
  colnames(merged_all_refs_tmp2) <- c("coord_id", "x", "y", "x1", "y1", "tri_id", "value")
  
  tri_ids <- unique(merged_all_refs_tmp2$tri_id)
  all_tri_means <- do.call(rbind, lapply(tri_ids, tri_means_just_scale, merged_all_refs_tmp2=merged_all_refs_tmp2))
  colnames(all_tri_means) <- c("tri_id", "mean")
  
  merged_all_refs_tmp3 <- merge(merged_all_refs_tmp2, all_tri_means, by="tri_id", all=TRUE, stringsAsFactors=FALSE)
  end_info2 <-  data.frame(merged_all_refs_tmp3$y, merged_all_refs_tmp3$x, merged_all_refs_tmp3$mean)
  
  p <- ggplot(data=end_info2)
  p + 
    geom_point(aes(x=merged_all_refs_tmp3.x, y=merged_all_refs_tmp3.y, col=merged_all_refs_tmp3.mean), na.rm=TRUE) + 
    geom_point(data=as.data.frame(verticesdf), aes(x=V1, y=V2), color='red')  +
    theme_bw()
  return(list(end_info2=end_info2, verticesdf=verticesdf, all_refs= all_refs, p=p))
}

just_tessalte_and_plot <- function(image_file, center_f, triangle_num=5, fovea_height, fovea_coords, n){
  model <- run_generate_model_and_find_fovea_dip_poisition(image_file)
  e <- registar(image_file, debug=FALSE)
  amodel <- center_align(model, center_f)
  #################
  #image_file <- all_params_odd[[x,1]]
  #fovea_height <- all_params_odd[[x,5]]
  fovea_height_scale <- fovea_height*(dim(amodel)[1])
  #fovea_coordx <- all_params_odd[[x,3]]
  fovea_coordx_scale <- fovea_coords[2]*(dim(amodel)[2])
  #fovea_coordy <- all_params_odd[[x,4]]
  fovea_coordy_scale <- fovea_coords[1]*(dim(amodel)[1])
  fovea_coords <- c(fovea_coordy_scale, fovea_coordx_scale)
  #################
  #smoothed_ellipse <- smooth_ellipse_max_slope(e, 31)
  #rotate_value <- new_find_rotate(smoothed_ellipse)
  #scale_centre_ids <- scale_and_centre_rotate(amodel, fovea_height_scale, fovea_coords, n)
  #rotated_matrix <- rotate_matrix(amodel, rotate_value)
  #rotate_value <- 90
  #rotated_matrix <- rotate_matrix(scale_centre_ids$ref_table, rotate_value, amodel)
  
  tmp_mat <- create_matrix(0, ((dim(amodel)[2])-1), 0, ((dim(amodel)[1])-1))
  rownames(tmp_mat) <- paste(tmp_mat[,1], ".", tmp_mat[,2], sep="")
  tmp_mat[,1] <- tmp_mat[,1]/dim(amodel)[2]
  tmp_mat[,2] <- tmp_mat[,2]/dim(amodel)[1]
  
  indicies <- make_tri_grid(amodel,triangle_num)
  indicies[,2]<- round(as.numeric(indicies[,2]),3)
  sqr_inds <- remove_to_square(indicies)
  verticesdf <- get_triangle_vertices(sqr_inds,amodel)
  all_refs <- get_triangle_ref_coords(tmp_mat, verticesdf)
  all_refs <- all_refs[complete.cases(all_refs),]
  merged_all_refs <- merge(tmp_mat, all_refs, by=0, all=TRUE, stringsAsFactors=FALSE)
  
  test <- do.call(rbind, apply(merged_all_refs, 1, Ids2coordinates_rotate, original=amodel))
  merged_all_refs_tmp2 <- merge(merged_all_refs, test, by=1, all=TRUE, stringsAsFactors=FALSE)
  colnames(merged_all_refs_tmp2) <- c("coord_id", "x", "y", "x2", "x3", "tri_id", "value")
  
  tri_ids <- unique(merged_all_refs_tmp2$tri_id)
  all_tri_means <- do.call(rbind, lapply(tri_ids, tri_means_just_scale, merged_all_refs_tmp2=merged_all_refs_tmp2))
  colnames(all_tri_means) <- c("tri_id", "mean")
  
  merged_all_refs_tmp3 <- merge(merged_all_refs_tmp2, all_tri_means, by="tri_id", all=TRUE, stringsAsFactors=FALSE)
  end_info2 <-  data.frame(merged_all_refs_tmp3$y, merged_all_refs_tmp3$x, merged_all_refs_tmp3$mean)
  
  p <- ggplot(data=end_info2)
  p + 
    geom_point(aes(x=merged_all_refs_tmp3.x, y=merged_all_refs_tmp3.y, col=merged_all_refs_tmp3.mean), na.rm=TRUE) + 
    geom_point(data=as.data.frame(verticesdf), aes(x=V1, y=V2), color='red')  +
    theme_bw()
  return(list(end_info2=end_info2, verticesdf=verticesdf, all_refs= all_refs, p=p))
}

add_standard_triangle <- function(end_info2, verticesdf){
  triangle_groups <- seq(from=1, to=dim(verticesdf)[1], by=3)
  matrix <- as.matrix(end_info2[,1:2])
  test_mctest <- lapply(triangle_groups, function(x){
    one_triangle_bnd <- verticesdf[x:(x+2),]
    inouts <- in.out(one_triangle_bnd, matrix)
    TRUTHS <- which(inouts)
    triangle_friends = NULL
    if(length(TRUTHS)>0) {
      triangle_friends <- end_info2[TRUTHS,]
      triangle_friends$tri_number <- (x-1)/3
    }
    return(triangle_friends)
  })
  test_mctest2 <- do.call(rbind, test_mctest)
  return(test_mctest2)
}

get_reduced_table <- function(x){
  tris <- unique(x$tri_number)
  table1 <- x
  per_tri <- lapply(tris, function(x){
    one_t <- filter(table1, table1$tri_number==x)
    val <- mean(one_t$merged_all_refs_tmp3.mean, na.rm=TRUE)
    row <- c(x,val)
    return(row)
  })
  per_tri_bound <- do.call(rbind, per_tri)
  return(per_tri_bound)
}

#Get the reduced table for a specified triangle iD
get_reduced_table_onetri <- function(x, tri_num){
  tris<- c(tri_num)
  table1 <- x
  per_tri <- lapply(tris,function(x){
    one_t <- filter(table1, table1$tri_number==x)
    val <- mean(one_t$merged_all_refs_tmp3.mean, na.rm=TRUE)
    row <- c(x,val)
    return(row)
  })
  per_tri_bound <- do.call(rbind, per_tri)
  return(per_tri_bound)
}

#Compare two different processed sets of image data
comp_per_tri <- function(x, comp1, comp2){
  evens <- comp1[[x]]
  odds <- comp2[[x]]
  comp_table <- data.frame(x,evens[,3], evens[,2], odds[,2])
  return(comp_table)  
}

median_normalise_intensity <- function(x, datax){
  tmp3 <- datax[[x]]$end_info2
  medi <- median(tmp3[,3], na.rm=TRUE)
  tmp3[,3] <- tmp3[,3]-medi
  return(tmp3)
}

add_standard_triangle_all <- function(x, datax_norm, datax){
  end_info <- datax_norm[[x]]
  verts <- datax[[x]]$verticesdf
  plus_tris <- add_standard_triangle(end_info, verts)
  return(plus_tris)
}

get_whole_reduced_table <- function(x, stdtritab){
  tmp <- lapply(stdtritab, get_reduced_table_onetri, tri_num=x)
  tmp <- do.call(rbind, tmp)
  tmp <- cbind(tmp, 1:length(stdtritab))
  return(tmp)
}

get_correls <- function(x, comp_table){
  table <- comp_table[[x]]
  correl <- cor(table[,3], table[,4], na.rm=TRUE)
  output <- data.frame(x, correl)
  return(output)
}

get_triplot_table <- function(x, ref_points, corr_tab){
  num_pos <- which(ref_points[,3]==x)
  small_tab <- ref_points[num_pos,]
  small_tab <- cbind(small_tab, corr_tab[x,2])
  return(small_tab)
}

making_mean_lorr <- function(x, left_index, right_index, odd_redtab, even_redtab){
  odd_tab <- odd_redtab[[x]]
  even_tab <- even_redtab[[x]]
  odd_subset_left <- odd_tab[left_index,]
  even_subset_left <- even_tab[left_index,]
  lefts_joined <- rbind(odd_subset_left, even_subset_left)
  lefts_mean <- colMeans(lefts_joined)
  odd_subset_right <- odd_tab[right_index,]
  even_subset_right <- even_tab[right_index,]
  rights_joined <- rbind(odd_subset_right, even_subset_right)
  rights_mean <- colMeans(rights_joined)
  return(list(lefts_mean = lefts_mean, rights_mean = rights_mean))
}

###################################
###############DATA################
###################################

image_file = "/Users/currant/OCTeyes/OCT_30_replicates_20161024/1145517/OCT/292422.dcm"
image_file = all_eye_paths[[1]]
center_f = 0.8
fovea_coords<- c(50,200)
fovea_height <- 35
smooth_f <- 30
pathname <- "/Users/currant/OCTeyes/OCT_30_replicates_20161024/"
n=20

###########DATA_PRODUCED BY THIS SCRIPT AND SAVED############
#saveRDS(all_params_even, paste("/Users/currant/OCTeyes/all_params_even",  strftime(Sys.time(), "%Y%m%d"), sep=""))
all_params_even <- readRDS("/Users/currant/OCTeyes/R_objects/all_params_even20170425")

#saveRDS(all_params_odd, paste("/Users/currant/OCTeyes/all_params_odd",  strftime(Sys.time(), "%Y%m%d"), sep=""))
all_params_odd <- readRDS("/Users/currant/OCTeyes/R_objects/all_params_odd20170425")

#saveRDS(biggy_table, paste("/Users/currant/OCTeyes/biggy_table",  strftime(Sys.time(), "%Y%m%d"), sep=""))
biggy_table <- readRDS("/Users/currant/OCTeyes/R_objects/biggy_table20170425")

#saveRDS(turning_outcome, paste("/Users/currant/OCTeyes/turning_outcome",  strftime(Sys.time(), "%Y%m%d"), sep=""))
turning_outcome <- readRDS("/Users/currant/OCTeyes/R_objects/turning_outcome20170425")

###################################
#############ANALYSIS##############
###################################

##########GENERATE PARAMETER TABLES FOR ~60 EYE SAMPLES###########
all_eye_paths <- get_eye_paths(pathname)
all_eye_params10 <- get_needed_params(all_eye_paths[1:10])
all_eye_paths_odd <- get_eye_pathsrepeats(pathname, 1)
all_params_odd <- get_needed_params(all_eye_paths_odd)
all_eye_paths_even <- get_eye_pathsrepeats(pathname, 2)
all_params_even <- get_needed_params(all_eye_paths_even)


##########REMOVE THOS THAT DO NOT WORK - REWORK WITH QUALITY CONTROL#######
all_params_odd_minus <- all_params_odd[-c(6),]
all_params_odd_minus <- all_params_odd_minus[-c(24),]
all_params_even_minus <- all_params_even[-c(6),]
all_params_even_minus <- all_params_even_minus[-c(24),]

############Generate diffrent combinations of analysis and save into R objects##############
#######ROTATE_AND_SCALE###############
odds <- lapply(1:dim(all_params_odd_minus)[1], function(x){
  image_file <- all_params_odd_minus[[x,1]]
  fovea_height <- all_params_odd_minus[[x,5]]
  fovea_coordx <- all_params_odd_minus[[x,3]]
  fovea_coordy <- all_params_odd_minus[[x,4]]
  fovea_coords <- c(fovea_coordy, fovea_coordx)
  outputz <- rotate_and_plot(image_file, 0.8, 5, fovea_height, fovea_coords, 20)
  return(outputz)
})
#saveRDS(odds, paste("/Users/currant/OCTeyes/odds_fix_final",  strftime(Sys.time(), "%Y%m%d"), sep=""))
odds <- readRDS("/Users/currant/OCTeyes/odds_fix_final20170523")

evens <- lapply(1:dim(all_params_even_minus)[1], function(x){
  image_file <- all_params_even_minus[[x,1]]
  fovea_height <- all_params_even_minus[[x,5]]
  fovea_coordx <- all_params_even_minus[[x,3]]
  fovea_coordy <- all_params_even_minus[[x,4]]
  fovea_coords <- c(fovea_coordy, fovea_coordx)
  outputz <- rotate_and_plot(image_file, 0.8, 5, fovea_height, fovea_coords, 20)
  return(outputz)
})
#saveRDS(evens, paste("/Users/currant/OCTeyes/evens_fix_final",  strftime(Sys.time(), "%Y%m%d"), sep=""))
evens <- readRDS("/Users/currant/OCTeyes/evens_fix_final20170523")

##########JUST SCALE#############
odds_scale <- lapply(1:dim(all_params_odd_minus)[1], function(x){
  image_file <- all_params_odd_minus[[x,1]]
  fovea_height <- all_params_odd_minus[[x,5]]
  fovea_coordx <- all_params_odd_minus[[x,3]]
  fovea_coordy <- all_params_odd_minus[[x,4]]
  fovea_coords <- c(fovea_coordy, fovea_coordx)
  outputz <- just_scale_and_plot(image_file, 0.8, 5, fovea_height, fovea_coords, 20)
  return(outputz)
})
#saveRDS(odds_scale, paste("/Users/currant/OCTeyes/odds_scale_fix_final",  strftime(Sys.time(), "%Y%m%d"), sep=""))
odds_scale <- readRDS("/Users/currant/OCTeyes/odds_scale_fix_final20170523")

evens_scale <- lapply(1:dim(all_params_even_minus)[1], function(x){
  image_file <- all_params_even_minus[[x,1]]
  fovea_height <- all_params_even_minus[[x,5]]
  fovea_coordx <- all_params_even_minus[[x,3]]
  fovea_coordy <- all_params_even_minus[[x,4]]
  fovea_coords <- c(fovea_coordy, fovea_coordx)
  outputz <- just_scale_and_plot(image_file, 0.8, 5, fovea_height, fovea_coords, 20)
  return(outputz)
})
#saveRDS(evens_scale, paste("/Users/currant/OCTeyes/evens_scale_fix_final",  strftime(Sys.time(), "%Y%m%d"), sep=""))
evens_scale <- readRDS("/Users/currant/OCTeyes/evens_scale_fix_final20170523")

#############JUST_ROTATE##############
odds_rotate <- lapply(1:dim(all_params_odd_minus)[1], function(x){
  image_file <- all_params_odd_minus[[x,1]]
  fovea_height <- all_params_odd_minus[[x,5]]
  fovea_coordx <- all_params_odd_minus[[x,3]]
  fovea_coordy <- all_params_odd_minus[[x,4]]
  fovea_coords <- c(fovea_coordy, fovea_coordx)
  outputz <- rotate_and_plot_no_scale(image_file, 0.8, 5, fovea_height, fovea_coords, 20)
  return(outputz)
})
#saveRDS(odds_rotate, paste("/Users/currant/OCTeyes/odds_rotate_fix_final",  strftime(Sys.time(), "%Y%m%d"), sep=""))
odds_rotate <- readRDS("/Users/currant/OCTeyes/odds_rotate_fix_final20170523")

evens_rotate <- lapply(1:dim(all_params_even_minus)[1], function(x){
  image_file <- all_params_even_minus[[x,1]]
  fovea_height <- all_params_even_minus[[x,5]]
  fovea_coordx <- all_params_even_minus[[x,3]]
  fovea_coordy <- all_params_even_minus[[x,4]]
  fovea_coords <- c(fovea_coordy, fovea_coordx)
  outputz <- rotate_and_plot_no_scale(image_file, 0.8, 5, fovea_height, fovea_coords, 20)
  return(outputz)
})
#saveRDS(evens_rotate, paste("/Users/currant/OCTeyes/evens_rotate_fix_final",  strftime(Sys.time(), "%Y%m%d"), sep=""))
evens_rotate <- readRDS("/Users/currant/OCTeyes/evens_rotate_fix_final")

#############JUST_TESSELATE##############
odds_tesselate <- lapply(1:dim(all_params_odd_minus)[1], function(x){
  image_file <- all_params_odd_minus[[x,1]]
  fovea_height <- all_params_odd_minus[[x,5]]
  fovea_coordx <- all_params_odd_minus[[x,3]]
  fovea_coordy <- all_params_odd_minus[[x,4]]
  fovea_coords <- c(fovea_coordy, fovea_coordx)
  outputz <- just_tessalte_and_plot(image_file, 0.8, 5, fovea_height, fovea_coords, 20)
  return(outputz)
})
#saveRDS(odds_tesselate, paste("/Users/currant/OCTeyes/odds_tesselate_fix_final",  strftime(Sys.time(), "%Y%m%d"), sep=""))
odds_tesselate <- readRDS("/Users/currant/O")

evens_tesselate <- lapply(1:dim(all_params_even_minus)[1], function(x){
  image_file <- all_params_even_minus[[x,1]]
  fovea_height <- all_params_even_minus[[x,5]]
  fovea_coordx <- all_params_even_minus[[x,3]]
  fovea_coordy <- all_params_even_minus[[x,4]]
  fovea_coords <- c(fovea_coordy, fovea_coordx)
  outputz <- just_tessalte_and_plot(image_file, 0.8, 5, fovea_height, fovea_coords, 20)
  return(outputz)
})
#saveRDS(evens_tesselate, paste("/Users/currant/OCTeyes/evens_tesselate_fix_final",  strftime(Sys.time(), "%Y%m%d"), sep=""))


odds_minus <- readRDS("/Users/currant/OCTeyes/R_objects/R_comp_objects_final/odds_fix_final20170518")
evens_minus <- readRDS("/Users/currant/OCTeyes/R_objects/R_comp_objects_final/evens_fix_final20170518")
odds_scale <- readRDS("/Users/currant/OCTeyes/odds_scale_fix_final20170523")
evens_scale <- readRDS("/Users/currant/OCTeyes/evens_scale_fix_final20170523")
odds_rotate <- readRDS("/Users/currant/OCTeyes/R_objects/R_comp_objects_final/odds_rotate_fix_final20170518")
evens_rotate <- readRDS("/Users/currant/OCTeyes/R_objects/R_comp_objects_final/evens_rotate_fix_final20170518")
odds_tesselate <- readRDS("/Users/currant/OCTeyes/R_objects/R_comp_objects_final/odds_tesselate_fix_final20170518")
evens_tesselate <- readRDS("/Users/currant/OCTeyes/R_objects/R_comp_objects_final/evens_tesselate_fix_final20170518")

#######ROTATE_AND_SCALE###############
odds_mean_rotate <- lapply(1:dim(all_params_odd_minus)[1], function(x){
  image_file <- all_params_odd_minus[[x,1]]
  fovea_height <- all_params_odd_minus[[x,5]]
  fovea_coordx <- all_params_odd_minus[[x,3]]
  fovea_coordy <- all_params_odd_minus[[x,4]]
  fovea_coords <- c(fovea_coordy, fovea_coordx)
  outputz <- rotate_and_plot_mean_rotate(image_file, 0.8, 5, fovea_height, fovea_coords, 20)
  return(outputz)
})
#saveRDS(odds_mean_rotate, paste("/Users/currant/OCTeyes/odds_mean_rotate_fix_final",  strftime(Sys.time(), "%Y%m%d"), sep=""))
odds_mean_rotate <- readRDS("/Users/currant/OCTeyes/odds_mean_rotate_fix_final20170626")

evens_mean_rotate <- lapply(1:dim(all_params_even_minus)[1], function(x){
  image_file <- all_params_even_minus[[x,1]]
  fovea_height <- all_params_even_minus[[x,5]]
  fovea_coordx <- all_params_even_minus[[x,3]]
  fovea_coordy <- all_params_even_minus[[x,4]]
  fovea_coords <- c(fovea_coordy, fovea_coordx)
  outputz <- rotate_and_plot_mean_rotate(image_file, 0.8, 5, fovea_height, fovea_coords, 20)
  return(outputz)
})
#saveRDS(evens_mean_rotate, paste("/Users/currant/OCTeyes/evens_mean_rotate_fix_final",  strftime(Sys.time(), "%Y%m%d"), sep=""))
evens_mean_rotate <- readRDS("/Users/currant/OCTeyes/evens_mean_rotate_fix_final20170626")

odds_rotate_mean_rotate <- lapply(1:dim(all_params_odd_minus)[1], function(x){
  image_file <- all_params_odd_minus[[x,1]]
  fovea_height <- all_params_odd_minus[[x,5]]
  fovea_coordx <- all_params_odd_minus[[x,3]]
  fovea_coordy <- all_params_odd_minus[[x,4]]
  fovea_coords <- c(fovea_coordy, fovea_coordx)
  outputz <- rotate_and_plot_no_scale_mean_rotate(image_file, 0.8, 5, fovea_height, fovea_coords, 20)
  return(outputz)
})
#saveRDS(odds_rotate_mean_rotate, paste("/Users/currant/OCTeyes/odds_rotate_mean_rotate_fix_final",  strftime(Sys.time(), "%Y%m%d"), sep=""))
odds_rotate_mean_rotate <- readRDS("/Users/currant/OCTeyes/odds_rotate_mean_rotate_fix_final20170626")

evens_rotate_mean_rotate <- lapply(1:dim(all_params_even_minus)[1], function(x){
  image_file <- all_params_even_minus[[x,1]]
  fovea_height <- all_params_even_minus[[x,5]]
  fovea_coordx <- all_params_even_minus[[x,3]]
  fovea_coordy <- all_params_even_minus[[x,4]]
  fovea_coords <- c(fovea_coordy, fovea_coordx)
  outputz <- rotate_and_plot_no_scale_mean_rotate(image_file, 0.8, 5, fovea_height, fovea_coords, 20)
  return(outputz)
})
#saveRDS(evens_rotate_mean_rotate, paste("/Users/currant/OCTeyes/evens_rotate_mean_rotate_fix_final",  strftime(Sys.time(), "%Y%m%d"), sep=""))
evens_rotate_mean_rotate <- readRDS("/Users/currant/OCTeyes/evens_rotate_mean_rotate_fix_final20170626")



########Median normalise all of the values##############
odds_minus_normalised <- lapply(1:length(odds_minus), median_normalise_intensity, datax=odds_minus)
evens_minus_normalised <- lapply(1:length(evens_minus), median_normalise_intensity, datax=evens_minus)
odds_scale_normalised <- lapply(1:length(odds_scale), median_normalise_intensity, datax=odds_scale)
evens_scale_normalised <- lapply(1:length(evens_scale), median_normalise_intensity, datax=evens_scale)
odds_rotate_normalised <- lapply(1:length(odds_rotate), median_normalise_intensity, datax=odds_rotate)
evens_rotate_normalised <- lapply(1:length(evens_rotate), median_normalise_intensity, datax=evens_rotate)
odds_tesselate_normalised <- lapply(1:length(odds_tesselate), median_normalise_intensity, datax=odds_tesselate)
evens_tesselate_normalised <- lapply(1:length(evens_tesselate), median_normalise_intensity, datax=evens_tesselate)

#########add standard triangles to all processed odd images########
odds_minus_norm_stdtri <- lapply(1:length(odds_minus_normalised), add_standard_triangle_all, datax_norm = odds_minus_normalised, datax=odds_minus)
evens_minus_norm_stdtri <- lapply(1:length(evens_minus_normalised), add_standard_triangle_all, datax_norm = evens_minus_normalised, datax=evens_minus)
odds_scale_norm_stdtri <- lapply(1:length(odds_scale_normalised), add_standard_triangle_all, datax_norm = odds_scale_normalised, datax=odds_scale)
evens_scale_norm_stdtri <- lapply(1:length(evens_scale_normalised), add_standard_triangle_all, datax_norm = evens_scale_normalised, datax=evens_scale)
odds_rotate_norm_stdtri <- lapply(1:length(odds_rotate_normalised), add_standard_triangle_all, datax_norm = odds_rotate_normalised, datax = odds_rotate)
evens_rotate_norm_stdtri <- lapply(1:length(evens_rotate_normalised), add_standard_triangle_all, datax_norm = evens_rotate_normalised, datax=evens_rotate)
odds_tesselate_norm_stdtri <- lapply(1:length(odds_tesselate_normalised), add_standard_triangle_all, datax_norm = odds_tesselate_normalised, datax=odds_tesselate)
evens_tesselate_norm_stdtri <- lapply(1:length(evens_tesselate_normalised), add_standard_triangle_all, datax_norm = evens_tesselate_normalised, datax=evens_tesselate)


get_max_min_tri_num <- function(stdtri_table_index, stdrtri_table_list){
  tris <- stdrtri_table_list[[stdtri_table_index]][,4]
  un_tris <- unique(tris)
  min <- min(un_tris)
  max <- max(un_tris)
  return(data.frame(min, max))
}

mandm_odds_minus <- lapply(1:length(odds_minus_norm_stdtri), get_max_min_tri_num, odds_minus_norm_stdtri)
mandm_odds_minus <- do.call(rbind, mandm_odds_minus)
final_min <- min(mandm_odds_minus[,1])
final_max <- max(mandm_odds_minus[,2])


#########get reduced version of table#############
redtab_oddsminus_norm <- lapply(final_min:final_max, get_whole_reduced_table, stdtritab = odds_minus_norm_stdtri)
redtab_evensminus_norm <- lapply(final_min:final_max, get_whole_reduced_table, stdtritab = evens_minus_norm_stdtri)
redtab_oddsscale_norm <- lapply(final_min:final_max, get_whole_reduced_table, stdtritab = odds_scale_norm_stdtri)
redtab_evensscale_norm <- lapply(final_min:final_max, get_whole_reduced_table, stdtritab = evens_scale_norm_stdtri)
redtab_oddsrotate_norm <- lapply(final_min:final_max, get_whole_reduced_table, stdtritab = odds_rotate_norm_stdtri)
redtab_evensrotate_norm <- lapply(final_min:final_max, get_whole_reduced_table, stdtritab = evens_rotate_norm_stdtri)
redtab_oddstesselate_norm <- lapply(final_min:final_max, get_whole_reduced_table, stdtritab = odds_tesselate_norm_stdtri)
redtab_evenstesselate_norm <- lapply(final_min:final_max, get_whole_reduced_table, stdtritab = evens_tesselate_norm_stdtri)


make_ewan_table <- function(list_index, listy_table){
  list_item <- listy_table[[list_index]]
  column_test <- data.frame(list_item[,2])
  colnames(column_test) <- list_item[1,1]
  return(column_test)
}

odd_names <- lapply(1:dim(all_params_odd_minus)[1], function(x, all_params_odd_minus){
  full_name <- all_params_odd_minus[[x,1]]
  split_name <- strsplit(full_name, "/")
  new_name <- paste(unlist(split_name)[[6]], ".", unlist(split_name)[[8]])
  return(new_name)
}, all_params_odd_minus)

even_names <- lapply(1:dim(all_params_even_minus)[1], function(x, all_params_even_minus){
  full_name <- all_params_even_minus[[x,1]]
  split_name <- strsplit(full_name, "/")
  new_name <- paste(unlist(split_name)[[6]], ".", unlist(split_name)[[8]])
  return(new_name)
}, all_params_even_minus)


redtab_oddsminus_ewan <- lapply(1:length(redtab_oddsminus_norm), make_ewan_table, listy_table = redtab_oddsminus_norm)
redtab_oddsminus_ewan <- do.call(cbind,redtab_oddsminus_ewan)
rownames(redtab_oddsminus_ewan) <- odd_names
write.table(redtab_oddsminus_ewan, "/Users/currant/OCTeyes/tables_for_ewan/odds_rotate_and_scale")
redtab_evensminus_ewan <- lapply(1:length(redtab_evensminus_norm), make_ewan_table, listy_table = redtab_evensminus_norm)
redtab_evensminus_ewan <- do.call(cbind, redtab_evensminus_ewan)
rownames(redtab_evensminus_ewan) <- even_names
write.table(redtab_evensminus_ewan, "/Users/currant/OCTeyes/tables_for_ewan/evens_rotate_and_scale")
redtab_oddsscale_ewan <- lapply(1:length(redtab_oddsscale_norm), make_ewan_table, listy_table = redtab_oddsscale_norm)
redtab_oddsscale_ewan <- do.call(cbind, redtab_oddsscale_ewan)
rownames(redtab_oddsscale_ewan) <- odd_names
write.table(redtab_oddsscale_ewan, "/Users/currant/OCTeyes/tables_for_ewan/odds_scale")
redtab_evensscale_ewan <- lapply(1:length(redtab_evensscale_norm), make_ewan_table, listy_table = redtab_evensscale_norm)
redtab_evensscale_ewan <- do.call(cbind, redtab_evensscale_ewan)
rownames(redtab_evensscale_ewan) <- even_names
write.table(redtab_evensscale_ewan, "/Users/currant/OCTeyes/tables_for_ewan/evens_scale")
redtab_oddsrotate_ewan <- lapply(1:length(redtab_oddsrotate_norm), make_ewan_table, listy_table = redtab_oddsrotate_norm)
redtab_oddsrotate_ewan <- do.call(cbind, redtab_oddsrotate_ewan)
rownames(redtab_oddsrotate_ewan) <- odd_names
write.table(redtab_oddsrotate_ewan, "/Users/currant/OCTeyes/tables_for_ewan/odds_rotate")
redtab_evensrotate_ewan <- lapply(1:length(redtab_evensrotate_norm), make_ewan_table, listy_table = redtab_evensrotate_norm)
redtab_evensrotate_ewan <- do.call(cbind, redtab_evensrotate_ewan)
rownames(redtab_evensrotate_ewan) <- even_names
write.table(redtab_evensrotate_ewan, "/Users/currant/OCTeyes/tables_for_ewan/evens_rotate")
redtab_oddstesselate_ewan <- lapply(1:length(redtab_oddstesselate_norm), make_ewan_table, listy_table = redtab_oddstesselate_norm)
redtab_oddstesselate_ewan <- do.call(cbind, redtab_oddstesselate_ewan)
rownames(redtab_oddstesselate_ewan) <- odd_names
write.table(redtab_oddstesselate_ewan, "/Users/currant/OCTeyes/tables_for_ewan/odds_tesselate")
redtab_evenstesselate_ewan <- lapply(1:length(redtab_evenstesselate_norm), make_ewan_table, listy_table = redtab_evenstesselate_norm)
redtab_evenstesselate_ewan <- do.call(cbind, redtab_evenstesselate_ewan)
rownames(redtab_evenstesselate_ewan) <- even_names
write.table(redtab_evenstesselate_ewan, "/Users/currant/OCTeyes/tables_for_ewan/evens_tesselate")

#######Create tables for comparison#########
tess_v_tess <- lapply((final_min+1):(final_max+1), comp_per_tri, comp1=redtab_evenstesselate_norm, comp2=redtab_oddstesselate_norm)
all_v_all <- lapply((final_min+1):(final_max+1), comp_per_tri, comp1=redtab_oddsminus_norm, comp2=redtab_evensminus_norm)
scale_v_scale <- lapply((final_min+1):(final_max+1), comp_per_tri, comp1=redtab_oddsscale_norm, comp2=redtab_evensscale_norm)
rotate_v_rotate <- lapply((final_min+1):(final_max+1), comp_per_tri, comp1=redtab_oddsrotate_norm, comp2=redtab_evensrotate_norm)

######Calculate comparison correlations##########
tess_v_tess_norm_corr <- lapply((final_min+1):(final_max+1), get_correls, comp_table = tess_v_tess)
tess_v_tess_corr <- do.call(rbind, tess_v_tess_norm_corr)

all_v_all_norm_corr <- lapply((final_min+1):(final_max+1), get_correls, comp_table = all_v_all)
all_v_all_corr <- do.call(rbind, all_v_all_norm_corr)

rotate_v_rotate_norm_corr <- lapply((final_min+1):(final_max+1), get_correls, comp_table = rotate_v_rotate)
rotate_v_rotate_corr <- do.call(rbind, rotate_v_rotate_norm_corr)

scale_v_scale_norm_corr <- lapply((final_min+1):(final_max+1), get_correls, comp_table = scale_v_scale)
scale_v_scale_corr <- do.call(rbind, scale_v_scale_norm_corr)

#############Defines the refernce plot ready for plotting back to triangles###############
ref_points <- evens_p[[1]]$all_refs

############Gets table of info to be able to plot back to triangles################
tess_v_tess_cor_triplot <- lapply(1:dim(tess_v_tess_corr)[1], get_triplot_table, ref_points = ref_points, corr_tab = tess_v_tess_corr)
tess_v_tess_cor_triplot <- do.call(rbind, tess_v_tess_cor_triplot)
tess_v_tess_cor_triplot <- as.data.frame(tess_v_tess_cor_triplot)

scale_v_scale_cor_triplot <- lapply(1:dim(scale_v_scale_corr)[1], get_triplot_table, ref_points = ref_points, corr_tab = scale_v_scale_corr)
scale_v_scale_cor_triplot <- do.call(rbind, scale_v_scale_cor_triplot)
scale_v_scale_cor_triplot <- as.data.frame(scale_v_scale_cor_triplot)

rotate_v_rotate_cor_triplot <- lapply(1:dim(rotate_v_rotate_corr)[1], get_triplot_table, ref_points = ref_points, corr_tab = rotate_v_rotate_corr)
rotate_v_rotate_cor_triplot <- do.call(rbind, rotate_v_rotate_cor_triplot)
rotate_v_rotate_cor_triplot <- as.data.frame(rotate_v_rotate_cor_triplot)

all_v_all_cor_triplot <- lapply(1:dim(all_v_all_corr)[1], get_triplot_table, ref_points = ref_points, corr_tab = all_v_all_corr)
all_v_all_cor_triplot <- do.call(rbind, all_v_all_cor_triplot)
all_v_all_cor_triplot <- as.data.frame(all_v_all_cor_triplot)



########################################################
##################MEAN ROTATE ANALYSIS##################
########################################################

########Median normalise all of the values##############
odds_minus_mean_normalised <- lapply(1:length(odds_mean_rotate), median_normalise_intensity, datax=odds_mean_rotate)
evens_minus_mean_normalised <- lapply(1:length(evens_mean_rotate), median_normalise_intensity, datax=evens_mean_rotate)
odds_rotate_mean_normalised <- lapply(1:length(odds_rotate_mean_rotate), median_normalise_intensity, datax=odds_rotate_mean_rotate)
evens_rotate_mean_normalised <- lapply(1:length(evens_rotate_mean_rotate), median_normalise_intensity, datax=evens_rotate_mean_rotate)

#########add standard triangles to all processed odd images########
odds_minus_mean_norm_stdtri <- lapply(1:length(odds_minus_mean_normalised), add_standard_triangle_all, datax_norm = odds_minus_mean_normalised, datax=odds_mean_rotate)
evens_minus_mean_norm_stdtri <- lapply(1:length(evens_minus_mean_normalised), add_standard_triangle_all, datax_norm = evens_minus_mean_normalised, datax=evens_mean_rotate)
odds_rotate_mean_norm_stdtri <- lapply(1:length(odds_rotate_mean_normalised), add_standard_triangle_all, datax_norm = odds_rotate_mean_normalised, datax = odds_rotate_mean_rotate)
evens_rotate_mean_norm_stdtri <- lapply(1:length(evens_rotate_mean_normalised), add_standard_triangle_all, datax_norm = evens_rotate_mean_normalised, datax=evens_rotate_mean_rotate)

get_max_min_tri_num <- function(stdtri_table_index, stdrtri_table_list){
  tris <- stdrtri_table_list[[stdtri_table_index]][,4]
  un_tris <- unique(tris)
  min <- min(un_tris)
  max <- max(un_tris)
  return(data.frame(min, max))
}

mandm_odds_minus <- lapply(1:length(odds_minus_mean_norm_stdtri), get_max_min_tri_num, odds_minus_mean_norm_stdtri)
mandm_odds_minus <- do.call(rbind, mandm_odds_minus)
final_min <- min(mandm_odds_minus[,1])
final_max <- max(mandm_odds_minus[,2])


#########get reduced version of table#############
redtab_oddsminus_mean_norm <- lapply(final_min:final_max, get_whole_reduced_table, stdtritab = odds_minus_mean_norm_stdtri)
redtab_evensminus_mean_norm <- lapply(final_min:final_max, get_whole_reduced_table, stdtritab = evens_minus_mean_norm_stdtri)
redtab_oddsrotate_mean_norm <- lapply(final_min:final_max, get_whole_reduced_table, stdtritab = odds_rotate_mean_norm_stdtri)
redtab_evensrotate_mean_norm <- lapply(final_min:final_max, get_whole_reduced_table, stdtritab = evens_rotate_mean_norm_stdtri)

#######Create tables for comparison#########
all_v_all_mean <- lapply((final_min+1):(final_max+1), comp_per_tri, comp1=redtab_oddsminus_mean_norm, comp2=redtab_evensminus_mean_norm)
rotate_v_rotate_mean <- lapply((final_min+1):(final_max+1), comp_per_tri, comp1=redtab_oddsrotate_mean_norm, comp2=redtab_evensrotate_mean_norm)

######Calculate comparison correlations##########
all_v_all_mean_norm_corr <- lapply((final_min+1):(final_max+1), get_correls, comp_table = all_v_all_mean)
all_v_all_mean_corr <- do.call(rbind, all_v_all_mean_norm_corr)

rotate_v_rotate_mean_norm_corr <- lapply((final_min+1):(final_max+1), get_correls, comp_table = rotate_v_rotate_mean)
rotate_v_rotate_mean_corr <- do.call(rbind, rotate_v_rotate_mean_norm_corr)

#############Defines the refernce plot ready for plotting back to triangles###############
ref_points <- evens_mean_rotate[[1]]$all_refs

############Gets table of info to be able to plot back to triangles################
rotate_v_rotate_mean_cor_triplot <- lapply(1:dim(rotate_v_rotate_mean_corr)[1], get_triplot_table, ref_points = ref_points, corr_tab = rotate_v_rotate_mean_corr)
rotate_v_rotate_mean_cor_triplot <- do.call(rbind, rotate_v_rotate_mean_cor_triplot)
rotate_v_rotate_mean_cor_triplot <- as.data.frame(rotate_v_rotate_mean_cor_triplot)

all_v_all_mean_cor_triplot <- lapply(1:dim(all_v_all_mean_corr)[1], get_triplot_table, ref_points = ref_points, corr_tab = all_v_all_mean_corr)
all_v_all_mean_cor_triplot <- do.call(rbind, all_v_all_mean_cor_triplot)
all_v_all_mean_cor_triplot <- as.data.frame(all_v_all_mean_cor_triplot)


p <- ggplot(data=rotate_v_rotate_mean_cor_triplot)
p + 
  geom_point(aes(x=V1, y=V2, col=V4), na.rm=TRUE) + 
  theme_bw()

p <- ggplot(data=all_v_all_mean_cor_triplot)
p + 
  geom_point(aes(x=V1, y=V2, col=V4), na.rm=TRUE) + 
  theme_bw()

all_v_all_mean_corr$type <- 'mean' 
all_v_all_corr$type <- 'norm'

rotate_v_rotate_mean_corr$type <-'mean'
rotate_v_rotate_corr$type <- 'norm'


tess_v_tess_cor_triplot$type <- 'tess'
scale_v_scale_cor_triplot



all_v_all_corr_mean_comp <- rbind(all_v_all_corr, all_v_all_mean_corr)

ggplot(all_v_all_corr_mean_comp, aes(correl, fill=type))+ geom_histogram(alpha = 0.5, aes(y=..density..), position = "identity")
ggplot(all_v_all_corr_mean_comp, aes(correl, fill=type))+ geom_density(alpha=0.2)

rotate_v_rotate_corr_mean_comp <- rbind(rotate_v_rotate_corr, rotate_v_rotate_mean_corr)
ggplot(rotate_v_rotate_corr_mean_comp, aes(correl, fill=type))+ geom_histogram(alpha = 0.5, aes(y=..density..), position = "identity")
ggplot(rotate_v_rotate_corr_mean_comp, aes(correl, fill=type))+ geom_density(alpha=0.2)

ggplot(all_v_all_corr_qc_comp, aes(correl, fill=qc_type))+ geom_histogram(alpha=0.5, aes(y=..density..), position = "identity")
ggplot(all_v_all_corr_qc_comp, aes(correl, fill=qc_type))+geom_density(alpha=0.2)


##################MAKING GRID OF FOUR CORR TRIPLOTS
all_v_all_cor_triplot$type <- 'all'
tess_v_tess_cor_triplot$type <- 'tess'
rotate_v_rotate_cor_triplot$type <- 'rotate'
scale_v_scale_cor_triplot$type <- 'scale'


all_v_all_corr$type <- 'all'
tess_v_tess_corr$type <- 'tess'
scale_v_scale_corr$type <- 'scale'
rotate_v_rotate_corr$type <- 'rotate'

df_all <- rbind(tess_v_tess_cor_triplot, scale_v_scale_cor_triplot, rotate_v_rotate_cor_triplot, all_v_all_cor_triplot)

df_transforms <- rbind(tess_v_tess_corr, scale_v_scale_corr, rotate_v_rotate_corr, all_v_all_corr)

ggplot(df_transforms, aes(correl, fill=type))+ geom_histogram(alpha=0.5, aes(y=..density..), position = "identity")
ggplot(df_transforms, aes(correl, fill=type))+geom_density(alpha=0.2) + labs(x = 'Correlation', y='Density')

p <- ggplot(df_all, aes(x=V1, y=V2,, col=V4), na.rm=TRUE) 
p+  
  facet_wrap(.~type) +
  theme_bw() +
  scale_color_gradient(low="snow", high="deeppink4")



p <- ggplot(data=scale_v_scale_cor_triplot)
p + 
  geom_point(aes(x=V1, y=V2, col=V4), na.rm=TRUE) + 
  theme_bw() +
  scale_color_gradient(low="palevioletred3", high="cadetblue")

#################################################
####################PLOTTING#####################
#################################################

###########Plotting scatter plots of correlation vs triangle############

pdf("/Users/currant/OCTeyes/scale_and_rotate_graphs/final_plots_we_hope/per_triangle_mean_norm_correlation_fix.pdf")
plot(tess_v_tess_corr, xlab="Triangle ID", ylab="Correlation", main="Tesselation vs tesselation per triangle")
plot(all_v_all_corr, xlab="Triangle ID", ylab="Correlation", main="All vs all per triangle")
plot(scale_v_scale_corr, xlab="Triangle ID", ylab="Correlation", main="Scale vs scale per triangle")
plot(rotate_v_rotate_corr, xlab="Triangle ID", ylab="Correlation", main="Rotate vs rotate per triangle")
dev.off()

############PLotting distribution of correlations ############
pdf("/Users/currant/OCTeyes/scale_and_rotate_graphs/final_plots_we_hope/per_triangle_mean_norm_correlation_histogram_fix.pdf")
hist(tess_v_tess_corr[,2], xlab="Correlation", main="Tesselation vs Tesselation")
hist(all_v_all_corr[,2], xlab="Correlation", main="All vs All")
hist(scale_v_scale_corr[,2], xlab="Correlation", main="Scale vs Scale")
hist(rotate_v_rotate_corr[,2], xlab="Correlation", main="Rotate vs Rotate")
dev.off()

#########Plotting back to triangle############
#MAKING TRANSFORMATION GRID
tess_v_tess_cor_triplot$type <- 'col1'
scale_v_scale_cor_triplot$type <- 'col2'
rotate_v_rotate_cor_triplot$type <- 'col1'
all_v_all_cor_triplot$type <- 'col2'
tess_v_tess_cor_triplot$type2 <- 'row1'
scale_v_scale_cor_triplot$type2 <- 'row1'
rotate_v_rotate_cor_triplot$type2 <- 'row2'
all_v_all_cor_triplot$type2 <- 'row2'

all_transform_df <- rbind(tess_v_tess_cor_triplot, scale_v_scale_cor_triplot, rotate_v_rotate_cor_triplot, all_v_all_cor_triplot)

p_tess <- ggplot(data = tess_v_tess_cor_triplot)
p_tess <- p_tess +
  geom_point(aes(x=V1, y=V2, col=V4), na.rm=TRUE) +
  theme_bw() +
  scale_color_gradient(low="snow", high="violetred4")
  
p_scale <- ggplot(data = scale_v_scale_cor_triplot)
p_scale <- p_scale +
  geom_point(aes(x=V1, y=V2, col=V4), na.rm=TRUE) +
  theme_bw() +
  scale_color_gradient(low="snow", high="violetred4")

p_rotate <- ggplot(data = rotate_v_rotate_cor_triplot)
p_rotate <- p_rotate +
  geom_point(aes(x=V1, y=V2, col=V4), na.rm=TRUE) +
  theme_bw() +
  scale_color_gradient(low="snow", high="violetred4")

p_all <- ggplot(data = all_v_all_cor_triplot)
p_all <- p_all +
  geom_point(aes(x=V1, y=V2, col=V4), na.rm=TRUE) +
  theme_bw() +
  scale_color_gradient(low="snow", high="violetred4")


grid.arrange(p_tess, p_scale, p_rotate, p_all, ncol=2)
plot_grid(p_tess, p_scale, p_rotate, p_all, labels = 'AUTO', ncol=2)


p_df <- ggplot(data = all_transform_df)
p_df +
  geom_point(aes(x=V1, y=V2, col=V4), na.rm=TRUE) +
  theme_bw() +
  scale_color_gradient(low="snow", high="violetred4") +
  facet_grid(type2~type) +
  theme(strip.background = element_blank(), strip.text.x = element_blank(), strip.text.y = element_blank()) +
  labs(colour = "Correlation \nbetween \nreplicates", x = "", y="") +
  theme(legend.key.size = unit(0.5, "in")) + 
  theme(legend.text = element_text(size = 12)) + 
  theme(legend.title = element_text(size=12)) +
  theme(title = element_text(size=18))

#MAKE COMP GRAPHS
pdf("/Users/currant/OCTeyes/scale_and_rotate_graphs/final_plots_we_hope/back_to_tris.pdf")
p <- ggplot(data=tess_v_tess_cor_triplot)
p + 
  geom_point(aes(x=V1, y=V2, col=V4), na.rm=TRUE) + 
  theme_bw() +
  scale_color_gradient(low="palevioletred3", high="cadetblue")


#INTENSITY PLOT
p <- ggplot(data=as.data.frame(redtab_evenstesselate_norm[[1]]))
p + 
  geom_point(aes(x=merged_all_refs_tmp3.x, y=merged_all_refs_tmp3.y, col=merged_all_refs_tmp3.mean), na.rm=TRUE, ) +
  #geom_point(data=as.data.frame(verticesdf), aes(x=V1, y=V2), color='red')  +
  scale_color_gradient(low = 'darkcyan', high = 'lightcyan') +
  labs(colour = "Normalised \nDepth", x = "", y="") +
  theme_bw() +
  theme(legend.key.size = unit(0.5, "in")) + 
  theme(legend.text = element_text(size = 12)) + 
  theme(legend.title = element_text(size=12)) +
  theme(title = element_text(size=18))


#Crosssection with bands
eye <- function (image_data, zpos, bands, name) {
  img = image_data$img
  image(t(img[, , zpos]), breaks = c(0, 50, 200), col = c("black", 
                                                          "cyan4"), main = name)
  for(x in 1:ncol(bands)) {
    matplot(1:dim(img)[2]/dim(img)[2], bands[, x]/dim(img)[1], 
            col = "hotpink", type = "l", add = T, lwd=1.5)
  }
  #matplot(1:dim(img)[2]/dim(img)[2], bands[, 2]/dim(img)[1], 
  #    col = "red", type = "l", add = T)
}


#####################################################
################LEFT AND RIGHT EYES##################
#####################################################
left_index <- which(all_params_odd_qc[,6]=="L")
right_index <- which(all_params_even_qc[,6]=="R")


mean_lorr <- lapply(1:length(redtab_oddstesselate_norm), making_mean_lorr, left_index = left_index, right_index = right_index, odd_redtab = redtab_oddstesselate_norm, even_redtab = redtab_evenstesselate_norm)
right_means <- lapply(1:length(mean_lorr), function(x){
  line <- mean_lorr[[x]]$rights_mean
  return(line)
})
right_means_table <- do.call(rbind, right_means)

left_means <- lapply(1:length(mean_lorr), function(x){
  line <- mean_lorr[[x]]$lefts_mean
  return(line)
})
left_means_table <- do.call(rbind, left_means)

right_means_triplot <- lapply(1:dim(right_means_table)[1], get_triplot_table, ref_points = ref_points, corr_tab = right_means_table)
right_means_triplot <- do.call(rbind, right_means_triplot)
right_means_triplot <- as.data.frame(right_means_triplot)

right_means_triplot_plot <- ggplot(data = right_means_triplot)
right_means_triplot_plot <- right_means_triplot_plot +
  geom_point(aes(x=V1, y=V2, col=V4), na.rm=TRUE) +
  theme_bw() +
  scale_color_gradient(low="lightgoldenrod1", high="violetred4")
right_means_triplot_plot

left_means_triplot <- lapply(1:dim(left_means_table)[1], get_triplot_table, ref_points = ref_points, corr_tab = left_means_table)
left_means_triplot <- do.call(rbind, left_means_triplot)
left_means_triplot <- as.data.frame(left_means_triplot)

left_means_triplot_plot <- ggplot(data = left_means_triplot)
left_means_triplot_plot <- left_means_triplot_plot +
  geom_point(aes(x=V1, y=V2, col=V4), na.rm=TRUE) +
  theme_bw() +
  scale_colour_gradient(low="lightgoldenrod1", high="violetred4")
left_means_triplot_plot

right_means_triplot$lorr <- "right eye"
left_means_triplot$lorr <- "left eye"
joint_lorr_table <- rbind(right_means_triplot, left_means_triplot)
joint_lorr_table$lorr <- factor(joint_lorr_table$lorr, levels = c("right eye", "left eye"))

joint_lorr_triplot_plot <- ggplot(data = joint_lorr_table)
joint_lorr_triplot_plot <- joint_lorr_triplot_plot +
  geom_point(aes(x=V2, y=V1, col=V4), na.rm=TRUE) +
  theme_bw() +
  scale_colour_gradient(low="lightgoldenrod1", high="violetred4") +
  facet_grid(.~lorr) +
  labs(colour = "Mean\nNormalised\nIntensity", x="", y="") +
  theme(strip.text.x = element_text(size= 14)) +
  theme(legend.text = element_text(size=14)) +
  theme(strip.background = element_blank()) +
  theme(axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank()) +
  coord_fixed(ratio=1)
joint_lorr_triplot_plot

joint_lorr_table_test <- factor(joint_lorr_table$lorr, levels = c("right", "left"))


####################################


ggplot(biggest, aes(correl, fill=type))+geom_density(alpha=0.2) + 
  facet_grid(qc~.) +
  labs(x ='Correlation', y='Density', fill="Transformation") +
  theme_bw() +
  theme(legend.key.size = unit(0.3, "in")) + 
  theme(legend.text = element_text(size = 12)) + 
  theme(legend.title = element_text(size=12)) +
  theme(strip.text.x = element_text(size = 12)) +
  coord_fixed(ratio = 0.05) +
  theme(strip.background = element_blank())



image_file = "/Users/currant/OCTeyes/OCT_30_replicates_20161024/1145517/OCT/292422.dcm"
image_data <- dicomInfo(image_file)
bands <- fit_major_bands(image_data, 60)
bands_i <- interpolate_outliers(bands)
test <- eye(image_data, 60, bands, "")



evens_tesselate_normalised_red <- evens_tesselate_normalised_red/(length(evens_tesselate_normalised))
evens_tesselate_normalised_red <- Reduce('+', evens_tesselate_normalised)


p <- ggplot(data=scale_v_scale_cor_triplot)
p + 
  geom_point(aes(x=V1, y=V2, col=V4), na.rm=TRUE) + 
  theme_bw()


p <- ggplot(data=rotate_v_rotate_cor_triplot)
p + 
  geom_point(aes(x=V1, y=V2, col=V4), na.rm=TRUE) + 
  theme_bw()

p <- ggplot(data=all_v_all_cor_triplot)
p + 
  geom_point(aes(x=V1, y=V2, col=V4), na.rm=TRUE) + 
  theme_bw()
dev.off()



###########################################
gen_ellipse <- function(x){
e <- registar(image_file, debug=FALSE)
smoothed_ellipse1 <- smooth_ellipse_max_slope(e, 31)
smoothed_ellipse2 <- smooth_ellipse_plateu(e,31)
return(list(max_slope_e = smoothed_ellipse1, plateu_e = smoothed_ellipse2))
}


all_even_ellipse <- lapply(all_params_even_qc[,1], gen_ellipse)
all_odds_ellipse <- lapply(all_params_odd_qc[,1], gen_ellipse)
odd_max_el_x <- lapply(1:length(all_odds_ellipse), function(x){
  two_cols <- all_odds_ellipse[[x]]$max_slope_e[,1]
  return(two_cols)
})
odd_max_el_x <- do.call(cbind, odd_max_el_x)
dim(odd_max_el_x)
mean_odd_max_el_x <- rowMeans(odd_max_el_x)
mean_odd_max_el_x

odd_max_el_y <- lapply(1:length(all_odds_ellipse), function(x){
  two_cols <- all_odds_ellipse[[x]]$max_slope_e[,2]
  return(two_cols)
})
odd_max_el_y <- do.call(cbind, odd_max_el_y)
mean_odd_max_el_y <- rowMeans(odd_max_el_y)
mean_odd_max_el <- cbind(mean_odd_max_el_x, mean_odd_max_el_y)


even_max_el_x <- lapply(1:length(all_even_ellipse), function(x){
  two_cols <- all_even_ellipse[[x]]$max_slope_e[,1]
  return(two_cols)
})
even_max_el_x <- do.call(cbind, even_max_el_x)
mean_even_max_el_x <- rowMeans(even_max_el_x)

even_max_el_y <- lapply(1:length(all_even_ellipse), function(x){
  two_cols <- all_even_ellipse[[x]]$max_slope_e[,2]
  return(two_cols)
})
even_max_el_y <- do.call(cbind, even_max_el_y)
mean_even_max_el_y <- rowMeans(even_max_el_y)

mean_max_x <- cbind(odd_max_el_x, even_max_el_x)
mean_max_x <- rowMeans(mean_max_x)

mean_max_y <- cbind(odd_max_el_y, even_max_el_y)
mean_max_y <- rowMeans(mean_max_y)

mean_max_ellipse <- cbind(mean_max_x, mean_max_y)


odd_plat_el_x <- lapply(1:length(all_odds_ellipse), function(x){
  two_cols <- all_odds_ellipse[[x]]$plateu_e[,1]
  return(two_cols)
})
odd_plat_el_x <- do.call(cbind, odd_plat_el_x)
mean_odd_plat_el_x <- rowMeans(odd_plat_el_x)

odd_plat_el_y <- lapply(1:length(all_odds_ellipse), function(x){
  two_cols <- all_odds_ellipse[[x]]$plateu_e[,2]
  return(two_cols)
})
odd_plat_el_y <- do.call(cbind, odd_plat_el_y)
mean_odd_plat_el_y <- rowMeans(odd_plat_el_y)

even_plat_el_x <- lapply(1:length(all_even_ellipse), function(x){
  two_cols <- all_even_ellipse[[x]]$plateu_e[,1]
  return(two_cols)
})
even_plat_el_x <- do.call(cbind, even_plat_el_x)
mean_even_plat_el_x <- rowMeans(even_plat_el_x)

even_plat_el_y <- lapply(1:length(all_even_ellipse), function(x){
  two_cols <- all_even_ellipse[[x]]$plateu_e[,2]
  return(two_cols)
})
even_plat_el_y <- do.call(cbind, even_plat_el_y)
mean_even_plat_el_y <- rowMeans(even_plat_el_y)

mean_plat_x <- cbind(odd_plat_el_x, even_plat_el_x)
mean_plat_x <- rowMeans(mean_plat_x)

mean_plat_y <- cbind(odd_plat_el_y, even_plat_el_y)
mean_plat_y <- rowMeans(mean_plat_y)

mean_plat_ellipse <- cbind(mean_plat_x, mean_plat_y)
plot(mean_max_ellipse)
mean_max_ellipse_df <- as.data.frame(mean_max_ellipse)
mean_plat_ellipse_df <- as.data.frame(mean_plat_ellipse)
mean_max_ellipse_df$ellipsetype <- 'max'
mean_plat_ellipse_df$ellipsetype <- 'plat'
colnames(mean_max_ellipse_df) <- c("x", "y", "ellipse_type")
colnames(mean_plat_ellipse_df) <- c("x", "y", "ellipse_type")
mean_ellipse <- rbind(mean_max_ellipse_df, mean_plat_ellipse_df)
mean_ellipse
plot(mean_ellipse[,1:2], colour = "ellipse_type")
mean_ellipse$x <- mean_ellipse$x/512
mean_ellipse$y <- mean_ellipse$y/512
image(t(amodelsss_even_mean[nrow(amodelsss_even_mean):1,]), col=heat_hcl(10), xaxt='n', yaxt='n')
matplot(mean_ellipse$x, mean_ellipse$y, add=TRUE, pch="•")



amodelsss_even <- lapply(all_params_even_qc[,1], function(x){
  model <- run_generate_model_and_find_fovea_dip_poisition(x)
  amodel <- center_align(model, 0.8)
  return(amodel)
})

amodelsss_odd <- lapply(all_params_odd_qc[,1], function(x){
  model <- run_generate_model_and_find_fovea_dip_poisition(x)
  amodel <- center_align(model, 0.8)
  return(amodel)
})

amodelsss_even_mean <- Reduce('+', amodelsss_even)/length(amodelsss_even)
amodelss
amodelsss_odd_mean <- Reduce('+', amodelsss_odd)



image_data = dicomInfo(image_file)
fac=0.15
bulls_eye_x = sapply(1:dim(image_data$img)[3], function(i) apply(image_data$img[,,i],1, mean))
aband = apply(bulls_eye_x, 2, mean)
l1 = lowess(aband, f=fac)
center = get_fovea_slice_robust(l1)

plot(apply(image_data$img[,,center],1, mean), xlab="Index", ylab="Intensity")
image(bulls_eye_x[nrow(bulls_eye_x):1,], col=rev(rainbow(100)), xlab="Mean Intensity", ylab="Slice Index")
plot(aband, xlab="Index", ylab="Slice Depth Estimate")
matplot(l1$x, l1$y, col="magenta", type="l", add=T)
abline(v=center)


