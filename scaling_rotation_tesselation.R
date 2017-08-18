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

Add_standard_triangle <- function(verticesdf, end_info2){
  triangle_groups <- seq(from=1, to=dim(verticesdf)[1], by=3)
  test_mctest <- lapply(triangle_groups, function(x){
    one_triangle_bnd <- verticesdf[x:(x+2),]
    matrix <- as.matrix(end_info2[,1:2])
    inouts <- in.out(one_triangle_bnd, matrix)
    TRUTHS <- which(inouts)
    triangle_friends <- end_info2[TRUTHS,]
    triangle_friends$tri_number <- x
    return(triangle_friends)
  })
  test_mctest2 <- do.call(rbind, test_mctest)
  return(test_mctest2)
}

#Function that adds_standards triangles to the outputs of rotate, plot, scale or tesselation
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

#Iterate through 
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
  medi <- median(tmp3, na.rm=TRUE)
  tmp3[,3] <- tmp3[,3]-medi
  return(tmp3)
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

#saveRDS(all_params_even, paste("/Users/currant/OCTeyes/all_params_even",  strftime(Sys.time(), "%Y%m%d"), sep=""))
all_params_even <- readRDS("/Users/currant/OCTeyes/R_objects/all_params_even20170425")

#saveRDS(all_params_odd, paste("/Users/currant/OCTeyes/all_params_odd",  strftime(Sys.time(), "%Y%m%d"), sep=""))
all_params_odd <- readRDS("/Users/currant/OCTeyes/R_objects/all_params_odd20170425")

#saveRDS(biggy_table, paste("/Users/currant/OCTeyes/biggy_table",  strftime(Sys.time(), "%Y%m%d"), sep=""))
biggy_table <- readRDS("/Users/currant/OCTeyes/R_objects/biggy_table20170425")

#saveRDS(turning_outcome, paste("/Users/currant/OCTeyes/turning_outcome",  strftime(Sys.time(), "%Y%m%d"), sep=""))
turning_outcome <- readRDS("/Users/currant/OCTeyes/R_objects/turning_outcome20170425")


#########NEW UPDATED READ IN ODDS ETC###############
odds_minus <- readRDS("/Users/currant/OCTeyes/R_objects/R_comp_objects_final/odds_fix_final20170518")
evens_mins <- readRDS("/Users/currant/OCTeyes/R_objects/R_comp_objects_final/evens_fix_final20170518")
odds_scale <- readRDS("/Users/currant/OCTeyes/R_objects/R_comp_objects_final/odds_scale_fix_final20170518")
evens_scale <- readRDS("/Users/currant/OCTeyes/R_objects/R_comp_objects_final/evens_scale_fix_final20170518")
odds_rotate <- readRDS("/Users/currant/OCTeyes/R_objects/R_comp_objects_final/odds_rotate_fix_final20170518")
evens_rotate <- readRDS("/Users/currant/OCTeyes/R_objects/R_comp_objects_final/evens_rotate_fix_final20170518")
odds_tesselate <- readRDS("/Users/currant/OCTeyes/R_objects/R_comp_objects_final/odds_tesselate_fix_final20170518")
evens_tesselate <- readRDS("/Users/currant/OCTeyes/R_objects/R_comp_objects_final/evens_tesselate_fix_final20170518")



#saveRDS(odds_minus, paste("/Users/currant/OCTeyes/odds_minus",  strftime(Sys.time(), "%Y%m%d"), sep=""))
odds_minus <- readRDS("/Users/currant/OCTeyes/R_comp_objects/odds_minus20170427")

#saveRDS(evens, paste("/Users/currant/OCTeyes/evens",  strftime(Sys.time(), "%Y%m%d"), sep=""))
evens_minus <- readRDS("/Users/currant/OCTeyes/R_comp_objects/evens20170427")

#saveRDS(odds_scale, paste("/Users/currant/OCTeyes/odds_scale",  strftime(Sys.time(), "%Y%m%d"), sep=""))
odds_scale <- readRDS("/Users/currant/OCTeyes/R_comp_objects/odds_scale20170427")

#saveRDS(evens_scale, paste("/Users/currant/OCTeyes/evens_scale",  strftime(Sys.time(), "%Y%m%d"), sep=""))
evens_scale <- readRDS("/Users/currant/OCTeyes/R_comp_objects/evens_scale20170428")

#saveRDS(odds_rotate, paste("/Users/currant/OCTeyes/odds_rotate",  strftime(Sys.time(), "%Y%m%d"), sep=""))
odds_rotate <- readRDS("/Users/currant/OCTeyes/R_comp_objects/odds_rotate20170428")

#saveRDS(evens_rotate, paste("/Users/currant/OCTeyes/evens_rotate",  strftime(Sys.time(), "%Y%m%d"), sep=""))
evens_rotate <- readRDS("/Users/currant/OCTeyes/R_comp_objects/evens_rotate20170502")

#saveRDS(odds_tesselate, paste("/Users/currant/OCTeyes/odds_tesselate",  strftime(Sys.time(), "%Y%m%d"), sep=""))
odds_tesselate <- readRDS("/Users/currant/OCTeyes/R_comp_objects/odds_tesselate20170502")

#saveRDS(evens_tesselate, paste("/Users/currant/OCTeyes/evens_tesselate",  strftime(Sys.time(), "%Y%m%d"), sep=""))
evens_tesselate <- readRDS("/Users/currant/OCTeyes/R_comp_objects/evens_tesselate20170502")


###################################
#############ANALYSIS##############
###################################

image_file = "/Users/currant/OCTeyes/OCT_30_replicates_20161024/1331287/OCT/306418.dcm"
new_im <- rotate_and_plot_byangle(0, image_file, centre_f, 5, fovea_height, fovea_coords, n)

p <- ggplot(data=tstststst)
p + 
  geom_point(aes(x=merged_all_refs_tmp3.x, y=merged_all_refs_tmp3.y, col=merged_all_refs_tmp3.mean), na.rm=TRUE) + 
  geom_point(data=as.data.frame(odds[[7]]$verticesdf, aes(x=V1, y=V2), color='red')  +
  theme_bw()

lapply(turning_outcome, function(x){
  outcome <- add_s
})

############Get paramaters for following analysis#############

all_eye_paths <- get_eye_paths(pathname)
all_eye_params10 <- get_needed_params(all_eye_paths[1:10])
all_eye_paths_odd <- get_eye_pathsrepeats(pathname, 1)
all_params_odd <- get_needed_params(all_eye_paths_odd)
all_eye_paths_even <- get_eye_pathsrepeats(pathname, 2)
all_params_even <- get_needed_params(all_eye_paths_even)

###########Adjust the paramaters to exclude some non-working scans, adjust so that this is done using quality control at some point############

all_params_odd_minus <- all_params_odd[-c(6),]
all_params_odd_minus <- all_params_odd_minus[-c(24),]
all_params_even_minus <- all_params_even[-c(6),]
all_params_even_minus <- all_params_even_minus[-c(24),]

###########Practice analysis on rotation outcomes #################
#Add standard triangle to all of the rotated images
all_stand_tri_ids <- lapply(turning_outcome, function(x){
  end_info <- x$end_info2
  verts <- x$verticesdf
  plus_tris <- add_standard_triangle(end_info, verts)
  return(plus_tris)
})

# crates a reduced table for rotated_triangles
reduced_tables <- lapply(all_stand_tri_ids, function(x){
  tris <- unique(x$tri_number)
  table1 <- x
  per_tri <- lapply(tris, function(x){
    one_t <- filter(table1, table1$tri_number==x)
    #print(str(one_t))
    #one_t <- complete.cases(one_t)
    val <- mean(one_t$merged_all_refs_tmp3.mean, na.rm=TRUE)
    row <- c(x,val)
    return(row)
  })
  per_tri_bound <- do.call(rbind, per_tri)
  return(per_tri_bound)
})


############Generate diffrent combinations of analysis##############
#######ROTATE_AND_SCALE###############
odds <- lapply(1:dim(all_params_odd_minus)[1], function(x){
  image_file <- all_params_odd_minus[[x,1]]
  fovea_height <- all_params_odd_minus[[x,5]]
  #fovea_height_scale <- fovea_height*(dim(amodel)[1])
  fovea_coordx <- all_params_odd_minus[[x,3]]
  #fovea_coordx_scale <- fovea_coordx*(dim(amodel)[2])
  fovea_coordy <- all_params_odd_minus[[x,4]]
  #fovea_coordy_scale <- fovea_coordy*(dim(amodel)[1])
  fovea_coords <- c(fovea_coordy, fovea_coordx)
  #fovea_coords <- c(fovea_coordy_scale, fovea_coordx_scale)
  outputz <- rotate_and_plot(image_file, 0.8, 5, fovea_height, fovea_coords, 20)
  return(outputz)
})



#saveRDS(odds_minus, paste("/Users/currant/OCTeyes/odds_minus",  strftime(Sys.time(), "%Y%m%d"), sep=""))
saveRDS(odds, paste("/Users/currant/OCTeyes/odds_fix_final",  strftime(Sys.time(), "%Y%m%d"), sep=""))


evens <- lapply(1:dim(all_params_even_minus)[1], function(x){
  image_file <- all_params_even_minus[[x,1]]
  fovea_height <- all_params_even_minus[[x,5]]
  #fovea_height_scale <- fovea_height*(dim(amodel)[1])
  fovea_coordx <- all_params_even_minus[[x,3]]
  #fovea_coordx_scale <- fovea_coordx*(dim(amodel)[2])
  fovea_coordy <- all_params_even_minus[[x,4]]
  #fovea_coordy_scale <- fovea_coordy*(dim(amodel)[1])
  fovea_coords <- c(fovea_coordy, fovea_coordx)
  #fovea_coords <- c(fovea_coordy_scale, fovea_coordx_scale)
  outputz <- rotate_and_plot(image_file, 0.8, 5, fovea_height, fovea_coords, 20)
  return(outputz)
})

#saveRDS(evens, paste("/Users/currant/OCTeyes/evens",  strftime(Sys.time(), "%Y%m%d"), sep=""))
saveRDS(evens, paste("/Users/currant/OCTeyes/evens_fix_final",  strftime(Sys.time(), "%Y%m%d"), sep=""))

##########JUST SCALE#############
odds_scale <- lapply(1:dim(all_params_odd_minus)[1], function(x){
  image_file <- all_params_odd_minus[[x,1]]
  fovea_height <- all_params_odd_minus[[x,5]]
  #fovea_height_scale <- fovea_height*(dim(amodel)[1])
  fovea_coordx <- all_params_odd_minus[[x,3]]
  #fovea_coordx_scale <- fovea_coordx*(dim(amodel)[2])
  fovea_coordy <- all_params_odd_minus[[x,4]]
  #fovea_coordy_scale <- fovea_coordy*(dim(amodel)[1])
  fovea_coords <- c(fovea_coordy, fovea_coordx)
  #fovea_coords <- c(fovea_coordy_scale, fovea_coordx_scale)
  outputz <- just_scale_and_plot(image_file, 0.8, 5, fovea_height, fovea_coords, 20)
  return(outputz)
})

#saveRDS(odds_scale, paste("/Users/currant/OCTeyes/odds_scale",  strftime(Sys.time(), "%Y%m%d"), sep=""))
saveRDS(odds_scale, paste("/Users/currant/OCTeyes/odds_scale_fix_final",  strftime(Sys.time(), "%Y%m%d"), sep=""))

x=10


evens_scale <- lapply(1:dim(all_params_even_minus)[1], function(x){
  image_file <- all_params_even_minus[[x,1]]
  fovea_height <- all_params_even_minus[[x,5]]
  #fovea_height_scale <- fovea_height*(dim(amodel)[1])
  fovea_coordx <- all_params_even_minus[[x,3]]
  #fovea_coordx_scale <- fovea_coordx*(dim(amodel)[2])
  fovea_coordy <- all_params_even_minus[[x,4]]
  #fovea_coordy_scale <- fovea_coordy*(dim(amodel)[1])
  fovea_coords <- c(fovea_coordy, fovea_coordx)
  #fovea_coords <- c(fovea_coordy_scale, fovea_coordx_scale)
  outputz <- just_scale_and_plot(image_file, 0.8, 5, fovea_height, fovea_coords, 20)
  return(outputz)
})

#saveRDS(evens_scale, paste("/Users/currant/OCTeyes/evens_scale",  strftime(Sys.time(), "%Y%m%d"), sep=""))
saveRDS(evens_scale, paste("/Users/currant/OCTeyes/evens_scale_fix_final",  strftime(Sys.time(), "%Y%m%d"), sep=""))


#############JUST_ROTATE##############
odds_rotate <- lapply(1:dim(all_params_odd_minus)[1], function(x){
  image_file <- all_params_odd_minus[[x,1]]
  fovea_height <- all_params_odd_minus[[x,5]]
  #fovea_height_scale <- fovea_height*(dim(amodel)[1])
  fovea_coordx <- all_params_odd_minus[[x,3]]
  #fovea_coordx_scale <- fovea_coordx*(dim(amodel)[2])
  fovea_coordy <- all_params_odd_minus[[x,4]]
  #fovea_coordy_scale <- fovea_coordy*(dim(amodel)[1])
  fovea_coords <- c(fovea_coordy, fovea_coordx)
  #fovea_coords <- c(fovea_coordy_scale, fovea_coordx_scale)
  outputz <- rotate_and_plot_no_scale(image_file, 0.8, 5, fovea_height, fovea_coords, 20)
  return(outputz)
})

#saveRDS(odds_rotate, paste("/Users/currant/OCTeyes/odds_rotate",  strftime(Sys.time(), "%Y%m%d"), sep=""))
saveRDS(odds_rotate, paste("/Users/currant/OCTeyes/odds_rotate_fix_final",  strftime(Sys.time(), "%Y%m%d"), sep=""))

evens_rotate <- lapply(1:dim(all_params_even_minus)[1], function(x){
  image_file <- all_params_even_minus[[x,1]]
  fovea_height <- all_params_even_minus[[x,5]]
  #fovea_height_scale <- fovea_height*(dim(amodel)[1])
  fovea_coordx <- all_params_even_minus[[x,3]]
  #fovea_coordx_scale <- fovea_coordx*(dim(amodel)[2])
  fovea_coordy <- all_params_even_minus[[x,4]]
  #fovea_coordy_scale <- fovea_coordy*(dim(amodel)[1])
  fovea_coords <- c(fovea_coordy, fovea_coordx)
  #fovea_coords <- c(fovea_coordy_scale, fovea_coordx_scale)
  outputz <- rotate_and_plot_no_scale(image_file, 0.8, 5, fovea_height, fovea_coords, 20)
  return(outputz)
})

#saveRDS(evens_rotate, paste("/Users/currant/OCTeyes/evens_rotate",  strftime(Sys.time(), "%Y%m%d"), sep=""))
saveRDS(evens_rotate, paste("/Users/currant/OCTeyes/evens_rotate_fix_final",  strftime(Sys.time(), "%Y%m%d"), sep=""))

#############JUST_TESSELATE##############
odds_tesselate <- lapply(1:dim(all_params_odd_minus)[1], function(x){
  image_file <- all_params_odd_minus[[x,1]]
  fovea_height <- all_params_odd_minus[[x,5]]
  #fovea_height_scale <- fovea_height*(dim(amodel)[1])
  fovea_coordx <- all_params_odd_minus[[x,3]]
  #fovea_coordx_scale <- fovea_coordx*(dim(amodel)[2])
  fovea_coordy <- all_params_odd_minus[[x,4]]
  #fovea_coordy_scale <- fovea_coordy*(dim(amodel)[1])
  fovea_coords <- c(fovea_coordy, fovea_coordx)
  #fovea_coords <- c(fovea_coordy_scale, fovea_coordx_scale)
  outputz <- just_tessalte_and_plot(image_file, 0.8, 5, fovea_height, fovea_coords, 20)
  return(outputz)
})

#saveRDS(odds_tesselate, paste("/Users/currant/OCTeyes/odds_tesselate",  strftime(Sys.time(), "%Y%m%d"), sep=""))
saveRDS(odds_tesselate, paste("/Users/currant/OCTeyes/odds_tesselate_fix_final",  strftime(Sys.time(), "%Y%m%d"), sep=""))

evens_tesselate <- lapply(1:dim(all_params_even_minus)[1], function(x){
  image_file <- all_params_even_minus[[x,1]]
  fovea_height <- all_params_even_minus[[x,5]]
  #fovea_height_scale <- fovea_height*(dim(amodel)[1])
  fovea_coordx <- all_params_even_minus[[x,3]]
  #fovea_coordx_scale <- fovea_coordx*(dim(amodel)[2])
  fovea_coordy <- all_params_even_minus[[x,4]]
  #fovea_coordy_scale <- fovea_coordy*(dim(amodel)[1])
  fovea_coords <- c(fovea_coordy, fovea_coordx)
  #fovea_coords <- c(fovea_coordy_scale, fovea_coordx_scale)
  outputz <- just_tessalte_and_plot(image_file, 0.8, 5, fovea_height, fovea_coords, 20)
  return(outputz)
})

#saveRDS(evens_tesselate, paste("/Users/currant/OCTeyes/evens_tesselate",  strftime(Sys.time(), "%Y%m%d"), sep=""))
saveRDS(evens_tesselate, paste("/Users/currant/OCTeyes/evens_tesselate_fix_final",  strftime(Sys.time(), "%Y%m%d"), sep=""))


############ADD STANDARD TRIANGLES TO DIFFERENTLY PROCESSED IMAGE DATA##############
#add_standard triangles to all processed odd images
oods_minus_stdtri <- lapply(odds, function(x){
  end_info <- x$end_info2
  verts <- x$verticesdf
  plus_tris <- add_standard_triangle(end_info, verts)
  return(plus_tris)
})

# add standard triangles to all processed even images
evens_minus_stdtri <- lapply(evens, function(x){
  end_info <- x$end_info2
  verts <- x$verticesdf
  plus_tris <- add_standard_triangle(end_info, verts)
  return(plus_tris)
})

#add standard triangles to scaled odd images
odds_scale_stdtri <- lapply(odds_scale, function(x){
  end_info <- x$end_info2
  verts <- x$verticesdf
  plus_tris <- add_standard_triangle(end_info, verts)
  return(plus_tris)
})

#add standard triangles to scaled even images
evens_scale_stdtri <- lapply(evens_scale, function(x){
  end_info <- x$end_info2
  verts <- x$verticesdf
  plus_tris <- add_standard_triangle(end_info, verts)
  return(plus_tris)
})

#add standard triangles to rotated odd images
odds_rotate_stdtri <- lapply(odds_rotate, function(x){
  end_info <- x$end_info2
  verts <- x$verticesdf
  plus_tris <- add_standard_triangle(end_info, verts)
  return(plus_tris)
})

#add standard triangles to rotated even images
evens_rotate_stdtri <- lapply(evens_rotate, function(x){
  end_info <- x$end_info2
  verts <- x$verticesdf
  plus_tris <- add_standard_triangle(end_info, verts)
  return(plus_tris)
})

#add standard traingles to plain tesselated odd images
odds_tesselate_stdtri <- lapply(odds_tesselate, function(x){
  end_info <- x$end_info2
  verts <- x$verticesdf
  plus_tris <- add_standard_triangle(end_info, verts)
  return(plus_tris)
})

#add standard triangles to plain tesselated even images
evens_tesselate_stdtri <- lapply(evens_tesselate, function(x){
  end_info <- x$end_info2
  verts <- x$verticesdf
  plus_tris <- add_standard_triangle(end_info, verts)
  return(plus_tris)
})

##############GETS THE REDUCED TABLE FOR DIFFERENTLY PROCESSED IMAGE DATA#############
#Get reduced table for completely processed odd images
all_redtab_oddsminus <- lapply(mintri:maxtri, function(x){
  tmp <- lapply(oods_minus_stdtri, get_reduced_table_onetri, tri_num=x)
  tmp<- do.call(rbind, tmp)
  tmp <- cbind(tmp, 1:28)
  return(tmp)
})

#Get reduced table for completely processed even images
all_redtab_evensminus <- lapply(mintri:maxtri, function(x){
  tmp <- lapply(evens_minus_stdtri, get_reduced_table_onetri, tri_num=x)
  tmp<- do.call(rbind, tmp)
  tmp <- cbind(tmp, 1:28)
  return(tmp)
})

#Get reduced table for scaled odd images
all_redtab_oddsscale <- lapply(mintri:maxtri, function(x){
  tmp <- lapply(odds_scale_stdtri, get_reduced_table_onetri, tri_num=x)
  tmp<- do.call(rbind, tmp)
  tmp <- cbind(tmp, 1:28)
  return(tmp)
})

#Get reduced table for scaled even images
all_redtab_evensscale <- lapply(mintri:maxtri, function(x){
  tmp <- lapply(evens_scale_stdtri, get_reduced_table_onetri, tri_num=x)
  tmp<- do.call(rbind, tmp)
  tmp <- cbind(tmp, 1:28)
  return(tmp)
})

#Get reduced table for rotated odd images
all_redtab_oddsrotate <- lapply(mintri:maxtri, function(x){
  tmp <- lapply(odds_rotate_stdtri, get_reduced_table_onetri, tri_num=x)
  tmp<- do.call(rbind, tmp)
  tmp <- cbind(tmp, 1:28)
  return(tmp)
})

#Get reduced table for rotated even images
all_redtab_evensrotate <- lapply(mintri:maxtri, function(x){
  tmp <- lapply(evens_rotate_stdtri, get_reduced_table_onetri, tri_num=x)
  tmp<- do.call(rbind, tmp)
  tmp <- cbind(tmp, 1:28)
  return(tmp)
})

#Get reduced table for tesselated odd images
all_redtab_oddstesselate <- lapply(mintri:maxtri, function(x){
  tmp <- lapply(odds_tesselate_stdtri, get_reduced_table_onetri, tri_num=x)
  tmp<- do.call(rbind, tmp)
  tmp <- cbind(tmp, 1:28)
  return(tmp)
})

#Get reduced table for tesselated even images
all_redtab_evenstesselate <- lapply(mintri:maxtri, function(x){
  tmp <- lapply(evens_tesselate_stdtri, get_reduced_table_onetri, tri_num=x)
  tmp<- do.call(rbind, tmp)
  tmp <- cbind(tmp, 1:28)
  return(tmp)
})



#Do comparisons of each processed odd set of images against similarly processed even images
all_tess_v_tess <- lapply(1:length(triangle_amount), comp_per_tri, comp1=all_redtab_evenstesselate, comp2=all_redtab_oddstesselate)
all_all_v_all <- lapply(1:length(triangle_amount), comp_per_tri, comp1=all_redtab_oddsminus, comp2=all_redtab_evensminus)
all_scale_v_scale <- lapply(1:length(triangle_amount), comp_per_tri, comp1=all_redtab_oddsscale, comp2=all_redtab_evensscale)
all_rotate_v_rotate <- lapply(1:length(triangle_amount), comp_per_tri, comp1=all_redtab_oddsrotate, comp2=all_redtab_evensrotate)


#Get the correlations per triangle for each comparison
all_tess_v_tess_corr <- lapply(1:length(triangle_amount), function(x){
  table <- all_tess_v_tess[[x]]
  correl <- cor(table[,3], table[,4])
  output <- c(x,correl)
})
all_tess_v_tess_corr <- do.call(rbind, all_tess_v_tess_corr)

all_all_v_all_corr <- lapply(1:length(triangle_amount), function(x){
  table <- all_all_v_all[[x]]
  correl <- cor(table[,3], table[,4])
  output <- c(x, correl)
})
all_all_v_all_corr <- do.call(rbind, all_all_v_all_corr)

all_rotate_v_rotate_corr <- lapply(1:length(triangle_amount), function(x){
  table <- all_rotate_v_rotate[[x]]
  correl <- cor(table[,3], table[,4])
  output <- c(x,correl)
})
all_rotate_v_rotate_corr <- do.call(rbind, all_rotate_v_rotate_corr)

all_scale_v_scale_corr <- lapply(1:length(triangle_amount), function(x){
  table <- all_scale_v_scale[[x]]
  correl <- cor(table[,3], table[,4])
  output<- c(x,correl)
})
all_scale_v_scale_corr <- do.call(rbind, all_scale_v_scale_corr)

###########Plotting scatter plots of correlation vs triangle############

pdf("/Users/currant/OCTeyes/graphics_and_diagrams/per_triangle_correlation.pdf")
plot(all_tess_v_tess_corr, xlab="Triangle ID", ylab="Correlation", main="Tesselation vs tesselation per triangle")
plot(all_all_v_all_corr, xlab="Triangle ID", ylab="Correlation", main="All vs all per triangle")
plot(all_scale_v_scale_corr, xlab="Triangle ID", ylab="Correlation", main="Scale vs scale per triangle")
plot(all_rotate_v_rotate_corr, xlab="Triangle ID", ylab="Correlation", main="Rotate vs rotate per triangle")
dev.off()

pdf("/Users/currant/OCTeyes/scale_and_rotate_graphs/per_triangle_correlation_fix.pdf")
plot(all_tess_v_tess_corr, xlab="Triangle ID", ylab="Correlation", main="Tesselation vs tesselation per triangle")
plot(all_all_v_all_corr, xlab="Triangle ID", ylab="Correlation", main="All vs all per triangle")
plot(all_scale_v_scale_corr, xlab="Triangle ID", ylab="Correlation", main="Scale vs scale per triangle")
plot(all_rotate_v_rotate_corr, xlab="Triangle ID", ylab="Correlation", main="Rotate vs rotate per triangle")
dev.off()

############PLotting distribution of correlations ############
pdf("/Users/currant/OCTeyes/per_triangle_correlation_histogram.pdf")
hist(all_tess_v_tess_corr[,2], xlab="Correlation", main="Tesselation vs Tesselation")
hist(all_all_v_all_corr[,2], xlab="Correlation", main="All vs All")
hist(all_scale_v_scale_corr[,2], xlab="Correlation", main="Scale vs Scale")
hist(all_rotate_v_rotate_corr[,2], xlab="Correlation", main="Rotate vs Rotate")
dev.off()

pdf("/Users/currant/OCTeyes/scale_and_rotate_graphs/per_triangle_correlation_histogram_fix.pdf")
hist(all_tess_v_tess_corr[,2], xlab="Correlation", main="Tesselation vs Tesselation")
hist(all_all_v_all_corr[,2], xlab="Correlation", main="All vs All")
hist(all_scale_v_scale_corr[,2], xlab="Correlation", main="Scale vs Scale")
hist(all_rotate_v_rotate_corr[,2], xlab="Correlation", main="Rotate vs Rotate")
dev.off()


############Generate a reduced table for each combination of analyis##############
oods_minus_red_tab <- lapply(oods_minus_stdtri, get_reduced_table)
evens_minus_red_tab <- lapply(evens_minus_stdtri, get_reduced_table)
odds_scale_red_tab <- lapply(odds_scale_stdtri, get_reduced_table)
evens_scale_red_tab <- lapply(evens_scale_stdtri,get_reduced_table)
odds_rotate_red_tab <- lapply(odds_rotate_stdtri, get_reduced_table)
evens_rotate_red_tab <- lapply(evens_rotate_stdtri, get_reduced_table)
odds_tesselate_red_tab <- lapply(odds_tesselate_stdtri, get_reduced_table)
evens_tesselate_red_tab <- lapply(evens_tesselate_stdtri, get_reduced_table)


############Generate a reduced table for each combination for triangle 87 only as a test##########
oods_minus_red_tab_87 <- lapply(oods_minus_stdtri, get_reduced_table_onetri, tri_num=87)
oods_minus_tab_87 <- do.call(rbind, oods_minus_red_tab_87)
oods_minus_tab_87 <- cbind(oods_minus_tab_87, 1:28)

evens_minus_red_tab_87 <- lapply(evens_minus_stdtri, get_reduced_table_onetri, tri_num=87)
evens_minus_tab_87 <- do.call(rbind, evens_minus_red_tab_87)
evens_minus_tab_87 <- cbind(evens_minus_tab_87, 1:28)

odds_scale_red_tab_87 <- lapply(odds_scale_stdtri, get_reduced_table_onetri, tri_num=87)
odds_scale_tab_87 <- do.call(rbind, odds_scale_red_tab_87)
odds_scale_tab_87 <- cbind(odds_scale_tab_87, 1:28)

evens_scale_red_tab_87 <- lapply(evens_scale_stdtri, get_reduced_table_onetri, tri_num=87)
evens_scale_tab_87 <- do.call(rbind, evens_scale_red_tab_87)
evens_scale_tab_87 <- cbind(evens_scale_tab_87, 1:28)

odds_rotate_red_tab_87 <- lapply(odds_rotate_stdtri, get_reduced_table_onetri, tri_num=87)
odds_rotate_tab_87 <- do.call(rbind, odds_rotate_red_tab_87)
odds_rotate_tab_87 <- cbind(odds_rotate_tab_87, 1:28)

evens_rotate_red_tab_87 <- lapply(evens_rotate_stdtri, get_reduced_table_onetri, tri_num=87)
evens_rotate_tab_87 <- do.call(rbind, evens_rotate_red_tab_87)
evens_rotate_tab_87 <- cbind(evens_rotate_tab_87, 1:28)

odds_tesselate_red_tab_87 <- lapply(odds_tesselate_stdtri, get_reduced_table_onetri, tri_num=87)
odds_tesselate_tab_87 <- do.call(rbind, odds_tesselate_red_tab_87)
odds_tesselate_tab_87 <- cbind(odds_tesselate_tab_87, 1:28)

evens_tesselate_red_tab_87 <- lapply(evens_tesselate_stdtri, get_reduced_table_onetri, tri_num=87)
evens_tesselate_tab_87 <- do.call(rbind, evens_tesselate_red_tab_87)
evens_tesselate_tab_87 <- cbind(evens_tesselate_tab_87, 1:28)


############Retroactively mean normailse############
odds_minus_retro <- lapply(1:length(odds_minus), function(x){
  tmp3 <- odds_minus[[x]]$end_info2
  medi <- median(tmp3[,3], na.rm=TRUE)
  tmp3[,3] <- tmp3[,3]-medi
  return(tmp3)
})

evens_minus_retro <- lapply(1:length(evens_minus), function(x){
  tmp3 <- evens_minus[[x]]$end_info2
  medi <- median(tmp3[,3], na.rm=TRUE)
  tmp3[,3] <- tmp3[,3]-medi
  return(tmp3)
})

odds_scale_retro <- lapply(1:length(odds_scale), function(x){
  tmp3 <- odds_scale[[x]]$end_info2
  medi <- median(tmp3[,3], na.rm=TRUE)
  tmp3[,3] <- tmp3[,3]-medi
  return(tmp3)
})

evens_scale_retro <- lapply(1:length(evens_scale), function(x){
  tmp3 <- evens_scale[[x]]$end_info2
  medi <- median(tmp3[,3], na.rm=TRUE)
  tmp3[,3] <- tmp3[,3]-medi
  return(tmp3)
})

odds_rotate_retro <- lapply(1:length(odds_rotate), function(x){
  tmp3 <- odds_rotate[[x]]$end_info2
  medi <- median(tmp3[,3], na.rm=TRUE)
  tmp3[,3] <- tmp3[,3]-medi
  return(tmp3)
})

evens_rotate_retro <- lapply(1:length(evens_rotate), function(x){
  tmp3 <- evens_rotate[[x]]$end_info2
  medi <- median(tmp3[,3], na.rm=TRUE)
  tmp3[,3] <- tmp3[,3]-medi
  return(tmp3)
})

odds_tesselate_retro <- lapply(1:length(odds_tesselate), function(x){
  tmp3 <- odds_tesselate[[x]]$end_info2
  medi <- median(tmp3[,3], na.rm=TRUE)
  tmp3[,3] <- tmp3[,3]-medi
  return(tmp3)
})

evens_tesselate_retro <- lapply(1:length(evens_tesselate), function(x){
  tmp3 <- evens_tesselate[[x]]$end_info2
  medi <- median(tmp3[,3], na.rm=TRUE)
  tmp3[,3] <- tmp3[,3]-medi
  return(tmp3)
})

#add_standard triangles to all processed odd images
oods_minus_retro_stdtri <- lapply(1:length(odds_minus_retro), function(x, odds_minus){
  end_info <- odds_minus_retro[[x]]
  verts <- odds_minus[[x]]$verticesdf
  plus_tris <- add_standard_triangle(end_info, verts)
  return(plus_tris)
}, odds_minus=odds_minus)

# add standard triangles to all processed even images
evens_minus_retro_stdtri <- lapply(1:length(evens_minus_retro), function(x, evens_minus){
  end_info <- evens_minus_retro[[x]]
  verts <- evens_minus[[x]]$verticesdf
  plus_tris <- add_standard_triangle(end_info, verts)
  return(plus_tris)
}, evens_minus=evens_minus)

#add standard triangles to scaled odd images
odds_scale_retro_stdtri <- lapply(1:length(odds_scale_retro), function(x, odds_scale){
  end_info <- odds_scale_retro[[x]]
  verts <- odds_scale[[x]]$verticesdf
  plus_tris <- add_standard_triangle(end_info, verts)
  return(plus_tris)
}, odds_scale=odds_scale)

#add standard triangles to scaled even images
evens_scale_retro_stdtri <- lapply(1:length(evens_scale_retro), function(x, evens_scale){
  end_info <- evens_scale_retro[[x]]
  verts <- evens_scale[[x]]$verticesdf
  plus_tris <- add_standard_triangle(end_info, verts)
  return(plus_tris)
}, evens_scale=evens_scale)

#add standard triangles to rotated odd images
odds_rotate_retro_stdtri <- lapply(1:length(odds_rotate_retro), function(x, odds_rotate){
  end_info <- odds_rotate_retro[[x]]
  verts <- odds_rotate[[x]]$verticesdf
  plus_tris <- add_standard_triangle(end_info, verts)
  return(plus_tris)
}, odds_rotate=odds_rotate)

#add standard triangles to rotated even images
evens_rotate_retro_stdtri <- lapply(1:length(evens_rotate_retro), function(x, evens_rotate){
  end_info <- evens_rotate_retro[[x]]
  verts <- evens_rotate[[x]]$verticesdf
  plus_tris <- add_standard_triangle(end_info, verts)
  return(plus_tris)
}, evens_rotate=evens_rotate)

#add standard traingles to plain tesselated odd images
odds_tesselate_retro_stdtri <- lapply(1:length(odds_tesselate_retro), function(x, odds_tesselate){
  end_info <- odds_tesselate_retro[[x]]
  verts <- odds_tesselate[[x]]$verticesdf
  plus_tris <- add_standard_triangle(end_info, verts)
  return(plus_tris)
}, odds_tesselate=odds_tesselate)

#add standard triangles to plain tesselated even images
evens_tesselate_retro_stdtri <- lapply(1:length(evens_tesselate_retro), function(x, evens_tesselate){
  end_info <- evens_tesselate_retro[[x]]
  verts <- evens_tesselate[[x]]$verticesdf
  plus_tris <- add_standard_triangle(end_info, verts)
  return(plus_tris)
}, evens_tesselate=evens_tesselate)


all_redtab_oddsminus_retro <- lapply(mintri:maxtri, function(x){
  tmp <- lapply(oods_minus_retro_stdtri, get_reduced_table_onetri, tri_num=x)
  tmp<- do.call(rbind, tmp)
  tmp <- cbind(tmp, 1:28)
  return(tmp)
})

all_redtab_evensminus_retro <- lapply(mintri:maxtri, function(x){
  tmp <- lapply(evens_minus_retro_stdtri, get_reduced_table_onetri, tri_num=x)
  tmp<- do.call(rbind, tmp)
  tmp <- cbind(tmp, 1:28)
  return(tmp)
})

all_redtab_oddsscale_retro <- lapply(mintri:maxtri, function(x){
  tmp <- lapply(odds_scale_retro_stdtri, get_reduced_table_onetri, tri_num=x)
  tmp<- do.call(rbind, tmp)
  tmp <- cbind(tmp, 1:28)
  return(tmp)
})

all_redtab_evensscale_retro <- lapply(mintri:maxtri, function(x){
  tmp <- lapply(evens_scale_retro_stdtri, get_reduced_table_onetri, tri_num=x)
  tmp<- do.call(rbind, tmp)
  tmp <- cbind(tmp, 1:28)
  return(tmp)
})

all_redtab_oddsrotate_retro <- lapply(mintri:maxtri, function(x){
  tmp <- lapply(odds_rotate_retro_stdtri, get_reduced_table_onetri, tri_num=x)
  tmp<- do.call(rbind, tmp)
  tmp <- cbind(tmp, 1:28)
  return(tmp)
})

all_redtab_evensrotate_retro <- lapply(mintri:maxtri, function(x){
  tmp <- lapply(evens_rotate_retro_stdtri, get_reduced_table_onetri, tri_num=x)
  tmp<- do.call(rbind, tmp)
  tmp <- cbind(tmp, 1:28)
  return(tmp)
})

all_redtab_oddstesselate_retro <- lapply(mintri:maxtri, function(x){
  tmp <- lapply(odds_tesselate_retro_stdtri, get_reduced_table_onetri, tri_num=x)
  tmp<- do.call(rbind, tmp)
  tmp <- cbind(tmp, 1:28)
  return(tmp)
})

 

all_redtab_evenstesselate_retro <- lapply(mintri:maxtri, function(x){
  tmp <- lapply(evens_tesselate_retro_stdtri, get_reduced_table_onetri, tri_num=x)
  tmp<- do.call(rbind, tmp)
  tmp <- cbind(tmp, 1:28)
  return(tmp)
})

all_tess_v_tess_retro <- lapply(1:length(triangle_amount), comp_per_tri, comp1=all_redtab_evenstesselate_retro, comp2=all_redtab_oddstesselate_retro)
all_all_v_all_retro <- lapply(1:length(triangle_amount), comp_per_tri, comp1=all_redtab_oddsminus_retro, comp2=all_redtab_evensminus_retro)
all_scale_v_scale_retro <- lapply(1:length(triangle_amount), comp_per_tri, comp1=all_redtab_oddsscale_retro, comp2=all_redtab_evensscale_retro)
all_rotate_v_rotate_retro <- lapply(1:length(triangle_amount), comp_per_tri, comp1=all_redtab_oddsrotate_retro, comp2=all_redtab_evensrotate_retro)

all_tess_v_tess_retro_corr <- lapply(1:length(triangle_amount), function(x){
  table <- all_tess_v_tess_retro[[x]]
  correl <- cor(table[,3], table[,4])
  output <- c(x,correl)
})
all_tess_v_tess_retro_corr <- do.call(rbind, all_tess_v_tess_retro_corr)

all_all_v_all_retro_corr <- lapply(1:length(triangle_amount), function(x){
  table <- all_all_v_all_retro[[x]]
  correl <- cor(table[,3], table[,4])
  output <- c(x, correl)
})
all_all_v_all_retro_corr <- do.call(rbind, all_all_v_all_retro_corr)

all_rotate_v_rotate_retro_corr <- lapply(1:length(triangle_amount), function(x){
  table <- all_rotate_v_rotate_retro[[x]]
  correl <- cor(table[,3], table[,4])
  output <- c(x,correl)
})
all_rotate_v_rotate_retro_corr <- do.call(rbind, all_rotate_v_rotate_retro_corr)

all_scale_v_scale_retro_corr <- lapply(1:length(triangle_amount), function(x){
  table <- all_scale_v_scale_retro[[x]]
  correl <- cor(table[,3], table[,4])
  output<- c(x,correl)
})
all_scale_v_scale_retro_corr <- do.call(rbind, all_scale_v_scale_retro_corr)

###########Plotting scatter plots of correlation vs triangle############

pdf("/Users/currant/OCTeyes/graphics_and_diagrams/per_triangle_mean_norm_correlation.pdf")
plot(all_tess_v_tess_retro_corr, xlab="Triangle ID", ylab="Correlation", main="Tesselation vs tesselation per triangle")
plot(all_all_v_all_retro_corr, xlab="Triangle ID", ylab="Correlation", main="All vs all per triangle")
plot(all_scale_v_scale_retro_corr, xlab="Triangle ID", ylab="Correlation", main="Scale vs scale per triangle")
plot(all_rotate_v_rotate_retro_corr, xlab="Triangle ID", ylab="Correlation", main="Rotate vs rotate per triangle")
dev.off()


############PLotting distribution of correlations ############
pdf("/Users/currant/OCTeyes/per_triangle_mean_norm_correlation_histogram.pdf")
hist(all_tess_v_tess_retro_corr[,2], xlab="Correlation", main="Tesselation vs Tesselation")
hist(all_all_v_all_retro_corr[,2], xlab="Correlation", main="All vs All")
hist(all_scale_v_scale_retro_corr[,2], xlab="Correlation", main="Scale vs Scale")
hist(all_rotate_v_rotate_retro_corr[,2], xlab="Correlation", main="Rotate vs Rotate")
dev.off()

###########Plotting scatter plots of correlation vs triangle############

pdf("/Users/currant/OCTeyes/scale_and_rotate_graphs/final_plots_we_hope/per_triangle_mean_norm_correlation_fix.pdf")
plot(all_tess_v_tess_retro_corr, xlab="Triangle ID", ylab="Correlation", main="Tesselation vs tesselation per triangle")
plot(all_all_v_all_retro_corr, xlab="Triangle ID", ylab="Correlation", main="All vs all per triangle")
plot(all_scale_v_scale_retro_corr, xlab="Triangle ID", ylab="Correlation", main="Scale vs scale per triangle")
plot(all_rotate_v_rotate_retro_corr, xlab="Triangle ID", ylab="Correlation", main="Rotate vs rotate per triangle")
dev.off()


############PLotting distribution of correlations ############
pdf("/Users/currant/OCTeyes/scale_and_rotate_graphs/final_plots_we_hope/per_triangle_mean_norm_correlation_histogram_fix.pdf")
hist(all_tess_v_tess_retro_corr[,2], xlab="Correlation", main="Tesselation vs Tesselation")
hist(all_all_v_all_retro_corr[,2], xlab="Correlation", main="All vs All")
hist(all_scale_v_scale_retro_corr[,2], xlab="Correlation", main="Scale vs Scale")
hist(all_rotate_v_rotate_retro_corr[,2], xlab="Correlation", main="Rotate vs Rotate")
dev.off()


############Generate a reduced table for each combination for triangle 87 only as a test##########
oods_minus_red_tab_retro_87 <- lapply(oods_minus_retro_stdtri, get_reduced_table_onetri, tri_num=87)
oods_minus_tab_retro_87 <- do.call(rbind, oods_minus_red_tab_retro_87)
oods_minus_tab_retro_87 <- cbind(oods_minus_tab_retro_87, 1:28)

evens_minus_red_tab_retro_87 <- lapply(evens_minus_retro_stdtri, get_reduced_table_onetri, tri_num=87)
evens_minus_tab_retro_87 <- do.call(rbind, evens_minus_red_tab_retro_87)
evens_minus_tab_retro_87 <- cbind(evens_minus_tab_retro_87, 1:28)

odds_scale_red_tab_retro_87 <- lapply(odds_scale_retro_stdtri, get_reduced_table_onetri, tri_num=87)
odds_scale_tab_retro_87 <- do.call(rbind, odds_scale_red_tab_retro_87)
odds_scale_tab_retro_87 <- cbind(odds_scale_tab_retro_87, 1:28)

evens_scale_red_tab_retro_87 <- lapply(evens_scale_retro_stdtri, get_reduced_table_onetri, tri_num=87)
evens_scale_tab_retro_87 <- do.call(rbind, evens_scale_red_tab_retro_87)
evens_scale_tab_retro_87 <- cbind(evens_scale_tab_retro_87, 1:28)

odds_rotate_red_tab_retro_87 <- lapply(odds_rotate_retro_stdtri, get_reduced_table_onetri, tri_num=87)
odds_rotate_tab_retro_87 <- do.call(rbind, odds_rotate_red_tab_retro_87)
odds_rotate_tab_retro_87 <- cbind(odds_rotate_tab_retro_87, 1:28)

evens_rotate_red_tab_retro_87 <- lapply(evens_rotate_retro_stdtri, get_reduced_table_onetri, tri_num=87)
evens_rotate_tab_retro_87 <- do.call(rbind, evens_rotate_red_tab_retro_87)
evens_rotate_tab_retro_87 <- cbind(evens_rotate_tab_retro_87, 1:28)

odds_tesselate_red_tab_retro_87 <- lapply(odds_tesselate_retro_stdtri, get_reduced_table_onetri, tri_num=87)
odds_tesselate_tab_retro_87 <- do.call(rbind, odds_tesselate_red_tab_retro_87)
odds_tesselate_tab_retro_87 <- cbind(odds_tesselate_tab_retro_87, 1:28)

evens_tesselate_red_tab_retro_87 <- lapply(evens_tesselate_retro_stdtri, get_reduced_table_onetri, tri_num=87)
evens_tesselate_tab_retro_87 <- do.call(rbind, evens_tesselate_red_tab_retro_87)
evens_tesselate_tab_retro_87 <- cbind(evens_tesselate_tab_retro_87, 1:28)

tess_v_tess_retro_87 <- data.frame(odds_tesselate_tab_retro_87[,3], odds_tesselate_tab_retro_87[,2], evens_tesselate_tab_retro_87[,2])
tess_v_tess_retro_cor_87 <- cor(tess_v_tess_retro_87[2:3])[1,2]
minus_v_minus_retro_87 <- data.frame(oods_minus_tab_retro_87[,3], oods_minus_tab_retro_87[,2], evens_minus_tab_retro_87[,2])
minus_v_minus_retro_cor_87 <- cor(minus_v_minus_retro_87[,2:3])[1,2]
scale_v_scale_retro_87 <- data.frame(odds_scale_tab_retro_87[,3], odds_scale_tab_retro_87[,2], evens_scale_tab_retro_87[,2])
scale_v_scale_retro_cors_87 <- cor(scale_v_scale_retro_87[,2:3])[1,2]
rotate_v_rotate_retro_87 <- data.frame(odds_rotate_tab_retro_87[,3], odds_rotate_tab_retro_87[,2], evens_rotate_tab_retro_87[,2])
rotate_v_rotate_retro_cor_87 <- cor(rotate_v_rotate_retro_87[,2:3])[1,2]

pdf("/Users/currant/OCTeyes/87_retro_correl.pdf")
plot(tess_v_tess_retro_87[,2:3], xlab='Image repeat 1', ylab='Image repeat 2', main='tesselation vs tesselation')
mtext(paste("correlation=", tess_v_tess_retro_cor_87[1,2], sep=""), 4)
abline(lm(tess_v_tess_retro_87[,2]~tess_v_tess_retro_87[,3]))

plot(scale_v_scale_retro_87[,2:3], xlab='Image repeat 1', ylab='Image repeat 2', main='scale vs scale')
mtext(paste("correlation=", scale_v_scale_retro_cor_87, sep=""), 4)
abline(lm(scale_v_scale_retro_87[,2]~scale_v_scale_retro_87[,3]))

plot(rotate_v_rotate_retro_87[,2:3], xlab='Image repeat 1', ylab='Image repeat 2', main='rotate vs rotate')
mtext(paste("correlation=", rotate_v_rotate_retro_cor_87, sep=""), 4)
abline(lm(rotate_v_rotate_retro_87[,2]~rotate_v_rotate_retro_87[,3]))

plot(minus_v_minus_retro_87[,2:3], xlab='Image repeat 1', ylab='Image repeat 2', main='all vs all')
mtext(paste("correlation=", minus_v_minus_retro_cor_87, sep=""),4)
abline(lm(minus_v_minus_retro_87[,2]~minus_v_minus_retro_87[,3]))


##########PLOT CORRELATIONS BACK TO TRIANGLES##################

ref_points <- evens_tesselate[[1]]$all_refs

get_triplot_table <- function(x, ref_points, corr_tab){
  num_pos <- which(ref_points[,3]==x)
  small_tab <- ref_points[num_pos,]
  small_tab <- cbind(small_tab, corr_tab[x,2])
  return(small_tab)
}
tess_v_tess_cor_triplot <- lapply(1:dim(all_tess_v_tess_retro_corr)[1], function(x){
  num_pos <- which(ref_points[,3]==x)
  small_tab <- ref_points[num_pos,]
  small_tab <- cbind(small_tab, all_tess_v_tess_retro_corr[x,2])
  return(small_tab)
})

tess_v_tess_cor_triplot <- do.call(rbind, tess_v_tess_cor_triplot)
tess_v_tess_cor_triplot <- as.data.frame(tess_v_tess_cor_triplot)

p <- ggplot(data=tess_v_tess_cor_triplot)
p + 
  geom_point(aes(x=V1, y=V2, col=V4), na.rm=TRUE) + 
  theme_bw()

scale_v_scale_cor_triplot <- lapply(1:dim(all_scale_v_scale_retro_corr)[1], function(x){
  num_pos <- which(ref_points[,3]==x)
  small_tab <- ref_points[num_pos,]
  small_tab <- cbind(small_tab, all_scale_v_scale_retro_corr[x,2])
  return(small_tab)
})

scale_v_scale_cor_triplot <- do.call(rbind, scale_v_scale_cor_triplot)
scale_v_scale_cor_triplot <- as.data.frame(scale_v_scale_cor_triplot)

p <- ggplot(data=scale_v_scale_cor_triplot)
p + 
  geom_point(aes(x=V1, y=V2, col=V4), na.rm=TRUE) + 
  theme_bw()

rotate_v_rotate_cor_triplot <- lapply(1:dim(all_rotate_v_rotate_retro_corr)[1], function(x){
  num_pos <- which(ref_points[,3]==x)
  small_tab <- ref_points[num_pos,]
  small_tab <- cbind(small_tab, all_rotate_v_rotate_retro_corr[x,2])
  return(small_tab)
})

rotate_v_rotate_cor_triplot <- do.call(rbind, rotate_v_rotate_cor_triplot)
rotate_v_rotate_cor_triplot <- as.data.frame(rotate_v_rotate_cor_triplot)

p <- ggplot(data=rotate_v_rotate_cor_triplot)
p + 
  geom_point(aes(x=V1, y=V2, col=V4), na.rm=TRUE) + 
  theme_bw()

all_v_all_cor_triplot <- lapply(1:dim(all_all_v_all_retro_corr)[1], function(x){
  num_pos <- which(ref_points[,3]==x)
  small_tab <- ref_points[num_pos,]
  small_tab <- cbind(small_tab, all_all_v_all_retro_corr[x,2])
  return(small_tab)
})

all_v_all_cor_triplot <- do.call(rbind, all_v_all_cor_triplot)
all_v_all_cor_triplot <- as.data.frame(all_v_all_cor_triplot)

p <- ggplot(data=all_v_all_cor_triplot)
p + 
  geom_point(aes(x=V1, y=V2, col=V4), na.rm=TRUE) + 
  theme_bw()

#############Gets the amount of triangles available, only works for 
min_and_max <- lapply(oods_minus_red_tab, function(x){
  minimumval <- min(x[,1])
  maximumval <- max(x[,1])
  maxandmin <- c(minimumval, maximumval)
  return(maxandmin)
})
max_and_min <- unlist(min_and_max)
maxtri <- max(max_and_min)
mintri <- min(max_and_min)
maxtri <- 159
mintri <- 0


triangle_amount <- 0:159



comp_function <- function(x, odd_tab, even_tab){
  odds <- odd_tab[[x]]
  evens <- even_tab[[x]]
  o_tris <-odds[,1]
  e_tris <- evens[,1]
  shared_tris <- o_tris[o_tris%in%e_tris]
  shared_df <- lapply(shared_tris, function(x){
    odd_val <- odds[which(odds[,1]==x),2]
    even_val <- evens[which(evens[,1]==x),2]
    return(c(x,odd_val,even_val))
  })
  shared_df <- do.call(rbind, shared_df)
  shared_df <- cbind(shared_df,x)
  return(shared_df)
}
  
tess_v_tess <- lapply(1:(length(odds_tesselate_red_tab)), function(x){
  odds <- odds_tesselate_red_tab[[x]]
  evens <- evens_tesselate_red_tab[[x]]
  o_tris <- odds[,1]
  e_tris <- evens[,1]
  shared_tris <- o_tris[o_tris%in%e_tris]
  shared_df <- lapply( shared_tris, function(x){
    odd_val <- odds[which(odds[,1]==x),2]
    even_val <- evens[which(evens[,1]==x),2]
    return(c(x,odd_val, even_val))
  })
  shared_df <- do.call(rbind, shared_df)
  shared_df <- cbind(shared_df, x)
  return(shared_df)
})
tess_v_tess <- do.call(rbind, tess_v_tess)

tess_v_tess_87 <- data.frame(odds_tesselate_tab_87[,3], odds_tesselate_tab_87[,2], evens_tesselate_tab_87[,2])
tess_v_tess_cor_87 <- cor(tess_v_tess_87[2:3])[1,2]
minus_v_minus_87 <- data.frame(oods_minus_tab_87[,3], oods_minus_tab_87[,2], evens_minus_tab_87[,2])
minus_v_minus_cor_87 <- cor(minus_v_minus_87[,2:3])[1,2]
scale_v_scale_87 <- data.frame(odds_scale_tab_87[,3], odds_scale_tab_87[,2], evens_scale_tab_87[,2])
scale_v_scale_cors_87 <- cor(scale_v_scale_87[,2:3])[1,2]
rotate_v_rotate_87 <- data.frame(odds_rotate_tab_87[,3], odds_rotate_tab_87[,2], evens_rotate_tab_87[,2])
rotate_v_rotate_cor_87 <- cor(rotate_v_rotate_87[,2:3])[1,2]


######PLOTTING 87
plot(tess_v_tess_87[,2:3], xlab='Image repeat 1', ylab='Image repeat 2', main='tesselation vs tesselation')
mtext(paste("correlation=", tess_v_tess_cor_87[1,2], sep=""), 4)
abline(lm(tess_v_tess_87[,2]~tess_v_tess_87[,3]))

plot(scale_v_scale_87[,2:3], xlab='Image repeat 1', ylab='Image repeat 2', main='scale vs scale')
mtext(paste("correlation=", scale_v_scale_cors_87, sep=""), 4)
abline(lm(scale_v_scale_87[,2]~scale_v_scale_87[,3]))

plot(rotate_v_rotate_87[,2:3], xlab='Image repeat 1', ylab='Image repeat 2', main='rotate vs rotate')
mtext(paste("correlation=", rotate_v_rotate_cor_87, sep=""), 4)
abline(lm(rotate_v_rotate_87[,2]~rotate_v_rotate_87[,3]))

plot(minus_v_minus_87[,2:3], xlab='Image repeat 1', ylab='Image repeat 2', main='all vs all')
mtext(paste("correlation=", minus_v_minus_cor_87, sep=""),4)
abline(lm(minus_v_minus_87[,2]~minus_v_minus_87[,3]))

splt_tess_v_tess <- split(as.data.frame(tess_v_tess), as.factor(tess_v_tess[,4]))
tess_v_tess_cors <- lapply(1:length(splt_tess_v_tess), function(x){
  data <- splt_tess_v_tess[[x]]
  correl <- cor(data[,2:3])[1,2]
  mean_val_one <- mean(data[,2])
  mean_val_two <- mean(data[,3])
  median_val_one <- median(data[,2])
  median_val_two <- median(data[,3])
  output <- c(x, correl, mean_val_one, mean_val_two, median_val_one, median_val_two)
  #output <- c(x, correl)
  return(output)
})
tess_v_tess_cors <- do.call(rbind, tess_v_tess_cors)

tess_v_scale <- lapply(1:(length(odds_tesselate_red_tab)), function(x){
  odds <- odds_tesselate_red_tab[[x]]
  evens <- evens_scale_red_tab[[x]]
  o_tris <- odds[,1]
  e_tris <- evens[,1]
  shared_tris <- o_tris[o_tris%in%e_tris]
  shared_df <- lapply( shared_tris, function(x){
    odd_val <- odds[which(odds[,1]==x),2]
    even_val <- evens[which(evens[,1]==x),2]
    return(c(x,odd_val, even_val))
  })
  shared_df <- do.call(rbind, shared_df)
  shared_df <- cbind(shared_df, x)
  return(shared_df)
})
tess_v_scale <- do.call(rbind, tess_v_scale)

splt_tess_v_scale <- split(as.data.frame(tess_v_scale), as.factor(tess_v_scale[,4]))
tess_v_scale_cors <- lapply(1:length(splt_tess_v_scale), function(x){
  data <- splt_tess_v_scale[[x]]
  correl <- cor(data[,2:3])[1,2]
  output <- c(x, correl)
  return(output)
})
tess_v_scale_cors <- do.call(rbind, tess_v_scale_cors)

tess_v_rotate <- lapply(1:(length(odds_tesselate_red_tab)), function(x){
  odds <- odds_tesselate_red_tab[[x]]
  evens <- evens_rotate_red_tab[[x]]
  o_tris <- odds[,1]
  e_tris <- evens[,1]
  shared_tris <- o_tris[o_tris%in%e_tris]
  shared_df <- lapply( shared_tris, function(x){
    odd_val <- odds[which(odds[,1]==x),2]
    even_val <- evens[which(evens[,1]==x),2]
    return(c(x,odd_val, even_val))
  })
  shared_df <- do.call(rbind, shared_df)
  shared_df <- cbind(shared_df, x)
  return(shared_df)
})
tess_v_rotate <- do.call(rbind, tess_v_rotate)

splt_tess_v_rotate <- split(as.data.frame(tess_v_rotate), as.factor(tess_v_rotate[,4]))
tess_v_rotate_cors <- lapply(1:length(splt_tess_v_rotate), function(x){
  data <- splt_tess_v_rotate[[x]]
  correl <- cor(data[,2:3])[1,2]
  output <- c(x, correl)
  return(output)
})
tess_v_rotate_cors <- do.call(rbind, tess_v_rotate_cors)

tess_v_minus <- lapply(1:(length(odds_tesselate_red_tab)), function(x){
  odds <- odds_tesselate_red_tab[[x]]
  evens <- evens_minus_red_tab[[x]]
  o_tris <- odds[,1]
  e_tris <- evens[,1]
  shared_tris <- o_tris[o_tris%in%e_tris]
  shared_df <- lapply( shared_tris, function(x){
    odd_val <- odds[which(odds[,1]==x),2]
    even_val <- evens[which(evens[,1]==x),2]
    return(c(x,odd_val, even_val))
  })
  shared_df <- do.call(rbind, shared_df)
  shared_df <- cbind(shared_df, x)
  return(shared_df)
})
tess_v_minus <- do.call(rbind, tess_v_minus)

splt_tess_v_minus <- split(as.data.frame(tess_v_minus), as.factor(tess_v_minus[,4]))
tess_v_minus_cors <- lapply(1:length(splt_tess_v_minus), function(x){
  data <- splt_tess_v_minus[[x]]
  correl <- cor(data[,2:3])[1,2]
  output <- c(x, correl)
  return(output)
})
tess_v_minus_cors <- do.call(rbind, tess_v_minus_cors)


rotate_v_rotate <- lapply(1:(length(odds_rotate_red_tab)), function(x){
  odds <- odds_rotate_red_tab[[x]]
  evens <- evens_rotate_red_tab[[x]]
  o_tris <- odds[,1]
  e_tris <- evens[,1]
  shared_tris <- o_tris[o_tris%in%e_tris]
  shared_df <- lapply( shared_tris, function(x){
    odd_val <- odds[which(odds[,1]==x),2]
    even_val <- evens[which(evens[,1]==x),2]
    return(c(x,odd_val, even_val))
  })
  shared_df <- do.call(rbind, shared_df)
  shared_df <- cbind(shared_df, x)
  return(shared_df)
})
rotate_v_rotate <- do.call(rbind, rotate_v_rotate)

splt_rotate_v_rotate <- split(as.data.frame(rotate_v_rotate), as.factor(rotate_v_rotate[,4]))
rotate_v_rotate_cors <- lapply(1:length(splt_rotate_v_rotate), function(x){
  data <- splt_rotate_v_rotate[[x]]
  correl <- cor(data[,2:3])[1,2]
  output <- c(x, correl)
  return(output)
})
rotate_v_rotate_cors <- do.call(rbind, rotate_v_rotate_cors)

scale_v_scale <- lapply(1:(length(odds_scale_red_tab)), function(x){
  odds <- odds_scale_red_tab[[x]]
  evens <- evens_scale_red_tab[[x]]
  o_tris <- odds[,1]
  e_tris <- evens[,1]
  shared_tris <- o_tris[o_tris%in%e_tris]
  shared_df <- lapply( shared_tris, function(x){
    odd_val <- odds[which(odds[,1]==x),2]
    even_val <- evens[which(evens[,1]==x),2]
    return(c(x,odd_val, even_val))
  })
  shared_df <- do.call(rbind, shared_df)
  shared_df <- cbind(shared_df, x)
  return(shared_df)
})
scale_v_scale <- do.call(rbind, scale_v_scale)

splt_scale_v_scale <- split(as.data.frame(scale_v_scale), as.factor(scale_v_scale[,4]))
scale_v_scale_cors <- lapply(1:length(splt_scale_v_scale), function(x){
  data <- splt_tess_v_minus[[x]]
  correl <- cor(data[,2:3])[1,2]
  output <- c(x, correl)
  return(output)
})
scale_v_scales_cors <- do.call(rbind, scale_v_scale_cors)


minus_v_minus <- lapply(1:(length(oods_minus_red_tab)), function(x){
  odds <- oods_minus_red_tab[[x]]
  evens <- evens_minus_red_tab[[x]]
  o_tris <- odds[,1]
  e_tris <- evens[,1]
  shared_tris <- o_tris[o_tris%in%e_tris]
  shared_df <- lapply( shared_tris, function(x){
    odd_val <- odds[which(odds[,1]==x),2]
    even_val <- evens[which(evens[,1]==x),2]
    return(c(x,odd_val, even_val))
  })
  shared_df <- do.call(rbind, shared_df)
  shared_df <- cbind(shared_df, x)
  return(shared_df)
})
minus_v_minus <- do.call(rbind, minus_v_minus)

splt_minus_v_minus <- split(as.data.frame(minus_v_minus), as.factor(minus_v_minus[,4]))
minus_v_minus_cors <- lapply(1:length(splt_minus_v_minus), function(x){
  data <- splt_minus_v_minus[[x]]
  correl <- cor(data[,2:3])[1,2]
  #mean_val_one <- mean(data[,2])
  #mean_val_two <- mean(data[,3])
  #output <- c(x, correl, mean_val_one, mean_val_two)
  output <- c(x, correl)
  return(output)
})
minus_v_minus_cors <- do.call(rbind, minus_v_minus_cors)

correl_comp <- data.frame(tess_v_tess_cors, tess_v_scale_cors[,2], tess_v_rotate_cors[,2], tess_v_minus_cors[,2], rotate_v_rotate_cors[,2], minus_v_minus_cors[,2])
correl_comp <- cbind(correl_comp, abs(tess_v_tess_cors[,5] - tess_v_tess_cors[,6]))
correl_comp_melt_median <- melt(correl_comp, id=c('X1', 'abs(tess_v_tess_cors[, 5] - tess_v_tess_cors[, 6])'))
ggplot(melted_correl_comps_pixel_val, aes(x=correl_comp_melt_median[,1], y=correl_comp_melt_median[,4], col=correl_comp_melt_median[,2])) + geom_point() + scale_color_gradientn(colours= rainbow(5))






############Quality control###############
##########################################

min_and_max <- lapply(reduced_tables, function(x){
  minimumval <- min(x[,1])
  maximumval <- max(x[,1])
  maxandmin <- c(minimumval, maximumval)
  return(maxandmin)
})
max_and_min <- unlist(min_and_max)
maxtri <- max(max_and_min)
mintri <- min(max_and_min)

triangle_amount <- mintri:maxtri

big_table_rotate <- lapply(reduced_tables, function(x){
  table <- x
  table <- as.data.frame(table)
  per_tri <- lapply(triangle_amount, function(x){
    tri_le_tri <- which(table[,1]==x)
    if (length(tri_le_tri)==0){
      tri_val <- NA
    }
    else{
      tri_val <- table[tri_le_tri,2]
    }
    #tri_le_tri <- filter(table, table[,1]==x)
    #tri_val <- tri_le_tri[,2]
    return(tri_val)
  })
  lineyline <- do.call(cbind, per_tri)
  return(lineyline)
})

biggy_table <- do.call(rbind, big_table_rotate)

############WORKING UP TO HERE###################

all_eye_paths <- get_eye_paths(pathname)
all_eye_params10 <- get_needed_params(all_eye_paths[1:10])
all_eye_paths_odd <- get_eye_pathsrepeats(pathname, 1)
all_params_odd <- get_needed_params(all_eye_paths_odd)
all_eye_paths_even <- get_eye_pathsrepeats(pathname, 2)
all_params_even <- get_needed_params(all_eye_paths_even)

all_params_odd_minus <- all_params_odd[-c(6),]
all_params_odd_minus <- all_params_odd_minus[-c(24),]
all_params_even_minus <- all_params_even[-c(6),]
all_params_even_minus <- all_params_even_minus[-c(24),]


##################

rotation_turns= list(0,20,40,60,80,100,120,140,160,180,200,220,240,260,280,300,320,340,360)
turning_outcome <- lapply(rotation_turns, rotate_and_plot, image_file=image_file, centre_f=0.8, triangle_num=5, fovea_height=fovea_height, fovea_coords=fovea_coords, n=n )

plot_all_rotations <- function(rotation){
  end_info2 <- rotation[[1]]
  verticesdf <- rotation[[2]]
  p <- ggplot(data=end_info2)
  p + 
    geom_point(aes(x=merged_all_refs_tmp3.x, y=merged_all_refs_tmp3.y, col=merged_all_refs_tmp3.mean), na.rm=TRUE) 
    geom_point(data=as.data.frame(verticesdf), aes(x=V1, y=V2), color='red')  +
    theme_bw()
}

pdf("~/Desktop/all_rotation_plot")
plotz <- lapply(turning_outcome,plot_all_rotations)
dev.off()

turning_plus_angle <- lapply(seq_along(turning_outcome), function(x){
  turning_outcome[[x]][[1]]$angle <- rotation_turns[[x]]
  return(turning_outcome[[x]][[1]])})

turning_plus_angle <- do.call(rbind, turning_plus_angle)
p <- ggplot(data=as.data.frame(turning_plus_angle))
p + 
  geom_point(aes(x=merged_all_refs_tmp3.x, y=merged_all_refs_tmp3.y, col=merged_all_refs_tmp3.mean), na.rm=TRUE) +
  facet_wrap( ~ as.factor(angle)) +
#geom_point(data=as.data.frame(verticesdf), aes(x=V1, y=V2), color='red')  +
  theme_bw()

diff_0_180 <- filter(turning_plus_angle, angle==0) - filter(turning_plus_angle, angle==180)
diff_0_180_complete <- diff_0_180[complete.cases(diff_0_180),]
sum(diff_0_180_complete[,3])

add_standard_tr_ids <- function(end_info2){
  tri_values <- unique(end_info2[,3])
  tri_values <- tri_values[!is.nan(tri_values)]
  new_end_info <- NULL
  end_info2_id <- lapply(1:length(tri_values), function(x){
    subgroup <- filter(end_info2, merged_all_refs_tmp3.mean==tri_values[x])
    subgroup$id <- x
    return(subgroup)
  })
}

end_info2_plus <- do.call(rbind,add_standard_tr_ids(end_info2))

end_info_tri_ids <- merge(end_info2_plus, all_refs, by=1)

triangle_groups <- seq(from=1, to=dim(verticesdf)[1], by=3)
triangle_ids <- do.call(rbind, lapply(1:length(triangle_groups), function(x){
  wee_matrix <- as.data.frame(verticesdf[triangle_groups[x]:((triangle_groups[x])+2), ])
  wee_matrix$value <- x
  return(wee_matrix)
}))


total_ids <- unique(triangle_ids$value)
testytest <-do.call(rbind, lapply(seq_along(total_ids), function(x){
  one_triangle <- filter(triangle_ids, triangle_ids$value==x)
  bnd <- as.matrix(data.frame(one_triangle[,2], one_triangle[,1]))
  m <- as.matrix(data.frame(end_info2_plus[,2], end_info2_plus[,1]))
  present <- in.out(bnd, m)
  TRUTHS <- which(present)
  end_info2_plus[TRUTHS, 5] <- x
  return(end_info2_plus)
}))


triangle_centres <- lapply(seq_along(triangle_groups), function(x){
  wee_matrix <- as.data.frame(verticesdf[x:(x+2),])
  wee_matrix$centrex <- (sum(wee_matrix[,1]))/3
  wee_matrix$centrey <- (sum(wee_matrix[,2]))/3
  return(wee_matrix)
})
triangle_centres <- do.call(rbind, triangle_centres)



rotating_amodel_test <- function(amodel, rotate_value, fovea_coords, fovea_height){
  tst_matrix <- create_matrix(0, (dim(amodel)[2])-1, 0, (dim(amodel)[1])-1)
  #tst_value <- apply(tst_matrix, 1, function(linez){
  #  x_coord <- ((as.numeric(linez[[1]]))+1)
  #  y_coord <- ((as.numeric(linez[[2]]))+1)
  #  val <- amodel[y_coord, x_coord]
  #  return(val)
    
  #})
  #tst_matrix <- as.data.frame(tst_matrix)
  #tst_matrix$val <- unlist(tst_value)
  rownames(tst_matrix) <- paste(tst_matrix[,1], tst_matrix[,2], sep = ".")
  id_mapping <- data.frame(paste(tst_matrix[,1], tst_matrix[,2], sep="."), x=tst_matrix[,1], y=tst_matrix[,2], stringsAsFactors = FALSE)
  rownames(id_mapping) <- paste(tst_matrix[,1], tst_matrix[,2], sep=".")
  d1 <- dim(amodel)[2]*(1/3)
  scale_co <- d1/fovea_height
  fovea_coords_scale <- fovea_coords*scale_co
  tmp_mat_scale <- tst_matrix*scale_co
  centre_mat <- c((dim(amodel)[1]*scale_co)/2, (dim(amodel)[2]*scale_co)/2)
  foveal_difference <- c(centre_mat[1]-fovea_coords_scale[1], centre_mat[2]-fovea_coords_scale[2])
  scale_and_shift <- tmp_mat_scale
  scale_and_shift[,1] <- scale_and_shift[,1] + foveal_difference[2]
  scale_and_shift[,2] <- scale_and_shift[,2] + foveal_difference[1]
  #scale_shift_small <- scale_and_shift
  #scale_shift_small[,1] <- scale_shift_small[,1]/dim(amodel)[2]
  #scale_shift_small[,2] <- scale_shift_small[,2]/dim(amodel)[1]
  #scale_shift_small <- scale_shift_small[complete.cases(scale_shift_small),]
  tst_matrix2 <- scale_and_shift
  tst_matrix2[,1] <- tst_matrix2[,1]-(dim(amodel)[2])%/%2
  tst_matrix2[,2] <- tst_matrix2[,2]-(dim(amodel)[1])%/%2
  rotate_value_clock <- rotate_value
  rotate_value_clock_rad <- deg2rad(rotate_value_clock)
  rotation_matrix <- matrix(c(cos(rotate_value_clock_rad), sin(rotate_value_clock_rad), -sin(rotate_value_clock_rad), cos(rotate_value_clock_rad)), nrow=2, ncol=2)
  rotated_amodel <- tst_matrix %*% rotation_matrix
  rotated_amodel[,1] <- rotated_amodel[,1] + (dim(amodel)[2])%/%2
  rotated_amodel[,2] <- rotated_amodel[,2] + (dim(amodel)[1])%/%2
  rotate_and_id <- merge(rotated_amodel, id_mapping, by=0, all=TRUE, stringsAsFactors=FALSE)
  assign_vals <- apply(rotate_and_id, 1, function(linez){
    x_coord <- ((as.numeric(linez[[5]]))+1)
    y_coord <- ((as.numeric(linez[[6]]))+1)
    val <- amodel[y_coord, x_coord]
    return(val)
  
    })
  rotate_and_id$val<- unlist(assign_vals) 
  needed_info <- data.frame(rotate_and_id$V1, rotate_and_id$V2, rotate_and_id$val)
  needed_info[,1] <- needed_info[,1]/dim(amodel)[2]
  needed_info[,2] <- needed_info[,2]/dim(amodel)[1]
  }
p <- ggplot(data=needed_info)
p + 
  geom_point(aes(x=rotate_and_id.V1, y=rotate_and_id.V2, col=rotate_and_id.val), na.rm=TRUE) + 
  theme_bw()

p <- ggplot(data=tst_matrix)
p + 
  geom_point(aes(x=V1, y=V2, col=val), na.rm=TRUE) + 
  theme_bw()
#################################
############DEVELOP##############
#################################
center_align_tmp <- function (model, f = 0.6) 
  model$ycenter = fit_deriv(model$model$top_band[model$zcenter, 
                                                 ])$midpoint
  model$zcenter = fit_deriv(model$model$top_band[, model$ycenter])$midpoint
  centers = (dim(model$model$top_band)/2) - c(model$zcenter, 
                                              model$ycenter)
  rangs = list()
  for (x in 1:length(centers)) {
    if (centers[x] < 0) {
      r = abs(centers[x] * 2):dim(model$model$top_band)[x]
    }
    else {
      r = 1:(dim(model$model$top_band)[x] - centers[x] * 
               2)
    }
    rangs[[x]] = r
  ns = dim(model$model$top_band) * f
  centered_band = model$model$top_band[rangs[[1]], rangs[[2]]]
  centered_band = align_model_tmp(centered_band, ns)
  return(centered_band)
}
  
  
align_model_tmp <- function (top_band_surface, ns = c(75, 306)) {
  ndata = NULL
  ds = dim(top_band_surface)
  st = ns/2
  indexa <- ((ds[1]/2) - st[1])
  indexb <- ((ds[1]/2) + st[1])
  indexc <- ((ds[2]/2) - st[2])
  indexd <- ((ds[2]/2) + st[2])
  if (indexa < 0 | indexb < 0 | indexc < 0 | indexd < 0){
    ndata = NULL
  }
  else {
  ndata = top_band_surface[indexa:indexb, indexc:indexd]
  }
  return(ndata)
}

all_refs <- NULL
triangle_groups <- seq(from=1, to=(dim(verticesdf)[1]), by=3)
for (x in triangle_groups){
  small_matrix <- verticesdf[x:(x+2),]
  inoutvalues <- in.out(small_matrix, rotated_matrix$rotated_amodel)
  TRUTHS <- which(inoutvalues)
  
  refcoords <- NULL
  for (y in 1:length(TRUTHS)){
    tmp <- NULL
    tmp <- rbind(tmp,(rotated_matrix$rotated_amodel[TRUTHS[y],]))
    tmp <- cbind(tmp, as.integer((x+2)/3))
    refcoords <- rbind(refcoords, tmp)
  }
  rownames(refcoords) <- rownames(rotated_matrix$rotated_amodel)[TRUTHS]
  #refcoors <- cbind(rownames(m)[TRUTHS], refcoords)
  all_refs<-rbind(all_refs,refcoords)
  
}




odds_rotate_and_scale_table <- lapply(1:dim(all_params_odd_minus)[1], function(x){
  image_file <- all_params_odd_minus[[x,1]]
  fovea_height <- all_params_odd_minus[[x,5]]
  centre_f <- 0.8
  model <- run_generate_model_and_find_fovea_dip_poisition(image_file)
  e <- registar(image_file, debug=FALSE)
  amodel <- center_align(model, centre_f)
  smoothed_ellipse <- smooth_ellipse_max_slope(e, 31)
  rotate_val <- new_find_rotate(smoothed_ellipse)
  #fovea_height <- fovea_height/dim(amodel)[1]
  d1 <- dim(amodel)[1]*(1/3)
  fovea_height_scale <- fovea_height*dim(amodel)[1]
  scale_co <- d1/fovea_height_scale
  return(data.frame(image_file=image_file, rotate_val=rotate_val, scale_co=scale_co))
})
odds_rotate_and_scale_tablex <- do.call(rbind, odds_rotate_and_scale_table)
saveRDS(odds_rotate_and_scale_tablex, paste("/Users/currant/OCTeyes/odds_rotate_and_scale_dataframe",  strftime(Sys.time(), "%Y%m%d"), sep=""))


evens_rotate_and_scale_table <- lapply(1:dim(all_params_even_minus)[1], function(x){
  image_file <- all_params_even_minus[[x,1]]
  fovea_height <- all_params_even_minus[[x,5]]
  centre_f <- 0.8
  model <- run_generate_model_and_find_fovea_dip_poisition(image_file)
  e <- registar(image_file, debug=FALSE)
  amodel <- center_align(model, centre_f)
  smoothed_ellipse <- smooth_ellipse_max_slope(e,31)
  rotate_val <- new_find_rotate(smoothed_ellipse)
  #fovea_height <- fovea_height/dim(amodel)[1]
  d1 <- dim(amodel)[1]*(1/3)
  fovea_height_scale <- fovea_height*dim(amodel)[1]
  scale_co <- d1/fovea_height_scale
  return(data.frame(image_file=image_file, rotate_val=rotate_val, scale_co=scale_co))
})
evens_rotate_and_scale_tablex <- do.call(rbind, evens_rotate_and_scale_table)

############FOR ONE IMAGE STRIAGHT##########


model = run_generate_model_and_find_fovea_dip_poisition(image_file)
e <- registar(image_file, debug=FALSE)

amodel = center_align(model, center_f)


scale_centre_ids <- scale_and_centre(amodel, fovea_height, fovea_coords, n)

indicies <- make_tri_grid(amodel,5)
indicies[,2]<- round(as.numeric(indicies[,2]),3)
sqr_inds <- remove_to_square(indicies)
verticesdf <- get_triangle_vertices(sqr_inds,amodel)
all_refs <- get_triangle_ref_coords(scale_centre_ids$ref_table, verticesdf)
all_refs <- all_refs[complete.cases(all_refs),]
merged_all_refs <- merge(scale_centre_ids, all_refs, by=0, all=TRUE, stringsAsFactors=FALSE)
#merged_all_refs_tmp <- merge(tmp_mat_scale_shift_small, all_refs, by=0, all=TRUE, stringsAsFactors=FALSE)





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


############FOR ONE IMAGE ROTATED##########
model = run_generate_model_and_find_fovea_dip_poisition(image_file)
e <- registar(image_file, debug=FALSE)
amodel = center_align(model, center_f)

smoothed_ellipse <- smooth_ellipse_max_slope(e, 31)
rotate_value <- new_find_rotate(smoothed_ellipse)
rotated_matrix <- rotate_matrix(amodel, rotate_value)


indicies <- make_tri_grid(amodel,5)
indicies[,2]<- round(as.numeric(indicies[,2]),3)
sqr_inds <- remove_to_square(indicies)
verticesdf <- get_triangle_vertices(sqr_inds,amodel)
all_refs <- get_triangle_ref_coords(rotated_matrix$rotated_amodel, verticesdf)
all_refs <- all_refs[complete.cases(all_refs),]
merged_all_refs <- merge(rotated_matrix, all_refs, by=0, all=TRUE, stringsAsFactors=FALSE)
#merged_all_refs <- merge(rotated_matrix$rotated_amodel, all_refs, by=0, all=TRUE, stringsAsFactors=FALSE)
#merged_all_refs_tmp <- merge(tmp_mat_scale_shift_small, all_refs, by=0, all=TRUE, stringsAsFactors=FALSE)




test <- do.call(rbind, apply(merged_all_refs, 1, Ids2coordinates_rotate, original=amodel))
merged_all_refs_tmp2 <- merge(merged_all_refs, test, by=1, all=TRUE, stringsAsFactors=FALSE)

colnames(merged_all_refs_tmp2) <- c("coord_id", "x", "y","id_mapping", "x1", "y1", "x2", "x3", "tri_id", "value")


tri_ids <- unique(merged_all_refs_tmp2$tri_id)
all_tri_means <- do.call(rbind, lapply(tri_ids, tri_means))
colnames(all_tri_means) <- c("tri_id", "mean")

merged_all_refs_tmp3 <- merge(merged_all_refs_tmp2, all_tri_means, by="tri_id", all=TRUE, stringsAsFactors=FALSE)



#end_info <- do.call(rbind, apply(rotated_matrix$rotated_amodel, 1, get_final_table_rotate, merged_all_refs_tmp3=merged_all_refs_tmp3))
#end_info <- data.frame(merged_all_refs_tmp3$y1, merged_all_refs_tmp3$x1, merged_all_refs_tmp3$mean)
end_info2 <-  data.frame(merged_all_refs_tmp3$y, merged_all_refs_tmp3$x, merged_all_refs_tmp3$mean)


#triangle_plot_m <- matrix(NA,dim(amodel)[1],dim(amodel)[2])
#for (x in 1:dim(end_info)[1]){
#triangle_plot_m[end_info[x,1], end_info[x,2]] <- end_info[x,3]
#}
#image(triangle_plot_m)

p <- ggplot(data=end_info2)
p + 
  geom_point(aes(x=merged_all_refs_tmp3.x, y=merged_all_refs_tmp3.y, col=merged_all_refs_tmp3.mean), na.rm=TRUE) + 
  geom_point(data=as.data.frame(verticesdf), aes(x=V1, y=V2), color='red')  +
  theme_bw()

View(end_info2)
end_info2[,1] <- end_info2[,1]*(dim(amodel)[1])
end_info2[,2] <- end_info2[,2]*(dim(amodel)[2])
View(end_info2)
end_info2_ceiling <- ceiling(end_info2[,1:2])
end_info2_ceiling <- cbind(end_info2_ceiling, end_info2[,3])
end_info2_floor <- floor(end_info2[,1:2])
end_info2_floor <- cbind(end_info2_floor, end_info2[,3])
end_info2x <- rbind(end_info2_ceiling, end_info2_floor)
end_info2[,1:2] <- round(end_info2[,1:2],0)
View(end_info2)
endo_info2 <- end_info2x
end_info2 <- end_info2x

end_info2_subset <- subset(end_info2, end_info2[,1] >0 & end_info2[,1] < dim(amodel)[1])
end_info2_subset <- subset(end_info2_subset, end_info2_subset[,2] >0 & end_info2_subset[,2] < dim(amodel)[2])
View(end_info2)
View(end_info2_subset)
f3 <- rasterFromXYZ(end_info2_subset)
f_m3 <- as.matrix(f3)
image(f_m3)
points(verticesdf[,1], verticesdf[,2], col='blue')