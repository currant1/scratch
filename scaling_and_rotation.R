################################
##########DEPENDENCIES##########
################################

library(foveafit)
library(colorspace)
library(NISTunits)

################################
###########FUNCTIONS############
################################

#Gets the full path of all eye dicom files within a given folder - currently only works for my organisation system obviously but very easily modified

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

#Gets the fovea centre and height of the fovea for use in scaling, using the registar output. Stores in table format.

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
    eye_data <- cbind(eye_data, eye_name, fovea_centre_x,fovea_centre_y, fovea_height, lorr, sex)
    all_eyes <- rbind( all_eyes, eye_data)
  }
  return(all_eyes)
}


#####################################
################DATA#################
#####################################


pathname <- "/Users/currant/OCTeyes/OCT_30_replicates_20161024/"



#####################################
#############ANALYSIS################
#####################################

all_eye_paths <- get_eye_paths(pathname)
all_eye_params10 <- get_needed_params(all_eye_paths[1:10])



######################################
###FOR TRIANGLES - Works with added 'next' loop control
######################################

find_rotation <- function(e){
  rotation_data <- NULL
  #rotation_data <- cbind(rotation_data, e$e2[1:180,])
  rotation_data <- cbind(rotation_data, e$fitted_center$plateu[1:180,])
  #rotation_data <- cbind(rotation_data, e$e2[181:360,])
  rotation_data <- cbind(rotation_data, e$fitted_center$plateu[181:360,])
  distance_value <- sqrt(((rotation_data[,1]-rotation_data[,3])^2) + (rotation_data[,2]-rotation_data[,4])^2)
  rotation_data <- cbind(rotation_data, distance_value)
  rotation_data_order <- rotation_data[order(rotation_data[,5]),]
  degree <- match(rotation_data_order[1,], rotation_data)
  degree <- degree[1]
  #degree <- which(rotation_data_order[1,]==rotation_data)
  degree <- degree - 90
  #degree <- degree + 60
  #new_radian <- NISTdegTOradian(degree)
  #return(new_radian)
  return(degree)
}

greatresult <- lapply(all_eye_paths[1], getFoveaSummary)
greatresult <- lapply("/Users/currant/OCTeyes/OCT_30_replicates_20161024/1145517/OCT/292422.dcm", getFoveaSummary)
getFoveaSummary <- function(image_file) {
  sf=8
  eye_name <- strsplit(image_file, "[/]")
  eye_name <- eye_name[[length(eye_name)]]
  batch <- eye_name[length(eye_name)-2]
  eye_name <- eye_name[length(eye_name)]
  model <- run_generate_model_and_find_fovea_dip_poisition(image_file)
  centre_f=0.8
  amodel <- center_align(model, centre_f)
  if (is.null(amodel)){
    return(NULL)
  }
  e <- registar(image_file)
  smoothed_ellipse <- smooth_ellipse_max_slope(e, 31)
  rotate_value <- new_find_rotate(smoothed_ellipse)
  fovea_coordsx <- all_eyes[[z,2]]
  fovea_coordsy <- all_eyes[[z,3]]
  fovea_coords <- c(fovea_coordsx, fovea_coordsy)
  fovea_height <- all_eyes[[z,4]]
  sf=8
  inds <- draw_tesselate_and_indices(fovea_coords,fovea_height,amodel, sf)
  inds_square <- remove_to_square(inds)
  verticesdf <- get_triangle_vertices(inds_square, amodel)
  all_refs <- get_triangle_ref_coords(amodel, verticesdf)
  all_triangle_means <- tesselate_mean(all_refs, amodel)
  triangle_plot_m <- plot_tesselation(all_triangle_means, amodel, all_refs)
  
  return(list(amodel=amodel, e= e, rotate_value=rotate_value, verticesdf=verticesdf, all_refs=all_refs, all_triangle_means=all_triangle_means, triangle_plot_m=triangle_plot_m))
}


rotate_matrix <- function(amodel, rotate_value){
  tst_matrix <- create_matrix(0, dim(amodel)[2], 0, dim(amodel)[1])
  rownames(tst_matrix) <- paste("x", tst_matrix[,1], "y", tst_matrix[,2], sep=".")
  tst_matrix[,1] <- tst_matrix[,1]-(dim(amodel)[2])%/%2
  tst_matrix[,2] <- tst_matrix[,2]-(dim(amodel)[1])%/%2
  rotate_value_clock <- rotate_value
  rotate_value_clock_rad <- deg2rad(rotate_value_clock)
  rotation_matrix_tst <- matrix(c(0,1,-1,0), nrow=2)
  rotation_matrix <- matrix(c(cos(rotate_value_clock_rad), sin(rotate_value_clock_rad), -sin(rotate_value_clock_rad), cos(rotate_value_clock_rad)), nrow=2, ncol=2)
  rotated_amodel <- tst_matrix %*% rotation_matrix
  return(rotated_amodel)
}
tst_matrix <- create_matrix(0, dim(amodel)[2], 0, dim(amodel)[1])
rownames(tst_matrix) <- paste("x", tst_matrix[,1], "y", tst_matrix[,2], sep=".")
tst_matrix[,1] <- tst_matrix[,1]-205
tst_matrix[,2] <- tst_matrix[,2]-52

rotate_value_clock <- rotate_value
rotate_value_clock_rad <- deg2rad(rotate_value_clock)
rotation_matrix_tst <- matrix(c(0,1,-1,0), nrow=2)
radninety <- deg2rad(45)
rotation_matrix <- matrix(c(cos(rotate_value_clock_rad), sin(rotate_value_clock_rad), -sin(rotate_value_clock_rad), cos(rotate_value_clock_rad)), nrow=2, ncol=2)
rotation_matrix <- matrix(c(cos(radninety), sin(radninety), -sin(radninety), cos(radninety)), nrow=2, ncol=2)
#rotated_amodel2 <- rotation_matrix %*% t(tst_matrix)
rotated_amodel <- tst_matrix %*% rotation_matrix

distancecalc <- function(rotated_amodel){
  
}
greatresult <- lapply(all_eye_paths[1:10], getFoveaSummary)
distancestry <- lapply(1:10, distancecalc)


distancevalstraight <- sqrt(((rotated_amodel[xval,1]- rotated_amodel[1,1])^2) + ((rotated_amodel[xval,2]-rotated_amodel[1,2])^2))
distancevalrotate <- sqrt(((rotated_amodel[dim(amodel)[1],1]- rotated_amodel[dim(amodel[1]),1])^2) + ((rotated_amodel[1,2]-rotated_amodel[1,2])^2))

for (z in 1:10){
#for (z in 1:dim(all_eyes)[1]){
  #image_file <- all_eye_paths[grep(all_eyes[z], all_eye_paths)]
  image_file <- all_eye_paths[[z]]
  
  eye_name <- strsplit(image_file, "[/]")
  eye_name <- eye_name[[length(eye_name)]]
  batch <- eye_name[length(eye_name)-2]
  eye_name <- eye_name[length(eye_name)]
  model <- run_generate_model_and_find_fovea_dip_poisition(image_file)
  centre_f=0.8
  amodel <- center_align(model, centre_f)
  if (is.null(amodel)){
    next
  }
  amodels <- append(amodels, amodel)
  fovea_coordsx <- all_eyes[[z,2]]
  fovea_coordsy <- all_eyes[[z,3]]
  fovea_coords <- c(fovea_coordsx, fovea_coordsy)
  fovea_height <- all_eyes[[z,4]]
  inds <- draw_tesselate_and_indices(fovea_coords,fovea_height,amodel, sf)
  inds_square <- remove_to_square(inds)
  points(fovea_coords[2], fovea_coords[1], col="blue", type="p", pch=16)
  xvals <- unlist(inds_square[,1])
  xvals <- round(xvals, 5)
  xvals_unique <- unique(xvals)
  xvals_unique_order <- sort(xvals_unique)
  xvals_above <- xvals_unique_order[xvals_unique_order > round(fovea_coords[2], 5) ]
  xvals_above_count <- length(xvals_above)
  xvals_below <- xvals_unique_order[xvals_unique_order < round(fovea_coords[2], 5) ]
  xvals_below_count <- length(xvals_below)
  triangles_across <- length(xvals_unique_order)-2
  yvals <- unlist(inds_square[,2])
  yvals <- round(yvals, 8)
  yvals_unique <- unique(yvals)
  yvals_unique_order <- sort(yvals_unique)
  triangles_down <- length(yvals_unique_order)-1
  verticesdf <- get_triangle_vertices(inds_square, amodel)
  all_refs <- get_triangle_ref_coords(amodel, verticesdf)
  all_triangle_means <- tesselate_mean(all_refs, amodel)
  triangle_plot_m <- plot_tesselation(all_triangle_means, amodel, all_refs)
  triangle_plots <- append(triangle_plots, triangle_plot_m)
  how_many_triangles <- length(unique(all_refs[,3]))
  how_many_triangles_row <- c(batch, eye_name, how_many_triangles, triangles_across, triangles_down)
  positional_row <- c(batch, eye_name, inds$horizontal_upper_count, inds$horizontal_lower_count, inds$posrange_upper_count, inds$posrange_lower_count, inds$negrange_upper_count, inds$negrange_lower_count, xvals_above_count, xvals_below_count)
  size_table_test <- rbind(size_table_test, how_many_triangles_row)
  #how_many_triangles <- how_many_triangleslength(how_many_triangles)
  positional2 <- rbind(positional2, positional_row)
 # row_new_table <- c(image_file, how_many_triangles)
  #size_table <- rbind(size_table, row_new_table)
}


size_table_test$total<- as.numeric(size_table_test[,4])*as.numeric(size_table_test[,5])
size_table_matrix <- matrix(unlist(size_table_test), nrow=10)
size_table_matrix[,6]<- as.numeric(size_table_matrix[,4])*as.numeric(size_table_matrix[,5])

smallest_x<- ((size_table_matrix)[order(size_table_matrix[,4])[1],])[4]
smallest_y<- ((size_table_matrix)[order(size_table_matrix[,4])[1],])[5]
size_table_numeric <- as.numeric(size_table_matrix[,3:6])
size_table_matrix<- cbind(size_table_matrix, as.numeric(size_table_matrix[,4]) - as.numeric(smallest_x))
size_table_matrix<- cbind(size_table_matrix, as.numeric(size_table_matrix[,5]) - as.numeric(smallest_y))

size_table_matrix
positional
belowx_order <- positional[order(positional[,10]),10] 
smallest_belowx <- belowx_order[1]
abovex_order <- positional[order(positional[,9]),9]
smallest_abovex <- abovex_order[1]
positional[,3:4] <- as.numeric(positional[,3:4])-1
belowy_order <- positional[order(positional[,4]),4]
smallest_belowy <- belowy_order[1]
abovey_order <- positional[order(positional[,3]),3]
smallest_belowx <- belowx_order[1]

#########SACLING TEST ON 1:10#################
triangle_amount_sort <- sort(size_table_matrix[,3])
smallest_number <- triangle_amount_sort[1]
smallest_grid_index <- which(size_table_matrix[,3]==smallest_number)
smallest_grid <- size_table_matrix[smallest_grid_index,]
size_table_matrix <- cbind(size_table_matrix, as.numeric(size_table_matrix[smallest_grid_index,3])-as.numeric(size_table_matrix[,3]))
size_table_matrix <- cbind(size_table_matrix, as.numeric(size_table_matrix[smallest_grid_index,4])-as.numeric(size_table_matrix[,4]))
size_table_matrix <- cbind(size_table_matrix, as.numeric(size_table_matrix[smallest_grid_index,5])-as.numeric(size_table_matrix[,5]))


###############TRYING TO SCALE TO ONE ANOTHER##################
test_subset <- size_table_matrix[5:10,]
test_subset[,3:6] <- as.numeric(test_subset[,3:6])
triangle_amount_sort <- sort(test_subset[,3])
smallest_number <- triangle_amount_sort[1]
smallest_grid_index <- which(test_subset[,3]==smallest_number)
smallest_grid <- test_subset[smallest_grid_index,]
test_subset <- cbind(test_subset, as.numeric(test_subset[smallest_grid_index,3])-as.numeric(test_subset[,3]))
test_subset <- cbind(test_subset, as.numeric(test_subset[smallest_grid_index,4])-as.numeric(test_subset[,4]))
test_subset <- cbind(test_subset, as.numeric(test_subset[smallest_grid_index,5])-as.numeric(test_subset[,5]))
for (x in 1){}






######################################
###FOR HEXAGONS
######################################
sf <- 8
amodels <- list()
hexagon_plots <- list()
size_table_test_hex <- NULL
#for (z in 1:dim(all_eyes)[1]){
for (z in 1:20){
  #image_file <- all_eye_paths[grep(all_eyes[z], all_eye_paths)]
  image_file <- all_eye_paths[[z]]
  eye_name <- strsplit(image_file, "[/]")
  eye_name <- eye_name[length(eye_name)]
  model <- run_generate_model_and_find_fovea_dip_poisition(image_file)
  centre_f=0.8
  amodel <- center_align(model, centre_f)
  amodels <- append(amodels, amodel)
  fovea_coordsx <- all_eyes[[z,2]]
  fovea_coordsy <- all_eyes[[z,3]]
  fovea_coords <- c(fovea_coordsx, fovea_coordsy)
  fovea_height <- all_eyes[[z,4]]
  inds <- draw_tesselate_and_indices(fovea_coords,fovea_height,amodel, sf)
  inds_square <- remove_to_square(inds)
  xvals <- unlist(inds_square[,1])
  xvals <- round(xvals, 8)
  xvals_unique <- unique(xvals)
  xvals_unique_order <- sort(xvals_unique)
  triangles_across <- length(xvals_unique_order)-2
  yvals <- unlist(inds_square[,2])
  yvals <- round(yvals, 8)
  yvals_unique <- unique(yvals)
  yvals_unique_order <- sort(yvals_unique)
  triangles_down <- length(yvals_unique_order)-1
  hexagon_indices <- get_hexagons(amodel, inds_square, fovea_coords)
  all_refs2 <- get_hexagon_ref_coords(amodel, hexagon_indices)
  all_hex_means <- tesselate_mean(all_refs2, amodel)
  hexagon_plot <- plot_tesselation_hex(all_hex_means, amodel,all_refs2)
  hexagon_plots <- append(hexagon_plots, triangle_plot_m)
  #how_many_triangles <- length(unique(all_refs[,3]))
  how_many_triangles_row <- c(eye_name, how_many_triangles)
  size_table_tes_hext <- rbind(size_table_test_hex, how_many_triangles_row)
  #how_many_triangles <- how_many_triangleslength(how_many_triangles)
  # row_new_table <- c(image_file, how_many_triangles)
  #size_table <- rbind(size_table, row_new_table)
}
names(amodel) <- all_eye_paths
names(triangle_plots) <- all_eye_paths



rotates <- NULL
for (x in all_eye_paths[1:20]){
  image_file <- x
  e <- registar(image_file, debug=FALSE)
  rotate_value <- find_rotation(e)
  rotates <- rbind(rotates, rotate_value)
}



for (x in all_eye_paths){
  image_file <- x
  model = run_generate_model_and_find_fovea_dip_poisition(image_file)
  center_f = 0.6
  amodel = center_align(model, center_f)
  fovea_coords <- c(0.5*dim(amodel)[1], 0.5*dim(amodel)[2])
  fovea_height <- which(all_eyes[])
}



sf <- 8
for (z in 1:10){
  #for (z in 1:dim(all_eyes)[1]){
  #image_file <- all_eye_paths[grep(all_eyes[z], all_eye_paths)]
  image_file <- all_eye_paths[[z]]
  eye_name <- strsplit(image_file, "[/]")
  eye_name <- eye_name[[length(eye_name)]]
  batch <- eye_name[length(eye_name)-2]
  eye_name <- eye_name[length(eye_name)]
  model <- run_generate_model_and_find_fovea_dip_poisition(image_file)
  centre_f=0.8
  amodel <- center_align(model, centre_f)
  if (is.null(amodel)){
    next
  }
  fovea_coordsx <- all_eyes[[z,2]]
  fovea_coordsy <- all_eyes[[z,3]]
  fovea_coords <- c(fovea_coordsx, fovea_coordsy)
  fovea_height <- all_eyes[[z,4]]
  inds <- draw_tesselate_and_indices(fovea_coords,fovea_height,amodel, sf)
  inds_square <- remove_to_square(inds)
  points(fovea_coords[2], fovea_coords[1], col="blue", type="p", pch=16)

}

x = "/Users/currant/OCTeyes/OCT_30_replicates_20161024/910147/OCT/516845.dcm"

