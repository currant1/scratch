#Tesselation script
#Takes in an image, scales and aligns it before performing a triangle tesselation and plotting means of each triangle.

###################################
############DEPENDECIES############
###################################

library(foveafit)
library(colorspace)

library(NISTunits)


###################################
############FUNCTIONS##############
###################################

#Gives intersect of positive and negative equations
get_intersect <- function(slope1, intercept1, slope2, intercept2) { 
  x_cor = (intercept2-intercept1)/ifelse(slope2<0, slope1-slope2, slope1+slope2)
  y_cor = (slope1*abs(x_cor))+intercept1
  return(list("x"= x_cor, "y"= y_cor))
}

#Calculates degree of rotation in radians based on shortest length of fovea and making it the y axis
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

#plots image with layer of tesselation and returns the vertice coordinates. triangles are scaled based on sf which is the desired number of triangles across the height of the foveal dip and are centred over the fovea.
draw_tesselate_and_indices_wrotate <- function(fovea_coords,fovea_height, amodel, sf, rotate_value){
  sf <- sf/2
  #fovea_coords[1] <- fovea_coords[1]/dim(amodel)[1]
  #fovea_coords[2] <- fovea_coords[2]/dim(amodel)[2]
  #fovea_height_square <- fovea_height/dim(amodel)[1]
  new_radian_sixty <- rotate_value + 60
  new_radian_no_s <- rotate_value
  new_radian_minus <- 60 - rotate_value
  new_radian_sixty <- NISTdegTOradian(new_radian_sixty)
  new_radian_no_s <- NISTdegTOradian(new_radian_no_s)
  new_radian_minus <- NISTdegTOradian(new_radian_minus)
  
  #step_value <- fovea_height_square/sf
  step_value <- fovea_height/sf
  pos_fovea_intersect <- fovea_coords[2]-((sin(new_radian_sixty))/((sin(new_radian_sixty))/(tan(new_radian_minus)))*fovea_coords[1])
  neg_fovea_intersect <- fovea_coords[2]+(((sin(new_radian_sixty))/(cos(new_radian_sixty)))*fovea_coords[1])
  posrange <- seq(from=pos_fovea_intersect, to= 1, by=step_value)
  posrange <- append(posrange, seq(from = pos_fovea_intersect, to=-((sin(new_radian_sixty))/(((sin(new_radian_sixty))/tan(new_radian_minus)))), by=-step_value))
  posrange <- unique(posrange)
  posrange <- sort(posrange)
  negrange <- seq(from=neg_fovea_intersect, to=1+((sin(new_radian_sixty))/(cos(new_radian_sixty))), by=step_value*2)
  negrange <- append(negrange, seq(from=neg_fovea_intersect, to=0, by=-step_value*2))
  negrange <- unique(negrange)
  negrange <- sort(negrange)
  fovea_coord_intercept <- fovea_coords[2]+(((sin(new_radian_no_s))/(cos(new_radian_no_s)))*fovea_coords[1])
  horizontal1 <- c(seq(from=fovea_coord_intercept, to=1+((sin(new_radian_no_s))/(cos(new_radian_no_s))), by=((step_value/2))))
  horizontal1 <- append(horizontal1, seq(from=fovea_coord_intercept, to=0, by=-((step_value/2))))
  horizontal1 <- unique(horizontal1)
  image(amodel)
  s = sapply(posrange, function(n) abline(a=n, b=(sin(new_radian_sixty))/((sin(new_radian_sixty))/(tan(new_radian_minus)))))
  s = sapply(negrange, function(n) abline(a=n, b=-((sin(new_radian_sixty))/cos(new_radian_sixty))))
  s = sapply(horizontal1, function(n) abline(a=n, b=-((sin(new_radian_no_s))/(cos(new_radian_no_s)))))
  indicies=NULL
  for (l in 1:length(posrange)){
    for (j in 1:length(negrange)){
      coords <- get_intersect(((sin(new_radian_sixty))/((sin(new_radian_sixty))/(tan(new_radian_minus)))), posrange[l], -((sin(new_radian_sixty))/(cos(new_radian_sixty))), negrange[j])
      indicies <- rbind(indicies, coords)
      if(sum(coords>=0)==2) {
        matplot(coords$x, coords$y, col="blue", pch=20, add=T)
      }
    }
  }
    return(indicies)
}




tmp_mat <- create_matrix(0, dim(amodel)[2], 0, dim(amodel)[1])
rownames(tmp_mat) <- paste("x", tmp_mat[,1], "y", tmp_mat[,2], sep="")
n <- 20
d1<-dim(amodel)[2]*(1/3)
fov_height_scale <- fovea_height*dim(amodel)[2]
fovea_coords_scale <- c(fovea_coords[1]*dim(amodel)[1], fovea_coords[2]*dim(amodel)[2])
scale_co <- d1/fov_height_scale
fovea_coords_scale2 <- fovea_coords_scale*scale_co
tmp_mat_scale <- tmp_mat*scale_co
centre_mat <- c((dim(amodel)[1]*scale_co)/2, (dim(amodel)[2]*scale_co)/2)
dist_b_fovs <- c(centre_mat[1]-fovea_coords_scale2[1], centre_mat[2]-fovea_coords_scale2[2])
tmp_mat_scale_shift <- tmp_mat_scale
tmp_mat_scale_shift[,1] <- tmp_mat_scale_shift[,1]+dist_b_fovs[1]
tmp_mat_scale_shift[,2] <- tmp_mat_scale_shift[,2]+dist_b_fovs[2]
tmp_mat_scale_shift_small <- tmp_mat_scale_shift
tmp_mat_scale_shift_small[,1]<- tmp_mat_scale_shift_small[,1]/dim(amodel)[2]
tmp_mat_scale_shift_small[,2]<- tmp_mat_scale_shift_small[,2]/dim(amodel)[1]



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

stretch_and_cut <- function(fovea_coords, fovea_height, amodel, sf, n){
  tmp_mat <- create_matrix(0, dim(amodel)[2], 0, dim(amodel)[1])
  fovea_coords_scale <- c(fovea_coords[1]*dim(amodel)[2], fovea_coords[2]*dim(amodel)[1])
  fovea_coords_scale_round <- round(fovea_coords_scale, 0)
  d1 <- (1/3)
  scale_co <- d1/fovea_height
  tmp_mat_scale <- tmp_mat*scale_co
  return(tmp_mat_scale)
}

#plots image with layer of tesselation and returns the vertice coordinates. triangles are scaled based on sf which is the desired number of triangles across the height of the foveal dip and are centred over the fovea.
draw_tesselate_and_indices <- function(fovea_coords,fovea_height, amodel, sf, new_radian=1.0472){
  sf <- sf/2
  #fovea_coords[1] <- fovea_coords[[1]]/(dim(amodel)[1])
  #fovea_coords[2] <- fovea_coords[[2]]/(dim(amodel)[2])
  #fovea_height <- fovea_height/dim(amodel)[1]
  #fovea_height_square <- fovea_height/dim(amodel)[1]
  step_value <- fovea_height/sf
  pos_fovea_intersect <- fovea_coords[[1]] - (1.732051*fovea_coords[[2]])
  neg_fovea_intersect <- fovea_coords[[1]] + (1.732051*fovea_coords[[2]])
  posrange_upper <- seq(from=pos_fovea_intersect, to= 1, by=step_value)
  posrange_upper_count <- length(posrange_upper)
  posrange_lower <- seq(from = posrange_upper[length(posrange_upper)], to=-1.732051, by=-step_value)
  posrange_lower <- unique(posrange_lower)
  posrange_lower_count <- length(posrange_lower)
  posrange <- append(posrange_upper, posrange_lower)
  posrange <- unique(posrange)
  posrange <- sort(posrange)
  negrange_upper <- seq(from=neg_fovea_intersect, to=make_tri_grid <- , by=step_value)
  negrange_upper <- unique(negrange_upper)
  negrange_upper_count <- length(negrange_upper)
  negrange_to1 <- negrange_upper[negrange_upper<1]
  negrange_to1 <- unique(negrange_to1)
  negrange_upper_to1_count <- length(negrange_to1)
  negrange_lower <- seq(from=neg_fovea_intersect, to=0, by=-step_value)
  negrange_lower <- unique(negrange_lower)
  negrange_lower_count <- length(negrange_lower)
  negrange <- append(negrange_upper, negrange_lower)
  negrange <- unique(negrange)
  negrange <- sort(negrange)
  horizontal_upper <- c(seq(from=fovea_coords[[1]], to=1, by=step_value/2))
  horizontal_upper <- unique(horizontal_upper)
  horizontal_upper_count <- length(horizontal_upper)
  horizontal_lower <- c(seq(from=fovea_coords[[1]], to=0, by=-step_value/2))
  horizontal_lower <- unique(horizontal_lower)
  horizontal_lower_count <- length(horizontal_lower)
  horizontal1 <- append(horizontal_upper, horizontal_lower)
  horizontal1 <- unique(horizontal1)
  image(amodel)
  s = sapply(posrange, function(n) abline(a=n, b=1.732051))
  s = sapply(negrange, function(n) abline(a=n, b=-1.732051))
  s = sapply(horizontal1, function(n) abline(a=n, b=0))
  indicies=NULL
  for (l in 1:length(posrange)){
    for (j in 1:length(negrange)){
      coords <- get_intersect(1.732051, posrange[l], -1.732051, negrange[j])
      indicies <- rbind(indicies, coords)
      if(sum(coords>=0)==2) {
        matplot(coords$x, coords$y, col="blue", pch=20, add=T)
      }
    }
  }
  return(list(indicies=indicies, horizontal_lower_count=horizontal_lower_count,horizontal_upper_count=horizontal_upper_count, posrange_upper_count=posrange_upper_count, posrange_lower_count=posrange_lower_count, negrange_upper_count=negrange_upper_count, negrange_lower_count=negrange_lower_count))
}

#plots image with layer of tesselation and returns the vertice coordinates, by tom
tesselate <- function(amodel, f = 10) {
  posrange = seq(from = 1, to = -1.732, length=f)
  negrange = seq(from = 0, to = 2.732, length=f)
  horizontal1 = c(posrange[posrange>0]/2, 1-(posrange[posrange>0]/2))
  image(amodel)
  s = sapply(posrange, function(n) abline(a=n, b='1.732051'))
  s = sapply(negrange, function(n) abline(a=n, b='-1.732051'))
  s = sapply(horizontal1, function(n) abline(a=n, b=0))
  indices = NULL
  for(i in 1:f) {
    for(j in 1:f) {
      coors = get_intersect(1.732051, posrange[i], -1.732051, negrange[j])
      indices = rbind(indices, coors)
      if(sum(coors>=0)==2) {
        matplot(coors$x, coors$y, col="blue", pch=20, add=T)
      }
    }
  }
  return(indices)
}


# aligns model
align_model <- function(top_band_surface, ns=c(75, 306)) {
  ndata = NULL
  ds = dim(top_band_surface)
  st = ns/2
  ndata = top_band_surface[((ds[1]/2)-st[1]): ((ds[1]/2)+st[1]), ((ds[2]/2)-st[2]): ((ds[2]/2)+st[2])]
  return(ndata)
} 

#centres and aligns based on a percentage of the picture as decided by factor f
center_align<- function(model, f=0.6) {
  model$ycenter = fit_deriv(model$model$top_band[model$zcenter,])$midpoint
  model$zcenter = fit_deriv(model$model$top_band[,model$ycenter])$midpoint
  centers =(dim(model$model$top_band)/2)-c(model$zcenter, model$ycenter)
  rangs = list()
  for(x in 1:length(centers)) {
    if(centers[x]<0) {
      r = abs(centers[x]*2):dim(model$model$top_band)[x]
    } else {
      r = 1:(dim(model$model$top_band)[x]-centers[x]*2)
    }
    rangs[[x]] = r
  }
  ns=dim(model$model$top_band) * f
  centered_band = model$model$top_band[rangs[[1]], rangs[[2]]]
  centered_band = align_model(centered_band, ns)
  return(centered_band)
} 

#cut out any partial triangles
remove_to_square <- function(inds){
  inds <- inds$indicies
  inds_square <- inds[inds [,1]>0,]
  inds_square <- inds_square[inds_square[,1]<1,]
  inds_square <- inds_square[inds_square[,2]>0,]
  inds_square <- inds_square[inds_square[,2]<1,]
  plot(inds_square)
  return(inds_square)
}

#cut out any partial triangles when rotated matrix
remove_to_square_w_rotate <- function(inds){
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
  vertices[,1] <- (unlist(vertices[,1]))*dim(amodel)[1]
  vertices[,2] <- (unlist(vertices[,2]))*dim(amodel)[2]
  verticesdf <- matrix(c(unlist(vertices)), nrow=dim(vertices), ncol=2)
  return(verticesdf)
}


#gets 6 indicesfor each hexagon
get_hexagons <- function(amodel, inds_square, fovea_coords){
  #fovea_coords[1] <- fovea_coords[1]/dim(amodel)[1]
  #fovea_coords[2] <- fovea_coords[2]/dim(amodel)[2] 
  sorted_square <- inds_square[order( unlist(inds_square[,2]), unlist(inds_square[,1])),]
  sorted_square <- lapply(sorted_square, round, 14)
  sorted_square <- matrix(unlist(sorted_square), ncol=2)
  row_values <- unique(sorted_square[,2])
  #row_values <- round(row_values, 8)
  line_by_line <- NULL
  for (a in row_values){
    a <- round(a, 6)
    line <- cbind(sorted_square[which(sorted_square[,2]==a)], a,match(a, row_values))
    line_by_line <- rbind(line_by_line, line)
  }
  fovea_yvalue <- (line_by_line[,3][which(fovea_coords[1]==line_by_line[,2])])[1]
  all_yvals <- unique(line_by_line[,3])
  fovea_row <- subset(line_by_line, line_by_line[,3]==fovea_yvalue)
  fovea_row <- (unique(fovea_row))
  fovea_row <- fovea_row[order(fovea_row[,1]),]
  non_fovea_row <- subset(line_by_line, line_by_line[,3]==(fovea_yvalue+1))
  non_fovea_row <- unique(non_fovea_row)
  non_fovea_row <- non_fovea_row[order(non_fovea_row[,1]),]
  fovea_xvalue <- which(fovea_row[,1]==fovea_coords[2])
  fovea_x_inds <- seq(from=fovea_xvalue, to=1, by=-3)
  fovea_x_inds <- append(fovea_x_inds, seq(from=fovea_xvalue, to=dim(fovea_row)[1], by=3))
  fovea_x_inds <- unique(fovea_x_inds)
  rows <- unique(line_by_line[,2:3])
  fovea_y_inds <- seq(from=fovea_yvalue, to =1, by=-2)
  fovea_y_inds <- append(fovea_y_inds, seq(from=fovea_yvalue, to=dim(rows)[1], by=2))
  fovea_y_inds <- unique(fovea_y_inds)
  other_yvalue <- fovea_yvalue +1
  if (length(fovea_row)%%2 != 0){
    other_xval <- fovea_xvalue -1
  } else {
    other_xval <- fovea_xvalue -2
  }
  #other_xval <- (which(abs(non_fovea_row - fovea_coords[2])==min(abs(non_fovea_row - fovea_coords[2])))) -1
  #other_xval <- (which(abs(non_fovea_row - fovea_coords[2])==min(abs(non_fovea_row - fovea_coords[2]))))-2
  other_x_inds <- seq(from=other_xval, to=1, by=-3)
  other_x_inds <- append(other_x_inds, seq(from=other_xval, to=dim(non_fovea_row)[1], by=3))
  other_x_inds <- unique(other_x_inds)
  #other_y_inds <- seq(from=(fovea_yvalue + 1), to=1, by=-2)
  other_y_inds <- seq(from=other_yvalue, to=1, by=-2)
  #other_y_inds <- append(other_y_inds, seq(from=(fovea_yvalue+1), to=dim(rows)[1], by=2))
  other_y_inds <- append(other_y_inds, seq(from=other_yvalue, to=dim(rows)[1], by=2))
  other_y_inds <- unique(other_y_inds)
  tmp1 <-lapply(fovea_y_inds, function(y) cbind(fovea_x_inds, rep(y,length(fovea_x_inds))))
  fovea_centre_coords <- do.call(rbind,tmp1)
  fovea_centre_coords <- unlist(fovea_centre_coords)
  fovea_centre_coords <- subset(fovea_centre_coords, fovea_centre_coords[,1] <= (dim(fovea_row)[1]-1))
  fovea_centre_coords <- subset(fovea_centre_coords, fovea_centre_coords[,2] <= (length(all_yvals)-1))
  fovea_centre_coords <- subset(fovea_centre_coords, fovea_centre_coords[,1] >= 2)
  fovea_centre_coords <- subset(fovea_centre_coords, fovea_centre_coords[,2] >= 2)
  tmp2 <-lapply(other_y_inds, function(y) cbind(other_x_inds, rep(y,length(other_x_inds))))
  other_centre_coords <- do.call(rbind,tmp2)
  other_centre_coords <- subset(other_centre_coords, other_centre_coords[,1] <= (dim(non_fovea_row)[1]-1))
  other_centre_coords <- subset(other_centre_coords, other_centre_coords[,2] <= (length(all_yvals)-1))
  other_centre_coords <- subset(other_centre_coords, other_centre_coords[,1] >= 2)
  other_centre_coords <- subset(other_centre_coords, other_centre_coords[,2] >= 2)
  hexagon_indices <- NULL
  for( c in 1:dim(fovea_centre_coords)[1]){
    print(c)
    centre_index <- fovea_centre_coords[c,]
    upper_lines <- subset(line_by_line, line_by_line[,3]==centre_index[2])
    upper_lines <- rbind(upper_lines, subset(line_by_line, line_by_line[,3]==(centre_index[2]+1)))
    sorted_upper_lines <- upper_lines[order(upper_lines[,1]),]
    sorted_upper_lines <- unique(sorted_upper_lines[,1:2])
    
    lower_lines <- subset(line_by_line, line_by_line[,3]==centre_index[2])
    lower_lines <- rbind(lower_lines, subset(line_by_line, line_by_line[,3]==(centre_index[2]-1)))
    sorted_lower_lines <- lower_lines[order(lower_lines[,1]),]
    sorted_lower_lines <- unique(sorted_lower_lines[,1:2])
    centre_in_two_lines <- fovea_row[centre_index[1],]
    centre_in_two_lines <- centre_in_two_lines[1]
    two_centre_index <- which(sorted_upper_lines[,1]==centre_in_two_lines)
    hexagon <- NULL
    hexagon <- rbind(sorted_upper_lines[two_centre_index-2,], sorted_upper_lines[two_centre_index-1,], sorted_upper_lines[two_centre_index+1,],sorted_upper_lines[two_centre_index+2,], sorted_lower_lines[two_centre_index-1,], sorted_lower_lines[two_centre_index+1,])
    hexagon[,3] <- c
    hexagon <- hexagon[order(hexagon[,2], hexagon[,1]),]
    hexagon_ordered <- NULL
    hexagon_ordered <- rbind(hexagon_ordered, hexagon[2,], hexagon[4,], hexagon[6,], hexagon[5,], hexagon[3,], hexagon[1,])
    hexagon_indices <- rbind(hexagon_indices, hexagon_ordered)
  }
  for( d in 1:dim(other_centre_coords)[1]){
    centre_index <- other_centre_coords[d,]
    upper_lines <- subset(line_by_line, line_by_line[,3]==centre_index[2])
    upper_lines <- rbind(upper_lines, subset(line_by_line, line_by_line[,3]==(centre_index[2]+1)))
    sorted_upper_lines <- upper_lines[order(upper_lines[,1]),]
    sorted_upper_lines <- unique(sorted_upper_lines)
    lower_lines <- subset(line_by_line, line_by_line[,3]==centre_index[2])
    lower_lines <- rbind(lower_lines, subset(line_by_line, line_by_line[,3]==(centre_index[2]-1)))
    sorted_lower_lines <- lower_lines[order(lower_lines[,1]),]
    sorted_lower_lines <- unique(sorted_lower_lines)
    centre_in_two_lines <- non_fovea_row[centre_index[1],]
    centre_in_two_lines <- centre_in_two_lines[1]
    two_centre_index <- which(sorted_upper_lines[,1]==centre_in_two_lines)
    hexagon <- NULL
    hexagon <- rbind(sorted_upper_lines[two_centre_index-2,], sorted_upper_lines[two_centre_index-1,], sorted_upper_lines[two_centre_index+1,],sorted_upper_lines[two_centre_index+2,], sorted_lower_lines[two_centre_index-1,], sorted_lower_lines[two_centre_index+1,])
    hexagon[,3] <- d + dim(fovea_centre_coords)[1]
    hexagon <- hexagon[order(hexagon[,2], hexagon[,1]),]
    hexagon_ordered <- NULL
    hexagon_ordered <- rbind(hexagon_ordered, hexagon[2,], hexagon[4,], hexagon[6,], hexagon[5,], hexagon[3,], hexagon[1,])
    hexagon_indices <- rbind(hexagon_indices, hexagon_ordered)
  }
  hexagon_indices[,1] <- (hexagon_indices[,1])*(dim(amodel)[1])
  hexagon_indices[,2] <- (hexagon_indices[,2])*(dim(amodel)[2])
  return(hexagon_indices)
}





# Creates a 2 column matrix used for reference
create_matrix <- function(xmin, xmax, ymin, ymax){
  tmp <-lapply(ymin:ymax, function(y) cbind(xmin:xmax, rep(y,length(xmin:xmax))))
  final <- do.call(rbind,tmp)
  return(final)
}

#gets the coordinates on the amodel for all points within each triangle
get_triangle_ref_coords <- function(amodel, verticesdf){
  m <- create_matrix(0,dim(amodel)[1], 0, dim(amodel)[2])
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
    all_refs<-rbind(all_refs,refcoords)
    
  }
  return(all_refs)
}

#gets the coordinates on the amodel for all points within each hexagon
get_hexagon_ref_coords <- function(amodel, hexagon_indices){
  m <- create_matrix(0, dim(amodel)[1], 0, dim(amodel)[2])
  all_refs2 <- NULL
  hexagon_groups <- seq(from=1, to=(dim(hexagon_indices)[1]), by=6)
  for (x in hexagon_groups){
    small_matrix <- hexagon_indices[x:(x+5),]
    inoutvalues <- in.out(small_matrix, m)
    TRUTHS <- which(inoutvalues)
    
    refcoords <- NULL
    for (y in 1:length(TRUTHS)){
      tmp <- NULL
      tmp <- rbind(tmp, (m[TRUTHS[y],]))
      tmp <- cbind(tmp, as.integer((x+2)/3))
      refcoords <- rbind(refcoords,tmp)
    }
    all_refs2 <- rbind(all_refs2, refcoords)
  }
  return(all_refs2)
}

#apply a mean across each triangle
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

#plot calculated mean to all points within triangle with that value and plot it back to a triangluar matrix
plot_tesselation <- function(all_triangle_means, amodel, all_refs){
  triangle_plot_m <- matrix(NA,dim(amodel)[1],dim(amodel)[2])
  for (x in 1:dim(all_triangle_means)[1]){
    triangle_value <- all_triangle_means[x,1]
    triangle_locations <- subset(all_refs, all_refs[,3]==x )
    for (y in 1:dim(triangle_locations)[1]){
      x_coord <- (triangle_locations[y,])[1]
      y_coord <- (triangle_locations[y,])[2]
      triangle_plot_m[x_coord, y_coord] <- triangle_value
      #triangle_plot_m[(triangle_locations[y,])[1], (triangle_locations[y,])[2]]<- triangle_value
    }
  }
  image(triangle_plot_m, col=heat_hcl(12))
  return(triangle_plot_m)
}


#plot calculated mean to all points within triangle with that value and plot it back to a triangluar matrix
plot_tesselation_hex <- function(all_triangle_means, amodel, all_refs){
  all_triangle_means <- all_hex_means
  all_refs <- all_refs2
  triangle_plot_m <- matrix(NA,dim(amodel)[1],dim(amodel)[2])
  for (x in 1:dim(all_triangle_means)[1]){
    triangle_value <- all_triangle_means[x,1]
    triangle_locations <- subset(all_refs, all_refs[,3]==all_triangle_means[x,2] )
    for (y in 1:dim(triangle_locations)[1]){
      x_coord <- (triangle_locations[y,])[1]
      y_coord <- (triangle_locations[y,])[2]
      triangle_plot_m[x_coord, y_coord] <- triangle_value
      #triangle_plot_m[(triangle_locations[y,])[1], (triangle_locations[y,])[2]]<- triangle_value
    }
  }
  image(triangle_plot_m, col=heat_hcl(12))
  return(triangle_plot_m)
}


#completes all functions and returns an image and the corresponding matrix for said image
complete_tesselation <- function(amodel,tf){
  inds <- tesselate(amodel,tf)
  inds_square <- remove_to_square(inds)
  verticesdf <- get_triangle_vertices(inds_square)
  all_refs<- get_triangle_ref_coords(amodel,verticesdf)
  all_triangle_means<- tesselate_mean(all_refs,amodel)
  triangle_plot_m<- plot_tesselation(all_triangle_means,amodel, all_refs)
  return(triangle_plot_m)
}


# gets all individual dicom paths from input directory and retrurns in list
get_eye_paths <- function(pathname){
  filenames <- list.files(path=pathname)
  all_eye_paths <- list()
  for (x in filenames){
    pathname1 <- paste(pathname, x, sep="")
    eye_files <- list.files(pathname1, pattern=".fda")
    for(y in eye_files){
      full_eye_path <- paste(pathname1, y, sep="")
      all_eye_paths <- append(all_eye_paths, full_eye_path)
    }
  } 
}

#Final fixed version of smooth ellipse function having fixed problem and lost location of original function
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


########################################
#################DATA###################
########################################



image_file = "/Users/currant/OCTeyes/OCT_30_replicates_20161024/1145517/OCT/292422.dcm"
image_file = "/Users/currant/OCTeyes/OCT_16_normal_scans_20160905/656057.dcm"    

########################################
###############ANALYSIS#################
########################################


model = run_generate_model_and_find_fovea_dip_poisition(image_file)
e <- registar(image_file, debug=FALSE)

enew_radian <- find_rotation(e)

center_f = 0.6
amodel = center_align(model, center_f)
inds = tesselate(amodel, 10)

inds_square <- remove_to_square(inds)


verticesdf <- get_triangle_vertices(inds_square, amodel)

hexagon_indices <- get_hexagons(amodel, inds_square, fovea_coords)
plot(hexagon_indices[,1:2])
all_refs2 <- get_hexagon_ref_coords(amodel, hexagon_indices)
all_hex_means <- tesselate_mean(all_refs2, amodel)
hexagon_plot <- plot_tesselation(all_hex_means, amodel,all_refs2)

all_refs <- get_triangle_ref_coords(amodel, verticesdf)


all_triangle_means<- tesselate_mean(all_refs, amodel)

triangle_plot_m <- plot_tesselation(all_triangle_means,amodel, all_refs)

tesselated <- complete_tesselation(amodel, 19)

fovea_coords<- c(38.5,154)
fovea_height <- 23.1

sf <- 8


inds <- draw_tesselate_and_indices(fovea_coords, fovea_height, amodel, sf, new_radian)
inds_square <- remove_to_square(inds)
hexagons <- get_hexagon_vertices(inds_square, fovea_coords,amodel, al_refs2)
image(amodel)
points(hexagons[,1:2])

pdf("/Users/currant/OCTeyes/getting_hexagons.pdf")

inds <- draw_tesselate_and_indices(fovea_coords, fovea_height, amodel, sf)
dev.off()


pathname <- "/Users/currant/OCTeyes/OCT_30_replicates_20161024/"

filenames <- list.files(path=pathname)
all_eye_paths <- list()
for (x in filenames){
  pathname1 <- paste(pathname, x, sep="")
  eye_files <- list.files(pathname1, pattern=".fda")
  for(y in eye_files){
    full_eye_path <- paste(pathname1, y, sep="")
    all_eye_paths <- append(all_eye_paths, full_eye_path)
  }
}

for (z in all_eye_paths){
  image_file <- z
  model <- run_generate_model_and_find_fovea_dip_poisition(image_file)
  center_f = 0.6
  amodel <- center_align(model, center_f)
  inds <- draw_tesselate_and_indices(fovea_coords, fovea_height, amodel, sf)
  inds_square <- remove_to_square(inds)
  verticesdf <- get_triangle_vertices(inds_square)
  all_refs <- get_triangle_ref_coords(amodel, verticesdf)
  all_triangle_means<- tesselate_mean(all_refs, amodel)
  triangle_plot_m <- plot_tesselation(all_triangle_means,amodel)
}




####################################

image(amodel, col=colorRampPalette(c('blue','green'))(100))
image(triangle_plot_m, add=TRUE, col=alpha(heat.colors(12),0.5))
image(triangle_plot_m)
points(inds_square)
