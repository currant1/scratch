# Development script to look at scaling between two different sized grids. How to chop up so they fit the same locationa nd can all be tracked back to one script




######################################
##############FUNCTIONS###############
######################################

# Updated fit_centered_surface function that fixes for get_needed_params
fit_centered_surface <-function (centered_band, debug = FALSE) 
  {
    v360 = extract_circle_params(centered_band)
    fitted_derivs <- lapply(lapply(v360$s.ss, fit_deriv, 
           debug), unlist)
    length_check <- lapply(fitted_derivs, length)
    index <- which(length_check==0)
    fitted_derivs[index] <- NA
    divs = do.call(rbind, fitted_derivs)
    xy.s = convert_coords(v360, divs)
    max_slope = cbind(c(xy.s$f.xy1[, 2], xy.s$f.xy1[, 3]), c(xy.s$f.xy2[, 
                                                                      2], xy.s$f.xy2[, 3]))
    plateu = cbind(c(xy.s$f.xy1[, 4], xy.s$f.xy1[, 5]), c(xy.s$f.xy2[, 
                                                                   4], xy.s$f.xy2[, 5]))
    return(list(circle_params = v360, max_slope = max_slope, 
              plateu = plateu, divs = divs))
}

# Updated registar to deal with any NA values, making get_needed params to work for all
registar <-function (image_file, debug = FALSE) 
  {
    model = run_generate_model_and_find_fovea_dip_poisition(image_file)
    if (abs(0.5 - c(model$zcenter)/dim(model$model$top_band)[1]) > 
        0.15) {
      model$zcenter = dim(model$model$top_band)[1]/2
    }
    if (abs(0.5 - c(model$ycenter)/dim(model$model$top_band)[2]) > 
        0.15) {
      model$ycenter = dim(model$model$top_band)[2]/2
    }
    centered_band = center_and_square(model)
    fitted_center = fit_centered_surface(centered_band, debug)
    fitted_center$max_slope[is.na(fitted_center$max_slope[,1]),1]= mean(fitted_center$max_slope, na.rm=TRUE)
    fitted_center$max_slope[is.na(fitted_center$max_slope[,2]),2]= mean(fitted_center$max_slope, na.rm=TRUE)
    e1 = get.ellipse(fit.ellipse(fitted_center$max_slope[, 1], 
                                 fitted_center$max_slope[, 2]))
    fitted_center$plateu[is.na(fitted_center$plateu[,1]),1]= mean(fitted_center$plateu, na.rm=TRUE)
    fitted_center$plateu[is.na(fitted_center$plateu[,2]),2]= mean(fitted_center$plateu, na.rm=TRUE)
    e2 = get.ellipse(fit.ellipse(fitted_center$plateu[, 1], fitted_center$plateu[, 
                                                                                 2]))
    if (debug) {
      image(t(centered_band - median(centered_band)), col = topo.colors(10))
      lines(e1/dim(centered_band)[1], lwd = 2, col = "blue")
      lines(e2/dim(centered_band)[1], lwd = 2, col = "blue")
      matplot(fitted_center$max_slope[, 1]/dim(centered_band), 
              fitted_center$max_slope[, 2]/dim(centered_band), 
              col = "red", add = T, pch = "+", cex = 0.5)
      matplot(fitted_center$plateu[, 1]/dim(centered_band), 
              fitted_center$plateu[, 2]/dim(centered_band), col = "red", 
              add = T, pch = "+", cex = 0.5)
    }
    return(list(e1 = e1, e2 = e2, model = model$model, ycenter = model$ycenter, 
                zcenter = model$zcenter, centered_band = centered_band, 
                fitted_center = fitted_center))
  }


align_model<-function(top_band_surface, ns=c(75, 306)) {
    ndata = NULL
    ds = dim(top_band_surface)
    st = ns/2
    ndata=NULL
    if (dim(top_band_surface)[1]>((ds[1]/2)-st[1]) 
        &((ds[1]/2)-st[1])>0 
        & dim(top_band_surface)[1]>((ds[1]/2)+st[1])
        &((ds[1]/2)+st[1])>0 
        & dim(top_band_surface)[2]>((ds[2]/2)-st[2])
        &((ds[2]/2)-st[2])>0 
        & dim(top_band_surface)[2]>((ds[2]/2)+st[2])
        &((ds[2]/2)+st[2])>0){
      ndata = top_band_surface[((ds[1]/2)-st[1]): ((ds[1]/2)+st[1]), ((ds[2]/2)-st[2]): ((ds[2]/2)+st[2])]
    }
    return(ndata)
}


# version of get_triangle_vertices that works for below to examples and gets up to run 10 when running through all params

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






#################NUMBER1###############


image_file1 = all_eye_paths[[3]]
e1 <- registar(image_file1, debug=FALSE)
smooth_ellipse1 <- smooth_ellipse_max_slope(e1, 30)
plot(smooth_ellipse1)
rotate_value1 <- find_rotation(smooth_ellipse1)
rotate_value1 <- rotate_value1 - 90
rotate_value1 <- rotate_value1 + 60
new_radian1 <- NISTdegTOradian(rotate_value1)


ordered_ring1 <- e1$e2[order(e1$e2[,1]),]
fovea_centre_x1 <- e1$zcenter/(dim(e1$model$top_band)[1])
fovea_centre_y1 <- e1$ycenter/(dim(e1$model$top_band)[2])
fovea_coords1 <- c(fovea_coordsx1, fovea_coordsy1)
fovea_height1 <- abs((ordered_ring1[1,1]) - (ordered_ring1[dim(ordered_ring1)[1],1]))
fovea_height1 <- fovea_height1/dim(e1$centered_band)[1]
eye_name1 <- strsplit(as.character(image_file1), "/")
eye_name1 <- eye_name1[[length(eye_name1)]]
eye_name1 <- eye_name1[[length(eye_name1)]]


model1 <- run_generate_model_and_find_fovea_dip_poisition(image_file1)
centre_f=0.8
amodel1 <- center_align(model1, centre_f)

fovea_coordsx1 <- fovea_centre_x1
fovea_coordsy1 <- fovea_centre_y1
fovea_coords1 <- c(fovea_coordsx1, fovea_coordsy1)

sf=8
inds1 <- draw_tesselate_and_indices(fovea_coords1,fovea_height1,amodel1, sf)


inds_square1 <- remove_to_square(inds1)
xvals1 <- unlist(inds_square1[,1])
xvals1 <- round(xvals1, 8)
xvals_unique1 <- unique(xvals1)
xvals_unique_order1 <- sort(xvals_unique1)
triangles_across1 <- length(xvals_unique_order1)-2
yvals1 <- unlist(inds_square1[,2])
yvals1 <- round(yvals1, 8)
yvals_unique1 <- unique(yvals1)
yvals_unique_order1 <- sort(yvals_unique1)
triangles_down1 <- length(yvals_unique_order1)-1
verticesdf1 <- get_triangle_vertices(inds_square1, amodel1)
all_refs1 <- get_triangle_ref_coords(amodel1, verticesdf1)
all_triangle_means1 <- tesselate_mean(all_refs1, amodel1)
triangle_plot_m1 <- plot_tesselation(all_triangle_means1, amodel1, all_refs1)
how_many_triangles1  <- length(unique(all_refs1[,3]))
how_many_triangles_row1 <- c(eye_name1, how_many_triangles1, triangles_across1, triangles_down1)

############NUMBER2###############

image_file2 = all_eye_paths[[5]]
e2 <- registar(image_file2, debug=FALSE)
smooth_ellipse2 <- smooth_ellipse_max_slope(e2, 30)
plot(smooth_ellipse2)
rotate_value2 <- find_rotation(smooth_ellipse2)
rotate_value2 <- rotate_value2 - 90
rotate_value2 <- rotate_value2 + 60
new_radian2 <- NISTdegTOradian(rotate_value2)


ordered_ring2 <- e2$e2[order(e2$e2[,1]),]
fovea_centre_x2 <- e2$zcenter/(dim(e2$model$top_band)[1])
fovea_centre_y2 <- e2$ycenter/(dim(e2$model$top_band)[2])
fovea_coords2 <- c(fovea_coordsx2, fovea_coordsy2)
fovea_height2 <- abs((ordered_ring1[1,1]) - (ordered_ring1[dim(ordered_ring1)[1],1]))
fovea_height2 <- fovea_height2/dim(e2$centered_band)[1]
eye_name2 <- strsplit(as.character(image_file2), "/")
eye_name2 <- eye_name2[[length(eye_name2)]]
eye_name2 <- eye_name2[[length(eye_name2)]]


model2 <- run_generate_model_and_find_fovea_dip_poisition(image_file2)
centre_f=0.8
amodel2 <- center_align(model2, centre_f)

fovea_coordsx2 <- fovea_centre_x2
fovea_coordsy2 <- fovea_centre_y2
fovea_coords2 <- c(fovea_coordsx2, fovea_coordsy2)

sf=8
inds2 <- draw_tesselate_and_indices(fovea_coords2,fovea_height2,amodel2, sf)


inds_square2 <- remove_to_square(inds2)
xvals2 <- unlist(inds_square2[,1])
xvals2 <- round(xvals2, 8)
xvals_unique2 <- unique(xvals2)
xvals_unique_order2 <- sort(xvals_unique2)
triangles_across2 <- length(xvals_unique_order2)-2
yvals2 <- unlist(inds_square2[,2])
yvals2 <- round(yvals2, 8)
yvals_unique2 <- unique(yvals2)
yvals_unique_order2 <- sort(yvals_unique2)
triangles_down2 <- length(yvals_unique_order2)-1
verticesdf2 <- get_triangle_vertices(inds_square2, amodel2)
all_refs2 <- get_triangle_ref_coords(amodel2, verticesdf2)
all_triangle_means2 <- tesselate_mean(all_refs2, amodel2)
triangle_plot_m2 <- plot_tesselation(all_triangle_means2, amodel2, all_refs2)
how_many_triangles2 <- length(unique(all_refs2[,3]))
how_many_triangles_row2 <- c(eye_name2, how_many_triangles2, triangles_across2, triangles_down2)


###############TRYING TO SCALE TO ONE ANOTHER##################
seq(from=1, to=dim(verticesdf)[1], by=3)






mean_point <- NULL
triangleids <- unique(all_refs[,3])
for (x in triangleids){
  tmp <- all_refs[which(all_refs[,3]==x),]
  meanrow <- colMeans(tmp)
  mean_point <- rbind(mean_point, meanrow)
}
mean_point_round <- round(mean_point, 0)


uniquex <- unique(verticesdf[,1])
uniquey <- unique(verticesdf[,2])


tabley <- NULL
triangle_ids <- unique(all_refs[,3])
for (x in triangle_ids){
  onetriangle = all_refs[all_refs[,3]==x,]
  orderyonetri <- onetriangle[order(onetriangle[,2]),]
  minonetri <- orderyonetri[1,2]
  tmp <- NULL
  tmp <- cbind(x, minonetri)
  tabley <- rbind(tabley, tmp)
}

for (x in 1:dim(positional)[1])
jumps <- seq(from=1, to=dim(tabley)[1], by=as.numeric(positional[x,9])+as.numeric(positional[x,10])-1)
tabley_lines <- NULL
for (x in 1:(length(jumps)-1)){
  tmp <- tabley[jumps[x]:(jumps[x+1]-1),]
  tmp <- cbind(tmp, x)
  tabley_lines <- rbind(tabley_lines, tmp)
}

all_refs_lines <- NULL
for (x in 1:length(tabley_lines[,1])){
  one_triangle <- all_refs[which(all_refs[,3]==tabley[x,1]),]
  one_triangle_rowid <- cbind(one_triangle, tabley_lines[x,3])
  all_refs_lines <- rbind(all_refs_lines, one_triangle_rowid)
}

all_refs_new_id_cols <- NULL
for (x in unique(all_refs_lines[,4])){
  one_line <- all_refs_lines[which(all_refs_lines[,3]==x),]
  cols <- seq(from=1, to=dim(one_line)[1], by=1)
  tmp <- cbind(one_line, cols)
  all_refs_new_id_cols <- rbind(all_refs_new_id_cols, tmp)
}

all_refs_newids <- NULL
for (x in 1:dim(all_refs_new_id_cols)[1]){
  tmp <- all_refs_new_id_cols[x,]
  tmp <- cbind(tmp, paste(tmp[4], tmp[5], sep=":"))
  all_refs_newids <- rbind(all_refs_newids, tmp)
}

# creates new ids for each triangle and puts them onto old tables to make values more comparable across runs
uniform_ids <- function(all_refs, all_triangle_means){
  triangle_ids <- unique(all_refs[,3])
  new_ids <- create_new_triangle_ids(triangle_ids)
  reduced_triangles_lowest <- NULL
  for (x in triangle_ids){
    one_triangle <- all_refs[which(all_refs[,3]==x),]
    one_triangle_sort <- one_triangle[order(one_triangle[,1], one_triangle[,2]),]
    lowest_point <- one_triangle_sort[1,]
    reduced_triangles_lowest <- rbind(reduced_triangles_lowest, lowest_point)
  }
  reduced_triangles_lowest_sort <- reduced_triangles_lowest[order(reduced_triangles_lowest[,2], reduced_triangles_lowest[,1]),]
  reduced_triangles_lowest_sort <- cbind(reduced_triangles_lowest_sort, new_ids)
  
  all_triangle_means_id <- NULL
  for (x in 1:length(reduced_triangles_lowest_sort[,3])){
    tmp <- all_triangle_means[which(all_triangle_means[,2]==reduced_triangles_lowest_sort[[x,3]]),]
    tmp <- c(tmp, reduced_triangles_lowest_sort[(reduced_triangles_lowest_sort[[x,3]]),4])
    tmp <- unlist(tmp)
    tmp[1] <- as.numeric(tmp[1])
    all_triangle_means_id <- rbind(all_triangle_means_id, tmp)
  }
  
  
  all_refs_id <- NULL
  for (x in 1:length(reduced_triangles_lowest_sort[,3])){
    tmp <- all_refs[which(all_refs[,3]==reduced_triangles_lowest_sort[[x,3]]),]
    id <- rep((reduced_triangles_lowest_sort[(reduced_triangles_lowest_sort[[x,3]]),4]), dim(tmp)[1])
    tmp <- cbind(tmp, unlist(id))
    all_refs_id <- rbind(all_refs_id, tmp)
  }
  output <- list(all_triangle_means_id=all_triangle_means_id, all_refs_id=all_refs_id)
  return(output)
}


#creates a list of new alphabetic ids that can be used to cross reference
create_new_triangle_ids <- function(triangle_ids){
  ids <- list()
  for (x in LETTERS){
    for (y in LETTERS){
      for (z in LETTERS){
        id <- paste(x, y, z, sep="")
        ids <- c(ids, id)
      }
    }
  }
  ids_new <- ids[1:length(triangle_ids)] 
  return(ids_new)
}





