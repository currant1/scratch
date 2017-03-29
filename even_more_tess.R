#Tesselation script
#Takes in an image, scales and aligns it before performing a triangle tesselation and plotting means of each triangle.

###################################
############DEPENDECIES############
###################################

library(foveafit)


###################################
############FUNCTIONS##############
###################################

#Gives intersect of positive and negative equations
get_intersect <- function(slope1, intercept1, slope2, intercept2) { 
  x_cor = (intercept2-intercept1)/ifelse(slope2<0, slope1-slope2, slope1+slope2)
  y_cor = (slope1*abs(x_cor))+intercept1
  return(list("x"= x_cor, "y"= y_cor))
}

#plots image with layer of tesselation and returns the vertice coordinates. triangles are scaled based on sf which is the desired number of triangles across the height og the foveal dip and are centred over the fovea.
draw_tesselate_and_indices <- function(fovea_coords,fovea_height, amodel, sf){
  sf <- sf/2
  fovea_coords[1] <- fovea_coords[1]/dim(amodel)[1]
  fovea_coords[2] <- fovea_coords[2]/dim(amodel)[2]
  fovea_height_square <- fovea_height/dim(amodel)[1]
  step_value <- fovea_height_square/sf
  pos_fovea_intersect <- fovea_coords[1] - (1.732051*fovea_coords[2])
  neg_fovea_intersect <- fovea_coords[1] + (1.732051*fovea_coords[2])
  posrange <- seq(from=pos_fovea_intersect, to= 1, by=step_value)
  posrange <- append(posrange, seq(from = posrange[length(posrange)], to=-1.732, by=-step_value))
  posrange <- unique(posrange)
  posrange <- sort(posrange)
  negrange <- seq(from=neg_fovea_intersect, to=2.73, by=step_value)
  negrange <- append(negrange, seq(from=neg_fovea_intersect, to=0, by=-step_value))
  negrange <- unique(negrange)
  negrange <- sort(negrange)
  horizontal1 <- c(seq(from=fovea_coords[1], to=1, by=step_value/2))
  horizontal1 <- append(horizontal1, seq(from=fovea_coords[1], to=0, by=-step_value/2))
  horizontal1 <- unique(horizontal1)
  image(amodel)
  s = sapply(posrange, function(n) abline(a=n, b='1.732051'))
  s = sapply(negrange, function(n) abline(a=n, b='-1.732051'))
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
  return(indicies)
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
  inds_square <- inds[inds [,1]>0,]
  inds_square <- inds_square[inds_square[,1]<1,]
  inds_square <- inds_square[inds_square[,2]>0,]
  inds_square <- inds_square[inds_square[,2]<1,]
  plot(inds_square)
  return(inds_square)
}

#get three vertices of each triangle and scale to values of model
get_triangle_vertices <- function(inds_square, amodel){
  sorted_square <- inds_square[order( unlist(inds_square[,2]), unlist(inds_square[,1])),]
  sorted_square <- lapply(sorted_square, round, 14)
  sorted_square <- matrix(unlist(sorted_square), ncol=2)
  yunique <- unique(sorted_square[,2])
  twolines_table_order <- NULL
  for (x in 1:(length(yunique)-1)) {
    twolines <- NULL
    twolines <- rbind(twolines, sorted_square[which(sorted_square[,2]==yunique[[x]]),])
    twolines <- rbind(twolines, sorted_square[which(sorted_square[,2]==yunique[[x+1]]),])  
    twolines[,1] <- round(twolines[,1], 8)
    twolines <- unique(twolines)
    twolines_order <- twolines[order(unlist(twolines[,1])),]
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
  return_items <- list("verticesdf"=verticesdf, "yunique"=yunique)
  return(return_items )
}


#gets 6 indicesfor each hexagon
get_hexagons <- function(amodel, inds_square, fovea_coords){
  fovea_coords[1] <- fovea_coords[1]/dim(amodel)[1]
  fovea_coords[2] <- fovea_coords[2]/dim(amodel)[2] 
  sorted_square <- inds_square[order( unlist(inds_square[,2]), unlist(inds_square[,1])),]
  sorted_square <- lapply(sorted_square, round, 14)
  sorted_square <- matrix(unlist(sorted_square), ncol=2)
  row_values <- unique(sorted_square[,2])
  line_by_line <- NULL
  for (a in row_values){
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
  other_xval <- (which(abs(non_fovea_row - fovea_coords[2])==min(abs(non_fovea_row - fovea_coords[2]))))-1
  other_x_inds <- seq(from=other_xval, to=1, by=-3)
  other_x_inds <- append(other_x_inds, seq(from=other_xval, to=dim(non_fovea_row)[1], by=3))
  other_x_inds <- unique(other_x_inds)
  other_y_inds <- seq(from=(fovea_yvalue + 1), to=1, by=-2)
  other_y_inds <- append(other_y_inds, seq(from=(fovea_yvalue+1), to=dim(rows)[1], by=2))
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
    centre_index <- fovea_centre_coords[c,]
    upper_lines <- subset(line_by_line, line_by_line[,3]==centre_index[2])
    upper_lines <- rbind(upper_lines, subset(line_by_line, line_by_line[,3]==(centre_index[2]+1)))
    sorted_upper_lines <- upper_lines[order(upper_lines[,1]),]
    sorted_upper_lines <- unique(sorted_upper_lines)
    lower_lines <- subset(line_by_line, line_by_line[,3]==centre_index[2])
    lower_lines <- rbind(lower_lines, subset(line_by_line, line_by_line[,3]==(centre_index[2]-1)))
    sorted_lower_lines <- lower_lines[order(lower_lines[,1]),]
    sorted_lower_lines <- unique(sorted_lower_lines)
    centre_in_two_lines <- fovea_row[centre_index[1],]
    centre_in_two_lines <- centre_in_two_lines[1]
    two_centre_index <- which(sorted_upper_lines[,1]==centre_in_two_lines)
    hexagon <- NULL
    hexagon <- rbind(sorted_upper_lines[two_centre_index-2,], sorted_upper_lines[two_centre_index-1,], sorted_upper_lines[two_centre_index+1,],sorted_upper_lines[two_centre_index+2,], sorted_lower_lines[two_centre_index-1,], sorted_lower_lines[two_centre_index+1,])
    hexagon[,3] <- c
    hexagon_indices <- rbind(hexagon_indices, hexagon)
  }
  ##################
  #Half works up until here
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
    hexagon_indices <- rbind(hexagon_indices, hexagon)
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
  image(triangle_plot_m)
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


########################################
#################DATA###################
########################################



image_file = "/Users/currant/OCTeyes/OCT_30_replicates_20161024/1145517/OCT/292422.dcm"


########################################
###############ANALYSIS#################
########################################


model = run_generate_model_and_find_fovea_dip_poisition(image_file)


center_f = 0.6
amodel = center_align(model, center_f)
inds = tesselate(amodel, 10)

inds_square <- remove_to_square(inds)


verticesdf <- get_triangle_vertices(inds_square)

hexagon_indices <- get_hexagons(amodel, inds_square, fovea_coords)
all_refs2 <- get_triangle_ref_coords(amodel, hexagon_indices)
all_hex_means <- tesselate_mean(all_refs2, amodel)
hexagon_plot <- plot_tesselation(all_hex_means, amodel, all_refs2)

all_refs <- get_triangle_ref_coords(amodel, verticesdf)


all_triangle_means<- tesselate_mean(all_refs, amodel)

triangle_plot_m <- plot_tesselation(all_triangle_means,amodel)

tesselated <- complete_tesselation(amodel, 19)

fovea_coords<- c(38.5,154)
fovea_height <- 23.1

sf <-4


inds <- draw_tesselate_and_indices(fovea_coords, fovea_height, amodel, sf)
inds_square <- remove_to_square(inds)
hexagons <- get_hexagon_vertices(inds_square, fovea_coords,amodel)
image(amodel)
points(hexagons[,1:2])

pdf("/Users/currant/OCTeyes/getting_hexagons.pdf")

inds <- draw_tesselate_and_indices(fovea_coords, fovea_height, amodel, sf)
dev.off()


pathname <- "/Users/currant/OCTeyes/OCT_30_replicates_20161024/"
test_path <- "/Users/currant/OCTeyes/OCT_30_replicates_20161024/1145517/OCT/292422.dcm"
image_file <- test_path

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

amodels <- list()
triangle_plots <- list()
size_table <- NULL
for (z in 1:length(all_eye_paths)){
  image_file <- all_eye_paths[[z]]
  model <- run_generate_model_and_find_fovea_dip_poisition(image_file)
  center_f = 0.6
  amodel <- center_align(model, center_f)
  amodels <- append(amodels, amodel)
  inds <- draw_tesselate_and_indices(fovea_coords, fovea_height, amodel, sf)
  inds_square <- remove_to_square(inds)
  verticesdf <- get_triangle_vertices(inds_square, amodel)
  all_refs <- get_triangle_ref_coords(amodel, verticesdf$verticesdf)
  all_triangle_means<- tesselate_mean(all_refs, amodel)
  triangle_plot_m <- plot_tesselation(all_triangle_means,amodel)
  triangle_plots <- append(triangle_plots, triangle_plot_m)
  how_many_triangles <- verticesdf$yunique
  how_many_triangles <- length(how_many_triangles)
  row_new_table <- c(z, how_many_triangles)
  size_table <- rbind(size_table, row_new_table)
}
names(amodel) <- all_eye_paths
names(triangle_plots) <- all_eye_paths 
####################################

image(amodel, col=colorRampPalette(c('blue','green'))(100))
image(triangle_plot_m, add=TRUE, col=alpha(heat.colors(12),0.5))
image(triangle_plot_m)
points(inds_square)