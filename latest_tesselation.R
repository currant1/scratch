get_intersect <- function(slope1, intercept1, slope2, intercept2) { 
  x_cor = (intercept2-intercept1)/ifelse(slope2<0, slope1-slope2, slope1+slope2)
  y_cor = (slope1*abs(x_cor))+intercept1
  return(list("x"= x_cor, "y"= y_cor))
}

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
        matplot(coors$x, coors$y, col="red", pch=20, add=T)
      }
    }
  }
  return(indices)
}


image_file = "/Users/currant/OCTeyes/OCT_30_replicates_20161024/1145517/OCT/292422.dcm"
model = run_generate_model_and_find_fovea_dip_poisition(image_file)


center_f = 0.6
amodel = center_align(model, center_f)

# RUN
inds = tesselate(amodel, 10)
inds[,1] <- (unlist(inds[,1]))*77
inds[,2] <- (unlist(inds[,2]))*308
inds_square <- inds[inds [,1]>0,]
inds_square <- inds_square[inds_square[,1]<77,]
inds_square <- inds_square[inds_square[,2]>0,]
inds_square <- inds_square[inds_square[,2]<308,]
plot(inds_square)

sorted_square <- inds_square[order( unlist(inds_square[,2]), unlist(inds_square[,1])),]
yunique <- unique(sorted_square[,2])

twolines_table <- NULL
for (x in 1:(length(yunique)-1)) {
  twolines <- NULL
  twolines <- rbind(twolines, sorted_square[which(sorted_square[,2]==yunique[[x]]),])
  twolines <- rbind(twolines, sorted_square[which(sorted_square[,2]==yunique[[x+1]]),])  
  twolines_table <- cbind2(twolines_table, twolines)
  twolines_table_order <- twolines_table[order(unlist(twolines_table[,1])),]
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

vertices[,1] <- (unlist(vertices[,1]))*77
vertices[,2] <- (unlist(vertices[,2]))*308









#vertices_newdim <- subset(vertices, vertices[,2]<77)

get_matrix_indices <- function(m) {
  inds = NULL
  s1 = 1:ncol(m)
  for(x in 1:nrow(m)) {
    inds = rbind(inds, (s1))
  }
  inds = cbind(as.vector(inds), as.vector(t(inds)))
  return(inds)
}

m <- matrix(0,77,308)
ref_matrix <- get_matrix_indices(m)

verticesdf <- matrix(c(unlist(vertices)), nrow=dim(vertices), ncol=2)

x=1

all_refs <- NULL
for (x in seq(from = 1, to = dim(verticesdf)[1], by=3)){
  #inoutmatrix <- verticesdf[x:(x+2),]
  #inoutmatrix[,2] <- (inoutmatrix[,1])/4
  #inoutvalues <- in.out(inoutmatrix, ref_matrix)
  inoutvalues <- in.out(verticesdf[x:(x+2),], ref_matrix)
  TRUTHS <- which(inoutvalues)
  refcoords <- NULL
  #ref_ids = NULL 
  for (y in 1:length(TRUTHS)){
    tmp <- NULL
    tmp <- rbind(tmp,(ref_matrix[TRUTHS[y],]))
    tmp <- cbind(tmp, ((x+2)/3))
    #tmp <- cbind(tmp, rep(((x+2)/3), nrow(tmp)))
    refcoords <- rbind(refcoords, tmp)
    #refcoords <- cbind(refcoords, rep(((x+2)/3), nrow(refcoords)) )
    #ref_ids = c(ref_ids,  rep(((x+2)/3), nrow(refcoords)))
    #refcoords$ID <- rep(((x+2)/3), nrow(refcoords))
  }
  #refcoords$ID <- ref_ids
  all_refs<-rbind(all_refs,refcoords)
  
}


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



triangle_plot_m <- matrix(0,77,308)
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


##########################################################


library(mgcv)
data(columb.polys)
bnd <- columb.polys[[2]]
plot(bnd,type="n")
polygon(bnd)
x <- seq(7.9,8.7,length=20)
y <- seq(13.7,14.3,length=20)
gr <- as.matrix(expand.grid(x,y))
inside <- in.out(bnd,gr)
points(gr,col=as.numeric(inside)+1)



triangle1inout <- in.out(verticesdf[1:3,], ref_matrix)
TRUTHS <- which(triangle1inout)
refcoords <- NULL
for (x in 1:length(TRUTHS)){
  refcoords <- rbind(refcoords, ref_matrix[TRUTHS[x],])
}

triangle1values <-NULL
for (x in 1:dim(refcoords)[1]){
  triangle1values <- rbind(triangle1values, amodel[(refcoords[x,])[1],(refcoords[x,])[2]])
}