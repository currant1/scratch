library(foveafit)
library(changepoint)
library(bcp)
library(reshape2)

image_data <- dicomInfo("~/OCTeyes/OCT_16_normal_scans_20160905/656057.dcm")
image_data <- image_data$img
image_slice <- image_data[,,70]


get_mean_of_matrices <- function(x,mean_factor, image_data){
  mean_factor_n <- mean_factor%/%2
  m <- matrix(0, dim(image_slice)[1], dim(image_slice)[2])
  for (y in ((x-mean_factor_n):(x+mean_factor_n))){
    todays_matrix <- image_data[,,y]
    m = m + todays_matrix
  }
  m <- m/5
  return(m)
}

mean_of_eye <- lapply((3:((dim(image_data)[3])-2)), get_mean_of_matrices, 5, image_data)

a_mean_of_eye <- do.call(abind, c(mean_of_eye, along = 3))


m <- matrix(0, dim(image_slice)[1], dim(image_slice)[2])
for (y in 68:72){
  todays_matrix <- image_data[,,y]
  m= m + todays_matrix
}

m <- m/5

mm <- mean(m)
m_new <- m - mm


image_slice_mean <- image_slice_mean/5
barplot(image_slice[1,])
bcp_trial <- bcp(image_slice[1,])


get_bcp <- function(x, image_slice){
  line <- image_slice[x,]
  bcp_trial <- bcp(line)
  return(bcp_trial$posterior.prob)
}

slice_bcp <- lapply(1:length(image_slice[1,]), get_bcp, image_slice = image_slice)

find_max_bcp <- function(x, slice_bcp){
  bcp_slice <- slice_bcp[[x]]
  max_bcp_slice <- max(bcp_slice, na.rm=TRUE)
  max_bcp_slice_id <- which(bcp_slice == max_bcp_slice)
  return(max_bcp_slice_id)
}

max_bcps <- lapply(1:length(slice_bcp), find_max_bcp, slice_bcp = slice_bcp){
  
}



meanintensity_bcp <- function(image_slice_loc, borrow_f, image_slice){
  borrow_n <- borrow_f%/%2
  image_slice_cross_group <- image_slice[,(image_slice_loc-borrow_n):(image_slice_loc+borrow_n)]
  image_slice_cross_mean <- rowMeans(image_slice_cross_group)
  bcp_across_mean <- bcp(image_slice_cross_mean)
  max_slice <- which(bcp_across_mean$posterior.prob == max(bcp_across_mean$posterior.prob, na.rm=TRUE))
  return(max_slice)
}

meanintensity_bcp_test5 <- lapply(3:((dim(image_slice)[2])-2), meanintensity_bcp, 5, image_slice)
meanintensity_bcp_test5_unlist <- melt(meanintensity_bcp_test5)
plot(meanintensity_bcp_test5_unlist[,2], meanintensity_bcp_test5_unlist[,1])

meanintensity_bcp_test9 <- lapply(5:((dim(image_slice)[2])-4), meanintensity_bcp, 9, image_slice)
meanintensity_bcp_test9_unlist <- melt(meanintensity_bcp_test9)
plot(meanintensity_bcp_test9_unlist[,2], meanintensity_bcp_test9_unlist[,1])

meanintensity_bcp_test3 <- lapply(2:((dim(image_slice)[2])-1), meanintensity_bcp, 3, image_slice)
meanintensity_bcp_test3_unlist <- melt(meanintensity_bcp_test3)
plot(meanintensity_bcp_test3_unlist[,2], meanintensity_bcp_test3_unlist[,1])


intensity_med_bcp <- function(image_slice_loc, med_factor, image_slice){
  image_slice_slice <- image_slice[,image_slice_loc]
  image_slice_med <- runmed(image_slice_slice, med_factor)
  image_slice_med <- as.matrix(image_slice_med)
  image_med_bcp <- bcp(image_slice_med)
  image_bcp_max <- which(image_med_bcp$posterior.prob == max(image_med_bcp$posterior.prob, na.rm = TRUE))
  return(image_bcp_max)
}

intensity_med_bcp_test3 <- lapply(1:dim(image_slice)[2], intensity_med_bcp, med_factor = 3, image_slice = image_slice)
intensity_med_bcp_test3_unlist <- melt(intensity_med_bcp_test3)
plot(intensity_med_bcp_test3_unlist[,2], intensity_med_bcp_test3_unlist[,1])

intensity_med_bcp_test21 <- lapply(1:dim(image_slice)[2], intensity_med_bcp, med_factor = 21, image_slice = image_slice)
intensity_med_bcp_test21_unlist <- melt(intensity_med_bcp_test21)
plot(intensity_med_bcp_test21_unlist[,2], intensity_med_bcp_test21_unlist[,1])


bcp_med <- function(image_slice_loc, med_factor, image_slice){
  image_slice_slice <- image_slice[,image_slice_loc]
  image_slice_bcp <- bcp(image_slice_slice)
  image_bcp_med <- runmed(na.omit(image_slice_bcp$posterior.prob), med_factor)
  image_bcp_med_max <- which(image_bcp_med ==max(image_bcp_med))
  return(image_bcp_med_max)
}

bcp_med_test3 <- lapply(1:dim(image_slice)[2], bcp_med, med_factor=3, image_slice = image_slice)
bcp_med_test3_unlist <- melt(bcp_med_test3)
plot(bcp_med_test3_unlist[,2], bcp_med_test3_unlist[,1])

bcp_med_test5 <- lapply(1:dim(image_slice)[2], bcp_med, med_factor=5, image_slice = image_slice)
bcp_med_test5_unlist <- melt(bcp_med_test5)
plot(bcp_med_test5_unlist[,2], bcp_med_test5_unlist[,1])

bcp_med_test11 <- lapply(1:dim(image_slice)[2], bcp_med, med_factor = 11, image_slice = image_slice)
bcp_med_test11_unlist <- melt(bcp_med_test11)
plot(bcp_med_test11_unlist[,2], bcp_med_test11_unlist[,1])

bcp_med_test21 <- lapply(1:dim(image_slice)[2], bcp_med, med_factor = 21, image_slice = image_slice)
bcp_med_test21_unlist <- melt(bcp_med_test21)
plot(bcp_med_test21_unlist[,2], bcp_med_test21_unlist[,1])

sumbcp_bcp <- function(image_slice, borrow_f){
  bcp_list <- lapply(1:dim(image_slice)[2], function(x){
    image_slice_slice <- image_slice[,x]
    bcp_slice <- bcp(image_slice_slice)
    return(bcp_slice)
  })
  borrow_n <- borrow_f%/%2
  bcp_sum <- lapply((1+borrow_n):((dim(image_slice)[2])-borrow_n), function(x){
    bcp_group_list <- bcp_list[(x-borrow_n):(x+borrow_n)]
    bcp_prob_list <- lapply(1:length(bcp_group_list), function(x, bcp_group_list){
      bcp_group_item <- bcp_group_list[[x]]
      output <- bcp_group_item$posterior.prob
      return(output)
    }, bcp_group_list = bcp_group_list)
    bcp_prob_table <- do.call(rbind, bcp_prob_list)
    bcp_group_sum <- colSums(bcp_prob_table)
    return(bcp_group_sum)
  })
  max_bcp <- lapply(1:length(bcp_sum), function(x, bcp_sum){
    bcp_item <- bcp_sum[[x]]
    max_bcp_item <- which(bcp_item == (max(bcp_item, na.rm = TRUE)))
    return(max_bcp_item)
  }, bcp_sum = bcp_sum)
  return(max_bcp)
}

mean_eye_sumbcp <- lapply(1:dim(a_mean_of_eye)[3], function(x){
  image_slice <- image_data[,,x]
  output <- sumbcp_bcp(image_slice,5)
  return(output)
}) 
  

symbcp_bcp_test5 <- sumbcp_bcp(image_slice, 5)
symbcp_bcp_test_melt5 <- melt(symbcp_bcp_test5)
symbcp_bcp_test_melt5[,1] <- runmed(symbcp_bcp_test_melt5[,1], 21)
plot(symbcp_bcp_test_melt5[,2], symbcp_bcp_test_melt5[,1])

symbcp_bcp_test11 <- sumbcp_bcp(image_slice, 11)
symbcp_bcp_test_melt11 <- melt(symbcp_bcp_test11)
plot(symbcp_bcp_test_melt11[,2], symbcp_bcp_test_melt11[,1])

symbcp_bcp_test3 <- sumbcp_bcp(image_slice, 3)
symbcp_bcp_test_melt3 <- melt(symbcp_bcp_test3)
plot(symbcp_bcp_test_melt3[,2], symbcp_bcp_test_melt3[,1])

symbcp_bcp_testm5 <- sumbcp_bcp(m, 5)
symbcp_bcp_test_meltm5 <- melt(symbcp_bcp_testm5)
plot(symbcp_bcp_test_meltm5[,2], symbcp_bcp_test_meltm5[,1])
