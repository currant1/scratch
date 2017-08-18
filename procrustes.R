######################################################
##################TRYING PROCRUSTES###################
######################################################

#######################################
#############DEPENDENCIES##############
#######################################

library(vegan)



#######################################
###############FUNCTIONS###############
#######################################



rotate_and_plot_p <- function(image_file, center_f, triangle_num=5, scale_f, fovea_coords,rotate_val, n){
  model <- run_generate_model_and_find_fovea_dip_poisition(image_file)
  e <- registar(image_file, debug=FALSE)
  amodel <- center_align(model, center_f)
  #################
  #image_file <- all_params_odd[[x,1]]
  #fovea_height <- all_params_odd[[x,5]]
  #fovea_coordx <- all_params_odd[[x,3]]
  fovea_coordx_scale <- fovea_coords[2]*(dim(amodel)[2])
  #fovea_coordy <- all_params_odd[[x,4]]
  fovea_coordy_scale <- fovea_coords[1]*(dim(amodel)[1])
  fovea_coords <- c(fovea_coordy_scale, fovea_coordx_scale)
  #################
  smoothed_ellipse <- smooth_ellipse_max_slope(e, 31)
  rotate_value <- rotate_val
  scale_centre_ids <- scale_and_centre_rotate_p(amodel, fovea_coords, n, scale_f)
  #rotated_matrix <- rotate_matrix(amodel, rotate_value)
  #rotate_value <- 90
  rotated_matrix <- rotate_matrix_p(scale_centre_ids$ref_table, rotate_val, amodel)
  
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

just_scale_and_plot_p <- function(image_file, center_f, triangle_num=5, scale_f, fovea_coords, n){
  model <- run_generate_model_and_find_fovea_dip_poisition(image_file)
  e <- registar(image_file, debug=FALSE)
  amodel <- center_align(model, center_f)
  #################
  #image_file <- all_params_odd[[x,1]]
  #fovea_height <- all_params_odd[[x,5]]
  #fovea_coordx <- all_params_odd[[x,3]]
  fovea_coordx_scale <- fovea_coords[2]*(dim(amodel)[2])
  #fovea_coordy <- all_params_odd[[x,4]]
  fovea_coordy_scale <- fovea_coords[1]*(dim(amodel)[1])
  fovea_coords <- c(fovea_coordy_scale, fovea_coordx_scale)
  #################
  smoothed_ellipse <- smooth_ellipse_max_slope(e, 31)
  #rotate_value <- new_find_rotate(smoothed_ellipse)
  scale_centre_ids <- scale_and_centre_rotate_p(amodel,fovea_coords, n, scale_f)
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

rotate_and_plot_no_scale_p <-function(image_file, center_f, triangle_num=5, fovea_coords, rotate_val, n){
  model <- run_generate_model_and_find_fovea_dip_poisition(image_file)
  e <- registar(image_file, debug=FALSE)
  amodel <- center_align(model, center_f)
  #################
  #image_file <- all_params_odd[[x,1]]
  #fovea_height <- all_params_odd[[x,5]]
  #fovea_coordx <- all_params_odd[[x,3]]
  fovea_coordx_scale <- fovea_coords[2]*(dim(amodel)[2])
  #fovea_coordy <- all_params_odd[[x,4]]
  fovea_coordy_scale <- fovea_coords[1]*(dim(amodel)[1])
  fovea_coords <- c(fovea_coordy_scale, fovea_coordx_scale)
  #################
  tmp_mat <- create_matrix(0, ((dim(amodel)[2])-1), 0, ((dim(amodel)[1])-1))
  rownames(tmp_mat) <- paste(tmp_mat[,1], ".", tmp_mat[,2], sep="")
  rotate_value <- rotate_val
  #scale_centre_ids <- scale_and_centre_rotate(amodel, fovea_height_scale, fovea_coords, n)
  #rotated_matrix <- rotate_matrix(amodel, rotate_value)
  #rotate_value <- 90
  rotated_matrix <- rotate_matrix_p(tmp_mat, rotate_val, amodel)
  
  
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

scale_and_centre_rotate_p <- function(amodel, fovea_coords, n, scale_f){
  fovea_coords[1] <- fovea_coords[1]/dim(amodel)[2]
  fovea_coords[2] <- fovea_coords[2]/dim(amodel)[1]
  tmp_mat <- create_matrix(0, ((dim(amodel)[2])-1), 0, ((dim(amodel)[1])-1))
  rownames(tmp_mat) <- paste(tmp_mat[,1], ".", tmp_mat[,2], sep="")
  id_mapping <- data.frame(id=paste(tmp_mat[,1], ".", tmp_mat[,2], sep=""), x=tmp_mat[,1], y=tmp_mat[,2], stringsAsFactors = FALSE)
  fovea_coords_big <- c(fovea_coords[1]*dim(amodel)[2], fovea_coords[2]*dim(amodel)[1])
  scale_co <- scale_f
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

rotate_matrix_p <- function(tst_matrix, rotate_val, amodel){
  #tst_matrix <- create_matrix(0, (dim(amodel)[2])-1, 0, (dim(amodel)[1])-1)
  #rownames(tst_matrix) <- paste( tst_matrix[,1], tst_matrix[,2], sep=".")
  id_mapping <- data.frame(id=paste(tst_matrix[,1], ".", tst_matrix[,2], sep=""), x=tst_matrix[,1], y=tst_matrix[,2], stringsAsFactors = FALSE)
  tst_matrix[,1] <- tst_matrix[,1]-(dim(amodel)[2])%/%2
  tst_matrix[,2] <- tst_matrix[,2]-(dim(amodel)[1])%/%2
  rotate_value_clock_rad <- rotate_val
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

procrustes_rotation_and_scale <- function(rotatey_eye_index, target_eye_es, smoothf = 30){
  rotatey_eye <- all_params_odd_qc[rotatey_eye_index,1][[1]]
  target_max_slope <- target_eye_es$max_slope
  target_plateau <- target_eye_es$plateau
  e <-registar(rotatey_eye)
  ellipse_max_slope <- smooth_ellipse_max_slope(e, smoothf)
  ellipse_plateu <- smooth_ellipse_plateu(e, smoothf)
  pro_max_slope <- procrustes(target_max_slope, ellipse_max_slope)
  pro_plateu <- procrustes(target_plateau, ellipse_plateu)
  max_slope_angle <- acos(pro_max_slope$rotation[1,1])
  plateu_angle <- acos(pro_plateu$rotation[1,1])
  return(data.frame(rotatey_eye = rotatey_eye, max_slope_angle=max_slope_angle, plateu_angle = plateu_angle, max_slope_scale = pro_max_slope$scale, plateu_scale = pro_plateu$scale))
}


######################################
###############ANALYSIS###############
######################################



target_eye1 <- registar(all_params_odd[1,1][[1]])
target_eye1 <- smooth_ellipse_max_slope(target_eye1, 30)
target_eye2 <- smooth_ellipse_plateu(target_eye1, 30)
target_eye_es <- list(max_slope=target_eye1, plateau=target_eye2)


all_params_even_qc <- all_params_even[quality_keep_joint,]
all_params_odd_qc <- all_params_odd[quality_keep_joint,]




procrustes_outcome_test <- lapply(2:dim(all_params_odd_qc)[1], procrustes_rotation_and_scale, target_eye_es=target_eye_es, smoothf=30)
procrustes_outcome_test_odds <- do.call(rbind, procrustes_outcome_test)
#procrustes_outcome_test_odds$max_slope_angle <- rad2deg(procrustes_outcome_test_odds$max_slope_angle)
#procrustes_outcome_test_odds$plateu_angle <- rad2deg(procrustes_outcome_test_odds$plateu_angle)
procrustes_outcome_test_oddsx <- data.frame(rotatey_eye = all_params_odd_qc[1,1], max_slope_angle = 0, plateu_angle = 0, max_slope_scale = 0, plateu_scale = 0)
colnames(procrustes_outcome_test_oddsx) <- c('rotatey_eye', 'max_slope_angle', 'plateu_angle', 'max_slope_scale', 'plateu_scale')
procrustes_outcome_test_oddsy <- rbind(procrustes_outcome_test_odds, procrustes_outcome_test_oddsx)

procrustes_outcome_test <- lapply(2:dim(all_params_even_qc)[1], procrustes_rotation_and_scale, target_eye_es=target_eye_es, smoothf=30)
procrustes_outcome_test_evens <- do.call(rbind, procrustes_outcome_test)
#procrustes_outcome_test_evens$max_slope_angle <- rad2deg(procrustes_outcome_test_evens$max_slope_angle)
#procrustes_outcome_test_evens$plateu_angle <- rad2deg(procrustes_outcome_test_evens$plateu)
procrustes_outcome_test_evensx <- data.frame(rotatey_eye = all_params_even_qc[1,1], max_slope_angle = 0, plateu_angle = 0, max_slope_scale = 0, plateu_scale = 0)
colnames(procrustes_outcome_test_evensx) <- c('rotatey_eye', 'max_slope_angle', 'plateu_angle', 'max_slope_scale', 'plateu_scale')
procrustes_outcome_test_evensy <- rbind(procrustes_outcome_test_evens, procrustes_outcome_test_evensx)


odds_p <- lapply(1:dim(all_params_odd_qc)[1], function(x){
  image_file <- all_params_odd_qc[[x,1]]
  scale_f <- procrustes_outcome_test_oddsy$max_slope_scale[[x]]
  fovea_coordx <- all_params_odd_qc[[x,3]]
  fovea_coordy <- all_params_odd_qc[[x,4]]
  fovea_coords <- c(fovea_coordy, fovea_coordx)
  rotate_val <- procrustes_outcome_test_oddsy$max_slope_angle[[x]]
  outputz <- rotate_and_plot_p(image_file, 0.8, 5, scale_f, fovea_coords, rotate_val, 20)
  return(outputz)
})
saveRDS(odds_p, paste("/Users/currant/OCTeyes/odds_fix_final_P",  strftime(Sys.time(), "%Y%m%d"), sep=""))
#odds_p <- readRDS("/Users/currant/OCTeyes/odds_fix_final20170523")

evens_p <- lapply(1:dim(all_params_even_qc)[1], function(x){
  image_file <- all_params_even_qc[[x,1]]
  scale_f <- procrustes_outcome_test_evensy$max_slope_scale[[x]]
  fovea_coordx <- all_params_even_qc[[x,3]]
  fovea_coordy <- all_params_even_qc[[x,4]]
  fovea_coords <- c(fovea_coordy, fovea_coordx)
  rotate_val <- procrustes_outcome_test_evensy$max_slope_angle[[x]]
  outputz <- rotate_and_plot_p(image_file, 0.8, 5, scale_f, fovea_coords, rotate_val, 20)
  return(outputz)
})
saveRDS(evens_p, paste("/Users/currant/OCTeyes/evens_fix_final_P", strftime(Sys.time(), "%Y%m%d"), sep=""))


odds_scale_p <- lapply(1:dim(all_params_odd_qc)[1], function(x){
  image_file <- all_params_odd_qc[[x,1]]
  scale_f <- procrustes_outcome_test_oddsy$max_slope_scale[[x]]
  fovea_coordx <- all_params_odd_qc[[x,3]]
  fovea_coordy <- all_params_odd_qc[[x,4]]
  fovea_coords <- c(fovea_coordy, fovea_coordx)
  outputz <- just_scale_and_plot_p(image_file, 0.8, 5, scale_f, fovea_coords, 20)
  return(outputz)
})
saveRDS(odds_scale_p, paste("/Users/currant/OCTeyes/odds_scale_fix_final_P",  strftime(Sys.time(), "%Y%m%d"), sep=""))
#odds_scale <- readRDS("/Users/currant/OCTeyes/odds_scale_fix_final20170523")

evens_scale_p <- lapply(1:dim(all_params_even_qc)[1], function(x){
  image_file <- all_params_even_qc[[x,1]]
  scale_f <- procrustes_outcome_test_evensy$max_slope_scale[[x]]
  fovea_coordx <- all_params_even_qc[[x,3]]
  fovea_coordy <- all_params_even_qc[[x,4]]
  fovea_coords <- c(fovea_coordy, fovea_coordx)
  outputz <- just_scale_and_plot_p(image_file, 0.8, 5, scale_f, fovea_coords, 20)
  return(outputz)
})
saveRDS(evens_scale_p, paste("/Users/currant/OCTeyes/evens_scale_fix_final_P",  strftime(Sys.time(), "%Y%m%d"), sep=""))
#evens_scale <- readRDS("/Users/currant/OCTeyes/evens_scale_fix_final20170523")

#############JUST_ROTATE##############
odds_rotate_p <- lapply(1:dim(all_params_odd_qc)[1], function(x){
  image_file <- all_params_odd_qc[[x,1]]
  fovea_coordx <- all_params_odd_qc[[x,3]]
  fovea_coordy <- all_params_odd_qc[[x,4]]
  fovea_coords <- c(fovea_coordy, fovea_coordx)
  rotate_val <- procrustes_outcome_test_oddsy$max_slope_angle[[x]]
  outputz <- rotate_and_plot_no_scale_p(image_file, 0.8, 5, fovea_coords, rotate_val, 20)
  return(outputz)
})
saveRDS(odds_rotate_p, paste("/Users/currant/OCTeyes/odds_rotate_fix_final_P",  strftime(Sys.time(), "%Y%m%d"), sep=""))
#odds_rotate <- readRDS("/Users/currant/OCTeyes/odds_rotate_fix_final20170523")

evens_rotate_p <- lapply(1:dim(all_params_even_qc)[1], function(x){
  image_file <- all_params_even_qc[[x,1]]
  fovea_coordx <- all_params_even_qc[[x,3]]
  fovea_coordy <- all_params_even_qc[[x,4]]
  fovea_coords <- c(fovea_coordy, fovea_coordx)
  rotate_val <- procrustes_outcome_test_evensy$max_slope_angle[[x]]
  outputz <- rotate_and_plot_no_scale_p(image_file, 0.8, 5, fovea_coords, rotate_val, 20)
  return(outputz)
})
saveRDS(evens_rotate_p, paste("/Users/currant/OCTeyes/evens_rotate_fix_final_P",  strftime(Sys.time(), "%Y%m%d"), sep=""))
#evens_rotate <- readRDS("/Users/currant/OCTeyes/evens_rotate_fix_final")



get_correls <- function(x, comp_table){
  table <- comp_table[[x]]
  table <- na.omit(table)
  correl <- cor(table[,3], table[,4])
  output <- data.frame(x, correl)
  return(output)
}

########Median normalise all of the values##############
odds_minus_normalised_p <- lapply(1:length(odds_p), median_normalise_intensity, datax=odds_p)
evens_minus_normalised_p <- lapply(1:length(evens_p), median_normalise_intensity, datax=evens_p)
odds_scale_normalised_p <- lapply(1:length(odds_scale_p), median_normalise_intensity, datax=odds_scale_p)
evens_scale_normalised_p <- lapply(1:length(evens_scale_p), median_normalise_intensity, datax=evens_scale_p)
odds_rotate_normalised_p <- lapply(1:length(odds_rotate_p), median_normalise_intensity, datax=odds_rotate_p)
evens_rotate_normalised_p <- lapply(1:length(evens_rotate_p), median_normalise_intensity, datax=evens_rotate_p)
odds_tesselate_normalised_p <- lapply(1:length(odds_tesselate), median_normalise_intensity, datax=odds_tesselate)
evens_tesselate_normalised_p <- lapply(1:length(evens_tesselate), median_normalise_intensity, datax=evens_tesselate)

#########add standard triangles to all processed odd images########
odds_minus_norm_stdtri_p <- lapply(1:length(odds_minus_normalised_p), add_standard_triangle_all, datax_norm = odds_minus_normalised_p, datax=odds_p)
evens_minus_norm_stdtri_p <- lapply(1:length(evens_minus_normalised_p), add_standard_triangle_all, datax_norm = evens_minus_normalised_p, datax=evens_p)
odds_scale_norm_stdtri_p <- lapply(1:length(odds_scale_normalised_p), add_standard_triangle_all, datax_norm = odds_scale_normalised_p, datax=odds_scale_p)
evens_scale_norm_stdtri_p <- lapply(1:length(evens_scale_normalised_p), add_standard_triangle_all, datax_norm = evens_scale_normalised_p, datax=evens_scale_p)
odds_rotate_norm_stdtri_p <- lapply(1:length(odds_rotate_normalised_p), add_standard_triangle_all, datax_norm = odds_rotate_normalised_p, datax = odds_rotate_p)
evens_rotate_norm_stdtri_p <- lapply(1:length(evens_rotate_normalised_p), add_standard_triangle_all, datax_norm = evens_rotate_normalised_p, datax=evens_rotate_p)
odds_tesselate_norm_stdtri_p <- lapply(1:length(odds_tesselate_normalised_p), add_standard_triangle_all, datax_norm = odds_tesselate_normalised_p, datax=odds_tesselate)
evens_tesselate_norm_stdtri_p <- lapply(1:length(evens_tesselate_normalised_p), add_standard_triangle_all, datax_norm = evens_tesselate_normalised_p, datax=evens_tesselate)


get_max_min_tri_num <- function(stdtri_table_index, stdrtri_table_list){
  tris <- stdrtri_table_list[[stdtri_table_index]][,4]
  un_tris <- unique(tris)
  min <- min(un_tris)
  max <- max(un_tris)
  return(data.frame(min, max))
}

mandm_odds_minus <- lapply(1:length(odds_minus_norm_stdtri_p), get_max_min_tri_num, odds_minus_norm_stdtri)
mandm_odds_minus <- do.call(rbind, mandm_odds_minus)
final_min <- min(mandm_odds_minus[,1])
final_max <- max(mandm_odds_minus[,2])


#########get reduced version of table#############
redtab_oddsminus_norm_p <- lapply(final_min:final_max, get_whole_reduced_table, stdtritab = odds_minus_norm_stdtri_p)
redtab_evensminus_norm_p <- lapply(final_min:final_max, get_whole_reduced_table, stdtritab = evens_minus_norm_stdtri_p)
redtab_oddsscale_norm_p <- lapply(final_min:final_max, get_whole_reduced_table, stdtritab = odds_scale_norm_stdtri_p)
redtab_evensscale_norm_p <- lapply(final_min:final_max, get_whole_reduced_table, stdtritab = evens_scale_norm_stdtri_p)
redtab_oddsrotate_norm_p <- lapply(final_min:final_max, get_whole_reduced_table, stdtritab = odds_rotate_norm_stdtri_p)
redtab_evensrotate_norm_p <- lapply(final_min:final_max, get_whole_reduced_table, stdtritab = evens_rotate_norm_stdtri_p)
redtab_oddstesselate_norm_p <- lapply(final_min:final_max, get_whole_reduced_table, stdtritab = odds_tesselate_norm_stdtri_p)
redtab_evenstesselate_norm_p <- lapply(final_min:final_max, get_whole_reduced_table, stdtritab = evens_tesselate_norm_stdtri_p)

tess_v_tess_p <- lapply((final_min+1):(final_max+1), comp_per_tri, comp1=redtab_evenstesselate_norm_p, comp2=redtab_oddstesselate_norm_p)
all_v_all_p <- lapply((final_min+1):(final_max+1), comp_per_tri, comp1=redtab_oddsminus_norm_p, comp2=redtab_evensminus_norm_p)
scale_v_scale_p <- lapply((final_min+1):(final_max+1), comp_per_tri, comp1=redtab_oddsscale_norm_p, comp2=redtab_evensscale_norm_p)
rotate_v_rotate_p <- lapply((final_min+1):(final_max+1), comp_per_tri, comp1=redtab_oddsrotate_norm_p, comp2=redtab_evensrotate_norm_p)

tess_v_tess_norm_corr_p <- lapply((final_min+1):(final_max+1), get_correls, comp_table = tess_v_tess_p)
tess_v_tess_corr_p <- do.call(rbind, tess_v_tess_norm_corr_p)

all_v_all_norm_corr_p <- lapply((final_min+1):(final_max+1), get_correls, comp_table = all_v_all_p)
all_v_all_corr_p <- do.call(rbind, all_v_all_norm_corr_p)

rotate_v_rotate_norm_corr_p <- lapply((final_min+1):(final_max+1), get_correls, comp_table = rotate_v_rotate_p)
rotate_v_rotate_corr_p <- do.call(rbind, rotate_v_rotate_norm_corr_p)

scale_v_scale_norm_corr_p <- lapply((final_min+1):(final_max+1), get_correls, comp_table = scale_v_scale_p)
scale_v_scale_corr_p <- do.call(rbind, scale_v_scale_norm_corr_p)

#############Defines the refernce plot ready for plotting back to triangles###############
ref_points <- evens_p[[1]]$all_refs

############Gets table of info to be able to plot back to triangles################
tess_v_tess_cor_triplot_p <- lapply(1:dim(tess_v_tess_corr_p)[1], get_triplot_table, ref_points = ref_points, corr_tab = tess_v_tess_corr_p)
tess_v_tess_cor_triplot_p <- do.call(rbind, tess_v_tess_cor_triplot_p)
tess_v_tess_cor_triplot_p <- as.data.frame(tess_v_tess_cor_triplot_p)

scale_v_scale_cor_triplot_p <- lapply(1:dim(scale_v_scale_corr_p)[1], get_triplot_table, ref_points = ref_points, corr_tab = scale_v_scale_corr_p)
scale_v_scale_cor_triplot_p <- do.call(rbind, scale_v_scale_cor_triplot_p)
scale_v_scale_cor_triplot_p <- as.data.frame(scale_v_scale_cor_triplot_p)

rotate_v_rotate_cor_triplot_p <- lapply(1:dim(rotate_v_rotate_corr_p)[1], get_triplot_table, ref_points = ref_points, corr_tab = rotate_v_rotate_corr_p)
rotate_v_rotate_cor_triplot_p <- do.call(rbind, rotate_v_rotate_cor_triplot_p)
rotate_v_rotate_cor_triplot_p <- as.data.frame(rotate_v_rotate_cor_triplot_p)

all_v_all_cor_triplot_p <- lapply(1:dim(all_v_all_corr_p)[1], get_triplot_table, ref_points = ref_points, corr_tab = all_v_all_corr_p)
all_v_all_cor_triplot_p <- do.call(rbind, all_v_all_cor_triplot_p)
all_v_all_cor_triplot_p <- as.data.frame(all_v_all_cor_triplot_p)

plot(tess_v_tess_corr_p, xlab="Triangle ID", ylab="Correlation", main = "Tess vs tess per triangle procrustes")
hist(tess_v_tess_corr_p[,2], xlab="Correlation", main="Tess vs tess procrustes")

plot(scale_v_scale_corr_p, xlab="Triangle ID", ylab="Correlation", main="Scale vs scale per triangle procrustes")
hist(scale_v_scale_corr_p[,2], xlab="Correlation", main="Scale vs Scale procrustes")

plot(rotate_v_rotate_corr_p, xlab="Triangle ID", ylab="Correlation", main = "Rotate vs rotate per triangle procrustes")
hist(rotate_v_rotate_corr_p[,2], xlab="Correlation", main="Rotate vs rotate procrustes")

plot(all_v_all_corr_p, xlab = "Triangle ID", ylab="Correlation", main="All vs all per triangle procrustes")
hist(all_v_all_corr_p[,2], xlab="Correlation", main = "All vs all procrustes")


all_v_all_corr_p$type <- 'all'
tess_v_tess_corr_p$type <- 'tess'
scale_v_scale_corr_p$type <- 'scale'
rotate_v_rotate_corr_p$type <- 'rotate'

df_transforms <- rbind(tess_v_tess_corr_p, scale_v_scale_corr_p, rotate_v_rotate_corr_p, all_v_all_corr_p)

ggplot(df_transforms, aes(correl, fill=type))+ geom_histogram(alpha=0.5, aes(y=..density..), position = "identity")
ggplot(df_transforms, aes(correl, fill=type))+geom_density(alpha=0.2) + labs(x = 'Correlation', y='Density')


procrustes_big_table <- rbind(all_v_all_corr_p, tess_v_tess_corr_p, scale_v_scale_corr_p, rotate_v_rotate_corr_p)
procrustes_big_table$qc <- "Procrustes Method"







##############################################################

smoothf <-30
e1 <- registar(all_params_odd[1,1][[1]])
Xellipse1 <- smooth_ellipse_max_slope(e1, smoothf)
Xellipse2 <- smooth_ellipse_plateu(e1, smoothf)


e2 <- registar(all_params_even[1,1][[1]])
Yellipse1 <- smooth_ellipse_max_slope(e2, smoothf)
Yellipse2 <- smooth_ellipse_plateu(e2, smoothf)

pro_max_slope <- procrustes(Xellipse1, Yellipse1)
pro_max_slope
summary(pro_max_slope)
plot(pro_max_slope)
plot(pro_max_slope, kind=2)
angle_max_slope1 <- acos(pro_max_slope$rotation[1,1])
scale_max_slope <- pro_max_slope$scale


pro_plateu <- procrustes(Xellipse2, Yellipse2)
pro_plateu
summary(pro_plateu)
plot(pro_plateu)
plot(pro_plateu, kind=2)
angle_plateau1 <- acos(pro_plateu$rotation[1,1])




X <- evens_tesselate[[2]]$end_info2[,1:2]
Y <- odds_tesselate[[2]]$end_info2[,1:2]
X[is.na(X)] <- 0
Y[is.na(Y)] <- 0
test_pro <- procrustes(X,Y)
 
summary(test_pro)
plot(test_pro)
test_pro$rotation
l <- asin(test_pro$rotation[1,2])
sin(l)
rad2deg(l)



vare.dist <- vegdist(wisconsin(varespec))
library(MASS)  ## isoMDS
mds.null <- isoMDS(vare.dist, tol=1e-7)
mds.alt <- isoMDS(vare.dist, initMDS(vare.dist), maxit=200, tol=1e-7)
vare.proc <- procrustes(mds.alt, mds.null)
vare.proc
summary(vare.proc)
plot(vare.proc)
plot(vare.proc, kind=2)
residuals(vare.proc)

