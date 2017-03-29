
# models the surface of the top and bottom major bands acorss oct scan slices
# @ surface modelling

#library(foveafit)

run_through <- function() {

	image_file = "/Users/tomas/projects/2016/eye_imaging/30_replicates/1853843/OCT/333879.dcm"
	image_data = dicomInfo(image_file)

	model_top_surface = NULL
	model_bottom_surface = NULL
	pb <- txtProgressBar(min = 0, max = dim(image_data$img)[3], style = 3);
	for(x in 1:dim(image_data$img)[3]) {
	    bands = fit_major_bands(image_data, x)
	    ibands = interpolate_outliers(bands)
	    model_bottom_surface = rbind(model_bottom_surface, ibands[,1])
	    model_top_surface = rbind(model_top_surface, ibands[,2])
	    setTxtProgressBar(pb, x)
	} 
	close(pb);

	majorbs = list("top_band"= model_top_surface, "bottom_band"=model_bottom_surface) 

}


fit.ellipse <- function (x, y = NULL) {
  # from:
  # http://r.789695.n4.nabble.com/Fitting-a-half-ellipse-curve-tp2719037p2720560.html
  #
  # Least squares fitting of an ellipse to point data
  # using the algorithm described in: 
  #   Radim Halir & Jan Flusser. 1998. 
  #   Numerically stable direct least squares fitting of ellipses. 
  #   Proceedings of the 6th International Conference in Central Europe 
  #   on Computer Graphics and Visualization. WSCG '98, p. 125-132 
  #
  # Adapted from the original Matlab code by Michael Bedward (2010)
  # michael.bedward@gmail.com
  #
  # Subsequently improved by John Minter (2012)
  # 
  # Arguments: 
  # x, y - x and y coordinates of the data points.
  #        If a single arg is provided it is assumed to be a
  #        two column matrix.
  #
  # Returns a list with the following elements: 
  #
  # coef - coefficients of the ellipse as described by the general 
  #        quadratic:  ax^2 + bxy + cy^2 + dx + ey + f = 0 
  #
  # center - center x and y
  #
  # major - major semi-axis length
  #
  # minor - minor semi-axis length
  #
  EPS <- 1.0e-8 
  dat <- xy.coords(x, y) 
  
  D1 <- cbind(dat$x * dat$x, dat$x * dat$y, dat$y * dat$y) 
  D2 <- cbind(dat$x, dat$y, 1) 
  S1 <- t(D1) %*% D1 
  S2 <- t(D1) %*% D2 
  S3 <- t(D2) %*% D2 
  T <- -solve(S3) %*% t(S2) 
  M <- S1 + S2 %*% T 
  M <- rbind(M[3,] / 2, -M[2,], M[1,] / 2) 
  evec <- eigen(M)$vec 
  cond <- 4 * evec[1,] * evec[3,] - evec[2,]^2 
  a1 <- evec[, which(cond > 0)] 
  f <- c(a1, T %*% a1) 
  names(f) <- letters[1:6] 
  
  # calculate the center and lengths of the semi-axes 
  #
  # see http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2288654/
  # J. R. Minter
  # for the center, linear algebra to the rescue
  # center is the solution to the pair of equations
  # 2ax +  by + d = 0
  # bx  + 2cy + e = 0
  # or
  # | 2a   b |   |x|   |-d|
  # |  b  2c | * |y| = |-e|
  # or
  # A x = b
  # or
  # x = Ainv b
  # or
  # x = solve(A) %*% b
  A <- matrix(c(2*f[1], f[2], f[2], 2*f[3]), nrow=2, ncol=2, byrow=T )
  b <- matrix(c(-f[4], -f[5]), nrow=2, ncol=1, byrow=T)
  soln <- solve(A) %*% b

  b2 <- f[2]^2 / 4
  
  center <- c(soln[1], soln[2]) 
  names(center) <- c("x", "y") 
  
  num  <- 2 * (f[1] * f[5]^2 / 4 + f[3] * f[4]^2 / 4 + f[6] * b2 - f[2]*f[4]*f[5]/4 - f[1]*f[3]*f[6]) 
  den1 <- (b2 - f[1]*f[3]) 
  den2 <- sqrt((f[1] - f[3])^2 + 4*b2) 
  den3 <- f[1] + f[3] 
  
  semi.axes <- sqrt(c( num / (den1 * (den2 - den3)),  num / (den1 * (-den2 - den3)) )) 
  
  # calculate the angle of rotation 
  term <- (f[1] - f[3]) / f[2] 
  angle <- atan(1 / term) / 2 
  
  list(coef=f, center = center, major = max(semi.axes), minor = min(semi.axes), angle = unname(angle)) 
}


get.ellipse <- function( fit, n=360 ) 
{
  # Calculate points on an ellipse described by 
  # the fit argument as returned by fit.ellipse 
  # 
  # n is the number of points to render 
  
  tt <- seq(0, 2*pi, length=n) 
  sa <- sin(fit$angle) 
  ca <- cos(fit$angle) 
  ct <- cos(tt) 
  st <- sin(tt) 
  
  x <- fit$center[1] + fit$maj * ct * ca - fit$min * st * sa 
  y <- fit$center[2] + fit$maj * ct * sa + fit$min * st * ca 
  
  cbind(x=x, y=y) 
}

Deriv <- function(mod, n = 200, eps = 1e-7, newdata, term) {
    if(inherits(mod, "gamm"))
        mod <- mod$gam
    m.terms <- attr(terms(mod), "term.labels")
    if(missing(newdata)) {
        newD <- sapply(model.frame(mod)[, m.terms, drop = FALSE],
                       function(x) seq(min(x), max(x), length = n))
        names(newD) <- m.terms
    } else {
        newD <- newdata
    }
    newDF <- data.frame(newD) ## needs to be a data frame for predict
    X0 <- predict(mod, newDF, type = "lpmatrix")
    #newDFeps_m <- newDF - eps
    newDF <- newDF + eps
    X1 <- predict(mod, newDF, type = "lpmatrix")
    #X_1 <- predict(mod, newDFeps_m, type = 'lpmatrix')
    Xp <- (X1 - X0) / eps
    #Xpp <- (X1 + X_1 - 2*X0)  / eps^2
    Xp.r <- NROW(Xp)
    Xp.c <- NCOL(Xp)
    ## dims of bs
    bs.dims <- sapply(mod$smooth, "[[", "bs.dim") - 1
    ## number of smooth terms
    t.labs <- attr(mod$terms, "term.labels")
    ## match the term with the the terms in the model
    if(!missing(term)) {
        want <- grep(term, t.labs)
        if(!identical(length(want), length(term)))
            stop("One or more 'term's not found in model!")
        t.labs <- t.labs[want]
    }
    nt <- length(t.labs)
    ## list to hold the derivatives
    lD <- vector(mode = "list", length = nt)
    names(lD) <- t.labs
    for(i in seq_len(nt)) {
        Xi <- Xp * 0
        #Xii <- Xp * 0
        want <- grep(t.labs[i], colnames(X1))
        Xi[, want] <- Xp[, want]
        df <- Xi %*% coef(mod)
        df.sd <- rowSums(Xi %*% mod$Vp * Xi)^.5

        #Xii[, want] <- Xpp[, want]
        #fd_d2 <- Xii %*% coef(mod)

        lD[[i]] <- list(deriv = df, se.deriv = df.sd)
    }
    class(lD) <- "Deriv"
    lD$gamModel <- mod
    lD$eps <- eps
    lD$eval <- newD - eps
    lD ##return
}

confint.Deriv <- function(object, term, alpha = 0.05, ...) {
    l <- length(object) - 3
    term.labs <- names(object[seq_len(l)])
    if(missing(term)) {
        term <- term.labs
    } else { ## how many attempts to get this right!?!?
        ##term <- match(term, term.labs)
        ##term <- term[match(term, term.labs)]
        term <- term.labs[match(term, term.labs)]
    }
    if(any(miss <- is.na(term)))
        stop(paste("'term'", term[miss], "not a valid model term."))
    res <- vector(mode = "list", length = length(term))
    names(res) <- term
    residual.df <- df.residual(object$gamModel)
    tVal <- qt(1 - (alpha/2), residual.df)
    ##for(i in term.labs[term]) {
    for(i in term) {
        upr <- object[[i]]$deriv + tVal * object[[i]]$se.deriv
        lwr <- object[[i]]$deriv - tVal * object[[i]]$se.deriv
        res[[i]] <- list(upper = drop(upr), lower = drop(lwr))
    }
    res$alpha = alpha
    res
}

signifD <- function(x, d, upper, lower, eval = 0) {
    miss <- upper > eval & lower < eval
    incr <- decr <- x
    want <- d > eval
    incr[!want | miss] <- NA
    want <- d < eval
    decr[!want | miss] <- NA
    list(incr = incr, decr = decr)
}

get_derivative <- function(data, eps= 1e-7) {
	time = 1:length(data)
	value = data
	z = data.frame(time, value)
	a <- gam(value~s(time),data=z)
	newDF <- with(z, data.frame(time = unique(time)))
	B <- predict(a,  newDF, type="response", se.fit=TRUE)
	X0 <- predict(a, newDF, type = 'lpmatrix')
	newDFeps_p <- newDF + eps
	X1 <- predict(a, newDFeps_p, type = 'lpmatrix')
	Xp <- (X0 - X1) / eps
	fd_d1 <- Xp %*% coef(a)
return(fd_d1)
}

rescale <- function (x, newrange) {
    if (missing(x) | missing(newrange)) {
        usage.string <- paste("Usage: rescale(x,newrange)\n", 
            "\twhere x is a numeric object and newrange is the new min and max\n", 
            sep = "", collapse = "")
        stop(usage.string)
    }
    if (is.numeric(x) && is.numeric(newrange)) {
        xna <- is.na(x)
        if (all(xna)) 
            return(x)
        if (any(xna)) 
            xrange <- range(x[!xna])
        else xrange <- range(x)
        if (xrange[1] == xrange[2]) 
            return(x)
        mfac <- (newrange[2] - newrange[1])/(xrange[2] - xrange[1])
        return(newrange[1] + (x - xrange[1]) * mfac)
    }
    else {
        warning("Only numeric objects can be rescaled")
        return(x)
    }
}

test_rescale <- function(centered_band) {
	image(t(centered_band), col=topo.colors(10))
 	fac = dim(centered_band)[1]/180
 	
	inds = x.xs[,2:3]/dim(centered_band)[1]
	inds = inds[nrow(inds):1,]


}


grad.des <- function() {
		# get data 
	rm(list = ls(all = TRUE)) # make sure previous work is clear
	ls()
	x0 <- c(1,1,1,1,1) # column of 1's
	x1 <- c(1,2,3,4,5) # original x-values
	 
	# create the x- matrix of explanatory variables
	 
	x <- as.matrix(cbind(x0,x1))
	 
	# create the y-matrix of dependent variables
	 
	y <- as.matrix(c(3,7,5,11,14))
	m <- nrow(y)
	 
	# implement feature scaling
	x.scaled <- x
	x.scaled[,2] <- (x[,2] - mean(x[,2]))/sd(x[,2])
	 
	# analytical results with matrix algebra
	 solve(t(x)%*%x)%*%t(x)%*%y # w/o feature scaling
	 solve(t(x.scaled)%*%x.scaled)%*%t(x.scaled)%*%y # w/ feature scaling
	 
	# results using canned lm function match results above
	summary(lm(y ~ x[, 2])) # w/o feature scaling
	summary(lm(y ~ x.scaled[, 2])) # w/feature scaling
	 
	# define the gradient function dJ/dtheata: 1/m * (h(x)-y))*x where h(x) = x*theta
	# in matrix form this is as follows:
	grad <- function(x, y, theta) {
	  gradient <- (1/m)* (t(x) %*% ((x %*% t(theta)) - y))
	  return(t(gradient))
	}
	 
	# define gradient descent update algorithm
	grad.descent <- function(x, maxit){
	    theta <- matrix(c(0, 0), nrow=1) # Initialize the parameters
	 
	    alpha = .05 # set learning rate
	    for (i in 1:maxit) {
	      theta <- theta - alpha  * grad(x, y, theta)   
	    }
	 return(theta)
	}
	 
	 
	# results without feature scaling
	print(grad.descent(x,1000))
	 
	# results with feature scaling
	print(grad.descent(x.scaled,1000))



}


index_circle <- function(mat, center, radius) {
  grid <- mat
  x.index <- center + radius * cos(seq(0, 2*pi, length = 360))
  y.index <- center + radius * sin(seq(0, 2*pi, length = 360))
 	for (i in seq(x.index)) {    
 		grid[x.index[i], y.index[i]] <- NA
 	}
  image(grid)
}  # end drawImage


dev_register <- function(file) {
	 model = run_generate_model_and_find_fovea_dip_poisition(file)
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
	 centered_band = model$model$top_band[rangs[[1]], rangs[[2]]]
	 centered_band = make_square(centered_band)
 	 grid <- centered_band
 	 diam <- unique(dim(centered_band))
 	 center <- radius <- diam/2
   coors = data.frame(center + radius * cos(seq(0, 2*pi, length = 360))
   										,  center + radius * sin(seq(0, 2*pi, length = 360))) 
   grand = diam/180
   x.xs = NULL
   for (x in 1:length(coors[,1])) {
   	x0 = floor(grand*x)
   	x1 = coors[x,1]
   	y0 = dim(centered_band)[1]-coors[x,2]
   	y1 = coors[x,2]
   	x.xs = rbind(x.xs, cbind(x0, x1, y0, y1))
   }
   x.xs = x.xs[x.xs[,1]<=diam,]
   s.ss = list()
   for(x in 1:length(x.xs[,1])) {
   	se1 = round(seq(x.xs[x,1], x.xs[x,2], length=diam))
   	se2 = round(seq(x.xs[x,3], x.xs[x,4], length=diam))
   	se1[se1<1] = 1; se1[se1>diam] = diam;
		se2[se2<1] = 1; se2[se2>diam] = diam;
		ss = vector()
		 for(y in 1:length(se1)) {
				ss[y] = centered_band[se2[y],se1[y]] 
			}
		s.ss[[x]] = ss	
	}
	p.ps = NULL
	for(x in 1:length(s.ss)) {
		ps = fit_deriv(s.ss[[x]])
		p.ps = rbind(p.ps, unlist(ps))
	}	
	xy.s = NULL
	xy.ss = NULL
	for(a in 1:dim(p.ps)[1]) {
		cx <- cy <- (dim(centered_band)[1]/2)
		r = p.ps[a,1] - p.ps[a,2]
		x = cx + r * cos(a)
		y = cy + r * sin(a)
		xy.s = rbind(xy.s, cbind(x, y))
		a = (180-a)+1
		r = p.ps[a,3] - p.ps[a,1]
		cx <- cy <- (dim(centered_band)[1]/2)
		x = cx + r * cos(a)
		y = cy + r * sin(a)
		xy.s = rbind(xy.s, cbind(x, y))

		r = abs (p.ps[a,1]-(diam/2))
		x = cx + r * cos(a)
		y = cy + r * sin(a)
		xy.ss = rbind(xy.ss, cbind(x, y))
	}

	e = get.ellipse(fit.ellipse(xy.s[,1], xy.s[,2]))
	image(t(centered_band), col=topo.colors(10))
	matplot(x.xs[,1]/dim(centered_band)[1], x.xs[,3]/dim(centered_band)[1], add=T, pch="*")
	matplot(x.xs[,1]/dim(centered_band)[1], 1-x.xs[,3]/dim(centered_band)[1], add=T, pch="*")
	matplot(x.xs[,3]/dim(centered_band)[1], x.xs[,1]/dim(centered_band)[1], add=T, pch="*")
	matplot( 1-x.xs[,3]/dim(centered_band)[1], x.xs[,1]/dim(centered_band)[1], add=T, pch="*")
	matplot(x.xs[,2]/dim(centered_band)[1], x.xs[,3]/dim(centered_band)[1], add=T, pch="*")
	matplot(x.xs[,2]/dim(centered_band)[1], 1-x.xs[,3]/dim(centered_band)[1], add=T, pch="*")
	abline(h=0.5)
	abline(v=0.5)
	points(xy.s[,1]/dim(centered_band)[1], xy.s[,2]/dim(centered_band)[1], col="red", pch="+", cex=0.8)
	lines(e/dim(centered_band)[1], lwd=2, col="blue")
}




create_tessalation <- function(m, mm, fac=10, main="", brks = 10) {
	
	inds = NULL
	s1 = seq(1, nrow(m), by=fac)
	for(x in 1:length(s1)) {
		inds = rbind(inds, (s1))
	}
	inds = cbind(as.vector(inds), as.vector(t(inds)))
	dxy1 <- deldir(inds[,1], inds[,2])
	delsgs <- dxy1$delsgs/nrow(m)
	dirsgs = dxy1$dirsgs/nrow(m)
        x1 <- delsgs[, 1]
        y1 <- delsgs[, 2]
        x2 <- delsgs[, 3]
        y2 <- delsgs[, 4]
        u1 <- dirsgs[, 1]
        v1 <- dirsgs[, 2]
        u2 <- dirsgs[, 3]
        v2 <- dirsgs[, 4]
   image((mm), col=topo.colors(brks), main=main)     
  # segments(u1, v1, u2, v2, col = "red", lty = 1, lwd=0.2)
	# segments(x1, y1, x2, y2, col = "red", lty = 1, lwd=0.2)
}

square_tes_extract <- function(m, fac=10) {
	values = NULL
	m = m-median(m)
	s1 = seq(1, dim(m)[1], by=fac)
	s2 = seq(1, dim(m)[2], by=fac)
	for(x in 2:length(s1)) {
		m1 = m[s1[x-1]:s1[x],]
		s = sapply(2:length(s2), function(i) mean(m1[,s2[i-1]:s2[i]]))
		values = rbind(values, s)
	}
	#image(values)
return(values)
}


tri_tes_extract <- function(m, fac=10, debug=FALSE) {
	triangles_from_square <- function(m1) {
			tri1 = NULL; tri2 = NULL
			for(y in 1:nrow(m1)) {
				tri1 = c(tri1, m1[y, 1:y])
				tri2 = c(tri2, m1[y,1:(nrow(m1)-y)])
			}
		return(list("t1" =tri1, "t2"=tri1))
	}
	values = NULL
	s1 = seq(1, dim(m)[1], by=fac)
	s2 = seq(1, dim(m)[2], by=fac)
	for(x in 2:length(s1)) {
		m1 = m[s1[x-1]:s1[x],]
		s = sapply(2:length(s2), function(i) unlist(lapply(triangles_from_square(m1[,s2[i-1]:s2[i]]), mean)))
		values = rbind(values, s)
	}
	if(debug) {
		image((values), col=topo.colors(1000))
		abline(v=seq(0,1, length=length(s1)))
		abline(h=seq(0,1, length=length(s2)))
		#abline(v=s1/dim(m)[1], col="lightgrey")
		#abline(h=s2/dim(m)[2], col="lightgrey")
		for(x in 2:length(s1)) {
			for(y in 2:length(s2)) {
				#segments(s1[x-1]/dim(m)[1], s2[y-1]/dim(m)[2], s1[x]/dim(m)[1], s2[y]/dim(m)[2])
			}
		} 
	}
return(values)
}

align_model <- function(data, ns=c(75, 306)) {
	ndata = NULL
	ds = dim(data)
	st = ns/2
	ndata = data[((ds[1]/2)-st[1]): ((ds[1]/2)+st[1]), ((ds[2]/2)-st[2]): ((ds[2]/2)+st[2])]
return(ndata)
} 

center_align<- function(model) {
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
	centered_band = model$model$top_band[rangs[[1]], rangs[[2]]]
	centered_band = align_model(centered_band)
return(centered_band)
} 

test <- function() {

	par(mfrow=c(1,2))
	data_folder = "/Users/tomas/projects/2016/eye_imaging/30_replicates"
	file_dirs = paste(data_folder, dir(data_folder), "OCT", sep="/")
	es = NULL
	vs = list()
	c = vector()
	for(x in 1:length(file_dirs)) {	
		files = paste(file_dirs[x], dir(file_dirs[x]), sep="/")
		image_file = files[1]
		model = run_generate_model_and_find_fovea_dip_poisition(image_file)
		amodel = center_align(model)
		v1 = square_tes_extract(amodel, 5)
		#tes = create_tessalation(amodel, 5)
		m1 = dim(amodel)
		image_file = files[2]
		model = run_generate_model_and_find_fovea_dip_poisition(image_file)
		amodel = center_align(model)
		v2 = square_tes_extract(amodel, 5)
		#tes = create_tessalation(amodel, 5)
		m2 = dim(amodel)
		es = rbind(es, c(m1, m2))

		vs[[x]] = list()
		vs[[x]]$rep1 = v1
		vs[[x]]$rep2 = v2
		c[x] = cor(c(v1), c(v2))
		print(x)	
	}

}

extract_overlay <- function() {
	lis = list()
	debug=F
	average_model = NULL
	for(x in 1:length(file_dirs)) {	
		if(x!=25) {
		files = paste(file_dirs[x], dir(file_dirs[x]), sep="/")
		image_file = files[1]
		model = run_generate_model_and_find_fovea_dip_poisition(image_file)
		amodel = center_align(model)
		centered_band = make_square(amodel)
		fitted_center = fit_centered_surface(centered_band, debug)
		e1 = get.ellipse(fit.ellipse(fitted_center$max_slope[,1], fitted_center$max_slope[,2]))
		e2 = get.ellipse(fit.ellipse(fitted_center$plateu[,1], fitted_center$plateu[,2]))
		image(t(centered_band-median(centered_band)), col=topo.colors(10))
		lines(e1/dim(centered_band)[1], lwd=2, col="blue")
		lines(e2/dim(centered_band)[1], lwd=2, col="blue")
		matplot(fitted_center$max_slope[,1]/dim(centered_band), fitted_center$max_slope[,2]/dim(centered_band), col="red", add=T, pch="+", cex=0.5)
		matplot(fitted_center$plateu[,1]/dim(centered_band), fitted_center$plateu[,2]/dim(centered_band), col="red", add=T, pch="+", cex=0.5)
		li = list("e1"=e1, "e2"=e2, "model"=model$model$top_band, "ycenter"  = model$ycenter, "zcenter" = model$zcenter
			, "smodel"=centered_band, "max_slope"=fitted_center$max_slope, "plateu"=fitted_center$plateu)
		lis[[x]] = li

		if(is.null(average_model)) {
			average_model = centered_band-median(centered_band)
		} else {
			average_model = average_model + (centered_band-median(centered_band))
		}	
	}
	}

	#image(t(lis[[1]]$smodel-median(lis[[1]]$smodel)), col=topo.colors(10))
	image(t(average_model), col=topo.colors(10))
	for(x in 1:length(lis)) {
		lines(lis[[x]]$e1/dim(lis[[x]]$smodel)[1], lwd=2, col="red")
		lines(lis[[x]]$e2/dim(lis[[1]]$smodel)[1], lwd=2, col="blue")
		#matplot(lis[[x]]$max_slope[,1]/dim(lis[[1]]$smodel), lis[[x]]$max_slope[,2]/dim(lis[[x]]$smodel), col="grey", add=T, pch="+", cex=0.4)
		#matplot(lis[[x]]$plateu[,1]/dim(lis[[x]]$smodel), lis[[x]]$plateu[,2]/dim(lis[[x]]$smodel), col="red", add=T, pch="+", cex=0.5)

	}


parms = NULL
for(x in 1:length(lis)) {

	x1 = lis[[x]]$e1[which.min(lis[[x]]$e1[,1]),1]
	x2 = lis[[x]]$e1[which.max(lis[[x]]$e1[,1]),1]

	y1 = lis[[x]]$e1[which.min(lis[[x]]$e1[,2]),2]
	y2 = lis[[x]]$e1[which.max(lis[[x]]$e1[,2]),2]

	parms = rbind(parms, c(x1, x2, y1, y2))
}


 plot(c(parms1[,c(2)]-parms1[,c(1)]), c(parms2[,c(2)]-parms2[,c(1)]), pch="+", xlab="replicate1", ylab=" replicate2", main="x-axis length" )
abline(lm(c(parms2[,c(2)]-parms2[,c(1)])~c(parms1[,c(2)]-parms1[,c(1)])))
legend("topleft", c(paste("R^2: ", round(summary(lm(c(parms2[,c(2)]-parms2[,c(1)])~c(parms1[,c(2)]-parms1[,c(1)])))$r.squared,2) )))
 plot(c(parms1[,c(4)]-parms1[,c(3)]), c(parms2[,c(4)]-parms2[,c(3)]), pch="+", xlab="replicate1", ylab=" replicate2", main="y-axis length" )
abline(lm(c(parms2[,c(4)]-parms2[,c(3)])~c(parms1[,c(4)]-parms1[,c(3)])))
legend("topleft", c(paste("R^2: ", round(summary(lm(c(parms2[,c(4)]-parms2[,c(3)])~c(parms1[,c(4)]-parms1[,c(3)])))$r.squared,2) )))

parms = NULL
for(x in 1:length(lis)) {
	y1 = lis[[x]]$e2[which.min(lis[[x]]$e2[,2] - lis[[x]]$e2[,1]),2]
	y2 = lis[[x]]$e2[which.max(lis[[x]]$e1[,2] - lis[[x]]$e2[,1]),2]

	x1 = lis[[x]]$e2[which.min(lis[[x]]$e2[,2] - lis[[x]]$e2[,1]),1]
	x2 = lis[[x]]$e2[which.max(lis[[x]]$e2[,2] - lis[[x]]$e2[,1]),1]
	parms = rbind(parms, c(x1, x2, y1, y2))
}



}


get_correlations <- function(f) {
	par(mfrow=c(1,2))
	data_folder = "/Users/tomas/projects/2016/eye_imaging/30_replicates"
	file_dirs = paste(data_folder, dir(data_folder), "OCT", sep="/")
	es = NULL
	vs = list()
	c = vector()
	for(x in 1:length(file_dirs)) {	
		files = paste(file_dirs[x], dir(file_dirs[x]), sep="/")
		image_file = files[1]
		model = run_generate_model_and_find_fovea_dip_poisition(image_file)
		amodel = center_align(model)
		v1 = square_tes_extract(amodel, f)
		#tes = create_tessalation(amodel, 5)
		m1 = dim(amodel)
		image_file = files[2]
		model = run_generate_model_and_find_fovea_dip_poisition(image_file)
		amodel = center_align(model)
		v2 = square_tes_extract(amodel, f)
		#tes = create_tessalation(amodel, 5)
		m2 = dim(amodel)
		es = rbind(es, c(m1, m2))

		vs[[x]] = list()
		vs[[x]]$rep1 = v1
		vs[[x]]$rep2 = v2
		c[x] = cor(c(v1), c(v2))
		print(x)	
	}

	randco = NULL
	for(x in 1:length(vs)) {
		for(y in 1:length(vs)) {
			if(x!=y) {
				c = cor(c(vs[[x]]$rep1), c(vs[[y]]$rep1 ))
				randco = rbind(randco, cbind(c, x, y))
				c = cor(c(vs[[x]]$rep2), c(vs[[y]]$rep2 ))
				randco = rbind(randco, cbind(c, x, y))
				c = cor(c(vs[[x]]$rep1), c(vs[[y]]$rep2 ))
				randco = rbind(randco, cbind(c, x, y))

			}
		}
	}
	
	fmat = NULL
	sample_cols = NULL
	mens = NULL
	sample_names = NULL
	for(x in 1:length(vs)) {
		fmat = cbind(fmat, cbind(c(vs[[x]]$rep1), c(vs[[x]]$rep2)))
		sample_cols = c(sample_cols, c(x, x))
		m1 = mean(max(c(vs[[x]]$rep1), max(c(vs[[x]]$rep2))))
		mens = c(mens, c(m1, m1))
		s1 = paste(x, "_1", sep="")
		s2 = paste(x, "_2", sep="")
		sample_names = c(sample_names, c(s1, s2))
	}
	heatmap.2(fmat, margins=c(10,10), tracecol=F, col=redgreen(75), ColSideColors=as.character(sample_cols))

	heatmap.2(fmat, margins=c(10,10), tracecol=F, Colv=F, dendrogram="none", col=redgreen(75), ColSideColors=as.character(sample_cols))
	a=  apply(fmat, 2, mean)
	heatmap.2(fmat[,order(a)], margins=c(10,10), tracecol=F, Colv=F, dendrogram="none", col=redgreen(75), ColSideColors=as.character(sample_cols[order(a)]))

	h = heatmap.2(t(fmat), margins=c(10,10), tracecol=F, col=redgreen(75), RowSideColors=as.character(sample_cols))
	ordered_names = labels(h$rowDendrogram)
	sample_key = substr(ordered_names, 1, nchar(ordered_names)-2)



return(c)
}


ewan_plots <- function(files) {
  plot.reg <- function(e, name="") {
    image(t(e$centered_band-median(e$centered_band)), col=topo.colors(10), main=name)
    lines(e$e1/dim(e$centered_band)[1], lwd=3, col="black")
    lines(e$e2/dim(e$centered_band)[1], lwd=3, col="black")
    matplot(e$fitted_center$max_slope[,1]/dim(e$centered_band), e$fitted_center$max_slope[,2]/dim(e$centered_band), col="red", add=T, pch="+", cex=0.5)
    matplot(e$fitted_center$plateu[,1]/dim(e$centered_band), e$fitted_center$plateu[,2]/dim(e$centered_band), col="blue", add=T, pch="+", cex=0.5)
    legend("bottomleft", c("Max slope", "Plateau"), col=c("red", "blue"), lty=1)
  }

  e1 = registar(files[1], debug=FALSE)
  e2 = registar(files[2], debug=FALSE)

  b1 = fit_four_bands(files[1], debug=FALSE)
  b2 = fit_four_bands(files[2], debug=FALSE)

  par(mfrow=c(2,2))
  plot.reg(e1, name="Top surface model - replicate 1")
  plot.reg(e2, name="Top surface model - replicate 2")


  eye(b1$image_data, b1$zpos, interpolate_outliers(b1$outer_band_positions[,c(1,2,3,5)], 1, 31), "Max dip z slice - replicate 1")
  legend("bottomleft", c("Fitted bands"), col="red", lty=1)

  abline(v=min(e1$e1[,1])/dim(e1$centered_band), col="red", lty="dashed")
  abline(v=max(e1$e1[,1])/dim(e1$centered_band), col="red", lty="dashed")
  abline(v=min(e1$e2[,1])/dim(e1$centered_band), col="blue", lty="dashed")
  abline(v=max(e1$e2[,1])/dim(e1$centered_band), col="blue", lty="dashed")


  b3 = b2
  sfac = median(b3$outer_band_positions[,1]) -  median(b1$outer_band_positions[,1])
  b3$image_data$img[,,b3$zpos]
  ei = (nrow(b3$image_data$img[,,b3$zpos]))
  temp1 =  b3$image_data$img[1:sfac,,b3$zpos]
  temp2 =  b3$image_data$img[(sfac+1):ei,,b3$zpos]
  t = rbind(temp2, temp1)
  b3$image_data$img[,,b3$zpos] = t

  eye(b3$image_data, b3$zpos, interpolate_outliers(b3$outer_band_positions[,c(1,2,3,5)]-sfac, 1, 31), "Max dip z slice - replicate 2")
  legend("bottomleft", c("Fitted bands"), col="red", lty=1)
  abline(v=min(e2$e1[,1])/dim(e2$centered_band), col="red", lty="dashed")
  abline(v=max(e2$e1[,1])/dim(e2$centered_band), col="red", lty="dashed")
  abline(v=min(e2$e2[,1])/dim(e2$centered_band), col="blue", lty="dashed")
  abline(v=max(e2$e2[,1])/dim(e2$centered_band), col="blue", lty="dashed")

}


