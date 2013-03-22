# -------------------------- LIST OF FUNCTIONS
# LL.to.ARBUTM()
# shape.prj.coords()
# log1()
# negf()
# invf()
# mginv()
# sizecolorgr()
# empirical.semivariogram()
# graph.semivariogram.trellis()
# emp.semivar.stream()
# emp.variogram.flow.connect()
# block.krige()
# graph.data.resid()
# graph.pred.data.krige()
# R2g()
# crossval()
# infldiagn()
# X.design()
# stream.graph.sizecolor()
# stream.graph.sizebw()
# graph.scatter()
# graph.histogram()
# graph.boxplot()
# graph.QQplot()
# fixed.effect.table()
# estimates.and.contrasts()
# hypoth.stream()
# pred.k()
# dist.geo.ani()
# dist.geo.iso()
# exp.Euclid()
# sph.Euclid()
# gau.Euclid()
# cau.Euclid()
# exp.tailup()
# sph.tailup()
# lws.tailup()
# mar.tailup()
# exp.taildown()
# sph.taildown()
# lws.taildown()
# mar.taildown()
# theta.ini()
# m2LL.stream()
# arg.error.check()
# make.cov.mat()
# untrans.parm()
# Info.crit.compare()
# stdX()
# unstdX.beta.hat()
# slm.stream()


# ------------------------------------------------------------------------------

LL.to.ARBUTM <-
function(cm,lat,lon)
# This function converts from Lat-Lon (decimal degrees) to the Universal
# Transverse Mercator Coordinates and returns the new coordinates in a
# 2-column matrix with x- and y- as columns.  In this program, the
# coordinates are calculated from a user supplied central meridian.
# Coordinates are returned in kilometers from the western-most
# longitude and the southern-most latitude observed in the data set.

{
# initialize some variables
	e2 <- 0.00676865799729
	a <- 6378206.4
	ep2 <- e2 / (1-e2)
	drc <- pi / 180
	sc <- 0.9996
	fe <- 500000
	ftm <- 0.30480371
#calculate some frequently used values
	lar <- lat * drc
	ls <- sin(lar)
	ls2 <- ls^2
	els2 <- ep2 * ls2
	lc <- cos(lar)
	lc2 <- lc^2
	lc3 <- lc^3
	lc5 <- lc^5
	elc2 <- ep2 * lc2
	lt2 <- tan(lar)^2
	lt4 <- lt2^2
# do the transformation
	v <- a/sqrt(1 - e2*ls2)
	p <- drc*(cm - lon)
	temp <- 5104.57388 - (lc2*(21.73607 - 0.11422*lc2))
	r1 <- 6367399.689*(lar - ls*lc*0.000001*temp)
	r2 <- (v*ls*lc*p^2)/2
	temp <- 5 - lt2 + 9*elc2 + (2*elc2)^2
	r3 <- (v*ls*lc3*p^4*temp)/24
	r4 <- v*lc*p
	temp <- 1 - lt2 + elc2
	r5 <- (v*lc3*p^3*temp)/6
	temp <- 61 - 58*lt2 + lt4 + 270*elc2 - 330*els2
	ra6 <- (v*ls*lc5*p^6*temp)/720
	temp <- 5 - 18*lt2 + lt4 + 14*elc2 - 58*els2
	rb5 <- (v*lc5*p^5*temp)/120
	northing <- sc*(r1 + r2 + r3 + ra6)
	easting <- -sc*(r4 + r5 + rb5)
	y<-(northing-min(northing))/1000
	x<-(easting-min(easting))/1000
	cbind(x,y)
}

#------------------- GET COORDINATES FROM ESRI SHAPE FILE

shape.prj.coords <- function(rshape){
	tmp <- matrix(0, nrow = length(rshape$att.data[,1]), ncol = 2)
	for (i in 1:length(rshape$att.data[,1])){
		tmp[i,] <- rshape$Shape[[i]]$verts
	}
	data.frame(Centroid.lon = tmp[,1], Centroid.lat = tmp[,2])
}


#------------------- SOME SIMPLE FUNCTIONS, USED WITH SIZE-COLOR GRAPHS

log1 <- function(x){log(x+1)}
negf <- function(x){-x}
invf <- function(x){1/x}


#------------------- GENERALIZED INVERSE OF A MATRIX

mginv <- function(X, tol = sqrt(.Machine$double.eps)) { 
	dnx <- dimnames(X) 
	if(is.null(dnx)) dnx <- vector("list", 2) 
	s <- svd(X) 
	nz <- s$d > tol * s$d[1] 
	structure( 
		if(any(nz)) s$v[, nz] %*% (t(s$u[, nz])/s$d[nz]) else X, 
		dimnames = dnx[2:1]) 
} 

#------------------- BASIC SIZE-COLOR GRAPH

size.color.gr <- function(data, 
	xcol = "x", ycol = "y", 
	sample.size.col, response.col, 
	sample.size.fun = eval, response.fun = eval,
	cex.start.sample.size = 1, cex.increment.sample.size = 0, 
	nclass.response = 10, nclass.sample = 10, 
	color.palette = NULL, breaks = NULL, ...)
{
	if(is.null(breaks)==TRUE) {
		lower.breaks <- matrix(0, nrow = nclass.response, ncol = 1)
		upper.breaks <- matrix(0, nrow = nclass.response, ncol = 1)
	} else {
		lower.breaks <- breaks[,1]
		upper.breaks <- breaks[,2]
	}
	if(length(color.palette) == 0)
		color.palette <- rainbow(nclass.response, start = .66, end = .99)
	col.stepi <- (max(sample.size.fun(data[,sample.size.col]),na.rm = TRUE) - 
		min(sample.size.fun(data[,sample.size.col]),na.rm = TRUE))/nclass.sample + 
		0.0001/(nclass.sample - 1)
	imax <- min(sample.size.fun(data[,sample.size.col]),na.rm = TRUE) - 0.0001
	col.stepj <- (max(response.fun(data[,response.col]),na.rm = TRUE) - 
		min(response.fun(data[,response.col]),na.rm = TRUE))/nclass.response + 
		0.0001/(nclass.response - 1)
	jmax1 <- min(response.fun(data[,response.col]),na.rm = TRUE) - 0.0001
	for (i in 1:nclass.sample){
		imin <- imax
		imax <- imax + col.stepi
		indi <- sample.size.fun(data[,sample.size.col]) >= imin & 
			sample.size.fun(data[,sample.size.col]) <= imax
		for (j in 1:nclass.response){
			if(j == 1) jmax <- jmax1
			jmin <- jmax
			if(is.null(breaks)) lower.breaks[j] <- jmin
			jmax <- jmax + col.stepj
			if(is.null(breaks)) upper.breaks[j] <- jmax
			indj <- response.fun(data[,response.col]) >= lower.breaks[j] & 
				response.fun(data[,response.col]) <= upper.breaks[j]
			points(data[indi & indj,xcol],data[indi & indj,ycol],
				col = color.palette[j], 
				pch = 19, 
				cex = cex.start.sample.size + (i-1)*cex.increment.sample.size)
		}
	}
	data.frame(lower.breaks = lower.breaks, upper.breaks = upper.breaks)
}

#--------- USUAL EMPIRICAL SEMIVARIOGRAM FOR EUCLIDEAN DISTANCE

empirical.semivariogram <-
function(data, x, y, var,
	nlag = 20, directions = c(0,45,90,135),
	tolerance = 22.5, inc = 0, maxlag = 1e32, nlagcutoff = 1,
	EmpVarMeth = "MethMoment")
# EMPIRICAL SEMIVARIOGRAM FUNCTION
# var1 is a matrix or data frame with x-coord in the first column
#                                     y-coord in the second column
#                                     z (response) in the third column
{
	n1 <- length(data[,1])
   # distance matrix among locations
	distance <- sqrt( ( matrix(data[,x],nrow=n1,ncol=1) %*%
		matrix(rep(1,times=n1),nrow=1,ncol=n1) -
		matrix(rep(1,times=n1),nrow=n1,ncol=1) %*%
		matrix(data[,x],nrow=1,ncol=n1) )^2 +
		( matrix(data[,y],nrow=n1,ncol=1) %*%
		matrix(rep(1,times=n1),nrow=1,ncol=n1) -
		matrix(rep(1,times=n1),nrow=n1,ncol=1) %*%
		matrix(data[,y],nrow=1,ncol=n1) )^2 )
	difx <- -(matrix(data[,y],nrow=n1,ncol=1) %*%
		matrix(rep(1,times=n1),nrow=1,ncol=n1) -
		matrix(rep(1,times=n1),nrow=n1,ncol=1) %*%
		matrix(data[,y],nrow=1,ncol=n1))
	signind <- -(matrix(data[,x],nrow=n1,ncol=1) %*%
		matrix(rep(1,times=n1),nrow=1,ncol=n1) -
		matrix(rep(1,times=n1),nrow=n1,ncol=1) %*%
		matrix(data[,x],nrow=1,ncol=n1)) < 0
	distance <- distance*1.0000000001
	theta.deg<-acos(difx/distance)*180/pi
   # matrix of degrees clockwise from north between locations
	theta.deg[signind] <- 360-theta.deg[signind]
	diff2 <- ( matrix(data[,var],nrow=n1,ncol=1) %*%
		matrix(rep(1,times=n1),nrow=1,ncol=n1) -
		matrix(rep(1,times=n1),nrow=n1,ncol=1) %*%
		matrix(data[,var],nrow=1,ncol=n1) )^2
	sqrtdiff <- sqrt(abs( matrix(data[,var],nrow=n1,ncol=1) %*%
		matrix(rep(1,times=n1),nrow=1,ncol=n1) -
		matrix(rep(1,times=n1),nrow=n1,ncol=1) %*%
		matrix(data[,var],nrow=1,ncol=n1) ) )
	if(EmpVarMeth == "CovMean") temp4cov <- data[,var] - mean(data[,var])
	else temp4cov <- data[,var]
	covprod <- (matrix(temp4cov,nrow=n1,ncol=1) %*%
		matrix(rep(1,times=n1),nrow=1,ncol=n1)) *
		(matrix(rep(1,times=n1),nrow=n1,ncol=1) %*%
		matrix(temp4cov,ncol=n1,nrow=1))
# convert to vectors
	distance <- matrix(distance, ncol = 1)
	theta.deg <- matrix(theta.deg, ncol = 1)
	diff2 <- matrix(diff2, ncol = 1)
	sqrtdiff <- matrix(sqrtdiff, ncol = 1)
	covprod <- matrix(covprod, ncol = 1)
# trim off values greater than maxlag
	indmax <- distance <= maxlag
	distance <- distance[indmax,]
	theta.deg <- theta.deg[indmax,]
	diff2 <- diff2[indmax,]
	sqrtdiff <- sqrtdiff[indmax,]
	covprod <- covprod[indmax,]

	maxd<-max(distance)
	if( inc <= 0) inc <- maxd/nlag
	ind <- distance==0
	ndir <- length(directions)
	store.results <- matrix(data = NA, ncol = 6,
		dimnames = list(NULL, c("distance", "gamma", "np", "azimuth", "hx", "hy")))
	for (j in 1:ndir) {
		for ( i in 1:nlag){
			if( (directions[j]-tolerance)<0 && (directions[j]+tolerance)>0 )
				ind1 <- theta.deg >= 360+directions[j]-tolerance |
					theta.deg < directions[j]+tolerance
			else if( (directions[j]+tolerance)>360 && (directions[j]-tolerance)<360 )
				ind1 <- theta.deg < directions[j]+tolerance-360 |
					theta.deg >= directions[j]-tolerance
			else
				ind1 <- theta.deg >= directions[j]-tolerance &
					theta.deg < directions[j]+tolerance
			ind<-distance>(i-1)*inc & distance<=i*inc &
				!is.na(theta.deg) & ind1
			nclass <- sum(ind)
			if(EmpVarMeth == "MethMoment") cv <- mean(diff2[ind])
			if(EmpVarMeth == "RobustMean") cv <- ((mean(sqrtdiff[ind]))^4)/(.457 + .494/sum(ind))
			if(EmpVarMeth == "RobustMedian") cv <- (median(sqrtdiff[ind]))^4/.457
			if(EmpVarMeth == "Covariance" | EmpVarMeth == "CovMean") cv <- mean(covprod[ind])
			mean.dis <- mean(distance[ind])
			if(nclass > 0) store.results<-rbind(store.results,
				c(mean.dis,cv,nclass,directions[j],0,0))
		}
	}
	store.results[,"hx"]<-store.results[,"distance"]*sin(store.results[,"azimuth"]*pi/180)
	store.results[,"hy"]<-store.results[,"distance"]*cos(store.results[,"azimuth"]*pi/180)
	store.results[,"gamma"]<-store.results[,"gamma"]/2
	ind <- store.results[,"np"] >= nlagcutoff
	store.results <- store.results[ind,]
	ind <- !is.na(store.results[,"hx"])
	store.results <- store.results[ind,]
	as.data.frame(store.results)
}

# ---------- TRELLIS GRAPH OF (DIRECTIONAL) EMPIRICAL SEMIVARIOGRAM VALUES

graph.semivariogram.trellis <- function(data, response.column, xcol, ycol, nlag = 20,
	directions = c(0, 45, 90, 135), tolerance = 22.5, inc = 0, maxlag = 1e32,
	nlagcutoff = 1, EmpVarMeth = "MethMoment", cex = 1)
{
  data <- data[!is.na(data[,response.column]),]
	emp.var.out <- empirical.semivariogram(data,
		x = xcol, y = ycol, var = response.column,
		nlag = nlag,
		directions = directions,
		tolerance = tolerance, inc = inc, nlagcutoff = nlagcutoff,
		maxlag = maxlag, EmpVarMeth = EmpVarMeth)
	win.graph()
		plot(c(1,1),type = "n")
		t.background <- trellis.par.get("background")
		t.background$col <- "white"
		trellis.par.set("background", t.background)
		xyplot(gamma ~ distance | azimuth, data=emp.var.out, pch = 19, col = "black",
			cex = cex, as.table = TRUE)
}

#--------- EMPIRICAL SEMIVARIOGRAM FOR STREAM DISTANCE

emp.semivar.stream <-
function(z, flow, distance,
	nlag = 20, inc = 0, maxlag = 1e32, nlagcutoff = 1,
	EmpVarMeth = "MethMoment")
# EMPIRICAL SEMIVARIOGRAM FUNCTION
# var1 is a matrix or data frame with x-coord in the first column
#                                     y-coord in the second column
#                                     z (response) in the third column
{
	n1 <- length(z)
	diff2 <- ( matrix(z,nrow=n1,ncol=1) %*% 
		matrix(rep(1,times=n1),nrow=1,ncol=n1) -
		matrix(rep(1,times=n1),nrow=n1,ncol=1) %*%
		matrix(z,nrow=1,ncol=n1) )^2 
	sqrtdiff <- sqrt(abs( matrix(z,nrow=n1,ncol=1) %*% 
		matrix(rep(1,times=n1),nrow=1,ncol=n1) -
		matrix(rep(1,times=n1),nrow=n1,ncol=1) %*%
		matrix(z,nrow=1,ncol=n1) ) )
	if(EmpVarMeth == "CovMean") temp4cov <- z - mean(z)
	else temp4cov <- z
	covprod <- (matrix(temp4cov,nrow=n1,ncol=1) %*% 
		matrix(rep(1,times=n1),nrow=1,ncol=n1)) *
		(matrix(rep(1,times=n1),nrow=n1,ncol=1) %*%
		matrix(temp4cov,ncol=n1,nrow=1))
# convert to vectors
	distance <- matrix(distance, ncol = 1)
	flow <- matrix(flow, ncol = 1)
	diff2 <- matrix(diff2, ncol = 1)
	sqrtdiff <- matrix(sqrtdiff, ncol = 1)
	covprod <- matrix(covprod, ncol = 1)
# trim off values greater than maxlag
	indmax <- distance <= maxlag
	distance <- distance[indmax,]
	flow <- flow[indmax,]
	diff2 <- diff2[indmax,]
	sqrtdiff <- sqrtdiff[indmax,]
	covprod <- covprod[indmax,]

	maxd<-max(distance)
	if( inc <= 0) inc <- maxd/nlag
	ind <- distance==0
	store.results <- matrix(data = NA, ncol = 6,
		dimnames = list(NULL, c("distance.connect", "gam.connect", "np.connect", 
		"distance.unconnect", "gam.unconnect","np.unconnect")))
	for ( i in 1:nlag){
		ind1 <- distance > (i-1)*inc & distance <= i*inc & flow
		ind2 <- distance > (i-1)*inc & distance <= i*inc & !flow
		nclass1 <- sum(ind1)
		nclass2 <- sum(ind2)
		if(EmpVarMeth == "MethMoment") {
			cv1 <- mean(diff2[ind1])
			cv2 <- mean(diff2[ind2])
		}
		if(EmpVarMeth == "RobustMean") {
			cv1 <- ((mean(sqrtdiff[ind1]))^4)/(.457 + .494/sum(ind1))
			cv2 <- ((mean(sqrtdiff[ind2]))^4)/(.457 + .494/sum(ind2))
		}
		if(EmpVarMeth == "RobustMedian") {
			cv1 <- (median(sqrtdiff[ind1]))^4/.457
			cv2 <- (median(sqrtdiff[ind2]))^4/.457
		}
		if(EmpVarMeth == "Covariance" | EmpVarMeth == "CovMean") {
			cv1 <- mean(covprod[ind1])
			cv2 <- mean(covprod[ind2])
		}
		mean.dis1 <- mean(distance[ind1])
		mean.dis2 <- mean(distance[ind2])
		if(nclass1 > 0 | nclass2 > 0) store.results <- rbind(store.results,
				c(mean.dis1, cv1, nclass1, mean.dis2, cv2, nclass2))
	}
	store.results[,"gam.connect"] <- store.results[,"gam.connect"]/2
	store.results[,"gam.unconnect"] <- store.results[,"gam.unconnect"]/2
	ind <- store.results[,"np.connect"] >= nlagcutoff |
		store.results[,"np.unconnect"] >= nlagcutoff
	store.results <- store.results[ind,]
	ind <- !is.na(store.results[,"np.connect"])
	store.results <- store.results[ind,]
	as.data.frame(store.results)
}

#---------- GRAPH EMPIRICAL SEMIVARIOGRAM FOR STREAM NETWORKS -------

emp.variogram.flow.connect <- function(data, response.col,
  flow.matrix, dist.junc,
  maxlag, nlag = 6, inc = 0, nlag.cutoff = 15,
  lower.data.trim = -1e50, upper.data.trim = 1e50,
  transformation = "none", min.cex = 1.5, max.cex = 6)
{
  zdata <- data[,response.col]
  ind <- !is.na(zdata) & zdata>lower.data.trim & zdata<upper.data.trim
  zdata <- zdata[ind]
  if(transformation == "sqrt") zdata <- sqrt(zdata)
  if(transformation == "log") zdata <- log(zdata + 1)

  flow.matrix <- flow.matrix + t(flow.matrix)
  diag(flow.matrix) <- 1
  stream.distance <- dist.junc + t(dist.junc)

  flow.matrix.obs <- flow.matrix[ind, ind]
  flow.matrix.01.obs <- flow.matrix.obs > 0
  stream.distance.obs <- stream.distance[ind, ind]

  # source("d:\\StreamNetworks\\slm_functions050310.txt")
  #  debug(empirical.semivariogram)
  ev <- emp.semivar.stream(zdata, flow.matrix.01.obs,
	  stream.distance.obs,
	  nlag = nlag, inc = inc, maxlag = maxlag, nlagcutoff = 1,
	  EmpVarMeth = "MethMoment")
	  ev[is.nan(ev[,"distance.unconnect"]),"distance.unconnect"] <- 0
	  ev[is.nan(ev[,"gam.unconnect"]),"gam.unconnect"] <- 0
	  ev[is.nan(ev[,"distance.connect"]),"distance.connect"] <- 0
	  ev[is.nan(ev[,"gam.connect"]),"gam.connect"] <- 0

  plot(c(0,
	  max(ev[,"distance.connect"],ev[,"distance.unconnect"])),
	  c(0,
	  max(ev[,"gam.connect"],ev[,"gam.unconnect"])), type = "n",
	  xlab = "Stream Distance",
	  ylab = "Semivariogram",
  	cex.lab = 1.5,
	  cex.axis = 1
  )
  nlag <- length(ev[,1])
  maxnp <- max(ev[,"np.connect"],ev[,"np.unconnect"])
  minnp <- min(ev[,"np.connect"],ev[,"np.unconnect"])
  np.range <- maxnp - minnp
  cex.inc <- max.cex - min.cex
  for(i in 1:nlag) {
	  if(ev[i,"np.connect"] > nlag.cutoff)
	  points(ev[i,"distance.connect"],
		  ev[i,"gam.connect"], pch = 19, col = "blue",
		  cex = min.cex + cex.inc*ev[i,"np.connect"]/maxnp)
	  if(ev[i,"np.unconnect"] > nlag.cutoff)
	  points(ev[i,"distance.unconnect"],
		  ev[i,"gam.unconnect"], pch = 19, col = "green",
		  cex = min.cex + cex.inc*ev[i,"np.unconnect"]/maxnp)
  }

}

#------------ MAKE MAP OF VARIOUS RESIDUALS

graph.data.resid <- function(rshape.outline, data, xcol, ycol, 
	sample.size.col, response.col, outfilen,
	sample.size.fun = eval, response.fun = eval, 
	nclass.response = 10, nclass.sample = 10,
	cex.start.sample.size = .5, cex.increment.sample.size = 0.3,
	leg.deci.dig = 2,
	color.palette = c("darkblue", 
		"blue",
		"cyan3", 
		"cyan",
		"lightgreen",
		"greenyellow",
		"yellow",
		"orange",
		"tomato",
		"red") )
{
	postscript(outfilen, onefile = TRUE, horizontal = FALSE)
		layout(matrix(c(1,2), ncol = 2), widths = c(2,1))
		par(mar=c(0,0,0,0))
		plot(rshape.outline, xaxt = "n", yaxt = "n", fg="white", xlab = "", ylab = "")
		breaks.j <- size.color.gr(data = data, 
			xcol = xcol, ycol = ycol, 
			sample.size.col = sample.size.col, response.col = response.col, 
			sample.size.fun = invf, response.fun = eval,
			nclass.response = nclass.response, nclass.sample = nclass.sample,
			cex.start.sample.size = cex.start.sample.size,
			cex.increment.sample.size = cex.increment.sample.size,
			color.palette = color.palette)
		plot(data.frame(x = 0, y = 0), type = "n", 
			xaxt = "n", yaxt = "n", bty = "n", xlab = "", ylab = "")
		nleg <- length(breaks.j[,1])
		dec.dig <- leg.deci.dig
		lowerb.cha <- as.character(as.numeric(as.integer(breaks.j[,"lower.breaks"]*10^dec.dig))/10^dec.dig)
		dash.j <- rep("to",times = nleg)
		upperb.cha <- as.character(as.numeric(as.integer(breaks.j[,"upper.breaks"]*10^dec.dig))/10^dec.dig)
		leg.lab <- paste(lowerb.cha,dash.j, upperb.cha)
		legend(x = -1, y = .25, pch = 19, cex = 1.1,
			legend = leg.lab,
			col = color.palette)
	if(graphic.device == "postscript") dev.off(dev.cur())
}

#------------ MAKE MAP OF KRIGING PREDICTIONS AND STANDARD ERRORS

graph.pred.data.krige <- function(rshape.outline, data, data.pred, xcol, ycol,
	sample.size.col, response.col, graphic.device = "windows", output.filename, 
	sample.size.fun = eval, response.fun = eval, 
	nclass.response = 10, nclass.sample = 10,
	cex.start.sample.size = .5, cex.increment.sample.size = 0.3,
	leg.deci.dig = 2,
	color.palette = c("darkblue", 
		"blue",
		"cyan3", 
		"cyan",
		"lightgreen",
		"greenyellow",
		"yellow",
		"orange",
		"tomato",
		"red"), ...
)
{
	if(graphic.device == "postscript") postscript(output.filename, onefile = TRUE, horizontal = FALSE)
	if(graphic.device == "windows") win.graph()
		layout(matrix(c(1,2), ncol = 2), widths = c(2,1))
		par(mar=c(0,0,0,0))
		plot(rshape.outline, xaxt = "n", yaxt = "n", fg="white", xlab = "", ylab = "")
		breaks.j <- size.color.gr(data = data.pred, 
			xcol = xcol, ycol = ycol, 
			sample.size.col = sample.size.col, response.col = response.col, 
			sample.size.fun = sample.size.fun, response.fun = response.fun,
			nclass.response = nclass.response, nclass.sample = nclass.sample,
			cex.start.sample.size = cex.start.sample.size,
			cex.increment.sample.size = cex.increment.sample.size,
			color.palette = color.palette)
		points(data[,xcol], data[,ycol], pch = 19, ...)
		plot(data.frame(x = 0, y = 0), type = "n", 
			xaxt = "n", yaxt = "n", bty = "n", xlab = "", ylab = "")
		nleg <- length(breaks.j[,1])
		dec.dig <- leg.deci.dig
		lowerb.cha <- as.character(as.numeric(as.integer(breaks.j[,"lower.breaks"]*10^dec.dig))/10^dec.dig)
		dash.j <- rep("to",times = nleg)
		upperb.cha <- as.character(as.numeric(as.integer(breaks.j[,"upper.breaks"]*10^dec.dig))/10^dec.dig)
		leg.lab <- paste(lowerb.cha,dash.j, upperb.cha)
		legend(x = -1, y = .25, pch = 19, cex = 1.1,
			legend = leg.lab,
			col = color.palette)
	if(graphic.device == "postscript") dev.off(dev.cur())
}

#------------ GENERALIZED R2 FOR SPATIAL MODELS

R2g <- function(z, Xdesign, bhat, Vi)
{
	SS.bhat <- t(z - Xdesign %*% bhat) %*% Vi %*% (z - Xdesign %*% bhat)
	mu.hat <- sum(Vi %*% z)/sum(Vi)
	SS.mu <- matrix(z - mu.hat, nrow = 1) %*% Vi %*% matrix(z - mu.hat, ncol = 1)
	1 - SS.bhat/SS.mu
}

#------------ CROSSVALIDATION DIAGNOSTICS

crossval <- function(z, X, V, Vi, n, p)
{
	cdd.out <- matrix(-999.9, nrow = n, ncol = 2)
	for(i in 1:n) {
		Vi.i <- Vi[(1:n) != i,(1:n) != i] - 
			matrix(Vi[(1:n) != i,i],ncol = 1) %*% 
			matrix(Vi[i,(1:n) != i],nrow = 1)/Vi[i,i]
		c.i <- matrix(V[(1:n) != i,i],ncol = 1)
		xi <- matrix(X[i,], ncol = 1)
		X.i <- X[(1:n) != i,]
		z.i <- matrix(z[(1:n) != i], ncol = 1)
		xxi <- xi - t(X.i) %*% Vi.i %*% c.i
		covb.i <- mginv(t(X.i) %*% Vi.i %*% X.i)
		si <- V[i,i]  - t(c.i) %*% Vi.i %*% c.i
		lam <- t(c.i + X.i %*% covb.i %*% xxi) %*% Vi.i
		cdd.out[i,1] <- lam %*% z.i
		cdd.out[i,2] <- sqrt(si + t(xxi) %*% covb.i %*% xxi)
	}
	cdd.out <- as.data.frame(cdd.out)
	names(cdd.out) <- c("cv.pred","cv.se")
	cdd.out
}


#------------ INFLUENCE DIAGNOSTICS

infldiagn <- function(z, V, X, n, p) 
{
	svd.out <- svd(V)
	V.5 <- svd.out$u %*% diag(1/sqrt(svd.out$d)) %*% t(svd.out$u)
	z.5 <- V.5 %*% z
	X.5 <- V.5 %*% X

	# Hat matrix, Montgomery and Peck, pg. 170
	Hat <- X.5 %*% mginv(t(X.5) %*% X.5) %*% t(X.5)
	# standardized residuals, Montgomery and Peck, pg. 170 
	e <- (diag(n) - Hat) %*% z.5
	# studentized residuals Montgomery and Peck, pg. 172
	r <- e/sqrt(diag((diag(n) - Hat)))
	# leverage - diagonal hat elements scaled by 2p/n, Montgomery and Peck, pg. 182
	hii <- diag(Hat)
	lever <- hii/(2*p/n)
	# Cook's D, Montgomery and Peck, pg. 182
	Di <- r^2*hii/((1-hii)*p)
	mat <- as.data.frame(cbind(e,r,lever,Di))
	colnames(mat) <- c("resid.stand", "resid.student", "leverage", "CooksD")
	mat
}

# ------------ CREATE ONE OF SEVERAL TYPES OF DESIGN MATRICES

X.design <- function(formula, data, Xmethod = "treatments")
{

	if(Xmethod == "treat.first.0") {
		options(contrasts=c("contr.treatment", "contr.treatment"))
		X <- model.matrix(formula, data = data)
		options(contrasts=c("contr.treatment", "contr.poly"))
	}

	if(Xmethod == "sum.to.0") {
		options(contrasts=c("contr.sum", "contr.sum"))
		X <- model.matrix(formula, data = data)
		options(contrasts=c("contr.treatment", "contr.poly"))
	}
	if(Xmethod == "over.par") {
		terms.list <- attr(terms(formula),"term.labels")
		X <- NULL
		if(attr(terms(formula),"intercept") == 1)
			X <- rep(1, times = length(data[,1]))
		for (i in 1:length(terms.list)) {
			form1 <- formula(paste("~ ", terms.list[i], " - 1"))
			X <- cbind(X,model.matrix(form1, data = data))
		}
	}
	X
}


#------------ CREATE A MAP WITH LOCATIONS COLORED ACCORDING TO THEIR VALUE

stream.graph.sizecolor <- function(obs.data, stream.shape = NULL, pred.loc = NULL, 
		xcol = xcol, ycol = ycol, 
		color.col, color.fun = eval, size.col, size.fun = eval, 
		cex.start.size = 1.5, cex.increment.size = .001, breaks = NULL,
		dec.dig, graphic.device = "windows", output.filename)
{
	if(graphic.device == "postscript") postscript(output.filename)
	if(graphic.device == "windows") win.graph()
		layout(matrix(c(1,2), ncol = 2), widths = c(2,1))
		par(mar=c(0,0,0,0))
		if(!is.null(stream.shape))
			plot(stream.shape, xaxt = "n", 
				yaxt = "n", fg="white", 
				xlab = "", ylab = "")
		else plot(obs.data[,c(xcol,ycol)], type = "n")
		breaks.j <- size.color.gr(data = obs.data, 
			xcol = xcol, ycol = ycol, 
			sample.size.col = size.col, response.col = color.col, 
			sample.size.fun = size.fun, response.fun = color.fun,
			nclass.response = 10, nclass.sample = 10,
			cex.start.sample.size = cex.start.size, 
			cex.increment.sample.size = cex.increment.size,
			breaks = breaks,
			color.palette = c("darkblue", 
				"blue",
				"cyan3", 
				"cyan",
				"lightgreen",
				"greenyellow",
				"yellow",
				"orange",
				"tomato",
				"red"))
		if(!is.null(pred.loc)) points(pred.loc, pch = 19, cex = 1, col = "black")
		plot(data.frame(x = 0, y = 0), type = "n", 
			xaxt = "n", yaxt = "n", bty = "n", xlab = "", ylab = "")
		nleg <- length(breaks.j[,1])
		dec.dig <- dec.dig
		lowerb.cha <- as.character(as.numeric(as.integer(
			breaks.j[,"lower.breaks"]*10^dec.dig))/10^dec.dig)
		dash.j <- rep("to",times = nleg)
		upperb.cha <- as.character(as.numeric(as.integer(breaks.j[,"upper.breaks"]*10^dec.dig))/10^dec.dig)
		leg.lab <- paste(lowerb.cha, dash.j, upperb.cha)
		legend(x = -1, y = .5, pch = 19, cex = 1.1,
			legend = leg.lab,
			col = c("darkblue", 
				"blue",
				"cyan3", 
				"cyan",
				"lightgreen",
				"greenyellow",
				"yellow",
				"orange",
				"tomato",
				"red"))
	if(graphic.device == "postscript") dev.off(dev.cur())
}

#------------ CREATE A MAP WITH LOCATIONS COLORED ACCORDING TO THEIR VALUE

stream.graph.sizebw <- function(obs.data, stream.shape = NULL, pred.loc = NULL, 
		xcol = xcol, ycol = ycol, 
		color.col, color.fun = eval, size.col, size.fun = eval, 
		cex.start.size = 1.5, cex.increment.size = .001, breaks = NULL,
		dec.dig, graphic.device = "windows", output.filename)
{
	if(graphic.device == "postscript") postscript(output.filename)
	if(graphic.device == "windows") win.graph()
		layout(matrix(c(1,2), ncol = 2), widths = c(2,1))
		par(mar=c(0,0,0,0))
		if(!is.null(stream.shape))
			plot(stream.shape, xaxt = "n", 
				yaxt = "n", fg="white", 
				xlab = "", ylab = "", lwd = 3)
		else plot(obs.data[,c(xcol,ycol)], type = "n", lwd = 2)
		breaks.j <- size.color.gr(data = obs.data, 
			xcol = xcol, ycol = ycol, 
			sample.size.col = size.col, response.col = color.col, 
			sample.size.fun = size.fun, response.fun = color.fun,
			nclass.response = 10, nclass.sample = 10,
			cex.start.sample.size = cex.start.size, 
			cex.increment.sample.size = cex.increment.size,
			breaks = breaks,
			color.palette = c("grey5", 
				"grey15",
				"grey25", 
				"grey35",
				"grey45",
				"grey55",
				"grey65",
				"grey75",
				"grey85",
				"grey95"))
		points(obs.data[,c(xcol,ycol)], cex = cex.start.size)		
		if(!is.null(pred.loc)) points(pred.loc, pch = 19, cex = 1, col = "black")
		plot(data.frame(x = 0, y = 0), type = "n", 
			xaxt = "n", yaxt = "n", bty = "n", xlab = "", ylab = "")
		nleg <- length(breaks.j[,1])
		dec.dig <- dec.dig
		lowerb.cha <- as.character(as.numeric(as.integer(
			breaks.j[,"lower.breaks"]*10^dec.dig))/10^dec.dig)
		dash.j <- rep("to",times = nleg)
		upperb.cha <- as.character(as.numeric(as.integer(breaks.j[,"upper.breaks"]*10^dec.dig))/10^dec.dig)
		leg.lab <- paste(lowerb.cha, dash.j, upperb.cha)
		legend(x = -1, y = .5, pch = 19, cex = 1.1,
			legend = leg.lab,
			col = c("darkblue", 
				"blue",
				"cyan3", 
				"cyan",
				"lightgreen",
				"greenyellow",
				"yellow",
				"orange",
				"tomato",
				"red"))
	if(graphic.device == "postscript") dev.off(dev.cur())
}


#------------ CREATE A (MULTIVARIATE) SCATTERPLOT

graph.scatter <- function(data, output.filename, graphic.device, ...)
{
	if(graphic.device == "postscript") postscript(output.filename)
	if(graphic.device == "windows") win.graph()
		plot(data, pch = 19, main = "Scatter Plots", ...)
	if(graphic.device == "postscript") dev.off(dev.cur())
}

#------------ CREATE A HISTOGRAM

graph.histogram <- function(data, response.column, output.filename, graphic.device)
{
	if(graphic.device == "postscript") postscript(output.filename)
	if(graphic.device == "windows") win.graph()
		hist(data[,response.column],lwd = 3, col = "blue", nclass = 15,
			xlim = c(min(data[,response.column]),max(data[,response.column])), 
			main = "Histogram", xlab = "Value Classes")
	if(graphic.device == "postscript") dev.off(dev.cur())
}

# ---------- CREATE A BOXPLOT 

graph.boxplot <- function(data, response.column, output.filename, graphic.device)
{
	if(graphic.device == "postscript") postscript(output.filename)
	if(graphic.device == "windows") win.graph()
		boxplot(data[,response.column], col = "blue", pch = 19, 
			cex.lab = 1.5, cex = 2, main = "Box Plot")
	if(graphic.device == "postscript") dev.off(dev.cur())
}

# ---------- CREATE A QQPLOT

graph.QQplot <- function(data, response.column, output.filename, graphic.device)
{
	if(graphic.device == "postscript") postscript(output.filename)
	if(graphic.device == "windows") win.graph()
		qqnorm(data[,response.column], pch = 19, cex = 1.2, cex.lab = 1.5, main = "QQ Plot")
	if(graphic.device == "postscript") dev.off(dev.cur())
}


# ---------- FORMAT FIXED EFFECTS TABLE TO INCLUDE FACTOR LEVELS SET TO 0

fixed.effect.table <- function(formula, data, fixed.effects.estimates)
{	
	fee <- fixed.effects.estimates
	form1 <- formula
	fee.i <- 0
	fee.xtra <- NULL
	levels.lab <- NULL
	if(attr(terms(form1),"intercept") == 1) { 
		levels.lab <- data.frame(
			factor = "intercept", 
			levels = "", 
			fee[1,]
			)
		fee.i <- 1
	}
	terms.table <- attr(terms(form1),"factors")
	if(length(attr(terms(form1),"term.labels")) > 0) {
		if(sum(terms.table[,1]) == 1 & attr(terms(form1),"intercept") == 0)
			terms.table[terms.table[,1] == 1,1] <- 2
		for (j in 1:length(terms.table[1,])) {
			ind.j <- terms.table[,j] > 0
			n.interact <- sum(ind.j)
			lab.i <- NULL
			for (i in 1:n.interact) {
				if(!any(rownames(terms.table)[ind.j][i] == colnames(data))) {levels.i <- matrix("")
					if(i == 1) lab.i <- matrix(levels.i, ncol = 1)
					if(i > 1) {
						lab.t <- NULL
						for (k in 2:i) lab.t <- cbind(lab.t,
							matrix(rep(lab.i[,k-1], times = length(levels.i)), ncol = 1))
						lab.tmp <- matrix(rep(levels.i, each = length(lab.i[,1])), ncol = 1)
						lab.i <- cbind(lab.t,lab.tmp)
					}
				}
				else if(is.factor(data[,rownames(terms.table)[ind.j][i]])) {
					levels.i <- levels(data[,rownames(terms.table)[ind.j][i]])
					if(i == 1) lab.i <- matrix(levels.i, ncol = 1)
						if(i > 1) {
						lab.t <- NULL
						for (k in 2:i) lab.t <- cbind(lab.t,
							matrix(rep(lab.i[,k-1], times = length(levels.i)), ncol = 1))
						lab.tmp <- matrix(rep(levels.i, each = length(lab.i[,1])), ncol = 1)
						lab.i <- cbind(lab.t,lab.tmp)
					}
				}
				else {levels.i <- matrix("")
					if(i == 1) lab.i <- matrix(levels.i, ncol = 1)
					if(i > 1) {
						lab.t <- NULL
						for (k in 2:i) lab.t <- cbind(lab.t,
							matrix(rep(lab.i[,k-1], times = length(levels.i)), ncol = 1))
						lab.tmp <- matrix(rep(levels.i, each = length(lab.i[,1])), ncol = 1)
						lab.i <- cbind(lab.t,lab.tmp)
					}
				}
			}
			cum.ind <- rep(TRUE, times = length(lab.i[,1]))
			levels.lab.i <- NULL
			for (i in 1:n.interact) {
				if(!any(rownames(terms.table)[ind.j][i] == colnames(data))) cum.ind <- cum.ind & TRUE
				else if(is.factor(data[,rownames(terms.table)[ind.j][i]])){
					if (terms.table[ind.j,j][i] == 1) {		
						cum.ind <- cum.ind & lab.i[,i] != 
							levels(data[,rownames(terms.table)[ind.j][i]])[1]
					}
				}
				else cum.ind <- cum.ind & TRUE
				if(i == 1) levels.lab.i <- lab.i[,1]
				if(i > 1) levels.lab.i <- paste(levels.lab.i, ",", lab.i[,i], sep = "")
			}
			n.fee <- sum(cum.ind)
			n.j <- length(cum.ind)
			fee.xtra <- matrix(c(0,NA,NA,NA,NA), nrow = n.j, ncol = 5, byrow = TRUE)
			fee.xtra[cum.ind,1:5] <- matrix(unlist(fee[(fee.i + 1):(fee.i + n.fee), 1:5]), ncol = 5)
			fee.xtra1 <- data.frame(factor = colnames(terms.table)[j], levels = levels.lab.i, 
				parameters = fee.xtra[,1], std.err = fee.xtra[,2], 
				df = fee.xtra[,3], t.value = fee.xtra[,4],
				prob.t = fee.xtra[,5])
			if(is.null(levels.lab)) levels.lab <- fee.xtra1
			else levels.lab <- rbind(levels.lab, fee.xtra1)		
			fee.i <- fee.i + n.fee
		}
	}
	rownames(levels.lab) <- 1:length(levels.lab[,1])
	colnames(levels.lab) <- c("effect", "levels", "estimate", "std.err", "df", "t.value", "prob.t")
	levels.lab
}



#------------ FUNCTION FOR ESTIMATES AND CONTRASTS, CHECKS FOR ESTIMABILITY


estimates.and.contrasts <- function(formula = formula, data = data, 
	est.con.mat = est.con.mat, est.con.lab = est.con.label, 
	fixed.effects.estimates = fixed.effects.estimates,
	cov.fix.eff = cov.fix.eff, p = p, n = n) 
{
		Xop <- X.design(formula = formula, data = data, Xmethod = "over.par")
		estimability <- matrix(FALSE, nrow = length(est.con.mat[,1]), ncol = 1)
		for(i in 1:length(est.con.mat[,1])) {
			p1 <- sum(svd(rbind(Xop, est.con.mat[i,]))$d > 1e-10)
			estimability[i,] <- p == p1
		}
		estimates.contrasts <- est.con.mat %*% matrix(fixed.effects.estimates[,3], ncol = 1)
		rownames(estimates.contrasts) <- 1:length(est.con.mat[,1])
		estimates.contrasts[!estimability] <- NA
		ec.var <- est.con.mat %*% cov.fix.eff %*% t(est.con.mat)
		std.err <- sqrt(diag(ec.var))
		std.err[!estimability] <- NA
		estimates.contrasts <- cbind(estimates.contrasts, std.err)
		np <- rep(n - p, times = length(est.con.mat[,1]))
		np[!estimability] <- NA
		estimates.contrasts <- cbind(estimates.contrasts, np)
		estimates.contrasts <- cbind(estimates.contrasts, 
			estimates.contrasts[,1]/estimates.contrasts[,2])
		estimates.contrasts <- cbind(estimates.contrasts,
			round(100000*(1- pt(abs(estimates.contrasts[,4]), df = n - p))*2)/100000)
		estimates.contrasts <- as.data.frame(estimates.contrasts)
		colnames(estimates.contrasts) <- c("estimate", "std.err", "df", "t.value", "prob.t")
		estimates.contrasts <- data.frame(label = est.con.lab, 
			estimable = estimability, estimates.contrasts)
		estimates.contrasts
}


#------------ FUNCTION FOR COMPUTING TYPE III GENERAL LINEAR HYPOTHESES

hypoth.stream <- function(formula, data, Vi)
{
	termlabs <- attr(terms(formula),"term.labels")
	if(length(termlabs) == 0) return(NULL)
	z <- data[,"work.col"]
	z.col <- "work.col"
	form1 <- termlabs[1]
	if(length(termlabs) > 1) for(i in 2:length(termlabs)) form1 <- paste(form1, "+", termlabs[i])
	formula <- as.formula(paste(z.col, "~", form1))
	X <- stdX(X.design(formula = formula, data = data, Xmethod = "sum.to.0"))
	terms.table <- attr(terms(formula),"factors")
	n <- length(z)
	p <- sum(svd(X)$d>1e-10)
	covb.sum0 <- mginv(t(X) %*% Vi %*% X)
	b.sum0 <- covb.sum0 %*% t(X) %*% Vi %*% z
	nterms <- length(terms.table[1,])
	cumcol <- 1
	storage <- matrix(0, ncol = 4, nrow = nterms)
	nfactor <- 0
	for(i in 1:nterms) {
		facts <- rownames(terms.table)[terms.table[,i] > 0]
		nfacts <- length(facts)
		nxcol <- 1
		for(j in 1:nfacts)
			if(is.factor(data[,facts[j]]))
				nxcol <- nxcol*(length(levels(data[,facts[j]])) - 1)
		K <- matrix(0, nrow = p, ncol = nxcol)
		K[(cumcol + 1):(cumcol + nxcol), 1:nxcol] <- diag(nxcol)
		cumcol <- cumcol + nxcol
		storage[i,1] <- nxcol
		storage[i,2] <- n - p
		storage[i,3] <- t(t(K) %*% b.sum0) %*% mginv(t(K) %*% covb.sum0 %*% K) %*%
			(t(K) %*% b.sum0)/sum(svd(K)$d>1e-10)
		storage[i,4] <- round(100000*(1 - pf(storage[i,3],storage[i,1], storage[i,2])))/100000
	}
	storage <- as.data.frame(storage)
	colnames(storage) <- c("Num.df", "Den.df", "F.value", "Prob.F")
	outpt <- data.frame(Effect = termlabs, storage)
	outpt
}


# -------------------- SCALED DISTANCE USING ANISOTROPY

dist.geo.ani <- function(xcoord, ycoord, range, minorp, rotate)
{
	if(length(xcoord) != length(ycoord)){
		return(list(errstate = 1,
			err.message = "number of x-coord not equal to number of y-coord",
			err.extra = data.frame(n.xcoord = length(xcoord), n.ycoord = length(ycoord))))
	}
	if(range < 0){
		return(list(errstate = 1,
			err.message = "range parameter less than 0",
			err.extra = data.frame(range.parm = range)))
	}
	if(rotate < 0 || rotate > 180){
		return(list(errstate = 1,
			err.message = "rotation parameter beyond 0 - 180 range",
			err.extra = data.frame(rotate.parm = rotate)))
	}
	if(minorp < 0 || minorp > 1){
		return(list(errstate = 1,
			err.message = "minor range proportion beyond 0 - 1 range",
			err.extra = dataframe(minorp.parm = minorp)))
	}
	# total number of observations
		n <- length(xcoord)
	# expand all x-coordinates
		sx <- matrix(xcoord,ncol = 1) %*%
			matrix(rep(1,times = n), nrow = 1)
	# find difference in x-coordinates between all pairwise locations
		sxdif <- t(sx) - sx
	# expand all y-coordinates
		sy <- matrix(ycoord,ncol = 1) %*%
			matrix(rep(1,times = n), nrow = 1)
	# find difference in x-coordinates between all pairwise locations
		sydif <- t(sy) - sy
	# rotate coordinates
		newx <- cos(rotate*.0174533)*sxdif - sin(rotate*.0174533)*sydif
		newy <- sin(rotate*.0174533)*sxdif + cos(rotate*.0174533)*sydif
	# scale coordinates by minor and major axes */
		newx <- newx/(range*minorp)
		newy <- newy/range
	# compute distance for the scaled and rotated coordinates */
		sqrt(newx^2 + newy^2)
}

# --------------- SCALED DISTANCE USING ISOTROPY

dist.geo.iso <- function(xcoord, ycoord, range)
{
	# total number of observations
		n <- length(xcoord)
	distance <- matrix(0, n, n)
	distance[lower.tri(distance)] <- dist(as.matrix(cbind(xcoord,ycoord)))
	(distance + t(distance))/range
}

#------------ MAKE COVARIANCE MATRIX

# ----------- EUCLIDEAN DISTANCE COVARIANCE MATRICES

exp.Euclid <- function(parsil, distance.matrix)
{
	parsil*exp(-3*distance.matrix)
}
sph.Euclid <- function(parsil, distance.matrix)
{
	CovMat <- (1 - 1.5*distance.matrix + 0.5*distance.matrix^3)
	CovMat[distance.matrix > 1] <- 0
	parsil*CovMat
}
gau.Euclid <- function(parsil, distance.matrix)
{
	parsil*exp(-distance.matrix^2)
}
cau.Euclid <- function(parsil, distance.matrix)
{
	parsil/(1 + distance.matrix^2)
}

# ----------- STREAM NETWORK COVARIANCE MATRICES

# ----------- Tail-up Models

exp.tailup <- function(dist.hydro, weight,
	parsil = parsil, range1 = range1)
{
	n <- length(dist.hydro[,1])
	parsil*exp(-3*dist.hydro/range1)*weight
}

sph.tailup <- function(dist.hydro, weight,
	parsil = parsil, range1 = range1)
{
	n <- length(dist.hydro[,1])
	parsil*(matrix(rep(1, times = n^2), nrow = n) - 1.5*(dist.hydro/range1) +
		0.5*(dist.hydro/range1)^3)*(dist.hydro < range1)*weight
}

lws.tailup <- function(dist.hydro, weight,
	parsil = parsil, range1 = range1)
{
	n <- length(dist.hydro[,1])
	parsil*(matrix(rep(1, times = n^2), nrow = n) - dist.hydro/range1)*
		(dist.hydro < range1)*weight
}

mar.tailup <- function(dist.hydro, weight,
	parsil = parsil, range1 = range1)
{
	n <- length(dist.hydro[,1])
	CovMat <- matrix(0, nrow = n, ncol = n)
	ind <- (1e-9 < dist.hydro/range1)
	CovMat[ind] <- (parsil*log(dist.hydro/range1 + 1)/(dist.hydro/range1)*weight)[ind]
	ind <- 1e-9 >= dist.hydro/range1
	CovMat[ind] <-	(parsil*weight)[ind]
	CovMat
}

# ----------- Tail-down Models

exp.taildown <- function(dist.hydro,
	parsil = parsil, range1 = range1)
{
	parsil*exp(-3*dist.hydro/range1)
}

sph.taildown <- function(dist.hydro, a.mat, b.mat,
	parsil = parsil, range1 = range1)
{
	flow.connect <- b.mat == 0
	n <- length(a.mat[,1])
	parsil*(matrix(rep(1, times = n^2), nrow = n) - 1.5*(dist.hydro/range1) +
		0.5*(dist.hydro/range1)^3)*(a.mat < range1)*flow.connect +
		parsil*(matrix(rep(1, times = n^2), nrow = n) - 1.5*b.mat/range1 +
		0.5*a.mat/range1)*(matrix(rep(1, times = n^2), nrow = n) - a.mat/range1)^2*
		(a.mat < range1)*(1 - flow.connect)
}

lws.taildown <- function(dist.hydro, a.mat, b.mat,
	parsil = parsil, range1 = range1)
{
	flow.connect <- b.mat == 0
	n <- length(a.mat[,1])
	parsil*(matrix(rep(1, times = n^2), nrow = n) - dist.hydro/range1)*
		(a.mat < range1)*flow.connect + 
		parsil*(matrix(rep(1, times = n^2), nrow = n) - a.mat/range1)*
		(a.mat < range1)*(1 - flow.connect)
}

mar.taildown <- function(dist.hydro, a.mat, b.mat,
	parsil = parsil, range1 = range1)
{
	n <- length(a.mat[,1])
  covmat <- matrix(0, ncol = n, nrow = n)
	flow.connect <- b.mat == 0
	covmat <- parsil*( (a.mat == 0)*1 +
		log(dist.hydro/range1 + 1)/(dist.hydro/range1)*flow.connect +
		(log(a.mat/range1 + 1) - log(b.mat/range1 + 1))/((a.mat - b.mat)/range1)*
			(a.mat > b.mat)*(1 - flow.connect) )
  covmat[is.nan(covmat)] <- 0
	covmat <- covmat +
      parsil/(a.mat/range1 + 1)*(a.mat == b.mat)*(1 - flow.connect)
	diag(covmat) <- parsil
	covmat	
}



#------------ FUNCTION TO PUT THETA BACK IN ORIGINAL PARAMETER SPACE

untrans.theta <- function(theta, CorModel1, CorModel2, CorModel3,
  use.anisotropy, use.nugget)
{
  # put theta back on untransformed scale
  ntheta.sofar <- 0
  if(!is.null(CorModel1)) {
  if(substr(unlist(strsplit(CorModel1,"\\."))[2],1,1) == "t")  {
    ntheta.sofar <- ntheta.sofar + 1
    theta[ntheta.sofar] <- exp(theta[ntheta.sofar])
    ntheta.sofar <- ntheta.sofar + 1
    theta[ntheta.sofar] <- exp(theta[ntheta.sofar])
  }
  else
    if(use.anisotropy == FALSE)  {
    ntheta.sofar <- ntheta.sofar + 1
    theta[ntheta.sofar] <- exp(theta[ntheta.sofar])
    ntheta.sofar <- ntheta.sofar + 1
    theta[ntheta.sofar] <- exp(theta[ntheta.sofar])
    }
    else   {
    ntheta.sofar <- ntheta.sofar + 1
    theta[ntheta.sofar] <- exp(theta[ntheta.sofar])
    ntheta.sofar <- ntheta.sofar + 1
    theta[ntheta.sofar] <- exp(theta[ntheta.sofar])
    ntheta.sofar <- ntheta.sofar + 1
    theta[ntheta.sofar] <- exp(theta[ntheta.sofar])/(1 +
      exp(theta[ntheta.sofar]))
    ntheta.sofar <- ntheta.sofar + 1
    theta[ntheta.sofar] <- 180*exp(theta[ntheta.sofar])/(1 +
      exp(theta[ntheta.sofar]))
    }
  }
  if(!is.null(CorModel2)) {
  if(substr(unlist(strsplit(CorModel2,"\\."))[2],1,1) == "t")  {
    ntheta.sofar <- ntheta.sofar + 1
    theta[ntheta.sofar] <- exp(theta[ntheta.sofar])
    ntheta.sofar <- ntheta.sofar + 1
    theta[ntheta.sofar] <- exp(theta[ntheta.sofar])
  }
  else
    if(use.anisotropy == FALSE)  {
    ntheta.sofar <- ntheta.sofar + 1
    theta[ntheta.sofar] <- exp(theta[ntheta.sofar])
    ntheta.sofar <- ntheta.sofar + 1
    theta[ntheta.sofar] <- exp(theta[ntheta.sofar])
    }
    else   {
    ntheta.sofar <- ntheta.sofar + 1
    theta[ntheta.sofar] <- exp(theta[ntheta.sofar])
    ntheta.sofar <- ntheta.sofar + 1
    theta[ntheta.sofar] <- exp(theta[ntheta.sofar])
    ntheta.sofar <- ntheta.sofar + 1
    theta[ntheta.sofar] <- exp(theta[ntheta.sofar])/(1 +
      exp(theta[ntheta.sofar]))
    ntheta.sofar <- ntheta.sofar + 1
    theta[ntheta.sofar] <- 180*exp(theta[ntheta.sofar])/(1 +
      exp(theta[ntheta.sofar]))
    }
  }
  if(!is.null(CorModel3)) {
  if(substr(unlist(strsplit(CorModel3,"\\."))[2],1,1) == "t")  {
    ntheta.sofar <- ntheta.sofar + 1
    theta[ntheta.sofar] <- exp(theta[ntheta.sofar])
    ntheta.sofar <- ntheta.sofar + 1
    theta[ntheta.sofar] <- exp(theta[ntheta.sofar])
  }
  else
    if(use.anisotropy == FALSE)  {
    ntheta.sofar <- ntheta.sofar + 1
    theta[ntheta.sofar] <- exp(theta[ntheta.sofar])
    ntheta.sofar <- ntheta.sofar + 1
    theta[ntheta.sofar] <- exp(theta[ntheta.sofar])
    }
    else   {
    ntheta.sofar <- ntheta.sofar + 1
    theta[ntheta.sofar] <- exp(theta[ntheta.sofar])
    ntheta.sofar <- ntheta.sofar + 1
    theta[ntheta.sofar] <- exp(theta[ntheta.sofar])
    ntheta.sofar <- ntheta.sofar + 1
    theta[ntheta.sofar] <- exp(theta[ntheta.sofar])/(1 +
      exp(theta[ntheta.sofar]))
    ntheta.sofar <- ntheta.sofar + 1
    theta[ntheta.sofar] <- 180*exp(theta[ntheta.sofar])/(1 +
      exp(theta[ntheta.sofar]))
    }
  }
	if(use.nugget == TRUE) {
    ntheta.sofar <- ntheta.sofar + 1
    theta[ntheta.sofar] <- exp(theta[ntheta.sofar])
  }
  theta
}

#------------ CREATE COVARIANCE MATRIX

make.cov.mat <- function(theta, dist.hydro, a.mat, b.mat, w.matrix,
    net.zero, xcoord, ycoord,
    CorModel1, CorModel2, CorModel3,
    use.nugget, use.anisotropy)
{

	# number of CorModels
	n.models <- 1 + !is.null(CorModel2) + !is.null(CorModel3)
	# dimension of covariance matrix
  if(!is.null(dist.hydro)) nn <- length(dist.hydro[,1])
  else if(!is.null(w.matrix)) nn <- length(w.matrix[,1])
  else if(!is.null(xcoord)) nn <- length(xcoord)

  V <- matrix(0, nrow = nn, ncol = nn)
  npar.sofar <- 0
  # CorModel1 possibilities
  if(!is.null(CorModel1)) {
  if(CorModel1 == "Exponential.taildown") {
    V <- V + exp.taildown(dist.hydro = dist.hydro,
      parsil = theta[1], range1 = theta[2])
    npar.sofar <- npar.sofar + 2
    if(!is.null(net.zero)) V <- V*net.zero
  }
  if(CorModel1 == "Spherical.taildown") {
    V <- V + sph.taildown(dist.hydro = dist.hydro,
      a.mat = a.mat, b.mat = b.mat,
      parsil = theta[1], range1 = theta[2])
    npar.sofar <- npar.sofar + 2
    if(!is.null(net.zero)) V <- V*net.zero
  }
  if(CorModel1 == "LinearSill.taildown") {
    V <- V + lws.taildown(dist.hydro = dist.hydro,
      a.mat = a.mat, b.mat = b.mat,
      parsil = theta[1], range1 = theta[2])
    npar.sofar <- npar.sofar + 2
    if(!is.null(net.zero)) V <- V*net.zero
  }
  if(CorModel1 == "Mariah.taildown") {
    V <- V + mar.taildown(dist.hydro = dist.hydro,
      a.mat = a.mat, b.mat = b.mat,
      parsil = theta[1], range1 = theta[2])
    npar.sofar <- npar.sofar + 2
    if(!is.null(net.zero)) V <- V*net.zero
  }
  if(CorModel1 == "Exponential.tailup") {
    V <- V + exp.tailup(dist.hydro = dist.hydro, weight = w.matrix,
      parsil = theta[1], range1 = theta[2])
    npar.sofar <- npar.sofar + 2
  }
  if(CorModel1 == "Spherical.tailup") {
    V <- V + sph.tailup(dist.hydro = dist.hydro, weight = w.matrix,
      parsil = theta[1], range1 = theta[2])
    npar.sofar <- npar.sofar + 2
  }
  if(CorModel1 == "LinearSill.tailup") {
    V <- V + lws.tailup(dist.hydro = dist.hydro, weight = w.matrix,
      parsil = theta[1], range1 = theta[2])
    npar.sofar <- npar.sofar + 2
  }
  if(CorModel1 == "Mariah.tailup") {
    V <- V + mar.tailup(dist.hydro = dist.hydro, weight = w.matrix,
      parsil = theta[1], range1 = theta[2])
    npar.sofar <- npar.sofar + 2
  }
  if(CorModel1 == "Exponential.Euclid") {
    if(use.anisotropy == FALSE) {
      dist.mat <- dist.geo.iso(xcoord, ycoord, theta[2])
      npar.sofar <- npar.sofar + 2
    }
    else {
      dist.mat <- dist.geo.ani(xcoord, ycoord, theta[2],
        theta[3], theta[4])
      npar.sofar <- npar.sofar + 4
    }
    V <- V + exp.Euclid(dist.mat,
      parsil = theta[1])
  }
  if(CorModel1 == "Spherical.Euclid") {
    if(use.anisotropy == FALSE) {
      dist.mat <- dist.geo.iso(xcoord, ycoord, theta[2])
      npar.sofar <- npar.sofar + 2
    }
    else {
      dist.mat <- dist.geo.ani(xcoord, ycoord, theta[2],
        theta[3], theta[4])
      npar.sofar <- npar.sofar + 4
    }
    V <- V + sph.Euclid(dist.mat,
      parsil = theta[1])
  }
  if(CorModel1 == "Gaussian.Euclid") {
    if(use.anisotropy == FALSE) {
      dist.mat <- dist.geo.iso(xcoord, ycoord, theta[2])
      npar.sofar <- npar.sofar + 2
    }
    else {
      dist.mat <- dist.geo.ani(xcoord, ycoord, theta[2],
        theta[3], theta[4])
      npar.sofar <- npar.sofar + 4
    }
    V <- V + gau.Euclid(dist.mat,
      parsil = theta[1])
  }
  if(CorModel1 == "Cauchy.Euclid") {
    if(use.anisotropy == FALSE) {
      dist.mat <- dist.geo.iso(xcoord, ycoord, theta[2])
      npar.sofar <- npar.sofar + 2
    }
    else {
      dist.mat <- dist.geo.ani(xcoord, ycoord, theta[2],
        theta[3], theta[4])
      npar.sofar <- npar.sofar + 4
    }
    V <- V + cau.Euclid(dist.mat,
      parsil = theta[1])
  }
  }
  # CorModel2 possibilities
  if(!is.null(CorModel2)) {
  if(CorModel2 == "Exponential.taildown") {
    V <- V + exp.taildown(dist.hydro = dist.hydro,
      parsil = theta[npar.sofar + 1], range1 = theta[npar.sofar + 2])
    npar.sofar <- npar.sofar + 2
    if(!is.null(net.zero)) V <- V*net.zero
  }
  if(CorModel2 == "Spherical.taildown") {
    V <- V + sph.taildown(dist.hydro = dist.hydro,
      a.mat = a.mat, b.mat = b.mat,
      parsil = theta[npar.sofar + 1], range1 = theta[npar.sofar + 2])
    npar.sofar <- npar.sofar + 2
    if(!is.null(net.zero)) V <- V*net.zero
  }
  if(CorModel2 == "LinearSill.taildown") {
    V <- V + lws.taildown(dist.hydro = dist.hydro,
      a.mat = a.mat, b.mat = b.mat,
      parsil = theta[npar.sofar + 1], range1 = theta[npar.sofar + 2])
    npar.sofar <- npar.sofar + 2
    if(!is.null(net.zero)) V <- V*net.zero
  }
  if(CorModel2 == "Mariah.taildown") {
    V <- V + mar.taildown(dist.hydro = dist.hydro,
      a.mat = a.mat, b.mat = b.mat,
      parsil = theta[npar.sofar + 1], range1 = theta[npar.sofar + 2])
    npar.sofar <- npar.sofar + 2
    if(!is.null(net.zero)) V <- V*net.zero
  }
  if(CorModel2 == "Exponential.tailup") {
    V <- V + exp.tailup(dist.hydro = dist.hydro, weight = w.matrix,
      parsil = theta[npar.sofar + 1], range1 = theta[npar.sofar + 2])
    npar.sofar <- npar.sofar + 2
  }
  if(CorModel2 == "Spherical.tailup") {
    V <- V + sph.tailup(dist.hydro = dist.hydro, weight = w.matrix,
      parsil = theta[npar.sofar + 1], range1 = theta[npar.sofar + 2])
    npar.sofar <- npar.sofar + 2
  }
  if(CorModel2 == "LinearSill.tailup") {
    V <- V + lws.tailup(dist.hydro = dist.hydro, weight = w.matrix,
      parsil = theta[npar.sofar + 1], range1 = theta[npar.sofar + 2])
    npar.sofar <- npar.sofar + 2
  }
  if(CorModel2 == "Mariah.tailup") {
    V <- V + mar.tailup(dist.hydro = dist.hydro, weight = w.matrix,
      parsil = theta[npar.sofar + 1], range1 = theta[npar.sofar + 2])
    npar.sofar <- npar.sofar + 2
  }
  if(CorModel2 == "Exponential.Euclid") {
    npar.parsil <- npar.sofar + 1
    if(use.anisotropy == FALSE) {
      dist.mat <- dist.geo.iso(xcoord, ycoord, theta[npar.sofar + 2])
      npar.sofar <- npar.sofar + 2
    }
    else {
      dist.mat <- dist.geo.ani(xcoord, ycoord, theta[npar.sofar + 2],
        theta[npar.sofar + 3], theta[npar.sofar + 4])
      npar.sofar <- npar.sofar + 4
    }
    V <- V + exp.Euclid(dist.mat,
      parsil = theta[npar.parsil])
  }
  if(CorModel2 == "Spherical.Euclid") {
    npar.parsil <- npar.sofar + 1
    if(use.anisotropy == FALSE) {
      dist.mat <- dist.geo.iso(xcoord, ycoord, theta[npar.sofar + 2])
      npar.sofar <- npar.sofar + 2
    }
    else {
      dist.mat <- dist.geo.ani(xcoord, ycoord, theta[npar.sofar + 2],
        theta[npar.sofar + 3], theta[npar.sofar + 4])
      npar.sofar <- npar.sofar + 4
    }
    V <- V + sph.Euclid(dist.mat,
      parsil = theta[npar.parsil])
  }
  if(CorModel2 == "Gaussian.Euclid") {
    npar.parsil <- npar.sofar + 1
    if(use.anisotropy == FALSE) {
      dist.mat <- dist.geo.iso(xcoord, ycoord, theta[npar.sofar + 2])
      npar.sofar <- npar.sofar + 2
    }
    else {
      dist.mat <- dist.geo.ani(xcoord, ycoord, theta[npar.sofar + 2],
        theta[npar.sofar + 3], theta[npar.sofar + 4])
      npar.sofar <- npar.sofar + 4
    }
    V <- V + gau.Euclid(dist.mat,
      parsil = theta[npar.parsil])
  }
  if(CorModel2 == "Cauchy.Euclid") {
    npar.parsil <- npar.sofar + 1
    if(use.anisotropy == FALSE) {
      dist.mat <- dist.geo.iso(xcoord, ycoord, theta[npar.sofar + 2])
      npar.sofar <- npar.sofar + 2
    }
    else {
      dist.mat <- dist.geo.ani(xcoord, ycoord, theta[npar.sofar + 2],
        theta[npar.sofar + 3], theta[npar.sofar + 4])
      npar.sofar <- npar.sofar + 4
    }
    V <- V + cau.Euclid(dist.mat,
      parsil = theta[npar.parsil])
  }
  }
  # CorModel3 possibilities
  if(!is.null(CorModel3)) {
  if(CorModel3 == "Exponential.taildown") {
    V <- V + exp.taildown(dist.hydro = dist.hydro,
      parsil = theta[npar.sofar + 1], range1 = theta[npar.sofar + 2])
    npar.sofar <- npar.sofar + 2
    if(!is.null(net.zero)) V <- V*net.zero
  }
  if(CorModel3 == "Spherical.taildown") {
    V <- V + sph.taildown(dist.hydro = dist.hydro,
      a.mat = a.mat, b.mat = b.mat,
      parsil = theta[npar.sofar + 1], range1 = theta[npar.sofar + 2])
    npar.sofar <- npar.sofar + 2
    if(!is.null(net.zero)) V <- V*net.zero
  }
  if(CorModel3 == "LinearSill.taildown") {
    V <- V + lws.taildown(dist.hydro = dist.hydro,
      a.mat = a.mat, b.mat = b.mat,
      parsil = theta[npar.sofar + 1], range1 = theta[npar.sofar + 2])
    npar.sofar <- npar.sofar + 2
    if(!is.null(net.zero)) V <- V*net.zero
  }
  if(CorModel3 == "Mariah.taildown") {
    V <- V + mar.taildown(dist.hydro = dist.hydro,
      a.mat = a.mat, b.mat = b.mat,
      parsil = theta[npar.sofar + 1], range1 = theta[npar.sofar + 2])
    npar.sofar <- npar.sofar + 2
    if(!is.null(net.zero)) V <- V*net.zero
  }
  if(CorModel3 == "Exponential.tailup") {
    V <- V + exp.tailup(dist.hydro = dist.hydro, weight = w.matrix,
      parsil = theta[npar.sofar + 1], range1 = theta[npar.sofar + 2])
    npar.sofar <- npar.sofar + 2
  }
  if(CorModel3 == "Spherical.tailup") {
    V <- V + sph.tailup(dist.hydro = dist.hydro, weight = w.matrix,
      parsil = theta[npar.sofar + 1], range1 = theta[npar.sofar + 2])
    npar.sofar <- npar.sofar + 2
  }
  if(CorModel3 == "LinearSill.tailup") {
    V <- V + lws.tailup(dist.hydro = dist.hydro, weight = w.matrix,
      parsil = theta[npar.sofar + 1], range1 = theta[npar.sofar + 2])
    npar.sofar <- npar.sofar + 2
  }
  if(CorModel3 == "Mariah.tailup") {
    V <- V + mar.tailup(dist.hydro = dist.hydro, weight = w.matrix,
      parsil = theta[npar.sofar + 1], range1 = theta[npar.sofar + 2])
    npar.sofar <- npar.sofar + 2
  }
  if(CorModel3 == "Exponential.Euclid") {
    npar.parsil <- npar.sofar + 1
    if(use.anisotropy == FALSE) {
      dist.mat <- dist.geo.iso(xcoord, ycoord, theta[npar.sofar + 2])
      npar.sofar <- npar.sofar + 2
    }
    else {
      dist.mat <- dist.geo.ani(xcoord, ycoord, theta[npar.sofar + 2],
        theta[npar.sofar + 3], theta[npar.sofar + 4])
      npar.sofar <- npar.sofar + 4
    }
    V <- V + exp.Euclid(dist.mat,
      parsil = theta[npar.sofar + 1])
  }
  if(CorModel3 == "Spherical.Euclid") {
    npar.parsil <- npar.sofar + 1
    if(use.anisotropy == FALSE) {
      dist.mat <- dist.geo.iso(xcoord, ycoord, theta[npar.sofar + 2])
      npar.sofar <- npar.sofar + 2
    }
    else {
      dist.mat <- dist.geo.ani(xcoord, ycoord, theta[npar.sofar + 2],
        theta[npar.sofar + 3], theta[npar.sofar + 4])
      npar.sofar <- npar.sofar + 4
    }
    V <- V + sph.Euclid(dist.mat,
      parsil = theta[npar.parsil])
  }
  if(CorModel3 == "Gaussian.Euclid") {
    npar.parsil <- npar.sofar + 1
    if(use.anisotropy == FALSE) {
      dist.mat <- dist.geo.iso(xcoord, ycoord, theta[npar.sofar + 2])
      npar.sofar <- npar.sofar + 2
    }
    else {
      dist.mat <- dist.geo.ani(xcoord, ycoord, theta[npar.sofar + 2],
        theta[npar.sofar + 3], theta[npar.sofar + 4])
      npar.sofar <- npar.sofar + 4
    }
    V <- V + gau.Euclid(dist.mat,
      parsil = theta[npar.parsil])
  }
  if(CorModel3 == "Cauchy.Euclid") {
    npar.parsil <- npar.sofar + 1
    if(use.anisotropy == FALSE) {
      dist.mat <- dist.geo.iso(xcoord, ycoord, theta[npar.sofar + 2])
      npar.sofar <- npar.sofar + 2
    }
    else {
      dist.mat <- dist.geo.ani(xcoord, ycoord, theta[npar.sofar + 2],
        theta[npar.sofar + 3], theta[npar.sofar + 4])
      npar.sofar <- npar.sofar + 4
    }
    V <- V + cau.Euclid(dist.mat,
      parsil = theta[npar.parsil])
  }
  }
  if(use.nugget == TRUE) {
    npar.sofar <- npar.sofar + 1
  	V <- V + diag(theta[npar.sofar], nrow = nn, ncol = nn)
  }
  V
}


#------------ CHECK ANY MIS-SPECIFICATIONS IN THE ARGUMENTS


arg.error.check <- function(formula, data,
  dist.junc = NULL,
  flow.matrix = NULL,
  xcol = NULL,
  ycol = NULL,
  CorModel1 = "Exponential.taildown",
  CorModel2 = NULL,
  CorModel3 = NULL,
  use.nugget = TRUE,
  use.anisotropy = FALSE,
  Distribution = "Gaussian",
	EstMeth = "REML",
	prediction.indices = NULL)
{
  if(!is.logical(use.nugget))
  return(list(Err = 1, message =
    "use.nugget argument must be TRUE or FALSE"))
  if(!is.logical(use.anisotropy))
  return(list(Err = 1, message =
    "use.anisotropy argument must be TRUE or FALSE"))
  if(EstMeth != "ML" & EstMeth != "REML")
  return(list(Err = 1, message =
    "EstMeth must be ML (Max. Likelihood) or REML (Restricted Max. Likelihood)"))
  if(!(Distribution == "Gaussian" | Distribution == "Binomial" |
    Distribution == "Poisson"))
  return(list(Err = 1, message =
    "distribution must be Gaussian or Binomial or Poisson"))

	# check to make sure we have necessary arguments for stream and
	# Euclidean distance models
  if(!is.null(CorModel1)) {
	if((CorModel1 == "Spherical.Euclid" | CorModel1 == "LinearSill.Euclid" |
		CorModel1 == "Exponential.Euclid" | CorModel1 == "Cauchy.Euclid" |
    CorModel1 == "Gaussian.Euclid") & (is.null(xcol) | is.null(ycol)) )
    return(list(Err = 1, message =
    "Need x- and y-coordinates [xcol and ycol arguments] for Euclidean distance models"))
	if((CorModel1 == "Spherical.taildown" | CorModel1 == "Exponential.taildown" |
    CorModel1 == "LinearSill.taildown" | CorModel1 == "Mariah.taildown" |
    CorModel1 == "Spherical.tailup" | CorModel1 == "Exponential.tailup" |
    CorModel1 == "LinearSill.tailup"  | CorModel1 == "Mariah.tailup")
    & is.null(dist.junc) )
    return(list(Err = 1, message =
    "Need distance-to-junction [dist.junc argument] matrix for stream distance models"))
	if((CorModel1 == "Spherical.tailup" | CorModel1 == "Exponential.tailup" |
    CorModel1 == "LinearSill.tailup"  | CorModel1 == "Mariah.tailup" )
    & is.null(flow.matrix) )
    return(list(Err = 1, message =
    "Need flow matrix [flow.matrix argument] for tail-up stream flow models"))
  }
  
  if(!is.null(CorModel2)) {
	if((CorModel2 == "Spherical.Euclid" | CorModel2 == "LinearSill.Euclid" |
		CorModel2 == "Exponential.Euclid" | CorModel2 == "Cauchy.Euclid" |
    CorModel2 == "Gaussian.Euclid") & (is.null(xcol) | is.null(ycol)) )
    return(list(Err = 1, message =
     "Need x- and y-coordinates [xcol and ycol arguments] for Euclidean distance models"))
	if((CorModel2 == "Spherical.taildown" | CorModel2 == "Exponential.taildown" |
    CorModel2 == "LinearSill.taildown" | CorModel2 == "Mariah.taildown" |
    CorModel2 == "Spherical.tailup" | CorModel2 == "Exponential.tailup" |
    CorModel2 == "LinearSill.tailup"  | CorModel2 == "Mariah.tailup")
    & is.null(dist.junc) )
    return(list(Err = 1, message =
    "Need distance-to-junction [dist.junc argument] matrix for stream distance models"))
	if((CorModel2 == "Spherical.tailup" | CorModel2 == "Exponential.tailup" |
    CorModel2 == "LinearSill.tailup"  | CorModel2 == "Mariah.tailup" )
    & is.null(flow.matrix) )
    return(list(Err = 1, message =
    "Need flow matrix [flow.matrix argument] for tail-up stream flow models"))
  }

  if(!is.null(CorModel3)) {
  if((CorModel3 == "Spherical.Euclid" | CorModel3 == "LinearSill.Euclid" |
		CorModel3 == "Exponential.Euclid" | CorModel3 == "Cauchy.Euclid" |
    CorModel3 == "Gaussian.Euclid") & (is.null(xcol) | is.null(ycol)) )
    return(list(Err = 1, message =
    "Need x- and y-coordinates [xcol and ycol arguments] for Euclidean distance models"))
	if((CorModel3 == "Spherical.taildown" | CorModel3 == "Exponential.taildown" |
    CorModel3 == "LinearSill.taildown" | CorModel3 == "Mariah.taildown" |
    CorModel3 == "Spherical.tailup" | CorModel3 == "Exponential.tailup" |
    CorModel3 == "LinearSill.tailup"  | CorModel3 == "Mariah.tailup")
    & is.null(dist.junc) )
    return(list(Err = 1, message =
    "Need distance-to-junction [dist.junc argument] matrix for stream distance models"))
	if((CorModel3 == "Spherical.tailup" | CorModel3 == "Exponential.tailup" |
    CorModel3 == "LinearSill.tailup"  | CorModel3 == "Mariah.tailup" )
    & is.null(flow.matrix) )
    return(list(Err = 1, message =
    "Need flow matrix [flow.matrix argument] for tail-up stream flow models"))
  }

  Err <- list(Err = 0)
  Err
}


#------------ MAKE PREDICTIONS (KRIGE) USING PREDICTION SET

pred.k <- function(z, Xdesign, X.na, ind, c.pred,
	covb, b.hat, Vi, sigma, trans.power, trans.shift)
{
	# covariance matrix between observed and prediction locations
	c.pred <- as.matrix(sigma[ind,!ind])

	# number of prediction locations
	np <- length(c.pred[1,])

	# number of observed data locations
	n <- length(c.pred[,1])

	# fitted variance
	sill <- sigma[1,1]

	# define some matrices to store prediction results
	storage.pred <- matrix(NA, nrow = np, ncol = 1)
	storage.se <- matrix(NA, nrow = np, ncol = 1)

	# cycle through prediction locations and use kriging to make predictions
	# see equations at bottom of page 123, Cressie, 1993
	for(i in 1:np) {
		mat1 <- matrix(as.matrix(X.na[i,]), ncol = 1) - t(Xdesign) %*% Vi %*% c.pred[,i]
		lam <- t(c.pred[,i] + Xdesign %*% covb %*% mat1) %*% Vi
		pred.var <- sill  - t(c.pred[,i]) %*% Vi %*% c.pred[,i] +
			t(mat1) %*% covb %*% mat1
		storage.pred[i] <- lam %*% z
		storage.se[i] <- sqrt(pred.var)
    if(!is.null(trans.power)) {
      ind.na <- !ind
		  if(trans.power == 0){
		  	m <- covb %*% (t(Xdesign) %*% Vi %*%
			  	c.pred[,i] - matrix(X.na[i,], ncol = 1))
			  C0 <- sill
			  muhat <- matrix(X.na[i,], nrow = 1) %*% b.hat
			  V.na <- sigma[ind.na, ind.na]
#			  storage.pred[i] <- exp(storage.pred[i]) +
#				  exp(muhat)*(storage.se[i]^2/2 + t(m) %*% matrix(X.na[i,], ncol = 1)) -
#				  trans.shift
			  storage.pred[i] <- exp(storage.pred[i] +
				  storage.se[i]^2/2 + t(m) %*% matrix(X.na[i,], ncol = 1)) -
				  trans.shift
        # prediction standard error
			  storage.se[i] <- exp(muhat)*storage.se[i]
		  }
		  if(trans.power > 0){
			  m <- covb %*% (t(Xdesign) %*% Vi %*%
				  c.pred[,i] - matrix(X.na[i,], ncol = 1))
			  C0 <- sill
			  muhat <- matrix(X.na[i,], nrow = 1) %*% b.hat
			  V.na <- sigma[ind.na, ind.na]
  		  storage.pred[i] <- (storage.pred[i])^(1/trans.power) +
				  (1/trans.power)*(1/trans.power - 1)*muhat^(1/trans.power - 2)*
          (storage.se[i]^2/2 + t(m) %*% matrix(X.na[i,], ncol = 1)) -
				  trans.shift
        # prediction standard error
			  storage.se[i] <- (1/trans.power)*muhat^(1/trans.power - 1)*storage.se[i]
		  }
		}
	}
 data.frame(predictions = storage.pred, pred.std.err = storage.se)
}

#------------ GET SENSIBLE INITIAL PARAMETER VALUES

theta.ini <- function(z, X, CorModel1, CorModel2, CorModel3, use.nugget,
  use.anisotropy, dist.hydro.data, dist.Euclid.data, n.models)
{
  var.resid <- mean((z - X %*% mginv(t(X) %*% X) %*% t(X) %*% z)^2)
	theta <- NULL
  if(!is.null(CorModel1)) {
  if(substr(unlist(strsplit(CorModel1,"\\."))[2],1,1) == "t")
	  theta <- rbind(theta,matrix(c(log(.9/n.models*var.resid),
                   log(mean(dist.hydro.data))),ncol = 1))
  else
    if(use.anisotropy == FALSE)
      theta <- rbind(theta,matrix(c(log(.9/n.models*var.resid),
                     log(mean(dist.Euclid.data))),ncol = 1))
    else theta <- rbind(theta,matrix(c(log(.9/n.models*var.resid),
                     log(mean(dist.Euclid.data)),0,0),ncol = 1))
  }
  if(!is.null(CorModel2)) {
    if(substr(unlist(strsplit(CorModel2,"\\."))[2],1,1) == "t")
	    theta <- rbind(theta,matrix(c(log(.9/n.models*var.resid),
                     log(mean(dist.hydro.data))),ncol = 1))
    else
      if(use.anisotropy == FALSE)
        theta <- rbind(theta,matrix(c(log(.9/n.models*var.resid),
                       log(mean(dist.Euclid.data))),ncol = 1))
      else theta <- rbind(theta,matrix(c(log(.9/n.models*var.resid),
                       log(mean(dist.Euclid.data)),0,0),ncol = 1))
  }
  if(!is.null(CorModel3)) {
    if(substr(unlist(strsplit(CorModel3,"\\."))[2],1,1) == "t")
	    theta <- rbind(theta,matrix(c(log(.9/n.models*var.resid),
                     log(mean(dist.hydro.data))),ncol = 1))
    else
      if(use.anisotropy == FALSE)
        theta <- rbind(theta,matrix(c(log(.9/n.models*var.resid),
                       log(mean(dist.Euclid.data))),ncol = 1))
      else theta <- rbind(theta,matrix(c(log(.9/n.models*var.resid),
                       log(mean(dist.Euclid.data)),0,0),ncol = 1))
  }
	if(use.nugget == TRUE) theta <- rbind(theta,log(.1*var.resid))

  theta
  
}

#------------ BLOCK KRIGING USING PREDICTION SET

block.krige <- function(z, pred.indices,
	covb, beta.hat, V, Vi, Sigma, X.all, ind.na)
{

	n.pred <- length(pred.indices)
	c.pred <- Sigma[!ind.na, pred.indices]
	sig.pp <- Sigma[pred.indices, pred.indices]
	sig.AA <- mean(sig.pp)
	Xdesign.pred <- as.matrix(X.all[pred.indices,])
	Xdesign <- as.matrix(X.all[!ind.na,])
	c.A <- matrix(rowMeans(c.pred), ncol = 1)
	x.A <- matrix(colMeans(Xdesign.pred), ncol = 1)
	lambda <- matrix(t(c.A) %*% Vi + (t(x.A) - t(c.A) %*% Vi %*% Xdesign) %*%
		covb %*% t(Xdesign) %*% Vi, ncol = 1)
	block.krige.est <- t(lambda) %*% z
	V <- Sigma[!ind.na, !ind.na]
	block.krige.se <- sqrt(sig.AA - 2*t(c.A) %*% lambda +
    t(lambda) %*% V %*% lambda)

	data.frame(block.prediction.estimate = block.krige.est,
    block.prediction.stderror = block.krige.se)
}


# ----------------------- STANDARDIZE THE DESIGN MATRIX --------------------------

stdX <- function(X) {
  if(length(X[1,]) == 1) return(X)
	X.std <- X
	for(i in 2:length(X[1,])) {
		X.std.means <- mean(X[,i])
		X.std.stdev <- sqrt(var(X[,i]))
    X.std[,i] <- (X[,i] - X.std.means)/X.std.stdev
	}
	X.std
}

# -------------- PUT FIXED EFFECT ESTIMATES BACK ON ORIGINAL SCALE --------------

unstdX.beta.hat <- function(X, beta.hat, covb) {
  if(length(X[1,]) == 1) return(list(beta.hat = beta.hat, covb = covb))
	bhat.se <- beta.hat*0
	B <- matrix(0, ncol = length(X[1,]), nrow =length(X[1,]))
	B[1,1] <- 1
	for(i in 2:length(X[1,])) {
		X.std.means <- mean(X[,i])
		X.std.stdev <- sqrt(var(X[,i]))
		B[i,i] <- 1/X.std.stdev
		B[1,i] <- -X.std.means/X.std.stdev
	}
	beta.hat1 <- B %*% beta.hat
	covb1 <- B %*% covb %*% t(B)
	junk <- 1
	list(beta.hat = beta.hat1, covb = covb1)
}

# -------- COMPARES MODELS BASED ON INFORMATION CRITERIA AND CROSS-VALIDATION

Info.crit.compare <- function(model.list)
{
   IC <- NULL
   for(i in 1:length(model.list)) {
     if(is.null(model.list[[i]]$CorModel1)) model1 <- NA
     else model1 <- model.list[[i]]$CorModel1
     if(is.null(model.list[[i]]$CorModel2)) model2 <- NA
     else model2 <- model.list[[i]]$CorModel2
     if(is.null(model.list[[i]]$CorModel3)) model3 <- NA
     else model3 <- model.list[[i]]$CorModel3
     IC <- rbind(IC, data.frame(
       formula = deparse(model.list[[i]]$formula, width.cutoff = 500),
       EstMethod = model.list[[i]]$EstMeth,
       model1 = model1,
       model2 = model2,
       model3 = model3,
       neg2LogL = model.list[[i]]$m2LL,
       AIC = model.list[[i]]$AIC,
       AICc = model.list[[i]]$AICc,
       BIC = model.list[[i]]$BIC,
       model.list[[i]]$cv.stats)
       )
   }
   IC
}

#------------ REML EQUATION TO MINIMIZE

m2LL.stream <- function(theta, m2LLdata, X, dist.hydro = NULL,
  a.mat = NULL, b.mat = NULL, net.zero = NULL,
	weight = NULL, x.dat = NULL, y.dat = NULL,
	Del.i = NULL, A.5 = NULL,
  CorModel1, CorModel2, CorModel3, use.nugget, use.anisotropy,
  EstMeth = "REML")
{
    if((max(theta) > 20) | (min(theta) < -20)) return(1e+32)
    theta1 <- untrans.theta(theta = theta,
      CorModel1 = CorModel1, CorModel2 = CorModel2, CorModel3 = CorModel3,
      use.nugget = use.nugget, use.anisotropy = use.anisotropy)
  	z <- m2LLdata
		n <- length(z)
		p <- length(X[1,])
    V <- make.cov.mat(theta1, dist.hydro = dist.hydro, w.matrix = weight,
      a.mat = a.mat, b.mat = b.mat, net.zero = net.zero,
      xcoord = x.dat, ycoord = y.dat,
      CorModel1 = CorModel1, CorModel2 = CorModel2, CorModel3 = CorModel3,
      use.nugget = use.nugget, use.anisotropy = use.anisotropy)
    if(!is.null(Del.i)) V <- diag(Del.i) %*% diag(A.5) %*% V %*%
      diag(A.5) %*% diag(Del.i)
    ViX <- solve(V,X)
		covbi <- t(X) %*% ViX
		covb <- mginv(covbi)
		b.hat <- covb %*% t(ViX) %*% z
		r <- z - X %*% b.hat
		f1 <- t(r) %*% solve(V,r) + sum(log(svd(V)$d))
		if(EstMeth == "REML") f1 <- f1 + sum(log(svd(covbi)$d))
		nmult <- (n - p)*log(2*pi)
		if(EstMeth == "ML") nmult <- n*log(2*pi)
		f1 + nmult
}


#------------ MASTER FUNCTION FOR SPATIAL LINEAR MODELS FOR STREAM NETWORKS

slm.stream <- function(formula, data,
  dist.junc = NULL,
  flow.matrix = NULL,
  net.zero = NULL,
  xcol = NULL,
  ycol = NULL,
  Distribution = "Gaussian",
  CorModel1 = "Exponential.tailup",
  CorModel2 = "LinearSill.taildown",
  CorModel3 = NULL,
  use.nugget = TRUE,
  use.anisotropy = FALSE,
  trialscol = NULL,
	EstMeth = "REML",
	trans.power = NULL,
  trans.shift = 0,
  prediction.indices = NULL)
{
  Err <- arg.error.check(formula = formula, data = data,
  dist.junc = dist.junc,  flow.matrix = flow.matrix,
  xcol = xcol, ycol = ycol,
  CorModel1 = CorModel1, CorModel2 = CorModel2, CorModel3 = CorModel3,
  use.nugget = use.nugget, use.anisotropy = use.anisotropy,
  Distribution = Distribution,
	EstMeth = EstMeth,
  prediction.indices = prediction.indices)
  
  if(Err$Err == 1) return(print(Err$message))
  
  # get a list of response and covariate names
  mod.names <- as.character(attr(terms(formula, data = data),"variables"))
  # get the number of names ( + 1, as the first is always "list")
  nc.tmp <- length(mod.names)
  # get the number of rows in the data set
  nr.tmp <- length(data[,1])
  # create a vector of all TRUE values
  ind.tmp <- rep(TRUE, times = nr.tmp)
  # if there are any covariates ...
  if(nc.tmp > 2) {
    # create a FALSE for a record with missing values of the covariates
    for(i in 3:nc.tmp) ind.tmp <- ind.tmp & !is.na(data[,mod.names[i]])
  }
  # remove records that had any missing values for any of the covariates
  data <- data[ind.tmp,]
  # name of the response variable
	response.col <- mod.names[2]
	# create a working column of the response variable
	data[,"work.col"] <- data[,response.col]
	# total number of observations
	n.all <- length(data[,1])
  # create a vector of all TRUE values
	ind <- rep(TRUE, times = n.all)
	# transformations on the working column
	if(!is.null(trans.power)) {
    if(Distribution != "Gaussian")
      return(print("Power transformation can only be used with Gaussian Distribution"))
    if(trans.power < 0) return(print("Power transformation must be > 0"))
    if(trans.power > .5) return(print("Power transformation must be < 0.5"))
    if(trans.power == 0)
      data[,"work.col"] <- log(data[,response.col] + trans.shift)
    else
      data[,"work.col"] <- (data[,response.col] + trans.shift)^trans.power
  }
	# set number of missing response values as 0
	n.na <- 0
	# create indicator vectors of observed and missing response values
	if(any(is.na(data[,"work.col"]))) {
		ind.na <- is.na(data[,"work.col"])
		ind <- !ind.na
		data.na <- data[ind.na,]
		n.na <- sum(ind.na)
		data.na[,response.col] <- 1
#		m.na <- model.frame(formula, data.na)
#		X.na1 <- model.matrix(formula, m.na)
#		X.na <- stdX(X.na1)
	}
	# grab all values of the response column
	z <- data[,response.col]
	# replace response column with all 1's for design matrix of all records
	data[,response.col] <- rep(1, times = n.all)
	# design matrix of all records
	Xdesign1 <- X.design(formula = formula,
      data = data, Xmethod = "treat.first.0")
  # standardize design matrix of all records
  Xdesign.all <- as.matrix(stdX(Xdesign1))
  # put response variable back in data set
	data[,response.col] <- z
	# data set of observed data only
  data1 <- data[ind,]
	if(is.factor(data1[,"work.col"]))
		data1[,"work.col"] <- as.numeric(as.character(data1[,"work.col"]))
	if(is.character(data1[,"work.col"]))
		data1[,"work.col"] <- as.numeric(data1[,"work.col"])
  # vector of observed response variable
	z <- data1[,"work.col"]
	n.obs <- length(z)
	# if response variable is binomial with n trials change z to proportion
	if(!is.null(trialscol) & Distribution == "Binomial"){
	  if(is.factor(data1[,"trialscol"]))
		  data1[,"trialscol"] <- as.numeric(as.character(data1[,"trialscol"]))
	  if(is.character(data1[,"trialscol"]))
		  data1[,"trialscol"] <- as.numeric(data1[,"trialscol"])
		trialsvec <- data1[ind,"trialscol"]
		z <- z/trialsvec
	}
	# else if response is Bernoulli, set trialsvec to all ones
	if(is.null(trialscol) & Distribution == "Binomial")
    trialsvec <- rep(1, times = n.obs)

	# design matrix of observed response records
  Xdesign <- as.matrix(Xdesign.all[ind,])
  # if any missing response values, design matrix of missing records
  if(n.na > 0) X.na <- matrix(Xdesign.all[ind.na,], nrow = n.na)

	# number of observed locations
	n <- length(z)
	# get the rank of the design matrix
	p <- sum(svd(Xdesign)$d>1e-10)

  a.mat <- NULL
  b.mat <- NULL
  a.mat.data <- NULL
  b.mat.data <- NULL
  dist.hydro <- NULL
  dist.hydro.data <- NULL
	# create any necessary matrices from distance and flow matrices
	if(!is.null(dist.junc) ) {
	    # maximum distance to common junction between two sites
			a.mat <- pmax(dist.junc,t(dist.junc))
			a.mat.data <- a.mat[ind,ind]
	    # minimum distance to common junction between two sites
			b.mat <- pmin(dist.junc,t(dist.junc))
			b.mat.data <- b.mat[ind,ind]
	    # hydrological distance
	    dist.hydro <- as.matrix(dist.junc + t(dist.junc))
	    # subset stream distance to observed locations only
	    dist.hydro.data <- dist.hydro[ind, ind]
	}
  w.matrix <- NULL
  w.matrix.data <- NULL
  if(max(flow.matrix) > 1) return(print("Flow matrix has values > 1"))
	if(!is.null(flow.matrix) ) {
			flow.sym.matrix <- flow.matrix + t(flow.matrix)
      if(max(flow.sym.matrix) > 1)
        return(print("Flow matrix has symmetry problems"))
			diag(flow.sym.matrix) <- 1
			# weight matrix for kernels if they split at a node
			w.matrix <- sqrt(flow.sym.matrix)
			# subset w.matrix to observed locations only
			w.matrix.data <- w.matrix[ind, ind]
	}
  net.zero.data <- net.zero[ind,ind]
  
  xcoord <- NULL
  ycoord <- NULL
  xcoord.data <- NULL
  ycoord.data <- NULL
  # create Euclidean distance matrix among observed data
	if(!is.null(xcol)){
    dist.Euclid.data <- dist.geo.iso(data[ind,xcol], data[ind,ycol], 1)
    xcoord <- data[,xcol]
    ycoord <- data[,ycol]
    xcoord.data <- xcoord[ind]
    ycoord.data <- ycoord[ind]
  }
  
	# number of CorModels
	n.models <- 1 + sum(!is.null(CorModel2)) + sum(!is.null(CorModel3))

  # Initial parameter estimates, fixed effects

	if(Distribution == "Binomial") {
	  beta.hat <- rep(0, times = length(Xdesign[1,]))
  	beta.hat[1] <- log((sum(z)/sum(trialsvec))/
		  (1 - sum(z)/sum(trialsvec)))
  }
	if(Distribution == "Poisson") {
	  beta.hat <- rep(0, times = length(Xdesign[1,]))
  	beta.hat[1] <- log(sum(z))
  }

	# ----------- START LOOPING HERE ---------------------------

  # create an indicator to stop looping
	stoploop <- 0
	# keep track of the number of iterations
	iter <- 0
	# keep track of number of inner iterations for beta
	inner.iter2 <- NULL
	# begin looping
	while(stoploop == 0) {
		if(Distribution == "Binomial") {
      beta.current <- beta.hat
		  eta.hat <- Xdesign %*% beta.hat
      #diagonal elements of Delta~^{-1} of my manuscript
		  Del.i <- as.vector((1 + exp(eta.hat))^2/exp(eta.hat))
      #diagonal elements of A^(1/2) of my manuscript
		  A.5 <- as.vector(sqrt(exp(eta.hat)/(1 +
        exp(eta.hat))^2/trialsvec))
		  #Binomial pseudo data
 		  zt <- Del.i*(z - exp(eta.hat)/(1 + exp(eta.hat))) + eta.hat
    }
		if(Distribution == "Poisson") {
      beta.current <- beta.hat
		  eta.hat <- Xdesign %*% beta.hat
      #diagonal elements of Delta~^{-1} of my manuscript
		  Del.i <- as.vector(1/exp(eta.hat))
      #diagonal elements of A^(1/2) of my manuscript
		  A.5 <- as.vector(sqrt(exp(eta.hat)))
		  #Poisson pseudo data
 		  zt <- Del.i*(z - exp(eta.hat)) + eta.hat
    }
    if(Distribution == "Gaussian"){
      A.5 <-  NULL
      Del.i <-  NULL
      zt <- z
    }
	  # initial parameter estimates
    theta <- theta.ini(z = zt, X = Xdesign,
      CorModel1 = CorModel1, CorModel2 = CorModel2, CorModel3 = CorModel3,
      use.nugget = use.nugget,  use.anisotropy = use.anisotropy,
      dist.hydro.data = dist.hydro.data, dist.Euclid.data = dist.Euclid.data,
      n.models = n.models)

	  # covariance parameter estimates using ML or REML
	  parmest1.out <- optim(theta, m2LL.stream, m2LLdata = zt, X = Xdesign,
		  dist.hydro = dist.hydro.data, weight = w.matrix.data,
		  net.zero = net.zero.data,
		  a.mat = a.mat.data, b.mat = b.mat.data,
		  x.dat = xcoord.data, y.dat = ycoord.data,
		  Del.i = Del.i, A.5 = A.5,
		  CorModel1 = CorModel1, CorModel2 = CorModel2, CorModel3 = CorModel3,
      use.nugget = use.nugget, use.anisotropy = use.anisotropy,
      EstMeth = EstMeth,
      method = "Nelder-Mead")
	  # covariance parameter estimates using ML or REML
	  parmest2.out <- optim(theta, m2LL.stream, m2LLdata = zt, X = Xdesign,
		  dist.hydro = dist.hydro.data, weight = w.matrix.data,
		  net.zero = net.zero.data,
		  a.mat = a.mat.data, b.mat = b.mat.data,
		  x.dat = xcoord.data, y.dat = ycoord.data,
		  Del.i = Del.i, A.5 = A.5,
		  CorModel1 = CorModel1, CorModel2 = CorModel2, CorModel3 = CorModel3,
      use.nugget = use.nugget, use.anisotropy = use.anisotropy,
      EstMeth = EstMeth,
      method = "BFGS")
    if(parmest1.out$value < parmest2.out$value) parmest.out <- parmest1.out
    else parmest.out <- parmest2.out
	  parmest <- parmest.out$par
	  
	  # go back to original scale for covariance parameters
    parmest <- untrans.theta(theta = parmest,
      CorModel1 = CorModel1, CorModel2 = CorModel2, CorModel3 = CorModel3,
      use.nugget = use.nugget, use.anisotropy = use.anisotropy)

	  V <- make.cov.mat(parmest, dist.hydro = dist.hydro.data,
      w.matrix = w.matrix.data, net.zero = net.zero.data,
      a.mat = a.mat.data, b.mat = b.mat.data,
      xcoord = xcoord.data, ycoord = ycoord.data,
      CorModel1 = CorModel1, CorModel2 = CorModel2, CorModel3 = CorModel3,
      use.nugget = use.nugget, use.anisotropy = use.anisotropy)
      
    if(!is.null(Del.i)) V <- diag(Del.i) %*% diag(A.5) %*% V %*%
      diag(A.5) %*% diag(Del.i)
    ViX <- solve(V,Xdesign)
		covbi <- t(Xdesign) %*% ViX
		covb <- mginv(covbi)
		if(Distribution != "Gaussian") beta.old <- beta.hat
		beta.hat <- covb %*% t(ViX) %*% zt
		if(Distribution == "Gaussian") beta.old <- beta.hat

    #convergence criteria on the fixed effect parameters
		if(all(abs(beta.hat - beta.old)/beta.old < 1e-5))
      stoploop <- 1
		if (iter >= 20) stoploop <- 1

		iter <- iter + 1
		
  }

	# ----------- DONE LOOPING HERE ---------------------------

	# covariance matrix between all locations (observed and prediction)
  sigma <- make.cov.mat(parmest, dist.hydro = dist.hydro, w.matrix = w.matrix,
      a.mat = a.mat, b.mat = b.mat, net.zero = net.zero,
      xcoord = xcoord, ycoord = ycoord,
      CorModel1 = CorModel1, CorModel2 = CorModel2, CorModel3 = CorModel3,
      use.nugget = use.nugget, use.anisotropy = use.anisotropy)

	# inverse covariance matrix between observed locations
	Vi <- solve(V)

	covbi <- t(Xdesign) %*% Vi %*% Xdesign
	covb <- solve(covbi)
	bhat.se <- sqrt(diag(covb))
	b.hat <- covb %*% t(Xdesign) %*% Vi %*% zt
  if(Distribution != "Gaussian") R2 <- NULL
  if(Distribution == "Gaussian") {
     R2 <- R2g(z = zt, Xdesign = Xdesign, bhat = b.hat, Vi = Vi)
  }
	bk.out <- NULL
  if(Distribution == "Gaussian") {
    if(any(is.na(data[,"work.col"]))) {
      #krige all missing values
      pk <- pred.k(z = zt, Xdesign = Xdesign, X.na = X.na, ind = ind, covb = covb,
        b.hat = b.hat, Vi = Vi, sigma = sigma, trans.power = trans.power,
        trans.shift = trans.shift)

      data[,"predictions"] <- rep(NA, times = n + n.na)
      data[,"pred.std.err"] <- rep(NA, times = n+ n.na)
      data[ind.na,"predictions"] <- pk[,"predictions"]
      data[ind.na,"pred.std.err"] <- pk[,"pred.std.err"]
    }

	  if(!is.null(prediction.indices)) {
	  	data.j <- data
	  	data.j[,response.col] <- 1
	  	m.all <- model.frame(formula, data.j)
	  	X.all <- stdX(model.matrix(formula, m.all))
      bk.out <- block.krige(z = z, pred.indices = prediction.indices,
	      covb = covb, beta.hat = beta.hat, V = V, Vi = Vi, Sigma = sigma,
        X.all = X.all, ind.na = ind.na)
	  }
  }

	# crossvalidation
	cv.out <- crossval(z = zt, X = Xdesign, V = V, Vi = Vi, n = n, p = p)
	cv.stats <- data.frame(bias = sum(cv.out[,1] - zt)/n,
		std.bias = sum((cv.out[,1] - zt)/sqrt(cv.out[,2]))/n,
		RMSPE = sqrt(sum((zt - cv.out[,1])^2)/n),
		RAV = sqrt(sum(cv.out[,2]^2)/n),
		std.MSPE = sum((cv.out[,1] - zt)^2/cv.out[,2]^2)/n,
		cov.80 = sum(abs((zt - cv.out[,1])/cv.out[,2]) < 1.281552)/n,
		cov.90 = sum(abs((zt - cv.out[,1])/cv.out[,2]) < 1.644854)/n,
		cov.95 = sum(abs((zt - cv.out[,1])/cv.out[,2]) < 1.959964)/n)
	infld <- infldiagn(z = zt, V = V, X = Xdesign, n = n, p = p)
	fit <- data.frame(data1, fit = Xdesign %*% b.hat,
		resid = zt - Xdesign %*% b.hat, infld, cv.out, cv.resid = zt - cv.out[,1],
		cv.stndr = (zt - cv.out[,1])/cv.out[,2])
	ubh <- unstdX.beta.hat(Xdesign1, b.hat, covb)
	b.hat <- ubh$beta.hat
  bhat.se <- sqrt(diag(ubh$covb))
	fixed.eff.est <- data.frame(parameters = b.hat,
			std.err = bhat.se,
			df = n - p,
			t.value = b.hat/bhat.se,
			prob.t = round(100000*(1 - pt(abs(b.hat/bhat.se), df = n - p))*2)/100000
			)
  data2 <- data1
  data2[,"Value"] <- zt
  data2[,"work.col"] <- zt
	fixed.effects.estimates <- fixed.effect.table(formula = formula, data = data2,
			fixed.effects.estimates = fixed.eff.est)
	typeIII.hypoth <- hypoth.stream(formula, data2, Vi)
	ind1 <- !is.na(fixed.effects.estimates[,4])
	cov.fix.eff <- matrix(0, nrow = length(fixed.effects.estimates[,4]),
		ncol = length(fixed.effects.estimates[,4]))
	cov.fix.eff[ind1,ind1] <- covb

	m2LL <- parmest.out$value
	if(EstMeth == "REML") {
		k <- length(parmest.out$par)
		n1 <- n - p
	}
	if(EstMeth == "ML") {
		k <- p + length(parmest.out$par)
		n1 <- n
	}
	AIC <- m2LL + 2*k
	AICc <- m2LL + 2*k*n1/(n1 - k - 1)
	BIC <- m2LL + log(n1)*k

  npar.sofar <- 0
  Cov.par1 <- NULL
  if(!is.null(CorModel1)) {
  if(substr(unlist(strsplit(CorModel1,"\\."))[2],1,1) == "t")  {
   Cov.par1 <- list(parsil1 = parmest[1], range1 = parmest[2])
   npar.sofar <- 2
  }
  else
    if(use.anisotropy == FALSE) {
       Cov.par1 <- list(parsil = parmest[1], range1 = parmest[2])
       npar.sofar <- 2
    }
    else {
       Cov.par1 <- list(parsil = parmest[1], range1 = parmest[2],
           minorp = parmest[3], rotate = parmest[4])
       npar.sofar <- 4
    }
  }
  Cov.par2 <- NULL
  if(!is.null(CorModel2)) {
  if(substr(unlist(strsplit(CorModel2,"\\."))[2],1,1) == "t")  {
   Cov.par2 <- list(parsil = parmest[npar.sofar + 1],
       range1 = parmest[npar.sofar + 2])
   npar.sofar <- npar.sofar + 2
  }
  else
    if(use.anisotropy == FALSE) {
       Cov.par2 <- list(parsil = parmest[npar.sofar + 1],
           range1 = parmest[npar.sofar + 2])
       npar.sofar <- npar.sofar + 2
    }
    else {
       Cov.par2 <- list(parsil = parmest[npar.sofar + 1],
           range1 = parmest[npar.sofar + 2], minorp = parmest[npar.sofar + 3],
           rotate = parmest[npar.sofar + 4])
       npar.sofar <- npar.sofar + 4
    }
  }
  Cov.par3 <- NULL
  if(!is.null(CorModel3)) {
  if(substr(unlist(strsplit(CorModel3,"\\."))[2],1,1) == "t")  {
   Cov.par3 <- list(parsil = parmest[npar.sofar + 1],
       range1 = parmest[npar.sofar + 2])
   npar.sofar <- npar.sofar + 2
  }
  else
    if(use.anisotropy == FALSE) {
       Cov.par3 <- list(parsil = parmest[npar.sofar + 1],
           range1 = parmest[npar.sofar + 2])
       npar.sofar <- npar.sofar + 2
    }
    else {
       Cov.par3 <- list(parsil = parmest[npar.sofar + 1],
           range1 = parmest[npar.sofar + 2], minorp = parmest[npar.sofar + 3],
           rotate = parmest[npar.sofar + 4])
       npar.sofar <- npar.sofar + 4
    }
  }
  nugget <- NULL
  if(use.nugget == TRUE)  nugget <- parmest[npar.sofar + 1]
  
	list(
    formula = formula,
		sample.size = n + n.na,
		obs.sample.size = n,
		missing.sample.size = n.na,
		Distribution = Distribution,
		CorModel1 = CorModel1,
		Covariance.Parameters1 = Cov.par1,
		CorModel2 = CorModel2,
		Covariance.Parameters2 = Cov.par2,
		CorModel3 = CorModel3,
		Covariance.Parameters3 = Cov.par3,
    nugget = nugget,
		EstMeth = EstMeth,
		fixed.effects.estimates = fixed.effects.estimates,
		typeIII.hypoth = typeIII.hypoth,
		R2g = R2,
		cov.fix.eff = cov.fix.eff,
		m2LL =	m2LL,
		AIC = AIC,
		AICc = AICc,
		BIC = BIC,
		cv.stats = cv.stats,
		data = data,
		fit = fit,
		block.krige = bk.out,
		V = V,
		Vi = Vi,
		X = Xdesign,
		min.eig.values = min(svd(sigma)$d)
	)
}


