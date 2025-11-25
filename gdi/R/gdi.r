##Function sellipse()
#'Estimate area of a superellipse. Assistant function for gdi.
#'
#' @param a First radius of the superellipse.
#' @param b Second radius of the superellipse.
#' @param k superellipse exponent.
#' @return A single number giving the area of the superellipse
#' @export sellipse
#' @examples
#' major_radius<-2
#' minor_radius<-3
#' exponent<-2.3
#' sellipse(major_radius, minor_radius, exponent)


sellipse <- function(a, b, k) {
  gamma1 <- gamma(1+1/k)
  gamma2 <- gamma(1+2/k)
  
  # Calculate area using formula
  area <- 4 * a * b * gamma1^2 / gamma2
  
  # Return the calculated area
  return(area)
}

##Function sellipse.coo()
#'calculate coordinates for plotting a superellipse for visualizing body cross-sections
#'
#' @param k superellipse exponent.
#' @param res the desired resolution
#' @return a data frame containing x and y coordinates for the outline of the (super)ellipse
#' @export sellipse.coo
#' @examples
#' sellipse.coo(2.0)->df #get coordinates for normal ellipse (exponent k=2)
#' plot(df$x,df$y,col="black", type="l") #plot normal ellipse
#' sellipse.coo(2.3)->df2 # get coordinates for superellipse with exponent 2.3
#' lines(df$x,df$y, col="blue") #plot superellipse

sellipse.coo <- function(k, res=100) {
  t <- seq(0, 2 * pi, length.out = res)
  x <- abs(cos(t))^(2 / k) * sign(cos(t))
  y <- abs(sin(t))^(2 / k) * sign(sin(t))
  
  df <- data.frame(x = x, y = y)
  return(df)
}


##Function cscorr()
#' Measure and analyze cross-sectional geometry supplied as an image.
#'
#' @param image_file Image to be read. Images can be jpeg or png files, or a previously read image saved as an array/matrix-type object in R.
#' @param threshold Reference value for color criterium after which pixels that are part of the cross-section are differentiated from the background.
#' @param channel Color channel to which to apply the threshold criterium. Default is 4 (alpha channel of rgba image). Channel setting needs to be adjusted depending on the color mode of the image used (e.g. there are two channels to choose from in a greyscale image, and 3 in an rgb image).
#' @param method Method for determining which pixels to count. Default "greater" counts pixels with value greater than threshold (e.g. higher opacity, in the case of an alpha channel). "less" counts pixels with a value less than the threshold. "not" counts all pixels not precisely matching threshold. Any other character string results in only pixels exactly matching the value given as threshold being counted.
#' @param return What value to return. Possible values are "area_corr" (Default, returns ratio between measured area and area of ellipse with same horizontal and vertical diameters), "aspect_ratio" (returns aspect ratio), "diameters" (returns diameters), "area" (returns area) and "rotI" (returns correction factors for rotational inertia calculations). Any other value for this parameter will prompt the function to return a vector containing all of these outputs.
#' @param k optional superellipse exponent for the (super)ellipse to which the measurements should be compared (for the "area_corr" setting for the parameter return).
#' @param scale Optional scale of the image (for raw area measurements).
#' @return Either a numeric of length 1 (depending on the input of the return parameter), defaulting to the area correction factor (if return=="area_corr"), or (if return is left empty or does not match any of the predefined settings) a numeric vector with 8 elements, containing all the possible outputs (x and y diameters, aspect ratio, area and area correction factor, correction factors representing ratios of rotational inertia in x, y and z planes relative an ellipse of equal diameters).
#' @import jpeg
#' @import png
#' @export cscorr
#' @examples
#' fdir <- system.file(package="gdi")
#' correction_factor <- cscorr(file.path(fdir,"exdata","cross_section.png"))


cscorr <- function(image_file, threshold=0.5, channel=4, method="greater", return="area_corr", k=2.0, scale=1) {
#load and save image data to variable img
if(!is.character(image_file)){
img<-image_file
}else if(grepl(".jpg",image_file)==TRUE | grepl(".jpeg",image_file)==TRUE | grepl(".JPG",image_file)==TRUE){
img <- jpeg::readJPEG(image_file)
#read image if it is jpg
}else if(grepl(".png",image_file)==TRUE | grepl(".PNG",image_file)==TRUE){
img <- png::readPNG(image_file)}#read image if it is png



# measure cross-section depth
nrows <- dim(img)[1]
ncols <- dim(img)[2]

if(!is.na(dim(img)[3])#if image has more than one channel, i.e. is an array, reduce the image to the channel to be analyzed
){
img<-img[,,channel]
}

##find maximum depth of cross-section
d <- numeric(0)
#loop through vertical lines of pixels, find all rows in which there are pixels that are part of the cross-section
for(x in 1:ncols){
    if(method=="greater"){
    d <- union(d,which(img[,x]>threshold))
    }else if(method=="less"){
    d <- union(d,which(img[,x]<threshold))
    }else if(method=="not"){
    d <- union(d,which(signif(img[,x],6)!=signif(threshold,6)))
    }else{
    d <- union(d,which(signif(img[,x],6)==signif(threshold,6)))
    }
}



vdiam <- length(d)/scale #calculate maximum depth in pixels

##find maximum width of cross-section
d <- numeric(0)
#loop through horizontal lines of pixels, find all columns in which there are pixels that are part of the cross-section

for(y in 1:nrows){
    if(method=="greater"){
    d <- union(d,which(img[y,]>threshold))
    }else if(method=="less"){
    d <- union(d,which(img[y,]<threshold))
    }else if(method=="not"){
    d <- union( d,which( signif(img[y,],6)!=signif(threshold,6)) )
    }else{
    d <- union( d,which( signif(img[y,],6)==signif(threshold,6)) )
    }
}



hdiam <- length(d)/scale #calculate maximum width in pixels


## calculate (super)elliptical cross-sectional area:
ecomp <- sellipse(vdiam/2, hdiam/2, k)

    ## calculate actual cross-sectional area
    if(method=="greater"){
    area <- sum(img[,]>threshold)
    }else if(method=="less"){
    area <- sum(img[,]<threshold)
    }else if(method=="not"){
    area <- sum( signif(img[,],6)!=signif(threshold,6) )
    }else{
    area <- sum( signif(img[,],6)==signif(threshold,6) )
    }
    area<-area/scale^2#set scale for area
    

    ## Return the calculated area or ratio
  if(return=="area"){return(area)
  }else if(return%in%c("area_corr","corr")){
  return(area/ecomp)
  }else if(return%in%c("asp","aspect_ratio")){
  return(hdiam/vdiam)
  }else if(return=="diameters"){
  diam<-c(x=hdiam, y=vdiam)
  return(diam)
  }else{
  
  
  ##rotational inertia contribution of 1 px², as rectangular plane
Ix<-(1^2)*1/12*1/(area*scale^2)#x direction
Iy<-(1^2)*1/12*1/(area*scale^2)#y direction

moments<-matrix(nrow=length(img), ncol=6)
colnames(moments)<-c("Ix","Iy","xpos", "ypos","Ix_","Iy_")
##loop through each pixel of image and save area moments and coordinates for each pixel belonging to the shape
i<-0#set index to start at 0
    for (x in 1:ncols){#loop through collumns
        for(y in 1:nrows){#loop through rows of collumn x
i<-i+1#increment the index
#select pixels belonging to shape depending on the method and threshold settings
if(method=="greater" & img[y,x]>threshold){
moments[i,1]<-Ix
moments[i,2]<-Iy
moments[i,3]<-x
moments[i,4]<-y
}else if(method=="less" & img[y,x]<threshold){
moments[i,1]<-Ix
moments[i,2]<-Iy
moments[i,3]<-x
moments[i,4]<-y
}else if(method=="not" & img[y,x]!=threshold){
moments[i,1]<-Ix
moments[i,2]<-Iy
moments[i,3]<-x
moments[i,4]<-y
}else if(img[y,x]==threshold){
moments[i,1]<-Ix
moments[i,2]<-Iy
moments[i,3]<-x
moments[i,4]<-y
}

        }
    }#end of loop

xcentroid<-mean(moments[,3], na.rm=T)#calculate centroid position for entire shape
ycentroid<-mean(moments[,4], na.rm=T)

#then use parallel axis theorem to convert pixel moments relative to overall centroid
moments[,5]<-Ix+1/(area*scale^2)*(moments[,4]-ycentroid)^2
moments[,6]<-Iy+1/(area*scale^2)*(moments[,3]-xcentroid)^2

#sum up all individual pixel moments for shape:
I_x_total<-sum(moments[,5],na.rm=T)
I_y_total<-sum(moments[,6],na.rm=T)

#calculate polar moment
I_z_total<-sum(I_x_total, I_y_total)

#rotational inertia for ellipse with same diameters and same mass
ellipse_y<-1/16*(hdiam*scale)^2
ellipse_x<-1/16*(vdiam*scale)^2
ellipse_polar<-ellipse_y+ellipse_x

I_corr<-c(I_x_total,I_y_total,I_z_total)/c(ellipse_x, ellipse_y, ellipse_polar)

  
  if(return=="rotI"){
  names(I_corr)<-c("I_corr_x_pitch","I_corr_y_yaw","I_corr_z_roll")
  return(I_corr)
  }else{
  full<-c(x=hdiam, y=vdiam, asp=hdiam/vdiam, area=area, area_corr=area/ecomp, I_corr_x_pitch=I_corr[1], I_corr_y_yaw=I_corr[2], I_corr_z_roll=I_corr[1])
  return(full)}
    }
}


##Function measuresil()
#'Take pixel-by-pixel measurements of a silhouette in jpeg or png format for use with the gdi function.
#'
#' @param image_file Image to be read. Images can be jpeg or png files, or a previously read image saved as an object in R.
#' @param threshold Reference value for color criterium after which pixels that are part of the silhouette are differentiated from the background.
#' @param channel Color channel to which to apply the threshold criterium. Default is 4 (alpha channel of rgba image). Channel setting needs to be adjusted depending on the color mode of the image used (e.g. there are two channels to choose from in a greyscale image with transparency, and 3 in an rgb image without transparency, or 4 in a full rgba image).
#' @param method Method for determining which pixels to count. Default "greater" counts pixels with value greater than threshold (e.g. higher opacity, in the case of an alpha channel). "less" counts pixels with a value less than the threshold. "not" counts all pixels not precisely matching threshold. Any other character string results in only pixels exactly matching the value given as threshold being counted.
#' @param align Indicate whether the silhouette long axis is aligned horizontally (setting "h", default), or vertically (any other parameter setting).
#' @param return Setting for what to return, default setting ("diameters") returns a single vector containing the diameters, any other setting returns a data frame containing centers and diameters.
#' @import jpeg
#' @import png
#' @export measuresil
#' @return A numeric vector giving the measurements of the silhouette
#' @examples
#' fdir <- system.file(package="gdi")
#' lat <- measuresil(file.path(fdir,"exdata","lat.png"))

measuresil<-function(image_file, threshold=0.5, channel=4, method="greater", align="h", return="diameters"){
#load and save image data to variable named img
if(!is.character(image_file)){
img<-image_file
}else if(grepl(".jpg",image_file)==TRUE | grepl(".jpeg",image_file)==TRUE | grepl(".JPG",image_file)==TRUE){
img <- jpeg::readJPEG(image_file)
#read image if it is jpg
}else if(grepl(".png",image_file)==TRUE | grepl(".PNG",image_file)==TRUE){
img <- png::readPNG(image_file)}#read image if it is png


#get dimensions
nrows <- dim(img)[1]
ncols <- dim(img)[2]

if(!is.na(dim(img)[3])#if image has more than one channel, i.e. is an array, reduce the image to the channel to be analyzed
){
img<-img[,,channel]
}

if(align=="h"){#if horizontally aligned, default
# Loop through each vertical line of pixels
depths <- rep(0, ncols)
centers <- rep(NA, ncols)



for (x in 1:ncols) {
        if(method=="greater"){
        s<-which(img[, x]>threshold)
        }else if(method=="less"){
        s<-which(img[, x]<threshold)
        }else if(method=="not"){
        s<-which(signif(img[, x],6)!=signif(threshold,6))
        }else{
        s<-which(signif(img[, x],6)==signif(threshold,6))
        }
        
        if(return!="diameters"){
        centers[x]<-nrows-mean(s)
        }
        depths[x]<-length(s)
        }

}else{#if vertically aligned
# Loop through each horizontal line of pixels
depths <- rep(0, nrows)
centers <- rep(NA, nrows)

for (y in 1:nrows) {
    if(method=="greater"){
    s<-which(img[y,]>threshold)
    }else if(method=="less"){
    s<-which(img[y,]<threshold)
    }else if(method=="not"){
    s<-which(signif(img[y,],6)!=signif(threshold,6))
    }else{
    s<-which(signif(img[y,],6)==signif(threshold,6))
    }
    
    if(return!="diameters"){
    centers[y]<-1+mean(s)
    }
    depths[y]<-length(s)
    }

}

if(return=="diameters"){
return(depths)
}else{
    depthscenters<-data.frame(diameter=depths,center=centers)
    attr(depthscenters, "x_resolution")<-ifelse(align=="h", ncols,nrows)#total resolution of image along the x axis
    attr(depthscenters, "y_resolution")<-ifelse(align=="h", nrows,ncols)

return(depthscenters)
}

}



## function inter_corr()
#' Helper function for interpolating cross-sectional metrics over the length of a silhouette
#' @param ... cross-sections between which to interpolate. Each object can be either a single numeric each giving the quantity to be interpolated between, a full output of cscorr() with option return="all", or a character vector giving filenames to use with cscorr().
#' @param x indices along silhouette to use as points of interpolation. If NULL, the function will treat indices as indices of the silhouette itself, not x, and attempt to prepend 1 and append sil (or its length) to it, to make the cross-section interpolate toward the default value at the first and last slice.
#' @param sil dataframe or vector (output of measuresil()) or single numeric giving the total number of slices of the silhouette for which to interpolate cross-sections
#' @param indices subset of x to which the elements in ... pertain, if ... does not have the same total number of elements as x
#' @param default default value to fill for elements of x not in indices, if length(x)>length(...)
#' @param default_k default superellipse exponent to use for elements of x not in indices, supersedes setting for default if not NULL and return parameter is set to "area_corr"
#' @param return which element of the cscorr() output to return in interpolated form (defaults to area_corr, can also be "asp"
#' @param v verbosity setting (logical)
#' @details This is a convenience function and wrapper around approx() and cscorr(), which creates a linear interpolation of cross-section geometry over the entire length of a silhouette. Input can in principle be any numeric value (but, is intended to primarily be cross-sectional geometry defined by "area_corr" and "asp" elements returned via cscorr().
#' @importFrom stats approx
#' @export inter_corr
#' @examples
#' fdir <- system.file(package="gdi")
#' cf <- cscorr(file.path(fdir,"exdata","cross_section.png"))
#' lat <- measuresil(file.path(fdir,"exdata","lat.png"))
#' inter_corr(cf,x=c(1,1000,2341), sil=lat, indices=2)

inter_corr<-function(...,x=NULL,sil,indices=NULL,default=1, default_k=2, return="area_corr",v=FALSE){
list(...)->dots
if(is.data.frame(sil)) sil<-nrow(sil)
if(is.numeric(sil) && length(sil)>1) sil<-length(sil)
if(is.null(x) && !is.null(indices)) c(1,indices,sil)->x
if(v) print(sil)
if(v) print(dots)

if(length(x)!=length(dots) && is.null(indices)) stop("x is not same length as supplied cross-sections and no i was provided to subset it")
#extract corrs
corrs<-numeric()

corr_temp<-list()
for(i in 1:length(dots)){

if(is.numeric(dots[[i]]) && length(dots[[i]])==1){
as.numeric(dots[[i]])[1]->corrs[i]

}else if(is.character(dots[[i]])){
if(dots[[i]]%in%names(corr_temp)){
corr_temp[[dots[[i]]]]->c0
}else{
cscorr(dots[[i]],return=return)->c0
c0->corr_temp[[dots[[i]]]]
}
c0->corrs[i]

}else if(is.numeric(dots[[i]]) && length(dots)>1){
dots[[i]]->c0
c0[return]->corrs[i]
}

}
if(v) print(corrs)
##interpolate corr factors
if(length(x)!=length(corrs) && !is.null(indices) && length(indices)==length(corrs)){
if(is.null(default_k) | return!="area_corr"){default->corr_default
}else{sellipse(1,1,default_k)/sellipse(1,1,2)->corr_default}


if(v) print(corr_default)
y_in<-rep(corr_default,length(x))
y_in[indices]<-corrs

approx(x=x, y=y_in, xout=c(1:sil))$y->out

}else{
corrs->y_in
approx(x=x, y=y_in, xout=c(1:sil))$y->out
}

#set information attributes
attr(out,"x")<-x
attr(out,"y_in")<-y_in

attr(out,"return")<-return

#output
return(out)
}##


##Function gdi()
#'Estimate volume using Graphic Double Integration.
#'
#' @param lat Measurements of diameter in lateral view/first of two orthogonal views to be used with the gdi. Can be either a numeric vector, a data.frame (output of measuresil(...,return="all") with a collumn named "diameter", or a text file with diameter measurements to be scanned.
#' @param dors Measurements of diameter in dorsal view/second of two orthogonal views to be used with the gdi. Can be either a numeric vector, a data.frame (output of measuresil(...,return="all") with a collumn named "diameter", or a text file with diameter measurements to be scanned. Must be the same length as lat.
#' @param indices Optional indices specifying a subset of the silhouette measurement vectors to be analyzed. Useful if separate segment calculations are desired.
#' @param scale Scale of the data, given in terms of how many units of the input data (e.g. pixels) are in one side of the desired unit of output volume. Defaults to 10.
#' @param sliceL Length of individual segments to be used in the GDI. Defaults to 1.
#' @param method Method to be used for the GDI. Default "raw" setting calculates each segment as an elliptical cylinder with volume = Area * SliceL. Any other string will result in volume being calculated as an elliptical frustum with base areas based on the measurements of segments i and i+1.
#' @param k Superellipse exponent to be used for the cross-sectional area. Defaults to 2.0 (normal ellipse).
#' @param corr Correction factor for area of cross-sections, calculated as the ratio between the actual cross-sectional area and that of a (super)ellipse (depending on the specified exponent k) with the same diameters. This setting enables the function to account for complex, non-elliptical cross-sections. Default value is 1, i.e. no correction. Can be either a single number, or a numeric vector of the same length as lat and dors (in the case of a changing cross-sectional geometry along the length of the body).
#' @param smooth.ends If method != "raw", specify whether first and last segments should be left raw, or taper to 0 (i.e. be approximated as cones). Only applies if there are no leading or following zeros in the measurement vectors.
#' @param return Determines whether to report the estimated total volume (if default/"total"), or a data.frame() with segment radii, areas and volumes (if left empty of any other character string).
#' @return Either a single number representing the total volume estimated (with names indicating the horizontal length of the silhouette in the unit determined by scale), or (if return!="total") a data.frame() containing columns with the radii in both dimensions, the estimated elliptical or superelliptical areas, and the segment volumes.
#' @export gdi
#' @examples
#' lateral <- rep(2,4) #generate example data
#' dorsal <- rep(2,4)
#' gdi(lat=lateral, dors=dorsal, scale=10, method="raw", k=2.0)
#' gdi(lat=lateral, dors=lateral/2, scale=10, method="smooth", k=2.3)


gdi<-function(lat, dors, indices=NULL, scale=10, sliceL=1, method="raw", k=2.0, corr=1, smooth.ends=FALSE, return="total"){
lat_<-lat
dors_<-dors

if(is.character(lat)){scan(lat)->lat_}
if(is.character(dors)){scan(dors)->dors_}#scan files, if inputs are filenames


if(is.data.frame(lat)){
        if(length(lat$diameter)!=nrow(lat)){
        stop("lat is a data.frame, but collumn \"diameter\" not found!")}
    lat_<-lat$diameter
}

if(is.data.frame(dors)){
        if(length(dors$diameter)!=nrow(dors)){
        stop("dors is a data.frame, but collumn \"diameter\" not found!")}
    dors_<-dors$diameter
}
#make data.frame
sil<-data.frame(ydiam_raw=lat_,zdiam_raw=dors_,ydiam_scaled=lat_/scale,zdiam_scaled=dors_/scale, slice_length=sliceL/scale, sliceL_raw=sliceL)#scale-adjust and save in dataframe
sliceL/scale->sliceL

if(return!="total"){
#save raw segment centroid positions on the x (anteroposterior) axis
##x
sil$x_center<-sil$sliceL_raw/2#half the raw sliceL
l<-cumsum(sil$x_center*2)#cumulative length of all slices
if(length(sil$x_center)>1){
for(i in 2:length(sil$x_center)){#for each segment, subtract half the segment length from the cumulative segment length to get the segment center
sil$x_center[i]<-l[i]-sil$x_center[i]
}}

##y
if(is.data.frame(lat)){
sil$y_center<-lat$center
}
##z
if(is.data.frame(dors)){
sil$z_center<-dors$center
}
}

#make subset based on indices
if(!is.null(indices)){
sil<-sil[indices,]}

#estimate cross-sections as (super)ellipses
sil$A<-sellipse(sil$ydiam_scaled/2, sil$zdiam_scaled/2, k)

if(method=="raw"){
    sil$V<-sil$A*sliceL#raw formula, approximating each segment as a cylinder
}else{
    for(i in 1:(nrow(sil)-1)){
    sil$V[i]<-(sil$A[i]+sil$A[i+1]+sqrt(sil$A[i]*sil$A[i+1]))/3*ifelse(length(sliceL)>1,sliceL[1],sliceL)
    #approximating segments as frusta, with A[i] and A[i+1] as base areas
    }
    if(smooth.ends==TRUE){
    sil$V[1]<-ifelse(length(sliceL)>1,sliceL[1],sliceL)*sil$A[1]/3#conical first segment (tapers to 0)
    sil$V[nrow(sil)]<-ifelse(length(sliceL)>1,sliceL[nrow(sil)],sliceL)*sil$A[nrow(sil)]/3#conical last segment (tapers to 0)
    }else{
    sil$V[nrow(sil)]<-ifelse(length(sliceL)>1,sliceL[nrow(sil)],sliceL)*sil$A[nrow(sil)]
    }
}

sil$V <- sil$V*corr
attr(sil, "x_resolution")<-attr(lat, "x_resolution")
attr(sil, "y_resolution")<-attr(lat, "y_resolution")
attr(sil, "z_resolution")<-attr(dors, "y_resolution")


sum(sil$V)->res#sum up segments
names(res) <- paste("x_dim",sum(sil$slice_length[sil$ydiam_raw!=0]), "units", sep="_")

if(return=="total"){
return(res)
}else{
return(sil)}
}



##Function fdetect()
#'Tool to help determine which threshold value and method to use with measuresil() or cscorr(). The function analyzes all pixels along the edges of the image to determine the background color, to help with deciding on appropriate settings and avoid errors introduced by inappropriate settings
#'
#' @param image_file Image to be read. Images can be jpeg or png files, or a previously read image saved as an object in R.
#' @param threshold Reference value for color criterium after which pixels that are part of the silhouette should be differentiated from the background.
#' @param channel Color channel to which to apply the threshold criterium. Default is 4 (alpha channel of rgba image). Channel setting needs to be adjusted depending on the color mode of the image used (e.g. there are two channels to choose from in a greyscale image, and 3 in an rgb image).
#' @param plot Whether to plot a histogram with the detected color values (if TRUE) or not (if FALSE, default).
#' @return A list()-object containing: $edgetable (a table of the different color values detected and their respective frequencies), $histogram (a histogram-object of the color values), $most_common (the most common color value found), $foreground (a character string, indicating whether the foreground color value is likely "greater" or "less" than the specified threshold), $result (a character string giving a summary of the results)
#' @export fdetect
#' @importFrom graphics hist
#' @importFrom graphics rug
#' @examples
#' fdir <- system.file(package="gdi")
#' fdetect(file.path(fdir,"exdata","lat.png"))




fdetect<-function(image_file, threshold=0.5, channel=4, plot=FALSE){
if(!is.character(image_file)){
img<-image_file
}else if(grepl(".jpg",image_file)==TRUE | grepl(".jpeg",image_file)==TRUE | grepl(".JPG",image_file)==TRUE){
img <- jpeg::readJPEG(image_file)
#read image if it is jpg
}else if(grepl(".png",image_file)==TRUE | grepl(".PNG",image_file)==TRUE){
img <- png::readPNG(image_file)}#read image if it is png

#get dimensions
nrows <- dim(img)[1]
ncols <- dim(img)[2]

if(!is.na(dim(img)[3])#if image has more than one channel, i.e. is an array, reduce the image to the channel to be analyzed
){
img<-img[,,channel]
}

#extract picture edges:
edges <- c(img[1,],img[nrows,],img[,1],img[,ncols])

out <- list()
out$edgetable <- table(edges)
out$histogram <- hist(edges, breaks=seq(0,1,0.05), plot=FALSE)

out$most_common <- signif(as.numeric(names(which.max(out$edgetable))),7)


#guess correct setting for method parameter in gdi()
if(out$most_common<threshold){
out$foreground <- "greater"
}else if(out$most_common>threshold){
out$foreground <- "less"
}else{
out$foreground <- NA
warning("Automatic background detection failed, please manually ascertain which color value corresponds to the foreground of your silhouette!")
}

#show histogram
if(plot==TRUE){
hist(edges, breaks=seq(0,1,0.05), xlab=paste("Edge pixel color values in channel", channel))
rug(edges)
}
out$result <- paste("Best guess is that the foreground has a channel", channel, "value", out$foreground, "than", threshold, "(background appears to be", out$most_common,")")
return(out)

} 




##Function imghist()
#'Simple histogram analysis for all color values in an input image. Can be used to help assess whether a chosen threshold value is appropriate for differentiating the silhouette from the background, or for general image analysis purposes.
#'
#' @param image_file Image to be read. Images can be jpeg or png files, or a previously read image saved as an object in R.
#' @param threshold Reference value for color criterium after which pixels that are part of the silhouette should be differentiated from the background.
#' @param channel color channel to which to apply the threshold criterium. Default is 4 (alpha channel of rgba image). Channel setting needs to be adjusted depending on the color mode of the image used (e.g. there are two channels to choose from in a greyscale image, and 3 in an rgb image).
#' @param breaks A vector of breaks for the histogram, defaults to a bin width of 0.05 between color values of 0 and 1.
#' @param plot Whether to plot a histogram, defaults to TRUE
#' @param unique Whether to return counts for unique color values, defaults to FALSE.
#' @return A plotted histogram (unless plot==FALSE), and a matrix containing the counts from the histogram (default) or the counts for unique color values (if unique==TRUE).
#' @export imghist
#' @importFrom graphics abline
#' @importFrom graphics hist
#' @importFrom graphics rug
#' @importFrom graphics lines
#' @importFrom stats density
#' @examples
#' fdir <- system.file(package="gdi")
#' imghist(file.path(fdir,"exdata","lat.png"))


imghist <- function(image_file, threshold=0.5, channel=4, breaks=seq(0,1,0.05), plot=TRUE, unique=FALSE){
if(!is.character(image_file)){
img<-image_file
}else if(grepl(".jpg",image_file)==TRUE | grepl(".jpeg",image_file)==TRUE | grepl(".JPG",image_file)==TRUE){
img <- jpeg::readJPEG(image_file)
#read image if it is jpg
}else if(grepl(".png",image_file)==TRUE | grepl(".PNG",image_file)==TRUE){
img <- png::readPNG(image_file)}#read image if it is png


if(!is.na(dim(img)[3])#if image has more than one channel, i.e. is an array, reduce the image to the channel to be analyzed
){
img<-img[,,channel]
}

if(plot==TRUE){#show histogram
hist(img[,], prob=TRUE, breaks=breaks, xlab=paste("All pixel color values in channel", channel))
lines(density(img[,]), col="grey", lwd=2)
abline(v=threshold)}

total <- length(img[,])

if(unique==FALSE){
h <- hist(img[,], breaks=breaks, plot=FALSE)
names<-h$mids
h <- h$counts
h <- cbind(h,h/total)
rownames(h) <- names
colnames(h) <- c("count","proportion")
return(as.data.frame(h))
}

if(unique==TRUE){
counts <- table(img[,])
sorted <- sort(counts, decreasing=TRUE)
names <- names(sorted)
sorted <- cbind(sorted, sorted/total)
rownames(sorted)<-names
colnames(sorted) <- c("count","proportion")
return(as.data.frame(sorted))
}

}


##Function hCOM()
#'Finds the horizontal (x axis, i.e. the axis vertical to the cross-sections) position of the center of mass (COM) of the volume. Experimental; only valid for "raw" gdi results with segment volumes approximated as elliptical prisms, or for manually supplied segment COMs. COM is calculated as a weighted mean of all segment COMs, with the segment mass as the weighting factor.
#'
#' @param x Either a data frame that is the output of gdi(..., return="all"), or a numeric vector of horizontal segment COM positions.
#' @param volumes An optional separate vector of volumes, required if x is not a data.frame containing volumes.
#' @param align alignment of the silhouette, if "h" (default) the silhouette is assumed to be horizontally aligned, if any other value (e.g. "v") then the silhouette is assumed to be vertically aligned.
#' @param subtract An optional separate vector of volumes, with length equal to the length or nrow() of x, to be subtracted from the volumes for the COM calculation.
#' @param densities An optional vector of segment densities, with length equal to the length or nrow() of x, to be multiplied with the volumes for the COM calculation. If both subtract and densities are supplied, the density is applied only to the "residual" volume that is left after subtraction.
#' @param scale Optional scale value (number of pixels to chosen unit of measurement)
#' @return An object of class numeric() containing the x coordinate of the center of mass of the shape, in pixels (or chosen units, if manually calculated)
#' @export hCOM
#' @importFrom stats weighted.mean
#' @examples
#' fdir <- system.file(package="gdi")
#' measuresil(file.path(fdir,"exdata","lat.png"), return="all")->lat_
#' measuresil(file.path(fdir,"exdata","dors.png"), return="all")->dors_
#' gdi(lat_, dors_, return="all")->gdiresults
#' hCOM(gdiresults)



hCOM<-function(x, volumes=NULL, align="h", subtract=NULL, densities=NULL, scale=1){
x_center<-x
if(is.data.frame(x)){
volumes<-x$V#look for collumn "V" containing segment volumes
x_center<-x$x_center#c(1:nrow(x))-0.5#save horizontal COM positions based on segment numbers
}
masses<-volumes
if(!is.null(subtract)){masses<-volumes-subtract}#subtract airspace volumes
if(!is.null(densities)){masses<-masses*densities}#and/or multiply by segment densities given as vector

if(align!="h"){
print("converted from top-down measurement")
    xmax<-attributes(x)$x_resolution
    return((xmax-weighted.mean(x_center, w=masses, na.rm=TRUE))/scale)
    #}
}else{
return(weighted.mean(x_center,w=masses, na.rm=TRUE)/scale)
}
}



##Function vCOM()
#'Finds the vertical (y axis, i.e. the axis parallel to the cross-section diameter) position of the center of mass (COM) of the volume. Experimental; only valid for "raw" gdi results with segment volumes approximated as elliptical prisms, or for manually supplied segment COMs. COM is calculated as a weighted mean of all segment COMs, with the segment mass as the weighting factor. Estimates have lower accuracy compared to hCOM, because cross-sectional geometry and variation in density throughout the cross-section is not taken into account.
#'
#' @param y A data.frame that is the output of gdi(..., return="all"), or a numeric vector containing vertical COM positions for segments
#' @param volumes An optional separate vector or data.frame (output of gdi(...,return="all") or vector of volumes.
#' @param subtract An optional separate vector of volumes, with length equal to the length or nrow() of x, to be subtracted from the volumes for the COM calculation.
#' @param densities An optional vector of segment densities, with length equal to the length or nrow() of x, to be multiplied with the volumes for the COM calculation. If both subtract and densities are supplied, the density is applied only to the "residual" volume that is left after subtraction.
#' @param scale Optional scale value (number of pixels to chosen unit of measurement)
#' @param from_top Whether the output coordinate should be measured from the top of the image (standard for image processing software), if TRUE, or from the bottom (standard for plotting in R (if FALSE, default). If TRUE, an attribute to y, containing the vertical dimension relative to which the measurement should be taken is required.
#' @return An object of class numeric() containing the y coordinate of the center of mass of the shape, in pixels (or chosen units, if manually calculated)
#' @export vCOM
#' @importFrom stats weighted.mean
#' @examples
#' fdir <- system.file(package="gdi")
#' measuresil(file.path(fdir,"exdata","lat.png"), return="all")->lat_
#' measuresil(file.path(fdir,"exdata","dors.png"), return="all")->dors_
#' gdi(lat_, dors_, return="all")->gdiresults
#' vCOM(gdiresults)


vCOM<-function(y,volumes=NULL,subtract=NULL, densities=NULL, scale=1, from_top=FALSE){
y_center<-y
    if(is.data.frame(y)){
    volumes<-y$V#look for collumn "center" containing vertical segment COM positions
    y_center<-y$y_center#look for collumn "center" containing vertical segment COM positions
    }

if(!is.null(volumes) & is.data.frame(volumes)){
volumes<-volumes$V#look for collumn "V" containing segment volumes
}


masses<-volumes
if(!is.null(subtract)){masses<-volumes-subtract}#subtract airspace volumes
if(!is.null(densities)){masses<-volumes*densities}#and/or multiply by segment densities given as vector

if(from_top==TRUE){
 ymax<-attributes(y)$y_resolution
 return((ymax-weighted.mean(y_center, w=masses, na.rm=TRUE))/scale)
    #}
}else{
return(weighted.mean(y_center, w=masses, na.rm=TRUE)/scale)
}
}


##Function expandConvexHull()
#' Estimate soft tissue expansion factors following Macaulay et al. 2023
#'
#' @param volume Numeric value containing volume of the convex hull (in liters). Defaults to NULL, in which case an isometric expansion factor is always returned. Setting for this parameter is ignored if method=="isometry"
#' @param method Method to use for deriving soft tissue expansion factor. Possible optios are "isometry" (default), "OLS" or "PGLS".
#' @param taxon Taxon whose specific regressions or expansion factors to use, can be "Croc_Lizard" for the expansions based on Crocodilians and Squamates, "Bird" for expansions based on birds of "All" for the pooled (allometric) or averaged (isometric) expansions.
#' @param segment Body segment for which expansion factor should be calculated. Possible options are: "Head", "Neck", "Trunk", "Tail", "Humerus", "Forearm", "Hand", "Thigh", "Shank", "MT" or "Pes"
#' @return A single numeric, containing the calculated expansion factor
#' @export expandConvexHull
#' @details This function returns a soft tissue expansion factor, calculated as the quotient of the skin volume and minimum convex hull volume of the selected body segment. This is based on expansion factors and regressions for birds and non-avian reptiles from Macaulay et al. 2023 (<10.1038/s41467-023-37317-y>). OLS and PGLS estimate the body volume allometrically following supplements 3, 4 and 5 of Macaulay et al., while isometry returns an isometric expansion factor following supplement 6. The body segment "Tail" should not be used in conjunction with the taxon "Bird", due to the absence of tails in this group – if this occurs nonethemess, an expansion factor of 1 is returned.
#' @examples
#' expandConvexHull(50,"isometry",taxon="Croc_Lizard",segment="Trunk")

expandConvexHull<-function(volume=NULL,method="isometry",taxon="Croc_Lizard",segment="Trunk"){

if(method%in%c("isometry","iso","Isometry","Iso") | is.null(volume)){
##isometric expansions
iso<-data.frame(Croc_Lizard=c(1.266, 8.646, 1.219, 3.369, 2.852, 2.866, 3.564, 5.102, 2.655, 3.494, 2.553),Bird=c(1.008, 3.825, 1.436, 1, 1.97, 1.736, 1.303, 4.538, 1.729, 0.792, 1.716))
iso$all<-c(1.137,6.235,1.328,3.369,2.411,2.301,2.434,4.820,2.192,2.143,2.135)
rownames(iso)<-c('Head', 'Neck', 'Trunk', 'Tail', 'Humerus', 'Forearm', 'Hand', 'Thigh', 'Shank', 'MT', 'Pes')

F_expansion<-iso[segment,taxon]
if(!(method%in%c("isometry","iso"))) warning("Selected method was NOT isometry, but no volume was provided; returning isometric expansion factor.")
}

if(method=="OLS"){
##allometric expansions OLS
if(taxon%in%c("Croc_Lizard","croc_lizard","Reptile","reptile")) ols_expansions<-data.frame(intercept=c(0.002, 0.73, 0.2, 0.387, 0.888, 0.677, 1.127, 0.951, 0.788, 0.908, 0.745),slope=c(0.978, 0.961, 1.039, 0.984, 1.09, 1.045, 1.106, 1.051, 1.075, 1.074, 1.067))
if(taxon%in%c("bird","Bird","Aves","avian","Avian")) ols_expansions<-data.frame(intercept=c(-0.085, 0.008, 0.213,0, 0.001, 0.017, 0.124, 0.487, 0.322, -0.402, 0.21),slope=c(0.982, 0.892, 1.018,1, 0.95, 0.963, 1.012, 0.975, 1.021, 0.946, 1.002))
if(taxon%in%c("all","All")) ols_expansions<-data.frame(intercept=c(-0.085, 0.008, 0.213, 0.001, 0.017, 0.124, 0.487, 0.322, -0.402, 0.21),slope=c(0.982, 0.892, 1.018,1, 0.95, 0.963, 1.012, 0.975, 1.021, 0.946, 1.002))
rownames(ols_expansions)<-c('Head', 'Neck', 'Trunk', 'Tail', 'Humerus', 'Forearm', 'Hand', 'Thigh', 'Shank', 'MT', 'Pes')

expansion_function<-function(x){exp(ols_expansions[segment,"slope"]*log(volume)+ols_expansions[segment,"intercept"])}

F_expansion<-expansion_function(volume)/volume
}

if(method=="PGLS"){
##allometric expansions PGLS
if(taxon%in%c("Croc_Lizard","croc_lizard","Reptile","reptile")) pgls_expansions<-data.frame(intercept=c(0.32, 0.711, 0.3, 0.27, 0.673, 0.795, 0.935, 1.013, 0.934, 0.941, 1.083),slope=c(1.052, 0.956, 1.076, 0.959, 1.047, 1.065, 1.072, 1.063, 1.107, 1.08, 1.134))
if(taxon%in%c("bird","Bird","Aves","avian","Avian")) pgls_expansions<-data.frame(intercept=c(-0.084, -0.048, 0.178,0, 0.03, 0.08, 0.128, 0.524, 0.33, -0.399, 0.228),slope=c(0.981, 0.879, 1.004,1, 0.951, 0.97, 1.01, 0.981, 1.023, 0.947, 1.003))
if(taxon%in%c("all","All")) pgls_expansions<-data.frame(intercept=c(-0.054, -0.063, 0.173,0, 0.006, 0.072, 0.115, 0.507, 0.324, -0.4, 0.219),slope=c(0.981, 0.876, 1.005,1, 0.95, 0.967, 1.011, 0.978, 1.021, 0.946, 1.009))

rownames(pgls_expansions)<-c('Head', 'Neck', 'Trunk', 'Tail', 'Humerus', 'Forearm', 'Hand', 'Thigh', 'Shank', 'MT', 'Pes')

expansion_function<-function(x){exp(pgls_expansions[segment,"slope"]*log(volume)+pgls_expansions[segment,"intercept"])}

F_expansion<-expansion_function(volume)/volume
}

return(F_expansion)

#Source: SuppData 4, 5 and 6 of Macaulay, S., Hoehfurtner, T., Cross, S.R.R., Marek, R.D., Hutchinson, J.R., Schachner, E.R., Maher, A.E., and Bates, K.T. 2023. Decoupling body shape and mass distribution in birds and their dinosaurian ancestors. Nature Communications: 14:1575. https://doi.org/10.1038/s41467-023-37317-y.
}
##


##Function plot_sil()
#' Plots a silhouette read by measuresil()
#'
#' @param sil A data frame that is the output of measuresil(..., return="all"), containing the center and the diameter of the silhouette at each value for x.
#' @param flip Logical indicating whether to flip axes (needed if measuresil() was performed using align="v", defaults to FALSE.
#' @param add Logical indicating whether to add silhoutte to an existing plot (defaults to FALSE)
#' @param expansion Numeric containing expansion factor to multiply with body diameters. Defaults to 1 (i.e. no expansion). XXX
#' @param xoffset Optional value by which to shift the silhouette on the x axis
#' @param yoffset Optional value by which to shift the silhouette on the y axis
#' @param alpha Opacity value for fill of polygon (defaults to 1)
#' @param col Fill color of polygon (defaults to "grey")
#' @param border Border color of polygon (defaults to "darkgrey")
#' @param scale Scale to use for plotting (given in pixels/unit). Defaults to 1.
#' @param xlab X axis label to use for plotting (if add=FALSE)
#' @param ylab Y axis label to use for plotting (if add=FALSE)
#' @param ... Other parameters to pass on to plot() or lines()
#' @return A plotted silhouette
#' @export plot_sil
#' @importFrom graphics lines
#' @importFrom graphics axis
#' @importFrom graphics polygon
#' @importFrom graphics par
#' @importFrom graphics plot.default
#' @importFrom grDevices col2rgb
#' @importFrom grDevices rgb
#' @examples
#' fdir <- system.file(package="gdi")
#' measuresil(file.path(fdir,"exdata","lat.png"), return="all")->lat_
#' plot_sil(lat_)

plot_sil<-function(sil, flip=FALSE, add=FALSE, expansion=1, xoffset=0, yoffset=0,alpha=1, col="grey", border="darkgrey", scale=1, xlab="", ylab="",...){

sil$diameter*expansion->sil$diameter#adjust diameters by expansion factor

if(flip==FALSE){
x<-c(1:nrow(sil),rev(1:nrow(sil)))+xoffset
y<-c(sil$center+sil$diameter/2, rev(sil$center-sil$diameter/2))+yoffset
}else{
y<--c(1:nrow(sil),rev(1:nrow(sil)))+nrow(sil)+yoffset
x<-c(sil$center+sil$diameter/2, rev(sil$center-sil$diameter/2))+xoffset
}

which(!is.na(y) & !is.na(x))->keep
x<-x[keep]
y<-y[keep]

x<-x/scale
y<-y/scale

aa<-function(col, alpha=0.5){
  apply(sapply(col, col2rgb)/255,2,function(x){rgb(x[1], x[2], x[3], alpha=alpha)})}
  
  if(alpha<1){col<-aa(col,alpha)}

if(add==TRUE){
polygon(y~x, border=border, col=col,...)

}else{

plot(y~x, type="n",axes=FALSE, xlab=xlab, ylab=ylab,...)
polygon(y~x, col=col, border=border)
axis(1)
}

}




##Function csI()
#' Calculates the second moment of area (=area moment of inertia, Ix and Iy) and polar moment of inertia (Iz or J) for a cross-section given as an image.
#'
#' @param image_file Image to be read. Images can be jpeg or png files, or a previously read image saved as an object in R.
#' @param threshold Reference value for color criterium after which pixels that are part of the silhouette should be differentiated from the background.
#' @param channel color channel to which to apply the threshold criterium. Default is 4 (alpha channel of rgba image). Channel setting needs to be adjusted depending on the color mode of the image used (e.g. there are two channels to choose from in a greyscale image, and 3 in an rgb image).
#' @param method Method for determining which pixels to count. Default "greater" counts pixels with value greater than threshold (e.g. higher opacity, in the case of an alpha channel). "less" counts pixels with a value less than the threshold. "not" counts all pixels not precisely matching threshold. Any other character string results in only pixels exactly matching the value given as threshold being counted.
#' @param scale Optional scale of the image (number of pixels per linear unit).
#' @param return What to return, defaults to returning both x and y second moments of area and polar moment of inertia for the entire shape (if return=="total"), otherwise returns raw data matrix for all pixels.
#' @return A numeric vector containing Ix, Iy and Iz for the shape (default), or a matrix containing area moments and coordinates for each pixel in the image, as well as area moments converted relative to the common centroid of the shape using parallel axis theorem.
#' @export csI
#' @examples
#' fdir <- system.file(package="gdi")
#' csI(file.path(fdir,"exdata","cross_section.png"))

csI<-function(image_file, threshold=0.5, channel=4, method="greater", scale=1, return="total"){
#load and save image data to variable named img
if(!is.character(image_file)){
img<-image_file
}else if(grepl(".jpg",image_file)==TRUE | grepl(".jpeg",image_file)==TRUE | grepl(".JPG",image_file)==TRUE){
img <- jpeg::readJPEG(image_file)
#read image if it is jpg
}else if(grepl(".png",image_file)==TRUE | grepl(".PNG",image_file)==TRUE){
img <- png::readPNG(image_file)}#read image if it is png


#get dimensions
nrows <- dim(img)[1]
ncols <- dim(img)[2]

if(!is.na(dim(img)[3])#if image has more than one channel, i.e. is an array, reduce the image to the channel to be analyzed
){
img<-img[,,channel]#img is now a matrix containing only the channel to be analyzed
}

width<-1/scale
height<-1/scale

#second moment of area for a square pixel
Ix<-(width*height^3)/12#x direction
Iy<-(width^3*height)/12#y direction

moments<-matrix(nrow=length(img), ncol=6)
colnames(moments)<-c("Ix","Iy","xpos", "ypos","Ix_","Iy_")

##loop through each pixel of image and save area moments and coordinates for each pixel belonging to the shape
i<-0#set index to start at 0
    for (x in 1:ncols){#loop through collumns
        for(y in 1:nrows){#loop through rows of collumn x
i<-i+1#increment the index
#select pixels belonging to shape depending on the method and threshold settings
if(method=="greater" & img[y,x]>threshold){
moments[i,1]<-Ix
moments[i,2]<-Iy
moments[i,3]<-x
moments[i,4]<-y
}else if(method=="less" & img[y,x]<threshold){
moments[i,1]<-Ix
moments[i,2]<-Iy
moments[i,3]<-x
moments[i,4]<-y
}else if(method=="not" & img[y,x]!=threshold){
moments[i,1]<-Ix
moments[i,2]<-Iy
moments[i,3]<-x
moments[i,4]<-y
}else if(img[y,x]==threshold){
moments[i,1]<-Ix
moments[i,2]<-Iy
moments[i,3]<-x
moments[i,4]<-y
}

        }
    }#end of loop

xcentroid<-mean(moments[,3], na.rm=T)#calculate centroid position for entire shape
ycentroid<-mean(moments[,4], na.rm=T)

#then use parallel axis theorem to convert pixel moments relative to overall centroid
moments[,5]<-Ix+width*height*((moments[,4]-ycentroid)/scale)^2
moments[,6]<-Iy+width*height*((moments[,3]-xcentroid)/scale)^2

#sum up all individual pixel moments for shape:
I_x_total<-sum(moments[,5],na.rm=T)
I_y_total<-sum(moments[,6],na.rm=T)

#calculate polar moment
I_z_total<-sum(I_x_total, I_y_total)

##return results
if(return=="total"){
total<-c(I_x=I_x_total, I_y=I_y_total, I_z=I_z_total)

attr(total, "centroid")<-c(x=xcentroid, y=ycentroid)
bbox<-rbind(c(min(moments[,3], na.rm=T),max(moments[,3], na.rm=T)),c(min(moments[,4], na.rm=T),max(moments[,4], na.rm=T)))
rownames(bbox)<-c("x","y")
colnames(bbox)<-c("min","max")
attr(total, "bounding_box")<-bbox


return(total)


}else(return(moments))

}




##Function rotI
#' Calculates the rotational inertia of a body. Only works with simple circular/elliptical and rectangular cross-sections, thus pixel-precise estimates are recommended.
#'
#' @param x Either a data frame that is structured like output of gdi(..., return="all") containing masses and diameters for pixel-wide segments, or a numeric vector of segment COM positions.
#' @param y An optional vector of vertical (dorsoventral) segment COM positions.
#' @param dors_diam An optional vector of transverse diameters of the silhouette, required if not contained in x.
#' @param lat_diam An optional vector of vertical diameters of the silhouette, required if not contained in x. Needed if inertia for "roll" or "pitch" should be calculated.
#' @param axis_coord An optional coordinate of the axis of rotation (in original units, i.e. pixels), defaults to the center of mass of the entire volume if not set.
#' @param axis Axis of rotation, defaults to "yaw" (i.e. rotation around vertical axis), can also be "pitch" (rotation around transverse axis) or "roll" (rotation around horizontal axis). For yaw rotation, the body is assumed to be bilaterally symmetrical, whereas for pitch rotation, dorsoventral variation in COM of segments is taken into account.
#' @param volumes An optional separate vector of volumes, required if x is not a data.frame containing them.
#' @param corr An optional correction factor for the cross-sectional shape, given as the ratio between the characteristic mass moment of inertia of a plane with the given shape (e.g. determined by cscorr()) and an elliptical plane with the same diameters and assigned mass. Allows the calculation of moments of inertia for bodies with arbitrary cross-sectional shapes.
#' @param densities An optional vector of segment densities, with length equal to the length or nrow() of x, to be multiplied with the volumes to calculate masses used in the inertial calculation.
#' @param scale Scale value, i.e. number of pixels representing 1 m
#' @return A numeric vector containing: The total mass (on the basis of gdi volumes and optional densities), the rotational inertia of the shape using a point mass approximation of each segment (yaw/pitch rotation only), rotational inertia using a cylindrical approximation for each segment, rotational inertia using a cuboidal approximation (note that this only changes the mass distribution, while segment masses are still assumed to correspond to gdi results multiplied by optional densities), and rotational inertia using a corrected cylindrical approximation based on value for corr.
#' @export rotI
#' @examples
#' fdir <- system.file(package="gdi")
#' measuresil(file.path(fdir,"exdata","lat.png"), return="all")->lat_
#' measuresil(file.path(fdir,"exdata","dors.png"), return="all")->dors_
#' gdi(lat_, dors_, return="all")->gdiresults
#' rotI(gdiresults)


rotI<-function(x,y=NULL, dors_diam=NULL, lat_diam=NULL, axis_coord=NULL, axis="yaw", volumes=NULL, densities=1, corr=1, scale=1){
x_center<-x

if(is.data.frame(x)){
volumes<-x$V#look for collumn "V" containing segment volumes

x_center<-x$x_center#c(1:nrow(x))-0.5#save horizontal COM positions based on segment numbers
if(is.numeric(x$y_center)){
y_center<-x$y_center#save vertical COM positions (perpendicular to cross-section) based on vertical segment COMS
}

if(is.null(dors_diam)){dors_diam<-x$zdiam_raw}#get dorsal and lateral diameters from data.frame, if not given elsewhere
if(is.null(lat_diam)){lat_diam<-x$ydiam_raw}

segL<-x$slice_length*(x$zdiam_raw/x$zdiam_scaled)
}else{
segL<-numeric(length(x_center))
for(i in 1:length(x_center)-1){
segL[i]<-x_center[i+1]-x_center[i]
}
segL[length(x_center)]<-segL[length(x_center)-1]
}

if(is.numeric(y)){
y_center<-y}
if(!exists("y_center")){
y_center<-rep(max(lat_diam, na.rm=T)/2, length(x_center))#define vertical COM positions as constant if y values provided

}

masses<-volumes*densities#multiply by segment densities

if(is.null(axis_coord)){
if(axis=="yaw" | axis== "y"){axis_coord<-gdi::hCOM(x=x_center, volumes=masses)#for yaw axis
}else if(axis=="pitch" | axis== "x" | axis=="roll" | axis=="x"){
axis_coord<-c(gdi::hCOM(x=x_center, volumes=masses), gdi::vCOM(y=y_center, volumes=masses))
}#for roll or pitch axis
}

#set scale
x_center/scale->x_center
if(is.numeric(y_center)){
y_center/scale->y_center}
axis_coord/scale->axis_coord
dors_diam/scale->dors_diam
lat_diam/scale->lat_diam
segL/scale->segL


    ##for yaw rotation
if(axis=="yaw" | axis=="y"){

##approximation assuming sections are point masses
    sum((x_center-axis_coord)^2*masses)->point_I

##approximation using cylindrical sections
segment_I<-1/12*masses*segL^2+1/4*masses*(dors_diam/2)^2 # formula for solid cylinder
segment_I_corr<-segment_I+(x_center-axis_coord)^2*masses#parallel axis theorem to convert for rotation around chosen axis or COM

    sum(segment_I_corr, na.rm=T)->exact_I_cylinder#sum up results

##approximation using cuboidal sections
segment_I2<-1/12*masses*(segL^2+dors_diam^2) #equivalent formula from https://en.wikiversity.org/wiki/PlanetPhysics/Rotational_Inertia_of_a_Solid_Cylinder
#segment_I2<-1/12*masses*segL^2+masses*(dors_diam/2)^2#cuboid formula
segment_I2_corr<-segment_I2+(x_center-axis_coord)^2*masses#use parallel axis theorem to convert around arbitrary axis or COM of entire shape

    sum(segment_I2_corr, na.rm=T)->exact_I_cuboid#sum up results

##approximation using cylindrical sections * correction factor
segment_I_<-segment_I*corr
segment_I_corr_<-segment_I_+(x_center-axis_coord)^2*masses#use parallel axis theorem to convert around arbitrary axis or COM of entire shape
    sum(segment_I_corr_, na.rm=T)->exact_I_cylinder_corrected#sum up results

    c(total_mass=sum(masses),I_point_masses=point_I, I_elliptical_cs=exact_I_cylinder,I_rectangular_cs=exact_I_cuboid, I_corrected_cs=exact_I_cylinder_corrected)->results_yaw

    
    
}else if(axis=="pitch" | axis=="z"){#for pitch rotation

##approximation assuming sections are point masses
    sum(((x_center-axis_coord[1])^2+(y_center-axis_coord[2])^2)*masses, na.rm=T)->point_I

##approximation using cylindrical sections
segment_I<-1/12*masses*segL^2+1/4*masses*(lat_diam/2)^2 # formula for solid cylinder
segment_I_corr<-segment_I+((x_center-axis_coord[1])^2+(y_center-axis_coord[2])^2)*masses#use parallel axis theorem to convert around arbitrary axis or COM of entire shape
    sum(segment_I_corr, na.rm=T)->exact_I_cylinder#sum up results

##approximation using cuboidal sections
segment_I2<-1/12*masses*(segL^2+lat_diam^2) #equivalent formula from https://en.wikiversity.org/wiki/PlanetPhysics/Rotational_Inertia_of_a_Solid_Cylinder
segment_I2_corr<-segment_I2+((x_center-axis_coord[1])^2+(y_center-axis_coord[2])^2)*masses#use parallel axis theorem to convert around arbitrary axis or COM of entire shape
    sum(segment_I2_corr, na.rm=T)->exact_I_cuboid#sum up results

##approximation using cylindrical sections * correction factor
segment_I_<-segment_I*corr
segment_I_corr_<-segment_I+((x_center-axis_coord[1])^2+(y_center-axis_coord[2])^2)*masses#use parallel axis theorem to convert around arbitrary axis or COM of entire shape
    sum(segment_I_corr_, na.rm=T)->exact_I_cylinder_corrected#sum up results

    c(total_mass=sum(masses),I_point_masses=point_I, I_elliptical_cs=exact_I_cylinder,I_rectangular_cs=exact_I_cuboid, I_corrected_cs=exact_I_cylinder_corrected)->results_pitch
}else if(axis=="roll" | axis=="x"){
##approximation using elliptical sections
segment_I<-1/4*masses*((lat_diam/2)^2+(dors_diam/2)^2)
segment_I_<-segment_I*corr

##approximation using cuboidal sections
segment_I2<-1/12*masses*(lat_diam^2+dors_diam^2)


if(is.numeric(y_center) & length(y_center)==length(x_center)){
segment_I_corr<-segment_I+(y_center-axis_coord[2])^2*masses#use parallel axis theorem to convert around arbitrary axis or COM of entire shape
segment_I2_corr<-segment_I2+(y_center-axis_coord[2])^2*masses

segment_I_corr_<-segment_I_+(y_center-axis_coord[2])^2*masses



}else{segment_I_corr<-segment_I
segment_I_corr_<-segment_I_
segment_I2_corr_<-segment_I2}

    sum(segment_I_corr, na.rm=T)->exact_I_ellipsoid#sum up results
    sum(segment_I_corr_, na.rm=T)->exact_I_ellipsoid_corrected#sum up results
    sum(segment_I2_corr, na.rm=T)->exact_I_cuboid#sum up results
    c(total_mass=sum(masses),I_elliptical_cs=exact_I_ellipsoid,I_rectangular_cs=exact_I_cuboid,I_corrected_cs=exact_I_ellipsoid_corrected)->results_roll

}

##return results
if(axis=="yaw" | axis=="y"){
return(results_yaw)}else if(axis=="pitch" | axis=="z"){
return(results_pitch)}else if(axis=="roll" | axis=="x"){
return(results_roll)
}


}##





##function transfer_ratio()
#' Transfer a vector of aspect ratios onto another body profile
#' @param lat diameters in lateral view, supply to estimate transverse diameters from dorsoventral ones
#' @param dors diameters in dorsal (or ventral) view, supply to estimate dorsoventral diameters from transverse ones
#' @param indices optional vector of indices to apply to subset lat or dors
#' @param lat0 (dorsoventral) diameters in lateral view for template
#' @param dors0 (transverse) diameters in dorsal/ventral view for template (vector is subsetted to same length as lat0 if lengths don’t match)
#' @param indices0 optional vector of indices to subset lat0 and dors0
#' @param smooth width of smoothing interval (as fraction of full silhouette length)
#' @param return what to return, either "diameters" to return diameters, or "ratios", to return ratios
#' @param ... additional arguments to pass on to paleoDiv::rmeana() for computation of running mean of diameter ratios
#' @return Either a numeric vector of interpolated transverse or dorsoventral diameters or width/depth or depth/width ratios depending on the chosen settings, or a data.frame of same format as dors0 or lat0 (with the "center" column copied)
#' @export transfer_ratio
#' @examples
#' fdir <- system.file(package="gdi")
#' lat <- measuresil(file.path(fdir,"exdata","lat.png"), return="all")
#' dors <- measuresil(file.path(fdir,"exdata","dors.png"), return="all")
#' gdi(lat,dors,scale=100) #real volume
#' dors_est <- transfer_ratio(lat=lat,lat0=lat,dors0=dors,smooth=0.005)
#' gdi(lat,dors_est,scale=100) #volume with dorsal view interpolated
#' lat_est <- transfer_ratio(dors=dors,lat0=lat,dors0=dors,smooth=0.005)
#' gdi(lat_est,dors,scale=100) #volume with lateral view interpolated
#' plot_sil(lat,asp=1) #visualize results
#' plot_sil(dors,add=TRUE)
#' plot_sil(lat_est,add=TRUE,col="blue",alpha=0.3)
#' plot_sil(dors_est,add=TRUE,col="blue",alpha=0.3)


transfer_ratio<-function(lat=NULL,dors=NULL,indices=NULL,lat0,dors0,indices0=NULL,smooth=0.005,return="diameters",...){

#define advanced running mean function
rmeana<-function(x0, y0, x1 = NULL, plusminus = 0.005,...){
    if (is.null(x1)) {x_ <- x0}else{x_ <- x1}
    y_ <- rep(NA, length(x_))
    
    for (i in 1:length(x_)) {
        indices <- which(abs(as.numeric(x0 - x_[i])) <= plusminus)
        if (length(indices) == 0) indices <- which.min(abs(x0 - x_[i]))
		y_[i] <- mean(y0[indices], na.rm = TRUE)
    }
    
    return(y_)} #end simplified rmeana function definition, see also package "paleoDiv"
    
#prepare template: equalize lengths
if(is.data.frame(lat0)){
if("center"%in%colnames(lat0)) lat0$center->lcenter
lat0$diameter->lat0
}
if(is.data.frame(dors0)){
if("center"%in%colnames(dors0)) dors0$center->dcenter
dors0$diameter->dors0
}
dors0[1:length(lat0)]->dors0

#subset template
if(is.null(indices0)) c(1:length(lat0))->indices0
lat0[indices0]->lat0
dors0[indices0]->dors0
if(exists("dcenter") && is.numeric(dcenter)) dcenter[indices0]->dcenter
if(exists("lcenter") && is.numeric(lcenter)) lcenter[indices0]->lcenter

#compute ratios
dors0/lat0->ratios0
c(1:length(lat0))->x0
x0_rel<-x0/length(lat0)
#

##do interpolation
if(!is.null(lat)){ ##lateral given, dorsal unknown
df<-FALSE
if(is.data.frame(lat)){
df<-TRUE
lat$diameter->lat
}
if(is.null(indices)) indices<-c(1:length(lat))
lat[indices]->lat

#
c(1:length(lat))->x
x_rel<-x/length(lat)

interpolated_ratios<-rmeana(x0=c(0,x0_rel), y0=c(ratios0[1],ratios0), x1=x_rel, plusminus=smooth/2,...)


if(any(is.na(interpolated_ratios))){ #remove nas resulting from interpolation, replace with closest non-na value
which(is.na(interpolated_ratios))->narat
which(!is.na(interpolated_ratios))->nnarat

for(i in narat){which(abs(nnarat-i)==min(abs(nnarat-i)))->j
interpolated_ratios[i]<-mean(interpolated_ratios[nnarat[j]],na.rm=TRUE)
}

}

if(any(is.infinite(interpolated_ratios))){ #remove infs resulting from interpolation, replace with closest non-inf value
which(is.infinite(interpolated_ratios))->narat
which(!is.infinite(interpolated_ratios))->nnarat

for(i in narat){which(abs(nnarat-i)==min(abs(nnarat-i)))->j
interpolated_ratios[i]<-mean(interpolated_ratios[nnarat[j]],na.rm=TRUE)
}

}


interpolated_ratios*lat->interpolated_diameters
if(df && exists("dcenter")){
mean(dcenter,na.rm=TRUE)->center
rep(center,length(interpolated_diameters))->center
}


}else if(!is.null(dors)){ ##dorsal given, lateral unknown
df<-FALSE
if(is.data.frame(dors)){
df<-TRUE
dors$diameter->dors
}
if(is.null(indices)) indices<-c(1:length(dors))
dors[indices]->dors

#
c(1:length(dors))->x
x_rel<-x/length(dors)

interpolated_ratios<-rmeana(x0=c(0,x0_rel), y0=c(ratios0[1],ratios0), x1=x_rel, plusminus=smooth/2,...)

if(any(is.na(interpolated_ratios))){ #remove nas resulting from interpolation, replace with closest non-na value
which(is.na(interpolated_ratios))->narat
which(!is.na(interpolated_ratios))->nnarat

for(i in narat){which(abs(nnarat-i)==min(abs(nnarat-i)))->j
interpolated_ratios[i]<-mean(interpolated_ratios[nnarat[j]],na.rm=TRUE)
}

}

if(any(is.infinite(interpolated_ratios))){ #remove infs resulting from interpolation, replace with closest non-inf value
which(is.infinite(interpolated_ratios))->narat
which(!is.infinite(interpolated_ratios))->nnarat

for(i in narat){which(abs(nnarat-i)==min(abs(nnarat-i)))->j
interpolated_ratios[i]<-mean(interpolated_ratios[nnarat[j]],na.rm=TRUE)
}

}

if( df && exists("lcenter") && is.numeric(lcenter) ){
rmeana(x0=x0_rel, y0=lcenter, x1=x_rel, plusminus=smooth/2,...)->center #interpolate lateral view center
}

dors/interpolated_ratios->interpolated_diameters
interpolated_ratios<-1/interpolated_ratios
}else{stop("lat or dors need to be supplied")}

interpolated_diameters[is.na(interpolated_diameters)]<-0
if(any(is.infinite(interpolated_diameters))) interpolated_diameters[is.infinite(interpolated_diameters)]<-0

if(any(is.infinite(interpolated_ratios))) interpolated_ratios[is.infinite(interpolated_ratios)]<-NA

#interpolated_ratios[is.na(interpolated_ratios)]<-0

if(exists("center") && is.numeric(center)){
interpolated_diameters<-data.frame(diameter=interpolated_diameters,center=center)
interpolated_diameters$center[which(interpolated_diameters$diameter==0)]<-NA
}
if(return%in%c("Diameters","diameters","diam","D")) return(interpolated_diameters)
if(return%in%c("ratios","Ratios","rat","R")) return(interpolated_ratios)

}##




