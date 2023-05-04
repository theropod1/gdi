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

##Function cscorr()
#'Measure and analyze cross-sectional geometry
#'
#' @param image_file Image to be read. Images can be jpeg or png files, or a previously read image saved as an object in R.
#' @param threshold Reference value for colour criterium after which pixels that are part of the cross-section are differentiated from the background.
#' @param channel Colour channel to which to apply the threshold criterium. Default is 4 (alpha channel of rgba image). Channel setting needs to be adjusted depending on the colour mode of the image used (e.g. there are two channels to choose from in a greyscale image, and 3 in an rgb image).
#' @param method Method for determining which pixels to count. Default "greater" counts pixels with value greater than threshold (e.g. higher opacity, in the case of an alpha channel). "less" counts pixels with a value less than the threshold. Any other character string results in only pixels exactly matching the value given as threshold being counted.
#' @param return What value to return. Possible values are "area_corr" (Default, returns ratio between measured area and area of ellipse with same horizontal and vertical diameters), "aspect_ratio" (returns aspect ratio), "diameters" (returns diameters) and "area" (returns area). Any other value for this parameter will prompt the function to return a vector containing all of these.
#' @param k optional superellipse exponent for the (super)ellipse to which the measurements should be compared (for the "area_corr" setting for the parameter return).
#' @param scale Optional scale of the image (for raw area measurements).
#' @return A single number, either the ratio between the measured cross-sectional area and that of an ellipse with the same vertical and horizontal diameters, or the measured area (if compare==FALSE).
#' @import jpeg
#' @import png
#' @export cscorr
#' @examples
#' fdir <- system.file(package="gdi")
#' correction_factor <- cscorr(file.path(fdir,"exdata","cross_section.png"))


cscorr <- function(image_file, threshold=0.5, channel=4, method="greater", return="area_corr", k=2.0, scale=1) {
#load and save image data to variable img
if(grepl(".jpg",image_file)==TRUE | grepl(".jpeg",image_file)==TRUE | grepl(".JPG",image_file)==TRUE){
img <- jpeg::readJPEG(image_file)}#read image if it is jpg

if(grepl(".png",image_file)==TRUE | grepl(".PNG",image_file)==TRUE){
img <- png::readPNG(image_file)}#read image if it is png

if(!is.character(image_file)){
img<-image_file
}

# measure cross-section depth
nrows <- dim(img)[1]
ncols <- dim(img)[2]

##find maximum depth of cross-section
d <- numeric(0)
#loop through vertical lines of pixels, find all rows in which there are pixels that are part of the cross-section
for(x in 1:ncols){
    if(method=="greater"){
    d <- union(d,which(img[,x,channel]>threshold))
    }else if(method=="less"){
    d <- union(d,which(img[,x,channel]<threshold))
    }else{
    d <- union(d,which(img[,x,channel]==threshold))
    }
}
vdiam <- length(d) #calculate maximum depth in pixels

##find maximum width of cross-section
d <- numeric(0)
#loop through horizontal lines of pixels, find all columns in which there are pixels that are part of the cross-section
for(y in 1:nrows){
    if(method=="greater"){
    d <- union(d,which(img[y,,channel]>threshold))
    }else if(method=="less"){
    d <- union(d,which(img[y,,channel]<threshold))
    }else{
    d <- union(d,which(img[y,,channel]==threshold))
    }
}
hdiam <- length(d) #calculate maximum width in pixels


# calculate (super)elliptical cross-sectional area:
ecomp <- sellipse(vdiam/2/scale, hdiam/2/scale, k)

# calculate actual cross-sectional area
    if(method=="greater"){
    area <- sum(img[,,channel]>threshold)/scale^2
    } else if(method=="less"){
    area <- sum(img[,,channel]<threshold)/scale^2
    }else{
    area <- sum(img[,,channel]==threshold)/scale^2
    }
  
  # Return the calculated area or ratio
  if(return=="area"){return(area)
  }else if(return=="area_corr"){
  return(area/ecomp)
  }else if(return=="aspect_ratio"){
  return(vdiam/hdiam)
  }else if(return=="diameters"){
  diam<-c(x=hdiam, y=vdiam)
  return(diam)
  }else{full<-c(x=hdiam, y=vdiam, asp=hdiam/vdiam, area=area, area_corr=area/ecomp)
  return(full)
    }
}


##Function measuresil()
#'Take pixel-by-pixel measurements of a silhouette in jpeg or png format for use with the gdi function.
#'
#' @param image_file Image to be read. Images can be jpeg or png files, or a previously read image saved as an object in R.
#' @param threshold Reference value for colour criterium after which pixels that are part of the silhouette are differentiated from the background.
#' @param channel Colour channel to which to apply the threshold criterium. Default is 4 (alpha channel of rgba image). Channel setting needs to be adjusted depending on the colour mode of the image used (e.g. there are two channels to choose from in a greyscale image, and 3 in an rgb image).
#' @param method Method for determining which pixels to count. Default "greater" counts pixels with value greater than threshold (e.g. higher opacity, in the case of an alpha channel). "less" counts pixels with a value less than the threshold. Any other character string results in only pixels exactly matching the value given as threshold being counted.
#' @import jpeg
#' @import png
#' @export measuresil
#' @return A numeric vector giving the measurements of the silhouette
#' @examples
#' fdir <- system.file(package="gdi")
#' lat <- measuresil(file.path(fdir,"exdata","lat.png"))

measuresil<-function(image_file, threshold=0.5, channel=4, method="greater"){
#load and save image data to variable named img
if(grepl(".jpg",image_file)==TRUE | grepl(".jpeg",image_file)==TRUE | grepl(".JPG",image_file)==TRUE){
img <- jpeg::readJPEG(image_file)}#read image if it is jpg

if(grepl(".png",image_file)==TRUE | grepl(".PNG",image_file)==TRUE){
img <- png::readPNG(image_file)}#read image if it is png

if(!is.character(image_file)){
img<-image_file
}

# Loop through each vertical line of pixels
nrows <- dim(img)[1]
ncols <- dim(img)[2]
depths <- rep(0, ncols)
for (x in 1:ncols) {
  depth <- 0
  for (y in 1:nrows) {
    # Get the color of the selected colour channel of the pixel
    color <- img[y, x, channel]
    
    #increment the depth if colour is greater than threshhold, e.g. transparency
    if(method=="greater"){
    if (color > threshold) {
      depth <- depth + 1
    }
    }else if(method=="less"){
    if (color < threshold) {
      depth <- depth + 1
    }
    }else{
    if (color == threshold) {
      depth <- depth + 1
    }}}
    
  
  depths[x] <- depth
}

return(depths)
}


##Function gdi()
#'Estimate volume using Graphic Double Integration.
#'
#' @param lat Measurements of diameter in lateral view/first of two orthogonal views to be used with the gdi. Can be either a numeric vector, or a text file to be scanned. Defaults to "lat.txt".
#' @param dors Measurements of diameter in dorsal view/second of two orthogonal views to be used with the gdi. Can be either a numeric vector, or a text file to be scanned. Must be the same length as lat. Defaults to "dors.txt".
#' @param scale Scale of the data in terms of how many units of the input data are in one side of the desired unit of output volume. Defaults to 10.
#' @param sliceL Length of individual segments to be used in the GDI. Defaults to 1/scale.
#' @param method Method to be used for the GDI. Default "raw" setting calculates each segment as an elliptical cylinder with volume = Area * SliceL. Any other string will result in volume being calculated as an elliptical frustum with base areas based on the measurements of segments i and i+1.
#' @param k Superellipse exponent to be used for the cross-sectional area. Defaults to 2.0 (normal ellipse).
#' @param corr Correction factor for area of cross-sections, calculated as the ratio between the actual cross-sectional area and that of a (super)ellipse (depending on the specified exponent k) with the same diameters. This setting enables the function to account for complex, non-elliptical cross-sections. Default value is 1, i.e. no correction. Can be either a single number, or a numeric vector of the same length as lat and dors (in the case of a changing cross-sectional geometry along the length of the body).
#' @param smooth.ends If method != "raw", specify whether first and last segments should be left raw, or taper to 0 (i.e. be approximated as cones). Only applies if there are no leading or following zeros in the measurement vectors.
#' @param return Determines whether to report the estimated total volume (if default/"total"), or a data.frame with segment radii, areas and volumes (if left empty of any other character string.
#' @return Either a single number representing the total volume estimated, or (if return!="total") a data.frame() containing columns with the radii in both dimensions, the estimated elliptical or superelliptical areas, and the segment volumes.
#' @export gdi
#' @examples
#' lateral <- rep(2,4) #generate example data
#' dorsal <- rep(2,4)
#' gdi(lat=lateral, dors=dorsal, scale=10, method="raw", k=2.0)
#' gdi(lat=lateral, dors=lateral/2, scale=10, method="smooth", k=2.3)


gdi<-function(lat, dors, scale=10, sliceL=1/scale, method="raw", k=2.0, corr=1, smooth.ends=FALSE, return="total"){

if(is.character(lat)){scan(lat)->lat}

if(is.character(dors)){scan(dors)->dors}

sil<-data.frame(lat/scale,dors/scale)#scale-adjust and save in dataframe

#estimate cross-sections as (super)ellipses
sil$A<-sellipse(sil$lat/2, sil$dors/2, k)

if(method=="raw"){
    sil$V<-sil$A*sliceL#raw formula, approximating each segment as a cylinder
}else{
    for(i in 1:(length(lat)-1)){
    sil$V[i]<-(sil$A[i]+sil$A[i+1]+sqrt(sil$A[i]*sil$A[i+1]))/3*ifelse(length(sliceL)>1,sliceL[1],sliceL)
    #approximating segments as frusta, with A[i] and A[i+1] as base areas
    }
    if(smooth.ends==TRUE){
    sil$V[1]<-ifelse(length(sliceL)>1,sliceL[1],sliceL)*sil$A[1]/3#conical first segment (tapers to 0)
    sil$V[length(lat)]<-ifelse(length(sliceL)>1,sliceL[length(lat)],sliceL)*sil$A[length(lat)]/3#conical last segment (tapers to 0)
    }else{
    sil$V[length(lat)]<-ifelse(length(sliceL)>1,sliceL[length(lat)],sliceL)*sil$A[length(lat)]
    }
}

sil$V <- sil$V*corr

sum(sil$V)->res#sum up segments

if(return=="total"){
return(res)
}else{
return(sil)}
}
