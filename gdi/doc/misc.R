## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(gdi)

## ---- fig.height=5, fig.width=7-----------------------------------------------
fdir <- system.file(package="gdi")

imghist(file.path(fdir,"exdata","femur.png"), channel=1, unique=TRUE)

## -----------------------------------------------------------------------------
cscorr(file.path(fdir,"exdata","femur.png"), channel=1, method="less", threshold=0.2, return="area")->compacta #measure compact bone (black pixels)

cscorr(file.path(fdir,"exdata","femur.png"), channel=1, method="greater", threshold=0.8, return="area")->medulla #measure medulla (white pixels)

medulla/(compacta+medulla) #calculate the airspace proportion

## -----------------------------------------------------------------------------
matrix(rep(1,100),nrow=10,ncol=10)->m #generate sample matrix representing the pixel color values of an image
csI(m) #return the second moments of area (Ix and Iy) and the polar moment of inertia (Iz) of the cross-section

## -----------------------------------------------------------------------------
csI(file.path(fdir,"exdata","femur.png"), channel=1, method="less", threshold=0.2)

## -----------------------------------------------------------------------------
csI(file.path(fdir,"exdata","femur.png"), channel=1, method="less", threshold=0.2, scale=676/11.6)

## -----------------------------------------------------------------------------
fdir <- system.file(package="gdi")

measurements_lateral <- measuresil(file.path(fdir,"exdata","lat.png"), return="all") # using the "return" parameter, we can make measuresil output a data.frame containing the centers on the y axis
measurements_dorsal <- measuresil(file.path(fdir,"exdata","dors.png")) #because our organism is bilaterally symmetric, we are not interested in the location of the COM on the transverse axis, and can therefore disregard recording centers for the dorsal view

gdi(measurements_lateral, measurements_dorsal, scale=100, return="all")->gdiresults # the parameter setting for "return" makes gdi() return a data.frame with all individual segment volumes, rather than only a single volume for the entire shape

## -----------------------------------------------------------------------------
rotI(gdiresults, scale=1000, axis="yaw", densities=1)

## -----------------------------------------------------------------------------
cscorr(file.path(fdir,"exdata","cross_section.png"),return="rotI")->c
c

## -----------------------------------------------------------------------------
rotI(gdiresults, scale=1000, corr=c[1])

