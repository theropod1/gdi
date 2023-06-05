## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(gdi)

## -----------------------------------------------------------------------------

fdir <- system.file(package="gdi")

measurements_lateral <- measuresil(file.path(fdir,"exdata","lat.png"))
measurements_dorsal <- measuresil(file.path(fdir,"exdata","dors.png"))

## -----------------------------------------------------------------------------
length(measurements_lateral)
length(measurements_dorsal)

## -----------------------------------------------------------------------------
gdi(measurements_lateral, measurements_dorsal, scale=100, method="raw")
gdi(measurements_lateral, measurements_dorsal, scale=100, method="smooth")

## -----------------------------------------------------------------------------
sil <- c(0,1)
gdi(sil, sil, method="raw", scale=1)
gdi(sil, sil, method="smooth", scale=1)


## -----------------------------------------------------------------------------
gdi(measurements_lateral,measurements_dorsal,k=2.3, scale=100)

## -----------------------------------------------------------------------------
sellipse.coo(2.0)->ellipse
plot(ellipse$x,ellipse$y,col="grey", type="l")
polygon(ellipse$x,ellipse$y,col="grey", border="grey")

sellipse.coo(2.3)->se2.3
lines(se2.3$x,se2.3$y, col="blue") #plot a superellipse with exponent 2.3

sellipse.coo(3)->se3
lines(se3$x,se3$y, col="red") #plot a superellipse with exponent 3

## -----------------------------------------------------------------------------
fdir <- system.file(package="gdi")
correction_factor <- cscorr(file.path(fdir,"exdata","cross_section.png"))
print(correction_factor)

## -----------------------------------------------------------------------------
gdi(measurements_lateral, measurements_lateral*0.9, scale=100)

## -----------------------------------------------------------------------------
aspect_ratio <- cscorr(file.path(fdir,"exdata","cross_section.png"))
print(aspect_ratio)

