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
gdi(measurements_lateral, measurements_dorsal, scale=100, method="smooth")

## -----------------------------------------------------------------------------
sil <- c(0,1)
gdi(sil, sil, method="raw", scale=1)
gdi(sil, sil, method="smooth", scale=1)


## -----------------------------------------------------------------------------
gdi(measurements_lateral,measurements_dorsal,k=2.3, scale=100)

## -----------------------------------------------------------------------------
fdir <- system.file(package="gdi")
correction_factor <- cscorr(file.path(fdir,"exdata","cross_section.png"))
print(correction_factor)

