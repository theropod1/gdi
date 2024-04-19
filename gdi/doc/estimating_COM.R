## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(gdi)

## -----------------------------------------------------------------------------
fdir <- system.file(package="gdi")

measurements_lateral <- measuresil(file.path(fdir,"exdata","lat.png"), return="all") # using the "return" parameter, we can make measuresil output a data.frame containing the centers on the y axis
measurements_dorsal <- measuresil(file.path(fdir,"exdata","dors.png")) #because our organism is bilaterally symmetric, we are not interested in the location of the COM on the transverse axis, and can therefore disregard recording centers for the dorsal view

## -----------------------------------------------------------------------------
knitr::kable(head(measurements_lateral, 10))

## -----------------------------------------------------------------------------
gdi(measurements_lateral, measurements_dorsal, scale=100, return="all")->gdiresults # the parameter setting for "return" makes gdi() return a data.frame with all individual segment volumes, rather than only a single volume for the entire shape

## -----------------------------------------------------------------------------
knitr::kable(head(gdiresults))

## -----------------------------------------------------------------------------
hCOM(gdiresults)
vCOM(gdiresults)

## ---- fig.cap="Plotted silhouette, with horizontal and vertical COM position marked.", fig.height=5, fig.width=7----
plot_sil(measurements_lateral, asp=1)
points(hCOM(gdiresults),vCOM(gdiresults))
abline(v=hCOM(gdiresults), lty=2)
abline(h=vCOM(gdiresults), lty=2)

## -----------------------------------------------------------------------------
hindlimb_lateral <- measuresil(file.path(fdir,"exdata","hl.png"),align="v", return="all")
forelimb_lateral <- measuresil(file.path(fdir,"exdata","fl.png"), align="v", return="all")

hl<-gdi(hindlimb_lateral,0.7*hindlimb_lateral,scale=100, return="all")
fl<-gdi(forelimb_lateral, 0.7*forelimb_lateral, scale=100, return="all")

## ---- fig.cap="Plotted silhouette, including limbs and individual segment COM positions", results=FALSE, message=FALSE, fig.height=5, fig.width=7----
plot_sil(measurements_lateral, asp=1)
plot_sil(forelimb_lateral, flip=T, add=T)#for silhouettes that were aligned, but digitized with align="v", we need to specify flip=TRUE here
plot_sil(hindlimb_lateral, flip=T, add=T)

points(hCOM(gdiresults),vCOM(gdiresults))#plot COM of axial segment
points(vCOM(fl),hCOM(fl,align="v"))#plot COM of forelimb segment
points(vCOM(hl),hCOM(hl,align="v"))#plot COM of forelimb segment

## -----------------------------------------------------------------------------
gdi(hindlimb_lateral,0.7*hindlimb_lateral,scale=100)->hindlimbV
gdi(forelimb_lateral, 0.7*forelimb_lateral, scale=100)->forelimbV
gdi(measurements_lateral, measurements_dorsal, scale=100, method="raw")->axialV

## -----------------------------------------------------------------------------
volumes<-c(axialV, forelimbV, hindlimbV)
x_COM<-c(hCOM(gdiresults), vCOM(fl), vCOM(hl))
y_COM<-c(vCOM(gdiresults), hCOM(fl, align="v"), hCOM(hl,align="v"))
knitr::kable(data.frame(volumes,x_COM,y_COM, row.names=c("axial_total","forelimb","hindlimb")))

## -----------------------------------------------------------------------------
weighted.mean(x_COM,w=volumes*c(1,2,2))#horizontal COM position
weighted.mean(y_COM,w=volumes*c(1,2,2))#vertical COM position

## ---- results=FALSE, message=FALSE, fig.cap="Plotted silhouette with overall COM position highlighted", fig.height=5, fig.width=7----
plot_sil(measurements_lateral, asp=1)
plot_sil(forelimb_lateral, flip=T, add=T)#for silhouettes that were aligned, but digitized with align="v", we need to specify flip=TRUE here
plot_sil(hindlimb_lateral, flip=T, add=T)

points(weighted.mean(x_COM,w=volumes*c(1,2,2)), weighted.mean(y_COM,w=volumes*c(1,2,2)), pch=16)
abline(v=weighted.mean(x_COM,w=volumes*c(1,2,2)), lty=2)
abline(h=weighted.mean(y_COM,w=volumes*c(1,2,2)), lty=2)

