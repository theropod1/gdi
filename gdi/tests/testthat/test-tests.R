test_that("fdetect() works", {

expect_equal(fdetect(system.file("exdata","lat.png", package="gdi"))$most_common, 0)

})



test_that("cscorr() works", {

#test that area correction factor is calculated correctly
expect_equal(signif(cscorr(system.file("exdata","cross_section.png", package="gdi")),7), signif(1.092215,7))

#test that area is calculated correctly
expect_equal(cscorr(system.file("exdata","cross_section.png", package="gdi"), return="area"), 119482)

#test that scale is applied correctly
expect_equal(cscorr(system.file("exdata","cross_section.png", package="gdi"), return="area", scale=10), 1194.82)
expect_equal(as.numeric(cscorr(system.file("exdata","cross_section.png", package="gdi"), return="diameters")[1]), 313)
expect_equal(as.numeric(cscorr(system.file("exdata","cross_section.png", package="gdi"), return="diameters", scale=10)[1]), 31.3)
expect_equal(as.numeric(signif(cscorr(system.file("exdata","cross_section.png", package="gdi"), return="diameters", scale=15)[1],4)), signif(20.8666667,4))

expect_equal(as.numeric(signif(cscorr(system.file("exdata","cross_section.png", package="gdi"), return="all", scale=10)[5],4)), as.numeric(signif(cscorr(system.file("exdata","cross_section.png", package="gdi"), return="area_corr")[1],4)))



  
})

test_that("measuresil() works", {
lat<-measuresil(system.file("exdata","lat.png", package="gdi"))
dors<-measuresil(system.file("exdata","dors.png", package="gdi"))

expect_equal(signif(lat[250],7), signif(68,7))
expect_equal(signif(lat[400],7), signif(82,7))

})

test_that("gdi() works", {

result <- gdi(c(1), c(1), scale=1)
exp <- 0.7853982
names(exp)<-paste("x_dim",sum(c(1)!=0)/1, "units", sep="_")

expect_equal(signif(result,7), signif(exp,7))

lat<-measuresil(system.file("exdata","lat.png", package="gdi"))
dors<-measuresil(system.file("exdata","dors.png", package="gdi"))

expect_equal(signif(as.numeric(gdi(lat,dors, scale=10)),7), signif(40256.77,7))



})

