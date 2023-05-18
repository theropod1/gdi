test_that("gdi() works", {

expect_equal(fdetect(system.file("exdata","lat.png", package="gdi"))$most_common, 0)

})



test_that("cscorr() works", {

#test that area correction factor is calculated correctly
expect_equal(signif(cscorr(system.file("exdata","cross_section.png", package="gdi")),7), signif(1.092215,7))

#test that area is calculated correctly
expect_equal(cscorr(system.file("exdata","cross_section.png", package="gdi"), return="area"), 119482)
  
})

test_that("measuresil() works", {
lat<-measuresil(system.file("exdata","lat.png", package="gdi"))

expect_equal(signif(lat[250],7), signif(69,7))
expect_equal(signif(lat[400],7), signif(83,7))

})

test_that("gdi() works", {

result <- gdi(c(1), c(1), scale=1)
exp <- 0.7853982
names(exp)<-paste("x_dim",sum(c(1)!=0)/1, "units", sep="_")

expect_equal(signif(result,7), signif(exp,7))

})

