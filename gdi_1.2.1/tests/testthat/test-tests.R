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
expect_equal(signif(gdi(c(1), c(1), scale=1),7), signif(0.7853982,7))

})

