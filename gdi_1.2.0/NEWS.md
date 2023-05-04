gdi v1.2.0 (Release date: 2023-05-04)
==============

Changes:
* added functions imghist() and fdetect() for aiding with setting appropriate parameters for measuresil() and cscorr()
* added function add.alpha() for adding transparency to a colour value when plotting.


gdi v1.1.2 (Release date: 2023-05-04)
==============

Changes:
* Fixed missing "," in description.
* Fixed errors in and improved vignette.
* Fixed error in cscorr()-function. Previous version would output the diameters of the silhouette in pixels, even if a scale was provided. New version converts them to units determined by scale.


gdi v1.1.1 (Release date: 2023-05-04)
==============

Changes:
* Fixed invalid doi in the description.


gdi v1.1.0 (Release date: 2023-05-03)
==============

Changes:
* Updated the description, vignettes and documentation.

* Added an additional parameter (corr) to the function gdi(). Similar to the parameter k, corr accounts for cross-sectional shapes other than simple ellipses. Differing from k, this is not done by calculating the area of a superellipse with a given exponent, but by applying a correction factor specifying the ratio in cross-sectional areas (such as the one calculated by cscorr()).

* Added the function cscorr(), which reads images of body cross-sections and computes various metrics, including a correction factor based on the ratio in area between the measured cross-section and an ellipse with the same diameters. This function enables the user to account for a wider variety of cross-sectional geometries than simple ellipses or superellipses, provided that a graphical representation of the cross-section is available.


gdi v1.0.0 (Release date: 2023-04-30)
==============

Changes:

* First version of this package.
 



