gdi v1.9.0 (Release date: 2025-09-15)
==============
* added transfer_ratio() function for interpolating missing dorsal or lateral view silhouettes from a template (=specimen with known dorsal and lateral view silhouette from which

gdi v1.8.1 (Release date: 2025-07-25)
==============
* fixed expandConvexHull() function; regressions are log-transformed

gdi v1.8.0 (Release date: 2025-05-19)
==============
Changes:
* Added function expandConvexHull, to generate soft tissue expansion factors based on regressions and data in Macaulay et al. 2023
* Added expansion parameter to plot_sil, to allow the visualization of soft tissue expansion factors.

gdi v1.7.0 (Release date: 2024-10-18)
==============
Changes:
* fixed handling of segment length input in gdi(): scale factor is now also applied to manual segment lengths internally, rather than requiring the segment lengths to already be given to the correct scale.


gdi v1.6.0 (Release date: 2023-09-23)
==============

Changes:
* cosmetic improvements to plot_sil function, which now plots polygons instead of lines and allows to specify a scale
* Increased size of plots in vignettes for better visibility

gdi v1.5.5 (Release date: 2023-09-23)
==============

Changes:
* csI() now outputs centroid and bounding box coordinates as attributes


gdi v1.5.4 (Release date: 2023-09-06)
==============

Changes:
* added new tests and vignette


gdi v1.5.3 (Release date: 2023-08-25)
==============

Changes:
* added collumn with x axis centroids to output of gdi(), for use in rotI() and hCOM() functions

gdi v1.5.2 (Release date: 2023-08-24)
==============

Changes:
* generalized rotI() function to calculate rotational inertia of shapes around any axis
* added function to cscorr to allow calculation of correction factors for moment of inertia.
* fixed issue in input to functions (rejection of previously read images supplied as matrix or array)


gdi v1.5.1 (Release date: 2023-08-23)
==============

Changes:
* added rotI() function to calculate rotational inertia of shapes


gdi v1.5.0 (Release date: 2023-08-19)
==============

Changes:
* added csI() function to calculate second moment of area for cross-sectional shapes



gdi v1.4.2 (Release date: 2023-08-08)
==============

Changes:
* fixed bug in cscorr (returning inverse aspect ratios depending on settings) – now aspect ratio is always max horizontal diameter/max vertical diameter


gdi v1.4.1 (Release date: 2023-07-28)
==============

Changes:
* modifications to image-processing functions to deal with single-channel images


gdi v1.4.0 (Release date: 2023-07-28)
==============

Changes:
* added functions vCOM and hCOM and expanded functionality of measuresil() and gdi() to allow estimating the position of the center of mass (COM)
* added function plot_sil() to plot silhouettes



gdi v1.3.1 (Release date: 2023-07-22)
==============

Changes:
* added parameter setting indices for gdi(), allowing to specify optional indices to specify a subset of the measurement vectors


gdi v1.3.0 (Release date: 2023-06-12)
==============

Changes:
* added parameter setting alignment for measuresil(), allowing the optional use of vertically aligned silhouettes if align!="h"
* updated unit tests to check that scale works correctly


gdi v1.2.3 (Release date: 2023-06-12)
==============

Changes:
* fixed cscorr() area output scaling (previous version didn’t set the scale for area resulting in unscaled area output even when scale is set)


gdi v1.2.2 (Release date: 2023-06-07)
==============

Changes:
* added function sellipse.coo() to visualize shapes of superelliptical cross-sections

gdi v1.2.1 (Release date: 2023-05-05)
==============

Changes:
* added tests using package 'testthat'
* gdi() now also reports the horizontal length of the silhouette as name of output (in the same unit cubed to calculate volume, i.e. decimeters if volume is intended in liters), intended to quickly spot scaling errors

gdi v1.2.0 (Release date: 2023-05-04)
==============

Changes:
* added functions imghist() and fdetect() for aiding with setting appropriate threshold and method parameters for measuresil() and cscorr()


gdi v1.1.2 (Release date: 2023-05-04)
==============

Changes:
* Fixed missing "," in description.
* Fixed errors in and improved vignette.
* Fixed error in cscorr()-function. Previous version would output the diameters of the silhouette in pixels, even if a scale was provided. New version converts them to units determined by scale parameter.


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
 



