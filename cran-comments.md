## R CMD check results

0 errors | 0 warnings | 0 note

* This is a resubmission.

## RHub CMD check results

* No errors, warnings or notes on platform "macos-highsierra-release-cran"

* 1 error on platform "macos-m1-bigsur-release":
Error: processing vignette 'Vignette.rmd' failed with diagnostics:
   polygon edge not found
   --- failed re-building ‘Vignette.rmd’
   
   SUMMARY: processing the following file failed:
     ‘Vignette.rmd’
     
* This error is reportedly due to a bug in the grid package which affects the
ggplot2 package on MacOS-arm64
