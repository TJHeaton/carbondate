# carbondate (development version)

* Added PPcalibrateMixedCurves to allow Poisson process summarisation of sets of sampels from a range of NH, SH and Marine environments
* Added marine calibration curves to data (Marine20, Marine13, Marine09, and Marine04)
* Added functions AddTextPlot(), AddLinePlot(), and AddShadingPlot() to enable users to overlay text and shading onto summary plots. This required edits to the summary plotting functions so they returned a list, with plotting par being the 2nd list component.  
* Fixed warning in find_predictive_density.cpp by using return std::move(retdata) where retdata was just named local variable
* Added function to create text at a specific point on the plots  

# carbondate 1.1.0

* Added `plot_lwd` argument to alter width of lines when plotting PP and DPMM. 
* Added `PlotRateIndividualRealisation()` function to plot individual posterior realisations of the Poisson Process rate.
* Edited `PlotPosteriorMeanRate()` function so users can plot posterior mean over time conditioned on a specified number of changes in the rate. We do not recommend this unless the user has independent knowledge of the number of changes. The conditioning is done by specifying the `n_changes` argument. The recommended default is NULL (no conditioning). 
* Updated website to explain can access posterior directly. 

# carbondate 1.0.1

* Fixing CRAN comments. 

# carbondate 1.0.0

* First release.

* Added a `NEWS.md` file to track changes to the package.
