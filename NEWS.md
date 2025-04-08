# carbondate (development version)

* Fixed warning in cpp using return std::move(retdata) where retdata is just named value

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
