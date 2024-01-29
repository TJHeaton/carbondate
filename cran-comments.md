## Resubmission
This is a resubmission. In this version I have:

* The writing, i.e., save(), commands in the tests/testthat/fixtures are in scripts that are not run as part of the tests. These were scripts to create static test fixtures where we aimed to follow the guidance on R packages (which said that the code should be available). We have commented out the save commands for safety.  

* Reset par back to user values when exiting ALL functions using suggested approach of:
oldpar <- par(no.readonly = TRUE)    # code line i
on.exit(par(oldpar))            # code line i + 1
The only instance where we have not done this is in the function .SetUpDensityPlot (in R/plot_helpers.R) which just calls par(new = TRUE). This function is a specific "helper" call to set up an overlay plot with a different scale. It is an internal/hidden function, and presumably anyone calling it would see no difference from if we had called any plot command in the function.

* Reset par back to user values in all vignettes/examples at end of each Rmarkdown chunk

* Altered documentation of functions (and vignettes) to incorporate above plotting changes.

* Fixed typos in vignettes and run MCMC in "Non-parametric-summed-density" for longer to show converged result


## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.
