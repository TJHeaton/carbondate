## ── R CMD check results ─────────────────────────────────────────── carbondate 1.0.1.9000 ────
Duration: 9m 34.9s

0 errors | 0 warnings | 1 note

❯ checking CRAN incoming feasibility ... [14s] NOTE
  Maintainer: 'Timothy J Heaton <T.Heaton@leeds.ac.uk>'
  
  Version contains large components (1.0.1.9000)

* This is default naming for development version


# Old versions

## Resubmission 1.0.1
This is a resubmission. In this version in response to your feedback, we have:

* The writing, i.e., save(), commands in the tests/testthat/fixtures are in scripts that are not run as part of the tests. These were scripts to create static test fixtures where we aimed to follow the guidance on R packages (which said that the code should be available). We have commented out the save commands for safety.  

* Reset par back to user values when exiting ALL functions using suggested approach of:
oldpar <- par(no.readonly = TRUE)    # code line i
on.exit(par(oldpar))            # code line i + 1
The only instance where we have not done this is in the function .SetUpDensityPlot (in R/plot_helpers.R) which just calls par(new = TRUE). This function is a specific "helper" call to set up an overlay plot with a different scale. It is an internal/hidden function, and presumably anyone calling it would see no difference from if we had called any plot command in the function.

* Reset par back to user values in all vignettes/examples at end of each Rmarkdown chunk

We have also:

* Altered documentation of functions (and vignettes) to incorporate above plotting changes (and fix a few typos to ensure parameter naming consistency across outputs).

* Fixed typos in vignettes

* Run MCMC in "Non-parametric-summed-density" for longer to show converged result, shortened other vignettes to get under 10min runtime

* Switched antialiaising back to default on png on vignettes.


## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.
