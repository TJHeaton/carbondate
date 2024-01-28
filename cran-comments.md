## Resubmission
This is a resubmission. In this version I have:

* Reset ALL par back to user values when exiting all functions using suggested approach of:
oldpar <- par(no.readonly = TRUE)    # code line i
on.exit(par(oldpar))            # code line i + 1
The only instance where we have not done this is in the function .SetUpDensityPlot (in R/plot_helpers.R) which just calls par(new). This function is a specific "helper" call to set up an overlay plot with a different scale. It is an internal/hidden function, and presumably anyone calling it would see no difference from if we had called any plot command in the function.

* Reset par back to user values in all vignettes/examples 





## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.
