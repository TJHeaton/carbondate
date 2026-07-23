# Package index

## Calibration of a Single Radiocarbon Sample

Function to independently calibrate a single radiocarbon determination
against a given calibration curve. The rest of the library concerns
calibration of multiple (related) radiocarbon samples.

- [`CalibrateSingleDetermination()`](https://tjheaton.github.io/carbondate/reference/CalibrateSingleDetermination.md)
  : Calibrate a Single Radiocarbon Determination

## Bayesian Non-Parametric Joint Density Estimation

Functions providing a rigorous Bayesian non-parametric alternative to
summed probability distributions (SPDs). These functions enable the
calibration and calendar age summarisation of multiple related
radiocarbon samples (and provide plotting of the results)

- [`PolyaUrnBivarDirichlet()`](https://tjheaton.github.io/carbondate/reference/PolyaUrnBivarDirichlet.md)
  : Calibrate and Summarise Multiple Radiocarbon Samples via a Bayesian
  Non-Parametric DPMM (with Polya Urn Updating)
- [`WalkerBivarDirichlet()`](https://tjheaton.github.io/carbondate/reference/WalkerBivarDirichlet.md)
  : Calibrate and Summarise Multiple Radiocarbon Samples via a Bayesian
  Non-Parametric DPMM (with Walker Updating)
- [`FindPredictiveCalendarAgeDensity()`](https://tjheaton.github.io/carbondate/reference/FindPredictiveCalendarAgeDensity.md)
  : Find Predictive Estimate of Shared Calendar Age Density from
  Bayesian Non-Parametric DPMM Output
- [`PlotPredictiveCalendarAgeDensity()`](https://tjheaton.github.io/carbondate/reference/PlotPredictiveCalendarAgeDensity.md)
  : Plot Predictive Estimate of Shared Calendar Age Density from
  Bayesian Non-Parametric DPMM Output
- [`PlotNumberOfClusters()`](https://tjheaton.github.io/carbondate/reference/PlotNumberOfClusters.md)
  : Plot Number of Calendar Age Clusters Estimated in Bayesian
  Non-Parametric DPMM Output
- [`PlotCalendarAgeDensityIndividualSample()`](https://tjheaton.github.io/carbondate/reference/PlotCalendarAgeDensityIndividualSample.md)
  : Plot Posterior Calendar Age Estimate for an Individual Determination
  after Joint Calibration
- [`PlotConvergenceData()`](https://tjheaton.github.io/carbondate/reference/PlotConvergenceData.md)
  : Plot KL Divergence of Predictive Density to Assess Convergence of
  Bayesian Non-Parametric DPMM Sampler

## Poisson Process Modelling

Functions for modelling the occurrence of radiocarbon samples as a
variable-rate (inhomogeneous) Poisson process. This is a further
(analogous) approach that enables the rigorous summarisation of calendar
age information from multiple related radiocarbon samples.

- [`PPcalibrate()`](https://tjheaton.github.io/carbondate/reference/PPcalibrate.md)
  : Model Occurrence of Multiple Radiocarbon Samples as a Variable-Rate
  Poisson Process
- [`PPcalibrateMixedCurves()`](https://tjheaton.github.io/carbondate/reference/PPcalibrateMixedCurves.md)
  : Model Occurrence of Multiple Radiocarbon Samples as a Variable-Rate
  Poisson Process when the Samples Come from a Variety of Environments
  and Multiple Calibration Curves are Needed
- [`FindPosteriorMeanRate()`](https://tjheaton.github.io/carbondate/reference/FindPosteriorMeanRate.md)
  : Find Posterior Mean Rate of Sample Occurrence for Poisson Process
  Model
- [`PlotPosteriorMeanRate()`](https://tjheaton.github.io/carbondate/reference/PlotPosteriorMeanRate.md)
  : Plot Posterior Mean Rate of Sample Occurrence for Poisson Process
  Model
- [`PlotNumberOfInternalChanges()`](https://tjheaton.github.io/carbondate/reference/PlotNumberOfInternalChanges.md)
  : Plot Number of Changepoints in Rate of Sample Occurrence for Poisson
  Process Model
- [`PlotPosteriorChangePoints()`](https://tjheaton.github.io/carbondate/reference/PlotPosteriorChangePoints.md)
  : Plot Calendar Ages of Changes in Rate of Sample Occurrence for
  Poisson Process Model
- [`PlotPosteriorHeights()`](https://tjheaton.github.io/carbondate/reference/PlotPosteriorHeights.md)
  : Plot Heights of Segments in Rate of Sample Occurrence for Poisson
  Process Model
- [`PlotCalendarAgeDensityIndividualSample()`](https://tjheaton.github.io/carbondate/reference/PlotCalendarAgeDensityIndividualSample.md)
  : Plot Posterior Calendar Age Estimate for an Individual Determination
  after Joint Calibration
- [`PlotRateIndividualRealisation()`](https://tjheaton.github.io/carbondate/reference/PlotRateIndividualRealisation.md)
  : Plot Individual Realisations of Posterior Rate of Sample Occurrence
  for Poisson Process Model

## Annotating the Summary Plots

Functions that can be used to annotate the various summary plots.

- [`AddTextPlot()`](https://tjheaton.github.io/carbondate/reference/AddTextPlot.md)
  : Add Text Annotation to Various Summary Plots
- [`AddLinePlot()`](https://tjheaton.github.io/carbondate/reference/AddLinePlot.md)
  : Add Straight Lines to Various Summary Plots
- [`AddShadingPlot()`](https://tjheaton.github.io/carbondate/reference/AddShadingPlot.md)
  : Add Shading to Various Summary Plots

## Common convergence plotting functions

Functions that can be used to plot information on convergence of the
MCMC for both the Bayesian non-parametric DPMM and the Poisson process
approaches to calibration and summarisation.

- [`PlotGelmanRubinDiagnosticSingleChain()`](https://tjheaton.github.io/carbondate/reference/PlotGelmanRubinDiagnosticSingleChain.md)
  : Plot Histogram of the Gelman-Rubin Convergence Diagnostic for a
  Single MCMC Chain
- [`PlotGelmanRubinDiagnosticMultiChain()`](https://tjheaton.github.io/carbondate/reference/PlotGelmanRubinDiagnosticMultiChain.md)
  : Plot Histogram of the Gelman-Rubin Convergence Diagnostic for
  Multiple Independent MCMC Chains

## Helper functions

General functions that may be of use. However, please do not use SPDs in
your inference (use either the Bayesian non-parametric or Poisson
Process approaches described above).

- [`FindSummedProbabilityDistribution()`](https://tjheaton.github.io/carbondate/reference/FindSummedProbabilityDistribution.md)
  : Find the summed probability distribution (SPD) for a set of
  radiocarbon observations
- [`InterpolateCalibrationCurve()`](https://tjheaton.github.io/carbondate/reference/InterpolateCalibrationCurve.md)
  : Interpolate a calibration curve at a set of calendar ages
- [`GenerateOxcalCode()`](https://tjheaton.github.io/carbondate/reference/GenerateOxcalCode.md)
  : Outputs code suitable for running in OxCal from a series of
  radiocarbon determinations

## Example datasets

Sets of radiocarbon determinations, either from real-life examples or
artificially-generated, that can be used as examples with the
calibration functions.

- [`alces`](https://tjheaton.github.io/carbondate/reference/alces.md) :
  Example real-life data - Alces in Yukon and Alaska
- [`armit`](https://tjheaton.github.io/carbondate/reference/armit.md) :
  Example real-life data - Population Decline in Iron Age Ireland
- [`bison`](https://tjheaton.github.io/carbondate/reference/bison.md) :
  Example real-life data - Bison in Yukon and Alaska
- [`buchanan`](https://tjheaton.github.io/carbondate/reference/buchanan.md)
  : Example real-life data - Palaeo-Indian demography
- [`cervus`](https://tjheaton.github.io/carbondate/reference/cervus.md)
  : Example real-life data - Cervus in Yukon and Alaska
- [`equus`](https://tjheaton.github.io/carbondate/reference/equus.md) :
  Example real-life data - Equus in Yukon and Alaska
- [`human`](https://tjheaton.github.io/carbondate/reference/human.md) :
  Example real-life data - Humans in Yukon and Alaska
- [`kerr`](https://tjheaton.github.io/carbondate/reference/kerr.md) :
  Example real-life data - Irish Rath
- [`mammuthus`](https://tjheaton.github.io/carbondate/reference/mammuthus.md)
  : Example real-life data - Mammuthus in Yukon and Alaska
- [`pp_uniform_phase`](https://tjheaton.github.io/carbondate/reference/pp_uniform_phase.md)
  : Example artificial data - Uniform Phase
- [`pp_uniform_phase_marine`](https://tjheaton.github.io/carbondate/reference/pp_uniform_phase_marine.md)
  : Example artificial marine data - Uniform Phase
- [`pp_uniform_phase_mixed`](https://tjheaton.github.io/carbondate/reference/pp_uniform_phase_mixed.md)
  : Example artificial data requiring multiple calibration curves -
  Uniform Phase
- [`two_normals`](https://tjheaton.github.io/carbondate/reference/two_normals.md)
  : Example artificial data - Mixture of Normal Phases
- [`two_normals_marine`](https://tjheaton.github.io/carbondate/reference/two_normals_marine.md)
  : Example artificial data - Mixture of Normal Phases Using Marine20
  calibration curve

## Calibration curves

- [`intcal04`](https://tjheaton.github.io/carbondate/reference/intcal04.md)
  : IntCal04 calibration curve
- [`intcal09`](https://tjheaton.github.io/carbondate/reference/intcal09.md)
  : IntCal09 calibration curve
- [`intcal13`](https://tjheaton.github.io/carbondate/reference/intcal13.md)
  : IntCal13 calibration curve
- [`intcal20`](https://tjheaton.github.io/carbondate/reference/intcal20.md)
  : IntCal20 calibration curve
- [`intcal98`](https://tjheaton.github.io/carbondate/reference/intcal98.md)
  : IntCal98 calibration curve
- [`shcal04`](https://tjheaton.github.io/carbondate/reference/shcal04.md)
  : SHCal04 calibration curve
- [`shcal13`](https://tjheaton.github.io/carbondate/reference/shcal13.md)
  : SHCal13 calibration curve
- [`shcal20`](https://tjheaton.github.io/carbondate/reference/shcal20.md)
  : SHCal20 calibration curve
- [`marine04`](https://tjheaton.github.io/carbondate/reference/marine04.md)
  : Marine04 calibration curve
- [`marine09`](https://tjheaton.github.io/carbondate/reference/marine09.md)
  : Marine09 calibration curve
- [`marine13`](https://tjheaton.github.io/carbondate/reference/marine13.md)
  : Marine13 calibration curve
- [`marine20`](https://tjheaton.github.io/carbondate/reference/marine20.md)
  : Marine20 calibration curve
