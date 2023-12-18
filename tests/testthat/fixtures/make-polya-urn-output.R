# 13th April 2022
# Code to implement NP Bayesian Calibration and Summaristion of Related 14C Determinations
# Original code written by TJHeaton and available at https://github.com/TJHeaton/NonparametricCalibration
# This code is adaped from NPCalibrationNeal to test that the R Package is giving the correct result

# NOTE:
# This version uses the Polya Urn method of sampling form the DP
# It is slower so I suggest using the Walker version

seednum <- 8
set.seed(seednum)

# Read in IntCal20 curve
calcurve <- read.table("tests/testthat/fixtures/helpers/intcal20.14c", sep = ",", header=FALSE, skip=11)
names(calcurve) <- c("calage", "c14age", "c14sig", "Delta14C", "DeltaSigma")

# Read in the necessary functions
# Read in the required slice updating functions
source("tests/testthat/fixtures/helpers/SliceUpdateFunsFinal.R")
source('tests/testthat/fixtures/helpers/WalkerDirichletMixtureUpdateFunsFinal.R') # This also reads in the slice sampling SliceUpdateFuns.R
source('tests/testthat/fixtures/helpers/NealDirichletMixtureMasterFunctionsFinal.R')
source("tests/testthat/fixtures/helpers/WalkerMasterFunctionFinal.R")
source("tests/testthat/fixtures/helpers/SimStudyFuncsFinal.R")

# Read in data
# x - c14ages
# xsig - corresponding 1 sigma
Kerr <- read.csv("tests/testthat/fixtures/helpers/kerr2014sss_sup.csv", header = FALSE, sep =  ",")
x <- Kerr[,3]
xsig <- Kerr[,4]

# Only choose the first 100 points for speed
x <- x[1:100]
xsig <- xsig[1:100]
#############################################################################
# Now choose hyperparameters
############################################################
# Prior on the concentration parameter
# Place  a gamma prior on alpha
# alpha ~ Gamma(alphaprshape, alphaprrate)
# A small alpha means more concentrated (i.e. few clusters)
# Large alpha not concentrated (many clusters)
cprshape <- alphaprshape <- 1
cprrate <- alphaprrate <- 1

#### Updated adaptive version
# Prior on mu theta for DP - very uninformative based on observed data
IntCalyrgrid <- FindCal(1:50000, calcurve$c14ag, calcurve$calage, calcurve$c14sig)
initprobs <- mapply(calibind, x, xsig, MoreArgs = list(calmu = IntCalyrgrid$mu, calsig = IntCalyrgrid$sigma))
inittheta <- apply(initprobs, 2, which.max)
# Choose A and B from range of theta
A <- median(inittheta)
B <- 1 / (max(inittheta) - min(inittheta))^2
maxrange <- max(inittheta) - min(inittheta)

# Parameters for sigma2 (sigma^2 ~ InvGamma(nu1, nu2))
# E[tau] = (1/100)^2 Var[tau] = (1/100)^4
# Interval for sigma2 is approx 1/ c(nu2/nu1 - 2*nu2^2/nu1, nu2/nu1 + 2*nu2^2/nu1)
tempspread <- 0.1 * mad(inittheta)
tempprec <- 1/(tempspread)^2
nu1 <- 0.25
nu2 <- nu1 / tempprec

# Setup the NP method
lambda <- (100/maxrange)^2  # Each muclust ~ N(mutheta, sigma2/lambda)


# Choose number of iterations for sampler
niter <- 1e5
nthin <- 5 # Don't choose too high, after burn-in we have (niter/nthin)/2 samples from posterior to potentially use
npostsum <- 500 # Current number of samples it will draw from this posterior to estimate fhat (possibly repeats)

mualpha <- NA
sigalpha <- NA
w <- max(1000, diff(range(x))/2)
m <- 10
nclusinit <- 10

save(x, xsig, lambda, nu1, nu2, A, B, mualpha, sigalpha, alphaprshape, alphaprrate, niter, nthin,
     inittheta, w, m, nclusinit, seednum, file="tests/testthat/fixtures/polya_urn_input.rda")

NealTemp <- BivarGibbsDirichletwithSlice(x = x, xsig = xsig,
                                         lambda = lambda, nu1 = nu1, nu2 = nu2,
                                         A = A, B = B,
                                         mualpha = mualpha, sigalpha = sigalpha,
                                         alphaprshape = alphaprshape, alphaprrate = alphaprrate,
                                         niter = niter, nthin = nthin, theta = inittheta,
                                         w = w, m = m,
                                         calcurve = calcurve, nclusinit = nclusinit)

##############################
# Also find the SPD estimate to plot alongside
##############################
# Find the independent calibration probabilities
yrange <- floor(range(NealTemp$theta))
yfromto <- seq(max(0,yrange[1]-400), min(50000, yrange[2]+400), by = 1)

# Find the calibration curve mean and sd over the yrange
CurveR <- FindCalCurve(yfromto, calcurve)

# Now we want to apply to each radiocarbon determination
# Matrix where each column represents the posterior probability of each theta in yfromto
indprobs <- mapply(calibind, x, xsig, MoreArgs = list(calmu = CurveR$curvemean, calsig = CurveR$curvesd))

# Find the SPD estimate (save as dataframe)
SPD <- data.frame(calage = yfromto, prob = apply(indprobs, 1, sum)/dim(indprobs)[2])

# To plot the predictive distribution then you run
seednum <- 14
set.seed(seednum)
source("tests/testthat/fixtures/helpers/NealPostProcessingFinal.R")

save(x, xsig, postdenCI, postden, tempx, file = "tests/testthat/fixtures/polya_urn_output.rda")

# To access the posterior calendar age estimate for individual determination then you can look at:
# NealTemp$theta[,10] # MCMC chain for 10th determination (will need to remove burn in)

# If we want to plot e.g. the posterior calendar age density against the curve then we can run the below
# ident is the determination you want to calibrate
ident <- 15
resolution <- 10
indpost <- plotindpost(NealTemp, ident = ident, y = x, er = xsig, calcurve = calcurve, resolution = resolution)
