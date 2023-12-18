# Find likelihood of calendar ages t for a single 14Cage observation
# Arguments
# x         - c14age of single observation
# xsig      - c14sig of single observation
# t         - vector of calendar ages tofind likelihood
# calcurve  - calibration curve
LegacyCalibrate <- function(x, xsig, t, calcurve) {
  calmu <- approx(calcurve$calage, calcurve$c14age, t)$y # Estimated calibration curve mean at calage t
  calsig <- approx(calcurve$calage, calcurve$c14sig, t)$y # Estimated calibratio curve sd at calage t
  lik <- dnorm(x, calmu, sqrt(calsig^2 + xsig^2)) # Likelihood of x at all given calendar ages
  # lik <- approx(calcurve$calage, dnorm(x, calcurve$c14age, sqrt(xsig^2 + calcurve$c14sig^2)), t, rule=2)$y
  return(lik)
}
