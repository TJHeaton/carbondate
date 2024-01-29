# In this file we will write code which
# finds the predictive dsitribution for a new theta given
# the existing thetas and mu

Temp <- NealTemp

# Find predictive density
npost <- dim(Temp$theta)[1]
nburn <- floor(npost/2)
SPDcol <- grey(0.1, alpha = 0.5)

# Create the range over which to plot the density (note we will have to rescale the pdfs to account for calculation grid)
ncalc <- 1001
tempx <- seq(floor(min(Temp$theta, na.rm = TRUE)), ceiling(max(Temp$theta, na.rm = TRUE)), length = ncalc)

#Now choose the ids of the posterior sample
sampid <- sample(x = nburn:npost, size = npostsum, replace = npostsum > (npost-nburn))

# Create a matrix where each column is the density for a particular sample id
# We can then find the mean along each row
postDmat <- apply(as.row(sampid), 2, function(i, out, x, lambda, nu1, nu2)
  NealFindpred(x, c = out$c[i,], phi = out$phi[[i]], tau = out$tau[[i]],
           alpha = out$alpha[i], muphi = out$muphi[i],
           lambda = lambda, nu1 = nu1, nu2 = nu2),
  out = Temp, x = tempx, lambda = lambda, nu1 = nu1, nu2 = nu2)

# Find CI and mean
postdenCI <- apply(postDmat, 1, quantile, probs = c(0.025, 0.975))
postden <- apply(postDmat, 1, mean)

# Create a layout which has 2/3 of the plot showing the predictive density and 1/3 showing the number of clusters
layout.matrix <- matrix(c(1, 2), nrow = 1, ncol = 2)
layout(mat = layout.matrix,
       heights = 1, # Heights of the two rows
       widths = c(10, 4.5)) # Widths of the two columns

# Plot the predictive joint density
par(mar = c(5, 4.5, 4, 2) + 0.1, las = 1)
plot(calcurve$calage, calcurve$c14age, col = "blue",
     ylim = range(x) + c(-2,2) * max(xsig), xlim = rev(range(tempx)),
     xlab = "Calendar Age (cal yr BP)", ylab = expression(paste(""^14, "C", " age (", ""^14, "C yr BP)")),
     type = "l", main = expression(paste(""^14,"C Calibration - Neal/P\u{F3}lya Urn DP")))
calcurve$ub <- calcurve$c14age + 1.96 * calcurve$c14sig
calcurve$lb <- calcurve$c14age - 1.96 * calcurve$c14sig
lines(calcurve$calage, calcurve$ub, lty = 3, col = "blue" )
lines(calcurve$calage, calcurve$lb, lty = 3, col = "blue")
polygon(c(rev(calcurve$calage), calcurve$calage), c(rev(calcurve$lb), calcurve$ub), col=rgb(0,0,1,.3), border=NA)
rug(x, side = 2)

legend("topright", legend = c("IntCal20", "NP Density Estimate - Neal", "95% Prob. Interval", "SPD Estimate"),
       lty = c(1,1,2, -1), pch = c(NA, NA, NA, 15),
       col = c("blue", "purple", "black", SPDcol), cex = 0.8, pt.cex = 2)

# Plot the density along the bottom
par(new = TRUE)
plot(tempx, postden, lty = 1, col = "purple", type = "n",
     ylim = c(0, 2*max(postdenCI)), xlim = rev(range(tempx)),
     axes = FALSE, xlab = NA, ylab = NA)
polygon(c(SPD$calage, rev(SPD$calage)), c(SPD$prob, rep(0, length(SPD$prob))),
        border = NA, col = SPDcol)
lines(tempx, postden, col = "purple")
lines(tempx, postdenCI[1,], col = "purple", lty = 3)
lines(tempx, postdenCI[2,], col = "purple", lty = 3)
#axis(side = 4, at = axTicks(4)/2, col = "purple", col.ticks = "purple")
mtext(paste0("(", letters[1], ")"), side = 3, adj = 0.05,
      line = -1.1)

# Plot the number of clusters
NealNClust <- apply(NealTemp$c, 1, max)
NealNClust <- NealNClust[nburn:npost]
hist(NealNClust, xlab = "Number of Clusters", main = "",
     probability = TRUE, breaks = seq(0.5, max(NealNClust)+1, by = 1))
mtext(paste0("(", letters[2], ")"), side = 3, adj = 1,
      line = -1.1)
