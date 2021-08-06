#' Wrapper function to calculate negative-sum of log-likelihoods
#' of simulated data against observed data
#'
#' @param obsData a wide-format data.table of observed species biomass or counts
#'  averaged/summed across a landscape (one line per year per repetition) or per \code{pixelIndex}.
#'  Should also have the columns \code{rep} and \code{year} and, optionally, \code{pixelIndex}.
#'  Note that usually observed data has the same values across repetitions.
#' @param simData as \code{obsData}, but with simulated/predicted data
#' @param reps vector of repetitions (i.e. \code{unique(simData$reps)})
#' @param years vector of years (i.e. \code{unique(simData$years)})
#' @param cols character. columns with data that will enter the log-likelihood simulation
#'  (usually species). Defaults to NULL and uses all columns that are not \code{rep}, \code{year}
#'  or \code{pixelIndex}.
#' @param varType what type of variable is being accessed? Only accepts "biomass" or "counts". See details.
#'
#' @return an \code{array} of negative-sum of log-likelihoods for each repetition and year.
#'
#' @details Negative sum of log-likelihoods are calculated by drawing from a density function of a known
#' distribution. For both "biomass" and "counts" (i.e. species presences), the only two cases implemented so far.
#' we use the multinomial distribution (\code{dmultinom}), as describes the probability of a given set of values (biomass or presences) being
#' attributes to different classes (the species). `simData` values are used as "probabilities" of observing
#' \code{x} (i.e. the \code{obsData} values) (see \code{?dmultinom}) and the output probabilities (i.e. likelihoods)
#' are logged. When probabilities are extremely small/large the calculated log-likelihood will be
#'  \code{abs(Inf)}. In this case, we take the smallest value of resulting log-likelihoods (across pixels only -
#'  at the landscape level it will be left as \code{Inf}) before calculating the \code{-sum} across pixels.

NegSumLogLikWrapper <- function(obsData, simData, reps, years, cols = NULL, varType = "biomass") {
  NegSumLogLik <- array(data = NA, dim = c(length(reps), length(years)),
                   dimnames = list(reps = reps, years = years))
  ## checks
  if (length(setdiff(names(obsData), names(simData)))) {
    stop("Observed and simulated data have different columns")
  }

  if (is.null(cols)) {
    cols <- setdiff(names(simData), c("rep", "year", "pixelIndex"))
  }

  cols2 <- if (any(grepl("pixelIndex", names(simData)))) {
    c("pixelIndex", cols)
  } else {
    cols
  }

  for (i in as.character(reps)) {
    for (j in as.character(years)) {
      obsDataSubset <- obsData[rep == i & year == j, .SD, .SDcols = cols2]
      simDataSubset <- simData[rep == i & year == j,  .SD, .SDcols = cols2]

      ## make sure order is the same for pixels and species
      if (identical(cols, cols2)) {
        obsDataSubset <- obsDataSubset[, ..cols]
      } else {
        obsDataSubset <- obsDataSubset[simDataSubset[, .(pixelIndex)], on = "pixelIndex"][, ..cols]
      }

      ## convert to matrix
      obsDataSubset <- as.matrix(obsDataSubset[, ..cols])
      simDataSubset <- as.matrix(simDataSubset[, ..cols])

      if (NROW(simDataSubset) != NROW(obsDataSubset)) {
        stop("Simulated and observed data differ in no. of rows")
      }

      if (varType == "biomass" | varType == "counts") {
        logLiks <- sapply(seq(NROW(obsDataSubset)),
                          FUN = function(x) {
                            dmultinom(prob = simDataSubset[x,], x = obsDataSubset[x,], log = TRUE)
                          })
      } else {
        stop("varType must be 'biomass' or 'counts'")
      }
      ## before summing make sure infinites get the lowest value
      logLiks[abs(logLiks) == Inf] <- min(logLiks[abs(logLiks) != Inf])
      NegSumLogLik[i, j] <- -sum(logLiks)
    }
  }

  return(NegSumLogLik)
}

