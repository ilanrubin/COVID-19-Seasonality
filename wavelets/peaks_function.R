
# function to find local maxima
findpeaks <- function(x, y, prom=1, w=1, ...) {
      n <- length(y)
      
      # find the maximum of a rolling window of w indices in each direction
      y.max <- zoo::rollapply(zoo(y), 2*w+1, max, 
                align="center")
      # find the difference between each of these rolling maximum and
      # every other time point excluding the beginning and end
      delta <- y.max - y[-c(1:w, n+1-1:w)]
      # find the indices of the times that are a local (in the rolling window) maximum
      i.max <- which(delta == 0) + w

      
  # Iteratively remove the least prominent peak until all remaining peaks
  # meet the prominence criterion.
  any.cut <- TRUE
  while (any.cut && length(i.max) > 0) {
    peak.lower.min <- rep(NA, length(i.max))
    peak.upper.min <- rep(NA, length(i.max))
    peak.height <- y[i.max]

    # Find the minimum value to the left of the first peak
    peak.lower.min[1] <- min(y[1:i.max[1]])

    # Find the minimum value to the right of the last peak
    peak.upper.min[length(i.max)] <- min(y[i.max[length(i.max)]:n])

    # Find minimum values between peaks
    if (length(i.max) > 1) {
      for (i in 2:length(i.max)) {
        peak.lower.min[i] <- min(y[i.max[i - 1]:i.max[i]])
        peak.upper.min[i - 1] <- min(y[i.max[i - 1]:i.max[i]])
      }
    }

    # Determine which peaks meet the prominence criterion
    peaks.are.prom <- (peak.height - peak.lower.min >= prom) &
                      (peak.height - peak.upper.min >= prom)

    # Check if any peaks need to be removed
    any.cut <- !all(peaks.are.prom,na.rm=TRUE)

    # Only cut *one* peak at a time (the least prominent)
    if (any.cut) {
      # Find the indices of peaks that *don't* meet the prominence
      # criterion, then select the *least* prominent among those
      non_prom_indices <- which(!peaks.are.prom)
      if(length(non_prom_indices) > 0){ #catch errors when the above vectors are empty
          index_to_remove <- non_prom_indices[order(peak.height[non_prom_indices])[1]]
          i.max <- i.max[-index_to_remove]
      } 
    }
  }
  return(list(x=x[i.max], i=i.max, y=y))
}
