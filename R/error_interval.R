winkler_score <- function(lt, ut, actual, level) {
  alpha <- 1 - level / 100
  score <- ifelse(
    actual < lt,
    (ut - lt) + (2 / alpha) * (lt - actual),
    ifelse(
      actual > ut,
      (ut - lt) + (2 / alpha) * (actual - ut),
      ut - lt
    )
  )
  # Age-year-specific
  #return(score)
  # Age-specific
  #return(rowMeans(score, na.rm=TRUE))
  # Year-specific
  #return(colMeans(score, na.rm=TRUE))
  # Overall
  return(mean(score, na.rm = TRUE))
}

