# Function to compute cutpoints
get_cutpoints <- function(ref_LP, probs) {

  # Get the number of digits
  se     <- sd(ref_LP, na.rm = TRUE)
  if (se < 1e-12) se <- 1e-12
  prec   <- 1/(0.001*se)
  digits <- sum(prec >= 10^(1:16))
  cuts   <- quantile(ref_LP, probs = probs, names = FALSE)
  cuts   <- round_down(cuts, digits = digits) # new 9/1/15
  cuts   <- c(unique(cuts), Inf)# new 9/1/15
  if (length(cuts) < 2) stop("ERROR: computing cutpoints")

  cuts

} # END: get_cutpoints

round_down <-function(num, digits){
  
  floor( num*10^digits ) / (10^digits)
  
}

get_endpoints <- function(cutpoints, miss_LP) {

  ncuts     <- length(cutpoints)
  miss_perc <- sum(cutpoints <= as.numeric(miss_LP))
  if (miss_perc == ncuts) {
    miss_perc <- miss_perc - 1
  } else if (!miss_perc) {
    miss_perc <- 1
  }
  a <- cutpoints[miss_perc]
  b <- cutpoints[miss_perc+1]

  c(a, b)

} # END: get_endpoints

# Function to return missing percentile
get_miss_perc <- function(cutpoints, miss_LP) {

  ncuts     <- length(cutpoints)
  miss_perc <- sum(cutpoints <= as.numeric(miss_LP))
  if (miss_perc == ncuts) {
    miss_perc <- miss_perc - 1
  } else if (!miss_perc) {
    miss_perc <- 1
  }

  miss_perc

} # END: get_miss_perc

# Functon to return the endpoints
get_newEndpoints <- function(cutpoints, miss_LP, miss_perc, ref_LP) {

  # Remove cutpoint defined by miss_perc + 1
  ncuts <- length(cutpoints)
  left  <- miss_perc - 1
  right <- miss_perc + 2
  for (i in 1:ncuts) {
    if (left < 1) left <- 1
    if (right > ncuts) right <- ncuts
    cut1 <- cutpoints[left]
    cut2 <- cutpoints[right]
    temp <- (ref_LP >= cut1) & (ref_LP < cut2)
    temp[is.na(temp)] <- FALSE
    if (any(temp)) break
    left  <- left - 1
    right <- right + 1
  }

  c(cut1, cut2)

} # END: get_newEndpoints
