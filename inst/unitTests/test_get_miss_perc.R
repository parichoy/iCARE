test_get_miss_perc <- function() {

  cutpoints <- 1:10
  ref_LP    <- 0:100
  miss_LP   <- 5.5
  probs     <- c(0, 0.25, 0.5, 0.75, 1)
  miss_perc <- 5

  ret       <- iCARE:::get_miss_perc(cutpoints, miss_LP) 

  checkEqualsNumeric(ret, 5, tolerance=1.0e-8)
  
}
