test_get_endpoints <- function() {

  cutpoints <- 1:10
  ref_LP    <- 0:100
  miss_LP   <- 5.5
  probs     <- c(0, 0.25, 0.5, 0.75, 1)
  miss_perc <- 5

  ret       <- iCARE:::get_endpoints(cutpoints, miss_LP) 

  checkEquals(sum(ret == 5:6), 2)
  
}
