test_get_newEndpoints <- function() {

  cutpoints <- 1:10
  ref_LP    <- 0:100
  miss_LP   <- 5.5
  probs     <- c(0, 0.25, 0.5, 0.75, 1)
  miss_perc <- 5

  ret       <- iCARE:::get_newEndpoints(cutpoints, miss_LP, miss_perc, ref_LP) 

  checkEquals(sum(ret == c(4, 7)), 2)
  
}
