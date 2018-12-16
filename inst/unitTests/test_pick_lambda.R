test_pick_lambda <- function() {

  lambda <- matrix(1:20, nrow=10, byrow=FALSE)
  ret    <- iCARE:::pick_lambda(5, lambda)

  checkEqualsNumeric(ret, 15, tolerance=1.0e-8)
  
}

