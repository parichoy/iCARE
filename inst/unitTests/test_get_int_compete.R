test_get_int_compete <- function() {

  lambda <- matrix(1:20, nrow=10, byrow=FALSE)
  rates  <- lambda
  ret    <- iCARE:::get_int_compete(1,2,rates)

  checkEqualsNumeric(ret, 11, tolerance=1.0e-8)
  
}

