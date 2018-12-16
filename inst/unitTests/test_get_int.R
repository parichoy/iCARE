test_get_int <- function() {

  lambda <- matrix(1:20, nrow=10, byrow=FALSE)
  ZBETA  <- 0.5
  rates  <- lambda
  ret    <- iCARE:::get_int(1,2,lambda, NULL, NULL, rates, ZBETA=ZBETA)

  checkEqualsNumeric(ret, 16.5, tolerance=1.0e-8)
  
}

