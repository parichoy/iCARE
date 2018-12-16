test_call_c1 <- function() {

  final       <- rep(-9, 2)
  lps         <- rep(-9, 2)
  ref_LP      <- 0:100
  miss        <- 1
  ref_risks   <- 0:1000
  refmat      <- matrix(data=1, nrow=length(ref_risks), ncol=2)
  betavec     <- rep(1, 2)
  zmat        <- matrix(data=1, nrow=length(final), ncol=2)
  pop.weights <- rep(1, length(final)) 

  ret <- iCARE:::call_c1(final, lps, ref_LP, miss, ref_risks, refmat, betavec, zmat, pop.weights, ncuts=10) 
  x   <- ret$final[1]

  checkEqualsNumeric(x, 0.5, tolerance=1e-8)
  
}
