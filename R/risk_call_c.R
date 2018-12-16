call_c1 <- function(final, lps, ref_full_LP, miss, ref_risks, refmat, betavec, zmat, pop.weights, ncuts = 100, debug = 0) {


  DMISS            <- 999999999.9
  DMISS_TEST       <- 999999999.0
  zmat             <- zmat[, miss, drop = FALSE]
  temp             <- !is.finite(zmat)
  zmat[temp]       <- DMISS
  temp             <- !is.finite(ref_risks)
  ref_risks[temp]  <- DMISS
  temp             <- !is.finite(refmat)
  refmat[temp]     <- DMISS
  n_beta           <- length(betavec)
  probs            <- seq(0, 1, 1/ncuts)
  n_probs          <- length(probs)
  nr_z             <- ncol(zmat)
  nc_z             <- nrow(zmat)
  nr_ref           <- ncol(refmat)
  nc_ref           <- nrow(refmat)
  retvec           <- rep(DMISS, length(miss))
  retlps           <- rep(DMISS, length(miss))
  retflag          <- 1
  dim(ref_full_LP) <- NULL

  temp <- .C("ref_risk1", as.numeric(ref_risks), as.numeric(betavec), as.numeric(t(zmat)),
       as.numeric(t(refmat)), as.integer(n_beta), as.integer(nr_z), as.integer(nc_z),
       as.integer(nr_ref), as.integer(nc_ref), as.numeric(probs), as.integer(n_probs),
       as.integer(debug), as.numeric(pop.weights), as.double(ref_full_LP),
       retvec = as.numeric(retvec), retflag = as.integer(retflag), retlps=as.numeric(retlps))

  retflag <- temp$retflag
  if (retflag) stop("ERROR in c code")
  retvec  <- temp$retvec
  retlps  <- temp$retlps
  temp    <- retvec > DMISS_TEST
  if (any(temp)) retvec[temp] <- NA
  temp    <- retlps > DMISS_TEST
  if (any(temp)) retlps[temp] <- NA

  final[miss] <- retvec
  lps[miss]   <- retlps

  list(final=final, lps=lps)

} # END: call_c1

call_c2 <- function(final, lps, ref_full_LP, miss, refmat, betavec, zmat, age_new, age_int, popSubFromLP,
                    lambda_0, compRates0, pop.weights, ncuts = 100, debug = 0) {

  DMISS            <- 999999999.9
  DMISS_TEST       <- 999999999.0
  zmat             <- zmat[, miss, drop = FALSE]
  temp             <- !is.finite(zmat)
  zmat[temp]       <- DMISS
  temp             <- !is.finite(refmat)
  refmat[temp]     <- DMISS
  n_beta           <- length(betavec)
  probs            <- seq(0, 1, 1/ncuts)
  n_probs          <- length(probs)
  nr_z             <- ncol(zmat)
  nc_z             <- nrow(zmat)
  nr_ref           <- ncol(refmat)
  nc_ref           <- nrow(refmat)
  retvec           <- rep(DMISS, length(miss))
  retlps           <- rep(DMISS, length(miss))
  retflag          <- 1
  dim(ref_full_LP) <- NULL

  age_new         <- age_new[miss]
  age_int         <- age_int[miss]

  # Get all values of age
  temp                         <- c(age_new, age_new+age_int)
  maxa                         <- max(temp)
  n_lambda                     <- maxa + 1
  lambda                       <- rep(0, n_lambda)
  lambda[lambda_0[, 1]+1]      <- lambda_0[, 2]
  compRates                    <- rep(0, n_lambda)
  compRates[compRates0[, 1]+1] <- compRates0[, 2]
  temp                         <- is.na(compRates)
  compRates[temp]              <- DMISS

  temp <- .C("ref_risk2", as.numeric(betavec), as.numeric(t(zmat)), as.numeric(t(refmat)),
            as.integer(n_beta), as.integer(nr_z), as.integer(nc_z), as.integer(nr_ref),
            as.integer(nc_ref), as.numeric(probs), as.integer(n_probs), as.integer(debug),
            as.integer(age_new), as.integer(age_int), as.integer(n_lambda), as.numeric(popSubFromLP),
            as.numeric(lambda), as.numeric(compRates), as.numeric(pop.weights), as.double(ref_full_LP),
              retvec = as.numeric(retvec), retflag = as.integer(retflag), retlps=as.numeric(retlps))

  retflag <- temp$retflag
  if (retflag) stop("ERROR in c code")
  retvec  <- temp$retvec
  retlps  <- temp$retlps
  temp    <- retvec > DMISS_TEST
  if (any(temp)) retvec[temp] <- NA
  temp    <- retlps > DMISS_TEST
  if (any(temp)) retlps[temp] <- NA

  final[miss] <- retvec
  lps[miss]   <- retlps

  list(final=final, lps=lps)

} # END: call_c2
