# History: May 29, 2015 Added option use.c.code
#          Jun 02, 2015 Updated cutpoints to first round them and then 
#                       remove duplicates
#          Jun 05, 2015 Added update_cutpoints function to update the 
#                       cutpoints in case
#                       whe there are no subjects between 2 cutpoints
#          Jun 08, 2015 Reorganized to be more modular, tidied up formatting
#          Jul 02, 2015 Modified for weights

computeAbsoluteRisk <- function(model.formula = NULL, model.cov.info = NULL, 
                                  model.snp.info = NULL, model.log.RR = NULL,
                                  model.ref.dataset = NULL, 
                                  model.ref.dataset.weights = NULL,
                                  model.disease.incidence.rates,
                                  model.competing.incidence.rates = NULL,
                                  model.bin.fh.name = NA, n.imp = 5,
                                  apply.age.start, apply.age.interval.length,
                                  apply.cov.profile  = NULL,
                                  apply.snp.profile = NULL,
                                  use.c.code = 1,  return.lp = FALSE, 
                                  return.refs.risk = FALSE){

  decision      <- decide_if_SNP_only(apply.cov.profile, model.formula, 
                     model.log.RR, model.ref.dataset, model.cov.info, 
                     model.snp.info, apply.snp.profile, apply.age.start, 
                     apply.age.interval.length)
  covs_in_model <-  decision[[1]]
  if(length(decision)>1) apply.snp.profile   <-  decision[[2]]

  
  handle.snps = 1 - is.null(model.snp.info)

  if (covs_in_model) {
    temp         <- check_age_lengths(apply.age.start, apply.age.interval.length, 
                          apply.cov.profile , "apply.cov.profile ")
    apply.age.start      <- temp[[1]]
    apply.age.interval.length <- temp[[2]]

    temp         <- check_model_inputs(apply.cov.profile , model.log.RR, 
                       model.ref.dataset, model.ref.dataset.weights, 
                       model.cov.info, model.formula, n.imp)
    data_use     <- temp[[1]]
    model.log.RR   <- temp[[2]]
    model.ref.dataset.weights <- temp[[3]]

    design_covs <- make_design_matrix_covs(data_use, apply.cov.profile, 
                      model.formula)

    check_design_matrix(model.log.RR, design_covs)
    Z_new        <- t(design_covs); rm(design_covs); gc()
    beta_est     <- as.matrix(nrow = length(model.log.RR), ncol=1, model.log.RR)
    pop.dist.mat <- make_design_matrix_dist(data_use, model.ref.dataset, 
                       model.formula)
    pop.weights  <- model.ref.dataset.weights
  }

  if(handle.snps){

    check_SNP_info(model.snp.info)

    snps.names <- paste(model.snp.info[,"snp.name"])
    snps.betas <- log( model.snp.info[,"snp.odds.ratio"] )
    snps.freqs <- model.snp.info[,"snp.freq"]

    processed_info = process_SNP_info(covs_in_model, apply.snp.profile, 
                      model.bin.fh.name, apply.cov.profile , 
                      model.ref.dataset, model.snp.info)

    attenuate.fh  <-  processed_info[[1]]
    fh.pop        <-  processed_info[[2]]
    fh.cov        <-  processed_info[[3]]
    apply.snp.profile  <-  processed_info[[4]]
  }

  if (covs_in_model) {

    if(handle.snps){
      covariate_stack5 <- do.call("rbind", replicate(n.imp, pop.dist.mat, 
                                   simplify = FALSE))
      simulated_snps   <- sim_snps(snps.betas, snps.freqs, 
                              cbind( rep( fh.pop, n.imp)) )
      pop.dist.mat     <- cbind( simulated_snps, covariate_stack5 )
      pop.weights      <- rep( model.ref.dataset.weights, n.imp )

      if(attenuate.fh){
        var_prs <- sum((snps.betas^2)*2*snps.freqs*(1-snps.freqs))
        alpha   <- var_prs/2
        if(model.cov.info[[which(model.bin.fh.name==colnames(apply.cov.profile ))]]$type=="factor"){
          model.bin.fh.name = paste("as.factor(", model.bin.fh.name, ")1", sep="")
        }
        beta_est[model.bin.fh.name,1] <- beta_est[model.bin.fh.name,1] - alpha
      }
      beta_est <- as.matrix(nrow = length(c( snps.betas, beta_est)), ncol = 1, c( snps.betas, beta_est) )
      Z_new    <- rbind(t(apply.snp.profile), Z_new)
    }
  } else {
    temp         <- check_age_lengths(apply.age.start, apply.age.interval.length, 
                       apply.snp.profile, "apply.snp.profile")
    apply.age.start      <- temp[[1]]
    apply.age.interval.length <- temp[[2]]

    if(handle.snps){
      pop.dist.mat <- sim_snps(snps.betas, snps.freqs, cbind( rep( fh.pop, n.imp)) )
      pop.weights  <- rep( 1/nrow(pop.dist.mat), nrow(pop.dist.mat) )
      beta_est     <- as.matrix(nrow = length(c( snps.betas)), ncol = 1, c( snps.betas) )
      Z_new        <- t(apply.snp.profile)
    }
  }
  lambda                    <- model.disease.incidence.rates; 
  rm(model.disease.incidence.rates)
  res                       <- check_rates(model.competing.incidence.rates, lambda,  
                                         apply.age.start, apply.age.interval.length)
  lambda                    <- res[[1]]
  model.competing.incidence.rates <- res[[2]]; rm(res)
  approx_expectation_rr     <- weighted.mean(  exp( pop.dist.mat%*% beta_est), w = pop.weights, na.rm=TRUE)
  calculation               <- precise_lambda0(lambda, approx_expectation_rr, beta_est, pop.dist.mat, pop.weights)
  lambda_0                  <- calculation[[1]]
  precise_expectation_rr    <- calculation[[2]]
  pop.dist.mat              <- t(pop.dist.mat )

  ###### Compute A_j  ## only those without NA
  final_risks <- rep(NA, ncol(Z_new))
  lps <- t(Z_new)%*%beta_est
  these = which(colSums( is.na(Z_new))==0 )

  if( length(these) > 0 ){
      temp               <- subset(Z_new, select=these)
      final_risks[these] <- comp_Aj(temp, apply.age.start[these], 
             apply.age.interval.length[these], lambda_0, beta_est, 
             model.competing.incidence.rates)
  }
  miss    <- which(is.na(final_risks))
  present <- which(!is.na(final_risks))
  ref_pop <- pop.dist.mat;
  ncuts   <- 100
  t00     <- proc.time()
  
  res         <- handle_missing_data(use.c.code, apply.age.start, 
                   apply.age.interval.length, Z_new, miss, present, ncuts, final_risks, 
                   ref_pop, pop.weights, lambda_0, beta_est, 
                   model.competing.incidence.rates, lps)
  final_risks <- res[[1]]
  lps         <- res[[2]]
  
  these = which(colSums(!is.na(Z_new))==0) # nothing to match on w referent dataset
  if(length(these)>0){
    ref_risks = get_refs_risk(ref_pop, apply.age.start, apply.age.interval.length, 
          lambda_0, beta_est, model.competing.incidence.rates, handle.snps, n.imp )[[1]]
    final_risks[these] <- weighted.mean(ref_risks, w = pop.weights[1:length(ref_risks)], 
                           na.rm = TRUE) 
  }
  t11    <- proc.time() - t00;  print(t11)
  result <- package_results(final_risks, Z_new, covs_in_model, handle.snps, apply.age.start, 
              apply.age.interval.length, apply.cov.profile , model.log.RR, beta_est, 
              apply.snp.profile, snps.names, return.lp, lps)
  if(return.refs.risk) result$refs.risk <- get_refs_risk(ref_pop, apply.age.start, 
                         apply.age.interval.length, lambda_0, beta_est, 
                      model.competing.incidence.rates, handle.snps, n.imp )[[1]]
  class(result) <- "icare"

  result
}

#' Building and Applying an Absolute Risk Model:  Compute Risk over Interval Split in Two Parts
#'
#' This function is used to build an absolute risk model that incorporates different input parameters before and after a given time point. The model is then applied to estimate absolute risks.
#'
#' @param apply.age.start single integer or vector of integer ages for the start of the interval over which to compute absolute risk.
#' @param apply.age.interval.length single integer or vector of integer years over which absolute risk should be computed.
#' @param apply.cov.profile dataframe containing the covariate profiles for which absolute risk will be computed. Covariates must be in same order with same names as in \code{model.formula}.
#' @param model.formula an object of class \code{formula}: a symbolic description of the model to be fitted, e.g. Y~Parity+FamilyHistory.
#' @param model.disease.incidence.rates two column matrix [ integer ages, incidence rates] or three column matrix [start age, end age, rate] with incidence rate of disease. Must fully cover age interval for estimation.
#' @param model.log.RR vector with log odds ratios corresponding to the model params; no intercept; names must match design matrix arising from \code{model.formula} and \code{model.cov.info}; check names using function \code{check_design_matrix()}.
#' @param model.ref.dataset dataframe of risk factors for a sample of subjects representative of underlying population, no missing values. Variables must be in same order with same names as in \code{model.formula}.
#' @param model.ref.dataset.weights optional vector of sampling weights for \code{model.ref.dataset}.
#' @param model.cov.info contains information about the risk factors in the model ;  a main list containing a list for each covariate, which must have the fields:\cr
#'        \itemize{ \item \code{"name"} : a string with the covariate name, matching name in model.formula \cr
#'        \item \code{"type"} : a string that is either "continuous" or "factor".}
#'        If factor variable, then:\cr
#'        \itemize{\item \code{"levels"} : vector with strings of level names  \cr
#'        \item \code{"ref"} : optional field, string with name of referent level
#'        }
#' @param use.c.code binary indicator of whether to run the c program for fast computation.
#' @param model.competing.incidence.rates two column matrix [ integer ages, incidence rates] or three column matrix [start age, end age, rate] with incidence rate of competing events. Must fully cover age interval for estimation.
#' @param return.lp binary indicator of whether to return the linear predictor for each subject in apply.cov.profile.
#' @param apply.snp.profile data frame with observed SNP data (coded 0,1, 2, or NA). May have missing values.
#' @param model.snp.info dataframe with three columns [ rs number, odds ratio, allele frequency ]
#' @param model.bin.fh.name string name of family history variable, if in model. This must refer to a variable that only takes values 0,1, NA.
#' @param apply.snp.profile data frame with observed SNP data (coded 0,1, 2, or NA). May have missing values.
#' @param cut.time integer age for which to split computation into before and after
#' @param apply.cov.profile.2 see \code{apply.cov.profile}, to be used for estimation in ages after the cutpoint
#' @param model.formula.2 see \code{model.formula}, to be used for estimation in ages after the cutpoint
#' @param model.log.RR.2 see \code{model.log.RR}, to be used for estimation in ages after the cutpoint
#' @param model.ref.dataset.2 see \code{model.ref.dataset}, to be used for estimation in ages after the cutpoint
#' @param model.ref.dataset.weights.2 see \code{model.ref.dataset.weights}, to be used for estimation in ages after the cutpoint
#' @param model.cov.info.2 see \code{model.cov.info}, to be used for estimation in ages after the cutpoint
#' @param model.bin.fh.name.2 see \code{model.bin.fh.name}, to be used for estimation in ages after the cutpoint
#' @param n.imp integer value for number of imputations for handling missing SNPs.
#' @details Individualized Coherent Absolute Risk Estimators (iCARE) is a tool that allows researchers to quickly build models for absolute risk and apply them to estimate individuals' risk based on a set of user defined input parameters. The software gives users the flexibility to change or update models rapidly based on new risk factors or tailor models to different populations based on the specification of simply three input arguments: \itemize{ \item (1) a model for relative risk assumed to be externally derived \item (2) an age-specific disease incidence rate and \item (3) the distribution of risk factors for the population of interest.} The tool can handle missing information on risk factors for risk estimation using an approach where all estimates are derived from a single model through appropriate model averaging.
#' @return This function returns a list of results objects, including: \cr
#' \itemize{ \item \code{risk} : absolute risk estimates over the specified interval for subjects given by \code{apply.cov.profile} \cr
#' \item \code{details}: dataframe with the start of the interval, the end of the interval, the covariate profile, and the risk estimates for each individual \cr
#' \item \code{beta.used} : the log odds ratios used in the model \cr
#' \item \code{lps} : linear predictors for subjects in \code{model.cov.profile}, if requested by \code{return.lp} \cr
#' \item \code{refs.risk} : absolute risk estimates for subjects in \code{model.ref.dataset}, if requested by \code{return.refs.risk} }
#' @export
#' @examples
#' print("Hello world")
#' 
#' 
computeAbsoluteRiskSplitInterval <- function(apply.age.start, apply.age.interval.length, apply.cov.profile, model.formula, model.disease.incidence.rates,
                    model.log.RR, model.ref.dataset, model.ref.dataset.weights=NULL, model.cov.info, use.c.code=1, model.competing.incidence.rates = NULL,
                    return.lp = FALSE, apply.snp.profile = NULL, model.snp.info = NULL, model.bin.fh.name = NULL, cut.time = NULL, apply.cov.profile.2 = NULL,
                    model.formula.2 = NULL, model.log.RR.2 = NULL, model.ref.dataset.2 = NULL, model.ref.dataset.weights.2 = NULL, model.cov.info.2 = NULL, model.bin.fh.name.2 = NULL, n.imp = 5, return.refs.risk=FALSE){

  ## if any of the second arguments are specified, then two two version
  if( is.null( cut.time )*is.null( apply.cov.profile.2 )*is.null( model.formula.2 )*is.null( model.log.RR.2 )*is.null( model.ref.dataset.2 )*is.null( model.cov.info.2 ) ==0 ){

    if(is.null( cut.time )){
      return_and_print(paste("ERROR: If you wish to use different model inputs over parts of the age interval, must specify cut.time.", sep="") )
      stop()
    }
    if( sum( (cut.time < apply.age.start) + (cut.time > apply.age.start + apply.age.interval.length) ) >0 ){
      cat("\nNote: You provided cut.times outside the interval defined by apply.age.start and apply.age.interval.length.\n")
      cat("Will compute risk using either inputs1 or inputs2 only, depending on whether the interval is below or above the cutpoint.\n")
    }
    ## for the rest, if not specified use originals
    if(is.null( apply.cov.profile.2 )){
      apply.cov.profile.2 = apply.cov.profile
    }
    if(is.null( model.formula.2 )){
      model.formula.2 = model.formula
    }

    if(is.null( model.log.RR.2 )){
      model.log.RR.2 = model.log.RR
    }
    if(is.null( model.ref.dataset.2 )){
      model.ref.dataset.2 = model.ref.dataset
    }
    if(is.null( model.cov.info.2 )){
      model.cov.info.2 = model.cov.info
    }
    if(is.null( model.bin.fh.name.2 )){
      model.bin.fh.name.2 = model.bin.fh.name
    }

    apply.age.start.1 = apply.age.start
    apply.age.interval.length.1 = (cut.time - apply.age.start)*( (cut.time - apply.age.start) < apply.age.interval.length) +   (apply.age.interval.length)*( 1 - ( (cut.time - apply.age.start) < apply.age.interval.length))
    apply.age.interval.length.1[which(apply.age.interval.length.1<0)]=0

    apply.age.start.2 = apply.age.start.1 + apply.age.interval.length.1
    apply.age.interval.length.2 = apply.age.start + apply.age.interval.length - apply.age.start.2
    
    result1 = computeAbsoluteRisk(apply.age.start=apply.age.start.1, apply.age.interval.length=apply.age.interval.length.1, apply.cov.profile=apply.cov.profile  , model.formula=model.formula,   model.disease.incidence.rates=model.disease.incidence.rates, model.log.RR=model.log.RR,    model.ref.dataset=model.ref.dataset,   model.ref.dataset.weights=model.ref.dataset.weights,  model.cov.info=model.cov.info,   use.c.code=use.c.code, model.competing.incidence.rates=model.competing.incidence.rates, return.lp=return.lp, apply.snp.profile = apply.snp.profile, model.snp.info = model.snp.info,  model.bin.fh.name = model.bin.fh.name,  n.imp = n.imp, return.refs.risk=return.refs.risk )
    result2 = computeAbsoluteRisk(apply.age.start=apply.age.start.2, apply.age.interval.length=apply.age.interval.length.2, apply.cov.profile=apply.cov.profile.2, model.formula=model.formula.2, model.disease.incidence.rates=model.disease.incidence.rates, model.log.RR=model.log.RR.2, model.ref.dataset=model.ref.dataset.2, model.ref.dataset.weights=model.ref.dataset.weights.2, model.cov.info=model.cov.info.2, use.c.code=use.c.code, model.competing.incidence.rates=model.competing.incidence.rates, return.lp=return.lp, apply.snp.profile = apply.snp.profile, model.snp.info = model.snp.info,  model.bin.fh.name = model.bin.fh.name.2, n.imp = n.imp, return.refs.risk=return.refs.risk )

    result = combine_risks_two_intervals(result1$details[,"Int_Start"], result1$details[,"Int_End"] - result1$details[,"Int_Start"], result1$details[,"Risk_Estimate"], result2$details[,"Int_Start"], result2$details[,"Int_End"] - result2$details[,"Int_Start"], result2$details[,"Risk_Estimate"])
    result$details = cbind(apply.age.start.1, apply.age.start.2, apply.age.start.2 + apply.age.interval.length.2, apply.cov.profile , apply.cov.profile.2, result$risk)
    if(return.lp){
      result$lps.1 = result1$lps
      result$lps.2 = result2$lps
    }
    if(return.refs.risk){
      result$refs.risk.1 = result1$refs.risk
      result$refs.risk.2 = result2$refs.risk}
  }
  else{  # just run standard function
    result = computeAbsoluteRisk(apply.age.start=apply.age.start, apply.age.interval.length=apply.age.interval.length, apply.cov.profile=apply.cov.profile , model.formula=model.formula, model.disease.incidence.rates=model.disease.incidence.rates, model.log.RR=model.log.RR, model.ref.dataset=model.ref.dataset, model.ref.dataset.weights=model.ref.dataset.weights, model.cov.info=model.cov.info, use.c.code=use.c.code, model.competing.incidence.rates=model.competing.incidence.rates, return.lp=return.lp, apply.snp.profile = apply.snp.profile, model.snp.info = model.snp.info, model.bin.fh.name = model.bin.fh.name, n.imp = n.imp)
  }
  class(result) <- "icareSplit"

  result
}




