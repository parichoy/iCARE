# lambda should be matrix of two columns (age, incidence rates ) diff ages in the rows

pick_lambda <- function(t, lambda){

  a <- which(t==lambda[,1])
  lambda[a,2]
}

# c should be a two column matrix with first column times

# computes internal integral over u for estimate of A_j
get_int <- function(a,t,lambda, Z_new, beta_est, model.competing.incidence.rates, ZBETA = NULL){

  holder<-0
  if (is.null(ZBETA)) {
    ZBETA <- t(exp(t(Z_new)%*%beta_est))
  }

  for(u in min(a, na.rm = TRUE):max(t, na.rm = TRUE)){
    temp0  <- (u>=a)*(u<t)
    holder <- holder + temp0*( (pick_lambda(u,lambda))*ZBETA + model.competing.incidence.rates[which(u==model.competing.incidence.rates[,1]),2])
  }
  holder
}

get_beta_given_names <-function(model.cov.info, model.formula){

  check_triple_check(model.cov.info)
  model.cov.info = process_triple_check(model.cov.info)

  if(is.null(colnames(model.cov.info))){
    return_and_print("ERROR: model.cov.info must have same names and order as predictors in model.formula.")
    stop()
  }
  if(length(colnames(model.cov.info))!=length(all.vars(model.formula)[2:length(all.vars(model.formula))])){
    return_and_print("ERROR: model.cov.info must have same names and order as predictors in model.formula.")
    stop()
  }
  if( sum(colnames(model.cov.info)!=all.vars(model.formula)[2:length(all.vars(model.formula))])>0 ){
    return_and_print("ERROR: model.cov.info must have same names and order as predictors in model.formula.")
    stop()
  }
  variables <- unique(all.vars(model.formula))[-1]
  data_use  <- subset(model.cov.info, select = variables)
  data_use[,all.vars(model.formula)[1]] <- rep(0, nrow(data_use))
  predictors <- as.matrix(model.matrix(model.formula, data = as.data.frame(data_use))[,2:ncol(model.matrix(model.formula, data = as.data.frame(data_use)))])
  colnames(predictors)

}

###  make design matrix for apply.cov.profile
make_design_matrix_covs <- function(data_use, apply.cov.profile, model.formula){

  options(na.action='na.pass')
  init_n   <- nrow(data_use)
  jump_cov <- nrow(apply.cov.profile )

  data_use   <- rbind(data_use, apply.cov.profile )
  data_use[,all.vars(model.formula)[1]] <- rep(0, nrow(data_use))

  predictors <- as.matrix(model.matrix(model.formula, data = as.data.frame(data_use))[,2:ncol(model.matrix(model.formula, data = as.data.frame(data_use)))])
  set        <- rep(FALSE, nrow(predictors))
  set[(init_n+1):(init_n+jump_cov)]=TRUE
  apply.cov.profile     <- subset( predictors, subset = set)
  apply.cov.profile
}

###  make design matrix for model.ref.dataset
make_design_matrix_dist <- function(data_use, model.ref.dataset, model.formula){

  init_n   <- nrow(data_use)
  jump_cov <- nrow(model.ref.dataset)

  data_use <- rbind(data_use, model.ref.dataset)
  data_use <- cbind(data_use, rep(0, nrow(data_use)))
  colnames(data_use) <- c(colnames(data_use)[1:(length(colnames(data_use))-1)],  all.vars(model.formula)[1])

  predictors   <- as.matrix(model.matrix(model.formula, data = as.data.frame(data_use))[,2:ncol(model.matrix(model.formula, data = as.data.frame(data_use)))])
  pop.dist.mat <- predictors[(init_n+1):(init_n+jump_cov),]
  cbind(pop.dist.mat)
}

###### Compute Absolute Risk
comp_Aj <- function(Z_new, apply.age.start, apply.age.interval.length, lambda, beta_est, model.competing.incidence.rates){

    ZBETA <- t(exp(t(beta_est)%*%Z_new))
    avec  <- apply.age.start + apply.age.interval.length
    tvec  <- min(apply.age.start, na.rm = TRUE):max(avec, na.rm = TRUE)
    temp  <- match(tvec, lambda[,1])
    c2    <- t(Z_new)%*%beta_est
    c3    <- log(lambda[temp, 2])

    Aj_est <- rep(0, ncol(Z_new))
    for(i in 1:length(tvec)){
      t         <- tvec[i]
      temp      <- get_int(apply.age.start,t,lambda, Z_new, beta_est, model.competing.incidence.rates, ZBETA = ZBETA)
      this_year <- (t >= apply.age.start)*(t < avec)*exp(c3[i] + c2 - temp)

      Aj_est    <- Aj_est + this_year
    }
    Aj_est
}

### iterate to obtain lambda_0
precise_lambda0 <- function(lambda, approx_expectation_rr, beta_est, pop.dist.mat, pop.weights){
  diagnose                <- list()
  lambda_0                <- lambda
  precise_expectation_rr0 <- approx_expectation_rr-1
  precise_expectation_rr1 <- approx_expectation_rr
  iter                    <- 0
  while( sum(abs(precise_expectation_rr1 - precise_expectation_rr0 )) >0.001 ){

    iter = iter+1
    precise_expectation_rr0 = precise_expectation_rr1
    diagnose[[iter]]=precise_expectation_rr0
    # new expectation rr implies lambda0
    lambda_0[,2] = lambda[,2] / precise_expectation_rr0
    # that lambda0 implies new expectation rr
    # this should be Nobs x Ntimes
    this = survival_givenX(lambda_0, beta_est, pop.dist.mat)*matrix(nrow=nrow(pop.dist.mat), ncol=nrow(lambda_0), pop.weights, byrow=TRUE)
    denom = 1/colSums(this)
    probX_givenT = sweep( this, MARGIN = 2, denom, FUN="*" )

    # to compute next iteration
    precise_expectation_rr1 = colSums( sweep( probX_givenT, MARGIN = 1,  exp( pop.dist.mat %*% beta_est ), FUN="*") )
  }
  lambda_0[,2] = lambda[,2] / precise_expectation_rr1
  res = list();  res[[1]] = lambda_0;  res[[2]] = precise_expectation_rr1; res[[3]]=diagnose
  res
}

### getting lambda_0 iteratively
survival_givenX <- function(lambda_0, beta_est, pop.dist.mat){ # must produce matrix with Nobs x Ntimes

  ## cumulative up to but not including current time
  mult   <- -exp( pop.dist.mat %*% beta_est )
  cumLam <- cumsum( c(0,lambda_0[,2]) )[1:length(lambda_0[,2])]
  exp( mult %*% cumLam  )
}

## transform model.cov.info into a form easily usable in the function
process_triple_check <- function(model.cov.info){
  maximum_dim = 0
  for(i in 1:length(model.cov.info)){
    maximum_dim = max(maximum_dim, length(model.cov.info[[i]]$levels) )
  }
  matt = data.frame(starter_column = rep(0, maximum_dim))

  for(i in 1:length(model.cov.info)){

    if(model.cov.info[[i]]$type == "continuous"){
        matt[,model.cov.info[[i]]$name] = rep(0,maximum_dim)
    }
    if(model.cov.info[[i]]$type == "factor"){
        matt[,model.cov.info[[i]]$name] = factor( c(model.cov.info[[i]]$levels, rep(model.cov.info[[i]]$levels[1], maximum_dim - length(model.cov.info[[i]]$levels))), levels = unique(c(model.cov.info[[i]]$levels, rep(model.cov.info[[i]]$levels[1], maximum_dim - length(model.cov.info[[i]]$levels)))))
        ### if there is a referent level specified, recode that way
        if(is.null(model.cov.info[[i]]$ref)==FALSE){
            matt[,model.cov.info[[i]]$name] = relevel(matt[,model.cov.info[[i]]$name], ref = as.character(model.cov.info[[i]]$ref))
        }
    }
  }
  matt[,1] <- NULL
  matt
}

## helper function
## creates model.cov.info from a dataframe
## option to flag factor variables by including a vector of Variable names
## if factor_vars input is not provided, then considered continuous if more than 12 unique levels, else factor
make_triple_check <- function(frame, factor_vars = NULL){

  if( !is.data.frame(frame) ){
    return_and_print("ERROR: input must be a dataframe")
    stop()
  }

  the_names = colnames(frame)
  model.cov.info = list()
  length(model.cov.info)=length(the_names)

  if(is.null(factor_vars)==FALSE){
    if(is.vector(factor_vars)!=TRUE ){
      return_and_print(paste("ERROR: The optional 'flag' of factor variables must be specified as a vector.", sep="") )
      stop()
    }
    if( prod( is.element(factor_vars, the_names)) == 0 ){
      return_and_print(paste("ERROR: The optional 'flag' of factor variables must contain names that match the column names of the dataframe.", sep="") )
      stop()
    }
  }

  for(i in 1:ncol(frame)){

    temp = list()
    temp$name = the_names[i]
    ###
    if(is.null(factor_vars)==FALSE){
      toggle = !is.element(temp$name, factor_vars)
    }
    if(is.null(factor_vars)){
      # binary indicator of whether the variable is continuous
      toggle = length(c(unique(frame[,i])))>12
      if( is.factor(frame[,13]) ==TRUE ){
        toggle = FALSE
      }
    }###
    if(toggle){
      temp$type = "continuous"
    }
    if(toggle == FALSE){
      temp$type = "factor"
      temp$levels = c(unique(frame[,i]))
      temp$ref = temp$levels[1]
    }
    model.cov.info[[i]]=temp
  }
  model.cov.info
}

## Function for combining Absolute Risks for Adjacent Intervals
combine_risks_two_intervals <- function(age_start_1, age_interval_1, absolute_risks_1, age_start_2, age_interval_2,
                                        absolute_risks_2){

  if( sum(age_start_1 + age_interval_1 != age_start_2) >0 ){
    return_and_print(paste("ERROR: To combine risks for two intervals, they must be adjacent.", sep="") )
    stop()
  }

  overall_age_start = age_start_1
  overall_age_interval = age_interval_1 + age_interval_2
  overall_absolute_risks =  absolute_risks_1 + (1 - absolute_risks_1)*absolute_risks_2

  result = list()
  result$details = cbind(overall_age_start, overall_age_interval, overall_absolute_risks)
  result$risk = cbind(overall_absolute_risks)
  result
}

## helper function for computing absolute risk of competing mortality
comp_AR_compete <- function(apply.age.start, apply.age.interval.length, model.competing.incidence.rates){
  ###### Compute A_j

  Aj_est <- 0
  for(t in min(apply.age.start, na.rm = TRUE):max(apply.age.start+apply.age.interval.length, na.rm = TRUE)){

    this_year <- (t>=apply.age.start)*(t<apply.age.start+apply.age.interval.length)*exp( log(pick_lambda(t, model.competing.incidence.rates)) - get_int_compete(apply.age.start,t, model.competing.incidence.rates) )#ex[t-apply.age.start+1]=Aj_est
    Aj_est    <- Aj_est + this_year
  }
  Aj_est
}

## helper function for computing survival of competing risk
get_int_compete <- function(a,t,model.competing.incidence.rates){

  holder <- 0
  for(u in min(a, na.rm = TRUE):max(t, na.rm = TRUE)){
    holder <- holder+(u>=a)*(u<t)*model.competing.incidence.rates[which(u==model.competing.incidence.rates[,1]),2]
  }
  holder
}

sim_snps <- function(betas, freqs, fh_status){

  snps <- matrix( NA, ncol = length(betas), nrow = length(fh_status))

  prob012_fh_no = cbind( (1-freqs)^2, 2*freqs*(1-freqs), freqs^2)

  beta_mat = matrix( betas, nrow = length(betas), ncol = 3, byrow = FALSE)
  top = exp( beta_mat*matrix(c(0,1,2), nrow = length(betas), ncol = 3, byrow = TRUE)/2 )*matrix(prob012_fh_no, nrow = length(betas), ncol = 3, byrow = FALSE)
  bottom = matrix(rowSums(top), nrow =length(betas), ncol = 3, byrow = FALSE)

  prob012_fh_yes = top/bottom

  fh_no  <- which(fh_status==0)
  fh_yes <- which(fh_status==1)

  vals = matrix( runif(length(betas)*length(fh_status)), ncol = length(betas), nrow = length(fh_status))

  if(length(fh_no)>0){
    snps[fh_no,] =  ( vals[fh_no,] > matrix(prob012_fh_no[,1], nrow = length(fh_no), ncol = length(prob012_fh_no[,1]), byrow = TRUE ) ) + ( vals[fh_no,] > matrix(rowSums(subset(prob012_fh_no, select=1:2)), nrow = length(fh_no), ncol = length(prob012_fh_no[,1]), byrow = TRUE ) )
  }
  if(length(fh_yes)>0){
    snps[fh_yes,] =  ( vals[fh_yes,] > matrix(prob012_fh_yes[,1], nrow = length(fh_yes), ncol = length(prob012_fh_yes[,1]), byrow = TRUE ) ) + ( vals[fh_yes,] > matrix(rowSums(subset(prob012_fh_yes, select=1:2)), nrow = length(fh_yes), ncol = length(prob012_fh_yes[,1]), byrow = TRUE ) )
  }
  snps

}

## helper function for handling missing data
handle_missing_data <- function(use.c.code, apply.age.start, apply.age.interval.length, Z_new, miss, present, ncuts, final_risks,
                                ref_pop, pop.weights, lambda_0, beta_est, model.competing.incidence.rates, lps){
  
  ref_full_LP <- t(ref_pop)%*%beta_est
  ###### Handle Missing Data  ##### If All times are the same
  if( length(unique(apply.age.start[miss]))==1 & length(unique(apply.age.interval.length[miss]))==1){
    
    pop.apply.age.start = unique(apply.age.start)  # change to single values so don't have to worry about dimension of ref_pop
    pop.apply.age.interval.length = unique(apply.age.interval.length)

    ###### Compute A_j_pop for ref_risks

    ref_risks   <- comp_Aj(ref_pop, pop.apply.age.start, pop.apply.age.interval.length, lambda_0, beta_est, model.competing.incidence.rates)

    if (use.c.code) {
      temp <- call_c1(final_risks, lps, ref_full_LP, miss, ref_risks, ref_pop, beta_est[, 1],
                             Z_new, pop.weights, ncuts=ncuts, debug = 0)
      final_risks <- temp$final ## will be: final_risks  <- temp[[1]]
      lps         <- temp$lps      ## will be: lps <- temp[[2]]
    } else {
      BETA  <- t(beta_est)
      Z_NEW <- Z_new
      REF   <- ref_pop
      PROBS <- seq(0, 1, 1/ncuts)

      for(i in 1:length(miss)){
        missi <- miss[i]

        ## make sure LPs based on non-missing covariates for the observation with missing
        present     <- which(is.na(Z_NEW[, missi])!=TRUE)
        if(length(present)==0){
          final_risks[missi] <- 0 
        }else{
        BETAV       <- BETA[1, present, drop = FALSE]
        ref_LP      <- BETAV%*%REF[present, , drop = FALSE] #just LP from predictors?
        miss_LP     <- BETAV%*%Z_NEW[present, missi, drop = FALSE] #just LP from predictors?

        cutpoints           <- get_cutpoints(ref_LP, PROBS)
        ncuts               <- length(cutpoints)
        miss_perc           <- get_miss_perc(cutpoints, miss_LP)
        these               <- (ref_LP >= cutpoints[miss_perc]) & (ref_LP < cutpoints[miss_perc+1])
        these[is.na(these)] <- FALSE
        if (!sum(these)) {
          temp                <- get_newEndpoints(cutpoints, miss_LP, miss_perc, ref_LP)
          these               <- (ref_LP >= temp[1]) & (ref_LP < temp[2])
          these[is.na(these)] <- FALSE
        }
        final_risks[missi]  <- weighted.mean(ref_risks[these], w = pop.weights[these], na.rm = TRUE)
        lps[missi]          <- weighted.mean(ref_full_LP[these], w = pop.weights[these], na.rm = TRUE)
        }
      }
    }
  } # END: if( length(unique(apply.age.start[miss]))==1 & length(unique(apply.age.interval.length[miss]))==1)


  ###### Handle Missing Data  ##### If All times are different
  if( length(unique(apply.age.start[miss])) > 1 || length(unique(apply.age.interval.length[miss])) > 1){
    if (use.c.code) {
      popSubFromLP = rep(0, ncol(ref_pop))
      temp  <- call_c2(final_risks, lps, ref_full_LP, miss, ref_pop, beta_est[, 1], Z_new, apply.age.start,
                               apply.age.interval.length, popSubFromLP, lambda_0, model.competing.incidence.rates, pop.weights,
                               ncuts = ncuts, debug = 0)
      final_risks <- temp$final ## will be: final_risks  <- temp[[1]]
      lps         <- temp$lps      ## will be: lps <- temp[[2]]
    } else {

      BETA  <- t(beta_est)
      Z_NEW <- Z_new
      REF   <- ref_pop
      CVEC  <- 2:(ncuts-1)
      PROBS <- seq(0, 1, 1/ncuts)

      for(i in 1:length(miss)){
        missi <- miss[i]

        pop.apply.age.start = apply.age.start[missi]  # change to single values so don't have to worry about dimension of ref_pop
        pop.apply.age.interval.length = apply.age.interval.length[missi];

        ## make sure LPs based on non-missing covariates for the observation with missing
        present = which(is.na(Z_new[,missi])!=TRUE)
        if(length(present)==0){
          final_risks[missi] <- 0 
        }else{
        BETAV       <- BETA[1, present, drop = FALSE]
        ref_LP      <- BETAV%*%REF[present, , drop = FALSE] #just LP from predictors?
        miss_LP     <- BETAV%*%Z_NEW[present, missi, drop = FALSE] #just LP from predictors?

        cutpoints           <- get_cutpoints(ref_LP, PROBS)
        ncuts               <- length(cutpoints)
        miss_perc           <- get_miss_perc(cutpoints, miss_LP)
        these               <- (ref_LP >= cutpoints[miss_perc]) & (ref_LP < cutpoints[miss_perc+1])
        these[is.na(these)] <- FALSE
        if (!sum(these)) {
          temp                <- get_newEndpoints(cutpoints, miss_LP, miss_perc, ref_LP)
          these               <- (ref_LP >= temp[1]) & (ref_LP < temp[2])
          these[is.na(these)] <- FALSE
        }
        temp                <- REF[,these, drop = FALSE]
        ref_risks           <- comp_Aj(temp, pop.apply.age.start, pop.apply.age.interval.length,
                                       lambda_0, beta_est, model.competing.incidence.rates)
        final_risks[missi]  <- weighted.mean(ref_risks, w = pop.weights[these], na.rm = TRUE)
        lps[missi]          <- weighted.mean(ref_full_LP[these], w = pop.weights[these], na.rm = TRUE)
      }}
    }
  }
  ## end function
  res = list(); res[[1]] = final_risks; res[[2]] = lps; res
}

decide_if_SNP_only <- function(apply.cov.profile , model.formula, model.log.RR, model.ref.dataset, model.cov.info,
                               model.snp.info, apply.snp.profile, apply.age.start, apply.age.interval.length){

  if(is.null(apply.cov.profile ) & is.null(model.formula) & is.null(model.log.RR) & is.null(model.ref.dataset) & is.null(model.cov.info)){
    covs_in_model = 0
    if( is.null(model.snp.info) ){
      return_and_print("ERROR: You are fitting a SNP-only Model, and thus must provide the 'model.snp.info' object.")
      stop()
    }
    if(is.null(apply.snp.profile)){
      if(length(apply.age.start)==1 & length(apply.age.interval.length)==1){
        apply.snp.profile = matrix(ncol = nrow(model.snp.info), nrow = 10000, NA)
        cat("\nNote: You did not provide apply.snp.profile.  Will impute SNPs for 10000 people.\n")
        cat("If require more, please provide apply.snp.profile input.\n")
      }else{
        apply.snp.profile = matrix(ncol = nrow(model.snp.info), nrow = length(apply.age.start), NA)
        print(paste("Note: You did not provide apply.snp.profile.  Imputed SNPs for ", length(apply.age.start) , " individuals, matching number of age intervals specified.", sep=""))
      }
      }
    if(is.vector(apply.snp.profile)){
      apply.snp.profile = cbind(apply.snp.profile)
     }
  }else{

    if( is.null(apply.cov.profile ) + is.null(model.formula) + is.null(model.log.RR) + is.null(model.ref.dataset) + is.null(model.cov.info) >0  ){
      return_and_print("ERROR: If any of apply.cov.profile , model.formula, model.log.RR, model.ref.dataset, or model.cov.info are NULL then they must all be NULL, and will define a SNP-only Model.")
      stop()
    }
    if(is.vector(apply.snp.profile)){
      apply.snp.profile = cbind(apply.snp.profile)
    }
    covs_in_model = 1
  }
  res <- list(); res[[1]] = covs_in_model; res[[2]] = apply.snp.profile; res
}

process_SNP_info <- function(covs_in_model, apply.snp.profile, model.bin.fh.name, apply.cov.profile , model.ref.dataset, model.snp.info){
  if( covs_in_model){

    if(is.null(apply.snp.profile)){
      apply.snp.profile = matrix(ncol = nrow(model.snp.info), nrow = nrow(apply.cov.profile ), NA)
      print("Note: You included snp_info, but did not provide apply.snp.profile.  Will impute all SNPs. ")
    }
    if(nrow(apply.snp.profile)!=nrow(apply.cov.profile )){
      return_and_print("ERROR: apply.cov.profile  and apply.snp.profile must have same number of rows. ")
      stop()
    }
    if(!is.na(model.bin.fh.name)){
      if( !is.element( model.bin.fh.name, colnames(apply.cov.profile )) ){
        return_and_print("ERROR: model.bin.fh.name must contain the variable name of family history (matching a column name in apply.cov.profile ) if it is in the model, otherwise NA.")
        stop()
      }else{
        fh.pop = model.ref.dataset[,model.bin.fh.name]
        fh.cov = apply.cov.profile [,model.bin.fh.name]
        attenuate.fh = 1
        if( prod(is.element(fh.pop, c(0,1,"0", "1", NA)) ) == 0){
          return_and_print("ERROR: The family history must be binary when using snp_info functionality. Check input for model.ref.dataset.")
          stop()
        }
        if( prod(is.element(fh.cov, c(0,1,"0", "1", NA)) ) == 0){
          return_and_print("ERROR: The family history must be binary when using snp_info functionality. Check input for apply.cov.profile .")
          stop()
        }}
    }else{ # model.bin.fh.name = NA - family history not in the model
      fh.pop = rep(0, nrow( model.ref.dataset))
      fh.cov = rep(0, nrow( apply.cov.profile ))
      attenuate.fh = 0
      print("Note: As specified, the model does not adjust SNP imputations for family history, as model.bin.fh.name = NA.")
    }
  }else{
    fh.pop = rep(0, 10000)
    fh.cov = rep(0, nrow( apply.snp.profile))
    attenuate.fh = 0
    print("Note: As specified, the model does not adjust SNP imputations for family history.")
  }
  res <- list(); res[[1]] = attenuate.fh; res[[2]] = fh.pop; res[[3]] = fh.cov; res[[4]] = apply.snp.profile; res
}

#' Building and Applying an Absolute Risk Model
#'
#' This function is used to build absolute risk models and apply them to estimate absolute risks.
#'
#' @param ref.pop design matrix created from dataframe of risk factors for a sample of subjects representative of underlying population, no missing values.
#' @param apply.age.start single integer or vector of integer ages for the start of the interval over which to compute absolute risk.
#' @param apply.age.interval.length single integer or vector of integer years over which absolute risk should be computed.
#' @param lambda_0 two column matrix [ integer ages, incidence rates] with incidence rate of disease. Must fully cover age interval for estimation.
#' @param beta_est vector with log odds ratios corresponding to the model params; no intercept; names must match design matrix, \code{ref_pop}.
#' @param model.competing.incidence.rates two column matrix [ integer ages, incidence rates] with incidence rate of competing events. Must fully cover age interval for estimation.
#' @param handle.snps binary indicator of whether SNPs are in the design matrix or not.  If so, design matrix has been stacked 5 times.
#' @return This function returns risks for the referent risk factor dataset.
#' @export
#' @examples
#' print("Hello world")
get_refs_risk <- function(ref_pop, apply.age.start, apply.age.interval.length, lambda_0, beta_est, model.competing.incidence.rates, handle.snps, n.imp ){

  refs.risk = comp_Aj(ref_pop, apply.age.start[1], apply.age.interval.length[1], lambda_0, beta_est, model.competing.incidence.rates)
  refs.lps = t(ref_pop)%*%beta_est
  
  if(handle.snps){
    refs.risk = rowMeans( matrix(refs.risk, ncol=n.imp, nrow=length(refs.risk)/n.imp, byrow=FALSE))
    refs.lps  = rowMeans( matrix(refs.lps,  ncol=n.imp, nrow=length(refs.lps)/n.imp, byrow=FALSE))
  }
  res = list(); res[[1]] = refs.risk; res[[2]] = refs.lps ; res
}


return_and_print <- function(x){
  print(x)
  return(x)  
}

