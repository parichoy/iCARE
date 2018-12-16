ModelValidation = function(study.data,
                            total.followup.validation = FALSE,
                            predicted.risk = NULL, 
                            predicted.risk.interval = NULL, 
                            linear.predictor = NULL, 
                            iCARE.model.object = 
                              list(model.formula = NULL,
                                   model.cov.info = NULL,
                                   model.snp.info = NULL,
                                   model.log.RR = NULL,
                                   model.ref.dataset = NULL,
                                   model.ref.dataset.weights = NULL,
                                   model.disease.incidence.rates = NULL,
                                   model.competing.incidence.rates = NULL,
                                   model.bin.fh.name = NA,
                                   apply.cov.profile  = NULL,
                                   apply.snp.profile = NULL, 
                                   n.imp = 5, use.c.code = 1,
                                   return.lp = TRUE, 
                                   return.refs.risk = TRUE),
                            number.of.percentiles = 10,
                            reference.entry.age = NULL, 
                            reference.exit.age = NULL,
                            predicted.risk.ref = NULL,
                            linear.predictor.ref = NULL,
                            linear.predictor.cutoffs = NULL,
                            dataset = "Example Dataset", 
                            model.name = "Example Risk Prediction Model"){
  
  model.formula = iCARE.model.object$model.formula
  model.cov.info = iCARE.model.object$model.cov.info
  model.snp.info = iCARE.model.object$model.snp.info
  model.log.RR = iCARE.model.object$model.log.RR
  model.ref.dataset = iCARE.model.object$model.ref.dataset
  model.ref.dataset.weights = 
    iCARE.model.object$model.ref.dataset.weights
  model.disease.incidence.rates = 
    iCARE.model.object$model.disease.incidence.rates
  model.competing.incidence.rates = 
    iCARE.model.object$model.competing.incidence.rates
  model.bin.fh.name = iCARE.model.object$model.bin.fh.name
  apply.cov.profile = iCARE.model.object$apply.cov.profile
  apply.snp.profile = iCARE.model.object$apply.snp.profile
  n.imp = iCARE.model.object$n.imp
  use.c.code = iCARE.model.object$use.c.code
  return.lp = iCARE.model.object$return.lp
  return.refs.risk = iCARE.model.object$return.refs.risk
  
  ##########Some coding adjustment to create output##############
  mf = paste(model.formula)
  if(total.followup.validation){
    timeframe = "Observed Followup"
  } else timeframe = paste(predicted.risk.interval,"years")
  ###############################################################
  
  observed.outcome = study.data$observed.outcome
  study.entry.age = study.data$study.entry.age
  study.exit.age = study.data$study.exit.age
  time.of.onset = study.data$time.of.onset
  sampling.weights = study.data$sampling.weights
  
  if(any(study.entry.age >= study.exit.age)){
    print("Error: entry age should be lower than exit age")
    stop()
  }
  
  #Compute the observed followup
  observed.followup = study.exit.age - study.entry.age
  
  # cases with time to event bigger than observed followup (data issue!) are censored
  observed.outcome[observed.outcome == 1 & time.of.onset > observed.followup] = 0
  time.of.onset[observed.outcome == 1 & time.of.onset > observed.followup] = Inf
  
  #the followup variable is subsequently adjusted based on the risk prediction
  #interval
  followup = observed.followup
  
  #Setting up the risk prediction interval if total followup validation is
  #desired
  if(total.followup.validation){
    predicted.risk.interval = observed.followup
  }
  
  #If a scalar is supplied for predicted.risk.interval (e.g., 5), create a
  #vector
  if(length(predicted.risk.interval) == 1){
    predicted.risk.interval =
      rep(predicted.risk.interval,length(observed.outcome))
  }
  
  #if time of onset falls within the risk prediction interval, for model
  #validation we adjust the followup to be equal to the minimum of risk prediction
  #interval and observed followup
  
  followup[(time.of.onset <= predicted.risk.interval & predicted.risk.interval
            <= observed.followup)] = 
    predicted.risk.interval[(time.of.onset <= predicted.risk.interval &
                               predicted.risk.interval <= observed.followup)]                              
  
  #if the time of onset is larger than the length of risk prediction interval
  #we make the outcome censored at the right end point of risk prediction
  #interval
  observed.outcome[predicted.risk.interval < time.of.onset & time.of.onset <=
                     observed.followup] = 0
  followup[predicted.risk.interval < time.of.onset & time.of.onset <=
             observed.followup] =
    predicted.risk.interval[predicted.risk.interval < time.of.onset &
                              time.of.onset <= observed.followup]
  
  #if the observed followup is larger than the length of the risk prediction
  #interval and the subject has not developed disease by end of followup, she
  #is considered censored at the right endpoint of the risk prediction
  #interval
  followup[predicted.risk.interval <= observed.followup & observed.followup <
             time.of.onset] = predicted.risk.interval[predicted.risk.interval
                                                      <= observed.followup &
                                                        observed.followup <
                                                        time.of.onset]
  #if either predicted risk or linear predictors are not supplied estimate
  #both using iCARE
  if(is.null(predicted.risk) | is.null(linear.predictor)){
    
    pred.risk = computeAbsoluteRisk(model.formula = model.formula,
                                    model.cov.info = model.cov.info,
                                    model.snp.info = model.snp.info,
                                    model.log.RR = model.log.RR,
                                    model.ref.dataset = model.ref.dataset,
                                    model.ref.dataset.weights = NULL,
                                    model.disease.incidence.rates =
                                      model.disease.incidence.rates,
                                    model.competing.incidence.rates =
                                      model.competing.incidence.rates,
                                    model.bin.fh.name = model.bin.fh.name,
                                    apply.age.start = study.entry.age,
                                    apply.age.interval.length = followup,
                                    apply.cov.profile 
                                    = apply.cov.profile,
                                    apply.snp.profile = apply.snp.profile,
                                    use.c.code = use.c.code, return.lp =
                                      return.lp, return.refs.risk =
                                      return.refs.risk)
    
    predicted.risk = as.numeric(pred.risk$risk)
    linear.predictor = as.numeric(pred.risk$lps)
    
  }
  
  if((is.null(predicted.risk.ref) | is.null(linear.predictor.ref)) & 
     (!is.null(reference.entry.age)) & (!is.null(reference.exit.age))){
    
    reference.followup = reference.exit.age - reference.entry.age
    
    predicted.risk.ref = computeAbsoluteRisk(model.formula = model.formula,
                                             model.cov.info = model.cov.info,
                                             model.snp.info = model.snp.info,
                                             model.log.RR = model.log.RR,
                                             model.ref.dataset = 
                                               model.ref.dataset,
                                             model.ref.dataset.weights = NULL,
                                             model.disease.incidence.rates =
                                               model.disease.incidence.rates,
                                             model.competing.incidence.rates =
                                               model.competing.incidence.rates,
                                             model.bin.fh.name = 
                                               model.bin.fh.name,
                                             apply.age.start = 
                                               reference.entry.age,
                                             apply.age.interval.length 
                                             = reference.followup,
                                             apply.cov.profile = 
                                               model.ref.dataset,
                                             apply.snp.profile = 
                                               NULL,
                                             use.c.code = use.c.code, 
                                             return.lp = return.lp,
                                             return.refs.risk =
                                               return.refs.risk)
    
    linear.predictor.ref = as.numeric(predicted.risk.ref$lps)
    predicted.risk.ref = as.numeric(predicted.risk.ref$refs.risk)
  
  }
  
  linear.predictor.cases = linear.predictor[observed.outcome == 1]
  linear.predictor.controls = linear.predictor[observed.outcome == 0]
  
  #matrix of indicators of risk score in cases bigger than risk score in 
  #controls
  indicator = vapply(linear.predictor.cases, function(x) x > 
                       linear.predictor.controls,
                     logical(length(linear.predictor.controls)))
  
  #In nested case-control study setting there will be a vector of
  #"sampling.weights"; in the full cohort study setting, "sampling.weights"
  #is set to NULL
  
  if(!is.null(sampling.weights)){
    
    freq = 1/sampling.weights
    
    #Computing Incidence Rates
    ages = (min(study.entry.age) + 1):(max(study.exit.age) - 1)
    len_ages = length(ages)
    age.of.onset = study.entry.age + time.of.onset
    
    study_incidence = rep(0,len_ages)
    
    for (i in 1:len_ages) {
      
      a=sum((study.entry.age <= ages[i] - 1 & age.of.onset >= ages[i] & 
               age.of.onset < ages[i] + 1 & study.exit.age >= ages[i]) * 
              freq)
      b=sum((study.entry.age <= ages[i] - 1 & study.exit.age >= ages[i] & 
               age.of.onset >= ages[i]) * freq)
      
      study_incidence[i] = a/b
      
    }
    
    study_incidence = cbind(ages,study_incidence)
    
    sampling.weights.cases = sampling.weights[observed.outcome == 1]
    sampling.weights.controls = sampling.weights[observed.outcome == 0]
    
    #computing the weight matrix for the IPW estimate of AUC
    freq.cases = freq[observed.outcome == 1]
    freq.controls = freq[observed.outcome == 0]
    
    weight.mat = matrix(kronecker(freq.controls,freq.cases), 
                        nrow = length(freq.controls), byrow = TRUE)
    auc = sum(indicator * weight.mat)/sum(weight.mat)
    
    #compute E_S_0[I(S_1 > S_0)]
    mean_S0_indicator = apply(indicator,2,weighted.mean,w = freq.controls)
    
    #compute E_S_1[I(S_1 > S_0)]
    mean_S1_indicator = apply(indicator,1,weighted.mean,w = freq.cases)
    
    #compute the variance of AUC
    term1 = wtd.var(mean_S0_indicator, weights = freq.cases) + 
      weighted.mean((mean_S0_indicator - auc)^2 * 
                      (1 - sampling.weights.cases)/sampling.weights.cases, 
                    w = freq.cases)
    
    term2 = wtd.var(mean_S1_indicator, weights = freq.controls) + 
      weighted.mean((mean_S1_indicator - auc)^2 * 
                      (1 - sampling.weights.controls)/sampling.weights.controls, 
                    w = freq.controls)
    
    var_auc = term1/sum(freq.cases) + term2/sum(freq.controls)
    
    #Overall expected by observed ratio
    p_E = (sum(predicted.risk * freq))/(sum(freq))
    p_O = (sum(observed.outcome * freq))/(sum(freq))
    
    #variance of observed risk
    var_obsrisk = (p_O * (1 - p_O) + sum((observed.outcome - p_O)^2 * 
                  (1 - sampling.weights)/sampling.weights^2)/sum(freq))/sum(freq)
    
    #variance of log expected by observed
    var_log_exp_by_obs = var_obsrisk/p_O^2
    
    exp_by_obs = p_E/p_O
    
    CI_exp_by_obs = exp(c(log(exp_by_obs) - 1.96 * sqrt(var_log_exp_by_obs),
                        log(exp_by_obs) + 1.96 * sqrt(var_log_exp_by_obs)))
    
    #create categories of the risk score
    
    if(!is.null(linear.predictor.cutoffs)){
      
      linear.predictor.cutoffs = c(min(linear.predictor),
                                   linear.predictor.cutoffs,
                                   max(linear.predictor))
      
      linear.predictor.cat = cut(linear.predictor, 
                                 breaks = linear.predictor.cutoffs, 
                                 include.lowest = TRUE)
    
    } else {
      linear.predictor.cat = weighted.quantcut(linear.predictor, weights =
                                               freq, q =
                                               number.of.percentiles)
    }
    
    #compute sample size at each category using weights
    wtd.samp.size.cat = aggregate(freq ~ linear.predictor.cat, FUN = 
                                    sum)$freq
    
    #Observed and predicted probabilities at the categories
    obs.prob.wt = observed.outcome * freq
    
    obs.prob.wt.cat = aggregate(obs.prob.wt ~ linear.predictor.cat, FUN =
                                  sum)
    wt.cat = aggregate(freq ~ linear.predictor.cat, FUN = sum)
    
    observed.prob.cat = obs.prob.wt.cat$obs.prob.wt/wt.cat$freq
    
    pred.prob.wt = predicted.risk * freq
    
    pred.prob.wt.cat = aggregate(pred.prob.wt ~ linear.predictor.cat, FUN =
                                   sum)
    
    predicted.prob.cat = pred.prob.wt.cat$pred.prob.wt/wt.cat$freq
    
    #create a variable that assigns each subject to his category specific 
    #observed 
    #risk the variable is "observed.risk.cat" and is of same length as the 
    #number of subjects
    
    observed.risk.cat = linear.predictor.cat  
    levels(observed.risk.cat) = observed.prob.cat
    observed.risk.cat = as.numeric(as.character(observed.risk.cat))
    
    #compute the variance correction term of absolute risk for each 
    #category
    var.absrisk.corr = (observed.outcome - observed.risk.cat)^2 * 
      (1 - sampling.weights)/(sampling.weights)^2
    
    var.absrisk.corr.cat = (aggregate(var.absrisk.corr ~ 
                                        linear.predictor.cat, FUN = 
                                        sum)$var.absrisk.corr)/wt.cat$freq
    
    #compute the variance of absolute risk in each category
    variance.absrisk.cat = (observed.prob.cat * (1 - observed.prob.cat) +
                              var.absrisk.corr.cat)/wt.cat$freq
    
    #compute the standard deviation of absolute risk in each cateogry
    sd.absrisk.cat = sqrt(variance.absrisk.cat)
    
    #compute the Hosmer-Lemeshow statistic for absolute risk
    chisq.absrisk = sum(((observed.prob.cat - predicted.prob.cat)^2)
                        /variance.absrisk.cat)
    
    #sigma2_bar = mean(variance.absrisk.cat)
    mu_bar = mean(observed.prob.cat)
    
    #Producing relative risk estimates
    
    observed.RR.cat = observed.prob.cat/mu_bar
    predicted.RR.cat = predicted.prob.cat/mean(predicted.prob.cat)
    
    #Compute the variance-covariance matrix for the absolute risks
    variance.mat.absrisk = diag(variance.absrisk.cat)
    
    #Compute the derivative matrix for log relative risk
    derivative_offdiagonal = -1/(number.of.percentiles * mu_bar)
    derivative_diagonal = diag(1/observed.prob.cat + derivative_offdiagonal)
    
    derivative_diagonal[lower.tri(derivative_diagonal)] = derivative_offdiagonal
    derivative_diagonal[upper.tri(derivative_diagonal)] = derivative_offdiagonal
    
    derivative.mat = derivative_diagonal[-number.of.percentiles,]
    
    varcov.logRR.cat = derivative_diagonal %*% variance.mat.absrisk %*% 
      derivative_diagonal
    
    sd.logRR.cat = sqrt(diag(varcov.logRR.cat))
    
    sigmainv_logRR = solve(derivative.mat %*% variance.mat.absrisk %*% 
                             t(derivative.mat))
    
    diff.logRR = (log(observed.RR.cat) - 
                    log(predicted.RR.cat))[-number.of.percentiles]
    
    chisq.logRR = as.numeric(diff.logRR %*% sigmainv_logRR %*% t(t(diff.logRR)))
    
    #p-value for Hosmer Lemeshow test
    PVAL_absrisk = 1 - pchisq(chisq.absrisk, number.of.percentiles) 
    
    #p-value for relative risk test
    PVAL_logRR = 1 - pchisq(chisq.logRR, number.of.percentiles - 1) 
    
    #Compute 95% confidence interval for observed absolute risk in each 
    #category
    upper.limit.absrisk.cat = observed.prob.cat + 1.96 * sd.absrisk.cat
    lower.limit.absrisk.cat = observed.prob.cat - 1.96 * sd.absrisk.cat
    
    #Compute 95% confidence interval for relative risk in each category
    upper.limit.RR.cat = exp(log(observed.RR.cat) + 1.96 * sd.logRR.cat)
    lower.limit.RR.cat = exp(log(observed.RR.cat) - 1.96 * sd.logRR.cat)
    
    #Compute 95% confidence interval for the AUC estimate
    upper.limit.auc = auc + 1.96 * sqrt(var_auc)
    lower.limit.auc = auc - 1.96 * sqrt(var_auc)
    
  } else {
    
    print("No sampling weights provided: assuming each subject has sampling weight 1")
    
    #Compute incidence rates
    
    ages = (min(study.entry.age) + 1):(max(study.exit.age) - 1)
    len_ages = length(ages)
    age.of.onset = study.entry.age + time.of.onset
    
    study_incidence = rep(0,len_ages)
    
    for (i in 1:len_ages) {
      
      a=sum(study.entry.age <= ages[i] - 1 & age.of.onset >= ages[i] 
            & age.of.onset < ages[i] + 1 & study.exit.age >= ages[i])
      b=sum(study.entry.age <= ages[i] - 1 & study.exit.age >= ages[i] 
            & age.of.onset >= ages[i])
      
      study_incidence[i] = a/b
      
    }
    
    study_incidence = cbind(ages,study_incidence)
    
    #estimate the AUC in full cohort setting
    auc = mean(indicator)
    
    #compute E_S_0[I(S_1 > S_0)]
    mean_S0_indicator = apply(indicator,2,mean)
    
    #compute E_S_1[I(S_1 > S_0)]
    mean_S1_indicator = apply(indicator,1,mean)
    
    #compute variance of AUC
    var_auc = var(mean_S0_indicator)/length(mean_S0_indicator) + 
      var(mean_S1_indicator)/length(mean_S1_indicator)
    
    #Overall expected by observed ratio
    p_E = mean(predicted.risk)
    p_O = mean(observed.outcome)
    
    #variance of observed risk
    var_obsrisk = (p_O * (1 - p_O))/length(predicted.risk)
    
    #variance of log expected by observed
    var_log_exp_by_obs = var_obsrisk/p_O^2
    
    exp_by_obs = p_E/p_O
    
    CI_exp_by_obs = exp(c(log(exp_by_obs) - 1.96 * sqrt(var_log_exp_by_obs),
                        log(exp_by_obs) + 1.96 * sqrt(var_log_exp_by_obs)))
    
    #create categories of the risk score
    if(!is.null(linear.predictor.cutoffs)){
      
      linear.predictor.cutoffs = c(min(linear.predictor),
                                   linear.predictor.cutoffs,
                                   max(linear.predictor))
      
      linear.predictor.cat = cut(linear.predictor, 
                                 breaks = linear.predictor.cutoffs, 
                                 include.lowest = TRUE)
    
    } else {
      
      linear.predictor.cat = quantcut(linear.predictor, 
                                      q = number.of.percentiles)
    }
    
    #compute sample size at each category
    samp.size.cat = as.numeric(table(linear.predictor.cat))
    
    #Observed and predicted probabilities at the categories
    observed.prob.cat = aggregate(observed.outcome ~ linear.predictor.cat,
                                  FUN = mean)$observed.outcome
    predicted.prob.cat = aggregate(predicted.risk ~ linear.predictor.cat,
                                   FUN = mean)$predicted.risk
    
    #compute the variance in each category
    variance.absrisk.cat = (observed.prob.cat * 
                              (1 - observed.prob.cat))/samp.size.cat
    
    #compute the standard deviation in each cateogry
    sd.absrisk.cat = sqrt(variance.absrisk.cat)
    
    #compute the Hosmer-Lemeshow statistic
    chisq.absrisk = sum(((observed.prob.cat - predicted.prob.cat)^2)
                        /variance.absrisk.cat)
    
    #sigma2_bar = mean(variance.absrisk.cat)
    mu_bar = mean(observed.prob.cat)
    
    #Producing relative risk estimates
    
    observed.RR.cat = observed.prob.cat/mu_bar
    predicted.RR.cat = predicted.prob.cat/mean(predicted.prob.cat)
    
    #Compute the variance-covariance matrix for the absolute risks
    variance.mat.absrisk = diag(variance.absrisk.cat)
    
    #Compute the derivative matrix for log relative risk
    derivative_offdiagonal = -1/(number.of.percentiles * mu_bar)
    derivative_diagonal = diag(1/observed.prob.cat + derivative_offdiagonal)
    
    derivative_diagonal[lower.tri(derivative_diagonal)] = derivative_offdiagonal
    derivative_diagonal[upper.tri(derivative_diagonal)] = derivative_offdiagonal
    
    derivative.mat = derivative_diagonal[-number.of.percentiles,]
    
    varcov.logRR.cat = derivative_diagonal %*% variance.mat.absrisk %*% 
      derivative_diagonal
    
    sd.logRR.cat = sqrt(diag(varcov.logRR.cat))
    
    sigmainv_logRR = solve(derivative.mat %*% variance.mat.absrisk %*% 
                             t(derivative.mat))
    
    diff.logRR = (log(observed.RR.cat) - 
                    log(predicted.RR.cat))[-number.of.percentiles]
    
    chisq.logRR = as.numeric(diff.logRR %*% sigmainv_logRR %*% t(t(diff.logRR)))
    
    #p-value for Hosmer Lemeshow test
    PVAL_absrisk = 1 - pchisq(chisq.absrisk, number.of.percentiles)  
    
    #p-value for relative risk test
    PVAL_logRR = 1 - pchisq(chisq.logRR, number.of.percentiles - 1) 
    
    #Compute 95% confidence interval for observed absolute risk in each 
    #category
    upper.limit.absrisk.cat = observed.prob.cat + 1.96 * sd.absrisk.cat
    lower.limit.absrisk.cat = observed.prob.cat - 1.96 * sd.absrisk.cat
    
    #Compute 95% confidence interval for relative risk in each category
    upper.limit.RR.cat = exp(log(observed.RR.cat) + 1.96 * sd.logRR.cat)
    lower.limit.RR.cat = exp(log(observed.RR.cat) - 1.96 * sd.logRR.cat)
    
    #Compute 95% confidence interval for the AUC estimate
    upper.limit.auc = auc + 1.96 * sqrt(var_auc)
    lower.limit.auc = auc - 1.96 * sqrt(var_auc)
    
  }
  
  #Output results
  DNAME <- paste(deparse(substitute(observed.frequency)),
                 deparse(substitute(expected.frequency)), 
                 sep = ", ")
  METHOD = "Hosmer and Lemeshow goodness of fit (GOF) test for Absolute Risk"
  names(chisq.absrisk) = "Chisquare"
  names(number.of.percentiles) = "df"
  
  h1 = structure(list(statistic = chisq.absrisk, parameter = 
                       number.of.percentiles, 
                     p.value = PVAL_absrisk, method = METHOD, 
                     data.name = DNAME),class = "htest")
  
  #cat("\n")
  
  DNAME <- paste(deparse(substitute(observed.frequency)),
                 deparse(substitute(expected.frequency)), 
                 sep = ", ")
  METHOD = "Goodness of fit (GOF) test for Relative Risk"
  names(chisq.logRR) = "Chisquare"
  names(number.of.percentiles) = "df"
  
  h2 = structure(list(statistic = chisq.logRR, parameter = 
                       number.of.percentiles - 1, 
                     p.value = PVAL_logRR, method = METHOD, data.name = DNAME), 
                 class = "htest")
  
  #cat("\n")
  
  #print(paste("Dataset: ",dataset))
  #print(paste("Model Name: ",model.name))
  #if(is.null(model.formula)){
  #  print("Model formula: Likely an additive SNP-only model")
  #} else {
  #print(paste("Model Formula:",mf[2],mf[1],mf[3]))
  #}
  #print(paste("Risk Prediction Interval:",timeframe))
  #cat("\n")
  #print(paste("Number of study subjects: "
  #                               ,length(observed.outcome)))
  #cat("\n")
  #print(paste("Number of cases: ",sum(observed.outcome)))
  #cat("\n")
  #print(paste("Follow-up time (years) [mean,range]: [",
  #                    round(mean(followup),3),", (",
  #                    round(range(followup)[1],3),",",
  #                    round(range(followup)[2],3),")"," ]"))
  #cat("\n")
  #print(paste("Baseline age (years) [mean,range]: [", 
  #            round(mean(study.entry.age),3),", (",
  #            round(range(study.entry.age)[1],3),",",
  #            round(range(study.entry.age)[2],3),")"," ]"))
  #cat("\n")
  #print("Absolute Risk Calibration")
  #print(h1)
  #cat("\n")
  #print("Relative Risk Calibration")
  #print(h2)
  #cat("\n")
  #print("Model Discrimination")
  #print(paste("Estimate of AUC:",round(auc,3)))
  #print(paste("95% CI of AUC: (",round(lower.limit.auc,3),",",
  #            round(upper.limit.auc,3),")"))
  #cat("\n")
  #print("Overall Expected to Observed Ratio")
  #print(paste("Estimate:",round(exp_by_obs,3)))
  #print(paste("95% CI:","(",round(CI_exp_by_obs[1],3),
  #            ",",round(CI_exp_by_obs[2],3),")"))
  
  
  results.cat = data.frame(cbind(levels(linear.predictor.cat),
                                 observed.prob.cat, predicted.prob.cat,
                                 lower.limit.absrisk.cat, upper.limit.absrisk.cat,
                                 observed.RR.cat, predicted.RR.cat,
                                 lower.limit.RR.cat, upper.limit.RR.cat
                                 ))
  
  names(results.cat) = c("Categories","Observed_Absolute_Risk","Predicted_Absolute_Risk",
                         "CI_Absolute_Risk_Lower","CI_Absolute_Risk_Upper",
                         "Observed_Relative_Risk","Predicted_Relative_Risk",
                         "CI_Relative_Risk_Lower","CI_Relative_Risk_Upper")
  
  ret <- list(Subject_Specific_Observed_Outcome = observed.outcome,
              Risk_Prediction_Interval = timeframe,
              Adjusted_Followup = followup, 
              Subject_Specific_Predicted_Absolute_Risk = predicted.risk, 
              Reference_Absolute_Risk = predicted.risk.ref,
              Subject_Specific_Risk_Score = linear.predictor, 
              Reference_Risk_Score = linear.predictor.ref,
              Population_Incidence_Rate = model.disease.incidence.rates,
              Study_Incidence_Rate = study_incidence,
              Category_Results = results.cat, 
              Category_Specific_Observed_Absolute_Risk = observed.prob.cat,
              Category_Specific_Predicted_Absolute_Risk = 
                predicted.prob.cat,
              Category_Specific_Observed_Relative_Risk = observed.RR.cat,
              Category_Specific_Predicted_Relative_Risk = predicted.RR.cat,
              Variance_Matrix_Absolute_Risk = variance.mat.absrisk,
              Variance_Matrix_LogRelative_Risk = varcov.logRR.cat,
              Hosmer_Lemeshow_Results = h1, 
              HL_pvalue = PVAL_absrisk, RR_test_result = h2, 
              RR_test_pvalue = PVAL_logRR,
              AUC = auc, Variance_AUC = var_auc, 
              CI_AUC = c(lower.limit.auc,upper.limit.auc),
              Overall_Expected_to_Observed_Ratio = exp_by_obs,
              CI_Overall_Expected_to_Observed_Ratio = CI_exp_by_obs,
              input.args=list(model.formula=model.formula, study.data=study.data,
                 dataset=dataset, model.name=model.name))
  
  class(ret) <- "icareValid"
  ret
}

weighted.quantcut = function (x, weights = NULL, q = 10, na.rm = TRUE, ...){
    
    if (length(q) == 1){ 
        q <- seq(0, 1, length.out = q + 1)
    }
    
    quant <- wtd.quantile(x, weights, q, type = "i/n", normwt = FALSE, 
                          na.rm = na.rm)
    dups <- duplicated(quant)
    if (any(dups)) {
        flag <- x %in% unique(quant[dups])
        retval <- ifelse(flag, paste("[", as.character(x), "]", 
            sep = ""), NA)
        uniqs <- unique(quant)

        reposition <- function(cut) {
            flag <- x >= cut
            if (sum(flag) == 0) 
                return(cut)
            else return(min(x[flag], na.rm = na.rm))
        }
        newquant <- sapply(uniqs, reposition)
        retval[!flag] <- as.character(cut(x[!flag], breaks = newquant, 
            include.lowest = TRUE, ...))
        levs <- unique(retval[order(x)])
        retval <- factor(retval, levels = levs)
        mkpairs <- function(x) sapply(x, function(y) if (length(y) == 
            2) 
            y[c(2, 2)]
        else y[2:3])
        levs  <- levs[!is.na(levs)]
        pairs <- mkpairs(strsplit(levs, "[^0-9+\\.\\-]+"))
        rownames(pairs) <- c("lower.bound", "upper.bound")
        colnames(pairs) <- levs
        closed.lower <- rep(F, ncol(pairs))
        closed.upper <- rep(T, ncol(pairs))
        closed.lower[1] <- TRUE
        for (i in 2:ncol(pairs)) if (pairs[1, i] == pairs[1, 
            i - 1] && pairs[1, i] == pairs[2, i - 1]) 
            closed.lower[i] <- FALSE
        for (i in 1:(ncol(pairs) - 1)) if (pairs[2, i] == pairs[1, 
            i + 1] && pairs[2, i] == pairs[2, i + 1]) 
            closed.upper[i] <- FALSE
        levs <- ifelse(pairs[1, ] == pairs[2, ], pairs[1, ], 
            paste(ifelse(closed.lower, "[", "("), pairs[1, ], 
                ",", pairs[2, ], ifelse(closed.upper, "]", ")"), 
                sep = ""))
        levels(retval) <- levs
    }
    else retval <- cut(x, quant, include.lowest = TRUE, ...)
    return(retval)
}

plotModelValidation = function(study.data, validation.results,
                                 dataset = "Example Dataset",
                                 model.name = "Example Model",
                                 x.lim.absrisk = NULL,
                                 y.lim.absrisk = NULL, 
                                 x.lab.absrisk = "Expected Absolute Risk (%)", 
                                 y.lab.absrisk = "Observed Absolute Risk (%)", 
                                 x.lim.RR = NULL,
                                 y.lim.RR = NULL, 
                                 x.lab.RR = "Expected Relative Risk", 
                                 y.lab.RR = "Observed Relative Risk",
                                 risk.score.plot.kernel = "gaussian",
                                 risk.score.plot.bandwidth = "nrd0",
                                 risk.score.plot.percent.smooth = 50){
  
  Category.Results = validation.results$Category_Results
  followup = validation.results$Adjusted_Followup
  observed.outcome = validation.results$Subject_Specific_Observed_Outcome
  linear.predictor = validation.results$Subject_Specific_Risk_Score
  model.disease.incidence.rates = 
    validation.results$Population_Incidence_Rate
  timeframe = validation.results$Risk_Prediction_Interval
  ages = validation.results$Study_Incidence_Rate[,1]
  study_incidence = validation.results$Study_Incidence_Rate[,2]
  
  study.entry.age = study.data$study.entry.age
  sampling.weights = study.data$sampling.weights
  
  if(is.null(sampling.weights))
    sampling.weights = rep(1,dim(study.data)[1])
  
  freq = 1/sampling.weights
  
  freq.cases = freq[observed.outcome == 1]
  freq.controls = freq[observed.outcome == 0]
  
  linear.predictor.cases = 
    linear.predictor[observed.outcome == 1]
  
  linear.predictor.controls = 
    linear.predictor[observed.outcome == 0]
  
  observed.prob.cat = 
    as.numeric(as.character(Category.Results$Observed_Absolute_Risk))
  predicted.prob.cat = 
    as.numeric(as.character(Category.Results$Predicted_Absolute_Risk))
  
  observed.RR.cat = 
    as.numeric(as.character(Category.Results$Observed_Relative_Risk))
  predicted.RR.cat = 
    as.numeric(as.character(Category.Results$Predicted_Relative_Risk))
  
  upper.limit.absrisk.cat = 
    as.numeric(as.character(Category.Results$CI_Absolute_Risk_Upper))
  lower.limit.absrisk.cat = 
    as.numeric(as.character(Category.Results$CI_Absolute_Risk_Lower))
  
  upper.limit.RR.cat = 
    as.numeric(as.character(Category.Results$CI_Relative_Risk_Upper))
  lower.limit.RR.cat = 
    as.numeric(as.character(Category.Results$CI_Relative_Risk_Lower))
  
  chisq.absrisk = 
    as.numeric(validation.results$Hosmer_Lemeshow_Results[["statistic"]])
  number.of.percentiles = 
    as.numeric(validation.results$Hosmer_Lemeshow_Results[["parameter"]])
  PVAL_absrisk = 
    as.numeric(validation.results$Hosmer_Lemeshow_Results[["p.value"]])
  
  chisq.logRR = as.numeric(validation.results$RR_test_result[["statistic"]])
  PVAL_logRR = as.numeric(validation.results$RR_test_result[["p.value"]])
  
  auc = validation.results$AUC
  upper.limit.auc = validation.results$CI_AUC[2]
  lower.limit.auc = validation.results$CI_AUC[1]
  
  Overall.Expected.to.Observed.Ratio = 
    validation.results$Overall_Expected_to_Observed_Ratio
  CI.Overall.Expected.to.Observed.Ratio = 
    validation.results$CI_Overall_Expected_to_Observed_Ratio

  
  #adjust x and y axes to draw absolute risk calibration plot
    if(length(x.lim.absrisk) < 2){
      x.lim.absrisk = c(min(lower.limit.absrisk.cat, 
                            predicted.prob.cat)*100,
                        max(upper.limit.absrisk.cat, 
                            predicted.prob.cat)*100)
    }
    
    if(length(y.lim.absrisk) < 2){
      y.lim.absrisk = c(min(lower.limit.absrisk.cat, 
                            predicted.prob.cat)*100,
                        max(upper.limit.absrisk.cat, 
                            predicted.prob.cat)*100)
    }
    
    #adjust x and y axes to draw relative risk calibration plot
    if(length(x.lim.RR) < 2){
      x.lim.RR = c(min(lower.limit.RR.cat, predicted.RR.cat),
                   max(upper.limit.RR.cat, predicted.RR.cat))
    }
    
    if(length(y.lim.RR) < 2){
      y.lim.RR = c(min(lower.limit.RR.cat, predicted.RR.cat),
                   max(upper.limit.RR.cat, predicted.RR.cat))
    }
    
  m <- rbind(c(1,2),c(3,4),c(5,5))
  layout(m)
  
  oldpar =  par(mar = rep(4,4))
  plotCI(predicted.prob.cat*100, observed.prob.cat *100, ui =
           upper.limit.absrisk.cat * 100, 
         li = lower.limit.absrisk.cat * 100, 
         xlab = x.lab.absrisk, ylab = y.lab.absrisk, xlim = 
           x.lim.absrisk, 
         ylim = y.lim.absrisk, col = "red",
         pch = 16, pty = "s",cex.lab = 1.2)
  abline(0,1,lty=2,col="black")
  #mtext(paste("Hosmer-Lemeshow p-value:",PVAL),side = 3, line = 1)
  mtext("Absolute Risk Calibration",side = 3, line =
          1,font = 4)
  
  plotCI(predicted.RR.cat, observed.RR.cat, ui =
           upper.limit.RR.cat, li = lower.limit.RR.cat, xlab =
           x.lab.RR, ylab = y.lab.RR, xlim = x.lim.RR, ylim = y.lim.RR,
         col = "red",
         pch = 16, pty = "s",cex.lab = 1.2)
  abline(0,1,lty=2,col="black")
  #mtext(paste("Hosmer-Lemeshow p-value:",PVAL),side = 3, line = 1)
  mtext("Relative Risk Calibration",side = 3, line =
          1,font = 4)
  
  plot(density(linear.predictor.controls, 
               bw = risk.score.plot.bandwidth,
               kernel = risk.score.plot.kernel,
               weights = freq.controls/sum(freq.controls),
               n = (risk.score.plot.percent.smooth/100) * 
                 length(linear.predictor.controls)),
       xlim = c(min(density(linear.predictor.controls)$x,
                    density(linear.predictor.cases)$x),
                max(density(linear.predictor.controls)$x,
                    density(linear.predictor.cases)$x)),
       ylim = c(min(density(linear.predictor.controls)$y,
                    density(linear.predictor.cases)$y),
                max(density(linear.predictor.controls)$y,
                    density(linear.predictor.cases)$y)),
       main = "",xlab = "Risk Score", pty = "s",cex.lab = 1.2)
  lines(density(linear.predictor.cases,
                bw = risk.score.plot.bandwidth,
                kernel = risk.score.plot.kernel,
                weights = freq.cases/sum(freq.cases),
                n = (risk.score.plot.percent.smooth/100) * 
                 length(linear.predictor.controls)),col = "red")
  legend("topright",legend = c("Controls","Cases"),
         col = c("black","red"),
         lty=1,cex = 0.75,y.intersp = 2, xjust = 0.5, 
         yjust = 0.5, x.intersp = 3,adj = c(0.3,0.6))
  mtext("Discrimination",side = 3, line =
          1,font = 4)
  
  plot(model.disease.incidence.rates[,1],
       predict(loess(model.disease.incidence.rates[,2] 
                     ~ model.disease.incidence.rates[,1])),
       log = "y", type = "l",
       xlim = c(18,max(model.disease.incidence.rates[,1])),
       ylim = c(min(model.disease.incidence.rates[,2],study_incidence) + 
                  10^(-10),
                max(model.disease.incidence.rates[,2],study_incidence)),
       xlab = "Age (in years)",ylab = "Incidence Rates",pty = "s",
       cex.lab = 1.2)
  lines(ages,study_incidence,col = "red")
  legend("bottomright",legend = c("Population","Study"),
         col = c("black","red"),lty=1,
         cex = 0.9,y.intersp = 2, x.intersp = 3,adj = c(0.3,0.6))
  mtext("Incidence Rates",side = 3, line =
          1,font = 4)
  
  plot(c(0,2),c(0,2),xaxt='n',yaxt='n',bty='n',pch='',ylab='',xlab='',
       pty = "s")
  text(x = 1, y = 2, paste("Dataset: ",dataset))
  text(x = 1, y = 1.8, paste("Model Name: ",model.name))
  text(x = 1, y = 1.6, paste("Risk Prediction Interval:",timeframe))
  text(x = 1, y = 1.4, paste("Number of subjects (cases): "
                             ,length(observed.outcome),"(",
                             sum(observed.outcome),")"))
  text(x = 1, y = 1.2, paste("Follow-up time (years) [mean,range]: [",
                               round(mean(followup),3),", (",
                               round(range(followup)[1],3),",",
                               round(range(followup)[2],3),")"," ]"))
  text(x = 1, y = 1, paste("Baseline age (years) [mean,range]: [", 
                               round(mean(study.entry.age),3),", (",
                               round(range(study.entry.age)[1],3),",",
                               round(range(study.entry.age)[2],3),")"," ]"))
  text(x = 1, y = 0.8, paste("E/O [Estimate, 95% CI]:"
                             ,"[",round(Overall.Expected.to.Observed.Ratio,3),
                             ",","(",round(CI.Overall.Expected.to.Observed.Ratio[1],3),
                             ",",round(CI.Overall.Expected.to.Observed.Ratio[2],3),
                             ")","]"))
  text(x = 0.3, y = 0.5, "Absolute Risk Calibration", font = 2)
  text(x = 0.3, y = 0.3, paste("HL Test, df:",
                               round(chisq.absrisk,3),",",number.of.percentiles))
  #text(x = 0.1, y = 0.3, paste("HL Test df:",number.of.percentiles))
  text(x = 0.3, y = 0.1, paste("p-value:",
                               noquote(ifelse(PVAL_absrisk < 
                                                1e-100,"<1e-100",
                                              format.pval(PVAL_absrisk,digits=20, 
                                                          eps = 1e-100, 
                                                          scientific = TRUE)))))
  text(x = 1, y = 0.5, "Relative Risk Calibration", font = 2)
  text(x = 1, y = 0.3, paste("Test, df:",
                               round(chisq.logRR,3),",",number.of.percentiles - 1))
  #text(x = 0.5, y = 0.3, paste("Test df:",number.of.percentiles - 1))
  text(x = 1, y = 0.1, paste("p-value:",
                               noquote(ifelse(PVAL_logRR < 1e-100,"<1e-100",
                                              format.pval(PVAL_logRR,digits=20,
                                                         eps = 1e-100, 
                                                          scientific = 
                                                           TRUE)))))
  text(x = 1.7, y = 0.5, "Model Discrimination", font = 2)
  text(x = 1.7, y = 0.3, paste("AUC est:",round(auc,3)))
  text(x = 1.7, y = 0.1, 
       paste("95% CI: (",round(lower.limit.auc,3),",",
             round(upper.limit.auc,3),")"))
  
  par(oldpar)
  
  NULL
}
