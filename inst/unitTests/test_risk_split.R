test_risk_split <- function() {

  risks <- c(0.08589966, 0.07212156, 0.14428262)

  data(bc_data, package="iCARE")
  set.seed(1)
  
  results <- computeAbsoluteRiskSplitInterval(model.formula=bc_model_formula, 
                                         cut.time = 50,
                                         model.cov.info       = bc_model_cov_info,
                                         model.snp.info       = bc_72_snps,
                                         model.log.RR         = bc_model_log_or,
                                         model.log.RR.2       = bc_model_log_or_post_50,
                                         model.ref.dataset    = ref_cov_dat,
                                         model.ref.dataset.2  = ref_cov_dat_post_50,
                                         model.disease.incidence.rates   = bc_inc,
                                         model.competing.incidence.rates = mort_inc, 
                                         model.bin.fh.name = "famhist",
                                         apply.age.start    = 30, 
                                         apply.age.interval.length = 40,
                                         apply.cov.profile  = new_cov_prof,
                                         apply.snp.profile  = new_snp_prof, 
                                         return.refs.risk   = TRUE)
  ret <- as.numeric(results$risk)

  checkEqualsNumeric(ret, risks, tolerance=1.0e-4)
  
}
