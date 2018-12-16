test_risk <- function() {

  risks <- c(0.10250188, 0.09017017, 0.16893543)

  data(bc_data, package="iCARE")
  set.seed(1)
  results <- computeAbsoluteRisk(model.formula     = bc_model_formula, 
                                         model.cov.info    = bc_model_cov_info,
                                         model.snp.info    = bc_72_snps,
                                         model.log.RR      = bc_model_log_or,
                                         model.ref.dataset = ref_cov_dat,
                                         model.disease.incidence.rates   = bc_inc,
                                         model.competing.incidence.rates = mort_inc, 
                                         model.bin.fh.name = "famhist",
                                         apply.age.start    = 50, 
                                         apply.age.interval.length = 30,
                                         apply.cov.profile  = new_cov_prof,
                                         apply.snp.profile  = new_snp_prof, 
                                         return.refs.risk   = TRUE)
  ret <- as.numeric(results$risk)
  
  
  checkEqualsNumeric(ret, risks, tolerance=1.0e-4)
  
}
