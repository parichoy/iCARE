\name{ModelValidation}
\alias{ModelValidation}
\title{Model Validation}
\usage{
ModelValidation(study.data, 
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
                model.name = "Example Risk Prediction Model")            
}
\arguments{
  \item{study.data}{Data frame which includes the variables below.
    \itemize{
      \item \code{observed.outcome}: 1 if disease has occurred by the end of followup, 0 if censored
      \item \code{study.entry.age}: age (in years) of entering the cohort
      \item \code{study.exit.age}: age (in years) of last followup visit
      \item \code{time.of.onset}: time (in years) of onset of disease; note that all subjects are disease free at 
          the time of entry and for those who do not develop disease by end of followup it is Inf
      \item \code{sampling.weights}: for a case-control study nested within a cohort study, this is a vector
          of sampling weights for each subject, i.e., probability of inclusion into the sample
    }
  }
  \item{total.followup.validation}{logical; TRUE if risk validation is performed over the total followup, 
      for all other cases (e.g., 5 year or 10 year risk validation) it is FALSE}
  \item{predicted.risk}{vector of predicted risks; should be supplied if risk prediction is done by some 
       method other than that implemented in \code{iCARE}; default is NULL}
  \item{predicted.risk.interval}{scalar or vector denoting the number of years after entering 
       the study over which risk validation is desired (e.g., 5 for validating a model for 5 year risk) 
       if \code{total.followup.validation = FALSE}; 
        if \code{total.followup.validation = TRUE},
        it can be set to NULL}
  \item{linear.predictor}{vector of risk scores for each subject, i.e. \code{x*beta}, where \code{x} is the 
       vector of risk factors and \code{beta} is the vector of log relative risks; in the current version if 
       both the arguments \code{predicted.risk} and \code{linear.predictor} are supplied the function will use the 
       supplied estimates to perform model validation, otherwise the function will compute these estimates
       using the \code{computeAbsoluteRisk} function}
  \item{iCARE.model.object}{A named list containing the input arguments to the function \code{\link{computeAbsoluteRisk}}.
                   The names in this list must match the argument names.
                   See \code{\link{computeAbsoluteRisk}}
                   }
  \item{number.of.percentiles}{the number of percentiles of the risk score that determines the number 
      of strata over which the risk prediction model is to be validated, default = 10}
  \item{reference.entry.age}{age of entry to be specified for computing absolute risk of the reference population}
  \item{reference.exit.age}{age of exit to be specified for computing absolute risk of the reference population}
  \item{predicted.risk.ref}{predicted absolute risk in the reference population assuming the entry age to be 
       as specified in \code{reference.entry.age} and exit age to be as specified in \code{reference.exit.age}}
  \item{linear.predictor.ref}{vector of risk scores for the reference population}
  \item{linear.predictor.cutoffs}{user specified cut-points for the linear predictor to define categories for
       absolute risk calibration and relative risk calibration}
  \item{dataset}{name and type of dataset to be displayed in the output, e.g.,
        "PLCO Full Cohort" or "Full Cohort Simulation"}
  \item{model.name}{name of the model to be displayed in output, e.g., "Synthetic Model" or "Simulation Setting"}
}

\value{
This function returns a list of the following objects: \cr
  \itemize{
    \item \code{Subject_Specific_Observed_Outcome}: observed outcome after adjusting the observed  followup according 
        to the risk prediction interval: 1 if disease has occurred by the end of followup, 0 if censored 
    \item \code{Risk_Prediction_Interval}: Character object showing the interval of risk prediction 
              (e.g., 5 years). If the risk prediction is over the total followup of the study,
              this reads "Observed Followup" 
    \item \code{Adjusted.Followup}: followup time (in years) after adjusting the observed followup according to the risk 
           prediction interval 
    \item \code{Subject_Specific_Predicted_Absolute_Risk}: predicted absolute risk of disease for each subject 
    \item \code{Reference_Absolute_Risk}: predicted absolute risk in the reference population 
    \item \code{Subject_Specific_Risk_Score}: estimated risk score for each subject; the missing covariates are handled 
                internally using the imputation in \code{iCARE} 
    \item \code{Reference_Risk_Score}: risk score for the reference population 
    \item \code{Population_Incidence_Rate}: age specific disease incidence rate in the population 
    \item \code{Study_Incidence_Rate}: estimated age specific incidence rate in the study 
    \item \code{Category_Results}: observed and predicted absolute risks and observed and predicted relative risks 
          in each category defined by the risk score 
    \item \code{Category_Specific_Observed_Absolute_Risk}: Observed absolute risk in each category defined by the risk score 
    \item \code{Category_Specific_Predicted_Absolute_Risk}: Predicted absolute risk in each category defined by the risk score 
    \item \code{Category_Specific_Observed_Relative_Risk}: Observed relative risk in each category defined by the risk score 
    \item \code{Category_Specific_Predicted_Relative_Risk}: Predicted relative risk in each category defined by the risk score 
    \item \code{Variance_Matrix_Absolute_Risk}: Variance-covariance matrix of the vector of cateogry specific absolute risks 
    \item \code{Variance_Matrix_LogRelative_Risk}: Variance-covariance matrix of the vector of cateogry specific relative risks 
    \item \code{Hosmer_Lemeshow_Results}: results of the Hosmer-Lemeshow type chisquare test comparing the observed and predicted absolute risks 
    \item \code{HL_pvalue}: pvalue of the Hosmer-Lemeshow type chisquare test 
    \item \code{RR_test_result}: results of the chisquare test comparing the observed and predicted relative risks 
    \item \code{RR_test_pvalue}: pvalue of the chisquare test of relative risk 
    \item \code{AUC}: estimate of the Area Under the Curve (AUC) defined as the probability that for a randomly sampled 
            case-control pair the case has a higher risk score than the control; for the full cohort setting we compute 
            the empirical proportion and for the nested case-control setting we compute the inverse probability weighted estimator 
    \item \code{Variance_AUC}: estimate of the variance of Area Under the Curve (AUC): for the full cohort setting the regular
        asymptotic variance is estimated and for the nested case-control setting the influence function based variance 
        estimate of the inverse probability weighted variance estimator is computed 
    \item \code{CI_AUC}: 95 percent Wald based confidence interval of Area Under the Curve (AUC) using the asymptotic variance 
    \item \code{Overall_Expected_to_Observed_Ratio}: The overall ratio of the expected risk to the observed risk 
    \item \code{CI_Overall_Expected_to_Observed_Ratio}: 95 percent Wald based confidence interval of the overall ratio of the expected risk to the observed risk 
  }
}
\description{
This function is used to validate absolute risk models.
}
%\details{}
\seealso{\code{\link{computeAbsoluteRisk}}}
\examples{

data(bc_data, package="iCARE")
validation.cohort.data$inclusion = 0
subjects_included = intersect(validation.cohort.data$id, 
                              validation.nested.case.control.data$id)
validation.cohort.data$inclusion[subjects_included] = 1

validation.cohort.data$observed.followup = validation.cohort.data$study.exit.age - 
  validation.cohort.data$study.entry.age

selection.model = glm(inclusion ~ observed.outcome 
                      * (study.entry.age + observed.followup), 
                      data = validation.cohort.data, 
                      family = binomial(link = "logit"))

validation.nested.case.control.data$sampling.weights =
  selection.model$fitted.values[validation.cohort.data$inclusion == 1]

set.seed(50)

data = validation.nested.case.control.data

snpDat     = bc_72_snps
form       = diagnosis ~ famhist + as.factor(parity)
info       = list(bc_model_cov_info[[1]], bc_model_cov_info[[3]])
vars       = all.vars(form)[-1]
risk.model = list(model.formula = form,
                  model.cov.info = info,
                  model.snp.info = snpDat,
                  model.log.RR = bc_model_log_or[c(1, 8:11)],
                  model.ref.dataset = ref_cov_dat[, vars],
                  model.ref.dataset.weights = NULL,
                  model.disease.incidence.rates = bc_inc,
                  model.competing.incidence.rates = mort_inc,
                  model.bin.fh.name = "famhist",
                  apply.cov.profile = data[,vars],
                  apply.snp.profile = data[,snpDat$snp.name],
                  n.imp = 5, use.c.code = 1, return.lp = TRUE,
                  return.refs.risk = TRUE)

# Not run since it can take a few minutes
# output = ModelValidation(study.data = data, total.followup.validation = TRUE,
#      predicted.risk.interval = NULL, iCARE.model.object = risk.model,
#      number.of.percentiles = 10)
output

}

