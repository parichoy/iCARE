\name{plotModelValidation}
\alias{plotModelValidation}
\alias{plot}
\title{Model Validation Plot}
\usage{
plotModelValidation(study.data, validation.results,
                    dataset = "Example Dataset",
                    model.name = "Example Model",
                    x.lim.absrisk = "",
                    y.lim.absrisk = "", 
                    x.lab.absrisk = "Expected Absolute Risk (\%)", 
                    y.lab.absrisk = "Observed Absolute Risk (\%)", 
                    x.lim.RR = "",
                    y.lim.RR = "", x.lab.RR = "Expected Relative Risk", 
                    y.lab.RR = "Observed Relative Risk",
                    risk.score.plot.kernel = "gaussian",
                    risk.score.plot.bandwidth = "nrd0",
                    risk.score.plot.percent.smooth = 50)
}
\arguments{
  \item{study.data}{See \code{\link{ModelValidation}} }
  \item{validation.results}{List returned from \code{\link{ModelValidation}} }
  \item{dataset}{Name and type of dataset to be displayed in the output, 
           e.g., "PLCO Full Cohort" or "Full Cohort Simulation"}
  \item{model.name}{Name of the model to be displayed in output, 
                     e.g., "Synthetic Model" or "Simulation Setting"}
  \item{x.lim.absrisk}{ Vector of length two specifying the x-axes limits in the
                               absolute risk calibration plot. If not specified, then default
                               limits will be computed.}
  \item{y.lim.absrisk}{ Vector of length two specifying the y-axes limits in the
                               absolute risk calibration plot. If not specified, then default
                               limits will be computed.}
  \item{x.lab.absrisk}{ String specifying the x-axes label in the  
                               absolute risk calibration plot.
                               The default is "Expected Absolute Risk (\%)".}
  \item{y.lab.absrisk}{ String specifying the y-axes label in the  
                               absolute risk calibration plot.
                               The default is "Observed Absolute Risk (\%)."}
  \item{x.lim.RR}{ Vector of length two specifying the x-axes limits in the
                               relative risk calibration plot. If not specified, then default
                               limits will be computed.}
  \item{y.lim.RR}{ Vector of length two specifying the y-axes limits in the
                               relative risk calibration plot. If not specified, then default
                               limits will be computed.}
  \item{x.lab.RR}{ String specifying the x-axes label in the  
                               relative risk calibration plot.
                               The default is "Expected Relative Risk".}
  \item{y.lab.RR}{ String specifying the y-axes label in the  
                               relative risk calibration plot.
                               The default is "Observed Relative Risk".}
  \item{risk.score.plot.kernel}{Character string giving the smoothing kernel to be used by the density function used 
         internally to plot the density of the risk scores. It should be one of "gaussian", "rectangular", "triangular",
          "epanechnikov", "biweight", "cosine" or "optcosine", with default "gaussian".}
  \item{risk.score.plot.bandwidth}{The options for bandwidth selection used by the density function internally to plot 
        the density of the risk scores. The options are one of the following: "nrd0", "nrd", "ucv", "bcv", "SJ-ste", "SJ-dpi"
        with the default being "nrd0". More information on these different options is available in the help pages that can be 
       accessed from R using the command \code{?bw.nrd}.}
  \item{risk.score.plot.percent.smooth}{Percentage of the number of sample points used for determining the number of equally 
      spaced points at which the density of the risk score is to be estimated. This number supplies the input parameter "n" to the 
      density function used internally to plot the densities of the risk score. The default value is 50.}  
}

\value{
This function returns NULL \cr
}
\description{
This function is used to create plots for model calibration, model discrimination and incidence rates.
}
%\details{}
\seealso{\code{\link{ModelValidation}}}
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
#output = ModelValidation(study.data = data, total.followup.validation = TRUE,
#          predicted.risk.interval = NULL, iCARE.model.object = risk.model,
#          number.of.percentiles = 10)

plot(output)

}

