% --- Source file: bc_data.Rd ---
\name{bc_data}
\alias{bc_data}
\alias{bc_72_snps}
\alias{bc_inc}
\alias{bc_model_formula}
\alias{bc_model_cov_info}
\alias{bc_model_log_or}
\alias{bc_model_log_or_post_50}
\alias{mort_inc}
\alias{new_cov_prof}
\alias{ref_cov_dat}
\alias{new_snp_prof}
\alias{ref_cov_dat_post_50}
\alias{validation.nested.case.control.data}
\alias{validation.cohort.data}
\alias{output}


%\docType{package}
\title{
  Data for examples
}
\description{
 Example data for \code{\link{computeAbsoluteRisk}}, \code{\link{computeAbsoluteRiskSplitInterval}}, 

 \code{\link{ModelValidation}}, and \code{\link{plotModelValidation}}.
}
\details{ 
\itemize{
\item \code{bc_model_cov_info}: a main list containing information on family history, age at menarche (years), 
   parity, age at first birth (years), age at menopause (years), height (meters), Body Mass Index (kg/sq.m.),
  use of hormone replacement therapy, use of estrogen and progesterone combined therapy, 
  use of estrogen only therapy, current use of hormone replacement therapy, alcohol (drinks/week), smoking status.; 
   information on each risk factor is given as a list 
\item \code{bc_model_formula}: formula for the specification of the models with risk factors
\item \code{bc_72_snps}: contains published SNP information from reference:
  Michailidou K, Lindstrom S, Dennis J, Beesley J, Hui S, Kar S, Lemacon A, Soucy P, Glubb D, Rostamianfar A, 
 et al. (2017) Association analysis identifies 65 new breast cancer risk loci. Nature 551:92-94
\item \code{bc_model_log_or}: vector of log-odds ratios of family history, age at menarche (years), parity, 
age at first birth (years), age at menopause (years), height (meters), Body Mass Index (kg/sq.m.), 
use of hormone replacement therapy, use of estrogen and progesterone combined therapy, 
use of estrogen only therapy, current use of hormone replacement therapy, alcohol (drinks/week), smoking status. 
\item \code{bc_model_log_or_post_50}: vector of log-odds ratios of family history, age at menarche (years), parity, age at first birth (years), 
age at menopause (years), height (meters), Body Mass Index (kg/sq.m.), use of hormone replacement therapy, 
use of estrogen and progesterone combined therapy, use of estrogen only therapy,
 current use of hormone replacement therapy, alcohol (drinks/week), smoking status. for women 50 years or older 
\item \code{ref_cov_dat}: contains individual level reference dataset of risk factors representative of the underlying population imputed using reference (4) and (5) 
\item \code{ref_cov_dat_post_50}: contains individual level reference dataset of the risk factors for women aged 50 years or older 
\item \code{bc_inc}: contains age-specific incidence rates of breast cancer from reference (3) 
\item \code{mort_inc}: contains age-specific incidence rates of all-cause mortality from reference (1) below 
\item \code{new_cov_prof}: Information on family history, age at menarche (years), parity, age at first birth (years), age at menopause (years), height (meters), 
  Body Mass Index (kg/sq.m.), use of hormone replacement therapy, use of estrogen and progesterone combined therapy, 
  use of estrogen only therapy, current use of hormone replacement therapy, alcohol (drinks/week), 
  smoking status for three women (given for illustration of absolute risk prediction) 
\item \code{new_snp_prof}: Information on 72 breast cancer associated SNPs for three women (given for illustration of absolute risk prediction) 
\item \code{validation.cohort.data}: Simulated  full cohort dataset of 50,000 women for illustration 
of model validation. The variables are: 
 \itemize{
  \item \code{id}: Subject id 
  \item \code{famhist}: Family history; binary indicator of presence/absence of disease among first degree relatives 
  \item \code{parity}: number of child births categorized as nulliparous (ref), 1 births, 2 births, 3 births, 4+ births
  \item \code{menarche_dec}: categories of age at menarche (years) with levels: less than 11,11-11.5,11.5-12,12-13(ref),13-14,14-15, greater than 15
  \item \code{birth_dec}: categories of age at first birth (years) with levels: less than 19 (ref), 19-22,22-23,23-25,25-27,27-30,30-34,34-38, greater than 38
  \item \code{agemeno_dec}: categories of age at menopause (years) with levels: less than 40 (ref), 40-45, 45-47, 47-48, 48-50, 50-51, 51-52, 52-53, 53-55, greater than 55
  \item \code{height_dec}: categories of height (meters) with levels: less than 1.55, 1.55-1.57, 1.57-1.60, 1.60-1.61, 1.61-1.63, 1.63-1.65, 1.65-1.66, 1.66-1.68, 1.68-1.71
  \item \code{bmi_dec}: categories of body mass index (kg/sq.m.) with levels: less than 21.5 (ref), 21.5-23, 23-24.2, 24.2-25.3, 25.3-26.5, 26.5-27.8, 27.8-29.3, 29.3-31.4, 31.4-34.6
  \item \code{rd_menohrt}: use of hormone replacement therapy with levels: premenopausal (ref), postmenopausal and never HRT user, postmenopausal and ever HRT user
  \item \code{rd2_everhrt_c}: binary indicator for postmenopausal and ever user of estrogen and progesterone combined therapy
  \item \code{rd2_everhrt_e}: binary indicator for postmenopausal and ever user of estrogen only therapy
  \item \code{rd2_currhrt}: binary indicator of postmenopausal and current HRT user
  \item \code{alcoholdweek_dec}: alcohol in drinks per week categorized into levels: none (ref), 0-0.4, 0.4-0.8, 0.8-1.5, 1.5-3.2, 3.2-5.7, 5.7-9.8, >9.8
  \item \code{ever_smoke}: binary indicator for ever smoker 
  \item \code{study.entry.age}: age of study entry 
  \item \code{study.exit.age}: age of study exit 
  \item \code{observed.outcome}: binary indicator of disease status (yes/no) 
  \item \code{time.of.onset}: time (in years) since study entry to the development of disease; 
   for subjects who have not developed disease beyond the observed followup, it is set to Inf 
  \item \code{observed.followup}: number of years the subject is followed up in the study 
   (difference between the age of study exit and age of study entry)
 } % END: itemize
\item \code{validation.nested.case.control.data}: A simulated example of a case-control study of 5285 women, 
  nested within the full cohort.
 In addition to the variables given above, it has information on the 72 breast cancer associated SNPs with 
variable names being the rs-identifiers.
\item \code{output}: object returned from \code{\link{computeAbsoluteRisk}}
} % END itemize
} % END: details
\examples{

 temp <- data(bc_data, package="iCARE")

 # Display the object names
 temp

}
\keyword{data}

\references{
(1)  Centers for Disease Control and Prevention (CDC), National Center for Health Statistics (NCHS). Underlying Cause of Death 1999-2011 on CDC WONDER Online Database, released 2014. Data are from the Multiple Cause of Death Files, 1999-2011, as compiled from data provided by the 57 vital statistics jurisdictions through the Vital Statistics Cooperative Program. \cr
 Accessed at http://wonder.cdc.gov/ucd-icd10.html on Aug 26, 2014.

(2)  Michailidou K, Beesley J, Lindstrom S, et al. Genome-wide association analysis of more than 120,000 individuals identifies 15 new susceptibility loci for breast cancer. Nature genetics 2015;47:373-80.

(3)  Surveillance, Epidemiology, and End Results (SEER) Program SEER*Stat Database: Incidence - SEER 18 Regs Research Data, Nov 2011 Sub, Vintage 2009 Pops (2000-2009) <Katrina/Rita Population Adjustment> - Linked To County Attributes - Total U.S., 1969-2010 Counties. In: National Cancer Institute D, Surveillance Research Program, Surveillance Systems Branch, ed. SEER18 ed.

(4)  2010 National Health Interview Survey (NHIS) Public Use Data Release, NHIS Survey Description. 2011. \cr
(Accessed at ftp://ftp.cdc.gov/pub/Health_Statistics/NCHS/Dataset_Documentation/NHIS/2010/srvydesc.pdf.)

(5)  Centers for Disease Control and Prevention (CDC). National Center for Health Statistics (NCHS). National Health and Nutrition Examination Survey Questionnaire. Hyattsville, MD: U.S. Department of Health and Human Services, Centers for Disease Control and Prevention; 2010.

}


