% --- Source file: iCARE.Rd ---
\name{iCARE}
\alias{iCARE}
\docType{package}
\title{
A Tool for Individualized Coherent Absolute Risk Estimation (iCARE)
}
\description{
Individualized Coherent Absolute Risk Estimators (iCARE) is a tool that allows researchers to quickly 
build models for absolute risk and apply them to estimate individuals' risk based on a set of user defined 
input parameters. The software gives users the flexibility to change or update models rapidly based on new 
risk factors or tailor models to different populations based on the specification of simply three input arguments: 
(1) a model for relative risk assumed to be externally derived 
(2) an age-specific disease incidence rate and 
(3) the distribution of risk factors for the population of interest. 
The tool can handle missing information on risk factors for risk estimation using an approach where all estimates 
are derived from a single model through appropriate model averaging.
}
\details{
The main functions for building and applying an absolute risk model are \code{\link{computeAbsoluteRisk}} 
and \code{\link{computeAbsoluteRiskSplitInterval}}.  
The first of these computes absolute risks over the specified time interval using a single set of paramters.  
The second provides more advanced functionality and computes absolute risk over the interval in two parts.  

\code{\link{computeAbsoluteRiskSplitInterval}} allows the user compute absolute risk over the interval in two parts,
 incorporating two different sets of paramters before and after a specified cutpoint.  
This function allows a different cutpoint for each covariate profile if desired. 
The function for validating an absolute risk model is \code{\link{ModelValidation}}, and
   \code{\link{plotModelValidation}} can be called for producing plots 
   for model calibration, model discrimination and incidence rates.
}
%\references{ }
\author{Paige Maas, Parichoy Pal Choudhury, Nilanjan Chatterjee and William Wheeler <wheelerb@imsweb.com>}
\keyword{package}