\name{woeBinning}
\docType{package}
\alias{woeBinning}
\title{Package for Supervised WOE Binning}
\description{
This package generates, visualizes, tabulates and deploys a supervised weight of evidence (WOE) binning of variables.
}
\details{
The package \code{woeBinning} adapts and automates the process of fine and coarse
classing and bins numeric variables and factors with respect to a dichotomous
target variable. Additionally, it visualizes the realized binning solution,
tabulates it and deploys it to (new) data. All functions can be used with single
variables or an entire data frame.
}
\section{Binning Functions}{
\itemize{ 
   \item \code{\link{woe.binning}} generates a supervised binning of numeric variables and factors. 
   \item \code{\link{woe.binning.plot}} visualizes the binning solution generated and saved via \code{woe.binning}. 
   \item \code{\link{woe.binning.table}} tabulates the binning solution generated and saved via \code{woe.binning}. 
   \item \code{\link{woe.binning.deploy}} deploys the binning solution generated and saved via \code{woe.binning} to (new) data. 
}
}
\references{
Siddiqi, N. 2006: \emph{Credit Risk Scorecards: Developing and Implementing Intelligent Credit Scoring}. Hoboken, New Jersey: John Wiley & Sons.

Anderson, R. 2007: \emph{The Credit Scoring Toolkit: Theory and Practice for Retail Credit Risk Management and Decision Automation}. Oxford / New York: Oxford University Press.
}
\examples{
# Load German credit data and create subset
data(germancredit)
df <- germancredit[, c('creditability', 'credit.amount', 'duration.in.month',
                  'savings.account.and.bonds', 'purpose')]

# Bin all variables of the data frame (apart from the target variable)
# with default parameter settings
binning <- woe.binning(df, 'creditability', df)

# Plot the binned variables
woe.binning.plot(binning)

# Tabulate the binned variables
tabulate.binning <- woe.binning.table(binning)
tabulate.binning
						  
# Deploy the binning solution to the data frame
# (i.e. add binned variables and corresponding WOE variables)
df.with.binned.vars.added <- woe.binning.deploy(df, binning,
                                               add.woe.or.dum.var='woe')		

}