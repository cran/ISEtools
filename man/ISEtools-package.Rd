\name{ISEtools-package}
\alias{ISEtools-package}
\alias{ISEtools}
\docType{package}
\title{
Analysis tools for ion selective electrodes (ISEs)
}
\description{
Bayesian calibration for single or multiple ISEs using R and OpenBUGS.  Estimation of analyte activities using single ISEs or ISE arrays.
}
\details{
\tabular{ll}{
Package: \tab ISEtools\cr
Type: \tab Package\cr
Version: \tab 2.5.2\cr
Depends: R (>2.13.0), R2OpenBUGS (>3.0)\cr
Date: \tab 2018-10-15\cr
License: \tab GPL-2\cr
SystemRequirements: OpenBUGS (>3.0)\cr
}
The primary funtions are loadISEdata (which loads calibration and experimental data from tab-delimited text files), describeISE (uses Bayesian calibration to estimate ISE parameters from calibration data), and analyseISE (combines calibration data with experimental data in basic or standard addition format to estimate analyte concentrations).
}
\author{
Peter Dillingham [aut, cre],
	Christina McGraw [ctb],
	Aleksandar Radu [ctb],
	Basim Alsaedi [ctb]

Maintainer: Peter Dillingham <peter.dillingham@otago.ac.nz>
}
\references{
Dillingham, P.W., Radu, T., Diamond, D., Radu, A. and McGraw, C.M. (2012). Bayesian Methods for Ion Selective Electrodes. \emph{Electroanalysis}, \strong{24}, 316-324. <doi:10.1002/elan.201100510>

Dillingham, P.W., Alsaedi, B.S.O. and McGraw, C.M. (2017). Characterising uncertainty in instrumental limits of detection when sensor response is non-linear. \emph{2017 IEEE SENSORS}, Glasgow, United Kingdom, pp. 1-3. <doi:10.1109/ICSENS.2017.8233898>
}
%\seealso{
%\link[http://www.maths.otago.ac.nz/~dillingh/software.html]
%}
\keyword{ package }
\examples{
data(LeadStdAdd)
print(LeadStdAdd)
summary(LeadStdAdd)
plot(LeadStdAdd)
\donttest{
example1 = describeISE(LeadStdAdd, Z =2, temperature=21)
print(example1)
summary(example1)
plot(example1)
example2 = analyseISE(LeadStdAdd, Z =2, temperature=21)
print(example2)
summary(example2)
plot(example2, ylim = c(-7, -3), xlab = "ID (Sample)", 
	ylab = expression(paste(log[10], " ", Pb^{paste("2","+",sep="")} )))
}
}
