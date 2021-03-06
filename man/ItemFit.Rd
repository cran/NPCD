\name{ItemFit}
\alias{ItemFit}

\title{Compute item fit statistics for outputs generated by estimation functions in the package}

\description{
  This function computes item fit statistics for outputs generated by estimation functions in the package, including \code{\link{AlphaNP}}, \code{\link{AlphaMLE}}, \code{\link{ParMLE}}, and\code{\link{JMLE}}. The function currently provides the RMSEA and Chi-square item fit statistics.
}

\usage{ItemFit(x, model=NULL, par=NULL)
}

\arguments{
  \item{x}{The output from the function (The list of all outputs).}
  \item{model}{This needs to be additionally specified only when \code{x} is output from \code{AlphaNP}. Currently support five models: \code{"DINA"}, \code{"DINO"}, \code{"NIDA"}, \code{"GNIDA"}, and \code{"RRUM"}. 
}
  \item{par}{This needs to be additionally specified only when \code{x} is output from \code{AlphaNP}. A list of parameters.                                            
          DINA & DINO --- \code{par$slip}: a vector of slipping parameters for each item;       
                   \code{par$guess}: a vector of guessing parameters for each item.  
          NIDA --- \code{par$slip}: a vector of slipping parameters for each attribute; 
                   \code{par$guess}: a vector of guessing parameters for each attribute.
          GNIDA --- \code{par$slip}: a matrix (items by attributes) of slipping parameters; 
                   \code{par$guess}: a matrix (items by attributes) of guessing parameters.
          RRUM --- \code{par$pi}: a vector of pi parameters for each item;
                   \code{par$r}: a matrix (items by attributes) of r parameters.                   
}
}

\value{
  \item{RMSEA}{The model-based root mean square error of approximation (Kunina-Habenicht et al., 2012) of each item based on the estimated or given item parameter, Q-vector, and alpha matrix.}
  \item{Chisq}{The Q1 Chi-square statistic (Wang et al., 2015; Yen, 1981) of each item based on the estimated or given item parameter, Q-vector, and alpha matrix.}
  \item{Chisq.p}{The p-values for the Chi-square statistic for each item.}
  \item{Chisq.df}{The degrees of freedom for the Chi-square statistic for each item.}
}

\seealso{
\code{\link{AlphaNP}}, \code{\link{AlphaMLE}}, \code{\link{ParMLE}}, \code{\link{JMLE}}
}

\references{
  Kunina-Habenicht, O., Rupp, A. A., & Wilhelm, O. (2012). The Impact of Model Misspecification on Parameter Estimation and Item-Fit Assessment in Log0Linear Diagnostic Classification Models. \emph{Journal of Educational Measurement, 49}(1), 59-81.
  
  Wang, C., Shu, Z., Shagn, Z., & Xu, G. (2015). Assessing Item-Level Fit for the DINA Model. \emph{Applied Psychological Measurement}, 1-14.
  
  Yen, W. M. (1981). Using Simulation Results to Choose a Latent Trait Model. \emph{Applied Psychological Measurement, 5}, 245-262.
}

\examples{
# See examples in AlphaNP, AlphaMLE, ParMLE, and JMLE.
}
