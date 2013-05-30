\name{JMLE}
\alias{JMLE}

\title{Joint maximum likelihood estimation of item parameters and examinee attribute profiles}
\description{
  This function returns joint maximum likelihood estimates of item parameters and examinee attribute profiles in cognitive diagnostic models. The algorithm starts from the nonparametric estimation of attribute profiles, implemented by the \code{AlphaNP} function, and then iteratively estimates item parameters and attribute profiles using conditional maximum likelihood estimation until the algorithm converges. Currently supported models include the DINA model, the DINO model, and the NIDA model.}

\usage{
JMLE(Y, Q, model = c("DINA", "DINO", "NIDA"), NP.method = c("Weighted", 
"Hamming", "Penalized"), wg=1, ws=1, conv.crit.par = 0.001, 
conv.crit.att = 0.01, max.ite = 100)
}
\arguments{
  \item{Y}{A matrix of binary responses. Rows represent persons and columns represent items. 1=correct, 0=incorrect.
}
  \item{Q}{The Q-matrix of the test. Rows represent items and columns represent attributes. 1=attribute required by the item, 0=attribute not required by the item.
}
  \item{model}{Currently support three models: \code{"DINA"}, \code{"DINO"}, and \code{"NIDA"}. The default is \code{"DINA"}.
}
  \item{NP.method}{The method of the nonparametric estimation in the initial stage. \code{"Hamming"}: the plain Hamming distance method; \code{"Weighted"}: the Hamming distance weighted by inversed item variance; \code{"Penalized"}: the Hamming distance weighted by inversed item variance and specified penalizing weights for guess and slip. The default is \code{"Weighted"}.
}
  \item{wg}{Additional argument for the "penalized" NP.method. \code{wg} is the weight assigned to guessing in the DINA or DINO models. A large value of \code{wg} results in a stronger impact on Hamming distance (larger loss function values) caused by guessing.
}
  \item{ws}{Additional input for the "penalized" NP.method. \code{ws} is the weight assigned to slipping in the DINA or DINO models. A large value of \code{ws} results in la stronger impact on Hamming distance (larger loss function values) caused by slipping.
}
  \item{conv.crit.par}{The critical value for the maximum absolute change in all item parameters values to determine convergence.
}
  \item{conv.crit.att}{The critical value for the percentage of examinee attribute profiles that are changed to determine convergence.
}
  \item{max.ite}{The maximum number of iterations allowed.
}
}

\value{
  \item{alpha.est}{JMLE estimates of examinee attribute profiles. Rows represent persons and columns represent attributes. 1=examinee masters the attribute, 0=examinee does not master the attribute.}
  \item{par.est}{JMLE estimates of item parameters, including \code{par.est$slip}, \code{par.est$guess}, \code{par.est$se.slip}, and \code{par.est$se.guess}.}
  \item{n.tie}{Number of ties in the final log-likelihood among the candidate attribute profiles for each person. When we encounter ties, one of the tied attribute profiles is randomly chosen.}
  \item{undefined.flag}{A binary vector indicating whether the parameters of each item are undefined. 1=undefined, 0=defined.}
  \item{loglike}{The final overall log-likelihood value from the estimated item parameters and attribute profiles based on the specified model.}
  \item{convergence}{A message on whether the algorithm converged.}
  \item{n.ite}{Number of iterations performed.}
  \item{loglike.matrix}{The values for the log-likelihood function in the last iteration for each candidate attribute profile by each person. Rows represent candidate attribute profiles in the same order with the pattern matrix; columns represent different examinees.}
  \item{est.class}{The final class number (row index in \code{pattern}) for each person's attribute profile. It can also be used for locating the log-likelihood value in \code{loglike.matrix} for the estimated attribute profile for each person.}
  \item{NP.loss.matrix}{The values for the loss function of the nonparametric estimation of Alpha. Rows represent candidate attribute profiles in the same order with the pattern matrix; columns represent different examinees.}
  \item{NP.alpha.est}{The estimates of examinee attribute profiles from the initial nonparameteric estimation.}
  \item{NP.est.class}{The class number (row index in \code{pattern}) for each person's attribute profile from the initial nonparametric classification. It can also be used for locating the loss function value in \code{NP.loss.matrix} for the estimated attribute profile for each person.}
  \item{pattern}{All possible attribute profiles in the search space.}
}

\seealso{
\code{\link{AlphaMLE}}, \code{\link{AlphaNP}}, \code{\link{ParMLE}}, \code{\link{print.JMLE}},  \code{\link{plot.JMLE}}
}

\references{
Chiu, C. (2011). \emph{Flexible approaches to cognitive diagnosis: nonparametric methods and small sample techniques}. Invited session of cognitive diagnosis and item response theory at 2011 Joint Statistical Meeting.
Chiu, C. Y., & Douglas, J. A. (In press). A model-free approach to cognitive diagnosis by proximity
to ideal response patterns. \emph{Journal of Classification}.
}

\examples{
data("Data.DINA")
JMLE.result <- JMLE(Data.DINA$response, Data.DINA$Q, model="DINA", conv.crit.par=0.001, 
conv.crit.att=0.001, max.ite=100)
print(JMLE.result)  # Print the estimated item parameters, standard errors, 
#and examinee attribute profiles
plot(JMLE.result, nperson=1)   # Plot the sorted loss function of different 
#attribute profiles for this examinee
}
