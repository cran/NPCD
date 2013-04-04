\name{JMLE}
\alias{JMLE}

\title{Joint maximum likelihood estimation of item parameters and examinee attribute profiles}
\description{
  This function returns joint maximum likelihood estimates of item parameters and examinee attribute profiles in cognitive diagnostic models. The algorithm starts from the nonparametric estimation of attribute profiles, implemented by the \code{AlphaNP} function, and then iteratively estimates item parameters and attribute profiles using conditional maximum likelihood estimation until the algorithm converges. Currently supported models include the DINA model, the DINO model, and the NIDA model.}

\usage{
JMLE(Y, Q, model = "DINA", NP.method = "Weighted", conv.crit.par = 0.001, conv.crit.att = 0.01, max.ite = 100)
}
\arguments{
  \item{Y}{A matrix of binary responses. Rows represent persons and columns represent items. 1=correct, 0=incorrect.
}
  \item{Q}{The Q-matrix of the test. Rows represent items and columns represent attributes. 1=attribute required by the item, 0=attribute not required by the item.
}
  \item{model}{Currently support three models: \code{"DINA"}, \code{"DINO"}, and \code{"NIDA"}. The default is \code{"DINA"}.
}
  \item{NP.method}{The method of the nonparametric estimation in the initial stage. \code{"Hamming"}: the plain Hamming distance method; \code{"Weighted"}: the Hamming distance weighted by inversed item variance; \code{"Panelized"}: the Hamming distance weighted by inversed item variance and specified panelizing weights for guess and slip. The default is \code{"Weighted"}.
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
  \item{undefined.flag}{A binary vector indicating whether the parameters of each item are undefined. 1=undefined, 0=defined.}
  \item{convergence}{A message on whether the algorithm converged.}
  \item{n.ite}{Number of iterations performed.}
  \item{loglike.matrix}{The values for the log-likehood function in the last iteration. Rows represent candidate attribute profiles in the same order with the pattern matrix; columns represent different examinees.}
  \item{NPloss.matrix}{The values for the loss function of the nonparametric estimation of Alpha. Rows represent candidate attribute profiles in the same order with the pattern matrix; columns represent different examinees.}
  \item{alpha.est.NP}{The estimates of examinee attribute profiles from the initial nonparameteric estimation.}
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
# Generate item and examinee profiles

natt <- 3
nitem <- 4
nperson <- 5
Q <- rbind(c(1, 0, 0), c(0, 1, 0), c(0, 0, 1), c(1, 1, 1))
alpha <- rbind(c(0, 0, 0), c(1, 0, 0), c(0, 1, 0), c(0, 0, 1), c(1, 1, 1))

# Generate DINA model-based response data

slip <- c(0.1, 0.15, 0.2, 0.25)
guess <- c(0.1, 0.15, 0.2, 0.25)
my.par <- list(slip=slip, guess=guess)

data <- matrix(NA, nperson, nitem)
eta <- matrix(NA, nperson, nitem)

for (i in 1:nperson) {
  for (j in 1:nitem) {
  eta[i, j] <- prod(alpha[i,] ^ Q[j, ])
  P <- (1 - slip[j]) ^ eta[i, j] * guess[j] ^ (1 - eta[i, j])
  u <- runif(1)
  data[i, j] <- as.numeric(u < P)
  }
}

# Using the function to estimate item parameters and examinee attribute profiles

JMLE.result <- JMLE(data, Q, model="DINA", conv.crit.par=0.0001, conv.crit.att=0.001, max.ite=100)

nperson <- 1   # Choose an examinee to investigate
print(JMLE.result)  # Print the estimated item parameters, standard errors, and examinee attribute profiles
plot(JMLE.result, nperson)   # Plot the sorted loss function of different attribute profiles for this examinee
}
