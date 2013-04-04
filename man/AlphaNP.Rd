\name{AlphaNP}
\alias{AlphaNP}

\title{Nonparametric estimation of attribute profiles}
\description{
  This function estimates attribute profiles using nonparametric approaches for both the "AND gate" and the "OR gate" cognitive diagnostic models. These algorithms select the attribute profile with the smallest loss function value (plain, weighted, or panelized Hamming distance, see below for details) as the estimate. If more than one attribute profiles have the smallest loss function value, one of them is randomly chosen.
}

\usage{
AlphaNP(Y, Q, gate="AND", method = "Hamming", wg = 1, ws = 1)
}
\arguments{
  \item{Y}{A matrix of binary responses. Rows represent persons and columns represent items. 1=correct, 0=incorrect.
}
  \item{Q}{The Q-matrix of the test. Rows represent items and columns represent attributes. 1=attribute required by the item, 0=attribute not required by the item.
}
  \item{gate}{\code{"AND"}: the examinee needs to possess all required attributes of an item in order to answer it correctly; \code{"OR"}: the examinee needs to possess only one of the required attributes of an item in order to answer it correctly.
}
  \item{method}{The method of nonparametric estimation. \code{"Hamming"}: the plain Hamming distance method; \code{"Weighted"}: the Hamming distance weighted by inversed item variance; \code{"Panelized"}: the Hamming distance weighted by inversed item variance and specified panelizing weights for guess and slip.
}
  \item{wg}{Additional argument for the "panelized" method. \code{wg} is the weight assigned to guessing in the DINA or DINO models. A large value of \code{wg} results in a stronger impact on Hamming distance (larger loss function values) caused by guessing.
}
  \item{ws}{Additional input for the "panelized" method. \code{ws} is the weight assigned to slipping in the DINA or DINO models. A large value of \code{ws} results in la stronger impact on Hamming distance (larger loss function values) caused by slipping.
}
}

\value{
  \item{alpha.est}{Estimated attribute profiles. Rows represent persons and columns represent attributes. 1=examinee masters the attribute, 0=examinee does not master the attribute.}
  \item{est.ideal}{Estimated ideal response to all items by all examinees. Rows represent persons and columns represent items. 1=correct, 0=incorrect.}
  \item{n.tie}{Number of ties in the Hamming distance. When we encounter ties, one of the tied attribute profiles is randomly chosen.}
  \item{pattern}{All possible attribute profiles in the search space.}
  \item{loss.matrix}{The matrix of the values for the loss function (the plain, weighted, or panelized Hamming distance. Rows represent candidate attribute profiles in the same order with the pattern matrix; columns represent different examinees.}
}

\seealso{
\code{\link{AlphaMLE}}, \code{\link{JMLE}}, \code{\link{print.AlphaNP}}, \code{\link{plot.AlphaNP}}
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

# Using the function to estimate examinee attribute profile

alpha.est.NP.H <- AlphaNP(data, Q, gate="AND", method="Hamming")
alpha.est.NP.W <- AlphaNP(data, Q, gate="AND", method="Weighted")
alpha.est.NP.P <- AlphaNP(data, Q, gate="AND", method="Panelized", wg=2, ws=1)

nperson <- 1   # Choose an examinee to investigate
print(alpha.est.NP.H) # Print the estimated examinee attribute profiles
plot(alpha.est.NP.H, nperson)   # Plot the sorted loss function of different attribute profiles for this examinee
}
