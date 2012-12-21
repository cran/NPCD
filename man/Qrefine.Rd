\name{Qrefine}
\alias{Qrefine}

\title{Refine the Q-matrix by minimizing the residual sum of square (RSS)}
\description{Refine the Q-matrix by minimizing the residual sum of square (RSS) betweenn the real responses and ideal responses. Examinee attribute profiles are estimated using the nonparametric method (plain Hamming) implemented by \code{\link{AlphaNP}}.
}

\usage{
Qrefine(Y, Q, gate="AND", max.ite = 50)
}
\arguments{
  \item{Y}{A matrix of binary responses. Rows represent persons and columns represent items. 1=correct, 0=incorrect.
}
  \item{Q}{The Q-matrix of the test. Rows represent items and columns represent attributes. 1=attribute required by the item, 0=attribute not required by the item.
}
  \item{gate}{\code{"AND"}: the examinee needs to possess all attributes required by an item in order to answer it correctly; \code{"OR"}: the examinee needs to possess only one of the attributes required by an item in order to answer it correctly.
}
  \item{max.ite}{The maximum number of iterations allowed.}
}

\value{
  \item{patterns}{All possible attribute profiles. Rows represent different patterns of attribute profiles and columns represent attributes. 1=examinee masters the attribute, 0=examinee does not master the attribute.}
  \item{initial.Q}{The initial Q-matrix. Rows represent items and columns represent attributes. 1=attribute required by the item, 0=attribute not required by the item. This is the preliminary Q-matrix to be refined.}
  \item{initial.class}{The row indices of \code{patterns} in the initial estimation of examinee attribute profiles.}
  \item{terminal.class}{The The row indices on \code{patterns} in the terminal estimation of examinee attribute profiles after the Q-matrix has been refined.}
  \item{modified.Q}{The modified Q-matrix. Rows represent items and columns represent attributes. 1=attribute required by the item, 0=attribute not required by the item.}
  \item{modified.entries}{The modified q-entries. Column 1 is the item ID of the modified entry; column 2 is the attribute ID of the modified entry.}
}

\seealso{
\code{\link{AlphaNP}}, \code{\link{plot.Qrefine}}
}

\references{
Chiu, C. Y. (2012). Statistical Refinement of the Q-matrix in Cognitive Diagnosis. \emph{Unpublished
manuscript}.
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

# Generate misspecified Q-matrix

Q_mis <- Q
Q_mis[1:3, 1] <- 1 - Q_mis[1:3, 1]

# Run Qrefine and create diagnostic plots

Qrefine.out <- Qrefine(data, Q_mis, gate="AND", max.ite=50)
plot(Qrefine.out)
}
